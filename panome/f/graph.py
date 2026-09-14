"""Sparse training graph, Leiden states, inductive state and diffusion projection."""
import joblib
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import eigsh
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.preprocessing import StandardScaler
import igraph as ig
import leidenalg


def neighbor_query(z,query,k,a):
    if len(z)>5000:
        from pynndescent import NNDescent
        index=NNDescent(z,n_neighbors=max(k+1,30),random_state=a.seed,n_jobs=a.cores)
        index.prepare()
        ids,dist=index.query(query,k=k)
    else:
        index=NearestNeighbors(n_neighbors=k,n_jobs=a.cores).fit(z)
        dist,ids=index.kneighbors(query)
    return index,ids,dist


def build_graph(z,part,a,out):
    tr=np.flatnonzero(part=='train');hold=np.flatnonzero(part!='train')
    scaler=StandardScaler().fit(z[tr]);zz=scaler.transform(z).astype('float32');train=zz[tr]
    n=len(train);k=min(a.neighbors,n-1)
    index,ids0,d0=neighbor_query(train,train,k+1,a)
    ids=np.empty((n,k),int);dist=np.empty((n,k),float)
    for i in range(n):
        ok=ids0[i]!=i
        ids[i]=ids0[i,ok][:k];dist[i]=d0[i,ok][:k]
    bandwidth=np.maximum(dist[:,-1],1e-8)
    w=np.exp(-dist/ bandwidth[:,None])
    adjacency=sparse.csr_matrix((w.ravel(),(np.repeat(np.arange(n),k),ids.ravel())),shape=(n,n))
    adjacency=adjacency.maximum(adjacency.T);adjacency.setdiag(0);adjacency.eliminate_zeros()
    sparse.save_npz(out/'sample_graph.npz',adjacency)
    coo=sparse.triu(adjacency,k=1).tocoo()
    graph=ig.Graph(n=n,edges=list(zip(coo.row.tolist(),coo.col.tolist())),directed=False)
    weights=coo.data.tolist()
    candidates=[];labels_by_resolution={}
    for resolution in a.resolutions:
        labels=np.asarray(leidenalg.find_partition(graph,leidenalg.RBConfigurationVertexPartition,
                          weights=weights,resolution_parameter=resolution,seed=a.seed).membership)
        counts=np.bincount(labels)
        valid=1<len(counts)<=a.max_states and counts.min()>=max(5,int(a.min_state_fraction*n))
        score=float(silhouette_score(train,labels,sample_size=min(3000,n),random_state=a.seed)) if valid else float('-inf')
        candidates.append(dict(resolution=resolution,n_states=len(counts),min_size=int(counts.min()),silhouette=score,eligible=valid))
        labels_by_resolution[resolution]=labels
    pd.DataFrame(candidates).to_csv(out/'resolution_search.csv',index=False)
    viable=[c for c in candidates if c['eligible']]
    if not viable: raise ValueError('No stable-sized multi-state partition; inspect resolution_search.csv and specify a broader resolution grid. Do not force five states.')
    selected=max(viable,key=lambda c:c['silhouette']);labels=labels_by_resolution[selected['resolution']]
    order=np.argsort(-np.bincount(labels));mapping={int(old):new for new,old in enumerate(order)}
    labels=np.array([mapping[v] for v in labels]);n_states=len(order)
    stability=[]
    # Seed stability is explicitly not resampling stability; also perturb graph edges.
    rng=np.random.default_rng(a.seed)
    for b in range(a.stability_repeats):
        keep=rng.random(len(weights))>.1
        perturbed=ig.Graph(n=n,edges=list(zip(coo.row[keep].tolist(),coo.col[keep].tolist())),directed=False)
        alt=leidenalg.find_partition(perturbed,leidenalg.RBConfigurationVertexPartition,
            weights=coo.data[keep].tolist(),resolution_parameter=selected['resolution'],seed=a.seed+b+1)
        stability.append(dict(replicate=b+1,edge_dropout=.1,ARI=adjusted_rand_score(labels,alt.membership)))
    pd.DataFrame(stability).to_csv(out/'graph_perturbation_stability.csv',index=False)
    if hasattr(index,'query'): hi,hd=index.query(zz[hold],k=k)
    else: hd,hi=index.kneighbors(zz[hold],n_neighbors=k)
    all_ids=np.empty((len(z),k),int);all_dist=np.empty((len(z),k))
    all_ids[tr]=ids;all_dist[tr]=dist;all_ids[hold]=hi;all_dist[hold]=hd
    ww=np.exp(-all_dist/np.maximum(all_dist[:,-1:],1e-8));ww/=ww.sum(axis=1,keepdims=True)
    soft=np.zeros((len(z),n_states))
    for state in range(n_states): soft[:,state]=(ww*(labels[all_ids]==state)).sum(axis=1)
    assigned=soft.argmax(axis=1);assigned[tr]=labels
    novelty_threshold=float(np.quantile(dist[:,0],.99))
    degree=np.asarray(adjacency.sum(axis=1)).ravel()
    sym=sparse.diags(1/np.sqrt(degree))@adjacency@sparse.diags(1/np.sqrt(degree))
    vals,vecs=eigsh(sym,k=min(a.graph_dims+1,n-2),which='LA',v0=rng.normal(size=n),tol=1e-4)
    order=np.argsort(vals)[::-1][1:];order=order[vals[order]>1e-5]
    vals=vals[order];psi=vecs[:,order]/np.sqrt(degree[:,None])
    diffusion=np.zeros((len(z),len(vals)),dtype='float32');diffusion[tr]=psi*vals
    # Nyström random-walk extension at diffusion time 1; no held-out graph refit.
    diffusion[hold]=(ww[hold,:,None]*psi[hi]).sum(axis=1)
    np.savez_compressed(out/'graph_coordinates.npz',soft=soft,diffusion=diffusion,state=assigned,
                        novelty=all_dist[:,0]>novelty_threshold,neighbor_distance=all_dist[:,0])
    joblib.dump(dict(index=index,scaler=scaler,labels=labels,train_rows=tr,psi=psi,eigenvalues=vals,
                     novelty_threshold=novelty_threshold,neighbors=k),out/'graph_model.joblib')
    summary=dict(**selected,n_edges=graph.ecount(),n_components=len(graph.connected_components()),
                 graph_dims=len(vals),perturbation_ARI_mean=float(np.mean([s['ARI'] for s in stability])),
                 novelty_threshold=novelty_threshold)
    return summary
