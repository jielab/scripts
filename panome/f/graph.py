"""Inductive person atlas. Held-out people never attend to one another."""
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import eigsh
from scipy.sparse.csgraph import connected_components
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score, adjusted_rand_score
import igraph
import leidenalg
from common import log

class PersonAtlas:
    def fit_index(self, train, a):
        if len(train) > 5000:
            from pynndescent import NNDescent
            index = NNDescent(train, n_neighbors=min(len(train)-1, max(40, a.neighbors+1)),
                               n_jobs=a.cores, random_state=a.seed)
            index.prepare()
        else:
            index = NearestNeighbors(n_jobs=a.cores).fit(train)
        return index

    def nearest(self, z, query_ids=None, query_groups=None):
        count = len(self.train_ids)
        extra = 1
        if query_groups is not None and self.train_groups is not None:
            extra += int(pd.Series(self.train_groups).value_counts().max())
        request = min(count, self.k + extra)
        if hasattr(self.index, "query"):
            ids, distance = self.index.query(z, k=request)
        else:
            distance, ids = self.index.kneighbors(z, n_neighbors=request)
        out_ids = np.empty((len(z), self.k), dtype=int)
        out_dist = np.empty((len(z), self.k), dtype=float)
        for i in range(len(z)):
            keep = np.ones(request, dtype=bool)
            if query_ids is not None:
                keep &= self.train_ids[ids[i]] != str(query_ids[i])
            if query_groups is not None and self.train_groups is not None:
                keep &= self.train_groups[ids[i]] != str(query_groups[i])
            candidates = np.flatnonzero(keep)[:self.k]
            if len(candidates) < self.k:
                raise ValueError("Insufficient independent neighbors after self/family exclusion")
            out_ids[i], out_dist[i] = ids[i, candidates], distance[i, candidates]
        return out_ids, out_dist

    def fit(self, z, ids, groups, a, out):
        self.scaler = StandardScaler().fit(z)
        train = self.scaler.transform(z).astype("float32")
        self.train_ids = np.asarray(ids, dtype=str)
        self.train_groups = None if groups is None else np.asarray(groups, dtype=str)
        n = len(z)
        largest_group = 1 if groups is None else int(pd.Series(groups).value_counts().max())
        self.k = min(a.neighbors, n - largest_group)
        if self.k < 2:
            raise ValueError("Insufficient people/groups for a molecular neighborhood")
        self.index = self.fit_index(train, a)
        neighbors, dist = self.nearest(train, self.train_ids, self.train_groups)
        weight = np.exp(-dist / np.maximum(dist[:, -1:], 1e-8))
        adjacency = sparse.csr_matrix((weight.ravel(), (
            np.repeat(np.arange(n), self.k), neighbors.ravel())), shape=(n, n))
        adjacency = adjacency.maximum(adjacency.T)
        adjacency.setdiag(0)
        adjacency.eliminate_zeros()
        sparse.save_npz(out / "sample_graph.npz", adjacency)
        coo = sparse.triu(adjacency, k=1).tocoo()
        graph = igraph.Graph(n=n, edges=list(zip(coo.row.tolist(), coo.col.tolist())))
        candidates, partitions = [], {}
        for resolution in a.resolutions:
            labels = np.asarray(leidenalg.find_partition(graph,
                leidenalg.RBConfigurationVertexPartition, weights=coo.data.tolist(),
                resolution_parameter=resolution, seed=a.seed).membership)
            counts = np.bincount(labels)
            valid = 1 < len(counts) <= a.max_states and counts.min() >= max(5, int(a.min_state_fraction*n))
            score = np.nan
            if valid:
                try:
                    score = float(silhouette_score(train, labels,
                        sample_size=min(3000, n), random_state=a.seed))
                except ValueError:
                    valid = False
            candidates.append(dict(resolution=resolution, n_states=len(counts),
                min_size=int(counts.min()), silhouette=score, eligible=bool(valid)))
            partitions[resolution] = labels
        pd.DataFrame(candidates).to_csv(out / "resolution_search.csv", index=False)
        good = [v for v in candidates if v["eligible"]]
        self.resolution = None
        if good:
            chosen = max(good, key=lambda v: v["silhouette"])
            self.resolution = chosen["resolution"]
            labels = partitions[self.resolution]
            order = np.argsort(-np.bincount(labels), kind="stable")
            mapping = {old: new for new, old in enumerate(order)}
            labels = np.array([mapping[v] for v in labels])
            state_status = "discovered"
        else:
            labels = np.zeros(n, dtype=int)
            state_status = "no_eligible_discrete_partition"
            log("NOTE", "s4_graph", "No eligible multistate partition; continuous atlas remains available.")
        self.labels = labels
        self.n_states = int(labels.max()+1)
        rng = np.random.default_rng(a.seed)
        stability = []
        if self.resolution is not None:
            for repeat in range(a.stability_repeats):
                keep = rng.random(len(coo.data)) >= .1
                pg = igraph.Graph(n=n, edges=list(zip(coo.row[keep].tolist(), coo.col[keep].tolist())))
                other = leidenalg.find_partition(pg, leidenalg.RBConfigurationVertexPartition,
                    weights=coo.data[keep].tolist(), resolution_parameter=self.resolution,
                    seed=a.seed+repeat+1)
                stability.append(dict(replicate=repeat+1, edge_dropout=.1,
                    ARI=float(adjusted_rand_score(labels, other.membership))))
        pd.DataFrame(stability, columns=["replicate", "edge_dropout", "ARI"]).to_csv(
            out / "graph_perturbation_stability.csv", index=False)
        n_components = connected_components(adjacency, directed=False, return_labels=False)
        self.psi = np.empty((n, 0), dtype="float32")
        eigenvalues = []
        requested = min(n-2, a.graph_dims + n_components)
        if requested < 256:
            degree = np.asarray(adjacency.sum(1)).ravel()
            sym = sparse.diags(1/np.sqrt(degree)) @ adjacency @ sparse.diags(1/np.sqrt(degree))
            vals, vectors = eigsh(sym, k=requested, which="LA",
                                  v0=rng.normal(size=n), tol=1e-5)
            order = np.argsort(vals)[::-1]
            # Every component has a trivial eigenvalue of one, not just the first.
            order = [j for j in order if 1e-6 < vals[j] < 1-1e-5][:a.graph_dims]
            self.psi = (vectors[:, order] / np.sqrt(degree[:, None])).astype("float32")
            eigenvalues = vals[order].tolist()
        self.novelty_threshold = float(np.quantile(dist[:, 0], .99))
        self.train_z = np.asarray(z, dtype="float32")
        self.summary = dict(n_states=self.n_states, state_status=state_status,
            resolution=self.resolution, n_components=int(n_components),
            n_edges=int(graph.ecount()), diffusion_dimensions=self.psi.shape[1],
            eigenvalues=eigenvalues, novelty_threshold=self.novelty_threshold,
            perturbation_ARI_mean=float(np.mean([v["ARI"] for v in stability])) if stability else None,
            limitation="Edge dropout is not full-pipeline resampling stability.")
        return self

    def project(self, z, ids=None, groups=None):
        scaled = self.scaler.transform(z).astype("float32")
        neighbors, dist = self.nearest(scaled, ids, groups)
        weights = np.exp(-dist / np.maximum(dist[:, -1:], 1e-8))
        weights /= weights.sum(1, keepdims=True)
        soft = np.stack([(weights * (self.labels[neighbors] == s)).sum(1)
                         for s in range(self.n_states)], 1)
        return dict(neighbors=neighbors, weights=weights.astype("float32"),
            soft=soft.astype("float32"), state=soft.argmax(1),
            neighbor_distance=dist[:, 0], novelty=dist[:, 0] > self.novelty_threshold,
            effective_neighbors=1/np.sum(weights**2, axis=1),
            neighborhood=np.einsum("nk,nkd->nd", weights, self.train_z[neighbors]).astype("float32"),
            # Same one-step extension for training and new people; self excluded.
            diffusion=np.einsum("nk,nkd->nd", weights, self.psi[neighbors]).astype("float32"))
