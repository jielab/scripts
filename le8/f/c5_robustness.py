"""FINAL C5 consensus/robustness layer.

Consumes the already-computed C5 annotation/CIGMA outputs. It never infers a
cell of origin from plasma abundance and never combines P values across evidence
families. The goal is to show which cell-context signals survive independent
background/null definitions and which upstream 5C evidence domains contribute.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd


def _read(path):
    path=Path(path)
    if not path.is_file(): return pd.DataFrame()
    try: return pd.read_csv(path)
    except Exception: return pd.DataFrame()


def _num(d,c):
    return pd.to_numeric(d[c],errors="coerce") if c in d else pd.Series(np.nan,index=d.index)


def final_consensus(stage, trait, layer):
    stage=Path(stage)
    if layer!="prot":
        pd.DataFrame([dict(status="not_applicable",detail="Metabolites have no unique gene/cell mapping; C5 consensus is protein-only")]).to_csv(stage/"c5.final_consensus_status.csv",index=False)
        return

    primary=_read(stage/"c5.cell.enrichment.csv")
    annotated=_read(stage/"c5.cell.enrichment_annotated_background.csv")
    matched=_read(stage/"c5.cell.matched_enrichment.csv")
    contrasts=_read(stage/"c5.cell.panel_contrasts.csv")
    pann=_read(stage/"c5.cell.panel_annotation.csv")
    cigma_cell=_read(stage/"c5.CIGMA_cell_types.csv")
    evidence=_read(stage/"c5.evidence_cell_context.csv")

    keys=["model","cell_type"]
    frames=[]
    if len(primary) and set(keys)<=set(primary):
        q="FDR_all_panel_cell_tests"
        a=primary.copy();a["primary_q"]=_num(a,q);a["primary_support"]=a.primary_q<.05
        frames.append(a[keys+[c for c in ["hits","odds_ratio","primary_q","primary_support"] if c in a]])
    if len(annotated) and set(keys)<=set(annotated):
        a=annotated.copy();a["annotated_q"]=_num(a,"FDR_all_panel_cell_tests");a["annotated_support"]=a.annotated_q<.05
        frames.append(a[keys+["annotated_q","annotated_support"]])
    if len(matched) and set(keys)<=set(matched):
        a=matched.copy();a["matched_q"]=_num(a,"FDR_matched_family");a["matched_support"]=a.matched_q<.05
        frames.append(a[keys+[c for c in ["expected_matched_hits","matched_q","matched_support"] if c in a]])

    if frames:
        z=frames[0]
        for a in frames[1:]: z=z.merge(a,on=keys,how="outer")
    else:
        z=pd.DataFrame(columns=keys)

    # Explicit CIGMA cell-specific support can be linked only through panel genes.
    if len(cigma_cell) and len(pann) and {"gene","cell_type"}<=set(cigma_cell) and {"model","gene"}<=set(pann):
        cc=cigma_cell.copy();cc["cigma_q"]=_num(cc,"FDR_cell_manifest")
        cc=cc[cc.cigma_q<.05]
        pg=pann[["model","gene"]].dropna().drop_duplicates()
        cg=pg.merge(cc[["gene","cell_type","cigma_q"]],on="gene",how="inner")
        if len(cg):
            cg=(cg.groupby(keys,as_index=False).agg(cigma_supported_genes=("gene","nunique"),cigma_best_q=("cigma_q","min")))
            cg["cigma_support"]=cg.cigma_supported_genes>0
            z=z.merge(cg,on=keys,how="outer")

    for c in ["primary_support","annotated_support","matched_support","cigma_support"]:
        if c not in z: z[c]=False
        z[c]=z[c].fillna(False).astype(bool)
    if len(z):
        z["independent_support_count"]=z[["primary_support","annotated_support","matched_support","cigma_support"]].sum(axis=1)
        z["robustness_tier"]=np.select([
            z.independent_support_count>=3,
            z.independent_support_count==2,
            z.independent_support_count==1],
            ["A: >=3 support routes","B: 2 support routes","C: 1 support route"],default="D: no FDR support")
        z["interpretation"]="Cell-expression context robustness; not cell of secretion, abundance, aging, mediation or causal cell of action"
        z=z.sort_values(["independent_support_count","model","cell_type"],ascending=[False,True,True])
    z.to_csv(stage/"c5.cell_robustness_consensus.csv",index=False)

    # Panel contrasts remain separate; summarize only FDR-significant direction, never merge with enrichment P values.
    if len(contrasts):
        qcol=next((c for c in ("FDR_all_contrasts","FDR_all_panel_cell_tests","FDR") if c in contrasts),None)
        if qcol:
            contrasts=contrasts.copy();contrasts["contrast_support"]=_num(contrasts,qcol)<.05
        contrasts.to_csv(stage/"c5.final_panel_contrasts.csv",index=False)

    # Upstream 5C evidence x cell context matrix. This is a count of supported source rows,
    # not a vote score and not a combined P value.
    if len(evidence) and {"evidence_domain","cell_type"}<=set(evidence):
        e=evidence.copy();p=_num(e,"p");q=_num(e,"q");val=_num(e,"value")
        support=np.zeros(len(e),dtype=bool)
        domain=e.evidence_domain.fillna("").astype(str)
        support |= domain.eq("C1_association") & (q<.05)
        support |= domain.eq("C2_MR") & (q<.05)
        # C3 value stores H4 in the current ledger.
        support |= domain.eq("C3_coloc") & (val>=.70)
        support |= domain.eq("C4_proxy") & (q<.05)
        support |= domain.eq("C4_path") & (q<.05)
        e["source_row_support"]=support
        m=(e[e.cell_type.notna() & e.source_row_support]
           .groupby(["cell_type","evidence_domain"],as_index=False)
           .agg(supported_rows=("source_row_support","size"),supported_features=("feature","nunique")))
        m.to_csv(stage/"c5.cell_by_5C_evidence.csv",index=False)

    status=pd.DataFrame([dict(status="completed",trait=trait,layer=layer,
        cell_model_rows=len(z),
        note="Independent cell-context routes summarized without P-value combination; CIGMA used only when explicit cell-specific tests exist")])
    status.to_csv(stage/"c5.final_consensus_status.csv",index=False)
