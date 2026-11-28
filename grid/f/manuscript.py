"""Publish aggregate-only tables and figures for the requested GRID traits."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

from core2 import atomic_table
from manuscript_stats import CORE, PRIMARY, SELECT
from output_layout import project_path, trait_path, module_results, upstream_root
from paper_compare import METHODS, publish_comparison

TABLES = ("overall", "paired", "screening", "support_counts", "ancestry", "anchor_associations",
          "binary_discrimination", "literature", "model_settings", "discovery", "cohort", "anchor_training",
          "csx_weights", "score_coverage", "individual_error_models", "assessment", "regional_gain", "ld_maps", "twin_screening", "blocknet_ablation","blocknet_twins","ld_boundaries")
FIGURES = ("fig1_ancestry", "fig2_abm_comparisons", "fig3_screening")
LABELS = {"abm":"ABM", PRIMARY:"PRS-CSx", "disco":"DiscoDivas",
          "cs.EUR":"PRS-CS EUR", "csx.input":"PRS-CSx fixed input"}
COLORS = {"abm":"#00765E", PRIMARY:"#175E9C", "disco":"#B65626","cs.EUR":"#715187"}


def tsv(path):
    return pd.read_csv(path, sep="\t")


def js(path):
    return json.loads(path.read_text()) if path.is_file() else {}


def context_files(root, trait):
    return [project_path(root, "tmp/cohort/cohort.json"), project_path(root, "qc/cohort.tsv"),
            trait_path(root, trait, "inference/provenance.json"), trait_path(root, trait, "benchmark.json"),
            trait_path(root, trait, "disco/model.json"), trait_path(root, trait, "abm/model.json"),
            project_path(root, f"qc/{trait}.coefficients.tsv"), project_path(root, f"qc/{trait}.disco_anchors.tsv"),
            trait_path(root, trait, "inference/weights/snp_counts.tsv"), module_results(root,"csx")/trait/"score_counts.tsv",
            trait_path(root,trait,"regions/metadata.json"),
            trait_path(root,trait,"regions/region_counts.tsv")]


def extract_context(root, traits, tables):
    cohort = js(project_path(root, "tmp/cohort/cohort.json"))
    counts = project_path(root, "qc/cohort.tsv")
    if counts.is_file():
        tables["cohort"].append(tsv(counts)[["role","ancestry","n","families"]])
    for t in traits:
        _, _, provfile, benchfile, discofile, abmfile, coeff, anchors, counts, score_path, regional_meta, regional_counts = context_files(root,t)
        prov, bench, disco, abm = map(js, (provfile, benchfile, discofile, abmfile))
        cols = ("geometry","projection_dimensions","radius","k","alpha","min_count","min_group_ess",
                "donor_count","radius_k","radius_quantile","global_residual_alpha","support_threshold_role")
        meta = js(regional_meta)
        files = meta.get("block_files",{})
        if not files and "consensus" in meta:
            files = {"shared": meta["consensus"].get("files",[])}
        observed_layout = meta.get("layout","consensus" if "consensus" in meta else "unrecorded")
        row = {k:abm.get(k) for k in cols}
        row.update(ld_layout=observed_layout, representation=meta.get("representation"),
                   geometry_phi=meta.get("geometry_phi"), captured_training_variance=abm.get("captured_training_variance"))
        if regional_counts.is_file():
            counts_by_map = tsv(regional_counts)
            for channel, d in counts_by_map.groupby("map"):
                records = files.get(channel,files.get("shared",[]))
                tables["ld_maps"].append(pd.DataFrame([dict(trait=t, channel=channel,
                    layout=observed_layout, block_build=meta.get("block_build","GRCh38" if "consensus" in meta else None),
                    assignment_position=meta.get("assignment_position"), blocks=int(d.n_catalog_blocks.sum()), empty=int(d.n_empty_blocks.sum()),
                    weight_snps=int(d.n_weight_snps.sum()), assigned=int(d.n_assigned.sum()),
                    unassigned=int(d.n_unassigned.sum()), coverage=d.n_assigned.sum()/d.n_weight_snps.sum(),
                    bed_files="|".join(r["path"] for r in records),
                    sha256="|".join(r.get("sha256","unrecorded") for r in records))]))
        row.update(trait=t, disco_status=disco.get("status","unrecorded"), disco_reason=disco.get("reason"),
                   active_anchors=",".join((disco.get("full") or {}).get("active_anchors",[])),
                   disco_input_phi=disco.get("input_phi"), csx_selection=bench.get("selected",{}).get(PRIMARY),
                   cs_best=bench.get("selected",{}).get("cs.best"),
                   support20_ess_cutoff=abm.get("policies",{}).get("support",{}).get("0p20"),
                   phenotype=json.dumps(cohort.get("phenotype_definitions",{}).get(t,{}), ensure_ascii=False),
                   covariates=",".join(bench.get("covariates",[])),
                   validation_status=cohort.get("validation_status","unrecorded"),
                   reference_build=prov.get("reference_build"),target_build=prov.get("target_build"))
        details = js(project_path(root,f"grid.paper/{t}.evaluation.json"))
        row["optional_methods_included"] = ",".join(details.get("optional_methods_included",[]))
        tables["model_settings"].append(pd.DataFrame([row]))
        for pop, rec in prov.get("raw_gwas", {}).items():
            fields = {k:rec.get(k) for k in ("n_gwas","n_gwas_original","n_source","sample_overlap","source_build")}
            tables["discovery"].append(pd.DataFrame([dict(trait=t, population=pop, **fields)]))
        if coeff.is_file():
            d = tsv(coeff)
            d["method"]=d.method.str.replace(r"^prscs\.", "cs.", regex=True)
            d = d.loc[d.method.isin([PRIMARY,"csx.input","cs.best"])]
            tables["csx_weights"].append(d.assign(trait=t))
        if anchors.is_file():
            tables["anchor_training"].append(tsv(anchors).assign(trait=t))
        if counts.is_file():
            d = tsv(counts).groupby("score_name").agg(posterior_variants=("n_variants","sum"),
                                                     posterior_chromosomes=("chromosome","nunique")).reset_index()
        else:
            d = pd.DataFrame(columns=["score_name"])
        if score_path.is_file():
            c = tsv(score_path)
            if {"score_name","n_scored_variants","chromosome"} <= set(c):
                c = c.groupby("score_name").agg(scored_variants=("n_scored_variants","sum"),
                                                scored_chromosomes=("chromosome","nunique")).reset_index()
                d = d.merge(c,on="score_name",how="outer",validate="one_to_one")
        if len(d):
            tables["score_coverage"].append(d.assign(trait=t))
        tables["individual_error_models"].append(pd.DataFrame([dict(trait=t, **details.get("error_head_audit",{}))]))


def ci_label(point, low, high, prefer_low=True):
    if not np.isfinite(point): return "unavailable"
    if not (np.isfinite(low) and np.isfinite(high)): return "interval unavailable"
    if low == high == point == 0: return "identical on this metric"
    if high < 0: return "favors first" if prefer_low else "favors reference"
    if low > 0: return "favors reference" if prefer_low else "favors first"
    return "inconclusive"


def assessments(tables):
    d = tables["paired"]
    out = []
    primary = d[d.ancestry.eq("ALL") & d.comparison.eq("method_minus_reference") &
                d.method.eq("abm") & d.reference.isin([PRIMARY,"disco"]) & d.metric.isin(["MSE","Brier"])]
    for r in primary.itertuples():
        out.append(dict(trait=r.trait,question="ABM prediction",selection=r.selection,method=r.method,
                        reference=r.reference,metric=r.metric,difference=r.difference,ci_low=r.ci_low,ci_high=r.ci_high,
                        interpretation=ci_label(r.difference,r.ci_low,r.ci_high)))
    for r in tables["screening"].query("selection == @SELECT").itertuples():
        for key, val, low, high in [
                ("selected minus all",r.selected_minus_all_loss,r.selection_low,r.selection_high),
                ("selected minus rejected, overlap standardized",r.risk_standardized_difference,
                 r.risk_difference_low,r.risk_difference_high)]:
            out.append(dict(trait=r.trait,question="ABM screening",method=r.method,reference=key,selection=r.selection,
                            metric=r.loss_name,difference=val,ci_low=low,ci_high=high,
                            interpretation=ci_label(val,low,high)))
    for r in d[d.ancestry.eq("ALL") & d.comparison.eq("selected_minus_all") & d.metric.eq("partial_R2")].itertuples():
        out.append(dict(trait=r.trait,question="ABM screening",method=r.method,reference="selected minus all",
                        selection=r.selection,metric=r.metric,difference=r.difference,ci_low=r.ci_low,ci_high=r.ci_high,
                        interpretation=ci_label(r.difference,r.ci_low,r.ci_high,False)))
    binary=tables["binary_discrimination"]
    if "ancestry" in binary:
        overall=binary[binary.ancestry.eq("ALL")]
        for r in overall[overall.method.eq("abm") & overall.reference.isin([PRIMARY,"disco"])].itertuples():
            for metric in ("auroc","average_precision"):
                point,lo,hi=(getattr(r,f"delta_{metric}{suffix}") for suffix in ("","_low","_high"))
                out.append(dict(trait=r.trait,question="ABM prediction",method=r.method,reference=r.reference,
                                selection=r.selection,metric=metric,difference=point,ci_low=lo,ci_high=hi,
                                interpretation=ci_label(point,lo,hi,False)))
        screened=overall[overall.selection.eq(SELECT) & overall.method.isin(CORE)].drop_duplicates("method")
        for r in screened.itertuples():
            for metric in ("auroc","average_precision"):
                point=getattr(r,f"selected_minus_all_{metric}",np.nan)
                lo,hi=(getattr(r,f"selection_{metric}_{suffix}",np.nan) for suffix in ("low","high"))
                out.append(dict(trait=r.trait,question="ABM screening",method=r.method,reference="selected minus all",
                                selection=r.selection,metric=metric,difference=point,ci_low=lo,ci_high=hi,
                                interpretation=ci_label(point,lo,hi,False)))
    return pd.DataFrame(out)


def literature_table(tables):
    """Pair the local estimate with its named reference and document comparability."""
    rows = []
    for r in tables["paired"].itertuples():
        metric = "Brier" if r.trait=="t2dm" else "partial_R2"
        if (r.metric!=metric or r.selection!="all" or r.comparison!="method_minus_reference"
                or r.method not in (PRIMARY,"disco")):
            continue
        comparable = (r.method==PRIMARY and r.reference.startswith("cs.")) or (
                      r.method=="disco" and r.reference in (PRIMARY,"csx.input"))
        if not comparable: continue
        rows.append(dict(trait=r.trait, ancestry=r.ancestry, method=r.method, reference=r.reference,
                         metric=metric, n=r.n, difference=r.difference, ci_low=r.ci_low,ci_high=r.ci_high,
                         interpretation=ci_label(r.difference,r.ci_low,r.ci_high,metric=="Brier"),
                         context=("Multi- vs single-discovery or independent stacking; local GWAS and target calibration"
                                  if r.method==PRIMARY else "Interpolated vs pooled CSx; different comparator from AJHG Fig. 3/4")))
    anchors=tables["anchor_associations"]
    for (trait,pop), group in anchors.groupby(["trait","ancestry"]):
        if pop=="ALL" or "disco" not in set(group.score):continue
        metric="log_or_per_sd" if trait=="t2dm" else "partial_r2"
        value=group.set_index("score")[metric]
        for ref in dict.fromkeys(["component_"+pop,"component_EUR"]):
            v=value.get(ref,np.nan)
            delta=value["disco"]-v
            desc="matched fine-tuning ancestry" if ref=="component_"+pop else "EUR fine-tuning ancestry"
            rows.append(dict(trait=trait, ancestry=pop, method="disco.genetic_score",reference=ref,
                n=group.loc[group.score.eq("disco"),"n"].iloc[0], metric=metric,method_estimate=value["disco"],
                reference_estimate=v,difference=delta, relative_percent=100*delta/v if metric=="partial_r2" and v>0 else np.nan,
                interpretation="anchor unavailable" if ref not in value else "descriptive point estimates",
                context=desc+"; same local samples and covariates; no paired CI on genetic association difference"))
    return pd.DataFrame(rows)


def draw_figures(tables, traits, out):
    plt.rcParams.update({"font.size":9,"axes.spines.right":False,"axes.spines.top":False,
                         "savefig.dpi":300})
    d = tables["ancestry"]
    for panel, stem in enumerate(FIGURES):
        fig, axes = plt.subplots(1,len(traits),figsize=(5.1*len(traits),4.3), squeeze=False)
        for ax, trait in zip(axes[0],traits):
            binary = trait == "t2dm"
            metric = "AUROC" if binary else "Partial R²"
            ax.set_title(trait.upper() if trait!="height" else "Height",loc="left",fontweight="bold")
            if panel == 0:
                frame = d[d.trait.eq(trait) & d.selection.eq("all")]
                pops = sorted(frame.ancestry.unique())
                for j,method in enumerate(("cs.EUR",PRIMARY,"disco","abm")):
                    m = frame[frame.method.eq(method)].set_index("ancestry")
                    key = "auroc" if binary else "partial_r2_vs_covariate"
                    if m.empty or key not in m: continue
                    ax.scatter(np.arange(len(pops)) + (j-1.5)*.14, m[key].reindex(pops),
                               label=LABELS[method],color=COLORS[method],s=23)
                labels = pops
                if binary and "n_cases" in frame:
                    cases = frame.drop_duplicates("ancestry").set_index("ancestry").n_cases
                    labels = [f"{pop}\n{int(cases[pop])} cases" if pd.notna(cases[pop]) else pop for pop in pops]
                ax.set_xticks(range(len(pops)),labels,rotation=25)
                ax.set_ylabel(metric)
            elif panel == 1:
                frame = tables["paired"]
                frame = frame[frame.trait.eq(trait) & frame.ancestry.eq("ALL") & frame.method.eq("abm") &
                              frame.reference.isin([PRIMARY,"disco"]) & frame.metric.eq("Brier" if binary else "MSE") &
                              frame.comparison.eq("method_minus_reference")]
                ax.axvline(0,color="#89939A",lw=1)
                labels=[]
                for i,r in enumerate(frame.itertuples()):
                    labels.append(("All" if r.selection=="all" else "High support")+" / "+LABELS[r.reference])
                    if np.isfinite(r.difference):
                        ax.plot([r.ci_low,r.ci_high],[i,i],color=COLORS[r.reference],lw=2)
                        ax.scatter(r.difference,i,color=COLORS[r.reference],zorder=3)
                ax.set_yticks(range(len(labels)),labels); ax.invert_yaxis()
                ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
                ax.ticklabel_format(axis="x",style="sci",scilimits=(-2,3),useMathText=True)
                ax.set_xlabel("ABM minus comparator "+("Brier score" if binary else "MSE"))
            else:
                frame = tables["screening"]
                frame = frame[frame.trait.eq(trait)]
                ax.axhline(0,color="#89939A",lw=1)
                for method in CORE:
                    v = frame[frame.method.eq(method)].sort_values("actual_coverage")
                    if v.empty: continue
                    x = v.actual_coverage.to_numpy(float)
                    ax.plot(x, v.selected_minus_all_loss,"o-",ms=4,label=LABELS[method],color=COLORS[method])
                    ax.fill_between(x,v.selection_low.to_numpy(float),v.selection_high.to_numpy(float),color=COLORS[method],alpha=.10)
                ax.set_xlabel("Fraction of test participants selected"); ax.set_xlim(-.02,1.02)
                ax.set_ylabel("Selected minus full "+("Brier score" if binary else "MSE"))
            if panel!=1 and ax.get_legend_handles_labels()[0]:
                ax.legend(fontsize=8,frameon=False)
        fig.tight_layout()
        fig.savefig(out/f"{stem}.png",bbox_inches="tight")
        plt.close(fig)


PAPER_NOTE = """<p>PRS-CSx (Ruan et al., Nature Genetics 2022, DOI 10.1038/s41588-022-01054-7)
showed gains against single-discovery methods across quantitative traits in several non-European target settings.
Compare csx.population with cs.EUR, the matching single-ancestry PRS-CS, and independently fitted
cs.linear_combination within each ancestry. cs.best is chosen in TUNE. csx.input holds the four
Disco input scores fixed. Overall GRID performance includes age, sex and PCs. The frozen covariate-only predictor supplies the partial R² denominator.</p>
<p>DiscoDivas (Ruan et al., AJHG 2026, DOI 10.1016/j.ajhg.2026.05.006) used LDpred2 inputs.
Its Fig. 3 summarizes seven quantitative traits, including LDL, with clearer gains in YGBB OTH.
Height is an extension. Their LDL was medication-adjusted. Fig. 4 reports adjusted log OR per SD for
diabetes, using AoU tuning and YGBB/MGBB testing. GRID uses CSx inputs, and its default T2DM endpoint is baseline prevalence.
The anchor association table evaluates the interpolated genetic score and each available component on identical
complete cases with the recorded covariates. Binary log ORs use one within-stratum score SD and family-sandwich intervals.
Continuous partial R² uses residuals from covariate adjustment. These association estimates use TEST outcomes only for evaluation.
A local OTH-matched anchor is unavailable under the default four-anchor specification. Compare OTH with the named
fixed anchors or pooled CSx and retain this difference in reporting. OTH assignment alone does not establish individual admixture.</p>"""


def publish(cfg, root):
    root=Path(root); traits=cfg["traits"]; out=project_path(root,"grid.paper")
    tables = {k:[] for k in TABLES}
    for t in traits:
        p=tsv(out/t/f"{t}.performance.tsv")
        tables["overall"].append(p[p.selection.eq("all")])
        tables["screening"].append(p[p.method.isin(METHODS) & p.selection.str.fullmatch(r"selected_support_\d+p\d+")])
        tables["twin_screening"].append(p[p.method.isin(METHODS) & p.selection.isin([SELECT,"selected_twin_count_0p20"])])
        for key,suffix in (("paired","contrasts"),("ancestry","ancestry"),("support_counts","support"),
                           ("anchor_associations","anchor_associations")):
            tables[key].append(tsv(out/t/f"{t}.{suffix}.tsv"))
        if t=="t2dm":
            tables["binary_discrimination"].append(tsv(out/t/f"{t}.paired_discrimination.tsv"))
    if "ldl" in traits:
        tables["ld_boundaries"].append(tsv(project_path(root,"grid.paper/ldl.boundaries.tsv")))
        tables["blocknet_twins"].append(tsv(project_path(root,"grid.paper/ldl.blocknet_twins.tsv")))
    extract_context(root,traits,tables)
    tables = {k:pd.concat(frames,ignore_index=True) if frames else pd.DataFrame({"status":["unavailable or not applicable"]})
              for k,frames in tables.items()}
    paired = tables["paired"]
    from ldl_report import ABLATIONS
    ablation=paired.loc[paired.trait.eq("ldl") & paired.selection.eq("all") &
                        pd.MultiIndex.from_frame(paired[["method","reference"]]).isin(ABLATIONS)].copy()
    tables["blocknet_ablation"]=ablation
    gain = paired[paired.ancestry.eq("ALL") & paired.selection.eq("all") &
                  paired.reference.isin([PRIMARY,"direct.linear"]) &
                  paired.method.isin(["abm","block.glmnet","block.lightgbm","abm.deep","direct.blocknet"]) &
                  paired.metric.isin(["MSE","Brier","relative_loss_reduction_pct","predictive_R2","partial_R2"])].copy()
    tables["regional_gain"] = gain
    tables["assessment"] = assessments(tables)
    tables["literature"] = literature_table(tables)
    for name,frame in tables.items():
        atomic_table(frame,out/"tables"/f"table_{name}.tsv")
    (out/"figures").mkdir(parents=True,exist_ok=True)
    draw_figures(tables,traits,out/"figures")
    publish_comparison(out,tables,traits)
    if "ldl" in traits:
        from ldl_report import plot
        plot(tables,out/"figures/fig4_ldl_blocknet.png")
    index=out/"manuscript.html"
    links=" · ".join(f'<a href="tables/table_{name}.tsv">{name}</a>' for name in tables)
    cols=["trait","question","selection","method","reference","metric","difference","ci_low","ci_high","interpretation"]
    view=tables["assessment"][cols].to_html(index=False,na_rep="—",float_format=lambda x:f"{x:.5g}")
    literature=tables["literature"]
    lit= literature[literature.ancestry.isin(["OTH","EAS","AFR","SAS"]) & (
        (literature.method.eq(PRIMARY) & literature.reference.eq("cs.EUR")) |
        (literature.method.eq("disco") & literature.reference.eq("csx.input")))]
    litview=lit[["trait","ancestry","n","method","reference","metric","difference","ci_low","ci_high","interpretation"]].to_html(
        index=False,na_rep="—",float_format=lambda v:f"{v:.4g}")
    overview = pd.concat([tables["overall"], tables["screening"]],ignore_index=True)
    overview = overview.loc[overview.method.isin(METHODS) & overview.selection.isin(["all",SELECT]),
                            ["trait","method","selection","n_evaluable","actual_coverage","primary_loss"] +
                            [k for k in ("partial_r2_vs_covariate","auroc","auprc","n_cases","case_retention") if k in overview]]
    overview = overview.sort_values(["trait","method","selection"])
    overview_html = overview.to_html(index=False,na_rep="—",float_format=lambda v:f"{v:.5g}")
    models=tables["model_settings"][["trait","ld_layout","geometry","projection_dimensions","donor_count","radius","k","alpha","support20_ess_cutoff","disco_status","active_anchors","validation_status"]]
    figures='<figure><img src="figures/fig5_method_comparison.png" alt="Paired loss differences from PRS-CSx"><figcaption>Paired differences in held-out loss with family-bootstrap intervals.</figcaption></figure>'+"".join(f'<figure><img src="figures/{s}.png" alt="{caption}"><figcaption>{caption}.</figcaption></figure>'
                    for s,caption in zip(FIGURES,["Predictive performance by ancestry","Paired differences on identical test participants",
                                                   "Test error as a function of calibration-frozen support coverage"]))
    if "ldl" in traits:
        figures += '<figure><img src="figures/fig4_ldl_blocknet.png" alt="LDL predictive R squared"><figcaption>LDL ablation controls and prespecified precision screening on the same TEST roster.</figcaption></figure>'
    ldl_view=""
    if "ldl" in traits:
        from html import escape
        fields=["ancestry","policy","n_evaluable","coverage","mse","predictive_r2","selected_minus_all_mse",
                "ci_low","ci_high","selected_minus_risk_only_mse","risk_only_ci_low","risk_only_ci_high",
                "risk_standardized_difference","blocknet_minus_linear_r2"]
        twins=tables["blocknet_twins"]
        if "policy" in twins:
            frame=twins.loc[twins.ancestry=="ALL"].reindex(columns=fields)
            ldl_view="<h2>LDL twin screening</h2>"+frame.to_html(index=False,na_rep="—",float_format=lambda x:f"{x:.5g}")
        boundary=tables["ld_boundaries"]
        if "metric" in boundary:
            frame=boundary.loc[boundary.metric=="predictive_R2",["ancestry","n","difference","ci_low","ci_high"]]
            ldl_view+="<h2>BlockNet ΔR² · uniform minus ancestry BEDs</h2>"+frame.to_html(index=False,float_format=lambda x:f"{x:.5g}")
        else:
            ldl_view+="<p>Boundary comparison: "+escape(str(boundary.status.iloc[0]))+"</p>"
    index.write_text(f"""<!doctype html><html lang="en"><meta charset="utf-8"><title>GRID manuscript results</title>
<style>body{{font:15px system-ui;max-width:1350px;margin:36px auto;padding:0 22px;color:#19323b}}p{{max-width:1000px;line-height:1.55}}
td,th{{padding:5px 8px;text-align:left;border-bottom:1px solid #ddd}}table{{border-collapse:collapse;font-size:13px}}
figure{{margin:26px 0}}img{{max-width:100%}}a{{color:#165783}}</style><h1>GRID manuscript results</h1>
<p><a href="comparison.html">Explore performance by trait, ancestry and support</a> · {links}</p><h2>Full test cohort and high-support participants</h2>{overview_html}<h2>Paired differences</h2>
<p>Prediction differences use both methods on the same people. Negative MSE and Brier differences favor the first method.
The high-support threshold ranks family ESS in CALIBRATION after the support eligibility rule. Actual TEST coverage
and case retention appear in the screening table. Selected-minus-all compares nested groups. The selected-minus-rejected
contrast standardizes on overlapping ancestry × frozen CSx risk cells. Inspect the overlap fractions and valid draws.
CIs are pointwise family bootstrap intervals conditional on the fitted models and thresholds.
T2DM AUROC and average precision contrasts appear in binary_discrimination. Compare partial R² and
discrimination alongside loss, since groups can differ in outcome variance and prevalence.</p>{view}
<h2>Frozen model settings</h2>{models.to_html(index=False,na_rep="—")}
<p>n_support counts training donors within the frozen radius after excluding the query family. k limits residual
borrowing. Family ESS measures the concentration of donor families. Alpha of zero returns the CSx prediction.
n_donors_used counts donors admitted to the chosen residual correction. A zero alpha gives zero donors used even when the radius contains many matches. The support_counts table includes selected and rejected TEST groups, the distribution of matched donors, and the fraction whose ABM prediction changed from CSx.</p>
<h2>Regional gains on the full test cohort</h2>
<p>The LDL ablation compares global ridge, global MLP, regional ridge and the full BlockNet. Each uses identical families, outcomes, side features and calibration rules. The blocknet_twins table reports coverage, paired loss differences, and the incremental accuracy of predicted errors. Thresholds are frozen before TEST. Risk-standardized differences compare selected and rejected people in overlapping ancestry and BlockNet risk strata.</p>
<p>The regional_gain table gives each model's paired loss and R² differences from CSx. Positive relative_loss_reduction_pct means lower loss, calculated as 100 × (reference loss − candidate loss) / reference loss. The ld_maps table identifies each BED and reports SNP coverage.</p>
{figures}{ldl_view}<h2>Comparison with published work</h2>{PAPER_NOTE}{litview}
<p>Consult discovery for recorded discovery–target overlap. An unknown entry leaves external generalizability unresolved.
Intervals for multiple traits, ancestries and coverages are descriptive. Independent validation is needed after an inspected TEST split.</p>
</html>""",encoding="utf-8")
    from html import escape
    from urllib.parse import quote
    from file_lock import file_lock
    from tempfile import NamedTemporaryFile
    main=upstream_root(root)/"paper"/"main"
    with file_lock(main.parent/"logs"/"index.lock"):
        links=(['<li><a href="manuscript.html">Main report</a></li>'] if (main/"manuscript.html").is_file() else []) + [f'<li><a href="{quote(p.parent.name)}/manuscript.html">{escape(p.parent.name)}</a></li>'
               for p in sorted(main.glob("*/manuscript.html"))]
        with NamedTemporaryFile("w",encoding="utf-8",dir=main,delete=False) as f:
            f.write('<!doctype html><meta charset="utf-8"><title>GRID results</title><h1>GRID results</h1><ul>'+"".join(links)+"</ul>")
        Path(f.name).replace(main/"experiments.html")
    print(f"Manuscript results: {index}",flush=True)
