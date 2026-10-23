"""Exploratory, matched molecular comparisons and a local evolution report."""
from __future__ import annotations

from collections import Counter, defaultdict
import html
import json
import math
from pathlib import Path

from ancient_join import number, table_rows
from core import dump_json, fingerprint, tsv_writer


def _bh(values):
    """Benjamini-Hochberg within the explicitly supplied family."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    result, bound = [None] * len(values), 1.
    for rank in range(len(order), 0, -1):
        i = order[rank - 1]
        bound = min(bound, values[i] * len(order) / rank)
        result[i] = bound
    return result


def compare_geometry(geometry_path, out_dir, *, draws=10000, seed=20261021, min_pairs=10):
    """Within-scorer normalized distance contrast, paired sign flips and bootstrap.

    This is observational, conditional on candidate selection. Physical block
    thinning does not guarantee LD independence; uncertainty is exploratory.
    Never regard weak selection evidence as proof of neutrality.
    """
    import numpy as np
    if draws < 100 or min_pairs < 2:
        raise ValueError("At least 100 resampling draws and min_pairs >= 2 are required")
    out = Path(out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    grouped, qc = defaultdict(dict), Counter()
    fields = ["site_id", "group", "match_set", "scorer", "status", "closeness_to_increasing",
              "n_features", "axis_sha256", "locus_block"]
    handle, writer = tsv_writer(out / "site_closeness.tsv", fields)
    with handle:
        for row in table_rows(geometry_path, ("site_id", "group", "match_set", "atlas_scorer", "atlas_status")):
            qc[row["atlas_status"]] += 1
            up = number(row.get("atlas_distance_to_estimated_increasing"))
            down = number(row.get("atlas_distance_to_estimated_decreasing"))
            valid = row["atlas_status"] == "ok" and up is not None and down is not None and up + down > 0
            score = (down - up) / (down + up) if valid else None
            writer.writerow(dict(site_id=row["site_id"], group=row["group"], match_set=row["match_set"],
                                 scorer=row["atlas_scorer"], status=row["atlas_status"],
                                 closeness_to_increasing=score if score is not None else ".",
                                 n_features=row.get("atlas_n_features", "."), axis_sha256=row.get("atlas_axis_sha256", "."),
                                 locus_block=row.get("locus_block", ".")))
            if not valid or row["match_set"] in ("", ".") or row["group"] not in ("strong", "weak_signal"):
                continue
            key = (row["atlas_scorer"], row["match_set"])
            bucket=grouped[key].setdefault(row['group'],[])
            if any(x[2]==row['site_id'] for x in bucket) or (row['group']=='strong' and bucket):
                raise ValueError(f"Repeated case/site in match set/scorer: {key}")
            bucket.append((score,row.get('atlas_n_features'),row['site_id'],row.get('atlas_axis_sha256'),row.get('locus_block','.')))
    by_scorer = defaultdict(list)
    handle, writer = tsv_writer(out / "paired_closeness.tsv", ["scorer", "axis_sha256", "match_set", "strong_site", "weak_site", "strong_closeness", "weak_closeness", "difference", "status", "n_controls"])
    with handle:
        for (scorer, pair), groups in sorted(grouped.items()):
            if set(groups) != {"strong", "weak_signal"}:
                qc["incomplete_geometry_pairs"] += 1
                continue
            case, controls = groups['strong'][0], groups['weak_signal']
            status=('feature_count_mismatch' if any(case[1]!=c[1] for c in controls) else
                    'feature_axis_mismatch' if not case[3] or case[3] in ('.','NA') or any(case[3]!=c[3] for c in controls) else 'ok')
            control=(sum(c[0] for c in controls)/len(controls),case[1],','.join(c[2] for c in controls),case[3])
            difference = case[0] - control[0]
            writer.writerow(dict(scorer=scorer, axis_sha256=case[3] or ".", match_set=pair, strong_site=case[2], weak_site=control[2],
                                 strong_closeness=case[0], weak_closeness=control[0], difference=difference, status=status,n_controls=len(controls)))
            if status == "ok":
                by_scorer[(scorer, case[3])].append((difference,case[4] if case[4] not in (".","NA","") else pair))
            else:
                qc[status] += 1
    results, tested = [], []
    rng = np.random.default_rng(seed)
    for (scorer, axis), values in sorted(by_scorer.items()):
        arr = np.asarray([x[0] for x in values], dtype=float)
        clusters=defaultdict(list)
        for value,block in values: clusters[block].append(value)
        groups=list(clusters.values())
        row = dict(scorer=scorer, axis_sha256=axis, n_pairs=len(arr),n_clusters=len(groups), mean_difference=float(arr.mean()),
                   bootstrap_ci_low=".", bootstrap_ci_high=".", sign_flip_p=".", bh_q=".",
                   status="too_few_pairs", inference="exploratory_matched_observational")
        if len(arr) >= min_pairs and len(groups)>=2:
            means, permuted = [], []
            for start in range(0, draws, 100):
                size = min(100, draws - start)
                totals=np.asarray([sum(g) for g in groups]);counts=np.asarray([len(g) for g in groups])
                chosen=rng.integers(0,len(groups),size=(size,len(groups)))
                means.extend((totals[chosen].sum(axis=1)/counts[chosen].sum(axis=1)).tolist())
                signs=rng.choice((-1.,1.),size=(size,len(groups)))
                permuted.extend(((signs*totals).sum(axis=1)/len(arr)).tolist())
            p = (1 + sum(abs(x) >= abs(row["mean_difference"]) for x in permuted)) / (draws + 1)
            lo, hi = np.quantile(means, [.025, .975])
            row.update(bootstrap_ci_low=float(lo), bootstrap_ci_high=float(hi), sign_flip_p=p, status="ok")
            tested.append((len(results), p))
        results.append(row)
    for (index, _), q in zip(tested, _bh([p for _, p in tested])):
        results[index]["bh_q"] = q
    fields = ["scorer", "axis_sha256", "n_pairs", "n_clusters", "mean_difference", "bootstrap_ci_low", "bootstrap_ci_high",
              "sign_flip_p", "bh_q", "status", "inference"]
    handle, writer = tsv_writer(out / "matched_functional_summary.tsv", fields)
    with handle:
        writer.writerows(results)
    result = dict(input=fingerprint(geometry_path), draws=draws, seed=seed, min_pairs=min_pairs,
                  metric="(distance_to_decreasing-distance_to_increasing)/(distance_to_decreasing+distance_to_increasing)",
                  qc=dict(qc), resampling="Whole case LD/proximity clusters; each case compared with the mean of all available exact-feature-axis controls",
                  multiple_testing="BH over tested scorer x exact feature-panel comparisons in this output only",
                  limitation="No causal or balancing-selection claim; approximate physical blocks do not establish independent LD loci. "
                             "Exact biological feature axes are matched, but correlated tracks and heterogeneous biological effects remain.",
                  outputs=[str(out / name) for name in ("site_closeness.tsv", "paired_closeness.tsv", "matched_functional_summary.tsv")])
    dump_json(out / "comparison_manifest.json", result)
    return result


def write_report(out_dir):
    """Small self-contained HTML index of actual local stage outputs."""
    out = Path(out_dir).resolve()
    cards, sections = [], []
    manifests = [("counts_audit", "原始 check 计数审计", "manifest.json"),
                 ("metadata", "样本信息", "metadata_qc.json"), ("join", "现代—古代选择统计连接", "manifest.json"),
                 ("candidates", "候选与匹配", "manifest.json"), ("universes", "U0–U3 独立分析集", "manifest.json"),
                 ("ukb", "UKB MFI 精确匹配", "U0_mfi.tsv.manifest.json"), ("panel", "古代面板覆盖", "manifest.json"),
                 ("trajectories", "频率轨迹", "trajectory_manifest.json"),
                 ("atlas", "稀疏功能查询", "query.json"), ("geometry", "分子效应几何", "geometry.json")]
    for folder, title, filename in manifests:
        path = out / folder / filename
        if not path.exists():
            continue
        data = json.loads(path.read_text())
        values = data.get("qc", data.get("target_counts", data.get("counts", data.get("status_counts", {}))))
        if values:
            rows = "".join(f"<tr><th>{html.escape(str(k))}</th><td>{html.escape(str(v))}</td></tr>" for k, v in values.items())
            cards.append(f"<section><h2>{title}</h2><table>{rows}</table></section>")
    for relative, title in [("figures/figure_manifest.json", "Fig1–Fig9 实际产物与缺失状态"),
                             ("figures/integrated_evidence.tsv", "观测证据分量与描述性分数"),
                             ("universes/U2_all_strong.tsv", "全部合格强信号"),
                             ("universes/U2_clustered_bootstrap.tsv", "按 LD block 聚类的 bootstrap"),
                             ("shared_metag_census.json", "共享 metag 全基因组 ERV/LTR 覆盖率"),
                             ("pca/pca_time.html", "PC1–PC2–时间：交互式三维图"),
                             ("counts_audit/counts_review.html", "旧 check 计数：按染色体的描述性汇总"),
                             ("trajectories/trajectory_bins.tsv", "按年代与地区分层的轨迹表"),
                             ("trajectories/trajectory_standardized.tsv", "固定构成标准化轨迹"),
                             ("panel/site_variant_links.tsv", "古代面板是否覆盖额外碱基"),
                             ("candidates/candidates.tsv", "候选位点与匹配对照"),
                             ("comparison/matched_functional_summary.tsv", "探索性功能比较"),
                             ("pca/commands.sh", "PCA 计算脚本")]:
        if (out / relative).exists():
            sections.append(f'<li><a href="{relative}">{title}</a></li>')
    text = """<!doctype html><html lang="zh"><meta charset="utf-8"><title>allele3 evolution</title>
<style>body{font:16px system-ui,sans-serif;background:#f4f7fb;color:#172b46;max-width:1150px;margin:32px auto;padding:0 24px}
h1{font-size:32px}h2{font-size:20px}section{background:white;border:1px solid #dde4ef;border-radius:12px;padding:20px;margin:16px 0}
table{border-collapse:collapse;width:100%}th,td{text-align:left;padding:7px;border-bottom:1px solid #edf1f7}a{color:#135fb2}.note{border-left:4px solid #e4aa42;padding:12px 18px;background:#fff8e8}</style>
<h1>allele3：现代分布、古代轨迹与分子效应</h1>
<p>此报告只汇总本地已完成步骤。空缺的 Atlas 或 PCA 结果表示该步骤尚未执行，不代表零效应或没有结构。</p>
<div class="note">REF 不是祖先状态；未列入古代双等位插补面板的第三碱基属于未覆盖。
本流程的选择统计来自已发表的精确等位基因对比，不把它转赋给第三碱基；频率曲线与 PCA 不能单独证明选择。</div>
<h2>结果入口</h2><ul>""" + "".join(sections) + "</ul>" + "".join(cards) + "</html>"
    (out / "report.html").write_text(text, encoding="utf-8")
    return str(out / "report.html")
