"""Read-only arithmetic audit and descriptive report for legacy counts.tsv.

These are conditional counts of legacy PASS/positive-AC ALT position groups.
REF occurrence, callable sequence, release/build and evolutionary processes
cannot be recovered from this summary. No hypothesis tests are performed.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from decimal import Decimal, InvalidOperation, localcontext
import base64
import csv
import gzip
import hashlib
import html
import json
import math
import os
from pathlib import Path
import re
import tempfile


TOOL_ID = "allele3.counts_audit"
SCHEMA_VERSION = 1
CATEGORIES = ("singleton", "non_singleton_biallelic", "tri_allelic", "four_plus_allelic")
TYPES = ("SNP", "INDEL")
REQUIRED = ("CHR", "TYPE", "CATEGORY", "COUNT", "TOTAL_SITES", "PERCENT")
AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))
STANDARD_CHROMS = AUTOSOMES + ("chrX", "chrY")
OUTPUT_NAMES = ("counts_audit.json", "counts_summary.tsv", "counts_by_chromosome.tsv",
                "counts_review.html", "counts_review.png", "manifest.json")
DEFINITIONS = {
    "expected_producer": "Reviewed legacy check implementation; a matching TSV header alone does not prove its producer or settings.",
    "retained_alt": "ALT retained by the legacy producer after PASS and AC > 0 filters.",
    "singleton": "Exactly one retained ALT and its AC equals 1.",
    "non_singleton_biallelic": "Exactly one retained ALT and its AC exceeds 1; the name does not establish observed REF copies.",
    "tri_allelic": "Exactly two retained ALT entries/states at the legacy position group; REF is an annotation label, not evidence of an observed third state.",
    "four_plus_allelic": "At least three retained ALT entries/states. Canonical A/C/G/T SNPs have at most three ALT bases relative to REF; observed REF is not established.",
    "INDEL": "Legacy groups use the start coordinate and were not fully FASTA-left-aligned; category size is not a reliable count of normalized polymorphic states.",
    "multi_alt": "tri_allelic + four_plus_allelic, i.e. legacy groups with at least two retained ALTs.",
    "alt_3plus": "four_plus_allelic, i.e. legacy groups with at least three retained ALTs.",
    "denominator": "TOTAL_SITES within each CHR/TYPE; aggregate fractions use summed counts divided by summed totals.",
}
LIMITATIONS = [
    "The counts file contains no release, reference-build, sample-size, population, or complete producer-configuration metadata.",
    "Rates are conditional fractions among legacy position groups, not densities per callable genomic base.",
    "REF labels do not establish observed REF copies. Two retained ALTs do not by themselves prove three observed alleles.",
    "INDEL start-coordinate groups are not fully reference-normalized; do not interpret their categories as normalized polymorphic-state counts.",
    "Missing chromosome/type groups are unavailable, not zero. X and Y are reported separately from autosomes.",
    "SNP and INDEL totals can overlap genomic positions; their sum is not asserted to count unique coordinates.",
    "Leave-one-chromosome ratios are descriptive and unadjusted for callability, ancestry, allele frequency, sequence context, or QC.",
    "No P values, confidence intervals, selection estimates, pleiotropy effects, or ancient-DNA conclusions are inferred.",
]
SUMMARY_FIELDS = ["scope", "type", "status", "observed_chromosomes", "missing_chromosomes", "total_sites"]
for _category in CATEGORIES:
    SUMMARY_FIELDS += [_category + "_count", _category + "_percent"]
SUMMARY_FIELDS += ["multi_alt_count", "multi_alt_percent", "alt_3plus_count", "alt_3plus_percent", "alt_3plus_within_multi_percent"]
CHROM_FIELDS = ["chrom", "type", "status", "total_sites"] + SUMMARY_FIELDS[6:] + [
    "reference_scope", "reference_total_sites", "multi_alt_reference_percent",
    "multi_alt_ratio_to_reference", "alt_3plus_reference_percent", "alt_3plus_ratio_to_reference",
]


def _open(path):
    if str(path).endswith((".gz", ".bgz")):
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return open(path, "rt", encoding="utf-8-sig", newline="")


def _sha(path):
    hasher = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1048576), b""):
            hasher.update(block)
    return hasher.hexdigest()


def _fingerprint(path):
    path = Path(path).resolve()
    return {"path": str(path), "size": path.stat().st_size, "sha256": _sha(path)}


def _chrom(value):
    value = value.strip()
    short = value[3:] if value.lower().startswith("chr") else value
    if short.isdigit() and 1 <= int(short) <= 22:
        return "chr" + str(int(short))
    if short.upper() in {"X", "Y"}:
        return "chr" + short.upper()
    if not value or any(c.isspace() for c in value):
        raise ValueError(f"Invalid chromosome label: {value!r}")
    return value


def _order(chrom):
    if chrom in STANDARD_CHROMS:
        return (STANDARD_CHROMS.index(chrom), "")
    return (100, chrom)


def _integer(value, where):
    if not re.fullmatch(r"[0-9]+", value.strip()):
        raise ValueError(f"{where}: expected a nonnegative integer, got {value!r}")
    return int(value)


def _read_counts(path):
    grouped = defaultdict(list)
    with _open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)) or not set(REQUIRED).issubset(fields):
            raise ValueError("Counts table requires unique CHR,TYPE,CATEGORY,COUNT,TOTAL_SITES,PERCENT columns")
        for line, raw in enumerate(reader, 2):
            if None in raw or any(v is None for v in raw.values()):
                raise ValueError(f"Counts line {line}: wrong number of fields")
            row = {key: value.strip() for key, value in raw.items()}
            chrom, typ, category = _chrom(row["CHR"]), row["TYPE"], row["CATEGORY"]
            if typ not in TYPES or category not in CATEGORIES:
                raise ValueError(f"Counts line {line}: unknown TYPE/CATEGORY {typ!r}/{category!r}")
            count = _integer(row["COUNT"], f"line {line} COUNT")
            total = _integer(row["TOTAL_SITES"], f"line {line} TOTAL_SITES")
            percent_text = row["PERCENT"]
            if total == 0 and percent_text in {".", "NA", "N/A"}:
                percent = None
            else:
                if not re.fullmatch(r"[0-9]+(?:\.[0-9]+)?", percent_text):
                    raise ValueError(f"line {line} PERCENT: expected a finite plain decimal in [0,100]")
                try:
                    percent = Decimal(percent_text)
                except InvalidOperation as exc:
                    raise ValueError(f"line {line}: invalid PERCENT") from exc
                if not percent.is_finite() or not 0 <= percent <= 100:
                    raise ValueError(f"line {line} PERCENT: expected a finite decimal in [0,100]")
            grouped[(chrom, typ)].append(dict(line=line, category=category, count=count,
                                              total=total, percent=percent, percent_text=percent_text))
    if not grouped:
        raise ValueError("Counts table contains no data")
    groups, audits, errors = [], [], []
    for (chrom, typ), rows in sorted(grouped.items(), key=lambda item: (_order(item[0][0]), item[0][1])):
        label = f"{chrom}/{typ}"
        if Counter(r["category"] for r in rows) != Counter(CATEGORIES):
            errors.append(f"{label}: each of the four categories must occur exactly once")
            continue
        totals = {r["total"] for r in rows}
        if len(totals) != 1:
            errors.append(f"{label}: repeated TOTAL_SITES values disagree")
            continue
        total = next(iter(totals))
        if sum(r["count"] for r in rows) != total:
            errors.append(f"{label}: sum(COUNT) does not equal TOTAL_SITES")
            continue
        checks = []
        for row in rows:
            if total == 0:
                if row["percent"] not in (None, Decimal(0)):
                    errors.append(f"{label}/{row['category']}: zero denominator cannot have a nonzero PERCENT")
                exact, delta, tolerance = None, None, None
            else:
                places = max(0, -row["percent"].as_tuple().exponent)
                with localcontext() as ctx:
                    ctx.prec = max(60, len(str(total)) + 20)
                    exact = Decimal(row["count"]) * 100 / Decimal(total)
                    delta = row["percent"] - exact
                    tolerance = Decimal(1).scaleb(-places) / 2
                    if abs(delta) > tolerance:
                        errors.append(f"{label}/{row['category']}: PERCENT disagrees with COUNT/TOTAL beyond rounding tolerance")
            checks.append(dict(category=row["category"], line=row["line"], reported_percent=row["percent_text"],
                               exact_percent=str(exact) if exact is not None else None,
                               error_percentage_points=str(delta) if delta is not None else None,
                               rounding_tolerance_percentage_points=str(tolerance) if tolerance is not None else None))
        counts = {r["category"]: r["count"] for r in rows}
        group = dict(chrom=chrom, type=typ, total_sites=total, counts=counts)
        groups.append(group)
        audits.append(dict(chrom=chrom, type=typ, category_count=4, category_sum=sum(counts.values()),
                           total_sites=total, reported_percent_sum=str(sum((r["percent"] or Decimal(0)) for r in rows)),
                           zero_denominator=total == 0, percent_checks=checks))
    if errors:
        raise ValueError("Counts audit failed: " + "; ".join(errors[:12]) + (f"; and {len(errors)-12} more" if len(errors)>12 else ""))
    return groups, audits, fields


def _percent(count, total):
    return 100.0 * count / total if count is not None and total else None


def _values(total, counts):
    result = {"total_sites": total}
    for category in CATEGORIES:
        count = counts.get(category) if counts is not None else None
        result[category + "_count"] = count
        result[category + "_percent"] = _percent(count, total)
    multi = counts["tri_allelic"] + counts["four_plus_allelic"] if counts is not None else None
    third = counts["four_plus_allelic"] if counts is not None else None
    result.update(multi_alt_count=multi, multi_alt_percent=_percent(multi, total),
                  alt_3plus_count=third, alt_3plus_percent=_percent(third, total),
                  alt_3plus_within_multi_percent=_percent(third, multi))
    return result


def _aggregate(groups, chromosomes, typ, scope):
    selected = [r for r in groups if r["chrom"] in chromosomes and r["type"] == typ]
    observed = sorted({r["chrom"] for r in selected}, key=_order)
    missing = [c for c in chromosomes if c not in observed]
    if not selected:
        status, total, counts = "missing", None, None
    else:
        total = sum(r["total_sites"] for r in selected)
        counts = {c: sum(r["counts"][c] for r in selected) for c in CATEGORIES}
        status = "zero_denominator" if total == 0 else "partial_coverage" if missing else "observed"
    return dict(scope=scope, type=typ, status=status, observed_chromosomes=observed,
                missing_chromosomes=missing, **_values(total, counts))


def _summaries(groups):
    observed = sorted({r["chrom"] for r in groups}, key=_order)
    scopes = [("autosomes", list(AUTOSOMES)), ("X", ["chrX"]), ("Y", ["chrY"]), ("all_observed", observed)]
    other = [c for c in observed if c not in STANDARD_CHROMS]
    if other:
        scopes.append(("other_contigs", other))
    summaries = [_aggregate(groups, chroms, typ, scope) for scope, chroms in scopes for typ in TYPES]
    by_chromosome = []
    for chrom in sorted(set(STANDARD_CHROMS) | set(observed), key=_order):
        for typ in TYPES:
            row = _aggregate(groups, [chrom], typ, chrom)
            row.pop("scope")
            row["chrom"] = chrom
            if chrom in AUTOSOMES:
                ref_chroms = [c for c in AUTOSOMES if c != chrom]
                reference_scope = "other_observed_autosomes"
            elif chrom in {"chrX", "chrY"}:
                ref_chroms = list(AUTOSOMES)
                reference_scope = "observed_autosomes"
            else:
                ref_chroms, reference_scope = [], "not_defined_for_other_contigs"
            reference = _aggregate(groups, ref_chroms, typ, reference_scope)
            row.update(reference_scope=reference_scope, reference_total_sites=reference["total_sites"],
                       reference_missing_chromosomes=reference["missing_chromosomes"])
            comparisons = {}
            for metric in (*CATEGORIES, "multi_alt", "alt_3plus"):
                numerator = row[metric + "_percent"]
                denominator = reference[metric + "_percent"]
                ratio = numerator / denominator if numerator is not None and denominator else None
                comparisons[metric] = dict(percent=numerator, reference_percent=denominator, ratio=ratio,
                                           difference_percentage_points=numerator-denominator if numerator is not None and denominator is not None else None)
                if metric in {"multi_alt", "alt_3plus"}:
                    row[metric + "_reference_percent"] = denominator
                    row[metric + "_ratio_to_reference"] = ratio
            row["comparisons"] = comparisons
            by_chromosome.append(row)
    extremes = {}
    for typ in TYPES:
        eligible = [r for r in by_chromosome if r["type"] == typ and r["chrom"] in AUTOSOMES and r["total_sites"]]
        extremes[typ] = {}
        for metric in ("multi_alt_percent", "alt_3plus_percent", "multi_alt_count", "alt_3plus_count", "total_sites"):
            extremes[typ][metric] = ({side: {"chrom": row["chrom"], "value": row[metric]}
                                     for side, row in (("min", min(eligible, key=lambda r:r[metric])),
                                                       ("max", max(eligible, key=lambda r:r[metric])))} if eligible else None)
    return summaries, by_chromosome, extremes


def _plot_backend(make_plots):
    if not make_plots:
        return None, "not_requested", None
    try:
        import matplotlib
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg
    except ImportError:
        return None, "unavailable_matplotlib", None
    return (Figure, FigureCanvasAgg), "available", matplotlib.__version__


def _make_plot(path, by_chromosome, summaries, backend):
    Figure, FigureCanvasAgg = backend
    fig = Figure(figsize=(14, 9.2), facecolor="white")
    FigureCanvasAgg(fig)
    axes = fig.subplots(2, 2)
    lookup = {(r["chrom"], r["type"]):r for r in by_chromosome}
    baseline = {r["type"]:r for r in summaries if r["scope"] == "autosomes"}
    positions = {c:i for i,c in enumerate(AUTOSOMES)} | {"chrX":23, "chrY":24.5}
    tick_labels = [str(i) for i in range(1,23)] + ["X", "Y"]
    colors = {"SNP":"#2166ac", "INDEL":"#16856c"}
    for col, typ in enumerate(TYPES):
        for row, metric in enumerate(("multi_alt_percent", "alt_3plus_percent")):
            ax = axes[row,col]
            records = [lookup[(chrom,typ)] for chrom in STANDARD_CHROMS]
            vals = [r[metric] for r in records if r[metric] is not None]
            for chrom in AUTOSOMES:
                value = lookup[(chrom,typ)][metric]
                if value is not None:
                    ax.scatter(positions[chrom], value, s=39, color=colors[typ], zorder=3)
            for chrom in ("chrX","chrY"):
                value = lookup[(chrom,typ)][metric]
                if value is not None:
                    ax.scatter(positions[chrom], value, marker="D", s=46, color="#d87920", zorder=4)
                else:
                    ax.text(positions[chrom], .025, "NA", transform=ax.get_xaxis_transform(),
                            ha="center", va="bottom", color="#9ca3af", fontsize=8)
            value = baseline[typ][metric]
            if value is not None:
                ax.axhline(value, color="#64748b", lw=1.15, ls=(0,(5,3)), zorder=1)
                ax.text(.02,.97,f"Autosomes weighted: {value:.4f}%", transform=ax.transAxes,
                        va="top", fontsize=10, color="#4b5563")
            valid_auto = [r for r in records if r["chrom"] in AUTOSOMES and r[metric] is not None]
            if valid_auto:
                highest, lowest = max(valid_auto,key=lambda r:r[metric]), min(valid_auto,key=lambda r:r[metric])
                for target, dy in ((highest,11),(lowest,-18)):
                    if highest is lowest and dy < 0:
                        continue
                    label = target["chrom"].removeprefix("chr")
                    ax.annotate(f"{label}: {target[metric]:.4f}%", (positions[target["chrom"]],target[metric]),
                                xytext=(0,dy), textcoords="offset points", ha="center", fontsize=9,
                                color=colors[typ], bbox=dict(facecolor="white",edgecolor="none",alpha=.9,pad=1))
            ax.axvline(22, color="#d3dce5", lw=1)
            ax.set_xticks([positions[c] for c in STANDARD_CHROMS],tick_labels, fontsize=8.5)
            for tick, chrom in zip(ax.get_xticklabels(), STANDARD_CHROMS):
                if lookup[(chrom,typ)]["status"] == "missing":
                    tick.set_color("#a4adb7")
            ax.set_xlim(-.8,25.3)
            ymax = max(vals) * 1.28 if vals and max(vals)>0 else 1.0
            ax.set_ylim(0,ymax)
            ax.set_ylabel("Percent of legacy position groups", fontsize=10)
            ax.set_xlabel("Chromosome",fontsize=10)
            label = "≥2 retained ALT" if row==0 else "≥3 retained ALT"
            ax.set_title(f"{'ABCD'[row*2+col]}   {typ}  |  {label}",loc="left",fontsize=13,fontweight="bold",pad=10)
            ax.grid(axis="y",color="#e5ebf0",lw=.7)
            ax.set_axisbelow(True)
            ax.spines[["top","right"]].set_visible(False)
            ax.spines[["left","bottom"]].set_color("#c4cdd6")
    fig.suptitle("PASS, AC>0 ALT categories by chromosome",fontsize=19,fontweight="bold",x=.5,y=.982)
    fig.text(.5,.942,"Autosomes: panel colour     X/Y: orange diamond     Dashed line: observed-autosome weighted fraction",ha="center",fontsize=10,color="#475569")
    fig.text(.055,.057,"Legacy conditional fractions, not callable-base densities. REF occurrence is not established.",fontsize=10,color="#475569")
    fig.text(.055,.034,"INDELs share start-coordinate groups without full reference left-alignment; their categories are not normalized allele-state counts.",fontsize=9.5,color="#475569")
    fig.text(.055,.012,"NA = chromosome/type not provided or zero denominator. No P values, confidence intervals, or evolutionary conclusions are shown.",fontsize=9.5,color="#475569")
    fig.subplots_adjust(top=.88,bottom=.15,left=.06,right=.98,hspace=.40,wspace=.20)
    fig.savefig(path,dpi=175,facecolor="white")


def _fmt(value, percent=False):
    if value is None:
        return "NA"
    if isinstance(value,float):
        return f"{value:.4f}%" if percent else f"{value:.6g}"
    return f"{value:,}" if isinstance(value,int) else str(value)


def _table(headers, rows):
    return '<div class="table-scroll"><table><thead><tr>'+''.join('<th>'+html.escape(h)+'</th>' for h in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in row)+'</tr>' for row in rows)+'</tbody></table></div>'


def _report(audit, image_path, plot_status):
    summary = audit["summaries"]
    table_rows = []
    scope_names = {"autosomes":"常染色体（已提供）", "X":"X", "Y":"Y", "all_observed":"全部已提供染色体", "other_contigs":"其他 contig"}
    for r in summary:
        table_rows.append([scope_names[r["scope"]],r["type"],r["status"],_fmt(r["total_sites"]),
                           _fmt(r["singleton_count"]),_fmt(r["non_singleton_biallelic_count"]),
                           _fmt(r["tri_allelic_count"]),_fmt(r["four_plus_allelic_count"]),
                           _fmt(r["multi_alt_percent"],True),_fmt(r["alt_3plus_percent"],True)])
    by_rows = []
    for r in audit["by_chromosome"]:
        by_rows.append([r["chrom"],r["type"],r["status"],_fmt(r["total_sites"]),
                        _fmt(r["multi_alt_percent"],True),_fmt(r["alt_3plus_percent"],True),
                        _fmt(r["multi_alt_ratio_to_reference"]),_fmt(r["alt_3plus_ratio_to_reference"])])
    image_block = ('<img class="figure" alt="SNP and INDEL retained ALT category percentages by chromosome" src="data:image/png;base64,'+
                   base64.b64encode(Path(image_path).read_bytes()).decode()+'">' if image_path else
                   '<p class="note">图未生成：'+html.escape(plot_status)+'。下方表格和 JSON/TSV 汇总仍完整可用。</p>')
    qc = audit["qc"]
    missing = '; '.join(typ+': '+(', '.join(qc['missing_standard_chromosomes_by_type'][typ]) or '无') for typ in TYPES)
    status_counts = json.dumps(qc["reported_percent_sum_distribution"],ensure_ascii=False)
    document = '''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>counts.tsv 条件计数审计</title><style>
body{margin:0;background:#edf2f6;color:#183247;font:15px/1.65 system-ui,-apple-system,"Segoe UI",sans-serif}main{max-width:1440px;margin:auto;padding:24px}header{padding:24px 30px;background:#17374f;color:white;border-radius:12px}h1{margin:0 0 10px;font-size:29px}h2{font-size:21px;margin:0 0 14px}.card{background:white;border:1px solid #d7e1e9;border-radius:10px;padding:22px;margin:18px 0}.note{background:#fff5df;border-left:4px solid #d69b38;padding:14px 17px}.figure{display:block;width:100%;height:auto}.table-scroll{overflow:auto}table{width:100%;border-collapse:collapse;font-size:13px;font-variant-numeric:tabular-nums}th,td{padding:9px 10px;border-bottom:1px solid #e0e8ee;text-align:right;white-space:nowrap}th{background:#edf4f8;color:#24445d}th:first-child,td:first-child{text-align:left}code{font-family:ui-monospace,monospace;word-break:break-all}small{color:#577084}li{margin:7px 0}.pill{display:inline-block;border:1px solid #8fc9ba;background:#e6f7f1;color:#175d4d;border-radius:5px;padding:3px 9px}a{color:#2166ac}footer{font-size:13px;color:#526c7e;margin:20px 4px}
</style></head><body><main><header><h1>counts.tsv 条件计数审计</h1><p>按已核验 legacy check 定义：统计 PASS 且 AC&gt;0 的 ALT 数量。所有比例采用各类型自己的 TOTAL_SITES 分母。</p></header>'''
    document += '<section class="card"><h2>数值核对</h2><p><span class="pill">算术校验通过</span> '+f"{qc['rows']} 个数据行，{qc['groups']} 个 CHR/TYPE 组；四类唯一、同组 TOTAL 一致、COUNT 合计和百分比舍入全部通过。"+'</p>'
    document += '<p>缺失标准染色体：'+html.escape(missing)+'。缺失项保留 NA，不补零。</p><p>四类已报告百分比之和（值→组数）：<code>'+html.escape(status_counts)+'</code>；逐项按其显示小数精度检查舍入。</p>'
    document += '<p><small>输入：'+html.escape(audit['input']['path'])+'<br>SHA-256：<code>'+audit['input']['sha256']+'</code></small></p></section>'
    document += '<section class="card"><h2>四类究竟表示什么</h2><ul><li><b>singleton</b>：仅一个保留 ALT，且其 AC=1。</li><li><b>non_singleton_biallelic</b>：仅一个保留 ALT，且 AC&gt;1。</li><li><b>tri_allelic</b>：两个保留 ALT。</li><li><b>four_plus_allelic</b>：至少三个保留 ALT；标准 A/C/G/T SNP 相对 REF 最多只有三个 ALT 碱基。</li></ul><p class="note"><b>REF 只是标签，不证明样本中存在 REF 拷贝。</b> 因而两个 ALT 不能自动称为三个已观测等位状态。INDEL 是同起始坐标分组，未完整结合 FASTA 左对齐，类别数也不能直接当作规范化多态状态数。</p><p>以上是此模块所要求的已核验 legacy producer 语义；其他来源即使列名相同也需核对源码。输入表本身没有 release/build、样本量、群体或完整过滤参数元数据。</p></section>'
    document += '<section class="card"><h2>染色体分布</h2>'+image_block+'<p><small>≥2 ALT = tri_allelic + four_plus_allelic；≥3 ALT = four_plus_allelic。虚线为已提供常染色体按 TOTAL 加权的比例，X/Y 单独展示。图的分母不是可调用碱基数。</small></p></section>'
    document += '<section class="card"><h2>汇总计数与加权比例</h2>'+_table(['范围','TYPE','状态','TOTAL','singleton','单ALT, AC>1','两ALT','至少三ALT','≥2 ALT 比例','≥3 ALT 比例'],table_rows)+'<p>完整四类加权比例另见 counts_summary.tsv；各类型 totals 可能覆盖相同坐标，其相加不代表唯一基因组位置数。</p></section>'
    document += '<section class="card"><h2>逐染色体与描述性比例比</h2><p>常染色体的参考为其余已提供常染色体；X/Y 参考为全部已提供常染色体。比例比是未校正的描述量，不是选择、功能或富集检验。</p>'+_table(['CHR','TYPE','状态','TOTAL','≥2 ALT','≥3 ALT','≥2 ALT / 参考','≥3 ALT / 参考'],by_rows)+'</section>'
    document += '<section class="card"><h2>结果可以支持的范围</h2><p>这些数据支持分类内计数、条件比例、染色体描述性比较和候选计算规模估算。它们不能单独给出可调用密度、变异产生率、祖先状态、首次出现时间、平衡选择或 pleiotropy 关系；未计算 P 值、置信区间或古代 DNA 结论。</p><p>较高的条件比例可用于安排后续序列背景、原始位点、规范化及 QC 审计；不能仅凭本表判定其生物学原因。</p></section>'
    return document+'<footer>独立离线 HTML：图像以 base64 内嵌，无 CDN、脚本或外部字体请求。完整机器可读审计见 counts_audit.json 和 manifest.json。</footer></main></body></html>'


def _write_tsv(path, fields, rows):
    with open(path,"w",encoding="utf-8",newline="") as handle:
        writer = csv.DictWriter(handle,fields,delimiter="\t",lineterminator="\n")
        writer.writeheader()
        for row in rows:
            clean = {}
            for name in fields:
                value = row.get(name)
                clean[name] = "." if value is None else ",".join(value) if isinstance(value,list) else format(value,".12g") if isinstance(value,float) else value
            writer.writerow(clean)


def _existing_manifest(destination):
    manifest_path = destination / "manifest.json"
    if manifest_path.exists():
        try:
            old = json.loads(manifest_path.read_text())
        except (ValueError,OSError) as exc:
            raise ValueError("Existing manifest.json is not a valid owned counts-audit manifest") from exc
        if old.get("tool") != TOOL_ID or old.get("schema_version") != SCHEMA_VERSION or not isinstance(old.get("outputs"),dict):
            raise ValueError("Refusing to overwrite outputs without an owned counts-audit manifest")
        owned = set(old["outputs"]) | {"manifest.json"}
        if not owned.issubset(OUTPUT_NAMES):
            raise ValueError("Existing counts-audit manifest lists unexpected output names")
    else:
        old, owned = None, set()
    for name in OUTPUT_NAMES:
        if (destination/name).exists() and name not in owned:
            raise ValueError(f"Refusing to overwrite unowned output: {destination/name}")
    return old


def audit_counts(counts_path, out_dir, make_plots=True):
    """Validate, summarize and publish owned offline outputs; never edit input.

    Invalid input raises ValueError before any output is replaced. Identical
    inputs/settings reuse checksum-verified outputs. Changed runs replace only
    outputs declared by this module's manifest, using staged atomic file moves.
    Matplotlib is optional; TSV/JSON/HTML are still produced if it is absent.
    """
    if not isinstance(make_plots,bool):
        raise ValueError("make_plots must be True or False")
    source = Path(counts_path).expanduser().resolve()
    destination = Path(out_dir).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    if source in {(destination/name).resolve() for name in OUTPUT_NAMES}:
        raise ValueError("Counts output must not overwrite its input")
    source_info = _fingerprint(source)
    groups, group_audit, source_columns = _read_counts(source)
    if _sha(source) != source_info["sha256"]:
        raise ValueError("Counts input changed while it was being read")
    summaries, by_chromosome, extremes = _summaries(groups)
    backend, backend_status, version = _plot_backend(make_plots)
    settings = dict(make_plots=make_plots,plot_backend_status=backend_status,matplotlib_version=version)
    code_hash = _sha(Path(__file__))
    signature = hashlib.sha256(json.dumps(dict(input=source_info,settings=settings,code_sha256=code_hash),sort_keys=True).encode()).hexdigest()
    destination.mkdir(parents=True,exist_ok=True)
    old = _existing_manifest(destination)
    if old and old.get("signature") == signature and all((destination/name).is_file() and _sha(destination/name)==spec["sha256"] for name,spec in old["outputs"].items()):
        return old
    missing = {typ:[c for c in STANDARD_CHROMS if not any(r["chrom"]==c and r["type"]==typ for r in groups)] for typ in TYPES}
    other = sorted({r["chrom"] for r in groups if r["chrom"] not in STANDARD_CHROMS})
    errors = [abs(float(c["error_percentage_points"])) for r in group_audit for c in r["percent_checks"] if c["error_percentage_points"] is not None]
    qc = dict(status="passed",rows=len(groups)*4,groups=len(groups),four_categories_unique=True,
              totals_consistent=True,category_sums_match=True,percentages_within_display_rounding=True,
              max_absolute_percent_error_pp=max(errors) if errors else None,
              reported_percent_sum_distribution=dict(Counter(r["reported_percent_sum"] for r in group_audit)),
              zero_denominator_groups=sum(r["total_sites"]==0 for r in groups),
              missing_standard_chromosomes_by_type=missing,other_contigs=other,
              missing_types=[typ for typ in TYPES if not any(r["type"]==typ for r in groups)])
    audit = dict(tool=TOOL_ID,schema_version=SCHEMA_VERSION,input=source_info,source_columns=source_columns,
                 release=None,reference_build=None,producer_metadata_in_input=False,
                 definitions=DEFINITIONS,limitations=LIMITATIONS,qc=qc,group_audit=group_audit,
                 summaries=summaries,by_chromosome=by_chromosome,autosomal_extremes=extremes)
    output_files = [name for name in OUTPUT_NAMES if name not in {"manifest.json","counts_review.png"}]
    plot_status = backend_status
    with tempfile.TemporaryDirectory(prefix=".counts-audit-",dir=destination) as temporary:
        stage = Path(temporary)
        plot_path = None
        if backend is not None:
            plot_path = stage/"counts_review.png"
            _make_plot(plot_path,by_chromosome,summaries,backend)
            output_files.append("counts_review.png")
            plot_status = "generated"
        (stage/"counts_audit.json").write_text(json.dumps(audit,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
        _write_tsv(stage/"counts_summary.tsv",SUMMARY_FIELDS,summaries)
        _write_tsv(stage/"counts_by_chromosome.tsv",CHROM_FIELDS,by_chromosome)
        (stage/"counts_review.html").write_text(_report(audit,plot_path,plot_status),encoding="utf-8")
        if _sha(source) != source_info["sha256"]:
            raise ValueError("Counts input changed during report generation")
        outputs = {name:{"path":str(destination/name),"size":(stage/name).stat().st_size,"sha256":_sha(stage/name)} for name in output_files}
        manifest = dict(tool=TOOL_ID,schema_version=SCHEMA_VERSION,signature=signature,input=source_info,
                        code_sha256=code_hash,settings=settings,plot_status=plot_status,outputs=outputs,
                        manifest_path=str(destination/"manifest.json"),report_path=str(destination/"counts_review.html"),
                        audit_path=str(destination/"counts_audit.json"),summary_path=str(destination/"counts_summary.tsv"),
                        by_chromosome_path=str(destination/"counts_by_chromosome.tsv"),
                        plot_path=str(destination/"counts_review.png") if plot_path else None,
                        qc=qc,definitions=DEFINITIONS,limitations=LIMITATIONS)
        (stage/"manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
        for name in output_files:
            os.replace(stage/name,destination/name)
        if old:
            for name in set(old["outputs"])-set(output_files):
                (destination/name).unlink(missing_ok=True)
        os.replace(stage/"manifest.json",destination/"manifest.json")
    return manifest
