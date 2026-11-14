#!/usr/bin/env python3
"""Summarize an existing allele3 result snapshot without re-running discovery.

Only explicit table contents are summarized. Empty annotations, unmatched
windows and an incomplete GWAS report remain visible. No HERV fraction or
ancient-origin finding is inferred from missing files.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import csv
import html
import json
from pathlib import Path
import re
from analysis_integrity import identity, transactional


def read_table(path):
    if not path.is_file():
        return []
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_table(path, rows, fields):
    with path.open("w", newline="") as f:
        w=csv.DictWriter(f, fields, delimiter="\t", lineterminator="\n"); w.writeheader()
        for row in rows: w.writerow({k:row.get(k, "NA") for k in fields})


def html_table(rows, fields):
    return '<table><thead><tr>'+''.join('<th>'+html.escape(label)+'</th>' for key,label in fields)+'</tr></thead><tbody>'+''.join(
        '<tr>'+''.join('<td>'+html.escape(str(row.get(key, 'NA')))+'</td>' for key,label in fields)+'</tr>' for row in rows)+'</tbody></table>'


@transactional("qc/manifest.json", input_dirs=("analysis_dir",))
def run(args):
    root=Path(args.analysis_dir)
    files={
        'occupancy':root/'gnomad_scan/results/occupancy_thresholds.tsv',
        'singletons':root/'gnomad_scan/results/singleton_sensitivity.tsv',
        'nominated':root/'hotspots/results/windows.nominated_windows.tsv',
        'selected':root/'hotspots/results/characterize.selected_windows.tsv',
        'matches':root/'hotspots/results/characterize.window_matches.tsv',
        'annotations':root/'hotspots/results/characterize.annotation_contrasts.tsv',
        'atlas':root/'hotspots/results/impact.window_atlas_summary.tsv',
        'selection':root/'hotspots/results/selection.window_selection_summary.tsv',
        'assoc_report':root/'assoc/reports/index.html',
    }
    before={key:identity(path) for key,path in files.items()}
    data={key:read_table(path) for key,path in files.items() if key!='assoc_report'}
    if not data['occupancy'] and not data['nominated']:
        raise ValueError('No scan occupancy or nominated window results in this snapshot')
    primary=sorted([r for r in data['nominated'] if r.get('is_primary')=='1'], key=lambda r:-float(r['oe']))
    counts=Counter(r['threshold'] for r in data['nominated'])
    atlas={}
    for row in data['atlas']:
        if row['metric']=='AVI_PHRED' and row['allele_scope']=='all_three_possible':
            if row['window_id'] in atlas:
                raise ValueError('Duplicate window AVI summary')
            atlas[row['window_id']]=row
    controls=defaultdict(list); unmatched=set()
    for row in data['matches']:
        if row['status']=='unmatched': unmatched.add(row['case_window_id'])
        elif row['status'] in {'matched','partial_match'}:
            controls[row['case_window_id']].append(row['control_window_id'])
    contrasts=[]
    for case, members in controls.items():
        if len(set(members))!=len(members): raise ValueError('Duplicate matched control')
        usable=[float(atlas[c]['mean']) for c in members if c in atlas]
        if case not in atlas or len(usable)!=len(members): continue
        control=sum(usable)/len(usable); case_mean=float(atlas[case]['mean'])
        contrasts.append(dict(case_window=case,n_controls=len(usable),case_avi_mean=case_mean,
                              matched_control_avi_mean=control,difference=case_mean-control,
                              inference='descriptive_selected_multialt_targets_not_independent_selection_test'))
    singleton=next((r for r in data['singletons'] if r.get('chrom')=='autosomes_observed'),{})
    status=[dict(module='gnomAD scan',status='table_available' if data['occupancy'] else 'not_in_snapshot',
                 evidence='Nonreference ALT occupancy; not joint REF/ALT frequency ranks'),
            dict(module='hotspots',status=f"{len(primary)} primary windows; {len(data['nominated'])} AC-tier records",
                 evidence='Working-model nominations; not independent replication'),
            dict(module='HERV/regulatory contrasts',status=f"{len(data['annotations'])} contrast rows",
                 evidence='Zero rows means no reported enrichment result'),
            dict(module='HERV personal census',status='not_summarized_by_these_tables',
                 evidence='Requires individual genome_fraction.json and denominator audit'),
            dict(module='ancient trajectories',status='not_summarized_by_these_tables',
                 evidence='Published selection lookups are not individual trajectories')]
    assoc_text=files['assoc_report'].read_text() if files['assoc_report'].is_file() else ''
    if re.search(r'incomplete legacy run|not final results',assoc_text,re.I):
        assoc_status='incomplete_legacy_run'
    elif assoc_text:
        assoc_status='report_present_completion_not_inferred'
    else: assoc_status='not_in_snapshot'
    status.append(dict(module='UKB association',status=assoc_status,evidence='Imputed GWAS is distinct from a native WGS rank3/4 burden test'))
    facts=dict(nominated_by_ac=dict(counts),primary_windows=len(primary),
               primary_strict_2x=sum(r.get('strict_2x_supported')=='1' for r in primary),
               primary_masks= dict(Counter(r.get('mask_status','NA') for r in primary)),
               selected_window_rows=len(data['selected']),atlas_window_rows=len(atlas),
               matched_avi_cases=len(contrasts),unmatched_cases=len(unmatched),
               matched_avi_negative=sum(r['difference']<0 for r in contrasts),
               matched_avi_mean_difference=sum(r['difference'] for r in contrasts)/len(contrasts) if contrasts else None,
               annotation_contrast_rows=len(data['annotations']),singleton_sensitivity=singleton,
               association_status=assoc_status)
    out=Path(args.out)
    for folder in ('results','reports','qc'): (out/folder).mkdir(parents=True,exist_ok=True)
    write_table(out/'results/module_status.tsv',status,['module','status','evidence'])
    write_table(out/'results/matched_avi_descriptive.tsv',contrasts,['case_window','n_controls','case_avi_mean','matched_control_avi_mean','difference','inference'])
    winfields=['chrom','start0','end0','threshold','n_eligible','observed_multi','oe','q_value','strict_2x_supported','mask_status']
    write_table(out/'results/primary_hotspots.tsv',primary,winfields)
    display=[]
    for r in data['occupancy']:
        display.append(dict(ac=r['threshold'],eligible=f"{int(r['n_sites']):,}",multi=f"{int(r['n_multi']):,}",
                            three=f"{int(r['n_three']):,}",percent=f"{100*float(r['multi_fraction']):.3f}%"))
    peaks=[dict(region=f"{r['chrom']}:{int(r['start0']):,}–{int(r['end0']):,}",oe=f"{float(r['oe']):.3f}",
                q=f"{float(r['q_value']):.3g}",mask=r.get('mask_status','NA')) for r in primary[:10]]
    parts=['<!doctype html><meta charset="utf-8"><title>allele3 result review</title>',
           '<style>body{font:16px/1.6 Arial,"Noto Sans CJK SC",sans-serif;max-width:1100px;margin:40px auto;color:#1c3447;padding:0 25px}h1,h2{color:#173b58}table{border-collapse:collapse;width:100%;font-size:14px}th,td{border-bottom:1px solid #dbe3e8;text-align:left;padding:9px}th{background:#edf3f6}aside{background:#fff5e8;padding:18px;border-left:4px solid #d79531}code{background:#edf3f6;padding:2px 5px}</style>',
           '<h1>allele3：当前结果与可检验假设</h1>',
           '<p>这是对现有结果快照的重新汇总；没有运行全基因组扫描、个体 HERV 检测或 UKB 关联。坐标区间为 GRCh38、0-based half-open。</p>',
           f'<p>代码快照：<code>{html.escape(args.scripts_revision)}</code>；结果快照：<code>{html.escape(args.analysis_revision)}</code>。</p>',
           '<h2>1. 已完成与未完成</h2>',html_table(status,[('module','模块'),('status','状态'),('evidence','解释')]),
           '<h2>2. AC 阈值敏感性</h2>',
           '<p>分母是至少一个非 REF 碱基达到 AC 阈值的已观察 SNV 位点；multi 指至少两个非 REF 碱基达到阈值，three 指三个非 REF 碱基。不能据此断言 REF 仍被观察或把 ALT 次序称为 R1–R4。</p>',
           html_table(display,[('ac','AC 阈值'),('eligible','分母位点'),('multi','≥2 非 REF'),('three','3 非 REF'),('percent','≥2 非 REF 比例')]),
           '<h2>3. 主分析候选热点</h2>',
           f'<p>AC≥5 主分析 {len(primary)} 个；所有 AC 分层合计 {len(data["nominated"])} 条提名记录。同一数据的不同 AC 阈值不是独立复现。主分析严格“下界超过 2 倍”标记有 {facts["primary_strict_2x"]} 个。</p>',
           html_table(peaks,[('region','区间'),('oe','观察/预期'),('q','原报告 q'),('mask','技术掩膜')]),
           '<aside>极低 q 值仍依赖工作模型和输入质量。未应用技术掩膜时，应先检查低可比对、重复序列、片段重复、着丝粒邻近等问题，再讨论生物学富集。</aside>',
           '<h2>4. AlphaGenome 描述性线索</h2>',
           f'<p>有 {len(contrasts)} 个带匹配对照的 case 窗口，另有 {len(unmatched)} 个未匹配 case；{facts["matched_avi_negative"]}/{len(contrasts)} 组的 case 平均 AVI-PHRED 低于其对照。每组 case−对照均值后，再跨组平均为 {facts["matched_avi_mean_difference"] if contrasts else "NA"}。</p>',
           '<p>这一结果来自每窗口有限抽样的 multi-ALT 位点，只能作为待验证的功能容忍度线索。AVI 用 gnomAD 频率标签训练，并含保守性信息；不能作为多态性或选择的独立验证。18 维 SHAP 贡献也不是 18 种已测生物功能。</p>',
           '<h2>5. 科学问题的四个层次</h2><ol><li>产生：独立突变、复发突变或基因流何时带来某个碱基？</li><li>被观察：样本量、覆盖度和分型能力能否检测到它？</li><li>维持：漂变、净适应度、显隐性及频率依赖是否维持多态？</li><li>后果：分子功能与疾病关联是否在独立数据中得到支持？</li></ol>',
           '<p>“功能多”不是稳定平衡的充分条件。新版动态模拟分别比较固定净优势、正/负频率依赖、杂合优势与环境权衡；不能把模拟参数称为真实人类选择系数。</p>',
           '<h2>参考来源</h2><ul><li><a href="https://doi.org/10.64898/2026.09.16.26363192">AlphaGenome Atlas（预印本）</a></li><li><a href="https://doi.org/10.1038/s41586-026-10358-1">Akbari et al., Nature 2026</a></li><li><a href="https://doi.org/10.1534/genetics.109.110510">Human Triallelic Sites, Genetics 2010</a></li><li><a href="https://doi.org/10.1093/molbev/msag188">Beyond fixation, MBE 2026</a></li></ul>',
           '<p>结果数据来源：<a href="https://github.com/jielab/analysis/tree/'+html.escape(args.analysis_revision)+'/allele3">固定 GitHub 结果快照</a>。</p>']
    (out/'reports/study_review.html').write_text('\n'.join(parts),encoding='utf-8')
    if any(identity(files[k])!=v for k,v in before.items()): raise ValueError('Result snapshot changed during summary')
    manifest=dict(status='complete',settings=vars(args),source_files=before,facts=facts,
                  interpretation='Descriptive audit of supplied snapshot; no genomic rerun or confirmation of biological hypotheses')
    (out/'qc/manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n')
    return manifest


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--analysis-dir',required=True,help='Existing allele3 result directory')
    p.add_argument('--scripts-revision',default='unspecified')
    p.add_argument('--analysis-revision',default='main')
    p.add_argument('--out',required=True)
    return p


def main(argv=None):
    result=run(parser().parse_args(argv));print(json.dumps(result['facts'],ensure_ascii=False));return 0


if __name__=='__main__':raise SystemExit(main())
