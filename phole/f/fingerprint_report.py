#!/usr/bin/env python3
"""Render an offline review of frozen plans or executed regional GWAS checks.

This renderer does not turn descriptive comparisons into significance tests,
colocalization, independent validation, or an assay capture-identity verdict.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile


TABLES = {
    "hypotheses": "hypotheses.tsv",
    "assays": "assay-plan.tsv",
    "regions": "regions.tsv",
    "queries": "lead-queries.tsv",
    "summary": "pair-summary.tsv",
    "coverage": "region-assay-qc.tsv",
    "cis": "cis-window-checks.tsv",
    "cis_plan": "cis-windows.tsv",
    "comparisons": "region-comparisons.tsv",
    "controls": "matched-reference-controls.tsv",
    "leads": "lead-effects.tsv",
    "consistency": "source-summary-checks.tsv",
    "family": "frozen-followup-family.tsv",
    "exclusions": "plan-exclusions.tsv",
}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_text(path, content):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix="." + path.name, dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            stream.write(content)
        os.replace(name, path)
    finally:
        Path(name).unlink(missing_ok=True)


def load_bundle(folder):
    folder = Path(folder)
    receipt = folder / "fingerprint-check.summary.json"
    info = json.loads(receipt.read_text())
    if info.get("mode") not in {"plan_only", "full_regions"}:
        raise ValueError("Unrecognized regional follow-up mode")
    hashes = info.get("output_sha256", {})
    data = {"run": info, "tables": {}, "columns": {}, "files": {}, "plot": None}
    inputs = {receipt.name: sha(receipt)}
    for key, name in TABLES.items():
        path = folder / name
        if info["mode"] == "plan_only" and key not in {"hypotheses", "assays", "regions", "queries", "family", "exclusions", "cis_plan"}:
            continue  # Never display old execution artifacts after a plan rerun.
        if not path.is_file():
            if key in {"hypotheses", "assays", "regions"}:
                raise ValueError("Missing frozen-plan input: " + str(path))
            continue
        actual = sha(path)
        expected = hashes.get(name)
        if isinstance(expected, dict):
            expected = expected.get("sha256")
        if not expected:
            if key in {"hypotheses", "assays", "regions"}:
                raise ValueError("Frozen-plan input lacks a current checksum: " + name)
            continue
        if actual != expected:
            raise ValueError("Stale or modified regional result: " + str(path))
        with path.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            data["tables"][key] = list(reader)
            data["columns"][key] = reader.fieldnames or []
        data["files"][key] = name
        inputs[name] = actual
    plot = folder / "plot-data.json"
    if info["mode"] == "full_regions" and plot.is_file() and hashes.get(plot.name):
        actual = sha(plot)
        expected = hashes.get(plot.name)
        if isinstance(expected, dict):
            expected = expected.get("sha256")
        if expected and actual != expected:
            raise ValueError("Stale regional plot data")
        data["plot"] = json.loads(plot.read_text())
        inputs[plot.name] = actual
    return data, inputs


HTML = r'''<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>PHOLE · 完整区域 GWAS 复核</title>
<style>
:root{font:15px/1.65 system-ui,sans-serif;color:#172c40;background:#f2f5f8;--blue:#116783;--muted:#56697c;--line:#dce4eb}
*{box-sizing:border-box}body{margin:0}header{padding:30px max(24px,calc((100vw - 1400px)/2));background:#123d54;color:white}h1{margin:0;font-size:29px}header p{max-width:1000px;color:#d3e5ee;margin:8px 0 0}.wrap{max-width:1448px;margin:auto;padding:24px}section{background:#fff;border:1px solid var(--line);border-radius:10px;padding:22px;margin-bottom:18px}h2{font-size:21px;margin:0 0 12px}h3{font-size:17px}.note{background:#edf4f8;border-left:4px solid var(--blue);padding:12px 16px}.notice{background:#fff5dd;border-left:4px solid #c5870d}.synthetic{background:#ffe9e9;border-left:5px solid #b32434;font-weight:700}.muted{color:var(--muted)}.metrics{display:flex;flex-wrap:wrap;gap:12px;margin:18px 0}.metric{background:#fff;border:1px solid var(--line);border-radius:8px;padding:12px 20px;min-width:155px}.metric b{display:block;font-size:27px;color:var(--blue)}.toolbar{display:flex;gap:14px;flex-wrap:wrap;align-items:end;margin:14px 0}.toolbar label{display:flex;flex-direction:column;gap:5px;font-weight:600}.toolbar label:first-child{flex:2;min-width:280px}select,input,button{font:inherit;border:1px solid #b8c6d1;border-radius:5px;padding:7px 10px;background:white;color:#172c40}select{max-width:100%}button{cursor:pointer}button:disabled{opacity:.4;cursor:default}a{color:#086285}.table-wrap{overflow:auto;max-height:500px;border:1px solid var(--line);border-radius:4px}table{border-collapse:collapse;width:100%;font:13px/1.5 ui-monospace,SFMono-Regular,Consolas,monospace}th,td{border-bottom:1px solid #e3e9ef;padding:8px 10px;text-align:left;white-space:nowrap;vertical-align:top}th{position:sticky;top:0;background:#eaf0f5;z-index:1}tr:nth-child(even){background:#f8fafc}.pages{display:flex;gap:12px;align-items:center;margin:10px 0}.empty{padding:16px;color:var(--muted)}.legend{display:flex;gap:24px;flex-wrap:wrap}.chip{padding:2px 8px;background:#e7eff4;border-radius:20px;font-size:13px}details{margin:14px 0}summary{cursor:pointer;font-weight:600}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f5f7fa;padding:14px;max-height:450px;overflow:auto;font-size:12px}#plot svg{width:100%;height:auto;max-height:1600px}.plot-links a{display:inline-block;margin:4px 15px 4px 0}code{background:#edf2f6;padding:2px 4px;border-radius:3px}footer{color:var(--muted);font-size:13px;padding:8px 2px} @media(max-width:650px){.wrap{padding:12px}section{padding:14px}h1{font-size:23px}.toolbar{display:block}.toolbar label{margin:10px 0}}
</style></head><body><header><h1>完整区域 GWAS 复核</h1><p>把候选 assay 与同名靶标、替代靶标放在同一组区域中比较。保留不显著、缺测、无法定向与原始 P 值。</p></header>
<main class="wrap"><div id="mode" class="note"></div><div id="synthetic"></div><div id="metrics" class="metrics"></div>
<section><h2>01 · 冻结的病例与数据覆盖</h2><p class="muted">病例来自已有摘要筛查。完整区域的再次读取属于跟进分析；被用于提名的队列和位点不会被计作新的独立验证。</p>
<div class="toolbar"><label>病例<select id="case"></select></label><label>区域<select id="region"></select></label></div><div id="case-description" class="note"></div>
<h3>预期 assay 与实际输入状态</h3><div id="assay-table"></div><p class="muted">未在映射表测量、缺少完整文件、build 不匹配、无可用变异、有值但不显著，是不同的状态。</p></section>
<section><h2>02 · 区域信号图</h2><p class="muted">横轴为同一已声明版本区域的位置；纵轴为来源 P 的 −log₁₀(P)。极小 P 按原文保留；源 P=0 且无有效 LOG10P 时不赋予定量高度，在原始结果表中保留，不能猜测其精确强度。只对显示采样，完整区域变异在 extracted/ 中。</p><div id="plot"></div><div id="plot-note" class="muted"></div></section>
<section><h2>03 · 完整结果与对照</h2><div class="note">Z / P profile 相关、同号比例和 lead 重现是描述性证据。区域 SNP 之间存在 LD，同队列测量共享样本；这些指标不是独立 SNP 检验、正式共定位、或已确认错捕的概率。跨平台 BETA 的绝对大小可能使用不同表型尺度。</div>
<div class="toolbar"><label>结果表<select id="table"></select></label><label>筛选表内文本<input id="search" type="search" placeholder="例如 rs7895470、ckb、unresolved"></label></div><p id="table-note" class="muted"></p><div id="result-table"></div><div class="plot-links" id="downloads"></div></section>
<section><h2>04 · 解释边界与来源</h2><p>完整数据可以核查自身 cis 是否遗漏、对齐 EA/NEA 后的方向是否相符、以及替代靶标在相同变异覆盖下是否更接近候选。即使这些检查通过，也仍需区分直接错捕、混合捕获、真实 trans 调控、蛋白复合物、表位效应和数据标签问题。</p><p>若要作正式多信号共定位或条件分析，需要与人群匹配的 LD 和适合的统计模型。此报告不会用区域相关系数替代这些步骤。</p><details><summary>运行来源记录</summary><pre id="provenance"></pre></details><p><a href="fingerprint-check.summary.json">原始运行汇总 JSON</a> · <a href="fingerprint-report.receipt.json">本报告来源与哈希</a></p></section>
<footer>PHOLE v5 · 原始表可完整下载；筛选只影响当前视图。此页面可离线打开。</footer></main>
<script id="phole-data" type="application/json">__DATA__</script>
<script>
"use strict";
const D=JSON.parse(document.getElementById('phole-data').textContent),$=id=>document.getElementById(id);
const H=D.tables.hypotheses||[],R=D.tables.regions||[];
const labels={summary:'病例汇总',cis:'精确编码基因 cis 窗口覆盖与最小 P',cis_plan:'冻结的精确 cis 窗口',coverage:'合并查询区域逐 assay 覆盖与最小 P',controls:'同一参考队列与相同变异的目标对照',comparisons:'区域间有方向的 profile 比较',leads:'冻结 lead 查询的完整 GWAS 效应',consistency:'原摘要与完整数据的源一致性',queries:'冻结的 lead 并集',regions:'冻结的区域',assays:'预期 assay',hypotheses:'提名病例',family:'冻结的探索性查询家族',exclusions:'计划中无法评价的项目'};
function value(obj,...keys){for(const k of keys){if(obj&&obj[k]!==undefined&&obj[k]!==null&&obj[k]!=='.')return obj[k];}return '';}
function text(el,s){el.textContent=String(s);return el;}function el(tag,s){let x=document.createElement(tag);if(s!==undefined)text(x,s);return x;}
function hid(r){return value(r,'HYPOTHESIS_ID','hypothesis_id');}function rid(r){return value(r,'REGION_ID','region_id');}
function number(v){if(v===undefined||v===null||v===''||v==='.'||typeof v==='boolean')return null;const n=Number(v);return Number.isFinite(n)?n:null;}
function fmt(v){return v===null||v===undefined?'未提供':String(v);}
const isPlan=D.run.mode==='plan_only';
text($('mode'),isPlan?'计划状态：这里只读取了现有摘要并冻结复核区域，没有读取完整 GWAS。下方覆盖状态不能作为阳性或阴性分析结果。':'执行状态：已按冻结计划读取可用的完整 GWAS 区域。请结合逐文件状态判断覆盖；完成运行不等于全部 assay 均有可用数据。');
$('mode').classList.toggle('notice',isPlan);
if(D.run.synthetic_input===true||D.run.synthetic_input===1||String(D.run.synthetic_input).toLowerCase()==='true'){$('synthetic').className='note synthetic';text($('synthetic'),'SYNTHETIC / 合成测试数据：本页仅验证程序，不是这些蛋白的真实 GWAS 结果。');}
for(const [name,n] of [['提名病例',H.length],['规范蛋白',D.run.canonical_proteins],['映射 assay',D.run.mapped_assays],['已检查完整文件',D.run.full_files_inspected||0],['已确认捕获身份',0]]){const x=el('div');x.className='metric';x.append(el('b',fmt(n)),el('span',name));$('metrics').append(x);}
for(const h of H){const o=el('option',`${value(h,'PRIMARY','PRIMARY_DATASET')} ${value(h,'PROT')} → ${value(h,'ALTERNATIVE_PROT')} (${value(h,'REFERENCE','REFERENCE_DATASET')})`);o.value=hid(h);$('case').append(o);}
for(const k of Object.keys(labels)){if(D.tables[k]){const o=el('option',labels[k]);o.value=k;$('table').append(o);}}
text($('provenance'),JSON.stringify(D.run,null,2));
function filtered(rows,withRegion=true){const h=$('case').value,r=$('region').value;return rows.filter(x=>(!h||!hid(x)||hid(x)===h)&&(!withRegion||!r||!rid(x)||rid(x)===r));}
function table(root,rows,columns){root.replaceChildren();if(!rows.length){root.append(el('div',isPlan?'尚无执行结果，或当前筛选无计划条目。':'当前筛选没有结果；这不是零效应。'));root.firstChild.className='empty';return;}
 const size=75;let page=0;const pages=el('div');pages.className='pages';const prev=el('button','上一页'),next=el('button','下一页'),status=el('span');const wrap=el('div');wrap.className='table-wrap';pages.append(prev,status,next);root.append(pages,wrap);
 function draw(){wrap.replaceChildren();const t=el('table'),thead=el('thead'),tr=el('tr');for(const c of columns)tr.append(el('th',c));thead.append(tr);t.append(thead);const b=el('tbody');for(const row of rows.slice(page*size,(page+1)*size)){const line=el('tr');for(const c of columns)line.append(el('td',row[c]??'.'));b.append(line);}t.append(b);wrap.append(t);text(status,`${page*size+1}–${Math.min(rows.length,(page+1)*size)} / ${rows.length}`);prev.disabled=page===0;next.disabled=(page+1)*size>=rows.length;}
 prev.onclick=()=>{page--;draw();};next.onclick=()=>{page++;draw();};draw();
}
function result(){const key=$('table').value,query=$('search').value.toLowerCase();const rows=filtered(D.tables[key]||[],!['cis','cis_plan'].includes(key)).filter(r=>!query||Object.values(r).some(v=>String(v).toLowerCase().includes(query)));table($('result-table'),rows,D.columns[key]||[]);text($('table-note'),`${labels[key]||key}：当前筛选 ${rows.length} 行；下载链接提供未经当前筛选删减的完整表。${['cis','cis_plan'].includes(key)?'精确 cis 表按病例展示全部编码基因窗口，不以合并查询区域代替 cis 边界。':''}`);}
for(const [key,name] of Object.entries(D.files)){const a=el('a',`${labels[key]||key} TSV`);a.href=name;a.download=name;$('downloads').append(a);}
// Plot-data is a collection of tracks with raw-P strings retained on points.
function tracks(){let p=D.plot;if(!p)return [];if(Array.isArray(p))return p;if(Array.isArray(p.tracks))return p.tracks;if(Array.isArray(p.regions))return p.regions.flatMap(r=>(r.assays||r.tracks||[]).map(t=>({...r,...t,tracks:undefined,assays:undefined})));return [];}
function neglog(p){const s=String(p??'').trim();if(!s||s==='.')return null;const m=s.match(/^([0-9]*\.?[0-9]+)(?:[eE]([+-]?\d+))?$/);if(!m)return null;const a=Number(m[1]);if(a===0)return null;const n=-Math.log10(a)-(m[2]?Number(m[2]):0);return Number.isFinite(n)&&n>=0?{n,bound:false}:null;}
function pointInfo(p){const pos=number(value(p,'POS','pos')),source=number(value(p,'LOG10P','log10p','source_log10p'));const raw=value(p,'P','p_text','p');const n=source!==null&&source>=0?{n:source,bound:false}:neglog(raw);return pos!==null&&n?{p,pos,y:Math.min(n.n,350),bound:n.bound,capped:n.n>350,raw}:null;}
const NS='http://www.w3.org/2000/svg';function sv(tag,attrs,content){const n=document.createElementNS(NS,tag);for(const [k,v] of Object.entries(attrs||{}))n.setAttribute(k,v);if(content!==undefined)n.textContent=String(content);return n;}
function plot(){const host=$('plot');host.replaceChildren();$('plot-note').textContent='';if(isPlan){host.append(el('p','完整 GWAS 尚未读取。运行 fingerprint-check 后，这里显示区域内全部可用变异的概览。'));return;}
 const chosen=tracks().filter(t=>hid(t)===$('case').value&&(!$('region').value||rid(t)===$('region').value));if(!chosen.length){host.append(el('p','本区域没有可绘制的数据。请先查看区域覆盖状态；缺图不代表无关联。'));return;}
 const region=R.find(r=>hid(r)===$('case').value&&rid(r)===$('region').value)||{},left=250,right=1080,width=1120,step=118;
 const groups=chosen.map(t=>({t,points:(t.points||t.variants||[]).map(pointInfo).filter(Boolean)}));let positions=groups.flatMap(g=>g.points.map(p=>p.pos));let start=number(region.START)??Math.min(...positions),end=number(region.END)??Math.max(...positions);if(!Number.isFinite(start)||!Number.isFinite(end)){host.append(el('p','当前区域没有可定量绘制的坐标/P 记录；源 P=0 无 LOG10P 的记录仍保留在原始表中。'));return;}if(end<=start)end=start+1;
 let ymax=Math.max(10,...groups.flatMap(g=>g.points.map(p=>p.y)));const svg=sv('svg',{viewBox:`0 0 ${width} ${groups.length*step+65}`,role:'img','font-family':'Arial, Noto Sans CJK SC, sans-serif','aria-label':'按 assay 排列的完整区域关联概览'});
 const stateLabels={not_mapped:'未测量',manifest_assay_missing:'缺完整文件记录',full_file_missing:'完整文件不存在',build_mismatch:'坐标版本不符',manifest_build_unresolved:'版本未确定',region_without_variants:'区域无可用变异',region_rows_all_failed_QC:'区域记录均未通过 QC',observed:'已读到区域数据'};
 groups.forEach(({t,points},i)=>{const y0=i*step+32,base=y0+70,state=value(t,'STATUS','status'),label=`${value(t,'DATASET','dataset')} · ${value(t,'ASSAY','assay')||'未映射 assay'}`;const name=sv('text',{x:6,y:y0+8,'font-size':14,fill:'#123d54'},label.length>27?label.slice(0,25)+'…':label);name.append(sv('title',{},label));svg.append(name);svg.append(sv('text',{x:6,y:y0+29,'font-size':11,fill:'#617386'},`${value(t,'PROT','prot')} | 定量显示 ${points.length}`));const statusText=stateLabels[state]||state||'状态未提供';svg.append(sv('text',{x:6,y:y0+47,'font-size':11,fill:state==='observed'?'#617386':'#9b5a18'},statusText));
 if(state&&state!=='observed'){svg.append(sv('rect',{x:left,y:y0,width:right-left,height:70,fill:'#f4f6f8'}));svg.append(sv('text',{x:(left+right)/2,y:y0+40,'text-anchor':'middle','font-size':13,fill:'#7d654b'},statusText+'；不能据此判断无关联'));return;}
 if(!points.length){svg.append(sv('text',{x:(left+right)/2,y:y0+40,'text-anchor':'middle','font-size':13,fill:'#7d654b'},'无可定量绘制的 P / LOG10P；请查看原始记录'));return;}
 svg.append(sv('line',{x1:left,y1:base,x2:right,y2:base,stroke:'#b5c3ce'}));for(const tick of [0,ymax]){const y=base-tick/ymax*68;svg.append(sv('text',{x:left-8,y:y+4,'text-anchor':'end','font-size':10,fill:'#52677a'},tick.toFixed(tick===0?0:1)));}
 const threshold=-Math.log10(5e-8);svg.append(sv('line',{x1:left,y1:base-threshold/ymax*68,x2:right,y2:base-threshold/ymax*68,stroke:'#b5c3ce','stroke-dasharray':'4 4'}));
 for(const q of points){const color=String(value(t,'ROLE','role')).includes('alternative')?'#b95030':'#137b9c',dot=sv('circle',{cx:left+(q.pos-start)/(end-start)*(right-left),cy:base-q.y/ymax*68,r:q.bound?2.7:2,fill:color,'fill-opacity':.7});dot.append(sv('title',{},`${value(q.p,'SNP','snp')} | ${q.pos}\nP=${q.raw}\nEA/NEA=${value(q.p,'EA','ea')}/${value(q.p,'NEA','nea')}\nZ=${value(q.p,'Z','z')}\n${q.bound?'P=0，无精确 LOG10P；仅下界':q.capped?'纵轴截于350；原值保留在表中':''}`));svg.append(dot);}}
 );const bottom=groups.length*step+12;for(let i=0;i<5;i++){const pos=start+(end-start)*i/4;svg.append(sv('text',{x:left+(right-left)*i/4,y:bottom,'text-anchor':'middle','font-size':12,fill:'#52677a'},(pos/1e6).toFixed(3)));}svg.append(sv('text',{x:(left+right)/2,y:bottom+25,'text-anchor':'middle','font-size':13},`chr${region.CHR||''} (Mb; build ${region.BUILD||'declared'})`));host.append(svg);text($('plot-note'),'虚线：P=5×10⁻⁸，仅作参照。源 P=0 且无有效 LOG10P 不画定量高度；纵轴显示最高350，来源值仍保留。每条轨道按来源保留；图上 Z 未跨轨道自动定向，正式方向比较请查 lead-effects 和 region-comparisons。采样规则见 plot-data.json，完整区域统计与变异导出不截断。');}
function selection(changeCase){if(changeCase){$('region').replaceChildren();for(const r of R.filter(x=>hid(x)===$('case').value)){const o=el('option',`${r.CHR}:${Number(r.START).toLocaleString()}–${Number(r.END).toLocaleString()} · ${r.REGION_ROLES||''}`);o.value=rid(r);$('region').append(o);}}
 const h=H.find(x=>hid(x)===$('case').value)||{};text($('case-description'),`${value(h,'PRIMARY','PRIMARY_DATASET')} ${value(h,'PRIMARY_ASSAY')}（${value(h,'PROT')}） → ${value(h,'ALTERNATIVE_PROT')}；参考 ${value(h,'REFERENCE','REFERENCE_DATASET')}/${value(h,'REFERENCE_ASSAY')}。${isPlan?'当前为冻结计划。':'所有测定身份结论保持未确认。'}`);
 table($('assay-table'),filtered(D.tables.assays||[],false),D.columns.assays||[]);result();plot();}
$('case').onchange=()=>selection(true);$('region').onchange=()=>selection(false);$('table').onchange=result;$('search').oninput=result;selection(true);
window.PHOLE_FULL_REVIEW={data:D,selection,result,plot};
</script></body></html>'''


def render(folder, cross_dir=None):
    folder = Path(folder)
    data, inputs = load_bundle(folder)
    serialized = json.dumps(data, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    serialized = serialized.replace("&", "\\u0026").replace("<", "\\u003c").replace(">", "\\u003e")
    serialized = serialized.replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")
    output = folder / "fingerprint-check.html"
    atomic_text(output, HTML.replace("__DATA__", serialized))
    receipt = {"schema": "phole-fingerprint-report-v1", "mode": data["run"]["mode"],
               "synthetic_input": data["run"].get("synthetic_input"), "input_sha256": inputs,
               "renderer_sha256": sha(__file__), "output": output.name, "output_sha256": sha(output),
               "scientific_results_modified": False}
    atomic_text(folder / "fingerprint-report.receipt.json", json.dumps(receipt, indent=2) + "\n")
    if cross_dir:
        cross = Path(cross_dir)
        # A dedicated navigation page remains valid after any cross-check rerun.
        # Existing scientific TSV/JSON and their integrity receipts are untouched.
        if cross.is_dir():
            relative = Path(os.path.relpath(output, cross)).as_posix()
            import html
            title = "完整区域 GWAS 复核计划" if data["run"]["mode"] == "plan_only" else "完整区域 GWAS 复核"
            link = '<a href="' + html.escape(relative, quote=True) + '">' + title + '</a>'
            nav = '<!doctype html><meta charset="utf-8"><title>' + title + '</title><h1>' + title + '</h1><p>' + link + '</p>'
            nav += '<p>执行状态与缺失信息见报告首页。摘要交叉检查：<a href="cross-check.html">cross-check.html</a>。</p>'
            nav_name = "full-gwas-plan.html" if data["run"]["mode"] == "plan_only" else "full-gwas-review.html"
            atomic_text(cross / nav_name, nav)
    return output


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--check-dir", type=Path, required=True)
    p.add_argument("--cross-dir", type=Path)
    a = p.parse_args(argv)
    print("Regional review: " + str(render(a.check_dir, a.cross_dir)))


if __name__ == "__main__":
    main()
