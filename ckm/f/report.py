#!/usr/bin/env python3
"""Create an offline CKM report from aggregate JSON, using only the stdlib.

The renderer never reads individual-level RDS data. The recursive identifier
check catches obvious accidental individual-level exports; it is not a formal
statistical disclosure-control policy. All tables keep the complete aggregate
input and render only one small page at a time.
"""
from __future__ import annotations

import argparse
import html
import json
import os
from pathlib import Path
import re
import tempfile


FORBIDDEN_KEYS = frozenset({
    "eid", "eids", "iid", "iids", "f_eid", "participant_id",
    "participant_ids", "patient_id", "patient_ids", "individual_id",
    "individual_ids", "individual_records", "participant_records",
    "raw_events", "date_raw", "date_of_birth", "date_birth", "birth_date",
    "dob", "nhs_number", "sample_ids", "subject_id", "subject_ids",
})


def validate_summary(data: object) -> None:
    """Reject malformed aggregate payloads and obvious participant identifiers."""
    if not isinstance(data, dict):
        raise ValueError("Report input must be an aggregate JSON object")

    def visit(value: object, path: str = "summary") -> None:
        if isinstance(value, dict):
            for key, item in value.items():
                normalized = re.sub(r"[^a-z0-9]+", "_", str(key).lower()).strip("_")
                if normalized in FORBIDDEN_KEYS:
                    raise ValueError(f"Individual-level key is not permitted: {path}.{key}")
                visit(item, f"{path}.{key}")
        elif isinstance(value, list):
            for item in value:
                visit(item, path + "[]")

    visit(data)
    # Invalid numbers must be corrected upstream, never silently made into zero.
    json.dumps(data, ensure_ascii=False, allow_nan=False)


TEMPLATE = r'''<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="color-scheme" content="light">
<title>__TITLE__</title>
<style>
:root{--ink:#183349;--muted:#637782;--teal:#087f8c;--teal-soft:#e8f5f4;--navy:#163c59;--bg:#f6f8fa;--line:#dce5eb;--paper:#fff;--gold:#997544;--font:Inter,"Segoe UI","PingFang SC","Microsoft YaHei",Arial,sans-serif}
*{box-sizing:border-box}body{margin:0;color:var(--ink);background:var(--bg);font-family:var(--font);font-size:14px;line-height:1.65;font-variant-numeric:tabular-nums}button,input,select{font:inherit;color:inherit}button,select,input{border:1px solid var(--line);border-radius:7px;background:#fff;padding:8px 11px}button{cursor:pointer}button:hover{border-color:#91b6c0;background:#f0f7f8}button:disabled{opacity:.45;cursor:default}button:focus-visible,input:focus-visible,select:focus-visible,a:focus-visible,summary:focus-visible{outline:3px solid #84c6cf;outline-offset:3px}[hidden]{display:none!important}h1,h2,h3,p{margin:0}h1{font-size:clamp(27px,3vw,39px);font-weight:650;line-height:1.3;letter-spacing:-.7px}h2{font-size:24px;font-weight:620;letter-spacing:-.4px;line-height:1.4}h3{font-size:15px;font-weight:650}a{color:var(--teal)}code{font-family:ui-monospace,SFMono-Regular,Consolas,monospace;font-size:.92em}svg{display:block;width:100%;height:auto;overflow:visible}svg text{font-family:var(--font)}.wrap{max-width:1280px;margin:auto;padding:0 32px}.mast{padding-top:36px;padding-bottom:22px}.eyebrow{display:flex;align-items:center;gap:9px;font-size:11px;letter-spacing:1.7px;font-weight:700;color:var(--teal)}.mark{width:22px;height:22px;border-radius:6px;background:var(--navy);color:white;display:grid;place-items:center;font-size:12px;letter-spacing:0}.heading{display:flex;align-items:flex-end;justify-content:space-between;gap:24px;margin-top:16px}.subtitle{font-size:13px;color:var(--muted);margin-top:10px;max-width:820px}.stamp{font-size:11px;color:var(--muted);text-align:right;white-space:pre-line;min-width:145px}.scope{margin-top:19px;background:#edf5f8;border-left:3px solid var(--teal);padding:11px 15px;border-radius:5px;font-size:12px;color:#3d6273}.demo{margin-top:15px;padding:10px 14px;border:1px solid #e0cc9d;background:#fff8e7;border-radius:7px;color:#735d2b}.cards{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin-bottom:23px}.card{background:white;border:1px solid var(--line);border-radius:10px;padding:18px 20px}.card-label{font-size:11px;color:var(--muted)}.card-value{font-size:30px;line-height:1.4;font-weight:620;color:var(--navy);margin:4px 0;letter-spacing:-.6px}.card-note{font-size:10px;color:var(--muted)}.tabs{display:flex;gap:4px;flex-wrap:wrap;border-bottom:1px solid var(--line);padding-bottom:10px;margin-bottom:22px}.tab{border:0;background:transparent;padding:11px 18px;color:var(--muted);font-size:13px;font-weight:600}.tab span{font-size:10px;margin-right:8px;opacity:.65}.tab[aria-selected=true]{background:var(--teal-soft);color:#14566a}.panel{background:var(--paper);border:1px solid var(--line);border-radius:12px;padding:26px 28px;margin-bottom:22px;min-width:0}.panel-head{display:flex;justify-content:space-between;gap:18px;align-items:flex-start;margin-bottom:20px}.kicker{font-size:10px;color:var(--teal);letter-spacing:1.4px;text-transform:uppercase;font-weight:700;margin-bottom:6px}.lede{font-size:12px;color:var(--muted);margin-top:9px;max-width:900px}.badge{font-size:10px;padding:4px 9px;border:1px solid var(--line);border-radius:5px;background:#f7f9fb;color:var(--muted);white-space:nowrap}.two{display:grid;grid-template-columns:1fr 1fr;gap:22px}.two>.panel{margin-bottom:22px}.control,.table-controls{display:flex;align-items:center;gap:10px;flex-wrap:wrap;margin-bottom:15px;color:var(--muted);font-size:11px}.control label,.table-controls label{display:flex;align-items:center;gap:6px}.control select,.table-controls select{font-size:11px;max-width:250px}.table-search{flex:1;min-width:190px;max-width:450px;font-size:12px}.table-actions{display:flex;gap:8px;align-items:center;margin-left:auto}.table-actions button{font-size:11px}.scroll{overflow:auto;max-width:100%;scrollbar-width:thin}.table-scroll{overflow:auto;max-height:610px;scrollbar-width:thin;border:1px solid var(--line);border-radius:7px}table{border-collapse:separate;border-spacing:0;width:100%;font-size:11px}thead th{background:#edf3f6;position:sticky;top:0;z-index:1;color:#3c5c70;white-space:nowrap;text-align:left;padding:0;border-bottom:1px solid #cfdee5}th button{display:block;width:100%;text-align:left;border:0;background:transparent;border-radius:0;padding:10px 12px;font-size:11px;font-weight:600;white-space:nowrap}td{border-bottom:1px solid #edf1f4;padding:9px 12px;white-space:nowrap;max-width:420px;overflow:hidden;text-overflow:ellipsis}tbody tr:last-child td{border-bottom:0}tbody tr:hover{background:#f3f9fa}td.num{text-align:right}td.status{color:var(--teal);font-size:10px}.pager{display:flex;align-items:center;gap:9px;flex-wrap:wrap;margin-top:13px;font-size:11px;color:var(--muted)}.pager .count{margin-right:auto}.pager button{font-size:11px;padding:5px 10px}.pager select{font-size:11px;padding:4px 8px}.empty{padding:29px 23px;background:#fafcfd;border:1px dashed #d2dfe6;border-radius:8px;color:var(--muted);font-size:12px;text-align:center}.caption{font-size:10px;line-height:1.8;color:var(--muted);margin-top:10px}.notes{color:var(--muted);font-size:12px;line-height:1.9}.notes ul{padding-left:19px;margin:0}.notes li{margin-bottom:5px}.legend{display:flex;flex-wrap:wrap;gap:8px 15px;margin-top:8px;font-size:10px;color:var(--muted)}.dot{display:inline-block;width:9px;height:9px;margin-right:5px;border-radius:3px}.details{border:1px solid var(--line);border-radius:9px;background:white;margin-bottom:22px}.details>summary{padding:15px 21px;cursor:pointer;font-size:13px;font-weight:600}.details>div{padding:5px 21px 21px}.methods-grid{display:grid;grid-template-columns:1fr 1fr;gap:17px}.method{padding:18px 20px;background:#f7fafb;border:1px solid var(--line);border-radius:8px}.method h3{font-size:13px;margin-bottom:7px}.method p{font-size:12px;color:var(--muted)}.section-title{font-size:12px;font-weight:650;margin:20px 0 10px;color:#345a71}.foot{padding:4px 0 33px;font-size:10px;color:var(--muted);display:flex;justify-content:space-between;gap:20px}.loading{padding:28px;color:var(--muted)}.error{border-left:3px solid #b25b48;background:#fff4f0;color:#8a4739;padding:18px;border-radius:7px}.muted{color:var(--muted)}
@media(max-width:900px){.cards{grid-template-columns:1fr 1fr}.two{grid-template-columns:1fr}.methods-grid{grid-template-columns:1fr}.heading{display:block}.stamp{text-align:left;margin-top:13px}.table-actions{margin-left:0}}
@media(max-width:600px){.wrap{padding:0 16px}.mast{padding-top:24px}.panel{padding:20px 17px}.cards{gap:9px}.card{padding:14px}.card-value{font-size:25px}.tabs{gap:2px}.tab{padding:9px 11px;font-size:12px}.badge{display:none}.table-search{max-width:100%;width:100%}.panel-head{gap:8px}.foot{display:block}}
@media print{body{background:white}.wrap{max-width:none;padding:0}.tabs,.control,.table-controls,.pager,.foot{display:none}.tab-panel{display:block!important}.panel,.card{break-inside:avoid}.table-scroll{max-height:none;overflow:visible}.cards{grid-template-columns:repeat(4,1fr)}.details>div{display:block}thead th{position:static}}
</style>
</head>
<body>
<header class="wrap mast">
  <div class="eyebrow"><span class="mark">C</span>CKM RESEARCH</div>
  <div class="heading"><div><h1 id="title">__TITLE__</h1><p class="subtitle">连续分期与阶段差异：从表型构建、风险验证到生活方式和多组学关联。</p></div><p class="stamp" id="stamp"></p></div>
  <p class="scope">研究用连续分期。Stage 内百分位细化当前人群的相对位置；它不是新的临床分级，也不是进入下一 stage 所需的年数。</p>
  <div class="demo" id="demo" hidden><strong>模拟数据示例</strong>：本报告用于检查流程与版式，不代表 UK Biobank 的研究结果。</div>
</header>
<main class="wrap">
  <section class="cards" id="cards" aria-label="研究概览"></section>
  <nav class="tabs" role="tablist" aria-label="报告章节">
    <button class="tab" id="tab-phenotype" role="tab" aria-controls="panel-phenotype" aria-selected="true" tabindex="0" data-panel="phenotype"><span>01</span>分期与表型</button>
    <button class="tab" id="tab-clock" role="tab" aria-controls="panel-clock" aria-selected="false" tabindex="-1" data-panel="clock"><span>02</span>Staging clock</button>
    <button class="tab" id="tab-factors" role="tab" aria-controls="panel-factors" aria-selected="false" tabindex="-1" data-panel="factors"><span>03</span>Shared &amp; distinct</button>
    <button class="tab" id="tab-yang" role="tab" aria-controls="panel-yang" aria-selected="false" tabindex="-1" data-panel="yang"><span>04</span>Stage 4 · Yang</button>
    <button class="tab" id="tab-methods" role="tab" aria-controls="panel-methods" aria-selected="false" tabindex="-1" data-panel="methods"><span>05</span>方法与运行信息</button>
  </nav>
  <section class="tab-panel" id="panel-phenotype" role="tabpanel" aria-labelledby="tab-phenotype">
    <div class="panel"><div class="panel-head"><div><p class="kicker">Baseline phenotype</p><h2>分期构成</h2><p class="lede">分期、事件和死亡均来自输入的汇总数据。未测量的亚临床疾病不能自动视为阴性。</p></div><span class="badge" id="phase"></span></div><div class="control"><label>显示指标 <select id="stage-measure"><option value="n">基线人数</option><option value="events">随访事件数</option><option value="deaths">死亡数</option></select></label></div><div id="stage-bars"></div><div class="caption" id="stage-caption"></div><div id="stage-counts" style="margin-top:22px"></div></div>
    <div class="two"><div class="panel"><h3>信息完整性</h3><p class="lede" style="margin-bottom:16px">保留输入定义的完整性口径和分母。</p><div id="completeness"></div></div><div class="panel"><h3>关键变量缺失</h3><p class="lede" style="margin-bottom:16px">缺失不等于正常，不在图中补零。</p><div id="missingness"></div></div></div>
    <div class="panel"><h3>两种 stage 3 定义的对应关系</h3><p class="lede" style="margin-bottom:17px">SCORE2 为欧洲人群适用的改编定义；PREVENT base 的风险定义与预测结局不同。比较表显示实际重分类人数。</p><div id="stage-comparison"></div></div>
    <div class="panel" id="phenotype-notes-panel" hidden><h3 style="margin-bottom:12px">表型说明</h3><div id="phenotype-notes" class="notes"></div></div>
  </section>
  <section class="tab-panel" id="panel-clock" role="tabpanel" aria-labelledby="tab-clock" hidden>
    <div class="panel"><div class="panel-head"><div><p class="kicker">Validation</p><h2>连续分期是否提供了更多信息？</h2><p class="lede">先在相同验证人群内比较模型的区分度和校准，再判断 stage 内的连续位置是否有价值。没有完成的分析会保留未运行或失败状态。</p></div><span class="badge">汇总验证结果</span></div><div id="clock-metrics"></div></div>
    <div class="panel"><h3>预测风险与观测风险</h3><p class="lede" style="margin-bottom:16px">每个点对应输入中的一个校准组。虚线表示预测与观测一致；图形只连接同一模型、同一人群、同一 stage 与同一时间窗的组。</p><div class="control" id="calibration-controls"></div><div id="calibration-chart"></div><div class="caption" id="calibration-caption"></div></div>
    <div class="panel"><h3>Stage 内的连续分布</h3><p class="lede" style="margin-bottom:16px">横轴保留输入的分组或百分位标签，纵轴为人数。不同 stage 的百分位不表示相同绝对风险。</p><div class="control" id="distribution-controls"></div><div id="distribution-chart"></div></div>
    <div class="panel"><h3>同一 stage 里的 G1、G2、G3</h3><p class="lede" style="margin-bottom:17px">按该 stage 的训练参照百分位分为三组，查看验证人群的预测与观测风险是否形成梯度。例如 S2.G1 表示 stage 2 内相对位置较低的研究组；它不是经过验证的临床“年级”。</p><div id="clock-grades"></div></div>
    <details class="details"><summary>校准数据表</summary><div id="calibration-table"></div></details>
    <details class="details"><summary>同一验证人群的模型增益</summary><div><p class="lede" style="margin:0 0 15px">Reference 指明比较基准。ΔAUC 大于零、ΔBrier 小于零分别表示相应指标改善；bootstrap 区间以已拟合模型为条件。</p><div id="clock-comparison"></div></div></details>
    <details class="details"><summary>各训练折选择的组学特征</summary><div><p class="lede" style="margin:0 0 15px">特征由各自训练折选择，不能把此表中的选择频率当作独立验证的关联显著性。</p><div id="clock-selection"></div></div></details>
    <details class="details"><summary>拟合状态与可评估人数</summary><div><p class="lede" style="margin:0 0 15px">检查每折训练事件数、参数数、受支持与不受支持的验证人数；失败状态会保留原始原因。</p><div id="clock-fit-status"></div></div></details>
    <div class="panel" id="clock-notes-panel" hidden><h3 style="margin-bottom:12px">Clock 说明</h3><div id="clock-notes" class="notes"></div></div>
  </section>
  <section class="tab-panel" id="panel-factors" role="tabpanel" aria-labelledby="tab-factors" hidden>
    <div class="panel"><div class="panel-head"><div><p class="kicker">Effect modification</p><h2>哪些因素共有，哪些随阶段而变？</h2><p class="lede">阶段差异依据 exposure × stage 的交互检验。一个阶段显著而另一个不显著，不能单独证明差异；关联结果也不能直接证明因果驱动或保护作用。</p></div><span class="badge">Lifestyle · Proteomics · Metabolomics</span></div><div id="pooled-table"></div></div>
    <div class="panel"><h3>各 stage 的效应及不确定性</h3><p class="lede" style="margin-bottom:17px">同时查看 HR、置信区间、事件数和交互作用。可按 source、stage 或 feature 筛选；所有筛选均覆盖完整汇总表。</p><div id="stage-associations"></div></div>
    <details class="details"><summary>共同效应模型中的斜率与分类对比</summary><div><p class="lede" style="margin:0 0 15px">共同斜率模型的 HR 和分类对比。主表中的交互检验与 pattern 判断另行保留。</p><div id="pooled-slopes"></div></div></details>
    <details class="details"><summary>Stage 之间的直接效应对比</summary><div><p class="lede" style="margin:0 0 15px">检查指定阶段间的效应差异、HR 比值及等效性检验。置信区间的水平以列名为准。</p><div id="stage-contrasts"></div></div></details>
    <details class="details"><summary>敏感性分析</summary><div id="sensitivity-table"></div></details>
    <details class="details"><summary>敏感性分析的 stage 分层结果</summary><div id="sensitivity-stage-table"></div></details>
    <details class="details"><summary>横断面的 stage 标志物</summary><div><p class="lede" style="margin:0 0 15px">描述基线各 stage 的分子差异。此表的结局是基线 stage；它与随访临床事件的关联、预测增益分别解释。</p><div id="markers-table"></div></div></details>
    <details class="details"><summary>注释与通路富集</summary><div><p class="lede" style="margin:0 0 15px">按实际检测到的特征与预设背景集合解释富集结果。</p><div id="enrichment-table"></div></div></details>
  </section>
  <section class="tab-panel" id="panel-yang" role="tabpanel" aria-labelledby="tab-yang" hidden>
    <div class="panel"><div class="panel-head"><div><p class="kicker">Baseline stage 4</p><h2>保留已经患病的人群</h2><p class="lede">按 Endpoint 区分基线 stage 4 的横断面患病标志物（OR）与 stage 4 人群的后续结局（HR）。基线患病者不能混入首次发病队列；其生物标志物还可能反映患病后的治疗与行为变化。</p></div><span class="badge">Yang analysis</span></div><div id="yang-table"></div></div>
  </section>
  <section class="tab-panel" id="panel-methods" role="tabpanel" aria-labelledby="tab-methods" hidden>
    <div class="panel"><p class="kicker">Interpretation</p><h2 style="margin-bottom:21px">结果的含义与边界</h2><div class="methods-grid">
      <div class="method"><h3>连续位置与真实时间</h3><p>Stage + stage 内百分位是研究用的位置表达。缺少重复测量时，不能据此估计个体从 stage 1 转到 stage 2 的真实时间。</p></div>
      <div class="method"><h3>独立验证与循环关联</h3><p>训练、预处理和校准应留在训练数据内。用于构建 clock 的同一特征再与 clock 关联，需要作为构成性或描述性关系解释。</p></div>
      <div class="method"><h3>分期定义与预测端点</h3><p>历史 SCORE2 adapted 与 PREVENT base 分开保留。PREVENT 的 ASCVD + HF 端点并不涵盖所有临床 stage 4 组成疾病。</p></div>
      <div class="method"><h3>竞争风险</h3><p>死亡可能阻止首次临床事件发生。绝对风险、校准与 Brier 指标须按模块声明的竞争风险方法解释；原因别 HR 与累计发生风险并非同一个量。</p></div>
      <div class="method"><h3>共有与特异</h3><p>交互检验评估阶段间效应是否不同。未发现交互作用不等于已证明所有阶段的效应相同；稀少事件会限制判断。</p></div>
      <div class="method"><h3>Driver / protector</h3><p>当前观察性扫描识别风险相关或保护相关因素。因果解释仍取决于设计、混杂控制与额外证据，不能由一张关联表直接得到。</p></div>
    </div></div>
    <div class="panel"><h3 style="margin-bottom:14px">输入与分析说明</h3><div id="input-notes" class="notes"></div></div>
    <div class="panel"><h3 style="margin-bottom:14px">验证与运行记录</h3><div id="validation"></div></div>
  </section>
  <div id="fatal-error" class="error" hidden></div>
  <noscript><div class="error">请在浏览器中启用 JavaScript，以查看此离线报告的图形和汇总表。报告不需要网络连接。</div></noscript>
  <footer class="foot"><span>CKM · 群体汇总 · 自包含离线报告</span><span id="footer-meta"></span></footer>
</main>
<script id="ckm-data" type="application/json">__CKM_JSON__</script>
<script>
'use strict';
(() => {
const $ = id => document.getElementById(id);
const E = x => String(x == null ? '' : x).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const arr = x => Array.isArray(x) ? x : [];
const rows = x => arr(x).filter(r=>r && typeof r==='object' && !Array.isArray(r));
const num = x => x==null || x==='' || typeof x==='boolean' || !Number.isFinite(Number(x)) ? null : Number(x);
const N = x => num(x)===null ? '—' : Number(x).toLocaleString('en-US',{maximumFractionDigits:0});
const F = (x,d=3) => num(x)===null ? '—' : Number(x).toLocaleString('en-US',{maximumFractionDigits:d});
const P = x => num(x)===null ? '—' : F(Number(x)*100,2)+'%';
const sum = xs => xs.reduce((s,v)=>s+(num(v)||0),0);
const valueText = x => x==null ? '' : typeof x==='object' ? JSON.stringify(x) : String(x);
const colors=['#087f8c','#264e73','#c38c51','#718c6f','#96708c','#4e9da7','#7888a3','#bb7169'];
const empty = x => '<div class="empty">'+E(x)+'</div>';
const unique = xs => [...new Set(xs.filter(x=>x!=null && x!=='').map(String))].sort((a,b)=>a.localeCompare(b,undefined,{numeric:true}));
const svg = (w,h,body,label) => '<svg viewBox="0 0 '+w+' '+h+'" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="'+E(label)+'"><title>'+E(label)+'</title>'+body+'</svg>';
const labels={stage:'Stage',source:'Source',feature:'Feature',cohort:'Cohort',model:'Model',reference:'Reference',endpoint:'Endpoint',term:'Contrast',lag_years:'Lag years',n:'N',events:'Events',deaths:'Deaths',HR:'HR',hr:'HR',effect:'Effect',effect_measure:'Effect measure',lower:'95% CI lower',upper:'95% CI upper',ci_low:'95% CI lower',ci_high:'95% CI upper',ci90_low:'90% CI lower',ci90_high:'90% CI upper',p:'P',q:'P-FDR',q_shared:'Shared effect FDR',p_shared:'Shared effect P',q_interaction:'Interaction FDR',p_interaction:'Interaction P',classification:'Classification',pattern:'Pattern',direction:'Direction',horizon:'Years',auc:'AUC',brier:'Brier',cindex:'C-index',status:'Status',predicted:'Predicted risk',observed:'Observed risk',stage_prevent:'Stage · PREVENT',stage_score2:'Stage · SCORE2',missing:'Missing',percent:'Percent',name:'Measure',variable:'Variable',bin:'Bin'};
const priority=['source','feature','cohort','model','reference','stage','grade','term','q_shared','q_interaction','pattern','direction','HR','hr','lower','upper','ci_low','ci_high','effect','effect_measure','n','events','deaths','horizon','auc','brier','cindex','p','q','p_shared','p_interaction','classification','status','endpoint','lag_years','outcome','contrast'];
const orderColumns = rs => {const keys=[...new Set(rs.flatMap(r=>Object.keys(r)))];return [...priority.filter(k=>keys.includes(k)),...keys.filter(k=>!priority.includes(k))];};
const fieldFormat=(k,v)=>{
  if(v==null || v==='')return '—';
  if(typeof v==='boolean')return v?'TRUE':'FALSE';
  if(typeof v==='object')return JSON.stringify(v);
  const n=num(v);
  if(n===null || typeof v==='string' && ['feature','source','model','cohort','stage','outcome','status','classification','variable','name','bin'].includes(k))return String(v);
  if(/(^p$|^q$|pvalue|p_value|p_fdr|^q_|^p_|_p$|_q$|fdr)/i.test(k))return n!==0 && Math.abs(n)<.001?n.toExponential(2):F(n,5);
  if(['predicted','observed','risk','risk10','predicted_risk','observed_risk'].includes(k))return P(n);
  if(['n','events','deaths','missing','count','total','n_total','n_events'].includes(k))return N(n);
  return F(n,4);
};
function noteList(target,x,placeholder='未提供额外说明。') {
  const items=arr(x);target.innerHTML=items.length?'<ul>'+items.map(v=>'<li>'+E(valueText(v))+'</li>').join('')+'</ul>':empty(placeholder);
}
let S;
try {
S=JSON.parse($('ckm-data').textContent);
const ph=S.phenotype||{},cl=S.clock||{},as=S.associations||{};
const mainAssociations=[as.global,as.assoc,as.pooled].find(v=>rows(v).length)||[];
const phenotypeOnly=S.phase==='phenotype_only' || S.mode==='phenotype_only';
const notRun=phenotypeOnly?'当前仅生成表型；本分析尚未运行。':'没有提供此分析的汇总结果；可能尚未运行或不适用，请查看运行记录。';
$('title').textContent=S.title||'CKM 连续分期与阶段差异';
$('stamp').textContent=[S.generated||S.generated_at||'',S.version?'Version '+S.version:''].filter(Boolean).join('\n');
$('footer-meta').textContent=[S.version||'',S.phase||''].filter(Boolean).join(' · ');
$('phase').textContent=phenotypeOnly?'表型已生成':S.phase||'分析汇总';
$('demo').hidden=!(S.demo||/synthetic|合成数据|模拟数据/i.test([S.title||'',S.version||'',...arr(S.input_notes)].join(' ')));
const sc=rows(ph.stage_counts);
const stageKey=x=>String(x==null?'':x).replace(/^stage[ _-]*/i,'').trim();
const early=sc.filter(r=>['0','1','2','3'].includes(stageKey(r.stage)));
const s4=sc.filter(r=>['4','4a','4b'].includes(stageKey(r.stage)));
// Do not add overlapping stage 4 total and stage 4a/4b subtotals.
const s4total=s4.find(r=>stageKey(r.stage)==='4');
const earlyEvents=early.length&&early.every(r=>num(r.events)!==null)?sum(early.map(r=>r.events)):null;
const featureCount=new Set(rows(mainAssociations).map(r=>[r.source||'',r.feature||''].join('\t')).filter(x=>x!=='\t')).size;
const total=num(ph.n_total)!==null?ph.n_total:null;
const cards=[['分析输入人数',total,'来自表型汇总的总人数'],['基线 stage 0–3 随访事件',earlyEvents,'按各 stage 的已提供事件数合计'],['基线 stage 4',s4total?s4total.n:s4.length?sum(s4.map(r=>r.n)):null,'Yang 人群独立分析'],['已汇总的暴露特征',featureCount||null,'按 source 与 feature 去重']];
$('cards').innerHTML=cards.map(([label,v,caption])=>'<div class="card"><div class="card-label">'+E(label)+'</div><div class="card-value">'+N(v)+'</div><div class="card-note">'+E(caption)+'</div></div>').join('');

function table(target,input,{filename='ckm-table',emptyText=notRun,filters=true}={}) {
  const data=rows(input),root=typeof target==='string'?$(target):target;
  if(!data.length){root.innerHTML=empty(emptyText);return;}
  const columns=orderColumns(data),categoryKeys=filters?['source','cohort','model','stage','endpoint','lag_years','reference'].filter(k=>columns.includes(k)&&unique(data.map(r=>r[k])).length>1):[];
  root.innerHTML='<div class="table-controls"><input class="table-search" type="search" placeholder="搜索完整表：特征、通路、状态…" aria-label="搜索完整汇总表"><div class="table-actions"><button type="button" class="download">下载筛选表</button></div></div><div class="table-scroll"><table><thead></thead><tbody></tbody></table></div><div class="pager"><span class="count" aria-live="polite"></span><label>每页 <select aria-label="每页行数"><option>50</option><option>100</option><option>250</option><option>500</option></select></label><button type="button" class="previous" aria-label="上一页">上一页</button><span class="page"></span><button type="button" class="next" aria-label="下一页">下一页</button></div>';
  const controls=root.querySelector('.table-controls'),search=root.querySelector('.table-search'),filterNodes={};
  categoryKeys.forEach(key=>{const label=document.createElement('label');label.textContent=(labels[key]||key)+' ';const sel=document.createElement('select');sel.setAttribute('aria-label','筛选 '+key);sel.innerHTML='<option value="">全部</option>'+unique(data.map(r=>r[key])).map(v=>'<option value="'+E(v)+'">'+E(v)+'</option>').join('');label.append(sel);controls.insertBefore(label,root.querySelector('.table-actions'));filterNodes[key]=sel;});
  const searchText=data.map(r=>columns.map(k=>valueText(r[k])).join('\t').toLocaleLowerCase());
  let page=0,size=50,sortKey=null,direction=1,filtered=[];
  const head=root.querySelector('thead'),body=root.querySelector('tbody');
  head.innerHTML='<tr>'+columns.map(k=>'<th scope="col" aria-sort="none"><button type="button" data-column="'+E(k)+'">'+E(labels[k]||k)+' <span aria-hidden="true">↕</span></button></th>').join('')+'</tr>';
  function refresh() {
    const query=search.value.trim().toLocaleLowerCase();
    filtered=data.map((r,i)=>({r,i})).filter(({r,i})=>(!query||searchText[i].includes(query))&&categoryKeys.every(k=>!filterNodes[k].value||String(r[k])===filterNodes[k].value));
    if(sortKey)filtered.sort((a,b)=>{const av=a.r[sortKey],bv=b.r[sortKey],am=av==null||av==='',bm=bv==null||bv==='';if(am||bm)return am===bm?a.i-b.i:am?1:-1;const an=num(av),bn=num(bv);const c=an!==null&&bn!==null?an-bn:valueText(av).localeCompare(valueText(bv),undefined,{numeric:true});return c===0?a.i-b.i:c*direction;});
    const pages=Math.max(1,Math.ceil(filtered.length/size));page=Math.max(0,Math.min(page,pages-1));
    const visible=filtered.slice(page*size,(page+1)*size);
    body.innerHTML=visible.length?visible.map(({r})=>'<tr>'+columns.map(k=>{const v=r[k],display=fieldFormat(k,v),numeric=num(v)!==null&&!['stage','feature','bin'].includes(k);return '<td'+(numeric?' class="num"':['classification','status'].includes(k)?' class="status"':'')+' title="'+E(valueText(v))+'">'+E(display)+'</td>';}).join('')+'</tr>').join(''):'<tr><td colspan="'+columns.length+'">没有匹配的汇总记录。</td></tr>';
    head.querySelectorAll('th').forEach(th=>{const button=th.querySelector('button'),active=button.dataset.column===sortKey;th.setAttribute('aria-sort',active?(direction===1?'ascending':'descending'):'none');button.querySelector('span').textContent=active?(direction===1?'↑':'↓'):'↕';});
    root.querySelector('.count').textContent='匹配 '+N(filtered.length)+' / '+N(data.length)+' 条 · 检索全部记录';
    root.querySelector('.page').textContent=(page+1)+' / '+pages;
    root.querySelector('.previous').disabled=page===0;root.querySelector('.next').disabled=page>=pages-1;
  }
  head.addEventListener('click',event=>{const b=event.target.closest('button[data-column]');if(!b)return;const k=b.dataset.column;direction=sortKey===k?-direction:1;sortKey=k;page=0;refresh();});
  search.addEventListener('input',()=>{page=0;refresh();});
  Object.values(filterNodes).forEach(sel=>sel.addEventListener('change',()=>{page=0;refresh();}));
  root.querySelector('.pager select').addEventListener('change',event=>{size=Math.min(500,Math.max(1,Number(event.target.value)));page=0;refresh();});
  root.querySelector('.previous').addEventListener('click',()=>{page--;refresh();});root.querySelector('.next').addEventListener('click',()=>{page++;refresh();});
  root.querySelector('.download').addEventListener('click',()=>{
    const cell=x=>{let s=valueText(x);if(typeof x==='string'&&/^[=+@-]/.test(s))s="'"+s;return '"'+s.replace(/"/g,'""')+'"';};
    const content='\ufeff'+columns.map(cell).join('\t')+'\r\n'+filtered.map(({r})=>columns.map(k=>cell(r[k])).join('\t')).join('\r\n');
    const url=URL.createObjectURL(new Blob([content],{type:'text/tab-separated-values;charset=utf-8'}));const a=document.createElement('a');a.href=url;a.download=filename+'.tsv';document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000);
  });
  refresh();
}

function barChart(input,{title='人数分布',key='n',label=r=>r.stage}={}) {
  const data=rows(input).filter(r=>num(r[key])!==null&&num(r[key])>=0);
  if(!data.length)return empty('未提供可绘制的有效汇总数值。');
  const W=1000,H=300,lm=75,rm=25,tm=28,bm=64,pw=W-lm-rm,ph=H-tm-bm,max=Math.max(1,...data.map(r=>num(r[key]))),step=pw/data.length;
  let b='';
  [0,.25,.5,.75,1].forEach(f=>{const y=tm+ph*(1-f);b+='<line x1="'+lm+'" y1="'+y+'" x2="'+(W-rm)+'" y2="'+y+'" stroke="#e1e9ed"/><text x="'+(lm-10)+'" y="'+(y+4)+'" text-anchor="end" fill="#637782" font-size="11">'+N(max*f)+'</text>';});
  data.forEach((r,i)=>{const n=num(r[key]),hh=ph*n/max,w=Math.max(1,step*.66),x=lm+step*(i+.17),lab=String(label(r)),rawstage=num(r.stage),color=colors[(rawstage===null?i:Math.floor(rawstage))%colors.length];b+='<rect x="'+x+'" y="'+(tm+ph-hh)+'" width="'+w+'" height="'+hh+'" rx="3" fill="'+color+'"><title>'+E(lab)+': '+N(n)+'</title></rect>';if(data.length<=18)b+='<text x="'+(x+w/2)+'" y="'+(tm+ph-hh-8)+'" text-anchor="middle" font-size="11" fill="#345a71">'+N(n)+'</text>';if(data.length<=24||i%Math.ceil(data.length/16)===0)b+='<text x="'+(x+w/2)+'" y="'+(H-bm+23)+'" text-anchor="middle" font-size="11" fill="#637782">'+E(lab.length>17?lab.slice(0,16)+'…':lab)+'</text>';});
  return svg(W,H,b,title);
}
function drawStages(){const key=$('stage-measure').value;$('stage-bars').innerHTML=barChart(sc,{title:'CKM stage 分布：'+key,key,label:r=>'Stage '+(r.stage==null?'未知':r.stage)});const unavailable=sc.filter(r=>num(r[key])===null).length;$('stage-caption').textContent=unavailable?'有 '+unavailable+' 个阶段未提供该指标；这些阶段未绘制，未被当成零。':'柱高直接对应各阶段的汇总数值。';}
$('stage-measure').addEventListener('change',drawStages);drawStages();
table('stage-counts',sc,{filename:'ckm-stage-counts',filters:false,emptyText:'没有提供基线分期汇总。'});
table('completeness',ph.completeness,{filename:'ckm-completeness',filters:false,emptyText:'未提供完整性汇总。'});
table('missingness',ph.missingness,{filename:'ckm-missingness',filters:false,emptyText:'未提供变量缺失汇总。'});
table('stage-comparison',ph.stage_comparison,{filename:'ckm-stage-comparison',filters:false,emptyText:'尚未提供两种分期定义的交叉表。'});
table('clock-metrics',cl.metrics,{filename:'ckm-clock-metrics'});
table('calibration-table',cl.calibration,{filename:'ckm-calibration'});
table('clock-comparison',cl.comparison,{filename:'ckm-clock-comparison'});
table('clock-selection',cl.selection,{filename:'ckm-clock-selection'});
table('clock-fit-status',cl.fit_status,{filename:'ckm-clock-fit-status'});
table('clock-grades',cl.grades,{filename:'ckm-clock-grades'});
table('pooled-table',mainAssociations,{filename:'ckm-associations'});
table('pooled-slopes',as.pooled,{filename:'ckm-pooled-slopes'});
table('stage-associations',as.stage,{filename:'ckm-stage-associations'});
table('stage-contrasts',as.contrasts,{filename:'ckm-stage-contrasts'});
table('sensitivity-table',as.sensitivity,{filename:'ckm-sensitivity'});
table('sensitivity-stage-table',as.sensitivity_stage,{filename:'ckm-sensitivity-stage'});
table('markers-table',cl.markers,{filename:'ckm-markers'});
table('yang-table',as.yang,{filename:'ckm-yang'});
table('enrichment-table',Array.isArray(S.enrichment)?S.enrichment:S.enrichment&&S.enrichment.results,{filename:'ckm-enrichment'});

function selectionControls(target,data,keys,onChange,{allowAll=['model'],defaultAll=['model']}={}) {
  const state={},container=$(target);container.innerHTML='';
  keys.forEach(key=>{const vals=unique(data.map(r=>r[key]));if(!vals.length)return;const label=document.createElement('label');label.textContent=(labels[key]||key)+' ';const sel=document.createElement('select');sel.setAttribute('aria-label','图形筛选 '+key);if(allowAll.includes(key))sel.append(new Option('全部',''));vals.forEach(v=>sel.append(new Option(v,v)));const preferred=vals.find(v=>/^(all|overall|0[-–]3|stage.?0[-–]3)$/i.test(v));sel.value=defaultAll.includes(key)&&allowAll.includes(key)?'':preferred||vals[0];state[key]=sel.value;sel.addEventListener('change',()=>{state[key]=sel.value;onChange(state);});label.append(sel);container.append(label);});return state;
}
const cal=rows(cl.calibration);
function calibration(state) {
  const selected=cal.filter(r=>Object.entries(state).every(([k,v])=>!v||String(r[k])===v));
  const data=selected.filter(r=>num(r.predicted)!==null&&num(r.observed)!==null&&num(r.predicted)>=0&&num(r.predicted)<=1&&num(r.observed)>=0&&num(r.observed)<=1);
  if(!data.length){$('calibration-chart').innerHTML=empty(cal.length?'当前筛选没有有效的 0–1 风险数据。':notRun);$('calibration-caption').textContent='';return;}
  const groups=new Map();data.forEach(r=>{const key=JSON.stringify([r.cohort||'',r.model||'',r.stage==null?'':r.stage,r.horizon==null?'':r.horizon]);if(!groups.has(key))groups.set(key,[]);groups.get(key).push(r);});
  const W=950,H=420,lm=82,rm=35,tm=25,bm=70,pw=W-lm-rm,ph=H-tm-bm;
  const maxData=Math.max(...data.flatMap(r=>[num(r.predicted),num(r.observed)]));const top=Math.min(1,Math.max(.05,Math.ceil(maxData*1.12/.05)*.05));
  const x=v=>lm+pw*v/top,y=v=>tm+ph*(1-v/top);let b='';
  for(let i=0;i<=5;i++){const v=top*i/5;b+='<line x1="'+lm+'" y1="'+y(v)+'" x2="'+(W-rm)+'" y2="'+y(v)+'" stroke="#e2e9ef"/><text x="'+(lm-12)+'" y="'+(y(v)+4)+'" text-anchor="end" font-size="11" fill="#637782">'+P(v)+'</text><text x="'+x(v)+'" y="'+(H-bm+24)+'" text-anchor="middle" font-size="11" fill="#637782">'+P(v)+'</text>';}
  b+='<line x1="'+x(0)+'" y1="'+y(0)+'" x2="'+x(top)+'" y2="'+y(top)+'" stroke="#91a5b3" stroke-dasharray="5 5"/><text x="'+(lm+pw/2)+'" y="'+(H-13)+'" text-anchor="middle" fill="#3e6279" font-size="12">预测风险</text><text transform="translate(20 '+(tm+ph/2)+') rotate(-90)" text-anchor="middle" fill="#3e6279" font-size="12">观测风险</text>';
  let legend='';[...groups.entries()].forEach(([key,rs],i)=>{const color=colors[i%colors.length],names=JSON.parse(key),name=names.map((v,j)=>v===''?'':j===2?'stage '+v:j===3?v+' y':v).filter(Boolean).join(' · ');rs.sort((a,b)=>num(a.predicted)-num(b.predicted));if(rs.length>1)b+='<polyline fill="none" stroke="'+color+'" stroke-width="2" points="'+rs.map(r=>x(num(r.predicted))+','+y(num(r.observed))).join(' ')+'"/>';rs.forEach(r=>{b+='<circle cx="'+x(num(r.predicted))+'" cy="'+y(num(r.observed))+'" r="4.3" fill="'+color+'" stroke="white" stroke-width="1.3"><title>'+E(name)+' · bin '+E(r.bin==null?'':r.bin)+' · n='+N(r.n)+' · predicted='+P(r.predicted)+' · observed='+P(r.observed)+'</title></circle>';});legend+='<span><i class="dot" style="background:'+color+'"></i>'+E(name)+'</span>';});
  $('calibration-chart').innerHTML=svg(W,H,b,'预测风险与观测风险校准图')+'<div class="legend">'+legend+'</div>';
  $('calibration-caption').textContent='绘制 '+N(data.length)+' 个汇总组。观测风险的估计方法见运行记录。'+(selected.length>data.length?' 另有 '+N(selected.length-data.length)+' 个缺失或超出 0–1 范围的点未绘制。':'');
}
const calState=selectionControls('calibration-controls',cal,['cohort','stage','horizon','model'],calibration);calibration(calState);
const dist=rows(cl.distribution);
function distribution(state){const data=dist.filter(r=>Object.entries(state).every(([k,v])=>!v||String(r[k])===v)).sort((a,b)=>num(a.bin)!==null&&num(b.bin)!==null?num(a.bin)-num(b.bin):String(a.bin).localeCompare(String(b.bin),undefined,{numeric:true}));$('distribution-chart').innerHTML=data.length?barChart(data,{title:'Stage 内连续分布',label:r=>r.bin==null?'未标记':r.bin}):empty(notRun);}
const distState=selectionControls('distribution-controls',dist,['stage'],distribution,{allowAll:[],defaultAll:[]});distribution(distState);
[['phenotype',ph.notes],['clock',cl.notes]].forEach(([name,notes])=>{if(arr(notes).length){$(name+'-notes-panel').hidden=false;noteList($(name+'-notes'),notes);}});
noteList($('input-notes'),S.input_notes);

function generic(target,value,depth=0) {
  if(value==null || typeof value==='object'&&Object.keys(value).length===0){target.innerHTML=empty('未提供验证或运行记录。');return;}
  if(Array.isArray(value)){if(rows(value).length===value.length)table(target,value,{filename:'ckm-validation',filters:false});else noteList(target,value);return;}
  if(typeof value!=='object'){target.textContent=valueText(value);return;}
  const scalars=[],nested=[];Object.entries(value).forEach(([k,v])=>{if(v!==null&&typeof v==='object')nested.push([k,v]);else scalars.push({item:k,value:v});});
  if(scalars.length){const box=document.createElement('div');target.append(box);table(box,scalars,{filename:'ckm-run-information',filters:false});}
  nested.forEach(([k,v])=>{const title=document.createElement('h3');title.className='section-title';title.textContent=k;target.append(title);const box=document.createElement('div');target.append(box);if(depth<5)generic(box,v,depth+1);else box.textContent=JSON.stringify(v);});
}
generic($('validation'),S.validation||{});
const tabs=[...document.querySelectorAll('.tab')];
function activate(id,focus=false){tabs.forEach(t=>{const active=t.dataset.panel===id;t.setAttribute('aria-selected',String(active));t.tabIndex=active?0:-1;$('panel-'+t.dataset.panel).hidden=!active;if(active&&focus)t.focus();});}
tabs.forEach((tab,i)=>{tab.addEventListener('click',()=>activate(tab.dataset.panel));tab.addEventListener('keydown',event=>{let next=null;if(event.key==='ArrowRight')next=(i+1)%tabs.length;if(event.key==='ArrowLeft')next=(i+tabs.length-1)%tabs.length;if(event.key==='Home')next=0;if(event.key==='End')next=tabs.length-1;if(next!==null){event.preventDefault();activate(tabs[next].dataset.panel,true);}});});
} catch(error) {
  $('fatal-error').hidden=false;$('fatal-error').textContent='报告无法读取汇总数据：'+String(error&&error.message||error);
  console.error(error);
}
})();
</script>
</body>
</html>'''


def render(data: dict, output: str | Path) -> Path:
    """Validate aggregate data and atomically write a complete offline HTML file."""
    validate_summary(data)
    packed = json.dumps(data, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    packed = (packed.replace("&", "\\u0026").replace("<", "\\u003c")
              .replace(">", "\\u003e").replace("\u2028", "\\u2028")
              .replace("\u2029", "\\u2029"))
    title = html.escape(str(data.get("title") or "CKM 连续分期与阶段差异"), quote=True)
    substitutions = {"__TITLE__": title, "__CKM_JSON__": packed}
    document = re.sub(r"__TITLE__|__CKM_JSON__",
                      lambda match: substitutions[match.group(0)], TEMPLATE)
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=destination.parent,
                                     prefix=".ckm-report-", suffix=".tmp", delete=False) as handle:
        temporary = Path(handle.name)
        try:
            handle.write(document)
            handle.flush()
            os.fsync(handle.fileno())
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    try:
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Aggregate CKM summary JSON")
    parser.add_argument("--output", required=True, help="Self-contained ckm.html destination")
    args = parser.parse_args(argv)
    with open(args.input, encoding="utf-8") as handle:
        data = json.load(handle)
    destination = render(data, args.output)
    print(f"CKM report: {destination}")


if __name__ == "__main__":
    main()
