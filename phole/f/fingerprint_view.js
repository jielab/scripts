/* Offline view of reported-locus hypotheses. No captured identity is inferred. */
(function(global){
'use strict';
const D={ukb:'UKB',ckb:'CKB',decode:'deCODE',scot:'SCOT'};
const TIERS={T1_distinctive_heldout:'去除提名区后仍有双平台非热点支持',T2_specific_full_profile:'完整指纹优先项；严格复核未全部通过',T3_ambiguous_multilocus:'多区域线索；含靶标或热点歧义',T4_specific_anchor_only:'单个特定锚点',T5_regional_or_hotspot_context:'区域／热点背景'};
const STATES={cis_gwas_only:'常规显著，未达严格 cis',cis_search_censored:'旧摘要可能漏收 cis',no_cis_reported_uncertified:'未报告 cis，完整阴性未认证',no_gwas_cis_certified:'完整窗口未检出显著 cis',cis_unknown:'cis 未知'};
const esc=v=>String(v==null?'':v).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const yes=v=>v===true||v===1||v==='1';
const yn=v=>v==null||v==='.'?'未知':yes(v)?'是':'未满足';
const val=v=>v==null||v==='.'||v===''?'—':String(v);
const num=v=>v==null||v==='.'?'—':Number(v).toLocaleString('en-US');
const frac=v=>v==null||v==='.'?'—':Number(v).toFixed(3);
const ds=v=>D[v]||val(v);
function download(name,rows){
 if(!rows.length)return;
 const keys=[...new Set(rows.flatMap(r=>Object.keys(r)))];
 const q=v=>/[\t\r\n"]/.test(String(v))?'"'+String(v).replace(/"/g,'""')+'"':String(v);
 const text=[keys.map(q).join('\t'),...rows.map(r=>keys.map(k=>q(r[k]==null?'.':r[k])).join('\t'))].join('\n')+'\n';
 const blob=new Blob(['\ufeff'+text],{type:'text/tab-separated-values;charset=utf-8'}),url=URL.createObjectURL(blob),a=document.createElement('a');
 a.href=url;a.download=name.replace(/[^\w.\-]/g,'_');a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
}
function table(head,rows,empty='没有符合当前条件的记录。'){
 return '<div class="fp-table-wrap"><table><thead><tr>'+head.map(x=>'<th>'+esc(x)+'</th>').join('')+'</tr></thead><tbody>'+ (rows.length?rows.join(''):'<tr><td colspan="'+head.length+'" class="fp-empty">'+esc(empty)+'</td></tr>')+'</tbody></table></div>';
}
function mount(root,data,options={}){
 if(!root)return null;
 if(!data){root.innerHTML='<p class="fp-note">本页未附遗传指纹详细复核。运行 phole.sh cross-check 可生成。</p>';return null;}
 const pairs=data.pairs||[],hyps=data.hypotheses||[],matches=data.matched_loci||[],contexts=data.same_platform_context||[],signals=data.profile_signals||[],sens=data.sensitivity||[];
 const index=new Map(hyps.map(r=>[r.HYPOTHESIS_ID,r]));
 const prior=pairs.filter(r=>yes(r.PRIOR_PRIORITY_CASE)),strict=pairs.filter(r=>r.TIER==='T1_distinctive_heldout');
 let selected=prior.length?prior[0]:pairs[0],filtered=[],page=0;
 root.classList.add('fp-review');
 root.innerHTML='<div class="fp-summary"><div><span>定向蛋白对</span><strong>'+num(pairs.length)+'</strong></div><div><span>原完整指纹重点</span><strong>'+num(prior.length)+'</strong></div><div><span>双平台非热点严格层</span><strong>'+num(strict.length)+'</strong></div><div><span>已确证抓错蛋白</span><strong>0</strong></div></div>'+
 '<p class="fp-note">所有数量按明确筛查规则给出。完整指纹、去除提名区域后的证据、同队列相似性分开呈现。共有 lead 与 Jaccard 描述筛选后的位点，尚未验证方向、LD 或捕获身份。</p>'+
 '<div class="fp-controls"><label>数据 <select data-fp="dataset"><option value="all">四组全部</option>'+Object.entries(D).map(([k,v])=>'<option value="'+k+'">'+v+'</option>').join('')+'</select></label><label>证据 <select data-fp="tier"><option value="prior">完整指纹重点</option><option value="exact">精确 anchor + 多区域</option><option value="T1_distinctive_heldout">双平台非热点严格层</option><option value="multi">至少两个匹配区域</option><option value="weak">弱 cis 扩展</option><option value="all">所有定向蛋白对</option></select></label><label class="fp-search">检索 <input data-fp="search" placeholder="SEMA3G、NID1 或 assay"></label><button type="button" data-fp="export">导出当前蛋白对</button></div>'+
 '<div data-fp="table"></div><div class="fp-pager"><span data-fp="count"></span><div><button type="button" data-fp="prev">上一页</button><span data-fp="page"></span><button type="button" data-fp="next">下一页</button></div></div>'+
 '<article class="fp-detail"><div class="fp-detail-heading"><div><h3 data-fp="title"></h3><p data-fp="sub"></p></div><div><button type="button" data-fp="locate">查看主报告</button> <button type="button" data-fp="exportCase">导出病例证据</button></div></div><label class="fp-probe-select">检测／参考组合 <select data-fp="hypothesis"></select></label><div data-fp="case"></div></article>'+
 '<section class="fp-within"><h3>同队列高相似性补充筛查</h3><p class="fp-note">预设描述条件：至少 4 个匹配区域且 Jaccard ≥0.60；仍要求 primary 的自身 cis 不强且有异平台自身强 cis。相同队列内的比较用于检查测量对象与命名，不计作独立复制。点击任一行查看实际 SNP、P 和外部参照。</p><div class="fp-controls"><label>显示 <select data-fp="withinScope"><option value="all">全部补充病例</option><option value="other">替代目标非 TNFSF14 / HSPA1L / VARS</option><option value="hub">替代目标为 TNFSF14 / HSPA1L / VARS</option></select></label><span data-fp="withinCount"></span></div><div data-fp="within"></div><button type="button" data-fp="exportWithin">导出全部补充筛查</button><article class="fp-detail" data-fp="withinDetail"></article></section>';
 const $=name=>root.querySelector('[data-fp="'+name+'"]');
 if(!options.onTargetSelect)$('locate').hidden=true;
 function getHyp(){return index.get($('hypothesis').value)||selected;}
 function row(r){return '<tr data-pair="'+esc(r.PAIR_ID)+'" tabindex="0"'+(selected&&r.PAIR_ID===selected.PAIR_ID?' class="fp-selected"':'')+'><td><b>'+esc(r.PROT)+'</b><small>'+esc(ds(r.PRIMARY))+' · '+esc(r.PRIMARY_ASSAY)+'</small></td><td><b>'+esc(r.ALTERNATIVE_PROT)+'</b><small>'+esc(ds(r.REFERENCE))+' · '+esc(r.REFERENCE_ASSAY)+'</small></td><td>'+esc(TIERS[r.TIER]||r.TIER)+'<small>'+esc(r.SCREEN_STRATUM)+'</small></td><td>'+num(r.FULL_SHARED_REGIONS)+' / '+num(r.FULL_EXACT_RSID)+'<small>J='+frac(r.FULL_JACCARD)+'；自身 J='+frac(r.FULL_INTENDED_MAX_JACCARD)+'</small></td><td>'+num(r.REGION_HELDOUT_SHARED_REGIONS)+' 区域<small>双平台非热点 '+num(r.REGION_HELDOUT_NONHOTSPOT_BOTH)+'；替代优于自身 '+(yes(r.REGION_HELDOUT_PREFERRED_TO_ALL_INTENDED)?'是':'未满足')+'</small></td></tr>';}
 function refresh(reset=true){
  if(reset)page=0;
  const d=$('dataset').value,t=$('tier').value,q=$('search').value.trim().toLowerCase();
  const exactPairs=new Set(hyps.filter(r=>yes(r.EXACT_ANCHOR_MULTILOCUS)).map(r=>r.PAIR_ID));
  filtered=pairs.filter(r=>(d==='all'||r.PRIMARY===d)&&(t==='all'||t==='prior'&&yes(r.PRIOR_PRIORITY_CASE)||t==='exact'&&exactPairs.has(r.PAIR_ID)||t==='multi'&&Number(r.FULL_SHARED_REGIONS)>=2||t==='weak'&&String(r.SCREEN_STRATA||r.SCREEN_STRATUM).includes('weak')||r.TIER===t)&&(!q||[r.PROT,r.ALTERNATIVE_PROT,r.PRIMARY_ASSAY,r.REFERENCE_ASSAY].join(' ').toLowerCase().includes(q)));
  if(!selected||!filtered.some(r=>r.PAIR_ID===selected.PAIR_ID)){selected=filtered[0]||null;renderDetail();}
  const pages=Math.max(1,Math.ceil(filtered.length/25));page=Math.max(0,Math.min(page,pages-1));
  $('table').innerHTML=table(['被核查的读出','替代蛋白参考','证据层级','完整指纹：区域 / 相同 rsID','去除提名区域'],filtered.slice(page*25,page*25+25).map(row));
  $('count').textContent=filtered.length.toLocaleString()+' 个定向蛋白对';$('page').textContent=(page+1)+' / '+pages;$('prev').disabled=page===0;$('next').disabled=page+1>=pages;$('export').disabled=!filtered.length;
 }
 function renderDetail(){
  if(!selected){$('title').textContent='当前没有匹配病例';$('sub').textContent='';$('hypothesis').innerHTML='';$('case').innerHTML='';$('exportCase').disabled=true;$('locate').disabled=true;return;}
  $('title').textContent=ds(selected.PRIMARY)+' '+selected.PROT+' → '+selected.ALTERNATIVE_PROT;
  $('sub').textContent='定向蛋白对包含 '+val(selected.N_ASSAY_HYPOTHESES)+' 个检测／参考组合；主表显示其中的代表组合。';
  const hs=hyps.filter(r=>r.PAIR_ID===selected.PAIR_ID);
  $('hypothesis').innerHTML=(hs.length?hs:[selected]).map(r=>'<option value="'+esc(r.HYPOTHESIS_ID)+'">'+esc(ds(r.PRIMARY)+' '+r.PRIMARY_ASSAY+' / '+ds(r.REFERENCE)+' '+r.REFERENCE_ASSAY)+'</option>').join('');
  $('hypothesis').value=hs.some(r=>r.HYPOTHESIS_ID===selected.HYPOTHESIS_ID)||!hs.length?selected.HYPOTHESIS_ID:hs[0].HYPOTHESIS_ID;$('exportCase').disabled=false;$('locate').disabled=false;renderCase();
 }
 function profilePlot(rows){
  const axes=data.axes||[];if(!axes.length||!rows.length)return '';
  const order=['primary_intended','primary_alternative','reference_intended','reference_alternative','external_intended','external_alternative'];
  const groups=new Map();for(const r of rows){const key=[r.ROLE,r.DATASET,r.ASSAY].join('::');if(!groups.has(key))groups.set(key,[]);groups.get(key).push(r);}
  const tracks=[...groups.values()].sort((a,b)=>order.indexOf(a[0].ROLE)-order.indexOf(b[0].ROLE)||a[0].ASSAY.localeCompare(b[0].ASSAY));
  const w=1200,left=225,right=28,top=28,rowH=43,h=top+tracks.length*rowH+45,total=Math.max(...axes.map(a=>Number(a.offset)+Number(a.span))),byChr=new Map(axes.map(a=>[String(a.chr),a]));
  const x=p=>left+Number(p)/total*(w-left-right);
  let s='<svg class="fp-profile-plot" viewBox="0 0 '+w+' '+h+'" role="img" aria-label="被核查与替代蛋白的基因组位点指纹；红 cis，绿 trans，竖线等高">';
  axes.forEach((a,i)=>{const x1=x(a.offset),x2=x(Number(a.offset)+Number(a.span));s+='<rect x="'+x1+'" y="14" width="'+(x2-x1)+'" height="'+(h-46)+'" fill="'+(i%2?'#f4f7f9':'#ffffff')+'"/><text x="'+((x1+x2)/2)+'" y="'+(h-14)+'" text-anchor="middle" fill="#718391" font-size="10">'+esc(a.chr)+'</text>';});
  tracks.forEach((rs,i)=>{const r=rs[0],y=top+i*rowH+17,kept=rs.filter(v=>yes(v.IN_STRICT_PROFILE));s+='<text x="8" y="'+(y-4)+'" fill="#294a5e" font-size="12">'+esc(ds(r.DATASET)+' · '+r.ASSAY)+'</text><text x="8" y="'+(y+10)+'" fill="#7b8b96" font-size="9">'+esc(r.ROLE)+' · '+kept.length+' lead</text><line x1="'+left+'" x2="'+(w-right)+'" y1="'+y+'" y2="'+y+'" stroke="#d9e3e8"/>';
   for(const q of kept){const a=byChr.get(String(q.CHR));if(!a)continue;const pos=x(Number(a.offset)+Number(q.POS));s+='<line x1="'+pos+'" x2="'+pos+'" y1="'+(y-11)+'" y2="'+(y+11)+'" stroke="'+(q.CLASS==='cis'?'#c84658':'#188369')+'" stroke-width="2.4"><title>'+esc(q.SNP+' · chr'+q.CHR+':'+q.POS+' · P='+q.P+' · '+q.CLASS)+'</title></line>';}
  });
  return '<div class="fp-plot-wrap">'+s+'</svg></div><p class="fp-note">按基因组坐标绘制所有严格阈值且通过当前 QC 的已报告 lead；红 cis、绿 trans，竖线等高。悬停可看原始 P。图中 lead 数可能与合并后的区域数不同；缺失区域不表示效应为零。</p>';
 }
 function renderCase(){
  const h=getHyp();if(!h)return;
  const id=h.HYPOTHESIS_ID,detail=index.has(id),m=matches.filter(r=>r.HYPOTHESIS_ID===id),c=contexts.filter(r=>r.HYPOTHESIS_ID===id);
  if(!detail){$('case').innerHTML='<p class="fp-note">本 HTML 为全部 '+num(pairs.length)+' 个定向蛋白对提供摘要，逐位点详情内嵌 '+num(hyps.length)+' 个重点／精确 anchor 多区域组合。当前背景蛋白对的完整记录见 fingerprints/assay-evidence.tsv、pair-summary.tsv；未内嵌信息不表示阴性结果。</p>'+table(['完整重合区域','相同 rsID','替代 Jaccard','自身最佳 Jaccard'],['<tr><td>'+num(h.FULL_SHARED_REGIONS)+'</td><td>'+num(h.FULL_EXACT_RSID)+'</td><td>'+frac(h.FULL_JACCARD)+'</td><td>'+frac(h.FULL_INTENDED_MAX_JACCARD)+'</td></tr>']);return;}
  const controls=[['完整指纹','FULL'],['去除提名区域','REGION_HELDOUT'],['去除 5 Mb 分组','BLOCK5MB_HELDOUT'],['去除提名染色体','CHROMOSOME_HELDOUT']];
  let html='<p class="fp-case-intro">'+esc(STATES[h.PRIMARY_STATE]||h.PRIMARY_STATE)+'；primary 自身 cis P='+esc(val(h.PRIMARY_CIS_P))+'。同一参考队列的自身强 cis：'+(yes(h.SAME_REFERENCE_INTENDED_STRICT_CIS)?'有':'未确认')+'，P='+esc(val(h.REFERENCE_INTENDED_CIS_P))+'。</p>'+
   '<div class="fp-anchor"><b>提名 anchor</b><span>'+esc(h.ANCHOR_SNP)+' · primary P='+esc(val(h.ANCHOR_PRIMARY_P))+' · 替代蛋白 cis P='+esc(val(h.ANCHOR_REFERENCE_P))+'</span></div>';
  html+=table(['比较范围','与替代蛋白重合','相同 rsID','Jaccard','自身最佳 Jaccard','优于全部自身 probes','两边均非热点','共同覆盖蛋白敏感性'],controls.map(([label,k])=>'<tr><td>'+label+'</td><td>'+num(h[k+'_SHARED_REGIONS'])+'</td><td>'+num(h[k+'_EXACT_RSID'])+'</td><td>'+frac(h[k+'_JACCARD'])+'</td><td>'+frac(h[k+'_INTENDED_MAX_JACCARD'])+'</td><td>'+yn(h[k+'_PREFERRED_TO_ALL_INTENDED'])+'</td><td>'+num(h[k+'_NONHOTSPOT_BOTH'])+'</td><td>'+num(h[k+'_NONHOTSPOT_BOTH_COMMON_UNIVERSE'])+'</td></tr>'));
  html+='<p class="fp-note">“去除”在比较双方使用相同坐标规则，并重新比较所有可用的同名 reference probes。5 Mb 和染色体是敏感性分析，不代表已经确定 LD 独立性。非热点阈值为各平台当前摘要涉及 ≤10 个目标，覆盖范围的影响见完整 TSV。</p>';
  html+='<p class="fp-cautions"><b>未通过严格层的原因：</b>'+esc(val(h.STRICT_TIER_NOT_MET_REASONS))+'<br><b>其他提示：</b>'+esc(val(h.CAUTIONS))+'</p>';
  html+='<h4>检测与替代蛋白的位点指纹</h4>'+profilePlot(signals.filter(r=>r.HYPOTHESIS_ID===id));
  const full=m.filter(r=>r.MODE==='FULL'&&r.COMPARISON==='alternative');
  html+='<h4>完整指纹的逐区域匹配</h4>'+table(['位置与 primary SNP','Primary P / CLASS','参考 SNP','参考 P / CLASS','精确对应','热点程度'],full.map(r=>'<tr><td>chr'+esc(r.CHR)+':'+esc(r.PRIMARY_POS)+'<small>'+esc(r.PRIMARY_SNP)+'</small></td><td class="fp-p">'+esc(r.PRIMARY_P)+'<small>'+esc(r.PRIMARY_CLASS)+'</small></td><td>'+esc(r.REFERENCE_SNP)+'<small>chr'+esc(r.CHR)+':'+esc(r.REFERENCE_POS)+'</small></td><td class="fp-p">'+esc(r.REFERENCE_P)+'<small>'+esc(r.REFERENCE_CLASS)+'</small></td><td>'+(yes(r.EXACT_RSID)?'同 rsID / 位置':num(r.DISTANCE_BP)+' bp')+'</td><td>两边非热点：'+(yes(r.NONHOTSPOT_BOTH)?'是':'否')+'<small>Primary '+esc(val(r.PRIMARY_REGION_TARGETS||r.PRIMARY_REGION_TRANS_TARGETS))+' / Reference '+esc(val(r.REFERENCE_REGION_TARGETS||r.REFERENCE_REGION_TRANS_TARGETS))+'</small></td></tr>'));
  html+='<h4>同一批对象的四种比较</h4>'+table(['比较类型','左侧测定','右侧测定','匹配区域 / 相同 rsID','Jaccard'],c.map(r=>'<tr><td>'+esc(r.COMPARISON)+'</td><td>'+esc(ds(r.LEFT_DATASET)+' '+r.LEFT_PROT)+'<small>'+esc(r.LEFT_ASSAY)+'</small></td><td>'+esc(ds(r.RIGHT_DATASET)+' '+r.RIGHT_PROT)+'<small>'+esc(r.RIGHT_ASSAY)+'</small></td><td>'+num(r.SHARED_REGIONS)+' / '+num(r.EXACT_RSID)+'</td><td>'+frac(r.JACCARD)+'</td></tr>'))+'<p class="fp-note">同队列的高相似性也可能来自共调控、复合物、样本处理或标签问题，不计作新的独立队列支持。</p>';
  const sr=sens.filter(r=>r.HYPOTHESIS_ID===id&&r.MODE==='REGION_HELDOUT');
  html+='<details><summary>匹配距离敏感性</summary>'+table(['匹配跨度','去提名区后重合','相同 rsID','优于自身 probes'],sr.map(r=>'<tr><td>'+num(r.MATCH_SPAN_BP)+' bp</td><td>'+num(r.PROFILE_SHARED_REGIONS)+'</td><td>'+num(r.PROFILE_EXACT_RSID)+'</td><td>'+(yes(r.PROFILE_PREFERRED_TO_ALL_INTENDED)?'是':'未满足')+'</td></tr>'))+'</details>';
  const ps=signals.filter(r=>r.HYPOTHESIS_ID===id);
  html+='<details><summary>四组测定的原始已报告位点（'+num(ps.length)+' 行）</summary>'+table(['角色','测定','位置','SNP','原始 P','CLASS / QC'],ps.map(r=>'<tr><td>'+esc(r.ROLE)+'</td><td>'+esc(ds(r.DATASET)+' '+r.ASSAY)+'</td><td>chr'+esc(r.CHR)+':'+esc(r.POS)+'</td><td>'+esc(r.SNP)+'</td><td class="fp-p">'+esc(r.P)+'</td><td>'+esc(r.CLASS)+' / '+(yes(r.QC_USABLE)?'可用':'需核查')+'</td></tr>'))+'</details>';
  $('case').innerHTML=html;
 }
 function renderWithin(){
  const all=data.same_cohort_screen||[],scope=$('withinScope').value,hub=r=>['TNFSF14','HSPA1L','VARS'].includes(r.ALTERNATIVE_PROT),rows=all.filter(r=>scope==='all'||(scope==='hub'?hub(r):!hub(r)));
  $('withinCount').textContent=rows.length+' / '+all.length+' 个 assay 配对';
  $('within').innerHTML=table(['队列','被核查读出','同队列替代蛋白','匹配区域 / 相同 rsID','Jaccard','非热点区域'],rows.map(r=>'<tr data-within="'+esc(r.HYPOTHESIS_ID)+'" tabindex="0"><td>'+esc(ds(r.PRIMARY))+'</td><td>'+esc(r.PROT)+'<small>'+esc(r.PRIMARY_ASSAY)+'</small></td><td>'+esc(r.ALTERNATIVE_PROT)+'<small>'+esc(r.ALTERNATIVE_ASSAY)+(hub(r)?' · 常见热点参照':'')+'</small></td><td>'+num(r.SHARED_REGIONS)+' / '+num(r.EXACT_RSID)+'</td><td>'+frac(r.JACCARD)+'</td><td>'+num(r.NONHOTSPOT_BOTH)+'</td></tr>'));
  $('exportWithin').disabled=!rows.length;
  renderWithinDetail(rows[0]);
 }
 function renderWithinDetail(r){
  if(!r){$('withinDetail').innerHTML='<p class="fp-note">没有匹配记录。</p>';return;}
  const m=(data.same_cohort_screen_matches||[]).filter(x=>x.HYPOTHESIS_ID===r.HYPOTHESIS_ID&&x.MODE==='FULL');
  let html='<h3>'+esc(ds(r.PRIMARY)+' '+r.PRIMARY_ASSAY+' → '+r.ALTERNATIVE_ASSAY)+'</h3><p class="fp-note">'+esc(STATES[r.PRIMARY_STATE]||r.PRIMARY_STATE)+'。异平台自身强 cis：'+esc(r.ORTHOGONAL_INTENDED_CIS_COHORTS)+'；最佳 P='+esc(val(r.ORTHOGONAL_INTENDED_BEST_CIS_P))+'。同队列提名 anchor：'+esc(r.ANCHOR_SNP)+'，匹配方式：'+esc(r.ANCHOR_MATCH)+'。</p>';
  html+=profilePlot((data.same_cohort_profile_signals||[]).filter(x=>x.HYPOTHESIS_ID===r.HYPOTHESIS_ID));
  html+=table(['Primary SNP / 位置','Primary P','替代蛋白 SNP / 位置','替代蛋白 P','精确对应','两边非热点'],m.map(x=>'<tr><td>'+esc(x.PRIMARY_SNP)+'<small>chr'+esc(x.CHR)+':'+esc(x.PRIMARY_POS)+'</small></td><td class="fp-p">'+esc(x.PRIMARY_P)+'</td><td>'+esc(x.REFERENCE_SNP)+'<small>chr'+esc(x.CHR)+':'+esc(x.REFERENCE_POS)+'</small></td><td class="fp-p">'+esc(x.REFERENCE_P)+'</td><td>'+(yes(x.EXACT_RSID)?'同 rsID / 位置':num(x.DISTANCE_BP)+' bp')+'</td><td>'+yn(x.NONHOTSPOT_BOTH)+'</td></tr>'));
  html+='<h4>外部平台参考</h4>'+table(['对照','队列 / assay','共享区域 / rsID','Jaccard'],['INTENDED','ALTERNATIVE'].map(k=>'<tr><td>'+(k==='INTENDED'?'同名蛋白':'替代蛋白')+'</td><td>'+esc(ds(r['BEST_ORTHOGONAL_'+k+'_COHORT']))+' / '+esc(val(r['BEST_ORTHOGONAL_'+k+'_ASSAY']))+'</td><td>'+num(r['BEST_ORTHOGONAL_'+k+'_SHARED_REGIONS'])+' / '+num(r['BEST_ORTHOGONAL_'+k+'_EXACT_RSID'])+'</td><td>'+frac(r['BEST_ORTHOGONAL_'+k+'_JACCARD'])+'</td></tr>'))+'<p class="fp-note">此表分别展示最佳可用外部组合，可能来自不同队列；并不自动满足主筛查的同一 reference 对照要求。缺测是未知。共同调控、复合物、热点和 assay 标签问题仍需区分。</p><button type="button" data-within-export="'+esc(r.HYPOTHESIS_ID)+'">导出该配对逐位点证据</button>';
  $('withinDetail').innerHTML=html;
  $('withinDetail').querySelector('[data-within-export]').onclick=()=>download(r.PRIMARY_ASSAY+'.'+r.PRIMARY+'.'+r.ALTERNATIVE_ASSAY+'.within-matches.tsv',m);
 }
 for(const name of ['dataset','tier'])$(name).addEventListener('change',()=>refresh());
 $('search').addEventListener('input',()=>refresh());$('prev').onclick=()=>{page--;refresh(false);};$('next').onclick=()=>{page++;refresh(false);};
 $('table').addEventListener('click',e=>{const row=e.target.closest('tr[data-pair]');if(row){selected=pairs.find(r=>r.PAIR_ID===row.dataset.pair);renderDetail();refresh(false);}});
 $('table').addEventListener('keydown',e=>{if(e.key==='Enter'){const row=e.target.closest('tr[data-pair]');if(row)row.click();}});
 $('hypothesis').onchange=renderCase;$('export').onclick=()=>download('fingerprint-pairs.filtered.tsv',filtered);
 $('withinScope').onchange=renderWithin;
 $('within').addEventListener('click',e=>{const row=e.target.closest('[data-within]');if(row)renderWithinDetail((data.same_cohort_screen||[]).find(r=>r.HYPOTHESIS_ID===row.dataset.within));});
 $('within').addEventListener('keydown',e=>{if(e.key==='Enter'){const row=e.target.closest('[data-within]');if(row)row.click();}});
 $('exportCase').onclick=()=>{const h=getHyp();if(!h)return;download(h.PROT+'.'+h.PRIMARY+'.'+h.ALTERNATIVE_PROT+'.summary.tsv',[h]);download(h.PROT+'.'+h.PRIMARY+'.'+h.ALTERNATIVE_PROT+'.matched-loci.tsv',matches.filter(r=>r.HYPOTHESIS_ID===h.HYPOTHESIS_ID));download(h.PROT+'.'+h.PRIMARY+'.'+h.ALTERNATIVE_PROT+'.contexts.tsv',contexts.filter(r=>r.HYPOTHESIS_ID===h.HYPOTHESIS_ID));};
 $('locate').onclick=()=>{const h=getHyp();if(h&&options.onTargetSelect)options.onTargetSelect(h.PRIMARY,h.PROT);};
 $('exportWithin').onclick=()=>download('fingerprint-same-cohort-high-similarity.tsv',data.same_cohort_screen||[]);
 renderDetail();refresh();renderWithin();
 return {select:(primary,prot,alt)=>{const r=pairs.find(r=>r.PRIMARY===primary&&r.PROT===prot&&(!alt||r.ALTERNATIVE_PROT===alt));if(r){$('tier').value='all';$('dataset').value=primary;$('search').value=prot;selected=r;renderDetail();refresh();}},data};
}
global.PHOLEFingerprintReview={mount};
})(window);
