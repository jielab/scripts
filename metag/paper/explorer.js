'use strict';
const D=JSON.parse(document.getElementById('data').textContent);
const el=id=>document.getElementById(id);
const all='*', fmt=n=>Number(n).toLocaleString('en-US',{maximumFractionDigits:2});
const percent=(a,b)=>b?(100*a/b).toFixed(2)+'%':'NA';
const esc=x=>String(x).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const list=(select,items,first)=>{select.replaceChildren();if(first)select.add(new Option(first,all));items.forEach(x=>select.add(new Option(x,x)));};
const cohort=s=>s.sample.split('/')[0];
const unique=arr=>Array.from(new Set(arr)).sort();
const sum=(arr,f)=>arr.reduce((a,r)=>a+Number(f(r)),0);
const meta={
 status:['执行状态','每个样本一行。not_run 或 split_changed 表示该分析尚无可用结果。'],
 seal_status:['威德尔海豹复核运行记录','每样本一行，记录是否实际运行及参考面板版本。'],
 seal_bins:['威德尔海豹双端比对','9713 为原种级分箱，9712 为原属级分箱。每个 bin 的结局互斥，分母为该 bin 全部读对。'],
 seal_loci:['威德尔海豹参考位置','每个点位窗口 10 kb。比较参考 accession、染色体区段、重复片段端点和跨样本相同窗口。'],
 sources:['食源相容与人源线索','逐样本列出人、植物、哺乳动物及已列入食源词表的其他类群，保留种属层级和原分类证据。'],
 panel:['人和植物的竞争面板','比较对象是否齐全、不同核基因组窗口数及最终种级读对。窗口数为截顶值。'],
 support:['目标物种的分类证据','同一物种可有多种证据。此表只包含 split 预先记录的目标类群种级分配。'],
 blast:['双端 BLAST 抽查','sampled_pairs 是该分箱实际抽查的读对数。参考范围与支持、冲突、分辨率不足等结局共同解释。'],
 amr:['群落 AMR 标记','每百万输入读对归一化。基因家族信号不能确定宿主菌株或表型耐药。'],
 viral:['病毒标签','此表复述原分类器的种级标签。零表示未见相符标签，解释时需结合材料与建库类型。'],
 truth:['独立标签准确率','只针对提供的独立逐对真值。positives、called、tp 分别为真实阳性、判为阳性和真阳性。final 与 full_kraken 可按同一批标签比较。'],
 fastq:['双端序列质量','按样本和 mate 汇总。q30 和 gc 以碱基数加权，范围 0–1。']
};
list(el('dataset'),unique(D.samples.map(cohort)),'全部队列');
list(el('group'),unique(D.catalog.map(t=>t[2])),'全部类群');
list(el('rank'),unique(D.catalog.map(t=>t[3])),'全部层级');
if(D.catalog.some(t=>t[3]==='species'))el('rank').value='species';
Object.entries(meta).forEach(([k,v])=>el('measure').add(new Option(v[0],k)));
el('origin').textContent=D.preview?'Wxiamen 公开汇总表预览 · 所有筛选均在本机完成':'已发布样本 · 离线交互报告';
if(D.preview)el('static').hidden=true;
let cohortRows=[], selected=[], names=new Set(), nNonempty=0, total=0, page=0, taxRows=[],detailRows=[],detailHeaders=[];
const ns='http://www.w3.org/2000/svg';
function svg(id,height,body,title){
 const div=el(id);
 div.innerHTML='<svg xmlns="'+ns+'" viewBox="0 0 640 '+height+'" role="img" aria-label="'+esc(title)+'" style="font:13px sans-serif;background:white"><title>'+esc(title)+'</title>'+body+'</svg>';
}
const txt=(x,y,t,extra='')=>'<text x="'+x+'" y="'+y+'" fill="#23424a" '+extra+'>'+esc(t)+'</text>';
function rect(x,y,w,h,color,tip){
 return '<rect x="'+x+'" y="'+y+'" width="'+w+'" height="'+h+'" fill="'+color+'"><title>'+esc(tip)+'</title></rect>';
}
function bars(id,rows,title){
 const max=Math.max(0,...rows.map(r=>r[1]));
 let body='';
 rows.forEach(([name,value,label],i)=>{
  const y=i*29+9;
  body+=txt(3,y+12,name.length>31?name.slice(0,30)+'…':name)+rect(235,y,max?305*value/max:0,15,'#248678',name+' · '+label)+txt(548,y+12,label);
 });
 svg(id,Math.max(48,rows.length*29+14),body||txt(4,27,'此选区没有记录'),title);
}
function resolution(){
 const size=20, pages=Math.max(1,Math.ceil(selected.length/size));
 page=Math.min(page,pages-1);
 let body='';
 selected.slice(page*size,page*size+size).forEach(({s},i)=>{
  const y=i*23+10,n=Number(s.pairs);let x=225;
  body+=txt(3,y+11,s.sample.length>31?'…'+s.sample.slice(-30):s.sample);
  [['species_pairs','#087a70'],['higher_rank_pairs','#7398a2'],['unclassified_pairs','#dce4e1']].forEach(([k,color])=>{
   const w=n?330*Number(s[k])/n:0;
   body+=rect(x,y,w,13,color,s.sample+' · '+k+' '+fmt(s[k])+' / '+fmt(n));x+=w;
  });
  body+=txt(563,y+11,percent(+s.species_pairs,n));
 });
 svg('resolution',Math.max(48,Math.min(size,selected.length-page*size)*23+17),body,'样本种级、较高层级和未分类比例');
 el('page').textContent=(page+1)+' / '+pages;el('prev').disabled=page===0;el('next').disabled=page+1===pages;
}
function table(id,headers,rows){
 el(id).innerHTML='<table><thead><tr>'+headers.map(x=>'<th>'+esc(x)+'</th>').join('')+'</tr></thead><tbody>'+rows.slice(0,150).map(row=>'<tr>'+row.map(x=>'<td>'+esc(x??'NA')+'</td>').join('')+'</tr>').join('')+'</tbody></table>';
}
const taxHeaders=['taxid','name','group','rank','pairs','input_pairs','ppm','positive_samples','nonempty_samples','detected_pct','top_sample','top_sample_pairs','top_sample_pct'];
function taxa(){
 const g=el('group').value,r=el('rank').value,q=el('query').value.trim().toLowerCase(),sort=el('sort').value;
 const vals=new Map();
 for(const {s,i} of selected){
  const v=D.bins[i];
  for(let j=0;j<v.length;j+=2){
   const id=v[j],n=v[j+1],t=D.catalog[id];
   if(!n||(g!==all&&t[2]!==g)||(r!==all&&t[3]!==r)||(q&&!t[1].toLowerCase().includes(q)&&!t[0].includes(q)))continue;
   let a=vals.get(id);
   if(!a){a={t,pairs:0,seen:0,peakN:0,top:''};vals.set(id,a);}
   a.pairs+=n;a.seen++;
   if(n>a.peakN||(n===a.peakN&&s.sample<a.top)){a.peakN=n;a.top=s.sample;}
  }
 }
 const ordered=Array.from(vals.values());
 for(const a of ordered)a.peak=a.peakN/a.pairs;
 ordered.sort((a,b)=>b[sort]-a[sort]||b.pairs-a.pairs||a.t[0].localeCompare(b.t[0]));
 taxRows=ordered.map(a=>[...a.t,a.pairs,total,Number((1e6*a.pairs/total).toFixed(4)),a.seen,nNonempty,Number((100*a.seen/nNonempty).toFixed(2)),a.top,a.peakN,Number((100*a.peak).toFixed(2))]);
 table('taxTable',taxHeaders,taxRows);el('taxCount').textContent=fmt(taxRows.length)+' 项，屏幕最多显示 150 项。下载包括当前筛选的全部行。';
 bars('leaders',ordered.slice(0,12).map(a=>[a.t[1],a.pairs,fmt(a.pairs)]),'所选分类单元的读对数');
 let b=rect(52,10,555,185,'#eef5f3','');
 const cloud=new Map();
 for(const a of ordered){
  const x=Math.round(52+555*a.seen/nNonempty),y=Math.round(195-185*a.peak);
  const key=x+','+y, p=cloud.get(key)||{x,y,n:0,first:a.t[1]};
  p.n++;cloud.set(key,p);
 }
 for(const p of cloud.values())b+='<circle cx="'+p.x+'" cy="'+p.y+'" r="'+Math.min(8,2+Math.log2(p.n))+'" fill="#06756e" fill-opacity=".56"><title>'+esc(p.first+' · 该位置 '+p.n+' 个分类单元')+'</title></circle>';
 b+=txt(0,22,'100%')+txt(12,196,'0%')+txt(49,214,'0%')+txt(570,214,'100%')+txt(135,237,'横轴：样本检出率　纵轴：最高样本占比');
 svg('scatter',248,b,'各分类单元的样本检出率与单一样本贡献');
}
function evidence(){
 const k=el('measure').value,q=el('find').value.trim().toLowerCase();
 const data=D.detail[k]||[];
 const matches=data.filter(x=>names.has(x.sample)&&(!q||Object.values(x).some(v=>String(v).toLowerCase().includes(q))));
 detailHeaders=data.length?Object.keys(data[0]):[];
 detailRows=matches.map(x=>detailHeaders.map(h=>x[h]));
 el('meaning').textContent=meta[k][1];el('detailCount').textContent=matches.length?fmt(matches.length)+' 行，屏幕最多显示 150 行。':'此选区没有可展示记录。执行情况见上方及“执行状态”表。';
 table('detailTable',detailHeaders,detailRows);el('exportDetail').disabled=!detailRows.length;
}
const cohortCols=['dataset','samples','nonempty_samples','input_pairs','species_pairs','pooled_species_pct','median_sample_species_pct','higher_rank_pct','unclassified_pct'];
function cohortTable(ss){
 const groups=new Map();
 for(const s of ss){const d=cohort(s);if(!groups.has(d))groups.set(d,[]);groups.get(d).push(s);}
 cohortRows=Array.from(groups,([d,rows])=>{
  const n=sum(rows,s=>s.pairs), species=sum(rows,s=>s.species_pairs);
  const rates=rows.filter(s=>+s.pairs>0).map(s=>100*s.species_pairs/s.pairs).sort((a,b)=>a-b);
  const i=Math.floor(rates.length/2);
  const med=rates.length?(rates.length%2?rates[i]:(rates[i]+rates[i-1])/2).toFixed(3):'NA';
  const pct=x=>n?(100*x/n).toFixed(3):'NA';
  return [d,rows.length,rates.length,n,species,pct(species),med,pct(sum(rows,s=>s.higher_rank_pairs)),pct(sum(rows,s=>s.unclassified_pairs))];
 });
 table('cohortTable',['队列','样本','非空样本','输入读对','种级读对','合并种级 %','样本种级中位 %','较高层级 %','未分类 %'],cohortRows);
 const methods=unique(ss.map(s=>s.code)).length, dbs=unique(ss.map(s=>s.database)).length;
 el('versions').textContent='当前选区有 '+methods+' 个记录的方法签名和 '+dbs+' 个分类库签名。逐样本版本见 samples.tsv。';
}
function refresh(){
 selected=D.samples.map((s,i)=>({s,i})).filter(({s})=>(el('dataset').value===all||cohort(s)===el('dataset').value)&&(el('sample').value===all||s.sample===el('sample').value));
 names=new Set(selected.map(x=>x.s.sample));const ss=selected.map(x=>x.s);
 nNonempty=ss.filter(s=>+s.pairs>0).length;total=sum(ss,s=>s.pairs);
 const exact=sum(ss,s=>s.species_pairs),higher=sum(ss,s=>s.higher_rank_pairs),unc=sum(ss,s=>s.unclassified_pairs);
 el('kpis').innerHTML=[['已完成样本',fmt(ss.length)],['输入读对',fmt(total)],['种级读对',percent(exact,total)],['较高层级',percent(higher,total)],['未分类',percent(unc,total)]].map(([l,v])=>'<div>'+l+'<b>'+v+'</b></div>').join('');
 cohortTable(ss);resolution();const reasons=new Map();
 for(const x of D.detail.evidence)if(names.has(x.sample))reasons.set(x.evidence,(reasons.get(x.evidence)||0)+Number(x.pairs));
 bars('reasons',Array.from(reasons).sort((a,b)=>b[1]-a[1]).map(([k,v])=>[k,total?100*v/total:0,percent(v,total)]),'分类证据占全部输入读对比例');
 const status=D.detail.status.filter(x=>names.has(x.sample)&&x.state==='complete');
 const count=k=>status.filter(x=>x[k]==='complete').length;
 el('coverage').textContent='当前选区 '+ss.length+' 个样本。证据汇总完成 '+status.length+'，BLAST 完成 '+count('blast')+'，AMR 完成 '+count('amr')+'，有独立真值标签 '+ss.filter(x=>Number(x.truth_pairs)>0).length+'，逐对核查 FASTQ '+ss.filter(x=>['bins','raw_and_bins'].includes(x.integrity)).length+'。';
 if(!D.detail.status.length)el('coverage').textContent='本报告的输入未提供 refine 执行记录。BLAST、AMR 的覆盖情况未知。当前选区有独立真值标签的样本 '+ss.filter(x=>Number(x.truth_pairs)>0).length+' 个。';
 const seal=(D.detail.seal_bins||[]).filter(x=>names.has(x.sample)&&String(x.bin_taxid)==='9713');
 el('sealSection').hidden=!seal.length;
 if(seal.length){
  const panelCount=new Set(seal.map(x=>x.panel_id)).size;
  el('sealNote').textContent='按所选样本合计，每项单位为读对。种级分箱使用 taxid 9713。属级分箱 9712 单列在证据表中。参考面板版本数 '+panelCount+'。';
  const categories=[['全部种级分箱','tested_pairs'],['成对支持核染色体','target_nuclear'],['其他或未定位 contig','target_other_contig'],['两端支持其他分类单元','alternative_taxon'],['海豹科内共享','shared_phocid']];
  bars('sealChart',categories.map(([label,k])=>[label,sum(seal,r=>r[k]),fmt(sum(seal,r=>r[k]))]),'威德尔海豹种级标签序列复核');
 }
 taxa();evidence();
}
function chooseDataset(){
 const ds=el('dataset').value;
 list(el('sample'),D.samples.filter(s=>ds===all||cohort(s)===ds).map(s=>s.sample),'全部样本');
 page=0;refresh();
}
function download(content,type,file){
 const url=URL.createObjectURL(new Blob([content],{type})),link=document.createElement('a');
 link.href=url;link.download=file;link.click();setTimeout(()=>URL.revokeObjectURL(url),3000);
}
function exportTable(headers,rows,name){
 const field=v=>String(v??'NA').replace(/[\t\r\n]/g,' ');
 download([headers,...rows].map(r=>r.map(field).join('\t')).join('\n')+'\n','text/tab-separated-values;charset=utf-8',name+'.tsv');
}
el('dataset').onchange=chooseDataset;el('sample').onchange=()=>{page=0;refresh();};
for(const id of ['group','rank','sort'])el(id).onchange=taxa;
el('query').oninput=taxa;el('find').oninput=evidence;el('measure').onchange=evidence;
el('prev').onclick=()=>{page--;resolution();};el('next').onclick=()=>{page++;resolution();};
el('exportCohorts').onclick=()=>exportTable(cohortCols,cohortRows,'cohort_comparison');
el('exportTax').onclick=()=>exportTable(taxHeaders,taxRows,'taxa_selection');
el('exportDetail').onclick=()=>exportTable(detailHeaders,detailRows,el('measure').value+'_selection');
document.querySelectorAll('[data-save]').forEach(b=>b.onclick=()=>{const id=b.dataset.save;download(el(id).querySelector('svg').outerHTML,'image/svg+xml',id+'.svg');});
chooseDataset();
