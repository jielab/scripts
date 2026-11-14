/* Offline fixed-coordinate explorer. No network, PCA fitting or inference. */
(function () {
'use strict';
const MISSING = new Set(['', '.', 'NA', 'NAN', 'NULL', 'NONE']);
const COLORS = ['#49b8e8','#f6a453','#46c9a4','#bd91e1','#e6cb5f','#ee7f9c','#73a9fa','#96be77','#c9996e','#86cbd0'];
const GREY = '#748b9d';
function number(value) {
  if (value === null || value === undefined || typeof value === 'boolean') return null;
  if (typeof value === 'string' && MISSING.has(value.trim().toUpperCase())) return null;
  const n = Number(value); return Number.isFinite(n) ? n : null;
}
function category(value) {
  return value === null || value === undefined || (typeof value === 'string' && MISSING.has(value.trim().toUpperCase())) ? null : String(value);
}
function extent(values) {
  let lo = Infinity, hi = -Infinity;
  for (const value of values) if (typeof value === 'number' && Number.isFinite(value)) { lo = Math.min(lo, value); hi = Math.max(hi, value); }
  return lo === Infinity ? null : [lo, hi];
}
function makeTracks(data) {
  const points = data.points;
  const tracks = [
    {id:'region',label:'Broad geographic region',kind:'categorical',values:points.map(p=>category(p.region)),units:'region label',provenance:'Supplied sample metadata; not inferred genetic ancestry',palette:data.region_colors || {}},
    {id:'cohort',label:'Ancient / modern',kind:'categorical',values:points.map(p=>p.cohort),units:'cohort label',provenance:'Supplied cohort labels',palette:{ancient:'#49b8e8',modern:'#f6a453'}},
    {id:'source',label:'Data source',kind:'categorical',values:points.map(p=>category(p.source)),units:'source label',provenance:'Supplied sample metadata'},
    {id:'age',label:'Recorded age (BP)',kind:'numeric',values:points.map(p=>p.age_bp),units:'years BP',provenance:'Original recorded dates, unchanged'},
    {id:'layer',label:'Age-bin midpoint (display)',kind:'numeric',values:points.map(p=>p.cohort==='modern'?0:(Math.floor(p.age_bp/data.display.layer_years)+.5)*data.display.layer_years),units:'display midpoint BP',provenance:'Display-only age-bin midpoint; original dates remain unchanged'}
  ];
  for (const [key,label] of [['coverage','Target coverage'],['iqs','IQS']]) {
    const values = points.map(p=>number(p[key]));
    if (values.some(x=>x !== null)) tracks.push({id:key,label,kind:'numeric',values,units:'as supplied',provenance:'Supplied sample-quality metadata; no missing-value imputation'});
  }
  (data.features || []).forEach(f=>tracks.push({id:'feature:'+f.feature_id,label:f.label,kind:f.kind,values:f.values.slice(),units:f.units,provenance:f.provenance}));
  (data.variants || []).forEach((v,i)=>tracks.push({id:'variant:'+i,label:v.variant_id+' · '+v.ref+'→'+v.alt,kind:'numeric',values:v.values === null ? points.map(()=>null) : v.values.slice(),units:v.count_kind === 'hardcall' ? 'ALT copies (0/1/2)' : 'ALT dosage (0–2)',provenance:'Supplied '+v.count_kind+' · '+(v.status || 'status undeclared')+'. 0/1/2 is one REF/ALT comparison, not three nucleotides.',variant:v}));
  return tracks.map(t=>{
    const observed=t.values.filter(v=>v !== null && v !== undefined);
    const categories=t.kind==='categorical' ? [...new Set(observed.map(String))].sort((a,b)=>a<b?-1:a>b?1:0) : [];
    const colorMap = new Map(categories.map((v,i)=>[v,(t.palette && Object.prototype.hasOwnProperty.call(t.palette,v))?t.palette[v]:COLORS[i % COLORS.length]]));
    colorMap.set(null,GREY);
    return Object.assign({},t,{range:t.kind==='numeric'?extent(t.values):null,categories,colorMap,observed:observed.length});
  });
}
function filterIndices(points,state,track) {
  const result=[];
  for (let i=0;i<points.length;i++) {
    const p=points[i];
    if (p.cohort==='ancient' && !state.ancient || p.cohort==='modern' && !state.modern) continue;
    if (state.scope==='primary' && p.cohort==='ancient' && p.age_bp>state.primaryMax) continue;
    if (!(p.cohort==='modern' && state.pinModern) && (p.age_bp<state.young || p.age_bp>state.old)) continue;
    if (state.region!=='all' && category(p.region)!==state.region || state.source!=='all' && category(p.source)!==state.source) continue;
    if (track && track.kind==='categorical' && state.excludedCategories && state.excludedCategories.has(track.values[i] ?? null)) continue;
    result.push(i);
  }
  return result;
}
function sizes(indices,track,base) {
  if (!track || !track.range || track.range[0]===track.range[1]) return indices.map(()=>base);
  const [lo,hi]=track.range;
  return indices.map(i=>track.values[i]===null ? base : base*(.65+1.35*(track.values[i]-lo)/(hi-lo)));
}
function transformedAge(point,placement,layerYears,scale) {
  const age=placement==='layered' ? (point.cohort==='modern' ? 0 : (Math.floor(point.age_bp/layerYears)+.5)*layerYears) : point.age_bp;
  return scale==='linear' ? age/1000 : Math.log1p(age/1000)/Math.LN10;
}
function toTSV(data,indices) {
  const fields=[...new Set(data.points.flatMap(p=>Object.keys(p)))];
  const extra=[...(data.features||[]).map(f=>'feature:'+f.feature_id),...(data.variants||[]).map(v=>'variant:'+v.variant_id+':'+v.alt)];
  const clean=value=>value===null||value===undefined ? 'NA' : String(value).replace(/[\t\r\n]/g,' ');
  const lines=[fields.concat(extra).map(clean).join('\t')];
  for (const i of indices) {
    const values=fields.map(k=>data.points[i][k]);
    for (const f of data.features||[]) values.push(f.values[i]);
    for (const v of data.variants||[]) values.push(v.values===null?null:v.values[i]);
    lines.push(values.map(clean).join('\t'));
  }
  return lines.join('\n')+'\n';
}
const core={number,category,extent,makeTracks,filterIndices,sizes,transformedAge,toTSV};
if (typeof module !== 'undefined' && module.exports) module.exports=core;
if (typeof document === 'undefined') return;

const el=id=>document.getElementById(id), D=JSON.parse(el('allele3-data').textContent), P=D.points;
const esc=value=>String(value ?? '').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const tracks=makeTracks(D), trackMap=new Map(tracks.map(t=>[t.id,t]));
const ALL=P.map((_,i)=>i), primaryMax=D.display.max_age_bp, layerYears=D.display.layer_years;
const allMax=P.reduce((top,p)=>Math.max(top,p.age_bp),layerYears), roundedMax=Math.ceil(allMax/100)*100;
const byId=new Map(P.map((p,i)=>[p.sample_id,i]));
const selected=new Set();
let active=null, visible=[], excludedCategories=new Set(), pending=false, drawing=false, initialized=false;
let playTimer=null, playing=false, userCamera=null, activeTheme=D.display.theme || 'dark';
const cameras={
  perspective:{projection:{type:'orthographic'},eye:{x:1.55,y:-1.85,z:.9},up:{x:0,y:0,z:1}},
  top:{projection:{type:'orthographic'},eye:{x:0,y:0,z:2.5},up:{x:0,y:1,z:0}},
  side:{projection:{type:'orthographic'},eye:{x:2.7,y:0,z:.05},up:{x:0,y:0,z:1}}
};
userCamera=structuredClone(cameras.perspective);
const config={responsive:true,displaylogo:false,scrollZoom:true,toImageButtonOptions:{format:'png',scale:2},modeBarButtonsToRemove:['toImage','autoScale2d']};
function setOption(id,value,text){const o=document.createElement('option');o.value=value;o.textContent=text;el(id).appendChild(o);}
for(const track of tracks){setOption('colorBy',track.id,track.label);if(track.kind==='numeric')setOption('sizeBy',track.id,track.label);}
for(const [id,key] of [['regionFilter','region'],['sourceFilter','source']]){
  const values=[...new Set(P.map(p=>category(p[key])).filter(v=>v!==null))].sort((a,b)=>a<b?-1:a>b?1:0);
  for(const value of values)setOption(id,value,value);
}
el('colorBy').value=D.display.default_color || 'region';
el('agePlacement').value=D.display.default_view || 'continuous';
el('timeScale').value=D.display.time_scale;
el('pointSize').value=D.display.point_size || 2.7;
el('opacity').value=D.display.opacity || .75;
el('ageScope').options[0].textContent='0–'+(primaryMax/1000)+' ka BP + modern';
el('simulated').hidden=!D.simulated;
el('totalCount').textContent=P.length.toLocaleString();
el('olderCount').textContent=D.counts.ancient_older_than_primary.toLocaleString();
el('sourceDescription').textContent='Basis: '+D.basis_label+' · coordinate kind: '+D.coordinate_kind+'. Source status: '+D.provenance.plan_validation_status+'.';
el('sourceCounts').textContent=D.counts.ancient.toLocaleString()+' ancient + '+D.counts.modern.toLocaleString()+' modern people. Primary age range contains '+D.counts.ancient_primary.toLocaleString()+' ancient people. '+D.counts.ancient_older_than_primary.toLocaleString()+' older ancient people remain available in all ages and all-point downloads.';
el('sourceIdentity').textContent=Object.entries(D.provenance.inputs || {}).map(([name,info])=>name+': '+info.path+' · SHA256 '+info.sha256).join(' | ');
function scopeTop(){return el('ageScope').value==='primary'?primaryMax:Math.max(roundedMax,layerYears);}
function ageControls(young,old){
  const top=scopeTop(); young=Math.max(0,Math.min(top,number(young) ?? 0));old=Math.max(0,Math.min(top,number(old) ?? top));
  if(young>old)[young,old]=[old,young];
  for(const id of ['youngRange','oldRange','youngNumber','oldNumber'])el(id).max=top;
  el('youngRange').value=el('youngNumber').value=young;el('oldRange').value=el('oldNumber').value=old;
}
ageControls(0,primaryMax);
function theme(){
  document.body.dataset.theme=activeTheme;el('themeButton').textContent=activeTheme==='dark'?'Light theme':'Dark theme';
  return activeTheme==='dark'?{panel:'#101f2f',text:'#b8cddd',grid:'#263b50',muted:'#6c879d',accent:'#65d3d0',selected:'#ffd275'}:
    {panel:'#ffffff',text:'#33576e',grid:'#e2eaf0',muted:'#9cb0bf',accent:'#087e89',selected:'#a86000'};
}
function state(){return {
  ancient:el('ancientToggle').checked,modern:el('modernToggle').checked,pinModern:el('pinModern').checked,
  scope:el('ageScope').value,primaryMax,young:Number(el('youngNumber').value),old:Number(el('oldNumber').value),
  region:el('regionFilter').value,source:el('sourceFilter').value,excludedCategories,
  colorBy:el('colorBy').value,sizeBy:el('sizeBy').value,shapeBy:el('shapeBy').value,placement:el('agePlacement').value,
  timeScale:el('timeScale').value,pointSize:Number(el('pointSize').value),opacity:Number(el('opacity').value),
  context:el('contextToggle').checked,frames:el('framesToggle').checked,projection:el('projection').value,
  viewMode:el('viewMode').value,theme:activeTheme
};}
function hover(i,s){
  const p=P[i],track=trackMap.get(s.colorBy),value=track.values[i];
  return '<b>'+esc(p.sample_id)+'</b><br>'+esc(p.cohort)+' · '+esc(p.region || 'region unavailable')+
    '<br>Recorded age: '+p.age_bp.toLocaleString()+' BP<br>PC1: '+p.pc1+'<br>PC2: '+p.pc2+
    '<br>'+esc(track.label)+': '+(value===null?'unavailable':esc(value))+
    (s.placement==='layered'?'<br>Age placement uses bin midpoints; recorded date above.':'');
}
function symbol(i,s){
  if(s.shapeBy==='cohort')return P[i].cohort==='modern'?'diamond':'circle';
  const labels=trackMap.get('region').categories,index=labels.indexOf(category(P[i].region));
  return ['circle','diamond','square','cross','x'][Math.max(0,index)%5];
}
function axisForProjection(i,s){
  return s.projection==='pc1_pc2'?[P[i].pc1,P[i].pc2]:
    [s.projection==='pc1_age'?P[i].pc1:P[i].pc2,transformedAge(P[i],s.placement,layerYears,s.timeScale)];
}
function addPoints(traces,ids,s,type,withScale){
  const color=trackMap.get(s.colorBy),size=trackMap.get(s.sizeBy),numeric=color.kind==='numeric';
  let scaleUsed=false;
  for(const cohort of ['ancient','modern'])for(const unavailable of [false,true]){
    const subset=ids.filter(i=>P[i].cohort===cohort && (color.values[i]===null || color.values[i]===undefined)===unavailable);
    if(!subset.length)continue;
    const shape=subset.map(i=>symbol(i,s));
    const marker={size:sizes(subset,size,s.pointSize).map(n=>n*(type==='scattergl'?1.7:1)*(cohort==='modern'?1.35:1)),
      opacity:s.opacity,symbol:shape,line:{width:0},color:unavailable?GREY:numeric?subset.map(i=>color.values[i]):subset.map(i=>color.colorMap.get(color.values[i]) || GREY)};
    if(numeric && !unavailable && color.range){
      const range=color.variant?[0,2]:color.range;
      marker.cmin=range[0]===range[1]?range[0]-.5:range[0];marker.cmax=range[0]===range[1]?range[1]+.5:range[1];
      marker.colorscale=color.variant?[[0,'#4694d6'],[.5,'#e8e0c9'],[1,'#f59d53']]:'Viridis';
      marker.showscale=withScale && !scaleUsed;scaleUsed=true;
      marker.colorbar={title:{text:color.units,font:{size:10}},thickness:10,len:.55,x:.96,tickfont:{size:10},outlinewidth:0};
    }
    const xy=type==='scattergl'?subset.map(i=>axisForProjection(i,s)):null;
    traces.push({type,mode:'markers',name:cohort+(unavailable?' · unavailable':''),showlegend:false,
      x:subset.map((i,j)=>xy?xy[j][0]:P[i].pc1),y:subset.map((i,j)=>xy?xy[j][1]:P[i].pc2),
      ...(type==='scatter3d'?{z:subset.map(i=>transformedAge(P[i],s.placement,layerYears,s.timeScale))}:{}),
      marker,customdata:subset,text:subset.map(i=>hover(i,s)),hoverinfo:'text'});
  }
}
function selectedTrace(ids,s,type,t){
  const a=ids.filter(i=>selected.has(i));if(!a.length)return null;
  const xy=type==='scattergl'?a.map(i=>axisForProjection(i,s)):null;
  return {type,mode:'markers',showlegend:false,name:'Selection',x:a.map((i,j)=>xy?xy[j][0]:P[i].pc1),y:a.map((i,j)=>xy?xy[j][1]:P[i].pc2),
    ...(type==='scatter3d'?{z:a.map(i=>transformedAge(P[i],s.placement,layerYears,s.timeScale))}:{}),
    customdata:a,text:a.map(i=>hover(i,s)),hoverinfo:'text',marker:{size:sizes(a,trackMap.get(s.sizeBy),s.pointSize).map(n=>n*(type==='scattergl'?2.2:1.5)+2),
      symbol:'circle-open',color:t.selected,opacity:1,line:{width:1.5,color:t.selected}}};
}
function plots(ids,s,t){
  const a=[],b=[];
  // Pinned references keep supplied dates even if an unusual modern record is
  // older than the primary range; a counted visible point must not be clipped.
  const referenceTop=ids.reduce((top,i)=>P[i].cohort==='modern'?Math.max(top,P[i].age_bp):top,0);
  const topBP=Math.max(scopeTop(),referenceTop)+(s.placement==='layered'?layerYears:0);
  const transform=bp=>s.timeScale==='linear'?bp/1000:Math.log1p(bp/1000)/Math.LN10;
  const ztop=transform(Math.max(topBP,layerYears));
  if(s.frames){
    const step=Math.max(layerYears,Math.ceil(topBP/(30*layerYears))*layerYears);
    for(let age=0;age<=topBP;age+=step){const z=transform(age);a.push({type:'scatter3d',mode:'lines',showlegend:false,hoverinfo:'skip',
      x:[D.limits.x[0],D.limits.x[1],D.limits.x[1]],y:[D.limits.y[0],D.limits.y[0],D.limits.y[1]],z:[z,z,z],
      line:{color:t.grid,width:age===0?2:1},opacity:.7});}
  }
  if(s.context){
    const shown=new Set(ids),context=ALL.filter(i=>!shown.has(i)&&(s.scope==='all'||P[i].cohort==='modern'||P[i].age_bp<=primaryMax));
    a.push({type:'scatter3d',mode:'markers',showlegend:false,hoverinfo:'skip',x:context.map(i=>P[i].pc1),y:context.map(i=>P[i].pc2),
      z:context.map(i=>transformedAge(P[i],s.placement,layerYears,s.timeScale)),marker:{size:Math.max(1,s.pointSize*.7),color:t.muted,opacity:.12}});
    const xy=context.map(i=>axisForProjection(i,s));b.push({type:'scattergl',mode:'markers',showlegend:false,hoverinfo:'skip',x:xy.map(v=>v[0]),y:xy.map(v=>v[1]),marker:{size:Math.max(1,s.pointSize),color:t.muted,opacity:.12}});
  }
  addPoints(a,ids,s,'scatter3d',true);addPoints(b,ids,s,'scattergl',false);
  const sa=selectedTrace(ids,s,'scatter3d',t),sb=selectedTrace(ids,s,'scattergl',t);if(sa)a.push(sa);if(sb)b.push(sb);
  let tickAges=[];
  if(s.timeScale==='linear'){const step=topBP<=12000?2000:Math.ceil(topBP/6000)*1000;for(let x=0;x<=topBP;x+=step)tickAges.push(x);}
  else {tickAges=[0,500,1000,2000,4000,6000,10000,15000,20000,30000,50000,100000,200000,500000,1000000].filter(x=>x<=topBP);if(tickAges.length<2)tickAges=[0,topBP];}
  const axis=(title,range)=>({title:{text:title,font:{size:11}},range,autorange:false,showbackground:false,gridcolor:t.grid,zerolinecolor:t.grid,color:t.text,tickfont:{size:10}});
  const zaxis=Object.assign(axis('Age (ka BP)'+(s.timeScale==='log1p'?' · log1p spacing':''),[-.025*ztop,1.04*ztop]),{tickmode:'array',tickvals:tickAges.map(transform),ticktext:tickAges.map(x=>x/1000)});
  const common={paper_bgcolor:t.panel,plot_bgcolor:t.panel,font:{family:'system-ui',color:t.text},showlegend:false};
  const three=Object.assign({},common,{margin:{l:0,r:0,t:12,b:0},uirevision:'fixed-pca-3d',
    scene:{xaxis:axis('PC1',D.limits.x),yaxis:axis('PC2',D.limits.y),zaxis,aspectmode:'manual',aspectratio:{x:1,y:1,z:1.45},camera:structuredClone(userCamera)}});
  const two=Object.assign({},common,{margin:{l:48,r:15,t:15,b:42},dragmode:'lasso',uirevision:'fixed-pca-2d-'+s.projection+'-'+s.timeScale+'-'+s.placement,
    xaxis:axis(s.projection==='pc2_age'?'PC2':'PC1',s.projection==='pc2_age'?D.limits.y:D.limits.x),
    yaxis:s.projection==='pc1_pc2'?Object.assign(axis('PC2',D.limits.y),{scaleanchor:'x',scaleratio:1}):zaxis,
    hovermode:'closest',clickmode:'event'});
  return {a,b,three,two};
}
function legend(track){
  el('legend').replaceChildren();
  if(track.kind==='categorical'){
    const values=track.categories.slice();if(track.observed<P.length)values.push(null);
    for(const value of values){
      const button=document.createElement('button'),dot=document.createElement('i'),text=document.createElement('span');
      dot.style.background=track.colorMap.get(value)||GREY;text.textContent=value===null?'Unavailable':value;
      button.append(dot,text);button.classList.toggle('off',excludedCategories.has(value));button.setAttribute('aria-pressed',String(!excludedCategories.has(value)));
      button.onclick=()=>{if(excludedCategories.has(value))excludedCategories.delete(value);else excludedCategories.add(value);requestDraw();};
      el('legend').appendChild(button);
    }
    el('legendText').textContent='Click a category to toggle it. Colors are fixed across filters.';
  }else{
    const range=track.variant&&track.range?[0,2]:track.range;
    el('legendText').textContent=(range?'Fixed color range: '+range[0]+' to '+range[1]+' '+track.units+'. ':'No observed values. ')+'Grey = unavailable.';
  }
}
function inspector(){
  const visibleSet=new Set(visible),nVisible=[...selected].filter(i=>visibleSet.has(i)).length;
  el('selectedCount').textContent=selected.size.toLocaleString();el('selectionVisible').textContent=selected.size?nVisible+' visible':'';
  el('selectionStatus').textContent=selected.size.toLocaleString()+' selected; '+nVisible.toLocaleString()+' in current view. Exports retain exact source values.';
  el('exportSelected').disabled=!selected.size;el('clearSelection').disabled=!selected.size;
  if(active===null)return;
  const p=P[active],track=trackMap.get(el('colorBy').value),info=[['PC1',p.pc1],['PC2',p.pc2],['Recorded age BP',p.age_bp],['Region',p.region],['Source',p.source],['Visible now',visibleSet.has(active)?'Yes':'Outside filters'],[track.label,track.values[active]]];
  for(const [k,v] of Object.entries(p))if(!['sample_id','pc1','pc2','age_bp','region','source','cohort'].includes(k))info.push([k,v]);
  el('sampleDetail').className='';el('sampleDetail').innerHTML='<div class="detail-id">'+esc(p.sample_id)+'</div><span class="pill">'+esc(p.cohort)+'</span>'+(p.region?'<span class="pill">'+esc(p.region)+'</span>':'')+
    '<dl class="detail-grid">'+info.map(([k,v])=>'<dt>'+esc(k)+'</dt><dd>'+esc(v===null||v===undefined?'Unavailable':v)+'</dd>').join('')+'</dl>';
}
function histogram(ids,t){
  const binWidth=Math.max(layerYears,Math.ceil(roundedMax/(40*layerYears))*layerYears);
  const n=Math.max(1,Math.ceil((roundedMax+1)/binWidth)),all=Array(n).fill(0),shown=Array(n).fill(0);
  for(const p of P)if(p.cohort==='ancient')all[Math.min(n-1,Math.floor(p.age_bp/binWidth))]++;
  for(const i of ids)if(P[i].cohort==='ancient')shown[Math.min(n-1,Math.floor(P[i].age_bp/binWidth))]++;
  const max=Math.max(...all,1),w=760/n,gap=Math.min(3,w*.2);
  let svg='<svg viewBox="0 0 760 86" role="img" aria-label="Original and filtered ancient sample counts by age interval">';
  for(let k=0;k<n;k++){
    const x=k*w+gap/2,bh=all[k]/max*65,sh=shown[k]/max*65;
    svg+='<g data-bin="'+k+'" role="button" tabindex="0" aria-label="'+k*binWidth+' to '+(k+1)*binWidth+' BP; '+shown[k]+' filtered of '+all[k]+' original"><title>'+k*binWidth+'–'+(k+1)*binWidth+' BP: '+shown[k]+' / '+all[k]+' people</title>'+
      '<rect x="'+x+'" y="'+(70-bh)+'" width="'+Math.max(.2,w-gap)+'" height="'+Math.max(2,bh)+'" rx="1.5" fill="'+t.grid+'"/>'+
      '<rect x="'+x+'" y="'+(70-sh)+'" width="'+Math.max(.2,w-gap)+'" height="'+sh+'" rx="1.5" fill="'+t.accent+'"/>'+
      '<rect x="'+(k*w)+'" y="0" width="'+w+'" height="76" fill="transparent"/></g>';
  }
  svg+='</svg>';el('histogram').innerHTML=svg;
  for(const bar of el('histogram').querySelectorAll('[data-bin]')){
    const choose=()=>{stopPlay();const k=Number(bar.dataset.bin);if((k+1)*binWidth>scopeTop())el('ageScope').value='all';ageControls(k*binWidth,Math.min((k+1)*binWidth,scopeTop()));requestDraw();};
    bar.onclick=choose;bar.onkeydown=e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();choose();}};
  }
}
function counts(ids,s,t){
  const ancient=ids.filter(i=>P[i].cohort==='ancient').length;
  el('visibleCount').textContent=ids.length.toLocaleString();el('viewCohorts').textContent=ancient.toLocaleString()+' ancient / '+(ids.length-ancient).toLocaleString()+' modern';
  const text=(s.young/1000)+'–'+(s.old/1000)+' ka BP'+(s.pinModern&&s.modern?' + modern references':'');
  el('activeCaption').textContent=text;el('timeStatus').textContent=text;
  const track=trackMap.get(s.colorBy);el('colorCaption').textContent=track.label;
  el('featureDescription').textContent=track.provenance+' · '+track.observed.toLocaleString()+'/'+P.length.toLocaleString()+' available.';
  el('pointSizeOutput').textContent=s.pointSize.toFixed(1);el('opacityOutput').textContent=Math.round(s.opacity*100)+'%';
  el('pngButton').textContent=s.viewMode==='2d'?'Export 2D PNG':'Export 3D PNG';
  legend(track);histogram(ids,t);inspector();
}
function pick(i,add){
  if(!Number.isInteger(i)||i<0||i>=P.length)return;
  if(!add)selected.clear();selected.add(i);active=i;requestDraw();
}
function attachEvents(){
  for(const id of ['plot3d','plot2d']){
    el(id).on('plotly_click',event=>{const pt=event.points&&event.points.find(p=>Number.isInteger(p.customdata));if(pt)pick(pt.customdata,Boolean(event.event&&event.event.shiftKey));});
  }
  el('plot2d').on('plotly_selected',event=>{
    if(!event||!Array.isArray(event.points))return;
    selected.clear();for(const p of event.points)if(Number.isInteger(p.customdata))selected.add(p.customdata);
    active=selected.size?[...selected][0]:null;requestDraw();
  });
  el('plot3d').on('plotly_relayout',event=>{if(event['scene.camera'])userCamera=structuredClone(event['scene.camera']);});
  el('plot3d').on('plotly_webglcontextlost',()=>{el('plotError').hidden=false;el('plotError').textContent='The browser lost its WebGL context. Reload this offline file to restore the plots; all data remain in the downloads.';});
}
async function requestDraw(){
  pending=true;if(drawing)return;drawing=true;
  try{
    while(pending){
      pending=false;const s=state(),t=theme();el('workspace').dataset.view=s.viewMode;
      visible=filterIndices(P,s,trackMap.get(s.colorBy));counts(visible,s,t);
      const q=plots(visible,s,t);
      await Plotly.react(el('plot3d'),q.a,q.three,config);
      await Plotly.react(el('plot2d'),q.b,q.two,config);
      if(!initialized){attachEvents();initialized=true;}
      el('loading').hidden=true;el('plotError').hidden=true;
    }
  }catch(error){
    el('plotError').hidden=false;el('plotError').textContent='Plot unavailable: '+String(error.message||error);el('loading').hidden=true;
    console.error(error);
  }finally{drawing=false;}
}
function stopPlay(){playing=false;if(playTimer!==null)clearTimeout(playTimer);playTimer=null;el('playButton').textContent='▶ Play slices';}
async function advance(){
  if(!playing)return;const width=Number(el('sliceWidth').value),top=scopeTop(),current=Number(el('youngNumber').value);
  const young=current+width>=top?0:current+width;ageControls(young,Math.min(young+width,top));
  await requestDraw();if(playing)playTimer=setTimeout(advance,900);
}
function download(content,name,type){
  const url=URL.createObjectURL(new Blob([content],{type})),a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
}
function exportTSV(ids,name){download(toTSV(D,ids),name,'text/tab-separated-values;charset=utf-8');}
async function exportPNG(which){
  try{await Plotly.downloadImage(el(which),{format:'png',filename:'allele3_'+which+'_fixed_coordinates',width:1600,height:1200,scale:1.5});}
  catch(error){el('plotError').hidden=false;el('plotError').textContent='PNG export unavailable: '+String(error.message||error);}
}
function search(){
  const term=el('sampleSearch').value.trim().toLowerCase();el('searchResults').replaceChildren();if(!term)return;
  const exact=byId.has(el('sampleSearch').value.trim())?[byId.get(el('sampleSearch').value.trim())]:[];
  const matches=exact.length?exact:ALL.filter(i=>P[i].sample_id.toLowerCase().includes(term)).slice(0,8);
  for(const i of matches){const button=document.createElement('button');button.textContent=P[i].sample_id+' · '+P[i].age_bp+' BP';button.onclick=()=>{pick(i,false);el('sampleSearch').value=P[i].sample_id;el('searchResults').replaceChildren();};el('searchResults').appendChild(button);}
  if(!matches.length){const div=document.createElement('div');div.className='hint';div.textContent='No source sample matches this ID.';el('searchResults').appendChild(div);}
}
el('themeButton').onclick=()=>{activeTheme=activeTheme==='dark'?'light':'dark';requestDraw();};
el('colorBy').onchange=()=>{excludedCategories=new Set();requestDraw();};
for(const id of ['sizeBy','shapeBy','ancientToggle','modernToggle','regionFilter','sourceFilter','pinModern','agePlacement','timeScale','contextToggle','framesToggle','projection'])el(id).onchange=()=>requestDraw();
for(const id of ['pointSize','opacity'])el(id).oninput=()=>requestDraw();
el('ageScope').onchange=()=>{stopPlay();ageControls(0,scopeTop());requestDraw();};
for(const [id,partner] of [['youngRange','oldNumber'],['oldRange','youngNumber'],['youngNumber','oldNumber'],['oldNumber','youngNumber']]){
  el(id).addEventListener(id.endsWith('Range')?'input':'change',()=>{stopPlay();if(id.startsWith('young'))ageControls(el(id).value,el(partner).value);else ageControls(el(partner).value,el(id).value);requestDraw();});
}
el('viewMode').onchange=()=>{requestDraw();setTimeout(()=>{Plotly.Plots.resize(el('plot3d'));Plotly.Plots.resize(el('plot2d'));},50);};
for(const button of document.querySelectorAll('[data-camera]'))button.onclick=()=>{userCamera=structuredClone(cameras[button.dataset.camera]);requestDraw();};
el('playButton').onclick=()=>{if(playing){stopPlay();return;}playing=true;el('playButton').textContent='❚❚ Pause';const width=Number(el('sliceWidth').value);ageControls(0,Math.min(width,scopeTop()));requestDraw().then(()=>{if(playing)playTimer=setTimeout(advance,900);});};
el('sliceWidth').onchange=()=>{stopPlay();};
el('sampleSearch').oninput=search;
el('sampleSearch').onkeydown=event=>{if(event.key==='Enter'){event.preventDefault();const first=el('searchResults').querySelector('button');if(first)first.click();}};
el('clearSelection').onclick=()=>{selected.clear();active=null;el('sampleDetail').className='empty-detail';el('sampleDetail').textContent='Selection cleared. Click a point, search a sample ID, or select in the 2D projection.';requestDraw();};
el('exportSelected').onclick=()=>exportTSV([...selected].sort((a,b)=>a-b),'allele3_selected_people.tsv');
el('exportFiltered').onclick=()=>exportTSV(visible,'allele3_filtered_people.tsv');
el('exportAll').onclick=()=>exportTSV(ALL,'allele3_all_original_people.tsv');
el('pngButton').onclick=()=>exportPNG(el('viewMode').value==='2d'?'plot2d':'plot3d');
el('export2d').onclick=()=>exportPNG('plot2d');
el('exportState').onclick=()=>{const s=state();s.excludedCategories=[...s.excludedCategories];download(JSON.stringify({schema:'allele3-explorer-view-1.0',created_utc:new Date().toISOString(),source_provenance:D.provenance,source_counts:D.counts,display:s,camera:userCamera,selected_sample_ids:[...selected].map(i=>P[i].sample_id),coordinates_changed:false},null,2)+'\n','allele3_view_settings.json','application/json');};
el('resetButton').onclick=()=>{
  stopPlay();selected.clear();active=null;excludedCategories=new Set();el('ageScope').value='primary';ageControls(0,primaryMax);
  for(const id of ['ancientToggle','modernToggle','pinModern','framesToggle'])el(id).checked=true;
  el('contextToggle').checked=false;el('regionFilter').value=el('sourceFilter').value='all';el('colorBy').value=D.display.default_color||'region';el('sizeBy').value='constant';el('shapeBy').value='cohort';
  el('agePlacement').value='continuous';el('timeScale').value=D.display.time_scale;el('viewMode').value='both';el('projection').value='pc1_pc2';
  el('pointSize').value=D.display.point_size;el('opacity').value=D.display.opacity;el('sampleSearch').value='';el('searchResults').replaceChildren();
  el('sampleDetail').className='empty-detail';el('sampleDetail').textContent='Click a point or use the linked 2D selection to inspect the original records.';
  userCamera=structuredClone(cameras.perspective);requestDraw();
};
window.__allele3Explorer={version:'2.0',getState:()=>Object.assign({},state(),{excludedCategories:[...excludedCategories]}),getVisibleSampleIDs:()=>visible.map(i=>P[i].sample_id),getSelectedSampleIDs:()=>[...selected].map(i=>P[i].sample_id),sourcePointCount:P.length,isReady:()=>initialized&&!drawing};
theme();requestDraw();
}());
