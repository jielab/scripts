"""Manuscript figures and matched Excel source data from frozen panome results.
No model fitting, state selection, or outcome-driven feature screening is performed here.
"""
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from sksurv.nonparametric import kaplan_meier_estimator

MODEL_LABELS={
    'clinical':'Clinical', 'clinical_pwas':'Clinical + PWAS score',
    'clinical_pca':'Clinical + PCA', 'clinical_ae':'Clinical + AE',
    'clinical_state':'Clinical + hard state', 'clinical_soft_state':'Clinical + soft state',
    'clinical_diffusion':'Clinical + diffusion', 'clinical_panome':'Clinical + AE + soft state',
    'clinical_elasticnet':'Clinical + elastic-net'}
MODEL_COLORS={m:plt.get_cmap('tab10')(i) for i,m in enumerate(MODEL_LABELS)}
DPI=400


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def numeric_table(path):
    try:return pd.read_csv(path)
    except pd.errors.EmptyDataError:return pd.DataFrame()


def workbook(path,tables,caption,metadata):
    wb=Workbook();wb.remove(wb.active)
    data={'Readme':pd.DataFrame({'item':['Caption',*metadata.keys()], 'value':[caption,*[str(v) for v in metadata.values()]]}),**tables}
    for name,frame in data.items():
        if len(frame)>1048575 or len(frame.columns)>16384:raise ValueError(f'Excel limits exceeded: {name}')
        ws=wb.create_sheet(name[:31]);ws.append([str(c) for c in frame.columns])
        for row_number,row in enumerate(frame.itertuples(index=False,name=None),start=2):
            values=[]
            for v in row:
                if v is None or pd.isna(v):v=None
                elif isinstance(v,(float,np.floating)) and not np.isfinite(v):v=str(v)
                elif isinstance(v,np.generic):v=v.item()
                values.append(v)
            ws.append(values)
            for column_number,value in enumerate(values,start=1):
                if isinstance(value,str) and value.startswith('='):ws.cell(row_number,column_number).data_type='s'
        ws.freeze_panes='B2';ws.auto_filter.ref=ws.dimensions
        for cell in ws[1]:
            cell.fill=PatternFill('solid',fgColor='173F50');cell.font=Font(color='FFFFFF',bold=True)
            cell.alignment=Alignment(wrap_text=True)
        for col in ws.columns:
            first=col[0];ws.column_dimensions[first.column_letter].width=min(40,max(14,len(str(first.value))+3))
        for row in ws.iter_rows(min_row=2):
            for cell in row:
                if isinstance(cell.value,float):cell.number_format='0.0000'
        if name=='Readme':
            ws.column_dimensions['A'].width=28;ws.column_dimensions['B'].width=110
            ws['B2'].alignment=Alignment(wrap_text=True,vertical='top');ws.row_dimensions[2].height=110
    wb.save(path)


def panel(ax,label,title):
    ax.set_title(f'{label}  {title}',loc='left',fontweight='bold',pad=12)
    ax.spines[['top','right']].set_visible(False)


def unavailable(ax,text):
    ax.text(.5,.5,text,ha='center',va='center',transform=ax.transAxes,wrap=True);ax.set_axis_off()


class Publication:
    def __init__(self,root,out):
        self.root=root;self.out=out;self.items=[]
        self.manifest=json.loads((root/'manifest.json').read_text());self.cfg=self.manifest['config']
        self.audit=json.loads((root/'s2_preprocess/cohort_audit.json').read_text())
        self.graph=json.loads((root/'s4_graph/summary.json').read_text())
        self.people=pd.read_csv(root/'s6_report/person_molecular_states.csv',dtype={self.cfg['id_col']:str})
        self.states=sorted(self.people.state.unique())
        self.colors={s:plt.get_cmap('tab10' if len(self.states)<=10 else 'tab20')(i) for i,s in enumerate(self.states)}
        self.counts=numeric_table(root/'s6_report/state_counts.csv')
        self.metrics=numeric_table(root/'s5_predict/test_metrics.csv')
        self.meta={'Run':str(root),'Trait':self.cfg['trait'],'Omics':self.cfg['biom'],
                   'Analysis signature':self.manifest['signature'],'Figure DPI':DPI,
                   'Study type':'Synthetic test' if self.cfg['demo'] else ('Pilot' if self.cfg['max_samples'] else 'Full input cohort'),
                   'Scope':'Internal split validation; baseline-free for the configured endpoint, not necessarily all diseases.',
                   'Risk interpretation':'Death-censored net risk; not competing-risk cumulative incidence.',
                   'States':'Exploratory communities; graph edge-dropout stability is not full-pipeline replication.'}

    def read(self,relative):return numeric_table(self.root/relative)

    def save(self,stem,fig,tables,caption):
        fig.text(.01,.006,f"Panome | {self.cfg['biom']} | {self.cfg['trait']} | {self.cfg['run_name']}",fontsize=8,color='#555555')
        fig.savefig(self.out/(stem+'.png'),dpi=DPI,facecolor='white',bbox_inches='tight')
        plt.close(fig)
        workbook(self.out/(stem+'.xlsx'),tables,caption,self.meta)
        self.items.append(dict(figure=stem,png=stem+'.png',xlsx=stem+'.xlsx',caption=caption))
        print(f'Final: {stem}.png + {stem}.xlsx',flush=True)

    def cohort(self):
        a=self.audit;history=self.read('s3_representation/ae_history.csv')
        split=pd.DataFrame(a['split_summary']).T.reindex(['train','validation','test']).rename_axis('split').reset_index().rename(columns={'size':'n','sum':'events'})
        rows=[('Input participants',a['joined']),('Baseline endpoint disease',a['prevalent']),
              ('Invalid follow-up',a['invalid_followup']),('Eligible before molecular QC',a['eligible']),
              ('Excluded: molecular missingness',a['excluded_sample_missing']),('Final cohort',a['post_qc_n'])]
        flow=pd.DataFrame(rows,columns=['criterion','n'])
        fig,axes=plt.subplots(1,2,figsize=(13,7),gridspec_kw={'width_ratios':[1.05,1]},layout='constrained')
        ax=axes[0];panel(ax,'A','Cohort assembly and frozen split');ax.set_axis_off()
        boxes=[f"Input: {a['joined']:,} participants",f"Baseline eligibility: {a['eligible']:,}",
               f"Molecular QC: {a['post_qc_n']:,} participants\n{a['retained_features']:,} retained features"]
        for i,text in enumerate(boxes):
            y=.84-i*.24
            ax.add_patch(FancyBboxPatch((.07,y-.08),.86,.15,boxstyle='round,pad=0.012',facecolor='#E7F1F4',edgecolor='#347E91',transform=ax.transAxes))
            ax.text(.5,y,text,ha='center',va='center',transform=ax.transAxes,fontsize=12)
            if i<2:ax.annotate('',xy=(.5,y-.17),xytext=(.5,y-.09),xycoords='axes fraction',arrowprops={'arrowstyle':'->','color':'#347E91'})
        labels=[f"{r.split.title()}: n={r.n:,}, events={r.events:,}" for r in split.itertuples()]
        ax.text(.5,.08,'\n'.join(labels),ha='center',va='center',transform=ax.transAxes,fontsize=11)
        ax=axes[1];panel(ax,'B','Autoencoder reconstruction')
        for col,label in [('train_mse','Training'),('validation_mse','Validation')]:ax.plot(history.epoch,history[col],label=label,lw=2)
        best=history.loc[history.validation_mse.idxmin()]
        ax.axvline(best.epoch,color='#777777',ls='--',lw=1,label=f'Best epoch: {int(best.epoch)}')
        ax.set(xlabel='Epoch',ylabel='Observed-entry mean squared error');ax.legend(frameon=False)
        self.save('Fig1.cohort_representation',fig,{'cohort_flow':flow,'cohort_audit':pd.DataFrame({'metric':list(a),'value':[json.dumps(v) for v in a.values()]}),'split_counts':split,'ae_training':history},
                  'A, cohort assembly under configured baseline disease and molecular QC rules; exclusion categories in the audit can overlap. B, observed-entry AE reconstruction losses. Early stopping uses validation reconstruction only. Disease outcomes do not train the AE.')

    def landscape(self):
        profiles=self.read('s6_report/training_state_profiles.csv').set_index('state')
        features=list(dict.fromkeys(f for _,row in profiles.iterrows() for f in row.abs().nlargest(5).index))
        plotted=profiles[features];fig=plt.figure(figsize=(14,10),layout='constrained');grid=fig.add_gridspec(2,2,height_ratios=[1,1.1])
        tables={};rng=np.random.default_rng(self.cfg['seed'])
        limits={c:(self.people[c].min(),self.people[c].max()) for c in ['AE1','AE2']}
        for j,split in enumerate(['train','test']):
            ax=fig.add_subplot(grid[0,j]);panel(ax,chr(65+j),f'{split.title()} molecular coordinates')
            d=self.people[self.people.split==split]
            if len(d)>8000:d=d.iloc[np.sort(rng.choice(len(d),8000,replace=False))]
            for s in self.states:
                q=d[d.state==s];ax.scatter(q.AE1,q.AE2,s=4,alpha=.55,color=self.colors[s],label=f'S{s}',rasterized=True)
            ax.set(xlabel='AE coordinate 1',ylabel='AE coordinate 2',xlim=limits['AE1'],ylim=limits['AE2']);ax.legend(frameon=False,ncol=4,markerscale=2,fontsize=8)
            tables[split+'_plotted_points']=d[['AE1','AE2','state','split']].reset_index(drop=True).rename_axis('plot_row').reset_index()
        ax=fig.add_subplot(grid[1,:]);panel(ax,'C','Training molecular profiles of candidate states')
        im=ax.imshow(plotted.values,cmap='RdBu_r',vmin=-2,vmax=2,aspect='auto')
        ax.set_xticks(range(len(features)),features,rotation=60,ha='right',fontsize=8);ax.set_yticks(range(len(profiles)),profiles.index)
        fig.colorbar(im,ax=ax,shrink=.7,label='Mean standardized molecular value')
        tables.update(plotted_heatmap=plotted.reset_index(),all_state_profiles=profiles.T.rename_axis('feature').reset_index(),state_counts=self.counts)
        self.save('Fig2.molecular_landscape',fig,tables,
                  'A–B, first two AE coordinates; up to 8,000 people per split are sampled using the recorded seed. The two coordinates are not a complete map of all latent dimensions. State colors are consistent across figures. C, union of five largest absolute training-profile means per state. Values are descriptive, not independent differential-expression tests; color scale saturates at ±2 SD. States remain exploratory.')

    def survival(self):
        test=self.people[self.people.split=='test'];tau=min(max(self.cfg['horizons']),test.time.max())
        fig=plt.figure(figsize=(14,8),layout='constrained');grid=fig.add_gridspec(2,2,height_ratios=[3,1])
        ax=fig.add_subplot(grid[0,0]);panel(ax,'A',f'Test-set {"CAD" if self.cfg["trait"]=="cvd_cad" else self.cfg["trait"]} net risk')
        curves=[];risks=[];ticks=np.linspace(0,tau,6)
        for s in self.states:
            d=test[test.state==s]
            if d.empty:continue
            t,surv,ci=kaplan_meier_estimator(d.event.astype(bool),d.time,conf_type='log-log')
            t=np.r_[0,t];surv=np.r_[1,surv];ci=np.c_[np.ones(2),ci]
            keep=t<=tau
            # Include the right boundary of the plotted step curve explicitly.
            end=np.searchsorted(t,tau,side='right')-1
            tt=np.r_[t[keep],tau];ss=np.r_[surv[keep],surv[end]]
            low=np.r_[ci[0,keep],ci[0,end]];high=np.r_[ci[1,keep],ci[1,end]]
            ax.step(tt,1-ss,where='post',color=self.colors[s],label=f'S{s}',lw=1.6)
            ax.fill_between(tt,1-high,1-low,step='post',color=self.colors[s],alpha=.08)
            curves.append(pd.DataFrame({'state':s,'time':tt,'survival':ss,'lower95':low,'upper95':high,'net_risk':1-ss,'risk_lower95':1-high,'risk_upper95':1-low}))
            risks.extend(dict(state=s,time=float(t0),at_risk=int((d.time>=t0).sum())) for t0 in ticks)
        ax.set(xlim=(0,tau),xlabel='Years since baseline',ylabel='Net risk (1 − Kaplan–Meier survival)',ylim=(0,min(1,max(c.risk_upper95.max() for c in curves)*1.1)));ax.legend(frameon=False,ncol=4)
        axr=fig.add_subplot(grid[1,0]);axr.set_axis_off()
        risk=pd.DataFrame(risks).pivot(index='state',columns='time',values='at_risk')
        tab=axr.table(cellText=risk.values,rowLabels=[f'S{s}' for s in risk.index],colLabels=[f'{t:g} y' for t in ticks],loc='center',cellLoc='center')
        tab.auto_set_font_size(False);tab.set_fontsize(9);tab.scale(1,1.12);axr.set_title('Number at risk',loc='left',fontsize=10)
        ax=fig.add_subplot(grid[:,1]);panel(ax,'B','Exploratory adjusted state associations')
        primary=self.read('s5_predict/test_state_cox.csv');hr=primary[primary.term.str.match(r'^state\d+$')].copy()
        hr['state']=hr.term.str.extract(r'(\d+)$').astype(int)
        ref=next((s for s in sorted(test.state.unique()) if s not in hr.state.values),None)
        for i,s in enumerate(self.states):
            row=hr[hr.state==s]
            if s==ref:ax.scatter(1,i,color=self.colors[s],marker='s');ax.text(1.04,i,'Reference',va='center',fontsize=9);continue
            if row.empty:continue
            row=row.iloc[0];v,l,u=[row[c] for c in ['exp(coef)','lower .95','upper .95']]
            if np.isfinite([v,l,u]).all() and min(v,l,u)>0:
                ax.errorbar(v,i,xerr=[[v-l],[u-v]],fmt='o',color=self.colors[s],capsize=3)
            else:ax.text(1,i,'Not estimable',va='center',fontsize=9)
        ax.axvline(1,color='#888888',ls='--',lw=1);ax.set(xscale='log',xlabel='Adjusted hazard ratio (95% CI)')
        ax.set_yticks(range(len(self.states)),[f'S{s} (events={int(test.loc[test.state==s,"event"].sum())})' for s in self.states]);ax.invert_yaxis()
        tables={'km_curves':pd.concat(curves,ignore_index=True),'numbers_at_risk':pd.DataFrame(risks),'state_counts':self.counts,'state_cox':hr,'all_cox_terms':primary,'ph_diagnostics':self.read('s5_predict/test_state_ph_diagnostics.csv')}
        optional=self.root/'validation/center_stratified_state_cox.csv'
        if optional.exists():tables['center_stratified_sensitivity']=numeric_table(optional)
        self.save('Fig3.state_survival',fig,tables,
                  f'A, test-set net risk (1 minus Kaplan–Meier survival) with transformed pointwise 95% log-log intervals to {tau:g} years and numbers at risk. Death is censored. B, frozen training-defined state contrasts from exploratory clinically adjusted test-set Cox, 95% Wald intervals; reference is the baseline factor level. Nonfinite intervals are not plotted. Sparse center covariates can separate (inspect all_cox_terms); any saved center-stratified post-hoc sensitivity is provided separately. Neither associations nor nominal significance establishes reproducible biological subtypes.')

    def prediction(self):
        c=self.metrics[self.metrics.metric=='Harrell_C'].set_index('model').reindex(MODEL_LABELS).dropna(subset=['value']).reset_index()
        ci=self.read('s5_predict/test_c_intervals.csv') if (self.root/'s5_predict/test_c_intervals.csv').exists() else pd.DataFrame(columns=['model'])
        c=c.merge(ci,on='model',how='left');base=float(c.loc[c.model=='clinical','value'].iloc[0]);c['delta']=c.value-base
        fig,axes=plt.subplots(1,2,figsize=(14,6.8),layout='constrained')
        for j,(value,lo,hi,title) in enumerate([('value','c_0.025','c_0.975','Independent test discrimination'),('delta','delta_vs_clinical_0.025','delta_vs_clinical_0.975','Increment beyond clinical baseline')]):
            ax=axes[j];panel(ax,chr(65+j),title)
            for i,row in c.iterrows():
                color=MODEL_COLORS[row.model];ax.scatter(row[value],i,color=color,s=38,zorder=3)
                if lo in c and np.isfinite([row[lo],row[hi]]).all():ax.hlines(i,row[lo],row[hi],color=color,lw=2)
            ax.set_yticks(range(len(c)),[MODEL_LABELS[m] for m in c.model]);ax.invert_yaxis()
            ax.set_xlabel('Harrell C (95% bootstrap interval)' if j==0 else 'Paired change in Harrell C (95% interval)')
            if j==1:ax.axvline(0,color='#999999',ls='--',lw=1)
            ax.grid(axis='x',alpha=.15)
        self.save('Fig4.prediction_benchmarks',fig,{'plotted_C_and_delta':c,'all_test_metrics':self.metrics,'validation_tuning':self.read('s5_predict/validation_tuning.csv')},
                  'A, Harrell C for nine models evaluated on the same independent test set. B, paired C difference relative to the clinical baseline. Intervals are 2.5th–97.5th percentiles from the stored test-person bootstrap, conditional on the fitted models and this split; they do not capture model retraining uncertainty. Hyperparameters were chosen in validation. Model order is prespecified, not sorted by test performance.')

    def calibration(self):
        cal=self.read('s5_predict/test_calibration.csv');horizons=sorted(cal.horizon.unique()) if not cal.empty else []
        fig,axes=plt.subplots(2,max(1,len(horizons)),figsize=(max(7,6.5*len(horizons)),10),layout='constrained',squeeze=False)
        selected=['clinical','clinical_pca','clinical_ae','clinical_panome','clinical_elasticnet'];plotted=[]
        for j,h in enumerate(horizons):
            ax=axes[0,j];panel(ax,chr(65+j),f'{h:g}-year net-risk calibration')
            for m in selected:
                d=cal[(cal.horizon==h)&(cal.model==m)].copy();plotted.append(d)
                ax.plot(d.predicted,d.observed_net_risk,'o-',label=MODEL_LABELS[m],color=MODEL_COLORS[m],markersize=4)
            vals=cal.loc[(cal.horizon==h)&cal.model.isin(selected),['predicted','observed_net_risk']].to_numpy()
            upper=max(.05,float(np.nanmax(vals))*1.12)
            ax.plot([0,upper],[0,upper],'--',color='#777777',lw=1);ax.set(xlim=(0,upper),ylim=(0,upper),xlabel='Mean predicted net risk',ylabel='Observed KM net risk')
            ax.legend(frameon=False,fontsize=8)
            ax=axes[1,j];panel(ax,chr(65+len(horizons)+j),f'{h:g}-year IPCW Brier score')
            d=self.metrics[(self.metrics.metric=='Brier_IPCW')&(self.metrics.horizon==h)].set_index('model').reindex(MODEL_LABELS).dropna(subset=['value'])
            ax.barh([MODEL_LABELS[m] for m in d.index],d.value,color=[MODEL_COLORS[m] for m in d.index]);ax.invert_yaxis();ax.set_xlabel('Brier score (lower is better)')
        if not horizons:
            for ax in axes.ravel():unavailable(ax,'No supported calibration horizon; see metric limitations.')
        limitations=json.loads((self.root/'s5_predict/metric_limitations.json').read_text())
        self.save('Fig5.calibration',fig,{'plotted_calibration':pd.concat(plotted,ignore_index=True) if plotted else pd.DataFrame({'status':['Unavailable']}),'all_calibration':cal,'horizon_metrics':self.metrics[self.metrics.horizon>0],'metric_limitations':pd.DataFrame(limitations) if limitations else pd.DataFrame({'status':['No unsupported configured horizons']})},
                  'Top, test-set net-risk calibration by model-specific risk quintile for five prespecified comparators; these are descriptive points without interval estimates. The diagonal denotes perfect calibration. Bottom, IPCW Brier score using the stored training censoring estimator. Death-censored net risk is not actual cumulative incidence with competing mortality. Excel includes all models and time-dependent AUC/Uno C as well as Brier scores.')

    def mosaic(self):
        ap=self.root/'s6_report/person_ae_cox_attributions.csv'
        fig=plt.figure(figsize=(14,9),layout='constrained');grid=fig.add_gridspec(2,2,height_ratios=[1,1.8])
        if not ap.exists():
            unavailable(fig.add_subplot(grid[:,:]),'Individual attributions were disabled for this run.')
            self.save('Fig6.individual_molecular_profiles',fig,{'status':pd.DataFrame({'status':['Attributions unavailable']})},'Individual attribution output was disabled; no individual profiles were inferred.');return
        attrs=pd.read_csv(ap,dtype={'person_id':str});idcol=self.cfg['id_col']
        pred=pd.read_csv(self.root/'s5_predict/test_predictions.csv',dtype={idcol:str})
        ppl=self.people.copy();ppl[idcol]=ppl[idcol].astype(str)
        candidates=pred[pred[idcol].isin(attrs.person_id.unique())].sort_values('clinical_ae_log_hazard').reset_index(drop=True)
        pairs=[(abs(candidates.iloc[i].clinical_ae_log_hazard-candidates.iloc[i-1].clinical_ae_log_hazard),i-1,i) for i in range(1,len(candidates)) if candidates.iloc[i].state!=candidates.iloc[i-1].state]
        if not pairs:
            unavailable(fig.add_subplot(grid[:,:]),'No attributable pair with different candidate states is available.')
            self.save('Fig6.individual_molecular_profiles',fig,{'status':pd.DataFrame({'status':['No eligible cross-state pair']})},'No outcome-independent pair satisfying the displayed comparison was available.');return
        _,i,j=min(pairs);chosen=candidates.iloc[[i,j]];weightcols=[c for c in ppl if c.startswith('state_weight_')]
        weights=[];barrows=[];subjects=[]
        ax=fig.add_subplot(grid[0,0]);panel(ax,'A','Similar AE-Cox scores, different neighborhoods')
        for n,(_,row) in enumerate(chosen.iterrows()):
            label=f'Example {chr(65+n)}';person=ppl[ppl[idcol]==row[idcol]].iloc[0];left=0.
            subjects.append({'example':label,'candidate_state':int(row.state),'ae_log_hazard':row.clinical_ae_log_hazard})
            for s,col in zip(self.states,weightcols):
                value=float(person[col]);ax.barh(n,value,left=left,color=self.colors[s],label=f'S{s}' if n==0 else None);left+=value
                weights.append(dict(example=label,state=s,weight=value))
            q=attrs[attrs.person_id==row[idcol]].drop(columns='person_id').copy();q.insert(0,'example',label);barrows.append(q)
        ax.set_yticks([0,1],['Example A','Example B']);ax.set(xlim=(0,1),xlabel='Local neighbor state weight');ax.invert_yaxis();ax.set_ylim(1.6,-.9);ax.legend(ncol=7,frameon=False,fontsize=8,loc='upper left')
        ax=fig.add_subplot(grid[0,1]);panel(ax,'B','Frozen-model predictions');ax.set_axis_off()
        ax.text(.03,.75,'\n\n'.join(f"{s['example']}: state S{s['candidate_state']}\nAE-Cox log hazard = {s['ae_log_hazard']:.4f}" for s in subjects),transform=ax.transAxes,fontsize=12,va='top')
        bars=pd.concat(barrows,ignore_index=True)
        plotted=[];contribution_limit=bars.integrated_gradient.abs().max()*1.08
        for n,label in enumerate(['Example A','Example B']):
            q=bars[bars.example==label].copy();q=q.loc[q.integrated_gradient.abs().nlargest(10).index].sort_values('integrated_gradient')
            plotted.append(q);ax=fig.add_subplot(grid[1,n]);panel(ax,chr(67+n),f'{label}: strongest stored contributions')
            ax.barh(q.feature,q.integrated_gradient,color=np.where(q.integrated_gradient>=0,'#C05C50','#347E91'))
            ax.axvline(0,color='#777777',lw=.8);ax.set_xlabel('Integrated gradient contribution to AE-Cox log hazard');ax.set_xlim(-contribution_limit,contribution_limit)
        completeness=self.read('s6_report/attribution_completeness.csv')
        self.save('Fig6.individual_molecular_profiles',fig,{'examples':pd.DataFrame(subjects),'state_weights':pd.DataFrame(weights),'plotted_contributions':pd.concat(plotted,ignore_index=True),'stored_top20_contributions':bars,'integration_error_summary':completeness[['log_hazard_delta','attribution_sum','approximation_error']].describe().reset_index()},
                  'Illustrative pair chosen without observed disease outcomes: among test people with saved attributions, choose adjacent AE-Cox log-hazard scores in different candidate states with the smallest score gap. A, local neighbor weights, not biological pathway percentages. B, frozen AE-Cox scores. C–D, top ten absolute contributions within the stored top-20 per-person integrated gradients. Positive/negative denote model contributions, not harmful/protective causal effects. Integration baseline is zero in the preprocessed molecular space; clinical terms are fixed. Examples are anonymized and are not independent validation of disease mechanisms.')

    def stability(self):
        stability=self.read('s4_graph/graph_perturbation_stability.csv');resolution=self.read('s4_graph/resolution_search.csv')
        fig,axes=plt.subplots(1,3,figsize=(15,5),layout='constrained')
        ax=axes[0];panel(ax,'A','Graph perturbation stability')
        ax.plot(stability.replicate,stability.ARI,'o-',color='#347E91');ax.axhline(stability.ARI.mean(),color='#C05C50',ls='--',label=f'Mean ARI = {stability.ARI.mean():.3f}')
        ax.set(xlabel='Edge-dropout replicate',ylabel='Adjusted Rand index',ylim=(0,1));ax.legend(frameon=False)
        ax=axes[1];panel(ax,'B','Resolution search')
        ax.plot(resolution.resolution,resolution.n_states,'o-',color='#347E91');ax.axvline(self.graph['resolution'],color='#C05C50',ls='--',label='Selected resolution')
        ax.set(xlabel='Leiden resolution',ylabel='Number of communities');ax.legend(frameon=False)
        ax=axes[2];panel(ax,'C','Distance-based novelty by state')
        for i,split in enumerate(['train','validation','test']):
            d=self.counts[self.counts.split==split].set_index('state').reindex(self.states)
            ax.bar(np.arange(len(self.states))+(i-1)*.25,d.novel,width=.25,label=split.title())
        ax.set_xticks(range(len(self.states)),[f'S{s}' for s in self.states]);ax.set_ylabel('Fraction beyond training distance threshold');ax.legend(frameon=False,fontsize=8)
        self.save('FigS1.state_diagnostics',fig,{'edge_perturbation_ARI':stability,'resolution_search':resolution,'state_counts_and_novelty':self.counts,'graph_summary':pd.DataFrame({'metric':list(self.graph),'value':[str(v) for v in self.graph.values()]})},
                  'A, repeat Leiden after randomly deleting 10% of graph edges; AE is not retrained. ARI evaluates graph perturbation stability only. B, prespecified resolution search; eligibility also requires community size constraints (see workbook), not just community count. C, fraction above the training nearest-other-distance 99th percentile. Low stability or increased novelty limits biological interpretation of candidate states.')

    def landmark(self):
        d=self.read('s5_predict/landmark_sensitivity.csv');fig,ax=plt.subplots(figsize=(10,6),layout='constrained');panel(ax,'A','Fixed-score landmark sensitivity')
        if d.empty:unavailable(ax,'Insufficient eligible follow-up/events for landmark analysis.')
        else:
            c=self.metrics[self.metrics.metric=='Harrell_C'][['model','value']]
            plot=d.merge(c,on='model');yy=np.arange(len(plot))
            ax.scatter(plot.value,yy-.1,label='All test participants',color='#347E91')
            ax.scatter(plot.Harrell_C,yy+.1,label=f'{plot.lag_years.iloc[0]:g}-year landmark subset',color='#C05C50')
            ax.set_yticks(yy,[MODEL_LABELS.get(m,m) for m in plot.model]);ax.invert_yaxis();ax.set_xlabel('Harrell C');ax.legend(frameon=False)
        self.save('FigS2.landmark_sensitivity',fig,{'landmark_results':d,'all_test_C':self.metrics[self.metrics.metric=='Harrell_C']},
                  'Fixed baseline predictions evaluated among people still event-free and under follow-up after the configured landmark; follow-up is reset at the landmark. No models are retrained. The full test and landmark populations differ, so changes are descriptive, not a paired improvement estimate. No interval estimates were calculated for this sensitivity.')


def export_final(root):
    root=Path(root).resolve();out=root/'publication'
    if not (root/'s6_report/DONE.json').exists():raise ValueError('Publication export requires completed s6_report')
    source_paths=[root/'manifest.json']
    for stage in ['s2_preprocess','s3_representation','s4_graph','s5_predict','s6_report']:
        source_paths.extend(sorted(p for p in (root/stage).iterdir() if p.suffix in ('.csv','.json')))
    extra=root/'validation/center_stratified_state_cox.csv'
    if extra.exists():source_paths.append(extra)
    sources={str(p.relative_to(root)):sha(p) for p in source_paths}
    import openpyxl
    payload={'exporter_sha256':sha(Path(__file__)),'sources':sources,'dpi':DPI,
             'versions':{'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__,'openpyxl':openpyxl.__version__}}
    signature=hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()
    config=json.loads((root/'manifest.json').read_text())['config']
    destination=root.parent if config['run_name']=='main' else out
    out.mkdir(exist_ok=True)
    lock=out/'.lock'
    try:fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY)
    except FileExistsError:raise RuntimeError(f'Publication export already locked: {lock}')
    os.write(fd,str(os.getpid()).encode());os.close(fd)
    try:
        pm=out/'publication_manifest.json';previous=json.loads(pm.read_text()) if pm.exists() else {}
        cached=previous.get('signature')==signature and all((out/name).is_file() and sha(out/name)==h for name,h in previous.get('outputs',{}).items()) and bool(previous.get('outputs'))
        if not cached:
            with tempfile.TemporaryDirectory(prefix='.figures-',dir=root) as td:
                temp=Path(td)
                with plt.rc_context({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,'axes.labelsize':11,'savefig.facecolor':'white','axes.linewidth':.8}):
                    pub=Publication(root,temp)
                    for method in ['cohort','landscape','survival','prediction','calibration','mosaic','stability','landmark']:getattr(pub,method)()
                pd.DataFrame(pub.items).to_csv(temp/'Fig_index.csv',index=False)
                notes=['# Panome manuscript figures',f'\nSource run: `{root}`','\nEach PNG is paired with an Excel workbook containing plotted values, supporting tables, and a Readme sheet. Candidate states remain exploratory.']
                for item in pub.items:notes.extend([f'\n## {item["figure"]}',item['caption']])
                (temp/'Fig_legends.md').write_text('\n\n'.join(notes))
                outputs={p.name:sha(p) for p in temp.iterdir() if p.is_file()}
                new={'signature':signature,**payload,'run':str(root),'outputs':outputs}
                for p in temp.iterdir():os.replace(p,out/p.name)
                tmp=out/'publication_manifest.json.tmp';tmp.write_text(json.dumps(new,indent=2));os.replace(tmp,pm)
        else:print('Final figures: verified source data and exporter; cached.',flush=True)
        record=json.loads(pm.read_text())
        if destination!=out:
            for name,h in record['outputs'].items():
                target=destination/name
                if target.exists() and sha(target)==h:continue
                temp=target.with_name('.'+target.name+'.tmp');shutil.copyfile(out/name,temp);os.replace(temp,target)
            target=destination/'Fig_manifest.json';tmp=target.with_suffix('.json.tmp');tmp.write_text(json.dumps(record,indent=2));os.replace(tmp,target)
        print(f'Final output: {destination}',flush=True)
        return destination
    finally:lock.unlink(missing_ok=True)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run-dir',required=True,type=Path)
    export_final(p.parse_args().run_dir)
