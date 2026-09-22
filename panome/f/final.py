"""Publication figures with corresponding source-data workbooks. No model fitting."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import STAGES, log, dump

COLORS = ["#007F86","#E0874A","#6875B8","#9B6FA7","#719949","#C56E79"]

def read(root,relative):
    path = root/relative
    if not path.is_file():
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()

def label(ax,letter,title):
    ax.set_title(f"{letter}  {title}",loc="left",fontweight="bold",fontsize=10,pad=10)
    ax.spines[["top","right"]].set_visible(False)

def unavailable(ax,message):
    ax.text(.5,.5,message,ha="center",va="center",transform=ax.transAxes,wrap=True,color="#777777")
    ax.set_xticks([]);ax.set_yticks([])

def export_final(root):
    root = Path(root)
    manifest = json.loads((root/"manifest.json").read_text())
    if manifest.get("schema") != 3:
        raise ValueError("final requires a v3 run; legacy runs such as main remain readable with the backed-up legacy code. Use --run-name v3 for new results.")
    cfg = manifest["config"]
    # Summarizing requires complete stage artifacts, but not raw inputs/dependencies for training.
    from panome import completed
    for stage in STAGES:
        if not completed(root/stage,manifest["signature"]):
            raise ValueError(f"final requires completed {stage}")
    out = root/"publication"
    out.mkdir(exist_ok=True)
    status = "SYNTHETIC SOFTWARE TEST" if cfg["demo"] else ("REAL-DATA PILOT" if cfg["max_samples"] else "HELD-OUT RESULTS")
    title = f"Panome | {cfg['trait']} | {cfg['biom']} | {status}"
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,
                         "axes.spines.top":False,"axes.spines.right":False})

    person_aliases = {}

    def source_table(table):
        table = table.copy()
        for col in [cfg["id_col"], "person_id", "reference_person_id", "person_1", "person_2"]:
            if col in table:
                def alias(value):
                    if pd.isna(value):
                        return value
                    key = str(value)
                    if key not in person_aliases:
                        person_aliases[key] = f"Example_{len(person_aliases)+1:06d}"
                    return person_aliases[key]
                table[col] = table[col].map(alias)
        return table

    def save(fig,name,tables,caption):
        fig.suptitle(title,fontsize=13,fontweight="bold")
        fig.text(.02,.012,caption,fontsize=8,color="#555555")
        fig.tight_layout(rect=(.01,.045,.995,.95),h_pad=2.5,w_pad=2.5)
        fig.savefig(out/f"{name}.png",dpi=300)
        fig.savefig(out/f"{name}.pdf")
        plt.close(fig)
        with pd.ExcelWriter(out/f"{name}.xlsx",engine="openpyxl") as writer:
            pd.DataFrame({"note":[title,caption,"No UKB result is implied by a synthetic run.",
                "Signature: "+manifest["signature"]]}).to_excel(writer,sheet_name="Readme",index=False)
            for sheet,table in tables.items():
                if len(table)>1_000_000:
                    table = table.iloc[:1_000_000]
                source_table(table).to_excel(writer,sheet_name=sheet[:31],index=False)
            for sheet in writer.book.worksheets:
                sheet.freeze_panes="A2"
                sheet.auto_filter.ref=sheet.dimensions
                for col in sheet.columns:
                    width=min(55,max(12,max(len(str(cell.value or "")) for cell in list(col)[:100])+2))
                    sheet.column_dimensions[col[0].column_letter].width=width
        log("DONE",name)

    persons = read(root,"s6_report/person_molecular_states.csv")
    profiles = read(root,"s6_report/training_state_profiles.csv")
    counts = read(root,"s6_report/state_counts.csv")
    stability = read(root,"s4_graph/graph_perturbation_stability.csv")
    history = read(root,"s3_representation/ae_history.csv")
    fig,axs = plt.subplots(2,2,figsize=(12,9))
    train = persons[persons.split=="train"]
    points = train.sample(min(len(train),6000),random_state=cfg["seed"])
    axs[0,0].scatter(points.AE1,points.AE2,c=points.state,s=7,alpha=.55,cmap="tab20",rasterized=True)
    axs[0,0].set(xlabel="AE coordinate 1",ylabel="AE coordinate 2")
    label(axs[0,0],"A","Discovery molecular landscape")
    for split,color in zip(["train","validation","test"],COLORS):
        part=counts[counts.split==split]
        axs[0,1].plot(part.state,part.n/part.n.sum(),marker="o",label=split,color=color)
    axs[0,1].set(xlabel="Projected molecular state",ylabel="Fraction of people")
    axs[0,1].legend(frameon=False)
    label(axs[0,1],"B","State representation across splits")
    axs[1,0].plot(history.epoch,history.train_mse,label="Train",color=COLORS[0])
    axs[1,0].plot(history.epoch,history.validation_mse,label="Validation",color=COLORS[1])
    axs[1,0].set(xlabel="Epoch",ylabel="Observed-entry MSE");axs[1,0].legend(frameon=False)
    label(axs[1,0],"C","Unsupervised reconstruction")
    if not stability.empty:
        axs[1,1].plot(stability.replicate,stability.ARI,"o-",color=COLORS[2])
        axs[1,1].set(xlabel="Edge-dropout replicate",ylabel="Adjusted Rand index",ylim=(-.05,1.05))
    else:
        unavailable(axs[1,1],"Stability estimates unavailable; see graph summary")
    label(axs[1,1],"D","Graph perturbation stability")
    save(fig,"Fig1_Atlas",{"landscape":points,"state_counts":counts,"reconstruction":history,
         "stability":stability,"state_profiles":profiles},
         "State labels are descriptive. Edge dropout does not measure complete pipeline reproducibility.")

    metrics = read(root,"s5_predict/test_metrics.csv")
    calibration = read(root,"s5_predict/test_calibration.csv")
    contrasts = read(root,"s5_predict/paired_contrasts.csv")
    landmark = read(root,"s5_predict/landmark_sensitivity.csv")
    model_status = read(root,"s5_predict/model_status.csv")
    metric = "Harrell_C" if cfg["outcome_type"]=="survival" else "R2"
    values = metrics[metrics.metric==metric].copy()
    fig,axs = plt.subplots(2,2,figsize=(14,11))
    axs[0,0].barh(values.model,values.value,color=[COLORS[1] if m.startswith("panome") else COLORS[0] for m in values.model])
    axs[0,0].set_xlabel(metric)
    label(axs[0,0],"A","All held-out model results")
    wanted = ["clinical","clinical_elasticnet","transformer_column","panome_transformer"]
    horizon = cfg["primary_horizon"] if cfg["outcome_type"]=="survival" else 0
    selected = calibration[(calibration.horizon==horizon)&calibration.model.isin(wanted)]
    for model,color in zip(wanted,COLORS):
        group = selected[selected.model==model]
        if not group.empty:
            axs[0,1].plot(group.predicted,group.observed,"o-",label=model,color=color)
    if selected.empty or selected.observed.notna().sum()==0:
        unavailable(axs[0,1],"Prespecified horizon is not estimable")
    else:
        bounds = np.r_[selected.predicted,selected.observed.dropna()]
        lo,hi = float(np.min(bounds)),float(np.max(bounds))
        axs[0,1].plot([lo,hi],[lo,hi],"--",color="#999999")
        axs[0,1].legend(frameon=False,fontsize=7)
        axs[0,1].set(xlabel="Predicted",ylabel="Observed")
    label(axs[0,1],"B","Calibration at prespecified horizon")
    sel = contrasts[(contrasts.model=="panome_transformer")|
                    ((contrasts.model=="panome_varying")&(contrasts.reference=="clinical_ae"))]
    if not sel.empty:
        ticks = np.arange(len(sel))
        axs[1,0].scatter(sel.delta,ticks,color=COLORS[1])
        for row,pos in zip(sel.itertuples(),ticks):
            if np.isfinite(row.lower) and np.isfinite(row.upper):
                axs[1,0].plot([row.lower,row.upper],[pos,pos],color=COLORS[1])
        axs[1,0].set_yticks(ticks,[f"{m} vs {b}" for m,b in zip(sel.model,sel.reference)],fontsize=8)
        axs[1,0].axvline(0,color="#999999",ls="--")
        axs[1,0].set_xlabel(f"Delta {metric} (model minus reference)")
    else:
        unavailable(axs[1,0],"Conditional model comparison unavailable")
    label(axs[1,0],"C","Prespecified individual-model increments")
    if cfg["outcome_type"]=="survival" and not landmark.empty:
        for model,color in zip(wanted,COLORS):
            group = landmark[landmark.model==model]
            if not group.empty:
                axs[1,1].plot(group.landmark,group.Harrell_C,"o-",label=model,color=color)
        axs[1,1].set(xlabel="Landmark year",ylabel="Conditional Harrell C")
        axs[1,1].legend(frameon=False,fontsize=7)
    else:
        errors=metrics[metrics.metric=="RMSE"]
        if not errors.empty:
            axs[1,1].barh(errors.model,errors.value,color=COLORS[2])
            axs[1,1].set_xlabel("RMSE in outcome units")
        else:
            unavailable(axs[1,1],"Landmark sensitivity unavailable")
    label(axs[1,1],"D","Landmark sensitivity" if cfg["outcome_type"]=="survival" else "Prediction error")
    save(fig,"Fig2_Prediction",{"metrics":metrics,"calibration":calibration,"contrasts":contrasts,
        "landmark":landmark,"model_status":model_status},
        "Intervals condition on fitted models. Survival probabilities are net risks with death censored.")

    pairs=read(root,"s6_report/matched_risk_pairs.csv")
    pair_features=read(root,"s6_report/matched_pair_features.csv")
    attr=read(root,"s6_report/person_ae_attributions.csv")
    completeness=read(root,"s6_report/attribution_completeness.csv")
    fig,axs=plt.subplots(2,2,figsize=(13,9))
    matrix=profiles.pivot(index="state",columns="feature",values="mean_standardized")
    top=matrix.abs().max().nlargest(min(16,matrix.shape[1])).index
    hm=axs[0,0].imshow(matrix[top],aspect="auto",cmap="RdBu_r",vmin=-2,vmax=2)
    axs[0,0].set_xticks(range(len(top)),top,rotation=75,ha="right",fontsize=7)
    axs[0,0].set_yticks(range(len(matrix)),[f"S{j}" for j in matrix.index])
    fig.colorbar(hm,ax=axs[0,0],shrink=.8,label="Mean standardized abundance")
    label(axs[0,0],"A","Training state profiles")
    if not pairs.empty:
        first=pair_features[pair_features.pair==pairs.iloc[0].pair].head(10)
        axs[0,1].plot(first.feature,first.person_1_z,"o-",color=COLORS[0],label="Person 1")
        axs[0,1].plot(first.feature,first.person_2_z,"o-",color=COLORS[1],label="Person 2")
        axs[0,1].tick_params(axis="x",rotation=70,labelsize=7)
        axs[0,1].legend(frameon=False);axs[0,1].set_ylabel("Processed molecular value")
    else:
        unavailable(axs[0,1],"No pairs meet the fixed risk/state criterion")
    label(axs[0,1],"B","A prespecified matched-risk example")
    if not attr.empty:
        first=attr[attr.person_id==attr.person_id.iloc[0]].head(10).iloc[::-1]
        axs[1,0].barh(first.feature,first.contribution,color=np.where(first.contribution>=0,COLORS[1],COLORS[0]))
        axs[1,0].set_xlabel("Contribution to clinical_ae prediction")
    else:
        unavailable(axs[1,0],"AE attribution not requested/available")
    label(axs[1,0],"C","One person's AE-model attribution")
    axs[1,1].hist(persons.loc[persons.split=="test","effective_neighbors"],bins=20,color=COLORS[2],alpha=.8)
    axs[1,1].set(xlabel="Effective reference neighbors",ylabel="Test people")
    label(axs[1,1],"D","Individual neighborhood support")
    save(fig,"Fig3_Individuals",{"state_profiles":profiles,"pairs":pairs,"pair_features":pair_features,
        "AE_attributions":attr,"integration_error":completeness,
        "support":persons[[cfg["id_col"],"split","state","novel","effective_neighbors","observed_fraction"]]},
        "Profiles and pairs are descriptive. AE attribution and reference attention do not identify causal mechanisms.")
    dump(out/"figure_manifest.json",dict(signature=manifest["signature"],
        files=sorted(p.name for p in out.glob("Fig*")),run_type=status))
