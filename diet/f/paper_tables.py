"""Assemble numbered supplement tables from analysis outputs and fixed design definitions."""
import argparse
import csv
import json
import os
import math
import shutil
import tempfile
from pathlib import Path

from openpyxl import Workbook, load_workbook
from openpyxl.styles import Font, PatternFill
from openpyxl.utils import get_column_letter

CONFIG = Path(__file__).resolve().parents[1] / "config"
UK = "dose/ukb/"
SC = UK + "04score/"
ML = UK + "03ml/"
ME = UK + "05mechanism/"
PAT = SC + "pattern/"
SENS = SC + "IS_sensitivity/"
M = "final/maha/"
UWB = "maha/ukb/Fig2.prospective.out.xlsx"
NW = M + "Fig3.nhanes.mortality.out.xlsx"
CW = M + "Fig4.chns.mortality.out.xlsx"
UW = M + "Fig2.ukb.precision.out.xlsx"

# Each block has a source table, optional workbook sheet and optional row selector.
def src(path, sheet=None, subset=None):
    return (path, sheet, subset)

def definition(study,n):
    return src("@"+study+"_s"+str(n)+".tsv")

def blocks(*paths):
    return [src(p) for p in paths]

def mortality(row):
    return any(word in str(row.get("Outcome","")).lower() for word in ("mort","death"))

DOSE = [
    ("Food items and groups", [definition("dose",1),src(UK+"02foodwas/food_group_map.csv")]),
    ("ICD-10 definitions", [definition("dose",2)]),
    ("UKB covariate missingness", blocks(UK+"00data/covariate_missingness.csv")),
    ("UKB baseline by DOSE quartile", blocks("dose/outputs/ukb_table1/Table1_data.csv")),
    ("Food groups and incident ischemic stroke", blocks(UK+"02foodwas/group/Group_asso_cvd_stroke_i.csv")),
    ("ML feature prioritization", blocks(*[ML+f+".csv" for f in ("Food_Importance","SHAP_CrossModel_Ranks","RFE_Ranks","Selection_Methods_Summary")])),
    ("Cross-validated prediction by feature set", blocks(ML+"SelectionValidation_Summary.csv")),
    ("Dietary patterns and stroke", blocks(PAT+"Table_Pattern_Stroke_M1_M3.csv")),
    ("DOSE quartiles and stroke", blocks(PAT+"Table_DOSE_Quartiles_M1_M3.csv")),
    ("Ischemic stroke sensitivity analyses", blocks(SENS+"Table_S13_Sensitivity.csv")),
    ("Ischemic stroke subgroups", blocks(SENS+"Table_subgroups.csv")),
    ("Joint and stratified genetic risk analyses", blocks(SENS+"Table_S15_PRS_DOSE_joint.csv")),
    ("NHANES two-day 2003-2018 participants", blocks("dose/nhanes/NHANES_baseline.csv")),
    ("Other health outcomes", blocks(SC+"otherOutcomes/Table_DOSE_OtherOutcomes_M1_excluding_stroke.csv")),
    ("Metabolite signature and ischemic stroke", blocks(*[UK+"06signature/metabolite/association_mediation/metabolite_signature_IS_"+f+".csv" for f in ("continuous","quartiles")])),
    ("Metabolite associations", blocks(ME+"metabolite/M2/01_all_metabolite_associations.csv")),
    ("Metabolite mediation", blocks(ME+"metabolite/M2/Table_S17.csv")),
    ("Inflammatory marker associations", blocks(ME+"inflammation/M2/01_all_inflammation_associations.csv")),
    ("Inflammatory marker mediation", blocks(ME+"inflammation/M2/Table_S19.csv")),
    ("Plasma protein associations", blocks(ME+"protein/M2/01_all_protein_associations.csv")),
    ("Brain MRI associations", blocks(ME+"04mri/02_DOSE_MRI.csv")),
]
MAHA = [
    ("Cross-cohort design and primary sample sizes", [definition("maha",1),src(UWB,"Fig2_main_assoc","counts_mortality"),src(NW,"associations","counts"),src(CW,"primary_results","counts"),src(CW,"high_event_results","counts")]),
    ("Score construction", [definition("maha",2)]),
    ("UKB outcome definitions", [definition("maha",3)]),
    ("UKB baseline by adherence", [src(M+"Table1.ukb.out.xlsx")]),
    ("CHNS baseline by adherence", [src(M+"Table1.chns.out.xlsx","Table1")]),
    ("UKB six focus outcomes", [src(UWB,"Fig2_main_assoc")]),
    ("Cross-cohort mortality associations", [src(UWB,"Fig2_main_assoc","mortality"),src(NW,"associations"),src(CW,"primary_results"),src(CW,"high_event_results")]),
    ("Mutually adjusted mortality contrasts", [src(UW,"head2head","mortality"),src(NW,"head2head"),src(CW,"primary_head2head"),src(CW,"high_event_head2head")]),
    ("Equivalence, Bayes factors and precision", [src(UW,"MAHA_precision"),src(NW,"equivalence"),src(NW,"bayes"),src(NW,"power")]),
    ("CHNS high-event sensitivity", [src(CW,"high_event_results")]+[src(M+"FigS7.chns.high_event_sensitivity.out.xlsx",s) for s in ("model5_end_wave","model4_last_contact","model5_last_contact","maha_algorithms")]),
    ("Guideline definitions and implementation", [definition("maha",11)])]

def cell_value(x):
    if x in ("","NA","NaN"): return None
    if isinstance(x,str):
        try:
            return int(x) if x.lstrip("-").isdigit() and len(x)<16 else float(x)
        except ValueError:
            return x
    return x

def read_table(path, sheet, subset):
    if path.suffix == ".xlsx":
        wb = load_workbook(path, read_only=True, data_only=True)
        try:
            ws = wb[sheet] if sheet else wb.worksheets[0]
            rows = list(ws.values)
        finally:
            wb.close()
    else:
        with path.open(encoding="utf-8-sig",newline="") as inp:
            rows = [list(map(cell_value,row)) for row in csv.reader(inp, delimiter="\t" if path.suffix==".tsv" else ",")]
    if len(rows)<2: raise ValueError("No data rows")
    head,*body=rows
    if subset:
        if "mortality" in subset:
            body=[row for row in body if mortality(dict(zip(head,row)))]
        if "counts" in subset:
            # Keep per-score rows where sample sizes differ across model fits.
            keep=[i for i,h in enumerate(head) if str(h).lower() in {"outcome","diet","score","model","n","n_total","n_event","n_events","events","deaths","cases","n_deaths","n_complete","n_unweighted","event"}]
            names={str(head[i]).lower() for i in keep}
            if not names.intersection({"n","n_total","n_complete","n_unweighted"}) or not names.intersection({"deaths","events","cases","n_event","n_events","n_deaths","event"}):
                raise ValueError("Missing sample-size or event columns: "+str(head))
            body=list(dict.fromkeys(tuple(row[i] for i in keep) for row in body))
            head=[head[i] for i in keep]
    if not body: raise ValueError("No matching rows")
    return [head,*body]

def checks(root,study):
    records=[]
    def measure(label, expected, path, column, sheet=None, select=None, digits=0):
        actual=None
        try:
            head,*body=read_table(root/path,sheet,None)
            d=[dict(zip(head,row)) for row in body]
            cells={float(row[column]) for row in d if select is None or select(row)}
            if len(cells)!=1: raise ValueError("Expected a single distinct value")
            actual=cells.pop()
            if not math.isfinite(actual): raise ValueError("Non-finite value")
            status="matches" if round(actual,digits)==expected else "differs"
        except (OSError,KeyError,TypeError,ValueError):
            status="unavailable"
        records.append([label,expected,actual,status,path,sheet])
    if study=="dose":
        qc="dose/outputs/ukb_table1/Table1_cohort_QC.csv"
        total=lambda r:r["group"]=="Total"
        measure("UKB primary participants",191494,qc,"N",select=total)
        measure("UKB incident ischemic strokes",2585,qc,"Events",select=total)
        for column,value,contrast in [("HR",.71,"P90 versus P10"),("HR_Q4vsQ1",.69,"Q4 versus Q1")]:
            measure("UKB DOSE M2 HR, "+contrast,value,PAT+"DOSE_IS_M2_summary.csv",column,select=lambda r:r["diet"]=="dose",digits=2)
        nh="dose/nhanes/NHANES_baseline.csv"
        measure("NHANES primary participants",34861,nh,"n",select=lambda r:r["variable"]=="N")
        measure("NHANES prevalent stroke cases",1380,nh,"n",select=lambda r:r["variable"]=="stroke_prevalent" and str(r["level"])=="1")
    else:
        for label,n,events,path,sheet,ncol,ecol in [
            ("UKB",128562,9853,UWB,"Fig2_main_assoc","N_total","N_event"),
            ("NHANES",27478,2900,NW,"associations","N_total","N_event"),
            ("CHNS fixed wave",7387,219,CW,"primary_results","N","Events"),
            ("CHNS high-event",8981,781,CW,"high_event_results","N","Events")]:
            def main(r):
                return str(r.get("Diet","")).upper()=="MAHA" and (label!="UKB" or any(k in str(r["Outcome"]).lower() for k in ("death","mort")))
            measure(label+" participants",n,path,ncol,sheet,main)
            measure(label+" mortality events",events,path,ecol,sheet,main)
    return records


def style_header(ws, row):
    for c in ws[row]:
        c.fill=PatternFill("solid",fgColor="235267")
        c.font=Font(color="FFFFFF",bold=True)

def write_book(root,study,scope):
    directory=root/"final"/study
    directory.mkdir(parents=True,exist_ok=True)
    onto=directory/(study.upper()+".Supplementary.Tables"+("" if scope=="all" else "."+scope)+".xlsx")
    wb=Workbook()
    index=wb.active
    index.title="Contents"
    index.append(["Table","Title","Analysis source","Source sheet","Rows","Status"])
    issues=[]
    for i,(title,sources) in enumerate(DOSE if study=="dose" else MAHA,1):
        if scope!="all":
            sources=[b for b in sources if b[0].startswith("@") or scope in b[0].lower()]
            if not sources: continue
        ws=wb.create_sheet(f"Table S{i}")
        ws.append([title]); ws["A1"].font=Font(size=14,bold=True,color="20485A")
        ws.freeze_panes="B4"
        for loc,sheet,subset in sources:
            path=CONFIG/"manuscripts"/loc[1:] if loc.startswith("@") else root/loc
            status="OK"
            try:
                data=read_table(path,sheet,subset)
            except (OSError,KeyError,ValueError) as e:
                status=str(e)
                if scope=="all": issues.append(f"{ws.title}: {loc} {status}")
                index.append([ws.title,title,loc,sheet,0,status])
                ws.append(["Unavailable",loc,status]); continue
            ws.append([loc,sheet]); header_row=ws.max_row+1
            for row in data: ws.append(row)
            style_header(ws,header_row)
            index.append([ws.title,title,loc,sheet,len(data)-1,status])
        for col in ws.iter_cols():
            lengths=[len(str(c.value or "")) for c in col[:120] if c.row!=1]
            ws.column_dimensions[col[0].column_letter].width=min(65,max(14,max(lengths,default=10)+2))
    qa=wb.create_sheet("Manuscript checks")
    qa.append(["Metric","Published","Recomputed","Comparison","Source","Source sheet"])
    for record in checks(root,study):
        if scope=="all" or record[0].lower().startswith(scope): qa.append(record)
    style_header(qa,1)
    for i,width in enumerate((39,16,17,15,80,28),1): qa.column_dimensions[get_column_letter(i)].width=width
    style_header(index,1)
    for k,v in enumerate((14,46,76,28,12,45),1): index.column_dimensions[get_column_letter(k)].width=v
    index.freeze_panes="C2"
    index.auto_filter.ref=index.dimensions
    with tempfile.NamedTemporaryFile(dir=directory,suffix=".xlsx",delete=False) as f: temp=Path(f.name)
    try:
        wb.save(temp); os.replace(temp,onto)
    finally: temp.unlink(missing_ok=True)
    print(onto,flush=True)
    return issues

def figures(root,study,scope):
    if scope!="all": return []
    dest=root/"final"/study
    missing=[]
    for spec in json.loads((CONFIG/"papers.json").read_text())[study]:
        if not (dest/spec["file"]).is_file(): missing.append(spec["file"])
    if study=="maha":
        status=dest/"publication.status.tsv"
        if status.is_file():
            with status.open() as f:
                for row in csv.DictReader(f,delimiter="\t"):
                    if row["status"] != "complete": missing.append(f"{row['item']}: {row['status']}")
        else: missing.append(str(status))
        for i in range(1,11):
            if not list(dest.glob(f"FigS{i}.*.png")): missing.append(f"FigS{i}")
    else:
        status=dest/"composite_status.tsv"
        if status.is_file():
            with status.open() as f:
                for row in csv.DictReader(f,delimiter="\t"):
                    if row["status"]!="complete": missing.append(f"{row['figure']}: {row['status']} {row['missing_panels']}")
        else: missing.append(str(status))
        for stem,source in [("FigS1.NHANES.flow","dose/nhanes/Sfig1_nhanes_primary"),("FigS2.Diet.correlations",PAT+"figures/Fig_Pattern_Correlations_10patterns")]:
            found=False
            for suffix in (".png",".pdf"):
                path=root/(source+suffix)
                if path.is_file():
                    shutil.copyfile(path,dest/(stem+suffix)); found=True
            if not found: missing.append(source)
    return missing

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,required=True); p.add_argument("--study",choices=["maha","dose","all"],required=True)
    p.add_argument("--scope",default="all")
    a=p.parse_args()
    issues=[]
    for study in (["maha","dose"] if a.study=="all" else [a.study]):
        issues.extend(write_book(a.output,study,a.scope))
        issues.extend(figures(a.output,study,a.scope))
    if issues: raise SystemExit("Publication inputs pending:\n"+"\n".join(issues))
