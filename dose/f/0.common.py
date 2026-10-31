"""Shared, working-directory-independent paths for DOSE analysis scripts."""

# Shared PNG byte limit; PDF/SVG export settings are preserved.
import sys as _png_sys
from pathlib import Path as _PngPath
_png_sys.path.append(str(_PngPath(__file__).resolve().parents[2] / "0f"))
from png_output import save_png_figure

import os
from pathlib import Path


def drive():
    return Path("D:/") if os.name == "nt" else Path("/mnt/d")


def output_dir():
    return Path(os.getenv("DOSE_OUTDIR", os.getenv("DOSE_PROJECT", str(drive() / "analysis/dose"))))


def mri_atlas_dir():
    if os.getenv("DOSE_MRI_ATLAS_DIR"):
        return Path(os.environ["DOSE_MRI_ATLAS_DIR"])
    local = output_dir() / ".tools/mri_atlases"
    if local.is_dir():
        return local
    return Path(os.getenv("DOSE_RUNTIME_DIR", str(Path.home() / ".local/share/dose"))) / "mri_atlases"


def ukb_dir():
    return Path(os.getenv("DOSE_UKB", str(drive() / "data/ukb")))


def files_dir():
    return Path(os.getenv("DOSE_FILES_DIR", str(drive() / "files/diet")))


def definition_file(name):
    if name not in {"diet", "vip", "icd10_other"}:
        raise ValueError(f"Unknown DOSE dictionary: {name}")
    return ukb_dir() / "common" / f"DOSE.{name}.lst"


def vip_fields(variables, dictionary=None):
    """Resolve baseline raw field names from the study-specific VIP list."""
    import csv
    import re
    path = Path(dictionary) if dictionary is not None else definition_file("vip")
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = [row for row in csv.reader(stream, delimiter="\t") if row and row[0].strip()]
    result = []
    for variable in variables:
        matches = [row for row in rows if len(row) >= 3 and row[1].strip() == variable]
        if len(matches) != 1:
            raise ValueError(f"{path}: expected one definition for {variable}; found {len(matches)}")
        field, _, instance = [value.strip() for value in matches[0][:3]]
        if not re.fullmatch(r"p[0-9]+", field) or instance != "i0":
            raise ValueError(f"{path}: {variable} must specify a raw field and baseline instance i0")
        result.append(f"{field}_{instance}")
    return result


def mri_dictionary():
    return Path(os.getenv("DOSE_MRI_DICTIONARY", str(files_dir() / "MRI_analysis_dictionary.csv")))


def jobs():
    value = int(os.getenv("DOSE_THREADS", "4"))
    if value < 1:
        raise ValueError("DOSE_THREADS must be positive")
    return value


def replace_results():
    return os.getenv("DOSE_REPLACE", "FALSE").upper() == "TRUE"


def lasso_preprocessor(log1p=False):
    import numpy as np
    from sklearn.impute import SimpleImputer
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import FunctionTransformer, StandardScaler
    steps = [("imputer", SimpleImputer(strategy="median"))]
    if log1p:
        steps.append(("log1p", FunctionTransformer(np.log1p, validate=False)))
    return Pipeline(steps + [("scaler", StandardScaler())])


def fit_lasso_full(x, y, alphas, folds, max_iter, seed, log1p=False):
    """Return only fitted parameters/small results, never the transformed matrix."""
    import numpy as np
    from sklearn.linear_model import Lasso, LassoCV, lasso_path
    pre = lasso_preprocessor(log1p)
    xs = pre.fit_transform(x)
    cv = LassoCV(alphas=alphas, cv=folds, max_iter=max_iter, random_state=seed, n_jobs=1).fit(xs, y)
    model = Lasso(alpha=float(cv.alpha_), max_iter=max_iter).fit(xs, y)
    selected = np.flatnonzero(np.abs(model.coef_) > 1e-12)
    if not selected.size:
        raise RuntimeError("LASSO selected no feature; inspect the alpha grid")
    full = model.predict(xs)
    path_alphas, path_coefs, _ = lasso_path(xs, y, alphas=alphas)
    return dict(pre=pre, cv=cv, model=model, selected=selected, full=full,
                path_alphas=path_alphas, path_coefs=path_coefs, full_r2=float(model.score(xs, y)))


def predict_lasso_fold(x, y, train_idx, test_idx, alphas, folds, max_iter, seed, log1p=False):
    """One outer fold at a time; joblib persists its prediction immediately."""
    from sklearn.linear_model import LassoCV
    pipeline = lasso_preprocessor(log1p)
    pipeline.steps.append(("lasso", LassoCV(alphas=alphas, cv=folds, max_iter=max_iter,
                                          random_state=seed, n_jobs=1)))
    pipeline.fit(x[train_idx], y[train_idx])
    return pipeline.predict(x[test_idx]), float(pipeline.named_steps["lasso"].alpha_)


def lasso_checkpoint_cache(out):
    import hashlib
    import importlib.metadata
    from joblib import Memory
    fingerprint = hashlib.sha256(Path(__file__).read_bytes() + Path(__file__).with_name("signature.py").read_bytes())
    for package in ("numpy", "pandas", "scikit-learn", "joblib"):
        fingerprint.update((package + importlib.metadata.version(package)).encode())
    cache = Memory(Path(out) / ".checkpoints" / fingerprint.hexdigest()[:16], verbose=0)
    if replace_results():
        cache.clear(warn=False)
    return cache


def lasso_reference_score(x, selected, pre, model, log1p=False):
    import numpy as np
    # Avoid imputing every omics column for a score that uses only selected ones.
    # Match the column-major selected matrix in the original full imputation.
    reference = np.array(x[:, selected], order="F", copy=True)
    reference = np.where(np.isnan(reference), pre.named_steps["imputer"].statistics_[selected], reference)
    if log1p:
        reference = np.log1p(reference)
    reference = reference - reference.min(axis=0) + 0.01
    return np.dot(reference, model.coef_[selected])


def ml_checkpoint_cache(out, step):
    import hashlib
    import importlib.metadata
    from joblib import Memory
    fingerprint = hashlib.sha256(Path(__file__).read_bytes() + Path(__file__).with_name("ml.py").read_bytes())
    for package in ("numpy", "pandas", "scikit-learn", "joblib", "lightgbm", "xgboost", "catboost"):
        try:
            fingerprint.update((package + importlib.metadata.version(package)).encode())
        except importlib.metadata.PackageNotFoundError:
            pass
    cache = Memory(Path(out) / ".checkpoints" / step / fingerprint.hexdigest()[:16], verbose=0)
    if replace_results():
        cache.clear(warn=False)
    return cache


def fit_predict_ml_fold(model, x_train, y_train, x_test):
    from sklearn.base import clone
    fitted = clone(model).fit(x_train, y_train)
    return ml_tools().predict_score(fitted, x_test)


def fit_importance_ml_fold(model, x_train, y_train, x_test):
    import numpy as np
    from sklearn.base import clone
    fitted = clone(model).fit(x_train, y_train)
    gain = fitted.booster_.feature_importance(importance_type="gain")
    split = fitted.booster_.feature_importance(importance_type="split")
    contrib = fitted.booster_.predict(x_test, pred_contrib=True)
    return gain, split, np.abs(contrib[:, :-1]).mean(axis=0)


import sys
from functools import lru_cache
from types import SimpleNamespace



def extract_columns(source, dest, fields, rename=None):
    import gzip, json, os
    source, dest = Path(source), Path(dest)
    meta_path = dest.with_suffix('.source.json')
    meta = dict(source=str(source.resolve()), size=source.stat().st_size,
                mtime_ns=source.stat().st_mtime_ns, fields=fields, rename=rename)
    if not replace_results() and dest.exists() and meta_path.exists() and json.loads(meta_path.read_text()) == meta:
        print('Up-to-date:', dest, flush=True)
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    temp = dest.with_suffix('.tmp')
    opener = gzip.open if source.suffix == '.gz' else open
    with opener(source, 'rb') as src, temp.open('wb') as out:
        header = src.readline().rstrip(b'\r\n').decode().split('\t')
        missing = set(fields) - set(header)
        if missing:
            raise ValueError(f'{source}: missing required fields {missing}')
        indices = [header.index(f) for f in fields]
        out.write(('\t'.join(rename or fields)+'\n').encode())
        n = 0
        for n, line in enumerate(src, 1):
            parts = line.rstrip(b'\r\n').split(b'\t', max(indices)+1)
            out.write(b'\t'.join(parts[i] for i in indices)+b'\n')
            if n % 50000 == 0:
                print(source.name, 'rows:', n, flush=True)
    os.replace(temp, dest)
    meta_path.write_text(json.dumps(meta, indent=2), encoding='utf-8')
    print('Wrote:', dest, 'rows:', n, flush=True)



def prepare_family_history():
    """Extract baseline family-history fields absent from older prepared UKB files."""
    import argparse
    parser = argparse.ArgumentParser(description=prepare_family_history.__doc__)
    parser.add_argument("--source", default=str(ukb_dir() / "rap/raw/pheno.tab.gz"))
    parser.add_argument("--output", default=str(output_dir() / "ukb/00data/inputs/family_history.tsv"))
    parser.add_argument("--dictionary", default=str(definition_file("vip")))
    args = parser.parse_args()
    variables = ["fhist_father", "fhist_mother", "fhist_siblings"]
    extract_columns(args.source, args.output,
                    ["eid"] + vip_fields(variables, args.dictionary),
                    ["eid"] + variables)


def run_cli(steps):
    """dose.sh supplies DOSE_STEP without changing a step's own CLI arguments."""
    step = os.getenv("DOSE_STEP", "")
    if step not in steps:
        raise SystemExit("Run through dose.sh; available steps: " + ", ".join(steps))
    steps[step]()


# ===========================================================================
@lru_cache(maxsize=1)
def ml_tools():
    """Shared utilities for the final DOSE s2 machine-learning module."""

    from pathlib import Path
    pass  # Shared path functions are defined above.
    import os
    import numpy as np
    import pandas as pd

    from sklearn.metrics import (
        roc_auc_score, average_precision_score, roc_curve,
        confusion_matrix, log_loss, brier_score_loss
    )
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.svm import LinearSVC
    from lightgbm import LGBMClassifier

    try:
        from xgboost import XGBClassifier
        HAVE_XGB = True
    except Exception:
        XGBClassifier = None
        HAVE_XGB = False

    try:
        from catboost import CatBoostClassifier
        HAVE_CAT = True
    except Exception:
        CatBoostClassifier = None
        HAVE_CAT = False

    SEED = 12345
    # Prespecified primary development outcome; do not switch based on post hoc results.
    TARGET = "cvd_stroke_i.Yt2e"
    TIME = "cvd_stroke_i.t2e"
    META_COLS = {"eid", "region", "region_code", TARGET, TIME}

    PROJECT = output_dir()
    OUT = PROJECT / "ukb/03ml"
    INFILE = OUT / "ML_input.csv"
    THREADS = jobs()


    def load_ml():
        df = pd.read_csv(INFILE)
        features = [c for c in df.columns if c not in META_COLS]
        folds = sorted(df["region_code"].dropna().astype(int).unique().tolist())
        if folds != list(range(10)):
            raise RuntimeError(f"region_code must be exactly 0..9; observed {folds}")
        if not features:
            raise RuntimeError("No food-group features found in ML_input.csv")
        if df[TARGET].isna().any():
            raise RuntimeError("Target contains missing values")
        if df[features].isna().any().any():
            bad = df[features].isna().sum()
            raise RuntimeError(
                "Food features contain missing values: " +
                ", ".join(f"{k}={v}" for k, v in bad[bad > 0].items())
            )
        return df, features, folds


    def pos_weight(y):
        y = np.asarray(y, dtype=int)
        n1 = np.sum(y == 1)
        n0 = np.sum(y == 0)
        return float(n0 / n1) if n1 > 0 else 1.0


    def make_model(name, y_train):
        if name == "LightGBM":
            return LGBMClassifier(
                n_estimators=1000, learning_rate=0.01, max_depth=-1,
                num_leaves=31, subsample=0.8, subsample_freq=1,
                colsample_bytree=1.0, reg_alpha=0.1, reg_lambda=1.0,
                min_child_samples=20, objective="binary", is_unbalance=True,
                random_state=SEED, n_jobs=THREADS, verbosity=-1,
            )
        if name == "XGBoost":
            if not HAVE_XGB:
                raise ImportError("xgboost is not installed/importable")
            return XGBClassifier(
                n_estimators=1000, learning_rate=0.01, max_depth=6,
                min_child_weight=1, subsample=0.8, colsample_bytree=0.8,
                reg_alpha=0.0, reg_lambda=1.0, objective="binary:logistic",
                eval_metric="logloss", scale_pos_weight=pos_weight(y_train),
                random_state=SEED, n_jobs=THREADS, verbosity=0,
            )
        if name == "RandomForest":
            return RandomForestClassifier(
                n_estimators=500, max_features="sqrt", min_samples_leaf=1,
                class_weight="balanced_subsample", random_state=SEED,
                n_jobs=THREADS,
            )
        if name == "Logistic":
            return LogisticRegression(
                C=1.0, class_weight="balanced", solver="lbfgs",
                max_iter=5000, random_state=SEED,
            )
        if name == "SVM":
            return LinearSVC(
                C=1.0, class_weight="balanced", dual="auto",
                max_iter=10000, random_state=SEED,
            )
        if name == "CatBoost":
            if not HAVE_CAT:
                raise ImportError("catboost is not installed/importable")
            return CatBoostClassifier(
                iterations=1000, learning_rate=0.01, depth=6,
                loss_function="Logloss", eval_metric="AUC",
                auto_class_weights="Balanced", random_seed=SEED,
                thread_count=THREADS, verbose=False,
            )
        raise ValueError(f"Unknown model: {name}")


    def make_rfe_model(name, y):
        """Lighter estimators for secondary RFE robustness analysis."""
        if name == "LightGBM":
            return LGBMClassifier(
                n_estimators=300, learning_rate=0.03, num_leaves=31,
                subsample=0.8, subsample_freq=1, colsample_bytree=1.0,
                objective="binary", is_unbalance=True,
                random_state=SEED, n_jobs=THREADS, verbosity=-1,
            )
        if name == "XGBoost":
            if not HAVE_XGB:
                raise ImportError("xgboost is not installed/importable")
            return XGBClassifier(
                n_estimators=300, learning_rate=0.03, max_depth=6,
                subsample=0.8, colsample_bytree=0.8,
                objective="binary:logistic", eval_metric="logloss",
                scale_pos_weight=pos_weight(y), random_state=SEED,
                n_jobs=THREADS, verbosity=0,
            )
        if name == "RandomForest":
            return RandomForestClassifier(
                n_estimators=300, max_features="sqrt",
                class_weight="balanced_subsample", random_state=SEED,
                n_jobs=THREADS,
            )
        if name == "Logistic":
            return LogisticRegression(
                C=1.0, class_weight="balanced", solver="lbfgs",
                max_iter=5000, random_state=SEED,
            )
        if name == "SVM":
            return LinearSVC(
                C=1.0, class_weight="balanced", dual="auto",
                max_iter=10000, random_state=SEED,
            )
        raise ValueError(name)


    def model_names(include_svm=True, include_catboost=True):
        names = ["LightGBM"]
        if HAVE_XGB:
            names.append("XGBoost")
        names += ["RandomForest", "Logistic"]
        if include_svm:
            names.append("SVM")
        if include_catboost and HAVE_CAT:
            names.append("CatBoost")
        return names


    def predict_score(model, X):
        if hasattr(model, "predict_proba"):
            p = model.predict_proba(X)
            return np.asarray(p[:, 1], dtype=float), True
        if hasattr(model, "decision_function"):
            return np.asarray(model.decision_function(X), dtype=float), False
        return np.asarray(model.predict(X), dtype=float), False


    def calc_metrics(y, score, is_probability=False):
        y = np.asarray(y, dtype=int)
        score = np.asarray(score, dtype=float)
        fpr, tpr, thresholds = roc_curve(y, score)
        j = tpr - fpr
        i = int(np.nanargmax(j))
        th = float(thresholds[i])
        pred = (score >= th).astype(int)
        tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
        se = tp / (tp + fn) if (tp + fn) else np.nan
        sp = tn / (tn + fp) if (tn + fp) else np.nan
        acc = (tp + tn) / (tp + tn + fp + fn)
        youden = se + sp - 1
        out = {
            "AUC": roc_auc_score(y, score),
            "AUPRC": average_precision_score(y, score),
            "Threshold_Youden": th,
            "Sensitivity_Youden": se,
            "Specificity_Youden": sp,
            "Accuracy_Youden": acc,
            "Youden_Index": youden,
            "BalancedAccuracy_Youden": np.nanmean([se, sp]),
            "EventRate": float(np.mean(y)),
        }
        if is_probability:
            p = np.clip(score, 1e-15, 1 - 1e-15)
            out["LogLoss"] = log_loss(y, p, labels=[0, 1])
            out["Brier"] = brier_score_loss(y, p)
        else:
            out["LogLoss"] = np.nan
            out["Brier"] = np.nan
        return out


    def normalize_importance(x):
        x = np.asarray(x, dtype=float)
        x[~np.isfinite(x)] = 0
        s = x.sum()
        return np.zeros_like(x) if s <= 0 else x / s
    return SimpleNamespace(**locals(), __file__=__file__)


# ===========================================================================
@lru_cache(maxsize=1)
def mri_tracts():
    """True UKB group-average tract surfaces coloured by observed association beta.
    Requires numpy pandas scipy nibabel scikit-image matplotlib vtk (9.4+ for headless Linux).
    """
    from pathlib import Path
    import argparse, json, hashlib, urllib.request
    import numpy as np
    import pandas as pd
    import nibabel as nib
    from scipy.ndimage import gaussian_filter, binary_fill_holes, label
    from skimage.measure import marching_cubes
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize, TwoSlopeNorm
    KEYS='ar_l ar_r atr_l atr_r cgc_l cgc_r cgh_l cgh_r cst_l cst_r fma fmi ifo_l ifo_r ilf_l ilf_r mcp ml_l ml_r ptr_l ptr_r slf_l slf_r str_l str_r unc_l unc_r'.split()
    BASE='https://www.fmrib.ox.ac.uk/ukbiobank/group_means/ukb/'
    NAMES=['TractsSummed_Tracts1-27.nii.gz','T1_T2FLAIR_fMRI.nii.gz']
    MASK_URL='https://raw.githubusercontent.com/Washington-University/HCPpipelines/master/global/templates/MNI152_T1_1mm_brain_mask.nii.gz'

    def mesh(a,affine,level,step=2):
        vertices,faces,_,_=marching_cubes(a,level=level,step_size=step)
        vertices=nib.affines.apply_affine(affine,vertices)
        pts=vtk.vtkPoints();pts.SetData(numpy_to_vtk(vertices,deep=True))
        cells=vtk.vtkCellArray()
        cells.ImportLegacyFormat(numpy_to_vtk(np.c_[np.full(len(faces),3),faces].ravel().astype(np.int64),deep=True,array_type=vtk.VTK_ID_TYPE))
        poly=vtk.vtkPolyData();poly.SetPoints(pts);poly.SetPolys(cells)
        smooth=vtk.vtkWindowedSincPolyDataFilter();smooth.SetInputData(poly);smooth.SetNumberOfIterations(15)
        smooth.SetPassBand(.12);smooth.BoundarySmoothingOff();smooth.NormalizeCoordinatesOn();smooth.Update()
        normals=vtk.vtkPolyDataNormals();normals.SetInputConnection(smooth.GetOutputPort());normals.SplittingOff();normals.Update()
        result=vtk.vtkPolyData();result.DeepCopy(normals.GetOutput());return result

    def actor(poly,color,opacity=1):
        mapper=vtk.vtkPolyDataMapper();mapper.SetInputData(poly);mapper.ScalarVisibilityOff()
        a=vtk.vtkActor();a.SetMapper(mapper);p=a.GetProperty();p.SetColor(*color[:3]);p.SetOpacity(opacity)
        p.SetAmbient(.25);p.SetDiffuse(.75);p.SetSpecular(.16);p.SetSpecularPower(20)
        return a

    BETA_LIMIT=.02
    TRACT_CMAP = TRACT_NORM = None
    def render(brain,meshes,rows,direction,size=(1200,1000)):
        win=vtk.vtkRenderWindow();win.SetOffScreenRendering(1);win.SetSize(*size);win.SetMultiSamples(0);win.SetAlphaBitPlanes(1)
        ren=vtk.vtkRenderer();ren.SetBackground(1,1,1);win.AddRenderer(ren)
        ren.SetUseDepthPeeling(1);ren.SetMaximumNumberOfPeels(100);ren.SetOcclusionRatio(.1)
        ren.AddActor(actor(brain,(.76,.79,.82),.09))
        for _,r in rows.iterrows():
            color=TRACT_CMAP(TRACT_NORM(r.beta))
            ren.AddActor(actor(meshes[r.atlas_key],color,.98))
        cam=ren.GetActiveCamera();center=np.array([0,-18,10]);direction=np.asarray(direction,float)
        cam.SetFocalPoint(*center);cam.SetPosition(*(center+direction/np.linalg.norm(direction)*450))
        cam.SetViewUp(*( (0,1,0) if abs(direction[2])>abs(direction[0])+abs(direction[1]) else (0,0,1)))
        cam.ParallelProjectionOn();cam.SetParallelScale(105);ren.ResetCameraClippingRange()
        win.Render();f=vtk.vtkWindowToImageFilter();f.SetInput(win);f.SetInputBufferTypeToRGB();f.ReadFrontBufferOff();f.Update()
        im=f.GetOutput();w,h,_=im.GetDimensions();array=vtk_to_numpy(im.GetPointData().GetScalars()).reshape(h,w,3)[::-1].copy()
        win.Finalize();return array

    def save(fig,out,name):
        for ext in ['png','pdf']:save_png_figure(fig, out/(name+'.'+ext),dpi=300,bbox_inches='tight',facecolor='white')
        plt.close(fig)

    def generate(r,out,atlas,threshold=1000):
        from types import SimpleNamespace
        args=SimpleNamespace(threshold=threshold)
        out=Path(out);atlas=Path(atlas)
        maskfile=atlas/'MNI152_T1_1mm_brain_mask.nii.gz'
        r=r[r.exposure.eq('dose') & r.family.isin(['FA27','MD27'])].copy()
        expected=[line.split()[0] for line in (atlas/'structureList.source.txt').read_text().splitlines() if line.strip()]
        if expected!=KEYS:raise ValueError('AutoPtx source order differs from renderer')
        if r.duplicated(['family','atlas_key']).any():raise ValueError('Duplicate tract results')
        for fam in ['FA27','MD27']:
            h=r[r.family.eq(fam)]
            if len(h)!=27 or set(h.atlas_key)!=set(KEYS):raise ValueError('Require all 27 named tracts')
            from scipy.stats import false_discovery_control
            if not np.allclose(false_discovery_control(h.P),h.FDR_family,atol=1e-10):raise ValueError('FDR mismatch')
        selected=r[r.status.eq('ok') & r.FDR_family.lt(.05)]
        nonlocal BETA_LIMIT, TRACT_CMAP, TRACT_NORM
        BETA_LIMIT=max(.02,float(selected.beta.abs().max())) if len(selected) else .02
        mixed=(selected.beta<0).any()
        TRACT_CMAP=plt.get_cmap('RdBu_r' if mixed else 'Reds')
        TRACT_NORM=TwoSlopeNorm(vmin=-BETA_LIMIT,vcenter=0,vmax=BETA_LIMIT) if mixed else Normalize(0,BETA_LIMIT)
        selected.to_csv(out/'displayed_tract_results.csv',index=False)
        im=nib.load(atlas/NAMES[0]);bg=nib.load(atlas/NAMES[1])
        if im.shape!=(182,218,182,28) or not np.allclose(im.affine,bg.affine):raise ValueError('Unexpected atlas geometry')
        # Volume 0 is the sum, volumes 1..27 follow official AutoPtx order.
        mapping=[];meshes={};manifest=[]
        masks=out/'cache/selected_tract_masks';masks.mkdir(parents=True,exist_ok=True)
        for i,key in enumerate(KEYS,1):
            a=np.asarray(im.dataobj[...,i],dtype=np.float32)
            xyz=nib.affines.apply_affine(im.affine,np.argwhere(a>a.max()*.15)).mean(0)
            if key.endswith('_l') and xyz[0]>=0 or key.endswith('_r') and xyz[0]<=0:raise ValueError('Hemisphere check failed')
            mapping.append(dict(atlas_key=key,nifti_index_zero_based=i,centroid_x=xyz[0],centroid_y=xyz[1],centroid_z=xyz[2],intensity_max=float(a.max())))
            if key in set(selected.atlas_key):
                if not 0<args.threshold<a.max():raise ValueError('Invalid isosurface threshold')
                meshes[key]=mesh(gaussian_filter(a,.7),im.affine,args.threshold)
                scale=float(a.max());path=masks/(key+'.nii.gz');nib.save(nib.Nifti1Image(a/scale,im.affine),path)
                manifest.append(dict(atlas_key=key,mask_path=str(path.relative_to(out)),threshold=args.threshold/scale,space='MNI152',source=BASE+NAMES[0]))
        pd.DataFrame(mapping).to_csv(out/'atlas_mapping_audit.csv',index=False)
        pd.DataFrame(manifest).to_csv(out/'verified_tract_manifest.csv',index=False)
        t1=np.asarray(bg.dataobj[...,0],dtype=np.float32)
        bm=nib.load(maskfile)
        from nibabel.processing import resample_from_to
        bm=resample_from_to(bm,(bg.shape[:3],bg.affine),order=0)
        mask=(bm.get_fdata()>.5) & (t1>45)
        labs,n=label(mask);counts=np.bincount(labs.ravel());counts[0]=0;mask=labs==counts.argmax()
        brain=mesh(gaussian_filter(mask.astype(np.float32),1.0),bg.affine,.5)
        plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42})
        # Main panel b: no colour is shown for nonsignificant MD.
        fig,axs=plt.subplots(1,2,figsize=(13,6.3))
        for ax,fam in zip(axs,['FA27','MD27']):
            rows=selected[selected.family.eq(fam)]
            rendered=render(brain,meshes,rows,(-1,-.35,.18))
            cache=out/'cache/tract_renders';cache.mkdir(parents=True,exist_ok=True)
            plt.imsave(cache/(fam+'.png'),rendered)
            ax.imshow(rendered);ax.axis('off')
            ax.set_title(f'{fam[:2]} | {len(rows)} of 27 tracts',fontsize=16,pad=4)
            if rows.empty:
                ax.text(.5,.52,'No association survives\nFDR correction',transform=ax.transAxes,ha='center',va='center',fontsize=16,
                        bbox=dict(facecolor='white',alpha=.9,edgecolor='none',pad=8))
        fig.text(.025,.97,'b',fontsize=25,fontweight='bold')
        cax=fig.add_axes([.15,.15,.27,.024]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=TRACT_NORM,cmap=TRACT_CMAP),cax=cax,orientation='horizontal')
        cb.set_label('Standardized beta',fontsize=11);cb.set_ticks([-BETA_LIMIT,0,BETA_LIMIT] if mixed else [0,BETA_LIMIT/2,BETA_LIMIT])
        fig.text(.5,.035,'UKB group-average tract surfaces | Within-modality FDR < 0.05',ha='center',fontsize=11)
        fig.subplots_adjust(left=.02,right=.98,top=.9,bottom=.22,wspace=.02);save(fig,out,'b_tract_anatomy')
        # Additional viewpoints expose overlap hidden by one camera.
        fa=selected[selected.family.eq('FA27')]
        fig,axs=plt.subplots(1,3,figsize=(17,6))
        for ax,direction,title in zip(axs,[(-1,-.35,.18),(1,-.35,.18),(0,0,1)],['Left oblique','Right oblique','Superior']):
            ax.imshow(render(brain,meshes,fa,direction));ax.axis('off');ax.set_title(title,fontsize=15)
        fig.suptitle(f'DOSE and white-matter FA | {len(fa)} FDR-significant tracts',fontsize=18)
        fig.text(.5,.045,'Red: positive beta; blue: negative beta. Tract overlap varies with viewpoint.',ha='center')
        fig.subplots_adjust(left=0,right=1,bottom=.12,top=.87,wspace=.01);save(fig,out,'FA_three_views')
        # Separate groups make all six anatomical structures identifiable.
        groups=[(['fma'],'Forceps major'),(['ifo_l','ifo_r'],'Inferior fronto-occipital fasciculus'),
                (['ilf_l','ilf_r'],'Inferior longitudinal fasciculus'),(['ptr_l','ptr_r'],'Posterior thalamic radiation'),
                (['slf_l','slf_r'],'Superior longitudinal fasciculus'),(['unc_l','unc_r'],'Uncinate fasciculus')]
        fig,axs=plt.subplots(2,3,figsize=(15,10))
        for ax,(keys,title) in zip(axs.flat,groups):
            rows=fa[fa.atlas_key.isin(keys)];ax.imshow(render(brain,meshes,rows,(-1,-.45,.25)));ax.axis('off')
            ax.set_title(title,fontsize=12)
            text='; '.join(('L' if v.atlas_key.endswith('_l') else 'R' if v.atlas_key.endswith('_r') else '')+f' beta={v.beta:.4f}' for _,v in rows.iterrows())
            ax.text(.5,.06,text,transform=ax.transAxes,ha='center',fontsize=10)
        fig.suptitle('Significant FA associations by tract group',fontsize=18)
        fig.subplots_adjust(left=.01,right=.99,bottom=.03,top=.92,wspace=.01,hspace=.06);save(fig,out,'FA_tract_groups')
        info={'atlas_source':BASE+NAMES[0],'background_source':BASE+NAMES[1],'background_brain_mask':MASK_URL,
              'license':'CC BY 4.0; University of Oxford / UK Biobank group-average resources',
              'threshold_atlas_intensity':args.threshold,'smoothing_sigma_voxels':.7,
              'geometry':'Marching-cubes surfaces of real group-average tract maps, not individual streamlines',
              'FA_significant':int((selected.family=='FA27').sum()),'MD_significant':int((selected.family=='MD27').sum()),
              'sha256':{name:hashlib.sha256((atlas/name).read_bytes()).hexdigest() for name in NAMES}}
        (out/'tract_provenance.json').write_text(json.dumps(info,indent=2))
        print('Finished real atlas-based tract rendering.',flush=True)
    return SimpleNamespace(**locals(), __file__=__file__)


# ===========================================================================
@lru_cache(maxsize=1)
def prepare_enigma():
    """Allow ENIGMA's wrappers to skip VTK's newer Python-only mixin classes."""
    import numpy as np
    if not hasattr(np, "string_"): np.string_ = np.bytes_
    if not hasattr(np, "unicode_"): np.unicode_ = np.str_
    from enigmatoolbox.vtk_interface.wrappers import base

    def wrap(obj):
        if obj is None or base.is_wrapper(obj):
            return obj
        if isinstance(obj, type):
            obj = obj()
        if not base.is_vtk(obj):
            raise ValueError(f"Unknown VTK object type: {type(obj)}")
        for cls in type(obj).__mro__:
            name = getattr(cls, "__vtkname__", "")
            wrapper = base.BSVTKObjectWrapperMeta.entries.get(name[3:]) if name.startswith("vtk") else None
            if wrapper is not None:
                return wrapper(obj)
        return base.BSVTKObjectWrapper(obj)

    base.BSWrapVTKObject = wrap

    def get_attribute(self, name):
        # Recent VTK exposes properties alongside Get* methods. Only methods
        # should become call wrappers; properties retain their actual values.
        try:
            value = getattr(self.VTKObject, name)
        except AttributeError:
            return self._handle_call("get", name, None)
        return base.VTKMethodWrapper(value) if callable(value) and not base.is_vtk(value) else base._wrap_output_data(value)

    base.BSVTKObjectWrapper.__getattr__ = get_attribute

    import vtk
    if not hasattr(vtk.vtkRenderer, "AddActor2D"):
        from enigmatoolbox.vtk_interface.wrappers.renderer import BSViewport
        from enigmatoolbox.vtk_interface.wrappers.actor import BSActor2D, BSScalarBarActor, BSTextActor
        def add_view_prop(actor_type):
            def add(self, obj=None, **kwargs):
                actor = actor_type(vtkobject=obj, **kwargs)
                self.VTKObject.AddViewProp(actor.VTKObject)
                return actor
            return add
        for method, actor_type in (("AddActor2D", BSActor2D), ("AddScalarBarActor", BSScalarBarActor), ("AddTextActor", BSTextActor)):
            setattr(BSViewport, method, add_view_prop(actor_type))


@lru_cache(maxsize=1)
def mri_volume():
    """DOSE regional volume associations; retain the existing 139-test FDR family."""
    import argparse
    import json
    import hashlib
    import os
    import re
    import urllib.request
    import xml.etree.ElementTree as ET
    from pathlib import Path
    import numpy as np
    import pandas as pd
    import nibabel as nib
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from skimage.measure import marching_cubes

    ASSETS = {
     'Cortical': ('HarvardOxford-Cortical.xml', 'HarvardOxford/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz'),
     'Subcortical': ('HarvardOxford-Subcortical.xml', 'HarvardOxford/HarvardOxford-sub-maxprob-thr25-2mm.nii.gz'),
     'Cerebellar': ('Cerebellum_MNIfnirt.xml', 'Cerebellum/Cerebellum-MNIfnirt-maxprob-thr25-2mm.nii.gz')}
    MIRROR='https://raw.githubusercontent.com/fithisux/FSL/master/data/atlases/'

    def norm(s):
     return re.sub(r'[^a-z0-9]', '', str(s).lower())

    def target_name(row):
     s=re.sub(r'^Volume of (grey matter in )?', '', row.label)
     s=re.sub(r'\s*\((left|right|vermis)\)$', '', s, flags=re.I)
     if row.category=='Cortical': return s
     if row.category=='Brainstem': return 'Brain-Stem'
     if row.category=='Subcortical': return f'{row.side} {s}'
     if row.category=='Cerebellar': return f'{row.side} {s.replace(" Cerebellum", "")}'
     raise ValueError(f'Unknown category: {row.category}')

    def read_results(path, exposure, model, qcol):
     d=path.copy() if isinstance(path,pd.DataFrame) else pd.read_csv(path)
     required={'exposure','model','family','field_id','category','label','side','beta','P',qcol}
     if not required.issubset(d): raise ValueError(f'Missing columns: {required-set(d)}')
     d=d[(d.exposure==exposure)&(d.model==model)&(d.family=='Volume139')].copy()
     if len(d)!=139 or d.field_id.nunique()!=139:
      raise ValueError(f'Expected 139 unique Volume139 rows, found {len(d)}; do not silently select a subset.')
     for c in ['beta','P',qcol]: d[c]=pd.to_numeric(d[c],errors='raise')
     if not np.isfinite(d[['beta','P',qcol]].to_numpy()).all(): raise ValueError('Non-finite results; resolve failed models before plotting.')
     if not d[['P',qcol]].apply(lambda x:x.between(0,1).all()).all(): raise ValueError('Invalid P/q values.')
     # Audit only, never replace supplied q values.
     p=d.P.to_numpy(); order=np.argsort(p); q=np.empty(len(p)); q[order]=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1].clip(0,1)
     if qcol=='FDR_family' and not np.allclose(q,d[qcol],atol=1e-8,rtol=1e-5):
      raise ValueError('FDR_family differs from BH across these 139 tests. Check the analysis family.')
     return d

    def load_atlases(root, download):
     ans={}
     for group,(xml,nii) in ASSETS.items():
      for rel in (xml,nii):
       p=root/Path(rel).name
       if not p.exists():
        if not download: raise FileNotFoundError(f'{p}\nSet --atlas-dir to FSL data/atlases, or use --download-atlases (FSL mirror).')
        p.parent.mkdir(parents=True,exist_ok=True)
        tmp=p.with_name(p.name+'.partial')
        print('Download:',rel,flush=True)
        urllib.request.urlretrieve(MIRROR+rel,tmp); tmp.replace(p)
      tree=ET.parse(root/xml)
      offset=1 if tree.findtext('./header/type').strip()=='Probabilistic' else 0
      labels={norm(e.text):(int(e.attrib['index'])+offset,e.text.strip()) for e in tree.findall('./data/label')}
      img=nib.load(root/Path(nii).name); data=np.asarray(img.dataobj)
      if data.ndim!=3 or not np.allclose(data,np.round(data)): raise ValueError('Expected discrete 3D atlas')
      ans[group]=(img,data.astype(np.int16),labels)
     return ans

    def make_mapping(d, atlases):
     rows=[]
     for r in d.itertuples():
      group='Subcortical' if r.category=='Brainstem' else r.category
      img,data,labels=atlases[group]; name=target_name(r)
      # Explicit spelling aliases only; no fuzzy matching.
      key=norm(name)
      aliases={'leftaccumbens':'leftaccumbens','rightaccumbens':'rightaccumbens'}
      key=aliases.get(key,key)
      if key not in labels: raise ValueError(f'No exact atlas match: {r.field_id}: {name}; candidates={list(labels)}')
      code,matched=labels[key]; mask=data==code
      if group=='Cortical':
       ijk=np.argwhere(mask); x=nib.affines.apply_affine(img.affine,ijk)[:,0]
       keep=x<0 if r.side=='Left' else x>0
       mask[:]=False; ijk=ijk[keep]; mask[tuple(ijk.T)]=True
      if not mask.any(): raise ValueError(f'Empty ROI: {r.field_id}')
      centroid=nib.affines.apply_affine(img.affine,np.argwhere(mask)).mean(0)
      if r.side=='Left' and centroid[0]>=0 or r.side=='Right' and centroid[0]<=0:raise ValueError(f'Hemisphere mismatch: {r.field_id}')
      rows.append(dict(centroid_x=float(centroid[0]),field_id=r.field_id,atlas_group=group,atlas_label=matched,atlas_code=code,n_voxels=int(mask.sum()),side=r.side))
     return d.merge(pd.DataFrame(rows),on=['field_id','side'],validate='one_to_one')

    def roi_mesh(img, data, code):
     mask=data==code
     inds=np.argwhere(mask); lo=np.maximum(inds.min(0)-1,0); hi=np.minimum(inds.max(0)+2,data.shape)
     crop=mask[tuple(slice(a,b) for a,b in zip(lo,hi))]
     # Padding ensures closed meshes even if an ROI touches a crop boundary.
     v,f,_,_=marching_cubes(np.pad(crop,1).astype(float),.5)
     v=nib.affines.apply_affine(img.affine,v+lo-1)
     return v,f

    def setup_axis(ax,xyz,elev,azim,title):
     lo=xyz.min(0); hi=xyz.max(0); centre=(hi+lo)/2; radius=max(hi-lo)/2*1.06
     ax.set(xlim=(centre[0]-radius,centre[0]+radius),ylim=(centre[1]-radius,centre[1]+radius),zlim=(centre[2]-radius,centre[2]+radius))
     ax.set_box_aspect((1,1,1),zoom=1.5);ax.set_proj_type('ortho');ax.view_init(elev=elev,azim=azim);ax.set_axis_off();ax.set_title(title,fontsize=11,pad=0)

    def collection(v,f,colors):
     c=Poly3DCollection(v[f],facecolors=colors,edgecolors=(0,0,0,0),linewidth=0,rasterized=True,shade=True)
     return c

    def render(d,atlases,args,out):
     from nilearn import datasets,surface
     scale=args.vmax or float(np.ceil(np.abs(d.beta).max()*1000)/1000)
     if scale<=0 or scale<np.abs(d.beta).max(): raise ValueError('--vmax must include all 139 beta values.')
     cm=plt.get_cmap('RdBu_r'); nm=Normalize(-scale,scale)
     neutral=(.83,.83,.83,1); background=(.95,.95,.95,1)
     def color(r): return cm(nm(r.beta)) if args.show_all or getattr(r,args.q_col)<args.alpha else neutral
     fs=datasets.fetch_surf_fsaverage('fsaverage5',data_dir=args.cache_dir)
     fig=plt.figure(figsize=(15,11));gs=fig.add_gridspec(3,4,left=.01,right=.99,top=.92,bottom=.10,hspace=.14,wspace=.01)
     cortex=atlases['Cortical'][0]; projected=[]
     for col,(side,view,azim) in enumerate([('Left','lateral',180),('Left','medial',0),('Right','medial',180),('Right','lateral',0)]):
      hemi=side.lower();pial=surface.load_surf_mesh(fs['pial_'+hemi]);mesh=surface.load_surf_mesh(fs['infl_'+hemi])
      labels=surface.vol_to_surf(cortex,pial,inner_mesh=fs['white_'+hemi],interpolation='nearest_most_frequent',depth=[.1,.3,.5,.7,.9])
      labels=np.nan_to_num(labels).astype(int)
      face_labels=labels[mesh.faces];a,b,c=face_labels.T
      # Majority vertex label, stable lowest label for a three-way boundary.
      codes=np.where((a==b)|(a==c),a,np.where(b==c,b,np.min(face_labels,axis=1)))
      colors=np.tile(background,(len(codes),1))
      rows=d[(d.category=='Cortical')&(d.side==side)]
      for r in rows.itertuples():
       colors[codes==r.atlas_code]=color(r)
       if view=='lateral': projected.append({'field_id':r.field_id,'surface_faces':int((codes==r.atlas_code).sum())})
      ax=fig.add_subplot(gs[0,col],projection='3d');ax.add_collection3d(collection(mesh.coordinates,mesh.faces,colors));setup_axis(ax,mesh.coordinates,0,azim,f'{side} {view}')
     pd.DataFrame(projected).to_csv(out/'cortical_projection_coverage.csv',index=False)
     absent=pd.DataFrame(projected).query('surface_faces == 0')
     if len(absent): print('WARNING: cortical ROIs absent from projection:',absent.field_id.tolist(),flush=True)
     for row,group in [(1,'Subcortical'),(2,'Cerebellar')]:
      sub=d[d.atlas_group==group];img,arr,_=atlases[group]
      meshes={r.field_id:roi_mesh(img,arr,r.atlas_code) for r in sub.itertuples()}
      views=[('Left structures',0,180,'Left'),('Right structures',0,0,'Right'),('Posterior',0,-90,None),('Inferior',-90,-90,None)] if row==1 else [('Left lateral',0,180,None),('Right lateral',0,0,None),('Posterior',0,-90,None),('Inferior',-90,-90,None)]
      for col,(title,elev,azim,side) in enumerate(views):
       ax=fig.add_subplot(gs[row,col],projection='3d')
       chosen=sub if side is None else sub[sub.side==side]
       # One collection for all meshes allows depth sorting across structures.
       polys=[];colors=[];points=[]
       for r in chosen.itertuples():
        v,f=meshes[r.field_id];polys.extend(v[f]);colors.extend([color(r)]*len(f));points.append(v)
       coll=Poly3DCollection(polys,facecolors=colors,edgecolors=(0,0,0,0),linewidth=0,rasterized=True,shade=True)
       ax.add_collection3d(coll);setup_axis(ax,np.concatenate(points),elev,azim,title)
     fig.text(.025,.965,'a  Cortical volume',weight='bold',fontsize=14)
     fig.text(.025,.663,'b  Subcortical and brainstem volume',weight='bold',fontsize=14)
     fig.text(.025,.38,'c  Cerebellar volume',weight='bold',fontsize=14)
     cax=fig.add_axes([.35,.072,.30,.016]);fig.colorbar(ScalarMappable(norm=nm,cmap=cm),cax=cax,orientation='horizontal',ticks=np.linspace(-scale,scale,5),label='Standardized β (positive = red; negative = blue)')
     selection='All coefficients' if args.show_all else f'Coloured: {args.q_col} < {args.alpha:g}; grey: not significant'
     fig.text(.5,.008,selection+' | Existing FDR retained across all 139 volumes',ha='center',fontsize=10)
     for ext in ['png','pdf']: save_png_figure(fig, out/f'volume_brain_maps.{ext}',dpi=300)
     plt.close(fig)
     # A complete significant-ROI forest plot preserves effects hidden in 3D views.
     sig=d[d[args.q_col]<args.alpha].sort_values(['category','field_id'])
     if len(sig) and {'LCI','HCI'}.issubset(sig):
      fig,ax=plt.subplots(figsize=(10,max(3,.34*len(sig)+1)))
      for i,r in enumerate(sig.itertuples()):
       ax.plot([r.LCI,r.HCI],[i,i],color=cm(nm(r.beta)));ax.plot(r.beta,i,'o',color=cm(nm(r.beta)))
      ax.set_yticks(range(len(sig)),[re.sub(r'^Volume of (grey matter in )?','',s) for s in sig.label]);ax.invert_yaxis();ax.axvline(0,color='.5',lw=.8);ax.set_xlabel('Standardized β (95% CI)');ax.set_title(f'{args.exposure}: {len(sig)} / 139 volumes, FDR < {args.alpha:g}')
      fig.tight_layout()
      for ext in ['png','pdf']:save_png_figure(fig, out/f'volume_significant_forest.{ext}',dpi=300)
      plt.close(fig)

    def generate(results,out,atlas_dir):
     from types import SimpleNamespace
     out=Path(out);atlas_dir=Path(atlas_dir)
     args=SimpleNamespace(exposure='dose',model='M2',q_col='FDR_family',alpha=.05,
        vmax=None,show_all=False,cache_dir=str(out/'cache/nilearn'))
     d=read_results(results,args.exposure,args.model,args.q_col)
     atlases=load_atlases(atlas_dir,False);d=make_mapping(d,atlases)
     d.to_csv(out/'volume_mapping_and_results.csv',index=False)
     counts=d.assign(significant=d[args.q_col]<args.alpha).groupby('category').agg(tested=('field_id','size'),significant=('significant','sum'))
     counts.to_csv(out/'volume_summary.csv');print(counts.to_string(),flush=True)
     render(d,atlases,args,out)
     files=[atlas_dir/Path(name).name for pair in ASSETS.values() for name in pair]
     info=dict(vars(args),atlas_dir=str(atlas_dir),
        source='FSL Harvard-Oxford cortical/subcortical and Cerebellum MNI-FNIRT atlases',
        documentation='https://fsl.fmrib.ox.ac.uk/fsl/docs/other/datasets.html',
        subcortical_interpretation='Same-named reference anatomy, NOT participant FIRST segmentations',
        sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
     (out/'volume_run_settings.json').write_text(json.dumps(info,indent=2),encoding='utf-8')
     (out/'volume_figure_caption.txt').write_text(
        'DOSE and regional brain volumes. Cortical: Harvard-Oxford labels projected to fsaverage5 inflated surfaces; '
        'subcortical/brainstem and cerebellar: surfaces from standard atlas masks. '
        'FIRST-derived volumes use same-named Harvard-Oxford regions for anatomical reference, not individual segmentation boundaries. '
        'Colours are standardized regression beta; grey indicates FDR >= 0.05. '
        'BH correction retained jointly across 139 volumes. Colours represent regional associations, not shape changes or vertexwise tests. '
        'M2 plus baseline age squared and centre, with ICV adjustment. '    'Some internal regions may be occluded in 3D; refer to the full significant-volume forest plot.',encoding='utf-8')
    return SimpleNamespace(**locals(), __file__=__file__)


# ===========================================================================
@lru_cache(maxsize=1)
def mri_focus():
    """Four requested subcortical regions on reference anatomy, with projected labels."""
    import numpy as np
    import pandas as pd
    import nibabel as nib
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    from scipy.ndimage import gaussian_filter
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize,TwoSlopeNorm
    mesh = mri_tracts().mesh; actor = mri_tracts().actor

    def render_focus(d,atlas,out):
        load_atlases = mri_volume().load_atlases; make_mapping = mri_volume().make_mapping
        volumes=d[d.family.eq('Volume139')].copy()
        mapped=make_mapping(volumes,load_atlases(atlas,False))
        chosen=mapped[mapped.category.eq('Subcortical') & mapped.label.str.contains('hippocampus|thalamus',case=False,regex=True)].copy()
        if len(chosen)!=4:raise ValueError('Require bilateral hippocampus and thalamus')
        im=nib.load(atlas/'HarvardOxford-sub-maxprob-thr25-2mm.nii.gz');a=np.asarray(im.dataobj)
        bg=nib.load(atlas/'ukb_tracts27/MNI152_T1_1mm_brain_mask.nii.gz')
        brain=mesh(gaussian_filter(np.asarray(bg.dataobj,dtype=np.float32),1),bg.affine,.5,step=3)
        win=vtk.vtkRenderWindow();win.SetOffScreenRendering(1);win.SetSize(1400,1100);win.SetMultiSamples(0);win.SetAlphaBitPlanes(1)
        ren=vtk.vtkRenderer();ren.SetBackground(1,1,1);win.AddRenderer(ren)
        ren.SetUseDepthPeeling(1);ren.SetMaximumNumberOfPeels(100);ren.SetOcclusionRatio(.05)
        ren.AddActor(actor(brain,(.78,.78,.78),.055))
        limit=float(np.ceil(chosen.beta.abs().max()*1000)/1000)
        signed=(chosen.beta<0).any();nm=TwoSlopeNorm(vmin=-limit,vcenter=0,vmax=limit) if signed else Normalize(0,limit)
        cm=plt.get_cmap('RdBu_r' if signed else 'Reds');centres={}
        for v in chosen.itertuples():
            mask=(a==v.atlas_code).astype(np.float32)
            poly=mesh(gaussian_filter(mask,.6),im.affine,.5,step=1)
            ren.AddActor(actor(poly,cm(nm(v.beta))))
            centres[v.field_id]=np.array(poly.GetCenter())
        cam=ren.GetActiveCamera();cam.SetFocalPoint(0,-20,0);cam.SetPosition(0,-430,240);cam.SetViewUp(0,0,1)
        cam.ParallelProjectionOn();cam.SetParallelScale(112);ren.ResetCameraClippingRange();win.Render()
        audit=[]
        for v in chosen.itertuples():
            ren.SetWorldPoint(*centres[v.field_id],1);ren.WorldToDisplay();x,y,_=ren.GetDisplayPoint()
            name=f'{v.side} '+('hippocampus' if 'hippocampus' in v.label.lower() else 'thalamus')
            q=v.FDR_family;stars='***' if q<.001 else '**' if q<.01 else '*' if q<.05 else ''
            audit.append(dict(field_id=v.field_id,label=name,beta=v.beta,FDR_family=q,stars=stars,
                x=x,y=1100-y,text_x=130 if v.side=='Left' else 1270,
                text_y=290 if 'thalamus' in name else 780,atlas_code=v.atlas_code))
        capture=vtk.vtkWindowToImageFilter();capture.SetInput(win);capture.SetInputBufferTypeToRGB();capture.ReadFrontBufferOff();capture.Update()
        im=capture.GetOutput();w,h,_=im.GetDimensions()
        rgb=vtk_to_numpy(im.GetPointData().GetScalars()).reshape(h,w,3)[::-1].copy();win.Finalize()
        cache=out/'cache/subcortical_focus';cache.mkdir(parents=True,exist_ok=True);plt.imsave(cache/'posterior_superior.png',rgb)
        pd.DataFrame(audit).to_csv(out/'Figure4_subcortical_source.csv',index=False)
        return rgb,audit,limit
    return SimpleNamespace(**locals(), __file__=__file__)


# ===========================================================================
@lru_cache(maxsize=1)
def mri_publication():
    """Compact MRI-only figure for placement alongside clinical outcomes in Figure 4."""
    from pathlib import Path
    import sys
    pass  # Shared path functions are defined above.
    import json
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize, TwoSlopeNorm

    def locate(out,name):
        for folder in ('','audit','main','supplementary'):
            p=out/folder/name
            if p.exists():return p
        raise FileNotFoundError(name)

    def generate(r,out, font_scale=1.0, return_figure=False):
        from matplotlib import font_manager
        # Fail explicitly if Arial is unavailable instead of silently substituting it.
        font_manager.findfont('Arial',fallback_to_default=False)
        plt.rcParams.update({'font.family':'Arial','font.size':10,'pdf.fonttype':42,'ps.fonttype':42})
        out=Path(out)
        d=r[r.exposure.eq('dose') & r.status.eq('ok')].copy()
        mapping=pd.read_csv(locate(out,'panel_a_mapping.csv'))
        source=mapping[mapping.status.eq('ok')].set_index('atlas_key')
        current=d[d.family.eq('ThicknessDK')].set_index('atlas_key')
        if set(source.index)!=set(current.index) or not np.allclose(source.loc[current.index,['beta','FDR_family']],current[['beta','FDR_family']]):
            raise ValueError('Cortical render cache differs from results; run the full plot script first')
        source=pd.read_csv(locate(out,'displayed_tract_results.csv')).set_index(['family','atlas_key'])
        current=d[d.family.isin(['FA27','MD27']) & d.FDR_family.lt(.05)].set_index(['family','atlas_key'])
        if set(source.index)!=set(current.index) or not np.allclose(source.loc[current.index,['beta','FDR_family']],current[['beta','FDR_family']]):
            raise ValueError('Tract render cache differs from results; run the full plot script first')
        hits=d[d.FDR_family.lt(.05)]
        counts={f:int((hits.family==f).sum()) for f in ['ThicknessDK','FA27','MD27','Volume139','WMH']}
        annotations=pd.read_csv(locate(out,'panel_a_annotation_coordinates.csv'))
        fig=plt.figure(figsize=(17,10),facecolor='white')
        fig.text(.025,.978,'DOSE and brain structure',fontsize=20,weight='bold',va='top')
        fig.text(.025,.93,f'a  Cortical thickness  |  {counts["ThicknessDK"]}/66 regions with FDR < 0.05',fontsize=13,weight='bold')
        ax=fig.add_axes([.012,.555,.976,.36])
        cortex_image=plt.imread(out/'cache/surface_renders/a_thickness_FDR05.png')
        # Replace the cached VTK colour-bar lettering with vector Arial text.
        cortex_image=cortex_image.copy()
        cortex_image[:,1430:,:3]=1
        ax.imshow(cortex_image);ax.axis('off')
        fig.canvas.draw()
        cortex_scale=ax.get_position().width*17/cortex_image.shape[1]
        widths=[]
        for start in (0,356,712,1068):
            yy,xx=np.where(np.min(cortex_image[:,start:start+356,:3],axis=2)<.98)
            widths.append((xx.max()-xx.min()+1)*cortex_scale)
        target_width=float(np.median(widths))
        size_audit=[dict(panel='a',brain_width_inches=w) for w in widths]
        def brain_axes(image,cx,cy,panel):
            rgb01=image[:,:,:3].astype(float)
            if image.dtype.kind in 'ui':
                rgb01/=np.iinfo(image.dtype).max
            yy,xx=np.where(np.min(rgb01,axis=2)<.995)
            scale=target_width/(xx.max()-xx.min()+1)
            w=image.shape[1]*scale;h=image.shape[0]*scale
            size_audit.append(dict(panel=panel,brain_width_inches=target_width))
            axis=fig.add_axes([cx-w/34,cy-h/20,w/17,h/10]);axis.imshow(image);axis.axis('off');return axis
        for v in annotations.itertuples():
            label=v.display_label.replace('Precentral gyrus','Precentral\ngyrus') if font_scale>1.25 else v.display_label
            ax.annotate(label,xy=(v.target_x,v.target_y),xytext=(v.text_x,v.text_y),
                ha='center',va='center',fontsize=8.5,
                arrowprops=dict(arrowstyle='-',color='.25',lw=.8,shrinkA=3,shrinkB=0),
                bbox=dict(facecolor='white',edgecolor='none',pad=.5))
        for x,label in zip([178,534,890,1246],['Left lateral','Left medial','Right medial','Right lateral']):
            ax.text(x,489,label,ha='center',fontsize=9,color='.3')
        ax.text(1500,25,'Standardized β',ha='center',fontsize=9)
        cortex_limit=float(d.loc[d.family.eq('ThicknessDK'),'beta'].abs().max())
        cortical_cax=ax.inset_axes([.90,.25,.022,.50])
        cortical_cb=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(-cortex_limit,cortex_limit),cmap='RdBu_r'),cax=cortical_cax)
        cortical_cb.set_ticks([-cortex_limit,cortex_limit])
        cortical_cb.set_ticklabels([f'{-cortex_limit:.4f}',f'{cortex_limit:.4f}'])
        cortical_cb.ax.tick_params(labelsize=8,length=0,pad=2)
        cortical_cb.outline.set_visible(False)
        fig.text(.025,.53,f'b  White-matter FA  |  {counts["FA27"]}/27 tracts with FDR < 0.05',fontsize=13,weight='bold')
        ax=brain_axes(plt.imread(out/'cache/tract_renders/FA27.png'),.24,.345,'b')
        selected=hits[hits.family.isin(['FA27','MD27'])]
        limit=max(.02,float(selected.beta.abs().max())) if len(selected) else .02
        mixed=(selected.beta<0).any()
        norm=TwoSlopeNorm(vmin=-limit,vcenter=0,vmax=limit) if mixed else Normalize(0,limit)
        cax=fig.add_axes([.16,.211,.16,.010])
        cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='RdBu_r' if mixed else 'Reds'),cax=cax,orientation='horizontal')
        cb.set_label('Standardized β',fontsize=10)
        cb.set_ticks([-limit,0,limit] if mixed else [0,limit/2,limit])
        fig.text(.025,.151,f'MD: {counts["MD27"]}/27 tracts with FDR < 0.05.',fontsize=10)
        fig.text(.025,.131,'Standard atlas surfaces; not individual tractography.',fontsize=9,color='.35')

        render_focus = mri_focus().render_focus
        rgb,regions,colour_limit=render_focus(d,mri_atlas_dir(),out)
        fig.text(.49,.53,'c  Hippocampus and thalamus volumes',fontsize=13,weight='bold')
        ax=brain_axes(rgb,.735,.345,'c')
        for v in regions:
            ax.annotate(f"{v['label']}\nβ = {v['beta']:.4f}{v['stars']}",xy=(v['x'],v['y']),
                xytext=(v['text_x'],v['text_y']),ha='center',va='center',fontsize=9,
                arrowprops=dict(arrowstyle='-',color='.3',lw=.8),bbox=dict(facecolor='white',edgecolor='none',pad=1))
        fig.text(.50,.185,'Posterior-superior view; L/R denote anatomical sides.',fontsize=9,color='.35')
        fig.text(.50,.163,'Darker red: larger positive β. * q<0.05; ** q<0.01; *** q<0.001.',fontsize=9)
        fig.text(.50,.141,'Selected regions highlighted; complete volume results in supplements.',fontsize=9,color='.35')
        wmh=d[d.family.eq('WMH')].iloc[0]
        fig.text(.025,.103,'d  White-matter hyperintensities',fontsize=12,weight='bold')
        fig.text(.025,.066,'Total WMH volume',fontsize=11,va='center')
        ax=fig.add_axes([.145,.046,.18,.041])
        ax.errorbar(wmh.beta,0,xerr=[[wmh.beta-wmh.LCI],[wmh.HCI-wmh.beta]],fmt='o',
            color='#2166ac',ms=5,capsize=3,lw=1.5)
        ax.axvline(0,color='.65',lw=.9,ls='--');ax.set_ylim(-1,1);ax.set_yticks([])
        ax.spines['left'].set_visible(False)
        from matplotlib.ticker import MaxNLocator
        lo,hi=min(0.,float(wmh.LCI)),max(0.,float(wmh.HCI))
        pad=max((hi-lo)*.12,1e-6)
        ax.set_xlim(lo-pad,hi+pad);ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.tick_params(labelsize=10)
        pd.DataFrame([dict(xmin=lo-pad,xmax=hi+pad,beta=wmh.beta,LCI=wmh.LCI,HCI=wmh.HCI)]).to_csv(
            out/'Figure4_WMH_axis_audit.csv',index=False)
        ax.set_xlabel('Standardized β (95% CI)',fontsize=10,labelpad=2)
        fig.text(.35,.066,f'β = {wmh.beta:.5f}\n95% CI {wmh.LCI:.5f} to {wmh.HCI:.5f}\nP = {wmh.P:.6f}',fontsize=9.5,va='center',linespacing=1.35)
        if font_scale != 1.0:
            from matplotlib.text import Text
            fig.canvas.draw()
            for text in fig.findobj(match=Text):
                text.set_fontsize(text.get_fontsize()*font_scale)
        if return_figure:
            return fig
        pd.DataFrame(size_audit).to_csv(out/'Figure4_brain_display_sizes.csv',index=False)
        d[d.family.eq('WMH')].to_csv(out/'Figure4_WMH_point_source.csv',index=False)
        for ext in ('png','pdf'):
            save_png_figure(fig, out/('Figure4_MRI.'+ext),dpi=300,facecolor='white',bbox_inches='tight')
        from matplotlib.transforms import Bbox
        for name,(x0,y0,x1,y1) in {
            'c_subcortical_focus':(.485,.125,.99,.552),
            'd_WMH_point_CI':(.02,.002,.52,.119)}.items():
            box=Bbox.from_extents(x0*17,y0*10,x1*17,y1*10)
            for ext in ('png','pdf'):save_png_figure(fig, out/(name+'.'+ext),dpi=300,facecolor='white',bbox_inches=box)
        plt.close(fig)
        shown=d[['field_id','family','category','label','beta','FDR_family']].copy()
        shown['panel_a']=d.family.eq('ThicknessDK') & d.FDR_family.lt(.05)
        shown['panel_b']=d.family.eq('FA27') & d.FDR_family.lt(.05)
        shown['panel_c']=d.field_id.isin([v['field_id'] for v in regions])
        shown['panel_d']=d.family.eq('WMH')
        shown.to_csv(out/'Figure4_MRI_display_plan.csv',index=False)
        (out/'Figure4_MRI_caption.txt').write_text(
            'DOSE and brain structural phenotypes. a, significant cortical thickness associations among 66 DK regions; '
            'selected anatomical labels shown, bilateral temporal poles unmeasured. b, significant white-matter FA tracts mapped '
            'to Oxford/UK Biobank group-average surfaces (CC BY 4.0); MD results explicitly reported. '
            'c, bilateral hippocampus and thalamus highlighted by user selection; all other volume results fully reported in supplements. '
            'Harvard-Oxford reference masks are not individual FIRST boundaries. Posterior-superior view; labels indicate anatomical sides. '
            'Darker red denotes larger positive beta; stars use BH across 139 volumes. '
            f'd, total WMH volume: standardized beta and 95% CI from the original adjusted model (N={int(wmh.N):,}); dashed line denotes beta=0. '
            'Red denotes positive association and blue negative association. Displayed brain outlines in b/c are scaled to the median width of an a subpanel; this is layout normalization, not a common anatomical millimetre scale. '
            'Both exposure and outcome were standardized within each complete-case model. '
            'Models adjusted for M2 covariates, baseline age and age squared, baseline centre; volume and WMH models additionally adjusted for ICV. '
            'Alcohol status not separately adjusted because DOSE includes Wine. '
            'BH performed separately within each score and prespecified MRI family, not within display categories. '
            'MRI sample size varies by outcome; see full results table. Observational associations, not MRI mediation or longitudinal atrophy. '
            'Panel letters a-d identify the MRI-only composite and can be relabelled when combined with clinical outcomes.',encoding='utf-8')
    return SimpleNamespace(**locals(), __file__=__file__)


if __name__ == "__main__":
    run_cli({"family_history": prepare_family_history})
