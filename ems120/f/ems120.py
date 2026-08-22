#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Imports, paths, and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EMS120 FINAL Python helpers: unified TEST-generation phone/geo/dx + stratified OOF CV.
# The former f/ml_phone.py, f/ml_geo.py, and f/ml_dx.py wrappers are no longer needed.
import os, json, time, re, math, gc, shutil
from pathlib import Path
from collections import Counter, defaultdict
from contextlib import nullcontext

import pandas as pd
import numpy as np
import joblib
import torch
try:
    from transformers import (
        AutoTokenizer,
        AutoModelForSequenceClassification,
        Trainer,
        TrainingArguments,
        set_seed,
    )
except ImportError:  # lets keyword/phone utilities be inspected without the AI environment
    AutoTokenizer = AutoModelForSequenceClassification = Trainer = TrainingArguments = set_seed = None
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import ElasticNetCV
from sklearn.metrics import (
    mean_squared_error,
    accuracy_score,
    f1_score,
    precision_recall_fscore_support,
    confusion_matrix,
)
from sklearn.model_selection import KFold, StratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Paths and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DIR0 = Path("D:/" if os.name == "nt" else "/mnt/d")
ANALYSIS_DIR = Path(os.environ.get(
    "EMS120_ANALYSIS_DIR",
    "D:/analysis/ems120" if os.name == "nt" else "/mnt/d/analysis/ems120",
))


def path0(*parts):
    return DIR0.joinpath(*parts).as_posix()


def display_path(path):
    text = os.fspath(path).replace("\\", "/")
    if os.name != "nt":
        return text
    m = re.match(r"^/mnt/([A-Za-z])(?:/(.*))?$", text)
    if not m:
        m = re.match(r"^/([A-Za-z])(?:/(.*))?$", text)
    if m:
        drive = m.group(1).upper()
        rest = (m.group(2) or "").replace("/", "\\")
        return f"{drive}:\\{rest}" if rest else f"{drive}:\\"

    if re.match(r"^[A-Za-z]:(?:/|$)", text):
        text = text.replace("/", "\\")
        return text if text.endswith("\\") or len(text) > 2 else text + "\\"
    return text


BASE_MODEL_DIR = path0("data.BIG", "AI_model_files", "hfl")
OUTPUT_MODEL_DIR = ANALYSIS_DIR.joinpath("ml", "hfl").as_posix()
TRAIN_FILE = ANALYSIS_DIR.joinpath("dat", "2019.train_dx.xlsx").as_posix()

TRAINED_MODEL_FILE = Path(OUTPUT_MODEL_DIR, "model.safetensors").as_posix()
KEYWORD_MODEL_FILE = Path(OUTPUT_MODEL_DIR, "keyword_model.json").as_posix()
LABELS_FILE = Path(OUTPUT_MODEL_DIR, "labels.json").as_posix()
MODEL_META_FILE = Path(OUTPUT_MODEL_DIR, "model_meta.json").as_posix()
DX_MODEL_VERSION = "macbert_dx_v2_fieldtag_cv5_kw4_20260809"
DX_CV_DIR = ANALYSIS_DIR.joinpath("raw", "hfl_cv5").as_posix()
DX_CV_MANIFEST = Path(DX_CV_DIR, "cv_manifest.json").as_posix()
DX_CV_SUMMARY = Path(DX_CV_DIR, "cv_summary.csv").as_posix()
DX_CV_SEED = 12010
DX_CV_FOLDS = 5
DX_CV_BOOTSTRAP_N = 1000

DX_TEXT_COLS = [
    "性别",
    "年龄",
    "呼救原因",
    "病因",
    "伤病程度",
    "症状",
    "主诉",
    "病史",
    "初步诊断",
    "补充诊断",
]
LABEL_COL = "疾病类型"
GEO_TEXT_COL = "地址"
GEO_LABEL_COL = "地址类型"
GEO_OUTPUT_MODEL_DIR = ANALYSIS_DIR.joinpath("ml", "hfl_geo").as_posix()
GEO_TRAINED_MODEL_FILE = Path(GEO_OUTPUT_MODEL_DIR, "model.safetensors").as_posix()
GEO_LABELS_FILE = Path(GEO_OUTPUT_MODEL_DIR, "labels.json").as_posix()
PHONE_OUTPUT_MODEL_DIR = ANALYSIS_DIR.joinpath("ml", "phone_sco3").as_posix()
PHONE_SCO3_MODEL_FILE = Path(PHONE_OUTPUT_MODEL_DIR, "phone_sco3_best_model.pkl").as_posix()
PHONE_SCO3_FEATURE_FILE = Path(PHONE_OUTPUT_MODEL_DIR, "phone_sco3_feature_names.pkl").as_posix()
PHONE_SCO3_REPORT_FILE = Path(PHONE_OUTPUT_MODEL_DIR, "phone_sco3_model_report.csv").as_posix()
SCRIPT_DIR = Path(__file__).resolve().parent

TOKENIZER = None
MODEL = None
LABEL2ID = None
ID2LABEL = None
KEYWORD_MODEL = None  # transparent field-aware weighted keyword comparator

GEO_TOKENIZER = None
GEO_MODEL = None
GEO_LABEL2ID = None
GEO_ID2LABEL = None


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Logging
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def ems_log(msg):
    print(f"[EMS120 {time.strftime('%H:%M:%S')}] {msg}", flush=True)

def gpu_mem():
    if torch is not None and torch.cuda.is_available():
        ems_log(f"GPU memory allocated: {torch.cuda.memory_allocated()/1024**3:.2f} GB")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Text preparation and transparent keyword comparator
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DX_TEXT_COLS = [
    "性别", "年龄", "呼救原因", "病因", "伤病程度", "症状", "主诉", "病史", "初步诊断", "补充诊断"
]
KEYWORD_FIELDS = ["呼救原因", "症状", "初步诊断", "补充诊断"]
KEYWORD_MODEL_VERSION = "field_kw_v4_logodds"

RAW_TO_GROUP = {
    "创伤-暴力事件": "Violence",
    "创伤-跌倒": "Fall",
    "创伤-交通事故": "Traffic",
    "创伤-其他原因": "Trauma",
    "创伤-高处坠落": "Trauma",
    "呼吸系统疾病": "Respiratory",
    "精神病": "Psychiatric",
    "理化中毒": "Intoxication",
    "泌尿系统疾病": "Urinary",
    "内分泌系统疾病": "Endocrine",
    "其他-死亡": "Death",
    "神经系统疾病-脑卒中": "Neurological",
    "神经系统疾病-其他疾病": "Neurological",
    "消化系统疾病": "Digestive",
    "心血管系统疾病-其他疾病": "CVD",
    "心血管系统疾病-胸痛": "CVD",
    "其他-胸闷": "CVD",
    "妇产科": "Ob/Gyn",
    "儿科": "Pediatric",
    "其他-其他症状": "Other",
    "其他-昏迷": "Other",
}
ANALYTIC_GROUP_ORDER = [
    "Violence", "Fall", "Traffic", "Trauma", "Respiratory", "Psychiatric",
    "Intoxication", "Urinary", "Endocrine", "Death", "Neurological",
    "Digestive", "CVD", "Ob/Gyn", "Pediatric", "Other"
]


def _require_transformers():
    if AutoTokenizer is None or AutoModelForSequenceClassification is None or Trainer is None:
        raise ImportError(
            "transformers is required for MacBERT training/inference. Activate the configured AI conda environment first."
        )


def _clean_fragment(x):
    if x is None or pd.isna(x):
        return ""
    s = str(x).strip()
    if s.lower() in {"", "nan", "none", "null"}:
        return ""
    return re.sub(r"\s+", " ", s)


def build_texts_from_df(df: pd.DataFrame):
    """Create field-tagged Chinese dispatch text.

    Field tags preserve the meaning of short fragments such as “中毒”, “胸痛、胸闷”
    and “损伤” without changing the underlying information available to the model.
    """
    cols = [c for c in DX_TEXT_COLS if c in df.columns]
    ems_log(f"Using DX text columns: {cols}")
    if not cols:
        return [""] * len(df)
    out = []
    values = df[cols].itertuples(index=False, name=None)
    for row in values:
        parts = []
        for c, v in zip(cols, row):
            z = _clean_fragment(v)
            if z:
                parts.append(f"{c} {z}")
        out.append("；".join(parts))
    return out


def _normalize_keyword_value(x):
    s = _clean_fragment(x)
    if not s:
        return ""
    s = re.sub(r"[\t\r\n]+", " ", s)
    s = re.sub(r"\s+", "", s)
    return s[:80]


def keyword_features_from_values(values: dict):
    """Small, auditable feature set for the transparent keyword comparator.

    We intentionally avoid a second black-box classifier. Features are exact values or
    short diagnosis fragments from four structured/free-text dispatch fields.
    """
    feats = set()
    for field in KEYWORD_FIELDS:
        s = _normalize_keyword_value(values.get(field, ""))
        if not s:
            continue
        feats.add(f"{field}={s}")
        # Diagnosis fields often contain several semicolon/comma separated statements.
        if field in {"初步诊断", "补充诊断"}:
            for part in re.split(r"[，,；;、。/|]+", s):
                part = part.strip("-:：?？()（）[]【】 ")
                if 2 <= len(part) <= 30:
                    feats.add(f"{field}~{part}")
    return feats


def _record_text(values: dict):
    return " ".join(_normalize_keyword_value(values.get(c, "")) for c in DX_TEXT_COLS if c in values)


# High-specificity, clinically transparent rules. Generic words (e.g. “呼吸”, “精神”,
# “外伤”) are deliberately excluded because they caused over-classification in the
# previous keyword implementation. The learned field-aware layer handles ambiguous cases.
DX_RULE_VERSION = "kw_v4_high_specificity"
DETERMINISTIC_DX_RULES = [
    {"label": "其他-死亡", "priority": 5.0, "threshold": 8.0, "patterns": {
        "无生命体征": 12, "心跳呼吸停止": 12, "呼吸心跳停止": 12, "心跳停止": 10,
        "呼吸停止": 10, "无脉搏": 10, "尸体": 10, "死亡": 9, "猝死": 9,
        "瞳孔散大": 6, "心肺复苏": 6,
    }},
    {"label": "创伤-交通事故", "priority": 4.0, "threshold": 8.0, "patterns": {
        "交通事故": 11, "车祸": 11, "被车撞": 11, "撞车": 10, "车撞": 10,
        "追尾": 9, "翻车": 9, "交通伤": 11,
    }},
    {"label": "创伤-高处坠落", "priority": 4.0, "threshold": 8.0, "patterns": {
        "高处坠落": 12, "高空坠落": 12, "坠楼": 12, "高坠伤": 12,
    }},
    {"label": "创伤-暴力事件", "priority": 3.5, "threshold": 8.0, "patterns": {
        "斗殴": 11, "打架": 11, "被人打": 10, "殴打": 10, "刀伤": 10,
        "砍伤": 10, "枪伤": 12,
    }},
    {"label": "理化中毒", "priority": 3.5, "threshold": 8.0, "patterns": {
        "一氧化碳": 12, "煤气中毒": 12, "有机磷": 12, "农药中毒": 11,
        "药物中毒": 11, "药物过量": 11, "过量服药": 11, "安眠药": 9,
        "吸毒": 10, "毒品": 10, "食物中毒": 10, "酒精中毒": 11,
        "急性酒精中毒": 12, "醉酒": 9,
    }},
    {"label": "精神病", "priority": 3.0, "threshold": 8.0, "patterns": {
        "精神异常": 11, "精神分裂": 12, "精神病": 10, "躁狂": 9, "幻觉": 9,
        "妄想": 9, "情绪失控": 8, "自杀": 9, "自伤": 9, "割腕": 10,
    }},
    {"label": "神经系统疾病-脑卒中", "priority": 3.0, "threshold": 8.0, "patterns": {
        "脑卒中": 12, "中风": 12, "偏瘫": 10, "失语": 10, "口角歪": 10,
        "嘴歪": 10, "脑血管意外": 11,
    }},
    {"label": "心血管系统疾病-胸痛", "priority": 2.5, "threshold": 8.0, "patterns": {
        "心肌梗死": 12, "心梗": 12, "心绞痛": 10, "急性冠脉": 10,
        "胸痛查因": 9, "胸痛待查": 9,
    }},
    {"label": "呼吸系统疾病", "priority": 2.0, "threshold": 8.0, "patterns": {
        "哮喘": 10, "肺炎": 9, "窒息": 10, "气喘": 8, "喘息": 8,
    }},
    {"label": "泌尿系统疾病", "priority": 2.0, "threshold": 8.0, "patterns": {
        "尿潴留": 10, "肾绞痛": 10, "血尿": 8, "尿痛": 8,
    }},
    {"label": "内分泌系统疾病", "priority": 2.0, "threshold": 8.0, "patterns": {
        "低血糖": 12, "高血糖": 10, "糖尿病酮症": 12,
    }},
    {"label": "妇产科", "priority": 2.0, "threshold": 8.0, "patterns": {
        "临产": 12, "分娩": 10, "宫缩": 9,
    }},
    {"label": "儿科", "priority": 1.5, "threshold": 8.0, "patterns": {
        "高热惊厥": 12, "小儿惊厥": 11,
    }},
]


def deterministic_keyword_predict(text: str, max_reason: int = 6):
    if text is None or not str(text).strip():
        return None, None
    s = str(text)
    scored = []
    for rule in DETERMINISTIC_DX_RULES:
        hits = [(pat, float(w)) for pat, w in rule["patterns"].items() if pat in s]
        if not hits:
            continue
        score = sum(w for _, w in hits) + float(rule.get("priority", 0.0))
        score += max(len(pat) for pat, _ in hits) / 10.0
        scored.append((score, float(rule.get("priority", 0.0)), rule, hits))
    if not scored:
        return None, None
    scored.sort(key=lambda z: (z[0], z[1]), reverse=True)
    best_score, _, best_rule, best_hits = scored[0]
    second_score = scored[1][0] if len(scored) > 1 else -np.inf
    threshold = float(best_rule.get("threshold", 8.0))
    # Ambiguous near-ties fall through to the learned field-aware layer.
    if best_score < threshold or (np.isfinite(second_score) and best_score - second_score < 1.0):
        return None, None
    best_hits = sorted(best_hits, key=lambda x: x[1], reverse=True)
    reason = ",".join(p for p, _ in best_hits[:max_reason])
    return best_rule["label"], f"RuleV4:{reason};score={best_score:.1f}"


def learn_keyword_model(df: pd.DataFrame, topk: int = 140, min_support: int = 2, alpha: float = 0.5):
    """Learn transparent field-value weights using smoothed class-vs-rest log odds."""
    if LABEL_COL not in df.columns:
        raise ValueError(f"training data missing label column: {LABEL_COL}")
    labels = df[LABEL_COL].astype(str)
    label_counts = labels.value_counts().to_dict()
    N = len(df)
    feat_global = Counter()
    feat_by_label = defaultdict(Counter)
    cols = [c for c in KEYWORD_FIELDS if c in df.columns]
    for row, lab in zip(df[cols].itertuples(index=False, name=None), labels):
        values = dict(zip(cols, row))
        feats = keyword_features_from_values(values)
        feat_global.update(feats)
        feat_by_label[str(lab)].update(feats)

    weights = {}
    for lab, n_lab in label_counts.items():
        candidates = {}
        for feat, n11 in feat_by_label[lab].items():
            fg = feat_global[feat]
            if fg < min_support:
                continue
            n10 = n_lab - n11
            n01 = fg - n11
            n00 = (N - n_lab) - n01
            # Smoothed log odds of feature presence in class vs all other classes.
            log_odds = math.log((n11 + alpha) / (n10 + alpha)) - math.log((n01 + alpha) / (n00 + alpha))
            if log_odds <= 0.15:
                continue
            support_boost = min(2.0, 0.75 + math.log1p(n11) / 3.0)
            candidates[feat] = float(log_odds * support_boost)
        top = sorted(candidates.items(), key=lambda x: x[1], reverse=True)[:topk]
        weights[lab] = {k: round(v, 6) for k, v in top}
    return {
        "version": KEYWORD_MODEL_VERSION,
        "fields": cols,
        "label_counts": {k: int(v) for k, v in label_counts.items()},
        "weights": weights,
    }


def keyword_predict_values(values: dict, model=None, max_reason: int = 6):
    model = KEYWORD_MODEL if model is None else model
    text = _record_text(values)
    lab_rule, reason_rule = deterministic_keyword_predict(text, max_reason=max_reason)
    if lab_rule is not None:
        return lab_rule, reason_rule
    if not model or "weights" not in model:
        return "其他-其他症状", "NoKeywordModel"
    feats = keyword_features_from_values(values)
    if not feats:
        return "其他-其他症状", "EmptyKeywordFields"
    best_lab, best_score, best_hits = "其他-其他症状", 0.0, []
    for lab, tok_w in model["weights"].items():
        hits = [(feat, float(tok_w[feat])) for feat in feats if feat in tok_w]
        if not hits:
            continue
        hits = sorted(hits, key=lambda x: x[1], reverse=True)
        score = sum(w for _, w in hits[:4])
        if score > best_score:
            best_lab, best_score, best_hits = lab, score, hits
    if best_score <= 0:
        return "其他-其他症状", "NoMatch"
    reason = ",".join(feat for feat, _ in best_hits[:max_reason])
    return best_lab, f"FieldKW:{reason};score={best_score:.2f}"


def keyword_predict_df(df: pd.DataFrame, model=None):
    cols = [c for c in set(DX_TEXT_COLS + KEYWORD_FIELDS) if c in df.columns]
    out_lab, out_reason = [], []
    for row in df[cols].itertuples(index=False, name=None):
        values = dict(zip(cols, row))
        lab, reason = keyword_predict_values(values, model=model)
        out_lab.append(lab)
        out_reason.append(reason)
    return out_lab, out_reason


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Dataset helpers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class DxDataset(torch.utils.data.Dataset):
    def __init__(self, enc, lab):
        self.enc = enc
        self.lab = lab

    def __getitem__(self, i):
        item = {k: torch.tensor(v[i]) for k, v in self.enc.items()}
        item["labels"] = torch.tensor(self.lab[i])
        return item

    def __len__(self):
        return len(self.enc["input_ids"])

def _build_label_maps_from_trainfile():
    df = pd.read_excel(TRAIN_FILE)
    if LABEL_COL not in df.columns:
        raise ValueError(f"training data missing label column: {LABEL_COL}")
    uniq = sorted(df[LABEL_COL].astype(str).unique())
    l2i = {l: i for i, l in enumerate(uniq)}
    i2l = {i: l for l, i in l2i.items()}
    return l2i, i2l


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 MacBERT 5-fold validation, final training, and model loading
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def _validate_training_df(df: pd.DataFrame, require_cv=False):
    missing = [c for c in [LABEL_COL] + DX_TEXT_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"training data missing required columns: {missing}")
    if df[LABEL_COL].isna().any():
        raise ValueError("training data contain missing disease labels")
    labels = df[LABEL_COL].astype(str).str.strip()
    if (labels == "").any():
        raise ValueError("training data contain blank disease labels")
    if require_cv:
        counts = labels.value_counts()
        if counts.min() < DX_CV_FOLDS:
            raise ValueError(
                f"stratified {DX_CV_FOLDS}-fold CV requires at least {DX_CV_FOLDS} records per class; "
                f"minimum class count is {counts.min()}"
            )
    return df


def _group_labels(labels):
    return [RAW_TO_GROUP.get(str(x), "Other") for x in labels]


def _metric_row(y_true, y_pred, system, level, fold=None):
    return {
        "system": system,
        "level": level,
        "fold": fold,
        "n": len(y_true),
        "accuracy": accuracy_score(y_true, y_pred),
        "macro_f1": f1_score(y_true, y_pred, average="macro", zero_division=0),
        "weighted_f1": f1_score(y_true, y_pred, average="weighted", zero_division=0),
    }


def _per_class_table(y_true, y_pred, label_order, system, level):
    pr, rc, f1, sup = precision_recall_fscore_support(
        y_true, y_pred, labels=label_order, zero_division=0
    )
    return pd.DataFrame({
        "system": system,
        "level": level,
        "class": label_order,
        "precision": pr,
        "recall": rc,
        "f1": f1,
        "support": sup.astype(int),
    })


def _bootstrap_metric_table(y_true, y_pred, system, level, n_boot=DX_CV_BOOTSTRAP_N, seed=DX_CV_SEED):
    """Bootstrap 95% CIs from the complete out-of-fold prediction set."""
    y_true = np.asarray(y_true, dtype=object)
    y_pred = np.asarray(y_pred, dtype=object)
    rng = np.random.default_rng(int(seed))
    n = len(y_true)
    vals = {"accuracy": [], "macro_f1": [], "weighted_f1": []}
    for _ in range(int(n_boot)):
        idx = rng.integers(0, n, size=n)
        yt, yp = y_true[idx], y_pred[idx]
        vals["accuracy"].append(accuracy_score(yt, yp))
        vals["macro_f1"].append(f1_score(yt, yp, average="macro", zero_division=0))
        vals["weighted_f1"].append(f1_score(yt, yp, average="weighted", zero_division=0))
    point = _metric_row(y_true, y_pred, system, level)
    rows = []
    for metric in ("accuracy", "macro_f1", "weighted_f1"):
        a = np.asarray(vals[metric], dtype=float)
        rows.append({
            "system": system,
            "level": level,
            "metric": metric,
            "estimate": float(point[metric]),
            "lo": float(np.quantile(a, 0.025)),
            "hi": float(np.quantile(a, 0.975)),
            "n_boot": int(n_boot),
            "n": int(n),
        })
    return pd.DataFrame(rows)


def _confidence_calibration_table(y_true, y_pred, confidence, level="grouped", n_bins=10):
    """Summarize OOF max-softmax confidence against observed correctness."""
    y_true = np.asarray(y_true, dtype=object)
    y_pred = np.asarray(y_pred, dtype=object)
    confidence = np.asarray(confidence, dtype=float)
    ok = np.isfinite(confidence)
    if not np.any(ok):
        return pd.DataFrame()
    d = pd.DataFrame({"true": y_true[ok], "pred": y_pred[ok], "confidence": confidence[ok]})
    d["correct"] = (d["true"] == d["pred"]).astype(int)
    d["confidence_bin"] = pd.qcut(
        d["confidence"].rank(method="first"), q=int(n_bins), labels=False
    ) + 1
    out = d.groupby("confidence_bin", as_index=False).agg(
        n=("correct", "size"),
        mean_confidence=("confidence", "mean"),
        observed_accuracy=("correct", "mean"),
        confidence_min=("confidence", "min"),
        confidence_max=("confidence", "max"),
    )
    out["level"] = level
    return out


def _predict_texts(model, tokenizer, texts, batch_size=256):
    model.eval()
    preds_all, conf_all = [], []
    use_cuda = DEVICE == "cuda"
    for i in range(0, len(texts), int(batch_size)):
        batch = texts[i:i + int(batch_size)]
        enc = tokenizer(batch, truncation=True, padding=True, max_length=128, return_tensors="pt")
        enc = {k: v.to(DEVICE) for k, v in enc.items()}
        amp_ctx = torch.amp.autocast("cuda", dtype=torch.float16) if use_cuda else nullcontext()
        with torch.inference_mode(), amp_ctx:
            logits = model(**enc).logits
            probs = torch.softmax(logits, dim=1)
        preds_all.extend(torch.argmax(probs, dim=1).cpu().numpy().tolist())
        conf_all.extend(probs.max(dim=1).values.cpu().numpy().tolist())
        del enc, logits, probs
    if use_cuda:
        torch.cuda.empty_cache()
    return preds_all, conf_all


def _cv_outputs_current():
    required = [
        DX_CV_MANIFEST, DX_CV_SUMMARY,
        Path(DX_CV_DIR, "cv_fold_metrics.csv").as_posix(),
        Path(DX_CV_DIR, "cv_per_class_raw.csv").as_posix(),
        Path(DX_CV_DIR, "cv_per_class_grouped.csv").as_posix(),
        Path(DX_CV_DIR, "cv_confusion_raw.csv").as_posix(),
        Path(DX_CV_DIR, "cv_confusion_grouped.csv").as_posix(),
        Path(DX_CV_DIR, "cv_oof_predictions.csv").as_posix(),
        Path(DX_CV_DIR, "cv_bootstrap_metrics.csv").as_posix(),
        Path(DX_CV_DIR, "cv_confidence_calibration.csv").as_posix(),
    ]
    if not all(os.path.exists(x) for x in required):
        return False
    try:
        meta = json.loads(Path(DX_CV_MANIFEST).read_text(encoding="utf-8"))
        return meta.get("dx_model_version") == DX_MODEL_VERSION and int(meta.get("folds", 0)) == DX_CV_FOLDS
    except Exception:
        return False


def refresh_dx_cv_derived_outputs():
    """Create bootstrap CIs and confidence-calibration tables from existing OOF predictions.

    This is intentionally a post-processing step: it does not retrain MacBERT.  It is
    useful when an earlier valid 5-fold OOF run already exists and only reviewer-ready
    uncertainty/calibration summaries are missing.
    """
    oof_file = Path(DX_CV_DIR, "cv_oof_predictions.csv")
    if not oof_file.exists():
        return {"status": "missing_oof", "dir": DX_CV_DIR}
    oof = pd.read_csv(oof_file)
    required = [LABEL_COL, "macbert_pred_raw", "macbert_confidence", "keyword_pred_raw"]
    missing = [c for c in required if c not in oof.columns]
    if missing:
        raise ValueError(f"Existing OOF predictions are missing columns: {missing}")

    y_true = oof[LABEL_COL].astype(str).to_numpy(dtype=object)
    y_pred = oof["macbert_pred_raw"].astype(str).to_numpy(dtype=object)
    y_kw = oof["keyword_pred_raw"].astype(str).to_numpy(dtype=object)
    confidence = pd.to_numeric(oof["macbert_confidence"], errors="coerce").to_numpy(dtype=float)
    y_true_g = np.asarray(_group_labels(y_true), dtype=object)
    y_pred_g = np.asarray(_group_labels(y_pred), dtype=object)
    y_kw_g = np.asarray(_group_labels(y_kw), dtype=object)

    bootstrap = pd.concat([
        _bootstrap_metric_table(y_true, y_pred, "MacBERT", "raw", seed=DX_CV_SEED + 11),
        _bootstrap_metric_table(y_true_g, y_pred_g, "MacBERT", "grouped", seed=DX_CV_SEED + 12),
        _bootstrap_metric_table(y_true, y_kw, "Keyword", "raw", seed=DX_CV_SEED + 13),
        _bootstrap_metric_table(y_true_g, y_kw_g, "Keyword", "grouped", seed=DX_CV_SEED + 14),
    ], ignore_index=True)
    bootstrap.to_csv(Path(DX_CV_DIR, "cv_bootstrap_metrics.csv"), index=False, encoding="utf-8-sig")

    calibration = _confidence_calibration_table(
        y_true_g, y_pred_g, confidence, level="grouped", n_bins=10
    )
    calibration.to_csv(Path(DX_CV_DIR, "cv_confidence_calibration.csv"), index=False, encoding="utf-8-sig")

    manifest_path = Path(DX_CV_MANIFEST)
    if manifest_path.exists():
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except Exception:
            manifest = {}
    else:
        manifest = {}
    manifest["dx_model_version"] = DX_MODEL_VERSION
    manifest["folds"] = DX_CV_FOLDS
    manifest["bootstrap_replicates"] = DX_CV_BOOTSTRAP_N
    manifest["derived_outputs_refreshed_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    ems_log("Refreshed CV bootstrap CIs and confidence calibration from existing OOF predictions")
    return {
        "status": "refreshed",
        "dir": DX_CV_DIR,
        "bootstrap": Path(DX_CV_DIR, "cv_bootstrap_metrics.csv").as_posix(),
        "calibration": Path(DX_CV_DIR, "cv_confidence_calibration.csv").as_posix(),
    }


def run_dx_stratified_cv(force=False):
    """Stratified 5-fold OOF validation on the manually annotated 2019 workbook.

    Outputs both the original fine-grained labels and the broader analytic groups used
    in the epidemiologic analyses. A transparent keyword comparator is trained within
    each fold, so its OOF performance is not calculated on records used to estimate
    its field-value weights.
    """
    _require_transformers()
    # Geo inference leaves global autograd disabled.  CV is a new training phase,
    # so restore it explicitly before creating any fold model.
    torch.set_grad_enabled(True)
    force = bool(force)
    if not force and _cv_outputs_current():
        ems_log(f"MacBERT CV outputs already current: {display_path(DX_CV_DIR)}")
        return {"status": "cached", "summary": DX_CV_SUMMARY, "dir": DX_CV_DIR}

    Path(DX_CV_DIR).mkdir(parents=True, exist_ok=True)
    df = _validate_training_df(pd.read_excel(TRAIN_FILE), require_cv=True).reset_index(drop=True)
    labels = df[LABEL_COL].astype(str).tolist()
    label_order = sorted(pd.Series(labels, dtype="string").unique().tolist())
    label2id = {lab: i for i, lab in enumerate(label_order)}
    id2label = {i: lab for lab, i in label2id.items()}
    y_id = np.array([label2id[x] for x in labels], dtype=int)
    texts = build_texts_from_df(df)

    tokenizer = AutoTokenizer.from_pretrained(BASE_MODEL_DIR)
    skf = StratifiedKFold(n_splits=DX_CV_FOLDS, shuffle=True, random_state=DX_CV_SEED)
    oof_pred = np.empty(len(df), dtype=object)
    oof_conf = np.full(len(df), np.nan, dtype=float)
    oof_kw = np.empty(len(df), dtype=object)
    oof_kw_reason = np.empty(len(df), dtype=object)
    fold_id = np.zeros(len(df), dtype=int)
    fold_metrics = []

    for fold, (tr_idx, va_idx) in enumerate(skf.split(np.zeros(len(df)), labels), start=1):
        ems_log(f"MacBERT CV fold {fold}/{DX_CV_FOLDS}: train={len(tr_idx)}, validation={len(va_idx)}")
        if set_seed is not None:
            set_seed(DX_CV_SEED + fold)
        fold_dir = Path(DX_CV_DIR, f"fold_{fold}_tmp")
        if fold_dir.exists():
            shutil.rmtree(fold_dir, ignore_errors=True)
        model = AutoModelForSequenceClassification.from_pretrained(
            BASE_MODEL_DIR, num_labels=len(label_order)
        ).to(DEVICE)
        model.train()
        model.requires_grad_(True)
        train_enc = tokenizer(
            [texts[i] for i in tr_idx], truncation=True, padding=True, max_length=128
        )
        train_y = [int(y_id[i]) for i in tr_idx]
        args = TrainingArguments(
            output_dir=fold_dir.as_posix(),
            per_device_train_batch_size=16,
            per_device_eval_batch_size=64,
            num_train_epochs=2,
            learning_rate=5e-5,
            weight_decay=0.01,
            logging_steps=100,
            logging_strategy="steps",
            disable_tqdm=False,
            save_strategy="no",
            fp16=torch.cuda.is_available(),
            seed=DX_CV_SEED + fold,
            data_seed=DX_CV_SEED + fold,
            report_to=[],
        )
        Trainer(model=model, args=args, train_dataset=DxDataset(train_enc, train_y)).train()
        pred_id, pred_conf = _predict_texts(model, tokenizer, [texts[i] for i in va_idx], batch_size=256)
        pred_lab = [id2label[int(i)] for i in pred_id]
        oof_pred[va_idx] = pred_lab
        oof_conf[va_idx] = pred_conf
        fold_id[va_idx] = fold

        kw_model = learn_keyword_model(df.iloc[tr_idx].reset_index(drop=True))
        kw_lab, kw_reason = keyword_predict_df(df.iloc[va_idx].reset_index(drop=True), model=kw_model)
        oof_kw[va_idx] = kw_lab
        oof_kw_reason[va_idx] = kw_reason

        y_val = [labels[i] for i in va_idx]
        fold_metrics.append(_metric_row(y_val, pred_lab, "MacBERT", "raw", fold))
        fold_metrics.append(_metric_row(_group_labels(y_val), _group_labels(pred_lab), "MacBERT", "grouped", fold))
        fold_metrics.append(_metric_row(y_val, kw_lab, "Keyword", "raw", fold))
        fold_metrics.append(_metric_row(_group_labels(y_val), _group_labels(kw_lab), "Keyword", "grouped", fold))

        del model, train_enc
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        shutil.rmtree(fold_dir, ignore_errors=True)

    y_true = np.array(labels, dtype=object)
    y_pred = np.array(oof_pred, dtype=object)
    y_kw = np.array(oof_kw, dtype=object)
    y_true_g = np.array(_group_labels(y_true), dtype=object)
    y_pred_g = np.array(_group_labels(y_pred), dtype=object)
    y_kw_g = np.array(_group_labels(y_kw), dtype=object)

    summary_rows = [
        _metric_row(y_true, y_pred, "MacBERT", "raw"),
        _metric_row(y_true_g, y_pred_g, "MacBERT", "grouped"),
        _metric_row(y_true, y_kw, "Keyword", "raw"),
        _metric_row(y_true_g, y_kw_g, "Keyword", "grouped"),
    ]
    summary = pd.DataFrame(summary_rows)
    summary["macbert_keyword_agreement"] = np.nan
    summary.loc[(summary.system == "MacBERT") & (summary.level == "raw"), "macbert_keyword_agreement"] = np.mean(y_pred == y_kw)
    summary.loc[(summary.system == "MacBERT") & (summary.level == "grouped"), "macbert_keyword_agreement"] = np.mean(y_pred_g == y_kw_g)
    summary.to_csv(DX_CV_SUMMARY, index=False, encoding="utf-8-sig")
    pd.DataFrame(fold_metrics).to_csv(Path(DX_CV_DIR, "cv_fold_metrics.csv"), index=False, encoding="utf-8-sig")

    bootstrap_metrics = pd.concat([
        _bootstrap_metric_table(y_true, y_pred, "MacBERT", "raw", seed=DX_CV_SEED + 11),
        _bootstrap_metric_table(y_true_g, y_pred_g, "MacBERT", "grouped", seed=DX_CV_SEED + 12),
        _bootstrap_metric_table(y_true, y_kw, "Keyword", "raw", seed=DX_CV_SEED + 13),
        _bootstrap_metric_table(y_true_g, y_kw_g, "Keyword", "grouped", seed=DX_CV_SEED + 14),
    ], ignore_index=True)
    bootstrap_metrics.to_csv(
        Path(DX_CV_DIR, "cv_bootstrap_metrics.csv"), index=False, encoding="utf-8-sig"
    )

    confidence_calibration = _confidence_calibration_table(
        y_true_g, y_pred_g, oof_conf, level="grouped", n_bins=10
    )
    confidence_calibration.to_csv(
        Path(DX_CV_DIR, "cv_confidence_calibration.csv"), index=False, encoding="utf-8-sig"
    )

    per_raw = pd.concat([
        _per_class_table(y_true, y_pred, label_order, "MacBERT", "raw"),
        _per_class_table(y_true, y_kw, label_order, "Keyword", "raw"),
    ], ignore_index=True)
    per_group = pd.concat([
        _per_class_table(y_true_g, y_pred_g, ANALYTIC_GROUP_ORDER, "MacBERT", "grouped"),
        _per_class_table(y_true_g, y_kw_g, ANALYTIC_GROUP_ORDER, "Keyword", "grouped"),
    ], ignore_index=True)
    per_raw.to_csv(Path(DX_CV_DIR, "cv_per_class_raw.csv"), index=False, encoding="utf-8-sig")
    per_group.to_csv(Path(DX_CV_DIR, "cv_per_class_grouped.csv"), index=False, encoding="utf-8-sig")

    cm_raw = pd.DataFrame(confusion_matrix(y_true, y_pred, labels=label_order), index=label_order, columns=label_order)
    cm_raw.index.name = "true_label"
    cm_raw.to_csv(Path(DX_CV_DIR, "cv_confusion_raw.csv"), encoding="utf-8-sig")
    cm_group = pd.DataFrame(confusion_matrix(y_true_g, y_pred_g, labels=ANALYTIC_GROUP_ORDER), index=ANALYTIC_GROUP_ORDER, columns=ANALYTIC_GROUP_ORDER)
    cm_group.index.name = "true_group"
    cm_group.to_csv(Path(DX_CV_DIR, "cv_confusion_grouped.csv"), encoding="utf-8-sig")

    oof = df[[c for c in ["ID", LABEL_COL] if c in df.columns]].copy()
    oof["fold"] = fold_id
    oof["macbert_pred_raw"] = y_pred
    oof["macbert_confidence"] = oof_conf
    oof["keyword_pred_raw"] = y_kw
    oof["keyword_reason"] = oof_kw_reason
    oof["true_group"] = y_true_g
    oof["macbert_pred_group"] = y_pred_g
    oof["keyword_pred_group"] = y_kw_g
    oof.to_csv(Path(DX_CV_DIR, "cv_oof_predictions.csv"), index=False, encoding="utf-8-sig")

    manifest = {
        "dx_model_version": DX_MODEL_VERSION,
        "keyword_model_version": KEYWORD_MODEL_VERSION,
        "folds": DX_CV_FOLDS,
        "seed": DX_CV_SEED,
        "bootstrap_replicates": DX_CV_BOOTSTRAP_N,
        "n": int(len(df)),
        "raw_classes": label_order,
        "raw_class_counts": {k: int(v) for k, v in pd.Series(labels).value_counts().to_dict().items()},
        "analytic_groups": ANALYTIC_GROUP_ORDER,
        "base_model_dir": display_path(BASE_MODEL_DIR),
        "training_file": display_path(TRAIN_FILE),
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    Path(DX_CV_MANIFEST).write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    ems_log(f"MacBERT CV finished; summary written to {display_path(DX_CV_SUMMARY)}")
    return {"status": "completed", "summary": DX_CV_SUMMARY, "dir": DX_CV_DIR}


def _dx_model_is_current():
    if not os.path.exists(TRAINED_MODEL_FILE) or not os.path.exists(MODEL_META_FILE):
        return False
    try:
        meta = json.loads(Path(MODEL_META_FILE).read_text(encoding="utf-8"))
        return meta.get("dx_model_version") == DX_MODEL_VERSION
    except Exception:
        return False


def train_dx_model():
    global TOKENIZER, MODEL, LABEL2ID, ID2LABEL, KEYWORD_MODEL
    _require_transformers()
    torch.set_grad_enabled(True)
    Path(OUTPUT_MODEL_DIR).mkdir(parents=True, exist_ok=True)
    ems_log(f"training final MacBERT disease model from {display_path(TRAIN_FILE)}")
    df = _validate_training_df(pd.read_excel(TRAIN_FILE), require_cv=False).reset_index(drop=True)
    texts = build_texts_from_df(df)
    labels = df[LABEL_COL].astype(str)
    uniq = sorted(labels.unique())
    LABEL2ID = {l: i for i, l in enumerate(uniq)}
    ID2LABEL = {i: l for l, i in LABEL2ID.items()}
    KEYWORD_MODEL = learn_keyword_model(df)
    y = labels.map(LABEL2ID).tolist()

    if set_seed is not None:
        set_seed(DX_CV_SEED)
    TOKENIZER = AutoTokenizer.from_pretrained(BASE_MODEL_DIR)
    enc = TOKENIZER(texts, truncation=True, padding=True, max_length=128)
    MODEL = AutoModelForSequenceClassification.from_pretrained(
        BASE_MODEL_DIR, num_labels=len(uniq)
    ).to(DEVICE)
    MODEL.train()
    MODEL.requires_grad_(True)
    gpu_mem()
    args = TrainingArguments(
        output_dir=OUTPUT_MODEL_DIR,
        per_device_train_batch_size=16,
        num_train_epochs=2,
        learning_rate=5e-5,
        weight_decay=0.01,
        logging_steps=50,
        logging_strategy="steps",
        disable_tqdm=False,
        save_strategy="no",
        fp16=torch.cuda.is_available(),
        seed=DX_CV_SEED,
        data_seed=DX_CV_SEED,
        report_to=[],
    )
    Trainer(model=MODEL, args=args, train_dataset=DxDataset(enc, y)).train()
    MODEL.save_pretrained(OUTPUT_MODEL_DIR)
    TOKENIZER.save_pretrained(OUTPUT_MODEL_DIR)
    Path(KEYWORD_MODEL_FILE).write_text(json.dumps(KEYWORD_MODEL, ensure_ascii=False, indent=2), encoding="utf-8")
    Path(LABELS_FILE).write_text(json.dumps(
        {"label2id": LABEL2ID, "id2label": {str(k): v for k, v in ID2LABEL.items()}},
        ensure_ascii=False, indent=2
    ), encoding="utf-8")
    Path(MODEL_META_FILE).write_text(json.dumps({
        "dx_model_version": DX_MODEL_VERSION,
        "keyword_model_version": KEYWORD_MODEL_VERSION,
        "training_n": int(len(df)),
        "raw_classes": uniq,
        "base_model_dir": display_path(BASE_MODEL_DIR),
        "training_file": display_path(TRAIN_FILE),
        "trained_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    ems_log(f"final MacBERT disease model saved: {display_path(OUTPUT_MODEL_DIR)}")


def load_model_once():
    global TOKENIZER, MODEL, KEYWORD_MODEL, LABEL2ID, ID2LABEL
    _require_transformers()
    if MODEL is not None:
        return
    if not _dx_model_is_current():
        ems_log("existing disease model is absent or from an older pipeline version; retraining")
        train_dx_model()
    else:
        ems_log(f"current disease model found: {display_path(TRAINED_MODEL_FILE)}")

    TOKENIZER = AutoTokenizer.from_pretrained(OUTPUT_MODEL_DIR)
    load_dtype = torch.float16 if DEVICE == "cuda" else torch.float32
    MODEL = AutoModelForSequenceClassification.from_pretrained(
        OUTPUT_MODEL_DIR, dtype=load_dtype
    ).to(DEVICE)
    MODEL.eval()
    MODEL.requires_grad_(False)
    torch.set_grad_enabled(False)

    try:
        KEYWORD_MODEL = json.loads(Path(KEYWORD_MODEL_FILE).read_text(encoding="utf-8"))
        if KEYWORD_MODEL.get("version") != KEYWORD_MODEL_VERSION:
            raise ValueError("old keyword model version")
    except Exception:
        ems_log("keyword model absent/outdated; rebuilding from expert-reviewed workbook")
        df = _validate_training_df(pd.read_excel(TRAIN_FILE), require_cv=False)
        KEYWORD_MODEL = learn_keyword_model(df)
        Path(KEYWORD_MODEL_FILE).write_text(json.dumps(KEYWORD_MODEL, ensure_ascii=False, indent=2), encoding="utf-8")

    if os.path.exists(LABELS_FILE):
        obj = json.loads(Path(LABELS_FILE).read_text(encoding="utf-8"))
        LABEL2ID = obj["label2id"]
        ID2LABEL = {int(k): v for k, v in obj["id2label"].items()}
    else:
        LABEL2ID, ID2LABEL = _build_label_maps_from_trainfile()
        Path(LABELS_FILE).write_text(json.dumps(
            {"label2id": LABEL2ID, "id2label": {str(k): v for k, v in ID2LABEL.items()}},
            ensure_ascii=False, indent=2
        ), encoding="utf-8")
    gpu_mem()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Phone SCO3 feature model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Phone model files are controlled below in ems120.py
def normalize_phone_tail_ml(x):
    if pd.isna(x) or x is None:
        return ""
    s = str(x).strip()
    if s.endswith(".0"):
        s = s[:-2]
    s = "".join(ch for ch in s if ch.isdigit())
    if len(s) == 11:
        s = s[3:]
    return s

def longest_run(s, mode="same"):
    if not s:
        return 0
    best = cur = 1
    for i in range(1, len(s)):
        if mode == "same":
            ok = s[i] == s[i - 1]
        elif mode == "asc":
            ok = int(s[i]) - int(s[i - 1]) == 1
        else:
            ok = int(s[i]) - int(s[i - 1]) == -1
        cur = cur + 1 if ok else 1
        best = max(best, cur)
    return best

def digit_entropy(s):
    if not s:
        return 0.0
    counts = np.array([s.count(str(d)) for d in range(10)], dtype=float)
    p = counts[counts > 0] / len(s)
    return float(-(p * np.log(p)).sum() / math.log(10))

def make_phone_features(phone_numbers):
    tails = [normalize_phone_tail_ml(x) for x in phone_numbers]
    rows = []
    for s in tails:
        valid = len(s) == 8 and s.isdigit()
        z = s if valid else ""
        row = {"phone_tail": z, "valid_phone": int(valid)}
        for d in range(10):
            row[f"n{d}"] = z.count(str(d))
        row.update(
            has_4=int("4" in z),
            n8=z.count("8"),
            n6=z.count("6"),
            n9=z.count("9"),
            n1=z.count("1"),
            tail1=int(z[-1:]) if valid else -1,
            tail2=int(z[-2:]) if valid else -1,
            tail3=int(z[-3:]) if valid else -1,
            tail4=int(z[-4:]) if valid else -1,
            has_88=int("88" in z),
            has_888=int("888" in z),
            has_8888=int("8888" in z),
            has_66=int("66" in z),
            has_666=int("666" in z),
            has_99=int("99" in z),
            has_999=int("999" in z),
            has_168=int("168" in z),
            has_518=int("518" in z),
            has_618=int("618" in z),
            has_123=int("123" in z),
            has_1234=int("1234" in z),
            has_6789=int("6789" in z),
            has_1314=int("1314" in z),
            has_520=int("520" in z),
            has_14=int("14" in z),
            has_74=int("74" in z),
            has_44=int("44" in z),
            has_444=int("444" in z),
            ends_8=int(z.endswith("8")) if valid else 0,
            ends_88=int(z.endswith("88")) if valid else 0,
            ends_888=int(z.endswith("888")) if valid else 0,
            ends_8888=int(z.endswith("8888")) if valid else 0,
            ends_168=int(z.endswith("168")) if valid else 0,
            ends_518=int(z.endswith("518")) if valid else 0,
            ends_668=int(z.endswith("668")) if valid else 0,
            ends_688=int(z.endswith("688")) if valid else 0,
            ends_999=int(z.endswith("999")) if valid else 0,
            longest_same_run=longest_run(z, "same") if valid else 0,
            longest_ascending_run=longest_run(z, "asc") if valid else 0,
            longest_descending_run=longest_run(z, "desc") if valid else 0,
            palindrome_tail4=int(valid and z[-4:] == z[-4:][::-1] and len(set(z[-4:])) > 1),
            unique_digit_count=len(set(z)) if valid else 0,
            entropy=digit_entropy(z) if valid else 0.0,
        )
        rows.append(row)
    return pd.DataFrame(rows)

def make_reason(s):
    s = normalize_phone_tail_ml(s)
    reasons = []
    if not (len(s) == 8 and s.isdigit()):
        return ""
    if "4" in s:
        reasons.append("contains digit 4, which strongly lowers the score")
    if s.count("8") >= 3:
        reasons.append("contains multiple 8s")
    if s.endswith("8888"):
        reasons.append("ends with 8888")
    elif s.endswith("888"):
        reasons.append("ends with 888")
    elif s.endswith("88"):
        reasons.append("ends with 88")
    if "168" in s:
        reasons.append("contains 168")
    if "518" in s:
        reasons.append("contains 518")
    if "618" in s:
        reasons.append("contains 618")
    if "666" in s:
        reasons.append("contains 666")
    if "999" in s:
        reasons.append("contains 999")
    if "1234" in s:
        reasons.append("contains ascending pattern 1234")
    if len(reasons) == 0:
        reasons.append("no strong lucky or unlucky pattern detected")
    return "; ".join(reasons)

def _candidate_models(random_state=12010):
    models = {
        "elastic_net": Pipeline(
            [
                ("scale", StandardScaler()),
                ("model", ElasticNetCV(l1_ratio=[0.1, 0.5, 0.9], cv=5, random_state=random_state, max_iter=20000)),
            ]
        ),
        "random_forest": RandomForestRegressor(
            n_estimators=500,
            min_samples_leaf=5,
            random_state=random_state,
            n_jobs=-1,
        ),
    }
    try:
        from lightgbm import LGBMRegressor

        models["lightgbm"] = LGBMRegressor(
            n_estimators=600,
            learning_rate=0.03,
            num_leaves=31,
            random_state=random_state,
            n_jobs=-1,
            verbose=-1,
        )
    except Exception:
        pass
    try:
        from xgboost import XGBRegressor

        models["xgboost"] = XGBRegressor(
            n_estimators=600,
            learning_rate=0.03,
            max_depth=4,
            subsample=0.9,
            colsample_bytree=0.9,
            objective="reg:squarederror",
            random_state=random_state,
            n_jobs=-1,
        )
    except Exception:
        pass
    return models

def train_phone_sco3_ml(train_file, model_file, feature_file, report_file=None, min_rows=100):
    df = pd.read_excel(train_file)
    phone_col = "phone" if "phone" in df.columns else "电话" if "电话" in df.columns else None
    if phone_col is None or "phone.sco" not in df.columns:
        raise ValueError("training file must contain phone/电话 and phone.sco")

    y = pd.to_numeric(df["phone.sco"], errors="coerce")
    feat = make_phone_features(df[phone_col])
    ok = (feat["valid_phone"] == 1) & y.between(0, 10)
    feat = feat.loc[ok].reset_index(drop=True)
    y = y.loc[ok].astype(float).reset_index(drop=True)
    if len(feat) < min_rows:
        raise ValueError(f"not enough usable expert-reviewed phone rows: {len(feat)}")

    X = feat.drop(columns=["phone_tail"])
    feature_names = X.columns.tolist()
    cv = KFold(n_splits=5, shuffle=True, random_state=12010)
    rows = []
    best_name, best_model, best_rmse = None, None, float("inf")
    for name, model in _candidate_models().items():
        try:
            scores = cross_val_score(model, X, y, cv=cv, scoring="neg_mean_squared_error", n_jobs=1)
            rmse_scores = np.sqrt(-scores)
            rmse = float(rmse_scores.mean())
            rows.append({"model": name, "cv_rmse": rmse, "cv_rmse_sd": float(rmse_scores.std()), "error": ""})
        except Exception as e:
            rows.append({"model": name, "cv_rmse": np.inf, "cv_rmse_sd": np.nan, "error": str(e)})
            ems_log(f"phone sco3 candidate skipped ({name}): {e}")
            continue
        if rmse < best_rmse:
            best_name, best_model, best_rmse = name, model, rmse

    if best_model is None:
        raise ValueError("no phone sco3 candidate model could be trained")

    best_model.fit(X, y)
    pred = np.clip(best_model.predict(X), 0, 10)
    train_rmse = float(np.sqrt(mean_squared_error(y, pred)))

    model_path = Path(model_file)
    model_path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(best_model, model_file)
    joblib.dump(feature_names, feature_file)
    report = pd.DataFrame(rows).sort_values("cv_rmse")
    report["selected"] = report["model"] == best_name
    report["train_rmse_selected"] = train_rmse
    if report_file:
        report.to_csv(report_file, index=False, encoding="utf-8-sig")
    return {"best_model": best_name, "cv_rmse": best_rmse, "train_rmse": train_rmse, "n": len(X)}

def predict_phone_sco3_ml(phone_numbers, model_file, feature_file):
    model = joblib.load(model_file)
    feature_names = joblib.load(feature_file)
    feat = make_phone_features(phone_numbers)
    X = feat[feature_names]
    pred = np.clip(model.predict(X), 0, 10)
    reasons = [make_reason(s) for s in phone_numbers]
    return pred.astype(float).tolist(), reasons


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Phone model helpers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# phone.sco is a 0-10 score: base + digit_luck + pattern_luck, then clamped to [min,max]
# phone.sco.reason is a compact "calculation string" to show how the score is formed.
# Digit-based luck (simple cultural weights)
RULE_SIMPLE = {
    "base": 4.0,
    "min": 0.0,
    "max": 10.0,
    # per-digit weights: contribution = count(digit) * weight
    "w": {"8": 1.0, "9": 1.0, "1": 0.5, "6": 0.5, "4": -4.0},
}

# Pattern-based luck (structured aesthetics) with a hard cap to avoid score saturation at 10
RULE_PATTERN = {
    "run_min": 3,     # 777 / 8888 ... (same digit run length threshold)
    "seq_min": 3,     # 123 / 987 ... (monotonic +/-1 length threshold)

    # weights for pattern contributions
    "w_run": 0.8,        # run contribution: len * w_run
    "w_seq": 0.6,        # sequence contribution: len * w_seq
    "w_repeat4": 3.5,    # 56785678 (repeat 4-digit block): fixed add
    "w_repeat2": 2.0,    # 23232323 (repeat 2-digit block): fixed add
    "w_mirror": 2.5,     # 23455432 (mirror): fixed add

    # cap the total pattern contribution
    "cap_pat": 6.0,
}

def _is_digit8(p: str) -> bool:
    return len(p) == 8 and p.isdigit()

def _repeat_block_score(p: str, block_len: int):
    """Detect repeated blocks: e.g., 56785678 (block_len=4) or 23232323 (block_len=2)."""
    if 8 % block_len != 0:
        return False, ""
    k = 8 // block_len
    blk = p[:block_len]
    if blk * k == p and len(set(blk)) > 1:
        return True, f"{blk}x{k}"
    return False, ""

def _mirror_score(p: str):
    """Detect mirror: 23455432 (a|a[::-1])."""
    a, b = p[:4], p[4:]
    if b == a[::-1] and len(set(a)) > 1:
        return True, f"{a}|{b}"
    return False, ""

def _seq_matches(p: str, seq_min: int):
    """Return monotonic +/-1 sequences as (substring, length, step)."""
    out = []
    start = 0
    cur_step = None
    cur_len = 1

    def flush(s, l, st):
        if l >= seq_min and st in (1, -1):
            out.append((p[s:s+l], l, st))

    for i in range(1, 8):
        diff = int(p[i]) - int(p[i - 1])
        if abs(diff) == 1:
            if cur_step is None:
                cur_step = diff
                cur_len = 2
                start = i - 1
            elif diff == cur_step:
                cur_len += 1
            else:
                flush(start, cur_len, cur_step)
                cur_step = diff
                cur_len = 2
                start = i - 1
        else:
            flush(start, cur_len, cur_step)
            cur_step = None
            cur_len = 1

    flush(start, cur_len, cur_step)
    return out

def _fmt_term(name: str, count: int, w: float):
    # e.g., "8(2*1.0)" or "4(1*-4.0)"
    return f"{name}({count}*{w:g})"

def _normalize_phone_tail(p):
    if pd.isna(p) or p is None:
        return ""
    p = str(p).strip()
    if p.endswith(".0"):
        p = p[:-2]
    p = re.sub(r"\D", "", p)
    if len(p) == 11:
        p = p[3:]
    return p

def score_phone_sco1(p: str):
    if not _is_digit8(p):
        return float("nan"), ""
    if "4" in p:
        return 0.0, "contains digit 4 => 0"
    n8, n9, n6, n1 = p.count("8"), p.count("9"), p.count("6"), p.count("1")
    raw = n8 * 1.0 + (n9 + n6 + n1) * 0.5
    score = min(10.0, raw / 8.0 * 10.0)
    return round(score, 1), f"raw=8({n8}*1)+9({n9}*0.5)+6({n6}*0.5)+1({n1}*0.5)={raw:.1f}; rescaled to {score:.1f}/10"

def _has_any(p, patterns):
    return [z for z in patterns if z in p]

def score_phone_sco2(p: str):
    """Rule score with a rescue mechanism.

    Unlike sco1, digit 4 is a penalty rather than an automatic zero. Strong lucky
    structure, such as 12345678, 888, 168, repeated blocks, or mirror patterns,
    can rescue the score into a reasonable range.
    """
    if not _is_digit8(p):
        return float("nan"), ""

    base = 3.0
    raw = base
    terms = ["base=3"]
    reasons = []

    n8, n9, n6, n1, n4 = p.count("8"), p.count("9"), p.count("6"), p.count("1"), p.count("4")
    digit_score = n8 * 1.0 + n9 * 0.7 + n6 * 0.6 + n1 * 0.4 - n4 * 1.5
    raw += digit_score
    terms.append(f"digit={digit_score:.1f}")

    if n4:
        raw -= 1.0
        terms.append("any4=-1")
        reasons.append("contains digit 4, but not forced to zero")
    for bad, penalty in {"44": 1.5, "444": 3.0, "14": 0.8, "74": 0.8}.items():
        if bad in p:
            raw -= penalty
            terms.append(f"{bad}=-{penalty:g}")

    rescue = 0.0
    tail_bonus = [("8888", 4.5), ("888", 3.3), ("168", 2.1), ("518", 1.8), ("668", 1.8), ("688", 2.0), ("999", 2.1), ("88", 1.6), ("8", 0.6)]
    for pat, bonus in tail_bonus:
        if p.endswith(pat):
            rescue += bonus
            terms.append(f"tail{pat}=+{bonus:g}")
            reasons.append(f"ends with {pat}")
            break

    for pat, bonus in {"12345678": 3.0, "23456789": 2.6, "1234": 1.2, "6789": 1.2, "168": 1.5, "198": 1.1, "518": 1.4, "618": 1.4, "668": 1.2, "688": 1.4, "520": 0.8, "1314": 0.6}.items():
        if pat in p:
            rescue += bonus
            terms.append(f"{pat}=+{bonus:g}")
            reasons.append(f"contains {pat}")

    for pat, bonus in {"88": 0.7, "99": 0.5, "66": 0.5, "888": 1.4, "999": 1.1, "666": 1.1}.items():
        if pat in p:
            rescue += bonus
            terms.append(f"{pat}=+{bonus:g}")
            if len(pat) >= 3:
                reasons.append(f"contains repeated {pat}")

    for sub, l, st in _seq_matches(p, 3):
        bonus = min(4.5, 0.70 * l)
        rescue += bonus
        terms.append(f"seq{sub}=+{bonus:g}")
        reasons.append(f"contains ascending sequence {sub}" if st == 1 else f"contains descending sequence {sub}")

    ok_rep4, rep4 = _repeat_block_score(p, 4)
    ok_rep2, rep2 = _repeat_block_score(p, 2)
    ok_mirror, mirror = _mirror_score(p)
    if ok_rep4:
        rescue += 2.5; terms.append(f"rep4 {rep4}=+2.5"); reasons.append(f"has repeated block {rep4}")
    if ok_rep2:
        rescue += 1.8; terms.append(f"rep2 {rep2}=+1.8"); reasons.append(f"has repeated block {rep2}")
    if ok_mirror:
        rescue += 1.8; terms.append(f"mirror {mirror}=+1.8"); reasons.append(f"has mirror pattern {mirror}")
    if p[-4:] == p[-4:][::-1] and len(set(p[-4:])) > 1:
        rescue += 1.2; terms.append("pal_tail4=+1.2"); reasons.append(f"has palindrome tail {p[-4:]}")

    rescue_capped = min(rescue, 6.5)
    raw += rescue_capped
    terms.append(f"rescue_cap={rescue_capped:.1f}")
    if n4 and rescue_capped >= 3.0:
        reasons.append("lucky structure rescues a number containing 4")
    if "4" not in p:
        raw += 0.8
        terms.append("no4=+0.8")
        reasons.append("contains no digit 4")

    final = max(0.0, min(10.0, round(raw, 1)))
    if not reasons:
        reasons.append("no strong lucky or unlucky pattern detected")
    reason = "; ".join(dict.fromkeys(reasons[:6]))
    return final, f"{reason}. ({' + '.join(terms)} = {final:g})"

def ensure_phone_sco3_model():
    if os.path.exists(PHONE_SCO3_MODEL_FILE) and os.path.exists(PHONE_SCO3_FEATURE_FILE):
        try:
            joblib.load(PHONE_SCO3_MODEL_FILE)
            joblib.load(PHONE_SCO3_FEATURE_FILE)
            return True
        except ModuleNotFoundError as e:
            ems_log(f"existing phone sco3 model requires missing module {e.name}; retraining")
        except Exception as e:
            ems_log(f"existing phone sco3 model cannot be loaded; retraining: {e}")
    try:
        info = train_phone_sco3_ml(
            TRAIN_FILE,
            PHONE_SCO3_MODEL_FILE,
            PHONE_SCO3_FEATURE_FILE,
            PHONE_SCO3_REPORT_FILE,
        )
        ems_log(f"phone sco3 model trained: {info}")
        return True
    except Exception as e:
        ems_log(f"phone sco3 training skipped: {e}")
        return False

def predict_phone_sco3(phone_tails, batch_size=2048):
    scores = [float("nan")] * len(phone_tails)
    reasons = [""] * len(phone_tails)
    idx = [i for i, p in enumerate(phone_tails) if _is_digit8(p)]
    if not idx:
        return scores, reasons
    if not ensure_phone_sco3_model():
        return scores, reasons
    for s in range(0, len(idx), batch_size):
        batch_idx = idx[s : s + batch_size]
        batch_phones = [phone_tails[i] for i in batch_idx]
        pred, rsn = predict_phone_sco3_ml(batch_phones, PHONE_SCO3_MODEL_FILE, PHONE_SCO3_FEATURE_FILE)
        for j, val, rr in zip(batch_idx, pred, rsn):
            scores[j] = round(float(val), 1)
            reasons[j] = rr
    return scores, reasons

def ml_phone(phone_list):
    """Batch scoring for 8-digit phone tail. Returns sco1/sco2/sco3 and matching reasons."""
    out = {"sco1": [], "sco2": [], "sco3": [], "reason1": [], "reason2": [], "reason3": []}
    if phone_list is None:
        out["sco"] = out["sco1"]
        out["reason"] = out["reason1"]
        return out

    tails = [_normalize_phone_tail(p) for p in list(phone_list)]
    for p in tails:
        s1, r1 = score_phone_sco1(p)
        s2, r2 = score_phone_sco2(p)
        out["sco1"].append(s1)
        out["reason1"].append(r1)
        out["sco2"].append(s2)
        out["reason2"].append(r2)

    s3, r3 = predict_phone_sco3(tails)
    out["sco3"] = s3
    out["reason3"] = r3
    out["sco"] = out["sco1"]
    out["reason"] = out["reason1"]
    return out


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 GEO model helpers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def _build_geo_train_df():
    df = pd.read_excel(TRAIN_FILE)
    missing = [c for c in [GEO_TEXT_COL, GEO_LABEL_COL] if c not in df.columns]
    if missing:
        raise ValueError(f"geo training data missing columns: {missing}")

    df = df[[GEO_TEXT_COL, GEO_LABEL_COL]].copy()
    df[GEO_TEXT_COL] = df[GEO_TEXT_COL].fillna("").astype(str).str.strip()
    df[GEO_LABEL_COL] = df[GEO_LABEL_COL].fillna("").astype(str).str.strip()
    df = df[(df[GEO_TEXT_COL] != "") & (df[GEO_LABEL_COL] != "")]
    if df.empty:
        raise ValueError("geo training data has no usable rows")
    return df

def _build_geo_label_maps_from_trainfile():
    df = _build_geo_train_df()
    uniq = sorted(df[GEO_LABEL_COL].unique())
    l2i = {l: i for i, l in enumerate(uniq)}
    i2l = {i: l for l, i in l2i.items()}
    return l2i, i2l

def train_geo_model():
    global GEO_TOKENIZER, GEO_MODEL, GEO_LABEL2ID, GEO_ID2LABEL
    torch.set_grad_enabled(True)

    ems_log(f"geo trained model file {display_path(GEO_TRAINED_MODEL_FILE)} not found")
    ems_log(f"train geo model using model files {display_path(BASE_MODEL_DIR)} and data file {display_path(TRAIN_FILE)}")

    df = _build_geo_train_df()
    texts = df[GEO_TEXT_COL].tolist()
    labels = df[GEO_LABEL_COL].astype(str)

    uniq = sorted(labels.unique())
    GEO_LABEL2ID = {l: i for i, l in enumerate(uniq)}
    GEO_ID2LABEL = {i: l for l, i in GEO_LABEL2ID.items()}
    y = labels.map(GEO_LABEL2ID).tolist()

    ems_log("Loading geo tokenizer...")
    GEO_TOKENIZER = AutoTokenizer.from_pretrained(BASE_MODEL_DIR)

    ems_log("Tokenizing geo training texts...")
    enc = GEO_TOKENIZER(texts, truncation=True, padding=True, max_length=96)

    ems_log("Loading geo base model...")
    GEO_MODEL = AutoModelForSequenceClassification.from_pretrained(
        BASE_MODEL_DIR, num_labels=len(uniq)
    ).to(DEVICE)

    GEO_MODEL.train()
    GEO_MODEL.requires_grad_(True)
    gpu_mem()

    args = TrainingArguments(
        output_dir=GEO_OUTPUT_MODEL_DIR,
        per_device_train_batch_size=16,
        num_train_epochs=2,
        logging_steps=50,
        logging_strategy="steps",
        disable_tqdm=False,
        save_strategy="no",
        fp16=torch.cuda.is_available(),
        report_to=[],
    )

    ems_log("Geo training started...")
    Trainer(model=GEO_MODEL, args=args, train_dataset=DxDataset(enc, y)).train()

    ems_log("Saving geo trained model...")
    GEO_MODEL.save_pretrained(GEO_OUTPUT_MODEL_DIR)
    GEO_TOKENIZER.save_pretrained(GEO_OUTPUT_MODEL_DIR)

    with open(GEO_LABELS_FILE, "w", encoding="utf-8") as f:
        json.dump(
            {"label2id": GEO_LABEL2ID, "id2label": {str(k): v for k, v in GEO_ID2LABEL.items()}},
            f,
            ensure_ascii=False,
        )

    ems_log(f"geo trained model saved: {display_path(GEO_TRAINED_MODEL_FILE)}")
    ems_log("geo training finished")

def load_geo_model_once():
    global GEO_TOKENIZER, GEO_MODEL, GEO_LABEL2ID, GEO_ID2LABEL
    if GEO_MODEL is not None:
        return

    if not os.path.exists(GEO_TRAINED_MODEL_FILE):
        train_geo_model()
    else:
        ems_log(f"geo trained model file {display_path(GEO_TRAINED_MODEL_FILE)} found")

    ems_log("Loading geo tokenizer/model into memory...")
    GEO_TOKENIZER = AutoTokenizer.from_pretrained(GEO_OUTPUT_MODEL_DIR)
    geo_load_dtype = torch.float16 if DEVICE == "cuda" else torch.float32
    GEO_MODEL = AutoModelForSequenceClassification.from_pretrained(
        GEO_OUTPUT_MODEL_DIR, dtype=geo_load_dtype
    ).to(DEVICE)

    GEO_MODEL.eval()
    GEO_MODEL.requires_grad_(False)
    torch.set_grad_enabled(False)

    if os.path.exists(GEO_LABELS_FILE):
        ems_log(f"geo label file {display_path(GEO_LABELS_FILE)} found; loading")
        with open(GEO_LABELS_FILE, "r", encoding="utf-8") as f:
            obj = json.load(f)
        GEO_LABEL2ID = obj["label2id"]
        GEO_ID2LABEL = {int(k): v for k, v in obj["id2label"].items()}
    else:
        ems_log(f"geo label file {display_path(GEO_LABELS_FILE)} not found; rebuilding from {display_path(TRAIN_FILE)}")
        GEO_LABEL2ID, GEO_ID2LABEL = _build_geo_label_maps_from_trainfile()
        with open(GEO_LABELS_FILE, "w", encoding="utf-8") as f:
            json.dump(
                {"label2id": GEO_LABEL2ID, "id2label": {str(k): v for k, v in GEO_ID2LABEL.items()}},
                f,
                ensure_ascii=False,
            )
        ems_log("geo label file recreated")

    gpu_mem()

def ml_geo(df, batch_size=512):
    batch_size = int(batch_size)
    load_geo_model_once()

    if GEO_TEXT_COL not in df.columns:
        raise ValueError(f"geo input data missing column: {GEO_TEXT_COL}")

    raw_texts = df[GEO_TEXT_COL].fillna("").astype(str).str.strip().tolist()
    out = [""] * len(raw_texts)
    infer_idx = [i for i, t in enumerate(raw_texts) if t]
    texts = [raw_texts[i] for i in infer_idx]

    ems_log(f"Starting geo inference: {len(texts)}/{len(raw_texts)} samples (batch={batch_size}, device={DEVICE})")
    if not texts:
        return out

    use_cuda = DEVICE == "cuda"
    amp_ctx = torch.amp.autocast("cuda", dtype=torch.float16) if use_cuda else nullcontext()

    preds_all = []
    for i in range(0, len(texts), batch_size):
        if i % 10000 == 0:
            ems_log(f"Geo inference progress: {i}/{len(texts)}")

        batch = texts[i : i + batch_size]
        enc = GEO_TOKENIZER(batch, truncation=True, padding=True, max_length=96, return_tensors="pt")
        enc = {k: v.to(DEVICE) for k, v in enc.items()}

        with torch.inference_mode(), amp_ctx:
            logits = GEO_MODEL(**enc).logits
            probs = torch.softmax(logits, dim=1)

        preds = torch.argmax(probs, dim=1).cpu().numpy()
        preds_all.extend(preds)

        del enc, logits, probs

    if use_cuda:
        torch.cuda.empty_cache()

    for row_i, pred_i in zip(infer_idx, preds_all):
        out[row_i] = GEO_ID2LABEL[int(pred_i)]

    ems_log("Geo inference finished")
    return out


def release_geo_model():
    """Release the geo transformer before disease CV/training."""
    global GEO_TOKENIZER, GEO_MODEL
    GEO_MODEL = None
    GEO_TOKENIZER = None
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    ems_log("Released geo model memory")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 DX model helpers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def ml_dx(df, data_name=None, batch_size=512):
    batch_size = int(batch_size)
    load_model_once()
    if data_name is None:
        data_name = f"rows={len(df)}"
    ems_log(f"now process disease classification: {data_name}")

    texts = build_texts_from_df(df)
    kw, kw_reason = keyword_predict_df(df, model=KEYWORD_MODEL)

    preds_all, conf_all = _predict_texts(MODEL, TOKENIZER, texts, batch_size=batch_size)
    ml = [ID2LABEL[int(i)] for i in preds_all]
    ml_reason = [f"{ID2LABEL[int(i)]} (p={c:.3f})" for i, c in zip(preds_all, conf_all)]
    ems_log("Disease inference finished")
    return {
        "kw": kw,
        "kw_reason": kw_reason,
        "ml": ml,
        "ml_reason": ml_reason,
        "ml_confidence": [float(x) for x in conf_all],
    }
