#!/usr/bin/env python3
from __future__ import annotations
import argparse, fcntl, hashlib, json, os, shutil, subprocess, tempfile
from pathlib import Path
import numpy as np
import pandas as pd

RELEASE="2026-10-17.c1-lockpaper-v1"
PRIMARY_SELECTOR = "risk_resid_error"
COVERAGES = (.2, .4, .5, .6, .8, 1.0)
SEED = 2026
PREDICTABILITY_TARGET_MODEL = "elasticnet_full"
DEFAULT_UKB_REL_FILE = os.getenv("LE8_UKB_REL_FILE", "/mnt/f/gen/ukb/37/typ/ukb_rel.dat")
KINSHIP_PROTOCOL = "2026-10-17.ukb-king-components-v1"
GROUPING_SOURCE = "UKB KING kinship connected components"
KINSHIP_ALGORITHM = ("union-find connected components on full released relatedness graph; "
                     "then subset to all.rds; nonrelated participants assigned singleton groups")
KINSHIP_SEMANTICS = "UKB KING close-relatedness file; all released rows used as graph edges"
GROUPING_AUDIT_FIELDS = (
    "relatedness_source_file", "relatedness_source_sha256", "group_map_sha256",
    "pair_rows", "related_participant_N", "singleton_participant_N", "component_N",
    "largest_component_size", "grouping_unit", "grouping_source",
    "family_disjoint_outer_split", "family_disjoint_ABM_roles",
)


def grouping_audit(record):
    """Allowlist for public tables: never include IDs or component membership."""
    return {k: record[k] for k in GROUPING_AUDIT_FIELDS if k in record}


def _protocol_path(path):
    return Path(str(path) + ".protocol.json")


def _validate_ids(values, context, unique=False):
    s = pd.Series(values, dtype="string")
    if s.empty or s.isna().any() or s.str.contains(r"\s|^$", regex=True).any():
        raise ValueError(f"Missing/invalid IDs in {context}")
    if unique and s.duplicated().any():
        raise ValueError(f"Duplicate eid in {context}")
    return s.astype(str)


def validate_related_pairs(pairs):
    required = {"ID1", "ID2", "HetHet", "IBS0", "Kinship"}
    if not required <= set(pairs):
        raise ValueError("UKB KING file requires ID1 ID2 HetHet IBS0 Kinship")
    d = pairs.copy()
    for c in ("ID1", "ID2"):
        d[c] = _validate_ids(d[c], "UKB KING source").to_numpy()
    if d.ID1.eq(d.ID2).any():
        raise ValueError("UKB KING source contains self pairs (ID1 == ID2)")
    canonical = pd.Series([tuple(sorted(p)) for p in zip(d.ID1, d.ID2)])
    if canonical.duplicated().any():
        raise ValueError("UKB KING source contains duplicated canonical pairs")
    for c in ("HetHet", "IBS0", "Kinship"):
        d[c] = pd.to_numeric(d[c], errors="raise")
        if not np.isfinite(d[c]).all():
            raise ValueError(f"UKB KING {c} must be finite")
    stats = dict(pair_rows=len(d), unique_related_ids_in_source=len(set(d.ID1) | set(d.ID2)),
                 kinship_min=float(d.Kinship.min()), kinship_median=float(d.Kinship.median()),
                 kinship_max=float(d.Kinship.max()),
                 kinship_below_0_0442_N=int((d.Kinship < .0442).sum()))
    # Conservative source-semantics guard, NEVER an edge-selection threshold.
    # Even a small unexpected tail requires source review; no rows are discarded.
    if stats["kinship_below_0_0442_N"]:
        raise ValueError("Unexpected UKB KING distribution below 0.0442; review source, "
                         "do not filter edges: " + json.dumps(stats))
    return d, stats


def build_kinship_groups(pairs, analysis_ids):
    """Build the FULL source graph before projecting onto the analytical cohort."""
    d, stats = validate_related_pairs(pairs)
    ids = _validate_ids(analysis_ids, "all.rds", unique=True)
    parent, size = {}, {}

    def find(x):
        if x not in parent:
            parent[x], size[x] = x, 1
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in zip(d.ID1, d.ID2):
        a, b = find(a), find(b)
        if a != b:
            if (size[a], a) < (size[b], b):
                a, b = b, a
            parent[b] = a
            size[a] += size[b]
    # Numeric EIDs sort numerically; lexical fallback supports synthetic IDs.
    def key(x):
        return (0, int(x), x) if x.isdecimal() else (1, x, x)

    minima = {}
    for x in parent:
        root = find(x)
        if root not in minima or key(x) < key(minima[root]):
            minima[root] = x
    membership = {x: "KIN_" + minima[find(x)] for x in parent}
    result = pd.DataFrame(dict(eid=ids, group=[membership.get(x, "IND_" + x) for x in ids]))
    result = result.sort_values("eid").reset_index(drop=True)
    lookup = result.set_index("eid").group
    both = d.ID1.isin(lookup.index) & d.ID2.isin(lookup.index)
    if not d.loc[both, "ID1"].map(lookup).eq(d.loc[both, "ID2"].map(lookup)).all():
        raise AssertionError("Related pair assigned to different groups")
    sizes = result.groupby("group").size()
    related = result.group.str.startswith("KIN_")
    stats.update(all_rds_N=len(result), all_rds_eid_hash=eid_hash(ids),
                 all_rds_related_N=int(related.sum()), all_rds_singleton_N=int((~related).sum()),
                 component_N=len(sizes), largest_analysis_component_size=int(sizes.max()),
                 connected_components_with_at_least_2_analysis_participants=int((sizes >= 2).sum()),
                 component_size_1_N=int((sizes == 1).sum()), component_size_2_N=int((sizes == 2).sum()),
                 component_size_3_N=int((sizes == 3).sum()), component_size_4plus_N=int((sizes >= 4).sum()))
    return result, stats


def ensure_kinship_groups(ukb_phe, rel_file, output, r_bin="Rscript", force=False):
    """Create/validate private protocol data; changed sources fail closed by default."""
    output, rel_file = Path(output), Path(rel_file)
    allrds = Path(ukb_phe) / "Rdata/all.rds"
    for source in (rel_file, allrds):
        if not source.is_file() or not source.stat().st_size:
            raise ValueError(f"Required nonempty kinship protocol source: {source}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(str(output) + ".lock", "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        meta = _protocol_path(output)
        requested = dict(release=RELEASE, kinship_protocol=KINSHIP_PROTOCOL,
                         source_file=str(rel_file.resolve()), source_sha256=sha(rel_file),
                         all_rds_file=str(allrds.resolve()), all_rds_sha256=sha(allrds),
                         algorithm=KINSHIP_ALGORITHM, source_semantics=KINSHIP_SEMANTICS)
        reuse = output.is_file() and meta.is_file() and not force
        if output.exists() != meta.exists() and not force:
            raise RuntimeError("Kinship group map or protocol sidecar is missing")
        if reuse:
            old = json.loads(meta.read_text())
            for k, v in requested.items():
                if old.get(k) != v:
                    raise RuntimeError(f"Kinship group map protocol mismatch for {k}; archive and regenerate explicitly")
            if old.get("group_map_sha256") != sha(output):
                raise RuntimeError("Kinship group map SHA256 mismatch")
            saved = pd.read_csv(output, sep="\t", dtype=str)
            if list(saved) != ["eid", "group"]:
                raise ValueError("Kinship group map must contain only eid,group")
            ids = saved.eid
        else:
            with tempfile.TemporaryDirectory() as td:
                ph = Path(td) / "phenotype.tsv"
                _rds_export(allrds, ph, r_bin)
                ids = pd.read_csv(ph, sep="\t", dtype=str).eid
        pairs = pd.read_csv(rel_file, sep=r"\s+", dtype={"ID1": str, "ID2": str})
        groups, stats = build_kinship_groups(pairs, ids)
        record = {**requested, **stats, "grouping_unit": "family", "grouping_source": GROUPING_SOURCE}
        if reuse:
            if not saved.equals(groups):
                raise RuntimeError("Kinship group map disagrees with full source graph")
            for k, v in record.items():
                if old.get(k) != v:
                    raise RuntimeError(f"Kinship group map audit mismatch for {k}")
        # Refuse a source that changed while the map was being built/validated.
        if sha(rel_file) != requested["source_sha256"] or sha(allrds) != requested["all_rds_sha256"]:
            raise RuntimeError("Kinship protocol source changed during validation")
        if not reuse:
            with tempfile.TemporaryDirectory(dir=output.parent) as td:
                staged = Path(td) / output.name
                groups.to_csv(staged, sep="\t", index=False)
                staged.chmod(0o600)
                record["group_map_sha256"] = sha(staged)
                sidecar = _protocol_path(staged)
                sidecar.write_text(json.dumps(record, indent=2) + "\n")
                os.replace(staged, output)
                os.replace(sidecar, meta)
        return output


def kinship_record(group_file, *, ukb_phe=None, r_bin="Rscript"):
    """Recognize generated maps even when supplied explicitly with --group-file."""
    if not group_file:
        return {}
    meta = _protocol_path(group_file)
    if not meta.is_file():
        if Path(group_file).name == "ukb_kinship_groups.tsv":
            raise RuntimeError("Kinship group map protocol sidecar is missing")
        return {}
    record = json.loads(meta.read_text())
    if "kinship_protocol" not in record:
        if Path(group_file).name == "ukb_kinship_groups.tsv":
            raise RuntimeError("Kinship group map protocol marker is missing")
        return {}
    if record.get("kinship_protocol") != KINSHIP_PROTOCOL or record.get("release") != RELEASE:
        raise RuntimeError("Kinship group map release/protocol mismatch")
    if ukb_phe is not None:
        ensure_kinship_groups(ukb_phe, record["source_file"], group_file, r_bin=r_bin)
    if record.get("group_map_sha256") != sha(group_file) or record.get("source_sha256") != sha(record["source_file"]):
        raise RuntimeError("Kinship group map/source SHA256 mismatch")
    return record


def _kinship_audit(record):
    if not record:
        return {}
    return dict(relatedness_source_file=record["source_file"], relatedness_source_sha256=record["source_sha256"],
                group_map_sha256=record["group_map_sha256"], pair_rows=record["pair_rows"],
                related_participant_N=record["all_rds_related_N"], singleton_participant_N=record["all_rds_singleton_N"],
                component_N=record["component_N"], largest_component_size=record["largest_analysis_component_size"],
                grouping_unit="family", grouping_source=GROUPING_SOURCE)


def validate_seed(seed):
    if int(seed) != SEED or float(seed) != SEED:
        raise ValueError("C1 LOCKPAPER fixes the analysis seed at 2026; "
                         "do not change the outer split after inspecting prior results.")


def eid_hash(ids):
    """Language-independent SHA256 of sorted unique IDs, one per line."""
    return hashlib.sha256("".join(x + "\n" for x in sorted(set(map(str, ids)))).encode()).hexdigest()


def roster_record(path):
    path = Path(path)
    meta = path.with_suffix(path.suffix + ".protocol.json")
    record = json.loads(meta.read_text())
    validate_seed(record.get("seed", -1))
    if record.get("release") != RELEASE or record.get("roster_sha256") != sha(path):
        raise RuntimeError("Frozen outer roster integrity/protocol mismatch: " + str(path))
    if record.get("group_file") and record.get("group_file_sha256") != sha(record["group_file"]):
        raise RuntimeError("Frozen outer roster group-file SHA256 mismatch")
    if record.get("relatedness_source_file"):
        current = kinship_record(record["group_file"])
        if not current or any(record.get(k) != v for k, v in _kinship_audit(current).items()):
            raise RuntimeError("Frozen outer roster kinship protocol mismatch")
    return record


def configure_roster(args, *, settings=None, dry_run=False):
    """Resolve shared CLI/environment settings before inspecting omics or Y."""
    validate_seed(args.seed)
    if settings is None:
        from le8_common import shared_analysis_settings
        settings, _, _ = shared_analysis_settings(args)
    for key, value in settings.items():
        setattr(args, key, value)
    validate_seed(args.seed)
    # Native result transactions use a scratch analysis root. Protocol inputs
    # always live in the published root and are never packaged as result tables.
    root = Path(os.getenv("LE8_PUBLISHED_ROOT") or args.analysis_root)
    roster = Path(args.outer_roster) if args.outer_roster else (
        root / "_protocol" / f"c1_outer_roster_seed{args.seed}.tsv"
    )
    phe = args.ukb_phe or os.getenv("UKB_PHE", "/mnt/d/data/ukb")
    if not args.group_file and not args.group_col:
        source = Path(getattr(args, "ukb_rel_file", None) or os.getenv("LE8_UKB_REL_FILE", DEFAULT_UKB_REL_FILE))
        # stat() propagates permission/I/O errors; only genuine absence permits fallback.
        try:
            source.stat()
        except FileNotFoundError:
            if source.is_symlink():
                raise RuntimeError(f"Broken UKB KING source symlink: {source}")
            print(f"[LE8 C1 LOCK] WARNING: UKB KING source does not exist: {source}; "
                  "grouping_unit=individual_unverified; close-relative leakage is unverified", flush=True)
        else:
            group_map = root / "_protocol" / "ukb_kinship_groups.tsv"
            if not dry_run:
                ensure_kinship_groups(phe, source, group_map, r_bin=args.r_bin)
            else:
                print(f"[LE8 C1 LOCK] planned UKB KING groups: {source} -> {group_map}", flush=True)
            args.group_file, args.group_col = str(group_map), ""
    # Native C1 and ABM children must inherit exactly this mapping.
    os.environ["LE8_GROUP_FILE"] = str(args.group_file or "")
    os.environ["LE8_GROUP_COLUMN"] = str(args.group_col or "")
    if not dry_run:
        ensure_roster(ukb_phe=phe,
                      output=roster, seed=args.seed, group_file=args.group_file,
                      group_col=args.group_col, r_bin=args.r_bin)
        print(f"[LE8 C1 LOCK] shared outer roster: {roster}; SHA256={sha(roster)}", flush=True)
        record = roster_record(roster)
        unit = "family" if record["group_file"] or record["group_col"] else "individual_unverified"
        print(f"[LE8 C1 LOCK] release={RELEASE}; seed={record['seed']}; "
              f"outer_test_N={record['test_N']}; test_eid_hash={record['test_eid_hash']}; "
              f"grouping_unit={unit}; independent_family_support_verified={unit == 'family'}", flush=True)
        if record.get("relatedness_source_file"):
            print("C1 LOCKPAPER relatedness:\n"
                  f"  source: {record['relatedness_source_file']}\n"
                  f"  source pairs: {record['pair_rows']}\n"
                  f"  source SHA256: {record['relatedness_source_sha256']}\n"
                  f"  grouping: {record['grouping_source']}\n"
                  f"  analysis N: {record['N']}\n"
                  f"  related analysis participants: {record['related_participant_N']}\n"
                  f"  singleton analysis participants: {record['singleton_participant_N']}\n"
                  f"  connected components (including singletons): {record['component_N']}\n"
                  f"  largest component: {record['largest_component_size']}\n"
                  f"  group-map SHA256: {record['group_map_sha256']}\n"
                  f"  seed: {record['seed']}\n"
                  f"  outer test fraction: {record['test_fraction']:.2f}\n"
                  "  family-disjoint split: TRUE\n  bootstrap unit: family", flush=True)
    else:
        print(f"[LE8 C1 LOCK] planned shared outer roster: {roster}", flush=True)
    args.outer_roster = str(roster)
    os.environ["LE8_OUTER_ROSTER"] = str(roster)
    return roster

def sha(path, block=1<<20):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(block),b""): h.update(b)
    return h.hexdigest()

def stable_hash(x,seed):
    return int.from_bytes(hashlib.sha256(f"{seed}:{x}".encode()).digest()[:8],"big")

def _rds_export(rds,out,r_bin,group_col=""):
    code = r'''
args<-commandArgs(trailingOnly=TRUE)
x<-readRDS(args[1])
if(!is.data.frame(x)||!"eid"%in%names(x))stop("all.rds must be a data.frame with eid")
g<-args[3]
z<-data.frame(eid=as.character(x$eid),stringsAsFactors=FALSE)
if(nzchar(g)){
 if(!g%in%names(x))stop("Group column missing: ",g)
 z$group<-as.character(x[[g]])
}
data.table::fwrite(z,args[2],sep="\t",quote=FALSE,na="NA")
'''
    subprocess.run([r_bin,"--vanilla","-e",code,str(rds),str(out),group_col],check=True)

def _read_group_file(path,r_bin):
    path=Path(path)
    if path.suffix.lower()==".rds":
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/"group.tsv"
            code = r'''
args<-commandArgs(trailingOnly=TRUE)
x<-readRDS(args[1]);if(!is.data.frame(x))stop("group RDS must be data.frame")
if(!all(c("eid","group")%in%names(x)))stop("group file requires eid,group")
data.table::fwrite(data.frame(eid=as.character(x$eid),group=as.character(x$group)),
                   args[2],sep="\t",quote=FALSE,na="NA")
'''
            subprocess.run([r_bin,"--vanilla","-e",code,str(path),str(out)],check=True)
            return pd.read_csv(out,sep="\t",dtype=str)
    sep="\t" if path.suffix.lower() in {".tsv",".txt"} else ","
    return pd.read_csv(path,sep=sep,dtype=str)

def ensure_roster(*,ukb_phe,output,seed=2026,test_fraction=.20,group_file="",group_col="",r_bin="Rscript",force=False):
    validate_seed(seed)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(str(output) + ".lock", "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        return _ensure_roster(ukb_phe=ukb_phe, output=output, seed=seed,
                              test_fraction=test_fraction, group_file=group_file,
                              group_col=group_col, r_bin=r_bin, force=force)


def _ensure_roster(*,ukb_phe,output,seed,test_fraction,group_file,group_col,r_bin,force):
    if test_fraction != .20:
        raise ValueError("C1 LOCKPAPER fixes the outer test fraction at 0.20")
    ukb_phe=Path(ukb_phe); output=Path(output)
    allrds=ukb_phe/"Rdata"/"all.rds"
    if not allrds.is_file(): raise FileNotFoundError(allrds)
    kin = kinship_record(group_file, ukb_phe=ukb_phe, r_bin=r_bin)
    meta=output.with_suffix(output.suffix+".protocol.json")
    requested=dict(release=RELEASE,seed=int(seed),test_fraction=float(test_fraction),
                   source=str(allrds.resolve()),source_sha256=sha(allrds),
                   group_file=str(Path(group_file).resolve()) if group_file else "",
                   group_file_sha256=sha(group_file) if group_file else "",
                   group_col=str(group_col or ""))
    requested.update(_kinship_audit(kin))
    # --replace (model fitting) and --force must never rerandomize a frozen roster.
    if force and (output.exists() or meta.exists()):
        raise RuntimeError("Frozen outer roster cannot be overwritten; archive the legacy individual protocol explicitly")
    if output.exists() != meta.exists() and not force:
        raise RuntimeError("Frozen outer roster or its protocol sidecar is missing")
    if output.is_file() and meta.is_file() and not force:
        old=roster_record(output)
        for k,v in requested.items():
            if old.get(k)!=v: raise RuntimeError(f"Frozen outer roster protocol mismatch for {k}; archive the pre-kinship roster explicitly before this one-time upgrade")
        if group_file:
            _assert_outer_groups(output, _read_group_file(group_file, r_bin))
        return output
    output.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as td:
        td=Path(td); ph=td/"phenotype.tsv"
        _rds_export(allrds,ph,r_bin,group_col)
        d=pd.read_csv(ph,sep="\t",dtype=str)
        if d.eid.isna().any() or d.eid.duplicated().any() or d.eid.str.strip().eq("").any(): raise ValueError("Invalid/duplicate eid in all.rds")
        if group_file:
            g=_read_group_file(group_file,r_bin)
            if not {"eid","group"}<=set(g): raise ValueError("Group file requires eid,group")
            if g.eid.isna().any() or g.eid.duplicated().any() or g.group.isna().any(): raise ValueError("Invalid group file")
            if "group" in d:
                supplied=d.eid.map(g.set_index("eid").group)
                if not d.group.eq(supplied).all(): raise ValueError("Conflicting family definitions in file and phenotype column")
                d=d.drop(columns="group")
            d=d.merge(g[["eid","group"]],on="eid",how="left",validate="one_to_one")
            if d.group.isna().any(): raise ValueError("Group file does not cover all phenotype IDs")
        elif "group" not in d:
            d["group"]=d.eid
        if d.group.isna().any() or d.group.astype(str).str.strip().eq("").any(): raise ValueError("Incomplete group mapping")
        groups=sorted(d.group.astype(str).unique(),key=lambda x:(stable_hash(x,seed),x))
        if len(groups)<2: raise ValueError("At least two independent groups are required")
        ntest=max(1,min(len(groups)-1,round(len(groups)*test_fraction)))
        test=set(groups[:ntest])
        d["role"]=d.group.astype(str).map(lambda g:"test" if g in test else "training")
        if d.groupby("group").role.nunique().max()!=1: raise AssertionError("Family crosses outer roles")
        outtmp=td/"outer.tsv"; d.sort_values("eid")[["eid","role"]].to_csv(outtmp,sep="\t",index=False); os.replace(outtmp,output)
    record={**requested,"N":int(len(d)),"groups":int(d.group.nunique()),
            "test_N":int((d.role=="test").sum()),"training_N":int((d.role=="training").sum()),
            "test_groups":int(len(test)),"roster_sha256":sha(output),
            "group_definition":group_col or ("external eid,group mapping" if group_file else "individual ID; family mapping unavailable"),
            "grouping_unit":"family" if group_file or group_col else "individual_unverified",
            "grouping_source":GROUPING_SOURCE if kin else (group_col or ("supplied eid,group mapping" if group_file else "individual IDs; kinship unverified")),
            "family_disjoint_outer_split":bool(group_file or group_col),
            "independent_family_support_verified":bool(group_file or group_col),
            "test_eid_hash":eid_hash(d.loc[d.role=="test","eid"]),
            "rule":"sort groups by SHA256(seed:group); exact prespecified group count assigned to test; outcome/omics blind"}
    tmp=Path(str(meta)+".tmp");tmp.write_text(json.dumps(record,indent=2),encoding="utf-8");os.replace(tmp,meta)
    return output


def _assert_outer_groups(roster, groups):
    d = pd.read_csv(roster, sep="\t", dtype=str)
    _validate_ids(d.eid, "outer roster", unique=True)
    _validate_ids(groups.eid, "group map", unique=True)
    if groups.group.isna().any() or groups.group.astype(str).str.strip().eq("").any():
        raise ValueError("Incomplete group mapping")
    if set(d.eid) - set(groups.eid) or set(d.role) != {"training", "test"}:
        raise ValueError("Invalid/incomplete outer roster")
    d = d.merge(groups[["eid", "group"]], on="eid", validate="one_to_one")
    if not d.groupby("group").role.nunique().eq(1).all():
        raise ValueError("Kinship component crosses outer roles")


def archive_pre_kinship_roster(output):
    """Explicit one-time migration only; never archive an already grouped roster."""
    output = Path(output)
    # Native transactions hold this shared analysis lock while fitting. Refuse
    # migration under a live run, whose already-loaded code/settings cannot upgrade.
    root_key = hashlib.sha256(str(output.parent.parent.resolve()).encode()).hexdigest()[:20]
    with open(f"/tmp/le8-{root_key}.lock", "a") as analysis_lock, open(str(output) + ".lock", "a") as lock:
        try:
            fcntl.flock(analysis_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("C1 analysis is running; stop the old run before archiving its outer roster") from exc
        fcntl.flock(lock, fcntl.LOCK_EX)
        record = roster_record(output)
        if record.get("grouping_unit") != "individual_unverified" or record.get("group_file") or record.get("group_col"):
            raise RuntimeError("Only the legacy individual_unverified roster can be archived by this migration")
        history = output.parent / "_history/pre_kinship_20261017"
        history.mkdir(parents=True, exist_ok=True)
        paths = (output, _protocol_path(output))
        if any((history / p.name).exists() for p in paths):
            raise FileExistsError("Pre-kinship archive already exists; do not repeat the protocol upgrade")
        # Copy and verify both before removing either original, retaining recovery
        # data even if an interrupted migration needs manual repair.
        for p in paths:
            shutil.copy2(p, history / p.name)
            if sha(p) != sha(history / p.name):
                raise RuntimeError("Pre-kinship archive verification failed")
        for p in paths:
            p.unlink()
        return history

def main(argv=None):
    p=argparse.ArgumentParser()
    p.add_argument("--ukb-phe",required=True);p.add_argument("--output",required=True)
    p.add_argument("--seed",type=int,default=2026);p.add_argument("--test-fraction",type=float,default=.20)
    p.add_argument("--group-file",default="");p.add_argument("--group-col",default="");p.add_argument("--r-bin",default="Rscript")
    p.add_argument("--ukb-rel-file", default=os.getenv("LE8_UKB_REL_FILE", DEFAULT_UKB_REL_FILE))
    p.add_argument("--archive-pre-kinship", action="store_true",
                   help="One-time upgrade: archive only an existing individual_unverified roster")
    a=p.parse_args(argv)
    validate_seed(a.seed)
    if a.test_fraction != .20:
        p.error("C1 LOCKPAPER fixes the outer test fraction at 0.20")
    a.analysis_root = Path(a.output).parent.parent
    a.outer_roster = a.output
    if a.archive_pre_kinship:
        if a.group_file or a.group_col:
            p.error("Pre-kinship migration uses the UKB KING source")
        # Validate/build the new source before touching the previous frozen split.
        ensure_kinship_groups(a.ukb_phe, a.ukb_rel_file,
                              a.analysis_root / "_protocol/ukb_kinship_groups.tsv", r_bin=a.r_bin)
        print("Archived legacy protocol:", archive_pre_kinship_roster(a.output), flush=True)
    configure_roster(a, settings={})
if __name__=="__main__": main()
