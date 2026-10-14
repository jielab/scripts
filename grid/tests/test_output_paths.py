"""Exercise stable workspaces across reuse, changed inputs and failed rebuilds."""

import os
from pathlib import Path
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[1]


SHELL = r'''
ROOT=$TEST_ROOT
GRID_OUTPUT_ROOT=$TEST_OUT
GRID_TRAIT=height
GRID_SCORE_DIR=$TEST_OUT/scores
GRID_PCA_FILE=$TEST_OUT/pca.tsv
GRID_MED_FILE=$TEST_OUT/centers.tsv
GRID_PCA_WEIGHT=$TEST_OUT/weights.tsv
GRID_PHE_FILE=$TEST_OUT/phe.rds
GRID_ANCESTRY_FILE=$TEST_OUT/ancestry.tsv
GRID_IMP_DIR=$TEST_OUT/genotypes
GRID_COV_PCS=10
GRID_DISTANCE_PCS=10
GRID_DISCO_A=1,1,1,1
GRID_REGRESS_PCA=TRUE
GRID_REMOVE=''
GRID_CHECK=FALSE
GRID_DRY_RUN=FALSE
GRID_REPLACE=${TEST_REPLACE:-FALSE}
GRID_JOBS=4
GRID_ETHNICITY_COL=ethnicity
GRID_ANCESTRY_PROB_MIN=0.90
GRID_ANCHOR_MAX_PER_GROUP=10000
suffix=${TEST_SUBSET:-}
work=$TEST_OUT/$TEST_METHOD/work
logdir=$TEST_OUT/$TEST_METHOD/log
POPS=(AFR EAS EUR SAS)
need() { [[ -s $1 ]]; }
_grid_die() { echo "$*" >&2; exit 1; }
join_comma() { local IFS=,; echo "$*"; }
publish() { cp -- "$1" "$2"; }
grid_select_r() { GRID_R=(fake_r); }
grid_run() { "$@"; }
grid_run_logged() { local dest=$1; shift; "$@" >"$dest" 2>&1; }
python3() {
    if [[ $1 == */0.common.py ]]; then
        case $2 in
            cache-path) echo "$TEST_OUT/score-lock" ;;
            signature) echo "$TEST_SIGNATURE" ;;
            disco-inputs)
                for name in centers pca AFR EAS EUR SAS; do printf 'input\n' >"$6/$name.tsv"; done
                echo 'Disco prepared'
                echo prepare >>"$TEST_OUT/events"
                ;;
            disco-output) need "$3" ;;
            publish) cp -- "$4" "$5" ;;
            *) return 99 ;;
        esac
    elif [[ $1 == */0.pca.py ]]; then
        need "$3"
    else
        command python3 "$@"
    fi
}
fake_r() {
    if [[ $1 == */2.disco.R ]]; then
        [[ ${TEST_FAIL:-FALSE} == FALSE ]] || return 17
        while [[ $1 != -o ]]; do shift; done
        printf 'scores\n' >"$2.tsv.gz"
        printf 'coefficients\n' >"$2.coef.tsv.gz"
    elif [[ $1 == */results.R ]]; then
        cp -- "$3" "$4"
    else
        # PCA aggregation/QC are stubbed; chromosome cache decisions are real.
        echo aggregate >>"$TEST_OUT/events"
    fi
}
plink2() {
    while [[ $1 != --out ]]; do shift; done
    printf 'scores\n' >"$2.sscore"
    printf 'variants\n' >"$2.sscore.vars"
    echo project >>"$TEST_OUT/events"
}
if [[ $TEST_METHOD == disco ]]; then disco_run; else pca_run; fi
'''


def run_script(tmp_path, method, signature="input-A", **settings):
    script = ROOT / ("2.disco.sh" if method == "disco" else "0.pca.sh")
    definitions = script.read_text().split('\nsource "$ROOT/f/0.common.sh"')[0]
    env = dict(os.environ, TEST_ROOT=str(ROOT), TEST_OUT=str(tmp_path),
               TEST_METHOD=method, TEST_SIGNATURE=signature, **settings)
    return subprocess.run(["bash", "-c", definitions + "\n" + SHELL],
                          env=env, capture_output=True, text=True)


@pytest.fixture
def inputs(tmp_path):
    for name in ("pca.tsv", "centers.tsv", "weights.tsv", "phe.rds"):
        (tmp_path / name).write_text("input\n")
    for subset in ("", "chr1,2"):
        scores = tmp_path / "scores/height" / subset
        scores.mkdir(parents=True, exist_ok=True)
        (scores / "1csx.scores.rds").write_text("input\n")
    genotypes = tmp_path / "genotypes"
    genotypes.mkdir()
    for chrom in range(1, 23):
        for ext in ("pgen", "pvar", "psam"):
            (genotypes / f"chr{chrom}.{ext}").write_text("input\n")
    return tmp_path


def assert_success(result):
    assert result.returncode == 0, result.stdout + result.stderr


def test_disco_reuses_flat_workspace_and_isolates_subsets(inputs):
    assert_success(run_script(inputs, "disco"))
    run = inputs / "disco/work/height"
    assert (run / "centers.tsv").is_file()
    assert not any(p.is_dir() for p in run.iterdir())
    assert_success(run_script(inputs, "disco"))
    assert inputs.joinpath("events").read_text().splitlines() == ["prepare"]
    assert_success(run_script(inputs, "disco", TEST_SUBSET=".chr1,2"))
    assert (run / "chr1,2/centers.tsv").is_file()
    assert (inputs / "disco/log/height/chr1,2/prepare.log").is_file()
    assert (run / "disco.signature").read_text().strip() == "input-A"


def test_disco_failed_replace_cannot_reuse_old_completion(inputs):
    assert_success(run_script(inputs, "disco"))
    failed = run_script(inputs, "disco", TEST_REPLACE="TRUE", TEST_FAIL="TRUE")
    assert failed.returncode == 17
    assert not (inputs / "disco/work/height/disco.signature").exists()
    assert_success(run_script(inputs, "disco"))
    assert inputs.joinpath("events").read_text().splitlines() == ["prepare"] * 3
    assert_success(run_script(inputs, "disco", signature="input-B"))
    assert (inputs / "disco/work/height/disco.signature").read_text().strip() == "input-B"


def test_pca_reuses_only_chromosomes_with_matching_inputs(inputs):
    assert_success(run_script(inputs, "pca"))
    run = inputs / "pca/work/projection"
    assert (run / "chr15.sscore").is_file()
    assert inputs.joinpath("events").read_text().splitlines().count("project") == 22
    # An interrupted aggregate must reuse matching completed chromosomes.
    (inputs / "pca/projection.signature").unlink()
    assert_success(run_script(inputs, "pca"))
    assert inputs.joinpath("events").read_text().splitlines().count("project") == 22
    assert_success(run_script(inputs, "pca", signature="input-B"))
    assert inputs.joinpath("events").read_text().splitlines().count("project") == 44
    assert (run / "projection.signature").read_text().strip() == "input-B"
