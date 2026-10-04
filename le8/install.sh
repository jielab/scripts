#!/usr/bin/env bash


# 🚩 install
# Explicit installation only. le8.sh never sources or runs this installer.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
profile=report
python_bin=${LE8_REPORT_PYTHON:-${PYTHON_BIN:-${HOME}/anaconda3/envs/le8/bin/python3}}
r_bin=${R_BIN:-Rscript}
install_r=TRUE
while (($#)); do
	case "$1" in
	--report)
		profile=report
		shift
		;;
	--abm)
		profile=abm
		shift
		;;
	--all)
		profile=tabicl
		shift
		;;
	--python)
		python_bin=${2:?--python requires FILE}
		shift 2
		;;
	--r-bin)
		r_bin=${2:?--r-bin requires FILE}
		shift 2
		;;
	--no-r)
		install_r=FALSE
		shift
		;;
	-h | --help)
		cat <<'HELP'
Usage: ./install.sh [--report|--abm|--all] [--python FILE] [--r-bin FILE] [--no-r]
Default: report Python dependencies + Shiny and workbook R dependencies.
--abm adds reference ABM; --all also adds pinned TabICLv2.
All Python profiles come from the single requirements.txt file.
External model weights and research data are never downloaded by this script.
Conda environment recipes are in environment*.yml in the project root.
HELP
		exit 0
		;;
	*)
		echo "ERROR: Unknown installation option: $1" >&2
		exit 2
		;;
	esac
done
if ! command -v "$python_bin" >/dev/null 2>&1; then python_bin=python3; fi
requirement_file=$(mktemp)
trap 'rm -f -- "$requirement_file"' EXIT
awk -v profile="$profile" '
  /^# \[report\]/ {section="report"; next}
  /^# \[abm\]/ {section="abm"; next}
  /^# \[tabicl\]/ {section="tabicl"; next}
  section=="report" || (profile!="report" && section=="abm") || profile=="tabicl" {print}
' "$ROOT/requirements.txt" >"$requirement_file"
"$python_bin" -m pip install -r "$requirement_file"
if [[ $install_r == TRUE ]]; then
	"$r_bin" - <<'R_INSTALL'
# Explicit dependency installation; app.R never installs software or runs analyses.
p <- c("shiny" , "DT" , "data.table" , "ggplot2" , "digest" , "openxlsx" , "jsonlite" , "R.utils")
m <- p[!vapply(p , requireNamespace , logical(1) , quietly = TRUE)]
if (length(m)) install.packages(m , repos = "https://cloud.r-project.org")
R_INSTALL
fi
