#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$$" > /mnt/d/scripts/gu/.audit/report-shiny.pid
export GU_SQLITE=/mnt/d/analysis/gu/gu.sqlite
export GU_PHYML_REPORT_DIR=/mnt/d/scripts/gu/.audit/final_integration/review
export GU_SHINY_PORT=3842
unset R_LIBS R_LIBS_USER
export R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null
exec /home/huangj/anaconda3/envs/gu/bin/Rscript --vanilla /mnt/d/scripts/gu/f/shiny.R
