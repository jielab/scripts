#!/usr/bin/env bash
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    case ${1:-} in
        -h|--help|help)
            cat <<'ALLELE_HELP'
Usage: bash f/allele3_bsub.sh [bsub OPTIONS] COMMAND [ARGS...]
内部 LSF 包装器，固定任务名 allele3，其他参数交给 bsub。
常规分析请使用顶层 .sh；它们管理依赖和提交记录。
独立小任务示例：bash f/allele3_bsub.sh -q 33 -n 10 -W 00:10 /bin/hostname
-h / --help 仅显示本帮助，不提交任务。
ALLELE_HELP
            exit 0 ;;
    esac
fi
# Scheduler display name for the public allele3 workflow.
set -euo pipefail
exec bsub -J allele3 "$@"
