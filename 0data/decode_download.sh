#!/usr/bin/env bash
# 下载 deCODE largescaleplasma-2023 bucket 中的全部共享文件。
# 默认目标：Windows E:\gwas\prot_decode（WSL 中为 /mnt/e/gwas/prot_decode）。
# 依赖 rclone；Ubuntu/WSL 可先运行：sudo apt update && sudo apt install -y rclone
# 使用：bash decode_download.sh
# 自定义目录：bash decode_download.sh /absolute/path/to/output
# 8 个文件并发下载；重复运行会跳过已完成且未变化的文件。
# 中断时尚未完成的文件可能重新下载，不保证跨进程的单文件断点续传。

set +x  # 避免 bash -x 将下方凭据打印到终端。
set -euo pipefail

if (( $# > 1 )); then
    printf '用法：bash %s [下载目录]\n' "$0" >&2
    exit 2
fi

if ! command -v rclone >/dev/null 2>&1; then
    cat >&2 <<'EOF'
未找到 rclone，请安装后重新运行本脚本。
Ubuntu/WSL：sudo apt update && sudo apt install -y rclone
Windows/Git Bash：安装 Windows 版 rclone 并将其加入 PATH。
官方下载：https://rclone.org/downloads/
EOF
    exit 127
fi

if (( $# == 1 )); then
    dest_dir="$1"
    if [[ -z "$dest_dir" ]]; then
        printf '错误：下载目录不能为空。\n' >&2
        exit 2
    fi
else
    case "$(uname -s)" in
        MINGW*|MSYS*|CYGWIN*)
            dest_dir='E:/gwas/prot_decode'
            if [[ ! -d 'E:/' ]]; then
                printf '错误：未找到 E 盘，请检查磁盘是否已连接。\n' >&2
                exit 1
            fi
            ;;
        Linux*)
            # 确认 E 盘已挂载，防止将大量文件误写到 WSL 系统盘。
            if ! mountpoint -q /mnt/e; then
                cat >&2 <<'EOF'
错误：E 盘未挂载到 /mnt/e。
请在 WSL 中先运行：
  sudo mkdir -p /mnt/e
  sudo mount -t drvfs E: /mnt/e
也可将已挂载的目标目录作为第一个参数传入。
EOF
                exit 1
            fi
            dest_dir='/mnt/e/gwas/prot_decode'
            ;;
        *)
            printf '请通过第一个参数指定下载目录。\n' >&2
            exit 2
            ;;
    esac
fi

# 使用进程环境配置 S3，无需运行 rclone config。
# 凭据只用于连接用户提供的 deCODE S3 服务。
export RCLONE_CONFIG_DECODE_PLASMA_TYPE='s3'
export RCLONE_CONFIG_DECODE_PLASMA_PROVIDER='Other'
export RCLONE_CONFIG_DECODE_PLASMA_ENV_AUTH='false'
export RCLONE_CONFIG_DECODE_PLASMA_ACCESS_KEY_ID='NLMINHZK0WJ0FJ8YHWHG'
export RCLONE_CONFIG_DECODE_PLASMA_SECRET_ACCESS_KEY='4yaXQ+4cEe6WYg01niHzIu1kIB3og+72ra6f4CHJ'
export RCLONE_CONFIG_DECODE_PLASMA_ENDPOINT='https://s3-ext.decode.is:443'
export RCLONE_CONFIG_DECODE_PLASMA_REGION='us-east-1'
export RCLONE_CONFIG_DECODE_PLASMA_FORCE_PATH_STYLE='true'

source_path='decode_plasma:largescaleplasma-2023'
mkdir -p -- "$dest_dir"
log_file="${dest_dir%/}/decode_download_$(date +%Y%m%d_%H%M%S)_$$.log"

printf '来源：%s\n目标：%s\n并发下载：8\n日志：%s\n' \
    "$source_path" "$dest_dir" "$log_file"

# copy 保留 bucket 内的目录结构。关闭单文件多流，限制为 8 个文件下载流。
if rclone copy "$source_path" "$dest_dir" \
    --transfers 8 \
    --checkers 8 \
    --multi-thread-streams 0 \
    --retries 10 \
    --low-level-retries 20 \
    --retries-sleep 10s \
    --contimeout 30s \
    --timeout 10m \
    --progress \
    --stats 10s \
    --log-level INFO \
    --log-file "$log_file"; then
    printf '\n下载完成：%s\n' "$dest_dir"
else
    exit_code=$?
    printf '\n下载未全部完成（退出码 %s）。请查看日志：%s\n' \
        "$exit_code" "$log_file" >&2
    printf '修复问题后重新运行相同命令，即可继续下载剩余文件。\n' >&2
    exit "$exit_code"
fi
