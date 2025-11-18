#!/bin/bash
# Utility functions for genome assembly pipeline
# 基因组组装流程工具函数库

set -o errexit
set -o pipefail

# ========== 颜色定义 ==========
readonly RED='\033[0;31m'
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly BLUE='\033[0;34m'
readonly NC='\033[0m' # No Color

# ========== 日志函数 ==========
log() {
    local level=$1
    shift
    local message="$@"
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    
    case $level in
        INFO)
            echo -e "${timestamp} ${GREEN}[ INFO]${NC} $message" | tee -a "${LOG_FILE:-/dev/null}"
            ;;
        WARN)
            echo -e "${timestamp} ${YELLOW}[ WARN]${NC} $message" | tee -a "${LOG_FILE:-/dev/null}"
            ;;
        ERROR)
            echo -e "${timestamp} ${RED}[ERROR]${NC} $message" | tee -a "${LOG_FILE:-/dev/null}" >&2
            ;;
        DEBUG)
            if [[ "${DEBUG:-false}" == "true" ]]; then
                echo -e "${timestamp} ${BLUE}[DEBUG]${NC} $message" | tee -a "${LOG_FILE:-/dev/null}"
            fi
            ;;
    esac
}

# ========== 状态记录函数 ==========
record_step() {
    local step_name=$1
    local status=$2  # START, FINISH, FAIL
    echo -e "${step_name}\t${status}\t$(date +"%Y-%m-%d %H:%M:%S")" >> "${STATUS_LOG:-status.log}"
}

# 检查步骤是否已完成
is_step_completed() {
    local step_name=$1
    if [[ -f "${STATUS_LOG}" ]]; then
        grep -q "^${step_name}\s*FINISH" "${STATUS_LOG}" && return 0
    fi
    return 1
}

# 获取上一个完成的步骤
get_last_completed_step() {
    if [[ -f "${STATUS_LOG}" ]]; then
        grep "FINISH" "${STATUS_LOG}" | tail -n 1 | awk '{print $1}'
    fi
}

# ========== 环境检查函数 ==========
check_command() {
    local cmd=$1
    if ! command -v "$cmd" &> /dev/null; then
        log ERROR "Command not found: $cmd"
        return 1
    fi
    return 0
}

check_file() {
    local file=$1
    if [[ ! -f "$file" ]]; then
        log ERROR "File not found: $file"
        return 1
    fi
    return 0
}

check_dir() {
    local dir=$1
    if [[ ! -d "$dir" ]]; then
        log ERROR "Directory not found: $dir"
        return 1
    fi
    return 0
}

# ========== 文件处理函数 ==========
create_dir() {
    local dir=$1
    if [[ ! -d "$dir" ]]; then
        mkdir -p "$dir"
        log INFO "Created directory: $dir"
    fi
}

# 检查fastq文件配对
check_fastq_pairs() {
    local dir=$1
    local type=$2
    local list_file="${dir}/${type}_file_list.txt"
    
    # 列出所有fastq文件
    find "$dir" -name "*.fastq.gz" -o -name "*.fq.gz" | sort > "$list_file"
    
    # 检查R1/R2配对
    local r1_count=$(grep -c "R1" "$list_file" || true)
    local r2_count=$(grep -c "R2" "$list_file" || true)
    
    if [[ $r1_count -ne $r2_count ]]; then
        log ERROR "Unmatched R1/R2 pairs in $dir (R1:$r1_count, R2:$r2_count)"
        return 1
    fi
    
    log INFO "Found $r1_count paired-end fastq files in $dir"
    return 0
}

# ========== 时间和资源统计 ==========
get_elapsed_time() {
    local start_time=$1
    local end_time=$(date +%s)
    local elapsed=$((end_time - start_time))
    
    local hours=$((elapsed / 3600))
    local minutes=$(((elapsed % 3600) / 60))
    local seconds=$((elapsed % 60))
    
    printf "%02d:%02d:%02d" $hours $minutes $seconds
}

# 记录步骤开始时间
start_timer() {
    echo $(date +%s)
}

# 打印步骤耗时
print_elapsed() {
    local start_time=$1
    local step_name=$2
    local elapsed=$(get_elapsed_time $start_time)
    log INFO "$step_name completed in $elapsed"
}

# ========== 数据验证函数 ==========
validate_fasta() {
    local fasta=$1
    if [[ ! -f "$fasta" ]]; then
        log ERROR "FASTA file not found: $fasta"
        return 1
    fi
    
    # 简单验证FASTA格式
    if ! grep -q "^>" "$fasta"; then
        log ERROR "Invalid FASTA format: $fasta"
        return 1
    fi
    
    local seq_count=$(grep -c "^>" "$fasta")
    log INFO "FASTA file contains $seq_count sequences"
    return 0
}

# ========== 清理函数 ==========
cleanup_intermediate() {
    local dir=$1
    if [[ "${KEEP_INTERMEDIATE:-true}" == "false" ]]; then
        log INFO "Cleaning up intermediate files in $dir"
        rm -rf "$dir"/*.tmp "$dir"/*.temp
    fi
}

# ========== 错误处理 ==========
error_exit() {
    local message=$1
    local step=$2
    log ERROR "$message"
    record_step "$step" "FAIL"
    exit 1
}

# ========== 配置文件解析 ==========
# 使用Python解析YAML (如果需要复杂解析)
parse_yaml() {
    local yaml_file=$1
    local prefix=$2
    python3 - <<EOF
import yaml
import sys

with open('$yaml_file', 'r') as f:
    config = yaml.safe_load(f)

def print_vars(d, prefix=''):
    for k, v in d.items():
        key = f"{prefix}{k}".upper()
        if isinstance(v, dict):
            print_vars(v, f"{prefix}{k}_")
        elif isinstance(v, list):
            print(f"{key}=({' '.join(map(str, v))})")
        else:
            print(f"{key}={v}")

print_vars(config, '$prefix')
EOF
}

# ========== 进度显示 ==========
show_progress() {
    local current=$1
    local total=$2
    local step_name=$3
    local percent=$((current * 100 / total))
    log INFO "Progress: [$current/$total] $percent% - $step_name"
}

# ========== 资源监控 ==========
check_disk_space() {
    local min_gb=$1
    local path=${2:-.}
    local available_gb=$(df -BG "$path" | awk 'NR==2 {print $4}' | sed 's/G//')
    
    if [[ $available_gb -lt $min_gb ]]; then
        log WARN "Low disk space: ${available_gb}GB available, ${min_gb}GB recommended"
        return 1
    fi
    log INFO "Disk space OK: ${available_gb}GB available"
    return 0
}

check_memory() {
    local min_gb=$1
    local available_gb=$(free -g | awk 'NR==2 {print $7}')
    
    if [[ $available_gb -lt $min_gb ]]; then
        log WARN "Low memory: ${available_gb}GB available, ${min_gb}GB recommended"
        return 1
    fi
    log INFO "Memory OK: ${available_gb}GB available"
    return 0
}

# ========== 导出所有函数 ==========
export -f log
export -f record_step
export -f is_step_completed
export -f get_last_completed_step
export -f check_command
export -f check_file
export -f check_dir
export -f create_dir
export -f check_fastq_pairs
export -f get_elapsed_time
export -f start_timer
export -f print_elapsed
export -f validate_fasta
export -f cleanup_intermediate
export -f error_exit
export -f show_progress
export -f check_disk_space
export -f check_memory
