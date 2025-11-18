#!/bin/bash
##############################################################################
# Genome Assembly Pipeline - Main Script
# 基因组组装流程 - 主脚本
# 
# Description: A modular pipeline for chromosome-level genome assembly
#              using PacBio HiFi and Hi-C data
# 
# Author: Genome Assembly Pipeline Team
# Version: 1.0.0
##############################################################################

set -o errexit
set -o pipefail
set -o nounset

# ========== 全局变量 ==========
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly WORK_DIR="${PWD}"
readonly CONFIG_FILE="${CONFIG_FILE:-${WORK_DIR}/config.yaml}"
readonly ENV_NAME="genome_assembly"  # 使用已有的虚拟环境名称

# 加载工具函数
source "${SCRIPT_DIR}/scripts/utils.sh"

# ========== 初始化 ==========
init_pipeline() {
    log INFO "=========================================="
    log INFO "Genome Assembly Pipeline Starting"
    log INFO "=========================================="
    log INFO "Work directory: ${WORK_DIR}"
    log INFO "Config file: ${CONFIG_FILE}"
    
    # 创建日志目录
    export LOG_DIR="${WORK_DIR}/logs"
    create_dir "$LOG_DIR"
    
    export LOG_FILE="${LOG_DIR}/assembly_$(date +%Y%m%d_%H%M%S).log"
    export STATUS_LOG="${LOG_DIR}/status.log"
    
    # 检查配置文件
    check_file "$CONFIG_FILE" || error_exit "Config file not found" "init"
    
    # 先激活环境（以便能使用Python解析YAML）
    activate_env_early
    
    # 解析配置文件
    parse_config
    
    log INFO "Initialization completed"
}

# ========== 提前激活环境 ==========
activate_env_early() {
    # 初始化conda
    if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
        source "$HOME/anaconda3/etc/profile.d/conda.sh"
    fi
    
    # 检查是否已经激活了正确的环境
    if [[ "${CONDA_DEFAULT_ENV:-}" != "${ENV_NAME}" ]]; then
        echo "ERROR: Please activate the genome_assembly environment first!"
        echo "Run the following command: conda activate genome_assembly"
        exit 1
    fi
}

# ========== 配置解析 ==========
parse_config() {
    log INFO "Parsing configuration..."
    
    # 使用Python解析YAML配置
    if command -v python3 &> /dev/null; then
        eval $(python3 - <<EOF
import yaml
import sys

try:
    with open('$CONFIG_FILE', 'r') as f:
        config = yaml.safe_load(f)
    
    print(f"export PROJECT_NAME='{config.get('project', {}).get('name', 'MyGenome')}'")
    print(f"export THREADS='{config.get('resources', {}).get('threads', 4)}'")
    print(f"export HIFI_DATA='{config.get('input', {}).get('hifi_data', '')}'")
    print(f"export HIC_DATA='{config.get('input', {}).get('hic_data', '')}'")
    print(f"export WGS_DATA='{config.get('input', {}).get('wgs_data', '')}'")
    print(f"export KRAKEN_DB='{config.get('databases', {}).get('kraken_db', '')}'")
    print(f"export BUSCO_DB='{config.get('databases', {}).get('busco_db', '')}'")
    print(f"export USE_CONDA='{config.get('software', {}).get('use_conda', 'true')}'")
    print(f"export CONDA_ENV='{config.get('software', {}).get('conda_env', 'genome_assembly')}'")
except Exception as e:
    print(f"echo 'Error parsing config: {e}' >&2", file=sys.stderr)
    sys.exit(1)
EOF
)
    else
        log WARN "Python3 not found, using default configuration"
        export PROJECT_NAME="${PROJECT_NAME:-TestGenome}"
        export THREADS="${THREADS:-4}"
        export HIFI_DATA="sequence_files/HiFi.fasta"
        export HIC_DATA="sequence_files"
    fi
    
    log INFO "Project name: ${PROJECT_NAME}"
    log INFO "Threads: ${THREADS}"
}

# ========== 环境准备 ==========
setup_environment() {
    local step_name="setup_environment"
    is_step_completed "$step_name" && log INFO "$step_name already completed, skipping" && return 0
    
    record_step "$step_name" "START"
    local start_time=$(start_timer)
    
    log INFO "Setting up environment..."
    
    # 环境已在init_pipeline中激活，这里只需验证
    if ! command -v hifiasm &> /dev/null; then
        log ERROR "hifiasm not found in activated environment"
        log ERROR "Current environment: ${CONDA_DEFAULT_ENV:-none}"
        log ERROR "Please run: bash setup_environment.sh"
        exit 1
    fi
    
    log INFO "Environment verified successfully"
    log INFO "Active environment: ${CONDA_DEFAULT_ENV}"
    log INFO "Hifiasm location: $(which hifiasm)"
    
    print_elapsed "$start_time" "$step_name"
    record_step "$step_name" "FINISH"
}

# ========== 步骤1: 数据质控预处理 ==========
run_data_qc() {
    local data_dir=$1
    local data_type=$2  # wgs, hic, rnaseq
    local step_name="qc_${data_type}"
    
    is_step_completed "$step_name" && log INFO "$step_name already completed, skipping" && return 0
    
    record_step "$step_name" "START"
    local start_time=$(start_timer)
    
    log INFO "Running quality control for ${data_type} data..."
    
    local qc_dir="${WORK_DIR}/01_qc/${data_type}"
    create_dir "$qc_dir"
    cd "$qc_dir"
    
    # 检查并列出所有fastq文件
    check_fastq_pairs "$data_dir" "$data_type"
    
    # 运行fastp进行质控
    local list_file="${data_dir}/${data_type}_file_list.txt"
    if [[ -f "$list_file" ]]; then
        # 提取R1/R2文件前缀
        grep "R1" "$list_file" | sed 's/R1.*//' | sort -u > "${data_type}_prefixes.txt"
        
        while read -r prefix; do
            local r1="${prefix}R1.fastq.gz"
            local r2="${prefix}R2.fastq.gz"
            local base=$(basename "$prefix")
            
            if [[ -f "$r1" ]] && [[ -f "$r2" ]]; then
                log INFO "Processing: $base"
                fastp -i "$r1" -I "$r2" \
                    -o "${base}clean_R1.fastq.gz" \
                    -O "${base}clean_R2.fastq.gz" \
                    -w "${THREADS}" \
                    -l 145 \
                    -h "${base}report.html" \
                    -j "${base}report.json" \
                    2>&1 | tee -a "${LOG_FILE}"
            fi
        done < "${data_type}_prefixes.txt"
    fi
    
    cd "$WORK_DIR"
    print_elapsed "$start_time" "$step_name"
    record_step "$step_name" "FINISH"
}

# ========== 步骤2: 基因组调查 ==========
run_genome_survey() {
    local step_name="genome_survey"
    is_step_completed "$step_name" && log INFO "$step_name already completed, skipping" && return 0
    
    record_step "$step_name" "START"
    local start_time=$(start_timer)
    
    log INFO "Running genome survey..."
    
    local survey_dir="${WORK_DIR}/02_genome_survey"
    create_dir "$survey_dir"
    cd "$survey_dir"
    
    # 使用WGS数据或HiFi数据进行基因组调查
    local input_data="${WGS_DATA:-${HIFI_DATA}}"
    
    # Jellyfish k-mer计数
    log INFO "Counting k-mers with Jellyfish..."
    jellyfish count -C -m 19 -s 1000000000 -t "${THREADS}" \
        -o genome.jf \
        <(find "${input_data}" -name "*.fastq.gz" -exec zcat {} \;) \
        2>&1 | tee -a "${LOG_FILE}"
    
    # 生成k-mer频率直方图
    jellyfish histo -t "${THREADS}" genome.jf > genome.histo
    
    # GenomeScope2分析
    if command -v genomescope2 &> /dev/null; then
        log INFO "Running GenomeScope2..."
        genomescope2 -i genome.histo -o ./ -k 19 -p 2 -n ${PROJECT_NAME} \
            2>&1 | tee -a "${LOG_FILE}"
    else
        log WARN "GenomeScope2 not found, skipping genome size estimation"
    fi
    
    cd "$WORK_DIR"
    print_elapsed "$start_time" "$step_name"
    record_step "$step_name" "FINISH"
}

# ========== 步骤3: Kraken2污染检测 ==========
run_kraken2() {
    local step_name="kraken2_contamination"
    is_step_completed "$step_name" && log INFO "$step_name already completed, skipping" && return 0
    
    record_step "$step_name" "START"
    local start_time=$(start_timer)
    
    log INFO "Running Kraken2 contamination detection..."
    
    local kraken_dir="${WORK_DIR}/03_kraken"
    create_dir "$kraken_dir"
    cd "$kraken_dir"
    
    # 验证HiFi数据
    validate_fasta "${HIFI_DATA}" || error_exit "Invalid HiFi data" "$step_name"
    
    # 运行Kraken2
    kraken2 --db "${KRAKEN_DB}" \
        --threads "${THREADS}" \
        --output "${PROJECT_NAME}_kraken.txt" \
        --report "${PROJECT_NAME}_report.txt" \
        --use-names \
        --confidence 0.5 \
        --minimum-hit-groups 300 \
        "${HIFI_DATA}" \
        2>&1 | tee -a "${LOG_FILE}"
    
    # 过滤保留目标物种序列
    if [[ -f "${MAMMAL_TAXIDS:-}" ]]; then
        log INFO "Filtering sequences by taxonomy..."
        awk -F "\t" '{if(NR==FNR){arr[$1]=0}else{if($1=="U"){print $2}else{match($3,/.*\(taxid (.*)\).*/,id);if(id[1] in arr){print $2}}}}' \
            "${MAMMAL_TAXIDS}" "${PROJECT_NAME}_kraken.txt" > keep_ids.txt
        
        # 使用seqkit提取序列
        seqkit grep -f keep_ids.txt "${HIFI_DATA}" > "${PROJECT_NAME}_filtered.fasta"
        
        # 更新HiFi数据路径
        export HIFI_DATA="${kraken_dir}/${PROJECT_NAME}_filtered.fasta"
        log INFO "Filtered HiFi data: ${HIFI_DATA}"
    fi
    
    cd "$WORK_DIR"
    print_elapsed "$start_time" "$step_name"
    record_step "$step_name" "FINISH"
}

# ========== 步骤4: Hifiasm基因组组装 ==========
run_hifiasm() {
    local step_name="hifiasm_assembly"
    is_step_completed "$step_name" && log INFO "$step_name already completed, skipping" && return 0
    
    record_step "$step_name" "START"
    local start_time=$(start_timer)
    
    log INFO "Running Hifiasm assembly..."
    
    local assembly_dir="${WORK_DIR}/04_hifiasm"
    create_dir "$assembly_dir"
    cd "$assembly_dir"
    
    # 转换为绝对路径
    local hifi_abs_path
    if [[ "${HIFI_DATA}" = /* ]]; then
        hifi_abs_path="${HIFI_DATA}"
    else
        hifi_abs_path="${WORK_DIR}/${HIFI_DATA}"
    fi
    
    # 验证HiFi数据是否存在
    if [[ ! -f "$hifi_abs_path" ]]; then
        log ERROR "HiFi data file not found: $hifi_abs_path"
        error_exit "HiFi data not found" "$step_name"
    fi
    
    # 构建Hi-C文件列表
    local hic1_files=""
    local hic2_files=""
    
    if [[ -d "${HIC_DATA}" ]]; then
        local hic_abs_path
        if [[ "${HIC_DATA}" = /* ]]; then
            hic_abs_path="${HIC_DATA}"
        else
            hic_abs_path="${WORK_DIR}/${HIC_DATA}"
        fi
        hic1_files=$(find "$hic_abs_path" -name "*R1*.fastq.gz" -o -name "*R1*.fastq" 2>/dev/null | tr '\n' ',' | sed 's/,$//')
        hic2_files=$(find "$hic_abs_path" -name "*R2*.fastq.gz" -o -name "*R2*.fastq" 2>/dev/null | tr '\n' ',' | sed 's/,$//')
    fi
    
    # 运行Hifiasm
    if [[ -n "$hic1_files" ]] && [[ -n "$hic2_files" ]]; then
        log INFO "Running Hifiasm with Hi-C data..."
        hifiasm -o "${PROJECT_NAME}" \
            -t 1 \
            --primary -l 0 \
            -m 1000 \
            --h1 "$hic1_files" \
            --h2 "$hic2_files" \
            "$hifi_abs_path" \
            2>&1 | tee -a "${LOG_FILE:-/dev/null}"
    else
        log INFO "Running Hifiasm without Hi-C data..."
        hifiasm -o "${PROJECT_NAME}" \
            -t 1 \
            --primary -l 0 \
            -m 1000 \
            "$hifi_abs_path" \
            2>&1 | tee -a "${LOG_FILE:-/dev/null}"
    fi
    
    # 转换GFA到FASTA
    if [[ -f "${PROJECT_NAME}.bp.p_ctg.gfa" ]]; then
        awk '/^S/{print ">"$2; print $3}' "${PROJECT_NAME}.bp.p_ctg.gfa" > "${PROJECT_NAME}.contigs.fa"
        export CONTIG_FILE="${assembly_dir}/${PROJECT_NAME}.contigs.fa"
    elif [[ -f "${PROJECT_NAME}.hic.p_ctg.gfa" ]]; then
        awk '/^S/{print ">"$2; print $3}' "${PROJECT_NAME}.hic.p_ctg.gfa" > "${PROJECT_NAME}.contigs.fa"
        export CONTIG_FILE="${assembly_dir}/${PROJECT_NAME}.contigs.fa"
    else
        log ERROR "No GFA output found from Hifiasm"
        error_exit "Hifiasm failed to produce output" "$step_name"
    fi
    
    log INFO "Contig file: ${CONTIG_FILE}"
    
    cd "$WORK_DIR"
    print_elapsed "$start_time" "$step_name"
    record_step "$step_name" "FINISH"
}

# ========== 步骤5: BUSCO质量评估 ==========
run_busco() {
    local step_name="busco_assessment"
    is_step_completed "$step_name" && log INFO "$step_name already completed, skipping" && return 0
    
    record_step "$step_name" "START"
    local start_time=$(start_timer)
    
    log INFO "Running BUSCO assessment..."
    
    local busco_dir="${WORK_DIR}/05_busco"
    create_dir "$busco_dir"
    cd "$busco_dir"
    
    # 运行BUSCO
    busco -i "${CONTIG_FILE}" \
        -o "${PROJECT_NAME}" \
        -l "${BUSCO_DB}" \
        -c "${THREADS}" \
        -m genome \
        --offline \
        2>&1 | tee -a "${LOG_FILE}"
    
    # 显示结果摘要
    if [[ -f "${PROJECT_NAME}/short_summary.txt" ]]; then
        log INFO "BUSCO Summary:"
        cat "${PROJECT_NAME}/short_summary.txt" | tee -a "${LOG_FILE}"
    fi
    
    cd "$WORK_DIR"
    print_elapsed "$start_time" "$step_name"
    record_step "$step_name" "FINISH"
}

# ========== 步骤6: QUAST评估 ==========
run_quast() {
    local step_name="quast_assessment"
    is_step_completed "$step_name" && log INFO "$step_name already completed, skipping" && return 0
    
    record_step "$step_name" "START"
    local start_time=$(start_timer)
    
    log INFO "Running QUAST assessment..."
    
    local quast_dir="${WORK_DIR}/06_quast"
    create_dir "$quast_dir"
    
    # 运行QUAST
    quast.py "${CONTIG_FILE}" \
        -o "$quast_dir" \
        -t "${THREADS}" \
        --eukaryote \
        --large \
        2>&1 | tee -a "${LOG_FILE}"
    
    cd "$WORK_DIR"
    print_elapsed "$start_time" "$step_name"
    record_step "$step_name" "FINISH"
}

# ========== 步骤7: Chromap Hi-C比对 ==========
run_chromap() {
    local step_name="chromap_alignment"
    is_step_completed "$step_name" && log INFO "$step_name already completed, skipping" && return 0
    
    record_step "$step_name" "START"
    local start_time=$(start_timer)
    
    log INFO "Running Chromap Hi-C alignment..."
    
    local chromap_dir="${WORK_DIR}/07_chromap"
    create_dir "$chromap_dir"
    cd "$chromap_dir"
    
    # 构建索引
    log INFO "Building Chromap index..."
    chromap -i -r "${CONTIG_FILE}" -o "${PROJECT_NAME}.index" \
        2>&1 | tee -a "${LOG_FILE}"
    
    # Hi-C比对
    local hic1_files=$(find "${HIC_DATA}" -name "*R1*.fastq.gz" | tr '\n' ',' | sed 's/,$//')
    local hic2_files=$(find "${HIC_DATA}" -name "*R2*.fastq.gz" | tr '\n' ',' | sed 's/,$//')
    
    log INFO "Aligning Hi-C reads..."
    chromap --preset hic \
        -r "${CONTIG_FILE}" \
        -x "${PROJECT_NAME}.index" \
        --remove-pcr-duplicates \
        -1 "$hic1_files" \
        -2 "$hic2_files" \
        --SAM \
        -o "${PROJECT_NAME}_hic.sam" \
        -t "${THREADS}" \
        2>&1 | tee -a "${LOG_FILE}"
    
    # 转换为BAM
    samtools view -@ "${THREADS}" -b "${PROJECT_NAME}_hic.sam" | \
        samtools sort -@ "${THREADS}" -o "${PROJECT_NAME}_hic.bam" -
    
    rm "${PROJECT_NAME}_hic.sam"
    
    export HIC_BAM="${chromap_dir}/${PROJECT_NAME}_hic.bam"
    log INFO "Hi-C BAM file: ${HIC_BAM}"
    
    cd "$WORK_DIR"
    print_elapsed "$start_time" "$step_name"
    record_step "$step_name" "FINISH"
}

# ========== 步骤8: YaHS染色体挂载 ==========
run_yahs() {
    local step_name="yahs_scaffolding"
    is_step_completed "$step_name" && log INFO "$step_name already completed, skipping" && return 0
    
    record_step "$step_name" "START"
    local start_time=$(start_timer)
    
    log INFO "Running YaHS scaffolding..."
    
    local yahs_dir="${WORK_DIR}/08_yahs"
    create_dir "$yahs_dir"
    cd "$yahs_dir"
    
    # 构建索引
    samtools faidx "${CONTIG_FILE}"
    
    # 运行YaHS
    yahs "${CONTIG_FILE}" "${HIC_BAM}" \
        -o "${PROJECT_NAME}" \
        2>&1 | tee -a "${LOG_FILE}"
    
    export SCAFFOLD_AGP="${yahs_dir}/${PROJECT_NAME}_scaffolds_final.agp"
    export HIC_BIN="${yahs_dir}/${PROJECT_NAME}.bin"
    
    log INFO "Scaffold AGP: ${SCAFFOLD_AGP}"
    log INFO "Hi-C bin file: ${HIC_BIN}"
    
    cd "$WORK_DIR"
    print_elapsed "$start_time" "$step_name"
    record_step "$step_name" "FINISH"
}

# ========== 步骤9: Juicer可视化 ==========
run_juicer() {
    local step_name="juicer_visualization"
    is_step_completed "$step_name" && log INFO "$step_name already completed, skipping" && return 0
    
    record_step "$step_name" "START"
    local start_time=$(start_timer)
    
    log INFO "Generating Juicer files..."
    
    local juicer_dir="${WORK_DIR}/09_juicer"
    create_dir "$juicer_dir"
    cd "$juicer_dir"
    
    # 运行juicer pre
    if command -v juicer &> /dev/null; then
        juicer pre -a -o "${PROJECT_NAME}" \
            "${HIC_BIN}" \
            "${SCAFFOLD_AGP}" \
            "${CONTIG_FILE}.fai" \
            > juicer_pre.log 2>&1
        
        # 生成.hic文件 (需要juicer_tools)
        if [[ -f "${JUICER_TOOLS_JAR:-}" ]]; then
            log INFO "Generating .hic file..."
            java -jar -Xmx500G "${JUICER_TOOLS_JAR}" pre \
                -j "${THREADS}" \
                --threads "${THREADS}" \
                "${PROJECT_NAME}.txt" \
                "${PROJECT_NAME}.hic" \
                <(grep "PRE_C_SIZE" juicer_pre.log | awk '{print $2" "$3}') \
                2>&1 | tee -a "${LOG_FILE}"
        fi
    else
        log WARN "Juicer not found, skipping visualization"
    fi
    
    cd "$WORK_DIR"
    print_elapsed "$start_time" "$step_name"
    record_step "$step_name" "FINISH"
}

# ========== 主流程 ==========
main() {
    local mode=${1:-"run"}
    
    case "$mode" in
        run)
            init_pipeline
            setup_environment
            
            # 读取步骤配置
            local RUN_QC=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE')).get('steps', {}).get('run_qc', False))" 2>/dev/null || echo "false")
            local RUN_GENOME_SURVEY=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE')).get('steps', {}).get('run_genome_survey', False))" 2>/dev/null || echo "false")
            local RUN_KRAKEN=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE')).get('steps', {}).get('run_kraken', False))" 2>/dev/null || echo "false")
            local RUN_HIFIASM=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE')).get('steps', {}).get('run_hifiasm', True))" 2>/dev/null || echo "true")
            local RUN_BUSCO=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE')).get('steps', {}).get('run_busco', False))" 2>/dev/null || echo "false")
            local RUN_QUAST=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE')).get('steps', {}).get('run_quast', True))" 2>/dev/null || echo "true")
            local RUN_CHROMAP=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE')).get('steps', {}).get('run_chromap', False))" 2>/dev/null || echo "false")
            local RUN_YAHS=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE')).get('steps', {}).get('run_yahs', False))" 2>/dev/null || echo "false")
            local RUN_JUICER=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE')).get('steps', {}).get('run_juicer', False))" 2>/dev/null || echo "false")
            
            # 数据质控
            if [[ "$RUN_QC" == "True" ]] && [[ -d "${WGS_DATA:-}" ]]; then
                run_data_qc "${WGS_DATA}" "wgs"
            fi
            if [[ "$RUN_QC" == "True" ]] && [[ -d "${HIC_DATA:-}" ]]; then
                run_data_qc "${HIC_DATA}" "hic"
            fi
            
            # 基因组调查与组装
            if [[ "$RUN_GENOME_SURVEY" == "True" ]]; then
                run_genome_survey
            fi
            if [[ "$RUN_KRAKEN" == "True" ]]; then
                run_kraken2
            fi
            if [[ "$RUN_HIFIASM" == "True" ]]; then
                run_hifiasm
            fi
            
            # 质量评估
            if [[ "$RUN_BUSCO" == "True" ]]; then
                run_busco
            fi
            if [[ "$RUN_QUAST" == "True" ]]; then
                run_quast
            fi
            
            # Hi-C挂载
            if [[ "$RUN_CHROMAP" == "True" ]]; then
                run_chromap
            fi
            if [[ "$RUN_YAHS" == "True" ]]; then
                run_yahs
            fi
            if [[ "$RUN_JUICER" == "True" ]]; then
                run_juicer
            fi
            
            log INFO "=========================================="
            log INFO "Pipeline completed successfully!"
            log INFO "=========================================="
            ;;
        
        resume)
            init_pipeline
            log INFO "Resuming from last completed step..."
            local last_step=$(get_last_completed_step)
            log INFO "Last completed step: ${last_step}"
            # 实现续跑逻辑...
            ;;
        
        clean)
            log INFO "Cleaning intermediate files..."
            rm -rf "${WORK_DIR}"/0*/{*.tmp,*.temp}
            log INFO "Cleanup completed"
            ;;
        
        *)
            echo "Usage: $0 {run|resume|clean}"
            echo "  run    - Run the complete pipeline"
            echo "  resume - Resume from last checkpoint"
            echo "  clean  - Clean intermediate files"
            exit 1
            ;;
    esac
}

# 执行主函数
main "$@"
