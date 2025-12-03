#!/bin/bash
# ==============================================================================
# 脚本名称: run_assembly.sh
# 功能: 基因组组装主流程
# ==============================================================================

set -e

# 获取脚本所在目录的绝对路径
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
CONFIG_FILE="${PROJECT_ROOT}/config/config.sh"

# 加载配置
if [ -f "$CONFIG_FILE" ]; then
    source "$CONFIG_FILE"
else
    echo "错误: 找不到配置文件 $CONFIG_FILE"
    exit 1
fi

# 激活 Conda 环境
# 注意：在脚本中激活 Conda 需要先 source conda.sh
if [ -f "${PROJECT_ROOT}/tools/miniconda3/bin/activate" ]; then
    source "${PROJECT_ROOT}/tools/miniconda3/bin/activate" "${CONDA_ENV_NAME}"
elif [ -f "$(conda info --base)/etc/profile.d/conda.sh" ]; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${CONDA_ENV_NAME}"
else
    echo "警告: 无法自动激活 Conda 环境，请确保您已手动激活 ${CONDA_ENV_NAME}"
fi

# 创建日志目录
mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/pipeline_$(date '+%Y%m%d_%H%M%S').log"

# 日志函数
log() {
    echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${LOG_FILE}"
}

log "========================================================"
log "开始运行基因组组装 Pipeline"
log "项目名称: ${PROJECT_NAME}"
log "工作目录: ${WORKDIR}"
log "========================================================"

# ------------------------------------------------------------------------------
# 1. 数据质控 (Fastp) - 针对 WGS 和 Hi-C 数据
# ------------------------------------------------------------------------------
run_fastp() {
    log "Step 1: 运行 Fastp 数据质控..."
    
    local OUT_DIR="${WORKDIR}/results/01_fastp"
    mkdir -p "${OUT_DIR}"

    # WGS 质控
    if [ -f "${WGS_R1}" ] && [ -f "${WGS_R2}" ]; then
        log "  正在处理 WGS 数据..."
        fastp -i "${WGS_R1}" -I "${WGS_R2}" \
              -o "${OUT_DIR}/wgs_clean_R1.fastq.gz" -O "${OUT_DIR}/wgs_clean_R2.fastq.gz" \
              -l "${FASTP_MIN_LEN}" \
              --thread "${THREADS}" \
              --html "${OUT_DIR}/wgs_fastp.html" \
              --json "${OUT_DIR}/wgs_fastp.json" 2>> "${LOG_FILE}"
    else
        log "  警告: 未找到 WGS 数据，跳过 WGS 质控。"
    fi

    # Hi-C 质控
    if [ -f "${HIC_R1}" ] && [ -f "${HIC_R2}" ]; then
        log "  正在处理 Hi-C 数据..."
        fastp -i "${HIC_R1}" -I "${HIC_R2}" \
              -o "${OUT_DIR}/hic_clean_R1.fastq.gz" -O "${OUT_DIR}/hic_clean_R2.fastq.gz" \
              -l "${FASTP_MIN_LEN}" \
              --thread "${THREADS}" \
              --html "${OUT_DIR}/hic_fastp.html" \
              --json "${OUT_DIR}/hic_fastp.json" 2>> "${LOG_FILE}"
    else
        log "  警告: 未找到 Hi-C 数据，跳过 Hi-C 质控。"
    fi
    
    log "Step 1 完成。"
}

# ------------------------------------------------------------------------------
# 2. 基因组调查 (Jellyfish + GenomeScope)
# ------------------------------------------------------------------------------
run_genome_survey() {
    log "Step 2: 运行基因组调查..."
    
    local IN_R1="${WORKDIR}/results/01_fastp/wgs_clean_R1.fastq.gz"
    local IN_R2="${WORKDIR}/results/01_fastp/wgs_clean_R2.fastq.gz"
    local OUT_DIR="${WORKDIR}/results/02_genome_survey"
    mkdir -p "${OUT_DIR}"

    if [ ! -f "${IN_R1}" ]; then
        log "  错误: 找不到质控后的 WGS 数据，跳过此步骤。"
        return
    fi

    # Jellyfish count
    log "  运行 Jellyfish count..."
    # 注意: <(zcat ...) 语法在某些 shell 中可能不兼容，这里使用 pigz -dc
    jellyfish count -C -m "${KMER_SIZE}" -s "${GENOME_SIZE_EST}" -t "${THREADS}" \
        -o "${OUT_DIR}/mer_counts.jf" \
        <(pigz -dc "${IN_R1}" "${IN_R2}") 2>> "${LOG_FILE}"

    # Jellyfish histo
    log "  运行 Jellyfish histo..."
    jellyfish histo -t "${THREADS}" "${OUT_DIR}/mer_counts.jf" > "${OUT_DIR}/mer_counts.histo"

    # GenomeScope
    log "  运行 GenomeScope..."
    genomescope2 -i "${OUT_DIR}/mer_counts.histo" -o "${OUT_DIR}" -k "${KMER_SIZE}" -p 2 2>> "${LOG_FILE}"

    log "Step 2 完成。"
}

# ------------------------------------------------------------------------------
# 3. 污染去除 (Kraken2) - 可选
# ------------------------------------------------------------------------------
run_kraken() {
    log "Step 3: 运行 Kraken2 污染检测..."
    
    if [ ! -d "${KRAKEN_DB}" ]; then
        log "  警告: Kraken2 数据库未找到 (${KRAKEN_DB})，跳过此步骤。"
        return
    fi

    local OUT_DIR="${WORKDIR}/results/03_kraken"
    mkdir -p "${OUT_DIR}"
    
    # 对 HiFi reads 进行检测
    kraken2 --db "${KRAKEN_DB}" --threads "${THREADS}" \
        --output "${OUT_DIR}/kraken_output.txt" \
        --report "${OUT_DIR}/kraken_report.txt" \
        --confidence "${KRAKEN_CONFIDENCE}" \
        "${HIFI_READS}" 2>> "${LOG_FILE}"
        
    # 注意：这里仅做检测和报告，自动过滤并提取 Clean Reads 逻辑较复杂，
    # 通常需要根据 TaxID 过滤。这里暂不自动执行序列提取，仅输出报告。
    log "  Kraken2 报告已生成，请检查报告以决定是否需要手动过滤。"
    log "Step 3 完成。"
}

# ------------------------------------------------------------------------------
# 4. 基因组组装 (Hifiasm)
# ------------------------------------------------------------------------------
run_assembly() {
    log "Step 4: 运行 Hifiasm 组装..."
    
    local OUT_DIR="${WORKDIR}/results/04_assembly"
    mkdir -p "${OUT_DIR}"
    
    cd "${OUT_DIR}"
    
    # Hifiasm 组装 (HiFi only 模式，如果有 Hi-C 可以加 --h1 --h2 参数进行 phasing)
    # 这里我们先做基础组装，后续用 Yahs 做 scaffolding
    log "  开始 Hifiasm..."
    
    # 直接运行 hifiasm
    hifiasm -o "${PROJECT_NAME}" -t "${THREADS}" "${HIFI_READS}" 2>> "${LOG_FILE}"
    
    # 转换 GFA 到 FASTA (Primary contigs)
    if [ -f "${PROJECT_NAME}.bp.p_ctg.gfa" ]; then
        awk '/^S/{print ">"$2;print $3}' "${PROJECT_NAME}.bp.p_ctg.gfa" > "${PROJECT_NAME}.p_ctg.fa"
        log "  组装完成，Primary Contigs: ${OUT_DIR}/${PROJECT_NAME}.p_ctg.fa"
    else
        log "  错误: 未生成 GFA 文件，组装可能失败，请检查日志。"
        exit 1
    fi
    
    log "Step 4 完成。"
}

# ------------------------------------------------------------------------------
# 5. 质量评估 (BUSCO) - Contig 级别
# ------------------------------------------------------------------------------
run_busco() {
    log "Step 5: 运行 BUSCO 评估 (Contig)..."
    
    if [ ! -d "${BUSCO_DB}" ]; then
        log "  警告: BUSCO 数据库未找到 (${BUSCO_DB})，跳过此步骤。"
        return
    fi

    local IN_FASTA="${WORKDIR}/results/04_assembly/${PROJECT_NAME}.p_ctg.fa"
    local OUT_DIR="${WORKDIR}/results/05_busco_contig"
    
    if [ ! -f "${IN_FASTA}" ]; then
        log "  错误: 找不到组装结果，跳过 BUSCO。"
        return
    fi

    busco -i "${IN_FASTA}" -l "${BUSCO_DB}" -o "${OUT_DIR}" -m genome -c "${THREADS}" --offline 2>> "${LOG_FILE}"
    
    log "Step 5 完成。"
}

# ------------------------------------------------------------------------------
# 6. Hi-C 比对 (Chromap)
# ------------------------------------------------------------------------------
run_chromap() {
    log "Step 6: 运行 Chromap Hi-C 比对..."
    
    local CONTIGS="${WORKDIR}/results/04_assembly/${PROJECT_NAME}.p_ctg.fa"
    local HIC_R1_CLEAN="${WORKDIR}/results/01_fastp/hic_clean_R1.fastq.gz"
    local HIC_R2_CLEAN="${WORKDIR}/results/01_fastp/hic_clean_R2.fastq.gz"
    local OUT_DIR="${WORKDIR}/results/06_chromap"
    mkdir -p "${OUT_DIR}"

    if [ ! -f "${HIC_R1_CLEAN}" ]; then
        log "  错误: 找不到 Clean Hi-C 数据，跳过此步骤。"
        return
    fi

    # 建立索引
    log "  建立索引..."
    chromap -i -r "${CONTIGS}" -o "${OUT_DIR}/contigs.index" 2>> "${LOG_FILE}"
    
    # 比对
    log "  开始比对..."
    chromap --preset "${CHROMAP_PRESET}" -r "${CONTIGS}" -x "${OUT_DIR}/contigs.index" \
        --remove-pcr-duplicates \
        -1 "${HIC_R1_CLEAN}" -2 "${HIC_R2_CLEAN}" \
        --SAM \
        -o "${OUT_DIR}/aligned.sam" -t "${THREADS}" 2>> "${LOG_FILE}"
        
    # 转 BAM 并排序
    log "  转换为 BAM..."
    samtools view -bS "${OUT_DIR}/aligned.sam" | samtools sort -@ "${THREADS}" -o "${OUT_DIR}/aligned.sorted.bam"
    rm "${OUT_DIR}/aligned.sam"
    
    log "Step 6 完成。"
}

# ------------------------------------------------------------------------------
# 7. Scaffolding (Yahs)
# ------------------------------------------------------------------------------
run_yahs() {
    log "Step 7: 运行 Yahs Scaffolding..."
    
    local CONTIGS="${WORKDIR}/results/04_assembly/${PROJECT_NAME}.p_ctg.fa"
    local BAM="${WORKDIR}/results/06_chromap/aligned.sorted.bam"
    local OUT_DIR="${WORKDIR}/results/07_yahs"
    mkdir -p "${OUT_DIR}"
    
    cd "${OUT_DIR}"
    
    # 建立 fai 索引
    samtools faidx "${CONTIGS}"
    
    # 运行 Yahs
    yahs "${CONTIGS}" "${BAM}" 2>> "${LOG_FILE}"
    
    log "  Scaffolding 完成，结果在 ${OUT_DIR}"
    log "Step 7 完成。"
}

# ------------------------------------------------------------------------------
# 8. 可视化 (Juicer)
# ------------------------------------------------------------------------------
run_juicer() {
    log "Step 8: 运行 Juicer Pre (生成 .hic 文件)..."
    
    local OUT_DIR="${WORKDIR}/results/08_juicer"
    mkdir -p "${OUT_DIR}"
    
    local BIN_FILE="${WORKDIR}/results/07_yahs/yahs.out.bin"
    local AGP_FILE="${WORKDIR}/results/07_yahs/yahs.out_scaffolds_final.agp"
    local FAI_FILE="${WORKDIR}/results/04_assembly/${PROJECT_NAME}.p_ctg.fa.fai"
    
    if [ ! -f "${JUICER_JAR}" ]; then
        log "  警告: Juicer Jar 未找到，跳过可视化步骤。"
        return
    fi

    # 使用 juicer pre
    # 注意: Yahs 输出的 bin 文件可以直接用于 juicer pre
    # 需要先生成 .txt 格式或者直接用 bin (取决于 juicer 版本，这里使用 yahs 推荐的流程)
    
    # Yahs 官方推荐流程:
    # juicer pre -a -o out_JBAT yahs.out.bin yahs.out_scaffolds_final.agp contigs.fa.fai >out_JBAT.log 2>&1
    # 然后用 java -jar juicer_tools pre ...
    
    # 这里简化处理，假设 juicer_tools 可用
    # 由于环境配置复杂，这里仅打印命令提示，防止因 Java 版本或内存问题导致整个流程崩溃
    log "  提示: 请手动运行以下命令生成 .hic 文件 (需要大量内存):"
    log "  java -Xmx${JUICER_JVM_MEM} -jar ${JUICER_JAR} pre ${BIN_FILE} ${AGP_FILE} ${FAI_FILE}"
    
    log "Step 8 (预处理) 完成。"
}

# ------------------------------------------------------------------------------
# 执行所有步骤
# ------------------------------------------------------------------------------
run_fastp
run_genome_survey
run_kraken
run_assembly
run_busco
run_chromap
run_yahs
run_juicer

log "========================================================"
log "所有流程执行完毕！"
log "日志文件: ${LOG_FILE}"
log "========================================================"
