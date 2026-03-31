#!/bin/bash
set -e
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "${SCRIPT_DIR}/../config/config.sh"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ${CONDA_ENV_NAME}

mkdir -p ${LOG_DIR} ${OUT_BASE}
exec > >(tee -a "${LOG_DIR}/pipeline.log") 2>&1
log() { echo -e "\033[34m[$(date '+%H:%M:%S')]\033[0m $1"; }

# Step 1: Fastp
log "Step 1: 质控..."
mkdir -p ${FASTP_DIR}
fastp -i ${WGS_R1} -I ${WGS_R2} -o ${FASTP_DIR}/clean_1.fq.gz -O ${FASTP_DIR}/clean_2.fq.gz -t ${THREADS}

# Step 4: Hifiasm
log "Step 4: 组装..."
mkdir -p ${HIFIASM_DIR} && cd ${HIFIASM_DIR}
hifiasm -o ${PROJECT_NAME} -t ${THREADS} --h1 ${HIC_R1} --h2 ${HIC_R2} ${HIFI_READS}
awk '/^S/{print ">"$2;print $3}' ${PROJECT_NAME}.hic.p_ctg.gfa > ${PROJECT_NAME}.p_ctg.fa
samtools faidx ${PROJECT_NAME}.p_ctg.fa

# Step 5: BUSCO (AI 自动判断匹配 odb)
log "Step 5: BUSCO 评估 (Auto-Lineage模式)..."
mkdir -p ${BUSCO_DIR}
if [ "${BUSCO_LINEAGE}" == "auto" ]; then
    BUSCO_MODE="--auto-lineage-euk"
else
    BUSCO_MODE="-l ${BUSCO_LINEAGE}"
fi
busco -i ${HIFIASM_DIR}/${PROJECT_NAME}.p_ctg.fa -o busco_res --out_path ${BUSCO_DIR} \
      -m genome -c ${THREADS} ${BUSCO_MODE} --download_path ${BUSCO_DB}

# Step 6 & 7: Chromap & YAHS
log "Step 7: Scaffolding..."
mkdir -p ${CHROMAP_DIR}
chromap --preset hic -r ${HIFIASM_DIR}/${PROJECT_NAME}.p_ctg.fa -1 ${HIC_R1} -2 ${HIC_R2} \
        -o ${CHROMAP_DIR}/aligned.bam -t ${THREADS} --output-format bam
mkdir -p ${YAHS_DIR} && cd ${YAHS_DIR}
yahs ${HIFIASM_DIR}/${PROJECT_NAME}.p_ctg.fa ${CHROMAP_DIR}/aligned.bam

# Step 8: Juicer Post
log "Step 8: 可视化处理..."
mkdir -p ${JUICER_DIR}
juicer pre -a -o ${JUICER_DIR}/${PROJECT_NAME} \
    ${YAHS_DIR}/yahs.out.bin ${YAHS_DIR}/yahs.out_scaffolds_final.agp ${HIFIASM_DIR}/${PROJECT_NAME}.p_ctg.fa.fai > ${JUICER_DIR}/pre.log 2>&1
cd ${JUICER_DIR}
SIZES=$(grep PRE_C_SIZE pre.log | awk '{print $2" "$3}')
java -Xmx${JUICER_JVM_MEM} -jar ${JUICER_JAR} pre -j ${THREADS} ${PROJECT_NAME}.txt ${PROJECT_NAME}.hic <(echo "$SIZES")

log "✅ 流程结束！"
