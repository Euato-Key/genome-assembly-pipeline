#!/bin/bash
set -o errexit

# author ： fei

readonly JUICER="juicer"
readonly MAMLIAN_TAXID="/hpcfile/home/zxf/references/Talpidae/Dpil/hifi_kraken/mammalia_taxids.txt"
JUICER_OUT=$2
LIFTOVER_AGP=$3
#NOW_TIME=`date +"%Y-%m-%d %H:%M:%S"`
#ont=""
#ont_RNASEQ=""
#pb_RNASEQ=""

KRAKEN_DB="/hpcfile/home/zxf/kraken_nt"
WORKDIR="$PWD"
WGS="/hpcfile/home/zxf/references/Talpidae/Utal/wgs_raw"
HIFI="/hpcfile/home/zxf/yanshu/Utal/hifi.fasta"  #must be fasta
RNASEQ="/hpcfile/home/zxf/references/Talpidae/Utal/rnaseq_raw"
HIC="/hpcfile/home/zxf/references/Talpidae/Utal/hic_raw" #dir of HIC
THREAD=12
DOCKER_NAME="assembly:0.6" 
NAME="UTAL" #project_NAME
BUSCO_DB="/hpcfile/home/zxf/busco_odb/laurasiatheria_odb10"
CONTIG=""
RERUN=""

usage() {
	echo -n -e \
	"  USAGE:  $(basename $0) [ main | rerun | post ]
	Example: $(basename $0) main    This will run main script.
	         $(basename $0) rerun   This will RERUN main script from last break point.
	         $(basename $0) post  JUICER_OUT  LIFTOVER_AGP
				     This will run post script, post script is used for 
				     final chromsome assembly.\n"
}

log() {
	local green='\033[0;32m'
    local red='\033[0;31m'
	local yellow='\033[1;33m'
    local nc='\033[0m'
    
    case $1 in
        INFO)  echo -e "$(date +"%Y-%m-%d %H:%M:%S") ${green}[ INFO]${nc} $2" ;;
        ERROR) echo -e "$(date +"%Y-%m-%d %H:%M:%S") ${red}[ERROR]${nc} $2" >&2 ;;
		WARN) echo -e "$(date +"%Y-%m-%d %H:%M:%S") ${yellow}[ WARN]${nc} $2" >&2 ;;
    esac | tee -a "$LOG_FILE"
}

slog() {
	echo -e "${FUNCNAME[0]}\tFINISH" >> ${WORKDIR}/log/status.log
}

log_check() {
	
	if [[ -f "${WORKDIR}/log/assembly.log" ]];then
		log INFO "Log info found! Make sure you are rerunning script!"
	else
		mkdir -p ${WORKDIR}/log
		declare -g LOG_FILE="${WORKDIR}/log/assembly.log"
		log INFO "Log dir created successfully!"
	fi	
}

status() {
	if [[ -f "${WORKDIR}/log/status.log" ]];then
		log WARN "Log info found! Make sure you are rerunning script!"
		RERUN="TRUE"
	else
		touch  ${WORKDIR}/log/status.log
		log INFO "Log created successfully!"
	fi
}

#run_data_pre [dir] [type(wgs,hic,rnaseq)]

run_data_pre() { 
	log INFO "${FUNCNAME[0]} starts running..."
	mkdir -p fastp && cd fastp

	local dir="${1}"
	local type="$2"

	ls -v ${dir}/*.fastq.gz > ${type}_need_fastp_lists.tmp
	
	# all fastq files must be XXXX.R1.fastq.gz XXXX.R2.fastq.gz, other files must delete
	
	for file in `cat ${type}_need_fastp_lists.tmp`;do
		if [[ $file =~ ^.*fastq\.gz$ ]];then
			echo "$file"
		else
			echo "${file} may not be a fastq file!"
			log ERROR "${file} may not be a fastq file!"
			exit 1
		fi	
	done

	## suffix must be R1/R2.fastq.gz
	for a in `cat ${type}_need_fastp_lists.tmp`;do
		match_dir=${a%/*}
		match_name=`basename $a .fastq.gz`
		match_name=${match_name%R*}
		echo "${match_dir}/${match_name}" >> ${type}_need_fastp_lists.tmp2	
	done

	sort -n -u ${type}_need_fastp_lists.tmp2 > ${type}_need_fastp_lists.tmp3
	
	NFL_line=`cat "${type}_need_fastp_lists.tmp3" | wc -l `
	#echo $NFL_line
	if [[ $((NFL_line%2)) == 0 ]];then
		log ERROR "Please check $a! It seem like having wrong number of fastq files."
		exit 1
	fi
	
	log INFO "${FUNCNAME[0]} ${type} finished."
	slog
	cd ${WORKDIR}
}

run_fastp2() { #run_fastp2 wgs
	log INFO "${FUNCNAME[0]} starts running..."
	#mkdir -p fastp && cd fastp
	[[ -d ${WORKDIR}/fastp ]] && cd ${WORKDIR}/fastp || (log ERROR "you need run data pre first!"; exit 1 )

	local type="$1"
	[[ -f ${WORKDIR}/fastp/${type}_need_fastp_lists.tmp3 ]] && mkdir -p ${type} || (log ERROR "you need run data pre first!";exit 1)

	#$((THREAD/4))
	docker run --rm -v $HOME:$HOME -w $PWD -u $(id -u):$(id -g) ${DOCKER_NAME} bash -c \
	"parallel --joblog ./parallel_${type}.log \
			 -j 1 \
			 -a ${type}_need_fastp_lists.tmp3 \
			 'cd wgs;fastp -w 16 -l 145 -i {}R1.fastq.gz -I {}R2.fastq.gz -o {/}clean_R1.fastq.gz -O {/}clean_R2.fastq.gz'"
			 
	if [[ $? != 0 ]];then
		log ERROR "run fastp of WGS error, check please!"
		#echo "${NOW_TIME} run fastp WGS error, check please!"
		exit 1
	fi

	log INFO "${FUNCNAME[0]} ${type} finish."
	slog
	cd ${WORKDIR}
}

#run_genome_survey type dir
run_genome_survey() {
	##genome survey
	log INFO "${FUNCNAME[0]} start running..."
	local type=$1
	local dir=${2:-"${WORKDIR}/fastp/${type}"}
	local genomescope2="genomescope.R"
	mkdir -p genome_survey && cd genome_survey

	#pigz -d -p 20 *.fq.gz
	docker run -u $(id -u):$(id -g) --rm -v $HOME:$HOME -w $PWD ${DOCKER_NAME} bash -c \
	"jellyfish count -C -m 19 -s 1000000000 -t ${THREAD} -o result.jf <(pigz -d -p ${THREAD} -c ${dir}/*.fastq.gz) && jellyfish histo -t $THREAD result.jf > result.histo && $genomescope2 -i result.histo -o ./ -k 19 -p 2 -n gs"

	if [[ $? != 0 ]];then
		log ERROR "${FUNCNAME[0]} error, check please!"
		#echo "${NOW_TIME} genome survey error, check please!"
		exit 1
	fi
	
	log INFO "${FUNCNAME[0]} finish."
	slog
	cd ${WORKDIR}
}

run_kraken2() {
	log INFO "${FUNCNAME[0]} start running..."
	local hifi=${1:-${HIFI}}
	mkdir -p kraken && cd kraken
	local confidence=0.5 #ngs seq set 0.4
	local mini_hit=300   #ngs seq set 3-5
	docker run -u $(id -u):$(id -g) -v $HOME:$HOME -w $PWD ${DOCKER_NAME} bash -c \
	"kraken2 -db $KRAKEN_DB --threads ${THREAD} \
		--output result_${NAME}.txt \
		--report report_${NAME}.txt \
		--use-names \
		--confidence ${confidence} \
		--minimum-hit-groups ${mini_hit} \
		--paired \
		$hifi"
		
	if [[ $? != 0 ]];then
		log ERROR "${FUNCNAME[0]} error,check please!"
		exit 1
	fi
	#MAMLIAN_TAXID="/hpcfile/home/zxf/references/Talpidae/Dpil/HIFI_kraken/mammalia_taxids.txt"

	awk -F "\t" '{if(NR==FNR){arr[$1]=0}else{if($1=="U"){print $2}else{match($3,/.*\(taxid (.*)\).*/,id);if(id[1] in arr){print $2}}}}' ${MAMLIAN_TAXID} result_${NAME}.txt > keep_lists.txt
	seqkit grep -f keep_lists.txt ${HIFI} > HIFI_kraken.fa

	[[ -f "${WORKDIR}/kraken/HIFI_kraken.fa" ]] && declare -g HIFI="${WORKDIR}/kraken/HIFI_kraken.fa"
	
	log INFO "${FUNCNAME[0]} finish."
	slog
	cd ${WORKDIR}
}

run_hifiasm() {
	log INFO "${FUNCNAME[0]} start running..."
	mkdir -p hifiasm && cd hifiasm
	local file="${1:+${workdir}/fastp/hic_need_fastp_lists.tmp3}"
	local hifi=${2:-$HIFI}
	local hic1=""
	local hic2=""
	for a in `cat ${file}`;do
		local fqname="${workdir}/fastq/hifi/`basename $a`"
		if [[ $hic1 = "" ]];then
			hic1="${fqname}R1.fastq.gz"
		else
			hic1="${hic1}','${fqname}R1.fastq.gz"
		fi
		
		if [[ $hic2 = "" ]];then
			hic2="${fqname}R2.fastq.gz"
		else
			hic2="${hic2}','${fqname}R2.fastq.gz"
		fi
	done
	
	docker run -u $(id -u):$(id -g) --rm -v $HOME:$HOME -w $PWD ${DOCKER_NAME} bash -c \
		"hifiasm -o ${NAME} -t ${THREAD} --primary -l 1 \
		${HIFI} \
		--h1 ${hic1} \
		--h2 ${hic2} && awk '/^S/{print ">"$2;print $3}' ${NAME}.hic.p_ctg.gfa > ${NAME}.hic.p_ctg.fa"
	
	declare -g CONTIG="${WORKDIR}/HIFIasm/${NAME}.HIC.p_ctg.fa"

	if [[ $? != 0 ]];then
		log ERROR "{FUNCNAME[0]} error, check please!"
		exit 1
	fi
	
	log INFO "${FUNCNAME[0]} finish."
	slog
	cd ${WORKDIR}
}

run_busco() {
	log INFO "${FUNCNAME[0]} start running..."
	
	local input=${1:-${WORKDIR}/hifiasm/${NAME}.hic.p_ctg.fa}
	mkdir busco && cd busco
	docker run -u $(id -u):$(id -g) --rm -v $HOME:$HOME -w $PWD ${DOCKER_NAME} bash -c \
	"busco -i ${input} \
		-o ${NAME} \
		--out_path ${NAME}_result \
		-l ${BUSCO_DB} \
		-c ${THREAD} \
		--offline \
		-m geno"
		
	log INFO "${FUNCNAME[0]} finish."
	slog
	cd ${WORKDIR}
	#echo "${NOW_TIME} BUSCO finished"
}

run_chromap() {
	log INFO "${FUNCNAME[0]} start running..."
	local type="hic"
	mkdir chromap && cd chromap
	local myhic1=""
	local myhic2=""
	for a in `cat ${workdir}/fastq/hic_need_fastp_lists.tmp3`;do
		myname="${workdir}/fastq/${type}/`basename $a`"
		if [[ $myhic1 = "" ]];then
			myhic1=${myname}R1.fastq.gz
		else
			myhic1=${myhic1}","${myname}R1.fastq.gz
		fi
		
		if [[ $myhic2 = "" ]];then
			myhic2=${myname}R2.fastq.gz
		else
			myhic2=${myhic2}","${myname}R2.fastq.gz
		fi
	done
	
	local sample=$NAME
	docker run -u $(id -u):$(id -g) --rm -v $HOME:$HOME -w $PWD ${DOCKER_NAME} bash -c \
	"chromap -i -r $CONTIG -o index_$sample && \
	chromap \
		--preset HIC \
		-r $CONTIG \
		-x index_$sample \
		--remove-pcr-duplicates \
		-1 $myhic1 \
		-2 $myhic2 \
		--SAM \
		-o ${sample}_aligned_HIC.sam \
		-t $THREAD \
		&& samtools view -b ${sample}_aligned_HIC.sam -- threads 8 -o ${sample}_aligned_HIC.bam && rm ${sample}_aligned_HIC.sam"
	declare -g ALIGNED_HIC="${PWD}/${sample}_aligned_HIC.bam"	
	#echo "${NOW_TIME} ${FUNCNAME[0]} finished"
	log INFO "${FUNCNAME[0]} finish."
	slog
	cd ${WORKDIR}
}

run_yahs() {
	log INFO "${FUNCNAME[0]} start running..."
	mkdir yahs && cd yahs
	#echo "${NOW_TIME} ${FUNCNAME[0]} start"
	#HIC="/hpcfile/home/zxf/references/Talpidae/Dpil/chromap/dpil_aligned_HIC.bam"
	docker run -u $(id -u):$(id -g) --rm -v $HOME:$HOME -w $PWD ${DOCKER_NAME} bash -c \
		"samtools faidx ${CONTIG} && yahs ${CONTIG} ${ALIGNED_HIC}"
	declare -g HIC_BIN="${PWD}/yahs.out.bin"
	declare -g AGP="${PWD}/yahs.out_scaffolds_final.agp"
	declare -g FAI="${CONTIG}.fai"
	#echo "${NOW_TIME} ${FUNCNAME[0]} finished"
	log INFO "${FUNCNAME[0]} finish."
	slog
	cd ${WORKDIR}
}

run_juicer() {
	log INFO "${FUNCNAME[0]} start running..."
	mkdir juicer && cd juicer
	#echo "${NOW_TIME} ${FUNCNAME[0]} start"
	 #
	jtool="/myscripts/soft/juicer_tools.3.0.0.jar" #must be juicer tools in juicer
	docker run -u $(id -u):$(id -g) --rm -v $HOME:$HOME -w $PWD ${DOCKER_NAME} bash -c \
	"$JUICER pre -a -o ${NAME} $HIC_BIN $AGP $FAI > juicer_pre.log && (java -jar -Xmx500G ${jtool} pre -j ${THREAD} --THREADs ${THREAD} ${NAME}.txt ${NAME}.HIC.part <(cat juicer_pre.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv ${NAME}.HIC.part ${NAME}.HIC)"
	log INFO "${FUNCNAME[0]} finish."
	slog
	cd ${WORKDIR}
}

run_juicer_post() {
	echo "${FUNCNAME[0]} start"
	docker run -u $(id -u):$(id -g) --rm -v $HOME:$HOME -w $PWD ${DOCKER_NAME} bash -c \
	"$JUICER post -o resout_${NAME} \
		$JUICER_OUT \
		$LIFTOVER_AGP \
		$CONTIG"
		
	if [[ $? != 0 ]];then
		echo "${NOW_TIME} ERROR: ${FUNCNAME[0]} exit, check log and RERUN."
		exit 1
	fi
	
	echo "${FUNCNAME[0]} finished"
}

main_run() {
	log_check
	status
	[[ $RERUN != "" ]] && log ERROR "There are logs in ${WORKDIR}/log, please delete log file. " && exit 2

	run_data_pre ${WGS} "wgs"
	run_fastp2 "wgs"

	run_data_pre ${HIC} "hic"
	run_fastp2 "hic"

	run_data_pre ${RNASEQ} "rnaseq"
	run_fastp2 "rnaseq"
	
	run_genome_survey
	run_kraken2
	run_hifiasm
	run_busco
	run_chromap
	run_yahs
	run_juicer
}

main_RERUN() {
	log_check
	status
	[[ $RERUN != "TRUE" ]] && log ERROR "There are no log file found, please run main!" && exit 2
	run=`tail -n1 ${WORKDIR}/log/status.log | awk -F "\t" '{print $1}'`
	
	case $run in
        run_fastp) 
			run_genome_survey
			run_kraken2
			run_hifiasm
			run_busco
			run_chromap
			run_yahs
			run_juicer
			;;
        run_genome_survey) 
			run_kraken2
			run_hifiasm
			run_busco
			run_chromap
			run_yahs
			run_juicer
			;;
		run_kraken2) 
			run_hifiasm
			run_busco
			run_chromap
			run_yahs
			run_juicer 
			;;
		run_HIFIasm) 
			run_busco
			run_chromap
			run_yahs
			run_juicer 
			;;
		run_busco) 
			run_chromap
			run_yahs
			run_juicer 
			;;
		run_chromap) 
			run_yahs
			run_juicer
			;;
		run_yahs) 
			run_juicer ;;
		run_juicer) log ERROR "Everything looks like fine. You don't need run script again! " && exit 3 ;;
		
		
		*) log ERROR "Something unexpected things happen, plesase check logs. " ;;
    esac 
}

main_post() {
	[[ $JUICER_OUT = "" || $LIFTOVER_AGP = "" ]] && (echo "need JUICER_OUT LIFTOVER_AGP"; exit 4 )
	run_juicer_post
}

main() {
	case $1 in
        main)  main_run ;;
        rerun) main_rerun ;;
		post) main_post ;;
		*) usage ;;
    esac 
}

#main "main"