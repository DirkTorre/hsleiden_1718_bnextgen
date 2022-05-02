#!/bin/bash 

# Additional information:
# =======================
#
# Information about the pipeline is given here.  
#

# Show usage information:
if [ "${1}" == "--h" ] || [ "${1}" == "--help" ] || [ "${1}" == "-h" ] || [ "${1}" == "-help" ]
then

    echo "" 
    echo "This pipeline automates variant detection in sequenced reads"
    echo "compared to a reference genome. Variants will be determinded."
    echo ""
    echo "The script can be executed from anywhere with:"
    echo "    Bash Start_Dirks_Pipeline.sh"
    echo ""
    echo "The script asks for a prefix for all output, including its folder."
    echo ""
    echo "A directory will be created in the same directory as the \"scripts\" directory."
    echo "It will be named: <prefix>_s_1089166_output_<YEAR>_<MONTH>_<DAY>_<HOUR>_<MINUTE>_<SECOND>"
    echo ""
    echo "The script containes internal variables which can be changed (comment out examples included):"
    echo "    READSFILE1=\"<string>\":  Full path to the first paired end read file."
    echo "    READSFILE2=\"<string>\":  Full path to the second paired end read file."
    echo "    MAX_BAD_PERCENTAGE=<integer>: Maximum percentage of bases under phred 20 score."
    echo "    MIN_READ_LEN=<integer>:   Mimimum read length that will be allowed to survive."
    echo "    REF_GENOME_FOLDER=\"<string>\":  Full path to the reference genome folder."
    echo "    REF_GENOME_FILE_NAME=\"<string>\":    Full name of the reference genome."
    echo ""
    echo "The script will give the path to the results directory."
    echo "This will contain:"
    echo "    \"<prefix>_consensus.fastq\":                   The fastq consensus sequence."
    echo "    \"<prefix>_indexed_genome_s1089166.fai.<etc>\": The indexed genome."
    echo "    \"<prefix>_my_var.bcf\":                        The bcf file."
    echo "    \"<prefix>_my_var.vcf\":                        The vcf file."
    echo "    \"<prefix>_output_sorted.bam\":                 The sorted bam file."
    echo "    \"<prefix>_output.sam\":                        The sam file."
    echo "    \"<prefix>_QC.txt\":                            Quality control file of reads before and after trimming."
    echo "    \"<prefix>_raw.bcf\":                           Mpileup converted to BCF."
    echo "    \"<prefix>_RunTime.txt\":                       States the start time and date of every script."
    echo "    \"<prefix>_s1089166_stats.txt\":                The file containing the amount of every mutation type."
    echo "    \"<prefix>_trimmed_paired_<1/2>.fastq\":        The trimmed read files."
    echo ""
    
    exit  
fi

# Checken op prefix.
while true
do
    read -p "Please enter a prefix that will be used for all output files: " PREFIX
    if [[ "${PREFIX}" =~ ^([A-Z]|[a-z]|[0-9]|_)+$ ]]
    then
        break
    fi
    echo "  Error: Only alphanumerical characters and underscores allowed."
done

ORIGINAL_PATH="$( pwd )"
SCRIPTS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$( echo ${SCRIPTS_DIR} | sed "s/\/scripts$/\//" )"
DATE=$(date +%Y_%m_%d_%H_%M_%S)
OUTPUT_DIR="${BASE_DIR}${PREFIX}_s1089166_${DATE}"

mkdir "${OUTPUT_DIR}"


MAX_BAD_PERCENTAGE=3
MIN_READ_LEN=60

# # laptop
# READSFILE1="/home/dirk/Documents/bnextgen/script/v8 deelopdracht 5/input/s_2_1_sequence.txt"
# READSFILE2="/home/dirk/Documents/bnextgen/script/v8 deelopdracht 5/input/s_2_2_sequence.txt"
# REF_GENOME_FOLDER="/home/dirk/Documents/bnextgen/script/assembly test/assembly"
# REF_GENOME_FILE_NAME="niet_infected_consensus.fasta"

# # server kleine files
# READSFILE1="/home/bnextgen/reads/testbestand1.fastq"
# READSFILE2="/home/bnextgen/reads/testbestand2.fastq"
# REF_GENOME_FOLDER="/home/bnextgen/refgenome"
# REF_GENOME_FILE_NAME="niet_infected_consensus.fasta"

# server grote files
READSFILE1="/home/bnextgen/reads/s_2_1_sequence.txt"
READSFILE2="/home/bnextgen/reads/s_2_2_sequence.txt"
REF_GENOME_FOLDER="/home/bnextgen/refgenome"
REF_GENOME_FILE_NAME="niet_infected_consensus.fasta"

# # desktop kleine reads files
# READSFILE1="/home/dirk/Documents/bioinfo/jaar4/bnextgen/opdracht/script/v8 deelopdracht 5/input/s_2_1_sequence.txt"
# READSFILE2="/home/dirk/Documents/bioinfo/jaar4/bnextgen/opdracht/script/v8 deelopdracht 5/input/s_2_2_sequence.txt"
# REF_GENOME_FOLDER="/home/dirk/Documents/bioinfo/jaar4/bnextgen/opdracht/script/assembly"
# REF_GENOME_FILE_NAME="niet_infected_consensus.fasta"

# # desktop grote reads files
# READSFILE1="/home/dirk/Documents/bioinfo/jaar4/bnextgen/opdracht/originele data/s_2_1_sequence_full.txt"
# READSFILE2="/home/dirk/Documents/bioinfo/jaar4/bnextgen/opdracht/originele data/s_2_2_sequence_full.txt"
# REF_GENOME_FOLDER="/home/dirk/Documents/bioinfo/jaar4/bnextgen/opdracht/script/assembly"
# REF_GENOME_FILE_NAME="niet_infected_consensus.fasta"


# Check presence of the script files.
onePipelineItemMissing=false

# All the necessary files and directories of the pipeline must be present.
# Their presences are checked.

for pipelineFile in \
        "DeelOpdracht1Majops.py" \
        "DeelOpdracht2Majops.py" \
        "DeelOpdracht4Majops.sh" \
        "DeelOpdracht5Majops.sh" \
        "DeelOpdracht6Majops.py" \
        "Dirks_Pipeline.sh"
# All the pipeline files and possible input files are listed here. 
do
	if ! [ -f "${SCRIPTS_DIR}/${pipelineFile}" ]
	then
		echo ${pipelineFile}" does not exist." 
		onePipelineItemMissing=true
	fi
done
#  
# Exit, if one or more pipeline items are missing.
if ${onePipelineItemMissing}
then
	echo "The NextGenSequencing Bash Pipeline is aborted."
	exit 	
fi

# The right command line arguments are generated out of the answers 
# of the questions here. 


# The Majops Skeleton Script Bash Pipeline starts.
printf "Script start:\n    " > "${OUTPUT_DIR}/RunTime.txt"
date >> "${OUTPUT_DIR}/RunTime.txt"


bash "${SCRIPTS_DIR}/Dirks_Pipeline.sh" "${SCRIPTS_DIR}" "${OUTPUT_DIR}" "${READSFILE1}" "${READSFILE2}" "${OUTPUT_DIR}/QC.txt" ${MAX_BAD_PERCENTAGE} ${MIN_READ_LEN} "${REF_GENOME_FOLDER}" "${REF_GENOME_FILE_NAME}" "${PREFIX}" "${ORIGINAL_PATH}"

printf "Script end:\n    " >> "${OUTPUT_DIR}/${PREFIX}_RunTime.txt"
date >> "${OUTPUT_DIR}/${PREFIX}_RunTime.txt"