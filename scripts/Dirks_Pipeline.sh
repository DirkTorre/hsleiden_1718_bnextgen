#!/bin/bash

# Show usage information:
if [ "${1}" == "--h" ] || [ "${1}" == "--help" ] || [ "${1}" == "-h" ] || [ "${1}" == "-help" ]
then
    echo ""
    echo "Takes arguments given by user in the Start script and executes the Next Gen Sequencing pipeline."
    echo ""
    echo "Way of usage:"
    echo ""
    echo "    This script must be executed by the starter script \"Start_Dirks_Pipeline.sh\""
    echo "    The script is executed in the following way:"
    echo "        bash Dirks_Pipeline <Argument 1> <Argument 2> etc..."
    echo ""
    echo "    The arguments are as following:"
    echo "        1)  Full path to directory of all scripts as string."
    echo "        2)  Full path to a exisiting output folder as string."
    echo "        3)  Full path to the first readfile as string."
    echo "        4)  Full path to the second readfile as string."
    echo "        5)  Full path to an output QC file"
    echo "        6)  A upper limit on the percentage bad bases in a read as a integer."
    echo "        7)  A lower limit on read length as a integer."
    echo "        8)  Full path to the reference genome folder as a string."
    echo "        9)  The reference genome file name as a string"
    echo "        10) A prefix for the output files."
    echo "        11) Full path to the current executing folder."
    echo ""
    echo "Example:"
    echo ""
    echo "    bash Dirks_Pipeline.sh \"/path/scripts/\" \"/path/output/\" \"/path/reads1.fastq\" \"/path/reads2.fastq\" \"/outputdir/QC.txt\" 15 20 \"/path/refgenomedir/\" \"refgenome.fasta\" \"run1\" \"/working_path/\""
    echo ""
    
    
    exit
fi

echo ""
echo "Pipeline owner: Dirk van der Torre - s1089166"
echo "The Next Gen Sequencing Pipeline has started." 
echo ""




# Real Bash script of the highway of the pipeline starts from here...

# Ophalen variabelen uit het start script.
SCRIPTS_DIR="${1}"
OUTPUT_DIR="${2}"
READSFILE1="${3}"
READSFILE2="${4}"
QCFILE="${5}"
MAX_BAD_PERCENTAGE=${6}
MIN_READ_LEN=${7}
REF_GENOME_FOLDER="${8}"
REF_GENOME_FILE_NAME="${9}"
PREFIX="${10}"
ORIGINAL_PATH="${11}"


# Quality control uitvoeren op paired-end reads.
echo "1/6: QC before trimming started"
python2 "${SCRIPTS_DIR}/DeelOpdracht1Majops.py" "${READSFILE1}" "${READSFILE2}" "${QCFILE}"
printf "1/6: QC before trimming done:\n    " >> "${OUTPUT_DIR}/RunTime.txt"

date >> "${OUTPUT_DIR}/RunTime.txt"


# Trimmen paired-end reads en verwijderen single reads.
echo "2/6: Trimming started"
python2 "${SCRIPTS_DIR}/DeelOpdracht2Majops.py" "${READSFILE1}" "${READSFILE2}" ${MAX_BAD_PERCENTAGE} ${MIN_READ_LEN} "${OUTPUT_DIR}/trimmed"
printf "2/6: Trimming done:\n    " >> "${OUTPUT_DIR}/RunTime.txt"
date >> "${OUTPUT_DIR}/RunTime.txt"


# Quality control uitvoeren op getrimde paired-end reads.
echo "3/6: QC after trimming started"
python2 "${SCRIPTS_DIR}/DeelOpdracht1Majops.py" "${OUTPUT_DIR}/trimmed_paired_1.fastq" "${OUTPUT_DIR}/trimmed_paired_2.fastq" "${QCFILE}"
printf "3/6: QC after trimming done:\n    " >> "${OUTPUT_DIR}/RunTime.txt"
date >> "${OUTPUT_DIR}/RunTime.txt"


# Indexeren van referentie genoom en alignen van getrimde reads daartegen.
echo "4/6: Indexing and alignment started"
bash "${SCRIPTS_DIR}/DeelOpdracht4Majops.sh" "${OUTPUT_DIR}/trimmed_paired_1.fastq" "${OUTPUT_DIR}/trimmed_paired_2.fastq" "${REF_GENOME_FOLDER}" "${REF_GENOME_FILE_NAME}" "${OUTPUT_DIR}"
printf "4/6: Indexing and alignment started:\n    " >> "${OUTPUT_DIR}/RunTime.txt"
date >> "${OUTPUT_DIR}/RunTime.txt"


# Maken variant lijst en consensus sequentie.
echo "5/6: Determining variant list and consensus seqence started"
bash "${SCRIPTS_DIR}/DeelOpdracht5Majops.sh" "${OUTPUT_DIR}" "${REF_GENOME_FOLDER}" "${REF_GENOME_FILE_NAME}"
printf "5/6: Determining variant list and consensus seqence done:\n    " >> "${OUTPUT_DIR}/RunTime.txt"
date >> "${OUTPUT_DIR}/RunTime.txt"


# Bepalen hoeveelheid inserties, deleties en soorten substituties.
echo "6/6: Counting insertions, deletions and substitutions started"
python2 "${SCRIPTS_DIR}/DeelOpdracht6Majops.py" "${OUTPUT_DIR}/my_var.vcf" "${OUTPUT_DIR}/s1089166_stats.txt"
printf "6/6: Counting insertions, deletions and substitutions done:\n    " >> "${OUTPUT_DIR}/RunTime.txt"
date >> "${OUTPUT_DIR}/RunTime.txt"

printf "eind tijd: " >> "${OUTPUT_DIR}/RunTime.txt\n    "
date >> "${OUTPUT_DIR}/RunTime.txt"

cd "${OUTPUT_DIR}"

for f in * ; do mv "$f" "${PREFIX}_$f" ; done

cd "${ORIGINAL_PATH}"

echo -e "\nAll generated data can be found at:\n ${OUTPUT_DIR}\n"

echo "Trimmed reads in:            ${PREFIX}_trimmed_paired_<1/2>.fastq"
echo "Quality control in:          ${PREFIX}_QC.txt"
echo "Indexed reference genome in: ${PREFIX}_indexed_genome_s1089166.fai.<number/rev.<number>>.bt2"
echo "SAM file in:                 ${PREFIX}_output.sam"
echo "BAM file in:                 ${PREFIX}_output_sorted.bam"
echo "Mpileup converted to bcf in: ${PREFIX}_raw.bcf"
echo "Consensus sequence in:       ${PREFIX}_consensus.fastq"
echo "Answers to assignment 6 in:  ${PREFIX}_s1089166_stats.txt"
echo "Run time in:                 ${PREFIX}_RunTime.txt"



echo ""
echo "The Next Gen Sequencing Pipeline has finished." 
echo ""
# 
# # Additional information:
# # =======================
# #
# # Remarks about the Majops Skeleton Script Bash Pipeline itself.
# # Description how it works.
# # Description which improvements can be done to improve the Majops Skeleton Script Bash Pipeline itself.
# #
