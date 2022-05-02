# Functie:
#   Het script indexeert het referentie genoom en aligned de getrimde reads daartegen.
# Input:
#   -Eerst wordt gekeken of de imputh om help vraagt. Zo niet dan wordt het script
#    normaal uitgevoerd.
#   -CONSENSUS_PATH: De volledige path naar de consensus sequentie.
#   -PAIRED1 & PAIRED2: De volledige path van de getrimde read files.
#   -FASTA_GENOME: De volledige path naar het genoom dat moet worden geindexeerd.
#   -OUTPUT-DIR: Hier wordt het geindexeerde referentie genoom neer gezet.
# Output:
#   -Geindexeerde genoom.
#   -Sam file met info over de alignment van de reads.

# Show usage information:
if [ "${1}" == "--h" ] || [ "${1}" == "--help" ] || [ "${1}" == "-h" ] || [ "${1}" == "-help" ]
then
    echo ""
    echo "This script indexes a reference genome and alignes reads against it with bowtie2."
    echo "Part 4 of Dirks Next Gen Sequencing Pipeline."
    echo ""
    echo "Way of usage:"
    echo ""
    echo "    bash DeelOpdracht4Majops.sh <Argument 1> <Argument 2> etc."
    echo ""
    echo "    The arguments are as following:"
    echo "        1) Path to rimmed read-pair 1 fastq file as string."
    echo "        2) Path to rimmed read-pair 2 fastq file as string."
    echo "        3) Path to reference genome FOLDER (not the file itself!) as string."
    echo "        4) Name of reference genome FILE (without path!) as string."
    echo "        5) Path to output directory as string."
    echo ""
    echo "Example:"
    echo ""
    echo "    bash DeelOpdracht4Majops.sh \"/path/Trimmed1.fastq\" \"/path/Trimmed2.fastq\" \"/folder/ref_genome_folder/\" \"HomoSapiensGenome.fasta\" \"/path/OutputDir/>\""
    echo ""
    
    exit  
fi

# Input van de variabelen

PAIRED1="${1}"
PAIRED2="${2}"
CONSENSUS_PATH="${3}"
FASTA_GENOME="${4}"
OUTPUT_DIR="${5}"

# Verplaasting naar de map met het referentie genoom om errors te voorkomen.
cd "${CONSENSUS_PATH}"

# Indexeren van het referentie genoom.
bowtie2-build -q "${FASTA_GENOME}" "${OUTPUT_DIR}/indexed_genome_s1089166.fai"

# Weer terug naar de output map.
cd "${OUTPUT_DIR}"

# Alignen van de getrimde reads met het referentie genoom.
bowtie2 --quiet -x indexed_genome_s1089166.fai -1 "${PAIRED1}" -2 "${PAIRED2}" -S output.sam


