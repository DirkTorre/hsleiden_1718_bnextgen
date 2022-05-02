# Functie:
#   Variant lijst en consensus sequentie maken.
# Input:
#   -OUTPUT_DIR: Volledige path voor de output bestanden.
#   -REF_GENOME_FOLDER: Volledige path naar het referentie genoom.
#   -REF_GENOME_FILE_NAME: De Naam van het referentie genoom.
# Output:
#   -my_var.vcf: VCF bestand (variant lijst) met alle varianten die zijn gevonden.
#   -my_var.bcf: De binary variant van de VCF.
#   -consensus.fastq: De sconsesus sequentie.
#   -pileup bestand.

# Show usage information:
if [ "${1}" == "--h" ] || [ "${1}" == "--help" ] || [ "${1}" == "-h" ] || [ "${1}" == "-help" ]
then
    echo ""
    echo "Creates variant list and consensus sequence with samtools."
    echo "Part 5 of Dirks Next Gen Sequencing Pipeline."
    echo ""
    echo "Way of usage:"
    echo ""
    echo "    bash DeelOpdracht5Majops.sh <Argument 1> <Argument 2> etc."
    echo ""
    echo "    The arguments are as following:"
    echo "        1) Path to output directory as string."
    echo "        2) Path to reference genome FOLDER (not the file itself!) as string."
    echo "        3) Reference genome file name as string."
    echo ""
    echo "Example:"
    echo ""
    echo "    bash DeelOpdracht5Majops.sh \"/path/OutputDir/\" \"/path/ref_genome_folder/\" \"HomoSapiensGenome.fasta\" "
    echo ""
    
    exit  
fi

OUTPUT_DIR="${1}"
REF_GENOME_FOLDER="${2}"
REF_GENOME_FILE_NAME="${3}"


### Converteren alignment.
# SAM naar BAM converteren.
samtools view -Sb "${OUTPUT_DIR}/output.sam" > "${OUTPUT_DIR}/output.bam"

# BAM sorteren (.bam gaat vanzelf in de output file).
samtools sort "${OUTPUT_DIR}/output.bam" "${OUTPUT_DIR}/output_sorted"

# Niet gesorteerde bam verwijderen.
rm "${OUTPUT_DIR}/output.bam"


### SNP calling.
# indexeren genoom voor pileup
samtools faidx "${REF_GENOME_FOLDER}/${REF_GENOME_FILE_NAME}"

# Pileup maken en omzetten naar bcf.
samtools mpileup -g -f "${REF_GENOME_FOLDER}/${REF_GENOME_FILE_NAME}" "${OUTPUT_DIR}/output_sorted.bam" > "${OUTPUT_DIR}/raw.bcf"

# Filteren van varianten.
bcftools view -bvcg "${OUTPUT_DIR}/raw.bcf" > "${OUTPUT_DIR}/my_var.bcf"

# BCF naar VCF omzetten.
bcftools view "${OUTPUT_DIR}/my_var.bcf" > "${OUTPUT_DIR}/my_var.vcf"


### Consensus sequentie maken.
# Om vcfutils.pl te vinden: locate samtools/vcfutils.pl
cat "${OUTPUT_DIR}/my_var.vcf" | /usr/share/samtools/vcfutils.pl vcf2fq > "${OUTPUT_DIR}/consensus.fastq"