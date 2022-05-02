import subprocess
import sys

"""
Dit script trimt paired-end reads en verwijderd single reads.
"""

def showUsageInformation():
    """
    Functie: Geeft info over het script en het gebruik er van.
    """
    print "" 
    print "This script trims fastq reads."
    print "Part 2 of Dirks Next Gen Sequencing Pipeline."
    print ""
    print "Way of usage:"
    print ""
    print "    python2 DeelOpdracht2Majops.py  <Argument 1> <Argument 2> etc..."
    print ""
    print "    The arguments are as following:"
    print "        1) Paired reads 1 file (with path)."
    print "        2) Paired reads 2 file (with path)."
    print "        3) Max percentage bad nucleotides allowed for read (as integer)."
    print "        4) Minimum acceptabele read length (as integer)."
    print ""
    print "Example:"
    print ""
    print "    python2 DeelOpdracht2Majops.py \"/path/reads1.fastq/\" \"/path/reads2.fastq\" 20 20"
    print ""
    print "Method used:"
    print ""
    print "    1) Trimming read by quality (substring determined by first and"
    print "       last base with quality above 19)."
    print "    2) Check if length of trimmed read is the same or higher as threshold."
    print "    3) Check if GC precentage of trimmed read is same or higher as threshold."
    print "    4) If both reads from the read pair come trough these checks the"
    print "       trimmed pair is written to the file"
    print ""
    
    sys.exit()


def trim(fastq1, fastq2, bestand_trimmed, max_percentage_slecht, min_read_len):
    """
    Functie: Trimt reads en schrijft ze naar files.
    Input:   -fastq1: volledige path van het eerste paired-end reads fastq bestand.
             -fastq2: volledige path van het tweede paired-end reads fastq bestand.
             -bestand_trimmed: path en prefix voor de trimmed reads.
             -max_percentage_slecht: maximale percentage slechte basen (als int).
             -min_read_len: minimale acceptabele lengte voor een read.
    Output:  -2 files met getrimde reads op opgegeven locatie met opgegeven prefix.
    Methode: 1) Reads worden getrimd (technisch gezien wordt de qualiteits string getrimd).
             2) De lengte van de resulterende read gaat door een check.
             3) Als de lengte goed is wordt gekeken of de read niet teveel
                lage kwaliteit nucleotiden heeft.
             4) Als het read pair door de checks komen dan worden ze weg geschreven.
                >>Beiden reads moeten dus door de checks komen!!!!!
    """
    
    # Maken en openen trimmed files
    paired_1_trimmed = open(bestand_trimmed+"_paired_1.fastq", 'w')
    paired_2_trimmed = open(bestand_trimmed+"_paired_2.fastq", 'w')
    
    # lezen read files
    reads1_file = open(str(fastq1), 'r')
    reads1 = reads1_file.readlines()
    reads1_file.close()
    reads2_file = open(str(fastq2), 'r')
    reads2 = reads2_file.readlines()
    reads2_file.close()
    
    # Door files loopen.
    for x in range(0,len(reads1)-3,4):
        # Reads per paar en per [ seq header, seq, qual header, qual ] ophalen.
        read1 = [ reads1[x].strip(), reads1[x+1].strip(), reads1[x+2].strip(),
                 reads1[x+3].strip() ]
        read2 = [ reads2[x].strip(), reads2[x+1].strip(), reads2[x+2].strip(),
                 reads2[x+3].strip() ]
        
        # Trimmen van de qualiteits scores.
        getrimde_indexen = [ list(trimEinden(reads1[x+3].strip())),
                            list(trimEinden(reads1[x+3].strip())) ]
        
        # Controleren of lengte na trimmen lang genoeg is.
        check_len = not (len(getrimde_indexen[0]) < min_read_len or 
                len(getrimde_indexen[1]) < min_read_len)
        if check_len:
            # Nu apart checken voor het geode nucleotide gehalte.
            # Dit kan ook nog aan de "if check_len" worden geplakt met een and,
            # maar dan moet je ook GC gehaltes bepalen van sequenties die
            # uberhaupt niet door de keuring gaan.
            check_geh = [ goedNucleotidenGehalte(
                reads1[3], getrimde_indexen[0], max_percentage_slecht) and 
            goedNucleotidenGehalte(reads2[3], getrimde_indexen[1], max_percentage_slecht) ]
            
            if check_geh:
                # Goede read pairs weg schrijven.
                begin1 = min(getrimde_indexen[0])
                einde1 = max(getrimde_indexen[0])
                begin2 = min(getrimde_indexen[1])
                einde2 = max(getrimde_indexen[1])
                read1_write = read1[0]+'\n'+read1[1][begin1:einde1+1]+'\n'
                read1_write += read1[2]+'\n'+read1[3][begin1:einde1+1]+'\n'
                read2_write = read2[0]+'\n'+read2[1][begin2:einde2+1]+'\n'
                read2_write += read2[2]+'\n'+read2[3][begin2:einde2+1]+'\n'
                paired_1_trimmed.write(read1_write)
                paired_2_trimmed.write(read2_write)

    # Sluiten van de getrimde read files.
    paired_1_trimmed.close()
    paired_2_trimmed.close()


def trimEinden(quality):
    """
    Functie: Trimt quality score string.
    Input:   -quality: De kwaliteits string van de read in phred 64.
    Output:  Een lijst van de indexen van een getrimde read.
    Methode: Aan beiden einden wordt gekeken waar de eerste goede phred score
             wordt gevonden. Het stuk tussen deze twee goede phreds zijn de
             goedgekeurde basen.
    """    
    
    # Eerste goede phred index in de lijst vinden.
    eerste_goede_phred = len(quality)-1
    #for x in range(len(quality)):
    for x in range(len(quality)-1,-1,-1):
        if ord(quality[x])-64 >= 20:
            eerste_goede_phred = x
    
    # Laatste goede phred index in de lijst vinden.
    laaste_goede_phred = 0
    for x in range(len(quality)):
        if ord(quality[x])-64 >= 20:
            laaste_goede_phred = x
    
    # De eerste goede phred index t/m de laatste goede phred index returnen.
    return range(eerste_goede_phred, laaste_goede_phred+1)


def goedNucleotidenGehalte(quality, getrimde_indexen, max_percentage_slecht):
    """
    Functie: Checkt voor slechte basen in getrimde read.
    Input:   -quality: De kwaliteits string van de read in phred 64.
             -getrimde_indexen: Een lijst met indexen van de read die aangeven
              welke basen over zijn gebelven na het rimmen.
             -max_percentage_slecht: Het maximale aantal slechte basen die
              over mogen blijven na trimmen.
    Output:  True als de read een laag genoeg gehalte aan slechte basen heeft.
    Methode: Eerst wordt de kleinste en grootste index bepaald. Dat is het gebied
             waarin moet worden gekeken. Daarna wordt in dat gebied het aantal
             basen bepaald met een Phred onder de 20. Ten slotte wordt het precentage
             bepaald van het aantal slechte basen met de lengte van de getrimde read.
    """   
    
    # Als het geen nucleotiden heeft met een min phred van 20, dan wordt de read 
    # afgekerud.
    if getrimde_indexen == []:
        return False
    
    slecht_nuc = 0
    begin = min(getrimde_indexen)
    einde = max(getrimde_indexen)
    
    # Door het Goede gebied loopen en slechte nucleotiden bepalen.
    for x in range(begin, einde+1):
        if ord(quality[x])-64 < 20:
            slecht_nuc += 1
    
    # Wanneer het percentage goede nucleotiden voldoende is True returnen, anders False.
    if slecht_nuc * 100 / len(getrimde_indexen) <= max_percentage_slecht:
        return True
    else:
        return False


def voldoendeLengte(getrimde_indexen, min_read_len):
    """
    Functie: Bepaald of de read wel lang genoeg is.
    Input:   -getrimde_indexen: Alleen de indexen die zijn overgebelven na trimmen.
             -min_read_len: De minimale lengte die een read moet hebben om door te mogen.
    Output:  Een True als read lang genoeg is, anders False.
    """
    if len(getrimde_indexen) >= min_read_len:
        return True
    else:
        return False
            

# functie wordt niet gebruikt, maar kan handig zijn ter handmatige controle
def markBadPhred(quality):
    """
    Functie: Zet slechte phred characters om in een '@', zodat snel kan worden
             gecontroleerd of een nieuwe trim methode werkt.
    Input:   -quality: de qualiteit van een read in een lijst van scores als
            de illumina phred characters.
    Output:  geeft een string van de phred characters waarbij alle characters die een lagere             
             phred representeren dan 20 als een '@' worden weergegeven.
    """
    goodAndBad = ""
    for x in range(0,len(quality)):
        if ord(quality[x])-64 < 20:
            goodAndBad+="@"
        else:
            goodAndBad+=quality[x]
    return goodAndBad


def main(argv):
    """
    Functie: Het uitvoeren van het script of hulp tonen voor het
             gebruik van het script.
    Input:   (Argumenten worden meegegeven via de command line)
             Een string argument "-h", voor help; of:
             -file1: reads fastq file 1, met path, van paired end reads.
             -file2: reads fastq file 2, met path, van paired end reads.
             -max_percentage_slecht_na_trim: Een maximaal aanvaardbaar
             percentage voor het GC gehalte.
             -min_read_len: De korste read lengte die aanvaardbaar is.
             -file_trimmed: Een path met basis voor de file namen b.v.: ~/map/run1
    Output:  Help voor het script (bij de "-h" optie) of twee getrimde read files.
    Gebruik: -Voor help: python2 DeelOpdracht2Majops.py -h
             -Voor uitvoering: python2 DeelOpdracht2Majops.py "path/readfile1" 
             "path/readfile2" <integer maximaal percentage slechte basen> 
             <integer minimale read lengte> <path/prefix_voor_getrimde_reads>
    """
    # Ophalen argumenten
    if len(argv) >= 2:
        if argv[1] == "-h" or argv[1] == "--h" or argv[1] == "-help" or argv[1] == "--help":
            showUsageInformation()
    
    file1 = argv[1]
    file2 = argv[2]
    max_percentage_slecht_na_trim = int(argv[3])
    min_read_len = int(argv[4])
    file_trimmed = str(argv[5])
    
    # Uitvoeren trimmen.
    trim(file1, file2, file_trimmed, max_percentage_slecht_na_trim, min_read_len)


if __name__ == "__main__":
    main(sys.argv)