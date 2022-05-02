import sys

def main(argv):
    """
    Functie: Het uitvoeren van het script of hulp tonen voor het
             gebruik van het script.
    Input:   Een string argument "-h", voor help; of een VCF file 
             met path en de output file met path.
    Output:  Een text file waaron staat hoeveel deleties/inserties
             en verschillende substituties er zijn.
    """
    if len(argv) >= 2:
        if argv[1] == "-h" or argv[1] == "--h" or argv[1] == "-help" or argv[1] == "--help":
            showUsageInformation()
    
    input_bestand = argv[1]
    output_bestand = argv[2]
    
    kolommen = fileInlezen(input_bestand)
    bepalenVeranderingen(kolommen, output_bestand)


def showUsageInformation():
    # Functie: Print informatie over het gebruik van het script.
    print "" 
    print "This script takes a VCF file and calculates the amount of mutation occurrences by type."
    print "Part 6 of Dirks Next Gen Sequencing Pipeline."
    print ""
    print "Way of usage:"
    print "" 
    print "    python2 DeelOpdracht6Majops.py <Argument 1> <Argument 2>"
    print ""
    print "    The arguments are as following:"
    print "        1) Full path to VCF file as string."
    print "        2) Full path to output file as string."
    print ""
    print "Example:"
    print ""
    print "    python2 DeelOpdracht6Majops.py \"/path/my_VCF.vcf\" \"/path/mutations.txt\""
    print ""
    print "Preview dummie output:"
    print ""
    print "    Deletions: 1"
    print "    Insertions: 0"
    print "    Deletions/Insertions ratio gives zero division."
    print "    A > T: 0"
    print "    A > C: 7"
    print "    A > G: 8"
    print "    T > A: 2"
    print "    T > C: 4"
    print "    T > G: 2"
    print "    C > A: 2"
    print "    C > G: 0"
    print "    C > T: 1"
    print "    G > A: 2"
    print "    G > C: 3"
    print "    G > T: 3"
    print ""
    
    sys.exit()



def fileInlezen(file_naam):
    # Functie:
    #   Inlezen van het vcf bestand.
    # Input:
    #   -file_naam: Volledige path van het vcf bestand.
    # Output:
    #   -data: twee kolommen. De eerste is voor de sequentie/nucleotide
    #          in het referentie genoom, de tweede bevat de
    #          sequentie/nucleotide voor de sample.
    openen = open(file_naam, 'r')
    data = openen.readlines()
    
    for x in range(len(data)):
        if data[x][0] == '#':
            data[x]=""
        else:
            data[x] = data[x].split("\t")
            data[x] = [ data[x][3], data[x][4] ]
    
    data = filter(None, data)
    
    return data
            
    
def bepalenVeranderingen(kolommen, output_bestand):
    # Functie:
    #   Bepalen van het aantal mutaties (per type).
    # Input:
    #   -kolommen: 2D array van strings. [["<A>", "<B>"],["<A>", "<B>"], ....]
    #              In de "A" kolom zit de sequentie/nucleotide van het referentie genoom,
    #              in de "B" kolom zit de sequentie/nucleotide van de sample.
    #   -output_bestand: Het bestand (met path) waar de statistieken naartoe worden geschreven.
    # Output:
    #   -data: twee kolommen. De eerste is voor de sequentie/nucleotide
    #          in het referentie genoom, de tweede bevat de
    #          sequentie/nucleotide voor de sample.
    deleties=0
    inserties=0
    # 0:AT 1:AC 2:AG ## 3:TA 4:TC 5:TG ## 6:CA 7:CG 8:CT ## 9:GA 10:GC 11:GT
    substituties = [0,0,0 ,0,0,0 ,0,0,0 ,0,0,0]
    # De opdracht houd geen rekening met substituties van elke mogelijke
    # nucleotide naar een specifieke (N -> G/A/T/C).
    mogelijke_mutaties = [ ["A", "T"], ["A", "C"], ["A", "G"],
                          ["T", "A"], ["T", "C"], ["T", "G"],
                          ["C", "A"], ["C", "G"], ["C", "T"],
                          ["G", "A"], ["G", "C"], ["G", "T"] ]
    
    # Tellen inserties, deleties en substituties.
    for x in range(len(kolommen)):
        if len(kolommen[x][0]) > len(kolommen[x][1]):
            deleties += 1
        elif len(kolommen[x][0]) < len(kolommen[x][1]):
            inserties += 1
        elif kolommen[x] in mogelijke_mutaties:
            index = mogelijke_mutaties.index(kolommen[x])
            substituties[index] += 1
    
    # Openen bestand.
    answerfile = open(output_bestand, 'w')
    
    # Schrijven van alle statistieken naar het opgegeven bestand.
    answerfile.write("Deletions: "+str(deleties)+"\n")
    answerfile.write("Insertions: "+str(inserties)+"\n")
    if inserties != 0:
        answerfile.write("Deletions/Insertions ratio: "+str(float(deleties)/float(inserties))+"\n")
    else:
        answerfile.write("Deletions/Insertions ratio gives zero division.\n")
    for x in range(len(mogelijke_mutaties)):
        answerfile.write(mogelijke_mutaties[x][0]+" > "+mogelijke_mutaties[x][1]+": "+str(substituties[x])+"\n")
        
    answerfile.close()
    

if __name__ == "__main__":
    main(sys.argv)