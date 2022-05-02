#!/usr/bin/python2
import os
import sys
import linecache
from multiprocessing import Process

"""
Dit script voert een quality control uit op paired-end reads.
"""

def showUsageInformation():
    """
    Functie: Geeft info over het script en het gebruik er van.
    """
    print "" 
    print "This script gives the QC from a paired-end read file couple."
    print "Part 1 and 3 of Dirks Next Gen Sequencing Pipeline."
    print ""
    print "Way of usage:"
    print "" 
    print "    python2 DeelOpdracht1Majops.py  <Argument 1> <Argument 2> etc..."
    print ""
    print "    The arguments are as following:"
    print "        1)  Full path to paired read file 1 as string."
    print "        2)  Full path to paired read file 2 as string."
    print "        3)  Full path to a QC output file as string."
    print ""
    print "Example:"
    print ""
    print "    python2 DeelOpdracht1Majops.py \"path/paired_reads_1.fastq\" \"path/paired_reads_2.fastq\" \"path/QC_output.txt\""
    print ""
    print "QC statistics to be written to file:"
    print ""
    print "    1) Reads files names and location."
    print "    2) Read count in file."
    print "    3) Average read length."
    print "    4) Smallest read length."
    print "    5) Largest read length."
    print "    6) GC content for read total."
    print "    7) GC content for all bases at a position."
    print ""

    sys.exit()
    
    
def multi(seqs, output_file):
    """
    Functie: Paralellisatie van het lezen van een fastq bestand en 
             schrijven van een QC.
    Input:   -seqs: De volledige path van een fastq bestand.
             -output_file: De file (met zijn path) waar de QC naartoe wordt geschreven.
    Output:  Geen uitvoer. De printAntwoorden() functie doet dat al.
    """
    # Statistieken voor de twee files worden bepaald en in een lijst gezet.
    statistieken = bepaalHoeveelheden(seqs)
    # De filenaam wordt bij de statistieken gezet.
    statistieken.insert(0, seqs)
    # De statistieken worden uitgeprint.
    printAntwoorden(statistieken, output_file)
    

def bepaalHoeveelheden(bestand):
    """
    Functie: Bepalen van de statistieken van de fastq file (zie output).
    Input:   -bestand: De volledige path van een fastq bestand.
    Output:  Een lijst met de elementen [aantal_regels, totaal_basen, min_lengte,
             max_lengte, totale_gc, per_base_lijst].
    """

    # Dit zou de snelste manier zijn om het aantal regels te bepalen van een
    # file groter dan 1 GB.
    aantal_regels = 0
    for line in open(bestand).xreadlines():
        aantal_regels += 1

    # Initialistaite van een paar variabelen.
    per_base_lijst = []
    totaal_basen = 0
    min_lengte = 1000000
    max_lengte = 0
    totale_gc = 0

    # De nodige data begint niet op de eerste regel, daarom een andere index.
    for x in range(2, aantal_regels, 4):
        # Hier wordt alleen een bepaalde regel uit het bestand gelezen.
        # Het volledige bestand wordt niet ingelezen.
        regel = linecache.getline(bestand, x).strip()
        totaal_basen += len(regel)
        if len(regel) < min_lengte:
            min_lengte = len(regel)
        if len(regel) > max_lengte:
            max_lengte = len(regel)
        for y in range(len(regel)):
            if regel[y].upper() in "GC":
                totale_gc += 1
            # Er wordt een lijst gemaakt met dictionary die de GC en ATN
            # gehaltes telt per sequentie per base. Om geheugen te besparen
            # wordt de lijst dynamisch aangepasst aan de langste sequentie.
            if len(per_base_lijst) < max_lengte:
                for z in range(max_lengte - len(per_base_lijst)):
                    per_base_lijst.append({"GC": 0, "rest": 0})
            if regel[y].upper() in "GC":
                per_base_lijst[y]["GC"] += 1
            elif regel[y].upper() in "ATN":
                per_base_lijst[y]["rest"] += 1

    # De per base lijst wordt meteen gerecyled om gc percentages op te slaan.
    # Dat scheelt geheugen.
    for base in range(len(per_base_lijst)):
        per_base_lijst[base] = round(float(per_base_lijst[base]["GC"]*100)/float(
            per_base_lijst[base]["GC"]+per_base_lijst[base]["rest"]),1)

    # Output van het script.
    return [aantal_regels, totaal_basen, min_lengte,
            max_lengte, totale_gc, per_base_lijst]


def printAntwoorden(QC, output):
    """
    Functie: De verkregen gegevens over de sequenties in de file wegschrijven
             naar een bestand opgegeven door de gebruiker.
    Input:   -QC: Alle statistieken van aan fastq file in een lijst. Inhoud is in
             de strings van onderstaande code af te leiden.
             -output: Volledige path naar de output file.
    Output:  Alle statistieken worden geschreven naar het bestand aangegeven in de
             variabele "output".
    """
    schrijven = open(output, 'a')
    schrijven.write("File:              "+QC[0]+"\n")
    schrijven.write("Aantal reads:      "+str(((QC[1]-2)/4)+1)+"\n")
    schrijven.write("Gemiddelde lengte: "+str(QC[2]/(((QC[1]-2)/4)+1))+"\n")
    schrijven.write("Minimum lengte:    "+str(QC[3])+"\n")
    schrijven.write("Maximum lengte:    "+str(QC[4])+"\n")
    schrijven.write("GC content:        "+str(round(float((QC[5]*100)/float(QC[2])),1))+"%"+"\n")
    schrijven.write("GC content per positie:"+"\n")
    for x in range(0, len(QC[6]), 5):
        for y in range(x, x+5):
            if y < len(QC[6]):
                #schrijven.write(str(y+1)+" "+QC[6][y])+"%\t\t")
                zin = "%5i %5.1f%s       " % (y+1, QC[6][y], "%")
                schrijven.write(zin)
        schrijven.write("\n")
    schrijven.write("\n")
    schrijven.close()


def main(argv):
    """
    Functie: Het uitvoeren van het script of hulp tonen voor het
             gebruik van het script.
    Input:   Een string argument "-h", voor help; of 2 read file namen
             met path en de output file met path.
    Output:  Help voor het script (bij de "-h" optie) of een file 
             met de QC van de reads.
    Gebruik: -Voor help: python2 DeelOpdracht1Majops.py -h
             -Voor uitvoering python2 DeelOpdracht1Majops.py "path/readfile1" "path/readfile2" 
             "path/QCoutputfile.txt"
    """
    if len(argv) >= 2:
        if argv[1] == "-h" or argv[1] == "--h" or argv[1] == "-help" or argv[1] == "--help":
            showUsageInformation()
    
    file1 = argv[1]
    file2 = argv[2]
    output_file = argv[3]
    
    p1 = Process(target = multi(file1, output_file))
    p1.start()
    p2 = Process(target = multi(file2, output_file))
    p2.start()
    
                               
if __name__ == "__main__":
    main(sys.argv)
    