import os
import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser("Aligns sequences to a reference and identifies variants.")
parser.add_argument("reference", help="FASTA file of reference sequence")
parser.add_argument("sequences", help="FASTA file of sequences to be aligned.")
parser.add_argument("site", help="Zero-indexed position of site of interest. Can be protein or nucleotide.",
                    type=int)
parser.add_argument("length", help="Number of sites to be considered, i.e. 3 for a codon.",
                    type=int)
args = parser.parse_args()

SeqReference = args.reference
SeqsToAlign = args.sequences
SiteOfInterest = args.site
SiteSize = args.length

# Given aligned reference and query sequences as strings,
# along with a ZERO-INDEXED site number (whether nucleotide or codon),
# returns the query genotype at that site.
def CompareAlignedSequences(ref, query, site, length):
    RefPosition = 0
    AlignPosition = 0
    while RefPosition < site:
        if ref[AlignPosition] != '-':
            RefPosition+=1
        AlignPosition+=1
    ReturnLength = 0
    ReturnSequence = ""
    try:
        while ReturnLength < length:
            if query[AlignPosition] != '-':
                ReturnSequence+=query[AlignPosition]
                ReturnLength +=1
            AlignPosition+=1
        return ReturnSequence
    except:
        return "---"



def main():

    FileStem = (SeqsToAlign.split('.')[0]).rsplit('/')[-1]

    #######################################################################
    # Align each sequence to the reference using needle.
    # Output the file as a FASTA format file in the directory alignments.
    # For n sequences, the alignment file will consist of n aligned pairs,
    # where the first sequence in each pair is the reference.
    #######################################################################

    SeqsAligned = "alignments/" + FileStem + "-aligned.fasta"
    if not os.path.exists(SeqsAligned):
        try:
            os.system('needle -asequence %s -bsequence %s '
                      '-gapopen 10.0 -gapextend 0.5 '
                      '-outfile %s -aformat fasta' %
                      (SeqReference, SeqsToAlign, SeqsAligned))
        except:
            print("Alignment error.")
    print("Query sequences aligned.")

    ########################################################################
    # Read in the aligned sequences, one pair at a time.
    # Determine the genotype of the query sequence at position 151.
    # Output the sequence names and their genotype at position 151.
    ########################################################################

    OutGenotypes = "out/" + FileStem + "-site" + str(SiteOfInterest+1) + "-genotypes.txt"

    # Read in the aligned sequences, one pair at a time.
    # Do not read the entire file into memory.
    with open(SeqsAligned,'r') as f:

        RefSeq = ""
        RefName = ""
        QuerySeq = ""
        QueryName = ""
        QueryNames = []
        QueryGenotypes = []

        line = f.readline()
        while line[0] == '>':

            # Store the reference sequence.
            RefName = line.strip('\n')
            nextline = next(f)
            while nextline[0] != '>':
                RefSeq += nextline.strip('\n')
                nextline = next(f)

            # Store the query sequence.
            QueryName = nextline.strip('\n')
            nextline = next(f)
            while nextline[0] != '>':
                QuerySeq += nextline.strip('\n')
                try:
                    nextline = next(f)
                except:
                    break

            # Analyze the genotype at site indicated for the reference and query sequences.
            QueryNames.append(QueryName)
            QueryGenotypes.append(CompareAlignedSequences(RefSeq, QuerySeq, SiteOfInterest, SiteSize))

            # Clear the current ref and query sequences.
            # Continue reading the file.
            RefSeq = ""
            QuerySeq = ""
            line = nextline
    # Output the genotype information into the OutGenotypes file.
    with open(OutGenotypes,'w') as f:
        for i in range(len(QueryNames)):
            f.write("%s\t%s\n" % (QueryGenotypes[i], QueryNames[i]))

    print("%d sequences analyzed." % (len(QueryNames)))


    #########################################################################
    # Parse the headers corresponding to the genotypes.
    # This parsing depends on the header format.
    # Parse out strains' year, country, and region and get rid of whitespace.
    #########################################################################

    OutGenotypesParsed = "out/" + FileStem + "-site" + str(SiteOfInterest+1) +  "-genotypes-parsed.txt"

    QueryYears = []
    QueryCountries = []
    QueryRegions = []
    QueryPassages = []

    for Name in QueryNames:
        if "Genbank" in FileStem:
            if "nucleotide" in FileStem:
                QueryRegions.append(((Name.split(" A/")[1]).split('/')[0]).replace(" ",""))
                QueryYear = (((Name.split(" A/")[1]).split('/')[2]).split(' ')[0]).replace(" ","")
                if len(QueryYear) != 4:
                    QueryYear = "20" + QueryYear
                QueryYears.append(QueryYear)
                QueryCountries.append(Name.split("(NA)")[1].split(" {")[0].replace(" ",""))
            else:
                QueryRegions.append(((Name.split(" A/")[1]).split('/')[0]).replace(" ",""))
                QueryYear = (((Name.split(" A/")[1]).split('/')[2]).split(' ')[0]).replace(" ","")
                if len(QueryYear) != 4:
                    QueryYear = "20" + QueryYear
                QueryYears.append(QueryYear)
                QueryCountries.append(Name.split(" NA ")[1].split(" gb|")[0].replace(" ",""))

        elif "GISAID" in FileStem:
            try:
                QueryRegions.append(((Name.split(">A_")[1]).split('_')[0]).split('/')[0].replace(" ",""))
            except:
                QueryRegions.append("0000")
            try:
                QueryYear = Name.split(" ")[-1][0:4]
            except:
                QueryYear = "0000"
            if len(QueryYear) == 2:
                QueryYear = "20" + QueryYear
            QueryYears.append(QueryYear)
            QueryCountries.append("None")

        if "Washington" or "Melbourne" in FileStem:
            QueryPassages.append(Name.split(' ')[1])
        else:
            QueryPassages.append("None")

    with open(OutGenotypesParsed,'w') as f:
        f.write("Genotype\tRegion\tYear\tCountry\tPassage\tFullName\n")
        for i in range(len(QueryNames)):
            f.write("%s\t%s\t%s\t%s\t%s\t%s\n" %
                    (QueryGenotypes[i],QueryRegions[i],
                     QueryYears[i],QueryCountries[i],
                     QueryPassages[i],QueryNames[i].replace(" ","")))



    print("Hello world!")

if __name__ == '__main__':
    main()