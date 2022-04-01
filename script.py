
#Couple of different ways to do this in this file.



#list number of times nucleotides occur and at what percentage in genome

seq_filename = "Ecoli_genome.txt"
total_length = 0
nucleotide = {}  #creates an empty dictionary

seq_file = open(seq_filename, "r") #r denotes "reading" of the file
for raw_line in seq_file:
    line = raw_line.rstrip("\r\n")
    length = len(line) # Python function to calculate the length of a string
    for nuc in line:
        if nucleotide.has_key(nuc):
            nucleotide[nuc] += 1
        else:
            nucleotide[nuc] = 1
    total_length += length
seq_file.close()
for n in nucleotide.keys():
    fraction = 100.0 * nucleotide[n] / total_length
    print "The nucleotide {0} occurs {1} times, or {2} %".format(n, nucleotide[n], fraction)



#reads codon squences for DNA/RNA pieces. Unsure for DNa but this works for mRNA
from sys import argv
script, filename = argv

def translate_dna(sequence):

    codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    proteinsequence = ''
    start = sequence.find('ATG')
    sequencestart = sequence[int(start):]
    stop = sequencestart.find('TAA')
    cds = str(sequencestart[:int(stop)+3])

    for n in range(0,len(cds),3):
        if cds[n:n+3] in codontable
            proteinsequence += codontable[cds[n:n+3]]
            print proteinsequence
        sequence = ''


header = ''
sequence = ''
for line in open(filename):
    if line[0] == ">":
        if header != '':
            print header
            translate_dna(sequence)

        header = line.strip()
        sequence = ''
    else:
        sequence += line.strip()

print header 
translate_dna(sequence)


#or this code will work as well:
from itertools import takewhile

def translate_dna(sequence, codontable, stop_codons = ('TAA', 'TGA', 'TAG')):       
    start = sequence.find('ATG')

    # Take sequence from the first start codon
    trimmed_sequence = sequence[start:]

    # Split it into triplets
    codons = [trimmed_sequence[i:i+3] for i in range(0, len(trimmed_sequence), 3)]
    print(len(codons))
    print(trimmed_sequence)
    print(codons)

    # Take all codons until first stop codon
    coding_sequence  =  takewhile(lambda x: x not in stop_codons and len(x) == 3 , codons)

    # Translate and join into string
    protein_sequence = ''.join([codontable[codon] for codon in coding_sequence])

    # This line assumes there is always stop codon in the sequence
    return "{0}_".format(protein_sequence)

