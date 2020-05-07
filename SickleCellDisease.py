# A program that translate each codon of DNA into the correct Amino Acid.

# Open and create normal DNA text file.
normal = open("normalDNA.txt", "w")
# Open and create mutated DNA text file.
mutated = open("mutatedDNA.txt", "w")

# Empty list variable for DNA SLC code.
slc = []
# Empty list variable for amino acids.
amino_acid = []

# Function for normal DNA that makes all letters in text file uppercase.
def normal_dna():
    # Open DNA text file.
    file = open("DNA.txt", "r+")
    # Read the file, save it as string stripping new lines and make all letters upper.
    dna = file.read().replace("\n", "").upper()
    # Convert string to list.
    dna_list = [dna[i:i + 3] for i in range(0, len(dna), 3)]
    # Close DNA file.
    file.close()
    # Return the DNA list.
    return dna_list

# Function for mutated DNA that replace lower "a" with upper "T".
def mutated_dna():
    # Open DNA text file.
    file = open("DNA.txt", "r+")
    # Read the file, save it as a string stripping new lines.
    dna = file.read().replace("\n", "")
    # Convert string to list.
    dna_list = [dna[i:i + 3] for i in range(0, len(dna), 3)]
    # Replace "a" in string with "T".
    dna_list = [s.replace('a', 'T') for s in dna_list]
    # Close DNA file.
    file.close()
    # Return DNA list.
    return dna_list

# Function that use created DNA list from DNA text file to translate.
def translate(dna_list):
    # List of SLC codes with their corresponding codons.
    I = ['ATT', 'ATC', 'ATA']
    L = ['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG']
    V = ['GTT', 'GTC', 'GTA', 'GTG']
    F = ['TTT', 'TTC']
    M = ['ATG']
    C = ['TGT', 'TGC']
    A = ['GCT', 'GCC', 'GCA', 'GCG']
    G = ['GGT', 'GGC', 'GGA', 'GGG']
    P = ['CCT', 'CCC', 'CCA', 'CCG']
    T = ['ACT', 'ACC', 'ACA', 'ACG']
    S = ['TCT', 'TCC', 'TCA', 'AGT', 'AGC']
    Y = ['TAT', 'TAC']
    W = ['TGG']
    Q = ['CAA', 'CAG']
    N = ['AAT', 'AAC']
    H = ['CAT', 'CAC']
    E = ['GAA', 'GAG']
    D = ['GAT', 'GAC']
    K = ['AAA', 'AAG']
    R = ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']
    Stop = ['TAA', 'TAG', 'TGA']

    # Go through DNA list and create SLC code list.
    for item in dna_list:
        if item in I:
            slc.append("I")
        elif item in L:
            slc.append("L")
        elif item in V:
            slc.append("V")
        elif item in F:
            slc.append("F")
        elif item in M:
            slc.append("M")
        elif item in C:
            slc.append("C")
        elif item in A:
            slc.append("A")
        elif item in G:
            slc.append("G")
        elif item in P:
            slc.append("P")
        elif item in T:
            slc.append("T")
        elif item in S:
            slc.append("S")
        elif item in Y:
            slc.append("Y")
        elif item in W:
            slc.append("W")
        elif item in Q:
            slc.append("Q")
        elif item in N:
            slc.append("N")
        elif item in H:
            slc.append("H")
        elif item in E:
            slc.append("E")
        elif item in D:
            slc.append("D")
        elif item in K:
            slc.append("K")
        elif item in R:
            slc.append("R")
        elif item in Stop:
            slc.append("Stop")
        else:
            slc.append("X")

    # Go through DNA list and create amino acids list.
    for item in dna_list:
        if item in I:
            amino_acid.append("Isoleucine")
        elif item in L:
            amino_acid.append("Leucine")
        elif item in V:
            amino_acid.append("Valine")
        elif item in F:
            amino_acid.append("Phenylananine")
        elif item in M:
            amino_acid.append("Methionine")
        elif item in C:
            amino_acid.append("Cysteine")
        elif item in A:
            amino_acid.append("Alanine")
        elif item in G:
            amino_acid.append("Glycine")
        elif item in P:
            amino_acid.append("Proline")
        elif item in T:
            amino_acid.append("Threonine")
        elif item in S:
            amino_acid.append("Serine")
        elif item in Y:
            amino_acid.append("Tyrosine")
        elif item in W:
            amino_acid.append("Tryptophan")
        elif item in Q:
            amino_acid.append("Glutamine")
        elif item in N:
            amino_acid.append("Asparagine")
        elif item in H:
            amino_acid.append("Histidine")
        elif item in E:
            amino_acid.append("Glutamic acid")
        elif item in D:
            amino_acid.append("Aspartic acid")
        elif item in K:
            amino_acid.append("Lysine")
        elif item in R:
            amino_acid.append("Arginine")
        elif item in Stop:
            amino_acid.append("Stop codon")
        else:
            amino_acid.append("X")

# Function to create normal DNA output.
def normal_output():
    # Calls the translate function using the normal DNA function.
    translate(normal_dna())
    # Append "Representing:" to the end of SLC list for the amino acids list.
    slc.append("\n\n(Representing: ")
    # Write SLC list to normal DNA output file.
    normal.write(str("".join(slc)))
    # Write amino acid list to normal DNA output file.
    normal.write(str(", ".join(amino_acid)))
    normal.write(")")
    # Deletes the SLC and amino acid list for re-use.
    del slc[:]
    del amino_acid[:]

# Function to create mutated DNA output.
def mutated_output():
    # Calls the translate function using the mutated DNA function.
    translate(mutated_dna())
    # Append "Representing:" to the end of SLC list for the amino acids list.
    slc.append("\n\n(Representing: ")
    # Write SLC list to normal DNA output file.
    mutated.write(str("".join(slc)))
    # Write amino acid list to normal DNA output file.
    mutated.write(str(", ".join(amino_acid)))
    mutated.write(")")
    # Deletes the SLC and amino acid list for re-use.
    del slc[:]
    del amino_acid[:]

# Call normal output function.
normal_output()
# Call mutated output function.
mutated_output()

# Close files.
normal.close()
mutated.close()
