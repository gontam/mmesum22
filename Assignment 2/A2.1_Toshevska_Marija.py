
# Reading the file
x = open(r'C:\Users\Marija\PycharmProjects\A2.1_Toshevska_Marija\BD137219.1.fasta','r')
a = x.read()
x.close()
#print(a)

dataset = {"TTT": "F", "TCT": "S" ,"TAT": "Y" ,"TGT": "C" ,"TTC": "F" ,"TCC": "S" ,"TAC": "Y" ,
           "TGC": "C" ,"TTA": "L" ,"TCA": "S" ,"TAA": "*" ,"TGA": "*" ,"TTG": "L" ,"TCG": "S" ,
           "TAG": "*" ,"TGG": "W" ,"CTT": "L" ,"CCT": "P" ,"CAT": "H" ,"CGT": "R" ,"CTC": "L" ,
           "CCC": "P" ,"CAC": "H" ,"CGC": "R" ,"CTA": "L" ,"CCA": "P" ,"CAA": "Q" ,"CGA": "R",
           "CTG": "L" ,"CCG": "P" ,"CAG": "Q","CGG": "R" ,"ATT": "I" ,"ACT": "T" ,"AAT": "N" ,
           "AGT": "S" ,"ATC": "I" ,"ACC": "T" ,"AAC": "N" ,"AGC": "S","ATA": "I" ,"ACA": "T" ,
           "AAA": "K" ,"AGA": "R" ,"ATG": "M" ,"ACG": "T","AAG": "K" ,"AGG": "R" ,"GTT": "V" ,
           "GCT": "A" ,"GAT": "D" ,"GGT": "G" ,"GTC": "V" ,"GCC": "A" ,"GAC": "D" ,"GGC": "G" ,
           "GTA": "V" ,"GCA": "A" ,"GAA": "E" ,"GGA": "G" ,"GTG": "V" ,"GCG": "A" ,"GAG": "E" ,"GGG": "G"
           }
dataset1 = {"TTT": "Phe", "TCT": "Ser" ,"TAT": "Tyr" ,"TGT": "Cys" ,"TTC": "Phe" ,"TCC": "Ser" ,"TAC": "Tyr" ,
           "TGC": "Cys" ,"TTA": "Leu" ,"TCA": "Ser" ,"TAA": "Ter" ,"TGA": "Ter" ,"TTG": "Leu" ,"TCG": "Ser" ,
           "TAG": "Ter" ,"TGG": "Trp" ,"CTT": "Leu" ,"CCT": "Pro" ,"CAT": "His" ,"CGT": "Arg" ,"CTC": "Leu" ,
           "CCC": "Pro" ,"CAC": "His" ,"CGC": "Arg" ,"CTA": "Leu" ,"CCA": "Pro" ,"CAA": "Gln" ,"CGA": "Arg",
           "CTG": "Leu" ,"CCG": "Pro" ,"CAG": "Gln","CGG": "Arg" ,"ATT": "Ile" ,"ACT": "Thr" ,"AAT": "Asn" ,
           "AGT": "Ser" ,"ATC": "Ile" ,"ACC": "Thr" ,"AAC": "Asn" ,"AGC": "Ser","ATA": "Ile" ,"ACA": "Thr" ,
           "AAA": "Lys" ,"AGA": "Arg" ,"ATG": "Met" ,"ACG": "Thr","AAG": "Lys" ,"AGG": "Arg" ,"GTT": "Val" ,
           "GCT": "Ala" ,"GAT": "Asp" ,"GGT": "Gly" ,"GTC": "Val" ,"GCC": "Ala" ,"GAC": "Asp" ,"GGC": "Gly" ,
           "GTA": "Val" ,"GCA": "Ala" ,"GAA": "Glu" ,"GGA": "Gly" ,"GTG": "Val" ,"GCG": "Ala" ,"GAG": "Glu" ,"GGG": "Gly"
           }

#user input
val = input("If you want one-letter coding please enter M. If you want three-letter coding, please enter Met: ")
protein = ""

#if user enters M
if val == "M":
 for i in range(0, len(a), 3):
    codon = a[i:i+3]
    if str(dataset.get(codon)) == "None":
        continue
    protein += str(dataset.get(codon))

 print(protein)

#if user enters Met
elif val == "Met":
    for i in range(0, len(a), 3):
        codon = a[i:i + 3]
        if str(dataset1.get(codon)) == "None":
            continue
        protein += str(dataset1.get(codon))

    print(protein)

with open("BD137219.1_translated_Toshevska.txt", "a") as o:
    o.write(protein)
