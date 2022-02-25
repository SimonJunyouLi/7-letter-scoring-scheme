def genopos (given):
    '''(num) -> str

    Return the the protein product name when given any genome
    position within the rnage of its gene.

    Restriction: The num must be within the range

    >>> genoposition(11052)
    '3C-like proteinse'
    '''

    pos = [540,2454,8289,9789,10707,11577,11826,12420,12759,13176,13203,15971,
           17774,19355,20393,21287,25127,25979,26257,26936,27128,27628,28008,29292,29409]
    name = ["leader protein","nsp2","nsp3","nsp4","3C-like proteinase","nsp6","nsp7",
            "nsp8","nsp9","nsp10","Segment 1 RNA-dependent RNA polymerase","Segment 2 RNA-dependent RNA polymerase",
            "helicase","3'-to-5' exonuclease","endoRNAse","2'-O-ribose methyltransferase",
            "surface glycoprotein","ORF3a protein","envelope protein","membrane glycoprotein",
            "ORF6 protein","ORF7a protein","ORF7b","ORF8 protein","nucleocapsid phosphoprotein",
            "ORF10 protein"]
    i = 0
    while(given > pos[i]):
        i += 1
    return name[i]

def position(given):
    pos = [(0, 540),(540, 2454),(2454, 8289),(8289, 9789),(9789, 10707),(10707, 11577),(11577, 11826),(11826, 12420),(12420, 12759),(12759, 13176),
           (13176, 13203),(13202, 15971),(15971, 17774),(17774, 19355),(19355, 20393),(20393, 21287),(21297, 25119),(25127, 25955),(25979, 26207),
           (26257, 26926),(26936, 27122),(27128, 27494),(27490, 27622),(27628, 27994),(28008, 29268),(29292, 29409)]
    name = ["leader protein","nsp2","nsp3","nsp4","3C-like proteinase","nsp6","nsp7",
            "nsp8","nsp9","nsp10","Segment 1 RNA-dependent RNA polymerase","Segment 2 RNA-dependent RNA polymerase",
            "helicase","3'-to-5' exonuclease","endoRNAse","2'-O-ribose methyltransferase",
            "surface glycoprotein","ORF3a protein","envelope protein","membrane glycoprotein","ORF6 protein",
            "ORF7a protein","ORF7b","ORF8 protein","nucleocapsid phosphoprotein","ORF10 protein"]
    i = 0
    while(given >= pos[i][1]):
        i += 1
    return (name[i],pos[i][0])

def amino_position(given):
    pos = [180,818,2763,3263,3569,3859,3942,4140,4253,4392,4401,5324,5925,6452,6798,7096,8370,8646,8722,8945,9007,9129,9173,9295,9715,9754]
    name = ["leader protein","nsp2","nsp3","nsp4","3C-like proteinase","nsp6","nsp7",
            "nsp8","nsp9","nsp10","Segment 1 RNA-dependent RNA polymerase","Segment 2 RNA-dependent RNA polymerase",
            "helicase","3'-to-5' exonuclease","endoRNAse","2'-O-ribose methyltransferase",
            "surface glycoprotein","ORF3a protein","envelope protein","membrane glycoprotein","ORF6 protein",
            "ORF7a protein","ORF7b","ORF8 protein","nucleocapsid phosphoprotein","ORF10 protein"]
    i = 0
    while(given > pos[i]-1):
        i += 1
    return (name[i],pos[i])
        
def aminoacid_name(triplet):
    tripletToAAname = {'TTT':'Phenylalanine','TTC':'Phenylalanine','TTA':'Leucine','TTG':'Leucine','CTT':'Leucine',
                       'CTC':'Leucine','CTA':'Leucine','CTG':'Leucine','ATT':'Isoleucine','ATC':'Isoleucine','ATA':'Isoleucine',
                       'ATG':'Methionine','GTT':'Valine','GTC':'Valine','GTA':'Valine','GTG':'Valine','TCT':'Serine','TCC':'Serine',
                       'TCA':'Serine','TCG':'Serine','CCT':'Prolinex','CCC':'Proline','CCA':'Proline','CCG':'Proline','ACT':'Threonine',
                       'ACC':'Threonine','ACA':'Threonine','ACG':'Threonine','GCT':'Alanine','GCC':'Alanine','GCA':'Alanine','GCG':'Alanine',
                       'TAT':'Tyrosine','TAC':'Tyrosine','TAA':'Stop','TAG':'Stop','CAT':'Histidine','CAC':'Histidine','CAA':'Glutamine',
                       'CAG':'Glutamine','AAT':'Asparagine','AAC':'Asparagine','AAA':'Lysine','AAG':'Lysine','GAT':'Aspartic acid','GAC':'Aspartic acid',
                       'GAA':'Glutamic acid','GAG':'Glutamic acid','TGT':'Cysteine','TGC':'Cysteine','TGA':'Stop','TGG':'Tryptophan',
                       'CGT':'Arginine','CGC':'Arginine','CGA':'Arginine','CGG':'Arginine','AGT':'Serine','AGC':'Serine','AGA':'Arginine','AGG':'Arginine',
                       'GGT':'Glycine','GGC':'Glycine','GGA':'Glycine','GGG':'Glycine'}
    return tripletToAAname[triplet]


def aminoacid_letter(triplet):
    tripletToAAname = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L',
                       'CTC':'L','CTA':'L','CTG':'L','ATT':'I','ATC':'I','ATA':'I',
                       'ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V','TCT':'S','TCC':'S',
                       'TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P','ACT':'T',
                       'ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
                       'TAT':'Y','TAC':'Y','TAA':'Ochre','TAG':'Amber','CAT':'H','CAC':'H','CAA':'Q',
                       'CAG':'Q','AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D',
                       'GAA':'E','GAG':'E','TGT':'C','TGC':'C','TGA':'Opal','TGG':'W',
                       'CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
                       'GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
    return tripletToAAname[triplet]


def specific_pos(num):
    protein = position(num)
    start_pos = protein[1]
    amino_num = num - start_pos + 1
    if  amino_num % 3 == 0:
        amino_pos = amino_num //3
    elif amino_num % 3 == 2:
        amino_pos = (amino_num - 2)//3
    else:
        amino_pos = (amino_num - 1)//3

    return (protein[0],amino_pos)

def alignment():
    '''
    Alignment, track mutations, sort, and determine the effects of mutations
    amino acids.
    '''
    f = open('Final.fasta', 'r')
    header = f.readline()
    reference = True
    mut_count = 0
    mut_total = 0

    mutation = []
    
    while True:
        
        line = f.readline()
        if not line: break
        
        genome = ''
        single = []
        
        while True:
            line = line.strip('\n')
            genome += line
            line = f.readline()
            if any(i in '>' for i in line):
                break
            elif not line: break

        if reference:
            refgeno = genome
            reference = False

        if refgeno != genome:
            index = 0
            for i in refgeno:
                if (i != genome[index]) & ((genome[index] == 'A')or(genome[index] == 'T')or(genome[index] == 'C')or(genome[index] == 'G')):
                    single.append((index,i +'->'+genome[index]))
                    print((index,i +'->'+genome[index]))
                    mut_total = mut_total + 1
                    if len(single)>100: #if it is a bad genome
                        single = []
                        mut_total = mut_total - 101
                        break
                index += 1
                 
            for l in single:
                if not mutation:
                    mutation.append((l[0],l[1],1))
                else:
                    temp = [item for item in mutation if (item[0] == l[0]) & (item[1] == l[1])] #find the index number and type of mutation in the mutation list that match
                    if not temp:
                        mutation.append((l[0],l[1],1))
                    else:
                        temp_pos = mutation.index(temp[0])
                        temp_tuple = temp[0]
                        count = temp_tuple[2] + 1
                        mutation[temp_pos] = (l[0],l[1],count)            
            
        line = ''
        
    f.close()
    mutation.sort(key=lambda x: x[2], reverse = True)#sorting from highest to lowest 
    print(mutation) 

    amino_acid_mutated = [] #list that contains the change in amino acid
    
    for i in mutation:
        geno_pos = i[0]
        triplet_pos = geno_pos % 3 #find out the pos of changed nucleoties in the triplets
        changed = i[1]
        arrow_pos = changed.index('->')+2
        changed = changed[arrow_pos:] #the mutated nucleotides
        print("changed: " + changed)
        if triplet_pos == 2: # when the last one is mutated
            codon = refgeno[geno_pos-2] + refgeno[geno_pos-1] + refgeno[geno_pos]
            codon_change = refgeno[geno_pos-2] + refgeno[geno_pos-1] + changed
        elif triplet_pos == 0: # when the first one is mutated
            codon = refgeno[geno_pos] + refgeno[geno_pos+1] + refgeno[geno_pos+2]
            codon_change = changed + refgeno[geno_pos+1] + refgeno[geno_pos+2]
        else: # when the second/middle one is mutated
            codon = refgeno[geno_pos-1] + refgeno[geno_pos] + refgeno[geno_pos+1]
            codon_change = refgeno[geno_pos-1] + changed + refgeno[geno_pos+1]
        print(codon)
        print(codon_change)
        ref_amino = aminoacid_name(codon) #find the corresponding amino acid
        mut_amino = aminoacid_name(codon_change)
        print(ref_amino)
        print(mut_amino)
        if(ref_amino != mut_amino): #if the amino acids are not the same (a mutation affecting the amino acid took place)
            amino_acid_mutated.append((i[0],i[1],i[2],ref_amino+'->'+mut_amino))
            mut_count = mut_count + 1
    amino_acid_mutated_top = amino_acid_mutated[:100]
    print(amino_acid_mutated_top)
    for i in amino_acid_mutated_top:
        print(i)
    print(mut_count / mut_total * 100)

def clean(name):
    '''(str)->file
    '''

    f = open(name+'.fasta', 'r')
    final = open('Final.fasta', 'w')
    header = f.readline()
    header1 = ''

    nsp1LeadChars = 'ATGGAGAGCC'

    while True:
         
        line = f.readline()
        if not line: break
        
        genome = ''
        
        while True:
            line = line.strip('\n')
            genome += line
            line = f.readline()
            if any(i in '>' for i in line):
                header1 = line
                break
            elif not line: break

        try:
            startingpos = genome.index(nsp1LeadChars)
            truncatedgenome = genome[startingpos:startingpos+29409]
            
            if len(truncatedgenome) == 29409:
                final.write(header)
                temp = ''
                state = True

                for x in truncatedgenome:
                    temp += x
                    state = True
                    if len(temp) == 60:
                        final.write(temp)
                        final.write('\n')
                        temp = ''
                        state = False
                if(state):
                    final.write(temp + '\n')
                    
                line = ''
                header = header1
        except:
            pass

    f.close()
    final.close()

def h_value():
    import math
    f = open('Final.fasta', 'r')
    header = f.readline()
    F=[0]*9754
    L=[0]*9754
    I=[0]*9754
    M=[0]*9754
    V=[0]*9754
    S=[0]*9754
    P=[0]*9754
    T=[0]*9754
    A=[0]*9754
    Y=[0]*9754
    H=[0]*9754
    Q=[0]*9754
    N=[0]*9754
    K=[0]*9754
    D=[0]*9754
    E=[0]*9754
    C=[0]*9754
    W=[0]*9754
    R=[0]*9754
    G=[0]*9754
    Amber=[0]*9754
    Opal=[0]*9754
    Ochre=[0]*9754
    total = [0]*9754
    h = [0]*9754
    
    pos = [(0, 540),(540, 2454),(2454, 8289),(8289, 9789),(9789, 10707),(10707, 11577),(11577, 11826),(11826, 12420),(12420, 12759),(12759, 13176),
           (13176, 13203),(13202, 15971),(15971, 17774),(17774, 19355),(19355, 20393),(20393, 21287),(21297, 25119),(25127, 25955),(25979, 26207),
           (26257, 26926),(26936, 27122),(27128, 27494),(27490, 27622),(27628, 27994),(28008, 29268),(29292, 29409)]
    name = ["leader protein","nsp2","nsp3","nsp4","3C-like proteinase","nsp6","nsp7",
            "nsp8","nsp9","nsp10","Segment 1 RNA-dependent RNA polymerase","Segment 2 RNA-dependent RNA polymerase",
            "helicase","3'-to-5' exonuclease","endoRNAse","2'-O-ribose methyltransferase",
            "surface glycoprotein","ORF3a protein","envelope protein","membrane glycoprotein","ORF6 protein",
            "ORF7a protein","ORF7b","ORF8 protein","nucleocapsid phosphoprotein","ORF10 protein"]
    import csv
    import random
    with open('data.csv', 'w') as p:
        writer = csv.writer(p)
        writer.writerow(['Overall Position']+['Protein Name']+['Amino Acid position within protein']+['H Value of AA Groups'])
        
        while True:
            
            line = f.readline()
            if not line: break
            
            genome = ''
            single = []
            position_count = 0
            
            while True:
                line = line.strip('\n')
                genome += line
                line = f.readline()
                if any(i in '>' for i in line):
                    break
                elif not line: break

            for k in pos:
                protein_section = genome[k[0]:k[1]]
               
                for i in range((k[1]-k[0])//3):
                    triplet = protein_section[3*i]+protein_section[3*i+1]+protein_section[3*i+2]
                   
                    if ((triplet[0] == 'A')or(triplet[0] == 'T')or(triplet[0] == 'C')or(triplet[0] == 'G'))&((triplet[1] == 'A')
                        or(triplet[1] == 'T')or(triplet[1] == 'C')or(triplet[1] == 'G'))&((triplet[2] == 'A')or(triplet[2] == 'T')or(triplet[2] == 'C')or(triplet[2] == 'G')):
                        name = aminoacid_letter(triplet)
                        count = vars()[name][position_count]+1
                        vars()[name][position_count] = count
                    else:
                        pass
                    position_count += 1
        print(Ochre)       
        for l in amino_acid_list:
            for k in range(0,9754):
                total[k] = total[k] + l[k]

        #print(total)

        for l in amino_acid_list:
            for k in range(0,9754):
                l[k] = l[k] / total[k]

        #print(F)
        
        for l in amino_acid_list:
            for k in range(0,9754):
                if l[k]!=0:
                    h[k] = h[k]+ -l[k]*math.log(l[k],2)
                    
        print(h)

        h_sorted = sorted(h,reverse = True)

        print(h_sorted[0],h_sorted[1],h_sorted[2])

        position_tuple  = ()
        protein = ()
        positional_count = 0
        print_count = 1

        for m in range(0,9754):
            protein = amino_position(m)
            if position_tuple == protein[1]:
                positional_count += 1
            else:
                positional_count = 1
                position_tuple = protein[1]
            print(protein[0],positional_count,h[m])
            writer.writerow([print_count]+[protein[0]]+[positional_count]+[h[m]])
            print_count += 1
            
def amino_acid_groups(given):
    amino_acid_to_group = {'A':'Aliphatic','V':'Aliphatic','L':'Aliphatic','I':'Aliphatic','M':'Aliphatic','C':'Aliphatic','F':'Aromatic',
                           'F':'Aromatic','W':'Aromatic','Y':'Aromatic','H':'Aromatic','S':'Polar','T':'Polar','N':'Polar','Q':'Polar',
                           'K':'Positive','R':'Positive','D':'Negative','E':'Negative','P':'P','G':'G','Ochre':'Stop','Opal':'Stop','Amber':'Stop'}
    return amino_acid_to_group[given]

def h_value_amino_groups():
    import math
    f = open('Final.fasta', 'r')
    header = f.readline()
    Aliphatic = [0]*9754
    Aromatic = [0]*9754
    Polar = [0]*9754
    Positive = [0]*9754
    Negative = [0]*9754
    P = [0]*9754
    G = [0]*9754
    Stop = [0]*9754
    total = [0]*9754
    h = [0]*9754
    amino_acid_list = [Aliphatic,Aromatic,Polar,Positive,Negative,P,G,Stop]
    amino_acid_list_name = ['Aliphatic','Aromatic','Polar','Positive','Negative','P','G','Stop']
    pos = [(0, 540),(540, 2454),(2454, 8289),(8289, 9789),(9789, 10707),(10707, 11577),(11577, 11826),(11826, 12420),(12420, 12759),(12759, 13176),
           (13176, 13203),(13202, 15971),(15971, 17774),(17774, 19355),(19355, 20393),(20393, 21287),(21297, 25119),(25127, 25955),(25979, 26207),
           (26257, 26926),(26936, 27122),(27128, 27494),(27490, 27622),(27628, 27994),(28008, 29268),(29292, 29409)]
    name = ["leader protein","nsp2","nsp3","nsp4","3C-like proteinase","nsp6","nsp7",
            "nsp8","nsp9","nsp10","Segment 1 RNA-dependent RNA polymerase","Segment 2 RNA-dependent RNA polymerase",
            "helicase","3'-to-5' exonuclease","endoRNAse","2'-O-ribose methyltransferase",
            "surface glycoprotein","ORF3a protein","envelope protein","membrane glycoprotein","ORF6 protein",
            "ORF7a protein","ORF7b","ORF8 protein","nucleocapsid phosphoprotein","ORF10 protein"]
    import csv
    import random
    with open('data.csv', 'w') as p:
        writer = csv.writer(p)
        writer.writerow(['Count']+['AA Group']+['Protein Name']+['Amino Acid position within protein']+['H Value'])
        while True:
            
            line = f.readline()
            if not line: break
            
            genome = ''
            single = []
            position_count = 0
            
            while True:
                line = line.strip('\n')
                genome += line
                line = f.readline()
                if any(i in '>' for i in line):
                    break
                elif not line: break

            for k in pos:
                protein_section = genome[k[0]:k[1]]
               
                for i in range((k[1]-k[0])//3):
                    triplet = protein_section[3*i]+protein_section[3*i+1]+protein_section[3*i+2]
                   
                    if ((triplet[0] == 'A')or(triplet[0] == 'T')or(triplet[0] == 'C')or(triplet[0] == 'G'))&((triplet[1] == 'A')
                        or(triplet[1] == 'T')or(triplet[1] == 'C')or(triplet[1] == 'G'))&((triplet[2] == 'A')or(triplet[2] == 'T')or(triplet[2] == 'C')or(triplet[2] == 'G')):
                        name = aminoacid_letter(triplet)
                        temp_amino_group = amino_acid_groups(name)
                        count = vars()[temp_amino_group][position_count]+1
                        vars()[temp_amino_group][position_count] = count
                    else:
                        pass
                    position_count += 1
               
        for l in amino_acid_list:
            for k in range(0,9754):
                total[k] = total[k] + l[k]

        #print(total)

        for l in amino_acid_list:
            for k in range(0,9754):
                l[k] = l[k] / total[k]

        #print(F)
        position_tuple  = ()
        protein = ()
        positional_count = 0
        
        
        for l in amino_acid_list:
            ind_pos = amino_acid_list.index(l)
            overall_count = 1
            for k in range(0,9754):
                if l[k]!=0:
                    l[k] = -l[k]*math.log(l[k],2)
                protein = amino_position(k)
                if position_tuple == protein[1]:
                    positional_count += 1
                else:
                    
                    positional_count = 1
                    position_tuple = protein[1]
                    
                print(amino_acid_list_name[ind_pos],protein[0],positional_count,l[k])
                writer.writerow([overall_count]+[amino_acid_list_name[ind_pos]]+[protein[0]]+[positional_count]+[l[k]])
                overall_count += 1
                    