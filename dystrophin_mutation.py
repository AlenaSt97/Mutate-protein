from gene_classes import *
import os.path

def write_to_file(text, filename, chars_per_line=70):
    outfile = open(filename, 'w')
    for i in range(0, len(text), chars_per_line):
        start = i
        end = start + chars_per_line
        outfile.write(text[start:end] + '\n')
    outfile.close()

def run_experiment(N,markov_chain):
    urlbase='https://raw.githubusercontent.com/AlenaSt97/Mutate-protein/main/'
    gene_file='dystrophin_v3.txt'
    exons_table='exon_positions_v4.txt'
    translation_table='nk_amino_translation.txt'
    protein=Protein((urlbase,gene_file),(urlbase,exons_table),\
                    (urlbase,translation_table),N,markov_chain)
    if N==0 and markov_chain is None:
        dys_prot=protein.get_product()
        return dys_prot
    else:
        count_norm=0
        count_mut=0
        for i in range(N):
            protein.introduce_mutations(markov_chain)
            dys_prot=protein.get_product()
            mut_len=len(dys_prot)
            if mut_len==norm_len:
                count_norm+=1
            else:
                count_mut+=1
                if markov_chain is None:
                    write_to_file(dys_prot,\
                                  f'mutate_proteins_equal/mutant_dystrophin_{i}.txt')
                else:
                    write_to_file(dys_prot,\
                                  f'mutate_proteins_non_equal/mutant_dystrophin_{i}.txt')
        print('normal:',count_norm)
        print('mutant:',count_mut)

if __name__ == '__main__':
    #create normal protein
    normal_prot=run_experiment(0,None)
    norm_len=len(normal_prot)
    if os.path.exists('normal_dystrophin.txt') is True:
        print('File already exists')
        pass
    else:
        write_to_file(normal_prot,'normal_dystrophin.txt')
        print('Protein sequence written to file')
    #introduce random SNP into exons
    print('Equal substitution probabilities')
    mutant_protein_v1=run_experiment(10000,None)
    
    markov_chain= {'A': {'A': 0.842,\
                         'C': 0.04,\
                         'G': 0.065,\
                         'T': 0.053},\
                   'C': {'A': 0.093,\
                         'C': 0.673,\
                         'G': 0.113,\
                         'T': 0.121},\
                   'G': {'A': 0.117,\
                         'C': 0.109,\
                         'G': 0.653,\
                         'T': 0.121},\
                   'T': {'A': 0.049,\
                         'C': 0.089,\
                         'G': 0.028,\
                         'T': 0.834}}

    print('Non-equal substitution probabilities')
    mutant_protein_v2=run_experiment(10000,markov_chain)
