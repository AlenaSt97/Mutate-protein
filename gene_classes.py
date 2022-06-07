import urllib.request, urllib.parse, urllib.error
import random
import ssl

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

#--------------------------------------------------------------------------
#Functions:

def read_dnafile(filename):
    lines = open(filename, 'r').readlines()
    dna = ''.join([line.strip() for line in lines])
    return dna

def read_exon_regions(filename):
    positions = []
    infile=open(filename,'r')
    for line in infile:
        start, end = line.split()
        start, end = int(start), int(end)
        positions.append((start, end))
    infile.close()
    return positions

def create_mRNA(exons):
    cDNA = ''.join(exons)
    mrna = cDNA.replace('T','U')
    return mrna

#random mutations, equal probabilities
def mutate_dna_v1(exons):
    cDNA = ''.join(exons)
    dna_list=list(cDNA)
    mutation_site=random.randint(0,len(dna_list)-1)
    dna_list[mutation_site]=random.choice(list('ATGC'))
    new_cDNA=''.join(dna_list)
    mrna = new_cDNA.replace('T','U')
    return mrna

#selection of a nucleotide taking into account
#the probability of a transition from a markov_chain
def draw(discrete_probdist):
    limit = 0
    r = random.random()
    for value in discrete_probdist:
        limit += discrete_probdist[value]
        if r < limit:
            return value

#unequal probabilities of nucleotide substitutions
def mutate_dna_v2(exons,markov_chain):
    cDNA = ''.join(exons)
    dna_list = list(cDNA)
    mutation_site = random.randint(0,len(dna_list)-1)
    from_base=cDNA[mutation_site]
    to_base=draw(markov_chain[from_base])
    dna_list[mutation_site]=to_base
    new_cDNA=''.join(dna_list)
    mrna = new_cDNA.replace('T','U')
    return mrna

def read_genetic_code(filename):
    infile = open(filename, 'r')
    genetic_code = {}
    for line in infile:
        columns = line.split()
        genetic_code[columns[0]] = columns[1]
    return genetic_code

def create_protein(mrna, genetic_code):
    protein = ''
    trans_start_pos = mrna.find('AUG')
    for i in range(len(mrna[trans_start_pos:])//3):
        start = trans_start_pos + i*3
        end = start + 3
        amino = genetic_code[mrna[start:end]]
        if amino == 'X':
            break
        protein += amino
    return protein

#---------------------------------------------------------------------------
#Classes:

class Gene:
    def __init__(self,dna,exon_regions):
        
        #determine dna
        if isinstance(dna,(list,tuple,str)) and len(dna) == 2:
            urllib.request.urlretrieve(dna[0]+dna[1], filename=dna[1])
            dna = read_dnafile(dna[1])
        elif isinstance(dna,str):
            pass
        else:
            raise TypeError('DNA file is not string',type(dna))
        self.dna=dna
        
        #determine exons
        if isinstance(exon_regions,(list,tuple,str)) and len(exon_regions)==2:
            urllib.request.urlretrieve(exon_regions[0]+exon_regions[1], filename=exon_regions[1])
            exon_regions=read_exon_regions(exon_regions[1])
        elif isinstance(exon_regions, (list,tuple,int)):
            pass
        else:
            raise TypeError('Exon regions file is not tuple', \
                            type(exon_regions))

        self.exon_regions=exon_regions
        exons=[]
        for start,stop in exon_regions:
            exons.append(dna[start-1:stop])
        self.exons=exons  

class Protein(Gene):
    def __init__(self,dna,exons,trans_code,N,markov_chain):
        Gene.__init__(self,dna,exons)

        #determine translation code
        if self.exons is not None:
            mrna=create_mRNA(self.exons)
            self.mrna=mrna
        else:
            raise ValueError('Exon regions not found')
        
        if isinstance(trans_code,(list,tuple,str)) and len(trans_code) == 2:
            urllib.request.urlretrieve(trans_code[0]+trans_code[1], filename=trans_code[1])
            genetic_code=read_genetic_code(trans_code[1])
        elif isinstance(trans_code,(list,tuple,str)) and len(trans_code)==5:
            pass
        else:
            raise TypeError('Translation code file is not tuple',type(trans_code))
        self.genetic_code=genetic_code

    #mutate cDNA
    def introduce_mutations(self,markov_chain):
        if markov_chain is None:
            self.mrna=mutate_dna_v1(self.exons)
        else:
            self.mrna=mutate_dna_v2(self.exons,markov_chain)

    def get_product(self):        
        return create_protein(self.mrna,self.genetic_code)










