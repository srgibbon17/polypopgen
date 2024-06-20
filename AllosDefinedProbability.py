import random
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import copy

"""Single Locus Model for Allopolyploids or Segmental Allopolyploids

This script allows for simulating a single locus in a population of allopolyploids or segmental allopolyploids.
It operates largely on a series of random choices for selection, mutation, and gamete formation. 
Specifically, this script uses derived gametic probabilities for differing numbers of homoeologous exchanges (HEs).

Basic demographic features are included such as population shifts or single bottlenecks. 
Currently, migration, population splits, additional bottlenecks, and other demography features are not included. 

A set of input parameters are described below using inline comments

Of note here is that selection is only modeled to be purifying and h1, h2, h3 are often set to be additive. 
In other words, h2 = 2h1 and h3 = 3h1 such that each allelic copy has an equivalent effect size.

Additionally, if the self_rate is set to zero, selfing will still occur with probability 1/N per individual per generation. 

For an allopolyploid with strict bivalent formation, set he0 = 1 and he1 = he2 = 0. 
^ This may be recoded in terms of r and b where r is the recombination rate and b is the probability of quadrivalent formation.
"""

# Let's define some parameters
pa = .4  # allele frequency for a (derived, deleterious allele) in subgenome a (0<x<1)  
qa = 1-pa  # allele frequency for A in subgenome a
pb = .7  # allele frequency for a in subgenome b (0<x<1) 
qb = 1-pb  # allele frequency for A in subgenome b
N = 1000  # population size
g = 100  # number of generations
g_bottleneck_start = 400  # generation at which bottleneck starts
g_bottleneck_length = 200  # number of generations for which the bottleneck lasts
N_bottleneck = 1000  # population size during the bottleneck
s = .01  # selection coefficient (0<x<1) 
h1 = .25  # dominance coefficient for G3 (0<x<1) 
h2 = .5  # dominance coefficient for G2 (0<x<1) 
h3 = .75  # dominance coefficient for G1 (0<x<1) 
mu = 0  # mutation rate (from 'A' to 'a') (0<x<1) 
nu = 0  # mutation rate (from 'a' to 'A') (0<x<1) 
self_rate = 0  # selfing probability (0<x<1) 
# note: he0, he1, and he2 can either be inputted as probabilities or relative frequencies
# e.g. probabilities [.5, .3, .2] are equivalent to relative frequencies [5, 3, 2]
he0 = 1  # probability of no HEs 
he1 = 0  # probability of 1 HE
he2 = 0  # probability of 2 HEs

plot_guidelines = False  # if true, plots expected genotype and gamete frequencies under 0 HEs


# Let's create some functions

# mutates pre-selected gametes
def mutation_bidirectional(maternal_gamete, paternal_gamete):

    """
    Mutates gametes both forwards ('A' to 'a') and backwards ('a' to 'A').
    The forwards rate is given as a variable 'mu' and the backwards rate is given as 'nu'

    Args:
        maternal_gamete: an unmutated gamete from the mother
        paternal_gamete: an unmutated gamete from the father

    Returns:
        maternal_gamete_mut: mutated gamete from mother
        paternal_gamete_mut: mutated gamete from father
    """

    global maternal_gamete_mut
    maternal_gamete_mut = list(copy.copy(maternal_gamete))
    global paternal_gamete_mut
    paternal_gamete_mut = list(copy.copy(paternal_gamete))

    # for each allele in the maternal gamete, mutates it with likelihoods mu and nu
    for l in range(len(maternal_gamete)):  
        mut_event_1 = random.choices([True, False], [mu, 1-mu], k=1)  # randomly sets mutation to occur (i.e. True) with likelihood mu
        mut_event_2 = random.choices([True, False], [nu, 1-nu], k=1)  # randomly sets mutation to occur (i.e. True) with likelihood nu
        if (mut_event_1[0] == True) and (maternal_gamete[l] == 'A'):  
            maternal_gamete_mut[l] = 'a'
        elif (mut_event_2[0] == True) and (maternal_gamete[l] == 'a'):
            maternal_gamete_mut[l] = 'A'
    # same as above, but for the paternal gamete
    for l in range(len(paternal_gamete)):  
        mut_event_1 = random.choices([True, False], [mu, 1-mu], k=1) 
        mut_event_2 = random.choices([True, False], [nu, 1-nu], k=1)   
        if (mut_event_1[0] == True) and (paternal_gamete[l] == 'A'):
            paternal_gamete_mut[l] = 'a'
        elif (mut_event_2[0] == True) and (paternal_gamete[l] == 'a'):
            paternal_gamete_mut[l] = 'A'


# counts relative gamete frequencies
def gamete_frequencies(gametes_j, mutation_boolean):
    
    """
    Counts relative gamete frequencies (g00, g01, g10, g11 scaled to N_j) for one generation

    Args:
        gamete_j: a list of all gametes used to form the jth generation
        mutation_boolean: True if mutation has occurred, False otherwise

    Returns:
        g00_freq: gamete freq. (pre or post mutation) of g00 
        g10_freq: gamete freq. (pre or post mutation) of g10
        g01_freq: gamete freq. (pre or post mutation) of g01
        g11_freq: gamete freq. (pre or post mutation) of g11
    """

    g00 = gametes_j.count(['a', 'a'])/(2*N_j)
    g10 = gametes_j.count(['A', 'a'])/(2*N_j)
    g01 = gametes_j.count(['a', 'A'])/(2*N_j)
    g11 = gametes_j.count(['A', 'A'])/(2*N_j)

    if mutation_boolean == True:
        g00_freq_post_mut.append(g00)
        g10_freq_post_mut.append(g10)
        g01_freq_post_mut.append(g01)
        g11_freq_post_mut.append(g11)
    elif mutation_boolean == False:
        g00_freq_pre_mut.append(g00)
        g10_freq_pre_mut.append(g10)
        g01_freq_pre_mut.append(g01)
        g11_freq_pre_mut.append(g11)

# counts relative genotype frequencies
def genotype_frequencies(gen_j_pre_sel):
   
    """
    Counts relative genotype frequencies 
        (G00, G01, G02,
         G10, G11, G12,
         G20, G21, G22)
    (all of which are scaled to N_j) for one generation

    Args:
        gen_j_pre_sel: a list of all individuals in a generation before selection

    Returns:
        G00_freq: genotype freq. of G00 
            etc. for G01, G02, G10, G11, G12, G20, G21, and G22
        
        G1_freq: genotype freq. of G1, equivalent to G01_freq + G10_freq
        G2_freq: genotype freq. of G2, equivalent to G02_freq + G11_freq + G20_freq
        G3_freq: genotype freq. of G3, equivalent to G021_freq + G12_freq
    """

    G00 = 0
    G01 = 0
    G02 = 0
    G10 = 0
    G11 = 0
    G12 = 0
    G20 = 0
    G21 = 0
    G22 = 0
    for i in range(len(gen_j_pre_sel)):
        allele_count_sga = gen_j_pre_sel[i][:2].count('A')
        allele_count_sgb = gen_j_pre_sel[i][2:].count('A')
        if (allele_count_sga == 0) and (allele_count_sgb == 0):
            G00 += 1
        elif (allele_count_sga == 0) and (allele_count_sgb == 1):
            G01 += 1
        elif (allele_count_sga == 0) and (allele_count_sgb == 2):
            G02 += 1
        elif (allele_count_sga == 1) and (allele_count_sgb == 0):
            G10 += 1
        elif (allele_count_sga == 1) and (allele_count_sgb == 1):
            G11 += 1
        elif (allele_count_sga == 1) and (allele_count_sgb == 2):
            G12 += 1
        elif (allele_count_sga == 2) and (allele_count_sgb == 0):
            G20 += 1
        elif (allele_count_sga == 2) and (allele_count_sgb == 1):
            G21 += 1
        elif (allele_count_sga == 2) and (allele_count_sgb == 2):
            G22 += 1
        
    G00_freq.append(G00/N_j)
    G01_freq.append(G01/N_j)
    G02_freq.append(G02/N_j)
    G10_freq.append(G10/N_j)
    G11_freq.append(G11/N_j)
    G12_freq.append(G12/N_j)
    G20_freq.append(G20/N_j)
    G21_freq.append(G21/N_j)
    G22_freq.append(G22/N_j)
    
    G1 = G01 + G10
    G2 = G20 + G11 + G02
    G3 = G21 + G12

    G1_freq.append(G1/N_j)
    G2_freq.append(G2/N_j)
    G3_freq.append(G3/N_j)

# forms gametes for next generation
def gamete_formation(genotype, parental_relation): 

    """
    Forms maternal and paternal gametes for the next generation given a parental genotype input. 

    Args:
        genotype: a given genotype from the list of possible genotypes 
            [G00, G01, G02, G10, G11, G12, G20, G21, G22]
        parental_relation: an input of the form 'maternal' or 'paternal' 

    Returns:
        
    """

    global maternal_gamete
    global paternal_gamete

    # gamete formation with probability distribution under no HEs
    def gamete_formation_HE0():
        global gamete
        if genotype[0] == 'G00':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 0, 0, 1], k=1)
        elif genotype[0] == 'G01':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 1, 0, 1], k=1)
        elif genotype[0] == 'G02':    
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 1, 0, 0], k=1)
        
        elif genotype[0] == 'G10':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 0, 1, 1], k=1)
        elif genotype[0] == 'G11':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 1, 1, 1], k=1)
        elif genotype[0] == 'G12':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 1, 0, 0], k=1)
        
        elif genotype[0] == 'G20':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 0, 1, 0], k=1)
        elif genotype[0] == 'G21':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 0, 1, 0], k=1)
        elif genotype[0] == 'G22':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 0, 0, 0], k=1)
            
    # gamete formation with probability distribution under one HE
    # special note for the HE1 case: an additional split is made across the type of recombination event  
    # (ex. Aa || Aa has 4 possible recombinations: A||A, a||a, A||a, and a||A and thus 4 distinct sets of gametic probabilities)
    def gamete_formation_HE1():
        global gamete

        if genotype[0] == 'G00':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 0, 0, 1], k=1)
        elif genotype[0] == 'G01':
            exchange_pattern = random.choices(['a||a', 'a||A'], [.5, .5], k=1)
            if exchange_pattern[0] == 'a||a':
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 1, 0, 1], k=1)
            elif exchange_pattern[0] == 'a||A':
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 3, 3, 9], k=1)
        elif genotype[0] == 'G02':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [3, 9, 1, 3], k=1)

        elif genotype[0] == 'G10':
            exchange_pattern = random.choices(['a||a', 'A||a'], [.5, .5], k=1)
            if exchange_pattern[0] == 'a||a':
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 0, 1, 1], k=1)
            elif exchange_pattern[0] == 'A||a':
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 3, 3, 9], k=1)
        elif genotype[0] == 'G11':
            exchange_pattern = random.choices(['A||A', 'a||a', 'A||a', 'a||A'],[.25, .25, .25, .25], k=1)
            if (exchange_pattern[0] == 'A||A') or (exchange_pattern[0] == 'a||a'):
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 1, 1, 1], k=1)
            elif exchange_pattern[0] == 'A||a':
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [3, 9, 1, 3], k=1)
            elif exchange_pattern[0] == 'a||A':
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [3, 1, 9, 3], k=1)
        elif genotype[0] == 'G12':
            exchange_pattern = random.choices(['A||A', 'a||A'], [.5, .5], k=1)
            if exchange_pattern[0] == 'A||A':
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 1, 0, 0], k=1)
            elif exchange_pattern[0] == 'a||A':
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [9, 3, 3, 1], k=1)

        elif genotype[0] == 'G20':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [3, 1, 9, 3], k=1)
        elif genotype[0] == 'G21':
            exchange_pattern = random.choices(['A||A', 'A||a'], [.5, .5], k=1)
            if exchange_pattern[0] == 'A||A':    
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 0, 1, 0], k=1)
            elif exchange_pattern[0] == 'A||a':
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [9, 3, 3, 1], k=1)
        elif genotype[0] == 'G22':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 0, 0, 0], k=1)

    # gamete formation with probability distribution under two HEs
    def gamete_formation_HE2():
        global gamete

        if genotype[0] == 'G00':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 0, 0, 1], k=1)
        elif genotype[0] == 'G01':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 3, 3, 9], k=1)
        elif genotype[0] == 'G02':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 1, 1, 1], k=1)

        elif genotype[0] == 'G10':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 3, 3, 9], k=1)
        elif genotype[0] == 'G11':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 1, 1, 1], k=1)
        elif genotype[0] == 'G12':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [9, 3, 3, 1], k=1)
        
        elif genotype[0] == 'G020':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 1, 1, 1], k=1)
        elif genotype[0] == 'G21':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [9, 3, 3, 1], k=1)
        elif genotype[0] == 'G22':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 0, 0, 0], k=1)
            
    
    # randomly sets the number of HEs to be 0, 1, or 2 based on the specified parameters
    number_HEs = random.choices([0, 1, 2], [he0, he1, he2], k=1)  
    
    if number_HEs[0] == 0:
        gamete_formation_HE0()
    elif number_HEs[0] == 1:
        gamete_formation_HE1()
    elif number_HEs[0] == 2:
        gamete_formation_HE2()


    if parental_relation == 'maternal':
        maternal_gamete = gamete
    elif parental_relation == 'paternal':
        paternal_gamete = gamete

def j_plus_1th_generation(): 
    for i in range(N_j):    

        selfing_event = random.choices([True, False], [self_rate, 1-self_rate], k=1)

        possible_genotypes = ['G00', 'G01', 'G02', 'G10', 'G11', 'G12', 'G20', 'G21', 'G22']

        if selfing_event[0] == True:
            # randomly samples one individual to be the mother; of form [['A','a','A','A']]
            parental_genotype = random.choices(possible_genotypes, pop_genotype_densities, k=1)  
            
            # for the maternal gamete (and later the paternal gamete), calls corresponding function given the number of HEs
            gamete_formation(parental_genotype, 'maternal')
            gamete_formation(parental_genotype, 'paternal')

        elif selfing_event[0] == False: 
            # randomly samples two parents (with replacement); parents are of form [['A','a','A','A'],['A','A','a','a']]
            maternal_genotype = random.choices(possible_genotypes, pop_genotype_densities, k=1)  
            paternal_genotype = random.choices(possible_genotypes, pop_genotype_densities, k=1)

            # for the maternal gamete (and later the paternal gamete), calls corresponding function given the number of HEs
            gamete_formation(maternal_genotype, 'maternal')
            gamete_formation(paternal_genotype, 'paternal')

        # append the maternal and paternal gametes to a list of all gametes BEFORE mutation for the jth generation
        gametes_j_pre_mut.append(maternal_gamete[0])  
        gametes_j_pre_mut.append(paternal_gamete[0])

        # runs mutation function to act on the ith set of gametes in the jth generation
        # the output (i.e. mutated gametes) are maternal/paternal_gamete_mut
        mutation_bidirectional(maternal_gamete[0], paternal_gamete[0])

        # append the maternal and paternal gametes to a list of all gametes AFTER mutation for the jth generation
        gametes_j_post_mut.append(maternal_gamete_mut)
        gametes_j_post_mut.append(paternal_gamete_mut)          
        
        # appends the alleles from the two lines above to be the ith individual in the jth generation
        """This is formatted as such because maternal_gamete_mut is (essentially) of the form ['A', 'b'] 
        where the first position indicates a dominant (A) or recessive (a) allele from the 'a' subgenome
        and the second position indicates a dominant (B) or recessive (b) allele from the 'b' subgenome.
        
        Although this is not precisely the case here, as A and a are used for both the 'a' and 'b' subgenomes,
        this distinction still holds significance for how the gamete sampling functions operate. 
        i.e. the gamete sampling operates differently on ['A','A','a','a'] than on ['A','a','A','a'] despite
        these genotypes having equivalent fitnesses and overall number of alleles. In other words, the distribution of
        alleles across subgenomes matters. 
        
        Thus, the order of the subgenomes is preserved with this approach, and would produce something of the form ['A','A','b','B']

        Note: this order is still accurate under HEs because double reduction cannot affect alleles near the centromere
        and here, the subgenome is defined based on the chromosomal region near the centromere 
        """
        individual = [maternal_gamete_mut[0], paternal_gamete_mut[0], maternal_gamete_mut[1], paternal_gamete_mut[1]]
        
        gen_j_pre_sel.append(individual)

# forms gametes
# this is composed of functions which are split across 0, 1, or 2 HEs 
# note: uses derived probability distributions for gametes under each number of HEs (not random processes)
def old_gamete_formation(a_allele_count, b_allele_count, k): 

    global maternal_gamete
    global paternal_gamete

    # gamete formation with probability distribution under no HEs
    def gamete_formation_HE0():
        global gamete
        if a_allele_count == 0:
            if b_allele_count == 0:  # G00
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 0, 0, 0], k=1)
            elif b_allele_count == 1:  # G01
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 0, 1, 0], k=1)
            elif b_allele_count == 2:  # G02
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 0, 1, 0], k=1)
        
        elif a_allele_count == 1:
            if b_allele_count == 0:  # G10
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 1, 0, 0], k=1)
            elif b_allele_count == 1:  # G11
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 1, 1, 1], k=1)
            elif b_allele_count == 2:  # G12
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 0, 1, 1], k=1)
        
        elif a_allele_count == 2:
            if b_allele_count == 0:  # G20
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 1, 0, 0], k=1)
            elif b_allele_count == 1:  #G21
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 1, 0, 1], k=1)
            elif b_allele_count == 2:  # G22
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 0, 0, 1], k=1)
            
    # gamete formation with probability distribution under one HE
    # special note for the HE1 case: an additional split is made across the type of recombination event  
    # (ex. Aa || Aa has 4 possible recombinations: A||A, a||a, A||a, and a||A and thus 4 distinct sets of gametic probabilities)
    def gamete_formation_HE1():
        global gamete

        if a_allele_count == 0:
            if b_allele_count == 0:  # G00
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 0, 0, 0], k=1)
            elif b_allele_count == 1:  # G01
                exchange_pattern = random.choices(['A||A', 'A||a'], [.5, .5], k=1)
                if exchange_pattern[0] == 'A||A':    
                    gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 0, 1, 0], k=1)
                elif exchange_pattern[0] == 'A||a':
                    gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [9, 3, 3, 1], k=1)
            elif b_allele_count == 2:  # G02
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [3, 1, 9, 3], k=1)

        elif a_allele_count == 1:
            if b_allele_count == 0:  # G10
                exchange_pattern = random.choices(['A||A', 'a||A'], [.5, .5], k=1)
                if exchange_pattern[0] == 'A||A':
                    gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 1, 0, 0], k=1)
                elif exchange_pattern[0] == 'a||A':
                    gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [9, 3, 3, 1], k=1)
            elif b_allele_count == 1:  # G11
                exchange_pattern = random.choices(['A||A', 'a||a', 'A||a', 'a||A'],[.25, .25, .25, .25], k=1)
                if (exchange_pattern[0] == 'A||A') or (exchange_pattern[0] == 'a||a'):
                    gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 1, 1, 1], k=1)
                elif exchange_pattern[0] == 'A||a':
                    gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [3, 9, 1, 3], k=1)
                elif exchange_pattern[0] == 'a||A':
                    gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [3, 1, 9, 3], k=1)
            elif b_allele_count == 2:  # G12
                exchange_pattern = random.choices(['a||a', 'A||a'], [.5, .5], k=1)
                if exchange_pattern[0] == 'a||a':
                    gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 0, 1, 1], k=1)
                elif exchange_pattern[0] == 'A||a':
                    gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 3, 3, 9], k=1)

        elif a_allele_count == 2:
            if b_allele_count == 0:  # G20
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [3, 9, 1, 3], k=1)
            elif b_allele_count == 1:  # G21
                exchange_pattern = random.choices(['a||a', 'a||A'], [.5, .5], k=1)
                if exchange_pattern[0] == 'a||a':
                    gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 1, 0, 1], k=1)
                elif exchange_pattern[0] == 'a||A':
                    gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 3, 3, 9], k=1)
            elif b_allele_count == 2:  # G22
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 0, 0, 1], k=1)

    # gamete formation with probability distribution under two HEs
    def gamete_formation_HE2():
        global gamete
        if a_allele_count == 0:
            if b_allele_count == 0: # G00
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 0, 0, 0], k=1)
            elif b_allele_count == 1: # G01
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [9, 3, 3, 1], k=1)
            elif b_allele_count == 2:  # G02
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 1, 1, 1], k=1)

        elif a_allele_count == 1:
            if b_allele_count == 0:  # G10
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [9, 3, 3, 1], k=1)
            elif b_allele_count == 1:  # G11
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 1, 1, 1], k=1)
            elif b_allele_count == 2:  # G12
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 3, 3, 9], k=1)
        
        elif a_allele_count == 2:
            if b_allele_count == 0:  # G20
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 1, 1, 1], k=1)
            elif b_allele_count == 1: # G21
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 3, 3, 9], k=1)
            elif b_allele_count == 2:  # G22
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 0, 0, 1], k=1)
            
    
    # randomly sets the number of HEs to be 0, 1, or 2 based on the specified parameters
    number_HEs = random.choices([0, 1, 2], [he0, he1, he2], k=1)  
    
    if number_HEs[0] == 0:
        gamete_formation_HE0()
    elif number_HEs[0] == 1:
        gamete_formation_HE1()
    elif number_HEs[0] == 2:
        gamete_formation_HE2()


    if k == 'maternal':
        maternal_gamete = gamete
    elif k == 'paternal':
        paternal_gamete = gamete

# creates the j+1th generation of individuals
# workflow: samples gametes from the parents, mutates gametes with mutation_bidirectional, and
# creates individuals from mutated gametes which are added to the j+1th generation
def old_j_plus_1th_generation(): 
    for i in range(N_j):    

        selfing_event = random.choices([True, False], [self_rate, 1-self_rate], k=1)

        if selfing_event[0] == True:
            # randomly samples one individual to be the mother; of form [['A','a','A','A']]
            mother = random.sample(gen_j_post_sel, 1)  
            
            # counts the number of 'A' alleles in each subgenome
            sga_allele_count = mother[0][:2].count('a')
            sgb_allele_count = mother[0][2:].count('a') 
            
            # for the maternal gamete (and later the paternal gamete), calls corresponding function given the number of HEs
            gamete_formation(sga_allele_count, sgb_allele_count, 'maternal')
            gamete_formation(sga_allele_count, sgb_allele_count, 'paternal')

        elif selfing_event[0] == False: 
            # randomly samples two parents (with replacement); parents are of form [['A','a','A','A'],['A','A','a','a']]
            parents = random.choices(gen_j_post_sel, k=2)  

            # counts the number of 'A' alleles in each subgenome of each parent
            maternal_sga_allele_count = parents[0][:2].count('a')  
            maternal_sgb_allele_count = parents[0][2:].count('a') 
            paternal_sga_allele_count = parents[1][:2].count('a') 
            paternal_sgb_allele_count = parents[1][2:].count('a') 

            # for the maternal gamete (and later the paternal gamete), calls corresponding function given the number of HEs
            gamete_formation(maternal_sga_allele_count, maternal_sgb_allele_count, 'maternal')
            gamete_formation(paternal_sga_allele_count, paternal_sgb_allele_count, 'paternal')

        # append the maternal and paternal gametes to a list of all gametes BEFORE mutation for the jth generation
        gametes_j_pre_mut.append(maternal_gamete[0])  
        gametes_j_pre_mut.append(paternal_gamete[0])

        # runs mutation function to act on the ith set of gametes in the jth generation
        # the output (i.e. mutated gametes) are maternal/paternal_gamete_mut
        mutation_bidirectional(maternal_gamete[0], paternal_gamete[0])

        # append the maternal and paternal gametes to a list of all gametes AFTER mutation for the jth generation
        gametes_j_post_mut.append(maternal_gamete_mut)
        gametes_j_post_mut.append(paternal_gamete_mut)          
        
        # appends the alleles from the two lines above to be the ith individual in the jth generation
        """This is formatted as such because maternal_gamete_mut is (essentially) of the form ['A', 'b'] 
        where the first position indicates a dominant (A) or recessive (a) allele from the 'a' subgenome
        and the second position indicates a dominant (B) or recessive (b) allele from the 'b' subgenome.
        
        Although this is not precisely the case here, as A and a are used for both the 'a' and 'b' subgenomes,
        this distinction still holds significance for how the gamete sampling functions operate. 
        i.e. the gamete sampling operates differently on ['A','A','a','a'] than on ['A','a','A','a'] despite
        these genotypes having equivalent fitnesses and overall number of alleles. In other words, the distribution of
        alleles across subgenomes matters. 
        
        Thus, the order of the subgenomes is preserved with this approach, and would produce something of the form ['A','A','b','B']

        Note: this order is still accurate under HEs because double reduction cannot affect alleles near the centromere
        and here, the subgenome is defined based on the chromosomal region near the centromere 
        """
        individual = [maternal_gamete_mut[0], paternal_gamete_mut[0], maternal_gamete_mut[1], paternal_gamete_mut[1]]
        
        gen_j_pre_sel.append(individual)



# purifying selection function 
# workflow: randomly chooses "True" or "False" for a selection event
# the respective frequencies for "True" are 0 (A, A), s*h1, s*h2, or s*h3 (heterozygotes), s (a, a)
# "True" represents that purifying selection occurs and (effectively) removes the individual
# "False" represents that the individual survives to reproduce and is kept in gen_j_post_sel
def selection_survivability(gen_j_pre_sel, gen_j_post_sel):
    for i in range(len(gen_j_pre_sel)):
        alleles = gen_j_pre_sel[i].count('A')
        # selection for recessive (G0) homozygotes - a2 a2 b2 b2 (G00) 
        if alleles == 0:
            sel_event_G0 = random.choices([True, False], [s, 1-s], k=1)
            if sel_event_G0[0] == False:
                gen_j_post_sel.append(gen_j_pre_sel[i])       
        # selection for G1 heterozygotes - a1 a2 a2 a2 (G10) or a2 a2 b1 b2 (G01) 
        elif alleles == 1:
            sel_event_G1 = random.choices([True, False], [h3*s, 1-h3*s], k=1)
            if sel_event_G1[0] == False:
                gen_j_post_sel.append(gen_j_pre_sel[i])
        # selection for G2 heterozygotes - a1 a1 b2 b2 (G20) or a2 a2 b1 b1 (G02) or a1 a2 b1 b2 (G11)
        elif alleles == 2:
            sel_event_G2 = random.choices([True, False], [h2*s, 1-h2*s], k=1)
            if sel_event_G2[0] == False:
                gen_j_post_sel.append(gen_j_pre_sel[i])
        # selection for G3 heterozygotes - a1 a1 a1 a2 (G21) or a1 a2 b1 b1 (G12)
        elif alleles == 3:
            sel_event_G3 = random.choices([True, False], [h1*s, 1-h1*s], k=1)
            if sel_event_G3[0] == False:
                gen_j_post_sel.append(gen_j_pre_sel[i])
        # selection for dominant (G4) homozygotes - a1 a1 a1 a1 (G22)
        elif alleles == 4:
            gen_j_post_sel.append(gen_j_pre_sel[i])

# purifying selection function
# 
def selection_fecundity(gen_j_pre_sel):
    
    global pop_genotype_densities

    G00 = 0
    G01 = 0
    G02 = 0
    G10 = 0
    G11 = 0
    G12 = 0
    G20 = 0
    G21 = 0
    G22 = 0

    for i in range(len(gen_j_pre_sel)):
        allele_count_sga = gen_j_pre_sel[i][:2].count('A')
        allele_count_sgb = gen_j_pre_sel[i][2:].count('A')
        if (allele_count_sga == 0) and (allele_count_sgb == 0):
            G00 += 1
        elif (allele_count_sga == 0) and (allele_count_sgb == 1):
            G01 += 1
        elif (allele_count_sga == 0) and (allele_count_sgb == 2):
            G02 += 1
        elif (allele_count_sga == 1) and (allele_count_sgb == 0):
            G10 += 1
        elif (allele_count_sga == 1) and (allele_count_sgb == 1):
            G11 += 1
        elif (allele_count_sga == 1) and (allele_count_sgb == 2):
            G12 += 1
        elif (allele_count_sga == 2) and (allele_count_sgb == 0):
            G20 += 1
        elif (allele_count_sga == 2) and (allele_count_sgb == 1):
            G21 += 1
        elif (allele_count_sga == 2) and (allele_count_sgb == 2):
            G22 += 1

    genotype_freq_gen_j = np.array([G00, G01, G02, G10, G11, G12, G20, G21, G22])
    genotype_proportions_gen_j = genotype_freq_gen_j/len(gen_j_pre_sel)

    absolute_fitnesses = np.array([1, 1-s*h1, 1-s*h2, 1-s*h1, 1-s*h2, 1-s*h3, 1-s*h2, 1-s*h3, 1-s])
    relative_fitnesses = absolute_fitnesses / (np.sum(np.multiply(genotype_proportions_gen_j, absolute_fitnesses)))

    pop_genotype_densities = np.multiply(relative_fitnesses, genotype_proportions_gen_j)

### Define a few variables and lists to store information
gen0 = []

# initializes variables to hold the gamete frequencies for all generations 
# these lists are delineated across both pre and post mutation as well as g00, g01, g10, and g11
# pre mutation
g00_freq_pre_mut = []
g10_freq_pre_mut = []
g01_freq_pre_mut = [] 
g11_freq_pre_mut = []  
# post mutation
g00_freq_post_mut = []
g10_freq_post_mut = []
g01_freq_post_mut = []  
g11_freq_post_mut = []  

# initializes variables to hold the genotype frequencies for all generations
G00_freq = []
G01_freq = []
G02_freq = []
G10_freq = []
G11_freq = []
G12_freq = []
G20_freq = []
G21_freq = []
G22_freq = []

G1_freq = []
G2_freq = []
G3_freq = []

# creates the initial generation of N individuals using random sampling with replacement
for i in range(N):  
    subgenome_a = random.choices(['A', 'a'], [pa, 1-pa], k=2)
    subgenome_b = random.choices(['A', 'a'], [pb, 1-pb], k=2)
    gen0.append(subgenome_a + subgenome_b)

# set the pre-selection generation to initially be the 0th generation created immediately above
gen_j_pre_sel = gen0  

for j in range(g):  
    
    # additional variables and lists

    # a post-selection generation (those that survive to reproduce)
    gen_j_post_sel = []  

    # a pre-mutation list of gametes
    gametes_j_pre_mut = [] 
    # a post-mutation list of gametes
    gametes_j_post_mut = []  

    # sets N_j to either be the population size (N) or the population size for the bottleneck (N_bottleneck)
    if j in range(g_bottleneck_start, g_bottleneck_start + g_bottleneck_length):  
        N_j = N_bottleneck
    else:
        N_j = N

    # runs selection function to act on the jth generation
    selection_fecundity(gen_j_pre_sel)
    # empties pre-selection genetation to be filled using the j+1th generation immediately below
    gen_j_pre_sel = []

    # creates the next generation (j+1th generation) of individuals of size N_j
    j_plus_1th_generation()          

    # counts and stores gamete frequencies pre-mutation
    gamete_frequencies(gametes_j_pre_mut, False)

    # counts and stores gamete frequencies post-mutation
    gamete_frequencies(gametes_j_post_mut, True) 
    
    # counts and stores genotype frequencies 
    genotype_frequencies(gen_j_pre_sel)

### Plotting Functions

# transitions lists of gamete frequencies to arrays so that the difference due to mutation can be calculated
def gamete_freq_as_arrays():
    global g00_freq_pre_mut
    global g10_freq_pre_mut
    global g01_freq_pre_mut
    global g11_freq_pre_mut

    global g00_freq_post_mut
    global g10_freq_post_mut
    global g01_freq_post_mut
    global g11_freq_post_mut

    global g00_freq_diff
    global g10_freq_diff
    global g01_freq_diff
    global g11_freq_diff

    g00_freq_pre_mut = np.array(g00_freq_pre_mut)
    g00_freq_post_mut = np.array(g00_freq_post_mut)
    g00_freq_diff = g00_freq_post_mut - g00_freq_pre_mut

    g10_freq_pre_mut = np.array(g10_freq_pre_mut)
    g10_freq_post_mut = np.array(g10_freq_post_mut)
    g10_freq_diff = g10_freq_post_mut - g10_freq_pre_mut

    g01_freq_pre_mut = np.array(g01_freq_pre_mut)
    g01_freq_post_mut = np.array(g01_freq_post_mut)
    g01_freq_diff = g01_freq_post_mut - g01_freq_pre_mut

    g11_freq_pre_mut = np.array(g11_freq_pre_mut)
    g11_freq_post_mut = np.array(g11_freq_post_mut)
    g11_freq_diff = g11_freq_post_mut - g11_freq_pre_mut


# plots both the gamete and genotype frequencies (each on its own axis)
# Sam's preferred plotting function/style
def plot_gamete_and_genotype_freq():
    
    x_index = range(g)  # creates a discrete time variable to plot against
    fig, axs = plt.subplots(2, 1, sharex = 'col', sharey = 'row')  

    fig.set_size_inches(14, 8) 
    axs[0].plot(x_index, g00_freq_pre_mut, color = '#008aff')
    axs[0].plot(x_index, g10_freq_pre_mut, color = '#8800da')
    axs[0].plot(x_index, g01_freq_pre_mut, color = '#00a50e')
    axs[0].plot(x_index, g11_freq_pre_mut, color = '#da1b00')
    axs[0].legend(['g00', 'g10', 'g01', 'g11'])
    #axs[0].add_artist(gamete_legend)
    axs[0].set_ylabel('Gamete Frequencies')

    axs[1].plot(x_index, G00_freq, color = '#008aff')
    axs[1].plot(x_index, G1_freq, color = '#8800da')
    axs[1].plot(x_index, G2_freq, color = '#eb7f00')
    axs[1].plot(x_index, G3_freq, color = '#00a50e')
    axs[1].plot(x_index, G22_freq, color = '#da1b00')
    axs[1].legend(['G0', 'G1', 'G2', 'G3', 'G4'])
    axs[1].set_ylabel('Genotype Frequencies')
    
    # plots initial expected gamete and genotype values
    if plot_guidelines == True:
        axs[0].axhline(y = ((pa+pb)/2)*((qa+qb)/2), color = '0', linestyle = 'dashed')
        axs[0].axhline(y = (9/11)*((pa+pb)/2)*((qa+qb)/2), color = '0', linestyle = 'dotted')
        axs[0].axhline(y = (2/3)*((pa+pb)/2)*((qa+qb)/2), color = '0', linestyle = (0, (1, 10)))
        #axs[0].legend(['HE0', 'HE1', 'HE2'])

        axs[0].axhline(y = qa*qb, color = '#008aff', linestyle = 'dashed')
        axs[0].axhline(y = pb*qa, color = '#00a50e', linestyle = 'dashed')
        axs[0].axhline(y = pa*qb, color = '#8800da', linestyle = 'dashed')
        axs[0].axhline(y = pa*pb, color = '#da1b00', linestyle = 'dashed')
        
        axs[1].axhline(y = qa**2*qb**2, color = '#008aff', linestyle = 'dashed')
        axs[1].axhline(y = 2*qa*qb*(pa*qb + pb*qa), color = '#8800da', linestyle = 'dashed')
        axs[1].axhline(y = pb**2*qa**2 + 4*pa*pb*qa*qb + pa**2*qb**2, color = '#eb7f00', linestyle = 'dashed')
        axs[1].axhline(y = 2*pa*pb*(pa*qb + pb*qa), color = '#00a50e', linestyle = 'dashed')
        axs[1].axhline(y = pa**2*pb**2, color = '#da1b00', linestyle = 'dashed')

    # labels
    fig.supxlabel('Time (in generations)')
    fig.suptitle('Gamete and Genotype Frequency Over Time')

    # creates a light gray shaded region for the population shift or bottleneck
    for ax in axs:
        shaded_region = Rectangle((g_bottleneck_start, 0), g_bottleneck_length, 1)
        shaded_region.set_color('0.85')
        ax.add_patch(shaded_region)
        ax.set_xlim(0, g-1)
        #ax.set_ylim(-.1, 1.1)

    plt.show()


# plots the three gamete frequencies (g00, g01/g10, and g11) 
# before mutation on three separate axes
def plot_gamete_freq():
    
    x_index = range(g)  # creates a discrete time variable to plot against
    fig, axs = plt.subplots(1, 4, sharex = 'col', sharey = 'row')  

    fig.set_size_inches(16, 6) 
    axs[0].plot(x_index, g00_freq_pre_mut)
    axs[0].set_xlabel('g00')
    axs[1].plot(x_index, g10_freq_pre_mut)
    axs[1].set_xlabel('g10')
    axs[2].plot(x_index, g01_freq_pre_mut)
    axs[2].set_xlabel('g01')
    axs[3].plot(x_index, g11_freq_pre_mut)
    axs[3].set_xlabel('g11')

    # labels
    fig.supxlabel('Time (in generations)')
    fig.suptitle('Gamete Frequency Over Time')
    fig.supylabel('Gamete frequency (pre-mut)')

    # creates a light gray shaded region for the population shift or bottleneck
    for ax in axs:
        shaded_region = Rectangle((g_bottleneck_start, 0), g_bottleneck_length, 1)
        shaded_region.set_color('0.85')
        ax.add_patch(shaded_region)
        ax.set_xlim(0, g-1)
        ax.set_ylim(-.1, 1.1)
    
    plt.show()
    

# plots the five genotype frequencies (G0, G1, G2, G3, and G4) on 5 axes in a 3 by 2 grid of 6 total axes
def plot_genotype_freq():

    x_index = range(g)  # creates a discrete time variable to plot against
    fig, axs = plt.subplots(2, 3, sharex = 'col', sharey = 'row')  

    fig.set_size_inches(16, 6)
    axs[0, 0].plot(x_index, G00_freq)
    axs[0, 0].set_xlabel('G0')
    axs[0, 1].plot(x_index, G1_freq)
    axs[0, 1].set_xlabel('G1')
    axs[0, 2].plot(x_index, G2_freq)
    axs[0, 2].set_xlabel('G2')
    axs[1, 0].plot(x_index, G3_freq)
    axs[1, 0].set_xlabel('G3')
    axs[1, 1].plot(x_index, G22_freq)
    axs[1, 1].set_xlabel('G4')

    # labels
    fig.supxlabel('Time (in generations)')
    fig.suptitle('Genotype Frequency Over Time')
    fig.supylabel('Genotype Frequency')

    # creates a light gray shaded region for the population shift or bottleneck
    for ax in axs[0]:
        shaded_region = Rectangle((g_bottleneck_start, 0), g_bottleneck_length, 1)
        shaded_region.set_color('0.85')
        ax.add_patch(shaded_region)
        ax.set_xlim(0, g-1)
    for ax in axs[1]:
        shaded_region = Rectangle((g_bottleneck_start, 0), g_bottleneck_length, 1)
        shaded_region.set_color('0.85')
        ax.add_patch(shaded_region)
        ax.set_xlim(0, g-1)
    for ax in axs[0]:
        ax.set_ylim(-.1, 1.1)
    for ax in axs[1]:
        ax.set_ylim(-.1, 1.1)

    plt.show()


# plots the three gamete frequencies (g00, g01/g10, and g11) 
# before mutation on the same axis
def plot_gamete_freq_single_axis():
    
    x_index = range(g)  # creates a discrete time variable to plot against
    fig, ax = plt.subplots(figsize=(16, 6))
    plt.axis([0, g-1, -.1, 1.1])  
 
    ax.plot(x_index, g00_freq_pre_mut, color = '#008aff')
    ax.plot(x_index, g10_freq_pre_mut, color = '#8800da')
    ax.plot(x_index, g01_freq_pre_mut, color = '#00a50e' )
    ax.plot(x_index, g11_freq_pre_mut, color = '#da1b00')
    ax.legend(['g00', 'g10', 'g01', 'g11'])
    
    # labels
    fig.supxlabel('Time (in generations)')
    fig.suptitle('Gamete Frequency Over Time')
    fig.supylabel('Gamete frequency (pre-mut)')

    # creates a light gray shaded region for the population shift or bottleneck
    shaded_region = Rectangle((g_bottleneck_start, 0), g_bottleneck_length, 1)
    shaded_region.set_color('0.85')
    ax.add_patch(shaded_region)

    plt.show()


# plots the five genotype frequencies (G0, G1, G2, G3, and G4) on a single axis
def plot_genotype_freq_single_axis():    
    
    x_index = range(g)  # creates a discrete time variable to plot against
    fig, ax = plt.subplots(figsize=(16, 6))
    plt.axis([0, g-1, -.1, 1.1])  
 
    ax.plot(x_index, G00_freq, color = '#008aff')
    ax.plot(x_index, G1_freq, color = '#00a50e')
    ax.plot(x_index, G2_freq, color = '#da1b00')
    ax.plot(x_index, G3_freq, color = '#8800da')
    ax.plot(x_index, G22_freq, color = '#eb7f00')
    ax.legend(['G0', 'G1', 'G2', 'G3', 'G4'])
    
    # labels
    fig.supxlabel('Time (in generations)')
    fig.suptitle('Genotype Frequency Over Time')
    fig.supylabel('Genotype frequency')

    # creates a light gray shaded region for the population shift or bottleneck
    shaded_region = Rectangle((g_bottleneck_start, 0), g_bottleneck_length, 1)
    shaded_region.set_color('0.85')
    ax.add_patch(shaded_region)

    plt.show()


# plots the gamete differential between before and after mutation on a single axis
# note: only useful for simulations with a 200 or fewer generations
def plot_gamete_diff():
    
    gamete_freq_as_arrays()

    x_index = range(g)  # creates a discrete time variable to plot against
    fig, ax = plt.subplots(figsize=(16, 6))  

    # labels
    fig.supxlabel('Time (in generations)')
    fig.suptitle('Gamete Differentials Over Time')
    fig.supylabel('Gamete Differential (post-mut - pre-mut)')

    ax.plot(x_index, g00_freq_diff, color = '#008aff')
    ax.plot(x_index, g10_freq_diff, color = '#8800da')
    ax.plot(x_index, g01_freq_diff, color = '#00a50e')
    ax.plot(x_index, g11_freq_diff, color = '#da1b00')
    ax.legend(['g00', 'g10', 'g01', 'g11'])

    # creates a light gray shaded region for the population shift or bottleneck
    shaded_region = Rectangle((g_bottleneck_start, -.25), g_bottleneck_length, .5)
    shaded_region.set_color('0.85')

    plt.show()


# plots all 9 genotypes as delineated across subgenomes... is very, very difficult to see 9 lines 
def plot_delineated_genotypes():
        
    x_index = range(g)  # creates a discrete time variable to plot against
    fig, ax = plt.subplots(figsize=(16, 6))
    plt.axis([0, g-1, -.1, 1.1])  
 
    ax.plot(x_index, G00_freq, color = '#008aff')
    ax.plot(x_index, G10_freq, color = '#00a545')
    ax.plot(x_index, G01_freq, color = '#53a500')
    ax.plot(x_index, G20_freq, color = '#da1b00')
    ax.plot(x_index, G11_freq, color = '#da0040')
    ax.plot(x_index, G02_freq, color = '#da009a')
    ax.plot(x_index, G21_freq, color = '#6d00da')
    ax.plot(x_index, G12_freq, color = '#a400da')
    ax.plot(x_index, G22_freq, color = '#eb7f00')
    ax.legend(['G00', 'G10', 'G01', 'G20', 'G11', 'G02', 'G21', 'G12', 'G22'])
    
    # labels
    fig.supxlabel('Time (in generations)')
    fig.suptitle('Genotype Frequency Over Time')
    fig.supylabel('Genotype frequency')

    # creates a light gray shaded region for the population shift or bottleneck
    shaded_region = Rectangle((g_bottleneck_start, 0), g_bottleneck_length, 1)
    shaded_region.set_color('0.85')
    ax.add_patch(shaded_region)

    plt.show()

plot_gamete_and_genotype_freq()
