import random
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import copy

# Let's define some parameters
pa = .5  # allele frequency for A in subgenome a (0<x<1)  
qa = 1-pa  # allele frequency for a in subgenome a
pb = .5  # allele frequency for A in subgenome b (0<x<1) 
qb = 1-pb  # allele frequency for a in subgenome b
N = 100  # population size
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

gen0 = []

for i in range(N):  
    subgenome_a = random.choices(['A', 'a'], [1-pa, pa], k=2)
    subgenome_b = random.choices(['A', 'a'], [1-pb, pb], k=2)
    gen0.append(subgenome_a + subgenome_b)


# forms gametes for the next generation
def gamete_sampling(genotype): 
    """
    Forms maternal and paternal gametes for the next generation given a parental genotype input. 

    Args:
        genotype: a given genotype from the list of possible genotypes 
            [G00, G01, G02, G10, G11, G12, G20, G21, G22]

    Returns:
        gamete: this will be either the maternal or paternal gamete, dependent upon input
            this is done using 3 sub-functions for 0 HEs, 1 HE, or 2 HEs
    """

    # gamete formation with probability distribution under no HEs
    def gamete_formation_HE0():
        """
        Samples a gamete under 0 HEs using pre-defined probability distributions

        Args:
            genotype (implicit): a given genotype from the list of possible genotypes 
                [G00, G01, G02, G10, G11, G12, G20, G21, G22]

        Returns:
            gamete: one of four possible gametes (i.e. [A, A], [A, a], [a, A], [a, a])
        """
        
        global gamete
        if genotype[0] == 'G00':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 0, 0, 0], k=1)
        elif genotype[0] == 'G01':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 0, 1, 0], k=1)
        elif genotype[0] == 'G02':    
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 0, 1, 0], k=1)
        
        elif genotype[0] == 'G10':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 1, 0, 0], k=1)
        elif genotype[0] == 'G11':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 1, 1, 1], k=1)
        elif genotype[0] == 'G12':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 0, 1, 1], k=1)
        
        elif genotype[0] == 'G20':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 1, 0, 0], k=1)
        elif genotype[0] == 'G21':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 1, 0, 1], k=1)
        elif genotype[0] == 'G22':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 0, 0, 1], k=1)
            
    # gamete formation with probability distribution under one HE
    # special note for the HE1 case: an additional split is made across the type of recombination event  
    # (ex. Aa || Aa has 4 possible recombinations: A||A, a||a, A||a, and a||A and thus 4 distinct sets of gametic probabilities)
    def gamete_formation_HE1():
        """
        Samples a gamete under 1 HE using pre-defined probability distributions. 
        Note: The probability distributions are also delineated across the type of recombination/HE event.

        Args:
            genotype (implicit): a given genotype from the list of possible genotypes 
                [G00, G01, G02, G10, G11, G12, G20, G21, G22]

        Returns:
            gamete: one of four possible gametes (i.e. [A, A], [A, a], [a, A], [a, a])
        """
        
        global gamete

        if genotype[0] == 'G00':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 0, 0, 0], k=1)
        elif genotype[0] == 'G01':
            exchange_pattern = random.choices(['A||A', 'A||a'], [.5, .5], k=1)
            if exchange_pattern[0] == 'A||A':    
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 0, 1, 0], k=1)
            elif exchange_pattern[0] == 'A||a':
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [9, 3, 3, 1], k=1)
        elif genotype[0] == 'G02':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [3, 1, 9, 3], k=1)

        elif genotype[0] == 'G10':
            exchange_pattern = random.choices(['A||A', 'a||A'], [.5, .5], k=1)
            if exchange_pattern[0] == 'A||A':
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 1, 0, 0], k=1)
            elif exchange_pattern[0] == 'a||A':
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [9, 3, 3, 1], k=1)
        elif genotype[0] == 'G11':
            exchange_pattern = random.choices(['A||A', 'a||a', 'A||a', 'a||A'],[.25, .25, .25, .25], k=1)
            if (exchange_pattern[0] == 'A||A') or (exchange_pattern[0] == 'a||a'):
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 1, 1, 1], k=1)
            elif exchange_pattern[0] == 'A||a':
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [3, 9, 1, 3], k=1)
            elif exchange_pattern[0] == 'a||A':
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [3, 1, 9, 3], k=1)
        elif genotype[0] == 'G12':
            exchange_pattern = random.choices(['a||a', 'A||a'], [.5, .5], k=1)
            if exchange_pattern[0] == 'a||a':
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 0, 1, 1], k=1)
            elif exchange_pattern[0] == 'A||a':
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 3, 3, 9], k=1)

        elif genotype[0] == 'G20':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [3, 9, 1, 3], k=1)
        elif genotype[0] == 'G21':
            exchange_pattern = random.choices(['a||a', 'a||A'], [.5, .5], k=1)
            if exchange_pattern[0] == 'a||a':
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 1, 0, 1], k=1)
            elif exchange_pattern[0] == 'a||A':
                gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 3, 3, 9], k=1)
        elif genotype[0] == 'G22':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 0, 0, 1], k=1)

    # gamete formation with probability distribution under two HEs
    def gamete_formation_HE2():
        """
        Samples a gamete under 2 HEs using pre-defined probability distributions

        Args:
            genotype (implicit): a given genotype from the list of possible genotypes 
                [G00, G01, G02, G10, G11, G12, G20, G21, G22]

        Returns:
            gamete: one of four possible gametes (i.e. [A, A], [A, a], [a, A], [a, a])
        """    

        global gamete

        if genotype[0] == 'G00':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 0, 0, 0], k=1)
        elif genotype[0] == 'G01':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [9, 3, 3, 1], k=1)
        elif genotype[0] == 'G02':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 1, 1, 1], k=1)

        elif genotype[0] == 'G10':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [9, 3, 3, 1], k=1)
        elif genotype[0] == 'G11':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 1, 1, 1], k=1)
        elif genotype[0] == 'G12':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 3, 3, 9], k=1)
        
        elif genotype[0] == 'G20':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 1, 1, 1], k=1)
        elif genotype[0] == 'G21':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [1, 3, 3, 9], k=1)
        elif genotype[0] == 'G22':
            gamete = random.choices([['A', 'A'],['a','A'],['A', 'a'],['a', 'a']], [0, 0, 0, 1], k=1)
            
    
    # randomly sets the number of HEs to be 0, 1, or 2 based on the specified parameters
    number_HEs = random.choices([0, 1, 2], [he0, he1, he2], k=1)  
    
    if number_HEs[0] == 0:
        gamete_formation_HE0()
    elif number_HEs[0] == 1:
        gamete_formation_HE1()
    elif number_HEs[0] == 2:
        gamete_formation_HE2()

    return gamete

# rescales offspring proportions based on fitness
def selection_fecundity(gen_j_pre_sel): 
    """
    Selection function which acts to create a weighting of offspring after selection acts. 
    Alternative version is selection_survivability (see below)

    Args:
        gen_j_pre_sel: a list of individuals before selection

    Returns:
        sel_weighted_genotype_densities: a fitness-scaled proportion of the number of individuals 
            that will reproduce in a given generation
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
        allele_count_sga = gen_j_pre_sel[i][:2].count('a')
        allele_count_sgb = gen_j_pre_sel[i][2:].count('a')
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

    sel_weighted_genotype_densities = np.multiply(relative_fitnesses, genotype_proportions_gen_j)

    return sel_weighted_genotype_densities

pop_densities = selection_fecundity(gen0)

def selfing_and_gamete_formation():
    """
    Gamete formation and selfing function which controls whether selfing occurs. 

    Args:
        self_rate: rate at which selfing is set to occur (implicit/parameter)
        pop_densities: selection scaled offspring production levels (defined elsewhere)

    Returns:
        gametes_list: a list of gametes for mutation to act on
    """
    
    # set a selfing_event to occur with probability = self_rate
    selfing_event = random.choices([True, False], [self_rate, 1-self_rate], k=1)

    possible_genotypes = ['G00', 'G01', 'G02', 'G10', 'G11', 'G12', 'G20', 'G21', 'G22']

    if selfing_event[0] == True:
        # randomly samples one individual to be the mother; of form [['A','a','A','A']]
        parental_genotype = random.choices(possible_genotypes, pop_densities, k=1)  
            
        # for the maternal gamete (and later the paternal gamete), calls corresponding function given the number of HEs
        maternal_gamete = gamete_sampling(parental_genotype)
        paternal_gamete = gamete_sampling(parental_genotype)

    elif selfing_event[0] == False: 
        # randomly samples two parents (with replacement); parents are of form [['A','a','A','A'],['A','A','a','a']]
        maternal_genotype = random.choices(possible_genotypes, pop_densities, k=1)  
        paternal_genotype = random.choices(possible_genotypes, pop_densities, k=1)

        # for the maternal gamete (and later the paternal gamete), calls corresponding function given the number of HEs
        maternal_gamete = gamete_sampling(maternal_genotype)
        paternal_gamete = gamete_sampling(paternal_genotype)

    gametes_list = [maternal_gamete[0], paternal_gamete[0]]
    
    return gametes_list

gametes = selfing_and_gamete_formation()
print(gametes)
