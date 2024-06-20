import random
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import copy

"""Single Locus Model for Autopolyploids under Strict Bivalent Formation

This script allows for simulating a single locus through time in a population of autopolyploids.
It operates largely on a series of random choices for selection, mutation, and gamete formation. 
Specifically, this script uses derived gametic probability frequencies under strict bivalent formation.

Basic demographic features are included such as population shifts or single bottlenecks. 
Currently, migration, population splits, additional bottlenecks, and other demography features are not included. 

A large set of input parameters are described below using inline comments

Of note here is that selection is only modeled to be purifying and h1, h2, h3 are often set to be additive. 
In other words, h2 = 2h1 and h3 = 3h1 such that each allelic copy has an equivalent effect size.

Additionally, if the self_rate is set to zero, selfing will still occur with probability 1/N per individual per generation.  
"""

# Let's define some parameters
p = .5  # allele frequency for A
q = 1-p  # allele frequency for a
N = 100  # population size
g = 100  # number of generations
g_bottleneck_start = 400  # generation at which bottleneck starts
g_bottleneck_length = 200  # number of generations for which the bottleneck lasts
N_bottleneck = 500  # population size during the bottleneck
s = .00001  # selection coefficient
h1 = .25  # dominance coefficient for G3
h2 = .5  # dominance coefficient for G2
h3 = .75  # dominance coefficient for G1
mu = 0  # mutation rate (from 'A' to 'a')
nu = 0  # mutation rate (from 'a' to 'A')
self_rate = .5  # selfing rate (what parameter/letter/variable to use for this?)

# Let's create some functions

# mutates pre-selected gametes
# mutation at rates mu (A to a) and nu (a to A) 
# gametes are copied to maintain a copy both before and after mutations
def mutation_bidirectional(maternal_gamete, paternal_gamete):

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


# counts gamete frequencies for a single generation
# note: gametes are stored in relative terms by dividing by 2N_j 
# these gametes are kept in a list of all gametes pre-mutation
def gamete_frequencies_pre_mut(gametes_j_pre_mut):
    
    g00 = gametes_j_pre_mut.count(['a', 'a'])/(2*N_j)
    g11 = gametes_j_pre_mut.count(['A', 'A'])/(2*N_j)
    g01 = 1 - g00 - g11

    g00_freq_pre_mut.append(g00)
    g01_freq_pre_mut.append(g01)
    g11_freq_pre_mut.append(g11)


# identical function to gamete_frequencies_pre_mut, but is for post-mutation
def gamete_frequencies_post_mut(gametes_j_post_mut):
    
    g00 = gametes_j_post_mut.count(['a', 'a'])/(2*N_j)
    g11 = gametes_j_post_mut.count(['A', 'A'])/(2*N_j)
    g01 = 1 - g00 - g11

    g00_freq_post_mut.append(g00)
    g01_freq_post_mut.append(g01)
    g11_freq_post_mut.append(g11)


# similar to gamete_frequencies, but counts genotype frequencies
# note: stores gamete frequencies both across subgenomes and for G0, G1, G2, G3, and G4
def genotype_frequencies(gen_j_pre_sel):
    G0 = 0
    G1 = 0
    G2 = 0
    G3 = 0
    G4 = 0
    for i in range(len(gen_j_pre_sel)):
        allele_count = gen_j_pre_sel[i].count('A')
        if allele_count == 0:
            G0 += 1
        elif allele_count == 1:
            G1 += 1
        elif allele_count == 2:
            G2 += 1
        elif allele_count == 3:
            G3 += 1
        else:
            G4 += 1
    G0_freq.append(G0/N_j)
    G1_freq.append(G1/N_j)
    G2_freq.append(G2/N_j)
    G3_freq.append(G3/N_j)
    G4_freq.append(G4/N_j)
    

# functions to form gametes 
# note: uses derived probability distributions for gametes under strict bivalent formation
def maternal_gamete_formation(allele_count):
    global maternal_gamete
    if allele_count == 0:
        maternal_gamete = random.choices([['A', 'A'],['A','a'],['a', 'a']], [0, 0, 1], k=1)
    elif allele_count == 1:
        maternal_gamete = random.choices([['A', 'A'],['A','a'],['a', 'a']], [0, 1, 1], k=1)      
    elif allele_count == 2:
        maternal_gamete = random.choices([['A', 'A'],['A','a'],['a', 'a']], [1, 4, 1], k=1)         
    elif allele_count == 3:
        maternal_gamete = random.choices([['A', 'A'],['A','a'],['a', 'a']], [1, 1, 0], k=1)        
    elif allele_count == 4:
        maternal_gamete = random.choices([['A', 'A'],['A','a'],['a', 'a']], [1, 0, 0], k=1)

def paternal_gamete_formation(allele_count):
    global paternal_gamete
    if allele_count == 0:
        paternal_gamete = random.choices([['A', 'A'],['A','a'],['a', 'a']], [0, 0, 1], k=1)
    elif allele_count == 1:
        paternal_gamete = random.choices([['A', 'A'],['A','a'],['a', 'a']], [0, 1, 1], k=1)      
    elif allele_count == 2:
        paternal_gamete = random.choices([['A', 'A'],['A','a'],['a', 'a']], [1, 4, 1], k=1)         
    elif allele_count == 3:
        paternal_gamete = random.choices([['A', 'A'],['A','a'],['a', 'a']], [1, 1, 0], k=1)        
    elif allele_count == 4:
        paternal_gamete = random.choices([['A', 'A'],['A','a'],['a', 'a']], [1, 0, 0], k=1)
                

# creates the j+1th generation of individuals
# workflow: samples gametes from the parents, mutates gametes with mutation_bidirectional, and
# creates individuals from mutated gametes which are added to the j+1th generation
def j_plus_1th_generation(): 
    for i in range(N_j):    
        selfing_event = random.choices([True, False], [self_rate, 1-self_rate], k=1)

        if selfing_event[0] == True:
            # randomly samples one individual to be the mother; of form [['A','a','A','A']]
            mother = random.choices(gen_j_post_sel, k=1)

            # counts the number of 'A' alleles in mother
            allele_count_mother = mother[0].count('A')

            # selects gametes under derived proability distributions for strict bivalent formation
            maternal_gamete_formation(allele_count_mother)
            paternal_gamete_formation(allele_count_mother)

        elif selfing_event[0] == False:
            # randomly samples two parents (with replacement); parents are of form [['A','a','A','A'],['A','A','a','a']]
            parents = random.choices(gen_j_post_sel, k=2) 

            # counts the number of 'A' alleles in each parent
            allele_count_maternal = parents[0].count('A')
            allele_count_paternal = parents[1].count('A')

            # selects gametes from each parent under strict bivalent formation
            maternal_gamete_formation(allele_count_maternal)
            paternal_gamete_formation(allele_count_paternal)
        
        # append the maternal and paternal gametes to a list of all gametes BEFORE mutation for the jth generation
        gametes_j_pre_mut.append(maternal_gamete[0])  
        gametes_j_pre_mut.append(paternal_gamete[0])

        # runs mutation function to act on the ith set of gametes in the jth generation
        mutation_bidirectional(maternal_gamete[0], paternal_gamete[0])

        # append the maternal and paternal gametes to a list of all gametes AFTER mutation for the jth generation
        gametes_j_post_mut.append(maternal_gamete_mut)
        gametes_j_post_mut.append(paternal_gamete_mut)          
        
        # appends the alleles from the two lines above as the next individual in the jth generation post-mutation, pre-selection
        gen_j_pre_sel.append(maternal_gamete_mut + paternal_gamete_mut)


# purifying selection function 
# workflow: randomly chooses "True" or "False" for a selection event
# the respective frequencies for "True" are 0 (A, A), s*h1, s*h2, or s*h3 (heterozygotes), s (a, a)
# "True" represents that purifying selection occurs and (effectively) removes the individual
# "False" represents that the individual survives to reproduce and is kept in gen_j_post_sel
def selection(gen_j_pre_sel, gen_j_post_sel):
    for i in range(len(gen_j_pre_sel)):
        alleles = gen_j_pre_sel[i].count('A')
        # selection for recessive homozygotes - a2 a2 a2 a2 or G0
        if alleles == 0:
            sel_event_G0 = random.choices([True, False], [s, 1-s], k=1)
            if sel_event_G0[0] == False:
                gen_j_post_sel.append(gen_j_pre_sel[i])       
        # selection for "heterozygotes" - a1 a2 a2 a2 or G1
        elif alleles == 1:
            sel_event_G1 = random.choices([True, False], [h3*s, 1-h3*s], k=1)
            if sel_event_G1[0] == False:
                gen_j_post_sel.append(gen_j_pre_sel[i])
        # selection for "heterozygotes" - a1 a1 a2 a2 or G2
        elif alleles == 2:
            sel_event_G2 = random.choices([True, False], [h2*s, 1-h2*s], k=1)
            if sel_event_G2[0] == False:
                gen_j_post_sel.append(gen_j_pre_sel[i])
        # selection for "heterozygotes" - a1 a1 a1 a2 or G3
        elif alleles == 3:
            sel_event_G3 = random.choices([True, False], [h1*s, 1-h1*s], k=1)
            if sel_event_G3[0] == False:
                gen_j_post_sel.append(gen_j_pre_sel[i])
        # selection for dominant homozygotes - a1 a1 a1 a1 or G4
        else:
            gen_j_post_sel.append(gen_j_pre_sel[i])


### Define a few variables and lists to store information
gen0 = []

# initializes variables to hold the gamete frequencies for all generations 
# these lists are delineated across both pre and post mutation as well as g00, g01/g10, and g11
# pre mutation
g00_freq_pre_mut = []
g01_freq_pre_mut = [] 
g11_freq_pre_mut = []  
# post mutation
g00_freq_post_mut = []
g01_freq_post_mut = []  
g11_freq_post_mut = []  

# initializes variables to hold the genotype frequencies for all generations
G0_freq = []
G1_freq = []
G2_freq = []
G3_freq = []
G4_freq = []

# creates the initial generation of N individuals using random sampling with replacement
for i in range(N):  
    individual = random.choices(['A', 'a'], [p, q], k=4)
    gen0.append(individual)

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
    if j in range(g_bottleneck_start, g_bottleneck_start+g_bottleneck_length):  
        N_j = N_bottleneck
    else:
        N_j = N

    # runs selection function to act on the jth generation
    selection(gen_j_pre_sel, gen_j_post_sel)
    # empties pre-selection genetation to be filled using the j+1th generation immediately below
    gen_j_pre_sel = [] 

    # creates the next generation (j+1th generation) of individuals of size N_j
    j_plus_1th_generation()          

    # counts and stores gamete frequencies pre-mutation
    gamete_frequencies_pre_mut(gametes_j_pre_mut)

    # counts and stores gamete frequencies post-mutation
    gamete_frequencies_post_mut(gametes_j_post_mut) 
    
    # counts and stores genotype frequencies
    genotype_frequencies(gen_j_pre_sel)

### Plotting Functions


# plots both the gamete and genotype frequencies (each on its own axis)
# Sam's preferred plotting function/style
def plot_gamete_and_genotype_freq():
    
    x_index = range(g)  # creates a discrete time variable to plot against
    fig, axs = plt.subplots(2, 1, sharex = 'col', sharey = 'row')  

    fig.set_size_inches(14, 8) 
    axs[0].plot(x_index, g00_freq_pre_mut, color = '#008aff')
    axs[0].plot(x_index, g01_freq_pre_mut, color = '#00a50e')
    axs[0].plot(x_index, g11_freq_pre_mut, color = '#da1b00')
    axs[0].legend(['g00', 'g01', 'g11'])
    axs[0].set_ylabel('Gamete Frequencies')

    axs[1].plot(x_index, G0_freq, color = '#008aff')
    axs[1].plot(x_index, G1_freq, color = '#eb7f00')
    axs[1].plot(x_index, G2_freq, color = '#00a50e')
    axs[1].plot(x_index, G3_freq, color = '#8800da')
    axs[1].plot(x_index, G4_freq, color = '#da1b00')
    axs[1].legend(['G0', 'G1', 'G2', 'G3', 'G4'])
    axs[1].set_ylabel('Genotype Frequencies')

    # labels
    fig.supxlabel('Time (in generations)')
    fig.suptitle('Gamete and Genotype Frequency Over Time')

    # creates a light gray shaded region for the population shift or bottleneck
    for ax in axs:
        shaded_region = Rectangle((g_bottleneck_start, 0), g_bottleneck_length, 1)
        shaded_region.set_color('0.85')
        ax.add_patch(shaded_region)
        ax.set_xlim(0, g-1)
        ax.set_ylim(-.1, 1.1)


# plots the three gamete frequencies (g00, g01/g10, and g11) 
# before mutation on three separate axes
def plot_gamete_freq():
    
    x_index = range(g)  # creates a discrete time variable to plot against
    fig, axs = plt.subplots(1, 3, sharex = 'col', sharey = 'row')  

    fig.set_size_inches(16, 6) 
    axs[0].plot(x_index, g00_freq_pre_mut)
    axs[0].set_xlabel('g00')
    axs[1].plot(x_index, g01_freq_pre_mut)
    axs[1].set_xlabel('g01')
    axs[2].plot(x_index, g11_freq_pre_mut)
    axs[2].set_xlabel('g11')

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
    

# plots the five genotype frequencies (G0, G1, G2, G3, and G4) on 5 axes in a 3 by 2 grid of 6 total axes
def plot_genotype_freq():

    x_index = range(g)  # creates a discrete time variable to plot against
    fig, axs = plt.subplots(2, 3, sharex = 'col', sharey = 'row')  

    fig.set_size_inches(16, 6)
    axs[0, 0].plot(x_index, G0_freq)
    axs[0, 0].set_xlabel('G0')
    axs[0, 1].plot(x_index, G1_freq)
    axs[0, 1].set_xlabel('G1')
    axs[0, 2].plot(x_index, G2_freq)
    axs[0, 2].set_xlabel('G2')
    axs[1, 0].plot(x_index, G3_freq)
    axs[1, 0].set_xlabel('G3')
    axs[1, 1].plot(x_index, G4_freq)
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


# plots the three gamete frequencies (g00, g01/g10, and g11) 
# before mutation on the same axis
def plot_gamete_freq_single_axis():
    
    x_index = range(g)  # creates a discrete time variable to plot against
    fig, ax = plt.subplots(figsize=(16, 6))
    plt.axis([0, g-1, -.1, 1.1])  
 
    ax.plot(x_index, g00_freq_pre_mut, color = '#008aff')
    ax.plot(x_index, g01_freq_pre_mut, color = '#00a50e' )
    ax.plot(x_index, g11_freq_pre_mut, color = '#da1b00')
    ax.legend(['g00', 'g01', 'g11'])
    
    # labels
    fig.supxlabel('Time (in generations)')
    fig.suptitle('Gamete Frequency Over Time')
    fig.supylabel('Gamete frequency (pre-mut)')

    # creates a light gray shaded region for the population shift or bottleneck
    shaded_region = Rectangle((g_bottleneck_start, 0), g_bottleneck_length, 1)
    shaded_region.set_color('0.85')
    ax.add_patch(shaded_region)


# plots the five genotype frequencies (G0, G1, G2, G3, and G4) on a single axis
def plot_genotype_freq_single_axis():    
    
    x_index = range(g)  # creates a discrete time variable to plot against
    fig, ax = plt.subplots(figsize=(16, 6))
    plt.axis([0, g-1, -.1, 1.1])  
 
    ax.plot(x_index, G0_freq, color = '#008aff')
    ax.plot(x_index, G1_freq, color = '#00a50e')
    ax.plot(x_index, G2_freq, color = '#da1b00')
    ax.plot(x_index, G3_freq, color = '#8800da')
    ax.plot(x_index, G4_freq, color = '#eb7f00')
    ax.legend(['G0', 'G1', 'G2', 'G3', 'G4'])
    
    # labels
    fig.supxlabel('Time (in generations)')
    fig.suptitle('Genotype Frequency Over Time')
    fig.supylabel('Genotype frequency')

    # creates a light gray shaded region for the population shift or bottleneck
    shaded_region = Rectangle((g_bottleneck_start, 0), g_bottleneck_length, 1)
    shaded_region.set_color('0.85')
    ax.add_patch(shaded_region)


# plots the gamete differential between before and after mutation on a single axis
# note: only useful for simulations with a 200 or fewer generations
def plot_gamete_diff_single_axis():
    global g00_freq_pre_mut
    global g01_freq_pre_mut
    global g11_freq_pre_mut

    global g00_freq_post_mut
    global g01_freq_post_mut
    global g11_freq_post_mut

    g00_freq_pre_mut = np.array(g00_freq_pre_mut)
    g00_freq_post_mut = np.array(g00_freq_post_mut)
    g00_freq_diff = g00_freq_post_mut - g00_freq_pre_mut

    g01_freq_pre_mut = np.array(g01_freq_pre_mut)
    g01_freq_post_mut = np.array(g01_freq_post_mut)
    g01_freq_diff = g01_freq_post_mut - g01_freq_pre_mut

    g11_freq_pre_mut = np.array(g11_freq_pre_mut)
    g11_freq_post_mut = np.array(g11_freq_post_mut)
    g11_freq_diff = g11_freq_post_mut - g11_freq_pre_mut

    x_index = range(g)  # creates a discrete time variable to plot against
    fig, ax = plt.subplots(figsize=(16, 6))  

    # labels
    fig.supxlabel('Time (in generations)')
    fig.suptitle('Gamete Differentials Over Time')
    fig.supylabel('Gamete Differential (post-mut - pre-mut)')

    ax.plot(x_index, g00_freq_diff, color = '#008aff')
    ax.plot(x_index, g01_freq_diff, color = '#00a50e')
    ax.plot(x_index, g11_freq_diff, color = '#da1b00')
    ax.legend(['g00', 'g01', 'g11'])

    # creates a light gray shaded region for the population shift or bottleneck
    shaded_region = Rectangle((g_bottleneck_start, -.25), g_bottleneck_length, .5)
    shaded_region.set_color('0.85')
        

plot_gamete_and_genotype_freq()

plt.show()
