import random
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import copy

# Let's define some parameters
pa = .5  # allele frequency for A
qa = 1-pa  # allele frequency for a
pb = .5  # allele frequency for B
qb = 1-pb  # allele frequency for b
N = 10000  # population size
g = 1000  # number of generations
g_bottleneck_start = 40000  # generation at which bottleneck starts
g_bottleneck_length = 200  # number of generations for which the bottleneck lasts
N_bottleneck = 100  # population size during the bottleneck
s = 0  # selection coefficient
h1 = .25  # dominance coefficient for G3
h2 = .5  # dominance coefficient for G2
h3 = .75  # dominance coefficient for G1
mu = 0  # mutation rate (from 'A' to 'a' and 'B' to 'b')
nu = 0  # mutation rate (from 'a' to 'A' and 'b' to 'B')
self_rate = 1  # selfing rate (what parameter/letter/variable to use for this?)

# Let's create some functions

# creates mutation in already selected gametes of the parents occuring at rates mu (A to a) and nu (a to A) 
# parents is a list of two individuals ex. [['A', 'a', 'b', 'b'], ['A', 'A', 'b', 'B']])
# gametes are defined both pre and post-mutation and are a random sample of any two alleles from the parents (disregarding sub-genomes)
def mutation_bidirectional(maternal_gamete, paternal_gamete):

    global maternal_gamete_mut
    maternal_gamete_mut = list(copy.copy(maternal_gamete))
    global paternal_gamete_mut
    paternal_gamete_mut = list(copy.copy(paternal_gamete))

    for l in range(len(maternal_gamete)):  # for each allele in the maternal gamete, will run mutation at random
        mut_event_1 = random.choices([True, False], [mu, 1-mu], k=1)  # randomly sets mutation to occur (i.e. True) with likelihood mu
        mut_event_2 = random.choices([True, False], [nu, 1-nu], k=1)  # randomly sets mutation to occur (i.e. True) with likelihood nu
        # if the mutation event is "True" and the allele at the lth position in the gamete is 'a1', changes 'a1' to 'a2'
        if (mut_event_1[0] == True) and (maternal_gamete[l] == 'A'):  
            maternal_gamete_mut[l] = 'a'
        elif (mut_event_2[0] == True) and (maternal_gamete[l] == 'a'):
            maternal_gamete_mut[l] = 'A'
        elif (mut_event_1[0] == True) and (maternal_gamete[l] == 'B'):
            maternal_gamete_mut[l] = 'b'
        elif (mut_event_2[0] == True) and (maternal_gamete[l] == 'b'):
            maternal_gamete_mut[l] = 'B'
    
    for l in range(len(paternal_gamete)):  # same as above, but for the paternal gamete
        mut_event_1 = random.choices([True, False], [mu, 1-mu], k=1) 
        mut_event_2 = random.choices([True, False], [nu, 1-nu], k=1)   
        if (mut_event_1[0] == True) and (paternal_gamete[l] == 'A'):
            paternal_gamete_mut[l] = 'a'
        elif (mut_event_2[0] == True) and (paternal_gamete[l] == 'a'):
            paternal_gamete_mut[l] = 'A'
        elif (mut_event_1[0] == True) and (paternal_gamete[l] == 'B'):
            paternal_gamete_mut[l] = 'b'
        elif (mut_event_2[0] == True) and (paternal_gamete[l] == 'b'):
            paternal_gamete_mut[l] = 'B'


# creates a function to count gamete frequencies for a single generation of the population (before mutation)
# then after counting, makes these gamete frequencies relative by dividing by 2N_j
# finally, appends these gamete frequencies to a list of gamete frequencies (before mutation) for all generations
# this is done once per generation and not for each sampling of maternal/paternal gametes
def gamete_frequencies_pre_mut(gametes_j_pre_mut):
    
    g00 = gametes_j_pre_mut.count(['a', 'b'])/(2*N_j)
    g10 = gametes_j_pre_mut.count(['A', 'b'])/(2*N_j)
    g01 = gametes_j_pre_mut.count(['a', 'B'])/(2*N_j)
    g11 = gametes_j_pre_mut.count(['A', 'B'])/(2*N_j)

    g00_freq_pre_mut.append(g00)
    g10_freq_pre_mut.append(g10)
    g01_freq_pre_mut.append(g01)
    g11_freq_pre_mut.append(g11)


# identical function to gamete_frequencies_pre_mut, but appends the gamete frequencies to a list after mutation
def gamete_frequencies_post_mut(gametes_j_post_mut):
    
    g00 = gametes_j_post_mut.count(['a', 'b'])/(2*N_j)
    g10 = gametes_j_post_mut.count(['A', 'b'])/(2*N_j)
    g01 = gametes_j_post_mut.count(['a', 'B'])/(2*N_j)
    g11 = gametes_j_post_mut.count(['A', 'B'])/(2*N_j)

    g00_freq_post_mut.append(g00)
    g10_freq_post_mut.append(g10)
    g01_freq_post_mut.append(g01)
    g11_freq_post_mut.append(g11)


# similar to gamete_frequencies, but counts genotype frequencies which are all appended to separate lists
def genotype_frequencies(gen_j_pre_sel):
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
        allele_count_sga = gen_j_pre_sel[i].count('A')
        allele_count_sgb = gen_j_pre_sel[i].count('B')
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
    
# creates a function to create the j+1th generation of individuals
# does so by sampling gametes from the parents, running the mutation_bidirectional function on these gametes,
# creating individuals for the next generation from the mutated gametes, 
# and appending these individuals to a list of all individuals for the j+1th generation
def j_plus_1th_generation(): 
    for i in range(N_j):    

        selfing_event = random.choices([True, False], [self_rate, 1-self_rate], k=1)

        if selfing_event[0] == True:
            mother = random.sample(gen_j_post_sel, 1)  # randomly samples one individual to be the mother
            
            mother_sga = mother[0][:2]
            mother_sgb = mother[0][2:]
            maternal_gamete = random.sample(mother_sga, 1) + random.sample(mother_sgb, 1)  # creates a maternal gamete under strict preferential pairing
            paternal_gamete =  random.sample(mother_sga, 1) + random.sample(mother_sgb, 1)  # creates a "paternal" (second maternal) gamete under strict preferential pairing

        elif selfing_event[0] == False: 
            parents = random.choices(gen_j_post_sel, k=2)
            maternal_gamete = random.sample(parents[0][:2], 1) + random.sample(parents[0][2:], 1)
            paternal_gamete = random.sample(parents[1][:2], 1) + random.sample(parents[1][2:], 1)
        
        # append the maternal and paternal gametes to a list of all gametes BEFORE mutation for the jth generation
        gametes_j_pre_mut.append(maternal_gamete)  
        gametes_j_pre_mut.append(paternal_gamete)

        # runs mutation function to act on the ith set of gametes in the jth generation
        mutation_bidirectional(maternal_gamete, paternal_gamete)

        # append the maternal and paternal gametes to a list of all gametes AFTER mutation for the jth generation
        gametes_j_post_mut.append(maternal_gamete_mut)
        gametes_j_post_mut.append(paternal_gamete_mut)          
        
        # appends the alleles from the two lines above to be the ith individual in the jth generation
        individual = [maternal_gamete_mut[0], paternal_gamete_mut[0], maternal_gamete_mut[1], paternal_gamete_mut[1]]
        
        gen_j_pre_sel.append(individual)

# creates selection on the jth generation using a random choice of a selection event being "True" or "False"
# the respective frequencies for "True" are 0 (a1, a1), s*h (heterozygote), s (a2, a2)
# "True" represents that a selection event does occur and the individual is effectively removed from the population
# "False" represents that a selection event did not occur and the individual is kept in the gen_post_sel
def selection(gen_j_pre_sel, gen_j_post_sel):
    for i in range(len(gen_j_pre_sel)):
        alleles_sga = gen_j_pre_sel[i].count('A')
        alleles_sgb = gen_j_pre_sel[i].count('B')
        alleles = alleles_sga + alleles_sgb
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


for i in range(N):  # creates the initial generation of N individuals using random sampling with replacement
    subgenome_a = random.choices(['A', 'a'], [pa, qa], k=2)
    subgenome_b = random.choices(['B', 'b'], [pb, qb], k=2)
    gen0.append(subgenome_a + subgenome_b)

gen_j_pre_sel = gen0  # initializes a list with the ancestral (0th) generation to include generations each with == N individuals before selection

for j in range(g):  # performs the following sequence of steps for g generations
    
    # additional variables and lists

    # use gen_j_pre/post_sel to 1) facilitate the selection function and 2) generate genotype frequencies over generations
    gen_j_post_sel = []  # initializes a list for the jth generation of individuals post selection (i.e. who survived to reproduce)

    # use gametes_j_pre/post_mut to 1) generate gamete frequencies over generations both before and after mutation and 
    # 2) measure the impact of mutation on gamete frequencies via comparison
    gametes_j_pre_mut = []  # initializes a list to hold the gametes for the jth generation before mutation
    gametes_j_post_mut = []  # initializes a list to hold the gametes for the jth generation after mutation 

    # sets N_j to either be the population size (N) or the population size for the bottleneck (N_bottleneck)
    if j in range(g_bottleneck_start, g_bottleneck_start + g_bottleneck_length):  
        N_j = N_bottleneck
    else:
        N_j = N

    # runs selection function to act on the jth generation
    selection(gen_j_pre_sel, gen_j_post_sel)
    gen_j_pre_sel = []  # empties the jth generation pre-selection list before filling it with the j+1th generation (see function)

    # creates the next generation (j+1th generation) of individuals of size N_j
    j_plus_1th_generation()          

    gamete_frequencies_pre_mut(gametes_j_pre_mut)

    gamete_frequencies_post_mut(gametes_j_post_mut) 
    
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
    axs[0].set_ylabel('Gamete Frequencies')

    axs[1].plot(x_index, G00_freq, color = '#008aff')
    axs[1].plot(x_index, G1_freq, color = '#8800da')
    axs[1].plot(x_index, G2_freq, color = '#eb7f00')
    axs[1].plot(x_index, G3_freq, color = '#00a50e')
    axs[1].plot(x_index, G22_freq, color = '#da1b00')
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


# plots the genotype frequencies delineated across subgenomes
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




