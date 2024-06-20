import random
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import copy

# Let's define some parameters
p = .5  # allele frequency for a1
q = 1-p  # allele frequency for a2
N = 100  # population size
g = 100  # number of generations
g_bottleneck_start = 400  # generation at which bottleneck starts
g_bottleneck_length = 200  # number of generations for which the bottleneck lasts
N_bottleneck = 500  # population size during the bottleneck
s = .00001  # selection coefficient
h1 = .25  # dominance coefficient for G3
h2 = .5  # dominance coefficient for G2
h3 = .75  # dominance coefficient for G1
mu = .000001  # mutation rate (from 'a1' to 'a2')
nu = .000001  # mutation rate (from 'a2' to 'a1')
self_rate = .8  # selfing rate (what parameter/letter/variable to use for this?)

# Let's create some functions

# creates mutation in already selected gametes of the parents occuring at rates mu (a1 to a2) and nu (a2 to a1) 
# parents is a list of two individuals ex. [['a1', 'a2', 'a2', 'a2'], ['a1', 'a1', 'a2', 'a1']])
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
        if (mut_event_1[0] == True) and (maternal_gamete[l] == 'a1'):  
            maternal_gamete_mut[l] = 'a2'
        elif (mut_event_2[0] == True) and (maternal_gamete[l] == 'a2'):
            maternal_gamete_mut[l] = 'a1'
    
    for l in range(len(paternal_gamete)):  # same as above, but for the paternal gamete
        mut_event_1 = random.choices([True, False], [mu, 1-mu], k=1) 
        mut_event_2 = random.choices([True, False], [nu, 1-nu], k=1)   
        if (mut_event_1[0] == True) and (paternal_gamete[l] == 'a1'):
            paternal_gamete_mut[l] = 'a2'
        elif (mut_event_2[0] == True) and (paternal_gamete[l] == 'a2'):
            paternal_gamete_mut[l] = 'a1'


# creates a function to count gamete frequencies for a single generation of the population (before mutation)
# then after counting, makes these gamete frequencies relative by dividing by 2N_j
# finally, appends these gamete frequencies to a list of gamete frequencies (before mutation) for all generations
# this is done once per generation and not for each sampling of maternal/paternal gametes
def gamete_frequencies_pre_mut(gametes_j_pre_mut):
    
    g00 = gametes_j_pre_mut.count(['a2', 'a2'])/(2*N_j)
    g11 = gametes_j_pre_mut.count(['a1', 'a1'])/(2*N_j)
    g01 = 1 - g00 - g11

    g00_freq_pre_mut.append(g00)
    g01_freq_pre_mut.append(g01)
    g11_freq_pre_mut.append(g11)


# identical function to gamete_frequencies_pre_mut, but appends the gamete frequencies to a list after mutation
def gamete_frequencies_post_mut(gametes_j_post_mut):
    
    g00 = gametes_j_post_mut.count(['a2', 'a2'])/(2*N_j)
    g11 = gametes_j_post_mut.count(['a1', 'a1'])/(2*N_j)
    g01 = 1 - g00 - g11

    g00_freq_post_mut.append(g00)
    g01_freq_post_mut.append(g01)
    g11_freq_post_mut.append(g11)


# similar to gamete_frequencies, but counts genotype frequencies which are all appended to separate lists
def genotype_frequencies(gen_j_pre_sel):
    G0 = 0
    G1 = 0
    G2 = 0
    G3 = 0
    G4 = 0
    for i in range(len(gen_j_pre_sel)):
        allele_count = gen_j_pre_sel[i].count('a1')
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
    

# creates a function to create the j+1th generation of individuals
# does so by sampling gametes from the parents, running the mutation_bidirectional function on these gametes,
# creating individuals for the next generation from the mutated gametes, 
# and appending these individuals to a list of all individuals for the j+1th generation
def j_plus_1th_generation(): 
    for i in range(N_j):    
        selfing_event = random.choices([True, False], [self_rate, 1-self_rate], k=1)

        if selfing_event[0] == True:
            mother = random.choices(gen_j_post_sel, k=1)  # randomly samples one individual from the prior (jth) generation after selection (without replacement)
            maternal_gamete = random.sample(mother[0], 2)
            paternal_gamete = random.sample(mother[0], 2)        
        elif selfing_event[0] == False:
            parents = random.choices(gen_j_post_sel, k=2) 
            maternal_gamete = random.sample(parents[0], 2)
            paternal_gamete = random.sample(parents[1], 2)
        
        # append the maternal and paternal gametes to a list of all gametes BEFORE mutation for the jth generation
        gametes_j_pre_mut.append(maternal_gamete)  
        gametes_j_pre_mut.append(paternal_gamete)

        # runs mutation function to act on the ith set of gametes in the jth generation
        mutation_bidirectional(maternal_gamete, paternal_gamete)

        # append the maternal and paternal gametes to a list of all gametes AFTER mutation for the jth generation
        gametes_j_post_mut.append(maternal_gamete_mut)
        gametes_j_post_mut.append(paternal_gamete_mut)          
        
        # appends the alleles from the two lines above as the next individual in the jth generation post-mutation, pre-selection
        gen_j_pre_sel.append(maternal_gamete_mut + paternal_gamete_mut)


# creates selection on the jth generation using a random choice of a selection event being "True" or "False"
# the respective frequencies for "True" are 0 (a1, a1), s*h (heterozygote), s (a2, a2)
# "True" represents that a selection event does occur and the individual is effectively removed from the population
# "False" represents that a selection event did not occur and the individual is kept in the gen_post_sel
def selection(gen_j_pre_sel, gen_j_post_sel):
    for i in range(len(gen_j_pre_sel)):
        alleles = gen_j_pre_sel[i].count('a1')
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

for i in range(N):  # creates the initial generation of N individuals using random sampling with replacement
    individual = random.choices(['a1', 'a2'], [p, q], k=4)
    gen0.append(individual)

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
    if j in range(g_bottleneck_start, g_bottleneck_start+g_bottleneck_length):  
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


