import random

# Let's define some parameters
p = .3  # allele frequency for A (derived allele)
q = 1-p  # allele frequency for a (ancestral allele)
g0_init_value = .3
g1_init_value = 0
g2_init_value = 1 - g0_init_value - g1_init_value
N = 10  # population size

gen0 = []

for i in range(N):
    individual = random.choices([['A', 'A'], ['A', 'a'], ['a', 'a']], [g0_init_value, g1_init_value, g2_init_value], k=2)
    gen0.append(individual[0] + individual[1])

print(gen0)
