from sympy import symbols, Eq, Rational, solve, factor, latex
from sympy import Rational

"""
Genotype and gamete equilibria for segmental allotetraploids with 1 HE. 

Outputs:
    Gamete equilibria (g00, g01, g10, g11)
    Delineated genotype equilibria (G00, G01, G02, G10, G11, G12, G20, G21, G22)
    Combined genotype equilibria (G0, G1, G2, G3, G4)
"""

# establish symbolic objects for computation
g00, g01, g10, g11, p, q = symbols('g00 g01 g10 g11 p q')

# establish symbolic objects for output
G00, G01, G10, G02, G11, G20, G12, G21, G22 = symbols('G00 G01 G10 G02 G11 G20 G12 G21 G22')
G0, G1, G2, G3, G4 = symbols('G0 G1 G2 G3 G4')

# determines which variable ('p' or 'q'), the equilibria are expressed in 
in_terms_of = 'q'

# gametic expressions derived from combined probabilities of gamete production after random mating of the previous generation
eqn1 = Eq(g00**2+Rational(17, 16)*g00*g01+Rational(17, 16)*g00*g10+Rational(3, 16)*g01**2+Rational(3, 16)*g10**2+Rational(7, 16)*g01*g10+Rational(7, 16)*g00*g11+Rational(1, 16)*g01*g11+Rational(1, 16)*g10*g11, g00)
eqn2 = Eq(Rational(3, 16)*g00*g01+Rational(11, 16)*g00*g10+Rational(1, 16)*g01**2+Rational(9, 16)*g10**2+Rational(9, 16)*g01*g10+Rational(9, 16)*g00*g11+Rational(3, 16)*g01*g11+Rational(11, 16)*g10*g11, g10)
eqn3 = Eq(Rational(11, 16)*g00*g01+Rational(3, 16)*g00*g10+Rational(9, 16)*g01**2+Rational(1, 16)*g10**2+Rational(9, 16)*g01*g10+Rational(9, 16)*g00*g11+Rational(11, 16)*g01*g11+Rational(3, 16)*g10*g11, g01)
eqn4 = Eq(Rational(1, 16)*g00*g01+Rational(1, 16)*g00*g10+Rational(3, 16)*g01**2+Rational(3, 16)*g10**2+Rational(7, 16)*g01*g10+Rational(7, 16)*g00*g11+Rational(17, 16)*g01*g11+Rational(17, 16)*g10*g11+g11**2, g11)

# simple identity that all gamete frequencies sum to 1
eqn5 = Eq(g00+g01+g10+g11, 1)

# determines if expressions for p or q will be used in solving the gamete equilibria
# note: here, g00+g01=g00+g10=q because HEs occur across subgenomes making qa and qb indistinguishable
if in_terms_of == 'q':
    eqn6 = Eq(g00+g01, q)   
    eqn7 = Eq(g00+g10, q)
elif in_terms_of == 'p':
    eqn6 = Eq(g10+g11, p)
    eqn7 = Eq(g01+g11, p)
else:
    raise Exception("Expressions can only be expressed in terms of p or q.")

# establishes the equations to be used in solving for equilibria and the variables to be solved for
eqn_system = [eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7]
variable_list = [g00, g01, g10, g11]

# solving function itself
gamete_soln_set = solve(eqn_system, variable_list, dict=True)

# sets the equilibrium gamete solutions equal to their respective variables
g00_eq = gamete_soln_set[0][g00]
g01_eq = gamete_soln_set[0][g01]
g10_eq = gamete_soln_set[0][g10]
g11_eq = gamete_soln_set[0][g11]

# creates a dictionary of factored expressions for the delineated genotypes (i.e. G00, G01, G02, etc.)
delineated_genotype_soln_set = {
    G00: factor(g00_eq**2),
    G01: factor(2*g00_eq*g01_eq),
    G10: factor(2*g00_eq*g10_eq),
    G02: factor(g01_eq**2),
    G11: factor(2*(g00_eq*g11_eq + g01_eq*g10_eq)),
    G20: factor(g10_eq**2),
    G12: factor(2*g01_eq*g11_eq),
    G21: factor(2*g10_eq*g11_eq),
    G22: factor(g11_eq**2),
}

# creates a dictionary of factored expressions for the combined genotypes (i.e. G0, G1, G2, G3, and G4)
combined_genotype_soln_set = {
    G0: factor(g00_eq**2),
    G1: factor(2*g00_eq*g01_eq + 2*g00_eq*g10_eq),
    G2: factor(g01_eq**2 + 2*(g00_eq*g11_eq + g01_eq*g10_eq) + g10_eq**2),
    G3: factor(2*g01_eq*g11_eq + 2*g10_eq*g11_eq),
    G4: factor(g11_eq**2),
}
