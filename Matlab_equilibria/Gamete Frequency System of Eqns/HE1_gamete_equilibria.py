from sympy import symbols, Eq, Rational, nonlinsolve, solve
from sympy import Rational

g00, g01, g10, g11, p, q = symbols('g00 g01 g10 g11 p q')

in_terms_of = 'q'

eqn1 = Eq(g00**2+Rational(17, 16)*g00*g01+Rational(17, 16)*g00*g10+Rational(3, 16)*g01**2+Rational(3, 16)*g10**2+Rational(7, 16)*g01*g10+Rational(7, 16)*g00*g11+Rational(1, 16)*g01*g11+Rational(1, 16)*g10*g11, g00)
eqn2 = Eq(Rational(3, 16)*g00*g01+Rational(11, 16)*g00*g10+Rational(1, 16)*g01**2+Rational(9, 16)*g10**2+Rational(9, 16)*g01*g10+Rational(9, 16)*g00*g11+Rational(3, 16)*g01*g11+Rational(11, 16)*g10*g11, g10)
eqn3 = Eq(Rational(11, 16)*g00*g01+Rational(3, 16)*g00*g10+Rational(9, 16)*g01**2+Rational(1, 16)*g10**2+Rational(9, 16)*g01*g10+Rational(9, 16)*g00*g11+Rational(11, 16)*g01*g11+Rational(3, 16)*g10*g11, g01)
eqn4 = Eq(Rational(1, 16)*g00*g01+Rational(1, 16)*g00*g10+Rational(3, 16)*g01**2+Rational(3, 16)*g10**2+Rational(7, 16)*g01*g10+Rational(7, 16)*g00*g11+Rational(17, 16)*g01*g11+Rational(17, 16)*g10*g11+g11**2, g11)

eqn5 = Eq(g00+g01+g10+g11, 1)

if in_terms_of == 'q':
    eqn6 = Eq(g00+g01, q)   
    eqn7 = Eq(g00+g10, q)
elif in_terms_of == 'p':
    eqn6 = Eq(g10+g11, p)
    eqn7 = Eq(g01+g11, p)
else:
    raise Exception("Expressions can only be expressed in terms of p or q.")

eqn_system = [eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7]
variable_list = [g00, g01, g10, g11]

soln_set = nonlinsolve(eqn_system, variable_list)

print(soln_set)

soln_set2 = solve(eqn_system, variable_list, set=True)

print(soln_set2)