import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as optimize
from sympy import *

def root_solve_init():
    g0, g1, s, h, mu, nu = symbols('g0 g1 s h mu nu')

    g1 = 1 - g0

    # equations to parameterize relative fitnesses
    w_bar = 1 - s*(2*h*g0*g1 + g1**2) # average fitness
    w0 = 1/w_bar # fitness of AA/ancestral homozygote
    w1 = (1-h*s)/w_bar # fitness of heterozygote
    w2 = (1-s)/w_bar # fitness of aa/derived homozygote

    # equations to model selection and meiosis
    sel_g0_eqn = w0*g0**2 + w1*g0*g1 
    sel_g1_eqn = w1*g0*g1 + w2*g1**2

    # equations to model mutation
    mut_g0_eqn = 1e7*(sel_g0_eqn*(1-mu) + sel_g1_eqn*nu - g0)
    mut_g1_eqn = sel_g0_eqn*mu + sel_g1_eqn*(1-nu) - g1

    g0_eqn_numeric = lambdify([(s, h, mu, nu)], mut_g0_eqn, 'scipy')

    return g0_eqn_numeric

def bifn_solve_init():
    g0, g1, s, h, mu, nu = symbols('g0 g1 s h mu nu')

    g1 = 1 - g0

    # equations to parameterize relative fitnesses
    w_bar = 1 - s*(2*h*g0*g1 + g1**2) # average fitness
    w0 = 1/w_bar # fitness of AA/ancestral homozygote
    w1 = (1-h*s)/w_bar # fitness of heterozygote
    w2 = (1-s)/w_bar # fitness of aa/derived homozygote

    # equations to model selection and meiosis
    sel_g0_eqn = w0*g0**2 + w1*g0*g1 
    sel_g1_eqn = w1*g0*g1 + w2*g1**2

    # equations to model mutation
    mut_g0_eqn = 1e7*(sel_g0_eqn*(1-mu) + sel_g1_eqn*nu - g0)
    mut_g1_eqn = sel_g0_eqn*mu + sel_g1_eqn*(1-nu) - g1

    g0_deriv_sym = diff(mut_g0_eqn, g0)

    g0_eqn = lambdify([(h, mu, nu)], mut_g0_eqn, 'scipy')
    g0_deriv = lambdify([(h, mu, nu)], g0_deriv_sym, 'scipy')

    return g0_eqn, g0_deriv

def cusp_solve_init():
    g0, g1, s, h, mu, nu = symbols('g0 g1 s h mu nu')

    g1 = 1 - g0

    # equations to parameterize relative fitnesses
    w_bar = 1 - s*(2*h*g0*g1 + g1**2) # average fitness
    w0 = 1/w_bar # fitness of AA/ancestral homozygote
    w1 = (1-h*s)/w_bar # fitness of heterozygote
    w2 = (1-s)/w_bar # fitness of aa/derived homozygote

    # equations to model selection and meiosis
    sel_g0_eqn = w0*g0**2 + w1*g0*g1 
    sel_g1_eqn = w1*g0*g1 + w2*g1**2

    # equations to model mutation
    mut_g0_eqn = 1e7*(sel_g0_eqn*(1-mu) + sel_g1_eqn*nu - g0)
    mut_g1_eqn = sel_g0_eqn*mu + sel_g1_eqn*(1-nu) - g1

    g0_deriv_sym = diff(mut_g0_eqn, g0)
    g0_sec_deriv_sym = diff(g0_deriv_sym, g0)

    g0_eqn = lambdify([(mu, nu)], mut_g0_eqn, 'scipy')
    g0_deriv = lambdify([(mu, nu)], g0_deriv_sym, 'scipy')
    g0_sec_deriv = lambdify([(mu, nu)], g0_sec_deriv_sym, 'scipy')

    return g0_eqn, g0_deriv, g0_sec_deriv


g0, s, h = symbols('g0 s h')

h_val = 1
mu_val = 5e-8
nu_val = 1e-9
s_val = 6.5e-7

#cusp equations
cusp_eqn_params_1, cusp_eqn_params_2, cusp_eqn_params_3 = cusp_solve_init()

cusp_eqn_1_eval = cusp_eqn_params_1([mu_val, nu_val])
cusp_eqn_2_eval = cusp_eqn_params_2([mu_val, nu_val])
cusp_eqn_3_eval = cusp_eqn_params_3([mu_val, nu_val])

cusp_set_func = lambdify([(g0, s, h)], [cusp_eqn_1_eval, cusp_eqn_2_eval, cusp_eqn_3_eval], 'scipy')

cusp_soln = optimize.root(cusp_set_func, [.5, 0, .75], method='hybr')

print(cusp_soln.x)

#bifn equations
bifn_eqn_params_1, bifn_eqn_params_2 = bifn_solve_init()

bifn_eqn_1_eval = bifn_eqn_params_1([h_val, mu_val, nu_val])
bifn_eqn_2_eval = bifn_eqn_params_2([h_val, mu_val, nu_val])

bifn_set_func = lambdify([(g0, s)], [bifn_eqn_1_eval, bifn_eqn_2_eval], 'scipy')

bifn_soln_1 = optimize.fsolve(bifn_set_func, [.5, .000001])
bifn_soln_2 = optimize.fsolve(bifn_set_func, [0, .00001])


#root/equilibria equations
g0_eqn_params = root_solve_init()

g0_eqn_eval = g0_eqn_params([s_val, h_val, mu_val, nu_val])

g0_eqn_func = lambdify(g0, g0_eqn_eval, 'scipy')

root_soln_1 = optimize.root_scalar(g0_eqn_func, bracket = [0, bifn_soln_2[0]], method ='toms748')
root_soln_2 = optimize.root_scalar(g0_eqn_func, bracket = [bifn_soln_2[0], bifn_soln_1[0]], method ='toms748')
root_soln_3 = optimize.root_scalar(g0_eqn_func, bracket = [bifn_soln_1[0], 1], method = 'toms748')

print(root_soln_1.root)
print(root_soln_2.root)
print(root_soln_3.root)

# iterations = 1000

# s_val_range = np.logspace(-6, -5, iterations)
# h_val = 1
# mu_val = 5e-8
# nu_val = 1e-9

# cusp_eqn_set = cusp_solve_init(mu_val, nu_val)
# cusp_soln = optimize.fsolve(cusp_eqn_set, [.5, .000001, .9])
# print(cusp_soln)

# bifn_eqn_set = bifn_solve_init(h_val, mu_val, nu_val)
# bifn_soln = optimize.fsolve(bifn_eqn_set, [.5, .000001])
# print(bifn_soln)
# bifn_soln_2 = optimize.fsolve(bifn_eqn_set, [0, .00001])
# print(bifn_soln_2)

# for i in range(iterations):
#     g0_eqn_numeric = root_eval_init(s_val_range[i], h_val, mu_val, nu_val)
#     soln = optimize.root_scalar(g0_eqn_numeric, bracket = [.5, 1], method='brentq')


