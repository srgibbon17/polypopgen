import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as optimize
from sympy import *

def root_eval_init(s_val, h_val, mu_val, nu_val):
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

    g0_eqn_sym = mut_g0_eqn.subs(s, s_val)
    g0_eqn_sym = g0_eqn_sym.subs(h, h_val)
    g0_eqn_sym = g0_eqn_sym.subs(mu, mu_val)
    g0_eqn_sym = g0_eqn_sym.subs(nu, nu_val)

    g0_eqn_numeric = lambdify(g0, g0_eqn_sym, 'scipy')

    return g0_eqn_numeric

def bifn_solve_init(h_val, mu_val, nu_val):
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

    g0_deriv = diff(mut_g0_eqn, g0)

    g0_eqn = mut_g0_eqn.subs(h, h_val)
    g0_eqn = g0_eqn.subs(mu, mu_val)
    g0_eqn = g0_eqn.subs(nu, nu_val)

    g0_deriv = g0_deriv.subs(h, h_val)
    g0_deriv = g0_deriv.subs(mu, mu_val)
    g0_deriv = g0_deriv.subs(nu, nu_val)

    g0_eqn_set = lambdify([(g0, s)], [g0_eqn, g0_deriv], 'scipy')

    return g0_eqn_set

def cusp_solve_init(mu_val, nu_val):
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

    g0_deriv = diff(mut_g0_eqn, g0)

    g0_sec_deriv = diff(g0_deriv, g0)

    g0_eqn = mut_g0_eqn.subs(mu, mu_val)
    g0_eqn = g0_eqn.subs(nu, nu_val)

    g0_deriv = g0_deriv.subs(mu, mu_val)
    g0_deriv = g0_deriv.subs(nu, nu_val)

    g0_sec_deriv = g0_sec_deriv.subs(mu, mu_val)
    g0_sec_deriv = g0_sec_deriv.subs(nu, nu_val)

    g0_eqn_set = lambdify([(g0, s, h)], [g0_eqn, g0_deriv, g0_sec_deriv], 'scipy')

    return g0_eqn_set




iterations = 1000

s_val_range = np.logspace(-6, -5, iterations)
h_val = 1
mu_val = 5e-8
nu_val = 1e-9

cusp_eqn_set = cusp_solve_init(mu_val, nu_val)
cusp_soln = optimize.fsolve(cusp_eqn_set, [.5, .000001, .9])
print(cusp_soln)

bifn_eqn_set = bifn_solve_init(h_val, mu_val, nu_val)
bifn_soln = optimize.fsolve(bifn_eqn_set, [.5, .000001])
print(bifn_soln)
bifn_soln_2 = optimize.fsolve(bifn_eqn_set, [0, .00001])
print(bifn_soln_2)

for i in range(iterations):
    g0_eqn_numeric = root_eval_init(s_val_range[i], h_val, mu_val, nu_val)
    soln = optimize.root_scalar(g0_eqn_numeric, bracket = [.5, 1], method='brentq')


