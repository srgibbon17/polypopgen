import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as optimize
from sympy import *

def func(x, *param_vals):
    # for diploids, define x as a vector of all possible parameters
    # x = [g0, s, h, mu, nu] such that:
        # x[0] = g0 = p = 1-g1 = 1-q
        # x[1] = s
        # x[2] = h
        # x[3] = mu
        # x[4] = nu
    # then the following sets of equations describe the change in gamete (allele) frequency
    # by modeling random mating, selection, meiosis, and mutation separately

    g0 = x
    g1 = 1-x
    s = param_vals[0][0]
    h = param_vals[0][1]
    mu = param_vals[0][2]
    nu = param_vals[0][3] 

    # equations to parameterize relative fitnesses
    w_bar = 1 - s*(2*h*g0*g1 + g1**2)# average fitness
    w0 = 1/w_bar # fitness of AA/ancestral homozygote
    w1 = (1-h*s)/w_bar # fitness of heterozygote
    w2 = (1-s)/w_bar # fitness of aa/derived homozygote

    # equations to model selection and meiosis
    sel_g0_eqn = w0*g0**2 + w1*g0*g1 
    sel_g1_eqn = w1*g0*g1 + w2*g1**2

    # equations to model mutation
    mut_g0_eqn = sel_g0_eqn*(1-mu) + sel_g1_eqn*nu - g0
    mut_g1_eqn = sel_g0_eqn*mu + sel_g1_eqn*(1-nu) - g1

    # because the dimensionality of the equations has been reduced using the linear relationship
    # g0 + g1 = 1 or g0 = 1 - g1, both mut_g0_eqn and mut_g1_eqn describe, in essence, the same dynamics
    # thus, only mut_g0_eqn is returned as the function

    return 1e7*(mut_g0_eqn)

def func_sym(x, *param_vals):
    g0, g1, s, h, mu, nu, x = symbols('g0 g1 s h mu nu x')

    # equations to parameterize relative fitnesses
    w_bar = 1 - s*(2*h*g0*g1 + g1**2) # average fitness
    w0 = 1/w_bar # fitness of AA/ancestral homozygote
    w1 = (1-h*s)/w_bar # fitness of heterozygote
    w2 = (1-s)/w_bar # fitness of aa/derived homozygote

    # equations to model selection and meiosis
    sel_g0_eqn = w0*g0**2 + w1*g0*g1 
    sel_g1_eqn = w1*g0*g1 + w2*g1**2

    # equations to model mutation
    mut_g0_eqn = sel_g0_eqn*(1-mu) + sel_g1_eqn*nu - g0
    mut_g1_eqn = sel_g0_eqn*mu + sel_g1_eqn*(1-nu) - g1

    mut_g0_eqn = mut_g0_eqn.subs(g1, 1-g0)
    g0_eqn_sym = mut_g0_eqn.subs(g0, x)
    g0_eqn_sym = g0_eqn_sym.subs(s, param_vals[0][0])
    g0_eqn_sym = g0_eqn_sym.subs(h, param_vals[0][1])
    g0_eqn_sym = g0_eqn_sym.subs(mu, param_vals[0][2])
    g0_eqn_sym = g0_eqn_sym.subs(nu, param_vals[0][3])

    print(g0_eqn_sym)

    g0_eqn_numeric = lambdify(x, 1e7*g0_eqn_sym, 'numpy')

    return g0_eqn_numeric

def func_fixed(x, *param_vals):
    g0 = x
    g1 = 1-x
    s = param_vals[0][0]
    h = param_vals[0][1]
    mu = param_vals[0][2]
    nu = param_vals[0][3]

    g0_eqn = g0**3*s*(2*h-1) + g0**2*s*(2-3*h+mu*h+nu*(1-h)) + g0*(-s*(1-h)+mu*(1-h*s)+nu*(1-2*s+h*s)) - nu*(1-s)

    return 1e7*g0_eqn

s_val = 2e-7
h_val = 1
mu_val = 5e-8
nu_val = 1e-9

param_vals = [s_val, h_val, mu_val, nu_val]

g0_0 = .4

soln = optimize.root_scalar(func_fixed, args=param_vals, x0=g0_0, method='secant')

print(soln)

#sol = optimize.fsolve(func, g0_0, args = param_vals)

#print(sol)

#sol_2 = optimize.fsolve(func, g0_0, args = param_vals)
#print(sol_2)

#g0_eqn_numeric = func_sym

#sol_3 = optimize.root_scalar(func_sym, args = param_vals, bracket = [0, 1], method='bisect')
#print(sol_3)

#x=0.8876250882518262717557968774926

#function = func(x, param_vals)

#print(function)