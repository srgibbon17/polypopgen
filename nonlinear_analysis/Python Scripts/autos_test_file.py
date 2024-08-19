import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as optimize
from sympy import *

def root_solve_init(scaling_coeff):
    """
    Initializes an equation set which can then be evaluated for a set of parameters
    and passed to scipy.optimize to evaluate roots of the equation

    Args:
        scaling_coeff: a value equal to 1/(mu_val) that scales the systems of equatons to
                       facilitate root-finding by offsetting the smallness of the parameter inputs

    Returns:
        g0_eqn_params and g1_eqn_params: a function which takes seven inputs:
            s_val: a value of a selection coefficient
            h1_val: a value of a simplex dominance coefficient
            h2_val: a value of a duplex dominance coefficient
            h3_val: a value of a triplex dominance coefficient
            mu_val: a value for forward mutation
            nu_val: a value for back mutation
            a_val: a value for double reduction rate (alpha)
        and, after evaluating at a set of parameters, can be passed to lambdify.
        This "re-lambdified" function is then passed to scipy.optimize to calculate a root.
    """

    # initialize the necessary symbolic variables
        # g0 denotes frequency of the homozygous ancestral gamete
        # g1 denotes frequency of the heterozygous gamete
        # g2 denotes frequency of the homozygous derived gamete
        # G0 denotes frequency of the homozygous ancestral genotype (nulliplex)
        # G1 denotes frequency of the simplex heterozygous genotype
        # G2 denotes frequency of the duplex heterozygous genotype
        # G3 denotes frequency of the triplex heterozygous genotype
        # G4 denotes frequency of the homozygous derived genotype (quadriplex)
        # s denotes selection coefficient defined in a Wright-Fisher model
        # h1 denotes the dominance coefficient for simplex heterozygotes
        # h2 denotes the dominance coefficient for duplex heterozygotes
        # h3 denotes the dominance coefficient for triplex heterozygotes
        # mu denotes forward mutation rate from ancestral to derived
        # nu denotes back mutation rate from derived to ancestral
        # a denotes the rate of double reduction (alpha) as defined by Fisher 
    g0, g1, g2, G0, G1, G2, G3, G4, s, h1, h2, h3, mu, nu, a = symbols('g0 g1 g2 G0 G1 G2 G3 G4 s h1 h2 h3 mu nu a')

    # a set of substitutions based on definitions and HWE
    g2 = 1 - g0 - g1

    #g1 = 2*(g0**Rational(1/2))*(1-(g0**Rational(1/2)))

    G0 = g0**2
    G1 = 2*g0*g1
    G2 = 2*g0*g2 + g1**2
    G3 = 2*g1*g2
    G4 = g2**2

    # equations to parameterize relative fitnesses
    w_bar = 1 - s*(G1*h1 + G2*h2 + G3*h3 + G4) # average fitness
    w0 = 1/w_bar # fitness of ancestral homozygote (nulliplex)
    w1 = (1-h1*s)/w_bar # fitness of simplex heterozygote
    w2 = (1-h2*s)/w_bar # fitness of duplex heterozygote
    w3 = (1-h3*s)/w_bar # fitness of triplex heterozygote
    w4 = (1-s)/w_bar # fitness of derived homozygote (quadriplex)

    # equations to model selection and meiosis
    sel_g0_eqn = G0*w0 + (1/2 + a/4)*G1*w1 + (1/6 + a/3)*G2*w2 + (a/4)*G3*w3
    sel_g1_eqn = (1/2 - a/2)*G1*w1 + (2/3 - (2*a/3))*G2*w2 + (1/2 - a/2)*G3*w3
    sel_g2_eqn = (a/4)*G1*w1 + (1/6 + a/3)*G2*w2 + (1/2 + a/4)*G3*w3 + G4*w4

    # equations to model mutation
    # note: the scaling coefficient of 1e7 is necessary for numerical analysis techniques 
    # to converge appropriately; it has no effect on the roots of the equation
    mut_g0_eqn = scaling_coeff*(sel_g0_eqn*(1-mu)**2 + sel_g1_eqn*(1-mu)*nu + sel_g2_eqn*nu**2)
    mut_g1_eqn = scaling_coeff*(2*sel_g0_eqn*(1-mu)*mu + sel_g1_eqn*(1-mu)*(1-nu) + 2*sel_g2_eqn*(1-nu)*nu)
    mut_g2_eqn = scaling_coeff*(sel_g0_eqn*mu**2 + sel_g1_eqn*mu*(1-nu) + sel_g2_eqn*(1-nu)**2)


    # converts symbolic expression to a scipy interpretable function
    g0_eqn_params = lambdify([(s, h1, h2, h3, mu, nu, a)], mut_g0_eqn, 'scipy')
    g1_eqn_params = lambdify([(s, h1, h2, h3, mu, nu, a)], mut_g1_eqn, 'scipy')

    return g0_eqn_params, g1_eqn_params

def bifn_solve_init():
    """
    Initializes a set of equations which can then be evaluated for a set of parameters
    and passed to scipy.optimize to evaluate bifurcations of the system

    Args:
        none

    Returns:
        g0_eqn_params: a function which takes three inputs:
            h_val: a value of a dominance coefficient
            mu_val: a value for forward mutation
            nu_val: a value for back mutation

        g0_deriv_params: a function which takes three inputs: 
            h_val: a value of a dominance coefficient
            mu_val: a value for forward mutation
            nu_val: a value for back mutation   

        Note: g0_deriv_params is the first derivative of the function g0_eqn_params
                
        After evaluating these for a set of parameters, te equations set can be passed to lambdify.
        This "lambdified" function is then passed to scipy.optimize to calculate a bifurcation point.
    """

    # initialize the necessary symbolic variables
        # g0 denotes gamete frequency of the ancestral allele (i.e. p)
        # g1 denotes gamete frequency of the derived allele (i.e. q)
        # s denotes selection coefficient defined in a Wright-Fisher model
        # h denotes dominance coefficient defined in a Wright-Fisher model
        # mu denotes forward mutation rate from ancestral to derived
        # nu denotes back mutation rate from derived to ancestral
    g0, g1, s, h, mu, nu = symbols('g0 g1 s h mu nu')

    # replaces g1 in the following expressions with 1 - g0 (equivalent to q = 1-p)
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

    # calculate the symbolic derivative of the differential equation
    g0_deriv_sym = diff(mut_g0_eqn, g0)

    # converts symbolic expressions to scipy interpretable functions
    g0_eqn_params = lambdify([(h, mu, nu)], mut_g0_eqn, 'scipy')
    g0_deriv_params = lambdify([(h, mu, nu)], g0_deriv_sym, 'scipy')

    return g0_eqn_params, g0_deriv_params

def return_jacobian_elements(scaling_coeff):
    """
    Initializes four functions which represent the Jacobian entries of the system

    Args:
        scaling_coeff: a value equal to 1/(mu_val) that scales the systems of equatons to
                       facilitate root-finding by offsetting the smallness of the parameter inputs

    Returns:
        g0_deriv_g0: g0 difference expression derived against g0 
        g0_deriv_g1: g0 difference expression derived against g1
        g1_deriv_g0: g1 difference expression derived against g0
        g1_deriv_g1: g1 difference expression derived against g1

        The following are inputs to each of the above returns
            s_val: a value of a selection coefficient
            h1_val: a value of a simplex dominance coefficient
            h2_val: a value of a duplex dominance coefficient
            h3_val: a value of a triplex dominance coefficient
            mu_val: a value for forward mutation
            nu_val: a value for back mutation
            a_val: a value for double reduction rate (alpha)
        and, after evaluating at a set of parameters, can be passed to lambdify.
        This "re-lambdified" function is then passed to scipy.optimize as the Jacobian.
    """

    # initialize the necessary symbolic variables
        # g0 denotes frequency of the homozygous ancestral gamete
        # g1 denotes frequency of the heterozygous gamete
        # g2 denotes frequency of the homozygous derived gamete
        # G0 denotes frequency of the homozygous ancestral genotype (nulliplex)
        # G1 denotes frequency of the simplex heterozygous genotype
        # G2 denotes frequency of the duplex heterozygous genotype
        # G3 denotes frequency of the triplex heterozygous genotype
        # G4 denotes frequency of the homozygous derived genotype (quadriplex)
        # s denotes selection coefficient defined in a Wright-Fisher model
        # h1 denotes the dominance coefficient for simplex heterozygotes
        # h2 denotes the dominance coefficient for duplex heterozygotes
        # h3 denotes the dominance coefficient for triplex heterozygotes
        # mu denotes forward mutation rate from ancestral to derived
        # nu denotes back mutation rate from derived to ancestral
        # a denotes the rate of double reduction (alpha) as defined by Fisher 
    g0, g1, g2, G0, G1, G2, G3, G4, s, h1, h2, h3, mu, nu, a = symbols('g0 g1 g2 G0 G1 G2 G3 G4 s h1 h2 h3 mu nu a')

    # a set of substitutions based on definitions and HWE
    g2 = 1 - g0 - g1

    #g1 = 2*(g0**Rational(1/2))*(1-(g0**Rational(1/2)))

    G0 = g0**2
    G1 = 2*g0*g1
    G2 = 2*g0*g2 + g1**2
    G3 = 2*g1*g2
    G4 = g2**2

    # equations to parameterize relative fitnesses
    w_bar = 1 - s*(G1*h1 + G2*h2 + G3*h3 + G4) # average fitness
    w0 = 1/w_bar # fitness of ancestral homozygote (nulliplex)
    w1 = (1-h1*s)/w_bar # fitness of simplex heterozygote
    w2 = (1-h2*s)/w_bar # fitness of duplex heterozygote
    w3 = (1-h3*s)/w_bar # fitness of triplex heterozygote
    w4 = (1-s)/w_bar # fitness of derived homozygote (quadriplex)

    # equations to model selection and meiosis
    sel_g0_eqn = G0*w0 + (1/2 + a/4)*G1*w1 + (1/6 + a/3)*G2*w2 + (a/4)*G3*w3
    sel_g1_eqn = (1/2 - a/2)*G1*w1 + (2/3 - (2*a/3))*G2*w2 + (1/2 - a/2)*G3*w3
    sel_g2_eqn = (a/4)*G1*w1 + (1/6 + a/3)*G2*w2 + (1/2 + a/4)*G3*w3 + G4*w4

    # equations to model mutation
    # note: the scaling coefficient of 1e7 is necessary for numerical analysis techniques 
    # to converge appropriately; it has no effect on the roots of the equation
    mut_g0_eqn = (sel_g0_eqn*(1-mu)**2 + sel_g1_eqn*(1-mu)*nu + sel_g2_eqn*nu**2 - g0)
    mut_g1_eqn = (2*sel_g0_eqn*(1-mu)*mu + sel_g1_eqn*(1-mu)*(1-nu) + 2*sel_g2_eqn*(1-nu)*nu - g1)
    mut_g2_eqn = (sel_g0_eqn*mu**2 + sel_g1_eqn*mu*(1-nu) + sel_g2_eqn*(1-nu)**2 - g2)

    g0_symderiv_g0 = diff(mut_g0_eqn, g0)
    g0_symderiv_g1 = diff(mut_g0_eqn, g1)
    g1_symderiv_g0 = diff(mut_g1_eqn, g0)
    g1_symderiv_g1 = diff(mut_g1_eqn, g1)

    g0_deriv_g0 = lambdify([(s, h1, h2, h3, mu, nu, a)], g0_symderiv_g0, 'scipy')
    g0_deriv_g1 = lambdify([(s, h1, h2, h3, mu, nu, a)], g0_symderiv_g1, 'scipy')
    g1_deriv_g0 = lambdify([(s, h1, h2, h3, mu, nu, a)], g1_symderiv_g0, 'scipy')
    g1_deriv_g1 = lambdify([(s, h1, h2, h3, mu, nu, a)], g1_symderiv_g1, 'scipy')

    return g0_deriv_g0, g0_deriv_g1, g1_deriv_g0, g1_deriv_g1

g0, g1, s = symbols('g0 g1 s')

s_val = 2e-6
h1_val = 1
h2_val = 1
h3_val = 1
mu_val = 5e-8
nu_val = 1e-9
a_val = 0

def eval_jacobian(s_val, h1_val, h2_val, h3_val, mu_val, nu_val, a_val):
    jac_elements = return_jacobian_elements(1/mu_val)

    g0_deriv_g0 = lambdify([(g0, g1)], jac_elements[0]([s_val, h1_val, h2_val, h3_val, mu_val, nu_val, a_val]), 'scipy')
    g0_deriv_g1 = lambdify([(g0, g1)], jac_elements[1]([s_val, h1_val, h2_val, h3_val, mu_val, nu_val, a_val]), 'scipy')
    g1_deriv_g0 = lambdify([(g0, g1)], jac_elements[2]([s_val, h1_val, h2_val, h3_val, mu_val, nu_val, a_val]), 'scipy')
    g1_deriv_g1 = lambdify([(g0, g1)], jac_elements[3]([s_val, h1_val, h2_val, h3_val, mu_val, nu_val, a_val]), 'scipy')

    jacobian_func = np.array([[g0_deriv_g0, g0_deriv_g1], [g1_deriv_g0, g1_deriv_g1]])

    return jacobian_func

jacobian_func = eval_jacobian(s_val, h1_val, h2_val, h3_val, mu_val, nu_val, a_val)

# root_eqn_set = root_solve_init(1/mu_val)

# root_eqn_1_eval = root_eqn_set([s_val, h1_val, h2_val, h3_val, mu_val, nu_val, a_val])

# root_set_func = lambdify(g0, root_eqn_1_eval, 'scipy')

# root_soln_1 = optimize.root(root_set_func, .0001, method='hybr')
# #root_soln_2 = optimize.root(root_set_func, [.25, .5], method='krylov')
# #root_soln_3 = optimize.root(root_set_func, [.95, .0494], method='krylov')

root_eqn_set = root_solve_init(1/mu_val)

root_eqn_1_eval = root_eqn_set[0]([s_val, h1_val, h2_val, h3_val, mu_val, nu_val, a_val])
root_eqn_2_eval = root_eqn_set[1]([s_val, h1_val, h2_val, h3_val, mu_val, nu_val, a_val])

root_set_func = lambdify([(g0, g1)], [root_eqn_1_eval, root_eqn_2_eval], 'scipy')

#options = {'xtol':1e-8, 'eps':1e-10}

root_soln = optimize.fixed_point(root_set_func, [.01, .01])

#root_soln_1 = optimize.root(root_set_func, [.0001, .0198], method='krylov')
#root_soln_2 = optimize.root(root_set_func, [.25, .5], method='hybr', options=options)
#root_soln_3 = optimize.root(root_set_func, [.99, .01], method='hybr')

#print(root_soln_1)
#print(root_soln_2)
#print(root_soln_3)

#optimize.show_options('root', 'hybr')

def nsolve_roots(s_val, h1_val, h2_val, h3_val, mu_val, nu_val, a_val):
    g0, g1, g2, G0, G1, G2, G3, G4, s, h1, h2, h3, mu, nu, a = symbols('g0 g1 g2 G0 G1 G2 G3 G4 s h1 h2 h3 mu nu a')

    # a set of substitutions based on definitions and HWE
    g2 = 1 - g0 - g1
    


    G0 = g0**2
    G1 = 2*g0*g1
    G2 = 2*g0*g2 + g1**2
    G3 = 2*g1*g2
    G4 = g2**2

    # equations to parameterize relative fitnesses
    w_bar = 1 - s*(G1*h1 + G2*h2 + G3*h3 + G4) # average fitness
    w0 = 1/w_bar # fitness of ancestral homozygote (nulliplex)
    w1 = (1-h1*s)/w_bar # fitness of simplex heterozygote
    w2 = (1-h2*s)/w_bar # fitness of duplex heterozygote
    w3 = (1-h3*s)/w_bar # fitness of triplex heterozygote
    w4 = (1-s)/w_bar # fitness of derived homozygote (quadriplex)

    # equations to model selection and meiosis
    sel_g0_eqn = G0*w0 + (1/2 + a/4)*G1*w1 + (1/6 + a/3)*G2*w2 + (a/4)*G3*w3
    sel_g1_eqn = (1/2 - a/2)*G1*w1 + (2/3 - (2*a/3))*G2*w2 + (1/2 - a/2)*G3*w3
    sel_g2_eqn = (a/4)*G1*w1 + (1/6 + a/3)*G2*w2 + (1/2 + a/4)*G3*w3 + G4*w4

    # equations to model mutation
    # note: the scaling coefficient of 1e7 is necessary for numerical analysis techniques 
    # to converge appropriately; it has no effect on the roots of the equation
    mut_g0_eqn = (1/mu_val)*(sel_g0_eqn*(1-mu)**2 + sel_g1_eqn*(1-mu)*nu + sel_g2_eqn*nu**2 - g0)
    mut_g1_eqn = (1/mu_val)*(2*sel_g0_eqn*(1-mu)*mu + sel_g1_eqn*(1-mu)*(1-nu) + 2*sel_g2_eqn*(1-nu)*nu - g1)
    mut_g2_eqn = (1/mu_val)*(sel_g0_eqn*mu**2 + sel_g1_eqn*mu*(1-nu) + sel_g2_eqn*(1-nu)**2 - g2)

    mut_g0_subbed = mut_g0_eqn.subs([(s, s_val), (h1, h1_val), (h2, h2_val), (h3, h3_val), (mu, mu_val), (nu, nu_val), (a, a_val), (g1, 2*(g0**Rational(1/2))*(1-(g0**Rational(1/2))))])

    #mut_g1_subbed = mut_g1_eqn.subs([(s, s_val), (h1, h1_val), (h2, h2_val), (h3, h3_val), (mu, mu_val), (nu, nu_val), (a, a_val)])


    root = nsolve([mut_g0_subbed], [g0], [.1])

    return root