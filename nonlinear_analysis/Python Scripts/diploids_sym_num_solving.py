import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as optimize
from sympy import *

def root_solve_init():
    """
    Initializes an equation which can then be evaluated for a set of parameters
    and passed to scipy.optimize to evaluate roots of the equation

    Args:
        none

    Returns:
        g0_eqn_params: a function which takes four inputs:
            s_val: a value of a selection coefficient
            h_val: a value of a dominance coefficient
            mu_val: a value for forward mutation
            nu_val: a value for back mutation
        and, after evaluating at a set of parameters, can be passed to lambdify.
        This "lambdified" function is then passed to scipy.optimize to calculate a root.
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
    w0 = 1/w_bar # fitness of ancestral homozygote
    w1 = (1-h*s)/w_bar # fitness of heterozygote
    w2 = (1-s)/w_bar # fitness of derived homozygote

    # equations to model selection and meiosis
    sel_g0_eqn = w0*g0**2 + w1*g0*g1 
    sel_g1_eqn = w1*g0*g1 + w2*g1**2

    # equations to model mutation
    # note: the scaling coefficient of 1e7 is necessary for numerical analysis techniques 
    # to converge appropriately; it has no effect on the roots of the equation
    mut_g0_eqn = 1e7*(sel_g0_eqn*(1-mu) + sel_g1_eqn*nu - g0)
    mut_g1_eqn = sel_g0_eqn*mu + sel_g1_eqn*(1-nu) - g1

    # converts symbolic expression to a scipy interpretable function
    g0_eqn_params = lambdify([(s, h, mu, nu)], mut_g0_eqn, 'scipy')

    return g0_eqn_params

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

def cusp_solve_init():
    """
    Initializes a set of equations which can then be evaluated for a set of parameters
    and passed to scipy.optimize to evaluate the cusp point of the system

    Args:
        none

    Returns:
        g0_eqn_params: a function which takes two inputs:
            mu_val: a value for forward mutation
            nu_val: a value for back mutation

        g0_deriv_params: a function which takes two inputs: 
            mu_val: a value for forward mutation
            nu_val: a value for back mutation  

        g0_sec_deriv_params: a function which takes two inputs: 
            mu_val: a value for forward mutation
            nu_val: a value for back mutation 

        Note: g0_deriv_params is the first derivative of the function g0_eqn_params and
              g0_sec_deriv_params is the second derivative of the function g0_eqn_params
                
        After evaluating these for a set of parameters, te equations set can be passed to lambdify.
        This "lambdified" function is then passed to scipy.optimize to calculate the cusp point.
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

    # calculate the symbolic derivatives of the differential equation
    g0_deriv_sym = diff(mut_g0_eqn, g0)
    g0_sec_deriv_sym = diff(g0_deriv_sym, g0)

    # converts symbolic expressions to scipy interpretable functions
    g0_eqn_params = lambdify([(mu, nu)], mut_g0_eqn, 'scipy')
    g0_deriv_params = lambdify([(mu, nu)], g0_deriv_sym, 'scipy')
    g0_sec_deriv_params = lambdify([(mu, nu)], g0_sec_deriv_sym, 'scipy')

    return g0_eqn_params, g0_deriv_params, g0_sec_deriv_params


def bifn_diagram_init(s_val_range, h_val, mu_val, nu_val):

    g0, s, h = symbols('g0 s h')

    # cusp equations
    cusp_eqn_params_1, cusp_eqn_params_2, cusp_eqn_params_3 = cusp_solve_init()

    # evaluates cusp equations at given mu and nu
    cusp_eqn_1_eval = cusp_eqn_params_1([mu_val, nu_val])
    cusp_eqn_2_eval = cusp_eqn_params_2([mu_val, nu_val])
    cusp_eqn_3_eval = cusp_eqn_params_3([mu_val, nu_val])

    # creates a vectorized set of equations with variables g0, s, and h callable by scipy.optimize
    cusp_set_func = lambdify([(g0, s, h)], [cusp_eqn_1_eval, cusp_eqn_2_eval, cusp_eqn_3_eval], 'scipy')

    # solves the set of equations for a point (g0, s, h)
    # initial guess is of the same form (g0, s, h)
    cusp_soln = optimize.root(cusp_set_func, [.5, 0, .75], method='hybr')

    # defines the critical h value at which the cusp occurs
    # for h_val > h_crit, two bifurcations occur
    # for h_val < h_crit, no bifurcations occur
    h_crit = cusp_soln.x[2]

    # if h>h_crit evaluates bifn equations
    if h_val > h_crit:
        # initializing bifn equations
        bifn_eqn_params_1, bifn_eqn_params_2 = bifn_solve_init()

        # evaluating at given h, mu, and nu
        bifn_eqn_1_eval = bifn_eqn_params_1([h_val, mu_val, nu_val])
        bifn_eqn_2_eval = bifn_eqn_params_2([h_val, mu_val, nu_val])

        # creates a vectorized equation set with variable g0 and s callable by scipy.optimize
        bifn_set_func = lambdify([(g0, s)], [bifn_eqn_1_eval, bifn_eqn_2_eval], 'scipy')

        # solves the set of equtions for two points of form (g0, s)
        # first initial condition finds bifn for higher g0 and lower s
        bifn_soln_1 = optimize.root(bifn_set_func, [.5, .000001], method='hybr')
        # second initial condition finds bifn for lower g0 and higher s
        bifn_soln_2 = optimize.root(bifn_set_func, [0, .00001], method='hybr')

    # initializing root/equilibria equation
    g0_eqn_params = root_solve_init()

    # numpy arrays to store values if no bifurcations occur
    eq_solns = np.zeros_like(s_val_range)
    s_vals = np.zeros_like(s_val_range)

    # numpy arrays to store stable eq. values under neutral evolution if bifurcations occur
    eq_solns_1 = np.zeros_like(s_val_range)
    s_vals_1 = np.zeros_like(s_val_range)

    # numpy arrays to store unstable eq. values if bifurcations occur
    eq_solns_2 = np.zeros_like(s_val_range)
    s_vals_2 = np.zeros_like(s_val_range)

    # numpy arrays to store stable eq. values under purifying selection if bifurcations occur
    eq_solns_3 = np.zeros_like(s_val_range)
    s_vals_3 = np.zeros_like(s_val_range)

    for i in range(len(s_val_range)):

        # evaluates the root equation for given s, h, mu, and nu
        g0_eqn_eval = g0_eqn_params([s_val_range[i], h_val, mu_val, nu_val])

        # creates a scalar equation with variable g0 callable by scipy.optimize
        g0_eqn_func = lambdify(g0, g0_eqn_eval, 'scipy')

        if h_val < h_crit:
            # implies that there are no bifurcations and will only be one (stable) equilibrium
            root_soln = optimize.root_scalar(g0_eqn_func, bracket = [0, 1], method ='toms748')
            # appends the root_soln and current s_val to respective arrays
            eq_solns[i] = root_soln.root
            s_vals[i] = s_val_range[i]

        else:
            # divides search range into three distinct "zones"
            # first zone is for small s wheere only neutral evolution occurs
            # second zone is for intermediate s where bistability occurs
            # third zone is for large s where only purifying selection occurs    
            if s_val_range[i] < bifn_soln_1.x[1]:
                # sets search range from zero up to the smaller g0 bifn value
                # root_soln_1 finds the neutral evolution, stable curve
                root_soln_1 = optimize.root_scalar(g0_eqn_func, bracket = [0, bifn_soln_2.x[0]], method ='toms748')
                # appends the root_soln and current s_val to respective arrays
                eq_solns_1[i] = root_soln_1.root
                s_vals_1[i] = s_val_range[i]

            elif s_val_range[i] > bifn_soln_1.x[1] and s_val_range[i] < bifn_soln_2.x[1]:
                # sets search range from zero up to the smaller g0 bifn value
                # root_soln_1 finds the neutral evolution, stable curve
                root_soln_1 = optimize.root_scalar(g0_eqn_func, bracket = [0, bifn_soln_2.x[0]], method ='toms748')
                # appends the root_soln and current s_val to respective arrays
                eq_solns_1[i] = root_soln_1.root
                s_vals_1[i] = s_val_range[i]

                # sets search range from smaller g0 bifn value to larger g0 bifn value
                # root_soln_2 finds the intermediate, unstable equilibria
                root_soln_2 = optimize.root_scalar(g0_eqn_func, bracket = [bifn_soln_2.x[0], bifn_soln_1.x[0]], method ='toms748')
                # appends the root_soln and current s_val to respective arrays
                eq_solns_2[i] = root_soln_2.root
                s_vals_2[i] = s_val_range[i]

                # sets search range from larger g0 bifn value to 1
                # root_soln_3 finds the purifying selection, stable curve
                root_soln_3 = optimize.root_scalar(g0_eqn_func, bracket = [bifn_soln_1.x[0], 1], method = 'toms748')
                # appends the root_soln and current s_val to respective arrays
                eq_solns_3[i] = root_soln_3.root
                s_vals_3[i] = s_val_range[i]

            elif s_val_range[i] > bifn_soln_2.x[1]:
                # sets search range from larger g0 bifn value to 1
                # root_soln_3 finds the purifying selection, stable curve
                root_soln_3 = optimize.root_scalar(g0_eqn_func, bracket = [bifn_soln_1.x[0], 1], method = 'toms748')
                # appends the root_soln and current s_val to respective arrays
                eq_solns_3[i] = root_soln_3.root
                s_vals_3[i] = s_val_range[i]

    eq_solns = np.trim_zeros(eq_solns)
    s_vals = np.trim_zeros(s_vals)
    eq_solns_1 = np.trim_zeros(eq_solns_1)
    s_vals_1 = np.trim_zeros(s_vals_1)
    eq_solns_2 = np.trim_zeros(eq_solns_2)
    s_vals_2 = np.trim_zeros(s_vals_2)
    eq_solns_3 = np.trim_zeros(eq_solns_3)
    s_vals_3 = np.trim_zeros(s_vals_3)
    

    if h_val < h_crit:
        return eq_solns, s_vals
    elif h_val > h_crit:
        return eq_solns_1, s_vals_1, eq_solns_2, s_vals_2, eq_solns_3, s_vals_3
        


iterations = 1000

h_val = 1
mu_val = 5e-8
nu_val = 1e-9
s_val_range = np.logspace(-7, -6, iterations)

outputs = bifn_diagram_init(s_val_range, h_val, mu_val, nu_val)

print(outputs)
