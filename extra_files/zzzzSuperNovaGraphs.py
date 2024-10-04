# Consequences of Fitting a Time-varying Equation of State Parameter with a Constant
# Author: Samuel Taylor

# https://arxiv.org/abs/0803.0982
# ^ to read
'''
pg 19: shows different values of w_o and w_a and 

future graphs compare values within linder model
'''


'''
NOTES: 
w is equation of state 
omega is the density parameter 
p is density

m for matter
Φ for energy
delta for dark energy
'''

from scipy.integrate import quad
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

HUBBLE_CONSTANT = .7


# Calculate list of distance modulus values from redshift 0 < z < 2
# given omega_m, omega_delta, and the equation of state parameter
def calculation(z_list, omega_m, omega_delta, w):
    dist_mod_values_list = []
    integration_val = 0
    luminosity_dist = 0
    distance_modulus = 0
    for z in z_list:
        f = lambda x: 1 / ((((omega_m) * ((1.0 + x) ** 3)) + ((omega_delta) * ((1.0 + x) ** (3.0 *(1.0+w))))) ** 0.5)
        integration_val = quad(f, 0, z)[0] # integrate fdx from 0 to z
        luminosity_dist = (1 + z) * float(integration_val) 
        distance_modulus = (5 * (math.log(luminosity_dist, 10))) + 42.38 - (5 * (math.log(HUBBLE_CONSTANT, 10)))
        # Create list of all distance modulus values
        dist_mod_values_list.append(distance_modulus)
    return dist_mod_values_list


def linder_calculate(z, omega_m, omega_delta, w_o, w_a):
    hubble_const = .7
    integration_val = 0
    luminosity_dist = 0
    f = lambda x: 1 / ((((omega_m) * ((1.0 + x) ** 3)) + (((omega_delta) * ((1.0 + x) ** (3.0 *(1.0+w_o+w_a))))) * pow(math.e, (3*(w_a*((1 / (x + 1)) -1)))) ) ** 0.5)
    # (1 / (x + 1)) is the scale factor
    integration_val = quad(f, 0, z)[0]

    luminosity_dist = (1 + z) * float(integration_val)
    distance_modulus = (5 * (math.log(luminosity_dist, 10))) + 42.38 - (5 * (math.log(hubble_const, 10)))
    return distance_modulus
    

def non_linear_model_calculate(z_list, omega_m, omega_delta,  w_o):
    # this is the K= 2 in the scherrer dutta Hilltop quintessence paper
    # equation 34.
    
    H_o = 68
    p_c = (H_o **2) * ((3 * ((3 * 10 ** 8)**2)) / (8 * math.pi * (6.67 * (10 **-11)))) # critical density
    
    dist_mod_values_list = []
    integration_val = 0
    luminosity_dist = 0
    distance_modulus = 0
    
    
    print("Non-linear model Table of Numbers")
    print("Constants used:")
    print(f" - Hubble constant: H_o={H_o}")
    print(f" - Critical density: p_c={p_c}")
    print(f" - Density parameter of matter: omega_m={omega_m}")
    print(f" - Density parameter of dark energy: omega_delta={omega_delta}")
    print(f" - Given equation of state: w_o={w_o}")
    
    print("[  z  ][ D_l  ][   mu   ]")
  
    for z in z_list:
        # https://lambda.gsfc.nasa.gov/education/graphic_history/hubb_const.html
        f = lambda x: 1 / ((((omega_m) * ((1.0 + x) ** 3)) + (( 0.7 * math.exp((1+w_o)*(1-(1/((1+x)**3)))) ))) ** 0.5)
        integration_val = quad(f, 0, z)[0] # integrate fdx from 0 to z
        luminosity_dist = (1 + z) * float(integration_val) 
        distance_modulus = (5 * (math.log(luminosity_dist, 10))) + 42.38 - (5 * (math.log(HUBBLE_CONSTANT, 10)))
        # Create list of all distance modulus values
        
        #printing all for the w_o case
        dist_mod_values_list.append(distance_modulus)
        
    return dist_mod_values_list


# Calculate list of distance modulus values using the Linder model
# given omega_m, omega_delta, and eq of state (VARYING WITH TIME)
def linder_distance_mod(z_list, omega_m, omega_delta, w_o, w_a):
    linder_values_list = []
    #z_list = np.arange(.02, 2, .02) # Cannot equal 0 --> ERROR
    for z in z_list:
        linder_values_list.append(linder_calculate(z, omega_m, omega_delta, w_o, w_a))
    return linder_values_list


# Adds the labels to the graphs
def add_labels(plt, title = ("mu VS z for Varying Equation of State Models")):
    plt.title(title)
    plt.xlabel("z (Redshift values)")
    plt.ylabel("mu (Distance Modulus)")


def chi_squared(actual, expected):
    chi2 = 0
    index=0
    while (index < len(actual) or index < len(expected)):
        chi2 += ( ((actual[index] - expected[index]) **2) * (.02))
        index += 1
    return chi2


def minimum_chi_squared(z_list, linder_dist_mods):
    omega_m = 0.3
    omega_delta = 0.7
    w = -2
    dist_mod_values = calculation(z_list, omega_m, omega_delta, w)
    X2 = chi_squared(linder_dist_mods, dist_mod_values)
    best_w = -2
    X2_list = []
    while (w <= 2):
        w += 0.01
        dist_mod_values = calculation(z_list, omega_m, omega_delta, w)
        new_X2 = chi_squared(linder_dist_mods, dist_mod_values)
        X2_list.append(new_X2)
        if (new_X2 < X2):
            best_w = w
        X2 = new_X2
    return best_w


def plot_linder_and_min_w(z_list, w_o, w_a):
    # 9/19: has a pivot a when you do this. You should get the same value of a every time.
    #       SOLVE AND PRINT FOR A. YOU SHOULD GET THE SAME VALUE FOR A EVERY TIME. 
    
    linder_dist_mod = linder_distance_mod(z_list, 0.3, 0.7, w_o, w_a) # fix whats making calculations off. should be close to w_o
    best_w = minimum_chi_squared(z_list, linder_dist_mod)
    print(f"The constant w that minimizes chi squared for the linder model with w_o = {w_o} and w_a = {w_a} is: {best_w:.2f}")
    
    a_piv = 1 - ((best_w-w_o)/w_a)
    print(f"the pivot a value for w_o = {w_o} and w_a = {w_a} w_* = {best_w:.2f} is {a_piv:.4f}")
    
    dist_mod_values_list_best =calculation(z_list, 0.3, 0.7, best_w)
    plt.figure()
    plt.plot(z_list, linder_dist_mod, label=f"wo = {w_o}\nwa = {w_a}")
    plt.plot(z_list, dist_mod_values_list_best, label=f"w = {best_w:.2f}")
    add_labels(plt, "Minimized Chi-Squared")
    plt.legend()
    #plt.show(block=True)


def plot_nonlin_and_min_w(z_list, omega_m, omega_delta, w_o):
    non_lin_dist_mod  = non_linear_model_calculate(z_list, 0.3, 0.7,  w_o) 
    best_w = minimum_chi_squared(z_list, non_lin_dist_mod)
    print(f"The constant w that minimizes chi squared for the non-linear model with w_o = {w_o} and omega_delta = {omega_delta}  is: {best_w:.2f}")
    
    dist_mod_values_list_best =calculation(z_list, 0.3, 0.7, best_w)
    #delete this later, this is just to print
    calculation_with_best_w(z_list, omega_m, omega_delta, best_w)
    plt.figure()
    plt.plot(z_list, non_lin_dist_mod, label=f"wo = {w_o}\nomega_delta = {omega_delta}")
    plt.plot(z_list, dist_mod_values_list_best, label=f"w = {best_w:.2f}")
    add_labels(plt, "Minimized Chi-Squared")
    plt.legend()
    plt.show(block=True)


    #delete this later this is just to print
def calculation_with_best_w(z_list, omega_m, omega_delta, w):
    dist_mod_values_list = []
    integration_val = 0
    luminosity_dist = 0
    distance_modulus = 0
    print("Constant w Table of Numbers")
    print("Constants used:")
    print(f" - Equation of state param: w={w:.2f}")
    print(f" - Hubble constant: H_o={HUBBLE_CONSTANT}")
    print(f" - Density parameter of matter: omega_m={omega_m}")
    print(f" - Density parameter of dark energy: omega_delta={omega_delta}")

    for z in z_list:
        f = lambda x: 1 / ((((omega_m) * ((1.0 + x) ** 3)) + ((omega_delta) * ((1.0 + x) ** (3.0 *(1.0+w))))) ** 0.5)
        integration_val = quad(f, 0, z)[0] # integrate fdx from 0 to z
        luminosity_dist = (1 + z) * float(integration_val) 
        distance_modulus = (5 * (math.log(luminosity_dist, 10))) + 42.38 - (5 * (math.log(HUBBLE_CONSTANT, 10)))
        # Create list of all distance modulus values
        dist_mod_values_list.append(distance_modulus)
    return


def all_best_w(z_list, w_o_start=-1, w_o_end=-.5, w_o_step=.1, w_a_start=-1, w_a_end=1, w_a_step=.1):
    w_o_iter = w_o_start
    w_a_iter = w_a_start
    calc_a = True
    print("Starting Calculation for best constant w fit for Linder Models")
    if calc_a:
        print("and calculating a in the following equation: w_o + (1-a)w_a = w")
    while (w_o_iter <= w_o_end):
        w_a_iter = w_a_start
        while (w_a_iter <= w_a_end):
            linder_dist_mod = linder_distance_mod(z_list, 0.3, 0.7, w_o_iter, w_a_iter)
            best_w = minimum_chi_squared(z_list, linder_dist_mod)
            print(f"Best-fit constant w for Linder model with w_o={w_o_iter:.1f} and w_a={w_a_iter:.1f} is: w={best_w:.2f}")
            if calc_a and (w_a_iter > (w_a_step/2) or w_a_iter < (w_a_step/2)): # making sure w_a is not 0 for the calculation
                a = 1 - ((best_w - w_o_iter) / w_a_iter)
                print(f"Solving for a using those variables: a={a:.2f}")
            w_a_iter += w_a_step
        w_o_iter += w_o_step
    print("Finished: All best values for w have been calculated")
    return


def plot_constant_w(z_list):
    plt.figure()
    dist_mod_values_list_1 = calculation(z_list, 0.3, 0.7, -1) #default one... makes exponent zero so just adding omega_delta
    plt.plot(z_list, dist_mod_values_list_1, label="w = -1")
    dist_mod_values_list8 =calculation(z_list, 0.3, 0.7, -0.8)
    plt.plot(z_list, dist_mod_values_list8, label="w = -.8")
    dist_mod_values_list85 =calculation(z_list, 0.3, 0.7, -0.85)
    plt.plot(z_list, dist_mod_values_list85, label="w = -.85")
    dist_mod_values_list9 =calculation(z_list, 0.3, 0.7, -0.9)
    plt.plot(z_list, dist_mod_values_list9, label="w = -.9")
    dist_mod_values_list95 =calculation(z_list, 0.3, 0.7, -0.95)
    plt.plot(z_list, dist_mod_values_list95, label="w = -.95")
    add_labels(plt)
    plt.legend()
    plt.show(block=True)
    
    
    
def plot_linder_model(z_list):
    plt.figure()
    plt.plot(z_list, linder_distance_mod(z_list, 0.3, 0.7, -0.8, -0.5), label="wo = -0.8\nwa = -0.5")
    plt.plot(z_list, linder_distance_mod(z_list, 0.3, 0.7, -0.8, 0.5), label="wo = -0.8\nwa = 0.5")
    plt.plot(z_list, linder_distance_mod(z_list, 0.3, 0.7, -1, -0.5), label="wo = -1\nwa = -0.5")
    plt.plot(z_list, linder_distance_mod(z_list, 0.3, 0.7, -1, 0.5), label="wo = -1\nwa = 0.5")
    plt.plot(z_list, linder_distance_mod(z_list, 0.3, 0.7, -1, 0), label="wo = -1\nwa = 0")
    plt.plot(z_list, linder_distance_mod(z_list, 0.3, 0.7, -1, 2), label="wo = -1\nwa = 2")
    add_labels(plt)
    plt.legend()
    plt.show(block=True)
    

def plot_constant_linder(z_list):
    # w = -1
    # wo = -.8
    # wa = -.5
    plt.figure()
    dist_mod_values_list_1 = calculation(z_list, 0.3, 0.7, -1) #default one... makes exponent zero so just adding omega_delta
    linder_dist_mod1 = linder_distance_mod(z_list, 0.3, 0.7, -0.8, -0.5)
    X2 = chi_squared(dist_mod_values_list_1, linder_dist_mod1)
    plt.plot(z_list, dist_mod_values_list_1, label="w = -1")
    plt.plot(z_list, linder_dist_mod1, label="wo = -0.8\nwa = -0.5")
    plt.plot([1], [36], color = "white", label= f"\nX^2 = {X2:.2e}")
    add_labels(plt)
    plt.legend()
    plt.show(block=True)
    print(f"Chi-squared for constant w = -1 and linder wo = -.8 , wa = -.5: {X2:.2e}")
    
    # w = -0.9
    # wo = -.8
    # wa = -.5
    plt.figure()
    dist_mod_values_list9 =calculation(z_list, 0.3, 0.7, -0.9)
    linder_dist_mod1 = linder_distance_mod(z_list, 0.3, 0.7, -0.8, -0.5)
    X2 = chi_squared(dist_mod_values_list9, linder_dist_mod1)
    plt.plot(z_list, dist_mod_values_list9, label="w = -0.9")
    plt.plot(z_list, linder_dist_mod1, label="wo = -0.8\nwa = -0.5")
    plt.plot([1], [36], color = "white", label= f"\nX^2 = {X2:.2e}")
    add_labels(plt)
    plt.legend()
    plt.show(block=True)
    print(f"Chi-squared for constant w = -0.9 and linder wo = -.8 , wa = -.5: {X2:.2e}")


def main():
    z_list = np.arange(.002, 2, .002)
    
    #plot_nonlin_and_min_w(z_list, omega_m=0.3, omega_delta=0.7,  w_o=-1)
    #plot_nonlin_and_min_w(z_list, omega_m=0.3, omega_delta=0.7,  w_o=-.9)
    #plot_nonlin_and_min_w(z_list, omega_m=0.3, omega_delta=0.7,  w_o=-.8)
    #plot_nonlin_and_min_w(z_list, omega_m=0.3, omega_delta=0.7,  w_o=-.7)
    
    #plot_constant_w(z_list)
    
    #plot_linder_model(z_list)

    plot_linder_and_min_w(z_list, -0.9, -.5)
    plot_linder_and_min_w(z_list, -0.8, -.5)
    plot_linder_and_min_w(z_list, -0.7, -.5)
    plot_linder_and_min_w(z_list, -0.6, -.5)
    
    #plot_nonlin_and_min_w(z_list, omega_m=0.3, omega_delta=0.7,  w_o=-1)
    #plot_nonlin_and_min_w(z_list, omega_m=0.3, omega_delta=0.7,  w_o=-.95)
    #plot_nonlin_and_min_w(z_list, omega_m=0.3, omega_delta=0.7,  w_o=-.9)
    #plot_nonlin_and_min_w(z_list, omega_m=0.3, omega_delta=0.7,  w_o=-.85)
    #plot_nonlin_and_min_w(z_list, omega_m=0.3, omega_delta=0.7,  w_o=-.8)
    #plot_nonlin_and_min_w(z_list, omega_m=0.3, omega_delta=0.7,  w_o=-.75)
    
    



    
    
    ### prints the best constant w for varying linder models
    #all_best_w(z_list, w_o_start=-1, w_o_end=-.5, w_o_step=.1, w_a_start=-1, w_a_end=1, w_a_step=.1)
   
    

    
    # w_o from -1 to -.5 w/ steps of .1
    # w_a + from -1 to 1 w/ steps of .1   <<find the minimum value for each number
    
    
    
    #send set of numbers


if __name__ == '__main__':
    main()


# 9/19. SCHERRER HILLTOP PAPER K=3. do the derivation and get omega phi as a function of a. Do derivation. Then once you can get this in a full w= then compare and use to find the best fit w.
# do the same to find pivot a for the cubic model. 