# This file creates the class for different models of dark energy
from scipy.integrate import quad
import math

HUBBLE_CONSTANT = .7
OMEGA_M = .3
OMEGA_PHI = .7

# Calculate list of distance modulus values from redshift 0 < z < 2
# given OMEGA_M, OMEGA_PHI, and the constant equation of state parameter
def calculation(z_list, w):
    dist_mod_values_list = []
    integration_val = 0
    luminosity_dist = 0
    distance_modulus = 0
    
    for z in z_list:
        f = lambda x: 1 / ((((OMEGA_M) * ((1.0 + x) ** 3)) + ((OMEGA_PHI) * ((1.0 + x) ** (3.0 *(1.0+w))))) ** 0.5)
        integration_val = quad(f, 0, z)[0] # integrate fdx from 0 to z
        luminosity_dist = (1 + z) * float(integration_val) 
        distance_modulus = (5 * (math.log(luminosity_dist, 10))) + 42.38 - (5 * (math.log(HUBBLE_CONSTANT, 10)))
        # Create list of all distance modulus values
        dist_mod_values_list.append(distance_modulus)
    return dist_mod_values_list



def chi_squared(actual, expected):
    chi2 = 0
    index=0
    while (index < len(actual) or index < len(expected)):
        chi2 += ( ((actual[index] - expected[index]) **2) * (.02))
        index += 1
    return chi2



            
def min_chi_2_best_w(z_list, dist_mod_list, w_val = 0, k_model = 0):
    lower_w = -1.2 #starting equation of state value to search from
    upper_w = -.6 
    
    dx = 0.000000001 #this is the incremental value to see if chi_X is inc or dec

    i = 0
    while i < 100:
        mid_w = (lower_w + upper_w)/2 
        mid_w_l = mid_w - dx
        mid_w_r = mid_w + dx
        
        dist_mod_w = calculation(z_list, mid_w)
        dist_mod_l = calculation(z_list, mid_w_l)
        dist_mod_r = calculation(z_list, mid_w_r)
        
        chi_x_w = chi_squared(dist_mod_list, dist_mod_w)
        chi_x_l = chi_squared(dist_mod_list, dist_mod_l)
        chi_x_r = chi_squared(dist_mod_list, dist_mod_r)
        
        if (chi_x_l < chi_x_w):
            upper_w = mid_w
        
        elif (chi_x_r < chi_x_w):
            lower_w = mid_w
        else:
            return mid_w  #returns this value if you found the valley peak before leaving the loop!
        
        i += 1

    best_w = mid_w       
    return best_w

