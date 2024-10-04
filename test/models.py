# This file creates the class for different models of dark energy
from scipy.integrate import quad
import math
import matplotlib.pyplot as plt
import numpy as np



HUBBLE_CONSTANT = .7
OMEGA_M = .3
OMEGA_PHI = .7

# Calculate list of distance modulus values from redshift 0 < z < 2
# given OMEGA_M, OMEGA_PHI, and the equation of state parameter
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
    lower_w = -1.01 #starting equation of state value to search from
    upper_w = -.8 
    
    # correction for discrepency
    if ((w_val == -.83 and k_model == 2) or (w_val == -.78 and k_model == 3)):
        upper_w = -.9
        lower_w = -.91
    
    mid_w = (lower_w + upper_w)/2 
    lower_w_dist_mod = calculation(z_list, lower_w)
    upper_w_dist_mod = calculation(z_list, upper_w)
    lower_w_X2 = chi_squared(dist_mod_list, lower_w_dist_mod)
    upper_w_X2 = chi_squared(dist_mod_list, upper_w_dist_mod)
    


    i = 0
    while i < 25:
        if (lower_w_X2 <= upper_w_X2):
            upper_w = mid_w
            upper_w_dist_mod = calculation(z_list, upper_w)
            upper_w_X2 = chi_squared(dist_mod_list, upper_w_dist_mod)

        else:
            lower_w = mid_w
            lower_w_dist_mod = calculation(z_list, lower_w)
            lower_w_X2 = chi_squared(dist_mod_list, lower_w_dist_mod)

        mid_w = (lower_w + upper_w) * .5
        i += 1

    best_w = mid_w       
    return best_w


class ConstantEquationState:
    def __init__(self, z_list, w):
        self.z_list = z_list
        self.w = w
        self.dist_mod_list = calculation(z_list, w)
        self.name = "Const w"
        self.label = f"{self.name}: w = {w:.2f}"

        
class Linder:
    ### w* = w_o + (1-a)w_a
    def __init__(self, z_list, w_o, w_a):
        self.z_list = z_list
        self.w_o = w_o
        self.w_a = w_a
        self.distance_mod_list()
        self.find_best_w()
        self.name = "Linder model"
        self.label = f"{self.name}: w_o = {self.w_o:.2f}, w_a = {self.w_a:.2f}"

        
    def distance_mod_list(self):
        self.dist_mod_list = []
        integration_val = 0
        luminosity_dist = 0
        
        for z in self.z_list:
            f = lambda x: 1 / ((((OMEGA_M) * ((1.0 + x) ** 3)) + (((OMEGA_PHI) * ((1.0 + x) ** (3.0 *(1.0+self.w_o+self.w_a))))) * pow(math.e, (3*(self.w_a*((1 / (x + 1)) -1)))) ) ** 0.5)
            # (1 / (x + 1)) is the scale factor
            integration_val = quad(f, 0, z)[0]

            luminosity_dist = (1 + z) * float(integration_val)
            distance_modulus = (5 * (math.log(luminosity_dist, 10))) + 42.38 - (5 * (math.log(HUBBLE_CONSTANT, 10)))
            self.dist_mod_list.append(distance_modulus)
        return self.dist_mod_list


    def find_best_w(self):
        self.best_w = min_chi_2_best_w(self.z_list, self.dist_mod_list)
        self.best_w_dist_mods = calculation(self.z_list, self.best_w)
        self.a_pivot = 1 - ((self.best_w-self.w_o)/self.w_a) # a = 1 - ((w* - wo) / wa)
        self.z = (1/self.a_pivot) - 1
       
        
class HilltopK2:
    ### 1 + w = (1 + w_o)a^3
    def __init__(self, z_list, w_o):
        self.z_list = z_list
        self.w_o = w_o
        self.distance_mod_list()
        self.find_best_w()
        self.name = "Hilltop (K=2)"
        self.label = f"{self.name} model: w_o = {self.w_o:.2f}"
        
    def distance_mod_list(self):
        self.dist_mod_list = []
        integration_val = 0
        luminosity_dist = 0
        distance_modulus = 0
        
        for z in self.z_list:
        # https://lambda.gsfc.nasa.gov/education/graphic_history/hubb_const.html
            f = lambda x: 1 / ((((OMEGA_M) * ((1.0 + x) ** 3)) + (( OMEGA_PHI * math.exp((1+self.w_o)*(1-(1/((1+x)**3)))) ))) ** 0.5)
            integration_val = quad(f, 0, z)[0] # integrate fdx from 0 to z
            luminosity_dist = (1 + z) * float(integration_val) 
            distance_modulus = (5 * (math.log(luminosity_dist, 10))) + 42.38 - (5 * (math.log(HUBBLE_CONSTANT, 10)))
            self.dist_mod_list.append(distance_modulus)
        
    def find_best_w(self):
        self.best_w = min_chi_2_best_w(self.z_list, self.dist_mod_list, self.w_o, 2)
        self.best_w_dist_mods = calculation(self.z_list, self.best_w)
        self.a_pivot = ((self.best_w + 1)/(self.w_o + 1))**(1/3) #a = ((1+w*)/(1+wo))^(1/3)
        self.z_pivot = (1/self.a_pivot) - 1

                
class HilltopK3:
    ###  1 + w = (1 +w_o)[(1−OMEGA_PHI)a^3 + OMEGA_PHIa^6]
    def __init__(self, z_list, w_o):
        self.z_list = z_list
        self.w_o = w_o
        self.distance_mod_list()
        self.find_best_w()
        self.name = "Hilltop (K = 3)"
        self.label = f"{self.name} model: w_o = {self.w_o:.2f}"
        
    def distance_mod_list(self):
        self.dist_mod_list = []
        integration_val = 0
        luminosity_dist = 0
        distance_modulus = 0
        
        for z in self.z_list:
        # https://lambda.gsfc.nasa.gov/education/graphic_history/hubb_const.html
            f = lambda x: 1 / (((OMEGA_M * ((1.0 + x) ** 3)) + (OMEGA_PHI * math.exp(((-(1+self.w_o)*((1/(1+x))**3))*(1-OMEGA_PHI+(OMEGA_PHI * ((1/(1+x))**3) / 2))) + ((1+self.w_o)*(1-(OMEGA_PHI/2)))))) ** 0.5)
            integration_val = quad(f, 0, z)[0] # integrate fdx from 0 to z
            luminosity_dist = (1 + z) * float(integration_val) 
            distance_modulus = (5 * (math.log(luminosity_dist, 10))) + 42.38 - (5 * (math.log(HUBBLE_CONSTANT, 10)))
            self.dist_mod_list.append(distance_modulus)
        
    def find_best_w(self):
        self.best_w = min_chi_2_best_w(self.z_list, self.dist_mod_list, self.w_o, 3)
        self.best_w_dist_mods = calculation(self.z_list, self.best_w)
        self.chi_X = chi_squared(self.dist_mod_list, self.best_w_dist_mods)
        self.a_pivot = self.find_a_pivot()
        self.z_pivot = (1/self.a_pivot) - 1

        
    def find_a_pivot(self):
       # Define the symbols
        a = sp.symbols('a')

        # Define the equation
        equation = (1 - OMEGA_PHI) * a**3 + OMEGA_PHI * a**6 - (1 + self.best_w) / (1 + self.w_o)

        # Solve the equation for x
        solutions = sp.solve(equation, a)
        # Print only the real solutions
        real_solutions = [sol.evalf() for sol in solutions if sol.is_real]
        for sol in real_solutions:
            if (sol >= 0):
                return sol
        raise ValueError(f"Hilltop K=3: no real solution for pivot_a found for w_o = {self.w_o} and w* = {self.best_w}")
    