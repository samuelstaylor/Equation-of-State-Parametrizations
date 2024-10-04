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



# Calculate list of distance modulus values from redshift 0 < z < 2
# given omega_m, omega_delta, and the equation of state parameter
def calculation(z_list, omega_m, omega_delta, w):
    return dist_mod_values_list


def linder_calculate(z, omega_m, omega_delta, w_o, w_a):
    return distance_modulus
    

def non_linear_model_calculate(z_list, omega_m, omega_delta,  w_o):    
    return dist_mod_values_list


# Calculate list of distance modulus values using the Linder model
# given omega_m, omega_delta, and eq of state (VARYING WITH TIME)
def linder_distance_mod(z_list, omega_m, omega_delta, w_o, w_a):
    return linder_values_list


# Adds the labels to the graphs
def add_labels(plt, title = ("mu VS z for Varying Equation of State Models")):
    plt.title(title)
    plt.xlabel("z (Redshift values)")
    plt.ylabel("mu (Distance Modulus)")


def chi_squared(actual, expected):
    return chi2


def minimum_chi_squared(z_list, linder_dist_mods):
    return best_w


def plot_linder_and_min_w(z_list, w_o, w_a):
    return


def plot_nonlin_and_min_w(z_list, omega_m, omega_delta, w_o):
    return

    #delete this later this is just to print
def calculation_with_best_w(z_list, omega_m, omega_delta, w):
    return


def all_best_w(z_list, w_o_start=-1, w_o_end=-.5, w_o_step=.1, w_a_start=-1, w_a_end=1, w_a_step=.1):
    return


def plot_constant_w(z_list):
    return
    
    
    
def plot_linder_model(z_list):
    return
    

def plot_constant_linder(z_list):
    return

def main():
    return

if __name__ == '__main__':
    main()

