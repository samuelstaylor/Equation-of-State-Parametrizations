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

from models import ConstantEquationState
from models import Linder
from models import HilltopK2
from models import HilltopK3

HUBBLE_CONSTANT = .7
OMEGA_M = .3
OMEGA_PHI = .7


# Prints a table of numbers for a given model, including redshift (z) and distance modulus (mu). Followed by w*
def print_data(model):
    print(f"{model.name} model Table of Numbers")
    print("[  z  ][   mu   ]")
    for i, z in enumerate(model.z_list):
        print(f"[{model.z_list[i]:.3F}][{model.dist_mod_list[i]:.5f}]")
    print(f"This is the corresponding best-fit const w for this model: {model.best_w:.2f}")


# Prints a table of numbers for a given model with specific redshift values (0.002, 0.5, 1, 1.5, 1.998).
def print_key_data(model):
    print(f"{model.name} w_o = {model.w_o} model Table of Numbers")
    print("[  z  ][   mu   ]")
    for i, z in enumerate(model.z_list):
        if (z==0.002 or z ==0.5 or z == 1 or z==1.5 or z == 1.998):
            print(f"[{model.z_list[i]:.3F}][{model.dist_mod_list[i]:.5f}]")
    print(f"This is the corresponding best-fit const w for this model: {model.best_w:.2f}")


# Plots distance modulus (mu) vs redshift (z) for multiple models. (looks like log(n))
def plot(z_list, models):
    plt.figure()
    for model in models:
        plt.plot(z_list, model.dist_mod_list, label=model.label)
    plt.title("mu VS z for Varying Equation of State Models")
    plt.xlabel("z (Redshift values)")
    plt.ylabel("mu (Distance Modulus)")
    plt.legend()
    plt.show(block=True)


# Plots the best-fit distance modulus (mu) with the original model's mu vs redshift (z).
def plot_best_fit(z_list, model):
    plt.figure()
    plt.plot(z_list, model.dist_mod_list, label=model.label)
    plt.plot(z_list, model.best_w_dist_mods, label=f'Best fit w: {model.best_w:.2f}')
    plt.title("Best fit w for Equation of State Models")
    plt.xlabel("z (Redshift values)")
    plt.ylabel("mu (Distance Modulus)")
    plt.legend()
    print(f"Best fit w for {model.label}: w={model.best_w:.2f}, Chi-square: {model.chi_X:.8f}")
    plt.show(block=True)


# Plots the model w_o vs the corresponding best-fit w (w*) for the same model with different w_o.
def plot_w_vs_const_w(models):
    model_w = []
    best_w = []
    
    for model in models:
        model_w.append(model.w_o)
        best_w.append(model.best_w)
        print(f"w_o = {model.w_o:.4f} : w* = {model.best_w:.4f}")
    
    plt.figure()
    plt.plot(model_w, best_w)
    plt.title(f"w_o and Corresponding Best-fit w for {models[0].name} Model")
    plt.xlabel("w_o")
    plt.ylabel("w_*")
    plt.legend()
    plt.show(block=True)
 
 
# Plots w_o vs corresponding best-fit w for models from two different sets.
def plot_w_vs_const_w_2models(models1,models2):
    model_w = []
    best_w_m1 = []
    best_w_m2 = []
    i = 0
    
    while (i < len(models1) or i < len(models2)):
        model_w.append(models1[i].w_o)
        best_w_m1.append(models1[i].best_w)
        #print(f"w_o = {models1[i].w_o:.4f} : w* = {models1[i].best_w:.4f}")
        #print(f"pivot z = {models1[i].z_pivot}")
        best_w_m2.append(models2[i].best_w)
        #print(f"w_o = {models2[i].w_o:.4f} : w* = {models2[i].best_w:.4f}")
        #print(f"pivot z = {models2[i].z_pivot}")
        i += 1
    
    plt.figure()
    plt.plot(model_w, best_w_m1, label=models1[0].name)
    plt.plot(model_w, best_w_m2, label=models2[0].name)
    plt.title(f"w_o and Corresponding Best-fit w for {models1[0].name} and {models2[0].name} Models")
    plt.xlabel("w_o")
    plt.ylabel("w_*")
    plt.legend()
    plt.show(block=True)
    

# Plots w_o vs corresponding pivot redshift (z_pivot) for models from two different sets.
def plot_z_pivot_w_2models(models1,models2):
    model_w = []
    z_pivot_m1 = []
    z_pivot_m2 = []
    i = 0
    
    while (i < len(models1) or i < len(models2)):
        model_w.append(models1[i].w_o)
        z_pivot_m1.append(models1[i].z_pivot)
        z_pivot_m2.append(models2[i].z_pivot)
        i += 1
    
    plt.figure()
    plt.plot(model_w, z_pivot_m1, label=models1[0].name)
    plt.plot(model_w, z_pivot_m2, label=models2[0].name)
    plt.title(f"w_o and Corresponding z_pivot for {models1[0].name} and {models2[0].name} Models")
    plt.xlabel("w_o")
    plt.ylabel("z")
    plt.legend()
    plt.show(block=True)
    

# Generates and returns HilltopK2 objects for a range of w values and given redshift values.
def hilltopk2plot(z_list):
    hilltopk2_plots = []
    hilltopk2_7 = HilltopK2(z_list, -.7)
    hilltopk2_plots.append(hilltopk2_7)
    hilltopk2_71 = HilltopK2(z_list, -.71)
    hilltopk2_plots.append(hilltopk2_71)
    hilltopk2_72 = HilltopK2(z_list, -.72)
    hilltopk2_plots.append(hilltopk2_72)
    hilltopk2_73 = HilltopK2(z_list, -.73)
    hilltopk2_plots.append(hilltopk2_73)
    hilltopk2_74 = HilltopK2(z_list, -.74)
    hilltopk2_plots.append(hilltopk2_74)
    hilltopk2_75 = HilltopK2(z_list, -.75)
    hilltopk2_plots.append(hilltopk2_75)
    hilltopk2_76 = HilltopK2(z_list, -.76)
    hilltopk2_plots.append(hilltopk2_76)
    hilltopk2_77 = HilltopK2(z_list, -.77)
    hilltopk2_plots.append(hilltopk2_77)
    hilltopk2_78 = HilltopK2(z_list, -.78)
    hilltopk2_plots.append(hilltopk2_78)
    hilltopk2_79 = HilltopK2(z_list, -.79)
    hilltopk2_plots.append(hilltopk2_79)
    hilltopk2_8 = HilltopK2(z_list, -.8)
    hilltopk2_plots.append(hilltopk2_8)
    hilltopk2_81 = HilltopK2(z_list, -.81)
    hilltopk2_plots.append(hilltopk2_81)
    hilltopk2_82 = HilltopK2(z_list, -.82)
    hilltopk2_plots.append(hilltopk2_82)
    hilltopk2_83 = HilltopK2(z_list, -.83)
    hilltopk2_plots.append(hilltopk2_83)
    hilltopk2_84 = HilltopK2(z_list, -.84)
    hilltopk2_plots.append(hilltopk2_84)
    hilltopk2_85 = HilltopK2(z_list, -.85)
    hilltopk2_plots.append(hilltopk2_85)
    hilltopk2_86 = HilltopK2(z_list, -.86)
    hilltopk2_plots.append(hilltopk2_86)
    hilltopk2_87 = HilltopK2(z_list, -.87)
    hilltopk2_plots.append(hilltopk2_87)
    hilltopk2_88 = HilltopK2(z_list, -.88)
    hilltopk2_plots.append(hilltopk2_88)
    hilltopk2_89 = HilltopK2(z_list, -.89)
    hilltopk2_plots.append(hilltopk2_89)
    hilltopk2_9 = HilltopK2(z_list, -.9)
    hilltopk2_plots.append(hilltopk2_9)
    hilltopk2_91 = HilltopK2(z_list, -.91)
    hilltopk2_plots.append(hilltopk2_91)
    hilltopk2_92 = HilltopK2(z_list, -.92)
    hilltopk2_plots.append(hilltopk2_92)
    hilltopk2_93 = HilltopK2(z_list, -.93)
    hilltopk2_plots.append(hilltopk2_93)
    hilltopk2_94 = HilltopK2(z_list, -.94)
    hilltopk2_plots.append(hilltopk2_94)
    hilltopk2_95 = HilltopK2(z_list, -.95)
    hilltopk2_plots.append(hilltopk2_95)
    hilltopk2_96 = HilltopK2(z_list, -.96)
    hilltopk2_plots.append(hilltopk2_96)
    hilltopk2_97 = HilltopK2(z_list, -.97)
    hilltopk2_plots.append(hilltopk2_97)
    hilltopk2_98 = HilltopK2(z_list, -.98)
    hilltopk2_plots.append(hilltopk2_98)
    hilltopk2_99 = HilltopK2(z_list, -.99)
    hilltopk2_plots.append(hilltopk2_99)
    return hilltopk2_plots


# Generates and returns HilltopK3 objects for a range of w values and given redshift values.
def hilltopk3plot(z_list):
    hilltopK3_plots = []
    hilltopK3_7 = HilltopK3(z_list, -.7)
    hilltopK3_plots.append(hilltopK3_7)
    hilltopK3_71 = HilltopK3(z_list, -.71)
    hilltopK3_plots.append(hilltopK3_71)
    hilltopK3_72 = HilltopK3(z_list, -.72)
    hilltopK3_plots.append(hilltopK3_72)
    hilltopK3_73 = HilltopK3(z_list, -.73)
    hilltopK3_plots.append(hilltopK3_73)
    hilltopK3_74 = HilltopK3(z_list, -.74)
    hilltopK3_plots.append(hilltopK3_74)
    hilltopK3_75 = HilltopK3(z_list, -.75)
    hilltopK3_plots.append(hilltopK3_75)
    hilltopK3_76 = HilltopK3(z_list, -.76)
    hilltopK3_plots.append(hilltopK3_76)
    hilltopK3_77 = HilltopK3(z_list, -.77)
    hilltopK3_plots.append(hilltopK3_77)
    hilltopK3_78 = HilltopK3(z_list, -.78)
    hilltopK3_plots.append(hilltopK3_78)
    hilltopK3_79 = HilltopK3(z_list, -.79)
    hilltopK3_plots.append(hilltopK3_79)
    hilltopK3_8 = HilltopK3(z_list, -.8)
    hilltopK3_plots.append(hilltopK3_8)
    hilltopK3_81 = HilltopK3(z_list, -.81)
    hilltopK3_plots.append(hilltopK3_81)
    hilltopK3_82 = HilltopK3(z_list, -.82)
    hilltopK3_plots.append(hilltopK3_82)
    hilltopK3_83 = HilltopK3(z_list, -.83)
    hilltopK3_plots.append(hilltopK3_83)
    hilltopK3_84 = HilltopK3(z_list, -.84)
    hilltopK3_plots.append(hilltopK3_84)
    hilltopK3_85 = HilltopK3(z_list, -.85)
    hilltopK3_plots.append(hilltopK3_85)
    hilltopK3_86 = HilltopK3(z_list, -.86)
    hilltopK3_plots.append(hilltopK3_86)
    hilltopK3_87 = HilltopK3(z_list, -.87)
    hilltopK3_plots.append(hilltopK3_87)
    hilltopK3_88 = HilltopK3(z_list, -.88)
    hilltopK3_plots.append(hilltopK3_88)
    hilltopK3_89 = HilltopK3(z_list, -.89)
    hilltopK3_plots.append(hilltopK3_89)
    hilltopK3_9 = HilltopK3(z_list, -.9)
    hilltopK3_plots.append(hilltopK3_9)
    hilltopK3_91 = HilltopK3(z_list, -.91)
    hilltopK3_plots.append(hilltopK3_91)
    hilltopK3_92 = HilltopK3(z_list, -.92)
    hilltopK3_plots.append(hilltopK3_92)
    hilltopK3_93 = HilltopK3(z_list, -.93)
    hilltopK3_plots.append(hilltopK3_93)
    hilltopK3_94 = HilltopK3(z_list, -.94)
    hilltopK3_plots.append(hilltopK3_94)
    hilltopK3_95 = HilltopK3(z_list, -.95)
    hilltopK3_plots.append(hilltopK3_95)
    hilltopK3_96 = HilltopK3(z_list, -.96)
    hilltopK3_plots.append(hilltopK3_96)
    hilltopK3_97 = HilltopK3(z_list, -.97)
    hilltopK3_plots.append(hilltopK3_97)
    hilltopK3_98 = HilltopK3(z_list, -.98)
    hilltopK3_plots.append(hilltopK3_98)
    hilltopK3_99 = HilltopK3(z_list, -.99)
    hilltopK3_plots.append(hilltopK3_99)
    return hilltopK3_plots


def main():
    z_list = np.arange(.002, 2, .002)
    
    print('Running... Fetching Hilltop K=2')
    hilltopK2_plots = hilltopk2plot(z_list)
    print('Running... Fetching Hilltop K=3')
    hilltopK3_plots = hilltopk3plot(z_list)
    print('Done with Hilltop Z-pivot calculations')
    
    
    plot_w_vs_const_w_2models(hilltopK2_plots, hilltopK3_plots)  
    plot_z_pivot_w_2models(hilltopK2_plots, hilltopK3_plots)  
    
if __name__ == "__main__":
    main()


# 9/19. SCHERRER HILLTOP PAPER K=3. do the derivation and get omega phi as a function of a. Do derivation. Then once you can get this in a full w= then compare and use to find the best fit w.
# do the same to find pivot a for the cubic model. 