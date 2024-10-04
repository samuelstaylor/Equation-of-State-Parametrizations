# Consequences of Fitting a Time-varying Equation of State Parameter Dark Energy Models with a Constant
# Author: Samuel Taylor


from models.const_w import ConstantEquationState
from models.linder import Linder
from models.hilltopK2 import HilltopK2
from models.hilltopK3 import HilltopK3
from models.hilltopK4 import HilltopK4

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.tri import Triangulation
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'

HUBBLE_CONSTANT = .7
OMEGA_M = .3
OMEGA_PHI = .7


# Prints a table of numbers for a given model, including redshift (z) and distance modulus (mu). Followed by w*
def print_data(model):
    print(f"{model.name} model Table of Numbers")
    print("[  z  ][   mu   ]")
    for i, z in enumerate(model.z_list):
        print(f"[{model.z_list[i]:.3F}][{model.dist_mod_list[i]:.5f}]")
    print(f"Corresponding best-fit const w for this model: {model.best_w:.4f}")
    print(f"Corresponding z-pivot for this model: {model.z_pivot:.4f}\n")



# Prints a table of numbers for a given model with specific redshift values (0.002, 0.5, 1, 1.5, 1.998).
def print_key_data(model):
    print(f"{model.name} w_o = {model.w_o} model Table of Numbers")
    print("[  z  ][   mu   ]")
    for i, z in enumerate(model.z_list):
        if (z==0.002 or z ==0.5 or z == 1 or z==1.5 or z == 1.998):
            print(f"[{model.z_list[i]:.3F}][{model.dist_mod_list[i]:.5f}]")
    print(f"Corresponding best-fit const w for this model: {model.best_w:.4f}")
    print(f"Corresponding z-pivot for this model: {model.z_pivot:.4f}\n")

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
    print(f"Best fit w for {model.label}: w={model.best_w:.2f}")
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
    
    plt.title(f"$w_{{o}}$ and Corresponding Best-fit $w_*$ for {models[0].name} Model")
    plt.xlabel(r"$w_{{o}}$", fontsize=12)
    plt.ylabel(r"$w_*$", fontsize=12)
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
        best_w_m2.append(models2[i].best_w)
        i += 1
    
    plt.figure()
    plt.plot(model_w, best_w_m1, label=models1[0].name)
    plt.plot(model_w, best_w_m2, label=models2[0].name)
    plt.title(f"w_o and Corresponding Best-fit w for {models1[0].name} and {models2[0].name} Models")
    plt.xlabel("w_o")
    plt.ylabel("w_*")
    plt.legend()
    plt.show(block=True)
    
    
# Plots w_o vs corresponding best-fit w for models from three different sets.
def plot_w_vs_const_w_3models(models1,models2,models3):
    model_w = []
    best_w_m1 = []
    best_w_m2 = []
    best_w_m3 = []

    i = 0
    
    while (i < len(models1) or i < len(models2) or i < len(models3)):
        model_w.append(models1[i].w_o)
        best_w_m1.append(models1[i].best_w)
        best_w_m2.append(models2[i].best_w)
        best_w_m3.append(models3[i].best_w)
        i += 1
    
    plt.figure()
    plt.plot(model_w, best_w_m1, label=models1[0].name)
    plt.plot(model_w, best_w_m2, label=models2[0].name)
    plt.plot(model_w, best_w_m3, label=models3[0].name)
    #plt.title(f"$w_{{o}}$ and Corresponding Best-fit $w_*$ for {models1[0].name}, {models2[0].name}, and {models3[0].name} Models", fontsize=14)
    plt.xlabel(r"$w_{{o}}$", fontsize=16,fontweight='bold')
    plt.ylabel(r"$w*$", fontsize=16,fontweight='bold')
    plt.legend()
    plt.savefig("images/hilltop_models_best_fit.png")
    plt.show(block=True)
    

def plot_linder_best_fit_w(models):
    model_w_o = []
    model_w_a = []
    best_w = []
    
    for model in models:
        model_w_o.append(model.w_o)
        model_w_a.append(model.w_a)
        best_w.append(model.best_w)
        print(f"w_o = {model.w_o:.4f}, w_a = {model.w_a:.4f}: w* = {model.best_w:.4f}")

    WO = np.linspace(min(model_w_o), max(model_w_o))
    WA = np.linspace(min(model_w_a), max(model_w_a))

    # Interpolate to grid
    W = griddata((model_w_o, model_w_a), best_w, (WO[None,:], WA[:,None]), method='cubic')
    W = np.clip(W, np.min(best_w), np.max(best_w))

    # Plot contour
    plt.figure()
    plt.contour(WO, WA, W, levels=15, cmap=plt.cm.jet)
    #plt.colorbar(label = r"$w*$",fontsize=16)
    cbar = plt.colorbar(label=r"$w*$")
    cbar.set_label(r"$w*$", fontsize=16)
    #plt.title(f"{models[0].name} model w*")
    plt.xlabel(r"$w_{{o}}$", fontsize=16,fontweight='bold')
    plt.ylabel(r"$w_{{a}}$", fontsize=16,fontweight='bold')
    plt.savefig("images/linder_best_w.png")
    plt.show()


def plot_linder_z_pivot(models):
    model_w_o = []
    model_w_a = []
    z_pivot = []
    
    for model in models:
        model_w_o.append(model.w_o)
        model_w_a.append(model.w_a)
        z_pivot.append(model.z_pivot)
        print(f"w_o = {model.w_o:.4f}, w_a = {model.w_a:.4f}: z_pivot = {model.z_pivot:.4f}")

    WO = np.linspace(min(model_w_o), max(model_w_o))
    WA = np.linspace(min(model_w_a), max(model_w_a))

    # Interpolate to grid
    Z = griddata((model_w_o, model_w_a), z_pivot, (WO[None,:], WA[:,None]), method='cubic')
    Z = np.clip(Z, np.min(z_pivot), np.max(z_pivot))

    # Plot contour
    plt.figure()
    plt.contour(WO, WA, Z, levels=15, cmap=plt.cm.jet)
    cbar = plt.colorbar(label=r"$z_{pivot}$")
    cbar.set_label(r"$z_{pivot}$", fontsize=16)
    #plt.title(f"{models[0].name} model z_pivot values")
    plt.xlabel(r'$w_{{o}}$', fontsize = 16)
    plt.ylabel(r'$w_{{a}}$', fontsize = 16)
    plt.savefig("images/linder_z_pivot")
    plt.show()


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
    #plt.title(f"w_o and Corresponding z_pivot for {models1[0].name} and {models2[0].name} Models")
    plt.xlabel("w_o")
    plt.ylabel("z")
    plt.legend()
    plt.show(block=True)
    
    
# Plots w_o vs corresponding pivot redshift (z_pivot) for models from two different sets.
def plot_z_pivot_w_3models(models1,models2,models3):
    model_w = []
    z_pivot_m1 = []
    z_pivot_m2 = []
    z_pivot_m3 = []
    i = 0
    
    while (i < len(models1) or i < len(models2) or i < len(models3)):
        model_w.append(models1[i].w_o)
        z_pivot_m1.append(models1[i].z_pivot)
        z_pivot_m2.append(models2[i].z_pivot)
        z_pivot_m3.append(models3[i].z_pivot)
        i += 1
    
    plt.figure()
    plt.plot(model_w, z_pivot_m1, label=models1[0].name)
    plt.plot(model_w, z_pivot_m2, label=models2[0].name)
    plt.plot(model_w, z_pivot_m3, label=models3[0].name)
    #plt.title(f"$w_{{o}}$ and Corresponding $z_{{pivot}}$ {models1[0].name}, {models2[0].name}, and {models3[0].name} Models", fontsize=14)
    plt.xlabel(r"$w_{{o}}$", fontsize=16, fontweight='bold')
    plt.ylabel(r"$z_{pivot}$", fontsize=16,fontweight='bold')
    plt.legend()
    plt.savefig("images/hilltop_models_pivot_z.png")
    plt.show(block=True)
    

# Generates and returns Linder modelsfor a range of w_o and w_a
def linder_plot(z_list):
    linder_plots = []
    w_o = -.99
    while (w_o <= -.7):
        w_a = -.20
        while (w_a <= .2):
            if (round(w_a,2) != 0):
                linder_plots.append(Linder(z_list, w_o, w_a))
            w_a += .01
        w_o += .01
    return linder_plots


# Generates and returns HilltopK2 objects for a range of w values and given redshift values.
def hilltopk2plot(z_list):
    hilltopk2_plots = []
    w_o = -.99
    while (w_o <= -.7):
        hilltopk2_plots.append(HilltopK2(z_list, w_o))
        w_o += .01
    return hilltopk2_plots


# Generates and returns HilltopK3 objects for a range of w values and given redshift values.
def hilltopk3plot(z_list):
    hilltopk3_plots = []
    w_o = -.99
    while (w_o <= -.7):
        hilltopk3_plots.append(HilltopK3(z_list, w_o))
        w_o += .01
    return hilltopk3_plots


# Generates and returns HilltopK4 objects for a range of w values and given redshift values.
def hilltopk4plot(z_list):
    hilltopk4_plots = []
    w_o = -.99
    while (w_o <= -.7):
        hilltopk4_plots.append(HilltopK4(z_list, w_o))
        w_o += .01
    return hilltopk4_plots


def main():
    print("Starting program: Consequences of Fitting a Time-varying Equation of State Parameter Dark Energy Models with a Constant\n")
    z_list = np.arange(.002, 2, .002)
    
    
    ### HILLTOP ###
    
    print('Running... Fetching Hilltop K=2')
    hilltopK2_plots = hilltopk2plot(z_list)
    print('Running... Fetching Hilltop K=3')
    hilltopK3_plots = hilltopk3plot(z_list)
    print('Running... Fetching Hilltop K=4')
    hilltopK4_plots = hilltopk4plot(z_list)
    print('Done with Hilltop calculations')
    plot_w_vs_const_w_3models(hilltopK2_plots, hilltopK3_plots, hilltopK4_plots)  
    plot_z_pivot_w_3models(hilltopK2_plots, hilltopK3_plots, hilltopK4_plots)  
    
    ###    
    
    
    ### LINDER ###
    print('Running... Fetching Linder')
    linder_plots = linder_plot(z_list)
    print("Done with Linder calculations")
    plot_linder_best_fit_w(linder_plots)
    plot_linder_z_pivot(linder_plots)
    ###
    
    
if __name__ == "__main__":
    main()


# 9/19. SCHERRER HILLTOP PAPER K=3. do the derivation and get omega phi as a function of a. Do derivation. Then once you can get this in a full w= then compare and use to find the best fit w.
# do the same to find pivot a for the cubic model. 