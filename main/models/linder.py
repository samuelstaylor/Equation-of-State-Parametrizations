from .model_methods import calculation, min_chi_2_best_w, HUBBLE_CONSTANT, OMEGA_M, OMEGA_PHI, math, quad


class Linder:
    ### w* = w_o + (1-a)w_a
    def __init__(self, z_list, w_o, w_a):
        self.z_list = z_list
        self.w_o = round(w_o, 2)
        self.w_a = round(w_a, 2) 
        if (self.w_a == 0):
            raise ValueError("w_a can not = 0 in Hilltop Model.")
        self.distance_mod_list()
        self.find_best_w()
        self.name = "Linder"
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
        # w_a can NOT = 0
        self.z_pivot = (1/self.a_pivot) - 1
       