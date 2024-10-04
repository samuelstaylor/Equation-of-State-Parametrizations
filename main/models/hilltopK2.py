from .model_methods import calculation, min_chi_2_best_w, HUBBLE_CONSTANT, OMEGA_M, OMEGA_PHI, math, quad


class HilltopK2:
    ### 1 + w = (1 + w_o)a^3
    def __init__(self, z_list, w_o):
        if (w_o == -1):
            raise ValueError("w_o can not = -1 in Hilltop Model.")
        self.z_list = z_list
        self.w_o = round(w_o, 2)
        self.distance_mod_list()
        self.find_best_w()
        self.name = "Hilltop (K = 2)"
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
        # Yes, this means w_o can NOT = -1
        self.z_pivot = (1/self.a_pivot) - 1
