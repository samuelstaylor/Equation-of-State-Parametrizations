from .model_methods import calculation

class ConstantEquationState:
    def __init__(self, z_list, w):
        self.z_list = z_list
        self.w = round(w, 2) # rounding the value to fix for float conversion.
        self.dist_mod_list = calculation(z_list, w)
        self.name = "Const w"
        self.label = f"{self.name}: w = {w:.2f}"
