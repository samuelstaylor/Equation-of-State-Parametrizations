import sympy as sp

# Define the symbols
x, C, y, o = sp.symbols('x C y o')

# Define the equation
equation = (1 - C) * x**3 + C * x**6 - (1 + y) / (1 + o)

# Solve the equation for x
solutions = sp.solve(equation, x)

# Print the solutions
print("Solutions for x:")
for sol in solutions:
    print(sol, "\n")


def create_data_file(z_list):
    cur_w = -1.05 #starting equation of state value to search from
    upper_w = -.75 
    
    with open("const_w_vals.dat", "w") as output_file:
        output_file.write("cur_w\tcalculation_result\n")  # Header line

        while (cur_w <= upper_w):
            dist_mod_list = calculation(z_list, cur_w)
            output_file.write(f"{cur_w:.5f}\t{dist_mod_list}\n")
            cur_w += .0001
            