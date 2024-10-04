import sympy as sp

OMEGA_PHI = .7
best_w = .7
w_o = -.7

# Define the symbols
a = sp.symbols('a')

# Define the equation
equation = (1 - OMEGA_PHI) * a**3 + OMEGA_PHI * a**6 - (1 + best_w) / (1 + w_o)

# Solve the equation for x
solutions = sp.solve(equation, a)
# Print only the real solutions
real_solutions = [sol.evalf() for sol in solutions if sol.is_real]
print("Real solutions for x:")
for sol in real_solutions:
    if (sol >= 0):
        print(sol)