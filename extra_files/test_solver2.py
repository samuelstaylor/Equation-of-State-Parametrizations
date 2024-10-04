import sympy as sp

# Define the symbols
x, C, y, o = sp.symbols('x C y o')

# Define the equation
equation = (1 - C) * x**3 + C * x**6 - (1 + y) / (1 + o)

# Substitute specific values for C, y, and o
C_value = 0.7
y_value = -.8
o_value = -.7

equation_substituted = equation.subs({C: C_value, y: y_value, o: o_value})

# Solve the equation for x
solutions = sp.solve(equation_substituted, x)

# Print only the real solutions
real_solutions = [sol.evalf() for sol in solutions if sol.is_real]
print("Real solutions for x:")
for sol in real_solutions:
    if (sol >= 0):
        print(sol)