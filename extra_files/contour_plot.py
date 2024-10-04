import numpy as np
import matplotlib.pyplot as plt

# Define the function f(x, y)
def f(x, y):
    return x**2 + y**2  # Replace this with your actual function

# Generate x and y values
x_values = np.linspace(-5, 5, 100)
y_values = np.linspace(-5, 5, 100)

# Create a meshgrid from x and y values
X, Y = np.meshgrid(x_values, y_values)

# Calculate the function values for each point in the meshgrid
Z = f(X, Y)

# Create a contour plot
plt.contour(X, Y, Z, cmap='viridis')

# Add labels and a colorbar
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Contour Plot of f(x, y)')
plt.colorbar(label='Function Value')

# Show the plot
plt.show()
