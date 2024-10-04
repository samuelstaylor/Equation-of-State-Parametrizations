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

# Create a filled contour plot with color mapping
contour = plt.contourf(X, Y, Z, cmap='viridis', levels=20)

# Add contour lines
contour_lines = plt.contour(X, Y, Z, colors='black', levels=20, linewidths=0.5)

# Add a colorbar
cbar = plt.colorbar(contour, label='Function Value')

# Add labels and title
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Contour Plot of f(x, y)')

# Add grid lines
plt.grid(True, linestyle='--', alpha=0.7)

# Customize the colorbar ticks
cbar.set_ticks(np.linspace(Z.min(), Z.max(), 5))

# Show the plot
plt.show()
