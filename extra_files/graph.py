import matplotlib.pyplot as plt

# Read data from the file
data_file = "chiX2.txt"
with open(data_file, 'r') as file:
    data = [line.split() for line in file.readlines()]

# Extract x and y values
x_values = [float(line[0]) for line in data]
y_values = [float(line[1]) for line in data]

# Plot the data
plt.plot(x_values, y_values, label="Data")

# Set labels and title
plt.xlabel("w_o")
plt.ylabel("chi_x")
plt.title("Graph of chi_x Data for hilltopk2 w_o = -.99")

# Show the plot
plt.legend()
plt.grid(True)
plt.show()
