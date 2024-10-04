X2 = chi_squared(dist_mod_values_list_1, dist_mod_values_list8)
plt.plot(z_list, dist_mod_values_list_1, label="w = -1")
plt.plot(z_list, dist_mod_values_list8, label="w = -.8")
 
# creating error
y_error = .5 #size of error
error_on_z = [0.5, 1, 1.5]
error_on_dist_mod =[]
index = -1
for val in error_on_z:
    index += 25
    error_on_dist_mod.append(dist_mod_values_list_1[index])


plt.errorbar(error_on_z, error_on_dist_mod,
             yerr=y_error,
             fmt='o')
 
print("Chi-squared value:", X2)
add_labels(plt)
plt.legend()
plt.show()
 