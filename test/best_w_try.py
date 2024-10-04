def min_chi_2_best_w(z_list, dist_mod_list):
    lower_w = -1.05 #starting equation of state value to search from
    upper_w = -.75 
    mid_w = (lower_w + upper_w)/2 
    lower_w_dist_mod = calculation(z_list, lower_w)
    upper_w_dist_mod = calculation(z_list, upper_w)
    lower_w_X2 = chi_squared(dist_mod_list, lower_w_dist_mod)
    upper_w_X2 = chi_squared(dist_mod_list, upper_w_dist_mod)

    i = 0
    while i < 1000:
        if (lower_w_X2 < upper_w_X2):
            upper_w = mid_w
            upper_w_dist_mod = calculation(z_list, upper_w)
            upper_w_X2 = chi_squared(dist_mod_list, upper_w_dist_mod)

        else:
            lower_w = mid_w
            lower_w_dist_mod = calculation(z_list, lower_w)
            lower_w_X2 = chi_squared(dist_mod_list, lower_w_dist_mod)

        mid_w = (lower_w + upper_w) / 2
        i += 1

    best_w = mid_w       
    return best_w



'''old chi_2 method'''

def min_chi_2_best_w(z_list, dist_mod_list):
    cur_w = -1.05 #starting equation of state value to search from
    upper_w = -.75 
    temp = calculation(z_list, cur_w)
    best_X2 = chi_squared(dist_mod_list,temp)
    best_w = cur_w
    
    while (cur_w <= upper_w):
        temp = calculation(z_list, cur_w)
        X2 = chi_squared(dist_mod_list,temp)
        if (X2 < best_X2):
            best_X2 = X2
            best_w = cur_w
        cur_w += .00001
    return best_w