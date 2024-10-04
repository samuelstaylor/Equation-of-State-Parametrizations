
def all_best_w(z_list, w_o_start=-1, w_o_end=-.5, w_o_step=.1, w_a_start=-1, w_a_end=1, w_a_step=.1):
    w_o_iter = w_o_start
    w_a_iter = w_a_start
    calc_a = True
    print("Starting Calculation for best constant w fit for Linder Models")
    if calc_a:
        print("and calculating a in the following equation: w_o + (1-a)w_a = w")
    while (w_o_iter <= w_o_end):
        w_a_iter = w_a_start
        while (w_a_iter <= w_a_end):
            linder_dist_mod = linder_distance_mod(z_list, 0.3, 0.7, w_o_iter, w_a_iter)
            best_w = minimum_chi_squared(z_list, linder_dist_mod)
            print(f"Best-fit constant w for Linder model with w_o={w_o_iter:.1f} and w_a={w_a_iter:.1f} is: w={best_w:.2f}")
            if calc_a and (w_a_iter > (w_a_step/2) or w_a_iter < (w_a_step/2)): # making sure w_a is not 0 for the calculation
                a = 1 - ((best_w - w_o_iter) / w_a_iter)
                print(f"Solving for a using those variables: a={a:.2f}")
            w_a_iter += w_a_step
        w_o_iter += w_o_step
    print("Finished: All best values for w have been calculated")
    return