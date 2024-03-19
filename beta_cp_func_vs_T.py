import numpy as np

# function of cp vs T in Kelvin
def cp_disc(x):
    if 273 <= x <= 600:
        result = 398.18776737 + -0.02021752*x + 0.00057432*x**2
    elif 600 < x <= 700:
        result = 1008.8601127 + -0.69177032*x
    else:
        result = 497.38204089 + -3.70561459e-03*x + 5.60806832e-05*x**2

    return result

# function of beta vs T in Kelvin
def beta_disc(x):
    if 273 <= x <= 631:
        return 1e-6*(46.29391086 + -1.80692616e-02*x + -6.04268663e-05*x**2 + 1.76344163e-07*x**3)
    elif 631 < x <= 688:
        return 1e-6*(125.72851444 + -0.11184634*x)
    else:
        return 1e-6*(48.50056087 + -1.60505423e-02*x + 2.33002672e-05*x**2)
    
# function (1 + beta)cp 
def f_beta_cp(x):
    if 273 <= x <= 600: #(1 + d)a
        return 1.01278e-16*x**5 - 3.82696e-14*x**4 + 6.10622e-11*x**3 + 0.000574323*x**2 - 0.0202257*x + 398.206
    elif 600 < x <= 631: #(1 + d)b
        return -1.2199e-13*x**4 + 2.19708e-10*x**3 - 4.84625e-8*x**2 - 0.691821*x + 1008.91
    elif 631 < x <= 688: #(1 + e)b
        return 7.7372e-8*x**2 - 0.69197*x + 1008.99
    elif 688 < x <= 700: #(1 + f)b
        return -1.61184e-11*x**3 + 3.461e-8*x**2 - 0.69182*x + 1008.91
    else: #(1 + f)c
        return 1.30669e-15*x**4 - 9.86467e-13*x**3 + 0.0000560951*x**2 - 0.00371378*x + 497.406
    

# integrated function (1 + beta)cp
def int_f(T_lower, T_upper):
    # integration in range 273 - 600K
    def int_a(x):
        return 1.01278e-16/6*x**6 - 3.82696e-14/5*x**5 + 6.10622e-11/4*x**4 + 0.000574323/3*x**3 - 0.0202257/2*x**2 + 398.206*x
    
    # integration in range 600 - 631K
    def int_b(x):
        return -1.2199e-13/5*x**5 + 2.19708e-10/4*x**4 - 4.84625e-8/3*x**3 - 0.691821/2*x**2 + 1008.91*x
    
    # integration in range 631K - 688K
    def int_c(x):
        return 7.7372e-8/3*x**3 - 0.69197/2*x**2 + 1008.99*x
    
    # integration in range 688K - 700K
    def int_d(x):
        return -1.61184e-11/4*x**4 + 3.461e-8/3*x**3 - 0.69182/2*x**2 + 1008.91*x
    
    # integration in range >700K
    def int_e(x):
        return 1.30669e-15/5*x**5 - 9.86467e-13/4*x**4 + 0.0000560951/3*x**3 - 0.00371378/2*x**2 + 497.406*x
    
    # --------------------------------------------------- #

    thresh = [273,600,631,688,700]
    if T_lower < thresh[0]:
        raise ValueError(f"Lower bound out of acceptable range; {T_lower} < 273K")
    
    T_intermediates = np.zeros(0)
    for i in reversed(range(len(thresh))):
        if T_lower >= thresh[i]:
            T_intermediates = np.append(T_intermediates, T_lower)
            T_intermediates = np.append(T_intermediates, thresh[i+1:])
            int_idx = i
            break
    T_intermediates = np.append(T_intermediates, T_upper)

    ints = [int_a, int_b, int_c, int_d, int_e]
    ints = ints[int_idx:]

    sum1 = 0
    for i in range(len(ints)):
        sum1 += ints[i](T_intermediates[i+1])

    sum2 = 0
    for i in range(len(ints)):
        sum2 += ints[i](T_intermediates[i])
    
    return sum1 - sum2     

    