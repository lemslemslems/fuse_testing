import numpy as np

# function for evaluating cp at specific temperature with T in Kelvin
def cp_disc(x):
    if 273 <= x <= 600:
        result = 398.18776737 + -0.02021752*x + 0.00057432*x**2
    elif 600 < x <= 700:
        result = 1008.8601127 + -0.69177032*x
    else:
        result = 497.38204089 + -3.70561459e-03*x + 5.60806832e-05*x**2

    return result

# integrated function of cp
def int_cp(T_lower,T_upper):

    # integration in range 273 - 600K
    def int_a(x):
        return 398.18776737*x + -0.02021752/2*x**2 + 0.00057432/3*x**3
    
    # integration in range 600 - 700K
    def int_b(x):
        return 1008.8601127*x + -0.69177032/2*x**2
    
    # integration in range >700K
    def int_c(x):
        return 497.38204089*x + -3.70561459e-03/2*x**2 + 5.60806832e-05/3*x**3
    
    thresh = [273,600,700]
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

    ints = [int_a, int_b, int_c]
    ints = ints[int_idx:]

    sum1 = 0
    for i in range(len(ints)):
        sum1 += ints[i](T_intermediates[i+1])

    sum2 = 0
    for i in range(len(ints)):
        sum2 += ints[i](T_intermediates[i])
    
    return sum1 - sum2      