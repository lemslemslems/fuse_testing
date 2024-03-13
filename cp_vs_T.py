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
    thresh = [273,600,700]

    for i in reversed(range(len(thresh))):
        if T_upper > thresh[i]:
            T_intermediates = np.zeros(i+1)
            break

    T_intermediates[0] = T_lower
    for i in range(len(T_intermediates)-1):
        T_intermediates[i] = thresh[i]
    T_intermediates[-1] = T_upper
    T_intermediates = np.append(T_intermediates, np.zeros(4 - len(T_intermediates)))

    int_T = np.zeros(4)
    for i in range(len(T_intermediates)):
        if 273 <= T_intermediates[i] <= 600:
            int_T[i] = 398.18776737*T_intermediates[i] + -0.02021752/2*T_intermediates[i]**2 + 0.00057432/3*T_intermediates**3
        elif 600 < T_intermediates[i] <= 700:
            int_T[i] = 1008.8601127*T_intermediates[i] + -0.69177032/2*T_intermediates[i]**2
        else:
            int_T[i] = 497.38204089*T_intermediates[i] + -3.70561459e-03/2*T_intermediates[i]**2 + 5.60806832e-05/2*T_intermediates[i]**2

    result = 0
    for i in range(len(int_T)):
        result += int_T[i] - int_T[i-1]

    return result       