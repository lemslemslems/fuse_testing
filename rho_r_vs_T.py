import numpy as np

# function of rho_r vs temperature in Kelvin
def rho_r_disc(x):
    if 273 <= x <= 630:
        result = -21.631256 + 1.65854877e-01*x + -3.42786503e-04*x**2 + 3.29913578e-07*x**3

    if x > 630:
        result = -3.57260069 + 6.24155480e-02*x + -1.68017495e-05*x**2 + -2.74364233e-10*x**3

    return result    

# integrated function of rho_r
def int_rho_r(T_lower, T_upper):

    # integration in range 273 - 630K
    def int_a(x):
        return -21.631256*x + 1.65854877e-01/2*x**2 + -3.42786503e-04/3*x**3 + 3.29913578e-07/4*x**4

    # integration in range > 630K
    def int_b(x):
        return -3.57260069*x + 6.24155480e-02/2*x**2 + -1.68017495e-05/3*x**3 + -2.74364233e-10/4*x**4


    thresh = [273,630]
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

    ints = [int_a, int_b]
    ints = ints[int_idx:]

    sum1 = 0
    for i in range(len(ints)):
        sum1 += ints[i](T_intermediates[i+1])

    sum2 = 0
    for i in range(len(ints)):
        sum2 += ints[i](T_intermediates[i])

    return sum1 - sum2  