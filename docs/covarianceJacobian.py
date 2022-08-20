# This script finds the jacobian required to calculate the covariance matrix under the coordinate transformation from the measured seed values to the state vector with theta and phi

import sympy as sym

if __name__ == "__main__":
    x1 = sym.Symbol('x^{(1)}')
    t1 = sym.Symbol('t^{(1)}')
    z1 = sym.Symbol('z^{(1)}')
    y1 = sym.Symbol('y^{(1)}')
    x2 = sym.Symbol('x^{(2)}')
    t2 = sym.Symbol('t^{(2)}')
    z2 = sym.Symbol('z^{(2)}')
    y2 = sym.Symbol('y^{(2)}')

    theta = sym.Symbol('theta')
    phi = sym.Symbol('phi')

    p_hat = [x1, t1, z1, y1, x2, t2, z2, y2]
    s_hat = [x1, t1, z1, sym.atan(sym.sqrt((x2-x1)**2 + (z2-z1)**2)/(y2 - y1)),
    sym.atan((z2-z1)/(x2-x1))]

    for i in range(len(s_hat)):
        for j in range(len(p_hat)):
            value = sym.diff(s_hat[i], p_hat[j])
            print("J_{" + str(i + 1) + str(j + 1) + "} &= " + str(value.simplify()) + " \\\\")

