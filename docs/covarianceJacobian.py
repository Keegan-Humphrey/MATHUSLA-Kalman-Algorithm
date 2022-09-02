# This script finds the jacobian required to calculate the covariance matrix under the coordinate transformation from the measured seed values to the state vector with theta and phi

import sympy as sym

if __name__ == "__main__":
    x1 = sym.Symbol('x_1')
    t1 = sym.Symbol('t_1')
    z1 = sym.Symbol('z_1')
    y1 = sym.Symbol('y_1')
    x2 = sym.Symbol('x_2')
    t2 = sym.Symbol('t_2')
    z2 = sym.Symbol('z_2')
    y2 = sym.Symbol('y_2')

    theta = sym.Symbol('theta')
    phi = sym.Symbol('phi')

    p_hat = [x1, t1, z1, y1, x2, t2, z2, y2]
    s_hat = [x1, t1, z1, sym.atan(sym.sqrt((x2-x1)**2 + (z2-z1)**2)/(y2 - y1)),
    sym.atan((z2-z1)/(x2-x1))]
    print("Old Reduced Coordinates")
    for i in range(len(s_hat)):
        for j in range(len(p_hat)):
            value = sym.diff(s_hat[i], p_hat[j])
            print("\tJ_{" + str(i + 1) + str(j + 1) + "} &= " + str(value.simplify()) + " \\\\")

    thx = sym.Symbol('thx')
    thz = sym.Symbol('thz')

    deltax = x2-x1
    deltay = y2-y1
    deltaz = z2-z1

    tanx = sym.Symbol('\\tan(\\theta_x)')
    tanz = sym.Symbol('\\tan(\\theta_z)')

    s_hat = [x1, t1, z1, deltax/deltay, deltaz/deltay]
    m_hat = [x1, t1, z1, y1, x2, t2, z2, y2]

    
    print("NEW Theta_x and Theta_z Coordinates")
    for i in range(len(s_hat)):
        for j in range(len(m_hat)):
            value = sym.diff(s_hat[i], m_hat[j])
            value.simplify
            value.subs(deltax/deltay, tanx)
            value.subs(deltaz/deltay, tanz)
            print("\tJ_{" + str(i + 1) + str(j + 1) + "} &= " + str(value) + " \\\\")

    
    x = sym.Symbol('x')
    y = sym.Symbol('y')
    z = sym.Symbol('z')
    t = sym.Symbol('t')

    tanx = sym.Symbol('\\tan(\\theta_x)')
    tanz = sym.Symbol('\\tan(\\theta_z)')

    beta = sym.Symbol('\\beta')

    N = sym.sqrt(1 + tanx**2 + tanz**2)
    c = sym.Symbol('c')
    
    s_hat = [x,t,z,tanx,tanz]
    v_hat = [x,y,z,beta*c*tanx/N,beta*c/N,beta*c*tanz/N,t]
    print("NEW State vector to xyz v t coordinates")
    for i in range(len(v_hat)):
        for j in range(len(s_hat)):
            value = sym.diff(v_hat[i], s_hat[j])
            print("\tJ_{" + str(i + 1) + str(j + 1) + "} &= " + str(value.simplify()) + " \\\\")
