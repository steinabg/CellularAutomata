import numpy as np


def T2calc_change_qcj(pt, Dj, q_cbj, Ej, gamma, q_th):
    '''
    This function calculates formula (49)
    :param pt: Time step of CA
    :param Dj: Deposition rate
    :param q_cbj: jth sediment sea bed density
    :param Ej: Erosion rate
    :param gamma: Porosity
    :param q_th: Turbidity current thickness
    :return: The change in the jth current sediment concentration
    '''
    diff = (Dj - q_cbj * Ej)
    factor = pt/(1-gamma)
    with np.errstate(divide='ignore', invalid='ignore'):
        res = np.nan_to_num(factor * diff / q_th[:, :, np.newaxis])[1:-1, 1:-1, :]

    return res

def T2_calc_change_qCBJ(pt, Dj, q_cbj, Ej, gamma, q_d):
    '''

    :param pt: Time step of CA
    :param Dj: Deposition rate
    :param q_cbj: jth sediment sea bed density
    :param Ej: Erosion rate
    :param gamma: Porosity
    :param q_d: Soft sediment thickness
    :return: The change in the jth bed sediment concentration
    '''

    factor = pt / (1 - gamma)
    res1 = np.nan_to_num(factor * np.sum(Dj - q_cbj * Ej, axis=2))

    firstDiff = (factor / q_d)[:, :, np.newaxis] * (Dj - q_cbj * Ej)

    secondDiff = (res1 / q_d)[:, :, np.newaxis] * q_cbj
    return np.nan_to_num(firstDiff - secondDiff)[1:-1, 1:-1, :]

def T2_calc_change_qd(pt, Dj, q_cbj, Ej, gamma):
    '''

    :param pt: Time step of CA
    :param Dj: Deposition rate
    :param q_cbj: jth sediment sea bed density
    :param Ej: Erosion rate
    :param gamma: Porosity
    :return: The change in soft sediment thickness
    '''
    factor = pt / (1 - gamma)
    var = factor * np.sum(Dj - q_cbj * Ej, axis=2)
    return np.nan_to_num(var)[1:-1,1:-1]


def calc_erotionRate(Z_mj):
    '''
    :type Z_mj: numpy.ndarray(Ny,Nx,Nj)
    :param Z_mj: jth value of Z as specified by eq. (38)

    :rtype: numpy.ndarray(Ny,Nx,Nj)
    :return: Deposition rate
    '''
    return 1.3e-07 * Z_mj ** 5 / (1 + 4.3e-07 * Z_mj ** 5)


def calc_Z_mj(kappa, Ustar, v_sj, f):
    '''
    :type kappa: double
    :param kappa: The value of kappa defined in eq. (39)

    :param Ustar: The value of U* (see. eq (46)).
    :type Ustar: numpy.ndarray(Ny,Nx)

    :param v_sj: Dimensionless sphere settling velocities
    :type v_sj: numpy.ndarray(Nj)

    :type f: numpy.ndarray(Nj)
    :param f: This function returns the value defined by eq.(40) TODO!

    :rtype: numpy.ndarray(Ny,Nx,Nj)
    :return: jth value of Z as specified by eq. (38)

    '''
    return kappa * (np.sqrt(Ustar ** 2)[:, :, np.newaxis]) * v_sj * f


def calc_kappa(D_s):  # TODO! SJEKK!
    '''
    This function computes kappa equation (39)

    :type D_s: numpy.ndarray(Nj)
    :param D_s: Array of sediment-particle diameters [m]

    :return: The value of kappa defined in eq. (39)
    :rtype: double

    '''

    phi = np.log2(D_s)
    sigma_phi = np.std(phi)  # Calculate standard deviation of phi
    return 1 - 0.288 * sigma_phi


def calc_Ustar(c_D, q_v):
    '''
    This function calculates the value of U*, eq. (46)

    :type c_D: double
    :param c_D: Bed drag coefficient [unit = 1]
    :type q_V: numpy.ndarray(Ny,Nx)
    :param q_V: Speed of turbidity current.

    :return: The value of U* (see. eq (46)).
    :rtype: numpy.ndarray(Ny,Nx)

    '''
    return c_D * q_v


def calc_fofR(R_pj):
    '''
    This function returns the value defined by eq.(40) TODO!

    :type: numpy.array(Nj)
    :param: Particle Reynolds number for particle type j.

    :rtype: numpy.ndarray(Nj)
    :return: This function returns the value defined by eq.(40) TODO!
    '''
    if (R_pj < 1):
        print("Undefined function value for R_pj<1 !")
    return np.where(R_pj >= 3.5, R_pj ** (0.6), 0.586 * R_pj ** (1.23))


def calc_Rpj(rho_j, rho, D_sj, nu, g=9.81):
    '''
    This function calculates and returns the particle Reynolds number as given by\
    equation (40). Assume rho = rho_ambient.

    :type rho_j: numpy.ndarray(Nj)
    :param rho_j: Density of sediment type no j. [kg/m^3]

    :type rho: double
    :param rho: density in equation (40). Assumed to be ambient density. [kg/m^3]. TODO!

    :type nu: double
    :param nu: Kinematic viscosity [m^2/s]

    :rtype: numpy.array(Nj)
    :return: Particle Reynolds number for particle type j.
    '''
    return np.sqrt(g * (rho_j - rho) * D_sj / rho) * D_sj / nu


def calc_nearBedConcentration_SusSed(D_sj, D_sg, q_cj):
    '''
    This function calculates the near-bed concentration of suspended sediment.\
    Equation (45).


    :type D_sj: numpy.ndarray(Nj)
    :param D_sj: jth sediment diameter

    :type D_sg: numpy.ndarray(Ny,Nx) TODO! Skal den være det?
    :param D_sg: Geometric mean size of suspended sediment mixture in cell

    :type q_cj: numpy.ndarray(Ny,Nx,Nj)
    :param q_cj: jth current sediment volume concentration in all cells

    :return: The near bed concentration
    :rtype: numpy.ndarray(Ny,Nx,Nj)
    '''
    # Nj = 2
    # Ny=Nx=4
    # D_sj= np.ones((Nj))
    # D_sg = np.arange(1,Ny*Nx+1).reshape(Ny,Nx)
    # res = D_sj/D_sg[:,:,np.newaxis]
    # res[:,:,0]
    # res[:,:,1]

    with np.errstate(divide='ignore', invalid='ignore'):
        res = (0.4 * (D_sj / D_sg[:, :, np.newaxis]) ** (1.64) + 1.64) * q_cj
    #         print("(D_sj/D_sg[:,:,np.newaxis])**(1.64).shape",((D_sj/D_sg[:,:,np.newaxis])**(1.64)).shape)
    return res


def calc_averageSedimentSize(q_cj, D_sj):
    '''
    :param q_cj: jth sediment volume concentration
    :type q_cj: numpy.ndarray(Ny,Nx,Nj)

    :param D_sj: jth sediment diameter
    :type D_sj: numpy.ndarray(Nj)

    :rtype: numpy.ndarray(Ny,Nx)
    :return: Geometric mean size of suspended sediment mixture in cell
    '''
    return np.sum(q_cj * D_sj, axis=2)


def calc_depositionRate(v_sj, c_nbj):
    '''
    TODO! Equation (36)

    :type v_sj: numpy.ndarray(Nj)
    :param v_sj: j'th sediment fall velocity. [unit = m/s]

    :param c_nbj: The near bed concentration. [unit = 1]
    :type c_nbh: numpy.ndarray(Ny,Nx,Nj)

    :rtype: numpy.ndarray(Ny,Nx,Nj)
    :return: jth deposition rate for all cells. [unit = m/s]

    '''
    # This array * cube multiplication is tested and should work.
    return v_sj * c_nbj


def calc_dimless_sphere_settlingVel(v_sj, g_reduced, nu):
    '''
    This function calculates the dimensionless sphere\
    settling velocity using scaling factor:\
    (g' * nu)^(-1/3).\

    TODO! Equation (37), but its altered!

    :type v_sj: numpy.ndarray(Nj)
    :param v_sj: j'th sediment fall velocity
    :type g_reduced: numpy.ndarray(Nj)
    :param g_reduced: Scaling factors
    :type nu: float/double
    :param nu: Kinematic viscosity

    :return: Dimensionless sphere settling velocities
    :rtype: numpy.ndarray(Nj)

    '''
    return v_sj * np.cbrt(1 / (g_reduced * nu))


def calc_g_reduced(rho_j, rho_a, g=9.81):
    '''
    This function is used in calculating the scaling for\
    the dimensionless sphere settling velocity.

    :type rho_j: numpy.ndarray(Nj)
    :param rho_j: Density of sediment type no j. [kg/m^3]

    :type rho_a: double
    :param rho_a: Ambient fluid density [kg/m^3]

    :return: List of reduced gravities for sediment type j.
    :rtype: numpy.ndarray(Nj)

    TODO! Assumed rho = rho_ambient
    '''

    return g * (rho_j - rho_a) / rho_a

