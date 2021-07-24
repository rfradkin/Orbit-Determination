# Method of Gauss functions, SSP New Mexico Tech 2021
# Name: Rom Fradkin

from os import remove
from numpy import array, cross, linalg, loadtxt, real, isreal, polynomial
from vpython import vector, rotate
from math import degrees, radians, sin, cos, asin, acos, atan2, pi

global k, μ, c_au
k = 0.0172020989484
μ = 1
c_au = 173.144643267

def conv_gaus_time(t1, t2, t3):

    # Convert to Gaussian time
    𝜏1 = k * (t1 - t2)
    𝜏3 = k * (t3 - t2)
    𝜏 = 𝜏3 - 𝜏1
    return array([𝜏1, 𝜏3, 𝜏])

def degr_amin_asec_to_radi(dd, am, asec): 

    # Convert from degrees, arcminutes, and arcseconds to radians
    deg_min_sec_deg = abs(dd) + am / 60 + asec / 3600
    if dd < 0:
        deg_min_sec_deg *= -1
    return radians(deg_min_sec_deg) 

def corr_time(𝜏0i, ρi):

    # Calculate the time correction (due to the time it takes for light to travel)
    return 𝜏0i - ρi / c_au

def ephe(a, e, i, Ω, ω, M0, t, t0, R, ε = radians(23.4374)):

    # Generate the ephemeris (right ascenscion and declination)
    n = k / a ** (3 / 2)
    M = n * (t - t0) + M0
    E = get_E(M, e)

    X = a * (cos(E) - e)
    Y = a * (1 - e ** 2) ** (1 / 2) * sin(E)
    Z = 0

    rota_matr = array([[cos(Ω) * cos(ω) - cos(i) * sin(Ω) * sin(ω), -cos(i) * cos(ω) * sin(Ω) - cos(Ω) * sin(ω),sin(Ω) * sin(i)], 
                          [sin(Ω) * cos(ω) + cos(i) * cos(Ω) * sin(ω), cos(i) * cos(ω) * cos(Ω)- sin(Ω) * sin(ω), -cos(Ω) * sin(i)], 
                          [sin(i) * sin(ω), sin(i) * cos(ω), cos(i)]])

    (x1, y1, z1) = rota_matr @ array([X, Y, Z]).reshape((3, 1))

    r = array([x1, y1, z1])
    ρ = R + r
    ρ_mag = linalg.norm(ρ)
    ρ_hat = ρ / ρ_mag
    equi_ρ = array([ρ_hat[0], 
                       ρ_hat[1] * cos(ε) - ρ_hat[2] * sin(ε), 
                       ρ_hat[1] * sin(ε) + ρ_hat[2] * cos(ε)])

    δ = asin(equi_ρ[2])
    cos_α = equi_ρ[0] / cos(δ)
    sin_α = equi_ρ[1] / cos(δ)
    α = atan2(sin_α, cos_α)
    return (degrees(α), degrees(δ))

def get_D0(ρ_hat1, ρ_hat2, ρ_hat3):

    # Calculate D0
    return ρ_hat1 @ cross(ρ_hat2, ρ_hat3)

def get_D1i(Ri, ρ_hat2, ρ_hat3):

    # Calculate D1
    return cross(Ri, ρ_hat2) @ ρ_hat3

def get_D2i(Ri, ρ_hat1, ρ_hat3):

    # Calculate D2
    return cross(ρ_hat1, Ri) @ ρ_hat3

def get_D3i(Ri, ρ_hat1, ρ_hat2):

    # Calculate D3
    return ρ_hat1 @ cross(ρ_hat2, Ri)

def get_E(M, e):

    # Calculate the eccentric anomoly
    E_guess = M
    M_guess = E_guess - e *  sin(E_guess)
    M_true = M
    while abs(M_guess - M_true) > 1e-4:
        M_guess = E_guess - e * sin(E_guess)
        E_guess = E_guess - (M_true - (E_guess - e * sin(E_guess))) / (e * cos(E_guess) - 1)
    return E_guess

def get_c1(f1, f3, g1, g3):

    # Calculate c1
    return g3 / (f1 * g3 - f3 * g1)

def get_c3(f1, f3, g1, g3):

    # Calculate c3
    return -g1 / (f1 * g3 - f3 * g1)

def get_d1(f1, f3, g1, g3):

    # Calculate d1
    return -f3 / (f1 * g3 - f3 * g1)

def get_d3(f1, f3, g1, g3):

    # Calculate d3
    return f1 / (f1 * g3 - f3 * g1)

def get_f_g(𝜏i, r2, r2_dot, flag):

    # Calculate the f and g series based using the method specified in the flag
    r2_mag = linalg.norm(r2)
    u = μ / (r2_mag ** 3)

    if r2_dot is not None:
        a = 1 / ((2 / r2_mag) - ((r2_dot @ r2_dot) / μ))
        e = (1 - (linalg.norm(cross(r2, r2_dot)) ** 2) / ((μ * a)) ** (1 / 2))
        n = ((μ / a ** 3) ** (1 / 2))
        z = r2 @ r2_dot / (r2_mag ** 2)
        q = ((r2_dot @ r2_dot) / (r2_mag ** 2)) - u

    if flag == 1 or type(flag) == str and flag.lower() == 'second':
        f = 1 - (𝜏i ** 2) / 2 * u
        g = 𝜏i - (𝜏i ** 3) / 6 * u
    elif flag == 2 or type(flag) == str and flag.lower() == 'third':
        f = 1 - (𝜏i ** 2) / 2 * u + (𝜏i ** 3) / 2 * u * z
        g = 𝜏i - (𝜏i ** 3) / 6 * u
    elif flag == 3 or type(flag) == str and flag.lower() == 'fourth':
        f = 1 - (𝜏i ** 2) / 2 * u + (𝜏i ** 3) / 2 * u * z + (𝜏i ** 4) / 24 * (3 * u * q - 15 * u * (z ** 2) + (u ** 2))
        g = 𝜏i - (𝜏i ** 3) / 6 * u + (𝜏i ** 4) / 4 * u * z
    elif flag == 4 or type(flag) == str and flag.lower() == 'newton_raphson':
        term = (r2 @ r2_dot) / (n * (a ** 2))
        if e < 0.1:
            sign = term * cos(n * 𝜏i - term) + (1 - r2_mag/ a) * sin(n * 𝜏i - term) > 0
            if sign:
                E1 = n * 𝜏i + 0.85 * e - term
            else:
                E1 = n * 𝜏i - 0.85 * e - term
        else:
            E1 = n * 𝜏i
        E0 = E1 + 1
        while abs(E1 - E0) > 1e-12:
            E0 = E1
            E1 = E0 - (E0- (1 - r2_mag / a) * sin(E0) + term * (1 - cos(E0)) - n * 𝜏i) / (1 - (1 - r2_mag / a) * cos(E0) + term * sin(E0))
        ΔE = E1
        f = 1 - (a / r2_mag) *  (1 - cos(ΔE))
        g = 𝜏i + (1 / n) * (sin(ΔE) - ΔE)
    else:
        raise ValueError('Flag not set properly.')

    return f, g

def get_jd(y, m, d, hh, mm, ss): 

    # Calculate the Julian date
    j_o = 367 * y - int(7 * (y + int((m + 9) / 12)) / 4) + int(275 * m / 9) + d + 1721013.5
    jd = j_o + (hh + mm / 60 + ss / 3600) / 24
    return jd

def get_orbi_elem(r, r_dot, t, t0):

    # Calculate the orbital elements
    r_mag = linalg.norm(r)
    r_cros_r_dot = h = cross(r, r_dot)
    mag_r_cros_r_dot = linalg.norm(r_cros_r_dot)
    
    a = 1 / (2 / r_mag - r_dot @ r_dot)
    e = (1 - (mag_r_cros_r_dot ** 2) / a) ** (1 / 2)
    i = acos(r_cros_r_dot[2] / mag_r_cros_r_dot)

    sin_Ω = r_cros_r_dot[0] / (mag_r_cros_r_dot * sin(i))
    cos_Ω = -r_cros_r_dot[1] / (mag_r_cros_r_dot * sin(i))

    Ω = atan2(sin_Ω, cos_Ω)
    
    sin_fw = r[2] / (r_mag * sin(i))
    cos_fw = 1 / (cos(Ω)) * (r[0] / r_mag + cos(i) * sin_fw * sin_Ω)

    fw = atan2(sin_fw, cos_fw)

    sin_f = (r @ r_dot) / (e * r_mag) * (a * (1 - e ** 2)) **  (1 / 2)
    cos_f = 1 / e * (a * (1 - e ** 2) / (r_mag) - 1)

    f = atan2(sin_f, cos_f)

    ω = fw - f

    n = k * (μ / (a ** 3)) ** (1 / 2)
    E2 = acos(1 / e * (1 - r_mag / a))

    f %= 2 * pi
    if pi < f < 2 * pi:
        E2 *= -1

    M2 = E2 - e * sin(E2)
    M0 = M2 + n * (t0 - t)

    Ω = degrees(Ω)
    i = degrees(i)
    ω = degrees(ω)
    M0 = degrees(M0)

    Ω %= 360
    i %= 360
    ω %= 360
    M0 %= 360

    return a, e, i, Ω, ω, M0

def get_r2_dot(d1, r1, d3, r3):

    # Calculate the r2_dot vector
    return  d1 * r1 + d3 * r3

def get_ρ(c1, c2, c3, D0, D1, D2, D3):

    # Calculate the rho vector
    ρ1 = (c1 * D1[0] + c2 * D1[1] + c3 * D1[2]) / (c1 * D0)
    ρ2 = (c1 * D2[0] + c2 * D2[1] + c3 * D2[2]) / (c2 * D0)
    ρ3 = (c1 * D3[0] + c2 * D3[1] + c3 * D3[2]) / (c3 * D0)
    return array([ρ1, ρ2, ρ3]).reshape((3,1))

def get_ρ_hat(ra, dec):

    # Calculate the rho hat vector
    return array([cos(ra) * cos(dec), sin(ra) * cos(dec), sin(dec)])

def hour_min_sec_to_radi(hh, mm, ss):

    # Convert from hours, minutes, and seconds to radians
    hour_min_sec_degr = hh * 15 + mm / 4 + ss / 240
    return radians(hour_min_sec_degr)

def load_file(file_name, path = ''):

    # Load the file
    with open(f'{path}{file_name}', 'r+') as data:
        _file = data.read()
        _file = _file.replace(':', ' ')

    upda_file = open(f'{path}updated_{file_name}', 'w')
    upda_file.write(_file)
    upda_file.close()

    file = loadtxt(f'{path}updated_{file_name}')
    remove(f'{path}updated_{file_name}')

    return file

def rota_equi_ecli(equi_vec, ε = radians(23.4374)):

    # Rotate from equitorial coordinates to ecliptic coordiantes
    rota_matr = array([ [1, 0, 0],
                        [0, cos(ε), sin(ε)],
                        [0, -sin(ε), cos(ε)]])
    return rota_matr @ equi_vec

def solv_SEL(𝜏, R2, ρ_hat2, D0, D2, mand_vali = 0): 

    # Solve the Scaler Equation of Lagrange
    A1 = 𝜏[1] / 𝜏[2]
    B1 = A1 / 6 * ((𝜏[2] ** 2) - (𝜏[1] ** 2))
    A3 = -𝜏[0] / 𝜏[2]
    B3 = A3 / 6 * ((𝜏[2] ** 2) - (𝜏[0] ** 2))

    A = (A1 * D2[0] - D2[1] + A3 * D2[2]) / -D0
    B = (B1 * D2[0] + B3 * D2[2]) / -D0
    E = -2 * ρ_hat2 @ R2
    F = linalg.norm(R2) ** 2

    a = -((A ** 2) + A * E + F)
    b = -μ * (2 * A * B + B * E)
    c = -(μ ** 2) * (B ** 2)

    roots = polynomial.polynomial.polyroots([c, 0, 0, b, 0, 0, a, 0, 1])

    vali_root = []
    for root in roots:
        if isreal(root):
            if root > 0:
                vali_root.append(real(root))   
    
    poss_ρ = []
    for root in vali_root:
        poss_ρ.append(A + μ * B / (root ** 3))

    vali_r2_mag = []
    for i in range(len(poss_ρ)):
        if poss_ρ[i] > 0 and vali_root[i] > 0:
            vali_r2_mag.append([vali_root[i], poss_ρ[i]])

    if mand_vali:
        return vali_r2_mag[mand_vali - 1][0]

    # If there are multiple possible r2 magnitudes, ask the user to choose 1
    if len(vali_r2_mag) > 1:
        print('Valid Roots (r2 magnitudes), Valid ρ')
        for i in range(len(vali_r2_mag)):
            print(f'{i + 1}: {vali_r2_mag[i]}')
        val_pair = int(input('There are multiple valid roots (r2 magnitudes) when solving the Scalar Equation of Lagrange. Which value pair should be used? ')) - 1
        return vali_r2_mag[val_pair][0]

    return vali_r2_mag[0][0]

def transform(E, a, e, i, ω, Ω):    

    # Transforms the position of the orbit    
    cartesian = vector(a * cos(E) - a * e, a * (1 - e ** 2) **  (1 / 2) * sin(E), 0)
    v1 = rotate(cartesian, angle = -ω, axis = vector(0, 0, 1))
    v2 = rotate(v1, angle = i, axis = vector(1, 0, 0))
    v3 = rotate(v2, angle = Ω, axis = vector(0, 0, 1))
    return v3




    
