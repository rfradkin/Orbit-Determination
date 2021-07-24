# Method of Gauss, SSP New Mexico Tech 2021
# Name: Rom Fradkin

from f import *
from numpy import array

def meth_of_gaus(file_name, path, t0_greg = '2021 7 24 07:00:00.000'):

    # Data preparation
    data = load_file(file_name, path = path)
    
    if len(data) != 3:
        print(f'There are {len(data)} possible observations. You must choose 3. To choose an observation, type yes or no.')
        obse = []
        for k in range(len(data)):
            if input(f'Would you like to use observation {k + 1}? ').lower() == 'yes':
                obse.append(k)
            if len(obse) == 3:
                break
        if len(obse) != 3:
            raise RuntimeError('You must choose only 3 observations.')
    else:
        obse = [0, 1, 2]

    y = array([data[obse[0], 0], data[obse[1], 0], data[obse[2], 0]])
    m = array([data[obse[0], 1], data[obse[1], 1], data[obse[2], 1]])
    d = array([data[obse[0], 2], data[obse[1], 2], data[obse[2], 2]])
    hh = array([data[obse[0], 3], data[obse[1], 3], data[obse[2], 3]])
    mm = array([data[obse[0], 4], data[obse[1], 4], data[obse[2], 4]])
    ss = array([data[obse[0], 5], data[obse[1], 5], data[obse[2], 5]])
    α_hour_min_sec = array([data[obse[0], 6:9], data[obse[1], 6:9], data[obse[2], 6:9]])
    δ_deg_amin_asec = array([data[obse[0], 9:12], data[obse[1], 9:12], data[obse[2], 9:12]])      
    R = array([data[obse[0], 12:15], data[obse[1], 12:15], data[obse[2], 12:15]])

    # Convert gregorian times to Julian dates
    t = array([get_jd(y[0], m[0], d[0], hh[0], mm[0], ss[0]), get_jd(y[1], m[1], d[1], hh[1], mm[1], ss[1]), get_jd(y[2], m[2], d[2], hh[2], mm[2], ss[2])])

    # Convert right ascencion and declination to degrees
    α = array([hour_min_sec_to_radi(α_hour_min_sec[0, 0], α_hour_min_sec[0, 1], α_hour_min_sec[0, 2]), hour_min_sec_to_radi(α_hour_min_sec[1, 0], \
        α_hour_min_sec[1, 1], α_hour_min_sec[1, 2]), hour_min_sec_to_radi(α_hour_min_sec[2, 0], α_hour_min_sec[2, 1], α_hour_min_sec[2, 2])])
    δ = array([degr_amin_asec_to_radi(δ_deg_amin_asec[0, 0], δ_deg_amin_asec[0, 1], δ_deg_amin_asec[0, 2]), degr_amin_asec_to_radi(δ_deg_amin_asec[1, 0], \
        δ_deg_amin_asec[1, 1], δ_deg_amin_asec[1, 2]), degr_amin_asec_to_radi(δ_deg_amin_asec[2, 0], δ_deg_amin_asec[2, 1], δ_deg_amin_asec[2, 2])])

    # Calculate the rho hat vector
    ρ_hat = array([get_ρ_hat(α[0], δ[0]), get_ρ_hat(α[1], δ[1]), get_ρ_hat(α[2], δ[2])])

    # Convert to Gaussian time
    𝜏 = conv_gaus_time(t[0], t[1], t[2])

    # Calculate the Ds
    D0 = get_D0(ρ_hat[0], ρ_hat[1], ρ_hat[2])
    D1 = array([get_D1i(R[0], ρ_hat[1], ρ_hat[2]), get_D1i(R[1], ρ_hat[1], ρ_hat[2]), get_D1i(R[2], ρ_hat[1], ρ_hat[2])])
    D2 = array([get_D2i(R[0], ρ_hat[0], ρ_hat[2]), get_D2i(R[1], ρ_hat[0], ρ_hat[2]), get_D2i(R[2], ρ_hat[0], ρ_hat[2])])
    D3 = array([get_D3i(R[0], ρ_hat[0], ρ_hat[1]), get_D3i(R[1], ρ_hat[0], ρ_hat[1]), get_D3i(R[2], ρ_hat[0], ρ_hat[1])])

    # Solve the Scaler Equation of Lagrange
    r2_mag = solv_SEL(𝜏, R[1], ρ_hat[1], D0, D2)

    # Calculate the f and g series as second degree series, as there is no r2_dot
    f1, g1 = get_f_g(𝜏[0], r2_mag, None, flag = 'second')
    f3, g3 = get_f_g(𝜏[1], r2_mag, None, flag = 'second')

    # Calculate the Cs 
    c1 = get_c1(f1, f3 ,g1, g3)
    c2 = -1
    c3 = get_c3(f1, f3 ,g1, g3)

    # Calculate the rho magnitudes
    ρ_mag = get_ρ(c1, c2, c3, D0, D1, D2, D3)

    # Calculate the Sun to asteroid vector
    r = (ρ_hat * ρ_mag) - R

    # Calculate the Ds
    d1 = get_d1(f1, f3, g1, g3)
    d3 = get_d3(f1, f3, g1, g3)

    # Calculate the r2_dot vector
    r2_dot = get_r2_dot(d1, r[0], d3, r[2])

    # Calculate the time correction (due to the time it takes for light to travel), then calculate the Gaussian times
    corr_t = array([corr_time(t[0], ρ_mag[0]), corr_time(t[1], ρ_mag[1]), corr_time(t[2], ρ_mag[2])])
    corr_𝜏 = conv_gaus_time(corr_t[0], corr_t[1], corr_t[2])

    # Initialize ρ2 magnitudes to a difference greater than 1e-12
    prev_ρ2 = 0
    ρ2 = 1

    itera = 0
    while abs(ρ2 - prev_ρ2) > 1e-12:

        prev_ρ2 = ρ_mag[1]

        # Calculate the f and g series using the Newton Raphson function
        f1, g1 = get_f_g(corr_𝜏[0], r[1], r2_dot, flag = 'Newton_Raphson')
        f3, g3 = get_f_g(corr_𝜏[1], r[1], r2_dot, flag = 'Newton_Raphson')  

        # Calculate the Cs
        c1 = get_c1(f1, f3 ,g1, g3)
        c2 = -1
        c3 = get_c3(f1, f3 ,g1, g3)
        ρ_mag = get_ρ(c1, c2, c3, D0, D1, D2, D3)

        # Calculate the Sun to asteroid vector
        r = (ρ_hat * ρ_mag) - R

        # Calculate the ds
        d1 = get_d1(f1, f3, g1, g3)
        d3 = get_d3(f1, f3, g1, g3)

        # Calculate the r2_dot vector
        r2_dot = get_r2_dot(d1, r[0], d3, r[2])

        # Calculate the time correction (due to the time it takes for light to travel), then calculate the Gaussian times
        corr_t = array([corr_time(t[0], ρ_mag[0]), corr_time(t[1], ρ_mag[1]), corr_time(t[2], ρ_mag[2])])
        corr_𝜏 = conv_gaus_time(corr_t[0], corr_t[1], corr_t[2])

        ρ2 = ρ_mag[1]

        if itera > 20000:
            raise RuntimeWarning('ρ2 taking too long to converge.')

        itera += 1

    r2 = r[1]
    r2_dot = r2_dot.reshape((3,1))

    # Rotate from equitorial to ecpliptic
    r2_ecli = rota_equi_ecli(r2).squeeze()
    r2_dot_ecli= rota_equi_ecli(r2_dot).squeeze()

    t0_greg = t0_greg.replace(':', ' ').split()
    for i in range(len(t0_greg)):
        t0_greg[i] = float(t0_greg[i])

    # Calculate the epoch time in Julian days
    y, m, d, hh, mm, ss = t0_greg
    t0 = get_jd(y, m, d, hh, mm, ss)

    # Return the orbital elements a, e, i, Ω, ω, M in degrees
    return get_orbi_elem(r2_ecli, r2_dot_ecli, t[1], t0)

print(meth_of_gaus('2003HA22.txt', '/Users/romfradkin/Desktop/SSP/'))

# 2003 HA22 Input File
# Y    D  M  h  m  s      α (Right Ascension) δ (Declination) R (Earth to Sun Vector) x, y, z
#
# 2021 06 29 03:47:27.000 15:32:47.94360 -13:29:50.0748 -1.298966056621338E-01 9.251322110560322E-01 4.010446264180329E-01
# 2021 07 04 03:11:41.090 15:46:18.88320 -14:59:26.4876 -2.129418402389331E-01 9.121786007644945E-01 3.953850295560716E-01
# 2021 07 11 02:05:38.080 16:09:54.91440 -17:18:35.5140 -3.263164629094907E-01 8.834531201812120E-01 3.829838974162144E-01
# 2021 07 17 02:06:21.0 16:34:36.15240 -19:24:12.0888 -4.205136842986850E-01 8.489850403484872E-01 3.680419069814028E-01

