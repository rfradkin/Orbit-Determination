# Monte Carlo Simulation for the Method of Gauss, SSP New Mexico Tech 2021
# Name: Rom Fradkin

from f import *
from numpy import array, vectorize, degrees, std
from numpy.random import normal
from matplotlib.pyplot import hist, xlabel, ylabel, title, show, axvline

def meth_of_gaus_mont_carl(α, δ, R, t, t0_greg = '2021 7 24 07:00:00.000'):

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
    r2_mag = solv_SEL(𝜏, R[1], ρ_hat[1], D0, D2, 1)

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

        if itera > 500:
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


def mont_carl(file_name, path, n = 2000):

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
    α_unce = array([data[obse[0], 15], data[obse[1], 15], data[obse[2], 15]])
    δ_unce = array([data[obse[0], 16], data[obse[1], 16], data[obse[2], 16]])  
    R = array([data[obse[0], 12:15], data[obse[1], 12:15], data[obse[2], 12:15]])

    # Convert gregorian times to Julian dates
    t = array([get_jd(y[0], m[0], d[0], hh[0], mm[0], ss[0]), get_jd(y[1], m[1], d[1], hh[1], mm[1], ss[1]), get_jd(y[2], m[2], d[2], hh[2], mm[2], ss[2])])

    # Convert right ascencion and declination to degrees
    α = array([hour_min_sec_to_radi(α_hour_min_sec[0, 0], α_hour_min_sec[0, 1], α_hour_min_sec[0, 2]), hour_min_sec_to_radi(α_hour_min_sec[1, 0], \
        α_hour_min_sec[1, 1], α_hour_min_sec[1, 2]), hour_min_sec_to_radi(α_hour_min_sec[2, 0], α_hour_min_sec[2, 1], α_hour_min_sec[2, 2])])
    δ = array([degr_amin_asec_to_radi(δ_deg_amin_asec[0, 0], δ_deg_amin_asec[0, 1], δ_deg_amin_asec[0, 2]), degr_amin_asec_to_radi(δ_deg_amin_asec[1, 0], \
        δ_deg_amin_asec[1, 1], δ_deg_amin_asec[1, 2]), degr_amin_asec_to_radi(δ_deg_amin_asec[2, 0], δ_deg_amin_asec[2, 1], δ_deg_amin_asec[2, 2])])

    # Create a normal distribution for each right ascencion and declination value
    ra_var1 = normal(α[0], α_unce[0], n)
    dec_var1 = normal(δ[0], δ_unce[0], n)
    ra_var2 = normal(α[1], α_unce[1], n)
    dec_var2 = normal(δ[1], δ_unce[1], n)
    ra_var3 = normal(α[2], α_unce[2], n)
    dec_var3 = normal(δ[2], δ_unce[2], n)

    # Create an array containing the varied right ascencion and declination values
    α_var = array([ra_var1, ra_var2, ra_var3])
    δ_var = array([dec_var1, dec_var2, dec_var3])

    # Create orbital element arrays containing the varied values
    a = []
    e = []
    i = []
    Ω = []
    ω = []
    M = []
    for j in range(n):
        try:
            elem = meth_of_gaus_mont_carl(α_var[:, j], δ_var[:, j], R, t)
            a.append(elem[0]) 
            e.append(elem[1]) 
            i.append(elem[2]) 
            Ω.append(elem[3]) 
            ω.append(elem[4]) 
            M.append(elem[5]) 
        except RuntimeWarning:
            pass

    # Convert pertinent values to degrees and reject the outliers
    a = array(a)
    e = array(e)
    i = degrees(array(i)) % 360
    Ω = degrees(array(Ω)) % 360
    ω = degrees(array(ω)) % 360
    M = degrees(array(M)) % 360

    # Create histograms for the varied orbital elements with a line showing the JPL Horizons predicted value
    hist(a, color = 'blue', bins = 80, histtype = 'step')
    hist(a, alpha = 0.3, color = 'blue', bins = 80)
    axvline(x = 1.876727658293008E+00, color = 'blue', linestyle = 'dotted')
    title('Sample Distribution of the Semi-Major Axis', fontname = 'Times New Roman', fontsize = 14)
    ylabel('Count', fontname = 'Times New Roman', fontsize = 14)
    xlabel('Semi-Major Axis (a) [AU]', fontname = 'Times New Roman', fontsize = 14)
    show()

    hist(e, color = 'green', bins = 80, histtype = 'step')
    hist(e, alpha = 0.3, color = 'green', bins = 80)
    axvline(x = 3.950197964098633E-01, color = 'green', linestyle = 'dotted')
    title('Sample Distribution of the Eccentricity', fontname = 'Times New Roman', fontsize = 14)
    ylabel('Count', fontname = 'Times New Roman', fontsize = 14)
    xlabel('Eccentricity (e)', fontname = 'Times New Roman', fontsize = 14)
    show()

    hist(i, color = 'red', bins = 80, histtype = 'step')
    hist(i, alpha = 0.3, color = 'red', bins = 80)
    axvline(x = 1.608955779163637E+00, color = 'red', linestyle = 'dotted')
    title('Sample Distribution of the Inclination', fontname = 'Times New Roman', fontsize = 14)
    ylabel('Count', fontname = 'Times New Roman', fontsize = 14)
    xlabel('Inclination (i)', fontname = 'Times New Roman', fontsize = 14)
    show()

    hist(Ω, color = 'turquoise', bins = 80, histtype = 'step')
    hist(Ω, alpha = 0.3, color = 'turquoise', bins = 80)
    axvline(x = 1.217877429178791E+02, color = 'turquoise', linestyle = 'dotted')
    title('Sample Distribution of the Longitude of the Ascending Node', fontname = 'Times New Roman', fontsize = 14)
    ylabel('Count', fontname = 'Times New Roman', fontsize = 14)
    xlabel('Longitude of the Ascending Node (Ω)', fontname = 'Times New Roman', fontsize = 14)
    show()

    hist(ω, color = 'magenta', bins = 80, histtype = 'step')
    hist(ω, alpha = 0.3, color = 'magenta', bins = 80)
    axvline(x = 1.635175564953411E+02, color = 'magenta', linestyle = 'dotted')
    title('Sample Distribution of the Longitude of the Argument of the Perihelion', fontname = 'Times New Roman', fontsize = 14)
    ylabel('Count', fontname = 'Times New Roman', fontsize = 14)
    xlabel('Argument of the Perihelion (ω)', fontname = 'Times New Roman', fontsize = 14)
    show()

    hist(M, color = 'orange', bins = 80, histtype = 'step')
    hist(M, alpha = 0.3, color = 'orange', bins = 80)
    axvline(x = 4.143177092212641E+00, color = 'orange', linestyle = 'dotted')
    title('Sample Distribution of the Mean Anomaly', fontname = 'Times New Roman', fontsize = 14)
    ylabel('Count', fontname = 'Times New Roman', fontsize = 14)
    xlabel('Argument of the Mean Anomaly (M)', fontname = 'Times New Roman', fontsize = 14)
    show()

    return std(a), std(e), std(i), std(Ω), std(ω), std(M)

print(mont_carl('2003HA22_Monte_Carlo.txt', '/Users/romfradkin/Desktop/SSP/', 2000))

# 2003 HA22 Monte Carlo Input File
# Y    D  M  h  m  s      α (Right Ascension) δ (Declination) R (Earth to Sun Vector) x, y, z, uncertainty α, uncertainty δ
#
# 2021 06 29 03:47:27.000 15:32:47.94360 -13:29:50.0748 -1.298966056621338E-01 9.251322110560322E-01 4.010446264180329E-01 2.0587612823419477e-06 1.9559651374123632e-06
# 2021 07 04 03:11:41.090 15:46:18.88320 -14:59:26.4876 -2.129418402389331E-01 9.121786007644945E-01 3.953850295560716E-01 1.4959616692067357e-06 3.1969512637175512e-06
# 2021 07 11 02:05:38.080 16:09:54.91440 -17:18:35.5140 -3.263164629094907E-01 8.834531201812120E-01 3.829838974162144E-01 8.488289466868082e-06 4.604510959183017e-06
# 2021 07 17 02:06:21.000 16:34:36.15240 -19:24:12.0888 -4.205136842986850E-01 8.489850403484872E-01 3.680419069814028E-01 1.6370097455364885e-07 3.011861142709904e-06