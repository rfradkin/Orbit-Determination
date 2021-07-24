# Orbit Visualization, SSP New Mexico Tech 2021
# Name: Rom Fradkin

# The radius' of the planets are not to scale
# The orbital periods are accurate
# The relative positions of the planets are accurate

from f import *
from vpython import vector, sphere, curve, canvas, rate, color, textures
from astroquery.jplhorizons import Horizons
from datetime import datetime

print('Input yes or no depending on your response.')

use_current_jd = input('Would you like to use the current Julian date? ')
if use_current_jd.lower() == 'yes':
    now = datetime.utcnow()
    t = get_jd(now.year, now.month, now.day, now.hour, now.minute, now.second)
else:
    t = input('Enter Julian date: ')

all_planets_AST = input('Would you like to see all the planets and the asteroid in our solar system? ')
if all_planets_AST != 'yes':
    show_mercury = input("Would you like to show Mercury's orbit? ")
    show_venus = input("Would you like to show Venus's orbit? ")
    show_earth = input("Would you like to show Earth's orbit? ")
    show_mars = input("Would you like to show Mars's orbit? ")
    show_jupiter = input("Would you like to show Jupiter's orbit? ")
    show_saturn = input("Would you like to show Saturn's orbit? ")
    show_uranus = input("Would you like to show Uranus's orbit? ")
    show_neptune = input("Would you like to show Neptune's orbit? ")
    show_AST = input("Would you like to show the asteroid's orbit? ")
    if show_AST == 'yes':
        HA22_2003 = input("Would you like to see 2003 HA22's orbit? ")
        if HA22_2003 == 'yes':
            ast_num = '302871'
        else:
            ast_num = input("Input the asteroid's JPL identification number: ")
        use_JPL_orbital = input("Would you like to use JPL Horizons's orbital elements for your asteroid? ")
        if use_JPL_orbital != 'yes':
            print('\nRelevant inputs must be in degreees.\n')
            a_AST = float(input('Input the semimajor-axis (a): '))
            e_AST = float(input('Input the eccentricity (e): '))
            i_AST = radians(float(input('Input the inclination (i): ')))
            Ω_AST = radians(float(input('Input the longitude of the ascending node (Ω): ')))
            ω_AST = radians(float(input('Input the argument of the perihelion (ω): ')))
            M_AST = radians(float(input('Input the mean anomoly (M): ')))
else:
    HA22_2003 = input("Would you like to see 2003 HA22's orbit? ")
    if HA22_2003 == 'yes':
        ast_num = '302871'
    else:
        ast_num = input("Input the asteroid's JPL identification number: ")
    use_JPL_orbital = input("Would you like to use JPL Horizons's orbital elements for your asteroid? ")
    if use_JPL_orbital != 'yes':
        print('\nRelevant inputs must be in degreees.\n')
        a_AST = float(input('Input the semimajor-axis (a): '))
        e_AST = float(input('Input the eccentricity (e): '))
        i_AST = radians(float(input('Input the inclination (i): ')))
        Ω_AST = radians(float(input('Input the longitude of the ascending node (Ω): ')))
        ω_AST = radians(float(input('Input the argument of the perihelion (ω): ')))
        M_AST = radians(float(input('Input the mean anomoly (M): ')))
                
print('\nPress control and click and drag to move around.')
print('Use two fingers to zoom in and out.')
print('Please wait while your orbits are loaded.')

# Create planet colors
mercury_color = vector(0.835,0.824,0.820)
venus_color = vector(1.00,0.271,0.000)
mars_color = vector(0.737, 0.152, 0.196)
jupiter_color = vector(0.645, 0.612, 0.525)
saturn_color = vector(0.671, 0.365, 0.291)
uranus_color = vector(0.675, 0.898, 0.933)
neptune_color = vector(0.129, 0.137, 0.329)

planets_AST_radius = 15
position_stretch = 350

# Create the background for the simulation
background = canvas(title = 'Solar System Orbits',
     width = 1400, height = 500,
     center = vector(0,0,0), background=color.white)

if all_planets_AST.lower() == 'yes' or show_mercury.lower() == 'yes':
    sun_mercury_jpl = Horizons(id = '199', location = '10', epochs = t, id_type = 'majorbody')
    sun_mercury_jpl_elem = sun_mercury_jpl.elements()
    #Defining Mecury
    a_mercury = sun_mercury_jpl_elem['a'][0]
    e_mercury = sun_mercury_jpl_elem['e'][0]
    M_mercury = radians(sun_mercury_jpl_elem['M'][0])
    Ω_mercury = radians(sun_mercury_jpl_elem['Omega'][0])
    i_mercury = radians(sun_mercury_jpl_elem['incl'][0])
    ω_mercury = radians(sun_mercury_jpl_elem['w'][0])
    p_mercury = (4 * pi ** (2) * a_mercury ** (3)) ** (1/2)
    mercury = sphere(pos = vector(0,0,0), canvas = background, radius = (planets_AST_radius), color=mercury_color)
    mercury.trail = curve (color=mercury_color)
    M_t_mercury = 2 * pi / p_mercury * 0 + M_mercury

if all_planets_AST.lower() == 'yes' or show_venus.lower() == 'yes':
    sun_venus_jpl = Horizons(id = '299', location = '10', epochs = t, id_type = 'majorbody')
    sun_venus_jpl_elem = sun_venus_jpl.elements()
    #Defining Venus
    a_venus = sun_venus_jpl_elem['a'][0]
    e_venus = sun_venus_jpl_elem['e'][0]
    M_venus = radians(sun_venus_jpl_elem['M'][0])
    Ω_venus = radians(sun_venus_jpl_elem['Omega'][0])
    i_venus = radians(sun_venus_jpl_elem['incl'][0])
    ω_venus = radians(sun_venus_jpl_elem['w'][0])
    p_venus = (4 * pi ** (2) * a_venus ** (3)) ** (1/2)
    venus = sphere(pos = vector(0,0,0), canvas = background, radius = (planets_AST_radius), color = venus_color)
    venus.trail = curve (color=venus_color)
    M_t_venus = 2 * pi / p_venus * 0 + M_venus

if all_planets_AST.lower() == 'yes' or show_earth.lower() == 'yes':
    sun_earth_jpl = Horizons(id = '399', location = '10', epochs = t, id_type = 'majorbody')
    sun_earth_jpl_elem = sun_earth_jpl.elements()
    #Defining Earth
    a_earth = sun_earth_jpl_elem['a'][0]
    e_earth = sun_earth_jpl_elem['e'][0]
    M_earth = radians(sun_earth_jpl_elem['M'][0])
    Ω_earth = radians(sun_earth_jpl_elem['Omega'][0])
    i_earth = radians(sun_earth_jpl_elem['incl'][0])
    ω_earth = radians(sun_earth_jpl_elem['w'][0])
    p_earth = (4 * pi ** (2) * a_earth ** (3)) ** (1/2)
    earth = sphere(pos = vector(0,0,0), canvas = background, radius = (planets_AST_radius), texture = textures.earth)
    earth.trail = curve (color = color.blue)
    M_t_earth = 2 * pi / p_earth * 0 + M_earth

if all_planets_AST.lower() == 'yes' or show_mars.lower() == 'yes':
    sun_mars_jpl = Horizons(id = '499', location = '10', epochs = t, id_type = 'majorbody')
    sun_mars_jpl_elem = sun_mars_jpl.elements()
    #Defining Mars
    a_mars = sun_mars_jpl_elem['a'][0]
    e_mars = sun_mars_jpl_elem['e'][0]
    M_mars = radians(sun_mars_jpl_elem['M'][0])
    Ω_mars = radians(sun_mars_jpl_elem['Omega'][0])
    i_mars = radians(sun_mars_jpl_elem['incl'][0])
    ω_mars = radians(sun_mars_jpl_elem['w'][0])
    p_mars = (4 * pi ** (2) * a_mars ** (3)) ** (1/2)
    mars = sphere(pos = vector(0,0,0), canvas = background, radius = (planets_AST_radius), color = mars_color)
    mars.trail = curve (color = mars_color)
    M_t_mars = 2 * pi / p_mars * 0 + M_mars

if all_planets_AST.lower() == 'yes' or show_jupiter.lower() == 'yes':
    sun_jupiter_jpl = Horizons(id = '599', location = '10', epochs = t, id_type = 'majorbody')
    sun_jupiter_jpl_elem = sun_jupiter_jpl.elements()
    #Defining Jupiter
    a_jupiter = sun_jupiter_jpl_elem['a'][0]
    e_jupiter = sun_jupiter_jpl_elem['e'][0]
    M_jupiter = radians(sun_jupiter_jpl_elem['M'][0])
    Ω_jupiter = radians(sun_jupiter_jpl_elem['Omega'][0])
    i_jupiter = radians(sun_jupiter_jpl_elem['incl'][0])
    ω_jupiter = radians(sun_jupiter_jpl_elem['w'][0])
    p_jupiter = (4 * pi ** (2) * a_jupiter ** (3)) ** (1/2)
    jupiter = sphere(pos = vector(0,0,0), canvas = background, radius = (planets_AST_radius), color = jupiter_color)
    jupiter.trail = curve (color=jupiter_color)
    M_t_jupiter = 2 * pi / p_jupiter * 0 + M_jupiter

if all_planets_AST.lower() == 'yes' or show_saturn.lower() == 'yes':
    sun_saturn_jpl = Horizons(id = '699', location = '10', epochs = t, id_type = 'majorbody')
    sun_saturn_jpl_elem = sun_saturn_jpl.elements()
    #Defining Saturn
    a_saturn = sun_saturn_jpl_elem['a'][0]
    e_saturn = sun_saturn_jpl_elem['e'][0]
    M_saturn = radians(sun_saturn_jpl_elem['M'][0])
    Ω_saturn = radians(sun_saturn_jpl_elem['Omega'][0])
    i_saturn = radians(sun_saturn_jpl_elem['incl'][0])
    ω_saturn = radians(sun_saturn_jpl_elem['w'][0])
    p_saturn = (4 * pi ** (2) * a_saturn ** (3)) ** (1/2)
    saturn = sphere(pos = vector(0,0,0), canvas = background, radius = (planets_AST_radius), color = saturn_color)
    saturn.trail = curve (color=saturn_color)
    M_t_saturn = 2 * pi / p_saturn * 0 + M_saturn

if all_planets_AST.lower() == 'yes' or show_uranus.lower() == 'yes':
    sun_uranus_jpl = Horizons(id = '799', location = '10', epochs = t, id_type = 'majorbody')
    sun_uranus_jpl_elem = sun_uranus_jpl.elements()
    #Defining Uranus
    a_uranus = sun_uranus_jpl_elem['a'][0]
    e_uranus = sun_uranus_jpl_elem['e'][0]
    M_uranus = radians(sun_uranus_jpl_elem['M'][0])
    Ω_uranus = radians(sun_uranus_jpl_elem['Omega'][0])
    i_uranus = radians(sun_uranus_jpl_elem['incl'][0])
    ω_uranus = radians(sun_uranus_jpl_elem['w'][0])
    p_uranus = (4 * pi ** (2) * a_uranus ** (3)) ** (1/2)
    uranus = sphere(pos = vector(0,0,0), canvas = background, radius = (planets_AST_radius), color = uranus_color)
    uranus.trail = curve (color=uranus_color)
    M_t_uranus = 2 * pi / p_uranus * 0 + M_uranus

if all_planets_AST.lower() == 'yes' or show_neptune.lower() == 'yes':
    sun_neptune_jpl = Horizons(id = '899', location = '10', epochs = t, id_type = 'majorbody')
    sun_neptune_jpl_elem = sun_neptune_jpl.elements()
    #Defining Neptune
    a_neptune = sun_neptune_jpl_elem['a'][0]
    e_neptune = sun_neptune_jpl_elem['e'][0]
    M_neptune = radians(sun_neptune_jpl_elem['M'][0])
    Ω_neptune = radians(sun_neptune_jpl_elem['Omega'][0])
    i_neptune = radians(sun_neptune_jpl_elem['incl'][0])
    ω_neptune = radians(sun_neptune_jpl_elem['w'][0])
    p_neptune = (4 * pi ** (2) * a_neptune ** (3)) ** (1/2)
    neptune = sphere(pos = vector(0,0,0), canvas = background, radius = (planets_AST_radius), color = neptune_color)
    neptune.trail = curve(color=neptune_color)
    M_t_neptune = 2 * pi / p_neptune * 0 + M_neptune

if all_planets_AST.lower() == 'yes' or show_AST.lower() == 'yes':
    if use_JPL_orbital == 'yes':
        sun_AST_jpl = Horizons(id = ast_num, location = '10', epochs = t, id_type = 'smallbody')
        sun_AST_jpl_elem = sun_AST_jpl.elements()
        #Defining AST
        a_AST = sun_AST_jpl_elem['a'][0]
        e_AST = sun_AST_jpl_elem['e'][0]
        M_AST = radians(sun_AST_jpl_elem['M'][0])
        Ω_AST = radians(sun_AST_jpl_elem['Omega'][0])
        i_AST = radians(sun_AST_jpl_elem['incl'][0])
        ω_AST = radians(sun_AST_jpl_elem['w'][0])
    p_AST = (4 * pi ** (2) * a_AST ** (3)) ** (1/2)
    AST = sphere(pos = vector(0,0,0), canvas = background, color = color.black, radius = (planets_AST_radius), texture = textures.gravel, \
                make_trail = True, trail_type = 'points', interval = 2)
    M_t_AST = 2 * pi / p_AST * 0 + M_AST

sun = sphere(position = vector(0,0,0), radius=(40), color = color.yellow)

if all_planets_AST.lower() == 'yes' or show_mercury.lower() == 'yes':
    #Placing of Mercury and labeling
    E_mercury = get_E(M_t_mercury, e_mercury)
    υ_mercury = transform(E_mercury, a_mercury, e_mercury, i_mercury, ω_mercury, Ω_mercury)
    mercury.pos = υ_mercury * position_stretch

if all_planets_AST.lower() == 'yes' or show_venus.lower() == 'yes':
    #Placing of Venus and labeling
    E_venus = get_E(M_t_venus, e_venus)
    υ_venus = transform(E_venus, a_venus, e_venus, i_venus, ω_venus, Ω_venus)
    venus.pos = υ_venus * position_stretch

if all_planets_AST.lower() == 'yes' or show_earth.lower() == 'yes':
    #Placing of Earth and labeling
    E_earth = get_E(M_t_earth, e_earth)
    υ_earth = transform(E_earth, a_earth, e_earth, i_earth, ω_earth, Ω_earth)
    earth.pos = υ_earth * position_stretch

if all_planets_AST.lower() == 'yes' or show_mars.lower() == 'yes':
    #Placing of Mars and labeling
    E_mars = get_E(M_t_mars, e_mars)
    υ_mars = transform(E_mars, a_mars, e_mars, i_mars, ω_mars, Ω_mars)
    mars.pos = υ_mars * position_stretch

if all_planets_AST.lower() == 'yes' or show_jupiter.lower() == 'yes':
    #Placing of Jupiter and labeling
    E_jupiter = get_E(M_t_jupiter, e_jupiter)
    υ_jupiter = transform(E_jupiter, a_jupiter, e_jupiter, i_jupiter, ω_jupiter, Ω_jupiter)
    jupiter.pos = υ_jupiter * position_stretch

if all_planets_AST.lower() == 'yes' or show_saturn.lower() == 'yes':
    #Placing of Saturn and labeling
    E_saturn = get_E(M_t_saturn, e_saturn)
    υ_saturn = transform(E_saturn, a_saturn, e_saturn, i_saturn, ω_saturn, Ω_saturn)
    saturn.pos = υ_saturn * position_stretch

if all_planets_AST.lower() == 'yes' or show_uranus.lower() == 'yes':
    #Placing of Uranus and labeling
    E_uranus = get_E(M_t_uranus, e_uranus)
    υ_uranus = transform(E_uranus, a_uranus, e_uranus, i_uranus, ω_uranus, Ω_uranus)
    uranus.pos = υ_uranus * position_stretch

if all_planets_AST.lower() == 'yes' or show_neptune.lower() == 'yes':
    #Placing of Neptune and labeling
    E_neptune = get_E(M_t_neptune, e_neptune)
    υ_neptune = transform(E_neptune, a_neptune, e_neptune, i_neptune, ω_neptune, Ω_neptune)
    neptune.pos = υ_neptune * position_stretch

if all_planets_AST.lower() == 'yes' or show_AST.lower() == 'yes':
    #Placing of AST and labeling
    E_AST = get_E(M_t_AST, e_AST)
    υ_AST = transform(E_AST, a_AST, e_AST, i_AST, ω_AST, Ω_AST)
    AST.pos = υ_AST * position_stretch

time = 0
Δt = .065

while True:

    rate(25)
    time += Δt

    if all_planets_AST.lower() == 'yes' or show_mercury.lower() == 'yes':
        #Update mercury position
        M_t_mercury = 2 * pi / p_mercury * time + M_mercury
        E_mercury = get_E(M_t_mercury, e_mercury)
        υ_mercury = transform(E_mercury, a_mercury, e_mercury, i_mercury, ω_mercury, Ω_mercury)
        mercury.pos = υ_mercury * position_stretch
        mercury.trail.append(pos = mercury.pos)

    if all_planets_AST.lower() == 'yes' or show_venus.lower() == 'yes':
        #Update Venus position
        M_t_venus = 2 * pi / p_venus * time + M_venus
        E_venus = get_E(M_t_venus, e_venus)
        υ_venus = transform(E_venus, a_venus, e_venus, i_venus, ω_venus, Ω_venus)
        venus.pos = υ_venus * position_stretch
        venus.trail.append(pos = venus.pos)

    if all_planets_AST.lower() == 'yes' or show_earth.lower() == 'yes':
        #Update Earth position
        M_t_earth = 2 * pi / p_earth * time + M_earth
        E_earth = get_E(M_t_earth, e_earth)
        υ_earth = transform(E_earth, a_earth, e_earth, i_earth, ω_earth, Ω_earth)
        earth.pos = υ_earth * position_stretch
        earth.trail.append(pos = earth.pos)

    if all_planets_AST.lower() == 'yes' or show_mars.lower() == 'yes':
        #Update Mars position
        M_t_mars = 2 * pi / p_mars * time + M_mars
        E_mars = get_E(M_t_mars, e_mars)
        υ_mars = transform(E_mars, a_mars, e_mars, i_mars, ω_mars, Ω_mars)
        mars.pos = υ_mars * position_stretch
        mars.trail.append(pos = mars.pos)

    if all_planets_AST.lower() == 'yes' or show_jupiter.lower() == 'yes':
        #Update Jupiter position
        M_t_jupiter = 2 * pi / p_jupiter * time + M_jupiter
        E_jupiter = get_E(M_t_jupiter, e_jupiter)
        υ_jupiter = transform(E_jupiter, a_jupiter, e_jupiter, i_jupiter, ω_jupiter, Ω_jupiter)
        jupiter.pos = υ_jupiter * position_stretch
        jupiter.trail.append(pos = jupiter.pos)

    if all_planets_AST.lower() == 'yes' or show_saturn.lower() == 'yes':
        #Update Saturn position
        M_t_saturn = 2 * pi / p_saturn * time + M_saturn
        E_saturn = get_E(M_t_saturn, e_saturn)
        υ_saturn = transform(E_saturn, a_saturn, e_saturn, i_saturn, ω_saturn, Ω_saturn)
        saturn.pos = υ_saturn * position_stretch
        saturn.trail.append(pos = saturn.pos)

    if all_planets_AST.lower() == 'yes' or show_uranus.lower() == 'yes':
        #Update Uranus position
        M_t_uranus = 2 * pi / p_uranus * time + M_uranus
        E_uranus = get_E(M_t_uranus, e_uranus)
        υ_uranus = transform(E_uranus, a_uranus, e_uranus, i_uranus, ω_uranus, Ω_uranus)
        uranus.pos = υ_uranus * position_stretch
        uranus.trail.append(pos = uranus.pos)

    if all_planets_AST.lower() == 'yes' or show_neptune.lower() == 'yes':
        #Update Neptune position
        M_t_neptune = 2 * pi / p_neptune * time + M_neptune
        E_neptune = get_E(M_t_neptune, e_neptune)
        υ_neptune = transform(E_neptune, a_neptune, e_neptune, i_neptune, ω_neptune, Ω_neptune)
        neptune.pos = υ_neptune * position_stretch
        neptune.trail.append(pos = neptune.pos)

    if all_planets_AST.lower() == 'yes' or show_AST.lower() == 'yes':
        #Update AST position
        M_t_AST = 2 * pi / p_AST * time + M_AST
        E_AST = get_E(M_t_AST, e_AST)
        υ_AST = transform(E_AST, a_AST, e_AST, i_AST, ω_AST, Ω_AST)
        AST.pos = υ_AST * position_stretch
    