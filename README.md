# Orbital-Determination
This program arrives at orbital elements that describe and can be used to model the motion of an asteroid after image reduction and processing using least-squares plate reduction astrometry. 

This code was developed as a part of the Summer Science Program in Astrophysics at the New Mexico Institute of Technology during the summer of 2021. 

## Output

### JPL HORIZONS Based Orbit of 2003 HA22 Visual Orbit
<img src="/Images/JPL_Orbit.jpg" alt="Default Login Screen" width="450"/>  

### Observation Based Orbit of 2003 HA22 Visual Orbit
<img src="/Images/Observational_Orbit.jpg" alt="Default Login Screen" width="450"/> 

### 2003 HA22 on June 29, 2021
<img src="/Images/2003_HA22_SLM.jpg" alt="Default Login Screen" width="450"/> 


| Orbital Element        |   Output           |  JPL    |
|:-------------:|:-------------:|:--------:|
| *a*  | 11.6 au | 1.876 au|
| *e*     | 0.893      |   0.395 |
| *i* | 2.723°       |  1.608°|
| *&Omega;* | 117.894°  |  121.804°|
| *&omega;* | 116.217°   |  164.439°|
| *M* | 0.221° |  294.957°|

## Files:

### R.Fradkin_Orbit_Determination.py  &  f.py

This is the main function. It imports functions from f.py. 

How to run:
1) Open R.Fradkin_OD.py and f.py
2) replace "input file here.txt" with the input file in R.Fradkin_OD.py/2003HA22.txt
3) Save and run R.Fradkin_OD.py
4) Debug :)

### R.Fradkin_Monte_Carlo.py  &  f.py

This is a modified Method of Gauss code that works with Monte Carlo and uses basic functions from f.py.
Replace the "input file here.txt" with the input file in R.Fradkin_Monte_Carlo.py/2003HA22_Monte_Carlo.txt to run. The input file can only have 3 observations. Uncertainties should be in the final two columns. 

### R.Fradkin_Orbit_Visualization.py  &  f.py

2003 HA22's orbit visualized with the Solar System.

