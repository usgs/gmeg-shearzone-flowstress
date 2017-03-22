# This code calculates strain rates, slip rates at a range of different temperatures incorporating changes 
# in fugacity as a function of temperature.

from __future__ import division # Import this to divide floats in python 2
import numpy as np 
import math

from scipy import optimize as opt 


# Calculate differential stress
d = 9.2 # measured grain size in microns



# Coefficients for calculating differential stress.
B = 2451 # 2451 (Holyoke 2010) #3631 (Stipp and Tullis, 2003)
m = -1.26 # Stress exponent, from Stipp and Tullis (2003).




#Create dictionary of flow law coefficients from literature.
"""
KT84: Kronenberg and Tullis, 1984
GT95wm: Gleason and Tullis, 1995 with melt
J84: Jaoul et al., 1984
K89: Koch, 1989
HC82: Hansen and Carter, 1982
LP92g: Luan and Patterson, 1992 with gel
LP92a: Luan and Patterson, 1992 with acid
H01: Hirth et al., 2001
RB04: Rutter and Brodie, 2004
"""
flow_laws = {
    "KT84": {"A": 2.2E-6, "n": 2.7, "Q": 1.2E5},
    'GT95wm': {"A": 1.8E-8, "n": 4, "Q": 1.37E5},
    'J84': {"A": 2.88E-3, "n": 1.8, "Q": 1.51E5},
    'K89': {"A": 1.1E-6, "n": 2.7, "Q": 1.34E5},
    'HC82': {"A": 1.99E-2, "n": 1.8, "Q": 1.67E5},
    'LP92g': {"A": 6.6E-8, "n": 3.1, "Q": 1.35E5},
    'LP92a': {"A": 3.98E-10, "n": 4, "Q": 1.35E5},
    "H01":  {"A": 6.3E-12, "n": 4, "Q": 1.35E5}, 
    "RB04": {"A": 1.2E-5, "n": 2.97, "Q": 2.42E5}
    }


#Define some basic functions

def temp_conv(T, direction = "C2K"):  #Convert from C to Kelvin or vice versa
    if direction == "C2K":
        Tk = T+273.15
        return Tk
    elif direction == "K2C":
        Tc = T-273.15
        return Tc
    else:
        print('ERROR: Not a valid direction of conversion')

def pres_conv(Pmpa): #Convert from MPa to Pa
    Ppa = round(Pmpa*1.0E6)
    return Ppa


#Defining PT Condition calculators and depth simulator, 
#alternatively use a range of temperatures at a given depth and pressure

def pressure_calculator(depth, density): #depth in km, density in g/cm3
    g = 9.8 #gravity
    h = depth*1000 #converting to SI
    row = density*1000 #converting to SI
    pressure = (row*g*h)/1.0E6 # pressure in MPa
    return pressure 

def temperature_calculator(depth, geothermal_gradient): #km, C/km
    temp = geothermal_gradient*depth 
    return temp
    
def depth_simulator(depth_range, density=2.7, geothermal_gradient=30): #depth range in km
    depth_conds = []
    for dep in range(depth_range[0], depth_range[1]):
        ptc = (round(pressure_calculator(dep, density)), temperature_calculator(geothermal_gradient, dep)) 
        depth_conds.append(ptc)
    return depth_conds


#Create the shear zone class

class Shearzone():
    def __init__(self, grain_size, constant_b, exponent):
        self.constant_b = constant_b 
        self.grain_size = grain_size 
        self.exponent = exponent
    
    def calculate_differential_stress(self): 
        part = (math.log(self.grain_size)-math.log(self.constant_b))/self.exponent
        s = math.exp(part)
        return s

#Create grain size simulator
def simulate_grain_size(range_min, range_max):
    B = 2451 #2451 (holyoke 2010) #3631 (stipp 2003)
    m = -1.26 # 
    diffs = []
    for x in range(range_min,range_max):
        part = (math.log(x)-math.log(B))/m
        s = math.exp(part)
        diffs.append(s)
    return diffs 




#Calculate fugacities.

#This section was adapted from PS.js by Tony Withers, see http://publish.uwo.ca/~awither5/fugacity/index.htm 
#for online calculator.

# Create a blank array
PScoeff = np.zeros([10,10])

# Input coefficients into array
PScoeff[0][2]=0.24657688*math.pow(10,6)
PScoeff[0][3]=0.51359951*math.pow(10,2) 
PScoeff[1][2]=0.58638965*math.pow(10,0) 
PScoeff[1][3]=-0.28646939*math.pow(10,-2) 
PScoeff[1][4]=0.31375577*math.pow(10,-4) 
PScoeff[2][2]=-0.62783840*math.pow(10,1) 
PScoeff[2][3]=0.14791599*math.pow(10,-1) 
PScoeff[2][4]=0.35779579*math.pow(10,-3) 
PScoeff[2][5]=0.15432925*math.pow(10,-7) 
PScoeff[3][3]=-0.42719875*math.pow(10,0) 
PScoeff[3][4]=-0.16325155*math.pow(10,-4) 
PScoeff[4][2]=0.56654978*math.pow(10,4) 
PScoeff[4][3]=-0.16580167*math.pow(10,2) 
PScoeff[4][4]=0.76560762*math.pow(10,-1) 
PScoeff[5][3]=0.10917883*math.pow(10,0) 
PScoeff[6][0]=0.38878656*math.pow(10,13) 
PScoeff[6][1]=-0.13494878*math.pow(10,9) 
PScoeff[6][2]=0.30916564*math.pow(10,6)
PScoeff[6][3]=0.75591105*math.pow(10,1) 
PScoeff[7][2]=-0.65537898*math.pow(10,5) 
PScoeff[7][3]=0.18810675*math.pow(10,3) 
PScoeff[8][0]=-0.14182435*math.pow(10,14) 
PScoeff[8][1]=0.18165390*math.pow(10,9) 
PScoeff[8][2]=-0.19769068*math.pow(10,6)
PScoeff[8][3]=-0.23530318*math.pow(10,2)
PScoeff[9][2]=0.92093375*math.pow(10,5)
PScoeff[9][3]=0.12246777*math.pow(10,3)


#Create function to make coefficient table for each T condition.
#Input coefficients in array from Eq 4 of Pitzer and Sterner (1994).

cs = np.zeros([10]) # Create blank array
def coeff_table():
    for i in range(0, 10):
        cs[i]=PScoeff[i][0]*math.pow(temperature,-4)+PScoeff[i][1]*math.pow(temperature,-2)\
        +PScoeff[i][2]*math.pow(temperature,-1)\
        +PScoeff[i][3]+PScoeff[i][4]*temperature+PScoeff[i][5]*math.pow(temperature,2)
        
#Call this function    
#coeff_table()


#Calculating equation of state and fugacity

#Solve Equation of state, Eq 2 of Pitzer and Sterner (1994)
#Returns pressure in Pa
def eos(T, V):
    den = 1/V
    R = 8314472
    var_num = cs[2]+2*cs[3]*den+3*cs[4]*math.pow(den,2)+4*cs[5]*math.pow(den,3)
    var_denom = math.pow((cs[1]+cs[2]*den+cs[3]*math.pow(den,2)+cs[4]*math.pow(den,3)+cs[5]*math.pow(den,4)),2)
    pressure=den+cs[0]*math.pow(den,2)-math.pow(den,2)*(var_num/var_denom)
    pressure= pressure + (cs[6]*math.pow(den,2)*math.exp(-cs[7]*den)+cs[8]*math.pow(den,2)*math.exp(-cs[9]*den))
    pressure = pressure*(R*T) #pressure in Pa
    return pressure

#Solve for fugacity, Eq 1 of Pitzer and Sterner (1994)
#Returns fugacity in MPa
def PSfug(P,T,V):
    den=1/V;
    R=8314472;
    quotient = cs[0]*den+(1/(cs[1]+cs[2]*den+cs[3]*math.pow(den,2)+cs[4]*math.pow(den,3)+cs[5]*math.pow(den,4))-1/cs[1])
    quotient-= cs[6]/cs[7]*(math.exp(-cs[7]*den)-1)
    quotient-= cs[8]/cs[9]*(math.exp(-cs[9]*den)-1)
    lnf=(math.log(den)+ quotient+P/(den*R*T))+math.log(R*T)-1
    return math.exp(lnf)/1e6 # fugacity in MPa




def fugacity_calculator():
    def fun(v):
        return eos(temperature, v)- pressure
    volume = opt.brentq(fun, 5, 30) #Volume in cc/mol

    #Calculate fugacity 
    fugacity = PSfug(pressure, temperature, volume)
    print(fugacity)
    
    return fugacity


#Calculating strain rate and slip rate, default is the Hirth et al (2001) flow law.

def calculate_strain_rate(s, flow_law = 'H01'): #s is differential stress, defaults to using Hirth et al flow law
    e = (flow_laws[flow_law]['A']*np.power(s,flow_laws[flow_law]['n'])*np.power(fugacity,1)*np.exp(-flow_laws[flow_law]['Q']/(8.3144598*temperature)))
    return e

def calculate_slip_rate(e, w): #w: width in m, output of mm/yr
    width = w*1000
    v = width*31536000*e
    return v







