from fugacity import FugacityCalculator 
import math
import numpy as np



#Constants from paleopiezometry equations of Stipp and Tullis (2003) with the Holyoke (2010) correction.

CONSTANT_B = 2451
EXPONENT = -1.26

#Create dictionary of flow law coefficients
#flow_laws = collections.OrderedDict()
FLOW_LAWS = {
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

def calculate_differential_stress(grain_size):
	differential_stress = []
	for grain in grain_size:
		part = (math.log(grain)-math.log(CONSTANT_B))/EXPONENT
		part2 = math.exp(part)
		differential_stress.append(part2)
	return differential_stress

def calculate_strain_rate(differential_stress, fugacity_calculations, flow_law=None): 

    strain_rate = []

    for temperature, pressure, fugacity in fugacity_calculations:
        for stress in differential_stress:
            sr = (FLOW_LAWS[flow_law]['A']*np.power(stress, FLOW_LAWS[flow_law]['n'])*np.power(fugacity,1)*np.exp(-FLOW_LAWS[flow_law]['Q']/(8.3144598*temperature)))
            srn = np.float64(sr).item()
            strain_rate.append(srn)
    return strain_rate

def calculate_slip_rate(strain_rate, width): #width in m, output of mm/yr
    
    velocities = []
    
    for w in width:   
        for strain in strain_rate:
            vel = w*1000*31536000*strain
            velocities.append(vel)
    return velocities






