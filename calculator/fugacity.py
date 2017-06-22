import numpy as np 
from scipy.constants.constants import C2K
import math
from scipy import optimize as opt # for optimization

import pdb


#Create coefficient tables to solve for equations of state

PS_COEFF = np.zeros([10,10])
PS_COEFF[0][2]=0.24657688*math.pow(10,6)
PS_COEFF[0][3]=0.51359951*math.pow(10,2) 
PS_COEFF[1][2]=0.58638965*math.pow(10,0) 
PS_COEFF[1][3]=-0.28646939*math.pow(10,-2) 
PS_COEFF[1][4]=0.31375577*math.pow(10,-4) 
PS_COEFF[2][2]=-0.62783840*math.pow(10,1) 
PS_COEFF[2][3]=0.14791599*math.pow(10,-1) 
PS_COEFF[2][4]=0.35779579*math.pow(10,-3) 
PS_COEFF[2][5]=0.15432925*math.pow(10,-7) 
PS_COEFF[3][3]=-0.42719875*math.pow(10,0) 
PS_COEFF[3][4]=-0.16325155*math.pow(10,-4) 
PS_COEFF[4][2]=0.56654978*math.pow(10,4) 
PS_COEFF[4][3]=-0.16580167*math.pow(10,2) 
PS_COEFF[4][4]=0.76560762*math.pow(10,-1) 
PS_COEFF[5][3]=0.10917883*math.pow(10,0) 
PS_COEFF[6][0]=0.38878656*math.pow(10,13) 
PS_COEFF[6][1]=-0.13494878*math.pow(10,9) 
PS_COEFF[6][2]=0.30916564*math.pow(10,6)
PS_COEFF[6][3]=0.75591105*math.pow(10,1) 
PS_COEFF[7][2]=-0.65537898*math.pow(10,5) 
PS_COEFF[7][3]=0.18810675*math.pow(10,3) 
PS_COEFF[8][0]=-0.14182435*math.pow(10,14) 
PS_COEFF[8][1]=0.18165390*math.pow(10,9) 
PS_COEFF[8][2]=-0.19769068*math.pow(10,6)
PS_COEFF[8][3]=-0.23530318*math.pow(10,2)
PS_COEFF[9][2]=0.92093375*math.pow(10,5)
PS_COEFF[9][3]=0.12246777*math.pow(10,3)

CS = np.zeros([10]) 


#Constants
#CONSTANT_B = 2451 #(Holyoke and Kronenberg, 2010) #3631 (Stipp and Tullis, 2003)
#EXPONENT =  -1.26
#1.45E4 -1.47 Twiss 1977 #Get more from literature
#339 -0.58 Koch 1983 # from Gleason and Tullis, 1993 (GRL)
PIEZOMETERS = {
    'ST03': {'Constant_B': 3631, 'Exponent': -1.26},
    'HK10': {'Constant_B': 2451, 'Exponent': -1.26},
    'K83': {'Constant_B': 339, 'Exponent': -0.58},
    'T77': {'Constant_B': 1.45E4, 'Exponent': -1.47}
}



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




def calculate_coefficient_table(temperature):
    for i in range(0, len(PS_COEFF)):
        CS[i]=PS_COEFF[i][0]*math.pow(temperature,-4)+PS_COEFF[i][1]*math.pow(temperature,-2)\
        +PS_COEFF[i][2]*math.pow(temperature,-1)\
        +PS_COEFF[i][3]+PS_COEFF[i][4]*temperature+PS_COEFF[i][5]*math.pow(temperature,2)



#CALCULATE EQUATIONS OF STATE AND FUGACITY
#Solve Equation of state, Eq 2 of Pitzer and Sterner (1994)
#Returns pressure in Pa
def eos(T, V):
    den = 1/V
    R = 8314472
    var_num = CS[2]+2*CS[3]*den+3*CS[4]*math.pow(den,2)+4*CS[5]*math.pow(den,3)
    var_denom = math.pow((CS[1]+CS[2]*den+CS[3]*math.pow(den,2)+CS[4]*math.pow(den,3)+CS[5]*math.pow(den,4)),2)
    pressure=den+CS[0]*math.pow(den,2)-math.pow(den,2)*(var_num/var_denom)
    pressure= pressure + (CS[6]*math.pow(den,2)*math.exp(-CS[7]*den)+CS[8]*math.pow(den,2)*math.exp(-CS[9]*den))
    pressure = pressure*(R*T) #pressure in Pa
    return pressure

#Solve for fugacity, Eq 1 of Pitzer and Sterner (1994)
#Returns fugacity in MPa
def PSfug(P,T,V):
    den=1/V;
    R=8314472;
    quotient = CS[0]*den+(1/(CS[1]+CS[2]*den+CS[3]*math.pow(den,2)+CS[4]*math.pow(den,3)+CS[5]*math.pow(den,4))-1/CS[1])
    quotient-= CS[6]/CS[7]*(math.exp(-CS[7]*den)-1)
    quotient-= CS[8]/CS[9]*(math.exp(-CS[9]*den)-1)
    lnf=(math.log(den)+ quotient+P/(den*R*T))+math.log(R*T)-1
    return math.exp(lnf)/1e6 # fugacity in MPa

#Optimizing equation to solve for volume
def fugacity_optimizer(temperature,pressure):

    def fun(v):
        return eos(temperature, v)- pressure
    volume = opt.brentq(fun, 5, 30) #Volume in cc/mo

    #Calculate fugacity 
    fugacity = PSfug(pressure, temperature, volume)
    
    return fugacity






    




    ##CALCULATE FUGACITY
#Takes imputs of pressure and temperature converts them from MPa and C to Pa and K 

class FlowStressCalculator():
    def __init__(self, temperature_values, pressure_values):
        self.temperature = C2K(np.array(temperature_values))
        self.pressure = np.array(pressure_values)*1.0E6
        self.grain_size = []
        self.fugacity = []
        self.differential_stress = []
        self.strain_rate = []
        self.slip_rate= []  

        
    def calculate_fugacity(self):
        
        for t in self.temperature:
            for p in self.pressure:
                calculate_coefficient_table(t)
                fug = fugacity_optimizer(t,p)
                self.fugacity.append(fug)
        return self.fugacity


    def calculate_differential_stress(self, grain_size, paleopiezometer='HK10'):
        
        for grain in grain_size:
            part = math.exp((math.log(grain)-math.log(PIEZOMETERS[paleopiezometer]['Constant_B']))/PIEZOMETERS[paleopiezometer]['Exponent'])
            self.differential_stress.append(part)
        return self.differential_stress


    def calculate_strain_rate(self, flow_law='H01'): 
        
        sr = []
        for num, stress in enumerate(self.differential_stress):
            for t, f in zip(self.temperature, self.fugacity):
                sr_i = (FLOW_LAWS[flow_law]['A']*np.power(stress, FLOW_LAWS[flow_law]['n'])*np.power(f,1)*np.exp(-FLOW_LAWS[flow_law]['Q']/(8.3144598*t)))
                self.strain_rate.append(sr_i)
        
        #st_rate = [sr[x:x+len(self.temperature)] for x in xrange(0, len(sr), len(self.temperature))]
        #self.strain_rate = st_rate

        return self.strain_rate
    
    
    def group_strain_rate(self):
        strain_rate_chunks = [self.strain_rate[x:x+len(self.temperature)] for x in xrange(0, len(self.strain_rate), len(self.temperature))]
        return strain_rate_chunks


    def calculate_slip_rate(self, width): #width in m, output of mm/yr
        
        for w in width: 
            for strain in self.strain_rate:
                vel = w*1000*31536000*strain
                self.slip_rate.append(vel)
            
        return self.slip_rate



#import matplotlib.pyplot as plt
#from matplotlib import cm
#from matplotlib.backends.backend_pdf import PdfPages



#class PlotFlowStress():



    


