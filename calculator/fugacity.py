import numpy as np 
from scipy.constants.constants import C2K
import math
from scipy import optimize as opt # for optimization


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

cs = np.zeros([10]) 

def calculate_coefficient_table(temperature):
    for i in range(0, len(PS_COEFF)):
        cs[i]=PS_COEFF[i][0]*math.pow(temperature,-4)+PS_COEFF[i][1]*math.pow(temperature,-2)\
        +PS_COEFF[i][2]*math.pow(temperature,-1)\
        +PS_COEFF[i][3]+PS_COEFF[i][4]*temperature+PS_COEFF[i][5]*math.pow(temperature,2)



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

#Optimizing equation to solve for volume
def fugacity_optimizer(temperature,pressure):

    def fun(v):
        return eos(temperature, v)- pressure
    volume = opt.brentq(fun, 5, 30) #Volume in cc/mol

    #Calculate fugacity 
    fugacity = PSfug(pressure, temperature, volume)
    
    return fugacity






    


#CALCULATE DIFFERENTIAL STRESS
#takes paleopiezometer constants from the literature and input of grain size in microns

#Constants
CONSTANT_B = 2451 #(Holyoke 2010) #3631 (Stipp and Tullis, 2003)
EXPONENT =  -1.26
#1.45E4 -1.47 Twiss 1977 #Get more from literature


#CALCULATE STRAIN RATE
#Takes flow law constants from the literature and input of differential stress and fugacity, pressure and temperature

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


    ##CALCULATE FUGACITY
#Takes imputs of pressure and temperature converts them from MPa and C to Pa and K 

class FlowStressCalculator():
    def __init__(self, temperature_values, pressure_values, grain_size):
        self.temperature = C2K(np.array(temperature_values))
        self.pressure = np.array(pressure_values)*1.0E6
        self.grain_size = grain_size
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


    def calculate_differential_stress(self):
        
        for grain in self.grain_size:
            part = math.exp((math.log(grain)-math.log(CONSTANT_B))/EXPONENT)
            self.differential_stress.append(part)
        return self.differential_stress


    def calculate_strain_rate(self, flow_law='H01'): 
          
        for t, f in zip(self.temperature, self.fugacity):
            for stress in self.differential_stress:
                sr = (FLOW_LAWS[flow_law]['A']*np.power(stress, FLOW_LAWS[flow_law]['n'])*np.power(f,1)*np.exp(-FLOW_LAWS[flow_law]['Q']/(8.3144598*t)))
                self.strain_rate.append(sr)
        
        return self.strain_rate

    def calculate_slip_rate(self, strain_rate, width): #width in m, output of mm/yr
        
        for w in width: 

            for strain in strain_rate:
            
                vel = w*1000*31536000*strain
                self.slip_rate.append(vel)
            
        return self.slip_rate

    


