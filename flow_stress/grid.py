import numpy as np
from scipy import optimize as opt 
from scipy.constants.constants import C2K, K2C



PS_COEFF = np.zeros([10,10])
PS_COEFF[0][2]=0.24657688*np.power(10,6)
PS_COEFF[0][3]=0.51359951*np.power(10,2) 
PS_COEFF[1][2]=0.58638965*np.power(10,0) 
PS_COEFF[1][3]=-0.28646939*np.power(10,-2) 
PS_COEFF[1][4]=0.31375577*np.power(10,-4) 
PS_COEFF[2][2]=-0.62783840*np.power(10,1) 
PS_COEFF[2][3]=0.14791599*np.power(10,-1) 
PS_COEFF[2][4]=0.35779579*np.power(10,-3) 
PS_COEFF[2][5]=0.15432925*np.power(10,-7) 
PS_COEFF[3][3]=-0.42719875*np.power(10,0) 
PS_COEFF[3][4]=-0.16325155*np.power(10,-4) 
PS_COEFF[4][2]=0.56654978*np.power(10,4) 
PS_COEFF[4][3]=-0.16580167*np.power(10,2) 
PS_COEFF[4][4]=0.76560762*np.power(10,-1) 
PS_COEFF[5][3]=0.10917883*np.power(10,0) 
PS_COEFF[6][0]=0.38878656*np.power(10,13) 
PS_COEFF[6][1]=-0.13494878*np.power(10,9) 
PS_COEFF[6][2]=0.30916564*np.power(10,6)
PS_COEFF[6][3]=0.75591105*np.power(10,1) 
PS_COEFF[7][2]=-0.65537898*np.power(10,5) 
PS_COEFF[7][3]=0.18810675*np.power(10,3) 
PS_COEFF[8][0]=-0.14182435*np.power(10,14) 
PS_COEFF[8][1]=0.18165390*np.power(10,9) 
PS_COEFF[8][2]=-0.19769068*np.power(10,6)
PS_COEFF[8][3]=-0.23530318*np.power(10,2)
PS_COEFF[9][2]=0.92093375*np.power(10,5)
PS_COEFF[9][3]=0.12246777*np.power(10,3)

CS = np.zeros([10]) 


def eos(T, V):
    den = 1/V
    R = 8314472
    var_num = CS[2]+2*CS[3]*den+3*CS[4]*np.power(den,2)+4*CS[5]*np.power(den,3)
    var_denom = np.power((CS[1]+CS[2]*den+CS[3]*np.power(den,2)+CS[4]*np.power(den,3)+CS[5]*np.power(den,4)),2)
    pressure=den+CS[0]*np.power(den,2)-np.power(den,2)*(var_num/var_denom)
    pressure= pressure + (CS[6]*np.power(den,2)*np.exp(-CS[7]*den)+CS[8]*np.power(den,2)*np.exp(-CS[9]*den))
    pressure = pressure*(R*T) #pressure in Pa
    return pressure

#Solve for fugacity, Eq 1 of Pitzer and Sterner (1994)
#Returns fugacity in MPa
def PSfug(P,T,V):
    den=1/V;
    R=8314472;
    quotient = CS[0]*den+(1/(CS[1]+CS[2]*den+CS[3]*np.power(den,2)+CS[4]*np.power(den,3)+CS[5]*np.power(den,4))-1/CS[1])
    quotient-= CS[6]/CS[7]*(np.exp(-CS[7]*den)-1)
    quotient-= CS[8]/CS[9]*(np.exp(-CS[9]*den)-1)
    lnf=(np.log(den)+ quotient+P/(den*R*T))+np.log(R*T)-1
    return np.exp(lnf)/1e6 # fugacity in MPa

#Optimizing equation to solve for volume
def fugacity_optimizer(temperature,pressure):

    def fun(v):
        return eos(temperature, v)- pressure
    volume = opt.brentq(fun, 5, 30) #Volume in cc/mo

    #Calculate fugacity 
    fugacity = PSfug(pressure, temperature, volume)
    
    return fugacity

def calculate_coefficient_table(temperature):
    for i in range(0, len(PS_COEFF)):
        CS[i] = PS_COEFF[i][0]*np.power(temperature,-4)+PS_COEFF[i][1]*np.power(temperature,-2)\
        +PS_COEFF[i][2]*np.power(temperature,-1)\
        +PS_COEFF[i][3]+PS_COEFF[i][4]*temperature+PS_COEFF[i][5]*np.power(temperature,2)

        return CS


            ##CALCULATE FLOW STRESS
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
        
        for t, p in zip(self.temperature, self.pressure):
            calculate_coefficient_table(t)
            fug = fugacity_optimizer(t,p)
            self.fugacity.append(fug)
        return self.fugacity


    def calculate_differential_stress(self, grain_size, paleopiezometer='HK10'):
        
        self.grain_size = grain_size
        for grain in self.grain_size:
            part = np.exp((np.log(grain)-np.log(PIEZOMETERS[paleopiezometer]['Constant_B']))/PIEZOMETERS[paleopiezometer]['Exponent'])
            self.differential_stress.append(part)
        return self.differential_stress


    def calculate_strain_rate(self, flow_law='H01'): 
        
        for stress in self.differential_stress:
            for t, f in zip(self.temperature, self.fugacity):
                sr_i = (FLOW_LAWS[flow_law]['A']*np.power(stress, FLOW_LAWS[flow_law]['n'])*np.power(f,1)*np.exp(-FLOW_LAWS[flow_law]['Q']/(8.3144598*t)))
                self.strain_rate.append(sr_i)

        return self.strain_rate
    

    def calculate_slip_rate(self, width): #width in m, output of mm/yr
        
        for w in width: 
            for strain in self.strain_rate:
                vel = w*1000*31536000*strain
                self.slip_rate.append(vel)
            
        return self.slip_rate

