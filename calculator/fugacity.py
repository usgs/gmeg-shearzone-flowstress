import numpy as np 
from scipy.constants.constants import C2K
import math
from scipy import optimize as opt # for optimization

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
def fugacity_calculator(temperature,pressure):

    def fun(v):
        return eos(temperature, v)- pressure
    volume = opt.brentq(fun, 5, 30) #Volume in cc/mol

    #Calculate fugacity 
    fugacity = PSfug(pressure, temperature, volume)
    
    return fugacity


PScoeff = np.zeros([10,10])

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

cs = np.zeros([10]) 

def calculate_coefficient_table(temperature):
    for i in range(0, len(PScoeff)):
        cs[i]=PScoeff[i][0]*math.pow(temperature,-4)+PScoeff[i][1]*math.pow(temperature,-2)\
        +PScoeff[i][2]*math.pow(temperature,-1)\
        +PScoeff[i][3]+PScoeff[i][4]*temperature+PScoeff[i][5]*math.pow(temperature,2)

class FugacityCalculator():
    def __init__(self, temperature_values, pressure_values):
        self.temperature = C2K(np.arange(temperature_values[0],temperature_values[1]+1))
        self.pressure = self.convert_mpa_to_pa(pressure_values)
        self.temperature_and_pressure = self.combine_temp_and_pressure_values()
    
    @staticmethod
    def convert_mpa_to_pa(pressure_values):
        pressures = np.arange(pressure_values[0],pressure_values[1]+1)
        for i,p in enumerate(pressures):
            Ppa = round(p*1.0E6)
            pressures[i] = Ppa
        return pressures 
    
    def combine_temp_and_pressure_values(self):
        all_values = []
        for i,temp in enumerate(self.temperature):
            temp_and_pressure = (temp,self.pressure[i])
            all_values.append(temp_and_pressure)
        return all_values
        
    def calculate_fugacity(self): 
        temperature_pressure_fugacity = []
        for t,p in self.temperature_and_pressure: 
            calculate_coefficient_table(t)
            fugacity = fugacity_calculator(t,p)
            values = {'fugacity':fugacity,'temperature':t, 'pressure':p}
            temperature_pressure_fugacity.append(values)
        return temperature_pressure_fugacity