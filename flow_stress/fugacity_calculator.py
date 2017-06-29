import numpy as np
from scipy import optimize as opt 
import math



class FugacityCalculator():

    def __init__(self, temperature, pressure):
        self.temperature = temperature
        self.pressure = pressure
        self.CS = np.zeros([10])
        self.PS_COEFF = np.zeros([10,10])
        self.PS_COEFF[0][2]=0.24657688*math.pow(10,6)
        self.PS_COEFF[0][3]=0.51359951*math.pow(10,2) 
        self.PS_COEFF[1][2]=0.58638965*math.pow(10,0) 
        self.PS_COEFF[1][3]=-0.28646939*math.pow(10,-2) 
        self.PS_COEFF[1][4]=0.31375577*math.pow(10,-4) 
        self.PS_COEFF[2][2]=-0.62783840*math.pow(10,1) 
        self.PS_COEFF[2][3]=0.14791599*math.pow(10,-1) 
        self.PS_COEFF[2][4]=0.35779579*math.pow(10,-3) 
        self.PS_COEFF[2][5]=0.15432925*math.pow(10,-7) 
        self.PS_COEFF[3][3]=-0.42719875*math.pow(10,0) 
        self.PS_COEFF[3][4]=-0.16325155*math.pow(10,-4) 
        self.PS_COEFF[4][2]=0.56654978*math.pow(10,4) 
        self.PS_COEFF[4][3]=-0.16580167*math.pow(10,2) 
        self.PS_COEFF[4][4]=0.76560762*math.pow(10,-1) 
        self.PS_COEFF[5][3]=0.10917883*math.pow(10,0) 
        self.PS_COEFF[6][0]=0.38878656*math.pow(10,13) 
        self.PS_COEFF[6][1]=-0.13494878*math.pow(10,9) 
        self.PS_COEFF[6][2]=0.30916564*math.pow(10,6)
        self.PS_COEFF[6][3]=0.75591105*math.pow(10,1) 
        self.PS_COEFF[7][2]=-0.65537898*math.pow(10,5) 
        self.PS_COEFF[7][3]=0.18810675*math.pow(10,3) 
        self.PS_COEFF[8][0]=-0.14182435*math.pow(10,14) 
        self.PS_COEFF[8][1]=0.18165390*math.pow(10,9) 
        self.PS_COEFF[8][2]=-0.19769068*math.pow(10,6)
        self.PS_COEFF[8][3]=-0.23530318*math.pow(10,2)
        self.PS_COEFF[9][2]=0.92093375*math.pow(10,5)
        self.PS_COEFF[9][3]=0.12246777*math.pow(10,3)




    def calculate_coefficient_table(self):
         
        for i in range(0, len(self.PS_COEFF)):
            self.CS[i]=self.PS_COEFF[i][0]*math.pow(self.temperature,-4)+self.PS_COEFF[i][1]*math.pow(self.temperature,-2)\
            +self.PS_COEFF[i][2]*math.pow(self.temperature,-1)\
            +self.PS_COEFF[i][3]+self.PS_COEFF[i][4]*self.temperature+self.PS_COEFF[i][5]*math.pow(self.temperature,2)
        #return self.CS

    #CALCULATE EQUATIONS OF STATE AND FUGACITY
    #Solve Equation of state, Eq 2 of Pitzer and Sterner (1994)
    #Returns pressure in Pa
    def eos(self, T, V):
        #self.V = V
        den = 1/V
        R = 8314472
        var_num = self.CS[2]+2*self.CS[3]*den+3*self.CS[4]*math.pow(den,2)+4*self.CS[5]*math.pow(den,3)
        var_denom = math.pow((self.CS[1]+self.CS[2]*den+self.CS[3]*math.pow(den,2)+self.CS[4]*math.pow(den,3)+self.CS[5]*math.pow(den,4)),2)
        pressure=den+self.CS[0]*math.pow(den,2)-math.pow(den,2)*(var_num/var_denom)
        pressure= pressure + (self.CS[6]*math.pow(den,2)*math.exp(-self.CS[7]*den)+self.CS[8]*math.pow(den,2)*math.exp(-self.CS[9]*den))
        pressure = pressure*(R*T) #pressure in Pa
        return pressure

    #Solve for fugacity, Eq 1 of Pitzer and Sterner (1994)
    #Returns fugacity in MPa
    def PSfug(self,P,T,V):
        den=1/V;
        R=8314472;
        quotient = self.CS[0]*den+(1/(self.CS[1]+self.CS[2]*den+self.CS[3]*math.pow(den,2)+self.CS[4]*math.pow(den,3)+self.CS[5]*math.pow(den,4))-1/self.CS[1])
        quotient-= self.CS[6]/self.CS[7]*(math.exp(-self.CS[7]*den)-1)
        quotient-= self.CS[8]/self.CS[9]*(math.exp(-self.CS[9]*den)-1)
        lnf=(math.log(den)+ quotient+P/(den*R*T))+math.log(R*T)-1
        return math.exp(lnf)/1e6 # fugacity in MPa

    def fugacity_optimizer(self):

        def fun(v):
            return self.eos(self.temperature, v)- self.pressure

        volume = opt.brentq(fun, 5, 30) #Volume in cc/mo

        #Calculate fugacity 
        fugacity = self.PSfug(self.pressure, self.temperature, volume)
        
        return fugacity





###same thing but not in class ...

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

	CS = np.zeros([10]) 
	
	def calculate_coefficient_table(temperature):
    
	    for i in range(0, len(PS_COEFF)):
	        CS[i]=PS_COEFF[i][0]*math.pow(temperature,-4)+PS_COEFF[i][1]*math.pow(temperature,-2)\
	        +PS_COEFF[i][2]*math.pow(temperature,-1)\
	        +PS_COEFF[i][3]+PS_COEFF[i][4]*temperature+PS_COEFF[i][5]*math.pow(temperature,2)
	    #return CS

	def fun(v):
		return eos(temperature, v)- pressure

	volume = opt.brentq(fun, 5, 30) #Volume in cc/mo

    #Calculate fugacity 
	fugacity = PSfug(pressure, temperature, volume)
    
	return fugacity



