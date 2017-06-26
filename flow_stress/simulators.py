#SIMULATORS - work in progress/future work

width=30
strain_rate=1.5E-11


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

def simulate_strain_rate(width_min, width_max): #constant stress condition
    v = 15
    e = []
    for x in range(width_min, width_max):
        strain_rate = (v/31536000)/(x*1000)
        e.append(strain_rate)
    return e

def simulate_slip_rate(width_min, width_max): #constant strain rate condition
    e = 1.0E-11
    v = []
    for x in range(width_min, width_max):
        slip_rate = (e*31536000)*(x*1000)
        v.append(slip_rate)
    return v

## Simulators, this is for future work/ work in progress.
from mpmath import mp

#define conversion from mm/yr to mm/sec
conv = float(mp.fdiv(1, 31536000))

#grain_size = [5, 6, 7, 8, 9, 10]
def simulate_width(slip_rate, strain_rate):
    #v=w*e
    width = []
    for slip in slip_rate:
        for strain in strain_rate:
            w = (slip/(strain*31536000))/1000
            if w < 30 or w == 30:
                
                width.append(w)
            else:
                width.append(str("ERROR: Too wide"))
            
    return width




def simulate_slip_and_strain(width, slip_rate=14.7, strain_rate=1.557027041956803e-11):
    slip_rate = []
    strain = []
    width = [width]
        
    print(len(slip_rate))
    
    if len(slip_rate) == 1: #constant stress condition
        for w in width:
            st = slip_rate*conv/(w*1000)
        return st
    
    elif len(strain_rate) == 1: #constant strain rate condition
        for w in width:
            sl = strain_rate*31536000*w*1000
        return sl




def simulator(width, strain_rate):
    v=15.
    
    w=(v/1000)/(strain_rate*31536000)
    
    e=(v/31536000)/(width*1000)

    return w, e, v



    
    