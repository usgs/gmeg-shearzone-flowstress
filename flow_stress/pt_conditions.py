#Defining PT Condition calculators and depth simulator, 
#alternatively use a range of temperatures as a given depth and pressure

def pressure_calculator(depth, density): #depth in km, density in g/cm3
    g = 9.8 #gravity
    h = depth*1000 #converting to SI
    row = density*1000 #converting to SI
    pressure = (row*g*h)/1.0E6 # pressure in MPa
    return pressure 

def temperature_calculator(depths, geothermal_gradients): #km, C/km
    #np.array()
    for depth in depth_and_geothermal_gradients: 
        temp = geothermal_gradient*depth 
    return np.array(temp)
    
def depth_simulator(depth_range, density=2.7, geothermal_gradient=30): #depth range in km
    depth_conds = []
    for dep in range(depth_range[0], depth_range[1]):
        ptc = (round(pressure_calculator(dep, density)), temperature_calculator(geothermal_gradient, dep)) 
        depth_conds.append(ptc)
    return depth_conds