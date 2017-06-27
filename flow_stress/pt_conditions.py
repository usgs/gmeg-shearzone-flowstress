#Defining PT Condition calculators and depth simulator, 
#alternatively use a range of temperatures as a given depth and pressure
class PTCalculator():

    def __init__(self, depth, density, geothermal_gradient):
        self.gravity = 9.8
        self.density = density*1000
        self.depth = depth
        self.geothermal_gradient = geothermal_gradient

    def pt_calculator(self): #depth in km, density in g/cm3

        pressure = [(self.density*self.gravity*d*1000)/1.0E6 for d in self.depth] #pressure in MPa
        temp = [self.geothermal_gradient*d for d in self.depth]
        
        return pressure, temp



    
    # def depth_simulator(self): #depth range in km
        
    #     depth_conds = []
    #     for dep in range(depth_range[0], depth_range[1]):
    #         ptc = (round(pressure_calculator(dep, density)), temperature_calculator(geothermal_gradient, dep)) 
    #         depth_conds.append(ptc)
        
    #     return depth_conds