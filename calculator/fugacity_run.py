from fugacity import FugacityCalculator

t = [450]
p = [400]

f = FugacityCalculator(t,p)
fugacity_calc = f.calculate_fugacity()

print(fugacity_calc)