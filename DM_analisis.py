import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as cte
from astropy import constants as u
from scipy.optimize import newton

def delta_t(f1,f2,DM, time):
    """f1, f2 in MHz; time in ms
    """
    cte = 4.15*10*6
    output = cte*(1./(f1**2)-1./(f2**2))*DM-time
    return output


f1 = 1000
f2 = 2000
tau = 2.**25/(62.5*10**3)
def test(DM):
    f1 = 1000
    f2 = 1200
    tau = 2.**25/(62.5*10**3)
    output = delta_t(f1,f2,DM,tau)
    return output


zero = newton(test, 300)
print(zero)


