import pandas as pd
import numpy as np
from astropy import constants as cte
from astropy import units as u



area = {'UTMOST':18*10**3,}    

#UTMOST 4.4*11.6 each module, they have 352 of each...

#data from FRB catalog
df = pd.read_csv('frbs.csv', delimiter=',')
data = df.get_values()

frb_name = data[:,0]
frb_bw = data[:,4]              #MHz
frb_center = data[:,5]     #MHz
frb_s = data[:,6]               #in Jansky
frb_flux = data[:,7]            #in Jy ms
frb_widht = data[:,8]           #in ms
frb_DM = data[:,9]              #cm**-3*pc

for i in range(len(frb_DM)):
    ind = np.char.find(frb_DM[i], '\xc2')
    frb_DM[i] = float(frb_DM[i][:ind])

frb_SNR = data[:,10]
frb_Tsys = data[:,11]

frb_init_frec = frb_center-frb_bw/2.
frb_end_frec = frb_center+frb_bw/2.


def delta_t(f1,f2,DM):
    """f1, f2 in MHz; the output is in ms
    """
    cte = 4.15*10*6
    output = cte*(1./(f1**2)-1./(f2**2))*DM
    return output










