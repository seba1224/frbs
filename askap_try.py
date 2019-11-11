import pandas as pd
import numpy as np
from astropy import constants as cte
from astropy import units as u


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


askap_index = 4
askap_v = 1*10**6*u.Hz


#askap info
aperture_eff =0.8 
element_area = 113
total_area = 113*36*aperture_eff*(u.m**2) #this should be the effective area of the array

tau = (2.*cte.k_B*frb_Tsys[askap_index]*u.K*frb_SNR[askap_index]/(frb_s[askap_index]*u.Jy*total_area))**2*1/askap_v
#print(tau.to(u.ms))



#second method

t = 4.15*10**6*(1./(frb_init_frec[askap_index])**2-1./(frb_end_frec[askap_index])**2)*frb_DM[askap_index]

#print(str(t)+' [ms]')

#print('\n'+'FRB170827:')


#frb170827 --- we have the image :) UTMOST
frb_index = 31
DM = frb_DM[frb_index]
f1 = frb_init_frec[frb_index]
f2 = frb_end_frec[frb_index]
T = 100*u.K
S = frb_s[frb_index]*u.Jy
Aeff = 18*10**3*(u.m)**2
SN = frb_SNR[frb_index]
nu = 0.4

t = 4.15*10**6*(f1**(-2)-f2**(-2))*DM 
print(str(t)+' [ms]  \t'+'DM equation'+'\t consistent with the image')

tau = (2*cte.k_B*T*SN/(S*Aeff))**2*1./(nu*10**6*u.Hz)
print(str(tau.to(u.ms))+'  \t'+'radiometer eq')



























