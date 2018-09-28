from JAC_Tools import error_propagation as ep
import numpy as np

def stellar_radius(variables):#variables=[L,teff]
    stbolt=5.670367*10**(-8)
    R=np.sqrt((variables[0])/(4*np.pi*stbolt*variables[1]**4))
    return R
def chromospheric_Activity_Age(logRHK):
    return 10**(-38.053-17.912*logRHK-1.6675*(logRHK)**2)


Lsun=3.828*10**(26)

truev_rad=[379.3149849736817*Lsun,4780]
er_rad=[(69.9059630128,-59.0274772996 ),(85,-85)]

truev_caa=[-6.0165]#,6025)
er_caa=[( 0.0800980536199,-0.0820910092786)]#,(96.24,-96.24)]

GS=ep.GridSearch(truev_caa,er_caa,chromospheric_Activity_Age)
lot=GS.errorpropagation(100)
print(lot)
GS=ep.GridSearch(truev_rad,er_rad,stellar_radius)
lot=GS.errorpropagation(100)
print(lot)
