import sys
import numpy as np
from scipy import odr
from scipy import constants as con
from pylab import *

q=con.elementary_charge
k=con.Bolzmann
T=293

ctr = 0
D = {}
V = list()
J = list()
#with open(sys.argv[1]) as f:
with open('/home/timo/Documents/Studium/master/minor/data/2013-06-24/IV/W1715d2.stf') as f:
    for line in f:
        ctr=ctr+1
        if ctr>18:
            a,b=line.split("\t")
            V.append(float(a.strip()))
            J.append(float(b.strip()))
        sepl=line.split('\t')
        if sepl[0].endswith(":"):
            D[sepl[0].strip().replace(":","")]=sepl[1].strip()

print D['Voc']

V = array(V)
J = array(J)

def ddiode(p,x):
    return(np.abs(p[4]) - np.abs(p[0]) * (np.exp(np.divide(q*(x[0]-x[1]*np.abs(p[2])) ,  k*T ))-1)
		         - np.abs(p[1]) * (np.exp(np.divide(q*(x[0]-x[1]*np.abs(p[2])),(2*k*T)))-1)
		         - np.divide(x[0]-x[1]*np.abs(p[2]),np.abs(p[3])) -x[1])
def jacd(p,x):
    return(np.vstack([],[]))

def jacb(p,x):
    e1 = np.exp(np.divide(q*(x[0]-x[1]*p[2]**2),  (k*T)))-1
    e2 = np.exp(np.divide(q*(x[0]-x[1]*p[2]**2),(2*k*T)))-1
    res=np.vstack([	-e1,
			-e2,
	      np.divide( x[1]*(2*p[0]*p[3]*q*(e1+1)+p[1]*p[3]*q*(e2+1)), 2*k*T*p[3] ),
	      np.divide( x[0]-x[1]*p[2] , np.power(p[3],3)  ),
			np.multiply(np.ones(x.shape[-1]),1) ])
    return res


mod = odr.Model(ddiode,implicit=1,fjacb=jacb)
dat = odr.Data([V,J],1) 
b = np.sqrt([np.power(10,-12),np.power(10,-8),1,12000,0.024])
imp_odr = odr.ODR(dat,mod,beta0=b)
out = imp_odr.run()
imp_odr = odr.ODR(dat,mod,beta0=out.beta)
out = imp_odr.run()
print out.beta
