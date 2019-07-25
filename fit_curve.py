import numpy as np
import symfit as sf
from symfit.core.minimizers import BFGS, NelderMead, SLSQP, Powell
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import argparse

def fit_data(data,guess=10.0,use_err=False):
    
    x_data = data[:,0]
    y_data = 1- data[:,1]
    if use_err:
       data_err= data[:,2]
    
    x, y = sf.variables('x, y')
    pKa, n = sf.parameters('pKa, n')
    model = sf.Model({y: (10**(n*(pKa-x))+1)**-1.0})
    pKa.value=guess
    n.value=1
    if use_err:
       fit = sf.Fit(model, x=x_data, y=y_data,sigma_y=data_err ,minimizer=Powell, absolute_sigma=True)
    else:
       fit = sf.Fit(model, x=x_data, y=y_data,minimizer=Powell)
    result = fit.execute()
    
    print("pKa.....................................", result.value(pKa),'+/-', result.stdev(pKa))
    print("n.......................................", result.value(n),'+/-', result.stdev(n))
    print("Regression coefficent:................", result.r_squared,'\n')
    
    x_out = np.arange(min(x_data),max(x_data),10**-3.0)
    y_out = fit.model(x=x_out, **result.params)[0]
    
    return x_out ,y_out, result.value(pKa), result.stdev(pKa), result.r_squared,result.value(n), result.stdev(n)


parser = argparse.ArgumentParser(description='Fit and plot titration curves.')
parser.add_argument('-f', dest='name', type=str,help='filename')
parser.add_argument('-n',dest='n', type=int,help='number of titratable sides')

args = parser.parse_args()


data=np.loadtxt(args.name)
data[:,1] = data[:,1]/args.n
data[:,2] = data[:,2]/args.n

x, y, pKa, err, R, n, nerr = fit_data(data,guess=10,use_err=True)

plt.plot(x,y,label=str(round(pKa,2)),c="#0039A6")

plt.errorbar(data[:,0],1-data[:,1],yerr=data[:,2],ls=' ',
                 lw=2,marker='>',c="#0039A6",markersize='7')

plt.legend()
plt.xlabel('pH',fontsize='12')
plt.ylabel('degree of deprotonation',fontsize='12')
plt.tight_layout()
plt.show()
