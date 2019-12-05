#exec(open('brute_curvefit.py').read())

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def bruteforce(func, x, y, bounds, ytol=1e-5, ntol = 10000, returnnfactor = 0.01):
    '''
    bounds: is an array of the form (min,max) where min is an array of minimum values possible for each of the free arguements.
    func: Any function with the first parameter as independednt parameter. The independednt parameter should be a range in which the independednt parameter can vary
    x: The independent parameter (list)
    y: The target. Should be over x
    ntol: Number of iterations
    returnnfactor: Factor of ntol to return. ie, how many list of parameters would be returned is ntol*returnnfactor
    '''
    arguements = func.__code__.co_varnames[1:]
    returnn = int(ntol*returnnfactor)
    bounds = np.array(bounds)
    paramlist = []
    errorlist = []
    for n in np.arange(ntol):
        currparam = []
        for i in np.arange(len(arguements)):
            currparam.append(np.random.rand(1)*(bounds[1,i]-bounds[0,i]) + bounds[0,i])
        # print(currparam)
        error = np.sum((func(x,*currparam)-y)**2)
        # print(error)
        paramlist.append(currparam)
        errorlist.append(error)
        print(i/ntol,end='\r')

    best_error_idx = np.argpartition(errorlist,returnn)[:returnn]
    best_params = np.array(paramlist)[best_error_idx]
    # best_error_idx = np.array(errorlist).argmin()
    # best_param = paramlist[best_error_idx]
    return [best_params, np.array(errorlist)[best_error_idx]]

def scipy_fit(func, x, y, bounds, p0list):
    '''
    bounds: is an array of the form (min,max) where min is an array of minimum values possible for each of the free arguements.
    func: Any function with the first parameter as independednt parameter. The independednt parameter should be a range in which the independednt parameter can vary
    x: The independent parameter (list)
    y: The target. Should be over x
    ntol: Number of iterations
    returnnfactor: Factor of ntol to return. ie, how many list of parameters would be returned is ntol*returnnfactor
    '''
    fitparams_list=[]
    error_list=[]
    for p0 in p0list:
        p0 = np.ravel(p0)
        fittedparam,cov = curve_fit(func, x, y, bounds=bounds, p0=p0)
        error = np.sum((func(x,*fittedparam)-y)**2)
        fitparams_list.append(fittedparam)
        error_list.append(error)
    best_error_idx = np.array(error_list).argmin()
    best_param = np.array(fitparams_list)[best_error_idx]
    return [best_param, np.array(error_list)[best_error_idx]]

def brute_then_scipy(func, x, y, bounds, ytol=1e-5, ntol = 10000, returnnfactor = 0.01):
    '''
    bounds: is an array of the form (min,max) where min is an array of minimum values possible for each of the free arguements.
    func: Any function with the first parameter as independednt parameter. The independednt parameter should be a range in which the independednt parameter can vary
    x: The independent parameter (list)
    y: The target. Should be over x
    ntol: Number of iterations
    returnnfactor: Factor of ntol to return. ie, how many list of parameters would be returned is ntol*returnnfactor
    '''
    paramsfitted,errors = bruteforce(func,x,y,bounds=bounds,ntol=ntol,returnnfactor=returnnfactor)
    paramfitted,error = scipy_fit(func,x,y,bounds=bounds, p0list=paramsfitted)
    return [paramfitted,error]

if __name__ == '__main__':
    def h(v, vhalf, k):
        return 1/(1+np.exp((v-vhalf)/-k))
    v = np.linspace(-0.100,0.100,3000)
    hinf = h(v,-0.050,-0.004)
    plt.plot(v,hinf, label='original')
    paramsfitted,errors = bruteforce(h,v,hinf,bounds=[[-1,-1],[1,1]])
    for param in paramsfitted:
        plt.plot(v,h(v,*param), label='fitted')
    plt.legend()
    plt.show()

    plt.plot(v,hinf, label='original')
    paramfitted,error = scipy_fit(h,v,hinf,bounds=[[-1,-1],[1,1]], p0list=paramsfitted)
    plt.plot(v,h(v,*paramfitted), label='fitted')
    plt.legend()
    plt.show()
