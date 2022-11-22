import pandas
from pandas import read_csv
from scipy.optimize import curve_fit
from matplotlib import pyplot
import numpy as np

data_list=[]

for i in range(3):
    data_list.append(read_csv('dynamics_{}.txt'.format(i)).values)
data = np.concatenate(data_list,axis=0)
print(data.shape)
data = np.delete(data,data[:,2]<-25,axis=0)
print(data.shape)

data = np.delete(data,data[:,2]>25,axis=0)
print(data.shape)

# define the true objective function
def objective(x, cm0, cm1,cm2,cm3):
    return cm1 * x[:,0] -cm2*x[:,1]-cm3*x[:,1]*x[:,1]-cm0
x, a = data[:, 0:2], data[:,2]
popt, _ = curve_fit(objective, x,a)
print(popt)

def objective1(x, cm0, cm1,cm2):
    return cm1 * x[:,0] -cm2*x[:,1]-cm0
popt, _ = curve_fit(objective1, x,a)
print(popt)


def objective2(x, cm0, cm1):
    return cm1 * x -cm0
popt, _ = curve_fit(objective2, x[:,0],a)
print(popt)


cm0,cm1 = popt
pyplot.scatter(x[:,0], a)
x_line = np.arange(min(x[:,0]), max(x[:,0]), 0.001)
a_fit = objective2(x_line, cm0, cm1)
pyplot.plot(x_line, a_fit, '--', color='red')
pyplot.show()