import fileinput
import numpy as np
import re
import os
import sys
import time
from matplotlib import pyplot as plt
# from matplotlib import rc # latex in matplotlib
# rc('text',usetex=True)


def set_l_max_scalar(new_value):
    for line in fileinput.input('test.ini',inplace=True):
        if line.startswith('l_max_scalar'):
            print('l_max_scalar = '+str(new_value))
        else:
            print(line.strip())


set_l_max_scalar(4400)
os.system('./camb test.ini')
test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
dlsn=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
print(dlsn.shape)
ls=np.arange(2501)
plt.figure(1)
plt.loglog(ls, dlsn[:2501,0])
plt.figure(2)
plt.semilogx(ls, dlsn[:2501,1])
plt.figure(3)
plt.semilogx(ls, dlsn[:2501,2])

set_l_max_scalar(4400)
os.system('./camb test.ini')
test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
dlsn=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
print(dlsn.shape)
plt.figure(1)
plt.loglog(ls, dlsn[:2501,0])
plt.figure(2)
plt.semilogx(ls, dlsn[:2501,1])
plt.figure(3)
plt.semilogx(ls, dlsn[:2501,2])
plt.show()

