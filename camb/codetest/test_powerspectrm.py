import numpy as np
from matplotlib import pyplot as plt

# test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
# insert_first_two_rows=np.array([[0,0,0],[0,0,0]])
# dls=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
# print(dls)
# ls=np.arange(len(dls[:,0]))
# plt.loglog(ls,dls[:,0])
# plt.show()



ali_noise_dls=np.loadtxt('Ali_noise.dat',usecols=(1,2))
print(ali_noise_dls)
insert_first_thirty_rows=np.zeros((30,2))
half_nls=np.insert(ali_noise_dls,0,insert_first_thirty_rows,axis=0)
insert_last_many_rows=np.zeros((2501-len(half_nls),2))
nls=np.insert(half_nls,len(half_nls),insert_last_many_rows,axis=0)


ls=np.arange(2501)
plt.plot(ls,nls[:,0])
plt.show()
