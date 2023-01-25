import numpy as np
mass_len=5
start_step=1e-5
end_step=5e3
mass_set=np.geomspace(start_step,end_step,mass_len)
print(mass_set)

# for i,mass in enumerate(mass_set):
#     print(mass)
#     np.save(f'sigma{i}.npy',mass)

sigma_set=np.zeros((5,1))
for i in range(len(sigma_set)):
    sigma_set[i]=2*np.load(f"sigma{i}.npy")

print(sigma_set)


