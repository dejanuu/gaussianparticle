import numpy as np

data = np.loadtxt("outputS.out")

theta_array = data[:, 0]
S11_array   = data[:, 1]

np.savez("target.npz", theta=theta_array, S11=S11_array)
