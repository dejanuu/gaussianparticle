import numpy as np
import matplotlib.pyplot as plt

# Assuming your file is named "data.txt"
import sys
filename = sys.argv[1]
data = np.loadtxt(filename)

# Let's plot the first column (time?) vs the second column (some measurement)
x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

plt.plot(x, (-1*z/y), marker='o')
plt.xlabel('Angle')
plt.ylabel('-S12/S11')
plt.title('My Data Plot')
plt.grid(True)
plt.show()
