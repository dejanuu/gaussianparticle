import numpy as np
import matplotlib.pyplot as plt

import sys
filename = sys.argv[1]
data = np.loadtxt(filename)

x = data[range(5, 175),  0]
y = data[range(5, 175), 1]
z = data[range(5, 175), 2]

ratio = -z / y

max_index = np.argmax(ratio)
max_ratio = ratio[max_index]
max_x = x[max_index]

print("Max (-z/y):", max_ratio)
print("Occurs at x:", max_x)

plt.plot(x, ratio, marker='o')
plt.xlabel('Angle')
plt.ylabel('-S12/S11')
plt.title('My Data Plot')
plt.grid(True)

plt.scatter([max_x], [max_ratio], color='red', s=80, label='Max value')
plt.axvline(x=max_x, color='red', linestyle='--', linewidth=1)
plt.axhline(y=max_ratio, color='red', linestyle='--', linewidth=1)
plt.savefig("testoutput1.jpg")
plt.show()


