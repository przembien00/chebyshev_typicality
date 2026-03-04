import numpy as np
import matplotlib.pyplot as plt


B_array = np.linspace(-4,4, 20)

plt.plot(B_array, np.tanh(B_array))
plt.savefig("Plots/tanh_plot.pdf", dpi=1000)