import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0, 10, 100)

fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 6))

axes[0].plot(x, np.sin(x))
axes[1].plot(x, np.cos(x))
axes[2].plot(x, np.tan(x))

# Remove vertical space
fig.subplots_adjust(hspace=0)

# Hide x tick labels on upper plots
for ax in axes[:-1]:
    ax.label_outer()

plt.show()
