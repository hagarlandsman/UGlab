import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("hist2d.csv", delimiter=",")
x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

# Create 2D histogram for plotting
plt.tricontourf(x, y, z, levels=100, cmap='viridis')
plt.colorbar(label="Counts")

plt.xlabel("Azimuth φ [deg]")
plt.ylabel("cos(θ)")
plt.title("Slant Depth")
plt.clim(300, 1400)  # Z-axis (color scale) range
cbar.set_label("Overburden [m]")
plt.show()
