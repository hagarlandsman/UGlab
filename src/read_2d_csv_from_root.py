import numpy as np
import matplotlib.pyplot as plt

density = 2.6
data = np.loadtxt("slant_2d_200520025_1730.csv", delimiter=",")
x = data[:, 0]
y = data[:, 1]
z = data[:, 2]
z = z * density
# Create 2D contour plot
contour = plt.tricontourf(x, y, z, levels=100, cmap='viridis')

# Set color scale range
contour.set_clim(300 * density , 1400 * density)

# Add colorbar and set label with LaTeX
cbar = plt.colorbar(contour)
cbar.set_label(r"Overburden [m.w.e.]")

# Axis labels and title with LaTeX
plt.xlabel(r"Azimuth $\phi$ [deg]")
plt.ylabel(r"$\cos(\theta)$")
# plt.title(r"Slant Depth")
plt.ylim(bottom=0.35)

plt.tight_layout()
plt.savefig("slant_depth_plot.png", dpi=300, bbox_inches='tight')

plt.show()
