import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

# Define the 2D Gaussian function
def gaussian_2d(x, y, x0, y0, sigma_x, sigma_y, amplitude, offset):
    return offset + amplitude * np.exp(
        -((x - x0)**2 / (2 * sigma_x**2) + (y - y0)**2 / (2 * sigma_y**2))
    )

# Create sample data
x = np.linspace(0, 10, 100)
y = np.linspace(0, 10, 100)
x, y = np.meshgrid(x, y)
data = gaussian_2d(x, y, x0=5, y0=5, sigma_x=1, sigma_y=1, amplitude=10, offset=1)

# Add some noise
data_noisy = data + 0.2 * np.random.normal(size=data.shape)

# Flatten the data for fitting
x_data = x.ravel()
y_data = y.ravel()
data_noisy_flat = data_noisy.ravel()

# Define the objective function for optimization
def objective(params):
    x0, y0, sigma_x, sigma_y, amplitude, offset = params
    return gaussian_2d(x_data, y_data, x0, y0, sigma_x, sigma_y, amplitude, offset) - data_noisy_flat

# Initial guess for the parameters
initial_guess = (5, 5, 1, 1, 10, 1)

# Perform the fit
params_opt, _ = opt.leastsq(objective, initial_guess)

# Extract the optimized parameters
x0_opt, y0_opt, sigma_x_opt, sigma_y_opt, amplitude_opt, offset_opt = params_opt

# Plot the results
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
ax[0].imshow(data_noisy, extent=(0, 10, 0, 10), origin='lower')
ax[0].set_title('Noisy Data')

fitted_data = gaussian_2d(x, y, x0_opt, y0_opt, sigma_x_opt, sigma_y_opt, amplitude_opt, offset_opt)
ax[1].imshow(fitted_data, extent=(0, 10, 0, 10), origin='lower')
ax[1].set_title('Fitted Data')

plt.show()
