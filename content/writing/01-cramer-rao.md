---
title: "Fisher Information Analysis of Dynamic Vision Sensors"
date: 2025-11-10
summary: "A verification of the Cramer-Rao lower bound derivation for 3D tracking using event cameras, based on 'CodedEvents'."
image: "/images/stats.jpg"
tags: ["Code", "Math", "Derivation"]
---

## 1. Optical Projection Model

The authors of [1] propose a **Cramer-Rao Lower Bound (CRLB)** for the precision of 3D object positioning from an event sensor, given a known **Point-Spread Function (PSF)**. 

### Ideal and Blurred Images
The ideal image $I_t$ at time $t$ is modeled as a Dirac delta function on the sensor plane $(u, v)$:
$$I_t(u, v) = \delta\left(u - \frac{f}{z}x(t), v - \frac{f}{z}y(t)\right) = \delta(u - u_0, v - v_0)$$

To account for optical blur, the real image $I^b_t$ is the convolution of the depth-variant PSF $h_z(t)$ with the ideal image:
$$I_t^b(u,v) = [h_z(t) * I_t](u, v)$$



### Coordinate Mapping
The mapping of 3D coordinates $(x, y, z)$ to sensor coordinates $(u, v)$ is defined by the focal length $f$ and the object depth $z + \Delta z(t)$:
$$u = f \frac{x(t)}{z + \Delta z(t)}, \quad v = f \frac{y(t)}{z + \Delta z(t)}$$

Assuming $f$ and $z$ are constant over small intervals, we use a scaling factor $S = \frac{f}{z}$ to approximate the projection:
$$(u, v) \approx S \cdot (x(t), y(t))$$

---

## 2. Event Measurement Model

Event cameras respond to changes in log-intensity. For accumulated events over a sufficiently long interval $\tau$, the measurement $O_t$ is approximately the log-difference of the blurred intensity:
$$O_{t} = \log(I_{t}^{b}) - \log(I_{t-\tau}^{b})$$

### Fisher Information Framework
The Fisher Information Matrix (FIM) $\mathcal{I}(\theta)$ measures the amount of information an observable random variable $X$ carries about an unknown parameter $\theta$:
$$\mathcal{I}(\theta)_{i,j} = \mathbb{E} \left[ \left( \frac{\partial}{\partial \theta_i} \log f(X; \theta) \right) \left( \frac{\partial}{\partial \theta_j} \log f(X; \theta) \right) \Bigg| \theta \right]$$

| Component | Physical Meaning | Role in Optimization |
| :--- | :--- | :--- |
| **Parameter** $\theta$ | Ground-truth 3D positions $\{x, y, z\}$ at $t$ and $t-\tau$. | The variables to be estimated. |
| **Ideal** $I_t$ | Perfect pin-hole projection (Dirac Delta). | The input to the optical model. |
| **PSF** $h_z$ | Depth-dependent blur pattern. | The engineered design variable. |
| **Image** $I_t^b$ | The "clean" blurred frame (mean $\lambda$). | The mean of the distribution. |
| **Measurement** $X$ | The noisy, binned event frame $O_t$. | The actual observed data. |
| **PDF** $f$ | Generative noise model (Normal). | Tool for calculating information. |

---

## 3. Statistical Approximation

Given a sufficiently large intensity $\lambda$, we approximate the Poisson distribution of photon counts as a Normal distribution:
$$X \sim \text{Poisson}(\lambda) \approx \mathcal{N}(\lambda, \lambda)$$

### The Delta Method for Ratios
Since the measurement involves the ratio of two independent Poisson variables $X = I^b_t$ and $Y = I^b_{t-\tau}$, we use the **Delta Method** (first-order Taylor expansion) to approximate the mean and variance:

$$\mathbb{E}\left[\frac{X}{Y}\right] \approx \frac{\mu_X}{\mu_Y}, \quad \text{Var}\left(\frac{X}{Y}\right) \approx \frac{\text{Var}(X)}{\mu_Y^2} + \frac{\mu_X^2 \text{Var}(Y)}{\mu_Y^4}$$

Substituting $\text{Var}(I) = \mathbb{E}[I] = \lambda$:
* **Mean** ($\mu$): $\frac{\lambda_t}{\lambda_{t-\tau}}$
* **Variance** ($\sigma^2$): $\frac{\lambda_t}{\lambda_{t-\tau}^2} + \frac{\lambda_t^2}{\lambda_{t-\tau}^3}$

This leads to the Normal approximation of the measurement ratio:
$$\frac{I^b_t}{I^b_{t-\tau}} \sim \mathcal{N}\left(\frac{\lambda_t}{\lambda_{t - \tau}}, \frac{\lambda_t}{\lambda_{t-\tau}^2} + \frac{\lambda_t^2}{\lambda_{t-\tau}^3}\right)$$

---

## 4. Analytical Derivation (SymPy)

To compute the FIM for the parameter set $\theta$, we apply the Gaussian Fisher Information identity. The code below calculates the coefficients $a, b,$ and $c$ representing the information contributions from $\mu$ (previous intensity) and $\nu$ (current intensity).

```python
import sympy as sp

# 1. Define variables
mu, nu = sp.symbols('mu nu') # mu = intensity at t-tau, nu = intensity at t
X = sp.symbols('X')          # Measurement ratio

# 2. Define mean and variance of the ratio
mean_X = nu / mu
var_X = (nu / mu**2) + (nu**2 / mu**3)

# 3. Partial derivatives w.r.t mu and nu
dm_dmu = sp.diff(mean_X, mu)
dm_dnu = sp.diff(mean_X, nu)
dv_dmu = sp.diff(var_X, mu)
dv_dnu = sp.diff(var_X, nu)

# 4. Fisher Information terms for a Gaussian distribution
# I(theta) = (1/sigma^2) * (d_mu/d_theta)^2 + (1/(2*sigma^4)) * (d_sigma^2/d_theta)^2
a = sp.factor(sp.simplify((1 / var_X) * (dm_dmu**2) + (1 / (2 * var_X**2)) * (dv_dmu**2)))
b = sp.factor(sp.simplify((1 / var_X) * (dm_dmu * dm_dnu) + (1 / (2 * var_X**2)) * (dv_dmu * dv_dnu)))
c = sp.factor(sp.simplify((1 / var_X) * (dm_dnu**2) + (1 / (2 * var_X**2)) * (dv_dnu**2)))
```

### Resulting Information Coefficients
The analytical solutions for the FIM components are:

$$
\begin{aligned}
a &= \frac{2 \mu^{2} \nu + 4 \mu^{2} + 2 \mu \nu^{2} + 12 \mu \nu + 9 \nu^{2}}{2 \mu^{2} (\mu + \nu)^{2}} \\
b &= - \frac{2 \mu^{2} \nu + 2 \mu^{2} + 2 \mu \nu^{2} + 7 \mu \nu + 6 \nu^{2}}{2 \mu \nu (\mu + \nu)^{2}} \\
c &= \frac{2 \mu^{2} \nu + \mu^{2} + 2 \mu \nu^{2} + 4 \mu \nu + 4 \nu^{2}}{2 \nu^{2} (\mu + \nu)^{2}}
\end{aligned}
$$

## 5. Numerical Solution for an Airy disk PSF.

```python
from scipy.special import j1


N_side = 60
pixel_size = 10e-6  # 10 microns
f = 0.05            # 50mm focal length
beta = 0.1          # Background photons
wavelength = 550e-9 
aperture = 0.01

# Create the sensor coordinate grid (u, v)
half_size = (N_side * pixel_size) / 2
u_vals = np.linspace(-half_size, half_size, N_side)
v_vals = np.linspace(-half_size, half_size, N_side)
u_grid, v_grid = np.meshgrid(u_vals, v_vals)
# Numerical Evaluation
a_func = sp.lambdify((mu, nu), a, "numpy")
b_func = sp.lambdify((mu, nu), b, "numpy")
c_func = sp.lambdify((mu, nu), c, "numpy")
# --- 3. Numerical Evaluation Setup ---
# Define test trajectories (Object moving 1mm in x and y)
x_tau, y_tau, z_tau = 0.0, 0.0, 1.0     # Position at t-tau
x_t, y_t, z_t = 0.001, 0.001, 1.0      # Position at t



def airy_disk_psf(u_grid, v_grid, x, y, z):
    # Mapping 3D (x,y,z) to 2D sensor coordinates (u0, v0)
    u0 = f * x / z
    v0 = f * y / z
    
    rho = np.sqrt((u_grid - u0)**2 + (v_grid - v0)**2)
    
    # Scale factor kappa changes with z (defocus approximation)
    w = (1.22 * wavelength * z) / aperture
    k = 3.8317 / w 
    
    val = k * rho
    val = np.where(val == 0, 1e-12, val)
    psf = (2 * j1(val) / val)**2
    
    return psf / np.sum(psf)

def get_full_grad(x, y, z, step=1e-6):
    """Calculates numerical gradients w.r.t x, y, and z."""
    # x gradient
    g_x = (airy_disk_psf(u_grid, v_grid, x + step, y, z) - 
           airy_disk_psf(u_grid, v_grid, x - step, y, z)) / (2 * step)
    # y gradient
    g_y = (airy_disk_psf(u_grid, v_grid, x, y + step, z) - 
           airy_disk_psf(u_grid, v_grid, x, y - step, z)) / (2 * step)
    # z gradient
    g_z = (airy_disk_psf(u_grid, v_grid, x, y, z + step) - 
           airy_disk_psf(u_grid, v_grid, x, y, z - step)) / (2 * step)
    
    return np.stack([g_x, g_y, g_z], axis=-1)


# Generate intensities (mu, nu) and their gradients
mu_grid = airy_disk_psf(u_grid, v_grid, x_tau, y_tau, z_tau) + beta
nu_grid = airy_disk_psf(u_grid, v_grid, x_t, y_t, z_t) + beta

grad_mu = get_full_grad(x_tau, y_tau, z_tau)  # Shape (60, 60, 3)
grad_nu = get_full_grad(x_t, y_t, z_t)        # Shape (60, 60, 3)

# --- 4. Fisher Information Calculation ---
# Evaluate a, b, c coefficients across the grid
a_grid = a_func(mu_grid, nu_grid)
b_grid = b_func(mu_grid, nu_grid)
c_grid = c_func(mu_grid, nu_grid)

# Vectorized summation using Einstein summation
# 'uv' indices are pixels, 'i'/'j' are the 3D params (x,y,z)
A_total = np.einsum('uv,uvi,uvj->ij', a_grid, grad_mu, grad_mu)
C_total = np.einsum('uv,uvi,uvj->ij', c_grid, grad_nu, grad_nu)
B_total = np.einsum('uv,uvi,uvj->ij', b_grid, grad_mu, grad_nu)

# Assemble 6x6 Matrix
I_total = np.block([
    [A_total,   B_total],
    [B_total.T, C_total]
])

# --- 5. Result: CRLB ---
# Invert to get the variance-covariance matrix
try:
    crlb_matrix = np.linalg.inv(I_total)
    std_devs = np.sqrt(np.diag(crlb_matrix))
    
    print("--- Tracking Precision (Standard Deviation) ---")
    print(f"Time t-tau: x={std_devs[0]:.2e}, y={std_devs[1]:.2e}, z={std_devs[2]:.2e}")
    print(f"Time t:     x={std_devs[3]:.2e}, y={std_devs[4]:.2e}, z={std_devs[5]:.2e}")
except np.linalg.LinAlgError:
    print("Matrix is singular! The PSF might be too symmetric or the signal-to-noise is too low.")



```

## Citations
[1] S. Shah et al., “CodedEvents: Optimal Point-Spread-Function Engineering for 3D-Tracking with Event Cameras”.


