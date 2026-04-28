---
title: "Wave Optics Cheat Sheet"
date: 2026-04-23
summary: "A cheat sheet for wave optics equations I frequently forget"
image: "/images/math-chalkboard.jpg"
tags: ["Code", "Math", "Derivation"]
---

In wave optics, light propagation through free space is modeled as a linear, shift-invariant (LTI) system, where lenses and phase masks act as spatial modulators.

### 1. Free-Space Propagation (The LTI System View)

Propagation between two parallel planes can be calculated either by spatial convolution (impulse response) or frequency domain multiplication (transfer function).

**Rayleigh-Sommerfeld Diffraction (Spatial Domain)**
Free space acts as a spatial filter. The field at distance $z$ is the 2D convolution of the initial field with the free-space impulse response $h(x,y,z)$.

$$U(x,y,z) = U(x,y,0) * h(x,y,z)$$

$$h(x,y,z) = \frac{z}{j\lambda} \frac{\exp(j k r)}{r^2}$$

* $r = \sqrt{x^2 + y^2 + z^2}$ is the distance between source and observation points.
* $k = \frac{2\pi}{\lambda}$ is the wavenumber.
* **When to use:** Theoretical formulations for exact scalar diffraction without paraxial restrictions. Rarely used directly in numerical computation due to the highly oscillatory integral.

**Angular Spectrum of Plane Waves (Frequency Domain)**

This is the spatial frequency equivalent of Rayleigh-Sommerfeld. The Fourier transform decomposes the wavefront into a set of plane waves traveling at different angles. The transfer function $H(f_X, f_Y)$ applies a phase delay to each spatial frequency component based on the propagation distance.

$$U(x,y,z) = \mathcal{F}^{-1}\left\{ \mathcal{F}\{U(x,y,0)\} \cdot H(f_X, f_Y) \right\}$$

$$H(f_X, f_Y) = \exp\left(j z \sqrt{k^2 - (2\pi f_X)^2 - (2\pi f_Y)^2}\right)$$

* $f_X, f_Y$ are spatial frequencies (cycles per unit length).
* **When to use:** Numerical simulations of propagation (like in `torchoptics`). It is exact for scalar fields and highly efficient due to the Fast Fourier Transform (FFT).

---

### 2. Paraxial Approximations

When the propagation distance $z$ is large compared to the transverse dimensions (the paraxial approximation), the exact distance $r$ is approximated using a binomial expansion. 

**Fresnel Diffraction (Near-Field)**

The impulse response simplifies to a quadratic phase factor (a spatial chirp). The field is the convolution of the source with this chirp.

$$U(x,y,z) = \frac{\exp(jkz)}{j\lambda z} \iint U(\xi,\eta,0) \exp\left( j\frac{k}{2z} [(x-\xi)^2 + (y-\eta)^2] \right) d\xi d\eta$$

* **When to use:** Analytical calculations in the near-field, or evaluating light immediately after it passes through an aperture or lens.

**Fraunhofer Diffraction (Far-Field)**
As $z$ becomes extremely large, the quadratic phase term across the aperture approaches unity. The diffraction pattern becomes the exact Fourier transform of the aperture, scaled by the propagation distance.

$$U(x,y,z) = \frac{\exp(jkz) \exp\left[j\frac{k}{2z}(x^2+y^2)\right]}{j\lambda z} \mathcal{F}\{U(\xi,\eta,0)\}\Bigg|_{f_X=\frac{x}{\lambda z}, f_Y=\frac{y}{\lambda z}}$$

* **When to use:** Analyzing far-field patterns (e.g., the Airy disk of a circular aperture) or propagation over macroscopic distances without lenses.

---

### 3. Optical Elements and Operations

**Thin Element Approximation**
A thin optical element modulates the complex amplitude of the incoming wave without spatially shifting the rays.

$$U_{out}(x,y) = U_{in}(x,y) \cdot A(x,y) \exp(j\phi(x,y))$$

* $A(x,y)$ represents amplitude attenuation (e.g., an amplitude grating or a binary mask).
* $\phi(x,y)$ represents phase modulation.
* **When to use:** Modeling Spatial Light Modulators (SLMs), coded apertures, or custom refractive optics where physical thickness $h(x,y)$ dictates the phase via $\phi(x,y) = \frac{2\pi}{\lambda} (n_{element} - n_{medium}) h(x,y)$.

**Thin Lens Amplitude Transmittance**
A spherical lens imparts a quadratic phase shift, effectively applying a spatial chirp to the wavefront that either converges or diverges the field.

$$t_{lens}(x,y) = \exp\left(-j \frac{k}{2f} (x^2 + y^2)\right)$$

* $f$ is the focal length. Positive $f$ indicates a converging (convex) lens.
* **When to use:** Representing ideal lenses in an optical path.

**The Fourier Transforming Property of a Lens**

If an object is placed at the front focal plane of a convex lens, the field at the back focal plane is the exact Fourier Transform of the object, without the residual quadratic phase curvature seen in Fraunhofer diffraction.

$$U_f(x_f,y_f) = \frac{1}{j\lambda f} \mathcal{F}\{U_0(x_0,y_0)\}\Bigg|_{f_X=\frac{x_f}{\lambda f}, f_Y=\frac{y_f}{\lambda f}}$$

* **When to use:** Designing 4f correlators, optical filters, or computational optical systems that manipulate light in the frequency domain.

---

### 4. Observables

**Intensity**
Optical sensors operate at timescales far slower than optical frequencies (approx. $10^{14}$ Hz). Therefore, phase is lost, and sensors integrate the absolute square of the complex amplitude over time.

$$I(x,y) = |U(x,y)|^2 = U(x,y) U^*(x,y)$$

* $U^*$ is the complex conjugate.
* **When to use:** The final step of any forward model before applying sensor noise, quantization, or thresholding constraints.

