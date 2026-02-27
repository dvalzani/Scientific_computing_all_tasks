### Task 06

 ```py
#!/usr/bin/env python3
import numpy as np
from scipy.fft import fft2, ifft2, rfft2, irfft2

# =============================
#   1) Create matrix A
# =============================
N = 1000
rng = np.random.default_rng(123)
A = rng.normal(loc=1.0, scale=1.0, size=(N, N))

# =============================
#   2) C2C FFT + inverse FFT
# =============================
C = fft2(A)
A_rec_c2c = np.real(ifft2(C))

# Errors
abs_err = np.abs(A_rec_c2c - A)
rel_err = abs_err / np.abs(A)

# RMS errors
mean_abs_rms    = np.sqrt(np.mean(abs_err**2))
median_abs_rms  = np.sqrt(np.median(abs_err**2))
mean_rel_rms    = np.sqrt(np.mean(rel_err**2))
median_rel_rms  = np.sqrt(np.median(rel_err**2))

print("\n=== C2C Reconstruction ===")
print("Mean ABS RMS   =", mean_abs_rms)
print("Median ABS RMS =", median_abs_rms)
print("Mean REL RMS   =", mean_rel_rms)
print("Median REL RMS =", median_rel_rms)

# =============================
#   3) R2C FFT + inverse (C2R)
# =============================
R = rfft2(A)
A_rec_r2c = irfft2(R, s=A.shape)

abs_err_r = np.abs(A_rec_r2c - A)
rel_err_r = abs_err_r / np.abs(A)

mean_abs_rms_r    = np.sqrt(np.mean(abs_err_r**2))
median_abs_rms_r  = np.sqrt(np.median(abs_err_r**2))
mean_rel_rms_r    = np.sqrt(np.mean(rel_err_r**2))
median_rel_rms_r  = np.sqrt(np.median(rel_err_r**2))

print("\n=== R2C Reconstruction ===")
print("Mean ABS RMS   =", mean_abs_rms_r)
print("Median ABS RMS =", median_abs_rms_r)
print("Mean REL RMS   =", mean_rel_rms_r)
print("Median REL RMS =", median_rel_rms_r)

# =============================
#   4) DC Components
# =============================
print("\nC[0,0] =", C[0,0])
print("R[0,0] =", R[0,0])
print("Sum(A) =", np.sum(A))

# Meaning:
# C[0,0] = R[0,0] = total sum of all entries in A

```
## The Output

```
=== C2C Reconstruction ===
Mean ABS RMS   = 5.604386163014069e-16
Median ABS RMS = 4.440892098500626e-16
Mean REL RMS   = 5.355488650613973e-13
Median REL RMS = 3.688878544893772e-16

=== R2C Reconstruction ===
Mean ABS RMS   = 5.992411392862008e-16
Median ABS RMS = 4.440892098500626e-16
Mean REL RMS   = 5.461079711699702e-13
Median REL RMS = 3.930980738975375e-16

C[0,0] = (1000545.7520460581-0j)
R[0,0] = (1000545.7520460581+0j)
Sum(A) = 1000545.7520460584

```


 ```
