# Choosing offline filters

## General form

## Exponential

## Butterworth (squared)

## Define your own

The filter coefficients are chosen as
```math
\begin{align}
    a_n &= \frac{\omega_c}{N}\sin{\frac{\pi}{2N}(2n-1)}\,, \\
    b_n &= \frac{\omega_c}{N}\cos{\frac{\pi}{2N}(2n-1)}\,, \\
    c_n &= \omega_c\sin{\frac{\pi}{2N}(2n-1)}\,, \\
    d_n &= \omega_c\cos{\frac{\pi}{2N}(2n-1)}\,.
\end{align}
```
This choice gives a weight function ``G(t)`` with frequency response 
```math
\begin{equation}
    \hat{G}(\omega) = \frac{1}{1 + \left(\omega/\omega_c\right)^{2N}}\,,
\end{equation}
```
which approaches a low-pass cutoff filter as ``N \rightarrow \infty``. 
