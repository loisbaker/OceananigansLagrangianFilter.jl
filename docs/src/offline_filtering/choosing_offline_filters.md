# Choosing offline filters
To be added. 

TLDR; set `N` and `freq_c` in [`OfflineFilterConfig`](@ref "OfflineFilterConfig") to get a filter with squared amplitude of a Butterworth order `N` filter with cutoff frequency `freq_c`.

If `N=1` then a double sided exponential filter is found. Otherwise, `N` should be even.

## General form

## Exponential

## Butterworth (squared)
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

## Define your own


