# Offline Filtering
We consider even filter kernels composed of sums of exponentials of the absolute value of ``t``:
```math
\begin{equation}
    G(t) = \sum_{n=1}^{N/2} e^{-c_n |t|} \left( a_n \cos(d_n |t|) + b_n \sin(d_n |t|) \right)\,.
\end{equation}
```
We define a set of ``N`` weight functions, for even ``N``. For ``k = 1,...,N/2`` we have
```math
\begin{align}
    G_{Ck}(t) &=e^{-c_k|t|}\cos d_k |t|\,, \\
    G_{Sk}(t) &=e^{-c_k|t|}\sin d_k |t|\,, \\
\end{align}
```
For ``t>0``, we have
```math
\begin{align} \label{forward_G_derivs}
    G_{Ck}'(t) &= - c_kG_{Ck}(t) - d_k G_{Sk}(t) \\
    G_{Sk}'(t) &= - c_kG_{Sk}(t) + d_k G_{Ck}(t) \\
\end{align}
```
and for ``t<0``
```math
\begin{align} \label{backward_G_derivs}
    G_{Ck}'(t) &= -G_{Ck}'(-t) = c_kG_{Ck}(t) + d_k G_{Sk}(t) \\
    G_{Sk}'(t) &= -G_{Sk}'(-t) = c_kG_{Sk}(t) - d_k G_{Ck}(t) \\
\end{align}
```
We then define a corresponding set of ``N`` filtered scalars
```math
\begin{align}\label{forwardgdef}
    g_{Ck}(\vb*{\varphi}(\vb*{a},t),t) &= \int_{-\infty}^t G_{Ck}(t-s)f(\vb*{\varphi}(\vb*{a},s),s)\, \mathrm{d} s\\
    g_{Sk}(\vb*{\varphi}(\vb*{a},t),t) &= \int_{-\infty}^t G_{Sk}(t-s)f(\vb*{\varphi}(\vb*{a},s),s)\, \mathrm{d} s\\
\end{align}
```
so that 
```math
\begin{equation}
    f_1^*(\vb*{x},t) = \sum_{n=1}^{N/2} a_n g_{Cn}(\vb*{x},t) + b_n g_{Sn}(\vb*{x},t)\,.
\end{equation}
```
We derive PDEs for the ``g_{Ck}`` and ``g_{Sk}`` by taking the time derivative of \eqref{forwardgdef}:
```math
\begin{align}
\frac{\partial g_{Ck}}{\partial t} + \vb*{u} \cdot \nabla g_{Ck} &= f - c_k g_{Ck} - d_k g_{Sk} \\
\frac{\partial g_{Sk}}{\partial t} + \vb*{u} \cdot \nabla g_{Sk} &=  - c_k g_{Sk} + d_k g_{Ck}
\end{align}
```
This system of equations can be solved with initial conditions ``g_{Ck}(\vb*{x},0) = g_{Sk}(\vb*{x},0)=0``.

If we want to find ``\bar{f}^{\mathrm{L}}``, we define map functions with which to interpolate ``f^*`` to ``\bar{f}^{\mathrm{L}}`` after the simulation. We define 
```math
\begin{align}
    \vb*{\Xi}_{Ck}(\vb*{\varphi}(\vb*{a},t),t) &= \int_{-\infty}^t G_{Ck}(t-s)\vb*{\varphi}(\vb*{a},s) \mathrm{d} s\\
    \vb*{\Xi}_{Sk}(\vb*{\varphi}(\vb*{a},t),t) &= \int_{-\infty}^t G_{Sk}(t-s)\vb*{\varphi}(\vb*{a},s) \mathrm{d} s\,,
\end{align}
```
so that
```math
\begin{equation}
\vb*{\Xi}(\vb*{x},t) = \sum_{n=1}^{N/2} a_n \vb*{\Xi}_{Ck}(\vb*{x},t) + b_n \vb*{\Xi}_{Sk}(\vb*{x},t)
\end{equation}
```
and
```math
\begin{align}
\frac{\partial \vb*{\Xi}_{Ck}}{\partial t} + \vb*{u} \cdot \nabla \vb*{\Xi}_{Ck} &= \vb*{x} - c_k \vb*{\Xi}_{Ck} - d_k \vb*{\Xi}_{Sk} \\
\frac{\partial \vb*{\Xi}_{Sk}}{\partial t} + \vb*{u} \cdot \nabla \vb*{\Xi}_{Sk} &=  - c_k \vb*{\Xi}_{Sk} + d_k \vb*{\Xi}_{Ck}\,.
\end{align}
```
To define perturbation equations, we set:
```math
\begin{align}
    \vb*{\Xi}_{Ck} &= \vb*{\xi}_{Ck} + \frac{c_k}{c_k^2 + d_k^2}\vb*{x} \\
    \vb*{\Xi}_{Sk} &= \vb*{\xi}_{Sk} + \frac{d_k}{c_k^2 + d_k^2}\vb*{x}\,,
\end{align}
```
where the coefficients of ``\vb*{x}`` are needed because each of the filters ``G_{Ck}`` and ``G_{Sk}`` are not individually normalised over the interval ``[-\infty,0]``.
The perturbation map equations are then given by
```math
\begin{align}
\frac{\partial \vb*{\xi}_{Ck}}{\partial t} + \vb*{u} \cdot \nabla \vb*{\xi}_{Ck} &= -\frac{c_k}{c_k^2 + d_k^2}\vb*{u} - c_k \vb*{\xi}_{Ck} - d_k \vb*{\xi}_{Sk} \\
\frac{\partial \vb*{\xi}_{Sk}}{\partial t} + \vb*{u} \cdot \nabla \vb*{\xi}_{Sk} &= -\frac{d_k}{c_k^2 + d_k^2}\vb*{u} - c_k \vb*{\xi}_{Sk} + d_k \vb*{\xi}_{Ck}\,.
\end{align}
```
Backward-pass equations of the same form are solved by time-reversing the velocity and field data, and changing the sign of the velocity. The final filtered field is then reconstructed by summing the forwards and backwards pass outputs at each time. 

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
