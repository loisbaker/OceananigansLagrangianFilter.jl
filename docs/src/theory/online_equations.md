# Online Lagrangian filtering equations

This page describes the Lagrangian filtering equations for the 'online' configuration of *OceananigansLagrangianFilter.jl*. They are very similar to the forward pass of the [offline configuration](@ref "Offline Lagrangian filtering equations").

We directly compute the Lagrangian mean of some scalar ``f`` as
```math
\begin{equation}\label{fstardef}
    f^*(\vb*{\varphi}(\vb*{a},t),t) = \int_{-\infty}^t G(t-s)f(\vb*{\varphi}(\vb*{a},s),s) \, \mathrm{d} s\,,
\end{equation}
```
and optionally compute
```math
\begin{equation}\label{Xidefonline}
\vb*{\Xi}(\vb*{\varphi}(\vb*{a},t),t) = \int_{-\infty}^t G(t-s)\vb*{\varphi}(\vb*{a},s)\,\mathrm{d} s\,,
\end{equation}
```
so that the generalised Lagrangian mean ``\bar{f}^{\mathrm{L}}`` (see definition in [Lagrangian averaging](@ref "Lagrangian averaging")) can be recovered by a post-processing interpolation step using
```math
\begin{equation}
\bar{f}^{\mathrm{L}}(\vb*{\Xi}(\vb*{x},t),t) = f^*(\vb*{x},t)\,.
\end{equation}
```

We consider filter kernels composed of sums of ``N`` exponentials, where ``N`` is even.
```math
\begin{equation}\label{onlinekernel}
    G(t) = \begin{cases}
    \sum_{n=1}^{N/2} e^{-c_n |t|} \left( a_n \cos(d_n |t|) + b_n \sin(d_n |t|) \right)\,, \hspace{1cm} t > 0\,,\\
    0\,,\hspace{1cm} t \leq 0\,.
    \end{cases}
\end{equation}
```
We require ``G`` to be normalised such that
```math
\begin{equation}
    \int_{-\infty}^\infty G(s) \, ds = 1\,,
\end{equation}
```
or equivalently, that ``\hat{G}(0) = 1``. This requires that
```math
\begin{equation}
    \sum_{n=1}^{N/2} \frac{a_nc_n + b_n d_n}{c_n^2 + d_n^2} = 1\,.
\end{equation}
```
This normalisation is only strictly required when we define a map that computes the trajectory mean position (`map_to_mean = true` in [`OnlineFilterConfig`](@ref "OnlineFilterConfig")) but we keep the requirement for now.

We define a set of ``N`` weight functions. For ``k = 1,...,N/2`` we have
```math
\begin{align}
    G_{Ck}(t) &= \begin{cases}
    e^{-c_k t}\cos d_k t\,,\hspace{1cm} &t > 0\,, \\
    0\,,\hspace{1cm} &t \leq 0\,,
    \end{cases}\\
    G_{Sk}(t) &= \begin{cases}
    e^{-c_k t}\sin d_k t\,, \hspace{1cm} &t > 0\,,\\
    0\,,\hspace{1cm} &t \leq 0\,.
    \end{cases}
\end{align}
```
For ``t>0``, we have
```math
\begin{align} \label{forward_G_derivs}
    G_{Ck}'(t) &= - c_kG_{Ck}(t) - d_k G_{Sk}(t)\,, \\
    G_{Sk}'(t) &= - c_kG_{Sk}(t) + d_k G_{Ck}(t)\,. \\
\end{align}
```

We then define a corresponding set of ``N`` filtered scalars
```math
\begin{align}\label{forwardgdef}
    g_{Ck}(\vb*{\varphi}(\vb*{a},t),t) &= \int_{-\infty}^t G_{Ck}(t-s)f(\vb*{\varphi}(\vb*{a},s),s)\, \mathrm{d} s\,,\\
    g_{Sk}(\vb*{\varphi}(\vb*{a},t),t) &= \int_{-\infty}^t G_{Sk}(t-s)f(\vb*{\varphi}(\vb*{a},s),s)\, \mathrm{d} s\,,\\
\end{align}
```
so that 
```math
\begin{equation}\label{reconstitutefstar}
    f^*(\vb*{x},t) = \sum_{n=1}^{N/2} a_n g_{Cn}(\vb*{x},t) + b_n g_{Sn}(\vb*{x},t)\,.
\end{equation}
```
We derive PDEs for the ``g_{Ck}`` and ``g_{Sk}`` by taking the time derivative of \eqref{forwardgdef}:
```math
\begin{align}
\frac{\partial g_{Ck}}{\partial t} + \vb*{u} \cdot \nabla g_{Ck} &= f - c_k g_{Ck} - d_k g_{Sk} \label{gCeqn}\\
\frac{\partial g_{Sk}}{\partial t} + \vb*{u} \cdot \nabla g_{Sk} &=  - c_k g_{Sk} + d_k g_{Ck} \label{gSeqn}
\end{align}
```
This system of equations can be solved with initial conditions ``g_{Ck}(\vb*{x},0) = g_{Sk}(\vb*{x},0)=0`` (TODO add more on ICs, spin-up)

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
\frac{\partial \vb*{\xi}_{Ck}}{\partial t} + \vb*{u} \cdot \nabla \vb*{\xi}_{Ck} &= -\frac{c_k}{c_k^2 + d_k^2}\vb*{u} - c_k \vb*{\xi}_{Ck} - d_k \vb*{\xi}_{Sk} \label{xiCeqn}\\
\frac{\partial \vb*{\xi}_{Sk}}{\partial t} + \vb*{u} \cdot \nabla \vb*{\xi}_{Sk} &= -\frac{d_k}{c_k^2 + d_k^2}\vb*{u} - c_k \vb*{\xi}_{Sk} + d_k \vb*{\xi}_{Ck}\,,\label{xiSeqn}
\end{align}
```
and solved with initial conditions ``\vb*{\xi}_{Ck}(\vb*{x},0) = \vb*{\xi}_{Sk}(\vb*{x},0)=0``.

