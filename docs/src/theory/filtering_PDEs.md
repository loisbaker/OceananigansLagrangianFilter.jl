# Background: PDEs for Lagrangian filtering
## Introduction 
Lagrangian averaging (or, equivalently, Lagrangian *filtering*) is an effective way to decompose complex multi-scale flows into wave and mean components, as it allows separation based on *intrinsic* frequency of processes, that is, frequency in the frame of the flow. 

Lagrangian means are usually found by seeding virtual particles in numerical simulations and keeping track of their positions. However, various methods for calculating Lagrangian means through the solution of PDEs have recently been developed [kafiabadComputingLagrangianMeans2023,bakerLagrangianFilteringWave2025,minzEfficientLagrangianAveraging2025](@citep), removing the need for particle tracking. These PDE-based methods are flexible, and allow various ways of computing the Lagrangian mean dependent on the use case. Some considerations include:

1. Should the filtering be performed *online* (at simulation time) or *offline* (after simulation time, using saved data)?
1. How important is the filter shape - can we get away with a less strict low-pass?
1. How often do we need to find the Lagrangian mean? Once per averaging interval, or at every time? 

These considerations will determine how the Lagrangian average should be performed, and how OceananigansLagrangianFilter is used. 

Before describing how to find Lagrangian averages with PDEs, we first describe the simpler case of finding Eulerian averages with ODEs. 

## ODEs for Eulerian time averages
We consider some scalar function ``f(t)`` (for now neglecting the spatial dimension), from which we would like to obtain ``\bar{f}``, its low-pass filter. We know ``f`` only at the current simulation time ``t``, and it can only depend on past times by causality, so we must have 
```math
\begin{equation}
    \bar{f}(t) = \int_{-\infty}^t K(t,s) f(s)\,\mathrm{d} s\,,
\end{equation}
```
for some weight function ``K(t,s)``, whose shape (a function of ``s``) could depend on the current time ``t``. 

We look for an evolution equation of ``\bar{f}`` and take the time derivative, finding:
```math
\begin{equation}\label{time_evol1}
    \dot{\bar{f}}(t) = K(t,t) f(t) + \int_{-\infty}^t \frac{\partial}{\partial t}K(t,s) \,\mathrm{d} s\,.
\end{equation}
```
To be able to close this equation so that it only depends on information available at the current time, we need either:
1. ``\frac{\partial}{\partial t}K(t,s) = \alpha(t)K(t,s)``, so that the integral in \eqref{time_evol1} can be expressed in terms of ``\bar{f}``.
1. ``\frac{\partial}{\partial t}K(t,s) = 0``, so that ``K(t,s) \equiv G(s)``, and the weight function ``G(s)`` does not change as ``t`` increases. 

These two options give rise to different schemes. 

### Temporal averaging with moving exponential weight functions

In this case, the most general weight function satisfies ``\frac{\partial}{\partial t}K(t,s) = \alpha(t)K(t,s)``. However, we additionally constrain the weight function by imposing that although the reference time of the weight function changes with the simulation time  ``t``, it's shape should not, so ``K(t,s) = G(t-s)``, and ``\bar{f}`` is a convolution between the impulse response ``G(t)`` and ``f``. This then implies that 
```math
\begin{equation}
    K(t,s) \equiv G(t -s) = \alpha e^{-\alpha (t -s)}
\end{equation}
```
for some constant ``\alpha``. With this special exponential weight function, the normalisation 
```math
\begin{equation}
    \int_{-\infty}^t G(t - s) \, ds = 1
\end{equation}
```
holds at all times so ``\bar{f}`` describes the exponential mean at all times (after some initial spin-up), and is given by
```math
\begin{equation}\label{exponentialODE}
    \dot{\bar{f}}(t) = \alpha(f(t) - \bar{f}(t))\,.
\end{equation}
```
Therefore, when the weight function ``G(t)`` is exponential (or composed of a small number of exponentials, as will be descibed later) equation \eqref{exponentialODE} can be solved alongside the governing equation for ``f(t)`` and the exponential mean found 'on-the-fly'. The methods that are currently implemented in OceananigansLagrangianFilter are based on sums of exponential weight functions.

### Temporal averaging with fixed arbitrary weight functions 
In this case, the weight function ``G(s)`` doesn't shift with the simulation time ``t``, so the low pass variable
```math
\begin{equation}
    \bar{f}(t) = \int_{-\infty}^t G(s) f(s)\,\mathrm{d} s
\end{equation}
```
only truly describes the low-pass filtered field as ``t \rightarrow \infty``. The weight function ``G(s)`` can take any form, but should satisfy the normalisation 
```math
\begin{equation}
    \int_{-\infty}^\infty G(s) \,\mathrm{d} s = 1\,.
\end{equation}
```

In reality, we consider some finite interval ``[t^* - T/2, t^* + T/2]`` for some reference time ``t^*``, and consider
```math
\begin{equation}
    \bar{f}(t,t^*) = \int_{t^* - T/2}^t G(t^* - s) f(s)\,\mathrm{d} s\,,
\end{equation}
```
where ``\bar{f}(t^* + T/2,t^*)`` is our desired output, ``G(t)`` is only non-zero on ``[-T/2,T/2]``, and ``\bar{f}`` satisfies
```math
\begin{equation}
    \dot{\bar{f}}(t,t^*) = G(t^* - t)f(t)\,.
\end{equation}
```
We can think of ``t^*_i = T/2 + iT``, where ``T`` is the averaging interval, and ``i \in \{0,1,2,....\}`` as defining a coarse time, and solve for the mean at each ``t^*_i`` by reinitialising ``\bar{f}`` after each time ``T`` so that ``\bar{f}(iT,t^*_i) = 0``.

Methods for filtering with arbitrary weight functions in this way are described in [bakerLagrangianFilteringWave2025](@citet). These are not yet implemented in OceananigansLagrangianFilter, but raise an issue on our [github](https://github.com/loisbaker/OceananigansLagrangianFilter.jl) if you're interested in using these methods. 

## Online exponential Lagrangian filtering 
Here, we demonstrate the how the (single) exponential Lagrangian mean can be found online. This is the basic idea behind OceananigansLagrangianFilter, and is a simplified version of the exponential Lagrangian averaging described in [minzEfficientLagrangianAveraging2025](@citep). We define (see [Lagrangian averaging](@ref "Lagrangian averaging") for general definitions):
```math
\begin{equation}\label{singleexpfstar}
    f^*(\vb*{\varphi}(\vb*{a},t),t) = \int_{-\infty}^t \alpha e^{-\alpha(t-s)}f(\vb*{\varphi}(\vb*{a},s),s)\,\mathrm{d} s\,,
\end{equation}
```
for some inverse timescale ``\alpha``. Taking the time derivative of \eqref{singleexpfstar} at fixed ``\vb*{a}``, and using the chain rule, gives
```math
\begin{equation}
    \frac{\partial f^*}{\partial t}(\vb*{\varphi}(\vb*{a},t),t) + \frac{\partial \vb*{\varphi}}{\partial t}(\vb*{a},t)\cdot \nabla f^*(\vb*{\varphi}(\vb*{a},t),t) = \alpha(f(\vb*{\varphi}(\vb*{a},t),t) - f^*(\vb*{\varphi}(\vb*{a},t),t))\,.
\end{equation}
```
Setting ``\vb*{\varphi}(\vb*{a},t) = \vb*{x}``, and noting that ``\frac{\partial \vb*{\varphi}}{\partial t} (\vb*{a},t)= \vb*{u}(\vb*{\varphi}(\vb*{a},t),t)`` by definition of the flow map, we have
```math
\begin{equation}\label{fstareqnsingleexp}
    \frac{\partial f^*}{\partial t}(\vb*{x},t) + \vb*{u} \cdot\nabla f^*(\vb*{x},t) = \alpha(f(\vb*{x},t) - f^*(\vb*{x},t))\,.
\end{equation}
```
This equation can then be solved alongside the dynamical equations of the simulation (which will determine ``f`` and ``\vb*{u}``) to find ``f^*`` at all times (after some suitable spin-up period). 

If we also want to find the generalised Lagrangian mean ``\bar{f}^{\mathrm{L}}`` (see definition in [Lagrangian averaging](@ref "Lagrangian averaging")), we define a map
```math
\begin{equation}\label{Xidefonline}
\vb*{\Xi}(\vb*{\varphi}(\vb*{a},t),t) = \int_{-\infty}^t \alpha e^{-\alpha(t-s)}\vb*{\varphi}(\vb*{a},s)\,\mathrm{d} s\,.
\end{equation}
```
Taking the time derivative, we find (c.f. \eqref{fstareqnsingleexp})
```math
\begin{equation}\label{Xieqnsingleexp}
    \frac{\partial \vb*{\Xi}}{\partial t}(\vb*{x},t) + \vb*{u} \cdot\nabla \vb*{\Xi}(\vb*{x},t) = \alpha(\vb*{x} - \vb*{\Xi}(\vb*{x},t))\,.
\end{equation}
```
Defining a perturbation ``\vb*{\xi}(\vb*{x},t) = \vb*{\Xi}(\vb*{x},t) - \vb*{x}``, we then have
```math
\begin{equation}\label{xieqnsingleexp}
    \frac{\partial \vb*{\xi}}{\partial t}(\vb*{x},t) + \vb*{u} \cdot\nabla \vb*{\xi}(\vb*{x},t) = - \vb*{u} - \alpha \vb*{\Xi}(\vb*{x},t)\,.
\end{equation}
```
After having solved \eqref{fstareqnsingleexp} for ``f^*`` and \eqref{xieqnsingleexp} for ``\vb*{\xi}``, the relation
```math
\begin{equation}
\bar{f}^{\mathrm{L}}(\vb*{\Xi}(\vb*{x},t),t) = f^*(\vb*{x},t)
\end{equation}
```
can be used to recover ``\bar{f}^{\mathrm{L}}`` by interpolation.


## Offline exponential Lagrangian filtering
While the exponential formulation is efficient as it finds the mean at all times, it limits the possible weight functions, and in particular limits us to causal averages, which depend at any time only on past data. For this reason, we also develop an approach to exponential filtering that combines a 'forward' exponential average (as in \eqref{fstareqnsingleexp}) with a 'backward' exponential pass. 

For \eqref{fstareqnsingleexp} to be run `backwards', the filtering must be performed offline and the data therefore saved at wave-resolving resolution. The forward and backward outputs are summed to give a total output with a more desirable effective filter shape, at the expense of the necessity of saving data. 

The goal is to find 
```math
\begin{equation}
    f^*(\vb*{\varphi}(\vb*{a},t), t) = \int_{-\infty}^{\infty} G(t-s) f(\vb*{\varphi}(\vb*{a},s), s) \, \mathrm{d} s\,,
\end{equation}
```
where ``G(t)`` is even such that ``G(t) \equiv G(|t|)``. 

This property means that the weight function will be centred on the reference time ``t``, and the frequency response of this filter will be real. Such filters have *linear phase shift* (or, in this case, *zero phase shift*,* since ``G`` is symmetric about ``t=0``). This means that the phases of frequencies in the pass-band are not modified, in contrast to filters like the single sided exponential used in [Online Lagrangian filtering](@ref "Online Lagrangian filtering equations"). See [Choosing online filters](@ref "Choosing online filters") for more explanation of weight functions.

For a single exponential (as opposed to sums of exponentials, to be introduced in [Online Lagrangian filtering equations](@ref "Online Lagrangian filtering equations") and [Offline Lagrangian filtering equations](@ref "Offline Lagrangian filtering equations")),
```math
\begin{equation}
G(t-s) = \frac{\alpha}{2} e^{-\alpha|t-s|}\,,
\end{equation}
```
i.e. a `double sided exponential'.

On the forward pass (which can be performed either 'online' at the same time as the simulation, or `offline' on saved data - our implementation does this offline since data must be saved for the backward pass anyway) we calculate
```math
\begin{equation}
    f^*_1(\vb*{\varphi}(\vb*{a},t), t) = \int_{-\infty}^{t} G(t-s) f(\vb*{\varphi}(\vb*{a},s), s) \, \mathrm{d} s\,,
\end{equation}
```
and on the backward pass (which must be performed offline) we calculate 
```math
\begin{equation}\label{f2star}
    f^*_2(\vb*{\varphi}(\vb*{a},t), t) = \int_{t}^{\infty} G(t-s) f(\vb*{\varphi}(\vb*{a},s), s) \, \mathrm{d} s\,,
\end{equation}
```
such that ``f^* = f^*_1 + f^*_2``.

The forward equation is identical to the online scheme \eqref{fstareqnsingleexp}, aside from a factor of two to ensure that the normalisation of the weight function still holds:
```math
\begin{equation}\label{f1stareqnsingleexp}
    \frac{\partial f^*_1}{\partial t}(\vb*{x},t) + \vb*{u} \cdot\nabla f^*_1(\vb*{x},t) = \frac{\alpha}{2}(f(\vb*{x},t) - f^*_1(\vb*{x},t))\,.
\end{equation}
```

The backward equation is found by first taking the time derivative of \eqref{f2star} to give
```math
\begin{equation}\label{f2stareqnsingleexp}
    \frac{\partial f^*_2}{\partial t}(\vb*{x},t) + \vb*{u} \cdot\nabla f^*_2(\vb*{x},t) = -\frac{\alpha}{2}(f(\vb*{x},t) - f^*(\vb*{x},t))\,,
\end{equation}
```
then setting ``\tilde{t} = T - t``, where ``T`` is the total simulation time (or end of the desired averaging interval), to give
```math
\begin{equation}\label{f2stareqnsingleexp_timereversed}
    \frac{\partial f^*_2}{\partial \tilde{t}}(\vb*{x},T - \tilde{t}) - \vb*{u} \cdot\nabla f^*_2(\vb*{x},T - \tilde{t}) = \frac{\alpha}{2}(f(\vb*{x},T - \tilde{t}) - f^*(\vb*{x},T - \tilde{t}))\,.
\end{equation}
```
Defining ``\tilde{f}^*_2(\vb*{x},\tilde{t}) = f^*_2(\vb*{x},T-\tilde{t})``, ``\tilde{f}(\vb*{x},\tilde{t}) = f(\vb*{x},T-\tilde{t})``, and ``\tilde{\vb*{u}}(\vb*{x},\tilde{t}) = \vb*{u}(\vb*{x},T-\tilde{t})``, we then have
```math
\begin{equation}\label{f2stareqnsingleexp_tilde}
    \frac{\partial \tilde{f}^*_2}{\partial \tilde{t}}(\vb*{x},\tilde{t}) + \tilde{\vb*{u}} \cdot\nabla \tilde{f}^*_2(\vb*{x},\tilde{t}) = \frac{\alpha}{2}(f(\vb*{x},\tilde{t}) - \tilde{f}^*(\vb*{x},\tilde{t}))\,,
\end{equation}
```
which is equivalent to \eqref{f1stareqnsingleexp} solved backwards with negated velocities. Equations for the forward and backward maps ``\vb*{\Xi}_1`` and ``\vb*{\Xi}_2`` (c.f. \eqref{Xidefonline}) can be found similarly by setting ``f`` to the identity in \eqref{fstareqnsingleexp} and \eqref{f2stareqnsingleexp_tilde}.

Having found ``f_1^*`` on the forward pass, and ``\tilde{f}_2^*`` on the backward pass, we then calculate 
```math
\begin{equation}
f^*(\vb*{x},t) = f_1^*(\vb*{x},t) + \tilde{f}_2^*(\vb*{x},T-t)\,.
\end{equation}
```