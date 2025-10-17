# Eulerian and Lagrangian averaging
## Eulerian averaging
The Eulerian weighted temporal mean of some scalar ``f`` is defined as
```math
\begin{equation}
\label{fbarEdef}
    \bar{f}^{\mathrm{E}}(\vb*{x},t) = \int_{-\infty}^{\infty} G(t-s)f(\vb*{x},s) \, \mathrm{d} s\,,
\end{equation}
```
where ``G(t)`` is some weight function, also referred to as a *filter kernel* or *impulse response*. The Fourier transform of ``G(t)`` is denoted ``\hat{G}(\omega)``, and is the *frequency response* of the weight function. Here, we will generally want ``G(t)`` to represent a low-pass filter, which retains the low frequencies and removes high frequencies. For example, two possible choices of weight function are 

1. Top-hat impulse response: ``G(t) = \begin{cases} 1/T \,, \hspace{1cm} -T/2 < t < T/2 \\ 0\,, \hspace{1cm} \mathrm{otherwise}\end{cases}\,,``
1. Top-hat frequency response:  ``\hat{G}(\omega) = \begin{cases} 1 \,, \hspace{1cm} -\omega_c < \omega < \omega_c \\ 0\,, \hspace{1cm} \mathrm{otherwise}\end{cases}\,.``


The Eulerian mean ``\bar{f}^{\mathrm{E}}(\vb*{x},t)`` is the field found by taking an average in time at a fixed spatial location ``\vb*{x}``. 

## Lagrangian averaging
In contrast, the Lagrangian mean finds the temporal average whilst moving with the flow on an (imaginary) fluid particle. We define the Lagrangian weighted temporal mean as 
```math
\begin{equation}
\label{fbardef}
    \bar{f}^{\mathrm{L}}(\bar{\vb*{\varphi}}(\vb*{a},t),t) = \int_{-\infty}^{\infty} G(t-s)f(\vb*{\varphi}(\vb*{a},s),s) \, \mathrm{d} s\,,
\end{equation}
```
where the flow map  ``\vb*{\varphi}(\vb*{a},t)`` is the position at time ``t`` of a particle with label ``\vb*{a}`` (which could be the initial position of the particle such that ``\vb*{\varphi}(\vb*{a},0) = \vb*{a}``). The mean flow map ``\bar{\vb*{\varphi}}`` is defined by

```math
\begin{equation}\label{phibardef}
\bar{\vb*{\varphi}}(\vb*{a} ,t) = \int_{-\infty}^{\infty} G(t-s)\vb*{\varphi}(\vb*{a},s)\,\mathrm{d} s\,.
\end{equation}
```

The definition \eqref{fbardef} ensures that ``\bar{f}^{\mathrm{L}}`` is the true generalised Lagrangian mean, in that (for strict band-pass filters) applying the same averaging procedure to the mean flow itself leaves it unchanged \citep{bakerLagrangianFilteringWave2025}. However, we also define an alternative Lagrangian mean, which is a rearrangement of ``\bar{f}^{\mathrm{L}}``:
```math
\begin{equation}\label{fstardef}
    f^*(\vb*{\varphi}(\vb*{a},t),t) = \int_{-\infty}^{\infty} G(t-s)f(\vb*{\varphi}(\vb*{a},s),s) \, \mathrm{d} s\,.
\end{equation}
```

While ``\bar{f}^{\mathrm{L}}(\vb*{x},t)`` describes the average along a particle trajectory whose mean position is ``\vb*{x}``, ``f^*(\vb*{x},t)`` defines the average along a particle trajectory whose position is ``\vb*{x}`` at time ``t``. It is often more desirable, or more convenient, to find ``f^*`` instead. If ``\bar{f}^{\mathrm{L}}`` is also needed, it can be found by a rearrangement of ``f^*`` using a map
```math
\begin{equation}
    \vb*{\Xi}(\vb*{\varphi}(\vb*{a},t),t) = \bar{\vb*{\varphi}}(\vb*{a},t)\,.
\end{equation}
```
