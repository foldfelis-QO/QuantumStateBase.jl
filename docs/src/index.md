```@meta
CurrentModule = QuantumStateBase
```

# QuantumStateBase

Documentation for [QuantumStateBase](https://github.com/foldfelis-QO/QuantumStateBase.jl).

## Installation

The package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
pkg> add QuantumStateBase
```

## Quick start

### Construct a squeezed thermal state and plot the Wigner function

```julia
julia> using QuantumStateBase, Plots

julia> state = SqueezedThermalState(0.5, 3π/2, 0.3, dim=35);

julia> w = wigner(state, LinRange(-3, 3, 101), LinRange(-3, 3, 101));

julia> heatmap(w.x_range, w.p_range,  w.𝐰_surface')
```

```@raw html
<img src="assets/squeezed_thermal_heatmap.png" width="50%"/>
```

## Wigner function

Wigner function is introduced by Eugene Winger, who describe quantum state using revised classical probability theory.

The quasi probability distribution is defined as:

```math
W_{mn}(x, p) = \frac{1}{2\pi} \int_{-\infty}^{\infty} dy \, e^{-ipy/h} \psi_m^*(x+\frac{y}{2}) \psi_n(x-\frac{y}{2})
```

Owing to the fact that the Moyal function is a generalized Wigner function. We can therefore implies that

```math
W(x, p) = \sum_{m, n} \rho_{m, n} W_{m, n}(x, p)
```

Here, ``\rho`` is the density matrix of the quantum state, defined as:

```math
\rho = \sum_{m, n, i} \, p_i \, | n \rangle \langle n | \hat{\rho}_i | m \rangle \langle m |
```
```math
\hat{\rho}_i = | \psi_i \rangle \langle \psi_i |
```
```math
\hat{\rho}_i \, \text{is a density operator of pure state.}
```

And, ``W_{m, n}(x, p)`` is the generalized Wigner function

```math
W_{m, n} = \{ \begin{array}{rcl}
\frac{1}{\pi} exp[-(x^2 + y^2)] (-1)^m  \sqrt{2^{n-m} \frac{m!}{n!}} (x-ip)^{n-m} L_m^{n-m} (2x^2 + 2p^2), \, n \geq m \\
\frac{1}{\pi} exp[-(x^2 + y^2)] (-1)^n  \sqrt{2^{m-n} \frac{n!}{m!}} (x+ip)^{m-n} L_n^{m-n} (2x^2 + 2p^2), \, n < m \\
\end{array}
```
