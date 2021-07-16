```@meta
CurrentModule = QuantumStateBase
```

# QuantumStateBase

Documentation for [QuantumStateBase](https://github.com/foldfelis-QO/QuantumStateBase.jl).

## Installation

The package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia-repl
pkg> add QuantumStateBase
```

## Quick start

### Construct a squeezed thermal state and plot the Wigner function

```julia-repl
julia> using QuantumStateBase, Plots

julia> state = SqueezedThermalState(ξ(0.5, 3π/2), 0.3);

julia> wf = WignerFunction(-10:0.1:10, -10:0.1:10);

julia> w = wf(state);

julia> heatmap(w.x_range, w.p_range,  w.𝐰_surface')
```

```@raw html
<img src="assets/squeezed_thermal_heatmap.png" width="50%"/>
```
