# Quick start guide

`QuantumStateBase` is a library to construct common quantum states for the study of quantum optics.

## Quantum state

### Fock state

Fock state or number state is the common basis to construct quantum space.
To construct a Fock state, a simple way is to use the constructor:

```julia
julia> FockState(1)
StateVector{ComplexF64}(dim=70, vec=[
 0.0 + 0.0im
 1.0 + 0.0im
     ⋮
 0.0 + 0.0im
 0.0 + 0.0im
])
```
```julia
julia> NumberState(1)
StateVector{ComplexF64}(dim=70, vec=[
 0.0 + 0.0im
 1.0 + 0.0im
     ⋮
 0.0 + 0.0im
 0.0 + 0.0im
])
```

### Coherent state

In quantum mechanics, coherent state is defined as an eigenstate of annihilation operator.

The simple way to construct a coherent state is to use the pre-defined constructure that provides:

```math
| \alpha \rangle = \hat{D}(\alpha) | 0 \rangle
```

```julia
julia> CoherentState(α(4.5, π/4))
StateVector{ComplexF64}(dim=70, vec=[
  4.006529739317876e-5 + 9.659393542585425e-17im
 0.0001274869956459149 - 0.00012748699564586142im
                       ⋮
 -6.121271897928516e-9 - 5.381634533850132e-23im
 -3.288576246443117e-9 + 3.288576246443058e-9im
])
```

It is also recommend to use the `displace!` operator to construct coherent state.

For example:

```math
| \psi \rangle = \hat{D}(\alpha) | 1 \rangle
```

```julia
julia> displace!(FockState(1), α(4.5, π/4))
StateVector{ComplexF64}(dim=70, vec=[
 -0.00012748699564590476 - 0.00012748699564585426im
  -0.0007712569748142598 + 4.4892457438879324e-17im
                         ⋮
   -4.558281860401069e-8 - 4.558281860401146e-8im
   -5.084710321874281e-8 - 4.331617601499496e-22im
])

```

### Squeezed state

Squeezed state is defined if its electric field strength for some phases
has a quantum uncertainty smaller than that of a coherent state.

The simple way to construct a squeezed state is to use the pre-defined constructure that provides:

```math
| \xi \rangle = \hat{S}(\xi) | 0 \rangle
```

```julia
julia> SqueezedState(ξ(0.5, π/4))
StateVector{ComplexF64}(dim=70, vec=[
     0.9417106158316756 + 3.209844670319301e-18im
                    0.0 + 0.0im
                        ⋮
 2.0592206061670102e-27 - 1.413495457663538e-12im
                    0.0 + 0.0im
])
```

It is also recommend to use the `squeeze!` operator to construct coherent state.

For example:

```math
| \psi \rangle = \hat{S}(\xi) | 1 \rangle
```

```julia
julia> squeeze!(SinglePhotonState(), ξ(0.3, π/8))
StateVector{ComplexF64}(dim=70, vec=[
                    0.0 + 0.0im
     0.9356524786986595 + 8.16994637882625e-19im
                        ⋮
                    0.0 + 0.0im
 1.1365395531228005e-18 - 1.136539553122705e-18im
])

```


### Thermal state

Thermal state is a mixed state with photon number distribution described by Bose-Einstein distribution.

```julia
julia> ThermalState(0.5)
StateMatrix{ComplexF64}(dim=70, 𝛒=[
 0.6666666666666666 + 0.0im                 0.0 + 0.0im  …                    0.0 + 0.0im
                0.0 + 0.0im  0.2222222222222222 + 0.0im                       0.0 + 0.0im
                0.0 + 0.0im                 0.0 + 0.0im                       0.0 + 0.0im
                0.0 + 0.0im                 0.0 + 0.0im                       0.0 + 0.0im
                0.0 + 0.0im                 0.0 + 0.0im                       0.0 + 0.0im
                0.0 + 0.0im                 0.0 + 0.0im  …                    0.0 + 0.0im
                0.0 + 0.0im                 0.0 + 0.0im                       0.0 + 0.0im
                0.0 + 0.0im                 0.0 + 0.0im                       0.0 + 0.0im
                    ⋮                                    ⋱
                0.0 + 0.0im                 0.0 + 0.0im                       0.0 + 0.0im
                0.0 + 0.0im                 0.0 + 0.0im                       0.0 + 0.0im
                0.0 + 0.0im                 0.0 + 0.0im  …                    0.0 + 0.0im
                0.0 + 0.0im                 0.0 + 0.0im                       0.0 + 0.0im
                0.0 + 0.0im                 0.0 + 0.0im                       0.0 + 0.0im
                0.0 + 0.0im                 0.0 + 0.0im                       0.0 + 0.0im
                0.0 + 0.0im                 0.0 + 0.0im     7.989915113185906e-34 + 0.0im
])
```

### Squeezed thermal state

Squeezed thermal state with `ξ = 0.3 exp(-im * π/4)` and `n̄ = 0.5` is equivalent to

```julia
squeeze!(ThermalState(0.5), ξ(0.3, π/4))
```

```julia
julia> SqueezedThermalState(ξ(0.3, π/4), 0.5)
StateMatrix{ComplexF64}(dim=70, 𝛒=[
      0.6407801270016841 + 8.058693321522632e-20im  …                      0.0 + 0.0im
                     0.0 + 0.0im                        4.9231584203398716e-33 + 8.397842435607661e-21im
    -0.08375298488546218 + 0.08375298488546218im                           0.0 + 0.0im
                     0.0 + 0.0im                         -6.19303323049067e-20 - 6.193033230493655e-20im
  1.4299632738856697e-18 - 0.026814344131536903im                          0.0 + 0.0im
                     0.0 + 0.0im                    …    5.067883911008702e-19 - 3.0814879110195774e-33im
    0.004524627317645765 + 0.004524627317645767im                          0.0 + 0.0im
                     0.0 + 0.0im                       -1.4496801073658087e-18 + 1.4496801073653357e-18im
                         ⋮                          ⋱
                     0.0 + 0.0im                         2.753897349626776e-18 - 2.7538973496267774e-18im
   4.567447222067013e-20 - 4.640681643033437e-33im                         0.0 + 0.0im
                     0.0 + 0.0im                    …    8.787055371266764e-34 + 2.405588755307263e-18im
  -8.341608280406275e-21 + 8.341608280408356e-21im                         0.0 + 0.0im
                     0.0 + 0.0im                        -8.741973447257788e-19 - 8.741973447257786e-19im
 -2.5492078509163883e-34 - 3.277032527478944e-21im                         0.0 + 0.0im
                     0.0 + 0.0im                        1.0347686817817904e-18 + 4.513898307157584e-36im
])
```
