var documenterSearchIndex = {"docs":
[{"location":"get_started/#Get-started","page":"Get started","title":"Get started","text":"","category":"section"},{"location":"get_started/","page":"Get started","title":"Get started","text":"QuantumStateBase is a library to construct common quantum states for the study of quantum optics.","category":"page"},{"location":"get_started/#Quantum-state","page":"Get started","title":"Quantum state","text":"","category":"section"},{"location":"get_started/#Fock-state","page":"Get started","title":"Fock state","text":"","category":"section"},{"location":"get_started/","page":"Get started","title":"Get started","text":"Fock state or number state is the common basis to construct quantum space. To construct a Fock state, a simple way is to use the constructor:","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"julia> FockState(1)\nStateVector{ComplexF64}(dim=70, vec=[\n 0.0 + 0.0im\n 1.0 + 0.0im\n     ⋮\n 0.0 + 0.0im\n 0.0 + 0.0im\n])","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"julia> NumberState(1)\nStateVector{ComplexF64}(dim=70, vec=[\n 0.0 + 0.0im\n 1.0 + 0.0im\n     ⋮\n 0.0 + 0.0im\n 0.0 + 0.0im\n])","category":"page"},{"location":"get_started/#Coherent-state","page":"Get started","title":"Coherent state","text":"","category":"section"},{"location":"get_started/","page":"Get started","title":"Get started","text":"In quantum mechanics, coherent state is defined as an eigenstate of annihilation operator.","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"The simple way to construct a coherent state is to use the pre-defined constructure that provides:","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":" alpha rangle = hatD(alpha)  0 rangle","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"julia> CoherentState(α(4.5, π/4))\nStateVector{ComplexF64}(dim=70, vec=[\n  4.006529739317876e-5 + 9.659393542585425e-17im\n 0.0001274869956459149 - 0.00012748699564586142im\n                       ⋮\n -6.121271897928516e-9 - 5.381634533850132e-23im\n -3.288576246443117e-9 + 3.288576246443058e-9im\n])","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"It is also recommend to use the displace! operator to construct coherent state.","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"For example:","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":" psi rangle = hatD(alpha)  1 rangle","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"julia> displace!(FockState(1), α(4.5, π/4))\nStateVector{ComplexF64}(dim=70, vec=[\n -0.00012748699564590476 - 0.00012748699564585426im\n  -0.0007712569748142598 + 4.4892457438879324e-17im\n                         ⋮\n   -4.558281860401069e-8 - 4.558281860401146e-8im\n   -5.084710321874281e-8 - 4.331617601499496e-22im\n])\n","category":"page"},{"location":"get_started/#Squeezed-state","page":"Get started","title":"Squeezed state","text":"","category":"section"},{"location":"get_started/","page":"Get started","title":"Get started","text":"Squeezed state is defined if its electric field strength for some phases has a quantum uncertainty smaller than that of a coherent state.","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"The simple way to construct a squeezed state is to use the pre-defined constructure that provides:","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":" xi rangle = hatS(xi)  0 rangle","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"julia> SqueezedState(ξ(0.5, π/4))\nStateVector{ComplexF64}(dim=70, vec=[\n     0.9417106158316756 + 3.209844670319301e-18im\n                    0.0 + 0.0im\n                        ⋮\n 2.0592206061670102e-27 - 1.413495457663538e-12im\n                    0.0 + 0.0im\n])","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"It is also recommend to use the squeeze! operator to construct coherent state.","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"For example:","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":" psi rangle = hatS(xi)  1 rangle","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"julia> squeeze!(SinglePhotonState(), ξ(0.3, π/8))\nStateVector{ComplexF64}(dim=70, vec=[\n                    0.0 + 0.0im\n     0.9356524786986595 + 8.16994637882625e-19im\n                        ⋮\n                    0.0 + 0.0im\n 1.1365395531228005e-18 - 1.136539553122705e-18im\n])\n","category":"page"},{"location":"get_started/#Thermal-state","page":"Get started","title":"Thermal state","text":"","category":"section"},{"location":"get_started/","page":"Get started","title":"Get started","text":"Thermal state is a mixed state with photon number distribution described by Bose-Einstein distribution.","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"julia> ThermalState(0.5)\nStateMatrix{ComplexF64}(dim=70, 𝛒=[\n 0.6666666666666666 + 0.0im                 0.0 + 0.0im  …                    0.0 + 0.0im\n                0.0 + 0.0im  0.2222222222222222 + 0.0im                       0.0 + 0.0im\n                0.0 + 0.0im                 0.0 + 0.0im                       0.0 + 0.0im\n                0.0 + 0.0im                 0.0 + 0.0im                       0.0 + 0.0im\n                0.0 + 0.0im                 0.0 + 0.0im                       0.0 + 0.0im\n                0.0 + 0.0im                 0.0 + 0.0im  …                    0.0 + 0.0im\n                0.0 + 0.0im                 0.0 + 0.0im                       0.0 + 0.0im\n                0.0 + 0.0im                 0.0 + 0.0im                       0.0 + 0.0im\n                    ⋮                                    ⋱\n                0.0 + 0.0im                 0.0 + 0.0im                       0.0 + 0.0im\n                0.0 + 0.0im                 0.0 + 0.0im                       0.0 + 0.0im\n                0.0 + 0.0im                 0.0 + 0.0im  …                    0.0 + 0.0im\n                0.0 + 0.0im                 0.0 + 0.0im                       0.0 + 0.0im\n                0.0 + 0.0im                 0.0 + 0.0im                       0.0 + 0.0im\n                0.0 + 0.0im                 0.0 + 0.0im                       0.0 + 0.0im\n                0.0 + 0.0im                 0.0 + 0.0im     7.989915113185906e-34 + 0.0im\n])","category":"page"},{"location":"get_started/#Squeezed-thermal-state","page":"Get started","title":"Squeezed thermal state","text":"","category":"section"},{"location":"get_started/","page":"Get started","title":"Get started","text":"Squeezed thermal state with ξ = 0.3 exp(-im * π/4) and n̄ = 0.5 is equivalent to","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"squeeze!(ThermalState(0.5), ξ(0.3, π/4))","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"julia> SqueezedThermalState(ξ(0.3, π/4), 0.5)\nStateMatrix{ComplexF64}(dim=70, 𝛒=[\n      0.6407801270016841 + 8.058693321522632e-20im  …                      0.0 + 0.0im\n                     0.0 + 0.0im                        4.9231584203398716e-33 + 8.397842435607661e-21im\n    -0.08375298488546218 + 0.08375298488546218im                           0.0 + 0.0im\n                     0.0 + 0.0im                         -6.19303323049067e-20 - 6.193033230493655e-20im\n  1.4299632738856697e-18 - 0.026814344131536903im                          0.0 + 0.0im\n                     0.0 + 0.0im                    …    5.067883911008702e-19 - 3.0814879110195774e-33im\n    0.004524627317645765 + 0.004524627317645767im                          0.0 + 0.0im\n                     0.0 + 0.0im                       -1.4496801073658087e-18 + 1.4496801073653357e-18im\n                         ⋮                          ⋱\n                     0.0 + 0.0im                         2.753897349626776e-18 - 2.7538973496267774e-18im\n   4.567447222067013e-20 - 4.640681643033437e-33im                         0.0 + 0.0im\n                     0.0 + 0.0im                    …    8.787055371266764e-34 + 2.405588755307263e-18im\n  -8.341608280406275e-21 + 8.341608280408356e-21im                         0.0 + 0.0im\n                     0.0 + 0.0im                        -8.741973447257788e-19 - 8.741973447257786e-19im\n -2.5492078509163883e-34 - 3.277032527478944e-21im                         0.0 + 0.0im\n                     0.0 + 0.0im                        1.0347686817817904e-18 + 4.513898307157584e-36im\n])","category":"page"},{"location":"get_started/#Wigner-function","page":"Get started","title":"Wigner function","text":"","category":"section"},{"location":"get_started/","page":"Get started","title":"Get started","text":"Wigner function is introduced by Eugene Winger, who describe quantum state using revised classical probability theory.","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"The quasi probability distribution is defined as:","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"W_mn(x p) = frac12pi int_-infty^infty dy  e^-ipyh psi_m^*(x+fracy2) psi_n(x-fracy2)","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"Owing to the fact that the Moyal function is a generalized Wigner function. We can therefore implies that","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"W(x p) = sum_m n rho_m n W_m n(x p)","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"Here, rho is the density matrix of the quantum state, defined as:","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"rho = sum_m n p_m n  m rangle langle n ","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"And, W_m n(x p) is the generalized Wigner function","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"W_m n =  beginarrayrcl\nfrac1pi exp-(x^2 + y^2) (-1)^m  sqrt2^n-m fracmn (x-ip)^n-m L_m^n-m (2x^2 + 2p^2)  n geq m \nfrac1pi exp-(x^2 + y^2) (-1)^n  sqrt2^m-n fracnm (x+ip)^m-n L_n^m-n (2x^2 + 2p^2)  n  m \nendarray","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"julia> state = SqueezedThermalState(ξ(0.5, 3π/2), 0.3);\n\njulia> wf = WignerFunction(-10:0.1:10, -10:0.1:10);\n\njulia> w = wf(state);\n\njulia> heatmap(w.x_range, w.p_range,  w.𝐰_surface')","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"<img src=\"assets/squeezed_thermal_heatmap.png\" width=\"50%\"/>","category":"page"},{"location":"get_started/#Quadrature-probability-density-function","page":"Get started","title":"Quadrature probability density function","text":"","category":"section"},{"location":"get_started/","page":"Get started","title":"Get started","text":"Usually, we describe a quantum state using two non-commuting observables X(position) and P(momentum) in phase space. The joint distribution is also known as Wigner function.","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"In experiments, we measure the E field of the light using the homodyne detector. The phase of the wave are the eigenvalues of the quadrature operator X_theta where X_theta = 0 = X and X_theta = pi2 = P","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"julia> θs = LinRange(0, 2π, 100);\n\njulia> xs = LinRange(-10, 10, 100);\n\njulia> ps = q_pdf(state, θs, xs);\n\njulia> heatmap(θs, xs, ps')","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"\n<img src=\"assets/squeezed_thermal_quad.png\" width=\"50%\"/>","category":"page"},{"location":"get_started/#Quantum-state-sampler","page":"Get started","title":"Quantum state sampler","text":"","category":"section"},{"location":"get_started/","page":"Get started","title":"Get started","text":"Here, we can sample points from quadrature probability density function of the quantum state. The sampler is implemented by special adaptive rejection method.","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"julia> points = rand(state, 4096);\n\njulia> scatter(points[1, :], points[2, :])","category":"page"},{"location":"get_started/","page":"Get started","title":"Get started","text":"<img src=\"assets/squeezed_thermal_sampled.png\" width=\"50%\"/>","category":"page"},{"location":"apis/#Index","page":"APIs","title":"Index","text":"","category":"section"},{"location":"apis/","page":"APIs","title":"APIs","text":"","category":"page"},{"location":"apis/#APIs","page":"APIs","title":"APIs","text":"","category":"section"},{"location":"apis/","page":"APIs","title":"APIs","text":"Modules = [QuantumStateBase]","category":"page"},{"location":"apis/#QuantumStateBase.ComplexVec","page":"APIs","title":"QuantumStateBase.ComplexVec","text":"ComplexVec{T <: Real}(r::T, θ::T)\n\nVector in polar coordinate for complex plane.\n\nv = r e^-itheta\n\n\n\n\n\n","category":"type"},{"location":"apis/#QuantumStateBase.StateMatrix","page":"APIs","title":"QuantumStateBase.StateMatrix","text":"StateMatrix{T <: Number} <: AbstractState\n\nDensity Matrix representation for pure and mixed quantum state. There are various constructures to construct different pure and mixed quantum states.\n\n\n\n\n\n","category":"type"},{"location":"apis/#QuantumStateBase.StateMatrix-Union{Tuple{StateVector{T}}, Tuple{T}} where T<:Number","page":"APIs","title":"QuantumStateBase.StateMatrix","text":"StateMatrix(state::StateVector{<:Number})\n\nConvert a StateVector to a StateMatrix.\n\nExamples\n\njulia> state = VacuumState();\n\njulia> StateMatrix(state);\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.StateVector","page":"APIs","title":"QuantumStateBase.StateVector","text":"StateVector{T <: Number} <: AbstractState\n\nVector representation for pure quantum state. There are various constructures to construct different pure quantum states.\n\n\n\n\n\n","category":"type"},{"location":"apis/#QuantumStateBase.WignerFunction-Tuple{AbstractRange, AbstractRange}","page":"APIs","title":"QuantumStateBase.WignerFunction","text":"WignerFunction(x_range::AbstractRange, p_range::AbstractRange; dim=DIM)\n\nDeclare a generalized Wigner function with x range and p range\n\nx_range: range of position.\np_range: range of momentum.\ndim: Maximum photon number for truncate, default is 70.\n\nExamples\n\njulia> state = VacuumState();\n\njulia> wf = WignerFunction(-10:0.1:10, -10:0.1:10);\n\njulia> w = wf(state);\n\n\n\n\n\n","category":"method"},{"location":"apis/#Base.copy-Union{Tuple{StateMatrix{T}}, Tuple{T}} where T<:Number","page":"APIs","title":"Base.copy","text":"Base.copy(state::StateMatrix{<:Number})\n\nReturn a new instance of a StateMatrix\n\nExamples\n\njulia> state = VacuumState(rep=StateMatrix);\n\njulia> new_state = copy(state);\n\n\n\n\n\n","category":"method"},{"location":"apis/#Base.copy-Union{Tuple{StateVector{T}}, Tuple{T}} where T<:Number","page":"APIs","title":"Base.copy","text":"Base.copy(state::StateVector{<:Number})\n\nReturn a new instance of a StateVector\n\nExamples\n\njulia> state = VacuumState();\n\njulia> new_state = copy(state);\n\n\n\n\n\n","category":"method"},{"location":"apis/#Base.rand-Tuple{StateMatrix, Integer, Type{IsGaussian}}","page":"APIs","title":"Base.rand","text":"Base.rand(state::AbstractState, n::Integer, ::Type{IsGaussian}; kwargs...)\n\nRandom points sampled from quadrature probability density function of Gaussian state.\n\nstate: Quantum state.\nn: n points.\nIsGaussian: To declare state is a Gaussian state.\nkwargs...: see gaussian_state_sampler\n\n\n\n\n\n","category":"method"},{"location":"apis/#Base.rand-Tuple{StateMatrix, Integer}","page":"APIs","title":"Base.rand","text":"Base.rand(state::AbstractState, n::Integer; kwargs...)\n\nRandom points sampled from quadrature probability density function of state.\n\nstate: Quantum state.\nn: n points.\nkwargs...: see state_sampler\n\n\n\n\n\n","category":"method"},{"location":"apis/#Base.rand-Tuple{StateMatrix, Type{IsGaussian}}","page":"APIs","title":"Base.rand","text":"Base.rand(state::AbstractState, ::Type{IsGaussian}; kwargs...)\n\nOne random point sampled from quadrature probability density function of Gaussian state.\n\nstate: Quantum state.\nIsGaussian: To declare state is a Gaussian state.\nkwargs...: see gaussian_state_sampler\n\n\n\n\n\n","category":"method"},{"location":"apis/#Base.rand-Tuple{StateMatrix}","page":"APIs","title":"Base.rand","text":"Base.rand(state::AbstractState, n::Integer; kwargs...)\n\nOne random point sampled from quadrature probability density function of state.\n\nstate: Quantum state.\nn: n points.\nkwargs...: see state_sampler\n\n\n\n\n\n","category":"method"},{"location":"apis/#Base.vec-Tuple{StateVector{var\"#s4\"} where var\"#s4\"<:Number}","page":"APIs","title":"Base.vec","text":"Base.vec(state::StateVector{<:Number})\n\nTo get the vector of a pure quantum state.\n\nExamples\n\njulia> state = FockState(1);\n\njulia> vec(state)\n70-element Vector{ComplexF64}:\n 0.0 + 0.0im\n 1.0 + 0.0im\n     ⋮\n 0.0 + 0.0im\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.Annihilation-Tuple{}","page":"APIs","title":"QuantumStateBase.Annihilation","text":"Annihilation(; dim=DIM)\n\nAnnihilation operator in matrix representation\n\nhata\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.CoherentState-Tuple{ComplexVec{var\"#s4\"} where var\"#s4\"<:Real}","page":"APIs","title":"QuantumStateBase.CoherentState","text":"CoherentState(α::ComplexVec{<:Real}; dim=DIM, rep=StateVector)\n\nCoherent state is defined as the eigenstate of annihilation operator.\n\nα: Eigenvalue of annihilation operator.\ndim: Maximum photon number for truncate, default is 70.\nrep: In which representation, default is StateVector.\n\nhata  alpha rangle = alpha  alpha rangle\n\nThis constructor will construct  alpha rangle = hatD(alpha)  0 rangle\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.Creation-Tuple{}","page":"APIs","title":"QuantumStateBase.Creation","text":"Creation(; dim=DIM)\n\nCreation operator in matrix representation\n\nhata^dagger\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.Displacement-Tuple{ComplexVec{var\"#s9\"} where var\"#s9\"<:Real}","page":"APIs","title":"QuantumStateBase.Displacement","text":"Displacement(α::ComplexVec{<:Real}; dim=DIM)\n\nDisplacement operator in matrix representation\n\nhatD(alpha) = exp(alpha hata^dagger - alpha^* hata)\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.FockState-Tuple{Any}","page":"APIs","title":"QuantumStateBase.FockState","text":"FockState(n; T=ComplexF64, dim=DIM, rep=StateVector)\n\nConstruct a Fock state in rep representation.\n\nn: Photon number of Fock state.\nT: Numeric data types, default is ComplexF64\ndim: Maximum photon number for truncate, default is 70\nrep: In which representation, default is StateVector.\n\nExamples\n\njulia> state = FockState(1);\n\njulia> state = FockState(1, dim=100);\n\njulia> state = FockState(1, rep=StateMatrix);\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.NumberState-Tuple{Any}","page":"APIs","title":"QuantumStateBase.NumberState","text":"NumberState(n; T=ComplexF64, dim=DIM, rep=StateVector)\n\nExactly the same with FockState.\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.SinglePhotonState-Tuple{}","page":"APIs","title":"QuantumStateBase.SinglePhotonState","text":"SinglePhotonState(; T=ComplexF64, dim=DIM, rep=StateVector)\n\nConstruct a single photon state in rep representation.\n\nT: Numeric data types, default is ComplexF64\ndim: Maximum photon number for truncate, default is 70\nrep: In which representation, default is StateVector.\n\nExamples\n\njulia> state = SinglePhotonState();\n\njulia> state = SinglePhotonState(dim=100);\n\njulia> state = SinglePhotonState(rep=StateMatrix);\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.SqueezedState-Tuple{ComplexVec{var\"#s6\"} where var\"#s6\"<:Real}","page":"APIs","title":"QuantumStateBase.SqueezedState","text":"SqueezedState(ξ::ComplexVec{<:Real}; dim=DIM, rep=StateVector)\n\nSqueezed state is defined if its electric field strength for some phases has a quantum uncertainty smaller than that of a coherent state.\n\nξ: Squeezing factor\ndim: Maximum photon number for truncate, default is 70.\nrep: In which representation, default is StateVector.\n\nThis constructor will construct  xi rangle = hatS(xi)  0 rangle\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.SqueezedThermalState-Tuple{ComplexVec{var\"#s6\"} where var\"#s6\"<:Real, Real}","page":"APIs","title":"QuantumStateBase.SqueezedThermalState","text":"SqueezedThermalState(ξ::ComplexVec{<:Real}, n̄::Real; dim=DIM)\n\nSqueezed state is defined if its electric field strength for some phases has a quantum uncertainty smaller than that of a coherent state.\n\nξ: Squeezing factor\nn̄: Average photon number at temperature T.\ndim: Maximum photon number for truncate, default is 70.\n\nThis constructor will construct rho = hatS(xi) rho_th hatS(xi)^T\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.Squeezing-Tuple{ComplexVec{var\"#s4\"} where var\"#s4\"<:Real}","page":"APIs","title":"QuantumStateBase.Squeezing","text":"Squeezing(ξ::ComplexVec{<:Real}; dim=DIM)\n\nSqueezing operator in matrix representation\n\nhatS(xi) = exp(frac12 (xi^* hata^2 - xi hata^dagger 2))\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.ThermalState-Tuple{Real}","page":"APIs","title":"QuantumStateBase.ThermalState","text":"ThermalState(n̄::Real; dim=DIM)\n\nThermal state is a mixed state with photon number distribution described by Bose-Einstein distribution.\n\nn̄: Average photon number at temperature T.\ndim: Maximum photon number for truncate, default is 70.\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.VacuumState-Tuple{}","page":"APIs","title":"QuantumStateBase.VacuumState","text":"VacuumState(; T=ComplexF64, dim=DIM, rep=StateVector)\n\nConstruct a vacuum state in rep representation.\n\nT: Numeric data types, default is ComplexF64\ndim: Maximum photon number for truncate, default is 70\nrep: In which representation, default is StateVector.\n\nExamples\n\njulia> state = VacuumState();\n\njulia> state = VacuumState(dim=100);\n\njulia> state = VacuumState(rep=StateMatrix);\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.annihilate!-Tuple{StateVector{var\"#s6\"} where var\"#s6\"<:Number}","page":"APIs","title":"QuantumStateBase.annihilate!","text":"annihilate!(state::AbstractState)\n\nApply annihilation operator on the quantum state.\n\nExamples\n\njulia> state = SinglePhotonState();\n\njulia> annihilate!(state);\n\njulia> vec(state) == vec(VacuumState())\ntrue\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.annihilate-Tuple{AbstractState}","page":"APIs","title":"QuantumStateBase.annihilate","text":"annihilate!(state::AbstractState)\n\nApply annihilation operator on the new instance of the quantum state.\n\nExamples\n\njulia> state = SinglePhotonState();\n\njulia> new_state = annihilate(state);\n\njulia> vec(state) == vec(SinglePhotonState())\ntrue\n\njulia> vec(new_state) == vec(VacuumState())\ntrue\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.create!-Tuple{StateVector{var\"#s6\"} where var\"#s6\"<:Number}","page":"APIs","title":"QuantumStateBase.create!","text":"create!(state::AbstractState)\n\nApply creation operator on the quantum state.\n\nExamples\n\njulia> state = VacuumState();\n\njulia> create!(state);\n\njulia> vec(state) == vec(SinglePhotonState())\ntrue\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.create-Tuple{AbstractState}","page":"APIs","title":"QuantumStateBase.create","text":"create(state::AbstractState)\n\nApply creation operator on the new instance of the quantum state.\n\nExamples\n\njulia> state = VacuumState();\n\njulia> new_state = create(state);\n\njulia> vec(state) == vec(VacuumState())\ntrue\n\njulia> vec(new_state) == vec(SinglePhotonState())\ntrue\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.displace!-Tuple{StateVector{var\"#s7\"} where var\"#s7\"<:Number, ComplexVec{var\"#s6\"} where var\"#s6\"<:Real}","page":"APIs","title":"QuantumStateBase.displace!","text":"displace!(state::AbstractState)\n\nApply displacement operator on the quantum state.\n\nExamples\n\njulia> state = VacuumState();\n\njulia> displace!(state, α(5., π/4));\n\njulia> vec(state) == vec(CoherentState(α(5., π/4)))\ntrue\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.gaussian_state_sampler-Tuple{AbstractState, Integer}","page":"APIs","title":"QuantumStateBase.gaussian_state_sampler","text":"gaussian_state_sampler(state::AbstractState, n::Integer; bias_phase=0.)\n\nRandom points sampled from quadrature probability density function of Gaussian state.\n\nstate: Quantum state.\nn: N points.\nbias_phase: The offset of the θ coordinate\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.laguerre-Tuple{Integer, Integer}","page":"APIs","title":"QuantumStateBase.laguerre","text":"laguerre(n::Integer, α::Integer)(x::Real)\n\nGeneralized Laguerre polynomials\n\nExamples\n\njulia> L = laguerre(1, 3);\n\njulia> L(5)\n-1\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.laguerre-Union{Tuple{T}, Tuple{Integer, Integer, T}} where T<:Real","page":"APIs","title":"QuantumStateBase.laguerre","text":"laguerre(n::Integer, α::Integer, x::Real)\n\nGeneralized Laguerre polynomials\n\nExamples\n\njulia> laguerre(0, 1, 5)\n1\n\njulia> laguerre(1, 3, 5)\n-1\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.purity-Tuple{StateMatrix{var\"#s1\"} where var\"#s1\"<:Number}","page":"APIs","title":"QuantumStateBase.purity","text":"purity(state::StateMatrix{<:Number})\n\nCalculate purity for quantum state in density matrix representation.\n\nExamples\n\njulia> state = VacuumState(rep=StateMatrix);\n\njulia> purity(state)\n1.0\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.purity-Tuple{StateVector{var\"#s1\"} where var\"#s1\"<:Number}","page":"APIs","title":"QuantumStateBase.purity","text":"purity(state::StateVector{<:Number})\n\nCalculate purity for quantum state in vector representation.\n\nExamples\n\njulia> state = VacuumState();\n\njulia> purity(state)\n1.0\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.q_pdf-Tuple{AbstractState, Any, Any}","page":"APIs","title":"QuantumStateBase.q_pdf","text":"q_pdf(state::AbstractState, θ::Real, x::Real; T=Float64)\n\nQuadrature prabability at points (θs, xs)\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.q_pdf-Tuple{AbstractState, Real, Real}","page":"APIs","title":"QuantumStateBase.q_pdf","text":"q_pdf(state::AbstractState, θ::Real, x::Real; T=Float64)\n\nQuadrature prabability at point (θ, x)\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.squeeze!-Tuple{StateVector{var\"#s7\"} where var\"#s7\"<:Number, ComplexVec{var\"#s6\"} where var\"#s6\"<:Real}","page":"APIs","title":"QuantumStateBase.squeeze!","text":"squeeze!(state::AbstractState)\n\nApply squeezing operator on the quantum state.\n\nExamples\n\njulia> state = VacuumState();\n\njulia> squeeze!(state, α(0.5, π/4));\n\njulia> vec(state) == vec(SqueezedState(ξ(0.5, π/4)))\ntrue\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.state_sampler-Tuple{AbstractState, Integer}","page":"APIs","title":"QuantumStateBase.state_sampler","text":"state_sampler(\n    state::AbstractState, n::Integer;\n    warm_up_n=128, batch_size=64, c=0.9, θ_range=(0., 2π), x_range=(-10., 10.),\n    show_log=false\n)\n\nRandom points sampled from quadrature probability density function of Gaussian state.\n\nstate: Quantum state.\nn: N points.\nwarm_up_n: N points sampled from uniform random and accepted by rejection method.\nbatch_size: Adapt g for every batch_size points.\nθ_range: Sampling range of θ.\nx_range: Sampling range of x.\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.α-Union{Tuple{T}, Tuple{T, T}} where T","page":"APIs","title":"QuantumStateBase.α","text":"α(r::Real, θ::Real)\n\nEigenvalue of annihilation operator.\n\nhata  alpha rangle = alpha  alpha rangle\n\nExamples\n\njulia> α(1.5, π/4)\nComplexVec{Float64}(1.5exp(-0.7853981633974483im))\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.ξ","page":"APIs","title":"QuantumStateBase.ξ","text":"ξ(r::Real, θ::Real)\n\nExamples\n\njulia> ξ(1.5, π/4)\nComplexVec{Float64}(1.5exp(-0.7853981633974483im))\n\n\n\n\n\n","category":"function"},{"location":"apis/#QuantumStateBase.𝛒-Tuple{StateMatrix{var\"#s1\"} where var\"#s1\"<:Number}","page":"APIs","title":"QuantumStateBase.𝛒","text":"𝛒(state::StateMatrix{<:Number})\n\nTo get the density matrix of a pure quantum state.\n\nExamples\n\njulia> state = FockState(1);\n\njulia> 𝛒(state)\n70×70 Matrix{ComplexF64}:\n 0.0+0.0im  0.0+0.0im  0.0+0.0im  …  0.0+0.0im  0.0+0.0im  0.0+0.0im\n 0.0+0.0im  1.0+0.0im  0.0+0.0im     0.0+0.0im  0.0+0.0im  0.0+0.0im\n    ⋮                             ⋱\n 0.0+0.0im  0.0+0.0im  0.0+0.0im     0.0+0.0im  0.0+0.0im  0.0+0.0im\n\n\n\n\n\n","category":"method"},{"location":"apis/#QuantumStateBase.𝛒-Tuple{StateVector{var\"#s1\"} where var\"#s1\"<:Number}","page":"APIs","title":"QuantumStateBase.𝛒","text":"𝛒(state::StateVector{<:Number})\n\nTo get the density matrix of a pure quantum state.\n\nExamples\n\njulia> state = FockState(1);\n\njulia> 𝛒(state)\n70×70 Matrix{ComplexF64}:\n 0.0+0.0im  0.0+0.0im  0.0+0.0im  …  0.0+0.0im  0.0+0.0im  0.0+0.0im\n 0.0+0.0im  1.0+0.0im  0.0+0.0im     0.0+0.0im  0.0+0.0im  0.0+0.0im\n    ⋮                             ⋱\n 0.0+0.0im  0.0+0.0im  0.0+0.0im     0.0+0.0im  0.0+0.0im  0.0+0.0im\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = QuantumStateBase","category":"page"},{"location":"#QuantumStateBase","page":"Home","title":"QuantumStateBase","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for QuantumStateBase.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package can be installed with the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run:","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add QuantumStateBase","category":"page"},{"location":"#Quick-start","page":"Home","title":"Quick start","text":"","category":"section"},{"location":"#Construct-a-squeezed-thermal-state-and-plot-the-Wigner-function","page":"Home","title":"Construct a squeezed thermal state and plot the Wigner function","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"julia> using QuantumStateBase, Plots\n\njulia> state = SqueezedThermalState(ξ(0.5, 3π/2), 0.3);\n\njulia> wf = WignerFunction(-10:0.1:10, -10:0.1:10);\n\njulia> w = wf(state);\n\njulia> heatmap(w.x_range, w.p_range,  w.𝐰_surface')","category":"page"},{"location":"","page":"Home","title":"Home","text":"<img src=\"assets/squeezed_thermal_heatmap.png\" width=\"50%\"/>","category":"page"},{"location":"#Plot-quadrature-probability-density-function-of-the-state","page":"Home","title":"Plot quadrature probability density function of the state","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"julia> θs = LinRange(0, 2π, 100);\n\njulia> xs = LinRange(-10, 10, 100);\n\njulia> ps = q_pdf(state, θs, xs);\n\njulia> heatmap(θs, xs, ps')","category":"page"},{"location":"","page":"Home","title":"Home","text":"\n<img src=\"assets/squeezed_thermal_quad.png\" width=\"50%\"/>","category":"page"},{"location":"#Sample-points-from-quadrature-PDF-of-the-state","page":"Home","title":"Sample points from quadrature PDF of the state","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"julia> points = rand(state, 4096);\n\njulia> scatter(points[1, :], points[2, :])","category":"page"},{"location":"","page":"Home","title":"Home","text":"<img src=\"assets/squeezed_thermal_sampled.png\" width=\"50%\"/>","category":"page"}]
}