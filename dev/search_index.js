var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = ImplicitIntegration","category":"page"},{"location":"#ImplicitIntegration","page":"Home","title":"ImplicitIntegration","text":"","category":"section"},{"location":"#overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package provides an integrate function to approximate volume and surface integrals over implicitly defined domains in arbitrary dimensions. More specifically, it allows for the computation of volume integrals of the form","category":"page"},{"location":"","page":"Home","title":"Home","text":"    int_phi(boldsymbolx)  0  cap  U  f(boldsymbolx)  mathrmdV","category":"page"},{"location":"","page":"Home","title":"Home","text":"and surface integrals of the form","category":"page"},{"location":"","page":"Home","title":"Home","text":"    int_phi(boldsymbolx) = 0  cap  U  f(boldsymbolx)  mathrmdS","category":"page"},{"location":"","page":"Home","title":"Home","text":"where phi  mathbbR^d to mathbbR is a level-set function implicitly defining the surface, and U = a_1 b_1 times ldots times a_d b_d is a bounding HyperRectangle.","category":"page"},{"location":"","page":"Home","title":"Home","text":"tip: Algorithm\nThe algorithm implemented is based on [1], and relies on the ability to bound functions and their derivatives over hyperrectangles. Reading the paper is a good idea if you want to understand the details and limitations of the method.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Since the package is not yet registered, you can install it by running","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg; Pkg.add(\"https://github.com/maltezfaria/ImplicitIntegration.jl\");","category":"page"},{"location":"#Basic-usage","page":"Home","title":"Basic usage","text":"","category":"section"},{"location":"#integrate","page":"Home","title":"integrate","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The main function provided by this package is integrate, which computes implicit integrals using an adaptive quadrature method. Here is how you can use it to compute the volume of a sphere:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using ImplicitIntegration, StaticArrays\nϕ = (x) -> x[1]^2 + x[2]^2 + x[3]^2 - 1\nf = (x) -> 1\nlc = (-1.1, -1.1, -1.1)\nhc = (1.1, 1.1, 1.1)\nint_volume  = integrate(f, ϕ, lc, hc)\nnothing # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"The integrate function returns a NamedTuple object containing the value of the integral in val, as well as a logger object containing information about the computation:","category":"page"},{"location":"","page":"Home","title":"Home","text":"println(\"Computed volume:  $(int_volume.val)\")\nprintln(\"Error of volume:  $(int_volume.val - 4π/3)\")\nprintln(int_volume.logger)","category":"page"},{"location":"","page":"Home","title":"Home","text":"To compute a surface integral isntead, simply set the surface keyword argument to true:","category":"page"},{"location":"","page":"Home","title":"Home","text":"int_surface = integrate(f, ϕ, lc, hc; surface = true)\nprintln(\"Computed surface: $(int_surface.val)\")\nprintln(\"Error of surface: $(int_surface.val - 4π)\")\nprintln(int_surface.logger)","category":"page"},{"location":"#quadgen","page":"Home","title":"quadgen","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"For situations where you need to compute the integrals of multiple functions over the same domain, you can use quadgen to generate a quadrature instead. It works similarly to integrate, but returns a Quadrature object instead of the integral value. Here is how to create a surface quadrature for a Cassini oval in 3D:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using ImplicitIntegration, StaticArrays\nusing GLMakie, LinearAlgebra\np1 = SVector(-1.0, 0, 0)\np2 = SVector(1.0, 0, 0)\nb = 1.1\nϕ = x -> (x - p1) ⋅ (x - p1) * (x - p2) ⋅ (x - p2) - b^2\nout = quadgen(ϕ, (-2, -2, -2), (2, 2, 2); order = 5, surface = true)\nxx = yy = zz = range(-1.5, 1.5, length = 200)\nvv = ϕ.(SVector.(Iterators.product(xx, yy, zz)))\nvolume(xx, yy, zz, vv, algorithm = :iso, transparency = true, alpha = 0.4, isovalue = 0)\nscatter!(out.quad.coords, markersize = 2, color = :red)\ncurrent_figure()","category":"page"},{"location":"","page":"Home","title":"Home","text":"See the quadgen docstrings for more information on the available options.","category":"page"},{"location":"#Interface-methods","page":"Home","title":"Interface methods","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"As alluded to in the overview section, the underlying algorithm relies on bounding the level-set function phi as well as its partial derivatives partial_x_i phi in order to find an appropriate height function. By default, ImplicitIntegration.jl uses a combination of IntervalArithmetic.jl and ForwardDiff.jl to perform these task. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"There are at least two situations, however, where it may be beneficial (or even required) to overload the defaults:","category":"page"},{"location":"","page":"Home","title":"Home","text":"ForwardDiff and/or IntervalArithmetic do not work with your function type\nThe bounds provided by the default methods are innacurate and/or slow","category":"page"},{"location":"","page":"Home","title":"Home","text":"In such cases, you can overload the following methods for your levelset ϕ:","category":"page"},{"location":"","page":"Home","title":"Home","text":"ϕ(x::SVector{N,<:Real}): evaluate the function\nϕ(I::SVector{N,<:Interval}): bound the function\nϕ(x::SVector{N,<:ForwardDiff.Dual{Tg,T<:Real,M}}), M <= N: evaluate the Jacobian vector product\nϕ(x::SVector{N,<:ForwardDiff.Dual{Tg,T<:Interval,M}}), M <= N: bound the Jacobian vector product","category":"page"},{"location":"#Common-issues","page":"Home","title":"Common issues","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Below is a list of known issues:","category":"page"},{"location":"#Bibliography","page":"Home","title":"Bibliography","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"R. I. Saye. High-order quadrature methods for implicitly defined surfaces and volumes in hyperrectangles. SIAM Journal on Scientific Computing 37, A993–A1019 (2015).\n\n\n\n","category":"page"},{"location":"docstrings/","page":"Docstrings","title":"Docstrings","text":"CurrentModule = ImplicitIntegration","category":"page"},{"location":"docstrings/#Docstrings","page":"Docstrings","title":"Docstrings","text":"","category":"section"},{"location":"docstrings/","page":"Docstrings","title":"Docstrings","text":"","category":"page"},{"location":"docstrings/","page":"Docstrings","title":"Docstrings","text":"Modules = [ImplicitIntegration]","category":"page"},{"location":"docstrings/#ImplicitIntegration.CellType","page":"Docstrings","title":"ImplicitIntegration.CellType","text":"@enum CellType\n\nEnumeration for different types of cells. Options are full_cell, empty_cell, and partial_cell.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#ImplicitIntegration.Config","page":"Docstrings","title":"ImplicitIntegration.Config","text":"struct Config\n\nThe Config struct represents the configuration for implicit integration, passed to the integrate function to customize the behavior of the algorithm. It contains the following fields:\n\nfind_zero a function with signature (f, a, b, tol) --> x such that f(x) ≈ 0, a ≤ x ≤ b. The tolerance tol is used to specify the absolule tolerance of the zero approximation (e.g. xatol in Roots).\nquad: a function with signature quad(f, a, b, tol) --> (I,E) such that I approximates the integral of f over [a,b] and E is the estimated error. a and b Tuple(s)/SVector(s) specifying the lower and upper bounds of the integration domain, and tol is the desired absolute tolerance.\nmin_vol: a function with signature (tol) --> Float64 used to specify the volume of boxes below which the spatial subdivision stops. low-order method is used to approximate the integral.\nmin_qual: a number between 0 and 1 used to specify the minimum quality factor for a height direction to be considered valid for recursion. The quality factor for a direction k is given by |∂ₖϕ| / |∇ϕ|.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#ImplicitIntegration.GaussLegendre-Tuple{}","page":"Docstrings","title":"ImplicitIntegration.GaussLegendre","text":"GaussLegendre(; order, T = Float64)\n\nConstruct a Gauss-Legendre quadrature rule of given order with nodes and weights of type T, callable as quad1d(f, a, b) where f is the function to integrate and a and b are the bounds of the integration interval.\n\nExamples\n\nquad1d = ImplicitIntegration.GaussLegendre(; order = 20)\nquad1d(cos, 0, 1) ≈ sin(1)\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#ImplicitIntegration.HyperRectangle","page":"Docstrings","title":"ImplicitIntegration.HyperRectangle","text":"struct HyperRectangle{N,T<:AbstractFloat}\n\nA struct representing a hyperrectangle in N-dimensional space.\n\nFields\n\nlc::SVector{N,T}: The lower corner of the hyperrectangle.\nhc::SVector{N,T}: The upper corner of the hyperrectangle.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#ImplicitIntegration.LogInfo","page":"Docstrings","title":"ImplicitIntegration.LogInfo","text":"struct LogInfo\n\nA structure to store logging information for integration processes.\n\nFields\n\nsubdivisions::Vector{Int}: A vector containing the subdivisions per dimension used during the integration process.\nloworder::Int: The number of times the low-order method was used.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#ImplicitIntegration.LogInfo-Tuple{Integer}","page":"Docstrings","title":"ImplicitIntegration.LogInfo","text":"LogInfo(d::Integer)\n\nInitialize a LogInfo object for integrating over ℝᵈ.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#ImplicitIntegration.Quadrature","page":"Docstrings","title":"ImplicitIntegration.Quadrature","text":"struct Quadrature{N,T}\n\nA collection of coords and weights for integration in N dimensions.\n\nQuadratures are the result of quadgen and can be used to integrate functions through integrate(f,::Quadrature).\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#ImplicitIntegration.Segment","page":"Docstrings","title":"ImplicitIntegration.Segment","text":"const Segment{T} = HyperRectangle{1,T}\n\nA one-dimensional hyperrectangle.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#ImplicitIntegration.TensorQuadrature","page":"Docstrings","title":"ImplicitIntegration.TensorQuadrature","text":"TensorQuadrature(quad1d)\n\nGiven a 1D quadrature rule quad1d with signature quad1d(f, a::Number, b::Number), return a tensor product quadrature rule callable through quadnd(f,a::SVector{N}, b::SVector{N}).\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#ImplicitIntegration._find_zeros!-Tuple{Any, Any, ImplicitIntegration.HyperRectangle{1, T} where T, Any, Any, Any}","page":"Docstrings","title":"ImplicitIntegration._find_zeros!","text":"_find_zeros!(roots, f, U::Segment, config, tol)\n\nReturn all zeros of the function f in the Segment U. f should be callable as f(x::SVector{1}).\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#ImplicitIntegration._integrand_eval-Union{Tuple{RET_TYPE}, Tuple{N}, Tuple{Any, Any, Any, ImplicitIntegration.HyperRectangle{N}, Int64, Any, Type{RET_TYPE}, Any, Any}} where {N, RET_TYPE}","page":"Docstrings","title":"ImplicitIntegration._integrand_eval","text":"_integrand_eval(f, phi_vec, s_vec, U, k, config, ::Type{RET_TYPE}, tol, logger)\n\nReturn a function f̃ : ℝᴺ⁻¹ -> ℝ that approximates the one-dimensional integral over I(x̃) ⊂ ℝ of the function t -> f(insert(x̃, k, t)), where the integration domain I is defined as I(x̃) = {t ∈ [a,b] : sᵢ*ϕᵢ(insert(̃x,k,t) ≥ 0 ∀ (ϕᵢ,sᵢ) ∈ V)}.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#ImplicitIntegration.bound-Tuple{Any, Any, Any}","page":"Docstrings","title":"ImplicitIntegration.bound","text":"bound(f, lc, hc) --> (lb, ub)\n\nReturn a lower and upper bound for the function f : U → ℝ valid for all x ∈ U.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#ImplicitIntegration.bound_gradient-Tuple{Any, Any, Any}","page":"Docstrings","title":"ImplicitIntegration.bound_gradient","text":"bound_gradient(f, lc, hc) --> bnds\n\nCompute a lower and upper bound for the gradient of a function f : U → ℝ valid for all x ∈ U in the sense that bnds[i][1] ≤ ∂f/∂xᵢ(x) ≤ bnds[i][2].\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#ImplicitIntegration.bounds-Tuple{ImplicitIntegration.HyperRectangle}","page":"Docstrings","title":"ImplicitIntegration.bounds","text":"bounds(rect::HyperRectangle)\n\nGet the lower and upper bounds of a HyperRectangle.\n\nArguments\n\nrect::HyperRectangle: The HyperRectangle object.\n\nReturns\n\nA tuple (lc, hc) representing the lower and upper bounds of the HyperRectangle.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#ImplicitIntegration.cell_type-NTuple{4, Any}","page":"Docstrings","title":"ImplicitIntegration.cell_type","text":"cell_type(ϕ, s, U, surface)\n\nCompute the CellType of a cell defined by the level-set function ϕ, a sign s, and the box U. If surface is true, then the cell is classified as per a surface integral.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#ImplicitIntegration.gradient-Tuple{Any, Any}","page":"Docstrings","title":"ImplicitIntegration.gradient","text":"gradient(f, x)\n\nCompute the gradient of a function f : ℝᵈ → ℝ at point x ∈ ℝᵈ.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#ImplicitIntegration.integrate-Union{Tuple{T}, Tuple{N}, Tuple{Any, Any, StaticArraysCore.SVector{N, T}, StaticArraysCore.SVector{N, T}}} where {N, T}","page":"Docstrings","title":"ImplicitIntegration.integrate","text":"integrate(f, ϕ, lc, hc; tol=1e-8, surface=false, log = false, config = Config()) -->\n(val, logger)\n\nIntegrate the function f over an implict domain defined by:\n\nΩ = {lc ≤ 𝐱 ≤ hc: ϕ(𝐱) < 0} if surface = false\nΓ = {lc ≤ 𝐱 ≤ hc: ϕ(𝐱) = 0} if surface = true\n\nwhere lc::NTuple and hc::NTuple denote the lower and upper corners of the bounding box.\n\ntol specifies the desired (absolute) tolerance for the approximation.\n\nThe function returns a tuple (val, logger) where val is the approximated value, and logger is a LogInfo object containing information about the integration process.\n\nFor a finer control over the integration process, pass a config object (see Config).\n\nNote that both f and ϕ must be callable with a single argument 𝐱 of type SVector. Furthemore, ϕ is expected to return a real value.\n\nSee also quadgen if you want to generate a quadrature instead of direcly computing the value of the integral.\n\nBy default, ImplicitIntegration uses ForwardDiff to compute gradients and IntervalArithmetic to compute bounds, both of which are needed for the algorithm to work. While these work reasonably well in most cases, you may want to overload\n\nInterface\n\nThe following methods should work for your input function ϕ:\n\nϕ(x::SVector{N,<:Real}) -> T to evaluate the level-set function at x.\nϕ(xI::SVector{N,<:Interval{<:Real}}) -> Interval{<:Real} to evaluate a bound on ϕ on the interval xI.\nϕ(xD::SVector{N,Dual{N,<:Real}}) -> Dual{N,<:Real} to evaluate the level-set function and its gradient at x.\nϕ(xDI::SVector{N,Dual{N,<:Interval{<:Real}}}) -> Dual{N,<:Interval{<:Real}} to evaluate a bound on ϕ and its gradient on the interval xDI.\n\nYou may need to overload the methods above if typeof(ϕ) is not supported by ForwardDiff and/or IntervalArithmetic, or if you have a better/faster implementation.\n\nExamples\n\nTo compute the area of a quarter of a disk of radius 1.0:\n\na, b = (0.0, 0.0), (1.5, 1.5)\nϕ = (x) -> x[1]^2 + x[2]^2 - 1\nf = (x) -> 1.0\nres = integrate(f, ϕ, a, b) # area of quarter of a disk\nres.val ≈ π / 4\n\nTo compute the perimeter of a quarter of a circle of radius 1.0:\n\na, b = (0.0, 0.0), (1.5, 1.5)\nϕ = (x) -> x[1]^2 + x[2]^2 - 1\nf = (x) -> 1.0\nres = integrate(x -> 1.0, ϕ, a, b; surface = true) # perimeter of quarter of a circle\nres.val ≈ 2π / 4\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#ImplicitIntegration.integrate-Union{Tuple{T}, Tuple{N}, Tuple{Any, ImplicitIntegration.Quadrature{N, T}}} where {N, T}","page":"Docstrings","title":"ImplicitIntegration.integrate","text":"integrate(f, Q::Quadrature)\n\nShorthand for ∑ᵢf(xᵢ)wᵢ, where xᵢ and wᵢ are the nodes and weights of the Quadrature.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#ImplicitIntegration.quadgen-Union{Tuple{T}, Tuple{N}, Tuple{Any, StaticArraysCore.SVector{N, T}, StaticArraysCore.SVector{N, T}}} where {N, T}","page":"Docstrings","title":"ImplicitIntegration.quadgen","text":"quadgen(ϕ, lc, hc; order, surface=false, config = nothing)\n\nReturn a Quadrature to integrate a function over an implict domain defined by:\n\nΩ = {lc ≤ 𝐱 ≤ hc: ϕ(𝐱) < 0} if surface = false\nΓ = {lc ≤ 𝐱 ≤ hc: ϕ(𝐱) = 0} if surface = true\n\nwhere lc::NTuple and hc::NTuple denote the lower and upper corners of the bounding box.\n\nThe order parameter specifies the degree of exactness of the quadrature rule; that is, the quadrature rule will integrate exactly polynomials of degree up to order, but not order+1. A GaussLegendre quadrature rule is used.\n\nFor a finer control of the integration process, the user can pass a config object to customize the behavior of various aspects of the algorithm (see Config for more details). In such cases, the order parameter is ignored and config.quad is used for the integration.\n\nNote that ϕ must be callable with a single argument 𝐱 of type SVector. Furthemore, ϕ is expected to return a real value.\n\nExamples\n\nTo compute the area of a quarter of a disk of radius 1.0:\n\na, b = (0.0, 0.0), (1.5, 1.5)\nϕ = (x) -> x[1]^2 + x[2]^2 - 1\nf = (x) -> 1.0\nout = quadgen(ϕ, a, b; order = 20)\nintegrate(f, out.quad) ≈ π / 4 # area of quarter of a disk\n\nTo compute the perimeter of a quarter of a circle of radius 1.0:\n\na, b = (0.0, 0.0), (1.5, 1.5)\nϕ = (x) -> x[1]^2 + x[2]^2 - 1\nf = (x) -> 1.0\nout = quadgen(ϕ, a, b; order = 20, surface = true)\nQ = out.quad\nintegrate(f, Q) ≈ 2π / 4 # perimeter of quarter of a circle\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#ImplicitIntegration.remove_dimension-Tuple{ImplicitIntegration.HyperRectangle, Any}","page":"Docstrings","title":"ImplicitIntegration.remove_dimension","text":"remove_dimension(rect::HyperRectangle, k)\n\nRemove a dimension from a HyperRectangle by deleting the k-th element from the lower and upper corners.\n\nArguments\n\nrect::HyperRectangle: The input hyperrectangle.\nk: The index of the dimension to be removed.\n\nReturns\n\nA new HyperRectangle with the k-th dimension removed.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#ImplicitIntegration.restrict-Tuple{Any, Any, Any}","page":"Docstrings","title":"ImplicitIntegration.restrict","text":" restrict(f, k, v)\n\nGiven a function f : ℝᵈ → ℝ, a value v ∈ ℝ and an integer 1 ≤ k ≤ d, return the function f̃ : ℝᵈ⁻¹ → ℝ defined by restricting f to the value v along dimension d; i.e. f̃(x) = f(x₁, ..., x_{k-1}, v, x_{k}, ..., x_d).\n\nnote: Note\nThe returned type should also implement the interface methods gradient, bound, and bound_gradient.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#ImplicitIntegration.sgn-Tuple{Any, Any, Bool, Any}","page":"Docstrings","title":"ImplicitIntegration.sgn","text":"sgn(m, s, S::Bool, σ)\n\nHelper function to compute the sign of lower and upper restrictions of a level-set function ϕᵢ in a box along a given height direction k. Here m is sign of ∂ₖϕᵢ, which is assume not to change throughout the box since k is assumed to be a height direction, s is the sign of ϕᵢ in the box, S is a flag to indicate whether we are in the special case of a surface integral.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#ImplicitIntegration.split-Tuple{ImplicitIntegration.HyperRectangle, Any}","page":"Docstrings","title":"ImplicitIntegration.split","text":"split(U::HyperRectangle, dir)\n\nSplit a hyperrectangle U along the specified direction dir.\n\nArguments\n\nU::HyperRectangle: The hyperrectangle to be split.\ndir: The direction along which to split the hyperrectangle.\n\nReturns\n\nUₗ: The left half of the split hyperrectangle.\nUᵣ: The right half of the split hyperrectangle.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#ImplicitIntegration.use_heuristic_bounds","page":"Docstrings","title":"ImplicitIntegration.use_heuristic_bounds","text":"use_heuristic_bounds(F::Type, n = 10)\n\nOverload interface methods for type F to use a sample-based heuristic when bounding functions of type F and its gradient.\n\nThe bounds are obtained by sampling the function (or its gradient) on an n × ... × n grid and taking the maximum/minimum of the attained values. This is obviously a heuristic, and may fail in practice.\n\n\n\n\n\n","category":"function"}]
}