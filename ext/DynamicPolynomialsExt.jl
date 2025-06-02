module DynamicPolynomialsExt

using DynamicPolynomials
import ImplicitIntegration

function __init__()
    @info "Loading DynamicPolynomials extension for ImplicitIntegration.jl"
end

function ImplicitIntegration.BernsteinPolynomial(p::DynamicPolynomials.Polynomial, lb, ub)
    T     = DynamicPolynomials.coefficienttype(p)
    sz    = mapreduce(m -> DynamicPolynomials.exponents(m), (a, b) -> max.(a, b), DynamicPolynomials.monomials(p)) .+ 1
    coefs = zeros(T, sz...)
    for t in terms(p)
        m = DynamicPolynomials.exponents(t)
        v = DynamicPolynomials.coefficient(t)
        I = m .+ 1
        coefs[Tuple(I)...] = v
    end
    b = ImplicitIntegration._monomial_to_bernstein(coefs, lb, ub)
    return ImplicitIntegration.BernsteinPolynomial(b, lb, ub)
end

function ImplicitIntegration.integrate(
    f,
    p::DynamicPolynomials.Polynomial,
    lb,
    ub,
    args...;
    kwargs...,
)
    b = ImplicitIntegration.BernsteinPolynomial(p, lb, ub)
    return ImplicitIntegration.integrate(f, b, lb, ub, args...; kwargs...)
end

function ImplicitIntegration.quadgen(
    p::DynamicPolynomials.Polynomial,
    lb,
    ub,
    args...;
    kwargs...,
)
    b = ImplicitIntegration.BernsteinPolynomial(p, lb, ub)
    return ImplicitIntegration.quadgen(b, lb, ub, args...; kwargs...)
end

end # module
