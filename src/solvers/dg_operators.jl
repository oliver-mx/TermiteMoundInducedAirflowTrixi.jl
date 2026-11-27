# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

@doc raw"""
    dg_operators(dg, cache)

DG-consistent operators on a 1D `TreeMesh`, all acting on nodal fields stored as
`(node, element)` matrices. Used by the "splitting" ([`GodunovSplittingCallback`](@ref)) and
the "coupled" ([`CoupledSourceTerms`](@ref)) methods.

# Returns
A `NamedTuple` with the fields
- `W`: quadrature weights (LGL weights times the element Jacobian), so that
  ``\int f \, dx \approx`` `sum(W .* f)`,
- `D`: reference derivative matrix and `inverse_jacobian` (physical derivative per element),
- `S`: reference integration matrix with the entries
  ``S_{ij} = \int_{-1}^{\xi_i} \ell_j(s) \, ds`` (see [`integration_matrix`](@ref)),
- `order`: element indices sorted by position (left to right),
- `volume`: length of the domain.
"""
function dg_operators(dg, cache)
    basis = dg.basis
    weights = basis.weights
    inverse_jacobian = cache.elements.inverse_jacobian
    n_elements = Trixi.nelements(dg, cache)
    W = [weights[i] / inverse_jacobian[e] for i in Trixi.eachnode(dg), e = 1:n_elements]
    order = sortperm([cache.elements.node_coordinates[1, 1, e] for e = 1:n_elements])
    return (;
        W,
        D = basis.derivative_matrix,
        S = integration_matrix(basis.nodes),
        inverse_jacobian,
        order,
        volume = sum(W),
    )
end

@doc raw"""
    integration_matrix(nodes)

Return the reference integration matrix `S` with the entries
``S_{ij} = \int_{-1}^{\xi_i} \ell_j(s) \, ds`` for the Lagrange polynomials ``\ell_j`` of the
reference `nodes` ``\xi_i``, computed exactly with Gauss quadrature.
"""
function integration_matrix(nodes)
    n = length(nodes)
    wbary = Trixi.barycentric_weights(nodes)
    gauss_nodes, gauss_weights = Trixi.gauss_nodes_weights(n)
    S = zeros(eltype(nodes), n, n)
    for i = 1:n
        half_length = (nodes[i] + 1) / 2
        for q = 1:n
            s = -1 + half_length * (gauss_nodes[q] + 1)
            ell = Trixi.lagrange_interpolating_polynomials(s, nodes, wbary)
            for j = 1:n
                S[i, j] += half_length * gauss_weights[q] * ell[j]
            end
        end
    end
    return S
end

@doc raw"""
    dg_derivative(f, ops)

Elementwise polynomial derivative ``\partial_x f`` of a nodal field `f` (`(node, element)` matrix).
"""
function dg_derivative(f, ops)
    df = similar(f)
    for e in axes(f, 2)
        df[:, e] = ops.inverse_jacobian[e] .* (ops.D * view(f, :, e))
    end
    return df
end

@doc raw"""
    dg_cumulative_integral(g, ops)

Cumulative integral ``\int_{x_{min}}^x g \, dy`` of a nodal field `g` (`(node, element)` matrix),
integrating the elementwise interpolation polynomial of `g` exactly. The result is zero at
the left end of the domain.
"""
function dg_cumulative_integral(g, ops)
    G = similar(g)
    offset = zero(eltype(g))
    for e in ops.order
        jacobian = inv(ops.inverse_jacobian[e])
        G[:, e] = offset .+ jacobian .* (ops.S * view(g, :, e))
        offset = G[end, e]
    end
    return G
end
end # @muladd
