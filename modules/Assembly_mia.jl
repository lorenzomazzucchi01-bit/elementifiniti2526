# Author: Ivan Bioli (https://github.com/IvanBioli)
# Inspired by code written by Jochen Hinz (https://github.com/JochenHinz) for MATH-451 @ EPFL

using Memoize
using SparseArrays

"""
    initialize_assembly!(mesh::Mesh)

Initialize the assembly process for the given mesh by computing the necessary geometric quantities.

# Arguments
- `mesh::Mesh`: The mesh object for which the assembly is initialized.
"""
function initialize_assembly!(mesh::Mesh)
    get_Bk!(mesh)
    get_detBk!(mesh)
    get_invBk!(mesh)
end

########################### GLOBAL ASSEMBLER ########################### 
"""
    assemble_global(mesh::Mesh, local_assembler!)

Assemble the global stiffness matrix and force vector for the given mesh using the provided local assembler function.

# Arguments
- `mesh::Mesh`: The mesh object.
- `local_assembler!`: A function that assembles the local stiffness matrix and force vector.

# Returns
- `K::SparseMatrixCSC`: The global stiffness matrix.
- `f::Vector`: The global force vector.
"""
function assemble_global(mesh::Mesh, local_assembler!)

    q = size(mesh.T, 2)
    n = size(mesh.p, 2)
    f = zeros(n)
    I = []
    J = []
    V = Float64[]

    for p = 1:q

        Tk = mesh.T[:, p]
        fe, Ke = local_assembler!(zeros(3,3), zeros(3), mesh, p)

        f[Tk] += fe
        
        for i = 1:3
            for j = 1:3
                push!(I, Tk[i])
                push!(J, Tk[j])
                push!(V, Ke[i, j])
            end
        end

    end

    A = sparse(I, J, V, n, n)

    return A, f

end


"""
    shapef_2DLFE(quadrule::TriQuad)

Compute the shape functions for the Poisson problem.

# Arguments
- `quadrule::TriQuad`: The quadrature rule.

# Returns
- `shapef`: The shape functions evaluated at the quadrature points.
"""
function shapef_2DLFE(quadrule::TriQuad)
    points = quadrule.points
    n = size(points,2)

    Phi = zeros(3, n)
    for i in 1:n
        point = points[:, i]
        Phi[1,i] = 1 - point[1] - point[2]
        Phi[2,i] = point[1]
        Phi[3,i] = point[2]
    end

    return Phi
end

"""
    ∇shapef_2DLFE(quadrule::TriQuad)

Compute the gradients of the shape functions for the Poisson problem.

# Arguments
- `quadrule::TriQuad`: The quadrature rule.

# Returns
- `∇shapef`: The gradients of the shape functions evaluated at the quadrature points.
"""
@memoize function ∇shapef_2DLFE(quadrule::TriQuad)
    q = size(quadrule.points,2)
    grads = [-1.0 1.0 0.0; -1.0 0.0 1.0]
    return repeat(reshape(grads, 2, 3, 1), 1, 1, q)
end

"""
    poisson_assemble_local!(Ke::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f)

Assemble the local stiffness matrix and force vector for the Poisson problem.

# Arguments
- `Ke::Matrix`: The local stiffness matrix to be assembled.
- `fe::Vector`: The local force vector to be assembled.
- `mesh::Mesh`: The mesh object.
- `cell_index::Integer`: The index of the current cell.
- `f`: The source term function.

# Returns
- `Ke`: The assembled local stiffness matrix.
- `fe`: The assembled local force vector.
"""
function poisson_assemble_local!(Ke::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f)
    
    quadrule = Q2_ref
    ak = mesh.ak[:, cell_index]
    Bk = mesh.Bk[:, :, cell_index]

    if mesh.invBk === nothing
        get_invBk!(mesh)
    end

    Bk_inv = mesh.invBk[:, :, cell_index]
    pe = Bk * quadrule.points + repeat(reshape(ak,2,1), 1, 3)
    Phi = shapef_2DLFE(quadrule)
    Grad = ∇shapef_2DLFE(quadrule)
    G_Phi = similar(Grad)
    for k in 1:size(Grad, 3)
        G_Phi[:, :, k] = Bk_inv' * Grad[:, :, k]
    end

    q = size(quadrule.points, 2)
    for p = 1:q
        wp = quadrule.weights[p] * abs( mesh.detBk[cell_index] )

        for i = 1:3
            v = Phi[i, p]
            G_v = G_Phi[:, i, p]
            fe[i] += wp * v * f(pe[:, p])

            for j = 1:3
                G_u = G_Phi[:, j, p]
                Ke[i, j] += wp * dot(G_v, G_u)
            end
        end
    end

    return fe, Ke
end