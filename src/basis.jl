using FastGaussQuadrature

"""
    lagrange_1d(points, i, x)

Evaluate the Lagrange interpolation polynomial defined over nodal `points` with index `i`
at point `x`.
"""
function lagrange_1d(points, i, x)
    l = 1.0
    x_i = points[i]
    for x_j in points
        if x_j != x_i
            l *= (x - x_j) / (x_i - x_j)
        end
    end
    return l
end

"""
    lagrange_diff(points, i, x)

Evaluate the derivative of the Lagrange interpolation polynomial defined over nodal `points`
with index `i` at point `x`.
"""
function lagrange_diff(points, i, x)
    dl = 0.0
    x_i = points[i]
    for x_j in points
        if x_j != x_i
            p = 1.0
            for x_m in points
                if (x_m != x_i) & (x_m != x_j)
                    p *= (x - x_m) / (x_i - x_m)
                end
            end
            dl += 1.0/(x_i - x_j) * p
        end
    end
    return dl
end

"""
    get_quadpoints(n)

Compute quadrature points and weights for Gaussian quadrature
of order `n`.
The points (and thus the corresponding weights) are normalized to the range
``[0.0, 1.0]``.

Return a tuple of `points, weights`.
"""
function get_quadpoints(n)
    x, w = gausslegendre(n)
    # They are scaled in [-1.0, 1.0]
    # and we need [0.0, 1.0]
    (x .+ 1)./2, w ./ 2
end

"""
    Basis

A standard 1-dimensional basis of `order::Integer`
with `quadpoints::Array{Float64,1}` and
corresponding `quadweights::Array{Float64,1}`
Is the basis (pun intended) for tensor-product
bases.

    Basis(order::Integer, dimensions)

Initialize a basis of `order::Integer` and
dimensions `dimensions`.
"""
struct Basis
    quadpoints::Array{Float64,1}
    quadweights::Array{Float64,1}
    order::Int64
    dimensions::Int64

    function Basis(order, dimensions)
        quadpoints, quadweights = get_quadpoints(order)
        new(quadpoints, quadweights, order, dimensions)
    end
end


"""
    Base.length(basis::Basis)

Return number of points for `basis` in n-dimensions.
"""
#Base.length(basis::Basis) = length(basis.quadpoints)*basis.dimensions
#Base.length(basis::Basis) = length(basis.quadpoints) ^ basis.dimensions
Base.length(basis::Basis) = basis.order ^ basis.dimensions

"""
    Base.size(basis::Basis)

Return number of points for `basis` for each dimensions as tuple.
"""
function Base.size(basis::Basis)
    ntuple(basis.dimensions) do x
       length(basis.quadpoints)
    end
end

"""
    Base.size(basis::Basis, dim)

Return number of points for `basis` for dimensions `dim`.

function Base.size(basis::Basis, dim::Integer)
    length(basis.quadpoints)
end
"""

function Base.size(basis::Basis, dim::Integer)
    if dim > basis.dimensions
        println("Error: Given Dimension exceeds the Dimension of the problem = ",basis.dimensions)
    else
        length(basis.quadpoints)
    end
end

"""
    evaluate_basis(basis::Basis, coeffs, x)

Evaluate the `basis` with coefficients
`coeffs` at point `x`.

# So the input is coeffs , X point, I need to give basis
#multiplication of the basis is used , indicate tensor product
"""
function evaluate_basis(basis::Basis, coeffs, x)
    #coeffs[1] #code given for order 1

    #To collect result
    basisvalue=0

    # To make linear indexes
    # collect(): makes a array of tuples of indexes, same like meshgrid in matlab
    # vec(): converts matrix into vector
    linearindex = vec(collect(Iterators.product(1:basis.order, 1:basis.order)))
    
    for i in 1:length(basis)
        basisvalue = basisvalue + ( coeffs[i] * (  lagrange_1d(basis.quadpoints, linearindex[i][1], x[1]) *  lagrange_1d(basis.quadpoints, linearindex[i][2], x[2]) ) )
    end

    return basisvalue
end

"""
    project_to_reference_basis(fun, basis::Basis, ndofs)

Project the result of the function `fun` to coefficients
of the basis built of a tensor-product of `basis`.
The function `fun(x,y)`  takes in the ``x, y``-coordinates
and returns a vector with size `ndofs`.
The corresponding coefficients are returned.

#The reshape() is an inbuilt function in julia which is used
to return an array with the same data as the specified array,
but with different specified dimension sizes.

# It returns the matrix of shape[length(basis),ndof]

"""
function project_to_reference_basis(fun, basis::Basis, ndofs)
    #reshape(fun(0.5, 0.5), (1,ndofs)) #made for constant polynomial of order 1
    
    # make a matrix of shape [length(basis),ndof]
    coeffarray = zeros( ( length(basis) , ndofs ) )

    # To make linear indexes
    # collect(): makes a array of tuples of indexes, same like meshgrid in matlab
    # vec(): converts matrix into vector

    linearindex = vec(collect(Iterators.product(1:basis.order, 1:basis.order)))

    for i in 1:length(basis)
        coeffarray[i,:] = fun(  basis.quadpoints[linearindex[i][1]]  ,   basis.quadpoints[linearindex[i][2]]  ) # coefficients are nothing but values of function at that point
    end

    return coeffarray # No need to reshaping ,as matrix is made fully and filled it
end


"""
    massmatrix(basis, dimensions)

Return the mass-matrix for a `dimensions`-dimensional
tensor-product basis built up from the 1d-basis `basis`.
"""
function massmatrix(basis, dimensions)
    #ones(1,1)
    mass_matrix = Matrix(1.0I, basis.order^2, basis.order^2)
    ij = 1
    for j = 1:basis.order
        for i = 1:basis.order                                  #not generalising to dimensions>2
            mass_matrix[ij,ij] = basis.quadweights[i] * basis.quadweights[j]                        #eq 4.12 and 4.19
            ij += 1
        end
    end
    mass_matrix

end

"""
    derivativematrix(basis)

Returns the 2-dimensional derivative matrix for `basis`.
Multiplying this with flux-coefficients of shape
`(dimensions * basissize_2d, ndofs)` returns the
coefficients of the corresponding derivative.
"""
function derivativematrix(basis)
    #zeros(1, basis.dimensions * 1)

    # Finding Derivative matrix Try3
    # Finding Dx only, Derivative of the basis function in the direction of x
    Dx = zeros(basis.order^2, basis.order^2)
    Dy = zeros(basis.order^2, basis.order^2)
    linearindex = vec(collect(Iterators.product(1:basis.order, 1:basis.order)))

    for j in 1:basis.order^2
        for p in 1:basis.order^2
            Dx[j,p] = lagrange_diff(basis.quadpoints,linearindex[j][1],basis.quadpoints[linearindex[p][1]]) * lagrange_1d(basis.quadpoints, linearindex[j][2], basis.quadpoints[linearindex[p][2]])
        end
    end

    for j in 1:basis.order^2
        for p in 1:basis.order^2
            Dy[j,p] = lagrange_diff(basis.quadpoints,linearindex[j][2],basis.quadpoints[linearindex[p][2]]) * lagrange_1d(basis.quadpoints, linearindex[j][1], basis.quadpoints[linearindex[p][1]])
        end
    end
    hcat(Dx,Dy)
end


"""
    get_face_quadpoints(basis::Basis, face)

Return the quadrature points at the face `face` for basis `basis`.
"""
function get_face_quadpoints(basis::Basis, face)    # Is there a more efficient/ general way to do this?
    #ones(1.0, 1.0)

    face_quadpoints = zeros(basis.order,2)
    if face == left
        face_quadpoints[:,1] .= 0
        face_quadpoints[:,2] = basis.quadpoints
    elseif face == top
        face_quadpoints[:,1] = basis.quadpoints
        face_quadpoints[:,2] .= 1
    elseif face == right
        face_quadpoints[:,1] .= 1
        face_quadpoints[:,2] = basis.quadpoints
    elseif face == bottom
        face_quadpoints[:,1] = basis.quadpoints
        face_quadpoints[:,2] .= 0
    end
    face_quadpoints
end

function face_projection_matrix(basis, face)
    #ones(1,1)
    face_quadpoints = get_face_quadpoints(basis,face)

    # Matrix of required size
    P_ij = zeros( basis.order  , (basis.order)^2 )

    # Linear indexing
    linearindex = vec(collect(Iterators.product(1:basis.order, 1:basis.order)))

    for j = 1:(basis.order)^2
        for i = 1:basis.order
            P_ij[i,j] = lagrange_1d( basis.quadpoints, linearindex[j][1] , face_quadpoints[i,1] ) *  lagrange_1d( basis.quadpoints, linearindex[j][2], face_quadpoints[i,2] )
        end
    end
    return P_ij
end


"""
    evaluate_m_to_n_vandermonde_basis(basis)

Return the Vandermonde matrix that converts between the 2D-modal
(normalized) Legendre-basis and the 2D-nodal tensor-product basis built
with `basis`.
"""
function evaluate_m_to_n_vandermonde_basis(basis)
    Array{Float64}(I, length(basis), length(basis))
end
