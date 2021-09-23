using LinearAlgebra

"""
    implementation of the Limiter
"""

# functions to be called when there is a need of Limiting the slopes
# data= dofs of particular cell
# k = varible number . which column in dofs to be modify
# basis =basis of a problem , node informaton
# a = a in linear eqaution, of that particular k
# b = b in linear eqaution, of that particular k
# c = c in linear eqaution, of that particular k
function update_solution(data, k, basis, a, b, c, cell)

    linearindex = vec(collect(Iterators.product(1:basis.order, 1:basis.order))) #linear indexing
    cellcenter = cell.center

    for i in 1:length(basis)
        x, y = globalposition(cell, (basis.quadpoints[linearindex[i][1]], basis.quadpoints[linearindex[i][2]]))
        data[i,k] = a + b * ( x - cellcenter[1] ) + c * ( y - cellcenter[2] )
    end

end

# Caluations required to decide limiting
function minmod_limiter(eq, scenario, dofs, grid, basis, ndofs, ndims)
    # Travelling each cell in a grid
    for i in eachindex(grid.cells)
        @views cell = grid.cells[i]
        facetypes = cell.facetypes
        @views data = dofs[:,:, cell.dataidx]
        weights = kron(basis.quadweights, basis.quadweights)

        # calculation of cell (weighted) averages of own cell and its neighbouring cells
        # calculation of a
        # a = ∑_n w_n u_n
        uown = vec(data' * weights)

        # For all neighbouring cells
        u_mean_neighbour = zeros(4,ndofs)
        for (i, neigh) in enumerate(cell.neighbors)
            @views dofsneigh = dofs[:,:,neigh.dataidx]
            u_mean_neighbour[i,:] =  dofsneigh' * weights
        end

        needs_limiting = BitArray(undef, ndofs) #initialise boolean for whether this cell will use limiting
        s_limit = Array{Float64}(undef, ndofs, ndims)

        for j in 1:ndims
            # points for (x - 1/2) or (y - 1/2)
            if j == 1
                points = repeat(basis.quadpoints, length(basis.quadpoints))
            else
                points = vcat(fill.(basis.quadpoints, length(basis.quadpoints))...)
            end

            #for x-dim:  s1 = 12* ∑_n w_n (x_n - 1/2) u_n
            s1 = 12 * vec(sum(weights .* (points .- 0.5) .* data, dims=1))
            s1 /= cell.size
            
            s2 = (uown - u_mean_neighbour[j, :]) / cell.size
            if facetypes[j] == 2  #at boundaries set slope to 0
                s2 = 0
            end
            s3 = (u_mean_neighbour[j + 2, :] - uown) / cell.size
            if facetypes[j + 2] == 2 #at boundaries set slope to 0
                s3 = 0
            end

            for k in 1:ndofs
                if s1[k]*s2[k] > 0 && s1[k]*s3[k] > 0  #if all have same sign
                    s_limit[k, j] = sign(s1[k])* min(abs(s1[k]),abs(s2[k]), abs(s3[k]))
                else
                    s_limit[k, j] = 0
                end

                if s_limit[k, j] != s1[k]
                    needs_limiting[k] = true
                end
            
            end

        end
        for k in 1:ndofs
            if needs_limiting[k]
                update_solution(data, k, basis, uown[k], s_limit[k,1], s_limit[k,2], cell)
            end
        end

    end
end
