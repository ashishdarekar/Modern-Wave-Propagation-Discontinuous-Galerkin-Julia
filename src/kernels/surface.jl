"""
    struct BuffersFaceIntegral

Stores all temporary buffers needed during face integrals.
Avoids costly re-allocation.
"""
struct BuffersFaceIntegral
    dofsface::Array{Float64,2}
    dofsfaceneigh::Array{Float64,2}
    fluxface::Array{Float64,2}
    fluxfaceneigh::Array{Float64,2}
    numericalflux::Array{Float64,2}
    numericalflux_scaled::Array{Float64,2}

    function BuffersFaceIntegral(basis, ndofs)
        @assert basis.dimensions == 2
        basissize_1d = size(basis, 1)
        basissize_nd = length(basis)
        dofsface = zeros(basissize_1d,ndofs)
        dofsfaceneigh = similar(dofsface)
        fluxface = zeros(basissize_1d * 2,ndofs)
        fluxfaceneigh = similar(fluxface)
        numericalflux = similar(dofsface)
        numericalflux_scaled = similar(numericalflux)

        new(dofsface, dofsfaceneigh, fluxface, fluxfaceneigh, numericalflux, numericalflux_scaled)
    end
end

"""
    rusanov(eq, dofs, dofsneigh, flux, fluxneigh, dx, normalidx, normalsign, numericalflux)

Computes the Rusanov (or local Lax-Friedrichs) numerical flux for degrees of freedom
`dofs`, degrees of freedoms of the neighbor `dofsneigh`, flux `flux`, flux of neighbor
`fluxneigh`, cellsize `dx`.
All quantities are assumed to be represented by a basis on the reference line.
The face is parametrized by a index `normalidx`, where 1 stands for a face in x-direction
and 2 for a face in y-direction.
The sign of the outer normal of the face is given by `normalsign`.

The numerical flux is stored in `numericalflux`.
Method also returns the maximal eigenvalue.
"""

function rusanov(eq, dofs, dofsneigh, flux, fluxneigh, dx, normalidx, normalsign, numericalflux)
    maxeigenval_left = max_eigenval(eq, dofs, normalidx)
    maxeigenval_right = max_eigenval(eq, dofsneigh, normalidx)
    maxeigenval = max(maxeigenval_left, maxeigenval_right)

    basissize_1d = size(numericalflux)[1]
    for i=1:size(dofs)[2]
        #check which value of flux need to use and also the direction of the normal
        for j = 1:basissize_1d
            if normalidx == 1 #For faces in the X direction that is left and right
            #@info dofs[i]
                numericalflux[j,i] = ( 1/2 * ( flux[j,i] + fluxneigh[j,i] )   ) +
                                                 ( maxeigenval* normalsign * 1/2 * (dofs[j,i] - dofsneigh[j,i]) )

            else #For faces in the Y direction that is Top and Bottom
                numericalflux[j,i] = ( 1/2 * ( flux[j+basissize_1d,i] + fluxneigh[j+basissize_1d,i] )   ) +
                                                 ( maxeigenval* normalsign * 1/2 * (dofs[j,i] - dofsneigh[j,i]) )
            end

            if eq == Acoustic()
                if i > 3
                    numericalflux[j, i] = 0  # For acoustic, material flux to 0
                end
            end
        end
    end

    #numericalflux =  (1/dx) * ( ( 1/2 * ( flux + fluxneigh ) )  +  ( maxeigenval* normalsign * (dofs - dofsneigh ) * 1/2 ) )
    maxeigenval # maximum eigenvalue
end

"""
    exactRiemann(eq, dofs, dofsneigh, flux, fluxneigh, dx, normalidx, normalsign, numericalflux)

Computes the exact flux for degrees of freedom `dofs`,
degrees of freedoms of the neighbor `dofsneigh`, flux `flux`,
flux of neighbor `fluxneigh`, cellsize `dx`.
All quantities are assumed to be represented by a basis on the reference line.
The face is parametrized by a index `normalidx`, where 1 stands for a face in x-direction
and 2 for a face in y-direction.
The sign of the outer normal of the face is given by `normalsign`.

The flux is stored in `numericalflux`.
Method also returns the maximal eigenvalue.
"""

function exactRiemann(eq, dofs, dofsneigh, flux, fluxneigh, normalidx, normalsign, numericalflux)
    maxeigenval_left = max_eigenval(eq, dofs, normalidx)
    maxeigenval_right = max_eigenval(eq, dofsneigh, normalidx)
    maxeigenval = max(maxeigenval_left, maxeigenval_right)

    # constants relevant to elastic equation
    λ = 2
    μ = 1
    ρ = 1
    cp = sqrt((λ + 2*μ)/ρ) # speed p-waves
    cs = sqrt(μ/ρ) #speed s-waves

    basissize_1d = size(numericalflux)[1]
    avr_flux = 1/2 * (flux + fluxneigh)
    dif_flux = 1/2 * (flux - fluxneigh)
    for j = 1:basissize_1d
        if normalidx == 1 #For faces in the X direction that is left and right
            numericalflux[j,1] = avr_flux[j,1] - cp*ρ * normalsign * dif_flux[j,4]
            numericalflux[j,2] = (avr_flux[j,2] - normalsign * dif_flux[j,2]
                                    - λ/(cp) * normalsign * dif_flux[j,4]
                                    + λ/(cp^2*ρ) * normalsign * dif_flux[j,1])
            numericalflux[j,3] = avr_flux[j,3] - cs*ρ * normalsign * dif_flux[j,5]
            numericalflux[j,4] = avr_flux[j,4] - 1/(cp*ρ) * normalsign * dif_flux[j,1]
            numericalflux[j,5] = avr_flux[j,5] - 1/(cs*ρ) * normalsign * dif_flux[j,3]
        else #For faces in the Y direction that is Top and Bottom
            numericalflux[j,1] = (avr_flux[j+basissize_1d,1] - normalsign * dif_flux[j+basissize_1d,1]
                                    - λ/cp * normalsign * dif_flux[j+basissize_1d,5]
                                    + λ/(cp^2*ρ) * normalsign * dif_flux[j+basissize_1d,2])
            numericalflux[j,2] = avr_flux[j+basissize_1d,2] - cp*ρ* normalsign * dif_flux[j+basissize_1d,5]
            numericalflux[j,3] = avr_flux[j+basissize_1d,3] - cs*ρ* normalsign * dif_flux[j+basissize_1d,4]
            numericalflux[j,4] = avr_flux[j+basissize_1d,4] - 1/(cs*ρ)* normalsign * dif_flux[j+basissize_1d,3]
            numericalflux[j,5] = avr_flux[j+basissize_1d,5] - 1/(cp*ρ)* normalsign * dif_flux[j+basissize_1d,2]
        end
    end

    maxeigenval # maximum eigenvalue
end

"""
    project_to_faces(globals, dofs, flux, dofsface, fluxface, face)

Projects degrees of freedom `dofs` and `flux` to `face`.
Result is stored in `dofsface` and `fluxface`.
Projection matrices are stored in `globals`.
"""

function project_to_faces(globals, dofs, flux, dofsface, fluxface, face)
    #dofsface .= dofs
    i = globals.basissize_1d
    dofsface .= globals.project_dofs_to_face[face] * dofs
    fluxface[1:i,:] = globals.project_flux_to_face[face] * flux[1:i*i,:]
    fluxface[i+1:end,:] = globals.project_flux_to_face[face] * flux[i*i+1:end,:]
    #dofs = 4x3  and flux = 8x3
    #face proj matrix = 2x4
    #order = 1
    #fluxface .= flux
end


"""
    evaluate_face_integral(eq, globals, buffers, cell, face, celldu)

Computes face integrals for equation `eq`, cell `cell` on face `face`.
Global matrices are passed in `globals`, buffers in `buffers`.
Result is stored in `celldu`. `solver` is specified in config file.
"""

function evaluate_face_integral(eq, globals, buffers, cell, face, celldu, solver)
    normalidx = globals.normalidxs[face]
    normalsign = globals.normalsigns[face]

    # Compute Riemann solver (in normal)
    buffers.numericalflux .= 0

    if solver == "exact"
        maxeigenval = exactRiemann(eq, buffers.dofsface, buffers.dofsfaceneigh, buffers.fluxface, buffers.fluxfaceneigh, normalidx, normalsign, buffers.numericalflux)
    else
        maxeigenval = rusanov(eq, buffers.dofsface, buffers.dofsfaceneigh, buffers.fluxface, buffers.fluxfaceneigh, cell.size[1], normalidx, normalsign, buffers.numericalflux)
    end

    #scale numerical flux with quadweights.
    buffers.numericalflux_scaled .= globals.quadweights_nd.*buffers.numericalflux*area(cell)

    celldu .= celldu .+ -normalsign*globals.project_dofs_from_face[face] * buffers.numericalflux_scaled

    maxeigenval
end
