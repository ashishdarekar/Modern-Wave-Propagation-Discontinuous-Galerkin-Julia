struct Elastic <: Equation end
@declare_dofs Elastic [:stxx, :styy, :stxy, :v, :w]

struct Elastic_waves <: Scenario end
struct GaussianPoint <: Scenario end
struct Elastic_source <: Scenario end
struct StiffInclusion <: Scenario end

function is_periodic_boundary(equation::Elastic, scenario::Elastic_waves)
    # Elastic_waves scenario with periodic boundary conditions
    true
end

function is_periodic_boundary(equation::Elastic, scenario::GaussianPoint)
    # Elastic_waves scenario with periodic boundary conditions
    true
end

function is_periodic_boundary(equation::Elastic, scenario::Elastic_source)
    # Elastic_waves scenario with periodic boundary conditions
    true
end

function is_periodic_boundary(equation::Elastic, scenario::StiffInclusion)
    # Elastic_waves scenario with periodic boundary conditions
    true
end
"""
# Refelcting Boundary conditions
function evaluate_boundary(eq::Elastic, scenario::Elastic_waves, face, normalidx, dofsface, dofsfaceneigh)
# dofsface and dofsfaceneigh have shape (num_2d_quadpoints, dofs)
# you need to set dofsfaceneigh
    dofsfaceneigh .= dofsface
    if normalidx == 1
        dofsfaceneigh[:, 1] = -dofsface[:, 1]
    else
        dofsfaceneigh[:, 2] = -dofsface[:, 2]
    end
end
"""

function get_initial_values(eq::Elastic, scenario::GaussianPoint, global_position; t=0.0)
    x, y = global_position
    if is_analytical_solution(eq, scenario)
        println("No analytical solution implemented")
    else
        # From Equation (22.67) from book of Leveque
        # TODO need to look fot actual values
        # some constants (Example 22.2 from book of Leveque only left side material) 
        λ = 2
        μ = 1
        ρ = 1
        
        #if (0.2 < x < 0.35) # As we have only domain fom (0 to 1)
        if (x <= 0.5 || y <= 0.5) || (x >= 0.5 && y >= 0.5) 
            #[λ+2*μ , λ , 0 , sqrt(5), 0]
            [0, 0, 0, exp(-1000 * (x - 0.5)^2 - 1000 *(y - 0.5)^2) + 1 , exp(-1000 * (x - 0.5)^2 - 1000 *(y - 0.5)^2) + 1 ]
        else
            [0 , 0 , 0 , 1 , 1]
        end

    end
end

function get_initial_values(eq::Elastic, scenario::Elastic_waves, global_position; t=0.0)
    x, y = global_position
    if is_analytical_solution(eq, scenario)
        println("No analytical solution implemented")
    else
        # From Equation (22.67) from book of Leveque
        # TODO need to look fot actual values
        # some constants (Example 22.2 from book of Leveque only left side material) 
        λ = 2
        μ = 1
        ρ = 1
        
        if (0.2 < x < 0.35) # As we have only domain fom (0 to 1)
            [λ+2*μ , λ , 0 , sqrt(5), 0]
        else
            [0 , 0 , 0 , 0 , 0]
        end

    end
end

<<<<<<< Updated upstream
=======
function get_initial_values(eq::Elastic, scenario::StiffInclusion, global_position; t=0.0)
    x, y = global_position
    if is_analytical_solution(eq, scenario)
        println("No analytical solution implemented")
    else
        # From Equation (22.67) from book of Leveque
        # TODO need to look fot actual values
        # some constants (Example 22.2 from book of Leveque only left side material) 
        λ = 4
        μ = 0.5
        ρ = 1
        
        if (0.6 < x < 0.65) # As we have only domain fom (0 to 1)
            [λ+2*μ , λ , 0 , sqrt(5), 0]
        else
            [0 , 0 , 0 , 0 , 0]
        end

    end
end

>>>>>>> Stashed changes
function get_initial_values(eq::Elastic, scenario::Elastic_source, global_position; t=0.0)
    x, y = global_position
    if is_analytical_solution(eq, scenario)
        #println("No analytical solution implemented")
        [( 1 * sin(6*pi*x/100) * sin(6*pi*y/100) * sin(6*pi*t/100) ) , ( 2 * sin(2*pi*x/50) * sin(2*pi*y/50) * sin(2*pi*t/50) ), ( 3 * sin(2*pi*x/100) * sin(2*pi*y/100) * sin(2*pi*t/100) ), ( 4 * sin(2*pi*x/50) * sin(2*pi*y/50) * sin(2*pi*t/10) ), ( 5 * sin(6*pi*x/100) * sin(6*pi*y/100) * sin(2*pi*t/5) ) ]
    else
        [0,0,0,0,0]
    end
    
end

<<<<<<< Updated upstream
function get_initial_values(eq::Elastic, scenario::StiffInclusion, global_position; t=0.0)
    x, y = global_position
    if is_analytical_solution(eq, scenario)
        println("No analytical solution implemented")
    else
        # From Equation (22.67) from book of Leveque
        # TODO need to look fot actual values
        # some constants (Example 22.2 from book of Leveque only left side material) 
        λ = 2
        μ = 1
        ρ = 1
        
        if (0.6 < x < 0.65) # As we have only domain fom (0 to 1)
            [λ+2*μ , λ , 0 , 2, 0]
        else
            [0 , 0 , 0 , 0 , 0]
        end
    end
end


=======
>>>>>>> Stashed changes

function is_analytical_solution(equation::Elastic, scenario::Elastic_waves)
    false
end

function is_analytical_solution(equation::Elastic, scenario::GaussianPoint)
    false
end

function is_analytical_solution(equation::Elastic, scenario::Elastic_source)
    true
end

function is_analytical_solution(equation::Elastic, scenario::StiffInclusion)
    false
end


# See equation (22.39) from book of Leveque for clarification 
function evaluate_flux(eq::Elastic, celldofs, cellflux)
    # size(cellflux)=(2*order*order, ndof)=(2,5)
    # size(celldofs)=(order*order, ndof)=(1,5)
    
    order_sq, ndof = size(celldofs)
    for j=1:order_sq
        # explicit names of dofs for clarity
        stxx = celldofs[j, 1]
        styy = celldofs[j, 2]
        stxy = celldofs[j, 3]
        v = celldofs[j, 4]
        w = celldofs[j, 5]
        
        #TODO need to look fot actual values
        #some constants (Example 22.2 from book of Leveque only left side material) 
        λ = 2
        μ = 1
        ρ = 1

        # here flux is computed elementwise
        cellflux[j,1]= (-(λ+2*μ)) * v
        cellflux[j + order_sq, 1] = (-λ) * w

        cellflux[j,2] = (-λ) * v
        cellflux[j + order_sq, 2] = (-(λ+2*μ)) * w

        cellflux[j,3] = (-μ) * w
        cellflux[j + order_sq, 3] = (-μ) * v

        cellflux[j,4] = (-1/ρ) * stxx
        cellflux[j + order_sq, 4] = (-1/ρ) * stxy

        cellflux[j,5] = (-1/ρ) * stxy
        cellflux[j + order_sq, 5] = (-1/ρ) * styy

    end

end

function evaluate_flux_with_stiffinclusion(eq::Elastic, scenario::StiffInclusion, celldofs, cellflux, global_position)
    # size(cellflux)=(2*order*order, ndof)=(2,5)
    # size(celldofs)=(order*order, ndof)=(1,5)
    x,y = global_position
    order_sq, ndof = size(celldofs)
    for j=1:order_sq
        # explicit names of dofs for clarity
        stxx = celldofs[j, 1]
        styy = celldofs[j, 2]
        stxy = celldofs[j, 3]
        v = celldofs[j, 4]
        w = celldofs[j, 5]
        
        if (x <=0.5 && y <= 0.5) ||  ( x<=0.5 && y <= ((-0.5*x) + 0.75) )
            λ = 4
            μ = 0.5
            ρ = 1
        else
            λ = 2
            μ = 1
            ρ = 1         
        end

        # here flux is computed elementwise
        cellflux[j,1]= (-(λ+2*μ)) * v
        cellflux[j + order_sq, 1] = (-λ) * w

        cellflux[j,2] = (-λ) * v
        cellflux[j + order_sq, 2] = (-(λ+2*μ)) * w

        cellflux[j,3] = (-μ) * w
        cellflux[j + order_sq, 3] = (-μ) * v

        cellflux[j,4] = (-1/ρ) * stxx
        cellflux[j + order_sq, 4] = (-1/ρ) * stxy

        cellflux[j,5] = (-1/ρ) * stxy
        cellflux[j + order_sq, 5] = (-1/ρ) * styy

    end

end


# ATTEMPT-1
"""
# Evaluation of flux with Source term 
# See equation (57-63) from paper of Martin Kaser for clarification 
function evaluate_flux_with_source(eq::Elastic, celldofs, cellflux, cell, basis)

    order_sq, ndof = size(celldofs)

    #TODO need to look fot actual values
    #some constants (Example 22.2 from book of Leveque only left side material) 
    λ = 4
    μ = 0.5
    ρ = 1

    up0 = [1,2,3,4,5]
    lamdapx = [42.9,50,100,50,42.9]
    lamdapy = [42.9,50,100,50,42.9]
    timep = [42.9,50,100,10,5]
    
    upsol = [
            (x,y,t) -> [up0[1] * sin((2*pi*x) / lamdapx[1]) * sin((2*pi*y) / lamdapy[1]) * sin((2*pi*t) / timep[1])],
            (x,y,t) -> [up0[2] * sin((2*pi*x) / lamdapx[2]) * sin((2*pi*y) / lamdapy[2]) * sin((2*pi*t) / timep[2])],
            (x,y,t) -> [up0[3] * sin((2*pi*x) / lamdapx[3]) * sin((2*pi*y) / lamdapy[3]) * sin((2*pi*t) / timep[3])],
            (x,y,t) -> [up0[4] * sin((2*pi*x) / lamdapx[4]) * sin((2*pi*y) / lamdapy[4]) * sin((2*pi*t) / timep[4])],
            (x,y,t) -> [up0[5] * sin((2*pi*x) / lamdapx[5]) * sin((2*pi*y) / lamdapy[5]) * sin((2*pi*t) / timep[5])]
            ]
    
    dudtp = [
            (x,y,t) -> [((2*pi) * timep[1]) * up0[1] * sin((2*pi*x) / lamdapx[1]) * sin((2*pi*y) / lamdapy[1]) * cos((2*pi*t) / timep[1])],
            (x,y,t) -> [((2*pi) * timep[2]) * up0[2] * sin((2*pi*x) / lamdapx[2]) * sin((2*pi*y) / lamdapy[2]) * cos((2*pi*t) / timep[2])],
            (x,y,t) -> [((2*pi) * timep[3]) * up0[3] * sin((2*pi*x) / lamdapx[3]) * sin((2*pi*y) / lamdapy[3]) * cos((2*pi*t) / timep[3])],
            (x,y,t) -> [((2*pi) * timep[4]) * up0[4] * sin((2*pi*x) / lamdapx[4]) * sin((2*pi*y) / lamdapy[4]) * cos((2*pi*t) / timep[4])],
            (x,y,t) -> [((2*pi) * timep[5]) * up0[5] * sin((2*pi*x) / lamdapx[5]) * sin((2*pi*y) / lamdapy[5]) * cos((2*pi*t) / timep[5])]
            ]
        
    dudxp = [
            (x,y,t) -> [((2*pi) / lamdapx[1]) * up0[1] * cos((2*pi*x) / lamdapx[1]) * sin((2*pi*y) / lamdapy[1]) * sin((2*pi*t) / timep[1])],
            (x,y,t) -> [((2*pi) / lamdapx[2]) * up0[2] * cos((2*pi*x) / lamdapx[2]) * sin((2*pi*y) / lamdapy[2]) * sin((2*pi*t) / timep[2])],
            (x,y,t) -> [((2*pi) / lamdapx[3]) * up0[3] * cos((2*pi*x) / lamdapx[3]) * sin((2*pi*y) / lamdapy[3]) * sin((2*pi*t) / timep[3])],
            (x,y,t) -> [((2*pi) / lamdapx[4]) * up0[4] * cos((2*pi*x) / lamdapx[4]) * sin((2*pi*y) / lamdapy[4]) * sin((2*pi*t) / timep[4])],
            (x,y,t) -> [((2*pi) / lamdapx[5]) * up0[5] * cos((2*pi*x) / lamdapx[5]) * sin((2*pi*y) / lamdapy[5]) * sin((2*pi*t) / timep[5])]
            ]

    dudyp = [
            (x,y,t) -> [((2*pi) / lamdapy[1]) * up0[1] * sin((2*pi*x) / lamdapx[1]) * cos((2*pi*y) / lamdapy[1]) * sin((2*pi*t) / timep[1])],
            (x,y,t) -> [((2*pi) / lamdapy[2]) * up0[2] * sin((2*pi*x) / lamdapx[2]) * cos((2*pi*y) / lamdapy[2]) * sin((2*pi*t) / timep[2])],
            (x,y,t) -> [((2*pi) / lamdapy[3]) * up0[3] * sin((2*pi*x) / lamdapx[3]) * cos((2*pi*y) / lamdapy[3]) * sin((2*pi*t) / timep[3])],
            (x,y,t) -> [((2*pi) / lamdapy[4]) * up0[4] * sin((2*pi*x) / lamdapx[4]) * cos((2*pi*y) / lamdapy[4]) * sin((2*pi*t) / timep[4])],
            (x,y,t) -> [((2*pi) / lamdapy[5]) * up0[5] * sin((2*pi*x) / lamdapx[5]) * cos((2*pi*y) / lamdapy[5]) * sin((2*pi*t) / timep[5])]
            ]

    #Hence source term is according to equation (2) from same paper size of soure term is: size is (5,1)
    sp = [
        (x,y,z) -> [dudtp[1](x,y,z) + (-(λ+2*μ)) * dudxp[1](x,y,z) + (-λ) * dudyp[1](x,y,z)],
        (x,y,z) -> [dudtp[2](x,y,z) + (-λ) * dudxp[2](x,y,z) + (-(λ+2*μ)) * dudyp[2](x,y,z)],
        (x,y,z) -> [dudtp[3](x,y,z) + (-μ) * dudxp[3](x,y,z) + (-μ) * dudyp[3](x,y,z)],
        (x,y,z) -> [dudtp[4](x,y,z) + (-1/ρ) * dudxp[4](x,y,z) + (-1/ρ)  * dudyp[4](x,y,z)],
        (x,y,z) -> [dudtp[5](x,y,z) + (-1/ρ) * dudxp[5](x,y,z) + (-1/ρ)  * dudyp[5](x,y,z)]
        ]
    
    #Updating Flux including source terms
    linearindex = vec(collect(Iterators.product(1:basis.order, 1:basis.order))) #linear indexing

    for j=1:order_sq
        # explicit names of dofs for clarity
        stxx = celldofs[j, 1]
        styy = celldofs[j, 2]
        stxy = celldofs[j, 3]
        v = celldofs[j, 4]
        w = celldofs[j, 5]
        
        # Finding global Coordinates of Cell
        localx,localy = globalposition(cell, (basis.quadpoints[linearindex[j][1]], basis.quadpoints[linearindex[j][2]]))
        #@info localx
        #@info localy
        localz = 0.0 # as problem is 2 Dimensional
        
        # here flux is computed elementwise, one source term is divided into 2 parts equally, one along X and other along Y
        cellflux[j,1]= ( (-(λ+2*μ)) * v ) - ( ((sp[1](localx,localy,localz))[1][1])/2 )
        cellflux[j + order_sq, 1] = ( (-λ) * w )- ( ((sp[1](localx,localy,localz))[1][1])/2 )
    
        cellflux[j,2] = ( (-λ) * v ) - ( ((sp[2](localx,localy,localz))[1][1])/2 )
        cellflux[j + order_sq, 2] = ( (-(λ+2*μ)) * w ) - ( ((sp[2](localx,localy,localz))[1][1])/2 )
    
        cellflux[j,3] = ( (-μ) * w ) - ( ((sp[3](localx,localy,localz))[1][1])/2 )
        cellflux[j + order_sq, 3] = ( (-μ) * v ) - ( ((sp[3](localx,localy,localz))[1][1])/2 )
    
        cellflux[j,4] = ( (-1/ρ) * stxx ) - ( ((sp[4](localx,localy,localz))[1][1])/2 )
        cellflux[j + order_sq, 4] = ( (-1/ρ) * stxy ) - ( ((sp[4](localx,localy,localz))[1][1])/2 )
    
        cellflux[j,5] = ( (-1/ρ) * stxy ) - ( ((sp[5](localx,localy,localz))[1][1])/2 )
        cellflux[j + order_sq, 5] = ( (-1/ρ) * styy ) - ( ((sp[5](localx,localy,localz))[1][1])/2 )
    end

end
"""

# ATTEMPT 2
# Evaluation of Source term
# See equation (57-63) from paper of Martin Kaser for clarification 
function evaluate_source(eq::Elastic, time, position)
<<<<<<< Updated upstream
    # Source term foe general problem (section 6.1 in same paper)
=======
    """
>>>>>>> Stashed changes
    λ = 2
    μ = 1
    ρ = 1
    
    up0 = [1,2,3,4,5]
    lamdapx = [42.9,50,100,50,42.9]
    lamdapy = [42.9,50,100,50,42.9]
    timep = [42.9,50,100,10,5]

    num_dofs = size(up0)
    upsol = zeros(num_dofs)
    dudtp = zeros(num_dofs)
    dudxp = zeros(num_dofs)
    dudyp = zeros(num_dofs)
    sp = zeros(num_dofs)
    
    upsol = [
            (x,y,t) -> [up0[1] * sin((2*pi*x) / lamdapx[1]) * sin((2*pi*y) / lamdapy[1]) * sin((2*pi*t) / timep[1])],
            (x,y,t) -> [up0[2] * sin((2*pi*x) / lamdapx[2]) * sin((2*pi*y) / lamdapy[2]) * sin((2*pi*t) / timep[2])],
            (x,y,t) -> [up0[3] * sin((2*pi*x) / lamdapx[3]) * sin((2*pi*y) / lamdapy[3]) * sin((2*pi*t) / timep[3])],
            (x,y,t) -> [up0[4] * sin((2*pi*x) / lamdapx[4]) * sin((2*pi*y) / lamdapy[4]) * sin((2*pi*t) / timep[4])],
            (x,y,t) -> [up0[5] * sin((2*pi*x) / lamdapx[5]) * sin((2*pi*y) / lamdapy[5]) * sin((2*pi*t) / timep[5])]
            ]
    for i in 1:
    dudtp = [
            (x,y,t) -> [((2*pi) * timep[1]) * up0[1] * sin((2*pi*x) / lamdapx[1]) * sin((2*pi*y) / lamdapy[1]) * cos((2*pi*t) / timep[1])],
            (x,y,t) -> [((2*pi) * timep[2]) * up0[2] * sin((2*pi*x) / lamdapx[2]) * sin((2*pi*y) / lamdapy[2]) * cos((2*pi*t) / timep[2])],
            (x,y,t) -> [((2*pi) * timep[3]) * up0[3] * sin((2*pi*x) / lamdapx[3]) * sin((2*pi*y) / lamdapy[3]) * cos((2*pi*t) / timep[3])],
            (x,y,t) -> [((2*pi) * timep[4]) * up0[4] * sin((2*pi*x) / lamdapx[4]) * sin((2*pi*y) / lamdapy[4]) * cos((2*pi*t) / timep[4])],
            (x,y,t) -> [((2*pi) * timep[5]) * up0[5] * sin((2*pi*x) / lamdapx[5]) * sin((2*pi*y) / lamdapy[5]) * cos((2*pi*t) / timep[5])]
            ]
        
    dudxp = [
            (x,y,t) -> [((2*pi) / lamdapx[1]) * up0[1] * cos((2*pi*x) / lamdapx[1]) * sin((2*pi*y) / lamdapy[1]) * sin((2*pi*t) / timep[1])],
            (x,y,t) -> [((2*pi) / lamdapx[2]) * up0[2] * cos((2*pi*x) / lamdapx[2]) * sin((2*pi*y) / lamdapy[2]) * sin((2*pi*t) / timep[2])],
            (x,y,t) -> [((2*pi) / lamdapx[3]) * up0[3] * cos((2*pi*x) / lamdapx[3]) * sin((2*pi*y) / lamdapy[3]) * sin((2*pi*t) / timep[3])],
            (x,y,t) -> [((2*pi) / lamdapx[4]) * up0[4] * cos((2*pi*x) / lamdapx[4]) * sin((2*pi*y) / lamdapy[4]) * sin((2*pi*t) / timep[4])],
            (x,y,t) -> [((2*pi) / lamdapx[5]) * up0[5] * cos((2*pi*x) / lamdapx[5]) * sin((2*pi*y) / lamdapy[5]) * sin((2*pi*t) / timep[5])]
            ]

    dudyp = [
            (x,y,t) -> [((2*pi) / lamdapy[1]) * up0[1] * sin((2*pi*x) / lamdapx[1]) * cos((2*pi*y) / lamdapy[1]) * sin((2*pi*t) / timep[1])],
            (x,y,t) -> [((2*pi) / lamdapy[2]) * up0[2] * sin((2*pi*x) / lamdapx[2]) * cos((2*pi*y) / lamdapy[2]) * sin((2*pi*t) / timep[2])],
            (x,y,t) -> [((2*pi) / lamdapy[3]) * up0[3] * sin((2*pi*x) / lamdapx[3]) * cos((2*pi*y) / lamdapy[3]) * sin((2*pi*t) / timep[3])],
            (x,y,t) -> [((2*pi) / lamdapy[4]) * up0[4] * sin((2*pi*x) / lamdapx[4]) * cos((2*pi*y) / lamdapy[4]) * sin((2*pi*t) / timep[4])],
            (x,y,t) -> [((2*pi) / lamdapy[5]) * up0[5] * sin((2*pi*x) / lamdapx[5]) * cos((2*pi*y) / lamdapy[5]) * sin((2*pi*t) / timep[5])]
            ]

    #Hence source term is according to equation (2) from same paper size of soure term is: size is (5,1)
    sp = [
        (x,y,t) -> [dudtp[1](x,y,t) + (-(λ+2*μ)) * dudxp[1](x,y,t) + (-λ) * dudyp[1](x,y,t)],
        (x,y,t) -> [dudtp[2](x,y,t) + (-λ) * dudxp[2](x,y,t) + (-(λ+2*μ)) * dudyp[2](x,y,t)],
        (x,y,t) -> [dudtp[3](x,y,t) + (-μ) * dudxp[3](x,y,t) + (-μ) * dudyp[3](x,y,t)],
        (x,y,t) -> [dudtp[4](x,y,t) + (-1/ρ) * dudxp[4](x,y,t) + (-1/ρ)  * dudyp[4](x,y,t)],
        (x,y,t) -> [dudtp[5](x,y,t) + (-1/ρ) * dudxp[5](x,y,t) + (-1/ρ)  * dudyp[5](x,y,t)]
        ]
    """
    
    # Source term from Lamb's problem:Ricker welvet (section 6.1 in same paper)
    rhosolid = 2200
    cp=3200
    cs=1847
    xs=[0.5,0.5] #Source application
    fc=14.5 #Ricker welvet central frequesncy
    a1=-2000
    a2=-(pi*pi*fc*fc)
    td=0.08 #source delay
    px,py = position
    y = [px-0.5, py-1]
    y_norm = y[1]^2 + y[2]^2
    dirac_delta = 5*(1/sqrt(pi))*exp(-25*(y_norm))

    g= a1 * (0.5 + (a2 * (time-td)^2 )) * (exp( a2 * (time -td)^2 ))

    sp = [(x,y,z) -> [0],
          (x,y,z) -> [0],
          (x,y,z) -> [0],
          (x,y,z) -> [0],
          (x,y,z) -> [(g/rhosolid)*dirac_delta]
        ]
<<<<<<< Updated upstream
    """
=======
    
>>>>>>> Stashed changes
    
    #@info "from function"
    #@info sp[1](1,2,3)    
    #Returning sp
    sp

end


# See equation (22.17 and 22.18) from book of Leveque
function max_eigenval(eq::Elastic, celldata, normalidx)
    maxeigenval = 0.

    #TODO need to look fot actual values
    #some constants (Example 22.2 from book of Leveque only left side material) 
    λ = 4
    μ = 0.5
    ρ = 1
    
    cp = sqrt((λ+2*μ)/ρ)
    cs = sqrt(μ/ρ)
    
    if (cp > maxeigenval && cp > cs)
        maxeigenval = cp
    elseif (cs > maxeigenval && cs > cp)
        maxeigenval = cs
    end
    maxeigenval
end
