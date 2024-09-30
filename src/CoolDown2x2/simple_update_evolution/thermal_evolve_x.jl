

include("diag_sqrt.jl")
include("invert_diag.jl")
include("truncate_ex.jl")


function thermal_evolve_x(Gamma::lattice2x1,lambdax::lattice2x1,lambday::lattice2x1,physical_legs::lattice_ind2x1,
    gx::lattice2x1,gt::Matrix{Symbol},D::Int64)

    N = size(gt)[1];
    f(x) = mod(x-1,N) + 1;
    sum_fx = 0

    for i = 1:1:N
        for j = 1:1:1

            lambda1 = getfield(lambdax,gt[f(i-1),f(j)])
            lambda2 = getfield(lambday,gt[f(i),f(j)])
            lambda3 = getfield(lambdax,gt[f(i),f(j)])
            lambda4 = getfield(lambday,gt[f(i),f(j-1)])
            lambda5 = getfield(lambday,gt[f(i+1),f(j)])
            lambda6 = getfield(lambday,gt[f(i+1),f(j-1)])
            lambda7 = getfield(lambdax,gt[f(i+1),f(j)])
            
            #      l4       l6
            # l1 - A - l3 - B - l7
            #      l2       l5

            gamma_a = getfield(Gamma,gt[f(i),f(j)])
            gamma_b = getfield(Gamma,gt[f(i+1),f(j)])

            A = (((gamma_a*lambda1)*lambda2)*lambda4)*diag_sqrt(lambda3); 
            B = (((gamma_b*diag_sqrt(lambda3))*lambda5)*lambda7)*lambda6; 

            ind_lambda3 = commonind(A,lambda3)
            ind_new = Index(ind_lambda3.space, dir = dir(ind_lambda3)) 
            xa1 = commonind(lambda3,gamma_a) 
            xa2 = commonind(lambda3,gamma_b)

            A = A*delta(ComplexF64,dag(xa2),ind_new)
            B = B*delta(ComplexF64,dag(ind_new),dag(xa1))


            AA,BB,l3,fx = truncate_ex(A,B,i,j,gx,gt,physical_legs,ind_new,D)
 
            AA = ((AA*invert_diag(lambda1))*invert_diag(lambda2))*invert_diag(lambda4)
            BB = ((BB*invert_diag(lambda5))*invert_diag(lambda7))*invert_diag(lambda6)
         
            setproperty!(lambdax,gt[f(i),f(j)],l3)
            setproperty!(Gamma,gt[f(i),f(j)],AA)
            setproperty!(Gamma,gt[f(i+1),f(j)],BB)

            sum_fx = sum_fx + fx
        end
    end

    return Gamma,lambdax, sum_fx


end