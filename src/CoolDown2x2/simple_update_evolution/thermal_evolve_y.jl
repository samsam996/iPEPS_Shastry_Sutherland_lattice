

include("diag_sqrt.jl")
include("invert_diag.jl")
include("truncate_ey.jl")


function thermal_evolve_y(Gamma::lattice2x1,lambdax::lattice2x1,lambday::lattice2x1,
    physical_legs::lattice_ind2x1,gy::lattice2x1,gt::Matrix{Symbol},D::Int64)
   
    N = size(gt)[1];
    f(x) = mod(x-1,N) + 1;
    sum_fy = 0
    
      for j = 1:1:N
        for i = 1:1:1

            lambda1 = getfield(lambdax,gt[f(i-1),f(j)]);
            lambda2 = getfield(lambday,gt[f(i),f(j)]);
            lambda3 = getfield(lambdax,gt[f(i),f(j)]);
            lambda4 = getfield(lambday,gt[f(i),f(j-1)]);
            lambda5 = getfield(lambdax,gt[f(i-1),f(j+1)]);
            lambda6 = getfield(lambdax,gt[f(i),f(j+1)]);
            lambda7 = getfield(lambday,gt[f(i),f(j+1)]);

            gamma_a = getfield(Gamma,gt[f(i),f(j)])
            gamma_c = getfield(Gamma,gt[f(i),f(j+1)])

            #       l4
            #  l1 - A(i,j) - l3
            #       l2
            #  l5 - C(i,j+1)- l6
            #       l7

            A = (((gamma_a*lambda1)*diag_sqrt(lambda2))*lambda3)*lambda4;
            C = (((gamma_c*lambda5)*lambda7)*lambda6)*diag_sqrt(lambda2);

            ya1 = commonind(lambda2,gamma_a); # ya1
            ya2 = commonind(lambda2,gamma_c); # ya2

            ind_lambda2 = commonind(A,lambda2)
            ind_new = Index(ind_lambda2.space, dir = dir(ind_lambda2))
            A = A*delta(dag(ya2),ind_new)
            C = C*delta(dag(ya1),dag(ind_new))

            C,A,l2,fy = truncate_ey(C,A,i,j,gy,gt,physical_legs,ind_new,D);

            A = ((A*invert_diag(lambda1))*invert_diag(lambda3))*invert_diag(lambda4);
            C = ((C*invert_diag(lambda5))*invert_diag(lambda7))*invert_diag(lambda6);
            
            setproperty!(lambday,gt[f(i),f(j)],l2);
            setproperty!(Gamma,gt[f(i),f(j)],A);
            setproperty!(Gamma,gt[f(i),f(j+1)],C);

            sum_fy = sum_fy + fy

        end
    end

    return Gamma,lambday, sum_fy

end