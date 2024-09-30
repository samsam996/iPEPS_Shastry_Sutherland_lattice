
include("thermal_evolve_x.jl")
include("thermal_evolve_y.jl")

function simpleupdateJ1(Gamma::lattice2x1,lambdax::lattice2x1,lambday::lattice2x1,
    physical_legs::lattice_ind2x1,gt::Matrix{Symbol},nsu::Float64,J1::Float64,J2::Float64,h::Float64,hs::Float64,D::Int64,FreeEnergy)

    N = size(gt)[1]

    gx,gy,mu,mu_hs = SSMHamiltonian(gt,physical_legs,nsu,J1,J2,h,hs)
    gx2,gy2,mu2,mu2_hs = SSMHamiltonian(gt,physical_legs,2*nsu,J1,J2,h,hs)

    Gamma,lambdax, ffx1 = thermal_evolve_x(Gamma,lambdax,lambday,physical_legs,gx2,gt,D)
    Gamma,lambday, ffy = thermal_evolve_y(Gamma,lambdax,lambday,physical_legs,gy,gt,D)
    Gamma,lambdax, ffx2 = thermal_evolve_x(Gamma,lambdax,lambday,physical_legs,gx2,gt,D)

    FreeEnergy = FreeEnergy + ffx1 + ffx2 + ffy; 

    for i = 1:N
        for j =1:1
            gamma = getproperty(Gamma,gt[i,j])
            mu_a = getproperty(mu,gt[i,j])
            gamma = gamma*mu_a
            gamma = noprime(gamma)
            setproperty!(Gamma,gt[i,j],gamma)
        end
    end


    return Gamma,lambdax,lambday, FreeEnergy

end