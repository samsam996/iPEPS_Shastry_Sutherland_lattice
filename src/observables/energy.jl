

include("horizontal_correlation.jl")
include("vertical_correlation.jl")
include("local_obs.jl")

function energy(C,T,tens_a,tens_A,gt::Matrix{Symbol},cxd,cyd,
    physical_legs,entanglement_legs,J1::Float64,J2::Float64,h::Float64,hs::Float64)

    ener_link = 0
    N = size(gt)[1]
    sx = 1/2*[0 1; 1 0]; sy = 1/2*[0 -1im; 1im 0]; sz = 1/2*[1 0; 0 -1]; id = [1 0; 0 1];

    # we don't want to parallelise the energy due to risk of memory shortage.
    for i = 1:N
        for j = 1:1
		
           println(i)
           @show horizontal = horizontal_correlation(C,T,tens_a,tens_A,gt,cxd,cyd,physical_legs,entanglement_legs,i,j)
           @show vertical = vertical_correlation(C,T,tens_a,tens_A,gt,cxd,cyd,physical_legs,entanglement_legs,i,j)
            ener_link = ener_link + horizontal + vertical             

        end
    end 

    # ener_link = sum(horizontal) + sum(vertical)

    obs1::Matrix{ComplexF64} = (kron(sz,sz) + kron(sx,sx) + kron(sy,sy)) 
    obs2::Matrix{ComplexF64} =  kron(sz,id) + kron(id,sz);
    expect_obs1 = 0
    expect_obs2 = 0


    mm = []

    for i = 1:N
        for j = 1:1
	    println(i)
            @show expect_obs1 = expect_obs1 + 
            local_obs(C,T,tens_a,tens_A,cxd,cyd,gt,physical_legs,entanglement_legs,i,j,obs1)
            @show expect_obs2 = expect_obs2 + 
            local_obs(C,T,tens_a,tens_A,cxd,cyd,gt,physical_legs,entanglement_legs,i,j,obs2)
            push!(mm,expect_obs2)
        end
    end

    if abs(hs) > 0
        energie = (J1*ener_link + J2*expect_obs1 - h*expect_obs2)/(N) - 1/N*(hs*(mm[1] - mm[2] + mm[4] - mm[5]))
    else
        energie = (J1*ener_link + J2*expect_obs1 - h*expect_obs2)/(N) 
    end

    return energie

end
