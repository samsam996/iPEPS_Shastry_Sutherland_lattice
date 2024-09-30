

# using MKL, FileIO, JLD2, KrylovKit, LinearAlgebra, ITensors

mutable struct lattice2x1
    A::ITensor
    B::ITensor
    gg::Matrix{Symbol}
    gt::Matrix{Symbol}
    
    function lattice2x1()
        name = fieldnames(lattice2x1)
        g1 = [name[1] name[2]; name[2] name[1]]
        g2 = [name[1] name[2]; name[2] name[1]]
        new(ITensor(),ITensor(),g1,g2)
    end
    
    function lattice2x1(a::ITensor,b::ITensor)

        name = fieldnames(lattice2x1)
        g1 = [name[1] name[2]; name[2] name[1]]
        g2 = [name[1] name[2]; name[2] name[1]]

        new(a,b,g1,g2)
            
    end

end

mutable struct lattice_ind2x1
    A::Index
    B::Index
    gg::Matrix{Symbol}
    gt::Matrix{Symbol}

    function lattice_ind2x1()
        
        name = fieldnames(lattice_ind2x1)
        g1 = [name[1] name[2]; name[2] name[1]]
        g2 = [name[1] name[2]; name[2] name[1]]

        new(ITensor(),ITensor(),g1,g2)
    end

    function lattice_ind2x1(a::Index,b::Index)

        name = fieldnames(lattice2x1)
        g1 = [name[1] name[2]; name[2] name[1]]
        g2 = [name[1] name[2]; name[2] name[1]]

        new(a,b,g1,g2)

    end

    function lattice_ind2x1(nature_of_the_legs::String)

        if nature_of_the_legs == "physical"
            new(
            Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"ia"),Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"ib"))
        elseif nature_of_the_legs == "ancilla"
            new(
            Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"sa"),Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"sb"))
        elseif nature_of_the_legs == "imaginary"
            new(
            Index([QN(0)=>1]),Index([QN(0)=>1]))
        end

    end

end


include("SSMHamiltonian.jl")
include("simple_update_evolution/simpleUpdateJ1.jl")
include("get_tens.jl")
include("initialisation.jl")
include("simpleupdate.jl")


function cool_system_down(D::Int64,J1::Float64,J2::Float64,h::Float64,hs::Float64,temperature) 
    
    final_temp = temperature[end]-1e-5
    N = 2
    dbetasu = 1e-2;
    temperature_save = simpleupdate(final_temp,temperature,dbetasu,J1,J2,h,hs,N,D)

    return temperature_save
  
end






