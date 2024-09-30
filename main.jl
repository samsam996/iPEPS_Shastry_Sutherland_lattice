
using MKL, ITensors, FileIO, JLD2, KrylovKit, LinearAlgebra, PyPlot

BLAS.set_num_threads(1)
ITensors.Strided.set_num_threads(1)
ITensors.disable_threaded_blocksparse()

mutable struct lattice
    A::ITensor
    B::ITensor
    C::ITensor
    D::ITensor
    E::ITensor
    F::ITensor
    gg::Matrix{Symbol}
    gt::Matrix{Symbol}
    
    function lattice()
        name = fieldnames(lattice)
        g1 = [
        name[1] name[2] name[3] name[4] name[5] name[6]; 
        name[6] name[1] name[2] name[3] name[4] name[5]; 
        name[5] name[6] name[1] name[2] name[3] name[4]; 
        name[4] name[5] name[6] name[1] name[2] name[3]; 
        name[3] name[4] name[5] name[6] name[1] name[2]; 
        name[2] name[3] name[4] name[5] name[6] name[1]];

        g2 = [
        name[1] name[6] name[5] name[4] name[3] name[2]; 
        name[2] name[1] name[6] name[5] name[4] name[3]; 
        name[3] name[2] name[1] name[6] name[5] name[4]; 
        name[4] name[3] name[2] name[1] name[6] name[5]; 
        name[5] name[4] name[3] name[2] name[1] name[6]; 
        name[6] name[5] name[4] name[3] name[2] name[1]]; 

        new(
        ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),g1,g2)
    end
    
    function lattice(a::ITensor,b::ITensor,c::ITensor,d::ITensor,e::ITensor,f::ITensor)

        name = fieldnames(lattice)
        g1 = [
        name[1] name[2] name[3] name[4] name[5] name[6]; 
        name[6] name[1] name[2] name[3] name[4] name[5]; 
        name[5] name[6] name[1] name[2] name[3] name[4]; 
        name[4] name[5] name[6] name[1] name[2] name[3]; 
        name[3] name[4] name[5] name[6] name[1] name[2]; 
        name[2] name[3] name[4] name[5] name[6] name[1]];

        g2 = [
        name[1] name[6] name[5] name[4] name[3] name[2]; 
        name[2] name[1] name[6] name[5] name[4] name[3]; 
        name[3] name[2] name[1] name[6] name[5] name[4]; 
        name[4] name[3] name[2] name[1] name[6] name[5]; 
        name[5] name[4] name[3] name[2] name[1] name[6]; 
        name[6] name[5] name[4] name[3] name[2] name[1]]; 

        new(a,b,c,d,e,f,g1,g2)
            
    end

end

mutable struct lattice_ind
    A::Index
    B::Index
    C::Index
    D::Index
    E::Index
    F::Index
    gg::Matrix{Symbol}
    gt::Matrix{Symbol}

    function lattice_ind()
        
        name = fieldnames(lattice_ind)
        g1 = [
        name[1] name[2] name[3] name[4] name[5] name[6]; 
        name[6] name[1] name[2] name[3] name[4] name[5]; 
        name[5] name[6] name[1] name[2] name[3] name[4]; 
        name[4] name[5] name[6] name[1] name[2] name[3]; 
        name[3] name[4] name[5] name[6] name[1] name[2]; 
        name[2] name[3] name[4] name[5] name[6] name[1]];

        g2 = [
        name[1] name[6] name[5] name[4] name[3] name[2]; 
        name[2] name[1] name[6] name[5] name[4] name[3]; 
        name[3] name[2] name[1] name[6] name[5] name[4]; 
        name[4] name[3] name[2] name[1] name[6] name[5]; 
        name[5] name[4] name[3] name[2] name[1] name[6]; 
        name[6] name[5] name[4] name[3] name[2] name[1]]; 

        new(ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),g1,g2)
    end

    function lattice_ind(a::Index,b::Index,c::Index,d::Index,e::Index,f::Index)

        name = fieldnames(lattice)
        g1 = [
        name[1] name[2] name[3] name[4] name[5] name[6]; 
        name[6] name[1] name[2] name[3] name[4] name[5]; 
        name[5] name[6] name[1] name[2] name[3] name[4]; 
        name[4] name[5] name[6] name[1] name[2] name[3]; 
        name[3] name[4] name[5] name[6] name[1] name[2]; 
        name[2] name[3] name[4] name[5] name[6] name[1]];

        g2 = [
        name[1] name[6] name[5] name[4] name[3] name[2]; 
        name[2] name[1] name[6] name[5] name[4] name[3]; 
        name[3] name[2] name[1] name[6] name[5] name[4]; 
        name[4] name[3] name[2] name[1] name[6] name[5]; 
        name[5] name[4] name[3] name[2] name[1] name[6]; 
        name[6] name[5] name[4] name[3] name[2] name[1]]; 
          
        new(a,b,c,d,e,f,g1,g2)

    end

    function lattice_ind(nature_of_the_legs::String)

        if nature_of_the_legs == "physical"
            new(
            Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"ia"),Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"ib"),
            Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"ic"),Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"id"),
            Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"ie"),Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"if"))
        elseif nature_of_the_legs == "ancilla"
            new(
            Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"sa"),Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"sb"),
            Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"sc"),Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"sd"),
            Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"se"),Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"sf"))
        elseif nature_of_the_legs == "imaginary"
            new(
            Index([QN(0)=>1]),Index([QN(0)=>1]),Index([QN(0)=>1]),Index([QN(0)=>1]),Index([QN(0)=>1]),Index([QN(0)=>1]))
        end

    end

end


function main(D, J1, J2, h, hs, temperature, unit_cell) 

    if unit_cell != "2x2" && unit_cell != "6x6"
        print("ERROR: your unit cell should either be 2x2 or 6x6")
        return
    end

    if abs(hs) != 0 
        print("ERROR: THERE IS NOT YET THE POSSIBILITY OF COOLING DOWN WITH A 6x6")
        return
    end

    include("src/observables/correlation_length.jl")
    include("src/observables/magnetisation.jl")
    include("src/observables/energy.jl")
    include("src/observables/z2order_parameter.jl")
    include("src/observables/PartitionPerSite.jl")
    include("src/CoolDown2x2/main_simpleupdate.jl")

    include("src/FreeEnergy"*unit_cell*"/compute_FreeEnergy"*unit_cell*".jl")

    temperature_save = cool_system_down(D,J1,J2,h,hs,temperature)

    file_names = []
    for temp in temperature_save

        file_name = string("Results/LocalTensorsD$(D)temp",string(round(temp,digits = 5)),"J$(J1)h",string(round(1e2*h,digits = 3)),".jld2")
        push!(file_names, file_name) 
    end

    k = get_free_energy(file_names)

    scatter(temperature_save, -temperature_save.*k,  color="blue")
    grid(true)
    show()
    savefig("figures/free_energy.png") 
end


let 
    D::Int64 = 3;  
    J1::Float64 = 0.63; 
    J2::Float64 = 1; 
    h::Float64 = 1.; 
    hs::Float64 = 0; 
    temperature = LinRange(5.22, 2.1, 2)
    unit_cell = "2x2"

    main(D,J1,J2,h,hs,temperature,unit_cell)

end

