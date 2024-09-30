
include("antidiagonal.jl")
include("diagonal.jl")


function z2order_parameter(C,T,gt,tens_a,tens_A,cxd,cyd,physical_legs,entanglement_legs)



    i = 1; j = 1; 

    diag = diagonal(C,T,gt,tens_a,tens_A,cxd,cyd,physical_legs,entanglement_legs,i,j)
    antidiag = antidiagonal(C,T,gt,tens_a,tens_A,cxd,cyd,physical_legs,entanglement_legs,i+1,j)

    # still missing the other one.
    @show diag
    @show antidiag

    return diag - antidiag

end