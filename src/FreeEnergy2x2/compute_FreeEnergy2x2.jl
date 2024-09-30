

include("get_tens2.jl")
include("ctmrg_evolution/ctm.jl")
include("ctmrg_evolution/OBC_PEPS.jl")



function get_free_energy(file_names)


file_name =[] 
magne = []
energie = []
tempe = []
invcorr1 = []
invcorr2 = []
invcorr3 = []
invcorr4 = []
wave_vec1 = []
wave_vec2 = []  
wave_vec3 = []
wave_vec4 = []
eig3 = []
eig4 = []
numb_iter = []
order_parameter_z2 = []
PPS = []



for file_name in file_names            

    D,J1,J2,h,hs,nsu,temp = load(file_name,"D","J1","J2","h","hs","nsu","temp")

    Gamma, lambdax, lambday, gt, FreeEnergy, physical_legs, ancilla_legs = load(file_name,"Gamma","lambdax","lambday","gt","FreeEnergy", "physical_legs","ancilla_legs")


    chi = D*D + 1; 
    chi = 10
    precision_ctm = 1e-7
    tens_a,tens_A,cxd,cyd = get_tens(Gamma,lambdax,lambday,physical_legs,ancilla_legs,gt)
    C,T = OBC_PEPS(tens_A,cxd,cyd,gt)

    C, T, it, err = ctm!(tens_a,tens_A,cxd,cyd,gt,physical_legs,ancilla_legs,chi,precision_ctm,C,T) 
    mm = magnetisation(C,T,tens_a,tens_A,gt,cxd,cyd,physical_legs,ancilla_legs)

    PPSf = 1/2*(PartitionPerSite(T,C,tens_a,gt,1,1) ) + 2*FreeEnergy

    omega0 = 1;
    omega1 = exp(2*pi*im/3)
    omega2 = exp(4*pi*im/3);

    ener = energy(C,T,tens_a,tens_A,gt,cxd,cyd,physical_legs,ancilla_legs,J1,J2,h,hs)

    order_parameter = 1/2*abs(omega0*mm[1,1] - mm[1,2])

    z2 = z2order_parameter(C,T,gt,tens_a,tens_A,cxd,cyd,physical_legs,ancilla_legs)
    @show push!(order_parameter_z2,z2)
    @show push!(magne,order_parameter)
    @show push!(tempe, temp)
    @show push!(PPS,PPSf)

    xi2, xi3, xi4, dq = correlation_length(C,T,gt)
                
    push!(invcorr1, xi2[1,1])
    push!(invcorr2, xi2[1,2])
    push!(invcorr3, xi2[2,1])
    push!(invcorr4, xi2[2,2])
    push!(wave_vec1, dq[1,1])
    push!(wave_vec2, dq[1,2])
    push!(wave_vec3, dq[2,1])
    push!(wave_vec4, dq[2,2])
    push!(eig3, xi3[1,1])
    push!(eig4, xi4[1,1])
    push!(numb_iter, it)



    name_data = string("Results/U1_Results2x2D",string(D),"temp",string(round(temp,digits = 5)),"_nsu_",string(nsu),"hs",string(round(1e4*hs,digits = 3)),"h",string(round(1e2*h,digits = 3)),".jld2")
    save(name_data,
    "D",D,
    "magne",order_parameter,
    "J2",J2,
    "h",h,
    "J1",J1,
    "hs",hs,
    "nsu",nsu,
    "energie",ener,
    "z2",z2,
    "tempe",temp,
    "invcorr1",xi2[1,1],
    "invcorr2",xi2[1,2],
    "invcorr3",xi2[2,1],
    "invcorr4",xi2[2,2],
    "wave_vec1",dq[1,1],
    "wave_vec2",dq[1,2],
    "wave_vec3",dq[2,1],
    "wave_vec4",dq[2,2],
    "magne_z", mm,
    "eig3",xi3[1,1],
    "eig4",xi4[1,1],
    "numb_iter",it,
    "chi",chi,
    "err_ctm",err,
    "C",C,
    "T",T,
    "tens_A",tens_A,
    "cxd",cxd,
    "cyd",cyd,
    "PPS",PPSf,
"gt",gt)


end


return PPS




 


end
