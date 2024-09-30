


include("ctmrg_evolution/OBC_PEPS.jl")
include("ctmrg_evolution/OBC_PEPS_ONES.jl")
include("ctmrg_evolution/ctm.jl")
include("get_tens.jl")



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
    PPS1x6 = []


    for file_name in file_names
        
        D,J1,J2,h,hs,nsu,temp = load(file_name,"D","J1","J2","h","hs","nsu","temp")

        Gamma_old, lambdax_old, lambday_old, gt_old, FreeEnergy, physical_legs_old, ancilla_legs_old = load(file_name,"Gamma","lambdax","lambday","gt","FreeEnergy", "physical_legs","ancilla_legs")

        gt = deepcopy(gt_old)

        Gamma = lattice(getproperty(Gamma_old,gt_old[1,1]),getproperty(Gamma_old,gt_old[1,2]),getproperty(Gamma_old,gt_old[1,1]),
        getproperty(Gamma_old,gt_old[1,2]),getproperty(Gamma_old,gt_old[1,1]),getproperty(Gamma_old,gt_old[1,2]))

        lambdax = lattice(getproperty(lambdax_old,gt_old[1,1]),getproperty(lambdax_old,gt_old[1,2]),getproperty(lambdax_old,gt_old[1,1]),
        getproperty(lambdax_old,gt_old[1,2]),getproperty(lambdax_old,gt_old[1,1]),getproperty(lambdax_old,gt_old[1,2]))

        lambday = lattice(getproperty(lambday_old,gt_old[1,1]),getproperty(lambday_old,gt_old[1,2]),getproperty(lambday_old,gt_old[1,1]),
        getproperty(lambday_old,gt_old[1,2]),getproperty(lambday_old,gt_old[1,1]),getproperty(lambday_old,gt_old[1,2]))

        physical_legs = lattice_ind(getproperty(physical_legs_old,gt[1,1]),getproperty(physical_legs_old,gt[1,2]),
        getproperty(physical_legs_old,gt_old[1,1]),getproperty(physical_legs_old,gt[1,2]),getproperty(physical_legs_old,gt[1,1]),
        getproperty(physical_legs_old,gt[1,2]))

        ancilla_legs = lattice_ind(getproperty(ancilla_legs_old,gt[1,1]),getproperty(ancilla_legs_old,gt[1,2]),
        getproperty(ancilla_legs_old,gt[1,1]),getproperty(ancilla_legs_old,gt[1,2]),getproperty(ancilla_legs_old,gt[1,1]),
        getproperty(ancilla_legs_old,gt[1,2]))


        gt = Gamma.gt
        chi = D*D + 1; 
        chi = 10
        precision_ctm = 1e-7
        tens_a,tens_A,cxd,cyd = get_tens2(Gamma,lambdax,lambday,physical_legs,ancilla_legs,gt)
        
        # C,T = OBC_PEPS_ONES(tens_A,cxd,cyd,gt,physical_legs)
        C,T = OBC_PEPS_ONES(tens_A,cxd,cyd,gt,physical_legs)

        C, T, it, err = ctm!(tens_a,tens_A,cxd,cyd,gt,physical_legs,ancilla_legs,chi,precision_ctm,C,T) 
        mm = magnetisation(C,T,tens_a,tens_A,gt,cxd,cyd,physical_legs,ancilla_legs)

        PPSf = 1/6*(PartitionPerSite(T,C,tens_a,gt,1,1) + PartitionPerSite(T,C,tens_a,gt,3,1) + PartitionPerSite(T,C,tens_a,gt,5,1)) + 
        2*FreeEnergy

        PPSf1x6 = 1/12*(PartitionPerSite(T,C,tens_a,gt,1,1) + PartitionPerSite(T,C,tens_a,gt,2,1) + PartitionPerSite(T,C,tens_a,gt,3,1) +
                PartitionPerSite(T,C,tens_a,gt,4,1) + PartitionPerSite(T,C,tens_a,gt,5,1) + PartitionPerSite(T,C,tens_a,gt,6,1)) + 2*FreeEnergy

        sx = 1/2*[0 1; 1 0]; sy = 1/2*[0 -1im; 1im 0]; sz = 1/2*[1 0; 0 -1]; id = [1 0; 0 1];
        obs1::Matrix{ComplexF64} = (kron(sz,sz) + kron(sx,sx) + kron(sy,sy)) 
        kk1 = local_obs(C,T,tens_a,tens_A,cxd,cyd,gt,physical_legs,ancilla_legs,1,1,obs1)
        kk2 = local_obs(C,T,tens_a,tens_A,cxd,cyd,gt,physical_legs,ancilla_legs,1,2,obs1)
        kk3 = local_obs(C,T,tens_a,tens_A,cxd,cyd,gt,physical_legs,ancilla_legs,1,3,obs1)
        kk4 = local_obs(C,T,tens_a,tens_A,cxd,cyd,gt,physical_legs,ancilla_legs,1,4,obs1)
        kk5 = local_obs(C,T,tens_a,tens_A,cxd,cyd,gt,physical_legs,ancilla_legs,1,5,obs1)
        kk6 = local_obs(C,T,tens_a,tens_A,cxd,cyd,gt,physical_legs,ancilla_legs,1,6,obs1)

        ener = energy(C,T,tens_a,tens_A,gt,cxd,cyd,physical_legs,ancilla_legs,J1,J2,h,hs)
        omega0 = 1;
        omega1 = exp(2*pi*im/3)
        omega2 = exp(4*pi*im/3);
        println([mm[1,1],mm[1,2],mm[1,3],mm[1,4],mm[1,5],mm[1,6]])

        order_parameter = 1/6*abs(omega0*mm[1,1]+omega1*mm[1,2]+omega2*mm[1,3]+omega0*mm[1,4]+omega1*mm[1,5]+omega2*mm[1,6])

        z2 = z2order_parameter(C,T,gt,tens_a,tens_A,cxd,cyd,physical_legs,ancilla_legs)
        @show push!(order_parameter_z2,z2)
        @show push!(magne,order_parameter)
        @show push!(tempe, temp)
        @show push!(PPS,PPSf)
        @show push!(PPS1x6,PPSf1x6)
        push!(energie, ener)
        @show ener
    
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


        name_data = string("Results/U1_Results6x6D",string(D),"temp",string(round(temp,digits = 5)),"_nsu_",string(nsu),"hs",string(round(1e4*hs,digits = 3)),"h",string(round(1e2*h,digits = 3)),".jld2")
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
        "PPS1x6",PPSf1x6,
        "ancilla_legs",ancilla_legs,
        "physical_legs",physical_legs,
        "kk1",kk1,
        "kk2",kk2,
        "kk3",kk3,
        "kk4",kk4,
        "kk5",kk5,
        "kk6",kk6,
        "gt",gt)

    end


    return PPS1x6


end
