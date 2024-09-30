

include("simple_update_evolution/diag_sqrt.jl")

function get_tens(Gamma,lambdax,lambday,physical_legs,ancilla_legs,gt)

    N = size(gt)[1]
    f(x) = mod(x-1,N) + 1

    tens_a = lattice2x1()
    cx = lattice2x1()
    cy = lattice2x1()

    tens_A = deepcopy(Gamma)
    
    for i = 1:N
        for j = 1
    
            Gamma0 = getproperty(tens_A,gt[f(i),f(j)])
            Gamma1 = getproperty(tens_A,gt[f(i+1),f(j)])

            lambdax1 = getproperty(lambdax,gt[f(i),f(j)])
            ind1 = commonind(Gamma0,lambdax1)
            ind2 = commonind(lambdax1,Gamma1)

            Gamma0 = Gamma0*diag_sqrt(lambdax1)
            Gamma0 = Gamma0*delta((ind1),dag(ind2))
            Gamma1 = diag_sqrt(lambdax1)*Gamma1
            setproperty!(tens_A,gt[f(i),f(j)],Gamma0)
            setproperty!(tens_A,gt[f(i+1),f(j)],Gamma1)

        end
    end

    for i = 1:N
        for j = 1
    
            Gamma0 = getproperty(tens_A,gt[f(i),f(j)])
            Gamma1 = getproperty(tens_A,gt[f(i),f(j+1)])

            lambday1 = getproperty(lambday,gt[f(i),f(j)])
            ind1 = commonind(Gamma0,lambday1)
            ind2 = commonind(lambday1,Gamma1)

            Gamma0 = Gamma0*diag_sqrt(lambday1)
            Gamma0 =Gamma0*delta((ind1),dag(ind2))
            Gamma1 = diag_sqrt(lambday1)*Gamma1
            setproperty!(tens_A,gt[f(i),f(j)],Gamma0)
            setproperty!(tens_A,gt[f(i),f(j+1)],Gamma1)

        end
    end

    for i = 1:N
        for j = 1

            ia = getproperty(physical_legs,gt[i,j])
            sa = getproperty(ancilla_legs, gt[i,j])
            A = getproperty(tens_A,gt[i,j])
            indA = inds(A)
            ind_dgtb = noncommoninds(indA, [ia,sa])
            A_prime = prime(A,ind_dgtb)
            a = A_prime*dag(A)
            setproperty!(tens_a,gt[i,j],a)

        end
    end


    i = 1; j = 1
    @show gt
    a1 = getproperty(tens_a,gt[f(i),f(j)])
    a2 = getproperty(tens_a,gt[f(i+1),f(j)])

    index1 = commoninds(a1,a2) 
    index2 = commoninds(a2,a1) 

    cx11 = combiner(index1[3],index1[7], dir = -dir(index1[3]), tags = "xa") # 1,1
    cx22 = combiner(index2[3],index2[7], dir = -dir(index2[3]), tags = "xb") # 2,1

    a1 = a1*cx11*dag(cx22)
    a2 = a2*dag(cx11)*(cx22)
    setproperty!(tens_a,gt[i,j],a1)
    setproperty!(tens_a,gt[f(i+1),f(j)],a2)
    setproperty!(cx,gt[2,j],cx22)
    setproperty!(cx,gt[1,j],cx11)

    a1 = getproperty(tens_a,gt[i,j])
    a2 = getproperty(tens_a,gt[f(i),f(j+1)])
    index1 = commoninds(a1,a2)
    index2 = commoninds(a2,a1)

    cy1 = combiner(index1[4],index1[6], dir = -dir(index1[4])) #1,1
    cy2 = combiner(index2[4],index2[6], dir = -dir(index2[4])) #1,2

    a1 = a1*(cy1)*dag(cy2)
    a2 = a2*dag(cy1)*(cy2)
    setproperty!(tens_a,gt[i,j],a1)
    setproperty!(tens_a,gt[f(i),f(j+1)],a2)
    setproperty!(cy,gt[i,1],cy1)
    setproperty!(cy,gt[i,2],cy2)

   

 
    return tens_a, tens_A, cx, cy

end