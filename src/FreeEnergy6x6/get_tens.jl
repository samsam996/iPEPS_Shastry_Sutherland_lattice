

# include("simple_update_evolution/diag_sqrt.jl")

function get_tens2(Gamma,lambdax,lambday,physical_legs,ancilla_legs,gt)

    N = size(gt)[1]
    f(x) = mod(x-1,N) + 1

    tens_a = lattice()
    cx = lattice()
    cy = lattice()

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
    aa1 = getproperty(tens_a,gt[f(i),f(j)])
    aa2 = getproperty(tens_a,gt[f(i+1),f(j)])
    index1 = commoninds(aa1,aa2) 
    index2 = commoninds(aa2,aa1) 

    for i = 1:N
        for j = 1:1
            # @show i
            a1 = getproperty(tens_a,gt[f(i),f(j)])
            a2 = getproperty(tens_a,gt[f(i+1),f(j)])
            cx11 = combiner(index1[3],index1[7], dir = -dir(index1[3]), tags = "xa") # 1,1
            cx22 = combiner(index2[3],index2[7], dir = -dir(index2[3]), tags = "xb") # 2,1

            if abs.(mod(i,2) - 1) < 1e-5
                a1 = a1*(cx11)#*dag(cx22)
                a2 = a2*dag(cx11)#*(cx22)
                setproperty!(cx,gt[i,j],cx11)
            else 
                a1 = a1*(cx22)
                a2 = a2*dag(cx22)
                setproperty!(cx,gt[i,j],cx22)
            end

            setproperty!(tens_a,gt[i,j],a1)
            setproperty!(tens_a,gt[f(i+1),f(j)],a2)
        
        end
    end
 

    for i = 1:1
        for j = 1:N

            a1 = getproperty(tens_a,gt[i,j])
            a2 = getproperty(tens_a,gt[f(i),f(j+1)])

            cy1 = combiner(index1[1],index1[5], dir = -dir(index1[1])) #1,1
            cy2 = combiner(index2[1],index2[5], dir = -dir(index2[1])) #1,2
        
            if abs.(mod(j,2) - 1) < 1e-5
            a1 = a1*(cy1)
            a2 = a2*dag(cy1)
            setproperty!(cy,gt[i,j],cy1)
            else 
            a1 = a1*(cy2)
            a2 = a2*dag(cy2)
            setproperty!(cy,gt[i,j],cy2)
            end
            setproperty!(tens_a,gt[i,j],a1)
            setproperty!(tens_a,gt[f(i),f(f(j+1))],a2)
    
        end
    end

 
    return tens_a, tens_A, cx, cy

end