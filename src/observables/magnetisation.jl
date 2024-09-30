

"""
    ind_x and ind_y are indices of tens_A 
    c_x_a and c_y_a are what maps the double indices into the indices of tens_a
    
"""


# function magnetisation(C::Vector{lattice},T::Vector{lattice},
#     tens_a::lattice,tens_A::lattice,gt::Matrix{Symbol},cxd::lattice,
#     cyd::lattice,physical_legs::lattice_ind,entanglement_legs::lattice_ind)
function magnetisation(C,T,tens_a,tens_A,gt::Matrix{Symbol},cxd,cyd,physical_legs,entanglement_legs)

    N = size(gt)[1]
    f(x) = mod(x-1,N) + 1

    mm = zeros(N,N)

    magne = lattice()

    for i = 1:N
        for j = 1:N

        Adown = getproperty(tens_A, gt[f(i),f(j)])

        sa = getproperty(entanglement_legs, gt[f(i),f(j)])
        ia = getproperty(physical_legs, gt[f(i),f(j)])
        
        all_ind_prime = noncommoninds(inds(Adown),[sa])
        Aup = prime(Adown, all_ind_prime)

        sz = 1/2*[1 0; 0 -1]; id = [1 0; 0 1]
        c1 = [-2 0 0 2]; c2 = [-2 0 0 2];

        obs = kron(sz,id) + kron(id,sz)
        szz = ITensor(dag(ia),ia')

        for i1 = 1:4
            for i2 = 1:4
                if abs(c1[i1] - c2[i2]) < 1e-9
                szz[ia=>i1,ia'=>i2] = obs[i1,i2]
                end
            end
        end

        AA = Aup*dag(szz)*dag(Adown)

        comb1 = getproperty(cxd, gt[f(i-1),f(j)])
        comb2 = getproperty(cyd, gt[f(i),f(j)])
        comb3 = getproperty(cxd, gt[f(i),f(j)])
        comb4 = getproperty(cyd, gt[f(i),f(j-1)])
    
        AAA = AA*dag(comb1)*comb2*comb3*dag(comb4)
    
        setproperty!(magne,gt[f(i),f(j)],AAA)
        
        end
    end


    for i = 1:N
        for j = 1:N

            c1 = getproperty(C[1],gt[f(i-1),f(j-1)])
            t1 = getproperty(T[1],gt[f(i),f(j-1)])
            c2 = getproperty(C[2],gt[f(i+1),f(j-1)])
            t4 = getproperty(T[4],gt[f(i-1),f(j)])
            aa = getproperty(tens_a,gt[f(i),f(j)])
            local_mm = getproperty(magne,gt[f(i),f(j)])
            t2 = getproperty(T[2],gt[f(i+1),f(j)])
            c4 = getproperty(C[4],gt[f(i-1),f(j+1)])
            t3 = getproperty(T[3],gt[f(i),f(j+1)])
            c3 = getproperty(C[3],gt[f(i+1),f(j+1)])

            mm[i,j] = (((((((((c1*t1)*t4)*local_mm)*c2)*t2)*c3)*t3)*c4)[1])/((((((((((c1*t1)*t4)*aa)*c2)*t2)*c3)*t3)*c4))[1])


        end
    end


    # then contract 

    return mm

end