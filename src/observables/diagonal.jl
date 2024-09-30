


function diagonal(C,T,gt,tens_a,tens_A,cxd,cyd,physical_legs,entanglement_legs,i,j)

    N = size(gt)[1]
    f(x) = mod(x-1,N) + 1

    szz = 1/2*[1 0; 0 -1]; id = [1 0; 0 1];
    obs::Matrix{Float64} =  kron(szz,id) + kron(id,szz);

    C1a = getproperty(C[1],gt[f(i-1),f(j-1)])
    T1b = getproperty(T[1],gt[f(i+0),f(j-1)])
    T1c = getproperty(T[1],gt[f(i+1),f(j-1)])
    C2d = getproperty(C[2],gt[f(i+2),f(j-1)])

    T4d = getproperty(T[4],gt[f(i-1),f(j)])
    a = getproperty(tens_a,gt[f(i),f(j)])
    b = getproperty(tens_a,gt[f(i+1),f(j)])
    T2c = getproperty(T[2],gt[f(i+2),f(j)])

    T4c = getproperty(T[4],gt[f(i-1),f(j+1)])
    d = getproperty(tens_a,gt[f(i+0),f(j+1)])
    T2b = getproperty(T[2],gt[f(i+2),f(j+1)])

    C4b = getproperty(C[4],gt[f(i-1),f(j+2)])
    T3c = getproperty(T[3],gt[f(i+0),f(j+2)])
    T3d = getproperty(T[3],gt[f(i+1),f(j+2)])
    C3a = getproperty(C[3],gt[f(i+2),f(j+2)])

    M1d = (((C1a*T1b)*T4d)*a)
    M2d = (((C2d*T1c)*T2c)*b)
    M3d = (((C3a*T2b)*T3d)*a)
    M4d = (((C4b*T3c)*T4c)*d)

    Z = ((M1d*M2d)*M4d)*M3d

    A1down = getproperty(tens_A, gt[f(i),f(j)])
    sa = getproperty(entanglement_legs, gt[f(i),f(j)])
    ia = getproperty(physical_legs, gt[f(i),f(j)])

    A2down = getproperty(tens_A, gt[f(i+1),f(j+1)])

    all_ind_prime = noncommoninds(inds(A1down),[sa])
    A1up = prime(A1down, all_ind_prime)

    all_ind_prime = noncommoninds(inds(A2down),[sa])
    A2up = prime(A2down, all_ind_prime)

    c1 = [-2 0 0 2]; c2 = [-2 0 0 2];
    sza = ITensor(dag(ia),ia')
    for i1 = 1:4
        for i2 = 1:4
            if abs(c1[i1] - c2[i2]) < 1e-9
            sza[ia=>i1,ia'=>i2] = obs[i1,i2]
            end
        end
    end

    AA1 = A1up*dag(sza)*dag(A1down)
    AA2 = A2up*dag(sza)*dag(A2down)

    # pour i,j
    # comb1 = getproperty(cxd, gt[f(i-1),f(j)])
    # comb2 = getproperty(cyd, gt[f(i),f(j)])
    # comb3 = getproperty(cxd, gt[f(i),f(j)])
    # comb4 = getproperty(cyd, gt[f(i),f(j-1)])

    comb1a = getproperty(cxd, gt[f(i-1),f(j)])
    comb2a = getproperty(cyd, gt[f(i),f(j)])
    comb3a = getproperty(cxd, gt[f(i),f(j)])
    comb4a = getproperty(cyd, gt[f(i),f(j-1)])

    Obsa1 = (((AA1*dag(comb1a))*comb2a)*comb3a)*dag(comb4a)
    Obsa2 = (((AA2*dag(comb1a))*comb2a)*comb3a)*dag(comb4a)

    M1u = (((C1a*T1b)*T4d)*Obsa1)
    M2u = (((C2d*T1c)*T2c)*b)
    M3u = (((C3a*T2b)*T3d)*Obsa2)
    M4u = (((C4b*T3c)*T4c)*d)


    Zup = ((M1u*M2u)*M4u)*M3u

    return Zup[1]/Z[1]


end