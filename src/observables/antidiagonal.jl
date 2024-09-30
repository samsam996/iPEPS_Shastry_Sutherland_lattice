

function antidiagonal(C,T,gt,tens_a,tens_A,cxd,cyd,physical_legs,entanglement_legs,i,j)

    N = size(gt)[1]
    f(x) = mod(x-1,N) + 1
    
    szz = 1/2*[1 0; 0 -1]; id = [1 0; 0 1];
    obs::Matrix{Float64} =  kron(szz,id) + kron(id,szz);


    C1a = getproperty(C[1],gt[f(i-1),f(j-1)])
    T1b = getproperty(T[1],gt[f(i+0),f(j-1)])
    T1c = getproperty(T[1],gt[f(i+1),f(j-1)])
    C2d = getproperty(C[2],gt[f(i+2),f(j-1)])

    T4d = getproperty(T[4],gt[f(i-1),f(j)])
    a = getproperty(tens_a,gt[f(i+0),f(j)])
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

    Bdown = getproperty(tens_A, gt[f(i+1),f(j)])
    sb = getproperty(entanglement_legs, gt[f(i+1),f(j)])
    ib = getproperty(physical_legs, gt[f(i+1),f(j)])

    Ddown = getproperty(tens_A, gt[f(i),f(j+1)])
    sd = getproperty(entanglement_legs, gt[f(i),f(j+1)])
    id = getproperty(physical_legs, gt[f(i),f(j+1)])

    all_ind_prime = noncommoninds(inds(Bdown),[sb])
    Bup = prime(Bdown, all_ind_prime)

    all_ind_prime = noncommoninds(inds(Ddown),[sd])
    Dup = prime(Ddown, all_ind_prime)

    c1 = [-2 0 0 2]; c2 = [-2 0 0 2];
    szb = ITensor(dag(ib),ib')
    szd = ITensor(dag(id),id')
    for i1 = 1:4
        for i2 = 1:4
            if abs(c1[i1] - c2[i2]) < 1e-9
            szb[ib=>i1,ib'=>i2] = obs[i1,i2]
            szd[id=>i1,id'=>i2] = obs[i1,i2]
            end
        end
    end

    BB = Bup*dag(szb)*dag(Bdown)
    DD = Dup*dag(szd)*dag(Ddown)

    comb1b = getproperty(cxd, gt[f(i),f(j)])
    comb2b = getproperty(cyd, gt[f(i+1),f(j)])
    comb3b = getproperty(cxd, gt[f(i+1),f(j)])
    comb4b = getproperty(cyd, gt[f(i+1),f(j-1)])
    Obsb = (((BB*dag(comb1b))*comb2b)*comb3b)*dag(comb4b)
    comb1d = getproperty(cxd, gt[f(i-1),f(j+1)])
    comb2d = getproperty(cyd, gt[f(i),f(j+1)])
    comb3d = getproperty(cxd, gt[f(i),f(j+1)])
    comb4d = getproperty(cyd, gt[f(i),f(j)])
    Obsd = (((DD*dag(comb1d))*comb2d)*comb3d)*dag(comb4d)

    M1u = (((C1a*T1b)*T4d)*a)
    M2u = (((C2d*T1c)*T2c)*Obsb)
    M3u = (((C3a*T2b)*T3d)*a)
    M4u = (((C4b*T3c)*T4c)*Obsd)


    Zup = ((M1u*M2u)*M4u)*M3u

    return Zup[1]/Z[1]



end