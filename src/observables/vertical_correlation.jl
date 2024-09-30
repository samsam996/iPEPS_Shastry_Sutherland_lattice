


function vertical_correlation(C,T,tens_a,tens_A,gt::Matrix{Symbol},cxd,cyd,
    physical_legs,entanglement_legs,i::Int64,j::Int64)

    N = size(gt)[1]
    f(x) = mod(x-1,N) + 1
    sx = 1/2*[0 1; 1 0]; sy = 1/2*[0 -1im; 1im 0]; sz = 1/2*[1 0; 0 -1]; id = [1 0; 0 1];
    ey = 0

    ss_xyz = zeros(ComplexF64,4,4,4)

    #       3|
    #        |
    #       4|
    # 
    #     -1----2-

    if mod(i+j,2) == 0 # pair -> up 
        ss_xyz = kron(kron(sx,id),kron(id,sx)) + kron(kron(sy,id),kron(id,sy)) + kron(kron(sz,id),kron(id,sz)) + 
                kron(kron(id,sx),kron(id,sx)) + kron(kron(id,sy),kron(id,sy)) + kron(kron(id,sz),kron(id,sz))
    end

    if mod(i+j,2) == 1 # impair -> down
        ss_xyz = kron(kron(sx,id),kron(sx,id)) + kron(kron(sy,id),kron(sy,id)) + kron(kron(sz,id),kron(sz,id)) + 
                kron(kron(sx,id),kron(id,sx)) + kron(kron(sy,id),kron(id,sy)) + kron(kron(sz,id),kron(id,sz))
    end
    ss_xyz = reshape(ss_xyz,(4,4,4,4))


    Adown = getproperty(tens_A, gt[f(i),f(j)])
    Cdown = getproperty(tens_A, gt[f(i),f(j+1)])

    sa = getproperty(entanglement_legs, gt[f(i),f(j)])
    sc = getproperty(entanglement_legs, gt[f(i),f(j+1)])

    ia = getproperty(physical_legs, gt[f(i),f(j)])
    ic = getproperty(physical_legs, gt[f(i),f(j+1)])

    prime_a = noncommoninds(inds(Adown),[sa])
    prime_c = noncommoninds(inds(Cdown),[sc])
     
    Aup = prime(Adown, prime_a)
    Cup = prime(Cdown, prime_c)

    c1 = [-2,0,0,2]
    Obs_cax = ITensor(dag(ic),dag(ia),ic',ia')

    for i1 = 1:4
        for i2 = 1:4
            for j1 = 1:4
                for j2 = 1:4
                    if abs(c1[i1] + c1[i2] - c1[j1] - c1[j2]) < 1e-9
                        Obs_cax[ic=>i1,ia=>i2,ic'=>j1,ia'=>j2] = ss_xyz[i1,i2,j1,j2]
                    end
                end
            end
        end
    end

    AA = Aup*dag(Adown)
    CC = Cup*dag(Cdown)
            
    comb1a = getproperty(cxd, gt[f(i-1),f(j)])
    comb2a = getproperty(cyd, gt[f(i),f(j)])
    comb3a = getproperty(cxd, gt[f(i),f(j)])
    comb4a = getproperty(cyd, gt[f(i),f(j-1)])

    comb1c = getproperty(cxd, gt[f(i-1),f(j+1)])
    comb2c = getproperty(cyd, gt[f(i),f(j+1)])
    comb3c = getproperty(cxd, gt[f(i),f(j+1)])
    comb4c = getproperty(cyd, gt[f(i),f(j)])

    aa = (((AA*dag(comb1a))*comb2a)*comb3a)*dag(comb4a)
    cc = (((CC*dag(comb1c))*comb2c)*comb3c)*dag(comb4c)
    #   C1a:i-1,j-1    T1b:i+0,j-1   C2a:i+1,j-1   
    #   T4c:i-1,j+0    A1d:i+0,j+0   T2c:i+1,j+0   
    #   T4a:i-1,j+1    A2b:i+0,j+1   T2a:i+1,j+1  
    #   C4c:i-1,j+2    T3d:i+0,j+2   C3c:i+1,j+2


    c1a = getproperty(C[1],gt[f(i-1),f(j-1)])
    t1b = getproperty(T[1],gt[f(i+0),f(j-1)])
    c2a = getproperty(C[2],gt[f(i+1),f(j-1)])

    t4c = getproperty(T[4],gt[f(i-1),f(j)])
    a1 = getproperty(tens_a,gt[f(i),f(j)])
    t2c = getproperty(T[2],gt[f(i+1),f(j)])

    t4a = getproperty(T[4],gt[f(i-1),f(j+1)])
    a3 = getproperty(tens_a,gt[f(i),f(j+1)])
    t2a = getproperty(T[2],gt[f(i+1),f(j+1)])

    c4c = getproperty(C[4],gt[f(i-1),f(j+2)])
    t3d = getproperty(T[3],gt[f(i),f(j+2)])
    c3c = getproperty(C[3],gt[f(i+1),f(j+2)])

    Left = (c1a*t1b)*c2a;
    Right = (c4c*t3d)*c3c;

    num = (((((((((Left)*t4c)*aa)*dag(Obs_cax))*t2c)*t4a)*cc)*t2a)*Right)[1]
    denom = ((((Left)*t4c*a1*t2c)*t4a*a3*t2a)*Right)[1]

    ey = num/denom

    return ey

end