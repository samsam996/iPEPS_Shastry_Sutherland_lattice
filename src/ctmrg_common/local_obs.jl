

function local_obs(C,T,tens_a,tens_A,cxd,
    cyd,gt::Matrix{Symbol},physical_legs,ancilla_legs,i::Int64,j::Int64,obs::Matrix{ComplexF64})

    N = size(gt)[1]
    f(x) = mod(x-1,N)+1;

    Adown = getproperty(tens_A, gt[f(i),f(j)])

    ia = getproperty(physical_legs, gt[f(i),f(j)])
    sa = getproperty(ancilla_legs, gt[f(i),f(j)])

    all_but_sa = noncommoninds(inds(Adown), [sa])
    Aup = prime(Adown,all_but_sa)

    # c1 = [-4 -2 -2 -2 -2 0 0 0 0 0 0 2 2 2 2 4]; 
    # c2 = [-4 -2 -2 -2 -2 0 0 0 0 0 0 2 2 2 2 4];

    c1 = [-2 0 0 2];
    c2 = [-2 0 0 2];

    ssz = ITensor(dag(ia),ia')

    for i1 = 1:4
        for i2 = 1:4
            if abs(c1[i1] - c2[i2]) < 1e-9
            ssz[ia=>i1,ia'=>i2] = obs[i1,i2]
            end
        end
    end

    AA = Aup*dag(ssz)*dag(Adown)

    comb1 = getproperty(cxd, gt[f(i-1),f(j)])
    comb2 = getproperty(cyd, gt[f(i),f(j)])
    comb3 = getproperty(cxd, gt[f(i),f(j)])
    comb4 = getproperty(cyd, gt[f(i),f(j-1)])

    Obs = (((AA*dag(comb1))*comb2)*comb3)*dag(comb4)

    c1 = getproperty(C[1],gt[f(i-1),f(j-1)])
    t1 = getproperty(T[1],gt[f(i),f(j-1)])
    c2 = getproperty(C[2],gt[f(i+1),f(j-1)])
    t4 = getproperty(T[4],gt[f(i-1),f(j)])
    aa = getproperty(tens_a,gt[f(i),f(j)])

    t2 = getproperty(T[2],gt[f(i+1),f(j)])
    c4 = getproperty(C[4],gt[f(i-1),f(j+1)])
    t3 = getproperty(T[3],gt[f(i),f(j+1)])
    c3 = getproperty(C[3],gt[f(i+1),f(j+1)])

    mm = ((((((((c1*t1)*t4)*Obs)*c2)*t2)*c3)*t3)*c4)[]/((((((((c1*t1)*t4)*aa)*c2)*t2)*c3)*t3)*c4)[]
   
    return mm

end