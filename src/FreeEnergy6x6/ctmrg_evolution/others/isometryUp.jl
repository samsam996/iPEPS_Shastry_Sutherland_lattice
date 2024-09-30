




function isometryUp(C,T,tens,gt::Matrix{Symbol},i::Int64,j::Int64,chi::Int64)

    N = size(gt)[1];
    f(x) = mod(x-1,N) + 1;
 
    c1a = (getfield(C[1],gt[f(i-2),f(j-1)]));
    t1b = (getfield(T[1],gt[f(i-1),f(j-1)]));
    t4c =  (getfield(T[4],gt[f(i-2),f(j)]));
    a4 = (getfield(tens,gt[f(i-1),f(j)]));
    t4a = (getfield(T[4],gt[f(i-2),f(j+1)])); 
    a2 = (getfield(tens,gt[f(i-1),f(j+1)]));
    c4c = (getfield(C[4],gt[f(i-2),f(j+2)]));
    t3d = (getfield(T[3],gt[f(i-1),f(j+2)]));

    M1 = ((c1a*t4c)*t1b)*a4
    M4 = ((c4c*t4a)*t3d)*a2
    M1 = M1/norm(((M1)))
    M4 = M4/norm(((M4)))    

    Rbar = M1*M4
    j1 = commonind(Rbar,t3d)
    j2 = commonind(Rbar,a2)
    c1 = combiner(j1; tags = "j1", dir = -dir(j1))
    c2 = combiner(j2; tags = "j2", dir = -dir(j2))
    Rbar = (Rbar*c1)*c2

    t1a = getfield(T[1],gt[f(i),f(j-1)]);
    c2b = getfield(C[2],gt[f(i+1),f(j-1)]);
    a3 = getfield(tens,gt[f(i),f(j)]);
    t2d = getfield(T[2],gt[f(i+1),f(j)]);
    a1 = getfield(tens,gt[f(i),f(j+1)]);
    t2b = getfield(T[2],gt[f(i+1),f(j+1)]);
    t3c = getfield(T[3],gt[f(i),f(j+2)]);
    c3d = getfield(C[3],gt[f(i+1),f(j+2)]);
 
    M2 = (((c2b*t1a)*t2d)*a3)
    M3 = (((c3d*t3c)*t2b)*a1)
    M2 = M2/norm(((M2)))
    M3 = M3/norm(((M3)))


    R = ((M2*M3)*dag(c1))*dag(c2)
    R = R/norm(((R)))
    Rbar = Rbar/norm(((Rbar)))

    xa = commonind(R,c1)
    cyc = commonind(R,c2)
    R = prime(R,xa)
    R = prime(R,cyc)

    Q2tilde = R*Rbar;
    Q2tilde = Q2tilde/norm(((Q2tilde)))

    u,s,v = svd(Q2tilde, (xa',cyc'), maxdim = chi)#, SVDMethod = "gesdd")#, cutoff = 1e-10);

    uindex = commonind(s,u);
    vindex = commonind(s,v);
 
    invs = invert_diag_sqrt(s);

    Ptilde = (Rbar*dag(v))*invs;
    P = (R*dag(u))*invs;
    P = P*delta(ComplexF64,vindex,uindex); # R vindex vindex uindex


    return Ptilde,P
 end
 