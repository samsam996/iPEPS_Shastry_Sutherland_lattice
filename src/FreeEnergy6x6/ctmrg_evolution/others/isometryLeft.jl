

function isometryLeft(C,T,tens,gt::Matrix{Symbol},
    i::Int64,j::Int64,chi::Int64)

    N = size(gt)[1];
    f(x) = mod(x-1,N)+1;
     
    c1d = getfield(C[1],gt[f(i-1),f(j-1)]);
    t1c = getfield(T[1],gt[f(i),f(j-1)]);
    t1d = getfield(T[1],gt[f(i+1),f(j-1)]);
    c2c = getfield(C[2],gt[f(i+2),f(j-1)]);
    t4b = getfield(T[4],gt[f(i-1),f(j)]);
    a1 = getfield(tens,gt[f(i),f(j)]);
    a2 = getfield(tens,gt[f(i+1),f(j)]);
    t2a = getfield(T[2],gt[f(i+2),f(j)]); 
 
    M1 = ((c1d*t1c)*t4b)*a1
    M2 = ((c2c*t2a)*t1d)*a2
    M1 = M1/norm(M1)
    M2 = M2/norm(M2)


    R = M1*M2
    j1 = commonind(R,t2a)
    j2 = commonind(R,a2)
    c1 = combiner(j1; tags = "j1", dir = -dir(j1))
    c2 = combiner(j2; tags = "j2", dir = -dir(j2))
    R = (R*c1)*c2

    # R : cyb ya yb' cya'
    yb = commonind(R,c1)
    cya = commonind(R,c2)
    R = prime(R,yb)
    R = prime(R,cya) # R : cyb ya yb' cya'


    t4d = getfield(T[4],gt[f(i-1),f(j+1)]); # ok
    a3 = getfield(tens,gt[f(i),f(j+1)]);
    a4 = getfield(tens,gt[f(i+1),f(j+1)]);
    t2c = getfield(T[2],gt[f(i+2),f(j+1)]);
    
    c4b = getfield(C[4],gt[f(i-1),f(j+2)]);
    t3a = getfield(T[3],gt[f(i),f(j+2)]);
    t3b = getfield(T[3],gt[f(i+1),f(j+2)]);
    c3a = getfield(C[3],gt[f(i+2),f(j+2)]);
         
    M4 = ((c4b*t4d)*t3a)*a3
    M3 = ((c3a*t3b)*t2c)*a4
    M4 = M4/norm(M4)
    M3 = M3/norm(M3)

    Rbar = ((M4*M3)*dag(c1))*dag(c2)
    Rbar = Rbar/norm(Rbar)
    R = R/norm(R)

    Q2tilde = R*Rbar; #  yb cya yb' cya'
    Q2tilde = Q2tilde/norm(Q2tilde)

    u,s,v = svd(Q2tilde,(cya',yb'), maxdim = chi)#, SVDMethod = "gesvd")#, cutoff = 1e-10)

    uindex = commonind(u,s);
    vindex = commonind(v,s);
    
    invs = invert_diag_sqrt(s);
    invs = invs/norm(invs)

    
    Ptilde = (Rbar*dag(v))*invs; 
    P = (R*dag(u))*invs; # R .. uindex uindex vindex 
    P = P*delta(dag(vindex),dag(uindex)); # R vindex vindex uindex
 

    return Ptilde, P 
 
 end
 
 