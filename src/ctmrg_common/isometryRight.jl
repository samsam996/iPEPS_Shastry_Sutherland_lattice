

include("invert_diag_sqrt.jl")

function isometryRight(C,T,tens,gt::Matrix{Symbol},i::Int64,j::Int64,chi::Int64)

   N = size(gt)[1];
   f(x) = mod(x-1,N)+1;

   c1a = getfield(C[1],gt[f(i-2),f(j-2)]);
   t1b = getfield(T[1],gt[f(i-1),f(j-2)]);
   t1a = getfield(T[1],gt[f(i),f(j-2)]);
   c2b = getfield(C[2],gt[f(i+1),f(j-2)]);
   t4c = getfield(T[4],gt[f(i-2),f(j-1)]);
   a1 = getfield(tens,gt[f(i-1),f(j-1)]);
   a2 = getfield(tens,gt[f(i),f(j-1)]);
   t2d = getfield(T[2],gt[f(i+1),f(j-1)]);
   
   M1 = ((c1a*t1b)*t4c)*a1
   M2 = ((c2b*t1a)*t2d)*a2

   M1 = M1/norm(((M1)))
   M2 = M2/norm(((M2)))

   Rbar = M1*M2
   jj1 = commoninds(Rbar,t4c)
   jj2 = commoninds(Rbar,a1)
   j1 = jj1[1]
   j2 = jj2[1]
   c1 = combiner(j1; tags = "j1", dir = -dir(j1))
   c2 = combiner(j2; tags = "j2", dir = -dir(j2))
   Rbar = Rbar*c1*c2

   yd = commonind(c1,Rbar)
   cyc = commonind(c2,Rbar)
   Rbar = prime(Rbar,yd)
   Rbar = prime(Rbar,cyc)

   t4a = getfield(T[4],gt[f(i-2),f(j)]);
   a3 = getfield(tens,gt[f(i-1),f(j)]);
   a4 = getfield(tens,gt[f(i),f(j)]);
   t2b = getfield(T[2],gt[f(i+1),f(j)]);
   c4c = getfield(C[4],gt[f(i-2),f(j+1)]);
   t3d = getfield(T[3],gt[f(i-1),f(j+1)]);
   t3c = getfield(T[3],gt[f(i),f(j+1)]);
   c3d = getfield(C[3],gt[f(i+1),f(j+1)]);

   M3 = ((c3d*t3c)*t2b)*a4
   M4 = ((c4c*t4a)*t3d)*a3

   M3 = M3/norm(((M3)))
   M4 = M4/norm(((M4)))

   R = ((M3*M4)*dag(c1))*dag(c2)

   Rbar = Rbar/norm(((Rbar)))
   R = R/norm(((R)))

   Q2tilde = (Rbar*R);     
   Q2tilde = Q2tilde/norm(((Q2tilde)))

   u,s,v = svd(Q2tilde,(cyc,yd), maxdim = chi)#, SVDMethod = "gesvd", cutoff = 1e-10);
   # u,s,v = svd(Q2tilde,(cyc,yd),maxdim = chi)#,cutoff = 1e-10);

   uindex = commonind(s,u);
   vindex = commonind(s,v);

   invs = invert_diag_sqrt(s);

   Ptilde = (Rbar*dag(v))*invs; 
   P = (R*dag(u))*invs; 
   P = P*(delta(vindex,uindex)); 

   return Ptilde, P

end