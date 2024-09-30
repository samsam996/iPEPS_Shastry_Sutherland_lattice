


function isometryDown(C,T,tens,gt::Matrix{Symbol},i::Int64,j::Int64,chi::Int64)
    
  N = size(gt)[1];
  f(x) = mod(x-1,N)+1;
 
  c1d = getfield(C[1],gt[f(i-1),f(j-2)]);
  t1c = getfield(T[1],gt[f(i),f(j-2)]);
  t4b = getfield(T[4],gt[f(i-1),f(j-1)]);
  a1 = getfield(tens,gt[f(i),f(j-1)]);
  t4d = getfield(T[4],gt[f(i-1),f(j)]);
  a3 = getfield(tens,gt[f(i),f(j)]);
  c4b = getfield(C[4],gt[f(i-1),f(j+1)]);
  t3a = getfield(T[3],gt[f(i),f(j+1)]);
     
  M1 = ((c1d*t1c)*t4b)*a1
  M4 = ((c4b*t3a)*t4d)*a3

  M1 = M1/norm(((M1)))
  M4 = M4/norm(((M4)))

  R = M1*M4
  j1 = commonind(R,t1c)
  j2 = commonind(R,a1)
  c1 = combiner(j1; tags = "j1", dir = -dir(j1))
  c2 = combiner(j2; tags = "j2", dir = -dir(j2))
  R = (R*c1)*c2

  t1d = getfield(T[1],gt[f(i+1),f(j-2)]);
  c2c = getfield(C[2],gt[f(i+2),f(j-2)]);
  a2 = getfield(tens,gt[f(i+1),f(j-1)]);
  t2a = getfield(T[2],gt[f(i+2),f(j-1)]);
  a4 = getfield(tens,gt[f(i+1),f(j)]);
  t2c = getfield(T[2],gt[f(i+2),f(j)]);
  t3b = getfield(T[3],gt[f(i+1),f(j+1)]);
  c3a = getfield(C[3],gt[f(i+2),f(j+1)]);
  
  M2 = ((c2c*t1d)*t2a)*a2
  M3 = ((c3a*t3b)*t2c)*a4

  M2 = M2/norm(((M2)))
  M3 = M3/norm(((M3)))

  Rbar = ((M2*M3)*dag(c1))*dag(c2)
  cxd = commonind(Rbar,dag(c1))
  xb = commonind(Rbar,dag(c2))
  Rbar = prime(Rbar,cxd)
  Rbar = prime(Rbar,xb)

  R = R/norm(((R)))
  Rbar = Rbar/norm(((Rbar)))

  Q2tilde = R*Rbar;
  Q2tilde = Q2tilde/norm(((Q2tilde)))

  u,s,v = svd(Q2tilde, (cxd,xb), maxdim = chi)#, SVDMethod = "gesvd")#,cutoff = 1e-10);

  uindex = commonind(s,u);
  vindex = commonind(s,v);
 
  invs = invert_diag_sqrt(s);
 
  Ptilde = ((Rbar*dag(v))*invs); 
  P = (R*dag(u))*invs;
  P = *(P,delta(vindex,uindex)); 
 
  return Ptilde,P

 end