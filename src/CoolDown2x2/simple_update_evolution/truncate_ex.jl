
function truncate_ex(A::ITensor,B::ITensor,i::Int,j::Int,gx::lattice2x1,gt::Matrix{Symbol},
  physical_legs::lattice_ind2x1,ind_new::Index,D::Int64)

  N = size(gt)[1];
  f(x) = mod(x-1,N) + 1;

  exphx = getfield(gx,gt[f(i),f(j)]);

  # A : xb1 ya2 xa yc1 ia sa
  # B : xa yb2 xb2 yd1 ib sb

  ia = getproperty(physical_legs,gt[f(i),f(j)]);
  ib = getproperty(physical_legs,gt[f(i+1),f(j)]);

  indsA = inds(A); # i1 i2 i3 i4 i5 i6
  qa = noncommoninds(indsA, [ind_new,ia]);

  ua,sa,va = svd(A,qa);
  ub,sb,vb = svd(B,[ind_new,ib]);

  uaindex = commonind(ua,sa);
  vbindex = commonind(sb,vb);
  va = sa*va;
  ub = ub*sb;

  M = (va*ub)*exphx;
  M = noprime(M)
  
  um,sm,vm = svd(M, [uaindex,ia], maxdim = D, cutoff = 1e-12);
  fx = log(norm(sm));
  sm = sm/norm(sm);

  umindex = commonind(sm,um)
  vmindex = commonind(sm,vm)

  AA = ua*um;
  BB = vm*vb;

  return AA,BB,sm, fx  

end