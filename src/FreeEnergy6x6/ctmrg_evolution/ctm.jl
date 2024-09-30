


include("update.jl")
include("../../ctmrg_common/local_obs.jl")



function ctm!(tens_a,tens_A,cxd,cyd,gt::Matrix{Symbol},physical_legs,ancilla_legs,chi::Int,precision::Float64,C,T)  
  
  sx = 1/2*[0 1; 1 0]; sy = 1/2*[0 -1im; 1im 0]; sz = 1/2*[1 0; 0 -1]; id = [1 0; 0 1];
  kk = 0;
  
  err = 1;
  tmp = Vector([2]); err = 1; iter = 0; 
  errtmp = 0;
  
  while err > precision 

    iter = iter + 1; 

    err_prev = err; 
    @showtime update!(C,T,tens_a,gt,chi)
  
    tmp = kk;

    obs1::Matrix{ComplexF64} = (kron(sz,sz) + kron(sx,sx) + kron(sy,sy)) 
    kk = local_obs(C,T,tens_a,tens_A,cxd,cyd,gt,physical_legs,ancilla_legs,1,1,obs1)

    err = norm(kk-tmp)
    println("err : ", err, "// it : ", iter)
    errtmp = deepcopy(err)

    if norm(err_prev - errtmp) < 1e-10
      err = 0
    end

    if iter > 200
      err = 0
    end

  end

  return C, T, iter, errtmp

end
