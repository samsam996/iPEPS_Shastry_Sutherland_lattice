

include("update.jl")
include("../../ctmrg_common/local_obs.jl")




function ctm!(tens_a,tens_A,cxd,cyd,gt::Matrix{Symbol},physical_legs,ancilla_legs,chi::Int,precision::Float64,C,T)  

  

  sx = 1/2*[0 1; 1 0]; sy = 1/2*[0 -1im; 1im 0]; sz = 1/2*[1 0; 0 -1]; id = [1 0; 0 1];
  kk = 0;
  
  tmp = Vector([2]); err = 1; iter = 0; 
  errtmp = 1

  while err > precision 

    iter = iter + 1; 

    @showtime update!(C,T,tens_a,gt,chi)
  
    tmp = kk;

    obs1::Matrix{ComplexF64} = zeros(4,4)
    obs1 = kron(sx,sx) + kron(sy,sy) + kron(sz,sz)
        
    kk = local_obs(C,T,tens_a,tens_A,cxd,cyd,gt,physical_legs,ancilla_legs,1,1,obs1) + 
    local_obs(C,T,tens_a,tens_A,cxd,cyd,gt,physical_legs,ancilla_legs,1,2,obs1)

    errtmp = norm(kk-tmp)
    err = norm(kk-tmp)
    println("err : ", err, "// it : ", iter)
                
    if iter > 300
      err = 0
    end

  end

  return C, T, iter, errtmp

end
