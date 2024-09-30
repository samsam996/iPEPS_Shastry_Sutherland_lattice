
include("LeftMove.jl")
include("RightMove.jl")
include("UpMove.jl")
include("DownMove.jl")
include("../../ctmrg_common/renormalisation.jl")


function update!(C,T,tens,gt::Matrix{Symbol},chi::Int64)

  renormalisation!(C,T,gt)  

  @showtime LeftMove!(C,T,tens,gt,chi,4)
  renormalisation!(C,T,gt)  

  @showtime RightMove!(C,T,tens,gt,chi,2)
  renormalisation!(C,T,gt)    

  @showtime UpMove!(C,T,tens,gt,chi,1)
  renormalisation!(C,T,gt)  
  
  @showtime DownMove!(C,T,tens,gt,chi,3)
  renormalisation!(C,T,gt)  
 
  nothing

end