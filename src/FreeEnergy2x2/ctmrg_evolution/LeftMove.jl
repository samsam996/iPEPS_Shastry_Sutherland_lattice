

include("../../ctmrg_common/isometryLeft.jl")
include("../../ctmrg_common/find_value.jl")


function LeftMove!(C,T,tens,gt::Matrix{Symbol},chi::Int64,theta::Int64)
        

    P = typeof(C[1])();
    Ptilde = typeof(C[1])();
 

    N = size(gt)[1]
    f(x) = mod(x-1,N) + 1
    
    mp = [2, 3, 4, 1]; mm = [4, 1, 2, 3];
    x = [1,0]; y = [0,1];
    dr = Array{Vector}(undef,4); dr[1] = x; dr[2] = y; dr[3] = -x; dr[4] = -y; 
    del = Array{Vector}(undef,4); del[1] = -y; del[2] = x; del[3] = y; del[4] = -x; 
    
    i = 1
    j = 1

    
    pt1,p1 = isometryLeft(C,T,tens,gt,i+dr[theta][1],j+dr[theta][2],chi)
    setproperty!(P, gt[f(i+dr[theta][1]),f(j+dr[theta][2])],p1)
    setproperty!(Ptilde, gt[f(i+dr[theta][1]),f(j+dr[theta][2])],pt1)

    pt2,p2 = isometryLeft(C,T,tens,gt,i,j,chi)
    setproperty!(P, gt[f(i),f(j)],p2)
    setproperty!(Ptilde, gt[f(i),f(j)],pt2)


    for i = 1:1
        
        # C_tmp = [lattice2x1() for k = 1:4]
        # T_tmp = [lattice2x1() for k = 1:4]
        C_tmp = [typeof(C[1])() for k = 1:4]
        T_tmp = [typeof(C[1])() for k = 1:4]
    

        for kk = 1:N
            setproperty!(C_tmp[1],gt[kk,1],getproperty(C[1],gt[kk,1]))
            setproperty!(C_tmp[2],gt[kk,1],getproperty(C[2],gt[kk,1]))
            setproperty!(C_tmp[3],gt[kk,1],getproperty(C[3],gt[kk,1]))
            setproperty!(C_tmp[4],gt[kk,1],getproperty(C[4],gt[kk,1]))
            setproperty!(T_tmp[1],gt[kk,1],getproperty(T[1],gt[kk,1]))
            setproperty!(T_tmp[2],gt[kk,1],getproperty(T[2],gt[kk,1]))
            setproperty!(T_tmp[3],gt[kk,1],getproperty(T[3],gt[kk,1]))
            setproperty!(T_tmp[4],gt[kk,1],getproperty(T[4],gt[kk,1]))
        end


        for j = 1:N
        
            value1,value2,value3 = find_value(C,T,P,Ptilde,gt,tens,del,dr,theta,mp,mm,i,j)

            setproperty!(C_tmp[mp[theta]],gt[f(i+dr[theta][1]),f(j+dr[theta][2])],value1)
            setproperty!(T_tmp[theta],gt[f(i),f(j)],value2)
            setproperty!(C_tmp[theta],gt[f(i-dr[theta][1]),f(j-dr[theta][2])],value3)
           
        end

       
        for kk = 1:N
            setproperty!(C[1],gt[kk,1],getproperty(C_tmp[1],gt[kk,1]))
            setproperty!(C[2],gt[kk,1],getproperty(C_tmp[2],gt[kk,1]))
            setproperty!(C[3],gt[kk,1],getproperty(C_tmp[3],gt[kk,1]))
            setproperty!(C[4],gt[kk,1],getproperty(C_tmp[4],gt[kk,1]))
            setproperty!(T[1],gt[kk,1],getproperty(T_tmp[1],gt[kk,1]))
            setproperty!(T[2],gt[kk,1],getproperty(T_tmp[2],gt[kk,1]))
            setproperty!(T[3],gt[kk,1],getproperty(T_tmp[3],gt[kk,1]))
            setproperty!(T[4],gt[kk,1],getproperty(T_tmp[4],gt[kk,1]))
        end




    end

    nothing

end
