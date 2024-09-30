


function find_value(C,T,P,Ptilde,gt::Matrix{Symbol},tens,
    del::Array{Vector},dr::Array{Vector},theta::Int64,mp::Vector{Int64},mm::Vector{Int64},i::Int64,j::Int64)
    
    N = size(gt)[1]
    f(x) = (mod(x-1,N)+1); 

    value1 = (getfield(C[mp[theta]],gt[f(i+dr[theta][1]+del[theta][1]),f(j+dr[theta][2]+del[theta][2])])*
    getfield(T[mp[theta]],gt[f(i+dr[theta][1]),f(j+dr[theta][2])]))*
    getfield(Ptilde,gt[f(i+dr[theta][1]),f(j+dr[theta][2])]);

    value2 = ((
    getfield(P,gt[f(i+dr[theta][1]),f(j+dr[theta][2])])*
    getfield(T[theta],gt[f(i+del[theta][1]),f(j+del[theta][2])]))*
    getfield(tens,gt[i,j]))*
    getfield(Ptilde,gt[i,j])

    value3 =  (
    getfield(C[theta],gt[f(i-dr[theta][1]+del[theta][1]),f(j-dr[theta][2]+del[theta][2])])*
    getfield(T[mm[theta]],gt[f(i-dr[theta][1]),f(j-dr[theta][2])]))*
    getfield(P,gt[f(i),f(j)])


    
    value1 = value1/norm((value1))
    value2 = value2/norm((value2))
    value3 = value3/norm((value3))
    

    return value1, value2, value3

end
