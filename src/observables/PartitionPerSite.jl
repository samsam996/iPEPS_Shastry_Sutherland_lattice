

function PartitionPerSite(T,C,tens,gt::Matrix{Symbol},i::Int64,j::Int64)


    N = size(gt)[1]
    f(x) = mod(x-1,N)+1;


    t1a = getfield(T[1], gt[i,j])
    t2a = getfield(T[2], gt[i,j])
    t3a = getfield(T[3], gt[i,j])
    t4a = getfield(T[4], gt[i,j])

    c3a = getfield(C[3],gt[f(i),f(j)])
    c2b = getfield(C[2],gt[f(i+1),f(j)])
    c4f = getfield(C[4],gt[f(i+5),f(j)])

    c1a = getfield(C[1],gt[f(i-1),f(j-1)])
    t1b = getfield(T[1],gt[f(i),f(j-1)])
    t1c = getfield(T[1],gt[f(i+1),f(j-1)])
    c2d = getfield(C[2],gt[f(i+2),f(j-1)])
    
    t4f = getfield(T[4],gt[f(i-1),f(j)])
    a1 = getfield(tens,gt[f(i),f(j)])
    a2 = getfield(tens,gt[f(i+1),f(j)])
    t2c = getfield(T[2],gt[f(i+2),f(j)])

    t4g = getfield(T[4],gt[f(i-1),f(j+1)])
    a3 = getfield(tens,gt[f(i),f(j+1)])
    a4 = getfield(tens,gt[f(i+1),f(j+1)])
    t2b = getfield(T[2],gt[f(i+2),f(j+1)])

    c4d = getfield(C[4],gt[f(i-1),f(j+2)])
    t3g = getfield(T[3],gt[f(i),f(j+2)])
    t3f = getfield(T[3],gt[f(i+1),f(j+2)])
    c3a = getfield(C[3],gt[f(i+2),f(j+2)])


    # a -- b -- c -- d
    # f -- a -- b -- c
    # g -- f -- a -- b
    # d -- g -- f -- a

    carree = (((c1a*c2b)*c3a)*c4f)[]
  
    # carree = carree*delta(dag(commonind(c2b,t2a)),commonind(t2b,c3a))*c3a
    # carree = carree*delta(dag(commonind(c4d,t3a)),(commonind(t3d,c3a)))*c4d
    # carree = carree*delta(dag(commonind(c1a,t4d)), commonind(t4a,c4d))



    C1 = *(c1a,t1b,t4f,a1)
    C2 = *(c2d,t1c,t2c,a2)
    C3 = *(c3a,t2b,t3f,a4)
    C4 = *(c4d,t4g,t3g,a3)
    Z = ((C1*C2)*(C3*C4))[]#*(c1d,t1c,t4b,a1,c2c,t2a,t1d,a2,c4b,t3a,t4d,a3,c3a,t2c,t3b,a4)

    # carre = (((c1d*c2c)*c3a)*c4b)[]
   


    # a -- b -- c -- d
    # f -- a -- b -- c
    # g -- f -- a -- b
    # d -- g -- f -- a

    c4f = getproperty(C[4],gt[f(i+5),j])
    c3c = getproperty(C[3],gt[f(i+2),j])
    t3b = getproperty(T[3],gt[f(i+1),j])
    t3a = getproperty(T[3],gt[i,j])

    c3g = getproperty(C[3],gt[f(i+4),f(j)])
    t2a = getproperty(T[2],gt[f(i),f(j)])
    t2f = getproperty(T[2],gt[f(i+5),f(j)])

    vert = (((((((c1a*c4f)*t1b)*t3a)*t1c)*t3b)*c2d)*c3c)[];
    
    hor = (((((((c1a*c2b)*t4f)*t2a)*t4g)*t2f)*c4d)*c3g)[];

    @show xx = (abs(Z*carree)/abs(vert*hor))

    return log(xx)

    
end