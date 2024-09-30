

function correlation_length(C,T,gt)

    N = size(gt)[1]
    f(x) = mod(x-1,N) + 1

    xi2 = zeros(N,N)
    xi3 = zeros(N,N)
    xi4 = zeros(N,N)
    dq = zeros(N,N)


 
    for i = 1:N
        for j = 1:1
        
            #  - T1a - T1b - x
            t1a = getproperty(T[1],gt[f(i),f(j)])
            t1b = getproperty(T[1],gt[f(i+1),f(j)])
            t1c = getproperty(T[1],gt[f(i+2),f(j)])
            t1d = getproperty(T[1],gt[f(i+3),f(j)])
            t1e = getproperty(T[1],gt[f(i+4),f(j)])
            t1f = getproperty(T[1],gt[f(i+5),f(j)])
            
            t3a = getproperty(T[3],gt[f(i),f(j+1)])
            t3b = getproperty(T[3],gt[f(i+1),f(j+1)])
            t3c = getproperty(T[3],gt[f(i+2),f(j+1)])
            t3d = getproperty(T[3],gt[f(i+3),f(j+1)])
            t3e = getproperty(T[3],gt[f(i+4),f(j+1)])
            t3f = getproperty(T[3],gt[f(i+5),f(j+1)])
            
            c2 = getproperty(C[2],gt[f(i+6),f(j)])
            c3 = getproperty(C[3],gt[f(i+6),f(j+1)])
            
            indexb = commonind(c2,t1f)
            indexd = commonind(c3,t3f)
            
            t1a = prime(t1a,indexb)
            t3c = prime(t3c,indexd)
            
            x0 = randomITensor(indexb,indexd)
            
            # add the other 4 t's
            g(x) = noprime((((((((((((t1f*x)*t3f)*t1e)*t3e)*t1d)*t3d)*t1c)*t3c)*t1b)*t3b)*t1a)*t3a)
            val,vec,info = eigsolve(g, x0, 9, :LM)
                        
            if length(val) > 2
                xi3[i,j] = -log(abs(val[3])/abs(val[1]));
            end
            if length(val) > 3
                xi4[i,j] = -log(abs(val[4])/abs(val[1]));
            end
            
            if length(val) > 1
            xi2[i,j] = -log(abs(val[2])/abs(val[1]));
                dq[i,j] = abs(angle(val[2]))
            end

        end
    end

    return xi2, xi3, xi4, dq

end