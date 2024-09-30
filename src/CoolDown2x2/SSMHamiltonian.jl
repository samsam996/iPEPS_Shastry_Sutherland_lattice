

function SSMHamiltonian(gt::Matrix{Symbol},physical_legs::lattice_ind2x1,nsu::Float64,J1::Float64,J2::Float64,h::Float64,hs::Float64)

    gx = lattice2x1();
    gy = lattice2x1();
    mu = lattice2x1();
    mu_hs = lattice2x1();

    N = size(gt)[1];
    f(x) = mod(x-1,N) + 1;

    sx = 1/2*[0 1; 1 0]; sy = 1/2*[0 -1im; 1im 0]; sz = 1/2*[1 0; 0 -1]; id = [1 0; 0 1];
    
    h2 = J2*(kron(sz,sz) + kron(sx,sx) + kron(sy,sy)) - h*kron(sz,id) - h*kron(id,sz);
    e2 = exp(-h2/nsu);
    gin = reshape(e2 ,(4, 4)); 

    h2_minus = - hs*kron(sz,id) - hs*kron(id,sz);
    e2 = exp(-h2_minus/nsu);
    ghs_minus = reshape(e2 ,(4, 4)); 

    h2_plus = + hs*kron(sz,id) + hs*kron(id,sz);
    e2 = exp(-h2_plus/nsu);
    ghs_plus = reshape(e2 ,(4, 4)); 

    h1 = J1*(kron(kron(sx,id),kron(sx,id)) + kron(kron(sy,id),kron(sy,id)) + kron(kron(sz,id),kron(sz,id))) + 
        J1*(kron(kron(id,sx),kron(sx,id)) + kron(kron(id,sy),kron(sy,id)) + kron(kron(id,sz),kron(sz,id))) 
    e1 = exp(-h1/nsu);
    gr = reshape(e1, (4,4,4,4));  # 21 43 21 43

    h2 = J1*(kron(kron(id,sx),kron(sx,id)) + kron(kron(id,sy),kron(sy,id)) + kron(kron(id,sz),kron(sz,id))) + 
        J1*(kron(kron(id,sx),kron(id,sx)) + kron(kron(id,sy),kron(id,sy)) + kron(kron(id,sz),kron(id,sz))) 
    e2 = exp(-h2/nsu);
    gl = reshape(e2, (4,4,4,4)); # 21 43 21 43
  
    h3 = J1*(kron(kron(sx,id),kron(id,sx)) + kron(kron(sy,id),kron(id,sy)) + kron(kron(sz,id),kron(id,sz))) + 
        J1*(kron(kron(id,sx),kron(id,sx)) + kron(kron(id,sy),kron(id,sy)) + kron(kron(id,sz),kron(id,sz))) 
    e3 = exp(-h3/nsu);
    gu = reshape(e3, (4,4,4,4)); # 21 43 21 43

    h4 = J1*(kron(kron(sx,id),kron(sx,id)) + kron(kron(sy,id),kron(sy,id)) + kron(kron(sz,id),kron(sz,id))) + 
        J1*(kron(kron(sx,id),kron(id,sx)) + kron(kron(sy,id),kron(id,sy)) + kron(kron(sz,id),kron(id,sz))) 
    e4 = exp(-h4/nsu);
    gd = reshape(e4, (4,4,4,4)); # 21 43 21 43

    for i = 1:1:N
        for j = 1:1:1

            ia = getfield(physical_legs,gt[f(i),f(j)]);
            ib = getfield(physical_legs,gt[f(i+1),f(j)]);
            ic = getfield(physical_legs,gt[f(i),f(j+1)]);
    
            tmpgx = ITensor(dag(ia),dag(ib),ia',ib');
            tmpgy = ITensor(dag(ic),dag(ia),ic',ia');
            tmpmu = ITensor(dag(ia),ia');
            tmphs = ITensor(dag(ia),ia')

            cin = [-2,0,0,2]
            cout = [-2,0,0,2]

            for i1 = 1:4
                for i2 = 1:4
                    for j1 = 1:4
                        for j2 = 1:4
                            if abs(mod(i+j,2) - 0) < 1e-9 # pair
                                if abs(cin[i1]+cin[i2]-cout[j1]-cout[j2]) < 1e-9                         
                                tmpgx[dag(ia)=>i1,dag(ib)=>i2,ia'=>j1,ib'=>j2] = gr[i1,i2,j1,j2];
                                end
                            end
                            if abs(mod(i+j,2) - 1) < 1e-9 # impair
                                if abs(cin[i1]+cin[i2]-cout[j1]-cout[j2]) <1e-9      
                                tmpgx[dag(ia)=>i1,dag(ib)=>i2,ia'=>j1,ib'=>j2] = gl[i1,i2,j1,j2];
                                end
                            end
                        end
                    end
                end
            end

            setproperty!(gx,gt[f(i),f(j)],tmpgx);

            for i1 = 1:4
                for i2 = 1:4
                    for j1 = 1:4
                        for j2 = 1:4
                            if abs(mod(i+j,2) - 0) < 1e-9 # pair
                                if abs(cin[i1]+cin[i2]-cout[j1]-cout[j2]) < 1e-9    
                                tmpgy[dag(ic)=>i1,dag(ia)=>i2,ic'=>j1,ia'=>j2] = gu[i1,i2,j1,j2];
                                end
                            end
                            if abs(mod(i+j,2)-1) < 1e-9 # impair
                                if abs(cin[i1]+cin[i2]-cout[j1]-cout[j2]) < 1e-9    
                                tmpgy[dag(ic)=>i1,dag(ia)=>i2,ic'=>j1,ia'=>j2] = gd[i1,i2,j1,j2];
                                end
                            end
                        end
                    end
                end
            end

            setproperty!(gy,gt[f(i),f(j)],tmpgy);

            for i1 = 1:4
                for j1 = 1:4
                    if abs(cin[i1]-cout[j1]) < 1e-9    
                    tmpmu[dag(ia)=>i1,ia'=>j1] = gin[i1,j1];
                    end
                end
            end

            setproperty!(mu,gt[f(i),f(j)],tmpmu);

            for i1 = 1:4
                for i2 = 1:4
                    if i == 1 || i == 4
                        if abs(cin[i1]-cout[i2]) < 1e-9
                            tmphs[dag(ia)=>i1,ia'=>i2] = ghs_minus[i1,i2]
                        end
                    elseif i == 2 || i == 5
                        if abs(cin[i1]-cout[i2]) < 1e-9
                            tmphs[dag(ia)=>i1,ia'=>i2] = ghs_plus[i1,i2]
                        end
                    end
                end
            end

            setproperty!(mu_hs,gt[f(i),f(j)],tmphs);

        end
    end

    return gx, gy, mu, mu_hs

end

