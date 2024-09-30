



function declare_indices(gt::Matrix{Symbol})

    N = size(gt)[1]
    f(x) = mod(x-1,N) + 1
    
    physical_legs = lattice_ind2x1("physical")
    ancilla_legs = lattice_ind2x1("ancilla")

    x1 = lattice_ind2x1("imaginary")
    x2 = lattice_ind2x1("imaginary")
    y1 = lattice_ind2x1("imaginary")
    y2 = lattice_ind2x1("imaginary")

    lambdax = lattice2x1()
    lambday = lattice2x1()
    Gamma = lattice2x1()

    for i = 1:N
        for j = 1:1

            indx1 = getproperty(x1,gt[i,j])
            indx2 = getproperty(x2,gt[i,j])
            indy1 = getproperty(y1,gt[i,j])
            indy2 = getproperty(y2,gt[i,j])

            lambdaxij = ITensor(dag(indx1),indx2)
            lambdayij = ITensor(dag(indy1),indy2)

            setproperty!(lambdax,gt[i,j],lambdaxij)
            setproperty!(lambday,gt[i,j],lambdayij)

        end
    end

    for i = 1:N
        for j = 1:1

            ind1 = getproperty(x2,gt[f(i-1),f(j)])
            ind2 = getproperty(y1,gt[f(i),f(j)])
            ind3 = getproperty(x1,gt[f(i),f(j)])
            ind4 = getproperty(y2,gt[f(i),f(j-1)])
            sa = getproperty(ancilla_legs,gt[f(i),f(j)])
            ia = getproperty(physical_legs,gt[f(i),f(j)])

            AA = ITensor(dag(ind1),ind2,ind3,dag(ind4),ia,dag(sa))
            setproperty!(Gamma,gt[i,j],AA)

        end
    end

    return Gamma, lambdax, lambday, physical_legs, ancilla_legs

end

function initialisation()

    tt = lattice2x1()
    gg = tt.gg;
    gt = tt.gt;

    N = size(gt)[1]
    f(x) = mod(x-1,N) + 1

    Gamma, lambdax, lambday, physical_legs, ancilla_legs = declare_indices(gt)

    for i = 1:N
        for j = 1:1
            
            tmp_gamma = ITensor(inds(getfield(Gamma,gt[f(i),f(j)])));
            leg1 = commonind(tmp_gamma,getfield(lambdax,gt[f(i-1),f(j)]))
            leg2 = commonind(tmp_gamma,getfield(lambday,gt[f(i),f(j)]))
            leg3 = commonind(tmp_gamma,getfield(lambdax,gt[f(i),f(j)]))
            leg4 = commonind(tmp_gamma,getfield(lambday,gt[f(i),f(j-1)]))
            leg5 = getproperty(physical_legs,gt[i,j])
            leg6 = getproperty(ancilla_legs,gt[i,j])

            tmp_gamma[leg1=>1,leg2=>1,leg3=>1,leg4=>1,leg5=>1,dag(leg6)=>1] = 1
            tmp_gamma[leg1=>1,leg2=>1,leg3=>1,leg4=>1,leg5=>2,dag(leg6)=>2] = 1
            tmp_gamma[leg1=>1,leg2=>1,leg3=>1,leg4=>1,leg5=>3,dag(leg6)=>3] = 1
            tmp_gamma[leg1=>1,leg2=>1,leg3=>1,leg4=>1,leg5=>4,dag(leg6)=>4] = 1
            
            setproperty!(Gamma,gt[f(i),f(j)],tmp_gamma);

            tmp_lx = ITensor(inds(getfield(lambdax,gt[f(i),f(j)])));
            tmp_ly = ITensor(inds(getfield(lambday,gt[f(i),f(j)])));

            xa11 = commonind(getfield(Gamma,gt[f(i),f(j)]),getfield(lambdax,gt[f(i),f(j)]))
            xa22 = commonind(getfield(Gamma,gt[f(i+1),f(j)]),getfield(lambdax,gt[f(i),f(j)]))

            ya11 = commonind(getfield(Gamma,gt[f(i),f(j)]),getfield(lambday,gt[f(i),f(j)]))
            ya22 = commonind(getfield(Gamma,gt[f(i),f(j+1)]),getfield(lambday,gt[f(i),f(j)]))

            for pp = 1:1                
                tmp_lx[xa11=>pp,xa22=>pp] = 1; 
                tmp_ly[ya11=>pp,ya22=>pp] = 1;
            end

            setproperty!(lambdax,gt[f(i),f(j)],tmp_lx);
            setproperty!(lambday,gt[f(i),f(j)],tmp_ly);

        end
    end

    return Gamma, lambdax, lambday, physical_legs, ancilla_legs, gt, gg
    
end