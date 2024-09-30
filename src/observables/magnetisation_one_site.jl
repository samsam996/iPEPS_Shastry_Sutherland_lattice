
"""
    ind_x and ind_y are indices of tens_A 
    c_x_a and c_y_a are what maps the double indices into the indices of tens_a

"""
function magnetisation_one_site(C::Vector{lattice},T::Vecto{lattice},tens_a_raw,tens_a,tens_A,gt,ind_x,ind_y,indx_a,indy_a,c_x_a,c_y_a,entanglement_legs,i,j)

    N = size(gt)[1]
    f(x) = mod(x-1,N) + 1

    magne = lattice()

    Adown = getproperty(tens_A, gt[f(i),f(j)])
    ind1 = getproperty(ind_x, gt[f(i),f(j)])
    ind2 = getproperty(ind_y, gt[f(i),f(j)])
    ind3 = getproperty(ind_x, gt[f(i-1),f(j)])
    ind4 = getproperty(ind_y, gt[f(i),f(j-1)])

    sa = getproperty(entanglement_legs, gt[f(i),f(j)])

    all_ind = inds(Adown)
    ia = noncommonind(all_ind,[ind1,ind2,ind3,ind4,sa])

    Aup = prime(Adown, ind1, ind2, ind3, ind4, ia)

    sz = 1/2*[1 0; 0 -1]; id = [1 0; 0 1]
    c1 = [-2 0 0 2]; c2 = [-2 0 0 2];

    obs = kron(sz,id) + kron(id,sz)
    sz = ITensor(dag(ia),ia')

    for i1 = 1:4
        for i2 = 1:4
            if abs(c1[i1] - c2[i2]) < 1e-9
            sz[ia=>i1,ia'=>i2] = obs[i1,i2]
            end
        end
    end

    AA = Aup*dag(sz)*dag(Adown)

    comb1 = getproperty(c_x_a, gt[f(i-1),f(j)])
    comb2 = getproperty(c_y_a, gt[f(i),f(j)])
    comb3 = getproperty(c_x_a, gt[f(i),f(j)])
    comb4 = getproperty(c_y_a, gt[f(i),f(j-1)])

    tens_a_Raw = getproperty(tens_a_raw,gt[f(i),f(j)]) 

    ind1 = commonind(comb1, tens_a_Raw)
    ind2 = commonind(comb2, tens_a_Raw)
    ind3 = commonind(comb3, tens_a_Raw)
    ind4 = commonind(comb4, tens_a_Raw)
   
    ind1_new = getproperty(indx_a, gt[f(i-1),f(j)])
    ind2_new = getproperty(indy_a, gt[f(i),f(j)])
    ind3_new = getproperty(indx_a, gt[f(i),f(j)])
    ind4_new = getproperty(indy_a, gt[f(i),f(j-1)])

    AAA = AA*comb1*comb2*comb3*comb4
    AAA = AAA*delta(ComplexF64,dag(ind1),dag(ind1_new))*delta(ComplexF64,dag(ind2),ind2_new)*delta(ComplexF64,dag(ind3),ind3_new)*delta(ComplexF64,dag(ind4),dag(ind4_new))

    setproperty!(magne,gt[f(i),f(j)],AAA)

    c1 = getproperty(C[1],gt[f(i-1),f(j-1)])
    t1 = getproperty(T[1],gt[f(i),f(j-1)])
    c2 = getproperty(C[2],gt[f(i+1),f(j-1)])
    t4 = getproperty(T[4],gt[f(i-1),f(j)])
    aa = getproperty(tens_a,gt[f(i),f(j)])
    local_mm = getproperty(magne,gt[f(i),f(j)])
    t2 = getproperty(T[2],gt[f(i+1),f(j)])
    c4 = getproperty(C[4],gt[f(i-1),f(j+1)])
    t3 = getproperty(T[3],gt[f(i),f(j+1)])
    c3 = getproperty(C[3],gt[f(i+1),f(j+1)])

    mm = ((((((((c1*t1)*t4)*local_mm)*c2)*t2)*c3)*t3)*c4)[1]/((((((((((c1*t1)*t4)*aa)*c2)*t2)*c3)*t3)*c4)))[1]
   
    return mm

end