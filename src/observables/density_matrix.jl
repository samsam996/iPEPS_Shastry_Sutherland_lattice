

function density_matrix(C::Vector{unit_cell.lattice},T::Vector{unit_cell.lattice},
    tens_a_raw::unit_cell.lattice,tens_a::unit_cell.lattice,tens_A::unit_cell.lattice,
    gt,ind_x,ind_y,indx_a,indy_a,c_x_a,c_y_a,entanglement_legs,i,j)

    N = size(gt)[1]
    f(x) = mod(x-1,N) + 1

    magne = twoxtwo.lattice()

    Adown = getproperty(tens_A, gt[f(i),f(j)])
    ind1 = getproperty(ind_x, gt[f(i),f(j)])
    ind2 = getproperty(ind_y, gt[f(i),f(j)])
    ind3 = getproperty(ind_x, gt[f(i-1),f(j)])
    ind4 = getproperty(ind_y, gt[f(i),f(j-1)])

    sa = getproperty(entanglement_legs, gt[f(i),f(j)])

    all_ind = inds(Adown)
    ia = noncommonind(all_ind,[ind1,ind2,ind3,ind4,sa])

    Aup = prime(Adown, ind1, ind2, ind3, ind4, ia)

    AA = Aup*dag(Adown)

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

    mm = ((c1*t1)*t4*local_mm*c2*t2*c3*t3*c4)/(((((((((c1*t1)*t4)*aa)*c2)*t2)*c3)*t3)*c4))
    mm = array(mm)
    lambda = eigvals(mm)

    return lambda

end