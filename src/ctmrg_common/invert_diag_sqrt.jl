

function invert_diag_sqrt(lambda::ITensor)

    index = inds(dag(lambda))
    inv_lambda_sqrt = ITensor(index)
    N = size(lambda)[1];

    for i = 1:N
        if lambda[index[1]=>i,index[2]=>i] > 1e-10
            inv_lambda_sqrt[index[1]=>i,index[2]=>i] = 1/sqrt(lambda[index[1]=>i,index[2]=>i]);
        end
    end

    return inv_lambda_sqrt
    
end
