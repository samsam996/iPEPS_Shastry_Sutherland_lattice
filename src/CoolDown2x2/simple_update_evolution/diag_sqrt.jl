

function diag_sqrt(lambda::ITensor)

    index = inds(lambda)
    lambda_sqrt = ITensor(index);
    N = size(lambda)[1]; 

    for i = 1:N
        lambda_sqrt[index[1]=>i,index[2]=>i] = sqrt(lambda[index[1]=>i,index[2]=>i]) 
    end

    return lambda_sqrt
    
end