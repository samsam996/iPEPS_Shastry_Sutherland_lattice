

function entropy_lambda(lambda::ITensor)

    s = 0;

    for i = 1:1:size(lambda)[1]
        if lambda[i,i] > 1e-12
            s = s - (lambda[i,i]*log(lambda[i,i]));
        end
    end

    return s
end