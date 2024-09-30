





function simpleupdate(final_temp::Float64,temperature,dbetasu::Float64,J1::Float64,J2::Float64,h::Float64,hs::Float64,N::Int64,D::Int64)

    f(x) = mod(x-1,N) + 1

    Gamma, lambdax, lambday, physical_legs, ancilla_legs, gt, gg = initialisation()

    beta = 0;
    temp =  Inf;
    nsu = 2/dbetasu;
    temperature_save = []

    relevant = 1
    FreeEnergy = 0;

    while final_temp < temp

        tempprev = 1/beta
        beta = beta + dbetasu; 
        temp = 1/beta;
        Gamma,lambdax,lambday, FreeEnergy = simpleupdateJ1(Gamma,lambdax,lambday,physical_legs,gt,nsu,J1,J2,h,hs,D, FreeEnergy)
   
        println(temp)

        if size(temperature)[1] >= relevant && size(temperature)[1] > 0 
            if (temperature[relevant] <= tempprev && temperature[relevant] > temp) 
                
                push!(temperature_save, temp)
                name_data2 = string("Results/LocalTensorsD",string(D),"temp",string(round(temp,digits = 5)),"J",string(J1),"h",string(round(1e2*h,digits = 3)),".jld2")
                save(name_data2,
                    "D",D,
                    "temp",temp,
                    "J1",J1,
                    "J2",J2,
                    "nsu",nsu,
                    "h",h,
                    "hs",hs,
                    "Gamma",Gamma,
                    "lambdax",lambdax,
                    "lambday",lambday,
                    "FreeEnergy", FreeEnergy,
                    "physical_legs",physical_legs,
                    "ancilla_legs",ancilla_legs,
                    "gt",gt,
                "gg",gg)

                relevant = relevant + 1

             
            end 

        end

    end

    return temperature_save  


end