Omega=0.05;
for Lambda=1.1
    for BathL=800
        for z=1
            for omegaratio=0.02:0.02:2
            gen_model(Omega,omegaratio,Lambda,BathL,z);
            end
        end
    end
end