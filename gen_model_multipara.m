Omega=0.05;
for Lambda=2
    for BathL=50
        for z=1
            for omegaratio=1:2
            gen_model(Omega,omegaratio,Lambda,BathL,z);
            end
        end
    end
end
