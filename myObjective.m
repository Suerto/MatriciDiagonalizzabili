%% Funzione myObjective
function [f, g] = myObjective(z, A)
    n = size(A, 1);
    B = A - z * eye(n);
    
    rankB = rank(B);
    if rankB < n - 1
        f = 0;
        g = 0;
        return;
    end

    [L, U, P] = lu(B);

    s = 0;
    copyP = P;
    for j = 1:n
        if copyP(j, j) ~= 1
            for k = j+1:n
                if copyP(k, j) == 1
                    copyP([j, k], :) = copyP([k, j], :);
                    s = s + 1;
                    break;
                end
            end
       end
    end
    detP = (-1)^s;
    
    detU = prod(diag(U));

    if abs(detU) < 1e-12
        detU = 1e-12;
    end
    
    f = detU / detP;
    
    if abs(f) < 1e-12
        f = 1e-12;
    end

    if rankB < n
        f = 1e-12;
        g = 0;
    else 
        InvB = U \ (L \ P);
    
        if abs(trace(InvB)) < 1e-12
            g = 1e3;
        else 
            g = 1 / trace(InvB);
        end
    end
end