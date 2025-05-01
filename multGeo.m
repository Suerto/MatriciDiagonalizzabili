%% Funzione multGeo
function k = multGeo(A, l, toll)
    [m, n] = size(A);
    if m ~= n
        error('La matrice A deve essere quadrata.');
    end
    
    if toll <= 0
        error('La tolleranza "toll" deve essere positiva.');
    end
    
    B = A - l * eye(n);
    
    [~, U, ~] = lu(B);

    k = sum(abs(diag(U)) < toll);
end