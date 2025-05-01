%% Funzione di generazione della matrice A = Q' * J * Q
function A = matriceTest(autovalori, geoMult, algMult) 
    verificaAutovalori(autovalori, geoMult, algMult);

    J = creaMatriceJordan(autovalori, geoMult, algMult);

    n = size(J, 1);
    [Q, ~] = qr(randn(n));
    
    if norm(Q' * Q - eye(n)) > 1e-10
        warning("Matrice Q non ortogonale");
    end
    A = Q' * J * Q;
end
%% Funzione helper : Verifica della validità dell'autovalore in input
function verificaAutovalori(autovalori, geoMult, algMult) 
    if ~isvector(autovalori) || ~isvector(geoMult) || ~isvector(algMult)
        error("Uno dei parametri forniti non è un vettore, programma terminato");
    end

    if length(autovalori) ~= length(geoMult) || length(autovalori) ~= length(algMult) || length(geoMult) ~= length(algMult)
        error("Uno dei vettori forniti ha una dimensione differente dagli altri, programma terminato");
    end

    if any(algMult < geoMult)
        error("La molteplicità algebrica deve essere maggiore o uguale alla molteplicità geometrica per ogni autovalore.");
    end
end
%% Funzione helper : Costruzione dei blocchi della matrice di Jordan
function bloccoJordan = creaBloccoJordan(autovalore, dim)
        if dim < 1
            error("La dimensione del blocco deve essere maggiore o uguale a 1.");
        end
        
        if dim == 1
            bloccoJordan = autovalore;
        else
            bloccoJordan = autovalore * eye(dim) + diag(ones(dim - 1, 1), 1);
        end
end
%% Funzione helper : Costruzione matrice di Jordan a blocchi
function matriceJordan = creaMatriceJordan(autovalori, geoMult, algMult)
    matriceJordan = [];
    
    for i = 1:length(autovalori)
        alg = algMult(i);
        lambda = autovalori(i);
        geo = geoMult(i);
        
        base_size = floor(alg / geo);
        remainder = mod(alg, geo);
        
        for j = 1:geo
            dim = base_size + (j <= remainder);
            bloccoJordan = creaBloccoJordan(lambda, dim);
            matriceJordan = blkdiag(matriceJordan, bloccoJordan);
        end
    end
end