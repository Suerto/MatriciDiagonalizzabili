%% function helper : calcolo della tolleranza e dello shift
function [tol, shift] = calcoloTolShift(A, autovalori, tipologia, AlgMult)
    mid = min(diff(sort(autovalori)));
    maxAv = max(AlgMult);

    switch tipologia
        case 'spaziati'
            tol = max(1e-8, norm(A) * eps);
            shift = mid / 10;
        case 'cluster'
            tol = 1e-8;
            if maxAv <= 3
                shift = mid / 3;
            elseif maxAv <= 5
                shift = mid / 5;
            else
                shift = mid / 10;
            end
    end
end
%% funzione di testing per autovalori spaziati
function testAutovaloriSpaziati()
    fprintf(" === AUTOVALORI SPAZIATI == \n\n");
    fprintf(" === TEST 1 : MATRICE NON DIAGONALIZZABILE | MOLTEPLICITA' BASSA === \n\n");

    AS1 = [25, 30, 35, 40];
    algMultAS1 = [3, 2, 2, 3];
    geoMultAS1 = [3, 1, 1, 2];

    MAS1 = matriceTest(AS1, geoMultAS1, algMultAS1);
    fprintf("Matrice di Testing Autovalori Spaziati : \n");
    disp(MAS1);

    testDiagonalizzabile(MAS1, AS1, algMultAS1, geoMultAS1, 'spaziati');
    
    fprintf("\n\n\n === TEST 2 : MATRICE DIAGONALIZZABILE | MOLTEPLICITA' BASSA === \n\n");

    AS2 = [10, 20, 30, 40];
    algMultAS2 = [1, 3, 2, 1];
    geoMultAS2 = [1, 3, 2, 1];

    MAS2 = matriceTest(AS2, geoMultAS2, algMultAS2);
    fprintf("Matrice di Testing Autovalori Spaziati : \n");
    disp(MAS2);

    testDiagonalizzabile(MAS2, AS2, algMultAS2, geoMultAS2, 'spaziati');
    
    fprintf("\n\n\n === TEST 3 : MATRICE NON DIAGONALIZZABILE | MOLTEPLICITA' ALTA === \n\n");
    
    AS3 = [12, 16, 20, 24];
    algMultAS3 = [5, 6, 7, 6];
    geoMultAS3 = [5, 6, 5, 4];

    MAS3 = matriceTest(AS3, geoMultAS3, algMultAS3);

    fprintf("Matrice di Testing Autovalori Spaziati : \n");
    disp(MAS3);

    testDiagonalizzabile(MAS3, AS3, algMultAS3, geoMultAS3, 'spaziati');

    fprintf("\n\n\n === TEST 4 : MATRICE DIAGONALIZZABILE | MOLTEPLICITA' ALTA === \n\n");
    
    AS4 = [14, 19, 24, 29];
    algMultAS4 = [6, 6, 5, 7];
    geoMultAS4 = [6, 6, 5, 7];

    MAS4 = matriceTest(AS4, geoMultAS4, algMultAS4);

    fprintf("Matrice di Testing Autovalori Spaziati : \n");
    disp(MAS4);

    testDiagonalizzabile(MAS4, AS4, algMultAS4, geoMultAS4, 'spaziati');
end
%% funzione di testing per autovalori clusterizzati
function testAutovaloriClusterizzati()
    fprintf(" ===  AUTOVALORI CLUSTER === \n\n");
    fprintf(" === TEST 1 : MATRICE NON DIAGONALIZZABILE | MOLTEPLICITA' BASSA === \n\n");

    AC1 = [4.103, 4.109, 4.111, 4.112];
    algMultAC1 = [3, 3, 1, 2];
    geoMultAC1 = [2, 2, 1, 2];

    MAC1 = matriceTest(AC1, geoMultAC1, algMultAC1);
    fprintf("Matrice di Testing Autovalori Clusterizzati : \n");
    disp(MAC1);
    testDiagonalizzabile(MAC1, AC1, algMultAC1, geoMultAC1, 'cluster');

    fprintf(" === TEST 2 : MATRICE DIAGONALIZZABILE | MOLTEPLICITA' BASSA === \n\n");

    AC2 = [5.261, 5.262, 5.264, 5.267];
    algMultAC2 = [1, 1, 3, 2];
    geoMultAC2 = [1, 1, 3, 2];

    MAC2 = matriceTest(AC2, geoMultAC2, algMultAC2);
    fprintf("Matrice di Testing Autovalori Clusterizzati : \n");
    disp(MAC2);
    testDiagonalizzabile(MAC2, AC2, algMultAC2, geoMultAC2, 'cluster');

    fprintf(" === TEST 3 : MATRICE NON DIAGONALIZZABILE | MOLTEPLICITA' ALTA === \n\n");

    AC3 = [1.724, 1.725, 1.726, 1.727];
    algMultAC3 = [4, 6, 7, 5];
    geoMultAC3 = [4, 4, 4, 4];

    MAC3 = matriceTest(AC3, geoMultAC3, algMultAC3);
    fprintf("Matrice di Testing Autovalori Clusterizzati : \n");
    disp(MAC3);
    testDiagonalizzabile(MAC3, AC3, algMultAC3, geoMultAC3, 'cluster');

    fprintf(" === TEST 4 : MATRICE DIAGONALIZZABILE | MOLTEPLICITA' ALTA === \n\n");

    AC4 = [8.675, 8.676, 8.678, 8.680];
    algMultAC4 = [5, 6, 7, 5];
    geoMultAC4 = [5, 6, 7, 5];

    MAC4 = matriceTest(AC4, geoMultAC4, algMultAC4);
    fprintf("Matrice di Testing Autovalori Clusterizzati : \n");
    disp(MAC4);
    testDiagonalizzabile(MAC4, AC4, algMultAC4, geoMultAC4, 'cluster');
end
%% funzione della verifica della diagnoalizzabilità
function testDiagonalizzabile(A, autovalori, algMult, geoMult, tipologia)
    fprintf("Verifica della diagonalizzabilità:\n\n");
    diag_ok = true;
    it = 15;
    maxit = 50;

    [tol, shift] = calcoloTolShift(A, autovalori, tipologia, algMult);
    for i = 1:length(autovalori)
        lambda = autovalori(i);

        [~, algCal, ~] = multAlg(A, lambda + shift, tol, it, maxit);
        geoCal = multGeo(A, lambda, tol);

        fprintf("Autovalore: %.3f\n", lambda);
        fprintf(" Algebrica attesa: %d | calcolata: %d\n", algMult(i), algCal);
        fprintf(" Geometrica attesa: %d | calcolata: %d\n", geoMult(i), geoCal);
        
        if abs(algCal - algMult(i)) > tol
            fprintf("Molteplicità algebrica non corrisponde!\n");
        end
        if abs(geoCal - geoMult(i)) > tol
            fprintf("Molteplicità geometrica non corrisponde!\n");
        end
        if abs(algCal - geoCal) > tol
            fprintf("Le molteplicità algebrica e geometrica sono differenti\n\n");
            diag_ok = false;
        else
            fprintf("Questo autovalore ha molteplicità coincidenti\n\n");
        end
    end

    if diag_ok
        fprintf("Matrice diagonalizzabile\n");
    else
        fprintf("Matrice NON diagonalizzabile\n");
    end
end
%% Testing
testAutovaloriSpaziati();
fprintf("\n\n\n");
testAutovaloriClusterizzati();