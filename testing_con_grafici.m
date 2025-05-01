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
    f = uifigure('Name', 'AUTOVALORI SPAZIATI');
    % Layout a griglia che riempie tutta la figura
    mainLayout = uigridlayout(f, [1, 1]);
    mainLayout.RowHeight = {'1x'};
    mainLayout.ColumnWidth = {'1x'};
    % Tab group che occupa tutto il layout
    tg = uitabgroup(mainLayout);

    fprintf(" === AUTOVALORI SPAZIATI == \n\n");
    fprintf(" === TEST 1 : MATRICE NON DIAGONALIZZABILE | MOLTEPLICITA' BASSA === \n\n");

    AS1 = [25, 30, 35, 40];
    algMultAS1 = [3, 2, 2, 3];
    geoMultAS1 = [3, 1, 1, 2];

    MAS1 = matriceTest(AS1, geoMultAS1, algMultAS1);
    fprintf("Matrice di Testing Autovalori Spaziati : \n");
    disp(MAS1);

    testDiagonalizzabile(MAS1, AS1, algMultAS1, geoMultAS1, 'spaziati', tg, 'Matrice S1');
    
    fprintf("\n\n\n === TEST 2 : MATRICE DIAGONALIZZABILE | MOLTEPLICITA' BASSA === \n\n");

    AS2 = [10, 20, 30, 40];
    algMultAS2 = [1, 3, 2, 1];
    geoMultAS2 = [1, 3, 2, 1];

    MAS2 = matriceTest(AS2, geoMultAS2, algMultAS2);
    fprintf("Matrice di Testing Autovalori Spaziati : \n");
    disp(MAS2);

    testDiagonalizzabile(MAS2, AS2, algMultAS2, geoMultAS2, 'spaziati', tg, 'Matrice S2');
    
    fprintf("\n\n\n === TEST 3 : MATRICE NON DIAGONALIZZABILE | MOLTEPLICITA' ALTA === \n\n");
    
    AS3 = [12, 16, 20, 24];
    algMultAS3 = [5, 6, 7, 6];
    geoMultAS3 = [5, 6, 5, 4];

    MAS3 = matriceTest(AS3, geoMultAS3, algMultAS3);

    fprintf("Matrice di Testing Autovalori Spaziati : \n");
    disp(MAS3);

    testDiagonalizzabile(MAS3, AS3, algMultAS3, geoMultAS3, 'spaziati', tg, 'Matrice S3');

    fprintf("\n\n\n === TEST 4 : MATRICE DIAGONALIZZABILE | MOLTEPLICITA' ALTA === \n\n");
    
    AS4 = [14, 19, 24, 29];
    algMultAS4 = [6, 6, 5, 7];
    geoMultAS4 = [6, 6, 5, 7];

    MAS4 = matriceTest(AS4, geoMultAS4, algMultAS4);

    fprintf("Matrice di Testing Autovalori Spaziati : \n");
    disp(MAS4);

    testDiagonalizzabile(MAS4, AS4, algMultAS4, geoMultAS4, 'spaziati', tg, 'Matrice S4');
end
%% funzione di testing per autovalori clusterizzati
function testAutovaloriClusterizzati()

    f = uifigure('Name', 'AUTOVALORI CLUSTER');
    % Layout a griglia che riempie tutta la figura
    mainLayout = uigridlayout(f, [1, 1]);
    mainLayout.RowHeight = {'1x'};
    mainLayout.ColumnWidth = {'1x'};
    % Tab group che occupa tutto il layout
    tg = uitabgroup(mainLayout);
    
    fprintf(" ===  AUTOVALORI CLUSTER === \n\n");
    fprintf(" === TEST 1 : MATRICE NON DIAGONALIZZABILE | MOLTEPLICITA' BASSA === \n\n");

    AC1 = [4.103, 4.109, 4.111, 4.112];
    algMultAC1 = [3, 3, 1, 2];
    geoMultAC1 = [2, 2, 1, 2];

    MAC1 = matriceTest(AC1, geoMultAC1, algMultAC1);
    fprintf("Matrice di Testing Autovalori Clusterizzati : \n");
    disp(MAC1);
    testDiagonalizzabile(MAC1, AC1, algMultAC1, geoMultAC1, 'cluster', tg, 'Matrice C1');

    fprintf(" === TEST 2 : MATRICE DIAGONALIZZABILE | MOLTEPLICITA' BASSA === \n\n");

    AC2 = [5.261, 5.262, 5.264, 5.267];
    algMultAC2 = [1, 1, 3, 2];
    geoMultAC2 = [1, 1, 3, 2];

    MAC2 = matriceTest(AC2, geoMultAC2, algMultAC2);
    fprintf("Matrice di Testing Autovalori Clusterizzati : \n");
    disp(MAC2);
    testDiagonalizzabile(MAC2, AC2, algMultAC2, geoMultAC2, 'cluster', tg, 'Matrice C2');

    fprintf(" === TEST 3 : MATRICE NON DIAGONALIZZABILE | MOLTEPLICITA' ALTA === \n\n");

    AC3 = [1.724, 1.725, 1.726, 1.727];
    algMultAC3 = [4, 6, 7, 5];
    geoMultAC3 = [4, 4, 4, 4];

    MAC3 = matriceTest(AC3, geoMultAC3, algMultAC3);
    fprintf("Matrice di Testing Autovalori Clusterizzati : \n");
    disp(MAC3);
    testDiagonalizzabile(MAC3, AC3, algMultAC3, geoMultAC3, 'cluster', tg, 'Matrice C3');

    fprintf(" === TEST 4 : MATRICE DIAGONALIZZABILE | MOLTEPLICITA' ALTA === \n\n");

    AC4 = [8.675, 8.676, 8.678, 8.680];
    algMultAC4 = [5, 6, 7, 5];
    geoMultAC4 = [5, 6, 7, 5];

    MAC4 = matriceTest(AC4, geoMultAC4, algMultAC4);
    fprintf("Matrice di Testing Autovalori Clusterizzati : \n");
    disp(MAC4);
    testDiagonalizzabile(MAC4, AC4, algMultAC4, geoMultAC4, 'cluster', tg, 'Matrice C4');
end
%% funzione della verifica della diagnoalizzabilità
function testDiagonalizzabile(A, autovalori, algMult, geoMult, tipologia, tabGroup, tabTitle)

    fprintf("Verifica della diagonalizzabilità:\n\n");
    diag_ok = true;
    it = 15;
    maxit = 50;

    [tol, shift] = calcoloTolShift(A, autovalori, tipologia, algMult);
    
    n = length(autovalori);
    autovalori_calcolati = zeros(1, n);
    alg_attesa = zeros(1, n);
    alg_calcolata = zeros(1, n);
    geo_attesa = zeros(1, n);
    geo_calcolata = zeros(1, n);

    for i = 1:n
        lambda = autovalori(i);

        [lamdacal, algCal, flag] = multAlg(A, lambda + shift, tol, it, maxit);
        geoCal = multGeo(A, lambda, tol);
        
        fprintf("Flag: %d\n", flag);
        fprintf("Autovalore: %.3f | calcolato: %.3f\n", lambda, lamdacal);
        fprintf(" Algebrica attesa: %d | calcolata: %d\n", algMult(i), algCal);
        fprintf(" Geometrica attesa: %d | calcolata: %d\n", geoMult(i), geoCal);
        
        autovalori_calcolati(i) = lamdacal;
        alg_attesa(i) = algMult(i);
        alg_calcolata(i) = algCal;
        geo_attesa(i) = geoMult(i);
        geo_calcolata(i) = geoCal;

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

    % Grafico   
    t = uitab(tabGroup, 'Title', tabTitle);
    gl = uigridlayout(t, [1,1]);
    ax = uiaxes(gl);
    hold(ax, 'on');

    x = 1:length(autovalori);
    plot(ax, x, alg_attesa, '-o', 'LineWidth', 2, 'DisplayName', 'Algebrica attesa');
    plot(ax, x, alg_calcolata, '--o', 'LineWidth', 2, 'DisplayName', 'Algebrica calcolata');
    plot(ax, x, geo_attesa, '-s', 'LineWidth', 2, 'DisplayName', 'Geometrica attesa');
    plot(ax, x, geo_calcolata, '--s', 'LineWidth', 2, 'DisplayName', 'Geometrica calcolata');

    title(ax, sprintf('Molteplicità per %s', tabTitle), 'FontSize', 18);
    xlabel(ax, 'Autovalore');
    ylabel(ax, 'Molteplicità');
    legend(ax, 'Location', 'best');
    grid(ax, 'on');

    ax.XTickLabel = arrayfun(@(i) sprintf('R:%.3f S:%.3f', autovalori(i), autovalori_calcolati(i)), 1:n, 'UniformOutput', false);
end
%% Testing
testAutovaloriSpaziati();
fprintf("\n\n\n");
testAutovaloriClusterizzati();