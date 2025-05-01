%% Funzione multAlg
function [l, m, flag] = multAlg(A, l0, toll, it, maxit)
    n = size(A, 1);
    l = l0;
    m = 1;
    total_calls = 0;
    history = zeros(it, 1);

    convergeInFase1 = false;
    %% Fase 1: Newton standard
    for k = 1:it
        [f, g] = myObjective(l, A);
        total_calls = total_calls + 1;

        history(k) = abs(g);
        l_new = l + g;

        if abs(l_new - l) < toll && abs(f) < toll && abs(g) < toll
            convergeInFase1 = true;
            break;
        end

        l = l_new;
    end
    %% Fase 2 : Stima di m
    m_est = 2;
    if it >= 2 && history(end-1) > toll && history(end) > 0
        ratio = abs(history(end)) / abs(history(end-1));
        if ratio < 1 && abs(1 - ratio) > eps
            m_est = max(2, round(1 / (1 - ratio)));
        end
    end

    fprintf('Stima di m: %d\n', m_est);

    if convergeInFase1
        l = l_new;
        m = 1;
        flag = 1;
        return;
    end

    %% Fase 3: Newton modificato
    flag = 0;
    for m_try = m_est:n
        l_current = l0;
        for step = 1:maxit
            if total_calls >= 10 * maxit
                flag = 0;
                return;
            end

            [~, g] = myObjective(l_current, A);
            total_calls = total_calls + 1;

            l_new = l_current + m_try * g;

            if abs(l_new - l_current) < toll
                l = l_new;
                m = m_try;
                flag = 1;
                return;
            end

            l_current = l_new;
        end
    end
end