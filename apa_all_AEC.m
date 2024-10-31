function [m, w] = apa_all_AEC(x, d, miu, ord, p, dlt, a, h1, Update, sel)
    %  Input:
    %       x   --> far-end signal
    %       d   --> reference signal
    %       miu --> step size
    %       ord --> length of the adaptive filter
    %       p   --> projection order
    %       dlt --> regularization constant
    %       a   --> IPNLMS parameter
    %       h1  --> true impulse response of the echo path
    %  Output:
    %       normalized misalignment in dB

    N = length(x);
    w = zeros(ord, 1);
    x1 = zeros(ord, 1);
    D = zeros(p, 1);
    X = zeros(ord, p);
    P = zeros(ord, p);
    m = zeros(N, 1);
    e = zeros(N, 1);
    S0 = dlt * eye(p, p);
    S = zeros(p, p);

    if Update == 1
        % AMIPAPA
        if sel == 1
            for n = 1:N
               
                x1 = [x(n); x1(1:ord-1)];
                X = [x1, X(:,1:p-1)];
                D = [d(n); D(1:p-1)];
                Y = X' * w;
                E = D - Y;
                kd = (1-a) / (2 * ord) + (1 + a) * abs(w) / (10^-8 + 2 * sum(abs(w)));

                P = [kd .* x1, P(:,1:p-1)];
                t = X' * P(:,1); S(1,:) = t';
                S(:,1) = t';

                if n < 100
                    S(1,1) = S(1,1) + 100 * dlt;
                else
                    S(1,1) = S(1,1) + dlt;
                end
                S(2:p, 2:p) = S0(1:p-1, 1:p-1);
                S0 = S;
                w = w + miu * P * inv(S) * E;
                m(n) = 20 * log10(norm(w - h1) / norm(h1));
                end
            

        % IPAPA
        elseif sel == 2
            for n = 1:N
                x1 = [x(n); x1(1:ord-1)];
                X = [x1, X(:,1:p-1)];
                D = [d(n); D(1:p-1)];
                Y = X' * w;
                E = D - Y;

                kd = (1-a) / (2 * ord) + (1 + a) * abs(w) / (10^-8 + 2 * sum(abs(w)));
                P = [kd .* x1, P(:,1:p-1)];

                S = dlt * eye(p) + X' * P;
                E2 = E' * inv(S);
                w = w + miu * P * E2';
                m(n) = 20 * log10(norm(w - h1) / norm(h1));
            end

        % APA
        elseif sel == 3
            for n = 1:N
                x1 = [x(n); x1(1:ord-1)];
                X = [x1, X(:,1:p-1)];
                D = [d(n); D(1:p-1)];
                E = D - X' * w;

                w = w + miu * X * inv(dlt * eye(p) + X' * X) * E;
                m(n) = 20 * log10(norm(w - h1) / norm(h1));
            end

        % NLMS
        elseif sel == 4
            for n = 1:N
                
                x1 = [x(n); x1(1:ord-1)];
                ep = d(n) - x1' * w;

                w = w + miu / (norm(x1)^2 + 10^-8) * x1 .* ep;
                m(n) = 20 * log10(norm(w - h1) / norm(h1));
            end

        % IPNLMS-l1 
        elseif sel == 5
            for n = 1:N
                x1 = [x(n); x1(1:ord-1)];
                ep = d(n) - x1' * w;

                kd = (1 - a) / (2 * ord) + (1 + a) * abs(w) / (10^-8 + 2 * sum(abs(w)));
                Kd = diag(kd);
                w = w + (miu / (x1' * Kd * x1 + 10^-8 * (1 - a) / (2 * ord))) * Kd * x1 .* ep;

                m(n) = 20 * log10(norm(w - h1) / norm(h1));
            end
        end
    end
end
