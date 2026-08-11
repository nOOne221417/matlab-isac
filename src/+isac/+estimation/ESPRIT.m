classdef ESPRIT
    methods(Static)
        function result = estimateAoA(Y, array, K)
            Mr = array.NumEs;
            lambda = array.WavelengthM;
            d = array.SpacingM;

            Xtheta = reshape(Y, Mr, []);  % 快拍
            Rx = Xtheta * Xtheta' / size(Xtheta, 2);  % 协方差矩阵
            [V, D] = eig(Rx);
            [~, idx] = sort(diag(D), 'descend');

            Es = V(:, idx(1:K));  % 信号子空间
            Es1 = Es(1:end-1, :);
            Es2 = Es(2:end, :);

            Psi = pinv(Es1) * Es2;
            eigVals  = eig(Psi);
            aoaDeg = asind(-lambda * angle(eigVals) / (2*pi*d));

            result = aoaDeg;
        end

        function result = estimateDoppler(Y, Bnp, waveform, K)
            P = waveform.NumSymbols;
            N = waveform.NumSubcarriers;

            Y_norm = Y ./ reshape(Bnp, 1, N, P);
            Xnu = reshape(permute(Y_norm, [3 2 1]), P, []);  % 快拍
            Rx = Xnu * Xnu' / size(Xnu, 2);  % 协方差矩阵
            [V, D] = eig(Rx);
            [~, idx] = sort(diag(D), 'descend');

            Es = V(:, idx(1:K));  % 信号子空间
            Es1 = Es(1:end-1, :);
            Es2 = Es(2:end, :);

            Psi = pinv(Es1) * Es2;
            eigVals  = eig(Psi);
            dopplerHz = angle(eigVals) / (2*pi*waveform.SymDur);

            result = dopplerHz;
        end

        function result = estimateRange(Y, Bnp, waveform, K)
            P = waveform.NumSymbols;
            N = waveform.NumSubcarriers;

            Y_norm = Y ./ reshape(Bnp, 1, N, P);
            Xtau = reshape(permute(Y_norm, [2 1 3]), N, []);  % 快拍
            Rx = Xtau * Xtau' / size(Xtau, 2);  % 协方差矩阵
            [V, D] = eig(Rx);
            [~, idx] = sort(diag(D), 'descend');

            Es = V(:, idx(1:K));  % 信号子空间
            Es1 = Es(1:end-1, :);
            Es2 = Es(2:end, :);

            Psi = pinv(Es1) * Es2;
            eigVals = eig(Psi);
            delayS = -angle(eigVals) / (2*pi*waveform.DeltaFHz);

            c = 299792458; % 光速
            rangeM = delayS * c / 2;  % 距离轴

            result = rangeM;
        end
    end
end