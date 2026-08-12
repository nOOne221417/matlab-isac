classdef ESPRIT
    methods(Static)
        function result = estimateAoA(Y, array, K, noisePower, noiseAlpha, modelOrderMethod)
            arguments
                Y {mustBeNumeric}
                array (1,1) isac.array.ULA
                K (1,1) double {mustBeInteger, mustBeNonnegative} = 0;
                noisePower (1,1) double {mustBeNonnegative} = 0;
                noiseAlpha (1,1) double {mustBePositive} = 1.2;
                modelOrderMethod (1,1) string {mustBeMember(modelOrderMethod, ...
                    ["noisePower", "mdl"])} = "mdl";
            end

            Mr = array.NumEs;
            lambda = array.WavelengthM;
            d = array.SpacingM;

            Xtheta = reshape(Y, Mr, []);  % 快拍
            Rx = Xtheta * Xtheta' / size(Xtheta, 2);  % 协方差矩阵
            [V, D] = eig(Rx);
            [eigenvalues, idx] = sort(real(diag(D)), 'descend');
            if K == 0
                K = ESPRIT.resolveModelOrder(eigenvalues, size(Xtheta, 2), ...
                    modelOrderMethod, noisePower, noiseAlpha, [], Mr - 1);
            end

            Es = V(:, idx(1:K));  % 信号子空间
            Es1 = Es(1:end-1, :);
            Es2 = Es(2:end, :);

            Psi = pinv(Es1) * Es2;
            eigVals  = eig(Psi);
            aoaDeg = asind(-lambda * angle(eigVals) / (2*pi*d));

            result = aoaDeg;
        end

        function result = estimateDoppler(Y, Bnp, waveform, K, noisePower, noiseAlpha, modelOrderMethod)
            arguments
                Y {mustBeNumeric}
                Bnp {mustBeNumeric}
                waveform (1,1) isac.waveform.OFDM
                K (1,1) double {mustBeInteger, mustBeNonnegative} = 0;
                noisePower (1,1) double {mustBeNonnegative} = 0;
                noiseAlpha (1,1) double {mustBePositive} = 1.2;
                modelOrderMethod (1,1) string {mustBeMember(modelOrderMethod, ...
                    ["noisePower", "mdl"])} = "mdl";
            end
            P = waveform.NumSymbols;
            N = waveform.NumSubcarriers;

            Y_norm = Y ./ reshape(Bnp, 1, N, P);
            Xnu = reshape(permute(Y_norm, [3 2 1]), P, []);  % 快拍
            Rx = Xnu * Xnu' / size(Xnu, 2);  % 协方差矩阵
            [V, D] = eig(Rx);
            [eigenvalues, idx] = sort(real(diag(D)), 'descend');
            if K == 0
                K = ESPRIT.resolveModelOrder(eigenvalues, size(Xnu, 2), ...
                    modelOrderMethod, noisePower, noiseAlpha, Bnp, P - 1);
            end

            Es = V(:, idx(1:K));  % 信号子空间
            Es1 = Es(1:end-1, :);
            Es2 = Es(2:end, :);

            Psi = pinv(Es1) * Es2;
            eigVals  = eig(Psi);
            dopplerHz = angle(eigVals) / (2*pi*waveform.SymDur);

            result = dopplerHz;
        end

        function result = estimateRange(Y, Bnp, waveform, K, noisePower, noiseAlpha, modelOrderMethod)
            arguments
                Y {mustBeNumeric}
                Bnp {mustBeNumeric}
                waveform (1,1) isac.waveform.OFDM
                K (1,1) double {mustBeInteger, mustBeNonnegative} = 0;
                noisePower (1,1) double {mustBeNonnegative} = 0;
                noiseAlpha (1,1) double {mustBePositive} = 1.2;
                modelOrderMethod (1,1) string {mustBeMember(modelOrderMethod, ...
                    ["noisePower", "mdl"])} = "mdl";
            end
            
            P = waveform.NumSymbols;
            N = waveform.NumSubcarriers;

            Y_norm = Y ./ reshape(Bnp, 1, N, P);
            Xtau = reshape(permute(Y_norm, [2 1 3]), N, []);  % 快拍
            Rx = Xtau * Xtau' / size(Xtau, 2);  % 协方差矩阵
            [V, D] = eig(Rx);
            [eigenvalues, idx] = sort(real(diag(D)), 'descend');
            if K == 0
                K = ESPRIT.resolveModelOrder(eigenvalues, size(Xtau, 2), ...
                    modelOrderMethod, noisePower, noiseAlpha, Bnp, N - 1);
            end

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

    methods(Static, Access = private)
        function K = resolveModelOrder(eigenvalues, numSnapshots, method, noisePower, noiseAlpha, Bnp, maxK)
            switch method
                case "noisePower"
                    effectiveNoisePower = ESPRIT.effectiveNoisePower(noisePower, Bnp);
                    K = ESPRIT.noisePowerMethod(eigenvalues, effectiveNoisePower, noiseAlpha, maxK);
                case "mdl"
                    K = ESPRIT.mdlMethod(eigenvalues, numSnapshots, maxK);
            end
        end

        function K = noisePowerMethod(eigenvalues, noisePower, noiseAlpha, maxK)
            arguments
                eigenvalues (:,1) double
                noisePower (1,1) double {mustBePositive}
                noiseAlpha (1,1) double {mustBePositive}
                maxK (1,1) double {mustBeInteger, mustBeNonnegative}
            end
            K = sum(eigenvalues(1:maxK) > noiseAlpha * noisePower);
        end

        function K = mdlMethod(eigenvalues, numSnapshots, maxK)
            dimension = numel(eigenvalues);
            mdl = zeros(maxK + 1, 1);
            for k = 0:maxK
                noiseEigenvalues = max(eigenvalues(k+1:dimension), eps);
                numNoiseEigenvalues = dimension - k;
                geometricMean = exp(mean(log(noiseEigenvalues)));
                logLikelihood = -numSnapshots * numNoiseEigenvalues * log( ...
                    geometricMean / mean(noiseEigenvalues));
                penalty = 0.5 * k * (2 * dimension - k) * log(numSnapshots);
                mdl(k+1) = logLikelihood + penalty;
            end
            [~, index] = min(mdl);
            K = index - 1;
        end

        function noisePower = effectiveNoisePower(noisePower, Bnp)
            arguments
                noisePower (1,1) double {mustBePositive}
                Bnp {mustBeNumeric}
            end
            if ~isempty(Bnp) && any(abs(Bnp(:)) == 0)
                error('isac:estimation:ZeroSymbol', ...
                    'Bnp must not contain zero-valued symbols.');
            end
            if isempty(Bnp)
                return;
            end
            noisePower = noisePower * mean(1 ./ abs(Bnp(:)).^2);
        end
    end
end
