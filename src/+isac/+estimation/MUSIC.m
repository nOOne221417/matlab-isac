classdef MUSIC
    methods (Static)
        function result = estimateAoA(Y, array, K, numAngleGrid, noisePower, noiseAlpha, modelOrderMethod)
            arguments
                Y {mustBeNumeric}
                array (1,1) isac.array.ULA
                K (1,1) double {mustBeInteger, mustBeNonnegative} = 0;
                numAngleGrid (1,1) double {mustBeInteger, mustBePositive} = 1024;
                noisePower (1,1) double {mustBeNonnegative} = 0;
                noiseAlpha (1,1) double {mustBePositive} = 1.2;
                modelOrderMethod (1,1) string {mustBeMember(modelOrderMethod, ...
                    ["noisePower", "mdl"])} = "mdl";
            end

            Mr = array.NumEs;

            Xtheta = reshape(Y, Mr, []);  % 快拍

            Rx = Xtheta * Xtheta' / size(Xtheta, 2); % 协方差矩阵 
            [V,D] = eig(Rx);  % EVD分解
            [eigenvalues, idx] = sort(real(diag(D)), 'descend');
            if K == 0
                K = MUSIC.resolveModelOrder(eigenvalues, size(Xtheta, 2), ...
                    modelOrderMethod, noisePower, noiseAlpha, [], Mr - 1);
            end
            En = V(:, idx(K+1:end));  % 噪声子空间
            
            spatialFrequency = (-numAngleGrid/2:numAngleGrid/2-1).' / numAngleGrid;
            angleAxis = asind(-spatialFrequency * array.WavelengthM / array.SpacingM);
            angleAxis = flipud(angleAxis);

            aR = array.steeringVector(deg2rad(angleAxis));  

            Q = En * En';
            denom = sum(conj(aR) .* (Q * aR), 1).';      
            spectrum = 1 ./ max(real(denom), eps);
            spectrumNorm = spectrum / max(spectrum);     % 归一化

            result.AngleDeg = angleAxis;
            result.Spectrum = spectrumNorm;
            result.PowerDb = 10 * log10(spectrumNorm + eps);
            result.NumSources = K;
        end

        function result = estimateDoppler(Y, Bnp, waveform, K, numDopplerGrid, noisePower, noiseAlpha, modelOrderMethod)
            arguments
                Y {mustBeNumeric}
                Bnp {mustBeNumeric}
                waveform (1,1) isac.waveform.OFDM
                K (1,1) double {mustBeInteger, mustBeNonnegative} = 0;
                numDopplerGrid (1,1) double {mustBeInteger, mustBePositive} = 1024;
                noisePower (1,1) double {mustBeNonnegative} = 0;
                noiseAlpha (1,1) double {mustBePositive} = 1.2;
                modelOrderMethod (1,1) string {mustBeMember(modelOrderMethod, ...
                    ["noisePower", "mdl"])} = "mdl";
            end

            P = waveform.NumSymbols;
            N = waveform.NumSubcarriers;
            Y_norm = Y ./ reshape(Bnp, 1, N, P);

            Xnu = reshape(permute(Y_norm, [3 2 1]), P, []);  % 快拍

            Rx = Xnu * Xnu' / size(Xnu, 2); % 协方差矩阵
            [V,D] = eig(Rx);  % EVD分解
            [eigenvalues, idx] = sort(real(diag(D)), 'descend');
            if K == 0
                K = MUSIC.resolveModelOrder(eigenvalues, size(Xnu, 2), ...
                    modelOrderMethod, noisePower, noiseAlpha, Bnp, P - 1);
            end
            En = V(:, idx(K+1:end));  % 噪声子空间

            dopplerAxis = (-numDopplerGrid/2 : numDopplerGrid/2-1).' / (numDopplerGrid * waveform.SymDur);   % Doppler轴
            aR = exp(1j * 2 * pi * ((0:P-1).') * waveform.SymDur * dopplerAxis.');  

            Q = En * En';
            denom = sum(conj(aR) .* (Q * aR), 1).';
            spectrum = 1 ./ max(real(denom), eps);
            spectrumNorm = spectrum / max(spectrum);     % 归一化

            result.DopplerHz = dopplerAxis;
            result.Spectrum = spectrumNorm;
            result.PowerDb = 10 * log10(spectrumNorm + eps);
            result.NumSources = K;
        end

        function result = estimateRange(Y, Bnp, waveform, K, numRangeGrid, noisePower, noiseAlpha, modelOrderMethod)
            arguments
                Y {mustBeNumeric}
                Bnp {mustBeNumeric}
                waveform (1,1) isac.waveform.OFDM
                K (1,1) double {mustBeInteger, mustBeNonnegative} = 0;
                numRangeGrid (1,1) double {mustBeInteger, mustBePositive} = 1024;
                noisePower (1,1) double {mustBeNonnegative} = 0;
                noiseAlpha (1,1) double {mustBePositive} = 1.2;
                modelOrderMethod (1,1) string {mustBeMember(modelOrderMethod, ...
                    ["noisePower", "mdl"])} = "mdl";
            end
            
            P = waveform.NumSymbols;
            N = waveform.NumSubcarriers;
            Y_norm = Y ./ reshape(Bnp, 1, N, P);

            Xtau = reshape(permute(Y_norm, [2 1 3]), N, []); % 快拍

            Rx = Xtau * Xtau' / size(Xtau, 2); % 协方差矩阵
            [V,D] = eig(Rx);  % EVD分解
            [eigenvalues, idx] = sort(real(diag(D)), 'descend');
            if K == 0
                K = MUSIC.resolveModelOrder(eigenvalues, size(Xtau, 2), ...
                    modelOrderMethod, noisePower, noiseAlpha, Bnp, N - 1);
            end
            En = V(:, idx(K+1:end));  % 噪声子空间

            rangeAxis = (0:numRangeGrid-1).' / (numRangeGrid * waveform.DeltaFHz);  % Range轴
            aR = exp(-1j * 2 * pi * ((0:N-1).') * waveform.DeltaFHz * rangeAxis.');

            Q = En * En';
            denom = sum(conj(aR) .* (Q * aR), 1).';
            spectrum = 1 ./ max(real(denom), eps);
            spectrumNorm = spectrum / max(spectrum);     % 归一化

            c = 299792458; % 光速

            result.RangeM = rangeAxis * c / 2;
            result.Spectrum = spectrumNorm;
            result.PowerDb = 10 * log10(spectrumNorm + eps);
            result.NumSources = K;
        end
    end

    methods(Static, Access = private)
        function K = resolveModelOrder(eigenvalues, numSnapshots, method, noisePower, noiseAlpha, Bnp, maxK)
            switch method
                case "noisePower"
                    effectiveNoisePower = MUSIC.effectiveNoisePower(noisePower, Bnp);
                    K = MUSIC.noisePowerMethod(eigenvalues, effectiveNoisePower, noiseAlpha, maxK);
                case "mdl"
                    K = MUSIC.mdlMethod(eigenvalues, numSnapshots, maxK);
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
