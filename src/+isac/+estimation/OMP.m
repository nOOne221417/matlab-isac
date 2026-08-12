classdef OMP
    methods(Static)
        function result = estimateAoA(Y, array, K, numAngleGrid, noisePower, noiseAlpha)
            arguments
                Y {mustBeNumeric}
                array (1,1) isac.array.ULA
                K (1,1) double {mustBeInteger, mustBeNonnegative} = 0;
                numAngleGrid (1,1) double {mustBeInteger, mustBePositive} = 1024;
                noisePower (1,1) double {mustBeNonnegative} = 0;
                noiseAlpha (1,1) double {mustBePositive} = 1.2;
            end

            Mr = array.NumEs;
            Xtheta = reshape(Y, Mr, []);  % 快拍

            spatialFrequency = (-numAngleGrid/2:numAngleGrid/2-1).' / numAngleGrid;
            angleAxis = asind(-spatialFrequency * array.WavelengthM / array.SpacingM);
            angleAxis = flipud(angleAxis);

            Adic = array.steeringVector(deg2rad(angleAxis)); 
            R = Xtheta; 
            support = [];

            if K == 0
                if noisePower <= 0
                    error('isac:estimation:NoisePowerRequired', ...
                        'noisePower must be positive when K is zero.');
                end
                threshold = noiseAlpha * sqrt(Mr * size(Xtheta, 2) * noisePower);
                while norm(R, 'fro') > threshold && numel(support) < Mr
                    energy = sum(abs(Adic' * R).^2, 2);
                    energy(support) = -Inf;

                    [~, idx] = max(energy);
                    support(end+1) = idx;

                    Aselected = Adic(:, support);
                    R = Xtheta - Aselected * pinv(Aselected) * Xtheta;
                end
            else
                for k = 1:K
                    energy = sum(abs(Adic' * R).^2, 2);
                    energy(support) = -Inf;

                    [~, idx] = max(energy);
                    support(end+1) = idx;

                    Aselected = Adic(:, support);
                    R = Xtheta - Aselected * pinv(Aselected) * Xtheta;
                end
            end
            
            result = angleAxis(support);
        end

        function result = estimateDoppler(Y, Bnp, waveform, K, numDopplerGrid, noisePower, noiseAlpha)
            arguments
                Y {mustBeNumeric}
                Bnp {mustBeNumeric}
                waveform (1,1) isac.waveform.OFDM
                K (1,1) double {mustBeInteger, mustBeNonnegative} = 0;
                numDopplerGrid (1,1) double {mustBeInteger, mustBePositive} = 1024;
                noisePower (1,1) double {mustBeNonnegative} = 0;
                noiseAlpha (1,1) double {mustBePositive} = 1.2;
            end

            P = waveform.NumSymbols;
            N = waveform.NumSubcarriers;
            Y_norm = Y ./ reshape(Bnp, 1, N, P);

            Xnu = reshape(permute(Y_norm, [3 2 1]), P, []);  % 快拍
            
            dopplerAxis = (-numDopplerGrid/2 : numDopplerGrid/2-1).' / (numDopplerGrid * waveform.SymDur);   % Doppler轴
            Adic = exp(1j * 2 * pi * ((0:P-1).') * waveform.SymDur * dopplerAxis.');  
            
            R = Xnu; 
            support = [];

            if K == 0
                effectiveNoisePower = OMP.normalizedNoisePower(noisePower, Bnp);
                threshold = noiseAlpha * sqrt(P * size(Xnu, 2) * effectiveNoisePower);
                while norm(R, 'fro') > threshold && numel(support) < P
                    energy = sum(abs(Adic' * R).^2, 2);
                    energy(support) = -Inf;

                    [~, idx] = max(energy);
                    support(end+1) = idx;

                    Aselected = Adic(:, support);
                    R = Xnu - Aselected * pinv(Aselected) * Xnu;
                end
            else
                for k = 1:K
                    energy = sum(abs(Adic' * R).^2, 2);
                    energy(support) = -Inf;

                    [~, idx] = max(energy);
                    support(end+1) = idx;

                    Aselected = Adic(:, support);
                    R = Xnu - Aselected * pinv(Aselected) * Xnu;
                end
            end
            
            result = dopplerAxis(support);
        end

        function result = estimateRange(Y, Bnp, waveform, K, numRangeGrid, noisePower, noiseAlpha)
            arguments
                Y {mustBeNumeric}
                Bnp {mustBeNumeric}
                waveform (1,1) isac.waveform.OFDM
                K (1,1) double {mustBeInteger, mustBeNonnegative} = 0;
                numRangeGrid (1,1) double {mustBeInteger, mustBePositive} = 1024;
                noisePower (1,1) double {mustBeNonnegative} = 0;
                noiseAlpha (1,1) double {mustBePositive} = 1.2;
            end

            P = waveform.NumSymbols;
            N = waveform.NumSubcarriers;
            Y_norm = Y ./ reshape(Bnp, 1, N, P);

            Xtau = reshape(permute(Y_norm, [2 1 3]), N, []);  % 快拍
            
            delayAxis = (0:numRangeGrid-1).' / (numRangeGrid * waveform.DeltaFHz);
            Adic = exp(-1j * 2 * pi * ((0:N-1).') * ...
                waveform.DeltaFHz * delayAxis.');
            
            R = Xtau; 
            support = [];

            if K == 0
                effectiveNoisePower = OMP.normalizedNoisePower(noisePower, Bnp);
                threshold = noiseAlpha * sqrt(N * size(Xtau, 2) * effectiveNoisePower);
                while norm(R, 'fro') > threshold && numel(support) < N
                    energy = sum(abs(Adic' * R).^2, 2);
                    energy(support) = -Inf;

                    [~, idx] = max(energy);
                    support(end+1) = idx;

                    Aselected = Adic(:, support);
                    R = Xtau - Aselected * pinv(Aselected) * Xtau;
                end
            else
                for k = 1:K
                    energy = sum(abs(Adic' * R).^2, 2);
                    energy(support) = -Inf;

                    [~, idx] = max(energy);
                    support(end+1) = idx;

                    Aselected = Adic(:, support);
                    R = Xtau - Aselected * pinv(Aselected) * Xtau;
                end
            end
            
            c = 299792458;
            result = delayAxis(support) * c / 2;
        end
    end

    methods(Static, Access = private)
        function effectiveNoisePower = normalizedNoisePower(noisePower, Bnp)
            if noisePower <= 0
                error('isac:estimation:NoisePowerRequired', ...
                    'noisePower must be positive when K is zero.');
            end
            if any(abs(Bnp(:)) == 0)
                error('isac:estimation:ZeroSymbol', ...
                    'Bnp must not contain zero-valued symbols.');
            end
            effectiveNoisePower = noisePower * mean(1 ./ abs(Bnp(:)).^2);
        end
    end
end
