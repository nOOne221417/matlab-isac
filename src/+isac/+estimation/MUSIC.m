classdef MUSIC
    methods (Static)
        function result = estimateAoA(Y, array, K, numfft)
            Mr = array.NumEs;

            Xtheta = reshape(Y, Mr, []);  % 快拍

            Rx = Xtheta * Xtheta' / size(Xtheta, 2); % 协方差矩阵 
            [V,D] = eig(Rx);  % EVD分解
            [~, idx] = sort(diag(D), 'descend');  % 从大到小排序
            En = V(:, idx(K+1:end));  % 噪声子空间
            
            spatialFrequency = (-numfft/2:numfft/2-1).' / numfft;
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
        end

        function result = estimateDoppler(Y, Bnp, waveform, K, numfft)
            P = waveform.NumSymbols;
            N = waveform.NumSubcarriers;
            Y_norm = Y ./ reshape(Bnp, 1, N, P);

            Xnu = reshape(permute(Y_norm, [3 2 1]), P, []);  % 快拍

            Rx = Xnu * Xnu' / size(Xnu, 2); % 协方差矩阵
            [V,D] = eig(Rx);  % EVD分解
            [~, idx] = sort(diag(D), 'descend');  % 从大到小排序
            En = V(:, idx(K+1:end));  % 噪声子空间

            dopplerAxis = (-numfft/2 : numfft/2-1).' / (numfft * waveform.SymDur);   % Doppler轴
            aR = exp(1j * 2 * pi * ((0:P-1).') * waveform.SymDur * dopplerAxis.');  

            Q = En * En';
            denom = sum(conj(aR) .* (Q * aR), 1).';
            spectrum = 1 ./ max(real(denom), eps);
            spectrumNorm = spectrum / max(spectrum);     % 归一化

            result.DopplerHz = dopplerAxis;
            result.Spectrum = spectrumNorm;
            result.PowerDb = 10 * log10(spectrumNorm + eps);
        end

        function result = estimateRange(Y, Bnp, waveform, K, numfft)
            P = waveform.NumSymbols;
            N = waveform.NumSubcarriers;
            Y_norm = Y ./ reshape(Bnp, 1, N, P);

            Xtau = reshape(permute(Y_norm, [2 1 3]), N, []); % 快拍

            Rx = Xtau * Xtau' / size(Xtau, 2); % 协方差矩阵
            [V,D] = eig(Rx);  % EVD分解
            [~, idx] = sort(diag(D), 'descend');  % 从大到小排序
            En = V(:, idx(K+1:end));  % 噪声子空间

            delayAxis = (0:numfft-1).' / (numfft * waveform.DeltaFHz);  % delay轴
            aR = exp(-1j * 2 * pi * ((0:N-1).') * waveform.DeltaFHz * delayAxis.');  

            Q = En * En';
            denom = sum(conj(aR) .* (Q * aR), 1).';
            spectrum = 1 ./ max(real(denom), eps);
            spectrumNorm = spectrum / max(spectrum);     % 归一化

            c = 299792458; % 光速

            result.RangeM = delayAxis * c / 2;  
            result.Spectrum = spectrumNorm;
            result.PowerDb = 10 * log10(spectrumNorm + eps);
        end
    end
end