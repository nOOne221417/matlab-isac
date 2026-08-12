classdef Periodogram
    methods (Static)
        function result = estimateAoA(Y, array, numfft)
            arguments
                Y {mustBeNumeric}
                array (1,1) isac.array.ULA
                numfft (1,1) double {mustBeInteger, mustBePositive} = 1024;
            end

            Mr = array.NumEs;

            Xtheta = reshape(Y, Mr, []);  % 快拍
            spectrum = mean(abs(fftshift(fft(Xtheta, numfft, 1), 1)).^2, 2); % fft
            spectrumNorm = spectrum / max(spectrum);  % 归一化

            spatialFrequency = (-numfft/2:numfft/2-1).' / numfft;
            angleAxis = asind(-spatialFrequency * array.WavelengthM / array.SpacingM);

           % arcsin 是奇函数，映射的时候因为有负号所以反过来
            angleAxis = flipud(angleAxis);
            spectrumNorm = flipud(spectrumNorm);

            result.AngleDeg = angleAxis;
            result.Spectrum = spectrumNorm;
            result.PowerDb = 10 * log10(spectrumNorm + eps);
        end

        function result = estimateDoppler(Y, Bnp, waveform, numfft)
            arguments
                Y {mustBeNumeric}
                Bnp {mustBeNumeric}
                waveform (1,1) isac.waveform.OFDM
                numfft (1,1) double {mustBeInteger, mustBePositive} = 1024;
            end

            P = waveform.NumSymbols;
            N = waveform.NumSubcarriers;
            Y_norm = Y ./ reshape(Bnp, 1, N, P);
            

            Xnu = reshape(permute(Y_norm, [3 2 1]), P, []);  % 快拍
            spectrum = mean(abs(fftshift(fft(Xnu, numfft, 1), 1)).^2, 2); % fft
            spectrumNorm = spectrum / max(spectrum);  % 归一化

            result.DopplerHz = (-numfft/2 : numfft/2-1).' ...
                / (numfft * waveform.SymDur);   % Doppler轴
            result.Spectrum = spectrumNorm;
            result.PowerDb = 10 * log10(spectrumNorm + eps);
        end

        function result = estimateRange(Y, Bnp, waveform, numfft)
            arguments
                Y {mustBeNumeric}
                Bnp {mustBeNumeric}
                waveform (1,1) isac.waveform.OFDM
                numfft (1,1) double {mustBeInteger, mustBePositive} = 1024;
            end
            
            P = waveform.NumSymbols;
            N = waveform.NumSubcarriers;
            Y_norm = Y ./ reshape(Bnp, 1, N, P);

            Xtau = reshape(permute(Y_norm, [2 1 3]), N, []); % 快拍
            spectrum = mean(abs(ifft(Xtau, numfft, 1)).^2, 2); % ifft
            spectrumNorm = spectrum / max(spectrum);  % 归一化
            delayAxis = (0:numfft-1).' / (numfft * waveform.DeltaFHz);  % 延迟轴

            c = 299792458; % 光速

            result.RangeM = delayAxis * c / 2;  % 距离轴
            result.Spectrum = spectrumNorm;
            result.PowerDb = 10 * log10(spectrumNorm + eps);
        end
    end
end
