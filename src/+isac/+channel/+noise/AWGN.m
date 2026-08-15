classdef AWGN < isac.channel.noise.NoiseModel
    methods
        function [Y, noisePower] = apply(~, cleanY, SNRdB)
            arguments
                ~
                cleanY {mustBeNumeric}
                SNRdB (1,1) double {mustBeReal}
            end

            signalPower = mean(abs(cleanY(:)).^2);
            noisePower = signalPower / (10^(SNRdB / 10));
            noise = sqrt(noisePower / 2) ...
                * (randn(size(cleanY)) + 1j * randn(size(cleanY)));
            Y = cleanY + noise;
        end
    end
end
