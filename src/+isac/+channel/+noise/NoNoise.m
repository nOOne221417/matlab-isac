classdef NoNoise < isac.channel.noise.NoiseModel
    methods
        function [Y, noisePower] = apply(~, cleanY, SNRdB)
            arguments
                ~
                cleanY {mustBeNumeric}
                SNRdB (1,1) double {mustBeReal}
            end
            if isnan(SNRdB)
                error('isac:channel:InvalidSNR', 'SNRdB must not be NaN.');
            end
            Y = cleanY;
            noisePower = 0;
        end
    end
end
