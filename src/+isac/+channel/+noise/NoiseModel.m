classdef (Abstract) NoiseModel
    methods (Abstract)
        [Y, noisePower] = apply(obj, cleanY, SNRdB)
    end
end
