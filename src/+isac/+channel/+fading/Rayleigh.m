classdef Rayleigh < isac.channel.fading.FadingModel
    methods
        function gains = coefficients(~, targets)
            arguments
                ~
                targets isac.scene.PointTarget
            end
            numPaths = numel(targets);
            gains = (randn(numPaths, 1) + 1j * randn(numPaths, 1)) / sqrt(2);
        end
    end
end
