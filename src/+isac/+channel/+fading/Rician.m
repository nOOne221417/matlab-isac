classdef Rician < isac.channel.fading.FadingModel
    properties (SetAccess = private)
        KFactor
    end

    methods
        function obj = Rician(kFactor)
            arguments
                kFactor (1,1) double {mustBeNonnegative} = 4
            end
            obj.KFactor = kFactor;
        end

        function gains = coefficients(obj, targets)
            arguments
                obj
                targets isac.scene.PointTarget
            end

            numPaths = numel(targets);
            diffuse = (randn(numPaths, 1) + 1j * randn(numPaths, 1)) / sqrt(2);
            specularScale = sqrt(obj.KFactor / (obj.KFactor + 1));
            diffuseScale = sqrt(1 / (obj.KFactor + 1));
            gains = specularScale + diffuseScale * diffuse;
        end
    end
end
