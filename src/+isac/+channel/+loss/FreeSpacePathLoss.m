classdef FreeSpacePathLoss < isac.channel.loss.LossModel
    methods
        function gains = coefficients(~, targets, responseModel)
            arguments
                ~
                targets isac.scene.PointTarget
                responseModel (1,1) isac.model.ResponseModel
            end
            distancesM = responseModel.pathDistancesM(targets);
            if any(distancesM(:) <= 0)
                error('isac:channel:InvalidPathDistance', ...
                    'All propagation distances must be positive.');
            end

            wavelengthM = responseModel.wavelengthM();
            legGains = wavelengthM ./ (4 * pi * distancesM);
            gains = prod(legGains, 2);
        end
    end
end
