classdef NoLoss < isac.channel.loss.LossModel
    methods
        function gains = coefficients(~, targets, responseModel)
            arguments
                ~
                targets isac.scene.PointTarget
                responseModel (1,1) isac.model.ResponseModel
            end
            responseModel.pathDistancesM(targets);
            gains = ones(numel(targets), 1);
        end
    end
end
