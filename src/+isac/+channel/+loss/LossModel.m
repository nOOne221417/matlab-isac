classdef (Abstract) LossModel
    methods (Abstract)
        gains = coefficients(obj, targets, responseModel)
    end
end
