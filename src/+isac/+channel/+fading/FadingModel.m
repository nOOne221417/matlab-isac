classdef (Abstract) FadingModel
    methods (Abstract)
        gains = coefficients(obj, targets)
    end
end
