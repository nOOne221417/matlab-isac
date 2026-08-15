classdef (Abstract) ResponseModel
    methods (Abstract)
        Yatom = atom(obj, target, Bnp)
        dimensions = outputSize(obj)
        distancesM = pathDistancesM(obj, targets)
        value = wavelengthM(obj)
    end
end
