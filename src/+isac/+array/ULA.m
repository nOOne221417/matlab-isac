classdef ULA
    properties (SetAccess = private)
        NumEs       % Number of elements
        SpacingM    % Element spacing (m)
        WavelengthM % Wavelength (m)
    end

    methods
        function obj = ULA(NumElements, SpacingM, WavelengthM)
            obj.NumEs = NumElements;
            obj.SpacingM = SpacingM;
            obj.WavelengthM = WavelengthM;
        end

        function a = steeringVector(obj, thetaRad)
            index = (0:obj.NumEs-1).';
            thetaRad = thetaRad(:).';
            a = exp(-1j * 2 * pi * obj.SpacingM /...
                 obj.WavelengthM * index * sin(thetaRad));
        end
    end
end
