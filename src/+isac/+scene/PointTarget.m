classdef PointTarget
    properties (SetAccess = private)
        AoADeg  % AoA (Deg)
        RangeM  % Range (m)
        velMps  % Velocity (m/s)
        Ref = 1 % Reflection coefficient
    end

    methods
        function obj = PointTarget(aoaDeg, rangeM, velMps, ref)
            obj.AoADeg = aoaDeg;
            obj.RangeM = rangeM;
            obj.velMps = velMps;
            if nargin > 3
                obj.Ref = ref;
            end
        end

        function thetaRed = angleRad(obj)
            thetaRed = deg2rad(obj.AoADeg);
        end

        function tau = delayS(obj)
           c = 299792458; % Speed of light (m/s)
           tau = 2 * obj.RangeM / c;
        end

        function nu = dopplerHz(obj, wavelengthM)
            nu = 2 * obj.velMps / wavelengthM;
        end
    end
end