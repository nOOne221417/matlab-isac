classdef PointTarget
    properties (SetAccess = private)
        AoADeg  % AoA (Deg)
        RangeM  % Range (m)
        velMps  % Velocity (m/s)
        Ref = 1 % Reflection coefficient
    end

    methods
        function obj = PointTarget(aoaDeg, rangeM, velMps, ref)
            arguments
                aoaDeg (1,1) double {mustBeReal}
                rangeM (1,1) double {mustBeNonnegative}
                velMps (1,1) double {mustBeReal}
                ref (1,1) double = 1
            end
            obj.AoADeg = aoaDeg;
            obj.RangeM = rangeM;
            obj.velMps = velMps;
            obj.Ref = ref;
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
