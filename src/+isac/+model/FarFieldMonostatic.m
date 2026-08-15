classdef FarFieldMonostatic < isac.model.ResponseModel
    properties (SetAccess = private)
        ReceiveArray    % Array (ULA)
        Waveform        % Waveform (OFDM)
    end

    methods
        function obj = FarFieldMonostatic(receiveArray, waveform)
            arguments
                receiveArray (1,1) isac.array.ULA
                waveform (1,1) isac.waveform.OFDM
            end

            obj.ReceiveArray = receiveArray;
            obj.Waveform = waveform;
        end

        function Yatom = atom(obj, target, Bnp)
            arguments
                obj
                target (1,1) isac.scene.PointTarget
                Bnp {mustBeNumeric}
            end
            validateattributes(Bnp, {'numeric'}, ...
                {'size', [obj.Waveform.NumSubcarriers, obj.Waveform.NumSymbols]});

            Mr = obj.ReceiveArray.NumEs;
            N = obj.Waveform.NumSubcarriers;
            P = obj.Waveform.NumSymbols;

            subcarrierIndex = 0:N-1;
            symbolIndex = 0:P-1;
            thetaRad = target.angleRad();
            tau = target.delayS();
            nu = target.dopplerHz(obj.ReceiveArray.WavelengthM);
            aR = obj.ReceiveArray.steeringVector(thetaRad);

            delay = exp(-1j * 2 * pi * subcarrierIndex ...
                * obj.Waveform.DeltaFHz * tau);
            doppler = exp(1j * 2 * pi * symbolIndex ...
                * obj.Waveform.SymDur * nu);

            delayDoppler = delay.' * doppler;
            Yatom = reshape(aR, Mr, 1, 1) ...
                .* reshape(Bnp .* delayDoppler, 1, N, P);
        end

        function dimensions = outputSize(obj)
            dimensions = [obj.ReceiveArray.NumEs, ...
                obj.Waveform.NumSubcarriers, obj.Waveform.NumSymbols];
        end

        function distancesM = pathDistancesM(~, targets)
            arguments
                ~
                targets isac.scene.PointTarget
            end
            rangesM = reshape([targets.RangeM], [], 1);
            distancesM = [rangesM, rangesM];
        end

        function value = wavelengthM(obj)
            value = obj.ReceiveArray.WavelengthM;
        end
    end
end
