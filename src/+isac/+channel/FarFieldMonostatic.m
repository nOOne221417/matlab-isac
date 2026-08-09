classdef FarFieldMonostatic
    properties (SetAccess = private)
        ReceiveArray    % Array (ULA)
        Waveform        % Waveform (OFDM)
    end

    methods 
        function obj = FarFieldMonostatic(ReceiveArray, waveform)
            obj.ReceiveArray = ReceiveArray;
            obj.Waveform = waveform;
        end

        function Y = createY(obj, targets, Bnp, SNRdB)
            Mr = obj.ReceiveArray.NumEs;
            N = obj.Waveform.NumSubcarriers;
            P = obj.Waveform.NumSymbols;

            Y_tamp = complex(zeros(Mr, N, P));

            subcarrierIndex = 0:N-1;
            symbolIndex = 0:P-1;

            for k = 1:numel(targets)
                target = targets(k);
                thetaRed = target.angleRad();

                tau = target.delayS();
                nu = target.dopplerHz(obj.ReceiveArray.WavelengthM);
                ak = target.Ref;
                aR = obj.ReceiveArray.steeringVector(thetaRed);

                delay = exp(-1j * 2 * pi * subcarrierIndex ...
                        * obj.Waveform.DeltaFHz * tau);
                doppler = exp(1j * 2 * pi * symbolIndex ...
                        * obj.Waveform.SymDur * nu);

                DaD = delay.' * doppler;
                Y_tamp = Y_tamp + ak * reshape(aR, Mr, 1, 1) ...
                    .* reshape(Bnp .* DaD, 1, N, P);
            end

            signalPower = mean(abs(Y_tamp(:)).^2);
            noisePower = signalPower / (10^(SNRdB/10));

            noise = sqrt(noisePower/2) * (randn(size(Y_tamp)) + 1j * randn(size(Y_tamp)));
            Y = Y_tamp + noise;
        end
    end
end
