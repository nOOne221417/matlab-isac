classdef OFDM
    properties (SetAccess = private)
        NumSubcarriers      % Number of subcarriers
        NumSymbols          % Number of symbols
        DeltaFHz            % Subcarrier spacing (Hz)
        SymDur              % Symbol duration (s)
        CPDur               % CP duration (s)
        BHz                 % Bandwidth (Hz)
    end

    methods
        function obj = OFDM(NumSubcarriers, NumSymbols, deltaFHz, cpRatio)
            obj.NumSubcarriers = NumSubcarriers;
            obj.NumSymbols = NumSymbols;
            obj.DeltaFHz = deltaFHz;
            time = 1 / deltaFHz;
            obj.CPDur = cpRatio * time;
            obj.SymDur = time + obj.CPDur;
            obj.BHz = NumSubcarriers * deltaFHz;
        end

        function Bnp = QPSK(obj)
            I = 2 * randi([0 1], obj.NumSubcarriers, obj.NumSymbols) - 1; % inPhase symbols
            Q = 2 * randi([0 1], obj.NumSubcarriers, obj.NumSymbols) - 1; % quadrature symbols
            Bnp = (I + 1j * Q) / sqrt(2);           % QPSK symbols
        end
    end
end