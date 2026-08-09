classdef MonostaticSimulator < handle
    properties (SetAccess = private)
        ReceiveArray    % aR (ULA)
        Waveform        % OFDM
        Targets         % 目标
        Channel         % FarFieldMonostatic
        SNRdB           % 信噪比
    end

    methods
        function obj = MonostaticSimulator(receiveArray, waveform, targets, snrdb)
            obj.ReceiveArray = receiveArray;
            obj.Waveform = waveform;
            obj.Targets = targets;
            obj.SNRdB = snrdb;

            obj.Channel = isac.channel.FarFieldMonostatic(receiveArray, waveform);
        
        end

        function data = run(obj)
            Bnp = obj.Waveform.QPSK();
            Y = obj.Channel.createY(obj.Targets, Bnp, obj.SNRdB);

            data.Y = Y;
            data.Symbols = Bnp;
            data.Array = obj.ReceiveArray;
            data.Waveform = obj.Waveform;
            data.Targets = obj.Targets;
        end
    end
end