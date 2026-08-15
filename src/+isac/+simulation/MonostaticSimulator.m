classdef MonostaticSimulator < handle
    properties (SetAccess = private)
        ReceiveArray    % aR (ULA)
        Waveform        % OFDM
        Targets         % 目标
        Channel         % SensingChannel
        SNRdB           % 信噪比
    end

    methods
        function obj = MonostaticSimulator(receiveArray, waveform, targets, snrdb, options)
            arguments
                receiveArray (1,1) isac.array.ULA
                waveform (1,1) isac.waveform.OFDM
                targets isac.scene.PointTarget
                snrdb (1,1) double {mustBeReal}
                options.Channel (1,1) isac.channel.SensingChannel = ...
                    isac.channel.SensingChannel( ...
                        isac.model.FarFieldMonostatic(receiveArray, waveform))
            end
            obj.ReceiveArray = receiveArray;
            obj.Waveform = waveform;
            obj.Targets = targets;
            obj.SNRdB = snrdb;
            obj.Channel = options.Channel;
        end

        function data = run(obj)
            Bnp = obj.Waveform.QPSK();
            [Y, channelInfo] = obj.Channel.createY(obj.Targets, Bnp, obj.SNRdB);

            data.Y = Y;
            data.noisePower = obj.Channel.noisePower;
            data.ChannelInfo = channelInfo;
            data.Symbols = Bnp;
            data.Array = obj.ReceiveArray;
            data.Waveform = obj.Waveform;
            data.Targets = obj.Targets;
        end
    end
end
