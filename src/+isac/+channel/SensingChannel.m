classdef SensingChannel < handle
    properties (SetAccess = private)
        ResponseModel
        FadingModel
        NoiseModel
        LossModel
        noisePower = 0
    end

    methods
        function obj = SensingChannel(responseModel, options)
            arguments
                responseModel (1,1) isac.model.ResponseModel
                options.FadingModel (1,1) isac.channel.fading.FadingModel = ...
                    isac.channel.fading.NoFading()
                options.LossModel (1,1) isac.channel.loss.LossModel = ...
                    isac.channel.loss.NoLoss()
                options.NoiseModel (1,1) isac.channel.noise.NoiseModel = ...
                    isac.channel.noise.AWGN()
            end
            obj.ResponseModel = responseModel;
            obj.FadingModel = options.FadingModel;
            obj.NoiseModel = options.NoiseModel;
            obj.LossModel = options.LossModel;
        end

        function [Y, info] = createY(obj, targets, Bnp, SNRdB)
            arguments
                obj
                targets isac.scene.PointTarget
                Bnp {mustBeNumeric}
                SNRdB (1,1) double {mustBeReal}
            end

            lossGains = obj.LossModel.coefficients(targets, obj.ResponseModel);
            fadingGains = obj.FadingModel.coefficients(targets);
            pathGains = fadingGains(:) .* lossGains(:);
            cleanY = complex(zeros(obj.ResponseModel.outputSize()));
            for k = 1:numel(targets)
                Yatom = obj.ResponseModel.atom(targets(k), Bnp);
                cleanY = cleanY + targets(k).Ref * pathGains(k) * Yatom;
            end
            [Y, obj.noisePower] = obj.NoiseModel.apply(cleanY, SNRdB);

            info.CleanY = cleanY;
            info.FadingGains = fadingGains;
            info.LossGains = lossGains;
            info.PathGains = pathGains;
            info.NoisePower = obj.noisePower;
            info.ResponseModel = obj.ResponseModel;
            info.FadingModel = obj.FadingModel;
            info.LossModel = obj.LossModel;
            info.NoiseModel = obj.NoiseModel;
        end
    end
end
