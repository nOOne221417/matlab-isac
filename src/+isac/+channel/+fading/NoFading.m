classdef NoFading < isac.channel.fading.FadingModel
    methods
        function gains = coefficients(~, targets)
            arguments
                ~
                targets isac.scene.PointTarget
            end
            gains = ones(numel(targets), 1);
        end
    end
end
