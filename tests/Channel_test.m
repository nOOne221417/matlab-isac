clear; clc;
rng(23);

c = 299792458;
fcHz = 28e9;
lambdaM = c / fcHz;

array = isac.array.ULA(4, lambdaM / 2, lambdaM);
waveform = isac.waveform.OFDM(8, 4, 120e3, 1/4);
target = isac.scene.PointTarget(10, 20, 5, 0.8);
Bnp = ones(waveform.NumSubcarriers, waveform.NumSymbols);

responseModel = isac.model.FarFieldMonostatic(array, waveform);
expectedY = target.Ref * responseModel.atom(target, Bnp);

channel = isac.channel.SensingChannel(responseModel);
[Y, info] = channel.createY(target, Bnp, Inf);
assert(isequal(Y, expectedY));
assert(info.NoisePower == 0);
assert(channel.noisePower == 0);
assert(info.FadingGains == 1);
assert(info.LossGains == 1);
assert(isa(info.ResponseModel, 'isac.model.FarFieldMonostatic'));
assert(isa(info.FadingModel, 'isac.channel.fading.NoFading'));
assert(isa(info.LossModel, 'isac.channel.loss.NoLoss'));
assert(isa(info.NoiseModel, 'isac.channel.noise.AWGN'));

rng(99);
[Y, info] = channel.createY(target, Bnp, 10);
rng(99);
expectedNoisePower = mean(abs(expectedY(:)).^2) / 10;
expectedNoise = sqrt(expectedNoisePower / 2) ...
    * (randn(size(expectedY)) + 1j * randn(size(expectedY)));
assert(isequal(Y, expectedY + expectedNoise));
assert(info.NoisePower == expectedNoisePower);

noiseFreeChannel = isac.channel.SensingChannel(responseModel, ...
    NoiseModel=isac.channel.noise.NoNoise());
[Y, info] = noiseFreeChannel.createY(target, Bnp, 10);
assert(isequal(Y, expectedY));
assert(info.NoisePower == 0);

rayleigh = isac.channel.fading.Rayleigh();
rayleighGains = rayleigh.coefficients(repmat(target, 3, 1));
assert(isequal(size(rayleighGains), [3, 1]));
assert(all(isfinite(rayleighGains)));

rician = isac.channel.fading.Rician(4);
ricianGains = rician.coefficients(repmat(target, 3, 1));
assert(isequal(size(ricianGains), [3, 1]));
assert(all(isfinite(ricianGains)));

ricianChannel = isac.channel.SensingChannel( ...
    responseModel, FadingModel=rician);
[Y, info] = ricianChannel.createY(target, Bnp, Inf);
assert(isequal(Y, target.Ref * info.PathGains * responseModel.atom(target, Bnp)));

loss = isac.channel.loss.FreeSpacePathLoss();
lossGains = loss.coefficients(target, responseModel);
expectedLoss = (lambdaM / (4 * pi * target.RangeM))^2;
assert(abs(lossGains - expectedLoss) < eps(expectedLoss) * 4);

lossChannel = isac.channel.SensingChannel(responseModel, LossModel=loss);
[Y, info] = lossChannel.createY(target, Bnp, Inf);
assert(isequal(Y, target.Ref * info.PathGains * responseModel.atom(target, Bnp)));

simulator = isac.simulation.MonostaticSimulator(array, waveform, target, 10);
data = simulator.run();
assert(data.noisePower > 0);
assert(data.ChannelInfo.NoisePower == data.noisePower);
