clear; clc; close all;
rng(23);

c = 299792458;
fcHz = 28e9;
lambdaM = c / fcHz;

array = isac.array.ULA(16, lambdaM / 2, lambdaM);
waveform = isac.waveform.OFDM(128, 64, 120e3, 1/4);

targets = [
    isac.scene.PointTarget(-20, 20, 8, 1.0)
    isac.scene.PointTarget( 10, 80, 12, 0.8)
    isac.scene.PointTarget( 45, 50, 20, 0.6)
];

simulator = isac.simulation.MonostaticSimulator(array, waveform, targets, 10);

data = simulator.run();

aoa = isac.estimation.MUSIC.estimateAoA(data.Y, array, 3, 1024);

range = isac.estimation.MUSIC.estimateRange(data.Y, data.Symbols, waveform, 3, 1024);

doppler = isac.estimation.MUSIC.estimateDoppler(data.Y, data.Symbols, waveform, 3, 1024);

figure;
plot(aoa.AngleDeg, aoa.PowerDb);
xlabel("AoA (deg)");
ylabel("Normalized power (dB)");

figure;
plot(range.RangeM, range.PowerDb);
xlabel("Range (m)");
ylabel("Normalized power (dB)");

figure;
plot(doppler.DopplerHz, doppler.PowerDb);
xlabel("Doppler (Hz)");
ylabel("Normalized power (dB)");

grid on;

