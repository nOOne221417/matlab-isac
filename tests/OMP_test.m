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

aoa = isac.estimation.OMP.estimateAoA(data.Y, array, 3);

range = isac.estimation.OMP.estimateRange(data.Y, data.Symbols, waveform, 3);

doppler = isac.estimation.OMP.estimateDoppler(data.Y, data.Symbols, waveform, 3);

expectedAoADeg = sort([targets.AoADeg].');
expectedRangeM = sort([targets.RangeM].');
expectedDopplerHz = sort(arrayfun(@(target) ...
    target.dopplerHz(lambdaM), targets(:)));

estimatedAoADeg = sort(real(aoa(:)));
estimatedRangeM = sort(real(range(:)));
estimatedDopplerHz = sort(real(doppler(:)));

aoaErrorDeg = estimatedAoADeg - expectedAoADeg;
rangeErrorM = estimatedRangeM - expectedRangeM;
dopplerErrorHz = estimatedDopplerHz - expectedDopplerHz;

disp('OMP AoA estimation results:');
disp(table(expectedAoADeg, estimatedAoADeg, aoaErrorDeg, ...
    'VariableNames', {'ExpectedDeg', 'EstimatedDeg', 'ErrorDeg'}));

disp('OMP range estimation results:');
disp(table(expectedRangeM, estimatedRangeM, rangeErrorM, ...
    'VariableNames', {'ExpectedM', 'EstimatedM', 'ErrorM'}));

disp('OMP Doppler estimation results:');
disp(table(expectedDopplerHz, estimatedDopplerHz, dopplerErrorHz, ...
    'VariableNames', {'ExpectedHz', 'EstimatedHz', 'ErrorHz'}));

assert(numel(estimatedAoADeg) == numel(expectedAoADeg), ...
    'OMP AoA estimator returned an unexpected number of targets.');
assert(numel(estimatedRangeM) == numel(expectedRangeM), ...
    'OMP range estimator returned an unexpected number of targets.');
assert(numel(estimatedDopplerHz) == numel(expectedDopplerHz), ...
    'OMP Doppler estimator returned an unexpected number of targets.');

assert(max(abs(aoaErrorDeg)) < 2.0, ...
    'OMP AoA error exceeds 2 degrees.');
assert(max(abs(rangeErrorM)) < 5.0, ...
    'OMP range error exceeds 5 meters.');
assert(max(abs(dopplerErrorHz)) < 100.0, ...
    'OMP Doppler error exceeds 100 Hz.');

fprintf('OMP test passed: max AoA error = %.3f deg, max range error = %.3f m, max Doppler error = %.3f Hz.\n', ...
    max(abs(aoaErrorDeg)), max(abs(rangeErrorM)), max(abs(dopplerErrorHz)));

grid on;
