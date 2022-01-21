%function plotAcquisition(acqResults)
Acq = load('GPSdata-DiscreteComponents-fs38_192-if9_55_cpp.mat');
%% Plot all results =======================================================
figure(101);

hAxes = newplot();

bar(hAxes, Acq.peakMetric);

title (hAxes, 'Acquisition results');
xlabel(hAxes, 'PRN number (no bar - SV is not in the acquisition list)');
ylabel(hAxes, 'Acquisition Metric');

oldAxis = axis(hAxes);
axis  (hAxes, [0, 33, 0, oldAxis(4)]);
set   (hAxes, 'XMinorTick', 'on');
set   (hAxes, 'YGrid', 'on');

%% Mark acquired signals ==================================================

acquiredSignals = Acq.peakMetric .* (Acq.carrFreq > 0);

hold(hAxes, 'on');
bar (hAxes, acquiredSignals, 'FaceColor', [0 0.8 0]);
hold(hAxes, 'off');

legend(hAxes, 'Not acquired signals', 'Acquired signals');
