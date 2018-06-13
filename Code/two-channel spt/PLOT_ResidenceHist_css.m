%PLOT_ResidenceHIST.m
%Anders Sejr Hansen, Oct 2015
clear; clc; close all;

%Title of the plots
PlotTitle = {'161119 H2B Cell#6 full data'};

pool_data = 0;
%Plot residence time histogram
%Path = '/Users/AndersSejrHansen/Dropbox/Data/Microscopy/Test_Data/20151009_U2OS_Clone32_500ms/';
%Path = '/Users/AndersSejrHansen/Dropbox/Data/Microscopy/Test_Data/20151009_mESC_Clone87_500ms/Using_10-6_Error/';
%load([Path, 'TrackedParticles.mat']);

load /Volumes/MyPassport/MFM_Data/161119/H2B_Cell6/SlowTrackingData/161119_H2B_Cell6_fulldata.mat
%load ./SlowTrackingData/20160125_C87_133pM_JF549/20160125_C87_133pM_L10-250_500ms_#4.rpt_tracked.mat

%Exposure Time
Exposure = 0.5;
StartFrameForFit = 5;

%Frames resident on DNA
ResidenceFrames = zeros(1,length(trackedPar));
for i=1:length(trackedPar)
    %ResidenceFrames(1,i) = length(trackedPar(i).TimeStamp);
    
    %Account for missed frames
    TempFrames = trackedPar(i).Frame;
    ResidenceFrames(1,i) = max(TempFrames)-min(TempFrames)+1;
end

%Empirical CDF

%Bin into a full histogram:
HistVec = 0:1:max(ResidenceFrames);
HistVecTime = Exposure.*HistVec;
ResidenceProb = histc(ResidenceFrames, HistVec)./length(ResidenceFrames);
ResidenceCDF = zeros(1,length(ResidenceProb));
for i=2:length(ResidenceProb)
    ResidenceCDF(1,i) = sum(ResidenceProb(1:i));
end
%ResidenceProb = ResidenceProb(2:end); %since single count frames is the second element
FractionBound = 1-ResidenceCDF;

%%%%%%%%%%%%%%%%%%%%%%%% CURVE FITTING %%%%%%%%%%%%%%%%%%%%%%%%

%Fit two exponentials to all of the data
yProb = FractionBound(StartFrameForFit:end);
xTime = HistVecTime(StartFrameForFit:end);

f = fittype('F*exp(-a*x) + (1-F)*exp(-b*x)');
[TwoExp_fit, TwoExp_param] = fit(xTime', yProb', f, 'Lower', [0 0.2 0], 'Upper', [1 5 0.2], 'StartPoint', [0.9 2 0.02]); 
TwoExp_CI = confint(TwoExp_fit);

xFit = 0:0.25:max(HistVec);
yFit = TwoExp_fit.F.*exp(-TwoExp_fit.a.*xFit) + (1-TwoExp_fit.F).*exp(-TwoExp_fit.b.*xFit);

Fit1_text(1) = {'Two Exponential fit: F*exp(-a*t) + (1-F)*exp(-b*t)'};
Fit1_text(2) = {['F = ', num2str(TwoExp_fit.F)]};
Fit1_text(3) = {['F (95% CI): [', num2str(TwoExp_CI(1,1)), ';', num2str(TwoExp_CI(2,1)), ']']};
Fit1_text(4) = {['a = ', num2str(TwoExp_fit.a), ' s^-1 or 1/a = ', num2str(1/TwoExp_fit.a), 's']};
Fit1_text(5) = {['a (95% CI): [', num2str(TwoExp_CI(1,2)), ';', num2str(TwoExp_CI(2,2)), ']']};
Fit1_text(6) = {['b = ', num2str(TwoExp_fit.b), ' s^-1 or 1/b = ', num2str(1/TwoExp_fit.b), 's']};
Fit1_text(7) = {['b (95% CI): [', num2str(TwoExp_CI(1,3)), ';', num2str(TwoExp_CI(2,3)), ']']};

%2 Power-Law fit
yPower = FractionBound;
xPower = HistVecTime;
p = fittype('F*x^-a + (1-F)*x^-b');
[TwoPow_fit, TwoPow_param] = fit(xTime', yProb', p, 'Lower', [0 0 0], 'Upper', [1 100 100], 'StartPoint', [0.9 2 0.02]); 
TwoPow_CI = confint(TwoPow_fit);

xPow = 0:0.25:max(HistVec);
yPow = TwoPow_fit.F.*xPow.^(-TwoPow_fit.a) + (1-TwoPow_fit.F).*xPow.^(-TwoPow_fit.b);



%Fit one exponential to only the long-lived component
LongThres = 30;
xLongProb = Exposure*(0:1:length(HistVec(LongThres:end)));
yLongProb = FractionBound(LongThres-1:end);
fLong = fittype('A*exp(-k*x)');
[oneExp_fit, oneExp_param] = fit(xLongProb', yLongProb', fLong, 'Lower', [0 0], 'Upper', [1 0.26], 'StartPoint', [0.05 0.02]); 
oneExp_CI = confint(oneExp_fit);

xFitLong = min(xLongProb):0.25:max(xLongProb);
yFitLong = oneExp_fit.A.*exp(-oneExp_fit.k.*xFitLong);

Fit2_text(1) = {'One Exponential fit: A*exp(-k*t)'};
Fit2_text(2) = {['A = ', num2str(oneExp_fit.A)]};
Fit2_text(3) = {['k = ', num2str(oneExp_fit.k), ' s^-1 or 1/k = ', num2str(1/oneExp_fit.k), 's']};
Fit2_text(4) = {['k (95% CI): [', num2str(TwoExp_CI(1,2)), ';', num2str(TwoExp_CI(2,2)), ']']};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PLOT the results
figure('position',[200 200 1200 250]); %[x y width height]
subplot(1,4,1);
hold on;
for i=2:length(HistVecTime)
    x1 = HistVecTime(1,i-1); x2 = HistVecTime(1,i);
    y = FractionBound(1,i-1);
    patch([x1 x1 x2 x2], [0 y y 0], [237/255, 28/255, 36/255],'LineStyle','none');
end
%[f,x ] = ecdf(ResidenceFrames);
%plot(Exposure*x, 1-f, 'b-');
plot(xFit+Exposure/2, yFit, 'k-', 'LineWidth', 2);
text(2,0.85*max(FractionBound),Fit1_text,'HorizontalAlignment','Left', 'FontSize',9, 'FontName', 'Helvetica');
axis([0 16.2 0 1.01*max(FractionBound)]);
title([PlotTitle, ' residence time'], 'FontSize',10, 'FontName', 'Helvetica');
ylabel('Fraction still bound', 'FontSize',10, 'FontName', 'Helvetica');
xlabel('Time (seconds)', 'FontSize',10, 'FontName', 'Helvetica');
hold off;

subplot(1,4,2);
hold on;
for i=2:length(HistVecTime)
    x1 = HistVecTime(1,i-1); x2 = HistVecTime(1,i);
    y = FractionBound(1,i-1);
    patch([x1 x1 x2 x2], [0 y y 0], [237/255, 28/255, 36/255],'LineStyle','none');
end
plot(xFit, yFit, 'k-', 'LineWidth', 2);
text(0.2*max(HistVecTime),0.065*max(FractionBound),[num2str(length(trackedPar)), ' trajectories'],'HorizontalAlignment','Left', 'FontSize',9, 'FontName', 'Helvetica');
axis([0 1.05*max(HistVecTime) 0 0.075*max(FractionBound)]);
title([PlotTitle, ' zoom-in on right tail'], 'FontSize',10, 'FontName', 'Helvetica');
ylabel('Fraction still bound', 'FontSize',10, 'FontName', 'Helvetica');
xlabel('Time (seconds)', 'FontSize',10, 'FontName', 'Helvetica');
hold off;

subplot(1,4,3);
hold on;
for i=LongThres:length(HistVecTime)
    x1 = HistVecTime(1,i-LongThres+1); x2 = HistVecTime(1,i-LongThres+2);
    y = FractionBound(1,i-1);
    patch([x1 x1 x2 x2], [0 y y 0], [237/255, 28/255, 36/255],'LineStyle','none');
end
plot(xFitLong, yFitLong, 'k-', 'LineWidth', 2);
text(1.2*LongThres*Exposure,0.95*FractionBound(LongThres),Fit2_text,'HorizontalAlignment','Left', 'FontSize',9, 'FontName', 'Helvetica');
axis([0 1.05*max(HistVecTime) 0 1.05*FractionBound(LongThres)]);
title([PlotTitle, ' long-lived molecules'], 'FontSize',10, 'FontName', 'Helvetica');
ylabel('Fraction still bound', 'FontSize', 10, 'FontName', 'Helvetica');
xlabel('Time (seconds)', 'FontSize', 10, 'FontName', 'Helvetica');
hold off;

subplot(1,4,4);
hold on;
%for i=2:length(HistVecTime)
%    x1 = HistVecTime(1,i-1); x2 = HistVecTime(1,i);
%    y = FractionBound(1,i-1);
%    patch([x1 x1 x2 x2], [0 y y 0], [237/255, 28/255, 36/255],'LineStyle','none');
%end
%[f,x ] = ecdf(ResidenceFrames);
%plot(Exposure*x, 1-f, 'b-');
plot(HistVecTime+Exposure/2, FractionBound, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
plot(xFit+Exposure/2, yFit, 'k-', 'LineWidth', 2);
plot(xPow+Exposure/2, yPow, 'k--', 'LineWidth', 2);
%text(2,0.85*max(FractionBound),Fit1_text,'HorizontalAlignment','Left', 'FontSize',9, 'FontName', 'Helvetica');
axis([0.4*Exposure 200 0.001 1.01*max(FractionBound)]);
title([PlotTitle, ' Power Law?'], 'FontSize',10, 'FontName', 'Helvetica');
ylabel('Fraction still bound', 'FontSize',10, 'FontName', 'Helvetica');
xlabel('Time (seconds)', 'FontSize',10, 'FontName', 'Helvetica');
legend('Data', 'Two Exp', 'Power Law');
set(gca,'xscale','log');
set(gca,'yscale','log');
hold off;


