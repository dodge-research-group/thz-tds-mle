function figgen
%FIGGEN  Makes figures for time-domain maximum likelihood paper

disp("Generate Monte Carlo data")
tic;
DataMC = datamc;
toc;

disp("Load experimental data")
tic;
DataExp = dataexp;
toc;

disp("ETFE fit of TF to Monte Carlo data")
tic;
CalcMC = calcmcfit(DataMC);
toc;

disp("MLE fit of TF to Monte Carlo data (takes a while)")
tic;
CalcTD = calctdfit(DataMC);
toc;

disp("MLE fit of noise model to experimental data")
tic;
CalcNoise = calcnoise(DataExp);
toc;

disp("MLE fits of TF models to experimental Data (takes a while)")
tic;
CalcTF = calctfexp(DataExp, CalcNoise);
toc;

Fig(1) = figsignalnoise(DataMC);
figNum(1) = 1;

Fig(2) = figetf(DataMC);
figNum(2) = 2;

Fig(3) = figetfefit(CalcMC);
figNum(3) = 3;

Fig(4) = figmlefit(CalcTD);
figNum(4) = 4;

Fig(5) = fignoise(DataExp, CalcNoise);
figNum(5) = 5;

Fig(6) = figtfexp(CalcTF, 11);
figNum(6) = 6;

nFig = length(figNum);

for iFig = 1:nFig
    print(Fig(iFig), '-depsc2', ['tdmle' num2str(iFig) '.eps']);
end

end