function Calc = calctdfit(varargin)
%CALCTDFIT time-domain fits to Monte Carlo simulation data
%

%% Parse inputs
Input = inputParser;

addOptional(Input, 'Data', struct([]), ...
    @(x) any([isstruct(x), isempty(x)]));

parse(Input, varargin{:});

Data = Input.Results.Data;
if isempty(Data)
    % Set default number of Monte Carlo iterations to 1000 pairs
    DefaultP.N = 256;
    DefaultP.A = 1;     % Amplitude [nA]
    DefaultP.T=0.05;    % Sampling time [ps]
    DefaultP.t0=2.5;	% Peak pulse time [ps]
    DefaultP.w=0.25;    % Pulse width [ps]
    
    DefaultP.sigmaAlpha=1e-4;	% Additive noise amplitude [rel. to peak]
    DefaultP.sigmaBeta=1e-2;	% Multiplicative noise amplitude [-]
    DefaultP.sigmaTau=1e-3;     % Time base noise amplitude [ps]
    
    DefaultP.Nmc = 500;
    DefaultP.seed = 0;
    
    Data = datamc('P',DefaultP);
end

t = Data.t;
N = Data.P.N;
T = Data.P.T;
ym = Data.ym;

ym1 = ym(:,1:2:end);
ym2 = ym(:,2:2:end);

Nmc = floor(Data.P.Nmc/2);

sigmaVec = [Data.P.sigmaAlpha; Data.P.sigmaBeta; Data.P.sigmaTau];

%%
% Transfer function definition and parameters

tfun = @(theta,w) theta(1)*exp(1i*theta(2)*w);

A0 = 1;          % amplitude ratio between pulses
eta0 = 0;           % delay between pulses [T]
theta0 = [A0;eta0]; % Initial parameter vector
Np = length(theta0);

f = fftfreq(N, T);
w = 2*pi*f;

%% Run fits

Fit.x0 = theta0;
Fit.lb = [];
Fit.ub = [];
Fit.solver = 'lsqnonlin';
Fit.options = optimoptions('lsqnonlin',...
    'Display','off',...
    'UseParallel',true);

pFit = zeros(Np,Nmc);
resnorm = zeros(1,Nmc);
residual = zeros(N,Nmc);

Init = cell(Nmc,1);
Diagnostic = struct('exitflag',Init,...
    'jacobian',Init);
cv = zeros(Np,Np,Nmc);

for jMC=1:Nmc
    
    sigmay1 = noiseamp(sigmaVec, ym1(:,jMC), T);
    sigmay2 = noiseamp(sigmaVec, ym2(:,jMC), T);
    
    Fit.objective = @(theta) ...
        costfunlsq(tfun,theta,ym1(:,jMC),ym2(:,jMC),...
        sigmay1,sigmay2,w);
    
    [p, resnorm(jMC), residual(:, jMC),...
        Diagnostic(jMC).exitflag,~,~,...
        Diagnostic(jMC).jacobian] = lsqnonlin(Fit);
    
    pFit(:, jMC) = p(:);
    
    [~, R] = qr(Diagnostic(jMC).jacobian);
    Rinv = R\eye(size(R));
    cv(:,:,jMC) = (Rinv*Rinv');
    
end

[resnorm, idxSort] = sort(resnorm);

Calc.t = t;
Calc.N = N;
Calc.resnorm = resnorm;
Calc.nu = N - Np;
Calc.pFit = pFit(:,idxSort);
Calc.residual = residual(:,idxSort);
Calc.cv = cv(:,:,idxSort);
Calc.idxSort = idxSort;

end