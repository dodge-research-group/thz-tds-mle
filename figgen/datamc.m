function Data = datamc(varargin)
%DATAMC generates Monte Carlo data for demonstration
%

%% Set defaults
DefaultP.N = 256;
DefaultP.A = 1;     % Amplitude [nA]
DefaultP.T=0.05;    % Sampling time [ps]
DefaultP.t0=2.5;      % Peak pulse time [ps]
DefaultP.w=0.25;    % Pulse width [ps]

DefaultP.sigmaAlpha=1e-4;   % Additive noise amplitude [relative to peak]
DefaultP.sigmaBeta=0.01;       % Multiplicative noise amplitude [-]
DefaultP.sigmaTau=1e-3;        % Time base noise amplitude [ps]

DefaultP.Nmc = 500;
DefaultP.seed = 0;

%% Parse inputs
Input = inputParser;
addOptional(Input, 'P' ,DefaultP, @(x) any([isstruct(x), isempty(x)]));
parse(Input, varargin{:});
P = Input.Results.P;

%% Set constants
if isempty(P)
    P = DefaultP;
end
N = P.N;
T = P.T;
t0 = P.t0;

sigmaAlpha = P.sigmaAlpha;
sigmaBeta = P.sigmaBeta;
sigmaTau = P.sigmaTau;
sigmaVec = [sigmaAlpha, sigmaBeta, sigmaTau];

Nmc = P.Nmc;
seed = P.seed;

%% Run simulation

rng(seed)

[y, t]=thzgen(N, T, t0, 'taur', 0.4);

sigmaT = noiseamp(sigmaVec, y, T);
ym = repmat(y, 1, Nmc) + repmat(sigmaT, 1, Nmc).*randn(N, Nmc);

f = fftfreq(N, T);
Nf = floor(N/2)+1;
Ym = fft(ym);
YmRatio = Ym(:,2:2:end)./Ym(:,1:2:end);

Vr = var(real(YmRatio), 0, 2);
Vi = var(imag(YmRatio), 0, 2);
V = (Vr + Vi);

Data.t = t;
Data.y0 = y;
Data.ym = ym;
Data.P = P;
Data.f = f;
Data.Ym = Ym;
Data.YmRatio = YmRatio;
Data.Vr = Vr;
Data.Vi = Vi;
Data.V = V;
Data.Nf = Nf;

end