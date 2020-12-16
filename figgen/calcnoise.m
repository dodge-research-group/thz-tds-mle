function Calc = calcnoise(varargin)
%CALCNOISE estimates noise parameters from air scans

%% Parse inputs
Input = inputParser;

addOptional(Input, 'Data', struct([]), ...
    @(x) any([isstruct(x), isempty(x)]));

parse(Input, varargin{:});

Data = Input.Results.Data;
if isempty(Data)
    Data = dataexp;
end

t = Data.AirScans.t(1:256);
X = Data.AirScans.X(1:256,:);

% Determine sampling time
dt = diff(t);
T = mean(dt(:));

% Compute derivative matrix
[N, M] = size(X);
fun = @(theta, w) -1i*w;
D = tdtf(fun, 0, N, T);

% Initialize parameter structure
iFit = 1;
P = struct('var', [], 'mu', [], 'A', [], 'eta', [], 'ts', []);

%% Fit for delay
% Assume constant noise, average signal, and constant amplitude

Fix = struct('logv', true, 'mu', true, 'A', true, 'eta', false);
Ignore = struct('A', true, 'eta', false);
v0 = mean(var(X, 1, 2))*[1;eps;eps];
mu0 = mean(X, 2);
Options = ...
    struct('v0', v0, 'mu0', mu0, 'ts', T, 'Fix', Fix, 'Ignore', Ignore);

P(iFit) = tdnoisefit(X, Options);
eta0 = P(iFit).eta;
iFit = iFit + 1;
    
%% Fit for amplitude
% Assume constant noise, average signal, and delays from previous fit

Fix = struct('logv', true, 'mu', true, 'A', false, 'eta', true);
Ignore = struct('A', false, 'eta', true);
Options = struct('v0', v0, 'mu0', mu0, 'eta0', eta0, 'ts', T, ...
    'Fix', Fix, 'Ignore', Ignore);

P(iFit) = tdnoisefit(X, Options);
A0 = P(iFit).A;
iFit = iFit + 1;

%% Revise mu0

Xadjusted = airscancorrect(X, P(end));
mu0 = mean(Xadjusted, 2);

%% Fit for var
% Assume constant signal, amplitude, and delays from previous fits

Fix = struct('logv', false, 'mu', true, 'A', true, 'eta', true);
Ignore = struct('A', false, 'eta', false);
Options = struct('v0', v0, 'mu0', mu0, 'A0', A0, 'eta0', eta0, 'ts', T, ...
    'Fix', Fix, 'Ignore', Ignore);

P(iFit) = tdnoisefit(X, Options);
v0 = P(iFit).var;
iFit = iFit + 1;

%% Fit for all parameters

Fix = struct('logv', false, 'mu', false, 'A', false, 'eta', false);
Ignore = struct('A', false, 'eta', false);
Options = struct('v0', v0, 'mu0', mu0, ...
    'A0', A0, 'eta0', eta0, 'ts', T, 'Fix', Fix, 'Ignore', Ignore);

[P(iFit), nllmin, Diagnostic] = tdnoisefit(X, Options);

%% Compare model to measurements

vEst = P(end).var;
muEst = P(end).mu;
AEst = P(end).A;
etaEst = P(end).eta;

vErr = Diagnostic.Err.var;
muErr = Diagnostic.Err.mu;
AErr = Diagnostic.Err.A;
etaErr = Diagnostic.Err.eta;

zeta = zeros(N, M);
S = zeros(N, N, M);
for m = 1:M
    S(:,:,m) = shiftmtx(P(end).eta(m), N, T);
    zeta(:,m) = P(end).A(m)*S(:,:,m)*P(end).mu;
end

Dmu = D*P(end).mu;
valpha = P(end).var(1);
vbeta = P(end).var(2)*P(end).mu.^2;
vtau = P(end).var(3)*(Dmu).^2;
vtot = valpha + vbeta + vtau;

delta = (X - zeta)./sqrt(vtot);

valphalow = P(end).var(1) - vErr(1);
vbetalow = (P(end).var(2) - vErr(2))*P(end).mu.^2;
vtaulow = (P(end).var(3) - vErr(3))*(Dmu).^2;
vtotlow = valphalow + vbetalow + vtaulow;

valphahigh = P(end).var(1) + vErr(1);
vbetahigh = (P(end).var(2) + vErr(2))*P(end).mu.^2;
vtauhigh = (P(end).var(3) + vErr(3))*(Dmu).^2;
vtothigh = valphahigh + vbetahigh + vtauhigh;

Xadjusted = airscancorrect(X, P(end));

Calc.P =  P;
Calc.nllmin = nllmin;
Calc.vEst = vEst;
Calc.vErr = vErr;
Calc.muEst = muEst;
Calc.muErr = muErr;
Calc.AEst = AEst;
Calc.AErr = AErr;
Calc.etaEst = etaEst;
Calc.etaErr = etaErr;
Calc.t = t;
Calc.X = X;
Calc.zeta = zeta;
Calc.delta = delta;
Calc.vtot = vtot;
Calc.vtotlow = vtotlow;
Calc.vtothigh = vtothigh;
Calc.Xadjusted = Xadjusted;

end