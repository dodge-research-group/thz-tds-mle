function Calc = calctfexp(varargin)
%CALCTFEXP estimate transfer function of experiment
%

%% Parse inputs
Input = inputParser;

addOptional(Input, 'Data', struct([]), ...
    @(x) any([isstruct(x), isempty(x)]));
addOptional(Input, 'Calc', struct([]), ...
    @(x) any([isstruct(x), isempty(x)]));
addOptional(Input, 'Nmax', 256, @(x) isnumeric(x));

parse(Input, varargin{:});

Data = Input.Results.Data;
if isempty(Data)
    Data = dataexp;
end

Calc = Input.Results.Calc;
if isempty(Calc)
    Calc = calcnoise(Data);
end

Nmax = Input.Results.Nmax;

t = Calc.t(1:Nmax);
X = Calc.X(1:Nmax,:);
[N, M] = size(X);
Npair = M-1;

T = mean(diff(t));
f = fftfreq(N, T);
w = 2*pi*f;

Nf = floor(N/2)+1;

Xf = conj(fft(X));  % Use -i*w*t convention for FFT
ETFE = Xf(:,2:end)./Xf(:,1:end-1);
Vr = var(real(ETFE), 0, 2);
Vi = var(imag(ETFE), 0, 2);
V = (Vr + Vi);

%%
% Transfer function definition and parameters

tfun = @(theta,w) theta(1)*exp(1i*theta(2)*w);
normresfunreim = @(H, theta) ...
    [real(H - tfun(theta, w))./sqrt(V);...
    imag(H - tfun(theta, w))./sqrt(V)];

tfunext = @(theta,w) (theta(1) + theta(3)*w.^2)...
    .*exp(1i*(theta(2) + theta(4)*w.^2).*w);
normresfunextreim = @(H, theta) ...
    [real(H - tfunext(theta, w))./sqrt(V);...
    imag(H - tfunext(theta, w))./sqrt(V)];

A0 = 1;          % amplitude ratio between pulses
eta0 = 0;           % delay between pulses [T]
theta0 = [A0;eta0]; % Initial parameter vector
Np = length(theta0);

theta0ext = [A0;eta0;0;0];
Npext = length(theta0ext);

%% Run FD fits

Fit.x0 = theta0;
Fit.lb = [];
Fit.ub = [];
Fit.solver = 'lsqnonlin';
Fit.options = optimoptions('lsqnonlin',...
    'Display','off',...
    'UseParallel',true);

pFitFD = zeros(Np,Npair);
pFitFDErr = zeros(Np,Npair);
resnormFD = zeros(1,Npair);
residualFD = zeros(2*N,Npair);
cvFD = zeros(Np,Np,Npair);

pFitFDext = zeros(Npext,Npair);
pFitFDextErr = zeros(Npext,Npair);
resnormFDext = zeros(1,Npair);
residualFDext = zeros(2*N,Npair);
cvFDext = zeros(Npext,Npext,Npair);

for j = 1:Npair
    Fit.objective = @(theta) normresfunreim(ETFE(:,j), theta);
    [pFitFD(:,j), resnormFD(j), residualFD(:,j), ~, ~, ~, jac] ...
        = lsqnonlin(Fit);
    [~, R] = qr(jac);
    Rinv = R\eye(size(R));
    cvFD(:,:,j) = (Rinv*Rinv');
    pFitFDErr(:,j) = sqrt(diag(squeeze(cvFD(:,:,j))));
end

Fit.x0 = theta0ext;
for j = 1:Npair
    Fit.objective = @(theta) normresfunextreim(ETFE(:,j), theta);
    [pFitFDext(:,j), resnormFDext(j), residualFDext(:,j), ~, ~, ~, jac] ...
        = lsqnonlin(Fit);
    [~, R] = qr(jac);
    Rinv = R\eye(size(R));
    cvFDext(:,:,j) = (Rinv*Rinv');
    pFitFDextErr(:,j) = sqrt(diag(squeeze(cvFDext(:,:,j))));
end

%% Run TD fits

ym1 = X(:,1:end-1);
ym2 = X(:,2:end);
vEst = Calc.vEst;
sigmaVec = sqrt(vEst);

pFitTD = zeros(Np,Npair);
pFitTDErr = zeros(Np,Npair);
resnormTD = zeros(1,Npair);
residualTD = zeros(N,Npair);
cvTD = zeros(Np,Np,Npair);

pFitTDext = zeros(Npext,Npair);
pFitTDextErr = zeros(Npext,Npair);
resnormTDext = zeros(1,Npair);
residualTDext = zeros(N,Npair);
cvTDext = zeros(Npext,Npext,Npair);

Fit.x0 = theta0;
for j=1:Npair
    sigmay1 = noiseamp(sigmaVec, ym1(:,j), T);
    sigmay2 = noiseamp(sigmaVec, ym2(:,j), T);
    Fit.objective = @(theta) ...
        costfunlsq(tfun,theta,ym1(:,j),ym2(:,j),...
        sigmay1,sigmay2,w);
    [p, resnormTD(j), residualTD(:, j),~,~,~,jac] = lsqnonlin(Fit);
    pFitTD(:, j) = p(:);
    [~, R] = qr(jac);
    Rinv = R\eye(size(R));
    cvTD(:,:,j) = (Rinv*Rinv');
    pFitTDErr(:,j) = sqrt(diag(squeeze(cvTD(:,:,j))));
end

Fit.x0 = theta0ext;
for j=1:Npair
    sigmay1 = noiseamp(sigmaVec, ym1(:,j), T);
    sigmay2 = noiseamp(sigmaVec, ym2(:,j), T);
    Fit.objective = @(theta) ...
        costfunlsq(tfunext,theta,ym1(:,j),ym2(:,j),...
        sigmay1,sigmay2,w);
    [p, resnormTDext(j), residualTDext(:, j),~,~,~,jac] ...
        = lsqnonlin(Fit);
    pFitTDext(:, j) = p(:);
    [~, R] = qr(jac);
    Rinv = R\eye(size(R));
    cvTDext(:,:,j) = (Rinv*Rinv');
    pFitTDextErr(:,j) = sqrt(diag(squeeze(cvTDext(:,:,j))));
end

%% Transfer variables to output structure
Calc.t = t;
Calc.X = X;
Calc.T = T;
Calc.f = f;
Calc.Xf = Xf;
Calc.ETFE = ETFE;
Calc.V = V;
Calc.N = N;
Calc.Nf = Nf;
Calc.M = M;
Calc.Np = Np;
Calc.nu = N - Np;
Calc.Npext = Npext;
Calc.nuext = N - Npext;

Calc.resnormFD = resnormFD;
Calc.AICFD = resnormFD + 2*Np;
Calc.pFitFD = pFitFD;
Calc.pFitFDErr = pFitFDErr;
Calc.residualFD = residualFD;
Calc.cvFD = cvFD;

Calc.resnormTD = resnormTD;
Calc.AICTD = resnormTD + 2*Np;
Calc.pFitTD = pFitTD;
Calc.pFitTDErr = pFitTDErr;
Calc.residualTD = residualTD;
Calc.cvTD = cvTD;

Calc.resnormFDext = resnormFDext;
Calc.AICFDext = resnormFDext + 2*Npext;
Calc.pFitFDext = pFitFDext;
Calc.pFitFDextErr = pFitFDextErr;
Calc.residualFDext = residualFDext;
Calc.cvFDext = cvFDext;

Calc.resnormTDext = resnormTDext;
Calc.AICTDext = resnormTDext + 2*Npext;
Calc.pFitTDext = pFitTDext;
Calc.pFitTDextErr = pFitTDextErr;
Calc.residualTDext = residualTDext;
Calc.cvTDext = cvTDext;

end