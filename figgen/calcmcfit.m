function Calc = calcmcfit(varargin)
%CALCMCFIT fits to Monte Carlo simulation data
%

%% Parse inputs
Input = inputParser;

addOptional(Input, 'Data', struct([]), ...
    @(x) any([isstruct(x), isempty(x)]));

parse(Input, varargin{:});

Data = Input.Results.Data;
if isempty(Data)
    Data = datamc;
end

%% Run fits

f = Data.f;
N = Data.P.N;
YmRatio = Data.YmRatio;

V = Data.V;

Nmc = Data.P.Nmc;

Fit.x0 = [1, 0];
Fit.lb = [];
Fit.ub = [];
Fit.solver = 'lsqnonlin';
Fit.options = optimoptions('lsqnonlin',...
    'Display', 'off');

w = 2*pi*f;
normresfunreim = @(H, A, tau) ...
    [real(H - A*exp(1i*w*tau))./sqrt(V);...
    imag(H - A*exp(1i*w*tau))./sqrt(V)];

NmcPair = Nmc/2;
pFit = zeros(2, NmcPair);
resnorm = zeros(NmcPair, 1);
residual = zeros(2*N, NmcPair);
cv = zeros(2, 2, NmcPair);

for j = 1:NmcPair
    Fit.objective = @(theta) ...
        normresfunreim(YmRatio(:,j), theta(1), theta(2));
    [pFit(:,j), resnorm(j), residual(:,j), ~, ~, ~, jac] = lsqnonlin(Fit);
    [~, R] = qr(jac);
    Rinv = R\eye(size(R));
    cv(:,:,j) = (Rinv*Rinv');
end

[resnorm, idxSort] = sort(resnorm);

Calc.f = f;
Calc.N = N;
Calc.resnorm = resnorm;
Calc.nu = N -2;
Calc.pFit = pFit(:,idxSort);
Calc.residual = residual(1:N,idxSort) + 1i*residual(N+1:end,idxSort);
Calc.cv = cv(:,:,idxSort);
Calc.idxSort = idxSort;

end