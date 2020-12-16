function varargout = plotmleecdf(varargin)
%PLOTMLEECDF plots ECDF for Monte Carlo simulations of MLE fits

%% Parse inputs
Input = inputParser;
addOptional(Input, 'Calc', struct([]), ...
    @(x) any([isstruct(x), isempty(x)]));
addOptional(Input, 'axes', gca, @ishandle);
addParameter(Input, 'index', 0);

parse(Input, varargin{:});

Calc = Input.Results.Calc;
if isempty(Calc)
    Calc = calctdfit;
end
Ax = Input.Results.axes;
idx = Input.Results.index;

resnorm = Calc.resnorm;
nu = Calc.nu;

%% Make plot
[cf,x]=ecdf(resnorm);
x = [0;x;1000];
cf = [0;cf;1];

Ax.NextPlot = 'add';
plot(Ax, x, cf, 'k-');
plot(Ax, x, chi2cdf(x, nu), '-.', 'Color', [0.5 0.5 0.5]);

if idx > 0
    plot(Ax, x(idx), cf(idx), 'x', 'Color', [0.8 0 0]);
end

chiChar = char(967);
xlabel(Ax, 'GOF statistic');
ylabel(Ax, 'Cumulative probability')
legend(Ax, 'S_{TLS}', [chiChar '^2'], 'Location', 'NE', 'Box', 'off')

if nargout > 0
    varargout = {gcf};
end

end