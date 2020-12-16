function varargout = plotresetfe(varargin)
%PLOTECDF plots ETFE residuals for Monte Carlo fit

%% Parse inputs
Input = inputParser;
addOptional(Input, 'Calc', struct([]), ...
    @(x) any([isstruct(x), isempty(x)]));
addOptional(Input, 'axes', gca, @ishandle);
addParameter(Input, 'index', 1);

parse(Input, varargin{:});

Calc = Input.Results.Calc;
if isempty(Calc)
    Calc = calcmcfit;
end
Ax = Input.Results.axes;
idx = Input.Results.index;

f = Calc.f;
N = length(f);
Nf = floor(N/2)+1;
f = f(1:Nf);
residual = Calc.residual(:,idx);

%% Make plot
Ax.NextPlot = 'add';
stem(Ax, f, real(residual(1:Nf)), 'o')
stem(Ax, f(2:N-Nf+1), imag(residual(2:N-Nf+1)), 'x')

xlabel(Ax, 'Frequency (THz)')
ylabel(Ax, 'Normed residual')
legend(Ax, 'Real', 'Imag', 'Location', 'NE', 'Box', 'off')

if nargout > 0
    varargout = {gcf};
end

end