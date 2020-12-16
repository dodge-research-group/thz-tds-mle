function varargout = plotresmle(varargin)
%PLOTRESMLE plots MLE residuals for Monte Carlo fit

%% Parse inputs
Input = inputParser;
addOptional(Input, 'Calc', struct([]), ...
    @(x) any([isstruct(x), isempty(x)]));
addOptional(Input, 'axes', gca, @ishandle);
addParameter(Input, 'index', 1);

parse(Input, varargin{:});

Calc = Input.Results.Calc;
if isempty(Calc)
    Calc = calctdfit;
end
Ax = Input.Results.axes;
idx = Input.Results.index;

t = Calc.t;
residual = Calc.residual(:,idx);

%% Make plot
Ax.NextPlot = 'add';
stem(Ax, t, residual, '.', 'Color', [0.4 0.4 0.4])

xlabel(Ax, 'Time (ps)')
ylabel(Ax, 'Normed residual')

if nargout > 0
    varargout = {gcf};
end

end