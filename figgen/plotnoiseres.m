function varargout = plotnoiseres(varargin)
%PLOTNOISERES plots the noise estimate residuals

%% Parse inputs
Input = inputParser;
addOptional(Input, 'Calc', struct([]), ...
    @(x) any([isstruct(x), isempty(x)]));
addOptional(Input, 'axes', gca, @ishandle);
parse(Input, varargin{:});
Calc = Input.Results.Calc;
Ax = Input.Results.axes;

if isempty(Calc)
    Calc = calcnoise;
end

t = Calc.t;
delta = Calc.delta;
M = size(Calc.X, 2);

%% Make plot
stem(Ax, t, delta(:,6)/sqrt(M/(M-1)), '.', 'Color', [0.4 0.4 0.4])

xlabel(Ax, 't (ps)');
ylabel(Ax, {'Normalized';'residual'})

ylim([-3 3])

%% Return figure handle if requested
if nargout > 0
    varargout = {gcf};
end

end