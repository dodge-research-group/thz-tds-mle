function varargout = plotmlenormprob(varargin)
%PLOTMLENORMPROB plots normal probability plot for MLE fit residuals

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

residual = Calc.residual(:,idx);

%% Make plot
h = probplot(Ax, residual, 'noref');
h(1).Marker = '.';
h(1).Color = 'k';
Ax.Title.String = '';
Ax.XGrid = 'on';
Ax.YGrid = 'on';
Ax.NextPlot = 'add';
plot(Ax, [-3.5 3.5], [-3.5 3.5], '-.', 'Color', [0.8 0.8 0.8])
Ax.XLim = [-3.5 3.5];
xlabel('Normed residual');
axis square

Ch = Ax.Children;
Ax.Children = flipud(Ch);

if nargout > 0
    varargout = {gcf};
end

end