function varargout = figmlefit(varargin)
%FIGMLEFIT plots ECDF and representative residuals for MLE fits
%
%   Fig = FIGMLEFIT returns figure handle

%% Parse inputs
Input = inputParser;

addOptional(Input, 'Calc', struct([]), ...
    @(x) any([isstruct(x), isempty(x)]));

parse(Input, varargin{:});

Calc = Input.Results.Calc;
if isempty(Calc)
    Calc = calctdfit;
end

%% Make figure and load default format specifications

Fig = figure(FigFormat);

panelLabelX = 0.05;
panelLabelY = 0.95;

insetScale = 0.4;
insetMarginRight = 0.0375;
insetMarginBottom = 0.2;

idx = round(length(Calc.resnorm)/2)+1;

% ECDF
Ax(1) = subplot(1, 2, 1);
plotmleecdf(Calc, Ax(1), 'index', idx);
Ax.XLim = [0 1000];
text(panelLabelX, panelLabelY, '(a)', 'Units', 'Normalized');

% Probability plot (inset)
pos = Ax(1).Position;
posInset = [pos(1) + (1 - insetScale)*pos(3) - insetMarginRight, ...
    pos(2) + insetMarginBottom, ...
    insetScale*pos(3), insetScale*pos(4)];
Inset = axes('Position',posInset);
plotmlenormprob(Calc, Inset, 'index', idx);

pTick = [0.01 0.1 0.5 0.9 0.99];
Inset.YLim = norminv([0.001 0.999]);
Inset.YTick = norminv(pTick);
Inset.YTickLabel = compose('%g', pTick);


% Residuals
Ax(2) = subplot(1, 2, 2);

plotresmle(Calc, Ax(2), 'index', idx);
xlim(Ax(2), [0 12.8])
ylim(Ax(2), [-3.5 3.5])
set(Ax(2), 'YTick', (-3:3))
text(panelLabelX, panelLabelY, '(b)', 'Units', 'Normalized');


if nargout > 0
    varargout{1} = Fig;
end

end