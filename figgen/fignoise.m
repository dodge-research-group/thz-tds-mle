function varargout = fignoise(varargin)
%FIGNOISE plots data, amplitude, and residuals for noise model
%
%   Fig = FIGNOISE returns figure handle

%% Parse inputs
Input = inputParser;

addOptional(Input, 'Data', struct([]), ...
    @(x) any([isstruct(x), isempty(x)]));
addOptional(Input, 'Calc', struct([]), ...
    @(x) any([isstruct(x), isempty(x)]));

parse(Input, varargin{:});

Data = Input.Results.Data;
if isempty(Data)
    Data = dataexp;
end

Calc = Input.Results.Calc;
if isempty(Calc)
    Calc = calcnoise(Data);
end

%% Make figure and load default format specifications

% Restrict fit to first N time points
N = 256;
Data.AirScans.t(N+1:end) = [];
Data.AirScans.X(N+1:end, :) = [];
tMax = Data.AirScans.t(N);

Fig = figure(FigFormat(1, 1.5));

% Rescale height for multiple panels
pos = get(Fig, 'defaultAxesPosition');
curpos = pos;
curpos(4) = pos(4)/3;

panelLabelX = 0.025;
panelLabelY = 0.9;

% Data
curpos(2) = pos(2) + 2*curpos(4);
Ax(1) = axes('Parent', Fig,...
    'Position', curpos,...
    'Layer', 'top');

plotnoisedata(Data, Ax(1));
Ax(1).XLim = [0 tMax];
Ax(1).XTick = (0:2:12);
Ax(1).XLabel = [];
Ax(1).XTickLabel = {};
Ax(1).YLim = [-250 700];
text(panelLabelX, panelLabelY, '(a)', 'Units', 'Normalized');

% Noise amplitude
curpos(2) = pos(2) + curpos(4);
Ax(2) = axes('Parent', Fig,...
    'Position', curpos,...
    'Layer', 'top');

plotnoiseest(Calc, Ax(2));
Ax(2).XLim = [0 tMax];
Ax(2).XTick = (0:2:12);
Ax(2).XLabel = [];
Ax(2).XTickLabel = {};
Ax(2).YLim = [0 11];
text(panelLabelX, panelLabelY, '(b)', 'Units', 'Normalized');

% Residual
curpos(2) = pos(2);
Ax(3) = axes('Parent', Fig,...
    'Position', curpos,...
    'Layer', 'top');

plotnoiseres(Calc, Ax(3));
Ax(3).XLim = [0 tMax];
Ax(3).XTick = (0:2:12);
Ax(3).YLim = [-3.5 3.5];
text(panelLabelX, panelLabelY, '(c)', 'Units', 'Normalized');


if nargout > 0
    varargout{1} = Fig;
end

end