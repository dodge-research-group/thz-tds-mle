function varargout = figetf(varargin)
%FIGETF plots real and imaginary parts of the ETFE
%
%   Fig = FIGETF returns figure handle

%% Parse inputs
Input = inputParser;

addOptional(Input, 'Data', struct([]), ...
    @(x) any([isstruct(x), isempty(x)]));

parse(Input, varargin{:});

Data = Input.Results.Data;
if isempty(Data)
    Data = datamc;
end

%% Make figure and load default format specifications

Fig = figure(FigFormat(1, 1.5));

% Rescale height for multiple panels
pos = get(Fig, 'defaultAxesPosition');
pos(4) = pos(4) - pos(2);
pos(2) = 2*pos(2);
hRescale = 1. + pos(4);
wRescale = 1. + pos(3) + pos(1);

panelLabelX = 0.05;
panelLabelYa = 0.95;
panelLabelYb = 0.95;

% Real part
curpos(1) = 1.5*pos(1)/wRescale;
curpos(2) = (pos(2)+pos(4))/hRescale;
curpos(3) = (pos(3) - 0.5*pos(1))/wRescale;
curpos(4) = pos(4)/hRescale;
Ax(1) = axes('Parent', Fig,...
    'Position', curpos,...
    'Layer', 'top');
plotetfre(Data, Ax(1));
Ax(1).XLim = [0 10];
Ax(1).XTick = (0:2:10);
Ax(1).XLabel = [];
Ax(1).XTickLabel = {};
Ax(1).YLim = [-1.75 2.75];
text(panelLabelX, panelLabelYa, '(a)', 'Units', 'Normalized');

% Imaginary part
curpos(2)=pos(2)/hRescale;
Ax(2) = axes('Parent', Fig,...
    'Position', curpos,...
    'Layer', 'top');

plotetfim(Data, Ax(2));
Ax(2).Layer = 'top';
Ax(2).XLim = [0 10];
Ax(2).XTick = (0:2:10);
Ax(2).YLim = [-2.25 2.25];
text(panelLabelX, panelLabelYb, '(b)', 'Units', 'Normalized');

% Real part, zoomed
curpos(1) = (3*pos(1) + pos(3))/wRescale;
curpos(2) = (pos(2) + pos(4))/hRescale;
Ax(3) = axes('Parent', Fig,...
    'Position', curpos,...
    'YAxisLocation', 'left',...
    'Layer', 'top');
plotetfre(Data, Ax(3));
Ax(3).XLim = [0 4];
Ax(3).XLabel = [];
Ax(3).XTickLabel = {};
Ax(3).YLim = [0.88 1.12];
Ax(3).YTick = (0.9:0.1:1.1);
text(panelLabelX, panelLabelYa, '(c)', 'Units', 'Normalized');

% Imaginary part, zoomed
curpos(2)=pos(2)/hRescale;
Ax(4) = axes('Parent', Fig,...
    'Position', curpos,...
    'YAxisLocation', 'left',...
    'Layer', 'top');

plotetfim(Data, Ax(4));
text(panelLabelX, panelLabelYb, '(d)', 'Units', 'Normalized');
Ax(4).XLim = [0 4];
Ax(4).YLim = [-0.12 0.12];
Ax(4).YTick = (-0.1:0.1:0.1);

if nargout > 0
    varargout{1} = Fig;
end

end