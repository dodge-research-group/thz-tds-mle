function varargout = figsignalnoise(varargin)
%FIGSIGNALNOISE plots simulated signal and noise
%
%   Fig = FIGSIGNALNOISE returns figure handle

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

Data.Ym = Data.Ym(:,1:10);
Fig = figure(FigFormat);

panelLabelX = 0.05;
panelLabelY = 0.95;

% Time traces
Ax(1) = subplot(1, 2, 1);
plotsignal(Data, Ax(1));
text(panelLabelX, panelLabelY, '(a)', 'Units', 'Normalized');

% Power spectrum
Ax(2) = subplot(1, 2, 2);

plotpsd(Data, Ax(2));
text(panelLabelX, panelLabelY, '(b)', 'Units', 'Normalized');


if nargout > 0
    varargout{1} = Fig;
end

end