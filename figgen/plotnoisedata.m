function varargout = plotnoisedata(varargin)
%PLOTNOISEDATA plots the raw data for the noise estimate

%% Parse inputs
Input = inputParser;
addOptional(Input, 'Data', struct([]), ...
    @(x) any([isstruct(x), isempty(x)]));
addOptional(Input, 'axes', gca, @ishandle);

parse(Input, varargin{:});

Data = Input.Results.Data;
if isempty(Data)
    Data = datamc;
end
Ax = Input.Results.axes;

t = Data.AirScans.t;
X = Data.AirScans.X;

%% Make plot
plot(Ax, t, X(:,1:end-1), 'Color', [0.2 0.2 0.2])
Ax.NextPlot = 'add';
plot(Ax, t, X(:,end), 'LineWidth', 1, 'Color', [0.8 0 0])

xlabel(Ax, 't (ps)');
ylabel(Ax, 'x (pA)')

ylim([-250 700])

%% Return figure handle if requested
if nargout > 0
    varargout = {gcf};
end

end