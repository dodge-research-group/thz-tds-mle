function varargout = plotsignal(varargin)
%PLOTSIGNAL plots Monte Carlo demo data in the time domain

%% Parse inputs
Input = inputParser;
addOptional(Input, 'Data', datamc, ...
    @(x) any([isstruct(x), isempty(x)]));
addOptional(Input, 'axes', gca, @ishandle);
parse(Input, varargin{:});
Data = Input.Results.Data;
Ax = Input.Results.axes;

t = Data.t;
y0 = Data.y0;

sigmaAlpha = Data.P.sigmaAlpha;
sigmaBeta = Data.P.sigmaBeta;
sigmaTau = Data.P.sigmaTau;
sigmaVec = [sigmaAlpha, sigmaBeta, sigmaTau];
T = Data.P.T;
sigmaT = noiseamp(sigmaVec, y0, T);

%% Make plot
sigmaChar = char(963);
muChar = char(956);
Ax.NextPlot = 'add';
PltNoise = plot(Ax, t, 30*sigmaT, '-.', 'Color', [0.5 0.5 0.5],...
    'DisplayName', ['Noise (30' sigmaChar ')']);
PltSignal = plot(Ax, t, y0, 'k', ...
    'DisplayName', ['Signal (' muChar ')']);

xlabel(Ax,'Time (ps)');
ylabel(Ax,['Amplitude (units of ' muChar '_p)'])
Ax.XLim = [0 12.8];
Ax.YLim = [-0.5 1.1];

legend([PltSignal, PltNoise], 'Location', 'NE', 'Box', 'off')

if nargout > 0
    varargout = {gcf};
end

end