function varargout = plotnoiseest(varargin)
%PLOTNOISEEST plots the noise estimate

%% Parse inputs
Input = inputParser;
addOptional(Input, 'Calc', struct([]), ...
    @(x) any([isstruct(x),isempty(x)]));
addOptional(Input, 'axes', gca, @ishandle);
parse(Input, varargin{:});
Calc = Input.Results.Calc;
if isempty(Calc)
    Data = dataexp;
    
    % Restrict fit to first N time points
    N = 256;
    Data.AirScans.t(N+1:end) = [];
    Data.AirScans.X(N+1:end, :) = [];
    Calc = calcnoise(Data);
end
Ax = Input.Results.axes;

t = Calc.t;
vtot = Calc.vtot;
Xadjusted = Calc.Xadjusted;
M = size(Xadjusted,2);

%% Make plot
Ax.NextPlot = 'add';
plot(Ax, t, std(Xadjusted, 0, 2), 'k.')
plot(Ax, t, sqrt(vtot)*M/(M-1), '-', 'Color', [0.8 0 0])

circumflexChar = char(770);
sigmaChar = char(963);
muChar = char(956);

xlabel(Ax, 't (ps)');
ylabel(Ax, [sigmaChar circumflexChar '*_' muChar ' (pA)'])

ylim([0 11])

%% Return figure handle if requested
if nargout > 0
    varargout = {gcf};
end

end