function varargout = plotetfre(varargin)
%PLOTETFRE plots the real part of the ETF

%% Parse inputs
Input = inputParser;
addOptional(Input, 'Data', datamc, ...
    @(x) any([isstruct(x), isempty(x)]));
addOptional(Input, 'axes', gca, @ishandle);
parse(Input, varargin{:});
Data = Input.Results.Data;
Ax = Input.Results.axes;

Nf = Data.Nf;
f = Data.f(1:Nf);
YmRatio = Data.YmRatio(1:Nf,:);

%% Make plot
HhatChar = char(292);
Ax.NextPlot = 'add';
plot(Ax, f, real(YmRatio), '.', 'Color', [0.8 0.8 0.8])
plot(Ax, f, real(YmRatio(:,1)), '.', 'Color', [0.8 0 0])
plot(Ax, f, real(YmRatio(:,1)), '-', 'Color', [0.8 0 0], 'LineWidth', 0.5)
plot(Ax, f, mean(real(YmRatio), 2), 'k-')

xlabel(Ax, 'Frequency (THz)');
ylabel(Ax, ['Re\{' HhatChar '_{ETFE}\}'])
ylim([-1 2])
if nargout > 0
    varargout = {gcf};
end

end