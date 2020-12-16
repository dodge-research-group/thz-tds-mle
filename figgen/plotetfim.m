function varargout = plotetfim(varargin)
%PLOTETFIM plots the imaginary part of the ETF

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
hatChar = char(292);
Ax.NextPlot = 'add';
plot(Ax, f, imag(YmRatio), '.', 'Color', [0.8 0.8 0.8])
plot(Ax, f, imag(YmRatio(:,1)), '.', 'Color', [0.8 0 0])
plot(Ax, f, imag(YmRatio(:,1)), '-', 'Color', [0.8 0 0], 'LineWidth', 0.5)
plot(Ax, f, mean(imag(YmRatio), 2), 'k-')

xlabel(Ax, 'Frequency (THz)');
ylabel(Ax, ['Im\{' hatChar '_{ETFE}\}'])
ylim([-2 2])
if nargout > 0
    varargout = {gcf};
end

end