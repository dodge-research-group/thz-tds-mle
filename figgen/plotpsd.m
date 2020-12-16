function varargout = plotpsd(varargin)
%PLOTPSD plots Monte Carlo demo data in the frequency domain

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
Ym = Data.Ym(1:Nf,:);

%% Make plot
Ax.NextPlot = 'add';
plot(Ax, f, 10*log10(abs(Ym(:,2:end)).^2/max(abs(Ym(:).^2))),...
    'Color', [0.5 0.5 0.5]);
plot(Ax, f, 10*log10(abs(Ym(:,1)).^2/max(abs(Ym(:).^2))), ...
    '-', 'Color', [0.8 0 0]);

xlabel(Ax, 'Frequency (THz)');
ylabel(Ax, 'Relative Power (dB)')
ylim([-70 10])

if nargout > 0
    varargout = {gcf};
end

end