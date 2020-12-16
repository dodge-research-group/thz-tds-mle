function varargout = figtfexp(CalcTF, idx)
%FIGTFEXP plots experimental transfer function fits
%
%   Fig = FIGTFEXP(CalcTF, idx) returns figure handle

% Make figure and load default format specifications
% Data = dataexp;
% CalcNoise = calcnoise(Data);
% CalcTF = calctfexp(CalcNoise);
t = CalcTF.t;
tMax = t(end);
X = CalcTF.X;

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

% Waveforms
curpos(1) = 1.5*pos(1)/wRescale;
curpos(2) = (pos(2)+pos(4))/hRescale;
curpos(3) = (pos(3) - 0.5*pos(1))/wRescale;
curpos(4) = pos(4)/hRescale;
Ax(1) = axes('Parent', Fig,...
    'Position', curpos,...
    'Layer', 'top');
plot(t, X(:,[idx , idx + 1]), 'Color', [0.2 0.2 0.2])
ylabel("x (pA)")
Ax(1).XLim = [0 tMax];
Ax(1).XTick = (0:2:tMax);
Ax(1).XLabel = [];
Ax(1).XTickLabel = {};
Ax(1).YLim = [-250 700];
text(panelLabelX, panelLabelYa, '(a)', 'Units', 'Normalized');

% Time-domain residuals
curpos(2)=pos(2)/hRescale;
Ax(2) = axes('Parent', Fig,...
    'Position', curpos,...
    'Layer', 'top');

stem(t, CalcTF.residualTD(:,idx), '.', 'Color', [0.4 0.4 0.4])
xlabel('Time (ps)')
ylabel(["Normalized"; "residual"])
Ax(2).Layer = 'top';
Ax(2).XLim = [0 tMax];
Ax(2).XTick = (0:2:tMax);
Ax(2).YLim = [-3.5 3.5];
text(panelLabelX, panelLabelYb, '(b)', 'Units', 'Normalized');

% ETFE
f = CalcTF.f;
Nf = CalcTF.Nf;
f = f(1:Nf);
w = 2*pi*f;
ETFE = CalcTF.ETFE(1:Nf,idx);
pFitFD = CalcTF.pFitFD(:,idx);
pFitFDext = CalcTF.pFitFDext(:,idx);
% pFitTD = CalcTF.pFitTD(:,i);
ETFEfitFD = pFitFD(1)*exp(1i*w*pFitFD(2));
ETFEfitFDext = (pFitFDext(1) + pFitFDext(3)*w.^2)...
    .*exp(1i*w.*(pFitFDext(2) + pFitFDext(4)*w.^2));
curpos(1) = (3*pos(1) + pos(3))/wRescale;
curpos(2) = (pos(2) + pos(4))/hRescale;
Ax(3) = axes('Parent', Fig,...
    'Position', curpos,...
    'YAxisLocation', 'left',...
    'Layer', 'top');

HhatChar = char(292);
plot(f, real(ETFE)-1, 'o', 'LineWidth', 0.5, 'MarkerSize', 4);
hold on
plot(f, imag(ETFE), 'x', 'LineWidth', 0.5, 'MarkerSize', 4);
% set(gca, 'ColorOrderIndex', 1);
% plot(f, [real(ETFE)-1, imag(ETFE)], '-', 'LineWidth', 0.25);
set(gca, 'ColorOrderIndex', 1);
plot(f, [real(ETFEfitFD)-1, imag(ETFEfitFD)], '-', 'LineWidth', 1.0);
set(gca, 'ColorOrderIndex', 1);
plot(f, [real(ETFEfitFDext)-1, imag(ETFEfitFDext)], ':', 'LineWidth', 1.0);
ylabel([HhatChar '_{ETFE} - 1'])
Ax(3).XLim = [0 4];
Ax(3).XLabel = [];
Ax(3).XTickLabel = {};
Ax(3).YLim = [-0.25 0.25];
Ax(3).YTick = (-0.2:0.1:0.2);
legend('Re', 'Im', 'Location', 'SW')
text(panelLabelX, panelLabelYa, '(c)', 'Units', 'Normalized');

% Frequency-domain residuals
curpos(2)=pos(2)/hRescale;
Ax(4) = axes('Parent', Fig,...
    'Position', curpos,...
    'YAxisLocation', 'left',...
    'Layer', 'top');

residualFD = CalcTF.residualFD(:,idx);
N = CalcTF.N;
stem(f, residualFD(1:Nf), 'o')
hold on
stem(f, residualFD(N + (1:Nf)), 'x')
xlabel('Frequency (THz)')
ylabel(["Normalized"; "residual"])
text(panelLabelX, panelLabelYb, '(d)', 'Units', 'Normalized');
Ax(4).XLim = [0 4];
Ax(4).YLim = [-1.75 1.75];
Ax(4).YTick = (-1:1);
legend('Re', 'Im', 'Location', 'SW')

if nargout > 0
    varargout{1} = Fig;
end

end