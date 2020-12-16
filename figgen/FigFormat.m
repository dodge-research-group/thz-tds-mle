classdef (InferiorClasses = {?matlab.graphics.figure.Figure}) FigFormat
%FigFormat Class to format figures

% Optics Express Style Guide:
% https://www.osapublishing.org/submit/style/oestyleguide.cfm
% The online guide refers to the style file for layout details. The maximum
% figure width is not specified explicitly, but Fig. 2 of Rehn et al., Opt.
% Express 25, 6712 (2017) is 4.64 in x 2 in, and uses 9 pt Times font for
% all of the labels. I changed Times to Helvetica.
% Printing area: 5.25 in x 8.5 in
% Justified caption width: 4.5 in
% 

properties
    Scale
    Figure
    Axes
    Line
    Stem
    Legend
    Text
end
properties (Constant,Access=private)
    DefaultFigureWidth = 4.5;
    DefaultFigureHeight = 2;
    DefaultViewScale = 2;
end
methods
    function obj = FigFormat(varargin)
        switch nargin
            case 0
                scaleWidth = 1;
                scaleHeight = 1;
                scaleView = FigFormat.DefaultViewScale;
            case 1
                scaleWidth = varargin{1};
                if ~isnumeric(scaleWidth)
                     error('Value must be numeric')
                end
                scaleHeight = 1;
                scaleView = FigFormat.DefaultViewScale;
            case 2
                scaleWidth = varargin{1};
                scaleHeight = varargin{2};
                if ~isnumeric([scaleWidth scaleHeight])
                     error('Values must be numeric')
                end
                scaleView = FigFormat.DefaultViewScale;
            case 3
                scaleWidth = varargin{1};
                scaleHeight = varargin{2};
                scaleView = varargin{3};
                if ~isnumeric([scaleWidth scaleHeight scaleView])
                     error('Values must be numeric')
                end
        end
        obj.Scale = scaleView;
        posMonitor = get(groot,'MonitorPositions');
        
        ppi = get(groot,'ScreenPixelsPerInch');
        figWidth = FigFormat.DefaultFigureWidth*scaleWidth;
        figHeight = FigFormat.DefaultFigureHeight*scaleHeight;
        
        posW = floor(obj.Scale*figWidth*ppi);
        posH = floor(obj.Scale*figHeight*ppi);
        posL = floor((posMonitor(1,3) - posW)/2);
        posB = floor(0.65*(posMonitor(1,4) - posH)/2);
        
        obj.Figure = struct('Position',[posL, posB , posW, posH]);
        obj.Axes = struct('FontSize',9*obj.Scale, ...
            'LabelFontSizeMultiplier', 10/9, ...
            'FontName', 'Helvetica', ...
            'LineWidth', obj.Scale/4, ...
            'Box', 'on');
        obj.Line = struct('LineWidth', obj.Scale, ...
            'MarkerSize', 4*obj.Scale);
        obj.Stem = struct('LineWidth', obj.Scale/4, ...
            'MarkerSize', 2*obj.Scale);
        obj.Legend = struct('FontSize', 8*obj.Scale);
        obj.Text = struct('FontSize', 9*obj.Scale, ...
            'FontName', 'Helvetica');
    end
    function varargout = figure(varargin)
        obj = varargin{1};
        if nargin == 1
            Fig = figure('Position',obj.Figure.Position);
        elseif nargin > 1
            Fig = figure('Position',obj.Figure.Position,varargin{2:end});
        end
        Type = fieldnames(obj);
        for iType = 1:length(Type)
            if ismember(Type(iType),{'Axes','Line','Stem','Legend','Text'})
                Prop = fieldnames(obj.(Type{iType}));
                for jProp = 1:length(Prop)
                    set(Fig,...
                        ['default' Type{iType} Prop{jProp}],...
                        obj.(Type{iType}).(Prop{jProp}))
                end
            end
        end
        if nargout > 0
            varargout{1} = Fig;
        end
    end
end
    
end

