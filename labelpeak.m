% This function adds a vertical line and labels a peak on a power spectrum
% in terms of that peak's period.
%
% IN:
% f: Frequency to label
% yloc: y-coordinate of label
% 'align': (default 'right') horizontal alignment of text ('left', 'center'
%   or 'right')
% 'stdf': (default 0) uncertainty in frequency, 1 standard deviation 
% 'fontsize': (default 9)
% 'plotline': (default false) plot a vertical line intersecting the peak
%
% OUT:
% htxt: handle to the label text
% hln: handle of line, if line is requested
%
% TO DO:
% - allow user to specify precision of labels
% - allow user to 
%
% Adrian Tasistro-Hart, adrianraph-at-gmail.com, 18.08.2018

function [htxt,hln] = labelpeak(f,yloc,varargin)

% parse inputs
parser = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(parser,'f',@isscalar)
addRequired(parser,'yloc',@isscalar)
addParameter(parser,'align','right',@ischar)
addParameter(parser,'stdf',0,validScalarPosNum)
addParameter(parser,'fontsize',9,validScalarPosNum)
addParameter(parser,'plotline',0,@(x) x==1 || x==0)

parse(parser,f,yloc,varargin{:})

f = parser.Results.f;
yloc = parser.Results.yloc;
aln = parser.Results.align;
stdf = parser.Results.stdf;
fontsize = parser.Results.fontsize;
plotline = parser.Results.plotline;

% validate alignment
aln = validatestring(aln,{'left','right','center'});

% formulate text
if stdf == 0
    txt = sprintf('%1.0f',1/f);
else
    txt = sprintf('%1.0f^{+%1.1f}_{-%1.1f}',...
        1/f, abs(1/f - 1/(f-2*stdf)), abs(1/f - 1/(f+2*stdf)));
end

if plotline
    hln = vline(f,'k-');
end

htxt = text(f,yloc,txt,'fontsize',fontsize,...
    'HorizontalAlignment',aln,'VerticalAlignment','bottom',...
    'BackgroundColor','none','Margin',1);

end