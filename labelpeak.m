% This function adds a vertical line and label a peak on a power spectrum.
% IN:
% f: Frequency to label
% yloc: y-coordinate of label
% 'align': (default 'right') horizontal alignment of text (left or right)
% 'stdf': (default 0) uncertainty in frequency, 1 standard deviation 
% 'fontsize': (default 9)
%
% OUT:
% h: handle of line
%
% TO DO:
% - allow user to specify precision of labels
%
% Adrian Tasistro-Hart, adrianraph-at-gmail.com, 18.08.2018

function h = labelpeak(f,yloc,varargin)

% parse inputs
parser = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(parser,'f',@isscalar)
addRequired(parser,'yloc',@isscalar)
addParameter(parser,'align','right',@ischar)
addParameter(parser,'stdf',0,validScalarPosNum)
addParameter(parser,'fontsize',9,validScalarPosNum)

parse(parser,f,yloc,varargin{:})

f = parser.Results.f;
yloc = parser.Results.yloc;
aln = parser.Results.align;
stdf = parser.Results.stdf;
fontsize = parser.Results.fontsize;

% validate alignment
aln = validatestring(aln,{'left','right'});

% formulate text
if stdf == 0
    txt = sprintf('%1.0f',1/f);
else
    txt = sprintf('%1.0f^{+%1.1f}_{-%1.1f}',...
        1/f, abs(1/f - 1/(f+2*stdf)), abs(1/f - 1/(f-2*stdf)));
end

h = vline(f,'k-');
text(f,yloc,txt,'fontsize',fontsize,...
    'HorizontalAlignment',aln,'VerticalAlignment','bottom',...
    'BackgroundColor','none','Margin',1)

end