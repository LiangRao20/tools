fprintf('\n startup.m entered \n');

% update toolboxcache
rehash toolboxcache

% unicode support
feature('DefaultCharacterSet', 'UTF8');

warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
warning('off','MATLAB:lang:cannotClearExecutingFunction');

%% Make good figures

% good colormap
set(groot,'DefaultFigureColormap',flipud(cbrewer('div', 'RdYlBu', 32)));
%set(groot, 'DefaultFigureColormap', parula); close;
fontName = 'Fira Sans';

% figure properties
%set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultFigureColor','w')
set(groot,'DefaultFigureRenderer','opengl')
set(groot,'DefaultFigurePaperPositionMode', 'auto');

set(groot,'DefaultTextFontName', fontName);
set(groot,'DefaultTextColor','k')
set(groot,'DefaultTextFontSize',24);

% line style and color cycling
set(groot,'DefaultAxesLineStyleOrder',{'-','--','-.'});
%set(groot,'DefaultAxesColorOrder', cbrewer('qual', 'Dark2', 8));
set(groot,'DefaultAxesColorOrder', linspecer(8, 'qualitative'));
%set(groot,'DefaultAxesColorOrder', distinguishable_colors(15));

% axes
set(groot,'DefaultAxesFontName',fontName)
set(groot, 'DefaultAxesTitleFontWeight','normal');
set(groot,'DefaultAxesFontWeight','normal')
set(groot,'DefaultAxesTickLength'  , [.01 .01]);
set(groot,'DefaultAxesFontSize',22)
set(groot,'DefaultAxesBox','on')
set(groot,'DefaultAxesTickDir','out')
%set(groot,'DefaultAxesXMinorTick','on')
%set(groot,'DefaultAxesYMinorTick','on')
%set(groot,'DefaultAxesZMinorTick','on')
set(groot,'DefaultAxesXColor',[1 1 1]*0.1)
set(groot,'DefaultAxesYColor',[1 1 1]*0.1)
set(groot,'DefaultAxesZColor',[1 1 1]*0.1)
set(groot,'DefaultAxesLineWidth',1)
set(groot,'DefaultLineLineWidth',2);

% multiple monitor support
%mm = get(groot,'MonitorPositions');
%if size(mm,1) == 2,
%    disp('2 monitors detected.');
%    set(groot,'DefaultFigurePosition',[0.725*(mm(1,3)+mm(2,3)) 0.5*mm(2,4) 560 420]);
%end


%% change to current working dir
if ~strcmpi(computer,'GLNXA64')
   cd('E:\Work\eddyshelf\');
else
   mach = evalc('system hostname');
   if strfind(mach,'poison')
     cd('/home/poison/deepak/ROMS/runs/eddyshelf/scripts/');
   else
     cd('/media/data/Work/eddyshelf/');
   end
end

