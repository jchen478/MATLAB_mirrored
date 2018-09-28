function [  ] = defaultFormat(  )


set ( 0 , 'DefaultAxesLineWidth' ,2 ,'DefaultAxesFontName' , 'CMU Serif',...
    'DefaultAxesFontSize' ,12 , 'DefaultAxesTickLength' , [ 0.025 0.025],...
    'DefaultLineLineWidth' ,2 , 'DefaultTextInterpreter' , 'Latex',...
    'DefaultLegendInterpreter' , 'Latex','defaultFigureColor','w')
set ( 0 , 'DefaultAxesTickLabelInterpreter' , 'Latex' )

colorOrder = [  rgb('DarkRed'); rgb('Crimson'); rgb('OrangeRed');
    rgb('Orange'); rgb('Gold'); rgb('Lime');
    rgb('Olive'); rgb('DarkGreen'); rgb('LightSkyBlue');
    rgb('MediumBlue'); rgb('Plum'); rgb('Purple')]

% colorOrder = [  rgb('Crimson'); rgb('OrangeRed');
%     rgb('Orange'); rgb('DarkGreen'); 
%     rgb('MediumBlue'); rgb('Purple')]
% 
% colorOrder = [ 
%     rgb('Orange'); rgb('DarkGreen'); 
%     rgb('MediumBlue'); rgb('Purple')]
% 
% colorOrder = [ rgb('Crimson');
%     rgb('Orange'); rgb('Gold');
%     rgb('DarkGreen'); rgb('LightSkyBlue');
%     rgb('MediumBlue'); rgb('Purple')]; 

colorOrder = [ rgb('Crimson'); rgb('Orange');
    rgb('MediumBlue'); rgb('DarkGreen'); rgb('Purple')];
% colorOrder = [rgb('Orange');
%     rgb('MediumBlue'); rgb('DarkGreen'); rgb('Purple')]; 

% colorOrder = [ rgb('Crimson'); rgb('MediumBlue')];
% colorOrder = [ rgb('Orange'); rgb('MediumBlue'); rgb('DarkGreen')];
set(groot,'DefaultAxesColorOrder',colorOrder, ...
    'DefaultAxesLineStyleOrder','-|--|:|-.');

end

