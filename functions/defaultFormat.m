function [  ] = defaultFormat(  )


set ( 0 , 'DefaultAxesLineWidth' ,2 ,'DefaultAxesFontName' , 'CMU Serif',...
'DefaultAxesFontSize' ,20 , 'DefaultAxesTickLength' , [ 0.025 0.015],...
'DefaultLineLineWidth' ,2 , 'DefaultTextInterpreter' , 'Latex',...
'DefaultLegendInterpreter' , 'Latex','defaultFigureColor','w')
set ( 0 , 'DefaultAxesTickLabelInterpreter' , 'Latex' )

colorOrder = [  rgb('DarkRed'); rgb('Crimson'); rgb('OrangeRed');
    rgb('Orange'); rgb('Gold'); rgb('Lime');
    rgb('Olive'); rgb('DarkGreen'); rgb('LightSkyBlue');
    rgb('MediumBlue'); rgb('Plum'); rgb('Purple')] 

set(groot,'DefaultAxesColorOrder',colorOrder, ...  
    'DefaultAxesLineStyleOrder','-|--|:|-.');

end

