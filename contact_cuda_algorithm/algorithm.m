clc;
clear;

tot = 6; 

fiberarrays = struct;
fiberarrays(1).list = [3]; 
fiberarrays(2).list = [3 5]; 
fiberarrays(3).list = [1 4 2]; 
fiberarrays(4).list = [3 5]; 
fiberarrays(5).list = [2 4]; 
fiberarrays(6).list = []; 
connections = 0; 
listofconnections = zeros(10,1); 
fiberstatus = zeros(tot,1); 
access = zeros(tot,1); 
for i=1:tot
    fiberstatus = Listgen(i,i,fiberarrays,fiberstatus);
end

fiberstatus

[listofconnections,connectionm,access] = compilelist(1,1,fiberarrays, fiberstatus,...
    listofconnections, connections, access);

listofconnections

