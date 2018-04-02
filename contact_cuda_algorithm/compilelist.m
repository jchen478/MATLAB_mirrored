function [listofconnections,connections,access]=compilelist(fiber,...
    thread,fiberarray,fiberstatus,listofconnections,connections,access)
   if fiberstatus(thread)==0
       access(thread) = 1;
       N = size(fiberarray(fiber).list,2);
       for n=1:N
          locfiber=fiberarray(fiber).list(n);
          %if locfiber~=lastfiber 
          if access(locfiber) == 0
              connections=connections+1;
              access(locfiber)=access(locfiber)+1
              listofconnections(connections)=locfiber;
              [listofconnections,connections,access]=compilelist(locfiber,...
                  thread,fiberarray,fiberstatus,...
                  listofconnections,connections,access);             
          end
          %lastfiber=fiber;
       end
   else
       disp('not a boss')
   end
end

