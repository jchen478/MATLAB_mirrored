classdef caseArr
    %CASEARR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        value
        legend
        dim
        ndata
    end
    
    methods
        function obj = caseArr(nameInput,valueInput,dimInput)
            if nargin > 0
                obj.dim = dimInput;
                obj.name = nameInput;
                obj.value = valueInput;
                obj.ndata = length(obj.value);
                obj.legend = cell(obj.ndata,1);
                for i=1:obj.ndata
                    obj.legend{i} = [obj.name,' = ', num2str(obj.value(i))];
                end
            end
        end
        function obj = pluck(obj,ind)
            obj.value(ind) = [];
            obj.legend(ind) = [];
            obj.ndata = obj.ndata - 1; 
        end
    end
    
end

