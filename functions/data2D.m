classdef data2D
    %CASEARR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        value
    end
    
    methods
        function obj = data2D(nameInput,valueInput)
            if nargin > 0                
                obj.name = nameInput;
                obj.value = valueInput;
            end
        end
    end
    
end

