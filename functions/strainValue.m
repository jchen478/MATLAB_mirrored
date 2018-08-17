classdef strainValue
    %strainValue Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        block_strain
        block_etarel
        name
        block_ndata
        block_detarel
    end
    
    methods
        function obj = strainValue()
            obj.block_ndata = 0; 
        end
    end
end
