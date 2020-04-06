classdef caseDistributionSystemMeasure < caseDistributionSystem
    % This is the class of distribution system. We assume all the
    % evaluations are conducted under practical measurements
    
    properties
        Property1
    end
    
    methods
        function obj = caseDistributionSystemMeasure(caseName, numSnap, range)
            % the construction function
            obj = obj@caseDistributionSystem(caseName, numSnap, range);
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 此处显示有关此方法的摘要
            %   此处显示详细说明
            outputArg = obj.Property1 + inputArg;
        end
    end
end

