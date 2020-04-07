classdef caseDistributionSystemMeasure < caseDistributionSystem
    % This is the class of distribution system. We assume all the
    % evaluations are conducted under practical measurements
    
    properties
        dataE               % the estimated data
    end
    
    methods
        function obj = caseDistributionSystemMeasure(caseName, numSnap, range)
            % the construction function
            obj = obj@caseDistributionSystem(caseName, numSnap, range);
        end
        
        function obj = preEvaluation(obj)
            % This method evaluate the parameters before approximating the
            % FIM. The evaluated value has low accuracy. We only use one
            % snapshot for the Vm and Va.
            
            % The first version is extremely simple
            
            % we first evaluate the vm
            obj.dataE.Vm = mean(obj.data.Vm_noised, 2);
            
            % We then evaluate the G and B. 
            obj.dataE.G = obj.data.G;
            
            obj.dataE.B = obj.data.B;
        end
    end
end

