classdef caseDistributionSystemMeasure < caseDistributionSystem
    % This is the class of distribution system. We assume all the
    % evaluations are conducted under practical measurements
    
    properties
        Property1
    end
    
    methods
        function obj = caseDistributionSystemMeasure(inputArg1,inputArg2)
            %UNTITLED ��������ʵ��
            %   �˴���ʾ��ϸ˵��
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 �˴���ʾ�йش˷�����ժҪ
            %   �˴���ʾ��ϸ˵��
            outputArg = obj.Property1 + inputArg;
        end
    end
end

