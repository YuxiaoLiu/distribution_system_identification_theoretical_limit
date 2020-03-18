classdef caseDistributionSystem
    % This is the class of distribution system
    
    properties
        caseName            % the name of the power system case
        mpc                 % the matpower struct of the power system
        numBus              % the number of bus
        numSnap             % the number of snapshot
        range               % % the deviation range
        
        addressLoadRaw      % the address of raw load data
        addressLoad         % the address of preprocessed load data
        
        loadP               % the active load of each bus
        loadQ               % the reactive load of each bus
        data                % the data struct contains operation data
    end
    
    methods
        function obj = caseDistributionSystem(caseName, numSnap, range)
            % the construction function
            obj.caseName = caseName;
            obj.addressLoadRaw = '.\data\file1csv.csv';
            obj.addressLoad = '.\data\dataLoad.csv';
            
            % load the distribution system
            obj.mpc = loadcase(caseName);
            obj.numBus = size(obj.mpc.bus, 1);
            obj.numSnap = numSnap;
            obj.range = range;
                       
        end
        
        function readLoadRaw(obj)
            % this method read and process the raw load data
            
            numDay = 20;            % we read 20 days of load data
            numCustomer = 979;      % the number of costomer of the first day
            loadRaw = xlsread(obj.addressLoadRaw);
            % read the data by rows
            load = zeros(numCustomer, numDay*48);
            idCustomer = loadRaw(1:numCustomer, 1);
            numRow = size(loadRaw, 1);
            idDay = 0;
            for i = 1:numRow
                if loadRaw(i, 2) > idDay % read the data of a new day
                    idDay = loadRaw(i, 2);
                end
                idRow = find(idCustomer == loadRaw(i, 1));
                if ~isempty(idRow)
                    rangeDay = (idDay-195)*48+1 : (idDay-194)*48;
                    load(idRow, rangeDay) = loadRaw(i, 3:end);
                end
            end
            % output the preprocessed load data
            xlswrite(obj.addressLoad, load);
        end
        
        function obj = readLoad(obj)
            % this method read the prepocessed load data, aggregate the
            % data, and cut the data into the appropriate size
            numAggregation = 5;     % aggregate serveral loads together
            loadRead = xlsread(obj.addressLoad);
            [numCust, numSnapRaw] = size(loadRead);
            numCustAggre = fix(numCust / numAggregation);
            load = zeros(numCustAggre, numSnapRaw);
            
            % aggregate and normalize the data
            idRow = 1;
            for i = 1:numCust
                if (mod(i,numAggregation) == 0)
                    custRange = i-numAggregation+1:i;
                    thisLoad = sum(loadRead(custRange,:));
                    load(idRow,:) = thisLoad/max(thisLoad);
                    idRow = idRow + 1;
                end
            end
            
            % cut the data
            load(obj.numBus:end,:) = []; % exclude the source bus
            load(:,obj.numSnap+1:end) = [];
            
            % rescale the data
            load = 1 - obj.range.P/2 + load*obj.range.P;
            obj.loadP = load;
            
            % generate the reactive load data
            rng(1);
            randQ = rand(size(load)) * obj.range.Q + 1 - obj.range.Q/2;
            obj.loadQ = load .* randQ;
        end
        
        function obj = genOperateData(obj)
            % this method generate the steady state operation data by
            % running power flow equations
            data_.P = zeros(obj.numBus, obj.numSnap);
            data_.Q = zeros(obj.numBus, obj.numSnap);
            data_.Vm = zeros(obj.numBus, obj.numSnap);
            data_.Va = zeros(obj.numBus, obj.numSnap);
            isSuccess = ones(obj.numSnap, 1);
            
            for i = 1:obj.numSnap
                mpcThis = obj.mpc;
                % update active and reactive load
                mpcThis.bus(2:end,3) = mpcThis.bus(2:end,3) .* obj.loadP(:, i);
                mpcThis.bus(2:end,4) = mpcThis.bus(2:end,4) .* obj.loadQ(:, i);
                % run power flow
                mpopt = mpoption('verbose',0,'out.all',0);
                mpcThis = runpf(mpcThis, mpopt);
                isSuccess(i, 1) = mpcThis.success;
                % output the data
                data_.P(:,i) = mpcThis.bus(:, 3);
                data_.Q(:,i) = mpcThis.bus(:, 4);
                data_.Vm(:,i) = mpcThis.bus(:, 8);
                data_.Va(:,i) = mpcThis.bus(:, 9)/180*pi;
            end
            
            % assert that all the power flow converge
            assert(isempty(find(isSuccess == 0, 1)));
            obj.data = data_;
        end
    end
end

