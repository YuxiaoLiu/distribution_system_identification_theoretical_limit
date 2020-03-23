classdef caseDistributionSystem
    % This is the class of distribution system
    
    properties
        caseName            % the name of the power system case
        mpc                 % the matpower struct of the power system
        numBus              % the number of bus
        numSnap             % the number of snapshot
        range               % the struct of deviation range
        numFIM              % the struct representing the size of FIM matrix
        
        addressLoadRaw      % the address of raw load data
        addressLoad         % the address of preprocessed load data
        
        loadP               % the active load of each bus
        loadQ               % the reactive load of each bus
        FIM                 % the fisher information matrix
        
        data                % the data struct contains operation data
        sigma               % the variance of each meansurement noise
        isMeasure           % whether we have a specific measurement device
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
                inject = makeSbus(mpcThis.baseMVA, mpcThis.bus, mpcThis.gen);
                data_.P(:,i) = real(inject);
                data_.Q(:,i) = imag(inject);
                data_.Vm(:,i) = mpcThis.bus(:, 8);
                data_.Va(:,i) = mpcThis.bus(:, 9)/180*pi;
            end
            
            % assert that all the power flow converge
            assert(isempty(find(isSuccess == 0, 1)));
            
            % generate the G and B matrix
            Y = makeYbus(obj.mpc);
            data_.G = real(full(Y));
            data_.B = imag(full(Y));
            obj.data = data_;
        end
        
        function obj = setAccuracy(obj)
            % This method set the accuracy of the measurement device and
            % generate the measurement noise. This method also set whether
            % we have the measurement of a certain state.
            
            % we first set the relative noise ratio, we assume the noise
            % ratio is the sigma/mean value
            ratio.P = 0.001;
            ratio.Q = 0.001;
            ratio.Vm = 0.00001;
            ratio.Va = 0.005;
            % we then have to set whether we have the measurement device
            obj.isMeasure.P = true(obj.numBus, 1);
            obj.isMeasure.Q = true(obj.numBus, 1);
            obj.isMeasure.Vm = true(obj.numBus, 1);
            obj.isMeasure.Va = false(obj.numBus, 1); % false
            obj.isMeasure.Vm(1) = false;
            obj.isMeasure.Va(1) = false;
            
            % We assume there is no noise in the source bus. We set the
            % enlarge ratio of each rows of measurement noise.
            obj.sigma.P = mean(abs(obj.data.P), 2) * ratio.P;
            obj.sigma.Q = mean(abs(obj.data.Q), 2) * ratio.Q;
            obj.sigma.Vm = mean(abs(obj.data.Vm), 2) * ratio.Vm;
            obj.sigma.Va = mean(abs(obj.data.Va), 2) * ratio.Va;
            obj.sigma.Vm(1) = 0;
            obj.sigma.Va(1) = 0;
            
            % we generate the measurement noise
            rng(1);
            obj.data.P_noise = randn(obj.numBus, obj.numSnap);
            obj.data.P_noise = bsxfun(@times, obj.data.P_noise, obj.sigma.P);
            rng(2);
            obj.data.Q_noise = randn(obj.numBus, obj.numSnap);
            obj.data.Q_noise = bsxfun(@times, obj.data.Q_noise, obj.sigma.Q);
            rng(3);
            obj.data.Vm_noise = randn(obj.numBus, obj.numSnap);
            obj.data.Vm_noise = bsxfun(@times, obj.data.Vm_noise, obj.sigma.Vm);
            rng(4);
            obj.data.Va_noise = randn(obj.numBus, obj.numSnap);
            obj.data.Va_noise = bsxfun(@times, obj.data.Va_noise, obj.sigma.Va);
            
            % the measurement data
            obj.data.P_noised = obj.data.P + obj.data.P_noise;
            obj.data.Q_noised = obj.data.Q + obj.data.Q_noise;
            obj.data.Vm_noised = obj.data.Vm + obj.data.Vm_noise;
            obj.data.Va_noised = obj.data.Va + obj.data.Va_noise;
        end
        
        function obj = buildFIM(obj)
            % This method build the fisher information matrix (FIM). We
            % build the FIM in the order of measurement device or
            % measurement functions.
            
            % initialize the FIM matrix
            obj.numFIM.G = (1 + obj.numBus) * obj.numBus / 2;
            obj.numFIM.B = (1 + obj.numBus) * obj.numBus / 2;
            obj.numFIM.Vm = obj.numSnap * (obj.numBus - 1); % exclude the source bus
            obj.numFIM.Va = obj.numSnap * (obj.numBus - 1);
            obj.numFIM.Sum = obj.numFIM.G + obj.numFIM.B + obj.numFIM.Vm + obj.numFIM.Va;
            obj.FIM = zeros(obj.numFIM.Sum, obj.numFIM.Sum);
            
            % calculate the sub-matrix of P of all snapshots and all buses
            for i = 1:obj.numBus
                if obj.isMeasure.P(i)
                    for j = 1:obj.numSnap
                        obj = buildFIMP(obj, i, j);
                    end
                end
            end
            % calculate the sub-matrix of Q of all snapshots and all buses
            for i = 1:obj.numBus
                if obj.isMeasure.Q(i)
                    for j = 1:obj.numSnap
                        obj = buildFIMQ(obj, i, j);
                    end
                end
            end
            % calculate the sub-matrix of Vm of all snapshots and all buses
            for i = 1:obj.numBus
                if obj.isMeasure.Vm(i)
                    for j = 1:obj.numSnap
                        obj = buildFIMVm(obj, i, j);
                    end
                end
            end
            % calculate the sub-matrix of Va of all snapshots and all buses
            for i = 1:obj.numBus
                if obj.isMeasure.Va(i)
                    for j = 1:obj.numSnap
                        obj = buildFIMVa(obj, i, j);
                    end
                end
            end
        end
        
        function obj = buildFIMP(obj, bus, snap)
            % This method build the P part of FIM a selected bus and a selected snapshot. 
            % We first build a matrix, then we reshape the matrix to a vector. At
            % last we add up the FIM matrix. We conduct both G and B
            % matrix. Note that the state variables of G and B form a half
            % triangle, while the measurement function forms a whole
            % matrix.
            
            h = zeros(obj.numFIM.Sum, 1);
            theta_ij = obj.data.Va(bus, snap) - obj.data.Va(:, snap);
            Theta_ij = repmat(obj.data.Va(:, snap), 1, obj.numBus) - repmat(obj.data.Va(:, snap)', obj.numBus, 1);
            % G_ij\cos(\Theta_ij)+B_ij\sin(\Theta_ij)
            GBThetaP = obj.data.G .* cos(Theta_ij) + obj.data.B .* sin(Theta_ij);
            % G_ij\sin(\Theta_ij)-B_ij\cos(\Theta_ij)
            GBThetaQ = obj.data.G .* sin(Theta_ij) - obj.data.B .* cos(Theta_ij);
            
%             % verify the PF calculation
%             P = (GBThetaP * obj.data.Vm(:, snap)) .* obj.data.Vm(:, snap);
%             deltaP = P - obj.data.P(:, snap);
%             assert (sum(abs(deltaP)) <= 1e-6 );
%             Q = (GBThetaQ * obj.data.Vm(:, snap)) .* obj.data.Vm(:, snap);
%             deltaQ = Q - obj.data.Q(:, snap);
%             assert (sum(abs(deltaQ)) <= 1e-6 );
            
            % G matrix
            H_G = zeros(obj.numBus, obj.numBus);
            H_G(bus, :) = obj.data.Vm(bus, snap) * obj.data.Vm(:, snap)' .* cos(theta_ij');
            h_G = matToCol(obj, H_G);
            assert (length(h_G) == obj.numFIM.G);
            h(1:obj.numFIM.G) = h_G;
            
            % B matrix
            H_B = zeros(obj.numBus, obj.numBus);
            H_B(bus, :) = obj.data.Vm(bus, snap) * obj.data.Vm(:, snap)' .* sin(theta_ij');
            h_B = matToCol(obj, H_B);
            assert (length(h_B) == obj.numFIM.B);
            h(obj.numFIM.G+1:obj.numFIM.G+obj.numFIM.B) = h_B;
            
            % Vm
            % the first order term of other Vm
            H_Vm = zeros(obj.numBus, obj.numSnap);
            h_Vm = obj.data.Vm(bus, snap) * GBThetaP(:, bus);
            % the second order term of Vm(bus)
            h_Vm(bus) = 2*obj.data.Vm(bus, snap) * GBThetaP(bus, bus);
            % the first order term of Vm(bus)
            fOrderVm = obj.data.Vm(:, snap) .* GBThetaP(:, bus);
            fOrderVm(bus) = 0;
            h_Vm(bus) = h_Vm(bus) + sum(fOrderVm);
            H_Vm(:, snap) = h_Vm;
            % remove the source bus whose magnitude is not the state variable
            H_Vm(1, :) = []; 
            h_VmLarge = reshape(H_Vm', [], 1);
            h(obj.numFIM.G+obj.numFIM.B+1:obj.numFIM.G+obj.numFIM.B+obj.numFIM.Vm) = h_VmLarge;
            
            % Va
            H_Va = zeros(obj.numBus, obj.numSnap);
            h_Va = obj.data.Vm(bus, snap) * obj.data.Vm(:, snap) .* GBThetaQ(:, bus);
            h_Va(bus) = h_Va(bus)-sum(obj.data.Vm(bus, snap) * obj.data.Vm(:, snap) .* GBThetaQ(:, bus));
            H_Va(:, snap) = h_Va;
            % remove the source bus whose magnitude is not the state variable
            H_Va(1, :) = []; 
            h_VaLarge = reshape(H_Va', [], 1);
            h(obj.numFIM.G+obj.numFIM.B+obj.numFIM.Vm+1:end) = h_VaLarge;
            
            % build FIMP
            h = h / obj.sigma.P(bus);
            FIMP = h * h';
            obj.FIM = obj.FIM + FIMP;
        end
        
        function obj = buildFIMQ(obj, bus, snap)
            % This method build the Q part of FIM a selected bus and a selected snapshot. 
            
            h = zeros(obj.numFIM.Sum, 1);
            theta_ij = obj.data.Va(bus, snap) - obj.data.Va(:, snap);
            Theta_ij = repmat(obj.data.Va(:, snap), 1, obj.numBus) - repmat(obj.data.Va(:, snap)', obj.numBus, 1);
            % G_ij\cos(\Theta_ij)+B_ij\sin(\Theta_ij)
            GBThetaP = obj.data.G .* cos(Theta_ij) + obj.data.B .* sin(Theta_ij);
            % G_ij\sin(\Theta_ij)-B_ij\cos(\Theta_ij)
            GBThetaQ = obj.data.G .* sin(Theta_ij) - obj.data.B .* cos(Theta_ij);
            
            % G matrix
            H_G = zeros(obj.numBus, obj.numBus);
            H_G(bus, :) = obj.data.Vm(bus, snap) * obj.data.Vm(:, snap)' .* sin(theta_ij');
            h_G = matToCol(obj, H_G);
            h(1:obj.numFIM.G) = h_G;
            
            % B matrix
            H_B = zeros(obj.numBus, obj.numBus);
            H_B(bus, :) = - obj.data.Vm(bus, snap) * obj.data.Vm(:, snap)' .* cos(theta_ij');
            h_B = matToCol(obj, H_B);
            h(obj.numFIM.G+1:obj.numFIM.G+obj.numFIM.B) = h_B;
            
            % Vm
            % the first order term of other Vm
            H_Vm = zeros(obj.numBus, obj.numSnap);
            h_Vm = obj.data.Vm(bus, snap) * GBThetaQ(:, bus);
            % the second order term of Vm(bus)
            h_Vm(bus) = 2*obj.data.Vm(bus, snap) * GBThetaQ(bus, bus);
            % the first order term of Vm(bus)
            fOrderVm = obj.data.Vm(:, snap) .* GBThetaQ(:, bus);
            fOrderVm(bus) = 0;
            h_Vm(bus) = h_Vm(bus) + sum(fOrderVm);
            H_Vm(:, snap) = h_Vm;
            % remove the source bus whose magnitude is not the state variable
            H_Vm(1, :) = []; 
            h_VmLarge = reshape(H_Vm', [], 1);
            h(obj.numFIM.G+obj.numFIM.B+1:obj.numFIM.G+obj.numFIM.B+obj.numFIM.Vm) = h_VmLarge;
            
            % Va
            H_Va = zeros(obj.numBus, obj.numSnap);
            h_Va = - obj.data.Vm(bus, snap) * obj.data.Vm(:, snap) .* GBThetaP(:, bus);
            h_Va(bus) = h_Va(bus)+sum(obj.data.Vm(bus, snap) * obj.data.Vm(:, snap) .* GBThetaP(:, bus));
            H_Va(:, snap) = h_Va;
            % remove the source bus whose magnitude is not the state variable
            H_Va(1, :) = []; 
            h_VaLarge = reshape(H_Va', [], 1);
            h(obj.numFIM.G+obj.numFIM.B+obj.numFIM.Vm+1:end) = h_VaLarge;
            
            % build FIMP
            h = h / obj.sigma.Q(bus);
            FIMQ = h * h';
            obj.FIM = obj.FIM + FIMQ;
        end
        
        function obj = buildFIMVm(obj, bus, snap)
            % This method build the Vm part of FIM a selected bus. 
            h = zeros(obj.numFIM.Sum, 1);
            H_Vm = zeros(obj.numBus, obj.numSnap);
            H_Vm(bus, snap) = 1 / obj.sigma.Vm(bus);
            % remove the source bus whose magnitude is not the state variable
            H_Vm(1, :) = []; 
            h_VmLarge = reshape(H_Vm', [], 1);
            h(obj.numFIM.G+obj.numFIM.B+1:obj.numFIM.G+obj.numFIM.B+obj.numFIM.Vm) = h_VmLarge;
            
            FIMVm = h * h';
            obj.FIM = obj.FIM + FIMVm;
        end
        
        function obj = buildFIMVa(obj, bus, snap)
            % This method build the Va part of FIM a selected bus. 
            h = zeros(obj.numFIM.Sum, 1);
            H_Va = zeros(obj.numBus, obj.numSnap);
            H_Va(bus, snap) = 1 / obj.sigma.Va(bus);
            % remove the source bus whose magnitude is not the state variable
            H_Va(1, :) = []; 
            h_VaLarge = reshape(H_Va', [], 1);
            h(obj.numFIM.G+obj.numFIM.B+obj.numFIM.Vm+1:end) = h_VaLarge;
            
            FIMVa = h * h';
            obj.FIM = obj.FIM + FIMVa;
        end
        
        function h = matToCol(~, H)
            % this method transform the matrix into the column of the half
            % triangle.
            H_up = tril(H, -1)'+triu(H);
            n = size(H, 1);
            N = (n + 1) * n / 2;
            h = zeros(N, 1);
            pt = 1;
            for i = 1:n
                h(pt:pt+n-i) = H_up(i, i:end);
                pt = pt+n-i+1;
            end
        end
    end
end

