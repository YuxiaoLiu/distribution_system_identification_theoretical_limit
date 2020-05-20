classdef caseDistributionSystemMeasure < caseDistributionSystem
    % This is the class of distribution system. We assume all the
    % evaluations are conducted under practical measurements
    
    properties
        dataE               % the estimated data
        boundA              % the approximated bound
        sigmaReal           % the deviation of the real state variables
        prior               % the prior assumptions of the G and B matrix
        
        A_FIM               % the approximated fisher information matrix
        A_FIMP              % the (sparse) FIM of active power injection
        A_FIMQ              % the (sparse) FIM of reactive power injection
        
        initPar             % the initial estimation of parameters and state variables
    end
    
    methods
        function obj = caseDistributionSystemMeasure(caseName, numSnap, range)
            % the construction function
            obj = obj@caseDistributionSystem(caseName, numSnap, range);
        end
        
        function obj = preEvaluation(obj, varargin)
            % This method evaluate the parameters before approximating the
            % FIM. The evaluated value has low accuracy. We only use one
            % snapshot for the Vm and Va.
            
            if nargin == 2
                obj.prior = varargin{1};
            elseif nargin == 1
                obj.prior.Gmin = 0.1;
                obj.prior.Bmin = 0.1;
                obj.prior.ratio = 0.05;
            end
            
            % we first evaluate the vm and the va
%             obj.dataE.Vm = obj.data.Vm;%_noised;
%             obj.dataE.Va = obj.data.Va;%_noised;

            obj.sigmaReal.Vm = cov(obj.data.Vm');
            mu = mean(obj.data.Vm, 2);
            rng(5);
            obj.dataE.Vm = mvnrnd(mu, obj.sigmaReal.Vm, obj.numSnap)';
            
            obj.sigmaReal.Va = cov(obj.data.Va');
            mu = mean(obj.data.Va, 2);
            rng(6);
            obj.dataE.Va = mvnrnd(mu, obj.sigmaReal.Va, obj.numSnap)';
            
%             obj.dataE.Vm = obj.data.Vm;%_noised;
            
            % We then evaluate the G and B. 
%             obj.dataE.G = obj.data.G;
%             obj.dataE.B = obj.data.B;
            
            obj = approximateY(obj);
            
            obj.dataE.Va = zeros(obj.numBus, obj.numSnap);
            obj.dataE.Va(2:end, :) = - (obj.dataE.G(2:end, 2:end)) \ obj.data.P_noised(2:end, :);
%             mu = mean(obj.data.Va, 2);
%             obj.sigmaReal.P = cov(obj.data.P');
%             obj.sigmaReal.Va = zeros(obj.numBus, obj.numBus);
%             obj.sigmaReal.Va(2:end, 2:end) = ...
%                 ((1.5*obj.dataE.G(2:end, 2:end)) \ obj.sigmaReal.P(2:end, 2:end)) / (1.5*obj.dataE.G(2:end, 2:end));
%             rng(7);
%             obj.dataE.Va = mvnrnd(mu, obj.sigmaReal.Va, obj.numSnap)';
%             obj.dataE.Va(2:end, :) = -obj.dataE.G(2:end, 2:end) \ obj.data.P_noised(2:end, :);
        end
        
        function obj = approximateFIM(obj, varargin)
            % This method approximate the fisher information matrix based
            % on the pre-evaluation results of the parameters.
            if nargin == 2
                obj.k = varargin{1};
            elseif nargin == 1
                obj.k.G = 5;
                obj.k.B = 10;
                obj.k.vm = 10;
                obj.k.va = 1000;
            end
            % initialize the A_FIM matrix
            obj.numFIM.G = (1 + obj.numBus) * obj.numBus / 2;
            obj.numFIM.B = (1 + obj.numBus) * obj.numBus / 2;
            obj.numFIM.Vm = obj.numSnap * (obj.numBus - 1); % exclude the source bus
            obj.numFIM.Va = obj.numSnap * (obj.numBus - 1);
            obj.numFIM.Sum = obj.numFIM.G + obj.numFIM.B + obj.numFIM.Vm + obj.numFIM.Va;
            
            obj.A_FIM = zeros(obj.numFIM.Sum, obj.numFIM.Sum);
            obj.A_FIMP = sparse(obj.numFIM.Sum, obj.numFIM.Sum);
            obj.A_FIMQ = sparse(obj.numFIM.Sum, obj.numFIM.Sum);
            obj.FIMVm = sparse(obj.numFIM.Sum, obj.numFIM.Sum);
            obj.FIMVa = sparse(obj.numFIM.Sum, obj.numFIM.Sum);
            
            % calculate the sub-matrix of P of all snapshots and all buses
            for i = 1:obj.numBus
                if obj.isMeasure.P(i)
                    for j = 1:obj.numSnap
                        obj = approximateFIMP(obj, i, j);
                    end
                end
            end
            obj.A_FIM = obj.A_FIM + full(obj.A_FIMP);
            % calculate the sub-matrix of Q of all snapshots and all buses
            for i = 1:obj.numBus
                if obj.isMeasure.Q(i)
                    for j = 1:obj.numSnap
                        obj = approximateFIMQ(obj, i, j);
                    end
                end
            end
            obj.A_FIM = obj.A_FIM + full(obj.A_FIMQ);
            % calculate the sub-matrix of Vm of all snapshots and all buses
            for i = 1:obj.numBus
                if obj.isMeasure.Vm(i)
                    for j = 1:obj.numSnap
                        obj = buildFIMVm(obj, i, j);
                    end
                end
            end
            obj.A_FIM = obj.A_FIM + full(obj.FIMVm);
            % calculate the sub-matrix of Va of all snapshots and all buses
            for i = 1:obj.numBus
                if obj.isMeasure.Va(i)
                    for j = 1:obj.numSnap
                        obj = buildFIMVa(obj, i, j);
                    end
                end
            end
            obj.A_FIM = obj.A_FIM + full(obj.FIMVa);
        end
        
        function obj = approximateFIMP(obj, bus, snap)
            % This method approximate the P part of FIM. We ignore the sin
            % part of the power flow equations.
            h = sparse(obj.numFIM.Sum, 1);
            theta_ij = obj.dataE.Va(bus, snap) - obj.dataE.Va(:, snap);
%             Theta_ij = repmat(obj.dataE.Va(:, snap), 1, obj.numBus) - repmat(obj.dataE.Va(:, snap)', obj.numBus, 1);
%             % G_ij\cos(\Theta_ij)+B_ij\sin(\Theta_ij)
%             GBThetaP = obj.dataE.G .* cos(Theta_ij) + obj.dataE.B .* sin(Theta_ij);
%             % G_ij\sin(\Theta_ij)-B_ij\cos(\Theta_ij)
%             GBThetaQ = obj.dataE.G .* sin(Theta_ij) - obj.dataE.B .* cos(Theta_ij);
            
            % G matrix
            H_G = zeros(obj.numBus, obj.numBus);
            H_G(bus, :) = obj.dataE.Vm(bus, snap) * obj.dataE.Vm(:, snap)' / obj.k.G; % .* cos(theta_ij')
            h_G = obj.matToCol(H_G);
            h(1:obj.numFIM.G) = h_G;
            
            % B matrix
            H_B = zeros(obj.numBus, obj.numBus);
            H_B(bus, :) = obj.dataE.Vm(bus, snap) * obj.dataE.Vm(:, snap)' .* sin(theta_ij') / obj.k.B;
            h_B = obj.matToCol(H_B);
            h(obj.numFIM.G+1:obj.numFIM.G+obj.numFIM.B) = h_B;
            
            % Vm
            % the first order term of other Vm
            H_Vm = zeros(obj.numBus, obj.numSnap);
            h_Vm = obj.dataE.Vm(bus, snap) * obj.dataE.G(:, bus) / obj.k.vm; % obj.dataE.G(:, bus)
            % the second order term of Vm(bus)
            h_Vm(bus) = 2*obj.dataE.Vm(bus, snap) * obj.dataE.G(bus, bus) / obj.k.vm; % obj.dataE.G(bus, bus)
            % the first order term of Vm(bus)
            fOrderVm = obj.dataE.Vm(:, snap) .* obj.dataE.G(:, bus) / obj.k.vm; % obj.dataE.G(:, bus)
            fOrderVm(bus) = 0;
            h_Vm(bus) = h_Vm(bus) + sum(fOrderVm);
            H_Vm(:, snap) = h_Vm;
            % remove the source bus whose magnitude is not the state variable
            H_Vm(1, :) = []; 
            h_VmLarge = reshape(H_Vm', [], 1);
            h(obj.numFIM.G+obj.numFIM.B+1:obj.numFIM.G+obj.numFIM.B+obj.numFIM.Vm) = h_VmLarge;
            
            % Va
            H_Va = zeros(obj.numBus, obj.numSnap);
            h_Va = obj.dataE.Vm(bus, snap) * obj.dataE.Vm(:, snap) .* (- obj.dataE.B(:, bus)) / obj.k.va; % (- obj.dataE.B(:, bus))
            h_Va(bus) = h_Va(bus)-sum(obj.dataE.Vm(bus, snap) * obj.dataE.Vm(:, snap) .* (- obj.dataE.B(:, bus))) / obj.k.va; % (- obj.dataE.B(:, bus)))
            H_Va(:, snap) = h_Va;
            % remove the source bus whose magnitude is not the state variable
            H_Va(1, :) = []; 
            h_VaLarge = reshape(H_Va', [], 1);
            h(obj.numFIM.G+obj.numFIM.B+obj.numFIM.Vm+1:end) = h_VaLarge;
            
            % build FIMP
            h = h / obj.sigma.P(bus);
            FIMPThis = h * h';
            obj.A_FIMP = obj.A_FIMP + FIMPThis;
        end
        
        function obj = approximateFIMQ(obj, bus, snap)
            % This method approximate the Q part of FIM. We ignore the sin
            % part of the power flow equations.
            h = sparse(obj.numFIM.Sum, 1);
            theta_ij = obj.dataE.Va(bus, snap) - obj.dataE.Va(:, snap);
%             Theta_ij = repmat(obj.dataE.Va(:, snap), 1, obj.numBus) - repmat(obj.dataE.Va(:, snap)', obj.numBus, 1);
%             % G_ij\cos(\Theta_ij)+B_ij\sin(\Theta_ij)
%             GBThetaP = obj.dataE.G .* cos(Theta_ij) + obj.dataE.B .* sin(Theta_ij);
%             % G_ij\sin(\Theta_ij)-B_ij\cos(\Theta_ij)
%             GBThetaQ = obj.dataE.G .* sin(Theta_ij) - obj.dataE.B .* cos(Theta_ij);
            
            % G matrix
            H_G = zeros(obj.numBus, obj.numBus);
            H_G(bus, :) = obj.dataE.Vm(bus, snap) * obj.dataE.Vm(:, snap)' .* sin(theta_ij') / obj.k.G;
            h_G = obj.matToCol(H_G);
            h(1:obj.numFIM.G) = h_G;
            
            % B matrix
            H_B = zeros(obj.numBus, obj.numBus);
            H_B(bus, :) = - obj.dataE.Vm(bus, snap) * obj.dataE.Vm(:, snap)' / obj.k.B; %  .* cos(theta_ij')
            h_B = obj.matToCol(H_B);
            h(obj.numFIM.G+1:obj.numFIM.G+obj.numFIM.B) = h_B;
            
            % Vm
            % the first order term of other Vm
            H_Vm = zeros(obj.numBus, obj.numSnap);
            h_Vm = obj.dataE.Vm(bus, snap) * (-obj.dataE.B(:, bus)) / obj.k.vm; % (-obj.dataE.B(:, bus))
            % the second order term of Vm(bus)
            h_Vm(bus) = 2*obj.dataE.Vm(bus, snap) * (-obj.dataE.B(bus, bus)) / obj.k.vm; % (-obj.dataE.B(bus, bus))
            % the first order term of Vm(bus)
            fOrderVm = obj.dataE.Vm(:, snap) .* (-obj.dataE.B(:, bus)) / obj.k.vm; % (-obj.dataE.B(:, bus))
            fOrderVm(bus) = 0;
            h_Vm(bus) = h_Vm(bus) + sum(fOrderVm);
            H_Vm(:, snap) = h_Vm;
            % remove the source bus whose magnitude is not the state variable
            H_Vm(1, :) = []; 
            h_VmLarge = reshape(H_Vm', [], 1);
            h(obj.numFIM.G+obj.numFIM.B+1:obj.numFIM.G+obj.numFIM.B+obj.numFIM.Vm) = h_VmLarge;
            
            % Va
            H_Va = zeros(obj.numBus, obj.numSnap);
            h_Va = - obj.dataE.Vm(bus, snap) * obj.dataE.Vm(:, snap) .* obj.dataE.G(:, bus) / obj.k.va; % obj.dataE.G(:, bus)
            h_Va(bus) = h_Va(bus)+sum(obj.dataE.Vm(bus, snap) * obj.dataE.Vm(:, snap) .* obj.dataE.G(:, bus)) / obj.k.va; % obj.dataE.G(:, bus))
            H_Va(:, snap) = h_Va;
            % remove the source bus whose magnitude is not the state variable
            H_Va(1, :) = []; 
            h_VaLarge = reshape(H_Va', [], 1);
            h(obj.numFIM.G+obj.numFIM.B+obj.numFIM.Vm+1:end) = h_VaLarge;
            
            % build FIMQ
            h = h / obj.sigma.Q(bus);
            FIMQThis = h * h';
            obj.A_FIMQ = obj.A_FIMQ + FIMQThis;
        end
        
        function obj = calABound(obj, varargin)
            % this method calculate the bound from the A_FIM matrix;
            
            if nargin == 3
                obj.admittanceOnly = varargin{1};
                obj.topoPrior = varargin{2};
            elseif nargin == 2
                obj.admittanceOnly = varargin{1};
                obj.topoPrior = false(obj.numBus, obj.numBus);
            elseif nargin == 1
                obj.admittanceOnly = false;
                obj.topoPrior = false(obj.numBus, obj.numBus);
            end
            
            % build the indexes we really care about
            delCols = [obj.matToCol(obj.topoPrior)>1e-4;obj.matToCol(obj.topoPrior)>1e-4];
            obj.numFIM.index = true(obj.numFIM.Sum, 1);
            obj.numFIM.index(delCols) = false;
            obj.numFIM.del = sum(delCols)/2;
            
            % for [A B; B' C], we calculate A-B/C*B'
            if obj.admittanceOnly
                obj.numFIM.index = obj.numFIM.index(1:obj.numFIM.G+obj.numFIM.B);
                A = obj.A_FIM(1:obj.numFIM.G+obj.numFIM.B, 1:obj.numFIM.G+obj.numFIM.B);
                B = obj.A_FIM(1:obj.numFIM.G+obj.numFIM.B, obj.numFIM.G+obj.numFIM.B+1:end);
                C = obj.A_FIM(obj.numFIM.G+obj.numFIM.B+1:end, obj.numFIM.G+obj.numFIM.B+1:end);
                obj.A_FIM = A - B/C*B';
                var = diag(obj.A_FIM(obj.numFIM.index, obj.numFIM.index)\eye(sum(obj.numFIM.index)));
            else
                var = diag(obj.A_FIM(obj.numFIM.index, obj.numFIM.index)\eye(sum(obj.numFIM.index)));
%                 % we construct a Hermitian matrix H and use Cholesky
%                 % decomposition to compute the inverse matrix
%                 FIM = obj.A_FIM(obj.numFIM.index, obj.numFIM.index);
%                 H = FIM * FIM';
%                 U = chol(H);
%                 Uinv = U \ eye(size(U));
%                 Cov = H' * (Uinv * Uinv');
            end
            if min(var) < 0
                var = abs(var);
                fprintf('We use the absolute value of the variance.\n');
            end
            
            obj.boundA.total = sqrt(var);
            
            boundG = zeros(obj.numFIM.G, 1);
            boundG(obj.numFIM.index(1:obj.numFIM.G)) = obj.boundA.total(1:obj.numFIM.G-obj.numFIM.del) / obj.k.G;
            obj.boundA.total(1:obj.numFIM.G-obj.numFIM.del) = obj.boundA.total(1:obj.numFIM.G-obj.numFIM.del) / obj.k.G;
            obj.boundA.G = obj.colToMat(boundG, obj.numBus);
            
            boundB = zeros(obj.numFIM.B, 1);
            boundB(obj.numFIM.index(1:obj.numFIM.G)) = ...
                obj.boundA.total(obj.numFIM.G+1-obj.numFIM.del:obj.numFIM.G+obj.numFIM.B-2*obj.numFIM.del) / obj.k.B;
            obj.boundA.total(obj.numFIM.G+1-obj.numFIM.del:obj.numFIM.G+obj.numFIM.B-2*obj.numFIM.del) = ...
                obj.boundA.total(obj.numFIM.G+1-obj.numFIM.del:obj.numFIM.G+obj.numFIM.B-2*obj.numFIM.del) / obj.k.B;
            obj.boundA.B = obj.colToMat(boundB, obj.numBus);
            
            obj.boundA.G_relative = abs(obj.boundA.G ./ repmat(diag(obj.data.G), 1, obj.numBus));
            obj.boundA.B_relative = abs(obj.boundA.B ./ repmat(diag(obj.data.B), 1, obj.numBus));
            obj.boundA.G_relative_col = reshape(obj.boundA.G_relative, [], 1);
            obj.boundA.B_relative_col = reshape(obj.boundA.B_relative, [], 1);
            
            if ~obj.admittanceOnly
                obj.boundA.Vm = ...
                    obj.boundA.total(obj.numFIM.G+obj.numFIM.B+1-2*obj.numFIM.del...
                    :obj.numFIM.G+obj.numFIM.B+obj.numFIM.Vm-2*obj.numFIM.del) / obj.k.vm;
                obj.boundA.total(obj.numFIM.G+obj.numFIM.B+1-2*obj.numFIM.del...
                    :obj.numFIM.G+obj.numFIM.B+obj.numFIM.Vm-2*obj.numFIM.del)...
                    = obj.boundA.Vm;
                obj.boundA.Va = ...
                    obj.boundA.total(obj.numFIM.G+obj.numFIM.B+obj.numFIM.Vm+1-2*obj.numFIM.del...
                    :obj.numFIM.Sum-2*obj.numFIM.del) / obj.k.va;
                obj.boundA.total(obj.numFIM.G+obj.numFIM.B+obj.numFIM.Vm+1-2*obj.numFIM.del...
                    :obj.numFIM.Sum-2*obj.numFIM.del)...
                    = obj.boundA.Va;
            end
        end
        
        function obj = approximateY(obj)
            % This method approximate the Y matrix by the measurements. We
            % use the simple Ohom's law to provide an initial value of Y.
            % We also assume the G/B ratio is a constant.
            rng(103);
            randG = 0.75 + 0.5 * randn(size(obj.data.G));
            rng(104);
            randB = 0.75 + 0.5 * randn(size(obj.data.B));
            obj.dataE.G = obj.data.G .* randG;
            obj.dataE.B = obj.data.B .* randB;
%             % The approximation of the diagonal elements
%             diagG = diag(obj.dataE.G);
%             diagB = diag(obj.dataE.B);
            
            % approximate the topology using Vm data only
            % ranking the Vm
            Vm = obj.data.Vm_noised;
            Topo = logical(eye(obj.numBus));
            VmMean = mean(Vm, 2);
            [~, VmOrder] = sort(VmMean,'descend');
            
            assert (VmOrder(1) == 1); % the first bus is the source bus
            
            Vm = movmean(Vm, floor(obj.numSnap/20)+1, 2);
            corr = corrcoef(Vm');
            corr(isnan(corr)) = 0; % one can also simulate some disturbance in the source bus voltage
            for i = 2:obj.numBus
                % iterate each bus
                [~, loc] = max(corr(VmOrder(i), VmOrder(1:i-1))); % the location of the connected bus
                Topo(VmOrder(i), VmOrder(loc)) = true;
                Topo(VmOrder(loc), VmOrder(i)) = true;
            end
            
%             Topo = obj.data.G~=0;

            % approximate the parameter
            IP = obj.data.IP_noised;
            IQ = obj.data.IQ_noised;
%             G = IP * obj.data.Vm_noised' / (obj.data.Vm_noised * obj.data.Vm_noised');
%             B = - IQ * obj.data.Vm_noised' / (obj.data.Vm_noised * obj.data.Vm_noised');
            G_ols = zeros(obj.numBus, obj.numBus);
            B_ols = zeros(obj.numBus, obj.numBus);
            for i = 1:obj.numBus
                j = VmOrder(i);
                filter = Topo(:, j);
                filter(j) = false;
                previous = VmOrder(1:i);
                previous = intersect(previous, find(filter));
                
                VmDelta = Vm(filter, :) - repmat(Vm(j, :), sum(filter), 1);
                yG = IP(i, :);
                yB = IQ(i, :);
                try
                    yG = yG - G_ols(previous, j) * VmDelta(filter(previous), :);
                    yB = yB + B_ols(previous, j) * VmDelta(filter(previous), :);
                catch
                    assert (i == 1);
                end
                
                rng(i);
                filter(previous) = false;
                VmDelta = Vm(filter, :) - repmat(Vm(j, :), sum(filter), 1);
%                 G_ols(j, filter) = obj.tls(VmDelta', yG');
                G_ols(j, filter) = yG * VmDelta' / (VmDelta * VmDelta');
                outlier = G_ols(j,:) > -obj.prior.Gmin;
                G_ols(j, filter & outlier') = - obj.prior.Gmin * (1+0.1*rand());
                G_ols(filter, j) = G_ols(j, filter);
                G_ols(j, j) = -sum(G_ols(j, :));
                
                B_ols(j, filter) = - yB * VmDelta' / (VmDelta * VmDelta');
                outlier = B_ols(j,:) < obj.prior.Bmin;
                B_ols(j, filter & outlier') = obj.prior.Bmin * (1+0.1*rand());
                B_ols(filter, j) = B_ols(j, filter);
                B_ols(j, j) = -sum(B_ols(j, :));
            end

            obj.dataE.G = G_ols;
            obj.dataE.B = B_ols;
            
            obj.dataE.G = obj.data.G;
            obj.dataE.B = obj.data.B;
        end
        
        function obj = initValue(obj)
            % This method provides the initial value (voltage angles?)
            
        end
        
        function obj = identifyMCMC(obj)
            % This method uses the Markov Chain Monte Carlo to sample the
            % distribution of the parameters and the topologies
            
            % Build the measurement function.
            % Currently, we assume we know all the P, Q, Vm, Va
            % measurements. We have to modify it later.
            data.Pn = obj.data.P_noised;
            data.Qn = obj.data.Q_noised;
            data.Vmn = obj.data.Vm_noised;
            data.Van = obj.data.Va_noised;
            
            par.G = obj.dataE.G;
            par.B = obj.dataE.B;
            par.Vm = obj.data.Vm_noised;
            par.Va = obj.data.Va_noised;
        end
    end
    
    methods (Static)
        function B = tls(xdata,ydata)
            % This method con
            SUM = sum(xdata,1);
            zero = find(SUM==0);
            xdata(:,zero)=[];

%             m       = length(ydata);       %number of x,y data pairs
            X       = xdata;
            Y       = ydata;
            n       = size(X,2);          % n is the width of X (X is m by n)
            Z       = [X Y];              % Z is X augmented with Y.
            [~, ~, V] = svd(Z,0);         % find the SVD of Z.
            VXY     = V(1:n,1+n:end);     % Take the block of V consisting of the first n rows and the n+1 to last column
            VYY     = V(1+n:end,1+n:end); % Take the bottom-right block of V.
            B       = -VXY/VYY;

            for i = zero
                B = [B(1:i-1); 0; B(i:end)];
            end
            end
    end
end

