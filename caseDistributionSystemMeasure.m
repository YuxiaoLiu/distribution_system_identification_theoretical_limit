classdef caseDistributionSystemMeasure < caseDistributionSystem
    % This is the class of distribution system. We assume all the
    % evaluations are conducted under practical measurements
    
    properties
        dataE               % the estimated data
        dataO               % the process data in the optimization iterations
        boundA              % the approximated bound
        sigmaReal           % the deviation of the real state variables
        prior               % the prior assumptions of the G and B matrix
        
        A_FIM               % the approximated fisher information matrix
        A_FIMP              % the (sparse) FIM of active power injection
        A_FIMQ              % the (sparse) FIM of reactive power injection
        
        initPar             % the initial estimation of parameters and state variables
        truePar             % the ground truth of the parameters
        
        grad                % the gradient vector
        gradChain           % the chain of the gradients
        gradP               % the gradient vector from the measurement of P
        gradQ               % the gradient vector from the measurement of Q
        gradVm              % the gradient vector from the measurement of Vm
        gradVa              % the gradient vector from the measurement of Va
        numGrad             % the number of the gradient elements
        loss                % the sum-of-squares loss function
        lossChain           % the chain of the loss functions
        parChain            % the chain of the parameters
        
        kZero               % the ratio that we set nondiagonal elements to zero
        maxIter             % the maximum iteration in the gradient-based methods
        step                % the step length of the iterations
        stepMin             % the minimum step length
        stepMax             % the maximum step length
        stepChain           % the chain of the step length
        iter                % the current iteration step
        updateStepFreq      % the frequency to update the step length
        
        momentRatio         % the part we maintain from the past gradient
        momentRatioMax      % the maximum momentRatio
        momentRatioMin      % the minimum momentRatio
        vmvaWeight          % the additional weight to the vm and va
        isConverge          % if the iteration concerges
        isGB                % whether to iterate the GB part
        
        H                   % the Hessian matrix
        HP                  % the P part
        HQ                  % the Q part
        HVm                 % the Vm part
        HVa                 % the Va part
        
        J                   % the Jacobian matrix
        Topo                % the topology matrix
        Tvec                % the topology vector
        thsTopo             % the threshold of Topology
        
        lambda              % the damping ratio of LM algorithm
        lambdaMax           % the maximum value
        lambdaMin           % the minimum value
        lambdaChain         % the chain of lambda ratio, the first order ratio
        gradOrigin          % the original gradient in LM algorithm
        gradPast            % the gradient of the last step
        lossMin             % the theoretical minimum loss
        momentLoss          % the moment of loss
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
                cov = obj.A_FIM(obj.numFIM.index, obj.numFIM.index)\eye(sum(obj.numFIM.index));
                var = diag(cov);
            else
                cov = obj.A_FIM(obj.numFIM.index, obj.numFIM.index)\eye(sum(obj.numFIM.index));
                var = diag(cov);
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
                cov = cov - diag(diag(cov)) + diag(var);
                fprintf('We use the absolute value of the variance.\n');
            end
            
            obj.boundA.total = sqrt(var);
            obj.boundA.cov = cov;
            
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
        
        function obj = identifyTopo(obj)
            % This method identifies the topology by the voltage magnitudes
            % measurements. The initial topology may not be accurate.
            T = obj.data.G~=0;
            Vm = movmean(obj.data.Vm_noised, floor(obj.numSnap/20)+1, 2);
            C = corrcoef(Vm');
            C(isnan(C)) = 0;
            C(1,1) = 1;
            CT = corrcoef(obj.data.Vm');
            CT(isnan(CT)) = 0;
            CT(1,1) = 1;
        end
        
        function obj = approximateY(obj)
            % This method approximates the Y matrix by the measurements. We
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
%             Vm = obj.data.Vm_noised;
            Vm = obj.data.Vm;
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
            
            T = obj.data.G~=0;

            % approximate the parameter
            IP = obj.data.IP_noised;
            IQ = obj.data.IQ_noised;
            G = IP * obj.data.Vm_noised' / (obj.data.Vm_noised * obj.data.Vm_noised');
            B = - IQ * obj.data.Vm_noised' / (obj.data.Vm_noised * obj.data.Vm_noised');
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
            
            
%             obj.dataE.G = (G+G')/2;
%             obj.dataE.B = (B+B')/2;
            
%             obj.dataE.G = obj.data.G;
%             obj.dataE.B = obj.data.B;
        end
        
        function obj = iterateY(obj)
            % This method iterate Y matrix considering the measurement
            % error from both inputs and outputs.
            
            % We first assume a flat diagonal element setting
            W = ones(obj.numBus*2, 1);
            obj = optimizeY(obj, W);
        end
        
        function [obj, Gopt, Bopt] = optimizeY(obj, W)
            % This method use some convex optimization method and provide
            % the G and B matrix
            
            % control variables
            G = sdpvar(obj.numBus, obj.numBus);
            B = sdpvar(obj.numBus, obj.numBus);
            % anxillary variables
            Pres = sdpvar(obj.numBus, obj.numSnap);
            Qres = sdpvar(obj.numBus, obj.numSnap);
            
            % constraints
            constP = Pres == G * obj.data.Vm_noised - obj.data.IP_noised;
            constQ = Qres == - B * obj.data.Vm_noised - obj.data.IQ_noised;
            constG = sum(G) == zeros(1, obj.numBus);
            constB = sum(B) == zeros(1, obj.numBus);
            constraints = [constP; constQ; constG; constB];
            for i = 1:obj.numBus
                for j = i+1:obj.numBus
                    constraints = [constraints; G(i,j)<=0];
                    constraints = [constraints; B(i,j)>=0];
                end
            end
            
            % objective function
            objective = sum(W(1:obj.numBus)' * (Pres .* Pres)...
                + W(1+obj.numBus:end)' * (Qres .* Qres));
            options = sdpsettings('solver','gurobi');
            sol = optimize(constraints,objective,options);
            
            Gopt = value(G);
            Bopt = value(B);
        end
        
        function obj = initValue(obj)
            % This method provides the initial value (voltage angles?)
            
        end
        
        function obj = identifyOptNLP(obj)
            % This method simply use the nonlinar programming techique to
            % solve the maximum identification problem
            
            % This version we simply assume we have all the measurements
            % We should bound all the control variables and all the
            % anxillary variables
            
            % control variables
            G = sdpvar(obj.numBus, obj.numBus);
            B = sdpvar(obj.numBus, obj.numBus);
            Pest = sdpvar(obj.numBus, obj.numSnap);
            Qest = sdpvar(obj.numBus, obj.numSnap);
            Vm = sdpvar(obj.numBus, obj.numSnap);
            Va = sdpvar(obj.numBus, obj.numSnap);
            % anxillary variables
            e_P = sdpvar(obj.numBus, obj.numSnap);
            e_Q = sdpvar(obj.numBus, obj.numSnap);
            e_Vm = sdpvar(obj.numBus, obj.numSnap);
            e_Va = sdpvar(obj.numBus, obj.numSnap);
            Theta_ij = sdpvar(obj.numBus, obj.numBus, obj.numSnap);
            GBThetaP = sdpvar(obj.numBus, obj.numBus, obj.numSnap);
            GBThetaQ = sdpvar(obj.numBus, obj.numBus, obj.numSnap);
            % some constaints
            maxGB = 1000;
            maxNoise = 10;
            
            % constraints
            Constraints = [];
            % the power flow equation, P and Q injections
            for snap = 1:obj.numSnap
                Theta_ij(:,:,snap) = repmat(Va(:, snap), 1, obj.numBus) - repmat(Va(:, snap)', obj.numBus, 1);
                % G_ij\cos(\Theta_ij)+B_ij\sin(\Theta_ij)
                GBThetaP(:,:,snap) = G .* cos(Theta_ij(:,:,snap)) + B .* sin(Theta_ij(:,:,snap));
                % G_ij\sin(\Theta_ij)-B_ij\cos(\Theta_ij)
                GBThetaQ(:,:,snap) = G .* sin(Theta_ij(:,:,snap)) - B .* cos(Theta_ij(:,:,snap));
                Constraints = [Constraints; Pest(:, snap) == (GBThetaP(:,:,snap) * Vm(:, snap)) .* Vm(:, snap)];
                Constraints = [Constraints; Qest(:, snap) == (GBThetaQ(:,:,snap) * Vm(:, snap)) .* Vm(:, snap)];
            end
            % the anxillary variable constraints
            Constraints = [Constraints; Pest + e_P == obj.data.P_noised];
            Constraints = [Constraints; Qest + e_Q == obj.data.Q_noised];
            Constraints = [Constraints; Vm + e_Vm == obj.data.Vm_noised];
            Constraints = [Constraints; Va + e_Va == obj.data.Va_noised];
            % zero noise for reference bus
            Constraints = [Constraints; e_Va(1,:) == zeros(1, obj.numSnap)];
            Constraints = [Constraints; e_Vm(1,:) == zeros(1, obj.numSnap)];
            % the sum of G and B
            Constraints = [Constraints; sum(G) == zeros(1, obj.numBus)];
            Constraints = [Constraints; sum(B) == zeros(1, obj.numBus)];
            % bound all the variables
%             for i = 1:obj.numBus
%                 for j = i+1:obj.numBus
%                     Constraints = [Constraints; -maxGB <= G(i,j) <= 0];
%                     Constraints = [Constraints; 0 <= B(i,j) <= maxGB];
%                 end
%             end
            Constraints = [Constraints; -obj.sigma.P*ones(1, obj.numSnap)*maxNoise <= e_P <= obj.sigma.P*ones(1, obj.numSnap)*maxNoise];
            Constraints = [Constraints; -obj.sigma.Q*ones(1, obj.numSnap)*maxNoise <= e_Q <= obj.sigma.Q*ones(1, obj.numSnap)*maxNoise];
            Constraints = [Constraints; -obj.sigma.Vm*ones(1, obj.numSnap)*maxNoise <= e_Vm <= obj.sigma.Vm*ones(1, obj.numSnap)*maxNoise];
            Constraints = [Constraints; -obj.sigma.Va*ones(1, obj.numSnap)*maxNoise <= e_Va <= obj.sigma.Va*ones(1, obj.numSnap)*maxNoise];
            
            % assign the initial value
            assign(G, obj.dataE.G);
            assign(B, obj.dataE.B);
            assign(Vm, obj.data.Vm_noised);
            assign(Va, obj.data.Va_noised);
            
            % objective function
            objective = sum((obj.sigma.P.^-2)' * (e_P.*e_P) ...
                + (obj.sigma.Q.^-2)' * (e_Q.*e_Q)...
                + (obj.sigma.Vm(2:end).^-2)' * (e_Vm(2:end,:).*e_Vm(2:end,:))...
                + (obj.sigma.Va(2:end).^-2)' * (e_Va(2:end,:).*e_Va(2:end,:)));
            options = sdpsettings('solver','ipopt','ipopt.max_iter',3000);
            sol = optimize(Constraints,objective,options);
            
            Gopt = value(G);
            Bopt = value(B);
            Pestopt = value(Pest);
            Qestopt = value(Qest);
            Vmopt = value(Vm);
            Vaopt = value(Va);
            e_Popt = value(e_P);
            e_Qopt = value(e_Q);
            e_Vmopt = value(e_Vm);
            e_Vaopt = value(e_Va);
        end
        
        function obj = identifyOptGradient(obj)
            % This method uses gradient-based method to solve the nonconvex
            % optimization problem.
            % Hopefully we could implement some power system domain
            % knowledge into the process because we know the ground truth
            % value.
            obj.maxIter = 2000;
            obj.step = 1;
            obj.stepMax = 2;
            obj.stepMin = 0.0001;
            obj.momentRatio = 0.9;
            obj.updateStepFreq = 20;
            obj.vmvaWeight = 1;
            obj.momentRatioMax = 0.9;
            obj.momentRatioMin = 0.9;
            obj.kZero = 0.0005;
            
            % we first initialize data
            obj.dataO.G = obj.dataE.G;
            obj.dataO.B = obj.dataE.B;
            % note that we should replace the Vm ro Va data to some
            % initialized data if we do not have the measurement devices
            obj.dataO.Vm = obj.data.Vm;
%             obj.dataO.Va = obj.data.Va;
            obj.dataO.Vm(2:end, :) = bsxfun(@times, obj.data.Vm_noised(2:end, :), obj.isMeasure.Vm(2:end));
            obj.dataO.Vm(obj.dataO.Vm == 0) = 1;
            obj.dataO.Va = bsxfun(@times, obj.data.Va_noised, obj.isMeasure.Va);
            
            % begin the iteration loop
            % initialize the gradient numbers
            obj.numGrad.G = (obj.numBus - 1) * obj.numBus / 2; % exclude the diagonal elements
            obj.numGrad.B = (obj.numBus - 1) * obj.numBus / 2;
            obj.numGrad.Vm = obj.numSnap * (obj.numBus - 1); % exclude the source bus
            obj.numGrad.Va = obj.numSnap * (obj.numBus - 1);
            obj.numGrad.Sum = obj.numGrad.G + obj.numGrad.B + obj.numGrad.Vm + obj.numGrad.Va;
            obj.iter = 1;
            obj.gradChain = zeros(obj.numGrad.Sum, obj.maxIter);
            obj.lossChain = zeros(5, obj.maxIter);
            obj.parChain = zeros(obj.numGrad.Sum, obj.maxIter);
            obj.stepChain = zeros(1, obj.maxIter);
            
            obj.isConverge = false;
            while (obj.iter <= obj.maxIter && ~obj.isConverge)
                % collect the paramter vector
                obj = collectPar(obj);
                % build the gradient
                obj = buildGradient(obj);
                % implement the re-weight techique.
                obj = tuneGradient(obj);
                % update the chains
                try
                    obj.gradChain(:, obj.iter) = obj.grad * (1-obj.momentRatio) + obj.gradChain(:, obj.iter-1) * obj.momentRatio;
                catch
                    obj.gradChain(:, obj.iter) = obj.grad;
                end
                obj.lossChain(:, obj.iter) = [obj.loss.total; obj.loss.P; obj.loss.Q; obj.loss.Vm; obj.loss.Va];
                % update the parameters
                obj = updatePar(obj);
                % if converge
                if mod(obj.iter, obj.updateStepFreq) == 0 %obj.iter > 10
                    if (mean(obj.lossChain(1, obj.iter-9:obj.iter-5)) < mean(obj.lossChain(1, obj.iter-4:obj.iter)))
                        obj.step = max(obj.step / 2, obj.stepMin);
                        obj.momentRatio = min(obj.momentRatio + 0.1, obj.momentRatioMax);
                    elseif (((obj.lossChain(1, obj.iter) - obj.lossChain(1, obj.iter-1)) < 0) ...
                        && ((obj.lossChain(1, obj.iter-1) - obj.lossChain(1, obj.iter-2)) < 0)...
                        && ((obj.lossChain(1, obj.iter-2) - obj.lossChain(1, obj.iter-3)) < 0))
%                         && ((obj.lossChain(1, obj.iter-3) - obj.lossChain(1, obj.iter-4)) < 0))
                        obj.step = min(obj.step * 1.2, obj.stepMax);
                        obj.momentRatio = max(obj.momentRatio - 0.1, obj.momentRatioMin);
                    elseif (((obj.lossChain(1, obj.iter) - obj.lossChain(1, obj.iter-1)) *...
                        (obj.lossChain(1, obj.iter-1) - obj.lossChain(1, obj.iter-2)) < 0)...
                        && ((obj.lossChain(1, obj.iter-7) - obj.lossChain(1, obj.iter-8)) *...
                        (obj.lossChain(1, obj.iter-8) - obj.lossChain(1, obj.iter-9)) < 0))
                        obj.step = max(obj.step / 2, obj.stepMin);
                        obj.momentRatio = min(obj.momentRatio + 0.1, obj.momentRatioMax);
                    end
%                     if (mean(obj.lossChain(1, obj.iter-9:obj.iter-5)) < mean(obj.lossChain(1, obj.iter-4:obj.iter)))
%                         isConverge = true;
%                     end
                end
                obj.stepChain(obj.iter) = obj.step;
                obj.iter = obj.iter + 1;
            end
        end
        
        function obj = buildGradient(obj)
            % This method build the gradient of the squared loss function
            
            % Initialize the gradient matrix
            
            obj.grad = zeros(obj.numGrad.Sum, 1);
            obj.gradP = zeros(obj.numGrad.Sum, 1);
            obj.gradQ = zeros(obj.numGrad.Sum, 1);
            obj.gradVm = zeros(obj.numGrad.Sum, 1);
            obj.gradVa = zeros(obj.numGrad.Sum, 1);
            
            obj.loss.total = 0;
            obj.loss.P = 0;
            obj.loss.Q = 0;
            obj.loss.Vm = 0;
            obj.loss.Va = 0;
            
            for i = 1:obj.numSnap
                % calculate some basic parameters at present state
                Theta_ij = repmat(obj.dataO.Va(:, i), 1, obj.numBus) - repmat(obj.dataO.Va(:, i)', obj.numBus, 1);
                % G_ij\cos(\Theta_ij)+B_ij\sin(\Theta_ij)
                GBThetaP = obj.dataO.G .* cos(Theta_ij) + obj.dataO.B .* sin(Theta_ij);
                % G_ij\sin(\Theta_ij)-B_ij\cos(\Theta_ij)
                GBThetaQ = obj.dataO.G .* sin(Theta_ij) - obj.dataO.B .* cos(Theta_ij);
                % P estimate
                Pest = (GBThetaP * obj.dataO.Vm(:, i)) .* obj.dataO.Vm(:, i);
                % Q estimate
                Qest = (GBThetaQ * obj.dataO.Vm(:, i)) .* obj.dataO.Vm(:, i);
                
                % calculate the sub-vector of P of all buses
                for j = 1:obj.numBus
                    if obj.isMeasure.P(j)
                        obj = buildGradientP(obj, i, j, GBThetaP, GBThetaQ, Pest);
                    end
                end
                
                % calculate the sub-vector of Q of all buses
                for j = 1:obj.numBus
                    if obj.isMeasure.Q(j)
                        obj = buildGradientQ(obj, i, j, GBThetaP, GBThetaQ, Qest);
                    end
                end
                
                % calculate the sub-vector of Vm of all buses
                for j = 1:obj.numBus
                    if obj.isMeasure.Vm(j)
                        obj = buildGradientVm(obj, i, j);
                    end
                end
                
                % calculate the sub-vector of Va of all buses
                for j = 1:obj.numBus
                    if obj.isMeasure.Va(j)
                        obj = buildGradientVa(obj, i, j);
                    end
                end
            end
            
            % collect the gradients and the losses
            obj.grad = obj.gradP + obj.gradQ + obj.gradVm + obj.gradVa;
            obj.loss.total = obj.loss.P + obj.loss.Q + obj.loss.Vm + obj.loss.Va;
        end
        
        function obj = buildGradientP(obj , snap, bus, GBThetaP, GBThetaQ, Pest)
            % This method builds the gradient from the measurement of P
            
            theta_ij = obj.dataO.Va(bus, snap) - obj.dataO.Va(:, snap);
            g = zeros(obj.numGrad.Sum, 1);
            
            % G matrix
            H_G = zeros(obj.numBus, obj.numBus);
            H_G(bus, :) = obj.dataO.Vm(bus, snap) * obj.dataO.Vm(:, snap)' .* cos(theta_ij');
            H_G(bus, :) = H_G(bus, :) - obj.dataO.Vm(bus, snap)^2; % the equivilance of diagonal elements
            h_G = obj.matToColDE(H_G);
            g(1:obj.numGrad.G) = h_G;
            
            % B matrix
            H_B = zeros(obj.numBus, obj.numBus);
            H_B(bus, :) = obj.dataO.Vm(bus, snap) * obj.dataO.Vm(:, snap)' .* sin(theta_ij');
            h_B = obj.matToColDE(H_B);
            g(obj.numGrad.G+1:obj.numGrad.G+obj.numGrad.B) = h_B;
            
            % Vm
            % the first order term of other Vm
            H_Vm = zeros(obj.numBus, obj.numSnap);
            h_Vm = obj.dataO.Vm(bus, snap) * GBThetaP(:, bus);
            % the second order term of Vm(bus)
            h_Vm(bus) = 2*obj.dataO.Vm(bus, snap) * GBThetaP(bus, bus);
            % the first order term of Vm(bus)
            fOrderVm = obj.dataO.Vm(:, snap) .* GBThetaP(:, bus);
            fOrderVm(bus) = 0;
            h_Vm(bus) = h_Vm(bus) + sum(fOrderVm);
            H_Vm(:, snap) = h_Vm;
            % remove the source bus whose magnitude is not the state variable
            H_Vm(1, :) = []; 
            h_VmLarge = reshape(H_Vm', [], 1);
            g(obj.numGrad.G+obj.numGrad.B+1:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) = h_VmLarge;
            
            % Va
            H_Va = zeros(obj.numBus, obj.numSnap);
            h_Va = obj.dataO.Vm(bus, snap) * obj.dataO.Vm(:, snap) .* GBThetaQ(:, bus);
            h_Va(bus) = - obj.dataO.Vm(bus, snap)^2 * obj.dataO.B(bus, bus)...
                       - obj.data.Q_noised(bus, snap); 
%             h_Va(bus) = h_Va(bus)-sum(GBThetaQ(bus, :) * obj.dataO.Vm(:, snap) * obj.dataO.Vm(bus, snap));
            H_Va(:, snap) = h_Va;
%             assert (H_Va(bus, snap) > 0)
            % remove the source bus whose magnitude is not the state variable
            H_Va(1, :) = []; 
            h_VaLarge = reshape(H_Va', [], 1);
            g(obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm+1:end) = h_VaLarge;
            
            % build GradientP and loss.P
            lossThis = (Pest(bus) - obj.data.P_noised(bus, snap));
            obj.loss.P = obj.loss.P + lossThis^2 * obj.sigma.P(bus).^(-2);
            gradPThis = obj.sigma.P(bus).^(-2) * lossThis * g;
            obj.gradP = obj.gradP + gradPThis;
        end
        
        function obj = buildGradientQ(obj , snap, bus, GBThetaP, GBThetaQ, Qest)
            % This method builds the gradient from the measurement of Q
            
            theta_ij = obj.dataO.Va(bus, snap) - obj.dataO.Va(:, snap);
            g = zeros(obj.numGrad.Sum, 1);
            
            % G matrix
            H_G = zeros(obj.numBus, obj.numBus);
            H_G(bus, :) = obj.dataO.Vm(bus, snap) * obj.dataO.Vm(:, snap)' .* sin(theta_ij');
            h_G = obj.matToColDE(H_G);
            g(1:obj.numGrad.G) = h_G;
            
            % B matrix
            H_B = zeros(obj.numBus, obj.numBus);
            H_B(bus, :) = - obj.dataO.Vm(bus, snap) * obj.dataO.Vm(:, snap)' .* cos(theta_ij');
            H_B(bus, :) = H_B(bus, :) + obj.dataO.Vm(bus, snap)^2; % the equivilance of diagonal elements
            h_B = obj.matToColDE(H_B);
            g(obj.numGrad.G+1:obj.numGrad.G+obj.numGrad.B) = h_B;
            
            % Vm
            % the first order term of other Vm
            H_Vm = zeros(obj.numBus, obj.numSnap);
            h_Vm = obj.dataO.Vm(bus, snap) * GBThetaQ(:, bus);
            % the second order term of Vm(bus)
            h_Vm(bus) = 2*obj.dataO.Vm(bus, snap) * GBThetaQ(bus, bus);
            % the first order term of Vm(bus)
            fOrderVm = obj.dataO.Vm(:, snap) .* GBThetaQ(:, bus);
            fOrderVm(bus) = 0;
            h_Vm(bus) = h_Vm(bus) + sum(fOrderVm);
            H_Vm(:, snap) = h_Vm;
            % remove the source bus whose magnitude is not the state variable
            H_Vm(1, :) = []; 
            h_VmLarge = reshape(H_Vm', [], 1);
            g(obj.numGrad.G+obj.numGrad.B+1:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) = h_VmLarge;
            
            % Va
            H_Va = zeros(obj.numBus, obj.numSnap);
            h_Va = - obj.dataO.Vm(bus, snap) * obj.dataO.Vm(:, snap) .* GBThetaP(:, bus);
            h_Va(bus) = - obj.dataO.Vm(bus, snap)^2 * obj.dataO.G(bus, bus) ...
                        + obj.data.P_noised(bus, snap);
%             h_Va(bus) = h_Va(bus)+sum(GBThetaP(bus, :) * obj.dataO.Vm(:, snap) * obj.dataO.Vm(bus, snap));
            H_Va(:, snap) = h_Va;
            % remove the source bus whose magnitude is not the state variable
            H_Va(1, :) = []; 
            h_VaLarge = reshape(H_Va', [], 1);
            g(obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm+1:end) = h_VaLarge;
            
            % build GradientQ and lossQ
            lossThis = (Qest(bus) - obj.data.Q_noised(bus, snap));
            obj.loss.Q = obj.loss.Q + lossThis^2 * obj.sigma.Q(bus).^(-2);
            gradQThis = obj.sigma.Q(bus).^(-2) * lossThis * g;
            obj.gradQ = obj.gradQ + gradQThis;
        end
        
        function obj = buildGradientVm(obj, snap, bus)
            % This method builds the gradient from the measurement of Vm
            g = zeros(obj.numGrad.Sum, 1);
            H_Vm = zeros(obj.numBus, obj.numSnap);
            H_Vm(bus, snap) = 1 / (obj.sigma.Vm(bus)^2);
            % remove the source bus whose magnitude is not the state variable
            H_Vm(1, :) = []; 
            h_VmLarge = reshape(H_Vm', [], 1);
            g(obj.numGrad.G+obj.numGrad.B+1:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) = h_VmLarge;
            
            % build GradientVm and lossVm
            lossThis = (obj.dataO.Vm(bus, snap) - obj.data.Vm_noised(bus, snap));
            obj.loss.Vm = obj.loss.Vm + lossThis^2 * obj.sigma.Vm(bus).^(-2);
            gradVmThis = lossThis * g;
            obj.gradVm = obj.gradVm + gradVmThis;
        end
        
        function obj = buildGradientVa(obj, snap, bus)
            % This method builds the gradient from the measurement of Va
            g = zeros(obj.numGrad.Sum, 1);
            H_Va = zeros(obj.numBus, obj.numSnap);
            H_Va(bus, snap) = 1 / (obj.sigma.Va(bus)^2);
            % remove the source bus whose angle is not the state variable
            H_Va(1, :) = []; 
            h_VaLarge = reshape(H_Va', [], 1);
            g(obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm+1:end) = h_VaLarge;
            
            % build GradientVa and lossVa
            lossThis = (obj.dataO.Va(bus, snap) - obj.data.Va_noised(bus, snap));
            obj.loss.Va = obj.loss.Va + lossThis^2 * obj.sigma.Va(bus).^(-2);
            gradVaThis = lossThis * g;
            obj.gradVa = obj.gradVa + gradVaThis;
        end
        
        function obj = tuneGradient(obj)
            % This method tunes the gradient according to the weights. In
            % this version we treat P and Q together.
            
            % The weight of the initial gradient
            wg.P.G = mean(abs(obj.gradP(1:obj.numGrad.G)));
            wg.P.B = mean(abs(obj.gradP(1+obj.numGrad.G:obj.numGrad.G+obj.numGrad.B)));
            wg.P.Vm = mean(abs(obj.gradP(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm)));
            wg.P.Va = mean(abs(obj.gradP(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end)));
            wg.Q.G = mean(abs(obj.gradQ(1:obj.numGrad.G)));
            wg.Q.B = mean(abs(obj.gradQ(1+obj.numGrad.G:obj.numGrad.G+obj.numGrad.B)));
            wg.Q.Vm = mean(abs(obj.gradQ(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm)));
            wg.Q.Va = mean(abs(obj.gradQ(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end)));
            
            wg.PQ_GB = mean([wg.P.G wg.P.B wg.Q.G  wg.Q.B]);
            wg.PQ_Vm = mean([wg.P.Vm wg.Q.Vm]);
            wg.PQ_Va = mean([wg.P.Va wg.Q.Va]);
            wg.Vm_Vm = mean(abs(obj.gradVm(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm)));
            wg.Va_Va = mean(abs(obj.gradVa(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end)));
            
            % The weight of the CRLB
            wb.GB = mean(obj.boundA.total(1:obj.numFIM.G+obj.numFIM.B));
            wb.Vm = mean(abs(obj.boundA.total(1+obj.numFIM.G+obj.numFIM.B:obj.numFIM.G+obj.numFIM.B+obj.numFIM.Vm)));
            wb.Va = mean(abs(obj.boundA.total(1+obj.numFIM.G+obj.numFIM.B+obj.numFIM.Vm:end)));
            % The weight of the loss function
            wl.total = sqrt(obj.loss.P + obj.loss.Q) + sqrt(obj.loss.Vm) + sqrt(obj.loss.Va);
            wl.PQ = sqrt(obj.loss.P + obj.loss.Q) / wl.total;
            wl.Vm = sqrt(obj.loss.Vm) * obj.vmvaWeight / wl.total; % one P measurements related to multiple Vm and Va, we should correct this.  * 2 * obj.numBus  
            wl.Va = sqrt(obj.loss.Va) * obj.vmvaWeight * 3 / wl.total; % * 2 * obj.numBus; the number five is the average degree of a distribution network times two  * 5
            % The conditional weights
            wl.Vm_PQ = wl.PQ / (wl.PQ + wl.Vm + 1e-9);
            wl.Vm_Vm = wl.Vm / (wl.PQ + wl.Vm + 1e-9);
            wl.Va_PQ = wl.PQ / (wl.PQ + wl.Va + 1e-9);
            wl.Va_Va = wl.Va / (wl.PQ + wl.Va + 1e-9);
            
            % tune the gradient vector. 
            % normalize the gradient weight by P Q Vm Va
            obj.gradP(1:obj.numGrad.G+obj.numGrad.B) = ...
                obj.gradP(1:obj.numGrad.G+obj.numGrad.B) / (wg.PQ_GB + 1e-9);
            obj.gradP(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) = ...
                obj.gradP(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) / (wg.PQ_Vm + 1e-9);
            obj.gradP(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) = ...
                obj.gradP(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) / (wg.PQ_Va + 1e-9);
            
            obj.gradQ(1:obj.numGrad.G+obj.numGrad.B) = ...
                obj.gradQ(1:obj.numGrad.G+obj.numGrad.B) / (wg.PQ_GB + 1e-9);
            obj.gradQ(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) = ...
                obj.gradQ(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) / (wg.PQ_Vm + 1e-9);
            obj.gradQ(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) = ...
                obj.gradQ(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) / (wg.PQ_Va + 1e-9);
            
            obj.gradVm(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) = ...
                obj.gradVm(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) / (wg.Vm_Vm + 1e-9);
            obj.gradVa(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) = ...
                obj.gradVa(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) / (wg.Va_Va + 1e-9);
            
            % we use the weight of the approximated CRLB and the weight
            % from the loss function
            
            obj.gradP(1:obj.numGrad.G+obj.numGrad.B) = ...
                obj.gradP(1:obj.numGrad.G+obj.numGrad.B) * wb.GB;
            obj.gradP(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) = ...
                obj.gradP(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) * wb.Vm * wl.Vm_PQ;
            obj.gradP(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) = ...
                obj.gradP(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) * wb.Va * wl.Va_PQ;
            
            obj.gradQ(1:obj.numGrad.G+obj.numGrad.B) = ...
                obj.gradQ(1:obj.numGrad.G+obj.numGrad.B) * wb.GB;
            obj.gradQ(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) = ...
                obj.gradQ(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) * wb.Vm * wl.Vm_PQ;
            obj.gradQ(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) = ...
                obj.gradQ(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) * wb.Va * wl.Va_PQ;
            
            obj.gradVm(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) = ...
                obj.gradVm(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) * wb.Vm * wl.Vm_Vm;
            obj.gradVa(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) = ...
                obj.gradVa(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) * wb.Va * wl.Va_Va;
            
            % collect the tuned gardient
            obj.grad = obj.gradP + obj.gradQ + obj.gradVm + obj.gradVa;
        end
        
        function obj = tuneGradientPQ(obj)
            % This method tunes the gradient according to the weights. In
            % this version we treat P and Q independently.
            
            % The weight of the initial gradient
            wg.P.G = mean(abs(obj.gradP(1:obj.numGrad.G)));
            wg.P.B = mean(abs(obj.gradP(1+obj.numGrad.G:obj.numGrad.G+obj.numGrad.B)));
            wg.P.Vm = mean(abs(obj.gradP(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm)));
            wg.P.Va = mean(abs(obj.gradP(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end)));
            wg.Q.G = mean(abs(obj.gradQ(1:obj.numGrad.G)));
            wg.Q.B = mean(abs(obj.gradQ(1+obj.numGrad.G:obj.numGrad.G+obj.numGrad.B)));
            wg.Q.Vm = mean(abs(obj.gradQ(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm)));
            wg.Q.Va = mean(abs(obj.gradQ(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end)));
            wg.Vm = mean(abs(obj.gradVm(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm)));
            wg.Va = mean(abs(obj.gradVa(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end)));
            % The weight of the CRLB
            wb.G = mean(obj.boundA.total(1:obj.numFIM.G));
            wb.B = mean(obj.boundA.total(1+obj.numFIM.G:obj.numFIM.G+obj.numFIM.B));
            wb.Vm = mean(abs(obj.boundA.total(1+obj.numFIM.G+obj.numFIM.B:obj.numFIM.G+obj.numFIM.B+obj.numFIM.Vm)));
            wb.Va = mean(abs(obj.boundA.total(1+obj.numFIM.G+obj.numFIM.B+obj.numFIM.Vm:end)));
            % The weight of the loss function
            wl.total = sqrt(obj.loss.P) + sqrt(obj.loss.Q) + sqrt(obj.loss.Vm) + sqrt(obj.loss.Va);
            wl.P = sqrt(obj.loss.P) / wl.total;
            wl.Q = sqrt(obj.loss.Q) / wl.total;
            wl.Vm = sqrt(obj.loss.Vm) * obj.vmvaWeight / wl.total; % one P measurements related to multiple Vm and Va, we should correct this.  * 2 * obj.numBus  
            wl.Va = sqrt(obj.loss.Va) * obj.vmvaWeight / wl.total; % * 2 * obj.numBus; the number five is the average degree of a distribution network times two  * 5
            % The conditional weights
            wl.GP = wl.P / (wl.P + wl.Q + 1e-9);
            wl.GQ = wl.Q / (wl.P + wl.Q + 1e-9);
            wl.BP = wl.P / (wl.P + wl.Q + 1e-9);
            wl.BQ = wl.Q / (wl.P + wl.Q + 1e-9);
            wl.VmP = wl.P / (wl.P + wl.Q + wl.Vm + 1e-9);
            wl.VmQ = wl.Q / (wl.P + wl.Q + wl.Vm + 1e-9);
            wl.VmVm = wl.Vm / (wl.P + wl.Q + wl.Vm + 1e-9);
            wl.VaP = wl.P / (wl.P + wl.Q + wl.Va + 1e-9);
            wl.VaQ = wl.Q / (wl.P + wl.Q + wl.Va + 1e-9);
            wl.VaVa = wl.Va / (wl.P + wl.Q + wl.Va + 1e-9);
            
            % tune the gradient vector. 
            % normalize the gradient weight by P Q Vm Va
            obj.gradP(1:obj.numGrad.G) = ...
                obj.gradP(1:obj.numGrad.G) / (wg.P.G + 1e-9);
            obj.gradP(1+obj.numGrad.G:obj.numGrad.G+obj.numGrad.B) = ...
                obj.gradP(1+obj.numGrad.G:obj.numGrad.G+obj.numGrad.B) / (wg.P.B + 1e-9);
            obj.gradP(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) = ...
                obj.gradP(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) / (wg.P.Vm + 1e-9);
            obj.gradP(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) = ...
                obj.gradP(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) / (wg.P.Va + 1e-9);
            
            obj.gradQ(1:obj.numGrad.G) = ...
                obj.gradQ(1:obj.numGrad.G) / (wg.Q.G + 1e-9);
            obj.gradQ(1+obj.numGrad.G:obj.numGrad.G+obj.numGrad.B) = ...
                obj.gradQ(1+obj.numGrad.G:obj.numGrad.G+obj.numGrad.B) / (wg.Q.B + 1e-9);
            obj.gradQ(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) = ...
                obj.gradQ(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) / (wg.Q.Vm + 1e-9);
            obj.gradQ(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) = ...
                obj.gradQ(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) / (wg.Q.Va + 1e-9);
            
            obj.gradVm(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) = ...
                obj.gradVm(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) / (wg.Vm + 1e-9);
            obj.gradVa(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) = ...
                obj.gradVa(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) / (wg.Va + 1e-9);
            
            % we use the weight of the approximated CRLB and the weight
            % from the loss function
            
            obj.gradP(1:obj.numGrad.G) = ...
                obj.gradP(1:obj.numGrad.G) * wb.G * wl.GP;
            obj.gradP(1+obj.numGrad.G:obj.numGrad.G+obj.numGrad.B) = ...
                obj.gradP(1+obj.numGrad.G:obj.numGrad.G+obj.numGrad.B) * wb.B * wl.BP;
            obj.gradP(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) = ...
                obj.gradP(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) * wb.Vm * wl.VmP;
            obj.gradP(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) = ...
                obj.gradP(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) * wb.Va * wl.VaP;
            
            obj.gradQ(1:obj.numGrad.G) = ...
                obj.gradQ(1:obj.numGrad.G) * wb.G * wl.GQ;
            obj.gradQ(1+obj.numGrad.G:obj.numGrad.G+obj.numGrad.B) = ...
                obj.gradQ(1+obj.numGrad.G:obj.numGrad.G+obj.numGrad.B) * wb.B * wl.BQ;
            obj.gradQ(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) = ...
                obj.gradQ(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) * wb.Vm * wl.VmQ;
            obj.gradQ(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) = ...
                obj.gradQ(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) * wb.Va * wl.VaQ;
            
            obj.gradVm(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) = ...
                obj.gradVm(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) * wb.Vm * wl.VmVm;
            obj.gradVa(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) = ...
                obj.gradVa(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) * wb.Va * wl.VaVa;
            
            % collect the tuned gardient
            obj.grad = obj.gradP + obj.gradQ + obj.gradVm + obj.gradVa;
        end
        
        function obj = collectPar(obj)
            % This method formulates the parameter vector
            par = zeros(obj.numGrad.Sum, 1);
            par(1:obj.numGrad.G) = obj.matOfColDE(obj.dataO.G);
            par(1+obj.numGrad.G:obj.numGrad.G+obj.numGrad.B) = obj.matOfColDE(obj.dataO.B);
            par(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) = ...
                reshape(obj.dataO.Vm(2:end,:)', [], 1); % we assume the value of the source bus is already known
            par(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end) = ...
                reshape(obj.dataO.Va(2:end,:)', [], 1);
            obj.parChain(:, obj.iter) = par;
        end
        
        function obj = updatePar(obj)
            % This method updates the parameters in the iteration process
            delta = obj.step * obj.gradChain(:, obj.iter);
            par = obj.parChain(:, obj.iter) - delta;
            % gather the par values
            G = par(1:obj.numGrad.G);
            B = par(1+obj.numGrad.G:obj.numGrad.G+obj.numGrad.B);
%             G(G>0) = 0;
%             B(B<0) = 0;
%             G(B==0) = 0; % we do not use it because it will cause sudden change
%             B(G==0) = 0;
            % we first do not assume any topologies, then we would add some
            % topology iteration techiques.
            obj.dataO.G = obj.colToMatDE(G, obj.numBus);
            obj.dataO.B = obj.colToMatDE(B, obj.numBus);
%             if mod(obj.iter, obj.updateVmVaFreq) == 0 
            Vm = par(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm);
            Va = par(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end);
            obj.dataO.Vm(2:end, :) = reshape(Vm, [], obj.numBus-1)'; % exclude the source bus
            obj.dataO.Va(2:end, :) = reshape(Va, [], obj.numBus-1)'; % exclude the source bus
%             end
        end
        
        function obj = identifyOptNewton(obj)
            % This method uses Newton method to update the parameters
            obj.maxIter = 500;
            obj.thsTopo = 0.01;
            obj.Topo = true(obj.numBus, obj.numBus);
            obj.Tvec = logical(obj.matOfColDE(obj.Topo));
            
            % we first initialize data
            obj.dataO.G = obj.dataE.G;
            obj.dataO.B = obj.dataE.B;
            % note that we should replace the Vm ro Va data to some
            % initialized data if we do not have the measurement devices
            obj.dataO.Vm = obj.data.Vm;
%             obj.dataO.Va = obj.data.Va;
            obj.dataO.Vm(2:end, :) = bsxfun(@times, obj.data.Vm_noised(2:end, :), obj.isMeasure.Vm(2:end));
            obj.dataO.Vm(obj.dataO.Vm == 0) = 1;
            obj.dataO.Va = bsxfun(@times, obj.data.Va_noised, obj.isMeasure.Va);
            
%             obj.dataO.P = obj.data.P_noised;
%             obj.dataO.Q = obj.data.Q_noised;
            obj.dataO.P = obj.data.P_noised;
            obj.dataO.Q = obj.data.Q_noised;
            
            % begin the iteration loop
            % initialize the gradient numbers
            obj.numGrad.G = (obj.numBus - 1) * obj.numBus / 2; % exclude the diagonal elements
            obj.numGrad.B = (obj.numBus - 1) * obj.numBus / 2;
            obj.numGrad.Vm = obj.numSnap * (obj.numBus - 1); % exclude the source bus
            obj.numGrad.Va = obj.numSnap * (obj.numBus - 1);
            obj.numGrad.Sum = obj.numGrad.G + obj.numGrad.B + obj.numGrad.Vm + obj.numGrad.Va;
            obj.iter = 1;
            obj.lossChain = zeros(5, obj.maxIter);
            obj.parChain = zeros(obj.numGrad.Sum, obj.maxIter);
            obj.isGB = false;
            
            obj.isConverge = 0;
            while (obj.iter <= obj.maxIter && obj.isConverge <= 2)
                % update Va by power flow calculation
%                 obj = updateParPF(obj);
                % collect the paramter vector
                obj = collectPar(obj);
                % build the Jacobian
                obj = buildJacobian(obj);
                % build the Hessian
                obj = buildHessian(obj);
                obj.lossChain(:, obj.iter) = ...
                    [obj.loss.total; obj.loss.P; obj.loss.Q; ...
                    obj.loss.Vm; obj.loss.Va];
                % update the parameters
                obj = updateParNewtonGV(obj);
%                 if (mod(obj.iter-1, 2) ~= 0)
%                     obj = updateParNewtonGV(obj, true); % update GB
%                 else
% %                     obj = updateParNewtonGV(obj, true); % update GB
%                     obj = updateParNewtonGV(obj, false);
%                 end
                obj.iter = obj.iter + 1;
            end
        end
        
        function obj = identifyOptLMPower(obj)
            % This function identify the topology and the parameters using
            % the LM-based strategy and the knowledge of power flow
            % equations
            obj.lambda = 1e3; % the proportion of first order gradient
            obj.lambdaMin = 1e-3;
            obj.lambdaMax = 1e3;
            
            obj.step = 1;
            obj.stepMin = 1e-4;
            obj.stepMax = 2;
            
            obj.maxIter = 3000;
            obj.thsTopo = 0.01;
            obj.Topo = true(obj.numBus, obj.numBus);
            obj.Tvec = logical(obj.matOfColDE(obj.Topo));
            obj.vmvaWeight = 1;
            
            obj.momentRatio = 0.9;
            
            % we first initialize the data
            obj.dataO.G = obj.dataE.G;
            obj.dataO.B = obj.dataE.B;
            
            obj.dataO.Vm = obj.data.Vm;
%             obj.dataO.Va = obj.data.Va;
            obj.dataO.Vm(2:end, :) = bsxfun(@times, obj.data.Vm_noised(2:end, :), obj.isMeasure.Vm(2:end));
            obj.dataO.Vm(obj.dataO.Vm == 0) = 1;
            obj.dataO.Va = bsxfun(@times, obj.data.Va_noised, obj.isMeasure.Va);
            
            obj.dataO.P = obj.data.P_noised;
            obj.dataO.Q = obj.data.Q_noised;
            
            % initialize the gradient numbers
            obj.numGrad.G = (obj.numBus - 1) * obj.numBus / 2; % exclude the diagonal elements
            obj.numGrad.B = (obj.numBus - 1) * obj.numBus / 2;
            obj.numGrad.Vm = obj.numSnap * (obj.numBus - 1); % exclude the source bus
            obj.numGrad.Va = obj.numSnap * (obj.numBus - 1);
            obj.numGrad.Sum = obj.numGrad.G + obj.numGrad.B + obj.numGrad.Vm + obj.numGrad.Va;
            
            % evaluate the minimum loss
            obj = evaluateLossMin(obj);
            
            % begin the iteration loop
            obj.iter = 1;
            obj.lossChain = zeros(5, obj.maxIter);
            obj.parChain = zeros(obj.numGrad.Sum, obj.maxIter);
            obj.stepChain = zeros(1, obj.maxIter);
            obj.lambdaChain = zeros(1, obj.maxIter);
            obj.isGB = false;
            obj.isConverge = 0;
            while (obj.iter <= obj.maxIter && obj.isConverge <= 2)
                % collect the parameter vector
                obj = collectPar(obj);
                % build the Hessian
                obj = buildHessian(obj);
                obj.lossChain(:, obj.iter) = ...
                    [obj.loss.total; obj.loss.P; obj.loss.Q; ...
                    obj.loss.Vm; obj.loss.Va];
                % implement the re-weight techique.
                obj.gradOrigin = obj.grad;
                obj = tuneGradient(obj);
                % update the parameters
                obj = updateParLMPower(obj);
                obj.iter = obj.iter + 1;
            end
            obj.lossChain(:, obj.iter:end) = [];
            obj.parChain(:, obj.iter:end) = [];
            obj.stepChain(:, obj.iter:end) = [];
            obj.lambdaChain(:, obj.iter:end) = [];
        end
        
        function obj = evaluateLossMin(obj)
            % This method evalutes the minimum loss
            obj.lossMin = obj.numSnap *...
                sum([obj.isMeasure.P;obj.isMeasure.Q;obj.isMeasure.Vm;obj.isMeasure.Va]);
        end
        
        function obj = updateParLMPower(obj)
            % This method updates the parameters using LM strategies and
            % iteratively updates the GB and VmVa
            % calculate the modified Hessian
            
            % we update all the parameters together
            id = 1:obj.numGrad.Sum;
            try
                moment = obj.gradPast * obj.momentRatio + obj.grad * (1-obj.momentRatio);
            catch
                moment = obj.grad;
            end
            obj.gradPast = moment;
            delta1 = obj.step * moment(id);
%             delta1 = moment(id);
            delta2 = obj.H(id, id) \ obj.gradOrigin(id);
            delta = delta1 * obj.lambda/(1+obj.lambda) + delta2 * 1/(1+obj.lambda);
%             delta = delta * obj.step;

%             H1 = moment * obj.lambda + obj.H;
%             delta = H1(id, id) \ obj.gradOrigin(id);
            par = zeros(obj.numGrad.Sum, 1);
            
            
            
            par(id) = obj.parChain(id, obj.iter) - delta(id);
            G = par(1:obj.numGrad.G);
            B = par(1+obj.numGrad.G:obj.numGrad.G+obj.numGrad.B);
%             G(G>0) = 0;
%             B(B<0) = 0;
            obj.dataO.G = obj.colToMatDE(G, obj.numBus);
            obj.dataO.B = obj.colToMatDE(B, obj.numBus);
            
            Vm = par(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm);
            Va = par(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end);
            obj.dataO.Vm(2:end, :) = reshape(Vm, [], obj.numBus-1)'; % exclude the source bus
            obj.dataO.Va(2:end, :) = reshape(Va, [], obj.numBus-1)'; % exclude the source bus
            
%             if obj.isGB
% %                 delta = zeros(obj.numGrad.G+obj.numGrad.B, 1);
% %                 id = [obj.Tvec; obj.Tvec];
%                 id = 1:obj.numGrad.G+obj.numGrad.B;
%                 delta = H1(id, id) \ obj.gradOrigin(id);
%                 % converge or not
%                 if (max(abs(delta)) < 1e-5)
%                     obj.isConverge = obj.isConverge + 1;
%                 else
%                     obj.isConverge = 0;
%                 end
%                 obj.isGB = false;
%                 % update the parameters
%                 par = zeros(obj.numGrad.Sum, 1);
%                 par(id) = obj.parChain(id, obj.iter) - delta(id);
%                 G = par(1:obj.numGrad.G);
%                 B = par(1+obj.numGrad.G:obj.numGrad.G+obj.numGrad.B);
%                 G(G>0) = 0;
%                 B(B<0) = 0;
%                 obj.dataO.G = obj.colToMatDE(G, obj.numBus);
%                 obj.dataO.B = obj.colToMatDE(B, obj.numBus);
%             else
%                 id = 1+obj.numGrad.B+obj.numGrad.G:obj.numGrad.Sum;
%                 delta = H1(id, id) \ obj.gradOrigin(id);
%                 % converge or not
%                 if (max(abs(delta)) < 1e-9)
%                     obj.isConverge = obj.isConverge + 1;
%                 else
%                     obj.isConverge = 0;
%                 end
%                 obj.isGB = true;
%                 % update the paramters
%                 par = zeros(obj.numGrad.Sum, 1);
%                 par(id) = obj.parChain(id, obj.iter) - delta;
%                 Vm = par(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm);
%                 Va = par(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end);
%                 obj.dataO.Vm(2:end, :) = reshape(Vm, [], obj.numBus-1)'; % exclude the source bus
%                 obj.dataO.Va(2:end, :) = reshape(Va, [], obj.numBus-1)'; % exclude the source bus
%             end
            % update the lambda and the step length
            try % I think we should use PD control here, currently it is only D control
                if obj.lossChain(1, obj.iter) < obj.lossChain(1, obj.iter-1)
%                     ratio = (obj.lossChain(1, obj.iter-1) / obj.lossChain(1, obj.iter));
                    obj.lambda = max(obj.lambda / 1.1, obj.lambdaMin);
                    obj.step = min(obj.step * 1.1, obj.stepMax);
                else
%                     ratio = (obj.lossChain(1, obj.iter) / obj.lossChain(1, obj.iter-1))^2;
                    ratio = log10(max(obj.loss.total, obj.lossMin * 10) / obj.lossMin);
                    dRatio = 10/ratio;
                    obj.lambda = min(obj.lambda * (1.1+dRatio), obj.lambdaMax);
%                     dRatio = 0.5;
                    obj.lambda = min(obj.lambda * (1.1+dRatio), obj.lambdaMax);
                    obj.step = max(obj.step / (1.1+dRatio), obj.stepMin);
                end
            catch
            end
%             if obj.iter > 1 % I think we should use PD control here, currently it is only D control
%                 if obj.loss.total < obj.momentLoss
%                     obj.lambda = max(obj.lambda / 1.2, obj.lambdaMin);
%                     obj.step = min(obj.step * 1.2, obj.stepMax);
%                 else
%                     obj.lambda = min(obj.lambda * 2, obj.lambdaMax);
%                     obj.step = max(obj.step / 2, obj.stepMin);
%                 end
%                 obj.momentLoss = obj.momentLoss * 0.1 + obj.loss.total * (1-0.1);
%             else
%                 obj.momentLoss = obj.loss.total;
%             end
            obj.lambdaChain(obj.iter) = obj.lambda;
            obj.stepChain(obj.iter) = obj.step;
%             obj.lambdaMax = log10(max(obj.loss.total, obj.lossMin * 10) / obj.lossMin) * 1000;
            % converge or not
            if obj.loss.total < obj.lossMin
                obj.isConverge = 3;
            end
        end
        
        function obj = buildJacobian(obj)
            % This method build the Jacobian matrix
            obj.J = zeros(2 * obj.numBus * obj.numSnap, obj.numGrad.G + obj.numGrad.B + obj.numGrad.Va);
            numP = obj.numBus * obj.numSnap;
            for snap = 1:obj.numSnap
                % calculate some basic parameters at present state
                Theta_ij = repmat(obj.dataO.Va(:, snap), 1, obj.numBus) - repmat(obj.dataO.Va(:, snap)', obj.numBus, 1);
                % G_ij\cos(\Theta_ij)+B_ij\sin(\Theta_ij)
                GBThetaP = obj.dataO.G .* cos(Theta_ij) + obj.dataO.B .* sin(Theta_ij);
                % G_ij\sin(\Theta_ij)-B_ij\cos(\Theta_ij)
                GBThetaQ = obj.dataO.G .* sin(Theta_ij) - obj.dataO.B .* cos(Theta_ij);
                for bus = 1:obj.numBus
                    theta_ij = obj.dataO.Va(bus, snap) - obj.dataO.Va(:, snap);
                    % P
                    hP = zeros(obj.numGrad.G + obj.numGrad.B + obj.numGrad.Va, 1);
                    % G matrix
                    H_G = zeros(obj.numBus, obj.numBus);
                    H_G(bus, :) = obj.dataO.Vm(bus, snap) * obj.dataO.Vm(:, snap)' .* cos(theta_ij');
                    H_G(bus, :) = H_G(bus, :) - obj.dataO.Vm(bus, snap)^2; % the equivilance of diagonal elements
                    h_G = obj.matToColDE(H_G);
                    hP(1:obj.numGrad.G) = h_G;

                    % B matrix
                    H_B = zeros(obj.numBus, obj.numBus);
                    H_B(bus, :) = obj.dataO.Vm(bus, snap) * obj.dataO.Vm(:, snap)' .* sin(theta_ij');
                    h_B = obj.matToColDE(H_B);
                    hP(obj.numGrad.G+1:obj.numGrad.G+obj.numGrad.B) = h_B;
                    
                    % Va
                    H_Va = zeros(obj.numBus, obj.numSnap);
                    h_Va = obj.dataO.Vm(bus, snap) * obj.dataO.Vm(:, snap) .* GBThetaQ(:, bus);
                    h_Va(bus) = - obj.dataO.Vm(bus, snap)^2 * obj.dataO.B(bus, bus)...
                       - obj.data.Q_noised(bus, snap); 
%                     h_Va(bus) = h_Va(bus)-sum(GBThetaQ(bus, :) * obj.dataO.Vm(:, snap) * obj.dataO.Vm(bus, snap));
                    H_Va(:, snap) = h_Va;
                    % remove the source bus whose magnitude is not the state variable
                    H_Va(1, :) = []; 
                    h_VaLarge = reshape(H_Va', [], 1);
                    hP(obj.numGrad.G+obj.numGrad.B+1:end) = h_VaLarge;
                    
                    % Q
                    hQ = zeros(obj.numGrad.G + obj.numGrad.B + obj.numGrad.Va, 1);
                    % G matrix
                    H_G = zeros(obj.numBus, obj.numBus);
                    H_G(bus, :) = obj.dataO.Vm(bus, snap) * obj.dataO.Vm(:, snap)' .* sin(theta_ij');
                    h_G = obj.matToColDE(H_G);
                    hQ(1:obj.numGrad.G) = h_G;

                    % B matrix
                    H_B = zeros(obj.numBus, obj.numBus);
                    H_B(bus, :) = - obj.dataO.Vm(bus, snap) * obj.dataO.Vm(:, snap)' .* cos(theta_ij');
                    H_B(bus, :) = H_B(bus, :) + obj.dataO.Vm(bus, snap)^2; % the equivilance of diagonal elements
                    h_B = obj.matToColDE(H_B);
                    hQ(obj.numGrad.G+1:obj.numGrad.G+obj.numGrad.B) = h_B;
                    
                    % Va
                    H_Va = zeros(obj.numBus, obj.numSnap);
                    h_Va = - obj.dataO.Vm(bus, snap) * obj.dataO.Vm(:, snap) .* GBThetaP(:, bus);
                    h_Va(bus) = - obj.dataO.Vm(bus, snap)^2 * obj.dataO.G(bus, bus) ...
                        + obj.data.P_noised(bus, snap);
%                     h_Va(bus) = h_Va(bus)+sum(GBThetaP(bus, :) * obj.dataO.Vm(:, snap) * obj.dataO.Vm(bus, snap));
                    H_Va(:, snap) = h_Va;
                    % remove the source bus whose magnitude is not the state variable
                    H_Va(1, :) = []; 
                    h_VaLarge = reshape(H_Va', [], 1);
                    hQ(obj.numGrad.G+obj.numGrad.B+1:end) = h_VaLarge;
                    
                    obj.J((snap-1)*obj.numBus+bus, :) = hP';
                    obj.J((snap-1)*obj.numBus+bus+numP, :) = hQ';
                end
            end
        end
        
        function obj = updateParPF(obj)
            % This method updates the voltage angles by power flow
            % calculation.
            % We first build the branch matrix
            GB = triu(obj.dataO.G - obj.dataO.B, 1);
            [fBus, tBus] = find(GB);
            branchLoc = find(GB);
            numBranch = length(fBus);
            branch = repmat(obj.mpc.branch(1,:), numBranch, 1);
            branch(:, 1:2) = [fBus, tBus];
            
            y = - obj.dataO.G(branchLoc) - 1j * obj.dataO.B(branchLoc);
            z = 1 ./ y;
            branch(:, 3) = real(z);
            branch(:, 4) = imag(z);
            
            % We then update the bus matrix can do the PF calculations
            mpcO = obj.mpc; 
            mpcO.branch = branch;
%             Y = makeYbus(mpcO);
%             Gtest = real(full(Y));
%             Btest = imag(full(Y));
            mpcO.bus(mpcO.bus(:, 2) == 2, 2) = 1; % change all PV buses to PQ buses
%             mpcO.bus(mpcO.bus(:, 2) == 1, 2) = 2; % change all PQ buses to PV buses
            mpopt = mpoption('verbose',0,'out.all',0);
            
            for snap = 1:obj.numSnap
                mpcO.bus(2:end, 3) = - obj.dataO.P(2:end, snap) * mpcO.baseMVA;
                mpcO.bus(2:end, 4) = - obj.dataO.Q(2:end, snap) * mpcO.baseMVA;
                mpcO.bus(:, 8) = obj.dataO.Vm(:, snap);
                mpcThis = runpf(mpcO, mpopt);
                obj.dataO.Va(:, snap) = mpcThis.bus(:,9)/180*pi;
            end
        end
        
        function obj = updateParNewton(obj)
            % This method updates the parameters using Newton kings of
            % strategies
%             delta = obj.H \ obj.grad;
%             par = obj.parChain(:, obj.iter) - delta;
            
            delta = obj.H(1:obj.numGrad.B+obj.numGrad.G, 1:obj.numGrad.B+obj.numGrad.G) \ obj.grad(1:obj.numGrad.B+obj.numGrad.G);
            par = zeros(obj.numGrad.Sum, 1);
            par(1:obj.numGrad.B+obj.numGrad.G) = obj.parChain(1:obj.numGrad.B+obj.numGrad.G, obj.iter) - delta;
            
%             delta = obj.H(1+obj.numGrad.B+obj.numGrad.G:end,1+obj.numGrad.B+obj.numGrad.G:end) \ obj.grad(1+obj.numGrad.B+obj.numGrad.G:end);
%             par = zeros(obj.numGrad.Sum, 1);
%             par(1+obj.numGrad.B+obj.numGrad.G:end) = obj.parChain(1+obj.numGrad.B+obj.numGrad.G:end, obj.iter) - delta;
            % gather the par values
            G = par(1:obj.numGrad.G);
            B = par(1+obj.numGrad.G:obj.numGrad.G+obj.numGrad.B);
%             G(G>0) = 0;
%             B(B<0) = 0;
            
            obj.dataO.G = obj.colToMatDE(G, obj.numBus);
            obj.dataO.B = obj.colToMatDE(B, obj.numBus);
            Vm = par(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm);
            Va = par(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end);
%             obj.dataO.Vm(2:end, :) = reshape(Vm, [], obj.numBus-1)'; % exclude the source bus
%             obj.dataO.Va(2:end, :) = reshape(Va, [], obj.numBus-1)'; % exclude the source bus
        end
        
        function obj = updateParNewtonGV(obj)
            % This method updates the parameters using Newton strategy by
            % iteratively updating GB and VmVa
            
            if obj.isGB % update GB value, and Va. We also update the topology
                
                % correct the topology
                
                diagEle = sum(abs(obj.dataO.G)) / 2;
                ratio1 = abs(bsxfun(@rdivide, obj.dataO.G, diagEle));
                ratio2 = abs(bsxfun(@rdivide, obj.dataO.G, diagEle'));
                ratio = max(ratio1, ratio2);
                obj.Topo = ratio > obj.thsTopo;
                obj.Tvec = logical(obj.matOfColDE(obj.Topo));
                obj.dataO.G(~obj.Topo) = 0;
                obj.dataO.B(~obj.Topo) = 0;
                
%                 % update GB only
%                 id = [obj.Tvec; obj.Tvec];
%                 delta = zeros(obj.numGrad.G+obj.numGrad.B, 1);
%                 delta(id) = obj.H(id, id) \ obj.grad(id);
                
%                 % update GBVa using state estimation
%                 id = [1:obj.numGrad.B+obj.numGrad.G obj.numGrad.B+obj.numGrad.G+obj.numGrad.Vm+1:obj.numGrad.Sum];
%                 delta = obj.H(id, id) \ obj.grad(id);

                % update GBVa using Pseudo power flow
                deltaP = obj.dataO.P - obj.data.P_noised;
                deltaQ = obj.dataO.Q - obj.data.Q_noised;
                deltaPQ = [reshape(deltaP, [], 1);reshape(deltaQ, [], 1)];
                id = true(obj.numGrad.B+obj.numGrad.G+obj.numGrad.Va, 1);
                id(1:obj.numGrad.G+obj.numGrad.B) = [obj.Tvec; obj.Tvec];
                delta = zeros(obj.numGrad.B+obj.numGrad.G+obj.numGrad.Va, 1);
                delta(id) = pinv(obj.J(:,id)) * deltaPQ;
                
%                 delta = pinv(obj.J) * deltaPQ;
                if (max(abs(delta)) < 1e-2)
                    obj.isGB = false;
                    obj.isConverge = obj.isConverge + 1;
                else
                    obj.isConverge = 0;
                end
%                 obj.isGB = false;
                
                par = zeros(obj.numGrad.Sum, 1);
                id = [obj.Tvec; obj.Tvec];
                par(id) = ...
                    obj.parChain(id, obj.iter) ...
                    - delta(id);
                G = par(1:obj.numGrad.G);
                B = par(1+obj.numGrad.G:obj.numGrad.G+obj.numGrad.B);
                G(G>0) = -max(G)*rand()*0.2;
                B(B<0) = max(B)*rand()*0.2;
                G(G<-300) = -300;
                B(B>300) = 300;
                obj.dataO.G = obj.colToMatDE(G, obj.numBus);
                obj.dataO.B = obj.colToMatDE(B, obj.numBus);


                par(obj.numGrad.B+obj.numGrad.G+obj.numGrad.Vm+1:end) = ...
                    obj.parChain(obj.numGrad.B+obj.numGrad.G+obj.numGrad.Vm+1:end, obj.iter)...
                    - delta(obj.numGrad.B+obj.numGrad.G+1:end);
                Va = par(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end);
                obj.dataO.Va(2:end, :) = reshape(Va, [], obj.numBus-1)'; % exclude the source bus
            else % update VmVa value
                delta = obj.H(1+obj.numGrad.B+obj.numGrad.G:end,1+obj.numGrad.B+obj.numGrad.G:end) \ obj.grad(1+obj.numGrad.B+obj.numGrad.G:end);
                if (max(abs(delta)) < 1e-5)
                    obj.isGB = true;
                    obj.isConverge = obj.isConverge + 1;
                else
                    obj.isConverge = 0;
                end
%                 obj.isGB = true;
                
                par = zeros(obj.numGrad.Sum, 1);
                par(1+obj.numGrad.B+obj.numGrad.G:end) = obj.parChain(1+obj.numGrad.B+obj.numGrad.G:end, obj.iter) - delta;
                Vm = par(1+obj.numGrad.G+obj.numGrad.B:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm);
                Va = par(1+obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm:end);
                obj.dataO.Vm(2:end, :) = reshape(Vm, [], obj.numBus-1)'; % exclude the source bus
                obj.dataO.Va(2:end, :) = reshape(Va, [], obj.numBus-1)'; % exclude the source bus
            end
            
            % update PQ
            for i = 1:obj.numSnap
                % calculate some basic parameters at present state
                Theta_ij = repmat(obj.dataO.Va(:, i), 1, obj.numBus) - repmat(obj.dataO.Va(:, i)', obj.numBus, 1);
                % G_ij\cos(\Theta_ij)+B_ij\sin(\Theta_ij)
                GBThetaP = obj.dataO.G .* cos(Theta_ij) + obj.dataO.B .* sin(Theta_ij);
                % G_ij\sin(\Theta_ij)-B_ij\cos(\Theta_ij)
                GBThetaQ = obj.dataO.G .* sin(Theta_ij) - obj.dataO.B .* cos(Theta_ij);
                % P estimate
                Pest = (GBThetaP * obj.dataO.Vm(:, i)) .* obj.dataO.Vm(:, i);
                obj.dataO.P(:, i) = Pest;
                % Q estimate
                Qest = (GBThetaQ * obj.dataO.Vm(:, i)) .* obj.dataO.Vm(:, i);
                obj.dataO.Q(:, i) = Qest;
            end
        end
        
        function obj = buildHessian(obj)
            % This method builds the Hessian matrix
            obj.H = zeros(obj.numGrad.Sum, obj.numGrad.Sum);
            obj.HP = sparse(obj.numGrad.Sum, obj.numGrad.Sum);
            obj.HQ = sparse(obj.numGrad.Sum, obj.numGrad.Sum);
            obj.HVm = sparse(obj.numGrad.Sum, obj.numGrad.Sum);
            obj.HVa = sparse(obj.numGrad.Sum, obj.numGrad.Sum);
            
            obj.grad = zeros(obj.numGrad.Sum, 1);
            obj.gradP = sparse(obj.numGrad.Sum, 1);
            obj.gradQ = sparse(obj.numGrad.Sum, 1);
            obj.gradVm = sparse(obj.numGrad.Sum, 1);
            obj.gradVa = sparse(obj.numGrad.Sum, 1);
            
            obj.loss.total = 0;
            obj.loss.P = 0;
            obj.loss.Q = 0;
            obj.loss.Vm = 0;
            obj.loss.Va = 0;
            
            for i = 1:obj.numSnap
                % calculate some basic parameters at present state
                Theta_ij = repmat(obj.dataO.Va(:, i), 1, obj.numBus) - repmat(obj.dataO.Va(:, i)', obj.numBus, 1);
                % G_ij\cos(\Theta_ij)+B_ij\sin(\Theta_ij)
                GBThetaP = obj.dataO.G .* cos(Theta_ij) + obj.dataO.B .* sin(Theta_ij);
                % G_ij\sin(\Theta_ij)-B_ij\cos(\Theta_ij)
                GBThetaQ = obj.dataO.G .* sin(Theta_ij) - obj.dataO.B .* cos(Theta_ij);
                % P estimate
                Pest = (GBThetaP * obj.dataO.Vm(:, i)) .* obj.dataO.Vm(:, i);
                obj.dataO.P(:, i) = Pest;
                % Q estimate
                Qest = (GBThetaQ * obj.dataO.Vm(:, i)) .* obj.dataO.Vm(:, i);
                obj.dataO.Q(:, i) = Qest;
                
                % calculate the sub-matrix of P of all buses
                for j = 1:obj.numBus
                    if obj.isMeasure.P(j)
                        obj = buildHessianP(obj, i, j, GBThetaP, GBThetaQ, Pest);
                    end
                end
                
                % calculate the sub-matrix of Q of all buses
                for j = 1:obj.numBus
                    if obj.isMeasure.Q(j)
                        obj = buildHessianQ(obj, i, j, GBThetaP, GBThetaQ, Qest);
                    end
                end
                
                % calculate the sub-matrix of Vm of all buses
                for j = 1:obj.numBus
                    if obj.isMeasure.Vm(j)
                        obj = buildHessianVm(obj, i, j);
                    end
                end
                
                % calculate the sub-matrix of Va of all buses
                for j = 1:obj.numBus
                    if obj.isMeasure.Va(j)
                        obj = buildHessianVa(obj, i, j);
                    end
                end
            end
            
            % collect the Hessians, gradients and the losses
            obj.H = full(obj.HP + obj.HQ + obj.HVm + obj.HVa);
            obj.gradP = full(obj.gradP);
            obj.gradQ = full(obj.gradQ);
            obj.gradVm = full(obj.gradVm);
            obj.gradVa = full(obj.gradVa);
            obj.grad = obj.gradP + obj.gradQ + obj.gradVm + obj.gradVa;
            obj.loss.total = obj.loss.P + obj.loss.Q + obj.loss.Vm + obj.loss.Va;
        end
        
        function obj = buildHessianP(obj, snap, bus, GBThetaP, GBThetaQ, Pest)
            % This method builds the Hessian from the measurement of P
            
            theta_ij = obj.dataO.Va(bus, snap) - obj.dataO.Va(:, snap);
            g = sparse(obj.numGrad.Sum, 1);
            
            % G matrix
            H_G = zeros(obj.numBus, obj.numBus);
            H_G(bus, :) = obj.dataO.Vm(bus, snap) * obj.dataO.Vm(:, snap)' .* cos(theta_ij');
            H_G(bus, :) = H_G(bus, :) - obj.dataO.Vm(bus, snap)^2; % the equivilance of diagonal elements
            h_G = obj.matToColDE(H_G);
            g(1:obj.numGrad.G) = h_G;
            
            % B matrix
            H_B = zeros(obj.numBus, obj.numBus);
            H_B(bus, :) = obj.dataO.Vm(bus, snap) * obj.dataO.Vm(:, snap)' .* sin(theta_ij');
            h_B = obj.matToColDE(H_B);
            g(obj.numGrad.G+1:obj.numGrad.G+obj.numGrad.B) = h_B;
            
            % Vm
            % the first order term of other Vm
            H_Vm = zeros(obj.numBus, obj.numSnap);
            h_Vm = obj.dataO.Vm(bus, snap) * GBThetaP(:, bus);
            % the second order term of Vm(bus)
            h_Vm(bus) = 2*obj.dataO.Vm(bus, snap) * GBThetaP(bus, bus);
            % the first order term of Vm(bus)
            fOrderVm = obj.dataO.Vm(:, snap) .* GBThetaP(:, bus);
            fOrderVm(bus) = 0;
            h_Vm(bus) = h_Vm(bus) + sum(fOrderVm);
            H_Vm(:, snap) = h_Vm;
            % remove the source bus whose magnitude is not the state variable
            H_Vm(1, :) = []; 
            h_VmLarge = reshape(H_Vm', [], 1);
            g(obj.numGrad.G+obj.numGrad.B+1:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) = h_VmLarge;
            
            % Va
            H_Va = zeros(obj.numBus, obj.numSnap);
            h_Va = obj.dataO.Vm(bus, snap) * obj.dataO.Vm(:, snap) .* GBThetaQ(:, bus);
            h_Va(bus) = - obj.dataO.Vm(bus, snap)^2 * obj.dataO.B(bus, bus)...
                       - obj.data.Q_noised(bus, snap); 
%             h_Va(bus) = h_Va(bus)-sum(GBThetaQ(bus, :) * obj.dataO.Vm(:, snap) * obj.dataO.Vm(bus, snap));
            H_Va(:, snap) = h_Va;
            % remove the source bus whose magnitude is not the state variable
            H_Va(1, :) = []; 
            h_VaLarge = reshape(H_Va', [], 1);
            g(obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm+1:end) = h_VaLarge;
            
            % build HP, gradP and loss.P
            lossThis = (Pest(bus) - obj.data.P_noised(bus, snap));
            obj.loss.P = obj.loss.P + lossThis^2 * obj.sigma.P(bus).^(-2);
            gradPThis = obj.sigma.P(bus).^(-2) * lossThis * g;
            obj.gradP = obj.gradP + gradPThis;
            HPThis = obj.sigma.P(bus).^(-2) * (g * g');
            obj.HP = obj.HP + HPThis;
        end
        
        function obj = buildHessianQ(obj , snap, bus, GBThetaP, GBThetaQ, Qest)
            % This method builds the Hessian from the measurement of Q
            
            theta_ij = obj.dataO.Va(bus, snap) - obj.dataO.Va(:, snap);
            g = sparse(obj.numGrad.Sum, 1);
            
            % G matrix
            H_G = zeros(obj.numBus, obj.numBus);
            H_G(bus, :) = obj.dataO.Vm(bus, snap) * obj.dataO.Vm(:, snap)' .* sin(theta_ij');
            h_G = obj.matToColDE(H_G);
            g(1:obj.numGrad.G) = h_G;
            
            % B matrix
            H_B = zeros(obj.numBus, obj.numBus);
            H_B(bus, :) = - obj.dataO.Vm(bus, snap) * obj.dataO.Vm(:, snap)' .* cos(theta_ij');
            H_B(bus, :) = H_B(bus, :) + obj.dataO.Vm(bus, snap)^2; % the equivilance of diagonal elements
            h_B = obj.matToColDE(H_B);
            g(obj.numGrad.G+1:obj.numGrad.G+obj.numGrad.B) = h_B;
            
            % Vm
            % the first order term of other Vm
            H_Vm = zeros(obj.numBus, obj.numSnap);
            h_Vm = obj.dataO.Vm(bus, snap) * GBThetaQ(:, bus);
            % the second order term of Vm(bus)
            h_Vm(bus) = 2*obj.dataO.Vm(bus, snap) * GBThetaQ(bus, bus);
            % the first order term of Vm(bus)
            fOrderVm = obj.dataO.Vm(:, snap) .* GBThetaQ(:, bus);
            fOrderVm(bus) = 0;
            h_Vm(bus) = h_Vm(bus) + sum(fOrderVm);
            H_Vm(:, snap) = h_Vm;
            % remove the source bus whose magnitude is not the state variable
            H_Vm(1, :) = []; 
            h_VmLarge = reshape(H_Vm', [], 1);
            g(obj.numGrad.G+obj.numGrad.B+1:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) = h_VmLarge;
            
            % Va
            H_Va = zeros(obj.numBus, obj.numSnap);
            h_Va = - obj.dataO.Vm(bus, snap) * obj.dataO.Vm(:, snap) .* GBThetaP(:, bus);
            h_Va(bus) = - obj.dataO.Vm(bus, snap)^2 * obj.dataO.G(bus, bus) ...
                        + obj.data.P_noised(bus, snap);
%             h_Va(bus) = h_Va(bus)+sum(GBThetaP(bus, :) * obj.dataO.Vm(:, snap) * obj.dataO.Vm(bus, snap));
            H_Va(:, snap) = h_Va;
            % remove the source bus whose magnitude is not the state variable
            H_Va(1, :) = []; 
            h_VaLarge = reshape(H_Va', [], 1);
            g(obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm+1:end) = h_VaLarge;
            
            % build HQ, GradientQ and lossQ
            lossThis = (Qest(bus) - obj.data.Q_noised(bus, snap));
            obj.loss.Q = obj.loss.Q + lossThis^2 * obj.sigma.Q(bus).^(-2);
            gradQThis = obj.sigma.Q(bus).^(-2) * lossThis * g;
            obj.gradQ = obj.gradQ + gradQThis;
            HQThis = obj.sigma.Q(bus).^(-2) * (g * g');
            obj.HQ = obj.HQ + HQThis;
        end
        
        function obj = buildHessianVm(obj, snap, bus)
            % This method builds the Hessian from the measurement of Vm
            g = sparse(obj.numGrad.Sum, 1);
            H_Vm = sparse(obj.numBus, obj.numSnap);
            H_Vm(bus, snap) = 1;
            % remove the source bus whose magnitude is not the state variable
            H_Vm(1, :) = []; 
            h_VmLarge = reshape(H_Vm', [], 1);
            g(obj.numGrad.G+obj.numGrad.B+1:obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm) = h_VmLarge;
            
            % build GradientVm and lossVm
            lossThis = (obj.dataO.Vm(bus, snap) - obj.data.Vm_noised(bus, snap));
            obj.loss.Vm = obj.loss.Vm + lossThis^2 * obj.sigma.Vm(bus).^(-2);
            gradVmThis = obj.sigma.Vm(bus).^(-2) * lossThis * g;
            obj.gradVm = obj.gradVm + gradVmThis;
            HVmThis = obj.sigma.Vm(bus).^(-2) * (g * g');
            obj.HVm = obj.HVm + HVmThis;
        end
        
        function obj = buildHessianVa(obj, snap, bus)
            % This method builds the Hessian from the measurement of Va
            g = sparse(obj.numGrad.Sum, 1);
            H_Va = sparse(obj.numBus, obj.numSnap);
            H_Va(bus, snap) = 1;
            % remove the source bus whose angle is not the state variable
            H_Va(1, :) = []; 
            h_VaLarge = reshape(H_Va', [], 1);
            g(obj.numGrad.G+obj.numGrad.B+obj.numGrad.Vm+1:end) = h_VaLarge;
            
            % build HVa, GradientVa and lossVa
            lossThis = (obj.dataO.Va(bus, snap) - obj.data.Va_noised(bus, snap));
            obj.loss.Va = obj.loss.Va + lossThis^2 * obj.sigma.Va(bus).^(-2);
            gradVaThis = obj.sigma.Va(bus).^(-2) * lossThis * g;
            obj.gradVa = obj.gradVa + gradVaThis;
            HVaThis = obj.sigma.Va(bus).^(-2) * (g * g');
            obj.HVa = obj.HVa + HVaThis;
        end
        
        function obj = identifyMCMCEIV(obj)
            % This method uses the Markov Chain Monte Carlo to sample the
            % distribution of the parameters and the topologies. We use the
            % error-in-variables(EIV) assumption.
            % We may have some errors regarding to the usage of reshape Vm
            % and Va
            
            % Build the measurement function.
            % Currently, we assume we know all the P, Q, Vm, Va
            % measurements. We have to modify it later.
            data.Pn = obj.data.P_noised;
            data.Qn = obj.data.Q_noised;
            data.Vmn = obj.data.Vm_noised;
            data.Van = obj.data.Va_noised;
            data.num = obj.numFIM;
            data.isMeasure = obj.isMeasure;
            data.sigma = obj.sigma;
            
            % build the parameters
            G = obj.matOfCol(obj.dataE.G);
            B = obj.matOfCol(obj.dataE.B);
            Vm = reshape(obj.data.Vm_noised(2:end,:)', [], 1); % we assume the value of the source bus is already known
            Va = reshape(obj.data.Va_noised(2:end,:)', [], 1);
            par = [G' B' Vm' Va'];
            assert (length(par) == obj.numFIM.Sum) % the number of total parameters
            
            % we also build the ground truth value of the parameters
            Gtr = obj.matOfCol(obj.data.G);
            Btr = obj.matOfCol(obj.data.B);
            Vmtr = reshape(obj.data.Vm(2:end,:), [], 1);
            Vatr = reshape(obj.data.Va(2:end,:), [], 1);
            obj.truePar = [Gtr' Btr' Vmtr' Vatr'];
            
            % build the params in the format of mcmc
            params = cell(1, data.num.Sum);
            for i = 1:data.num.G
                params{i} = ...
                    {sprintf('G_{%d}',i), G(i), -Inf, Inf, obj.boundA.total(i)};
            end
            for i = 1:data.num.B
                params{i+data.num.G} = ...
                    {sprintf('B_{%d}',i), B(i), -Inf, Inf, obj.boundA.total(i+data.num.G)};
            end
            for i = 1:data.num.Vm
                params{i+data.num.G+data.num.B} =...
                    {sprintf('Vm_{%d}',i), Vm(i), Vm(i)*0.99, Vm(i)*1.01, obj.boundA.total(i+data.num.G+data.num.B)};
            end
            for i = 1:data.num.Va
                params{i+data.num.G+data.num.B+data.num.Vm} =...
                    {sprintf('Vm_{%d}',i), Va(i), Va(i)-abs(Va(i))*0.01, Va(i)+abs(Va(i))*0.01, obj.boundA.total(i+data.num.G+data.num.B+data.num.Vm)};
            end
            
            % build the model
            model.ssfun = @sumOfSquaresEIV;
            % build the sigma2 (the sum of squares error of the measurements)
            sigma2P = obj.sigma.P(obj.isMeasure.P).^2';
            sigma2Q = obj.sigma.Q(obj.isMeasure.Q).^2';
            sigma2Vm = obj.sigma.Vm(obj.isMeasure.Vm).^2';
            sigma2Va = obj.sigma.Va(obj.isMeasure.Va).^2';
            model.sigma2 = [sigma2P sigma2Q sigma2Vm sigma2Va] * obj.numSnap.^2 *1000;
            model.S20 = model.sigma2;
            model.N = obj.numSnap;
            
            % build the options
            options.nsimu = 500000;
            options.qcov = obj.boundA.cov;
            numGB = obj.numFIM.G+obj.numFIM.B;
            options.qcov(1:numGB,1:numGB) = options.qcov(1:numGB,1:numGB) * 1;
%             options.updatesigma = 1;

            % run the mcmc simulation
            [res,chain,s2chain] = mcmcrun(model,data,params,options);
            
            % run the mcmc simulation and update the cov matrix iteratively
            options.nsimu = 2000;
            numIter = 10;
            Gs = cell(1, numIter);
            Bs = cell(1, numIter);
            chains = [];
            errorInit = sum(sumOfSquaresEIV(par, data) ./ model.sigma2);
            errorTrue = sum(sumOfSquaresEIV(obj.truePar, data) ./ model.sigma2);
            errorEval = sum(sumOfSquaresEIV(res.theta', data) ./ model.sigma2);
            for i = 1:10
                [res,chain,~] = mcmcrun(model,data,params,options);
                errorEval = sum(sumOfSquaresEIV(res.theta', data) ./ model.sigma2)
                Gs{i} = res.theta(1:data.num.G);
                Bs{i} = res.theta(1+data.num.G:data.num.G+data.num.B);
                chains = [chains;chain];
                % rebuild the FIM matrix and the cov matrix
                obj.dataE.G = obj.colToMat(Gs{i}, obj.numBus);
                obj.dataE.B = obj.colToMat(Bs{i}, obj.numBus);
                obj = approximateFIM(obj);
                obj = calABound(obj);
                options.qcov = obj.boundA.cov;
                numGB = obj.numFIM.G+obj.numFIM.B;
                options.qcov(1:numGB,1:numGB) = options.qcov(1:numGB,1:numGB)*100;
            end
            errorInit = sum(sumOfSquaresEIV(par, data) ./ model.sigma2);
            errorTrue = sum(sumOfSquaresEIV(obj.truePar, data) ./ model.sigma2);
            errorEval = sum(sumOfSquaresEIV(res.theta', data) ./ model.sigma2);
        end
        
        function obj = identifyMCMCEIO(obj)
            % This method uses the Markov Chain Monte Carlo to sample the
            % distribution of the parameters and the topologies. We use the
            % error-in-outputs(EIV) assumption.
            % Build the measurement function.
            
            % Currently, we assume we know all the P, Q, Vm, Va
            % measurements. We have to modify it later.
            data.Pn = obj.data.P_noised;
            data.Qn = obj.data.Q_noised;
            data.Vmn = obj.data.Vm_noised;
            data.Van = obj.data.Va_noised;
            data.num = obj.numFIM;
            data.isMeasure = obj.isMeasure;
            data.sigma = obj.sigma;
            
            % build the parameters
%             G = obj.matOfCol(obj.dataE.G);
%             B = obj.matOfCol(obj.dataE.B);
            G = obj.matOfCol(obj.data.G);
            B = obj.matOfCol(obj.data.B);
            par = [G' B'];
            
            % we also build the ground truth value of the parameters
            Gtr = obj.matOfCol(obj.data.G);
            Btr = obj.matOfCol(obj.data.B);
            obj.truePar = [Gtr' Btr'];
            
            % build the params in the format of mcmc
            params = cell(1, data.num.G+data.num.B);
            for i = 1:data.num.G
                params{i} = ...
                    {sprintf('G_{%d}',i), G(i), -Inf, Inf, obj.boundA.total(i)};
            end
            for i = 1:data.num.B
                params{i+data.num.G} = ...
                    {sprintf('B_{%d}',i), B(i), -Inf, Inf, obj.boundA.total(i+data.num.G)};
            end
            
            % build the model
            model.ssfun = @sumOfSquaresEIO;
            % build the sigma2 (the sum of squares error of the measurements)
            % we use the summation of G and B matrices to approximate the
            % first order of measurement noises
%             sumG = diag(obj.dataE.G);
%             sumB = diag(obj.dataE.B);
            sumG = diag(obj.data.G);
            sumB = diag(obj.data.B);
            sigma2P = sumG(obj.isMeasure.P).^2';
            sigma2Q = sumB(obj.isMeasure.Q).^2';
            model.sigma2 = [sigma2P sigma2Q] / 10000000;
            model.S20 = model.sigma2;
            model.N = obj.numSnap;
            
            % build the options
            options.nsimu = 50000;
            numGB = obj.numFIM.G+obj.numFIM.B;
            options.qcov = obj.boundA.cov(1:numGB, 1:numGB);
            
            % run the mcmc
            [res,chain,s2chain] = mcmcrun(model,data,params,options);
            errorInit = sum(sumOfSquaresEIO(par, data) ./ model.sigma2);
            errorTrue = sum(sumOfSquaresEIO(obj.truePar, data) ./ model.sigma2);
            errorEval = sum(sumOfSquaresEIO(res.theta', data) ./ model.sigma2);
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
        
        function h = matOfCol(H)
            % This method get the half triangle of a matrix
            H_up = tril(H, 1)';
            n = size(H, 1);
            N = (n + 1) * n / 2;
            h = zeros(N, 1);
            pt = 1;
            for i = 1:n
                h(pt:pt+n-i) = H_up(i, i:end);
                pt = pt+n-i+1;
            end
        end
        
        function h = matOfColDE(H)
            % This method get the half triangle of a matrix The name DE
            % denotes diagonal exclude, which means we consider the
            % diagonal elements as the negative summation of the rest elements.
%             H_up = tril(H, 1)';
            n = size(H, 1);
            N = (n - 1) * n / 2;
            h = zeros(N, 1);
            pt = 1;
            for i = 1:n
                h(pt:pt+n-i-1) = H(i, i+1:end);
                pt = pt+n-i;
            end
        end
    end
end

