function ss = sumOfSquaresEIO(par, data)
    % This method calculate the sum of squares error of the
    % likelihood function.This is an error-in-outputs formulation.
    
    % first we set a constant to reduce the value of the ss
%     kReduce = 100; % the reduce ratio is k^2
%     data.sigma.P = data.sigma.P * kReduce;
%     data.sigma.Q = data.sigma.Q * kReduce;
%     data.sigma.Vm = data.sigma.Vm * kReduce;
%     data.sigma.Va = data.sigma.Va * kReduce;
    
    [numBus, numSnap] = size(data.Pn);
    numMP = sum(data.isMeasure.P);
    numMQ = sum(data.isMeasure.Q);
    numMeasure = numMP + numMQ;
    % recover the parameters
    G = colToMat(par(1:data.num.G), numBus);
    B = colToMat(par(1+data.num.G:data.num.G+data.num.B), numBus);
    
    % define the estimation matrices
    est.P = zeros(numBus, numSnap);
    est.Q = zeros(numBus, numSnap);
    est.Vm = zeros(numBus, numSnap);
    est.Va = zeros(numBus, numSnap);
    
    for snap = 1:numSnap
        % the power flow equation, P and Q injections
        Theta_ij = repmat(data.Van(:, snap), 1, numBus) - repmat(data.Van(:, snap)', numBus, 1);
        % G_ij\cos(\Theta_ij)+B_ij\sin(\Theta_ij)
        GBThetaP = G .* cos(Theta_ij) + B .* sin(Theta_ij);
        % G_ij\sin(\Theta_ij)-B_ij\cos(\Theta_ij)
        GBThetaQ = G .* sin(Theta_ij) - B .* cos(Theta_ij);
        est.P(:, snap) = (GBThetaP * data.Vmn(:, snap)) .* data.Vmn(:, snap);
        est.Q(:, snap) = (GBThetaQ * data.Vmn(:, snap)) .* data.Vmn(:, snap);
    end
    
    ss = zeros(1, numMeasure);
    ss(1:numMP) = (sum((est.P - data.Pn) .^ 2, 2))'; % ./ data.sigma.P
    ss(1+numMP:numMP+numMQ) = (sum((est.Q - data.Qn) .^ 2, 2) )'; % ./ data.sigma.Q
end

function H = colToMat(h, n)
    % This method transform the column of half triangle to a
    % symmetric matrix
    H = zeros(n, n);
    pt = 1;
    for i = 1:n
        H(i, i:end) = h(pt:pt+n-i);
        pt = pt+n-i+1;
    end
    H = H + triu(H, 1)';
end

