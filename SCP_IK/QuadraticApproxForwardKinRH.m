% Generate a quadratic approximation of the forward kinematics of the right
% hand, in the trust region, with Xcurrent (at time T) as the center
% Example usage:
% [P,q,r, residuals] = QuadraticApproxForwardKinRH(TrustRegionMin, TrustRegionMax, numParticles, XcurrentAtT);
function [P,q,r, residual] = QuadraticApproxForwardKinRH(TrustRegionMin, TrustRegionMax, numParticles, XcurrentAtT, index)
    n = size(TrustRegionMin,1);
    Z = GenerateRandomX(TrustRegionMin, TrustRegionMax, numParticles);
    m = 6; % Problem specific to the ForwardKinRH function
    Y = zeros(m,numParticles);
    for i=1:numParticles
        Y(:,i) = ForwardKinRH(Z(:,i));
        %Y(:,i) = [Z(:,i);1].^2;
    end

    cvx_begin quiet
        variable P1(n,n)
        variable q1(n)
        variable r1
        f = 0;
        for i=1:numParticles
           f = f + abs((Z(:,i)-XcurrentAtT)'*P1*(Z(:,i)-XcurrentAtT) + q1'*(Z(:,i)-XcurrentAtT) + r1 - Y(index,i));
        end
        
        minimize f 
        subject to 
            P1 == semidefinite(n,n);
        
    cvx_end
    
    P = P1;
    q = q1;
    r = r1;
    residual = cvx_optval/(max(Y(index,:))-min(Y(index,:)))/numParticles;
    
%     test = zeros(size(Y(index,:)));
%     for i=1:numParticles
%        test(i) = (Z(:,i)-XcurrentAtT)'*P*(Z(:,i)-XcurrentAtT) + q'*(Z(:,i)-XcurrentAtT) + r;
%     end
%     timeAxis = 1:numParticles;
%     figure; plot(timeAxis, Y(index,:), 'b', timeAxis, test, 'r');
%     figure; plot(Y(index,:)-test);
%     errorPercent = norm(Y(index,:)-test,1)/(max(Y(index,:))-min(Y(index,:)))/numParticles;
%     disp(errorPercent);
    
end