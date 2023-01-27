function [outputWeights, outputIntercept] = gaussian_elimination_gcvHansen(trainInput, trainTeacher, stateCollection, esn)

% PSEUDOINVERSE computes the outputWeights using the standard pseudoinverse
% operation
%
% inputs:
% trainInput = input time series of size nInputUnits x (nTrainingPoints-nForgetPoints) 
% trainTeacher = teacher time series of size nOutputUnits x (nTrainingPoints-nForgetPoints) 
% stateCollectMat = matrix of size nInternalUnits x (nTrainingPoints-nForgetPoints)
%
% output:
% outputWeights = ifRegressInput = 0 then it is matrix of size nInternalUnits x nOutputUnits containing the learnt weights
%   ifRegressInput = 1 then it is matrix of size (nInternalUnits + nInputUnits) x nOutputUnits containing the learnt weights

% outputIntercept = vector of size nOutputUnits x 1 containing the learnt
%   intercept


% Created April 30, 2006, D. Popovici
% Copyright: Fraunhofer IAIS 2006 / Patent pending
% Revision 1, Feb 13, 2017, L. Grigoryeva and J.-P. Ortega: different
% conventions.
% Revision 2, March 13, 2017, L. Grigoryeva and J.-P. Ortega:  Generalized to multivariate inputs and outputs

if esn.ifRegressInput == 0
    x = stateCollection;
else
    x =[stateCollection; trainInput];
end
y = feval(esn.inverseOutputActivationFunction, trainTeacher);
[N, T] = size(x);

x1 = x;
y1 = y;
x = x - repmat(mean(x,2), 1, T);
y = y - mean(y);
x = x./repmat(std(x1, [], 2), 1, T);
% 
% mx = mean(x1, 2);
% stdx = std(x1,0,2);
% mx = mx';
% stdx = stdx';
% MX = mx(ones(T,1),:);
% STDX = stdx(ones(T,1),:);
% Z = (x1' - MX) ./ STDX;
% Gamma0 = (Z' * Z);
% invGamma0 = inv(Gamma0 + esn.ridgeLambda * eye(N));
% Cov_XY = Z' * y';
% outputWeights = (invGamma0 * Cov_XY);
% 
% muY = mean(y1, 2);
% muX = mean(x1, 2);
% 
%        outputWeights = outputWeights./stdx';
%        outputIntercept = muY- outputWeights' * muX;
% sum((y1' - outputIntercept *ones(1999,1) - x1'*outputWeights).^2)
% % get the mean of the reservoir output


% 
% % compute the GAMMA(0) of the reservoir output
% Gamma0 = (x1 * x1') - T * muX * muX';
% %Gamma0 = 1/T * (x * x') - muX * muX';
% 
% threshold = esn.ridgeLambda;
% min_eigval = min(eig(Gamma0, 'vector'));
% if abs(min_eigval) < threshold
%     %invGamma0 = invCov_corrected(Gamma0, 0, threshold);
%     invGamma0 = inv(Gamma0 + esn.ridgeLambda * eye(N));
% else
%     invGamma0 = inv(Gamma0);
% end
% %invGamma0 = inv(Gamma0 + esn.ridgeLambda * eye(N));
% % get the Cov(x,y) of the reservoir output and the teaching signal
% Cov_XY = (1/T) * x * y' - muX * muY';
% % get the reservoir output matrix and intercept
% %outputWeights = ((Gamma0 + esn.ridgeLambda * eye(N))) \ Cov_XY;
% outputWeights = (T * invGamma0 * Cov_XY);
% outputIntercept = muY- outputWeights' * muX;
% % 
% % 
% 
% l=0:0.0000001:0.0001
% b = ridge(y',x',l,0);

% fVal = @(lambda)gcv(lambda, x, y);
% options = optimoptions('fmincon', 'Display', 'off');
% lambda = fmincon(fVal, 1.e-6, [], [], [], [], 1.e-6, 1.e-1, [], options);
% fValx = gcv(lambda, x, y);
% i = 1;
% fValx1=[];
% for lambda = linspace(1.e-8, 1.e-4,500)
%     fValx1(i) =gcv(lambda, x, y);
%     i = i + 1;
% end
%try
[U, s, V] = csvd(x');
[reg_min, G, reg_param] = gcvHansen(U, s, y', 'Tikh'); % performs cross validation, while respecting the temporal properties
display(reg_min);
b = ridge(y1', x1', reg_min, 0);
outputWeights = b(2:end);
outputIntercept = b(1);
% catch
%     b = ridge(y1', x1', 1e-4, 0);
%     outputWeights = b(2:end);
%     outputIntercept = b(1);
% end

% er = 1/T*sum((y1' - outputIntercept *ones(T,1) - x1'*outputWeights).^2)



% hyperopts = struct('Optimizer','gridsearch','ShowPlots', 0, 'Verbose', 1);%, 'NumGridDivisions', 500);
% [Mdl, FitInfo] = fitrlinear(x1, y1, 'ObservationsIn','columns','Learner','leastsquares','Lambda', ...
%     reg_min/T);
% [Mdl, FitInfo] = fitrlinear(x1, y1, 'ObservationsIn','columns','Lambda', ...
%     10.^(-(10:-2:2)), 'Learner','leastsquares','OptimizeHyperParameters', {'Lambda'}, ...
%     'HyperparameterOptimizationOptions', hyperopts);
% 
% outputIntercept = Mdl.Bias;
% outputWeights = Mdl.Beta;
% hyperopts = struct('Optimizer','gridsearch','ShowPlots', 0, 'Verbose', 0);
% [Mdl] = fitrlinear(x1', y1','Learner','leastsquares','OptimizeHyperParameters', {'Lambda'}, ...
%     'HyperparameterOptimizationOptions', hyperopts);
% %FitInfo, HyperparameterOptimizationResults
% 
% YHat = predict(Mdl,x1');
% YHat = YHat';
% er = 1/T*sum((y1 - YHat).^2)
% 
% 
% % YHat1 = Mdl.Bias + x1'*Mdl.Beta;
% % YHat1 = YHat1';
% 
% outputIntercept = Mdl.Bias;
% outputWeights = Mdl.Beta;

%display(esn.ridgeLambda);
% 
% er = 1/T*sum((y1' - outputIntercept *ones(T,1) - x1'*outputWeights).^2)
% er = 1/T*sum((y1' - outputIntercept1 *ones(T,1) - x1'*outputWeights1).^2)
% figure
% plot(l,b,'LineWidth',2)
% ylim([-100 100])
% grid on
% xlabel('Ridge Parameter')
% ylabel('Standardized Coefficient')
% title('{\bf Ridge Trace}')
% legend('x1','x2','x3','x1x2','x1x3','x2x3')
% 
% 
% [N, T] = size(x);
% x=x-repmat(mean(x,2),1,T);
% x = x./repmat(std(x,[],2),1,T);
% y = y - mean(y);
% er=[];
% j = 1;
% for l=linspace(0, 5e-4,500)
% er(j) = 1/T * norm((eye(T) - x'*((x*x' + T*l*eye(N))\x)) * y')^2 / ...
%     (1/T * trace(eye(T) - x'*((x*x' + T*l*eye(N))\x) ))^2;
% j = j+1;
% end

