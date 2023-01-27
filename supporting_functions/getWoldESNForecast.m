function [esn1_t, esn2_t, esn3_t, Forecasts, epsilons, Vt, StoPart, trainError, testError, InSampleFit] = getWoldESNForecast(...
    esn1, esn2, esn3, sample, T, Horizon, nForecasts, bootstrapBlocksize, nForgetPoints)

% nForgetPoints is used in the training, computation of errors, and choice
% of epsilons for bootstrapping
if isempty(nForgetPoints)
    nForgetPoints = 0;
end

trainTeacher = sample(1, 1:T);
testTeacher = sample(1, T + 1:end);

% First reservoir

[esn1_t, ~] = ...
    train_esn([], trainTeacher, esn1, nForgetPoints, 'ifRegressInput', 0);

[stateCollectionTrain1, predictedTrainOutput1] = simulate_esn([], trainTeacher, esn1_t, 0);
predictedTrainOutput1 = predictedTrainOutput1(1,1:end-1);

trainError = compute_NRMSE(predictedTrainOutput1(nForgetPoints + 1:end), trainTeacher(:, nForgetPoints + 1:end));
disp(sprintf('train NRMSE = %s', num2str(trainError)))

if ~isempty(testTeacher)
    [~, predictedTestOutput1] = simulate_esn([], testTeacher, esn1_t, 0,...
        'startingState', stateCollectionTrain1(:, end-1), 'startingTeacher', trainTeacher(1, end));
    predictedTestOutput1 = predictedTestOutput1(1,1:end-1);
    
    testError = compute_NRMSE(predictedTestOutput1, testTeacher(:, 1:end));
    disp(sprintf('test NRMSE = %s', num2str(testError)))
else
    disp(sprintf('test not possible for testing sample is empty'));
    testError = [];
end

% compute innovations
epsilons = trainTeacher - predictedTrainOutput1;

% Second reservoir

trainInput2 = epsilons;
esn2.teacherMask = zeros(esn2.nInternalUnits, esn2.nInputUnits);
[esn2_t, ~] = ...
    train_esn(trainInput2, trainTeacher, esn2, nForgetPoints, 'ifRegressInput', 0);
[stateCollectionTrain2, predictedTrainOutput2] = simulate_esn(trainInput2, [], esn2_t, 0);

% compute predictable innovations and stochastic linear part in Wold decompostion
StoPart = predictedTrainOutput2;
Vt = trainTeacher - predictedTrainOutput2;

% third reservoir

trainTeacher3 = Vt;
[esn3_t, ~] = ...
    train_esn([], trainTeacher3, esn3, nForgetPoints, 'ifRegressInput', 0);
[stateCollectionTrain3, predictedTrainOutput3] = simulate_esn([], trainTeacher3, esn3_t, 0);

InSampleFit = predictedTrainOutput2 + predictedTrainOutput3(1,1:end-1);

% Forecasting

q = fix(Horizon/bootstrapBlocksize);
r = rem(Horizon,bootstrapBlocksize);
epsilons_forecast = zeros(1,Horizon);
Epsilons = epsilons(1, nForgetPoints+1:end);
Forecasts = zeros(nForecasts, Horizon);
for i=1:nForecasts
    [~, idx] = datasample(Epsilons(1:end-bootstrapBlocksize),q + 1);
    
    for j = 1:q
        for k=1:bootstrapBlocksize
            epsilons_forecast(1, (j-1)*bootstrapBlocksize + k) = Epsilons(1, idx(j) + k);
        end
    end
    for j=1:r
        epsilons_forecast(1, q*bootstrapBlocksize + j) = Epsilons(1, idx(q+1) + j);
    end
    
    
    [stateCollectionForecast2, ForecastOutput2] = simulate_esn(epsilons_forecast, [], esn2_t, 0,...
        'startingState', stateCollectionTrain2(:, end));
    
    
    [stateCollectionForecast3, ForecastOutput3] = simulate_esn([], [], esn3_t, 0,...
        'startingState', stateCollectionTrain3(:, end-1), 'startingTeacher', predictedTrainOutput3(:, end-1),...
        'lengthSimul', Horizon);
    
    Forecasts(i, :) = ForecastOutput2 + ForecastOutput3;
end
