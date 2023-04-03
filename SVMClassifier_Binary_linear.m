function [features, Observed, Shuffled] = SVMClassifier_Binary_linear(Data,varargin)
% Camden MacDowell 2019
%
% See MATLAB documentation for fitcsvm for more details.  
% Citation for decoding neural data using classifiers: 
% Pereira, Mitchell, Botvinick, 2008 NeuroImage. Machine Learning 
% Classifiers and fMRI: a tutorial overview
% 
% @SYNPOSIS
% Performs linear binary classification on matrix Data according to the 
% class identities in the last column of Data (see @Data). Optionally
% tunes classifier hyperparameters kernalscale and boxcontrain using pre-built 
% baysian optimization (baysopt fnc). Optionally performs feature selection
% (currently only supports ANOVA feature selection). 
% 
% @INPUTS
% Data: observations x features matrix. (i.e. trials x pixels for imaging data).
% Final column contains the response variable with the class identify of
% each row. 
%
% @varargin: see opts structure below
%
% @Optimized Hyperparameters:
% @Box contraint: is the regularization parameter, where
% increase = hard margin and thus high cost on missclassified points. 
%  
% @Kernal scale is 1/gamma. This is the radius of influence of the sample
% selected by the model as support vectors. If small then you end up with
% lots of support vectors that really only apply to those specific points
% (low radius of influence). If large than the opposite and the model is to
% 'rigid'/contrained to capture the shape fo the data.
%
% @Kernel
% Linear

%% Set optional paramters
opts.holdout = 0.25; %fraction hold out data for validation
opts.nshuf = 0; %Number of shuffles for shuffle test set to zero for no shuffle test
opts.featureselect  = 'none'; %Options: 'none', 'anova': 
opts.numberfeatures = 50; %Number of features to select
opts.tuneparam = 0; %Tune hyperparameters using baysian optimization
opts.pca = 0;  %First perform dimensionality reduction. Uses #PCs to explain 99% of variance. 
opts.solver = 1; %SVM solver to use: 1=SMO,2=ISDA,3=L1QP


%Process optional inputs
if mod(length(varargin), 2) ~= 0, error('Must pass key/value pairs for options.'); end
for i = 1:2:length(varargin)
    try
        opts.(varargin{i}) = varargin{i+1};
    catch
        error('Couldn''t set option ''%s''.', varargin{2*i-1});
    end
end

%%
%set appropriate solver;
solver = {'SMO','ISDA','L1QP'};

Observed = struct(); %Observed results structure
Shuffled = struct(); %Shuffled results structure
predictors = Data(:, 1:end-1);
response = Data(:,end);

%Hold out data for validation. This is not used in tuning of classifier
cvp = cvpartition(response, 'Holdout', opts.holdout);
trainingPredictors = predictors(cvp.training, :);
trainingResponse = response(cvp.training, :);
validationPredictors = predictors(cvp.test, :);
validationResponse = response(cvp.test, :);

%preform pca on the testing data, project validation data into pca space
if opts.pca
    fprintf('\nPerforming pca')
    [coef, score, ~, ~, explain, mu] = pca(trainingPredictors);  
    trainingPredictors = score(:,1:find(cumsum(explain)>99,1));  
    validationPredictors = (coef(:,1:find(cumsum(explain)>99,1))'*(validationPredictors-repmat(mu,size(validationPredictors,1),1))')';    
end

%Optional feature selection on training data; 
[features] = FeatureSelection(trainingPredictors,trainingResponse,opts.featureselect,opts.numberfeatures);
trainingPredictors(:,~features)=[];
trainingResponse(~features) = [];
validationPredictors(:,~features)=[];
validationResponse(~features) = [];

%Optional Tuning of classifier hyperparameters using cross validation: Only use the training data for this.
if opts.tuneparam
    c = cvpartition(trainingResponse, 'kfold', 10);
    optstune = struct('Optimizer','bayesopt','ShowPlots',false,'CVPartition',c,...
        'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',100);
    svmmod = fitcsvm(trainingPredictors,trainingResponse,'KernelFunction','linear',...
        'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',optstune,'Solver',solver{opts.solver});
    % Build a classifier using tuned parameters and all training data
    classificationSVM = fitcsvm(trainingPredictors,trainingResponse,'KernelFunction','linear','BoxConstraint',...
        svmmod.HyperparameterOptimizationResults.XAtMinObjective.BoxConstraint,...
        'KernelScale',svmmod.HyperparameterOptimizationResults.XAtMinObjective.KernelScale,'Solver',solver{opts.solver});
else
    classificationSVM = fitcsvm(trainingPredictors,trainingResponse,'KernelFunction','linear','BoxConstraint',1,...
        'KernelScale',1,'Solver',solver{opts.solver});
end

% Create the result struct with predict function
svmPredictFcn = @(x) predict(classificationSVM, x);
validationPredictFcn = @(x) svmPredictFcn(x);

% Compute validation predictions
[validationPredictions, validationScores] = validationPredictFcn(validationPredictors);
Outstats = [];
for i = 1:numel(unique(response))
    [Outstats(i).X,Outstats(i).Y,Outstats(i).T,Outstats(i).AUC,~,Outstats(i).SUBY] = perfcurve(validationResponse,validationScores(:,i),i);
end

% Compute validation accuracy
correctPredictions = (validationPredictions == validationResponse);
isMissing = isnan(validationResponse);
correctPredictions = correctPredictions(~isMissing);

%Save off all desired information : Combine since it's a binary classifier 
Observed.AUC = mean([Outstats(:).AUC]);
Observed.X = [cat(1,Outstats(:).X),cat(1,Outstats(:).X)];
Observed.Y = [cat(1,Outstats(:).Y),cat(1,Outstats(:).Y)];
Observed.T = [cat(1,Outstats(:).T),cat(1,Outstats(:).T)];
Observed.Accurary = sum(correctPredictions)/length(correctPredictions);
Observed.Predictions = validationPredictions;
Observed.CorrectResponse = validationResponse;
Observed.Scores = validationScores;
Observed.Classifier = classificationSVM;
Observed.trainingResponse = trainingResponse;
Observed.trainingPredictors = trainingPredictors;
Observed.validationPredictors = validationPredictors;
Observed.validationResponse = validationResponse;
Observed.betas = {classificationSVM.Beta};
Observed.bias = classificationSVM.Bias;


% Now randomly shuffle the held out labels and use the same classifier 
% Compute validation predictions
if opts.nshuf>=1
    for shuf = 1:opts.nshuf
        validationResponse_shuf = validationResponse(randperm(numel(validationResponse)));
        [validationPredictions, validationScores] = validationPredictFcn(validationPredictors);
        Outstats = [];

        for i = 1:numel(unique(response))
            [Outstats(i).X,Outstats(i).Y,Outstats(i).T,Outstats(i).AUC,~,~] =...
            perfcurve( validationResponse_shuf,validationScores(:,i),i);
        end

        % Compute validation accuracy
        correctPredictions = (validationPredictions ==  validationResponse_shuf);
        isMissing = isnan( validationResponse_shuf);
        correctPredictions = correctPredictions(~isMissing);
        validationAccuracy = sum(correctPredictions)/length(correctPredictions);

        Shuffled(shuf).AUC = mean([Outstats(:).AUC]);
        Shuffled(shuf).X = [cat(1,Outstats(:).X),cat(1,Outstats(:).X)];
        Shuffled(shuf).Y = [cat(1,Outstats(:).Y),cat(1,Outstats(:).Y)];
        Shuffled(shuf).T = [cat(1,Outstats(:).T),cat(1,Outstats(:).T)];
        Shuffled(shuf).Accurary = validationAccuracy;
        Shuffled(shuf).Predictions = validationPredictions;
        Shuffled(shuf).CorrectResponse =  validationResponse_shuf;
        Shuffled(shuf).Scores = validationScores;
    end
else
    Shuffled = struct();
end
end


%% Subfunctions
function [features] = FeatureSelection(trainingPredictors,trainingResponse,type,number)
    
    switch type
        case 'none'
%             fprintf('\nPerforming no feature selection'); 
            features = (1:1:numel(trainingResponse));
        case 'anova'
            fprintf('\nPerforming feature selection using ANOVA');
            %Compute the anova fstat between response for each feature
            stat = zeros(1,size(trainingPredictors,2));
            for i = 1:size(trainingPredictors,2)
                stat(i) = anovan(trainingPredictors(:,i),trainingResponse,'display','off');                
            end
            stat = -log10(stat);    
            if numel(stat>number)
                [~,features] = maxk(stat,number); %reduce number of features
            else
                features = (1:1:numel(stat));
            end
        otherwise
            fprintf('\nUnknown feature selection\n no selection performed')
            features = (1:1:numel(trainingResponse));
    end   
        
    
end









