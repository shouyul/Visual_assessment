function [train_acc, test_acc] = compute_train_test_acc(data_mat,data_label)

% get training and testing accuracies
% Initialize matrices to store training and testing accuracies
num_sets = max(data_label(:,2));

%%
train_acc = zeros(num_sets, 1);
test_acc = zeros(num_sets, 1);

% Iterate over each fold
for i = 1:num_sets
    % Extract training indices for the current fold
    trainIndices = ones(size(data_label,1),1);
    trainIndices = data_label(:,2)~=i;
    
    % Extract training data and labels for the current fold
    X_train = data_mat(trainIndices, :);
    y_train = data_label(trainIndices,1);
    
    % Create an SVM template
    t = templateSVM('Standardize',true','KernelFunction','linear',...
        'BoxConstraint',1);

    % Train a new model on the training data
    SVMModel_fold = fitcecoc(X_train, y_train,...
        'Learners',t);
    
    % Evaluate the model on the training data
    y_pred_train = predict(SVMModel_fold, X_train);
    
    % Compute the training accuracy for the current fold
    train_acc(i) = sum(y_pred_train == y_train) / numel(y_train);
    
    % Extract test indices for the current fold
    testIndices = ones(size(data_label,1),1);
    testIndices = data_label(:,2)==i;
    
    % Extract test data and labels for the current fold
    X_test = data_mat(testIndices, :);
    y_test = data_label(testIndices,1);
    
    % Use the trained model to predict labels for the test data
    y_pred_test = predict(SVMModel_fold, X_test);
    
    % Compute the testing accuracy for the current fold
    test_acc(i) = sum(y_pred_test == y_test) / numel(y_test);
end

train_acc = mean(train_acc);
test_acc = mean(test_acc);

end