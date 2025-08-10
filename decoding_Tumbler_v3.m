function [acc] = decoding_Tumbler_v3(X,y,crossval_method,kfold)
% 
% % X: chan x feat x trial
% % y: trial type
% 
X = reshape(permute(X,[3,1,2]),size(X,3),[]);

t = templateSVM('Standardize',true);
Mdl = fitcecoc(X,y,'Learners',t);

switch crossval_method
    case 'kfold'
        CVMdl = crossval(Mdl,'KFold',kfold);
    case 'leaveout'
        CVMdl = crossval(Mdl,'Leaveout','on');
end

% acc = 1-kfoldLoss(CVMdl);
[y_pred,scores] = kfoldPredict(CVMdl);

acc = mean(y==y_pred);

end