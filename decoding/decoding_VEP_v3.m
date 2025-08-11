function [acc,dprime] = decoding_VEP_v3(trials,Info,crossval_method,kfold)
% 
% % trials: chan x feat x trial
% % Info: trial type, block number, trial number
% 
X = reshape(permute(trials,[3,1,2]),size(trials,3),[]);
y = Info(:,1);

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

conf = confusionmat(y,y_pred);
acc(1) = conf(1,1)/sum(conf(1,:));
acc(2) = conf(2,2)/sum(conf(2,:));
acc(3) = conf(3,3)/sum(conf(3,:));
acc(4) = conf(4,4)/sum(conf(4,:));

dprime = sEEG_compute_dprime(y,y_pred);

end