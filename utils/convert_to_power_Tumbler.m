function [trial_power,dB_present,dB_absent,freq] = convert_to_power_Tumbler(trials,trialtype)

srate = 250;
[~,freq] = cwt(trials(1,:,1),srate);
trialtype = trialtype(1:size(trials,3));
% idx = (freq>=1&freq<=42);

for ch = 1:size(trials,1)
    
    temp_wavelet = nan(length(freq),size(trials,2),size(trials,3));
    
    for t = 1:size(trials,3)
        temp_wavelet(:,:,t) = cwt(trials(ch,:,t),srate);
    end

    temppower = abs(temp_wavelet);
    baseline = mean(mean(temppower,3),2);
    
    trial_power(ch,:,:) = squeeze(mean(temppower,2));

    temppower_present = temppower(:,:,trialtype==1);
    temppower_absent = temppower(:,:,trialtype==0);

    dB_present(ch,:,:) = 10*log10(mean(temppower_present,3)./baseline);
    dB_absent(ch,:,:) = 10*log10(mean(temppower_absent,3)./baseline);
end

end