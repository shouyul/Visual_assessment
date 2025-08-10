function [trial_power, power_time, dB_time, freq] = convert_to_power_SSVEP(trials)

srate = 250;
[~,freq] = cwt(trials(1,:,1),srate);

for ch = 1:size(trials,1)
    disp(ch)
    y = [];
    m = [];
    for t = 1:size(trials,3)
        y(:,t) = fft(trials(ch,:,t),500);
        m(:,t) = abs(y(:,t));
    end

    y(m<1e-6) = 0;
    f = (0:size(y,1)-1)*100/size(y,1);        % Frequency vector
end