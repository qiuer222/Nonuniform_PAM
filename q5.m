q = 5;
% n = floor((q-1)/2);
di = 0.01;
snr_y = -1:0.001:6;
sigma_n = sqrt(1./(10.^(snr_y/10)));
C = zeros(1,length(snr_y));
MI_goal = 0.5;
figure;
[~,now] = size(snr_y);
for delta1 = -1:di:0
    %     for delta2 = -1:di:0
    %         for delta3 = -1:di:0
    %             for delta4 = -1:di:0
    %             for delta4 = -1:di:1
    %                 for delta5 = -1:di:1
    %                     delta = [delta1 delta2 delta3 delta4 delta5]
    %                     delta = [delta1 delta2 delta3  delta4 -delta4 -delta3 -delta2 -delta1]
    delta = [-1 delta1 0 -delta1 1]
    delta = delta/sqrt(mean(delta.^2));
    [~,isu] = size(unique(delta))
    if(isu~=q)
        continue;
    end
    py = zeros(size(y));
    for i = 1:1:q
        py = 1/(q*sqrt(2*pi*sigma_n(now)^2))*exp(-(y-delta(i)).^2/(2*sigma_n(now)^2))+py;
    end
    C_now = -sum(py.*log2(py+eps))*0.01-0.5*log2(2*pi*exp(1)*sigma_n(now)^2);
    if C_now>MI_goal
        for j = now:-1:1
            py = zeros(size(y));
            for i = 1:1:q
                py = 1/(q*sqrt(2*pi*sigma_n(j)^2))*exp(-(y-delta(i)).^2/(2*sigma_n(j)^2))+py;
            end
            C(j) = -sum(py.*log2(py+eps))*0.01-0.5*log2(2*pi*exp(1)*sigma_n(j)^2);
            if(C(j)<MI_goal)
                C(j:now)
                plot(snr_y(j:now),C(j:now));hold on;
                now = j+1;
                delta_now = delta;
                break
            end
            
        end
%     else
%         break;
    end
    
    xlabel('SNR/dB');
    ylabel('MI');
    
end
%         end
%     end
%         end
%     end
% end
delta = delta_now
snr_y = -5:0.1:20;
sigma_n = sqrt((1)./(10.^(snr_y/10)));
C = zeros(1,length(snr_y));
for j = 1:1:length(snr_y)
    py = zeros(size(y));
    for i = 1:1:q
        py = 1/(q*sqrt(2*pi*sigma_n(j)^2))*exp(-(y-delta(i)).^2/(2*sigma_n(j)^2))+py;
    end
    C(j) = -sum(py.*log2(py+eps))*0.01-0.5*log2(2*pi*exp(1)*sigma_n(j)^2);
end
%  delta =  [  -1.5102   -0.4682         0    0.4682    1.5102]%MI=1

