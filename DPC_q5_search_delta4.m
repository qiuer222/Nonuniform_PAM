%norm sir fixed s,Px=1
q = 5;
s = [-1 1]
% n = floor((q-1)/2);
di = 0.001;
delta = zeros(1,q)
snr_y = -2:0.001:8;
dy = 0.001;
y = -5*q:dy:5*q;
sigma_n = sqrt((1)./(10.^(snr_y/10)));
C = zeros(1,length(snr_y));
MI_goal = 1;
figure;
[~,now] = size(snr_y);
ub = 0;
counts = 0;
for delta1 = -1:di:0
    
    %     for delta2 = delta1:di:-di
    %         for delta3 = -1:di:0
    %             for delta4 = -1:di:0
    %             for delta4 = -1:di:1
    %                 for delta5 = -1:di:1
    %                     delta = [delta1 delta2 delta3 delta4 delta5]
    %                     delta = [delta1 delta2 delta3  delta4 -delta4 -delta3 -delta2 -delta1]
    delta = [-1 delta1 0 -delta1 1]
    [~,isu] = size(unique(delta));
    if(isu~=q)
        continue;
    end
    %%
    delta = delta/sqrt(mean(delta.^2));
    for del_delta = (delta(end)-delta(1)):0.1:(delta(end)-delta(1))*(q+2)/(q-1)
        %norm Px
        for pp = 0.01:0.01:3
            delta_pp = pp*delta;
            sub_mat = delta_pp' - s;
            del_delta_pp = pp*del_delta;
            k_max = floor( (max(sub_mat,[],'all')-delta_pp(1))/del_delta_pp);
            k_min = floor( (min(sub_mat,[],'all')-delta_pp(1))/del_delta_pp);
            ec = delta_pp'+del_delta_pp*(k_min:k_max);
            lambda_set = del_delta_pp*(k_min:k_max);
            [~,Q_sub_mat_i] = min(abs(sub_mat(:)'-lambda_set(:)));
            x = sub_mat-reshape(lambda_set(Q_sub_mat_i),q,[]);
            Px = mean(x.^2,'all')
            if(abs(Px-1)<1e-3)
                break
            end
        end
        if pp==3
            counts = counts+1;
            continue
        end
        pp
        del_delta_pp = pp*del_delta;
        n_s = length(s);
        % c-s_q
        sub_mat = delta_pp' - s;
        k_max = floor( (max(sub_mat,[],'all')-delta_pp(1))/del_delta_pp);
        k_min = floor( (min(sub_mat,[],'all')-delta_pp(1))/del_delta_pp);
        ec = delta_pp'+del_delta_pp*(k_min:k_max);
        lambda_set = del_delta_pp*(k_min:k_max);
        [~,Q_sub_mat_i] = min(abs(sub_mat(:)'-lambda_set(:)));
        x = sub_mat-reshape(lambda_set(Q_sub_mat_i),q,[]);
        Px = mean(x.^2,'all');
        
        r = x+s;
        
        r_set = unique(r)';
        n_r = length(r_set);
        mat_rc = zeros(q,n_r);
        for i = 1:1:q
            mat_rc(i,:) = histcounts(r(i,:),[r_set,10]);
        end
        
        % counts
        p_r = sum(mat_rc,1)/(sum(mat_rc,'all')+eps);
        p_rc = mat_rc./(sum(mat_rc,2)+eps);
        C3 = zeros(1,length(snr_y));
        
        p_yr=1/(sqrt(2*pi*sigma_n(now)^2))*exp(-(y-r_set').^2/(2*sigma_n(now)^2));
        p_y = p_r*p_yr;
        p_yc = p_rc*p_yr;
        C_ub = -sum(p_y.*log2(p_y+eps))*dy+sum(sum(p_yc.*log2(p_yc+eps)))/q*dy
        if(C_ub>MI_goal)
            for j = now:-1:1
                p_yr=1/(sqrt(2*pi*sigma_n(j)^2))*exp(-(y-r_set').^2/(2*sigma_n(j)^2));
                p_y = p_r*p_yr;
                p_yc = p_rc*p_yr;
                C3(j) = -sum(p_y.*log2(p_y+eps))*dy+sum(sum(p_yc.*log2(p_yc+eps)))/q*dy;
                if(C3(j)<MI_goal)
                    plot(snr_y(j:now),C3(j:now));hold on;
                    now = j+1;
                    delta_now = delta_pp;
                    del_delta_now = del_delta_pp;
                    
                    break
                end
            end
        else
            break;
        end
    end
    %         end
    %     end
    %         end
    %     end
end

fprintf("snr_c: %f\n",snr_y(now));
% 0.405000

%%
delta = delta_now
del_delta = del_delta_now
snr_y = -5:0.1:20;
sigma_n = sqrt((1)./(10.^(snr_y/10)));

% s_q = delta(1:4);
n_s = length(s);
% c-s_q
sub_mat = delta' - s;
% mod delta
k_max = floor( (max(sub_mat,[],'all')-delta(1))/del_delta);
k_min = floor( (min(sub_mat,[],'all')-delta(1))/del_delta);
ec = delta'+del_delta*(k_min:k_max);
lambda_set = del_delta*(k_min:k_max);
[~,Q_sub_mat_i] = min(abs(sub_mat(:)'-lambda_set(:)));
x = sub_mat-reshape(lambda_set(Q_sub_mat_i),q,[]);
Px = mean(x.^2,'all')
r = x+s;

r_set = unique(r)';
n_r = length(r_set);

mat_rc = zeros(q,n_r);
for i = 1:1:q
    mat_rc(i,:) = histcounts(r(i,:),[r_set,10]);
end

% counts
p_r = sum(mat_rc,1)/(sum(mat_rc,'all')+eps);
p_rc = mat_rc./(sum(mat_rc,2)+eps);
C3 = zeros(1,length(snr_y));
for j = 1:1:length(snr_y)
    p_yr=1/(sqrt(2*pi*sigma_n(j)^2))*exp(-(y-r_set').^2/(2*sigma_n(j)^2));
    p_y = p_r*p_yr;
    p_yc = p_rc*p_yr;
    C3(j) = -sum(p_y.*log2(p_y+eps))*dy+sum(sum(p_yc.*log2(p_yc+eps)))/q*dy;
end