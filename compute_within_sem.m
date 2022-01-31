function [D_new,sem] = compute_within_sem(D,varargin)
we_want_ci = false;
if length(varargin)==1
    % we want a CI
    alpha = varargin{1};
    we_want_ci = true;
end   
D_new = D - mean(D,1,'omitnan');

D_new = D_new + mean(mean(D,1,'omitnan'),2,'omitnan');
% 
if we_want_ci == true
    % CI uses alpha
    for cond_i = size(D_new,2):-1:1
        n    = size(D_new(~isnan(D_new(:,cond_i))),1);
        tVal = tinv(1-alpha/2,n-1);

        sem(1,cond_i) = mean(D_new(:,cond_i),1,'omitnan') + tVal * std(D_new(:,cond_i),1,'omitnan')./sqrt(n);
        sem(2,cond_i) = mean(D_new(:,cond_i),1,'omitnan') - tVal * std(D_new(:,cond_i),1,'omitnan')./sqrt(n);
    end
else
    % lack of CI means standard error 
    for cond_i = size(D_new,2):-1:1
        n    = size(D_new(~isnan(D_new(:,cond_i))),1);
        sem(cond_i,1) = std(D_new(:,cond_i),1,'omitnan') / sqrt(n);
    end
end

end