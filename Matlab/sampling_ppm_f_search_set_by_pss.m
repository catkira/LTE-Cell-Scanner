% Xianjun Jiao (putaoshu@msn.com)
% Find out LTE PSS in the signal stream and correct sampling&carrier error.
% A script of project: https://github.com/JiaoXianjun/rtl-sdr-LTE

function [ppm, f_set, xc, fo_idx_set, pss_idx_set, fo_pss_idx_set, fo_with_all_pss_idx, extra_info] = sampling_ppm_f_search_set_by_pss(s, fo_search_set, td_pss, pss_fo_set, max_reserve, num_pss_period_try, combined_pss_peak_range, par_th, num_peak_th)
extra_info = [];

len_pss = size(pss_fo_set, 1);
num_fo_pss = size(pss_fo_set, 2);

len = length(s);
len_short = len - (len_pss-1);

% corr_store = zeros(len_short, num_fo_pss);
% 
% tic;
% for i=1:num_fo_pss
%     tmp_corr = abs( filter(pss_fo_set(end:-1:1, i), 1, s) ).^2;
%     tmp_corr = tmp_corr(len_pss:end);
%     corr_store(:,i) = tmp_corr;
% end
% cost_time1 = toc;
% disp(['PSS xcorr cost ' num2str(cost_time1)]);

tic;
[corr_store, fo_search_set] = fft_corr(s, td_pss, fo_search_set);
cost_time2 = toc;
disp(['PSS xcorr cost ' num2str(cost_time2)]);

ppm = inf;
f_set = fo_search_set;

fo_idx_set = 1:length(f_set);
n_f = length(f_set);
xc=zeros(3,len_short,n_f);
for foi=1:n_f
  for t=1:3
    col_idx = (t-1)*length(fo_search_set) + fo_idx_set(foi);
    xc(t,:,foi)=corr_store(:,col_idx);
  end
end

pss_idx_set = inf;
fo_pss_idx_set = inf;
fo_with_all_pss_idx = inf;
disp('Corr done.');
