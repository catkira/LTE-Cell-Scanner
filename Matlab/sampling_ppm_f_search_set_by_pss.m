% Xianjun Jiao (putaoshu@msn.com)
% Find out LTE PSS in the signal stream and correct sampling&carrier error.
% A script of project: https://github.com/JiaoXianjun/rtl-sdr-LTE

function xc = sampling_ppm_f_search_set_by_pss(s, fo_search_set)

[~, td_pss] = pss_gen(); % generate td zadoff chu sequences
len_pss = size(td_pss, 1);
len = length(s);
len_short = len - (len_pss-1);

if 1
    pss_fo_set = pss_fo_set_gen(td_pss, fo_search_set);
    num_fo_pss = size(pss_fo_set,2);
    corr_store = zeros(len_short, num_fo_pss);
    tic;
    for i=1:num_fo_pss
        tmp_corr = abs( filter(pss_fo_set(end:-1:1, i), 1, s) ).^2; % perform convolution and abs()^2
        tmp_corr = tmp_corr(len_pss:end); % remove len_pss items from beginning, because of non full overlap
        corr_store(:,i) = tmp_corr;
    end
    cost_time1 = toc;
    disp(['PSS xcorr cost ' num2str(cost_time1)]);
else
    tic;
    [corr_store, fo_search_set] = fft_corr(s, td_pss, fo_search_set);
    cost_time2 = toc;
    disp(['PSS xcorr cost ' num2str(cost_time2)]);
end

%fo_idx_set = 1:length(fo_search_set);
n_f = length(fo_search_set);
xc=zeros(3,len_short,n_f);
for foi=1:n_f
  for t=1:3
    col_idx = (t-1)*length(fo_search_set) + foi;
    xc(t,:,foi)=corr_store(:,col_idx);
  end
end

disp('Corr done.');
