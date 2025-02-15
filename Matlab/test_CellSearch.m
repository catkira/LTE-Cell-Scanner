% Xianjun Jiao (putaoshu@msn.com)
% CellSearch.m
% Improved LTE-Cell-Scanner (written by James Peroulas: https://github.com/Evrytania/LTE-Cell-Scanner).

% Some scripts are borrowed from:
% https://github.com/JiaoXianjun/rtl-sdr-LTE
% https://github.com/Evrytania/LTE-Cell-Scanner
% https://github.com/Evrytania/Matlab-Library
% https://github.com/JiaoXianjun/multi-rtl-sdr-calibration

% clear all;
% close all;
function test_CellSearch(test_sp, test_ep)
test_source_info = regression_test_source('../regression_test_signal_file');

test_sp = 2;
test_ep = 2;
%test_ep = length(test_source_info);
%test_sp = 10; test_ep = 10;

f_search_set = -140e3:5e3:140e3;

filename = ['CellSearch_test' num2str(test_sp) 'to' num2str(test_ep) '_fo' num2str(min(f_search_set)/1e3) 'to' num2str(max(f_search_set)/1e3) '.mat'];

for i = test_sp : test_ep
    disp(test_source_info(i).filename);
    coef_pbch = pbch_filter_coef_gen(test_source_info(i).fs);
    
    r_raw = get_signal_from_bin(test_source_info(i).filename, inf, test_source_info(i).dev);
    r_raw = r_raw - mean(r_raw); % remove DC

    r_pbch = filter_wo_tail(r_raw, coef_pbch, (30.72e6/16)/test_source_info(i).fs);
    [~, ~, ~] = CellSearch(r_pbch, r_raw, f_search_set, test_source_info(i).fc, test_source_info(i).fs);
    save(filename, 'test_source_info', 'test_sp', 'test_ep', 'f_search_set');
end

