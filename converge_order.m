clear;clc;
tic
fdm_obj=fdm_class();
%计算时间的阶
Nh=300;
Nt=[5,10,20,40,80];
num_Nt=length(Nt);
Nt_initial=Nt(1);
results_t=zeros(num_Nt,7);
results_t(:,1)=Nt';

u_fdm_t_1=fun_fdm(Nh,Nt_initial,fdm_obj);
error_t_1=error_compute(Nh,u_fdm_t_1,fdm_obj);
results_t(1,2:2:end)=error_t_1;

for i=2:num_Nt
    Nt_i=Nt(i);
    u_fdm_t_2=fun_fdm(Nh,Nt_i,fdm_obj);
    error_t_2=error_compute(Nh,u_fdm_t_2,fdm_obj);
    results_t(i,2:2:end)=error_t_2;

    %计算收敛阶
    results_t(i,3:2:end)=log(results_t(i-1,2:2:end)./results_t(i,2:2:end))./log(2);
end

disp('关于时间的误差与收敛阶')
disp('-----------------------------------------------------');
disp('网格数 Nt |    L1误差   | L1阶  |    L2误差   | L2阶  |   Linf误差  | Linf阶  |');
disp('-----------------------------------------------------');
disp(results_t);
% 将结果写入 CSV 文件（数值按 5 位有效数字）
filename = fullfile(pwd, 'time_error_convergence1.csv'); % 保存到当前工作目录

% 准备带表头的字符串表格，数值按 5 位有效数字格式化
headers = {'N','L1_error','L1_order','L2_error','L2_order','Linf_error','Linf_order'};
num_rows = size(results_t,1);
formatted = cell(num_rows, numel(headers));
for i = 1:num_rows
    formatted{i,1} = sprintf('%d', results_t(i,1)); % N 用整数格式
    for j = 2:size(results_t,2)
        formatted{i,j} = sprintf('%.5g', results_t(i,j)); % 5 位有效数字
    end
end

T = cell2table(formatted, 'VariableNames', headers);
writetable(T, filename);
fprintf('时间误差结果已写入 %s\n', filename);



%计算空间的阶
Nt=500;
Nh=[5,10,20,40,80];
num_Nh=length(Nh);
Nh_initial=Nh(1);
results_h=zeros(num_Nh,7);
results_h(:,1)=Nh';

u_fdm_h_1=fun_fdm(Nh_initial,Nt,fdm_obj);
error_h_1=error_compute(Nh_initial,u_fdm_h_1,fdm_obj);
results_h(1,2:2:end)=error_h_1;

for i=2:num_Nh
    Nh_i=Nh(i);
    u_fdm_h_2=fun_fdm(Nh_i,Nt,fdm_obj);
    error_h_2=error_compute(Nh_i,u_fdm_h_2,fdm_obj);
    results_h(i,2:2:end)=error_h_2;

    %计算收敛阶
    results_h(i,3:2:end)=log(results_h(i-1,2:2:end)./results_h(i,2:2:end))./log(2);
end

disp('关于空间的误差与收敛阶')
disp('-----------------------------------------------------');
disp('网格数 Nh |    L1误差   | L1阶  |    L2误差   | L2阶  |   Linf误差  | Linf阶  |');
disp('-----------------------------------------------------');
disp(results_h);

% 将结果写入 CSV 文件（数值按 5 位有效数字）
filename = fullfile(pwd, 'space_error_convergence1.csv'); % 保存到当前工作目录

% 准备带表头的字符串表格，数值按 5 位有效数字格式化
headers = {'N','L1_error','L1_order','L2_error','L2_order','Linf_error','Linf_order'};
num_rows = size(results_h,1);
formatted = cell(num_rows, numel(headers));
for i = 1:num_rows
    formatted{i,1} = sprintf('%d', results_h(i,1)); % N 用整数格式
    for j = 2:size(results_h,2)
        formatted{i,j} = sprintf('%.5g', results_h(i,j)); % 5 位有效数字    
    end
end

T = cell2table(formatted, 'VariableNames', headers);
writetable(T, filename);
fprintf('空间误差结果已写入 %s\n', filename);
toc

