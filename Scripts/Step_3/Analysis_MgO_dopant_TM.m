%%
% Code written by 
% Chun-Yen Liu
% Ph.D. Student
% Department of Chemical and Biomolecular Engineering
% Rice University, Houtson, TX 77005
%
% Adapted from:
% Interaction trends between single metal atoms and oxide supports identified with density functional theory and statistical learning
% Nolan J. O'Connor, A S M Jonayat, Michael J. Janik*, Thomas P. Senftle* 

% Distributed as part of the publication - 
% Using Statistical Learning to Predict Interactions Between Single Metal Atoms and Modified MgO(100) Supports
% Chun-Yen Liu, Shijia Zhang, Daniel Martinez, Meng Li*, and Thomas P. Senftle*

% Department of Chemical and Biomolecular Engineering and Department of Statistics, Rice University, Houston, TX 77005 (USA)
% *meng@rice.edu
% *tsenftle@rice.edu

%%
clear all;
clc

%% Inputs
file = '../Step_2/dopant/MgO_dopant_LS.txt'; % Specify the descriptor set
Dimension = 1:8;                             % Determine the dimension

%% load dataset generated in data generation from step 1 directory
load('../Step_1/dopant/data_set_dopant.mat')

%% Calculate RMSE and R^2 for SL-derived descriptors

%% Import data: the selected descriptors
fileID = fopen(file,'r');
read = fgetl(fileID);
Pair = zeros(1,2);
line = 0;
while ischar(read)
    line = line + 1;
    read2 = split(read, '"');
    read2 = read2(~cellfun('isempty',read2));
    read3 = split(read2, '*');
    if length(read3) == 1
        Pair(line,1) = str2double(read3);
    else
        Pair(line,:) = str2double(read3);
    end   
    if isnan(Pair(line,1)) | isnan(Pair(line,2))
        Pair(line,:) = [];
        line = line - 1;
    end
    read = fgetl(fileID);
end
fclose(fileID);

%%  Use the descriptors selected by dopant-modified MgO
Y = Data_d_train;
D_p = Data_d_p_train;
head_p = head_d_all;
Data = zeros(length(Y),size(Pair,1));
Head = strings(1,size(Pair,1));
for i = 1:size(Pair,1)
    if Pair(i,2) == 0
        Data(:,i) = D_p(:,Pair(i,1));
        Head(i)   = head_p(Pair(i,1));
    else
    Data_temp = [D_p(:,Pair(i,1)) D_p(:,Pair(i,2))];
    Head_temp = [head_p(Pair(i,1)) head_p(Pair(i,2))];
    Unit_temp = [1, 1];
    Operation = {'*'};
    [Data(:,i), Head(i), Unit_temp_1] = GenBasicFeature(Data_temp, Head_temp, Unit_temp, Operation, 0);
    end
end

%%  Use the descriptors selected by dopant-modified MgO for test set
Y_test = Data_d_test;
D_p_test = Data_d_p_test;
Data_t = zeros(length(Y_test),size(Pair,1));
Head = strings(1,size(Pair,1));
for i = 1:size(Pair,1)
    if Pair(i,2) == 0
        Data_t(:,i) = D_p_test(:,Pair(i,1));
        Head(i)   = head_p(Pair(i,1));
    else
    Data_temp = [D_p_test(:,Pair(i,1)) D_p_test(:,Pair(i,2))];
    Head_temp = [head_p(Pair(i,1)) head_p(Pair(i,2))];
    Unit_temp = [1, 1];
    Operation = {'*'};
    [Data_t(:,i), Head(i), Unit_temp_1] = GenBasicFeature(Data_temp, Head_temp, Unit_temp, Operation, 0);
    end
end

%% Refined with l0-norm method

RMSE_L_zero = zeros(size(Pair,1), size(Pair,1));
Rsq_L_zero = zeros(size(Pair,1), size(Pair,1));
Data_L_zero = Data(:,1:size(Pair,1));
min_RMSE = zeros(size(Pair,1), 1);
max_Rsq = zeros(size(Pair,1), 1);
RMSE_test = zeros(size(Pair,1), 1);
Rsq_test = zeros(size(Pair,1), 1);

for i = Dimension
    %tic
    num_des = 1:1:size(Pair,1);
    sz = nchoosek(size(Pair,1),i);
    index = int8(nchoosek(num_des,i));
    coeff = zeros(i+1, sz);
    parfor j = 1: sz
        Data_temp = ones(size(Y,1), i+1);
        for k = 1: i
            Data_temp(:,k) = Data(:,index(j,k));
        end
        coeff(:, j) = regress(Y, Data_temp);
        Predicted = Data_temp*coeff(:, j);
        RMSE_L_zero(j, i) = sqrt(sum((Y - Predicted).^2)/size(Y,1));
        Rsq_L_zero(j, i) = 1 - sum((Y - Predicted).^2)/sum((Y - mean(Y)).^2);
    end
    
    idx_RMSE = [];
    idx_Rsq = [];
    
    [min_RMSE(i), idx_RMSE] = min(RMSE_L_zero(1:sz,i));
    [max_Rsq(i), idx_Rsq] = max(Rsq_L_zero(1:sz,i));
    
    % Calculate the test error
    coeff_t = regress(Y, [Data(:,index(idx_RMSE, :)) ones(length(Y),1)]);
    Data_t_2 = Data_t(:,index(idx_RMSE, :));
    Predicted_t = [Data_t_2 ones(length(Y_test),1)]*coeff_t;
    RMSE_test(i) = sqrt(sum((Y_test - Predicted_t).^2)/size(Y_test,1));
    Rsq_test(i) = 1 - sum((Y_test - Predicted_t).^2)/sum((Y_test - mean(Y_test)).^2);
    
    disp('Descriptors for lowest RMSE')
    for j = 1: i
        disp(Head(index(idx_RMSE,j)))
    end
    for j = 1: i
        disp([ Pair(index(idx_RMSE,j),1), Pair(index(idx_RMSE,j),2)])
    end
    disp('Coefficients for lowest RMSE')
    disp(coeff(:, idx_RMSE)')
    disp('Index for lowest RMSE')
    disp(index(idx_RMSE, :))
    
    %toc
end

figure(1)
subplot(2,1,1);
plot(1:max(Dimension), min_RMSE(Dimension),'bo') 
hold on
plot(1:max(Dimension), min_RMSE(Dimension),'b-')
hold on
plot(1:max(Dimension), RMSE_test(Dimension),'go') 
hold on
plot(1:max(Dimension), RMSE_test(Dimension),'g-')
xlabel('Number of descriptors')
ylabel('RMSE (eV)')
xlim([0 max(Dimension)+1])
ylim([0 0.8])
subplot(2,1,2);
plot(1:max(Dimension), max_Rsq(Dimension),'ro')
hold on
plot(1:max(Dimension), max_Rsq(Dimension),'r-')
xlabel('Number of descriptors')
ylabel('R^2')
xlim([0 max(Dimension)+1])
ylim([0 1.2])

figure(2)
Data_p = Data(:,index(idx_RMSE, :));
coeff_p = regress(Y, [Data_p ones(length(Y),1)]);
Predicted_p = [Data_p ones(length(Y),1)]*coeff_p;
plot(Predicted_p, Y,'ro')
xlabel('Predicted Binding Energy Difference (eV)')
ylabel('DFT Binding Energy Difference (eV)')
hold on
Data_t_2 = Data_t(:,index(idx_RMSE, :));
Predicted_test = [Data_t_2 ones(length(Y_test),1)]*coeff_p;
plot(Predicted_test, Y_test,'go')
hold on
plot([1 -5], [1, -5],'k-')


