%
% Adapted from:
% Interaction trends between single metal atoms and oxide supports identified with density functional theory and statistical learning
% Nolan J. O'Connor, A S M Jonayat, Michael J. Janik*, Thomas P. Senftle* 

% Distributed as part of the publication - 
% Using Statistical Learning to Predict Interactions Between Single Metal Atoms and Modified MgO(100) Supports
% Chun-Yen Liu, Shijia Zhang, Daniel Martineza, Meng Li*, and Thomas P. Senftle*

% Department of Chemical and Biomolecular Engineering and Department of Statistics, Rice University, Houston, TX 77005 (USA)
% *meng@rice.edu
% *tsenftle@rice.edu

%%
clear all;
clc

%% Inputs
file = '../Step_2/MgO_dopant_LS.txt'; % Specify the descriptor set
Dimension = 1:3;                      % Determine the dimension

%% load dataset generated in data generation from step 1 directory
load('../Step_1/data_set.mat')

%% Calculate RMSE and R^2 for SL-derived descriptors

%% Import data: most relevant pairs
fileID = fopen(file,'r');
Pair = fscanf(fileID, [char(34) '%d' char(42) '%d' char(34) '\n'], [2 Inf])';
fclose(fileID);

%% Check whether a dopant-modified or adsorbate-modified data 
if isempty(strfind(file,'adsorbate'))
    Y = dBE_eV_d;
    D_p = D_p_d_all;
    head_p = head_d_all;
else
    Y = dBE_eV_A;
    D_p = D_p_A_all;
    head_p = head_A_all;
end

%%  Use the descriptors selected by all three methods
Data = zeros(length(dBE_eV)/2,size(Pair,1));
Head = strings(1,size(Pair,1));
for i = 1:size(Pair,1)
    Data_temp = [D_p(:,Pair(i,1)) D_p(:,Pair(i,2))];
    Head_temp = [head_p(Pair(i,1)) head_p(Pair(i,2))];
    Unit_temp = [1, 1];
    Operation = {'*'};
    [Data(:,i), Head(i), Unit_temp_1] = GenBasicFeature(Data_temp, Head_temp, Unit_temp, Operation, 0);
end

%% Check whether a dopant-modified or adsorbate-modified data 
if ~isempty(strfind(file,'adsorbate'))
    Y = [Y(1:26); Y(28:37); Y(39:end)];
    Data = [Data(1:26, :); Data(28:37, :); Data(39:end, :)];
end

%% Refined with l0-norm method

RMSE_L_zero = zeros(size(Pair,1), size(Pair,1));
Rsq_L_zero = zeros(size(Pair,1), size(Pair,1));
Data_L_zero = Data(:,1:size(Pair,1));
min_RMSE = zeros(size(Pair,1), 1);
max_Rsq = zeros(size(Pair,1), 1);

for i = Dimension
    num_des = 1:1:size(Pair,1);
    sz = nchoosek(size(Pair,1),i);
    index = nchoosek(num_des,i);
    coeff = zeros(i+1, sz);
    for j = 1: sz
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

end

figure(1)
subplot(2,1,1);
plot(1:max(Dimension), min_RMSE(Dimension),'bo') 
hold on
plot(1:max(Dimension), min_RMSE(Dimension),'b-')
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
plot([1 -5], [1, -5],'k-')


