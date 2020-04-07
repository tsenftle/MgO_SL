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

% Department of Chemical and Biomolecular Engineering, Rice University, Houston, TX 77005 (USA)
% *tsenftle@rice.edu

%%
clear all;
clc
disp('Importing Data');

%% Import data: binding energy
[~, ~, raw] = xlsread('BindingEnergy_adsorbate_V1.xlsx','binding_energy_difference');
raw = raw(2:74,:);
cellVectors = raw(:,[1,2]);
raw = raw(:,3);
data = reshape([raw{:}],size(raw));

%% Allocate imported array to column variable names
metal_name  =   cellVectors(:,1);      % Adamtom metal name
AdsDop      =   cellVectors(:,2);      % Adsorbate or Dopant name
dBE_eV_A    =   data;                  % Binding energy difference between ...
                                       % clean surface and modified surface (eV)

%% Clear temporary variables
clearvars raw cellVectors data;

%% Import data: atomic properties
[~, ~, raw] = xlsread('BindingEnergy_adsorbate_V1.xlsx','atomic_properties');
raw = raw(3:23,:);
cellVectors = raw(:,1);
raw = raw(1:21,2:17);
data = reshape([raw{:}], size(raw));

%% Allocate imported array to column variable names
atom_name = cellVectors(:,1);   % atom name
Z = data(:, 1);                 % Atomic number
EN_P = data(:,2);               % Electronegativity - Pauling scale
EN_MB = data(:,3);              % Electronegativity - Martynov-Batsanov
IE_1_eV = data(:,4);            % 1st ionization energy (eV)
IE_2_eV = data(:,5);            % 2nd ionization energy (eV)
EA_eV = data(:,6);              % Electron affinity (eV)
Hs_eV = data(:,7)/60.2/1.602;   % Heat of sublimination of adatom metal
Hf_eV = data(:,8)/60.2/1.602;   % Heat of oxide formation of most stable metal oxide
dHo_eV = data(:,9);             % Heat of oxide formation of adatom metal oxide
Zunger_rs = data(:,10);         % Zunger orbital radii, rs
Zunger_rp = data(:,11);         % Zunger orbital radii, rp
WC_rs = data(:,12);             % Waber and Cromer radii, rs
WC_rp = data(:,13);             % Waber and Cromer radii, rp
Val_elec = data(:,14);          % Number of valence electrons
nws13 = data(:,15);             % Miedema - 1st parameter
phi = data(:,16);               % Miedema - 2nd parameter

%% Clear temporary variables
clearvars raw cellVectors data;

%% Import data: adsorbate properties
[~, ~, raw] = xlsread('BindingEnergy_adsorbate_V1.xlsx','adsorbate_properties');
raw = raw(3:6,:);
cellVectors = raw(:,1);
raw = raw(1:4,2:10);
data = reshape([raw{:}], size(raw));

%% Allocate imported array to column variable names
adsorbate_name = cellVectors(:,1);   % Adsorbate name
A_IE_eV = data(:,1);                 % Ionization energy (eV)
A_EA_eV = data(:,2);                 % Electron affinity (eV)
A_Val_elec = data(:,3);              % Number of valence electrons
A_cn = data(:,4);                    % Coordination number of adsorbate on surface
A_DE_eV = data(:,5)/60.2/1.602;      % Bond disscociation energy (eV)
A_EN_P = data(:,6);                  % Electronegativity - Pauling scale of attached atom
A_EN_P_del = data(:,7);              % Electronegativity - Pauling scale - difference between adsorbate atoms
A_EN_MB = data(:,8);                 % Electronegativity - Martynov-Batsanov of attached atom
A_EN_MB_del = data(:,9);             % Electronegativity - Martynov-Batsanov - difference between adsorbate atoms

%% Clear temporary variables
clearvars raw cellVectors data;

% END OF IMPORT

%% Create D matrix

n_s  = size(dBE_eV_A,1);          % Number of samples

% Atomic properties
head_a = {'m_Z', 's_Z', 'o_Z', 'm_EN_P', 's_EN_P','o_EN_P', 'm_EN_MB', 's_EN_MB','o_EN_MB',...
          'm_IE_1_eV', 'm_IE_2_eV', 's_IE_1_eV', 's_IE_2_eV', 'o_IE_1_eV', 'o_IE_2_eV', ...
          'm_EA_eV', 's_EA_eV', 'o_EA_eV', 'm_Hs_eV', 's_Hs_eV','o_Hs_eV',...
          'm_Hf_eV', 's_Hf_eV', 'm_dHo_eV', 's_dHo_eV','m_Zunger_rs', 'm_Zunger_rp', ...
          's_Zunger_rs', 's_Zunger_rp', 'o_Zunger_rs', 'o_Zunger_rp',...
          'm_WC_rs', 'm_WC_rp', 's_WC_rs', 's_WC_rp', 'o_WC_rs', 'o_WC_rp',...
          'm_Val_elec', 's_Val_elec', 'o_Val_elec', 'm_nws13', 's_nws13', 'm_phi', 's_phi'};
% head_1: title of each primary descriptors from atomic properties
% m_*: properties of adsorbed metals
% s_*: properties of parent metal in oxide surface
% o_*: properties of oxygen
n_ap  = size(head_a, 2);         % Number of primary descriptors from atomic properties
D_p_a = zeros(n_s,n_ap);         % Create a matrix filled with atomic properties

% Adsorbate properties
head_A = {'a_IE_eV', 'a_EA_eV', 'a_Val_elec', 'a_cn', 'a_DE_eV',...
          'a_EN_P', 'a_EN_P_del', 'a_EN_MB', 'a_EN_MB_del'};
% head_3: title of each primary descriptors from adsorbate properties
% a_*: properties of dopants
n_Ap  = size(head_A, 2);         % Number of primary descriptors from adsorbate properties
D_p_A = zeros(n_s, n_Ap);        % Create a matrix filled with adsorbate properties

% Create primary descriptor matrix
% Atomic properties
for i=1:n_s
    idm=find(strcmp(atom_name, metal_name(i)));      % index of adsorbed metals
    ids=find(strcmp(atom_name, 'Mg'));               % index of parent metal in oxide surface
    ido=find(strcmp(atom_name, 'O'));                % index of oxygen 
    D_p_a(i,:)= [ Z(idm) Z(ids) Z(ido) EN_P(idm) EN_P(ids) EN_P(ido) ...
                 EN_MB(idm) EN_MB(ids) EN_MB(ido) IE_1_eV(idm) IE_2_eV(idm) ...
                 IE_1_eV(ids) IE_2_eV(ids) IE_1_eV(ido) IE_2_eV(ido) ...
                 EA_eV(idm) EA_eV(ids) EA_eV(ido) Hs_eV(idm) Hs_eV(ids) Hs_eV(ido) ...
                 Hf_eV(idm) Hf_eV(ids) dHo_eV(idm) dHo_eV(ids) Zunger_rs(idm) ...
                 Zunger_rp(idm) Zunger_rs(ids) Zunger_rp(ids) Zunger_rs(ido) Zunger_rp(ido) ...
                 WC_rs(idm) WC_rp(idm) WC_rs(ids) WC_rp(ids) WC_rs(ido) WC_rp(ido) ...
                 Val_elec(idm) Val_elec(ids) Val_elec(ido) nws13(idm) nws13(ids) ...
                 phi(idm) phi(ids)];
end

% Adsorbate properties
for i=1:n_s
    ida=find(strcmp(adsorbate_name, AdsDop(i)));     % index of adsorbate
    D_p_A(i,:)= [ A_IE_eV(ida) A_EA_eV(ida) A_Val_elec(ida) A_cn(ida)...
                  A_DE_eV(ida) A_EN_P(ida) A_EN_P_del(ida) A_EN_MB(ida) A_EN_MB_del(ida)];
end

% End of data import

%% Generate matrix with dopant properties
disp('Generating Feature Set of Atoms');

%% Baseline: simple secondary descriptors of Atoms

% Simple secondary descriptors for atoms
Data_temp = D_p_a;
Head_temp = head_a;
Unit_temp = ones(1, size(Data_temp,2));
Operation_1 = {'^r', '^I', '^2', '^3', 'abs'};
[Data_temp_1, Head_temp_1, Unit_temp_1] = GenSelfFeature(Data_temp, Head_temp, Unit_temp, Operation_1, 1);
D_p_a = Data_temp_1;
head_a = Head_temp_1;
clearvars Data_temp Head_temp Unit_temp Data_temp_1 Head_temp_1 Unit_temp_1 Operation_1

%% Generate matrix with adsorbate properties
disp('Generating Feature Set of Adsorbates');

% Simple secondary descriptors for dopants
Data_temp = D_p_A;
Head_temp = head_A;
Unit_temp = ones(1, size(Data_temp,2));
Operation_1 = {'^r', '^I', '^2', '^3', 'abs'};
[Data_temp_1, Head_temp_1, Unit_temp_1] = GenSelfFeature(Data_temp, Head_temp, Unit_temp, Operation_1, 1);
D_p_A = Data_temp_1;
head_A = Head_temp_1;
clearvars Data_temp Head_temp Unit_temp Data_temp_1 Head_temp_1 Unit_temp_1 Operation_1

%% Atomic number: metal, surface and oxygen
Data_temp = D_p_a(:,1:3);
Head_temp = head_a(1:3);
Unit_temp = ones(1, size(Data_temp,2));
Operation_1 = {'-', '-abs', '/abs','*abs'};
[Data_temp_1, Head_temp_1, Unit_temp_1] = GenBasicFeature(Data_temp, Head_temp, Unit_temp, Operation_1, 0);
Operation_2 = {'^r'; '^I'; '^2'};
[Data_temp_2, Head_temp_2, Unit_temp_2] = GenSelfFeature(Data_temp_1, Head_temp_1, Unit_temp_1, Operation_2, 1);
D_p_a_Z = Data_temp_2;
head_a_Z = Head_temp_2;
clearvars Data_temp Head_temp Unit_temp Data_temp_1 Head_temp_1 Unit_temp_1 ...
          Data_temp_2 Head_temp_2 Unit_temp_2 Operation_1 Operation_2

%% Electronegativity - Pauling scale: metal, surface, oxygen and adsorbate
Data_temp = [D_p_a(:,4:6) D_p_A(:,6:7)];
Head_temp = [head_a(4:6) head_A(6:7)];
Unit_temp = ones(1, size(Data_temp,2));
Operation_1 = {'-', '-abs'};
[Data_temp_1, Head_temp_1, Unit_temp_1] = GenBasicFeature(Data_temp, Head_temp, Unit_temp, Operation_1, 1);
Operation_2 = {'/abs'};
[Data_temp_2, Head_temp_2, Unit_temp_2] = GenBasicFeature(Data_temp_1, Head_temp_1, Unit_temp_1, Operation_2, 1);
Operation_3 = {'^r'; '^I'; '^2'};
[Data_temp_3, Head_temp_3, Unit_temp_3] = GenSelfFeature(Data_temp_2, Head_temp_2, Unit_temp_2, Operation_3, 1);
D_p_A_EN_P = Data_temp_3;
head_A_EN_P = Head_temp_3;
clearvars Data_temp Head_temp Unit_temp Data_temp_1 Head_temp_1 Unit_temp_1 ...
          Data_temp_2 Head_temp_2 Unit_temp_2 Data_temp_3 Head_temp_3 Unit_temp_3 ...
          Operation_1 Operation_2 Operation_3
      
%% Electronegativity - Martynov-Batsanov scale: metal, surface, oxygen and adsorbate
Data_temp = [D_p_a(:,7:9) D_p_A(:,8:9)];
Head_temp = [head_a(7:9) head_A(8:9)];
Unit_temp = ones(1, size(Data_temp,2));
Operation_1 = {'-', '-abs'};
[Data_temp_1, Head_temp_1, Unit_temp_1] = GenBasicFeature(Data_temp, Head_temp, Unit_temp, Operation_1, 1);
Operation_2 = {'/abs'};
[Data_temp_2, Head_temp_2, Unit_temp_2] = GenBasicFeature(Data_temp_1, Head_temp_1, Unit_temp_1, Operation_2, 1);
Operation_3 = {'^r'; '^I'; '^2'};
[Data_temp_3, Head_temp_3, Unit_temp_3] = GenSelfFeature(Data_temp_2, Head_temp_2, Unit_temp_2, Operation_3, 1);
D_p_A_EN_MB = Data_temp_3;
head_A_EN_MB = Head_temp_3;
clearvars Data_temp Head_temp Unit_temp Data_temp_1 Head_temp_1 Unit_temp_1 ...
          Data_temp_2 Head_temp_2 Unit_temp_2 Data_temp_3 Head_temp_3 Unit_temp_3 ...
          Operation_1 Operation_2 Operation_3
      
%% Ionization energies: metal, surface, oxygen and adsorbate
Data_temp = [D_p_a(:,10:15) D_p_A(:,1)];
Head_temp = [head_a(10:15) head_A(1)];
Unit_temp = ones(1, size(Data_temp,2));
Operation_1 = {'-', '-abs'};
[Data_temp_1, Head_temp_1, Unit_temp_1] = GenBasicFeature(Data_temp, Head_temp, Unit_temp, Operation_1, 0);
[Data_temp_2, Head_temp_2, Unit_temp_2] = RmVarZero(Data_temp, Head_temp, Unit_temp, 10^-10);
[Data_temp_3, Head_temp_3, Unit_temp_3] = RmVarZero(Data_temp_1, Head_temp_1, Unit_temp_1, 10^-10);
Data_temp_3 = [Data_temp_2 Data_temp_3];
Head_temp_3 = [Head_temp_2 Head_temp_3];
Unit_temp_3 = [Unit_temp_2 Unit_temp_3];
Operation_2 = {'/abs'};
[Data_temp_4, Head_temp_4, Unit_temp_4] = GenBasicFeature(Data_temp_3, Head_temp_3, Unit_temp_3, Operation_2, 1);
Data_temp_4 = [Data_temp Data_temp_4];
Head_temp_4 = [Head_temp Head_temp_4];
Unit_temp_4 = [Unit_temp Unit_temp_4];
Operation_3 = {'^I'};
[Data_temp_5, Head_temp_5, Unit_temp_5] = GenSelfFeature(Data_temp_4, Head_temp_4, Unit_temp_4, Operation_3, 1);
D_p_A_IE = Data_temp_5;
head_A_IE = Head_temp_5;
clearvars Data_temp Head_temp Unit_temp Data_temp_1 Head_temp_1 Unit_temp_1 ...
          Data_temp_2 Head_temp_2 Unit_temp_2 Data_temp_3 Head_temp_3 Unit_temp_3 ...
          Data_temp_4 Head_temp_4 Unit_temp_4 Data_temp_5 Head_temp_5 Unit_temp_5 ...
          Operation_1 Operation_2 Operation_3

%% Electron affinity: metal, surface, oxygen and adsorbate
Data_temp = [D_p_a(:,16:18) D_p_A(:,2)];
Head_temp = [head_a(16:18) head_A(2)];
Unit_temp = ones(1, size(Data_temp,2));
Operation_1 = {'-', '-abs'};
[Data_temp_1, Head_temp_1, Unit_temp_1] = GenBasicFeature(Data_temp, Head_temp, Unit_temp, Operation_1, 0);
[Data_temp_2, Head_temp_2, Unit_temp_2] = RmVarZero(Data_temp, Head_temp, Unit_temp, 10^-10);
[Data_temp_3, Head_temp_3, Unit_temp_3] = RmVarZero(Data_temp_1, Head_temp_1, Unit_temp_1, 10^-10);
Data_temp_3 = [Data_temp_2 Data_temp_3];
Head_temp_3 = [Head_temp_2 Head_temp_3];
Unit_temp_3 = [Unit_temp_2 Unit_temp_3];
Operation_2 = {'/abs'};
[Data_temp_4, Head_temp_4, Unit_temp_4] = GenBasicFeature(Data_temp_3, Head_temp_3, Unit_temp_3, Operation_2, 1);
Data_temp_4 = [Data_temp Data_temp_4];
Head_temp_4 = [Head_temp Head_temp_4];
Unit_temp_4 = [Unit_temp Unit_temp_4];
Operation_3 = {'^I'};
[Data_temp_5, Head_temp_5, Unit_temp_5] = GenSelfFeature(Data_temp_4, Head_temp_4, Unit_temp_4, Operation_3, 1);
D_p_A_EA = Data_temp_5;
head_A_EA = Head_temp_5;
clearvars Data_temp Head_temp Unit_temp Data_temp_1 Head_temp_1 Unit_temp_1 ...
          Data_temp_2 Head_temp_2 Unit_temp_2 Data_temp_3 Head_temp_3 Unit_temp_3 ...
          Data_temp_4 Head_temp_4 Unit_temp_4 Data_temp_5 Head_temp_5 Unit_temp_5 ...
          Operation_1 Operation_2 Operation_3

%% Heat of oxide formation of adatom metal oxide: metal and surface 
Data_temp = D_p_a(:,24:25);
Head_temp = head_a(24:25);
Unit_temp = ones(1, size(Data_temp,2));
Operation_1 = {'-', '-abs'};
[Data_temp_1, Head_temp_1, Unit_temp_1] = GenBasicFeature(Data_temp, Head_temp, Unit_temp, Operation_1, 1);
Operation_2 = {'/abs'};
[Data_temp_2, Head_temp_2, Unit_temp_2] = GenBasicFeature(Data_temp_1, Head_temp_1, Unit_temp_1, Operation_2, 1);
Operation_3 = {'^I'};
[Data_temp_3, Head_temp_3, Unit_temp_3] = GenSelfFeature(Data_temp_2, Head_temp_2, Unit_temp_2, Operation_3, 1);
D_p_a_dHo = Data_temp_3;
head_a_dHo = Head_temp_3;
clearvars Data_temp Head_temp Unit_temp Data_temp_1 Head_temp_1 Unit_temp_1 ...
          Data_temp_2 Head_temp_2 Unit_temp_2 Data_temp_3 Head_temp_3 Unit_temp_3 ...
          Operation_1 Operation_2 Operation_3

%% Zunger orbital radii, rs and rp: metal, surface and oxygen
Data_temp = D_p_a(:,26:31);
Head_temp = head_a(26:31);
Unit_temp = ones(1, size(Data_temp,2));
Operation_1 = {'-', '-abs'};
[Data_temp_1, Head_temp_1, Unit_temp_1] = GenBasicFeature(Data_temp, Head_temp, Unit_temp, Operation_1, 1);
[Data_temp_2, Head_temp_2, Unit_temp_2] = RmVarZero(Data_temp_1, Head_temp_1, Unit_temp_1, 10^-10);
Operation_2 = {'^r'; '^I'; '^2'};
[Data_temp_3, Head_temp_3, Unit_temp_3] = GenSelfFeature(Data_temp_2, Head_temp_2, Unit_temp_2, Operation_2, 1);
D_p_a_Z_r = Data_temp_3;
head_a_Z_r = Head_temp_3;
clearvars Data_temp Head_temp Unit_temp Data_temp_1 Head_temp_1 Unit_temp_1 ...
          Data_temp_2 Head_temp_2 Unit_temp_2 Data_temp_3 Head_temp_3 Unit_temp_3 ...
          Operation_1 Operation_2

%% Waber and Cromer orbital radii, rs and rp: metal, surface and oxygen
Data_temp = D_p_a(:,32:37);
Head_temp = head_a(32:37);
Unit_temp = ones(1, size(Data_temp,2));
Operation_1 = {'-', '-abs'};
[Data_temp_1, Head_temp_1, Unit_temp_1] = GenBasicFeature(Data_temp, Head_temp, Unit_temp, Operation_1, 1);
[Data_temp_2, Head_temp_2, Unit_temp_2] = RmVarZero(Data_temp_1, Head_temp_1, Unit_temp_1, 10^-10);
Operation_2 = {'^r'; '^I'; '^2'};
[Data_temp_3, Head_temp_3, Unit_temp_3] = GenSelfFeature(Data_temp_2, Head_temp_2, Unit_temp_2, Operation_2, 1);
D_p_a_WC_r = Data_temp_3;
head_a_WC_r = Head_temp_3;
clearvars Data_temp Head_temp Unit_temp Data_temp_1 Head_temp_1 Unit_temp_1 ...
          Data_temp_2 Head_temp_2 Unit_temp_2 Data_temp_3 Head_temp_3 Unit_temp_3 ...
          Operation_1 Operation_2

%% Number of valence electrons: metal, surface, oxygen and adsorbate
Data_temp = [D_p_a(:,38:40) D_p_A(:,3)];
Head_temp = [head_a(38:40) head_A(3)];
Unit_temp = ones(1, size(Data_temp,2));
Operation_1 = {'-', '-abs'};
[Data_temp_1, Head_temp_1, Unit_temp_1] = GenBasicFeature(Data_temp, Head_temp, Unit_temp, Operation_1, 1);
Operation_2 = {'^r'; '^I'; '^2'};
[Data_temp_2, Head_temp_2, Unit_temp_2] = GenSelfFeature(Data_temp_1, Head_temp_1, Unit_temp_1, Operation_2, 1);
D_p_A_NVal = Data_temp_2;
head_A_NVal = Head_temp_2;
clearvars Data_temp Head_temp Unit_temp Data_temp_1 Head_temp_1 Unit_temp_1 ...
          Data_temp_2 Head_temp_2 Operation_1 Operation_2 Operation_3

%% Miedema - 1st parameter: metal and surface
Data_temp = D_p_a(:,41:42);
Head_temp = head_a(41:42);
Unit_temp = ones(1, size(Data_temp,2));
Operation_1 = {'-', '-abs'};
[Data_temp_1, Head_temp_1, Unit_temp_1] = GenBasicFeature(Data_temp, Head_temp, Unit_temp, Operation_1, 1);
Operation_2 = {'^r'; '^I'; '^2'};
[Data_temp_2, Head_temp_2, Unit_temp_2] = GenSelfFeature(Data_temp_1, Head_temp_1, Unit_temp_1, Operation_2, 1);
D_p_a_M_1 = Data_temp_2;
head_a_M_1 = Head_temp_2;
clearvars Data_temp Head_temp Unit_temp Data_temp_1 Head_temp_1 Unit_temp_1 ...
          Data_temp_2 Head_temp_2 Operation_1 Operation_2 Operation_3
      
%% Miedema - 2nd parameter: metal and surface
Data_temp = D_p_a(:,43:44);
Head_temp = head_a(43:44);
Unit_temp = ones(1, size(Data_temp,2));
Operation_1 = {'-', '-abs'};
[Data_temp_1, Head_temp_1, Unit_temp_1] = GenBasicFeature(Data_temp, Head_temp, Unit_temp, Operation_1, 1);
Operation_2 = {'^r'; '^I'; '^2'};
[Data_temp_2, Head_temp_2, Unit_temp_2] = GenSelfFeature(Data_temp_1, Head_temp_1, Unit_temp_1, Operation_2, 1);
D_p_a_M_2 = Data_temp_2;
head_a_M_2 = Head_temp_2;
clearvars Data_temp Head_temp Unit_temp Data_temp_1 Head_temp_1 Unit_temp_1 ...
          Data_temp_2 Head_temp_2 Operation_1 Operation_2 Operation_3
      
%% Absolute Electronegativity: metal, surface, oxygen and adsorbate
Data_temp = [D_p_a(:,[10 12 14 16:18]) D_p_A(:,1:2)];
Head_temp = [head_a([10 12 14 16:18]) head_A(1:2)];
Unit_temp = [1 2 3 1 2 3 4 4];
Operation_1 = {'+'};
[Data_temp_1, Head_temp_1, Unit_temp_1] = GenBasicFeature(Data_temp, Head_temp, Unit_temp, Operation_1, 0);
Unit_temp_1 = ones(1, size(Data_temp_1,2));
Operation_2 = {'-', '-abs'};
[Data_temp_2, Head_temp_2, Unit_temp_2] = GenBasicFeature(Data_temp_1, Head_temp_1, Unit_temp_1, Operation_2, 1);
Operation_3 = {'/abs'};
[Data_temp_3, Head_temp_3, Unit_temp_3] = GenBasicFeature(Data_temp_2, Head_temp_2, Unit_temp_2, Operation_3, 1);
Operation_4 = {'^I'};
[Data_temp_4, Head_temp_4, Unit_temp_4] = GenSelfFeature(Data_temp_3, Head_temp_3, Unit_temp_3, Operation_4, 1);
D_p_a_AEN = Data_temp_4;
head_a_AEN = Head_temp_4;
clearvars Data_temp Head_temp Unit_temp Data_temp_1 Head_temp_1 Unit_temp_1 ...
          Data_temp_2 Head_temp_2 Unit_temp_2 Data_temp_3 Head_temp_3 Unit_temp_3 ...
          Data_temp_4 Head_temp_4 Unit_temp_4 Operation_1 Operation_2 Operation_3 ...
          Operation_4

%% Bond Dissociation Energy and Coordination Number
Data_temp = D_p_A(:,4:5);
Head_temp = head_A(4:5);
Unit_temp = ones(1, size(Data_temp,2));
Operation_1 = {'*'};
[Data_temp_1, Head_temp_1, Unit_temp_1] = GenBasicFeature(Data_temp, Head_temp, Unit_temp, Operation_1, 1);
D_p_A_CNDE = Data_temp_1;
head_A_CNDE = Head_temp_1;
clearvars Data_temp Head_temp Unit_temp Data_temp_1 Head_temp_1 Unit_temp_1

%% Add up matrices of adsorbate

D_p_A_all = [ D_p_a D_p_A D_p_a_Z D_p_A_EN_P D_p_A_EN_MB D_p_A_IE D_p_A_EA D_p_a_dHo ...
              D_p_a_Z_r D_p_a_WC_r D_p_A_NVal D_p_a_M_1 D_p_a_M_2 D_p_a_AEN D_p_A_CNDE];
head_A_all = [ head_a head_A head_a_Z head_A_EN_P head_A_EN_MB head_A_IE head_A_EA head_a_dHo ...
               head_a_Z_r head_a_WC_r head_A_NVal head_a_M_1 head_a_M_2 head_a_AEN head_A_CNDE];

%% Split the matrix into test set and training set

rng(0);
idx = randperm(n_s);
Data_A_BE = zeros(length(dBE_eV_A),1);
Data_A_all = D_p_A_all*0;

for i = 1:n_s
    Data_A_BE(i) = dBE_eV_A(idx(i));
    Data_A_all(i,:) = D_p_A_all(idx(i),:);
end

Data_A_train = Data_A_BE(1:52);
Data_A_p_train = Data_A_all(1:52,:);
Data_A_test = Data_A_BE(53:end);
Data_A_p_test = Data_A_all(53:end,:);

save('data_set_adsorbate.mat')



