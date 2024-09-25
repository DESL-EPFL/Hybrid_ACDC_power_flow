% main script for linear power system state estimation

% Willem lambrichts
% General and Unified Model of the Power Flow Problem in Multiterminal AC/DC Networks
% IEEE Transactions on Power Systems 
% 10.1109/TPWRS.2024.3378926

clear all;
close all;
clc;

addpath(genpath(pwd))


%% initialize parameters
nb_phases = 3;
[Grid_para,idx,Filter_para,simulation_para] = initialize(nb_phases);


%% Get the EMTP measurements

%OPTION 1 (prefered) - Load data of nodal voltage and power injections from the .mat file
% load('data_balanced.mat')           %balanced grid
% load('data_unbalanced_light.mat')   %slightly unbalanced grid
% load('data_unbalanced_strong.mat')  %strongly unbalanced grid
% load('data_unbalanced_strong_wlosses.mat')

%OPTION 2 - Get data directly from EMPT simulation 
% emtp_data = load('ACDC_balanced.mat');  %balanced grid
% emtp_data = load('ACDC_unbalanced_light.mat'); %slightly unbalanced grid
% emtp_data = load('ACDC_unbalanced_strong.mat'); %strongly unbalanced grid
% emtp_data = load('ACDC_unbalanced_strong_wlosses.mat'); %strongly unbalanced grid + filter
% data = process_EMTP_data(emtp_data);


%% Balanced - load flow 

%get data
load('data_balanced.mat')
E_star = data.E_star;
S_star = data.S_star;

%solve load flow
[E,~,~] = NR_rectangularACDC_3ph_general(Grid_para,Filter_para,S_star,E_star,idx,simulation_para);

% compute error
E_delta_balanced = E - E_star;
S_delta_balanced = E.*conj( Grid_para.YY*E) - S_star;

%% Unbalanced - load flow 

%get data
load('data_unbalanced_strong.mat')
E_star = data.E_star;
S_star = data.S_star;

%solve load flow
[E,~,~] = NR_rectangularACDC_3ph_general(Grid_para,Filter_para,S_star,E_star,idx,simulation_para);

% compute error
E_delta_unbalanced = E - E_star;
S_delta_unbalanced = E.*conj( Grid_para.YY*E) - S_star;


%% plot (figure 5 in paper)

f = make_plot_Fig5(E_delta_balanced,E_delta_unbalanced,Grid_para);

folder = './Figures';
saveas(f,[folder filesep() 'jpg' filesep() 'Fig5-LoadFlow_Error_Microgrid'],'jpg');
saveas(f,[folder filesep() 'eps' filesep() 'Fig5-LoadFlow_Error_Microgrid'],'epsc');
