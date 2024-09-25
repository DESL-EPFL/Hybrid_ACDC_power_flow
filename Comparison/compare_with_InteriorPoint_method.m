% script for linear power system state estimation
% Coparison of the proposed Newton Raphson based load flow
% to an interior point method

% Willem lambrichts
% General and Unified Model of the Power Flow Problem in Multiterminal AC/DC Networks
% IEEE Transactions on Power Systems 
% 10.1109/TPWRS.2024.3378926


% time NR: 0.05 s
% time IP: 2.5 s
% accuracy in order of 10e-8

clear all;
close all;
clc;

addpath(genpath(pwd))


%% initialize parameters
nb_phases = 1;
[Grid_para,idx,Filter_para,simulation_para] = initialize_InteriorPoint(nb_phases);

%% Get the EMTP measurements

% Load data of nodal voltage and power injections from the .mat file
load('data_balanced.mat') %data_balanced.mat OR data_unbalanced_light.mat OR data_unbalanced_strong.mat OR data_unbalanced_strong_wlosses.mat
E_star = data.E_star;
S_star = data.S_star;

if Grid_para.n_ph == 1
    id_v = [1: 3 : Grid_para.n_ac*3,Grid_para.n_ac*3+1:Grid_para.n_ac*3 + Grid_para.n_dc ];

    E_star = E_star(id_v);
    S_star = [transpose(sum(reshape(S_star(1:Grid_para.n_ac*3),3,Grid_para.n_ac))/3); S_star(Grid_para.n_ac*3+1:end)];
end


%% Newton-Raphson based load flow

[E_PF,J,n_iter,time_NR] = NR_rectangularACDC_1ph_general(Grid_para,Filter_para,S_star,E_star,idx,simulation_para);
S_PF = E_PF.*conj(Grid_para.YY*E_PF);


%% Interior point base load flow
E_abs = sdpvar(length(simulation_para.E_0),1, 'full');
E_arg = sdpvar(length(simulation_para.E_0),1, 'full');
E = E_abs .* exp(1j*E_arg);

temp = E.*conj(Grid_para.YY*E);

time_IP = tic;
con = [];
%PQ nodes
con = [con (real(temp(idx.pqac)) == real(S_star(idx.pqac)))     ];
con = [con (imag(temp(idx.pqac)) == imag(S_star(idx.pqac)))     ];
con = [con (real(temp(idx.pdc)) == real(S_star(idx.pdc)))     ];
con = [con (imag(temp(idx.pdc)) == imag(S_star(idx.pdc)))     ];

%IC nodes
con = [con (real(temp(idx.vscac_vq)) == -real(temp(idx.vscdc_vq))) ];
con = [con (imag(temp(idx.vscac_vq)) == imag(S_star(idx.vscac_vq))) ];


%slack node
con = [con (E_abs(1) == abs(E_star(1)) )    ];
con = [con (E_arg(1) == angle(E_star(1)) )    ];

%DC nodes
con = [con (E_abs(idx.vscdc_vq) == real(E_star(idx.vscdc_vq)) )    ];
con = [con (E_arg(idx.vscdc_vq) == 0 )    ];
con = [con (E_arg(idx.pdc) == 0)];


con = [con (E_abs <= 1.05)&(E_abs >= 0.90)]; 

options = sdpsettings('verbose', 0);

sol = optimize(con, 1, options);

time_IP = toc(time_IP);

disp(['Newton-Raphson based load flow: ', num2str(time_NR), ' secondes'])
disp(['Interior point base load flow: ', num2str(time_IP), ' secondes'])

% mean(double(E) - E_PF)
% max(double(E) - E_PF)
% 
% mean(double(E.*conj(Grid_para.YY*E)) - S_PF)
% max(double(E.*conj(Grid_para.YY*E)) - S_PF)
