% script for linear power system state estimation
% Loadability assesment

% Willem lambrichts
% General and Unified Model of the Power Flow Problem in Multiterminal AC/DC Networks
% IEEE Transactions on Power Systems 
% 10.1109/TPWRS.2024.3378926

clear all;
close all;
clc;

addpath(genpath(pwd))

%% Base values
A_b = 1e5;
V_b= 400;
Y_b = A_b/V_b^2; 

Vdc_b = 800;
Adc_b = A_b;
Ydc_b = Adc_b/Vdc_b^2; 

%% Set the Grid parameters
Grid_para.n_dc = 8;
Grid_para.n_ac = 18;
Grid_para.n_ph = 1;
Grid_para.n_nodes = Grid_para.n_ac*Grid_para.n_ph + Grid_para.n_dc;
Grid_para.V_b = V_b;
Grid_para.Y_b = Y_b;

[Yac, YYLac, YLac, YTac, YYTac, I_bac, Ampacitiesac, y_lx, y_tx, Aac, linedata_ac]  = Ymatrix('linedata_AC.txt',A_b,V_b,[]);
[Ydc, YYLdc, YLdc, YTdc, YYTdc, I_bdc, Ampacitiesdc, y_ihdc, y_idc, Adc, linedata_dc]  = Ymatrix('linedata_DC.txt',Adc_b,Vdc_b,[]);

Ydc = Ydc(19:26,19:26)/2;
YYLdc = YYLdc(19:26,19:26)/2;
YYTdc = YYTdc(19:26,19:26)/2;
Yac = cell2mat(arrayfun(@(x) x*eye(Grid_para.n_ph),Yac,'UniformOutput',false));
YYLac = cell2mat(arrayfun(@(x) x*eye(Grid_para.n_ph),YYLac,'UniformOutput',false));
YYTac = cell2mat(arrayfun(@(x) x*eye(Grid_para.n_ph),YYTac,'UniformOutput',false));
YY = blkdiag(Yac,Ydc);
Grid_para.G = real(YY);
Grid_para.B = imag(YY);
Grid_para.YY = YY;

%% Set the nodes types
idx1.slack = 1;
idx1.pqac = [2:14]';
idx1.pvac = []';

idx1.pdc = [23:26]';
idx1.vdc = []';

idx1.vscac_pq = [16:18]';
idx1.vscac_vq = [15]';

idx1.vscdc_pq = [20:22]';
idx1.vscdc_vq = [19]';

% idx1.vscac_pq = [17,15]';
% idx1.vscac_vq = [18,16]';
% 
% idx1.vscdc_pq = [21,19]';
% idx1.vscdc_vq = [22,20]';


idx = Get_multiphase_Node_indices(idx1,Grid_para);
linedata = [linedata_ac;linedata_dc];
Grid_para = Get_Converter_para(idx1,linedata,Grid_para);


%% Set the filter parameters
Filter_para.R = 0.008*Y_b; %checked
Filter_para.X = 0.04*Y_b;  %checked
Filter_para.IGBT_piecewise = [  0                   0
                                0.04926559627563	0.7
                                2.30625399327864	0.8
                                15.7793399043317	0.85
                                107.547461516782	0.8999
                                735.837403888342	0.9499
                                1588.01477341768	0.9699];
Filter_para.Exclude_losses = 1;


%% Initialize
if Grid_para.n_ph == 3
    simulation_para.E_0 = [repmat([1; exp(-1i*pi*2/3);  exp(1i*pi*2/3) ], Grid_para.n_ac,1 ) ; ones(Grid_para.n_dc,1)];
elseif Grid_para.n_ph == 1
    simulation_para.E_0 = [repmat(1, Grid_para.n_ac,1 ) ; ones(Grid_para.n_dc,1)];
end
    
simulation_para.tol = 1e-7;
simulation_para.n_max = 100;

%% Get the EMTP measurements

% Load data of nodal voltage and power injections from the .mat file
load('data_balanced.mat') %data_balanced.mat OR data_unbalanced_light.mat OR data_unbalanced_strong.mat OR data_unbalanced_strong_wlosses.mat


folder = './Figures';

%% B26_p
E_all = [];
S_all = [];
for i = 0:0.01:1.5
    
    E_star = data.E_star([1:3:54,55:62]);
    S_star = data.S_star([1:3:54,55:62])/10;
    S_star(26) = complex(-i , imag(S_star(26)));

    [E,J,n_iter] = NR_rectangularACDC_1ph_general(Grid_para,Filter_para,S_star,E_star,idx,simulation_para);

    if n_iter < 99
        E_all = [E_all, E];
        S_all = [S_all, E.*conj(YY*E)];
    end

end
f1 = plot_loadability(S_all,E_all,1,26,'real');
saveas(f1,[folder filesep() 'jpg' filesep() 'Fig8(e)-Loadability'],'jpg');
saveas(f1,[folder filesep() 'eps' filesep() 'Fig8(e)-Loadability'],'epsc');


%% BIC2_p
E_all = [];
S_all = [];
for i = 0:0.01:1.2

    E_star = data.E_star([1:3:54,55:62]);
    S_star = data.S_star([1:3:54,55:62])/10;
    S_star(16) = complex(-i , imag(S_star(16)));

    [E,J,n_iter] = NR_rectangularACDC_1ph_general(Grid_para,Filter_para,S_star,E_star,idx,simulation_para);

    if n_iter < 99
        E_all = [E_all, E];
        S_all = [S_all, E.*conj(YY*E)];
    end

end
f2 = plot_loadability(S_all,E_all,1,16,'real');
saveas(f2,[folder filesep() 'jpg' filesep() 'Fig8(a)-Loadability'],'jpg');
saveas(f2,[folder filesep() 'eps' filesep() 'Fig8(a)-Loadability'],'epsc');

%% BIC2_q
E_all = [];
S_all = [];
for i = -3.5:0.01:2.5

    E_star = data.E_star([1:3:54,55:62]);
    S_star = data.S_star([1:3:54,55:62])/10;
    S_star(16) = complex(real(S_star(16)), -i);

    [E,J,n_iter] = NR_rectangularACDC_1ph_general(Grid_para,Filter_para,S_star,E_star,idx,simulation_para);

    if n_iter < 99
        E_all = [E_all, E];
        S_all = [S_all, E.*conj(YY*E)];
    end

end
f3 = plot_loadability(S_all,E_all,1,16,'imag');
saveas(f3,[folder filesep() 'jpg' filesep() 'Fig8(b)-Loadability'],'jpg');
saveas(f3,[folder filesep() 'eps' filesep() 'Fig8(b)-Loadability'],'epsc');


%% B14_p
E_all = [];
S_all = [];
for i = 0:0.01:2

    E_star = data.E_star([1:3:54,55:62]);
    S_star = data.S_star([1:3:54,55:62])/10;
    S_star(14) = complex(-i , imag(S_star(14)));

    [E,J,n_iter] = NR_rectangularACDC_1ph_general(Grid_para,Filter_para,S_star,E_star,idx,simulation_para);

    if n_iter < 99
        E_all = [E_all, E];
        S_all = [S_all, E.*conj(YY*E)];
    end

end
f4 = plot_loadability(S_all,E_all,1,14,'real');
saveas(f4,[folder filesep() 'jpg' filesep() 'Fig8(c)-Loadability'],'jpg');
saveas(f4,[folder filesep() 'eps' filesep() 'Fig8(c)-Loadability'],'epsc');


%% B14_q
E_all = [];
S_all = [];
for i = -4.5:0.01:3

    E_star = data.E_star([1:3:54,55:62]);
    S_star = data.S_star([1:3:54,55:62])/10;
    S_star(14) = complex(real(S_star(14)), -i);

    [E,J,n_iter] = NR_rectangularACDC_1ph_general(Grid_para,Filter_para,S_star,E_star,idx,simulation_para);

    if n_iter < 99
        E_all = [E_all, E];
        S_all = [S_all, E.*conj(YY*E)];
    end

end
f5 = plot_loadability(S_all,E_all,1,14,'imag');

saveas(f5,[folder filesep() 'jpg' filesep() 'Fig8(d)-Loadability'],'jpg');
saveas(f5,[folder filesep() 'eps' filesep() 'Fig8(d)-Loadability'],'epsc');



