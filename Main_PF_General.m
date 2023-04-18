% main script for linear power system state estimation

clear all;
close all;
clc;

%% Simulation options
addpath(genpath('/acdc-power-flow'))

%% Import paths

% addpath '/Users/willem/Documents/phd/Micro grid specs/'
% addpath('/Users/willem/Documents/phd/State_estimation')
% VSI_name = '/Users/willem/Documents/phd/Micro grid specs/DCdata.txt';
% DC_grid_name = '/Users/willem/Documents/phd/Micro grid specs/DCgrid_data2.txt';
% linedata_DCgrid = '/Users/willem/Documents/phd/Micro grid specs/linedata_DCgrid.txt';
% folder= '/Users/willem/Documents/phd/State_estimation/matlab_Full_SE_3phase/bad data rejection';
% addpath '/Users/willem/Documents/phd/Optimal_power_flow/Sensitivity_codes'

   
    
%% Generalized LF algorithm
% clear all;
% close all;
% clc;

% Base values
A_b = 1e5;
V_b= 400;
Y_b = A_b/V_b^2; 

Vdc_b = 800;
Adc_b = A_b;
Ydc_b = Adc_b/Vdc_b^2; 



% Set the Grid parameters
Grid_para.n_dc = 8;
Grid_para.n_ac = 18;
Grid_para.n_ph = 3;
Grid_para.n_nodes = Grid_para.n_ac*Grid_para.n_ph + Grid_para.n_dc;
Grid_para.V_b = V_b;
Grid_para.Y_b = Y_b;

[Yac, YYL, YL, YT, YYT, I_b, Ampacities, y_lx, y_tx, A, linedata_ac]  = Ymatrix('linedata_AC.txt',A_b,V_b,[]);
[Ydc, YYLdc, YLdc, YT_dc, YYTdc, I_bdc, Ampacitiesdc, y_ihdc, y_idc, Adc, linedata_dc]  = Ymatrix('linedata_DC.txt',Adc_b,Vdc_b,[]);

Ydc = Ydc(19:26,19:26)/2;
Yac = cell2mat(arrayfun(@(x) x*eye(Grid_para.n_ph),Yac,'UniformOutput',false));
YY = blkdiag(Yac,Ydc);
Grid_para.G = real(YY);
Grid_para.B = imag(YY);


% Set the nodes types
idx1.slack = 1;
idx1.pqac = [2:14]';
idx1.pvac = []';

idx1.pdc = [23:26]';
idx1.vdc = []';

idx1.vscac_pq = [17,15]';
idx1.vscac_vq = [18,16]';

idx1.vscdc_pq = [21,19]';
idx1.vscdc_vq = [22,20]';

idx = Get_multiphase_Node_indices(idx1,Grid_para);
linedata = [linedata_ac;linedata_dc];
Grid_para = Get_Converter_para(idx1,linedata,Grid_para);


% Get the EMTP measurements
repeat=1;
ZIN_polyphase = [];
[Nodal_V_mag,Nodal_V_angle, Nodal_I_mag, Nodal_I_angle, Flow_I_mag, Flow_I_angle, Idc_flow, Idc_inj, Vdc_LF, n_timesteps, Mreal, Mimag, Mabs,Nodal_P,Nodal_Q,Pdc_inj,M_LF] = GetEMTPdata_LF(A_b,V_b,Adc_b, Vdc_b,repeat,ZIN_polyphase,Grid_para.n_ph);

Nodal_V_mag = Nodal_V_mag(end,:);
Nodal_V_angle =  Nodal_V_angle(end,:);

Nodal_P =  Nodal_P(end,:);
Nodal_Q =  Nodal_Q(end,:);
Pdc_inj =  Pdc_inj(end,:);

V_complex_LF = transpose(complex(Nodal_V_mag.*cos(Nodal_V_angle), Nodal_V_mag.*sin(Nodal_V_angle)));
Vdc_LF = transpose(Vdc_LF);

E_star = [V_complex_LF(1:end,1); Vdc_LF(:,1)];
S_star = E_star.*conj(YY*E_star); %[transpose(complex(Nodal_P, Nodal_Q)); transpose(complex(Pdc_inj,0))];

% Set the filter parameters
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


% Initialize
E_0 = [repmat([1; exp(-1i*pi*2/3);  exp(1i*pi*2/3) ], Grid_para.n_ac,1 ) ; ones(Grid_para.n_dc,1)];

% Solve the PF
tol = 1e-7;
n_max = 100;
[E,J,n_iter] = NR_rectangularACDC_3ph_general(Grid_para,Filter_para,S_star,E_star,E_0,idx,tol,n_max);

% Visualize the solution

E_delta = E - E_star;
S_delta = E.*conj(YY*E) - S_star;

figure (1)
subplot(2,1,1)
scatter(1:62,real(E_delta))
ylabel('Delta E - real')

subplot(2,1,2)
scatter(1:62,imag(E_delta))
ylabel('Delta E - imaginary')

% figure (2)
% subplot(2,1,1)
% scatter(1:62,real(S_delta))
% ylabel('Delta S - real')
% 
% subplot(2,1,2)
% scatter(1:62,imag(S_delta))
% ylabel('Delta S - imaginary')

