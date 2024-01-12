
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

[Yac, YYL, YL, YT, YYT, I_b, Ampacities, y_lx, y_tx, A, linedata_ac]  = Ymatrix('linedata_AC.txt',A_b,V_b,[]);
[Ydc, YYLdc, YLdc, YT_dc, YYTdc, I_bdc, Ampacitiesdc, y_ihdc, y_idc, Adc, linedata_dc]  = Ymatrix('linedata_DC.txt',Adc_b,Vdc_b,[]);

Ydc = Ydc(19:26,19:26)/2;
Yac = cell2mat(arrayfun(@(x) x*eye(Grid_para.n_ph),Yac,'UniformOutput',false));
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

idx1.vscac_pq = []';
idx1.vscac_vq = [15:18]';

idx1.vscdc_pq = []';
idx1.vscdc_vq = [19:22]';

idx = Get_multiphase_Node_indices(idx1,Grid_para);
linedata = [linedata_ac;linedata_dc];
Grid_para = Get_Converter_para(idx1,linedata,Grid_para);

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

% Get data directly from EMPT simulation
% 
% repeat=1;
% ZIN_polyphase = [];
% [Nodal_V_mag,Nodal_V_angle, Nodal_I_mag, Nodal_I_angle, Flow_I_mag, Flow_I_angle, Idc_flow, Idc_inj, Vdc_LF, n_timesteps, Mreal, Mimag, Mabs,Nodal_P,Nodal_Q,Pdc_inj,M_LF] = GetEMTPdata_LF(A_b,V_b,Adc_b, Vdc_b,repeat,ZIN_polyphase,Grid_para.n_ph);
% 
% Nodal_V_mag = Nodal_V_mag(end,:);
% Nodal_V_angle =  Nodal_V_angle(end,:);
% 
% Nodal_P =  Nodal_P(end,:);
% Nodal_Q =  Nodal_Q(end,:);
% Pdc_inj =  Pdc_inj(end,:);
% 
% V_complex_LF = transpose(complex(Nodal_V_mag.*cos(Nodal_V_angle), Nodal_V_mag.*sin(Nodal_V_angle)));
% Vdc_LF = transpose(Vdc_LF);
% 
% E_star = [V_complex_LF(1:end,1); Vdc_LF(:,1)];
% S_star = [transpose(complex(Nodal_P, Nodal_Q)); transpose(complex(Pdc_inj,0))];


%% Set the filter parameters
Filter_para.R = 1e-5 + 0*0.008*Y_b; %checked
Filter_para.X = 1e-5 + 0*0.04*Y_b;  %checked
Filter_para.IGBT_piecewise = [  0                   0
                                0.04926559627563	0.7
                                2.30625399327864	0.8
                                15.7793399043317	0.85
                                107.547461516782	0.8999
                                735.837403888342	0.9499
                                1588.01477341768	0.9699];
Filter_para.Exclude_losses = 1;

%% Initialize
E_0 = [repmat([1], Grid_para.n_ac,1 ) ; ones(Grid_para.n_dc,1)];

if Grid_para.n_ph == 3
    E_0 = [repmat([1; exp(-1i*pi*2/3);  exp(1i*pi*2/3) ], Grid_para.n_ac,1 ) ; ones(Grid_para.n_dc,1)];
end
%% Solve the PF
tol = 1e-7;
n_max = 100;
[E_PF,J,n_iter,time_NR] = NR_rectangularACDC_1ph_general(Grid_para,Filter_para,S_star,E_star,E_0,idx,tol,n_max);
S_PF = E_PF.*conj(YY*E_PF);

time_NR

%% IP method
time_IP_arr = [];

for i = 1:1
    i;

E_abs = sdpvar(length(E_0),1, 'full');
E_arg = sdpvar(length(E_0),1, 'full');

time_IP = tic;
% E = E_re + 1j*E_im;
E = E_abs .* exp(1j*E_arg);

temp = E.*conj(YY*E);

time_IP = tic;
con = [];
%PQ
con = [con (real(temp(idx.pqac)) == real(S_star(idx.pqac)))     ];
con = [con (imag(temp(idx.pqac)) == imag(S_star(idx.pqac)))     ];
con = [con (real(temp(idx.pdc)) == real(S_star(idx.pdc)))     ];
con = [con (imag(temp(idx.pdc)) == imag(S_star(idx.pdc)))     ];

%IC

% idx_vscac_vq_res = reshape(idx.vscac_vq,Grid_para.n_ph,length(idx.vscac_vq)/Grid_para.n_ph);
% con = [con (transpose(real(sum(temp(idx_vscac_vq_res))/3)) == -real(temp(idx.vscdc_vq))) ];

con = [con (real(temp(idx.vscac_vq)) == -real(temp(idx.vscdc_vq))) ];
con = [con (imag(temp(idx.vscac_vq)) == imag(S_star(idx.vscac_vq))) ];


%slack
con = [con (E_abs(1) == abs(E_star(1)) )    ];
con = [con (E_arg(1) == angle(E_star(1)) )    ];

%DC
con = [con (E_abs(idx.vscdc_vq) == real(E_star(idx.vscdc_vq)) )    ];
con = [con (E_arg(idx.vscdc_vq) == 0 )    ];
con = [con (E_arg(idx.pdc) == 0)];


con = [con (E_abs <= 1.05)&(E_abs >= 0.90)]; %B21

options = sdpsettings('verbose', 2);
% options.DualReductions=0;
% options.InfUnbdInfo=1;


sol = optimize(con, 1, options);

time_IP = toc(time_IP)
time_IP_arr = [time_IP_arr, time_IP];

end 

%figure;
%plot(time_IP_arr);
mean(double(E) - E_PF);
max(double(E) - E_PF);

mean(double(E.*conj(YY*E)) - S_PF);
max(double(E.*conj(YY*E)) - S_PF);
%%
% fun = 1;
% 
% lb = 0.95*ones(length(E_star));
% ub = 1.05*ones(length(E_star));
% A = [];
% b = [];
% Aeq = [];
% beq = [];
% 
% x = fmincon(fun,E_0,A,b,Aeq,beq,lb,ub,power(E,S_star,idx,YY))
% 
% function [c] = power(E,S_star,idx,YY)
% temp = E.*conj(YY*E);
% 
% c1 = temp(idx.pqac) - S_star((idx.pqac));
% c2 = temp(idx.pqdc) - S_star((idx.pqdc));
% c3 = temp(idx.vscac_vq) - temp(idx.vscdc_vq);
% 
% c = [c1;c2;c3];
% end



