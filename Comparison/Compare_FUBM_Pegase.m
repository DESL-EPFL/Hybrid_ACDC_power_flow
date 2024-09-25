% script for linear power system state estimation
% Comparison with the FUBM-based matpower model 

% Willem lambrichts
% General and Unified Model of the Power Flow Problem in Multiterminal AC/DC Networks
% IEEE Transactions on Power Systems 
% 10.1109/TPWRS.2024.3378926

clear all
close all
clc

addpath(genpath(pwd))

%% FUBM for hybrid microgrid
[MVAbase, bus,gen, branch, success, et, time_fubm, Ybus, iter_fubm] = runpf('IEEE_Pegase_HVDC'); 

E_fubm = bus(:,8) .* exp(1i*bus(:,9)*pi/180);
S_fubm = E_fubm.*conj(Ybus*E_fubm);

%% Proposed method (EPFL)

mpc = loadcase('IEEE_Pegase_HVDC'); 

%set bus types
bus_idx_ac = [find(mpc.bus(:,7)== 0)];
bus_idx_dc = [find(mpc.bus(:,7)== 3); find(mpc.bus(:,7)== 4) ];
bus_idx_vsc_dc = [find(mpc.bus(:,7)== 5); find(mpc.bus(:,7)== 6) ];
bus_idx_vsc_ac = [find(mpc.bus(:,7)== 7); find(mpc.bus(:,7)== 8) ];



%% Get base values
A_b = mpc.baseMVA*1e6;
V_b= mpc.bus(1,10)*1e3;
Y_b = A_b/V_b^2; 

Vdc_b = mpc.bus(end,10)*1e3;
Adc_b = A_b;
Ydc_b = Adc_b/Vdc_b^2; 


%% Get power setpoints
S_star = complex(mpc.bus(:,3), mpc.bus(:,4));

%Get Q setpoint from vsc
ind_branch_Q = find(mpc.branch(:,17));

ind_bus_from_Q = get_index(mpc.bus(:,1), mpc.branch(ind_branch_Q,1));
ind_bus_to_Q   = get_index(mpc.bus(:,1), mpc.branch(ind_branch_Q,2));

S_star(ind_bus_to_Q) = S_star(ind_bus_to_Q) + complex(-mpc.branch(ind_branch_Q,14),mpc.branch(ind_branch_Q,17));
S_star(ind_bus_from_Q) = S_star(ind_bus_from_Q) + mpc.branch(ind_branch_Q,14);
S_star = -S_star/mpc.baseMVA; % matpower uses a different convention


%% Get voltage setpoints 
E_star = mpc.bus(:,8) .* exp(1i* mpc.bus(:,9));

%Get Vdc setpoint from vsc
ind_branch_Vdc = find(mpc.branch(:,22));
mpc.branch(ind_branch_Vdc,1);
mpc.branch(ind_branch_Q,22);

ind_bus_Vdc = [];
for i = 1:length(ind_branch_Vdc)
ind_bus_Vdc = [ind_bus_Vdc ; find(mpc.bus(:,1)== mpc.branch(ind_branch_Vdc(i),1))];
end

E_star(ind_bus_Vdc) = mpc.branch(ind_branch_Vdc,22);


bus_ac = mpc.bus(bus_idx_ac,1);
bus_dc = mpc.bus(bus_idx_dc,1);
bus_vsc_dc = mpc.bus(bus_idx_vsc_dc,1);
bus_vsc_ac = mpc.bus(bus_idx_vsc_ac,1);

branch_idx_ac = [];
branch_idx_dc = [];
branch_idx_vsc_dc = [];
branch_idx_vsc_ac = [];

for i = 1:length(mpc.branch(:,1))
    
    if any(bus_ac(:) == mpc.branch(i,1)) && any(bus_ac(:) == mpc.branch(i,2)) && mpc.branch(i,11)
        branch_idx_ac = [branch_idx_ac; i];
        
    elseif any(bus_dc(:) == mpc.branch(i,1)) && any(bus_dc(:) == mpc.branch(i,2))  && mpc.branch(i,11)
        branch_idx_dc = [branch_idx_dc; i];

    elseif any(bus_dc(:) == mpc.branch(i,1)) && any(bus_vsc_dc(:) == mpc.branch(i,2))  && mpc.branch(i,11)
        branch_idx_vsc_dc = [branch_idx_vsc_dc; i];
        
    elseif any(bus_ac(:) == mpc.branch(i,1)) && any(bus_vsc_ac(:) == mpc.branch(i,2))  && mpc.branch(i,11)
        branch_idx_vsc_ac = [branch_idx_vsc_ac; i];
    else
        mpc.branch(i,1:2);
    end
end

branch_vsc_dc = mpc.branch(branch_idx_vsc_dc,[1,2]);
branch_vsc_dc = branch_vsc_dc(:);

branch_vsc_ac = mpc.branch(branch_idx_vsc_ac,[1,2]);
branch_vsc_ac = branch_vsc_ac(:);



branch_vscVQ = mpc.branch(find(mpc.branch(:,22)),[1,2]);
branch_vscPQ = mpc.branch(find(mpc.branch(:,14)),[1,2]);

%% Set the Grid parameters
Grid_para.n_dc = length(bus_dc) + length(bus_vsc_dc);
Grid_para.n_ac = length(bus_ac) + length(bus_vsc_ac);
Grid_para.n_ph = 1;
Grid_para.n_nodes = Grid_para.n_ac*Grid_para.n_ph + Grid_para.n_dc;
Grid_para.V_b = V_b;
Grid_para.Y_b = Y_b;


[YYac, linedata_ac, YYLac, YYTac]  =  Ymatrix_matpower(mpc.branch([branch_idx_ac;branch_idx_vsc_ac ],1:5));
[YYdc, linedata_dc, YYLdc, YYTdc]  =  Ymatrix_matpower(mpc.branch([branch_idx_dc;branch_idx_vsc_dc ],1:5));
YYac = YYac(sort([bus_ac;bus_vsc_ac]),sort([bus_ac;bus_vsc_ac]));
YYdc = YYdc(sort([bus_dc;bus_vsc_dc]),sort([bus_dc;bus_vsc_dc]));

YYLac = YYLac(sort([bus_ac;bus_vsc_ac]),sort([bus_ac;bus_vsc_ac]));
YYLdc = YYLdc(sort([bus_dc;bus_vsc_dc]),sort([bus_dc;bus_vsc_dc]));

YYTac = YYTac(sort([bus_ac;bus_vsc_ac]),sort([bus_ac;bus_vsc_ac]));
YYTdc = YYTdc(sort([bus_dc;bus_vsc_dc]),sort([bus_dc;bus_vsc_dc]));

Grid_para.YY = blkdiag(YYac,YYdc);
Grid_para.G = real(Grid_para.YY);
Grid_para.B = imag(Grid_para.YY);



% indices
bus_index = mpc.bus(:,1);
idx1.slack = intersect(mpc.areas(:,2), bus_ac); %
idx1.pqac = [setdiff(bus_ac,[bus_vsc_ac;idx1.slack;mpc.gen(:,1)])];
idx1.pvac = [setdiff(mpc.gen(:,1),idx1.slack)];

idx1.pdc = setdiff(bus_dc,bus_vsc_dc);
idx1.vdc = []';

idx1.vscac_pq = intersect(branch_vscPQ,bus_vsc_ac);
idx1.vscac_vq = intersect(branch_vscVQ,bus_vsc_ac);

idx1.vscdc_pq = intersect(branch_vscPQ,bus_vsc_dc);
idx1.vscdc_vq = intersect(branch_vscVQ,bus_vsc_dc);

linedata = [linedata_ac;linedata_dc];
Grid_para = Get_Converter_para_FUBM(idx1,linedata,Grid_para,mpc);


idx1.slack = get_index(bus_index, idx1.slack);
idx1.pvac = get_index(bus_index, idx1.pvac);
idx1.pqac = get_index(bus_index, idx1.pqac);
idx1.vscac_pq = get_index(bus_index, idx1.vscac_pq);
idx1.vscac_vq = get_index(bus_index, idx1.vscac_vq);
idx1.vscdc_pq = get_index(bus_index, idx1.vscdc_pq);
idx1.vscdc_vq = get_index(bus_index, idx1.vscdc_vq);
idx1.pdc = get_index(bus_index, idx1.pdc);


idx = idx1;
Grid_para.pos_ac3(:,1) = get_index(bus_index, Grid_para.pos_ac3(:,1));
Grid_para.pos_ac3(:,2) = get_index(bus_index, Grid_para.pos_ac3(:,2));


Grid_para.pos_dc3(:,1) = get_index(bus_index, Grid_para.pos_dc3(:,1));
Grid_para.pos_dc3(:,2) = get_index(bus_index, Grid_para.pos_dc3(:,2));



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
    E_0 = [repmat([1; exp(-1i*pi*2/3);  exp(1i*pi*2/3) ], Grid_para.n_ac,1 ) ; ones(Grid_para.n_dc,1)];
elseif Grid_para.n_ph == 1
    E_0 = [ones(Grid_para.n_ac + Grid_para.n_dc,1)];
end

%% Solve the PF
simulation_para.tol = 1e-6;
simulation_para.n_max = 100;
simulation_para.E_0 = E_0;

time_array = [];
for i = 1:10
    
    % Solve the LF without a quadratic loss term
    % [E,J,n_iter,time] = NR_rectangularACDC_1ph_general(Grid_para,Filter_para,S_star,E_star,idx,simulation_para);

    % Solve the LF with a quadratic loss term
    [E_epfl,J,n_iter,time] = NR_rectangularACDC_1ph_general_quadratic_loss(Grid_para,Filter_para,S_star,E_star,idx,simulation_para);

    time_array = [time_array, time];
end
S_epfl = E_epfl.*conj(Grid_para.YY*E_epfl);
time_epfl = mean(time_array);
iter_epfl = n_iter;

disp(['FUBM-based method - time: ', num2str(time_fubm), ' - nb iterations: ', num2str(iter_fubm)])
disp(['Proposed method - time: ', num2str(time_epfl), ' - nb iterations: ', num2str(iter_epfl)])
