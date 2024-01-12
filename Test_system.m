%% main script for linear power system state estimation

clear all;
close all;
clc;

addpath(genpath(pwd))

% t = tic;
% mpc = loadcase('FUBM_test_grid');
mpc = loadcase('fubm_case_57_14_2MTDC_ctrls_EPFL_4_2'); %FUBM_m_grid, fubm_case_57_14_2MTDC_ctrls_EPFL, fubm_case_30_2MTDC_ctrls_vt2_pf_EPFL, fubm_case1354pegase_2MTDC_ctrls_pf_qt_dp_EPFL


%% Get power setpoints
S_star = complex(mpc.bus(:,3), mpc.bus(:,4));

%Get Q setpoint from vsc
ind_branch_Q = find(mpc.branch(:,17));
mpc.branch(ind_branch_Q,2);
mpc.branch(ind_branch_Q,17);

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



%% Base values
A_b = mpc.baseMVA*1e6;
V_b= mpc.bus(1,10)*1e3;
Y_b = A_b/V_b^2; 

Vdc_b = mpc.bus(end,10)*1e3;
Adc_b = A_b;
Ydc_b = Adc_b/Vdc_b^2; 




%% Set the nodes types



% [Yac, YYL, YL, YT, YYT, I_b, Ampacities, y_lx, y_tx, A, linedata_ac]  = Ymatrix('linedata_AC_test.txt',A_b,V_b,[]);
% [Ydc, YYLdc, YLdc, YT_dc, YYTdc, I_bdc, Ampacitiesdc, y_ihdc, y_idc, Adc, linedata_dc]  = Ymatrix('linedata_DC_test.txt',Adc_b,Vdc_b,[]);
% 
% Ydc = Ydc(Grid_para.n_ac+1:end,Grid_para.n_ac+1:end)/2;
% Yac = cell2mat(arrayfun(@(x) x*eye(Grid_para.n_ph),Yac,'UniformOutput',false));

%% µG 
% bus_idx_ac = [find(mpc.bus(:,7)== 1)];
% bus_idx_dc = [find(mpc.bus(:,7)== 2); ];
% bus_idx_vsc_dc = [find(mpc.bus(:,7)== 3);  ];
% bus_idx_vsc_ac = [find(mpc.bus(:,7)== 4); ];

%% 30 node
% bus_idx_ac = [find(mpc.bus(:,7)== 1)];
% bus_idx_dc = [find(mpc.bus(:,7)== 2); find(mpc.bus(:,7)== 3) ];
% bus_idx_vsc_dc = [find(mpc.bus(:,7)== 6); find(mpc.bus(:,7)== 7) ];
% bus_idx_vsc_ac = [find(mpc.bus(:,7)== 4); find(mpc.bus(:,7)== 5) ];

%% 57 node
bus_idx_ac = [find(mpc.bus(:,7)== 1); find(mpc.bus(:,7)== 2) ];
bus_idx_dc = [find(mpc.bus(:,7)== 3); find(mpc.bus(:,7)== 4) ];
bus_idx_vsc_dc = [find(mpc.bus(:,7)== 5); find(mpc.bus(:,7)== 6) ];
bus_idx_vsc_ac = [find(mpc.bus(:,7)== 7); find(mpc.bus(:,7)== 8) ];

%% 1354 node
% bus_idx_ac = [find(mpc.bus(:,7)== 0)];
% bus_idx_dc = [find(mpc.bus(:,7)== 3); find(mpc.bus(:,7)== 4) ];
% bus_idx_vsc_dc = [find(mpc.bus(:,7)== 5); find(mpc.bus(:,7)== 6) ];
% bus_idx_vsc_ac = [find(mpc.bus(:,7)== 7); find(mpc.bus(:,7)== 8) ];




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



% id_ac_i = ([(find(mpc.bus(:,7)== 1));(find(mpc.bus(:,7)== 7)); (find(mpc.bus(:,7)== 2));(find(mpc.bus(:,7)== 8))]);
% id_dc_i = ([(find(mpc.bus(:,7)== 3));(find(mpc.bus(:,7)== 5)); (find(mpc.bus(:,7)== 4));(find(mpc.bus(:,7)== 6))]);
% 
% id_ac = ([mpc.bus(find(mpc.bus(:,7)== 1));mpc.bus(find(mpc.bus(:,7)== 7)); mpc.bus(find(mpc.bus(:,7)== 2));mpc.bus(find(mpc.bus(:,7)== 8))]);
% id_dc = ([mpc.bus(find(mpc.bus(:,7)== 3));mpc.bus(find(mpc.bus(:,7)== 5)); mpc.bus(find(mpc.bus(:,7)== 4));mpc.bus(find(mpc.bus(:,7)== 6))]);
% 
% YYac = YYac(id_ac,id_ac);
% YYdc = YYdc(id_dc,id_dc);


bus_index = mpc.bus(:,1);

YY = blkdiag(YYac,YYdc);
Grid_para.G = real(YY);
Grid_para.B = imag(YY);



% indices
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
%% Ground truth

% Load data of nodal voltage and power injections from the .mat file
% load('data_balanced.mat') %data_balanced.mat OR data_unbalanced_light.mat OR data_unbalanced_strong.mat OR data_unbalanced_strong_wlosses.mat
% E_star_true = data.E_star;
% S_star_true = data.S_star;
% 
% E_star_true = E_star_true([1:3:Grid_para.n_ac*3,Grid_para.n_ac*3+1:end]);
% S_star_true = S_star_true([1:3:Grid_para.n_ac*3,Grid_para.n_ac*3+1:end]);
% S_star_true - E_star_true.*conj(YY*E_star_true);
% 
% -S_star*0.1

E_star_true = [0.999999999541630 - 1.04057185227647e-09i;0.988162425456039 - 0.0117867269709593i;0.968998299781958 - 0.0123571332798927i;0.984095999030410 - 0.0168957249977548i;0.978314789055129 - 0.0165779665794751i;0.982583693380744 - 0.0352843639758142i;0.982151593348496 - 0.0405380861746744i;0.974606425845804 - 0.0272541337196129i;0.974378898935833 - 0.0279640779474437i;0.946794288082694 - 0.0483317187156929i;0.935591981657944 - 0.0567359819951060i;0.942091637992704 - 0.0514353228869956i;0.938060790708089 - 0.0540954742917767i;0.982378146978744 - 0.0428993479789816i;0.974683408536357 - 0.0284573876656165i;0.938438952993418 - 0.0545902471027461i;0.935701734682472 - 0.0573302735330433i;0.982194867397862 - 0.0411311674317191i;0.999993428599251 + 0.00000000000000i;0.999993432300937 + 0.00000000000000i;0.999993434701429 + 0.00000000000000i;0.999993422962788 + 0.00000000000000i;1.00105968576700 + 0.00000000000000i;1.00191438906190 + 0.00000000000000i;0.998919257446391 + 0.00000000000000i;0.998061576904007 + 0.00000000000000i];
S_star_true = [1.19858842534728 - 0.464427350534137i;9.40942338545655e-08 - 2.25137583800892e-06i;-0.300003441130744 - 9.29922099683583e-05i;2.79505538174175e-07 - 5.57122986090230e-07i;-0.374982912657203 - 0.0751166925054122i;2.48776573643737e-08 + 3.05986437111428e-09i;-3.24599359925611e-09 - 2.28391606785539e-08i;-3.00279939926335e-07 + 3.17079741887731e-07i;-0.100539148166344 + 0.000175112822269798i;-2.12288144408844e-09 + 1.86552180316074e-09i;-0.150011014113213 + 0.0299495950049383i;8.74133425162435e-09 + 1.71419284697175e-08i;-0.202033363862924 + 0.000442433739335894i;-0.00254075300730539 + 0.124281463779260i;0.0170343712117197 + 0.100053535517447i;0.0307715450787493 + 0.100497607422195i;-0.0171862673417051 + 0.100256304939446i;-0.0309069461434549 + 0.100303513410919i;-0.0170600025751033 - 1.68527898446222e-09i;-0.0307351063151647 - 1.75365451500785e-09i;0.0171867232438888 - 1.51404537489410e-09i;0.0309093336473239 - 1.44543228531351e-09i;0.0499996805218759 - 4.65894443326827e-09i;0.0499996864745647 - 3.01988200214097e-09i;-0.0500003237304102 - 4.94113590997407e-09i;-0.0500003236411775 - 3.38015562695775e-09i];

-S_star_true*0.1
% Get data directS_star_truely from EMPT simulation
% 
% repeat=1;
% ZIN_polyphase = [];
% [Nodal_V_mag,Nodal_V_angle, Nodal_I_mag, Nodal_I_angle, Flow_I_mag, Flow_I_angle, Idc_flow, Idc_inj, Vdc_LF, n_timesteps, Mreal, Mimag, Mabs,Nodal_P,Nodal_Q,Pdc_inj,M_LF] = GetEMTPdata_LF(A_b,V_b,Adc_b, Vdc_b,repeat,ZIN_polyphase,Grid_para.n_ph);
% 
% [Nodal_V_mag,Nodal_V_angle, Nodal_I_mag, Nodal_I_angle, Flow_I_mag, Flow_I_angle, Idc_flow, Idc_inj, Vdc_LF, n_timesteps, M_real, M_imag, Mabs,Nodal_P,Nodal_Q,Pdc_inj,M_LF] = GetEMTPdata_LF_simplified3(A_b,V_b,Adc_b, Vdc_b,1,[],Grid_para.n_ph); 
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
% E_star_true = [V_complex_LF(1:end,1); Vdc_LF(:,1)];
% S_star_true = [transpose(complex(Nodal_P, Nodal_Q))/3; transpose(complex(Pdc_inj,0))];



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
tol = 1e-8;
n_max = 20;

time_array = [];
for t = 1:1
    
% [E,J,n_iter] = NR_rectangularACDC_3ph(Grid_para,Filter_para,S_star,E_star,E_0,idx,tol,n_max);
[E,J,n_iter,time,K] = NR_rectangularACDC_1ph_general_slack(Grid_para,Filter_para,S_star,E_star,E_0,idx,tol,n_max);

time_array = [time_array ; time];


end

disp({'time: ', mean(time_array), ' - nb iterations: ', n_iter, ' - nb states: ', size(J,1)})

% toc(t)
S = E.*conj(YY*E);

E_res_me = [abs(E) , angle(E)*180/pi];
E_res_true = [abs(E_star_true) , angle(E_star_true)*180/pi];


% figure
% scatter(1:length(E),abs(E))

%orderdc = [94,98,97,96,95,92,93,91];
% orderac = [78,75,74,79,80,76,77,73];

orderdc = [94,95,96,97,98,91,92,93];
orderac = [80,78,79,74,75,73,76,77];
[real(S(orderac)), imag(S(orderac)), abs(E(orderac)), real(E(orderdc)) ];



%% validate quadratic
real(S(orderac))
S(orderdc)

A_ic = 11.033e-3;
B_ic = 3.464e-3;
C_ic = 11e-3;
I_mag = abs(YY * E);
P_quad_ic = A_ic + B_ic*I_mag + C_ic*I_mag.^2;

P_quad_ic(orderac)


results = [real(S(orderac)), imag(S(orderac)), abs(E(orderac)), real(E(orderdc)), real(P_quad_ic(orderac)) ];


%% show AC results
Sac = S(1:Grid_para.n_ac)
Sdc = S(1+Grid_para.n_ac:end)


sum(Sac)
sum(Sac.*(Sac>0))
sum(Sac.*(Sac<0))

sum(Sdc)
sum(Sdc.*(Sdc>0))
sum(Sdc.*(Sdc<0))


Grid_para.n_dc = length(bus_dc) + length(bus_vsc_dc);
Grid_para.n_ac = length(bus_ac) + length(bus_vsc_ac);



bus_idx_ac57 = [find(mpc.bus(:,7)== 1); find(mpc.bus(:,7)== 7) ];
% bus_idx_ac57 = bus_idx_ac57(2:end);
bus_idx_ac14 = [find(mpc.bus(:,7)== 2); find(mpc.bus(:,7)== 8) ];
% bus_idx_ac14 = bus_idx_ac14(2:end);

bus_idx_dc2 = [find(mpc.bus(:,7)== 3); find(mpc.bus(:,7)== 5) ];
bus_idx_dc1 = [find(mpc.bus(:,7)== 4); find(mpc.bus(:,7)== 6) ];

Sac_57 = S(bus_idx_ac57);
Sac_14 = S(bus_idx_ac14);
Sdc_1 = real(S(bus_idx_dc1));
Sdc_2 = real(S(bus_idx_dc2));




%P ac
[sum(real(Sac_57(1:end)).*(real(Sac_57(1:end))>0)) , sum(real(Sac_57(1:end)).*(real(Sac_57(1:end))<0)), sum(real(Sac_57));
 sum(real(Sac_14(1:end)).*(real(Sac_14(1:end))>0)) , sum(real(Sac_14(1:end)).*(real(Sac_14(1:end))<0)), sum(real(Sac_14))]

%Q ac
[sum(imag(Sac_57(1:end)).*(imag(Sac_57(1:end))>0)) , sum(imag(Sac_57(1:end)).*(imag(Sac_57(1:end))<0)), sum(imag(Sac_57));
 sum(imag(Sac_14(1:end)).*(imag(Sac_14(1:end))>0)) , sum(imag(Sac_14(1:end)).*(imag(Sac_14(1:end))<0)), sum(imag(Sac_14))]


%P dc
[sum(Sdc_1.*(Sdc_1>0)), sum(Sdc_1.*(Sdc_1<0)),sum(Sdc_1);
 sum(Sdc_2.*(Sdc_2>0)), sum(Sdc_2.*(Sdc_2<0)),sum(Sdc_2)]

% table paper
[sum(real(Sac_57(1:end)).*(real(Sac_57(1:end))>0)) , sum(imag(Sac_57(1:end)).*(imag(Sac_57(1:end))>0)) ,sum(real(Sac_57(1:end)).*(real(Sac_57(1:end))<0)),  sum(imag(Sac_57(1:end)).*(imag(Sac_57(1:end))<0)),sum(real(Sac_57)), sum(imag(Sac_57));
 sum(real(Sac_14(1:end)).*(real(Sac_14(1:end))>0)) , sum(imag(Sac_14(1:end)).*(imag(Sac_14(1:end))>0)) ,sum(real(Sac_14(1:end)).*(real(Sac_14(1:end))<0)),  sum(imag(Sac_14(1:end)).*(imag(Sac_14(1:end))<0)),sum(real(Sac_14)), sum(imag(Sac_14))]



%% losses computation
I = YY*E;
YYL = blkdiag(YYLac,YYLdc);
YYT = blkdiag(YYTac,YYTdc);

E_57 = E(bus_idx_ac57);
E_14 = E(bus_idx_ac14);
E_1 = real(E(bus_idx_dc1));
E_2 = real(E(bus_idx_dc2));

YYL_57 = YYL(bus_idx_ac57,bus_idx_ac57);
YYL_14 = YYL(bus_idx_ac14,bus_idx_ac14);
YYL_1 = real(YYL(bus_idx_dc1,bus_idx_dc1));
YYL_2 = real(YYL(bus_idx_dc2,bus_idx_dc2));

YYT_57 = YYT(bus_idx_ac57,bus_idx_ac57);
YYT_14 = YYT(bus_idx_ac14,bus_idx_ac14);
YYT_1 = real(YYT(bus_idx_dc1,bus_idx_dc1));
YYT_2 = real(YYT(bus_idx_dc2,bus_idx_dc2));

Iflow_57 = zeros(size(YYL_57));
SLflow_57 = zeros(size(YYL_57));
STflow_57 = zeros(size(YYL_57));
for i = 1:size(YYL_57,1)
    for j = 1:size(YYL_57,2)
        Iflow_57(i,j) = YYL_57(i,j)*(E_57(j) - E_57(i)) + 0*YYT_57(i,j)*E_57(i);
        SLflow_57(i,j) = (E_57(j) - E_57(i)) .* conj(YYL_57(i,j)*(E_57(j) - E_57(i)));
        STflow_57(i,j) = E_57(i) .* conj(YYT_57(i,j) * E_57(i));
    end
end

Iflow_14 = zeros(size(YYL_14));
SLflow_14 = zeros(size(YYL_14));
STflow_14 = zeros(size(YYL_14));
for i = 1:size(YYL_14,1)
    for j = 1:size(YYL_14,2)
        Iflow_14(i,j) = YYL_14(i,j)*(E_14(j) - E_14(i)) + 0*YYT_14(i,j)*E_14(i);
        SLflow_14(i,j) = (E_14(j) - E_14(i)) .* conj(YYL_14(i,j)*(E_14(j) - E_14(i)));
        STflow_14(i,j) = E_14(i) .* conj(YYT_14(i,j) * E_14(i));
    end
end

Iflow_1 = zeros(size(YYL_1));
SLflow_1 = zeros(size(YYL_1));
STflow_1 = zeros(size(YYL_1));
for i = 1:size(YYL_1,1)
    for j = 1:size(YYL_1,2)
        Iflow_1(i,j) = YYL_1(i,j)*(E_1(j) - E_1(i)) + 0*YYT_1(i,j)*E_1(i);
        SLflow_1(i,j) = (E_1(j) - E_1(i)) .* conj(YYL_1(i,j)*(E_1(j) - E_1(i)));
        STflow_1(i,j) = E_1(i) .* conj(YYT_1(i,j) * E_1(i));
    end
end

Iflow_2 = zeros(size(YYL_2));
SLflow_2 = zeros(size(YYL_2));
STflow_2 = zeros(size(YYL_2));
for i = 1:size(YYL_2,1)
    for j = 1:size(YYL_2,2)
        Iflow_2(i,j) = YYL_2(i,j)*(E_2(j) - E_2(i)) + 0*YYT_2(i,j)*E_2(i);
        SLflow_2(i,j) = (E_2(j) - E_2(i)) .* conj(YYL_2(i,j)*(E_2(j) - E_2(i)));
        STflow_2(i,j) = E_2(i) .* conj(YYT_2(i,j) * E_2(i));
    end
end

pos_57 = find(triu(YYL_57) ~= 0);
pos_14 = find(triu(YYL_14) ~= 0);
pos_1 = find(triu(YYL_1) ~= 0);
pos_2 = find(triu(YYL_2) ~= 0);

sum(SLflow_57(pos_57))
sum(STflow_57(pos_57))

sum(SLflow_14(pos_14))
sum(STflow_14(pos_14))

sum(SLflow_1(pos_1))
sum(STflow_1(pos_1))

sum(SLflow_2(pos_2))
sum(STflow_2(pos_2))


%% Visualize the solution
% E_delta = E - E_star;
% S_delta = E.*conj(YY*E) - S_star;
% % 
% figure
% subplot(2,1,1)
% hold on
% scatter(1:length(E_res_true), E_res_true(:,1) - E_res_me(:,1),'filled')
% scatter(1:length(E_res_true), E_res_true(:,1) - E_res_fubm(:,1),'filled')
% ylabel('Delta E - mag')
% legend({'EPFL','FUBM'})
% 
% subplot(2,1,2)
% hold on
% scatter(1:length(E_res_true), E_res_true(:,2) - E_res_me(:,2),'filled')
% scatter(1:length(E_res_true), E_res_true(:,2) - E_res_fubm(:,2),'filled')
% ylabel('Delta E - angle (degree)')
% legend({'EPFL','FUBM'})
% 
