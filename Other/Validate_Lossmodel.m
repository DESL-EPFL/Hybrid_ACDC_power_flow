% main script for linear power system state estimation

clear all;
close all;
% clc;

set(groot,'defaultFigureVisible','on')
% set(groot,'defaultFigureVisible','off')

%% Simulation options
n = 2; %0 for for M_abs 
       %1 for M real + j M_imag
       %2 for symmetric components (only 3ph)
n_phases = 3; % 1 (single-phase equivalent)
adaptive = 0; % 0 for DFK
              % 1 for PECE
              % 2 for PECE Diag 
repeat = 1; %adds the measurements profiles n times after eachother
DC_state = 1; %0 = V, 1 = I %Choose the state variable
threshold = 4.5;
add_BD = 0; %0=no, 1=PMU level mag, 2=measurement level, 3=PMU level rad
BD_detection = 0;
N = 50;


%% Import paths
addpath('/Users/willem/Documents/phd/State_estimation/mosek/9.2/toolbox/r2015a');
addpath('/Users/willem/Documents/phd/State_estimation/maxdet');
addpath(genpath('/Users/willem/Documents/phd/State_estimation/YALMIP-master'));
addpath('/Users/willem/Documents/phd/State_estimation/sdpt3-master/')
addpath('/Users/willem/Documents/tools/aboxplot/')
run('install_sdpt3.m')
addpath '/Users/willem/Documents/phd/Micro grid specs/'
addpath('/Users/willem/Documents/phd/State_estimation')
VSI_name = '/Users/willem/Documents/phd/Micro grid specs/DCdata.txt';
DC_grid_name = '/Users/willem/Documents/phd/Micro grid specs/DCgrid_data2.txt';
linedata_DCgrid = '/Users/willem/Documents/phd/Micro grid specs/linedata_DCgrid.txt';
folder= '/Users/willem/Documents/phd/State_estimation/matlab_Full_SE_3phase/bad data rejection';
addpath '/Users/willem/Documents/phd/Optimal_power_flow/Sensitivity_codes'

%%  Step 1: Load nodal admittance matrix and base values
%   The data need to be saved in the folder 'Data_SE'
%   in the same directory as this script.



modess = {'P5';'E4';'Q3';'P2'}
% for zzz= 1:4
    zzz = 1;
% ! Indicate the number of phases !
n_dclinks = 4;

% AC side
A_b = 1e5;
V_b= 400;
Y_b = A_b/V_b^2; % base admittance

%DC side
Vdc_b = 800;
Adc_b = A_b;
Ydc_b = Adc_b/Vdc_b^2; % base admittance

%Linkage AC-DC
Kv = Vdc_b/V_b;
Ka = Adc_b/A_b;
Ky = Ydc_b/Y_b;
% 
% Kroner reduction
ZIN = []; %[6,12]; %6,12
ZIN_polyphase = polyphase_indices(ZIN,n_phases);
% % Get admittance matrix
% [YY, YYL, YL, YT, YYT, I_b, Ampacities, y_ih, y_i, A, linedata]  = Ymatrix('linedata_mG_dc.txt',A_b,V_b,ZIN);
% YY = cell2mat(arrayfun(@(x) x*eye(n_phases),YY,'UniformOutput',false));
% YYL = cell2mat(arrayfun(@(x) x*eye(n_phases),YYL,'UniformOutput',false));
% YYT = cell2mat(arrayfun(@(x) x*eye(n_phases),YYT,'UniformOutput',false));
% 
% 
% %%  Step 2: Construct the reduced nodal admittance matrix
% %   Since the voltage of the slack node is fixed,
% %   it does not need be estimated.
% 
% Y_SE = YY;
% n_nodes = size(Y_SE,1)/n_phases;
% n_nodesac = n_nodes;
% n_nodesdc = 8;
% %%  Step 4: Configure the state estimator
% 
% % ! Select which nodes to estimate (node indices) !
% sel_SE = transpose(1:n_nodes); % estimate entire network
% % Create polyphase indices
% idx_SE = polyphase_indices(sel_SE,n_phases);
% 
% n_states = 2 * length(idx_SE) + n_dclinks;
% 
% %%  Step 5: Construct the measurement model
% 
% % AC side
% % ! Select voltage measurement locations (node incices) !
% sel_PMU_V = [1,3,5,9,11,13,14];% transpose(1:n_nodes-4);%
% sel_PMU_V = kron_reduce_PMU(sel_PMU_V,ZIN);
% 
% % ! Select current injection measurement locations (node indices) !
% sel_PMU_I = [1,3,5,9,11,13,14];%transpose(1:n_nodes-4);%
% sel_PMU_I = kron_reduce_PMU(sel_PMU_I,ZIN);
% 
% 
% % ! Select current flow measurement locations (line indices)
% loc_PMU_Iflow = [ 9,15; %DC1
%                  13,16; %DC2
%                  11,17; %DC3
%                   7,18]; %DC4
% loc_PMU_Iflow = kron_reduce_PMU(loc_PMU_Iflow,ZIN);

% % ! Indicate standard deviation of voltage measurement noise !
% sd_PMU_V = 1e-3;
% 
% % ! Indicate standard deviation of current measurement noise !
% sd_PMU_I = 1e-3;

%% DC side
% ! Select voltage measurement locations (node incices) !
% sel_PMU_Vdc = transpose(1:n_dclinks*2); %
% 
% % ! Select current injection measurement locations (node indices) !
% sel_PMU_Idc = transpose(1:n_dclinks); %

% % ! Indicate standard deviation of dc current measurement noise !
% sd_PMU_Idc = 1e-3;
 
% % ! Indicate standard deviation of dc voltage measurement noise !
% sd_PMU_Vdc = 1e-3;

% Create polyphase indices
% idx_PMU_V = polyphase_indices(sel_PMU_V,n_phases);
% idx_PMU_I = polyphase_indices(sel_PMU_I,n_phases);
% idx_PMU_Iflow = polyphase_indices(loc_PMU_Iflow,n_phases);
% idx_PMU_Vdc = sel_PMU_Vdc;%sel_PMU_Vdc;
% idx_PMU_Idc = sel_PMU_Idc;
% 
% n_measurements = 2*(numel(idx_PMU_V)+numel(idx_PMU_I)+length(idx_PMU_Iflow))+length(idx_PMU_Vdc)+length(idx_PMU_Idc);
% 
% in_memory = 1;
% [Nodal_V_mag,Nodal_V_angle, Nodal_I_mag, Nodal_I_angle, Flow_I_mag, Flow_I_angle, Idc_flow, Idc_inj, Vdc_LF, n_timesteps, Mreal, Mimag, Mabs] = GetEMTPdataV6_IpDKF(A_b,V_b,Adc_b, Vdc_b,repeat,ZIN_polyphase,n_phases,in_memory) ;
% Mreal = Mreal/2;
% Mimag = Mimag/2;
% % % 
% [Nodal_V_mag2,Nodal_V_angle2, Nodal_I_mag2, Nodal_I_angle2, Flow_I_mag2, Flow_I_angle2, Idc_flow2, Idc_inj2, Vdc_LF2, n_timesteps2, Mreal2, Mimag2, Mabs2] = GetEMTPdataV6(A_b,V_b,Adc_b, Vdc_b,repeat,ZIN_polyphase,n_phases);

% [Nodal_V_mag,Nodal_V_angle, Nodal_I_mag, Nodal_I_angle, Flow_I_mag, Flow_I_angle, Idc_flow, Idc_inj, Vdc_LF, n_timesteps, Mreal, Mimag, Mabs] = GetEMTPdataV6(A_b,V_b,Adc_b, Vdc_b,repeat,ZIN_polyphase,n_phases);
[Nodal_V_mag,Nodal_V_angle, Nodal_I_mag, Nodal_I_angle, Flow_I_mag, Flow_I_angle, Idc_flow, Idc_inj, Vdc_LF, n_timesteps, Mreal, Mimag, Mabs,Nodal_P,Nodal_Q,Pdc_inj,M_LF] = GetEMTPdata_LF_simplified3(A_b,V_b,Adc_b, Vdc_b,repeat,ZIN_polyphase,n_phases);

%[Nodal_V_mag2,Nodal_V_angle2, Nodal_I_mag2, Nodal_I_angle2, Flow_I_mag2, Flow_I_angle2, Idc_flow2, Idc_inj2, Vdc_LF2,n_timesteps2, Mreal2, Mimag2, Mabs2] = GetEMTPdataV6(A_b,V_b,Adc_b,Vdc_b,repeat,ZIN_polyphase,n_phases);
% [Nodal_V_mag,Nodal_V_angle, Nodal_I_mag, Nodal_I_angle, Flow_I_mag, Flow_I_angle, Idc_flow, Idc_inj, Vdc_LF,n_timesteps, Mreal, Mimag, Mabs] = GetEMTPdataV5(A_b,V_b,Adc_b,Vdc_b,repeat,ZIN_polyphase,n_phases);


% Compute current injections using the admittance matrix
Nodal_I_mag2 = Nodal_I_mag;
Nodal_I_angle2 = Nodal_I_angle;

Flow_I_mag2 = Nodal_I_mag2; 
Flow_I_angle2 = Nodal_I_angle2; 


% Make LF data
V_complex_LF = transpose(complex(Nodal_V_mag.*cos(Nodal_V_angle), Nodal_V_mag.*sin(Nodal_V_angle)));
I_complex_LF = transpose(complex(Nodal_I_mag2.*cos(Nodal_I_angle2), Nodal_I_mag2.*sin(Nodal_I_angle2)));
If_complex_LF = transpose(complex(Flow_I_mag2.*cos(Flow_I_angle2), Flow_I_mag2.*sin(Flow_I_angle2))); %I_LF(15:18,:); %
Vdc_LF = transpose(Vdc_LF);
Idc_flow = transpose(Idc_flow);
Idc_inj = transpose(Idc_inj);

Pac = real(V_complex_LF(7:9,:).*conj(I_complex_LF(7:9,:)))/3;
Pdc = -Vdc_LF(1,:).*Idc_inj(1,:)/2;

figure
hold on
plot(sum(Pac))
plot(Pdc)

figure
plot(sum(Pac) - Pdc)

Y_b = 0.625;
R = 0.008*Y_b; %checked
X = 0.04*Y_b;  %checked
Zf = complex(R,X);

        a = 5.96e-13;
        b = 36.38;
        Vj0 = 0.82/V_b;
        Rj0 = 1/900*Y_b;
        Imag = abs(I_complex_LF(7:9,:));
        dio = @(I) 1/b*log(I/a);
        R_eq_ctu = arrayfun(dio,Imag)/V_b./Imag*4/pi;
        R_eq_lin = (Vj0*4/pi./Imag + Rj0);
        IGBT_piecewise = [
        0	0
        0.04926559627563	0.7
        2.30625399327864	0.8
        15.7793399043317	0.85
        107.547461516782	0.8999
        735.837403888342	0.9499
        1588.01477341768	0.9699];

        I_b =Y_b*V_b;
        R_eq_ctu = interp1(IGBT_piecewise(:,1),IGBT_piecewise(:,2),Imag*I_b)./V_b./Imag;
        R_eq = R_eq_ctu*4/pi; 
        R_eq(isnan(R_eq))=0;
        R_eq(isinf(R_eq))=0;

        
Pfilter = -real(Zf*I_complex_LF(7:9,:).*conj(I_complex_LF(7:9,:)))/3;
Pvsc = -real(R_eq.*I_complex_LF(7:9,:).*conj(I_complex_LF(7:9,:)))/3;
Plossdc = 0.001*Idc_inj(1,:).*Idc_inj(1,:);

Pother1 = -0.0065*abs(Pac);
u = 11.033
v = 3.464
w = 6.667e-3
Pother2 = u + v*abs(I_complex_LF(7:9,:)) + w*abs(I_complex_LF(7:9,:)).^2;


figure
hold on
plot(sum(Pac) - Pdc)
plot(sum(Pvsc)+sum(Pfilter))
legend('Pac - Pdc','Ploss theory')

plot(sum(Pother1))
plot(-sum(Pother2))

figure
plot(sum(Pvsc)+sum(Pfilter) - (sum(Pac) - Pdc))

% 
mode = char(modess(zzz));%'P5' %'P5'E4' 'Q3'
z1 = 85;
switch mode
    case 'P5'
        z2 = 200;
    case 'E4'
        z2 = 382;
    case 'Q3'
        z2 = 549;
    case 'P2'
        z2 = 700;
    otherwise
        disp('NOPE')
end

%% Two port model
% AVG model
    Rf = 1e-2/2*ones(3,1);
    Xf = 5e-2/2*ones(3,1);
%Switching model
%     Rf = 1e-2/2;
%     Xf = 0.3/2;
% AVG model big R
%     Rf = 1e-1/2*ones(6,1);
%     Xf = 5e-2/2*ones(6,1);

Zf = Rf + 1i*Xf;

% M = Mreal + 1i*Mimag;
% M = M(z,1+1);
% 
% t  = 0:0.0001:0.04;
% M_t = abs(M).*sin(2*pi*50*t+ angle(M));
% figure
% plot(t,M_t)
% 
% Vacg = V_complex_LF(37+1,z)
% Vacg_t = abs(Vacg).*sin(2*pi*50*t+ angle(Vacg));
% 
% 
% Vac = V_complex_LF(37+1,z) + If_complex_LF(1+1,z).*complex(Rf,Xf); 
% Vac_t = abs(Vac).*sin(2*pi*50*t+ angle(Vac));
% 
% Iac = If_complex_LF(1+1,z);
% Iac_t = abs(Iac).*sin(2*pi*50*t+ angle(Iac));
% 
% Pac = Nodal_P(z,43+1);
% Qac = Nodal_Q(z,43+1);
% 
% Vdc = Vdc_LF(1+0,z);
% Vdc_t = Vdc*ones(1,length(t));
% 
% Idc = Idc_flow(1+0,z)/2;
% Idc_t = Idc*ones(1,length(t));
% 
% Pdc = Pdc_inj(z,1+0)
% 
% coup = []

M = conj(M_LF(z1,:))';
Vacg = V_complex_LF(:,z1);
% Vac = V_complex_LF(7:9,z1) + If_complex_LF(:,z1).*complex(Rf,Xf);%.*[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]'; 
Iac = If_complex_LF(:,z1);%.*[1 1 1 1]';
Pac = Nodal_P;
Qac = Nodal_Q;
Vdc = reshape(repmat(Vdc_LF(:,z1)',[1,3]),6,1);
Idc = -reshape(repmat(Idc_flow(:,z1)/2,[1,3])',6,1); 

%__________________
       a = 5.96e-13;
        b = 36.38;
        Vj0 = 0.82/V_b;
        Rj0 = 1/900*Y_b;
        Imag = abs(Iac);
        dio = @(I) 1/b*log(I/a);
        R_eq_ctu = 2*arrayfun(dio,Imag)/V_b./Imag*4/pi;
        R_eq_lin = 2*(Vj0*4/pi./Imag + Rj0);
        IGBT_piecewise = [
        0	0
        0.04926559627563	0.7
        2.30625399327864	0.8
        15.7793399043317	0.85
        107.547461516782	0.8999
        735.837403888342	0.9499
        1588.01477341768	0.9699];

        I_b =Y_b*V_b;
        R_eq_ctu = interp1(IGBT_piecewise(:,1),IGBT_piecewise(:,2),Imag*I_b)./V_b./Imag;
        R_eq = R_eq_ctu*4/pi; % ???
        R_eq(isnan(R_eq))=0;
        R_eq(isinf(R_eq))=0;
        
        %_________________

%% QacVdc nodes

En = 0.995;
Eac = 1;
Edc = 1;
Em = 1.01;

ZINac = [];
[Ydc, YYL, YL, YT, YYT, I_b, Ampacities, y_ih, y_i, A, linedata_dc]  = Ymatrix([1,2,0.2/2,	1E-09,	1E-09,1,100,0,0],A_b,V_b,ZIN);
[Yac, YYL, YL, YT, YYT, I_b, Ampacities, y_ih, y_i, A, linedata_ac]  = Ymatrix([1,2,0.00816,0.003581416,0.0000030159,1,100,0,5;2,3,0.00816,0.003581416,0.0000030159,1,100,0,5],A_b,V_b,ZINac);
Yac3 = cell2mat(arrayfun(@(x) x*eye(n_phases),Yac,'UniformOutput',false));

[Yf, YYL, YL, YT, YYT, I_b, Ampacities, y_ih, y_i, A, linedata]  = Ymatrix([1,2,Rf(1)/Y_b,Xf(1)/Y_b,1E-09,1,100,0,5],A_b,V_b,ZIN);


idx_slack = 1;

% idx1phac.slack = idx_slack;
% idx1phac.pq = [1]'; %[2:26];
% idx1phac.pv = []';%[];idx_qv;
% idx1phac.qv = []';
% idx3phac.slack = idx_slack;
% idx3phac.pq = idx1phac.pq;
% idx3phac.pv = idx1phac.pv;
% idx3phac.qv = idx1phac.qv;
% idx1phdc.slack = idx_slack;
% idx1phdc.pq = [2]'; %[2:26];
% idx1phdc.pv = []';%[];idx_qv;
% idx1phdc.qv = [1]';
% idx3phdc.slack = idx_slack;
% idx3phdc.pq = idx1phdc.pq;
% idx3phdc.pv = idx1phdc.pv;
% idx3phdc.qv = idx1phdc.qv;
% idxCtrlac = [1:2];
% idxCtrldc = [1:2];
% vdep.alpha = repmat([1 0 0],2,1);
% vdep.beta = repmat([1 0 0],2,1);
% vdep.lambda = repmat([0 0 0],2,1);
% vdep.omega = repmat([0 0 0],2,1);
% nph = 1;
% 
% 
% S0ac = [0.995+1i*0.8;1+1i*0.8];
% S0dc = [0.995+1i*0.8;1+1i*0.8];
% Eac = Vacg([1,4]);
% Edc = Vdc([1,4]);


% [Eac(2)*(conj(Yac(2,1))*conj(Eac(1)) + conj(Yac(2,2))*conj(Eac(2))) + real(Iac(4)*conj(Zf(4)*Iac(4)));
% -(Edc(1)*(conj(Ydc(1,2))*conj(Edc(2)) + conj(Ydc(1,1))*conj(Edc(1))) - 1i*Qac(4))]
% 
% 
% Edc(1)*(conj(Yf(1,2))*conj(Edc(2)) + conj(Yf(1,1))*conj(Edc(1)))
% 
% 
% abs(Yac(2,1))^2 * abs(Eac(2)-Eac(1))^2 * conj(Zf)
% 
% 
% 
% Vac(1)*conj() - Eac(2)
% 
% Vac(1)*conj(Iac(4))
% Eac(2)*conj(Iac(4))
% Edc(1)*conj(Idc(4)) - 1i*Qac(4)
% 
% Iac(4)*conj(Zf*Iac(4))
% 
% Iac(4)^2*(Zf)
% 
% conj(conj(Iac(4))*(Zf*Iac(4)))
% Sf = Iac(4)*conj(Zf(4)*Iac(4));
% 
% conj(Zf(4))*conj(Iac(4))*Iac(4)
% abs(Iac(4))^2

% (Eac(2)-Eac(1))/complex(0.00816*Y_b+0*Rf(1),(0.00358141*Y_b+0*Xf(1)))

% [K, Time] = SC_Voltage_simplified(Yac,Ydc,S0ac,S0dc,Eac,Edc,idx1phac,idx3phac,idx1phdc,idx3phdc,idxCtrlac,idxCtrldc,nph,vdep)

M = conj(M(:))';
Vacg = V_complex_LF(:,z1);
% Vac = V_complex_LF(:,z1) + If_complex_LF(:,z1).*complex(Rf,Xf);%.*[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]'; 
Iac = If_complex_LF(:,z1);%.*[1 1 1 1]';
Pac = Nodal_P;
Qac = Nodal_Q;
Vdc = reshape(repmat(Vdc_LF(:,z1)',[1,3]),6,1);
Idc = -reshape(repmat(Idc_flow(:,z1)/2,[1,3])',6,1); 

% 
% alpha = ((conj(Yac(2,1))*conj(Eac(1)) + conj(Yac(2,2))*conj(Eac(2)))*conj(Zf)*Yac(2,2)*Eac(2) +  ...
%         ((Yac(2,1))*(Eac(1)) + (Yac(2,2))*(Eac(2)))*conj(Zf)*conj(Yac(2,2))*conj(Eac(2)) )/abs(Eac(2))
% 
% alpha2 = conj(Zf)*abs(Yac(1,2))^2*(Eac(2)/abs(Eac(2))*(conj(Eac(2))-conj(Eac(1))) + Eac(2)/abs(Eac(2))*(Eac(2)-Eac(1)))
%     
%     
% beta = ((conj(Yac(2,1))*conj(Eac(1)) + conj(Yac(2,2))*conj(Eac(2)))*conj(Zf)*Yac(2,2)*Eac(2) -  ...
%         ((Yac(2,1))*(Eac(1)) + (Yac(2,2))*(Eac(2)))*conj(Zf)*conj(Yac(2,2))*conj(Eac(2)) )
% beta2 = conj(Zf)*abs(Yac(1,2))^2*(Eac(2)/abs(Eac(2))*(conj(Eac(2))-conj(Eac(1))) - Eac(2)/abs(Eac(2))*(Eac(2)-Eac(1)))


S_star = [transpose(complex(Nodal_P(z1,1:3:end), Nodal_Q(z1,1:3:end))); transpose(complex(Pdc_inj(z1,:),0))];
E_star = [V_complex_LF(1:3:end,z1); Vdc_LF(:,z1)];
S_star(3) = (S_star(3));
S_star(2) = (S_star(2));
S_star(1) = (S_star(1));
S_star2 = [transpose(complex(Nodal_P(z2,1:3:end), Nodal_Q(z2,1:3:end))); transpose(complex(Pdc_inj(z2,:),0))];
E_star2 = [V_complex_LF(1:3:end,z2); Vdc_LF(:,z2)];
S_star2(3) = (S_star2(3));
S_star2(2) = (S_star2(2));
S_star2(1) = (S_star2(1));
S0ac = S_star(1:3);
S0dc = S_star(4:5);
Eac = E_star(1:3);
Edc = E_star(4:5);
E = E_star;

S0ac2 = S_star2(1:3);
S0dc2 = S_star2(4:5);
Eac2 = E_star2(1:3);
Edc2 = E_star2(4:5);

E_0 = ones(5,1);
Y = blkdiag(Yac,Ydc);
tol = 1e-7;
n_max = 100;


idx1ph.slack = 1;
idx1ph.pqac = [2]'; %[2:26];
idx1ph.pqdc = [5]'; %[2:26];
idx1ph.vscac = [3]';%[];idx_qv;
idx1ph.vscdc = [4]';%[];idx_qv;

n_ac = 3
n_dc = 2
pos_ac = [3 2];
pos_dc = [4 5];
pos_dc3 = pos_dc +n_ac*2;


S_star(4) = S_star(4);
S_star(5) = S_star(5);
S_star2(4) = S_star2(4);
S_star2(5) = S_star2(5);

% [E,J,n_iter] = NR_rectangularACDC(Y,S_star,E_star,E_0,idx1ph,tol,n_max);
% E_delta = E - E_star
% S_delta = E.*conj(Y*E) - S_star 
% 
% 
% [E2,J2,n_iter2] = NR_rectangularACDC_1ph(Y,S_star,E_star,E_0,idx1ph,tol,n_max,pos_ac,pos_dc);
% E2_delta = E2 - E_star
% S2_delta = E2.*conj(Y*E2) - S_star 


%% 3 phase


Y3 = blkdiag(Yac3,Ydc);
E_star3 = [V_complex_LF(:,z1); Vdc_LF(:,z1)];
S_star3 = [transpose(complex(Nodal_P(z1,:), Nodal_Q(z1,:))); transpose(Pdc_inj(z1,:))];
E_03 = [repmat([1; exp(-1i*pi*2/3);  exp(1i*pi*2/3) ], n_ac,1 ) ; ones(n_dc,1)];
pos_ac3 = polyphase_indices(pos_ac,n_phases);

alp = exp(2*pi/3*1i);
A = 1/3*[1 1     1; 
         1 alp   alp^2; 
         1 alp^2 alp];
ACell =  repmat({A}, 1, n_ac);
ICell =  repmat({1}, 1, n_dc); 
Atot = blkdiag(ACell{:},ICell{:});
     

idx3ph.slack = polyphase_indices(idx_slack,n_phases);
idx3ph.pqac = polyphase_indices(idx1ph.pqac,n_phases);
idx3ph.pqdc =  n_ac*(n_phases -1) + idx1ph.pqdc;
idx3ph.vscac = polyphase_indices(idx1ph.vscac,n_phases);
idx3ph.vscdc = n_ac*(n_phases -1) + idx1ph.vscdc;

[E3,J3,n_iter] = NR_rectangularACDC_3ph(Y3,S_star3,E_star3,E_03,idx3ph,tol,n_max,pos_ac3,pos_dc3);
E_delta = E3 - E_star3
S_delta = E3.*conj(Y3*E3) - S_star3

Atot*E_star3.*conj(Atot*Y3*E_star3)

     real(A*E_star3(7:9,:).*conj(A*Y3(7:9,:)*E_star3) )
     -[0; E_star3(10) * conj(Y3(10,:)*E_star3); 0]
     
     1/3*(1*E_star3(7)+1*E_star3(8)+1*E_star3(9))* 1/3*conj(1*(Y3(7,7)*E_star3(7)+Y3(4,4)*(V2a_re + 1i*V2a_im)) +1*(Y3(8,8)*E_star3(8)+Y3(5,5)*(V2b_re + 1i*V2b_im)) + 1*(Y3(9,9)*E_star3(9)+Y3(6,6)*(V2c_re + 1i*V2c_im)))
     1/3*(1*E_star3(7)+alp*E_star3(8)+alp^2*E_star3(9))* 1/3*conj(1*(Y3(7,7)*E_star3(7)+Y3(7,4)*(V2a_re + 1i*V2a_im)) +alp*(Y3(8,8)*E_star3(8)+Y3(8,5)*(V2b_re + 1i*V2b_im)) + alp^2*(Y3(9,9)*E_star3(9)+Y3(9,6)*(V2c_re + 1i*V2c_im)))
     1/3*(1*E_star3(7)+alp^2*E_star3(8)+alp*E_star3(9))* 1/3*conj(1*(Y3(7,7)*E_star3(7)+Y3(7,4)*(V2a_re + 1i*V2a_im)) +alp^2*(Y3(8,8)*E_star3(8)+Y3(8,5)*(V2b_re + 1i*V2b_im)) + alp*(Y3(9,9)*E_star3(9)+Y3(9,6)*(V2c_re + 1i*V2c_im)))


     1/3*(E_star3(7)+alp^2*E_star3(8)+alp*E_star3(9))* 1/3*conj(Y3(7,7)*E_star3(7)+alp^2*Y3(8,8)*E_star3(8)+alp*Y3(9,9)*E_star3(9))
     S_star3(1:3,:)
     
     
G3 = real(Y3);
B3 = imag(Y3);
alpha = (G3(10,11)* real(E_star3(11,:))).^2 - 4 * G3(11,11) * real(A*E_star3(7:9,:).*conj(A*Y3(7:9,:)*E_star3));

-Y3(10,11)/(2*Y3(10,10)).* E_star3(11,:) + sqrt(alpha)/(2*Y3(10,10))

P2 == real(E_star3(4:6,:) .* conj(Y3(4:6,1:3)*E_star3(1:3,:) + Y3(4:6,4:6)*E_star3(4:6,:) + Y3(4:6,7:9)*E_star3(7:9,:)))
Q2 == imag(E_star3(4:6,:) .* conj(Y3(4:6,1:3)*E_star3(1:3,:) + Y3(4:6,4:6)*E_star3(4:6,:) + Y3(4:6,7:9)*E_star3(7:9,:)))
Q3 == imag(E_star3(7:9,:) .* conj(Y3(7:9,1:3)*E_star3(1:3,:) + Y3(7:9,4:6)*E_star3(4:6,:) + Y3(7:9,7:9)*E_star3(7:9,:)))

P5 == real(E_star3(11,:) .* conj(Y3(11,10)*E_star3(10,:) + Y3(11,11)*E_star3(11,:) ))

V4 == -Y3(10,11)/(2*Y3(10,10)).* E_star3(11,:) + sqrt(alpha)/(2*Y3(10,10))

     
x = 0.8:0.001:1.2;
figure
hold on
scatter(x,abs(E2(2,:))')
scatter(x,abs(E2(3,:))')
scatter(x,abs(E2(5,:))')

E_star.*conj(Y*E_star) 


Zf = (0.008+0.04*1i)*Y_b;
Pf = Iac(7)*Zf*Zf

(E_star(2) - E_star(3))*Y(3,2)*Zf*Zf
real(S_star(3) + S_star(4) - Pf)
 

real(E_star(3).*conj(Y(3,:)*E_star) + E_star(4).*conj(Y(4,:)*E_star) + Zf*Zf*Y(3,:)*E_star)
beta = real(E_star(3).*conj(Y(3,:)*E_star) + Zf*Zf*Y(3,:)*E_star)

real((-Y(4,5)*E_star(5) + sqrt((Y(4,5)*E_star(5))^2 - 4*Y(4,4)*beta))/2/Y(4,4))

G = real(Y);
B = imag(Y);
1/(2*G(4,4)) * (1/(2*sqrt((Y(4,5)*E_star(5))^2 - 4*Y(4,4)*beta))) * (-4*G(4,4)) * (real(E_star(3)) * G(3,2) - imag(E_star(3)) * B(3,2)) 


J2 = imag((S_star(3) - S_star2(3)))./(E_star - E_star2);
J7 = [real(J2(2)) real(J2(3)) J2(5) imag(J2(2)) imag(J2(3))]

J2r = real((S_star(2) - S_star2(2)))./real(E_star - E_star2);
J2i = real((S_star(2) - S_star2(2)))./imag(E_star - E_star2);
J7 = [J2r(2) J2r(3) J2r(5) J2i(2) J2i(3)]

J2r = real((E_star(4) - E_star2(4)))./real(E_star - E_star2);
J2i = real((E_star(4) - E_star2(4)))./imag(E_star - E_star2);
J7 = [J2r(2) J2r(3) J2r(5) J2i(2) J2i(3)]

J2r = real((E_star(4) - E_star2(4)))./(E_star - E_star2);
J7 = [real(J2(2)) real(J2(3)) J2(5) imag(J2(2)) imag(J2(3))]


(E_star(4) - E_star2(4)) / real(E_star(2) - E_star2(2))


%% second start

E_s = E
S_s = E.*conj(Y*E);
alpha = (G(4,5)* real(E_s(5,:))).^2 - 4 * G(5,5) * real(S_s(3,:));
-Y(4,5)/(2*Y(4,4)).* E_s(5,:) + sqrt(alpha)/(2*Y(4,4))


P2 == real(E_s(2,:) .* conj(Y(2,1)*E_s(1,:) + Y(2,2)*E_s(2,:) + Y(2,3)*E_s(3,:)))
Q2 == imag(E_s(2,:) .* conj(Y(2,1)*E_s(1,:) + Y(2,2)*E_s(2,:) + Y(2,3)*E_s(3,:)))
Q3 == imag(E_s(3,:) .* conj(Y(3,1)*E_s(1,:) + Y(3,2)*E_s(2,:) + Y(3,3)*E_s(3,:)))

P5 == real(E_s(5,:) .* conj(Y(5,4)*E_s(4,:) + Y(5,5)*E_s(5,:) ))

V4 == -Y(4,5)/(2*Y(4,4)).* E_s(5,:) + sqrt(alpha)/(2*Y(4,4))

G = real(Y);
B = imag(Y); 
E_re = real(E_star);
E_im = imag(E_star);

J11 = G(2,1)*E_re(1) - B(2,1)*E_im(1) +2*G(2,2)*E_re(2) + G(2,3)*E_re(3) - B(2,3)*E_im(3)
J12 = G(2,3)*E_re(2) + B(2,3)*E_im(2)
J13 = 0
J14 = B(2,1)*E_re(1) + G(2,1)*E_im(1) + 2*G(2,2)*E_im(2) + B(2,3)*E_re(3) + B(2,3)*E_im(3)
J15 = -B(2,3)*E_re(2) + G(2,3)*E_im(2)

%% SC comparison
Fl = [2,3]; 
idxCtrl = [1:5];
vdep.alpha = repmat([1 0 0],5,1);
vdep.beta = repmat([1 0 0],5,1);
vdep.lambda = repmat([0 0 0],5,1);
vdep.omega = repmat([0 0 0],5,1);
nph = 1;
filter = 0; unblanced_3ph = 0;
[K, Time] = SC_Voltage_V5(Yac,S_star(1:3),E_star(1:3),Ydc,S_star(4:5),E_star(4:5),idx1ph,idx1ph,idxCtrl,nph,vdep,Zf,Fl,unblanced_3ph,filter);

n_nodes = 5;
J_PR = zeros(n_nodes);
J_PX = zeros(n_nodes);
J_QR = zeros(n_nodes);
J_QX = zeros(n_nodes);
J_VR = zeros(n_nodes);
J_VX = zeros(n_nodes);
for k = 1:size(K,1)
    if( sum( K{k,1} == idx1ph.pqac ))
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
        J_QR(:,k) = real(K{k,2}{2,1});
        J_QX(:,k) = imag(K{k,2}{2,1});
    elseif( sum( K{k,1} == idx1ph.vscac ))
        J_QR(:,k) = real(K{k,2}{1,1});
        J_QX(:,k) = imag(K{k,2}{1,1});
    elseif( sum( K{k,1} == idx1ph.vscdc ))
        J_VR(:,k) = real(K{k,2}{1,1});
        J_VX(:,k) = imag(K{k,2}{1,1});
        J_QR(:,k) = real(K{k,2}{2,1});
        J_QX(:,k) = imag(K{k,2}{2,1});
    elseif( sum( K{k,1} == idx1ph.pqdc ) )
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
    end
end

J_PR(:,setdiff(1:n_nodes,[idx1ph.pqac; idx1ph.pqdc])) = [];
J_PX(:,setdiff(1:n_nodes,[idx1ph.pqac; idx1ph.pqdc])) = [];
J_QR(:,setdiff(1:n_nodes,[idx1ph.pqac; idx1ph.vscac])) = [];
J_QX(:,setdiff(1:n_nodes,[idx1ph.pqac; idx1ph.vscac])) = [];
J_VR(:,setdiff(1:n_nodes,[idx1ph.vscdc])) = [];
J_VX(:,setdiff(1:n_nodes,[idx1ph.vscdc])) = [];

J_PR(idx1ph.slack,:) = [];
J_PX(idx1ph.slack,:) = [];
J_QR(idx1ph.slack,:) = [];
J_QX(idx1ph.slack,:) = [];
J_VR(idx1ph.slack,:) = [];
J_VX(idx1ph.slack,:) = [];


%% Gurobi
G = real(Y);
B = imag(Y);
E_re = real(E_star);
E_im = imag(E_star);
V2_re = sdpvar(1,1, 'full'); %11
V2_im = sdpvar(1,1, 'full'); %11
V3_re = sdpvar(1,1, 'full'); %11
V3_im = sdpvar(1,1, 'full'); %11
V5_re = sdpvar(1,1, 'full'); %11
% alpha = sdpvar(1,1, 'full'); %11

P2 = real(S_star(2));
Q2 = imag(S_star(2));
Q3 = imag(S_star(3));
V4_re = real(E_star(4));
P5 = real(S_star(5));

% syms V2_re V2_im V3_re V3_im V5_re

% V2_re = real(E_star(2)); %11
% V2_im = imag(E_star(2)); %11
% V3_re = real(E_star(3)); %11
% V3_im = imag(E_star(3)); %11
% V5_re = real(E_star(5)); %11



con1 = [P2 == V2_re*( G(2,1)*E_re(1) + G(2,2)*V2_re + G(2,3)*V3_re - B(2,1)*E_im(1) - B(2,2)*V2_im - B(2,3)*V3_im ) + ...
                V2_im*(B(2,1)*E_re(1) + B(2,2)*V2_re + B(2,3)*V3_re + G(2,1)*E_im(1) + G(2,2)*V2_im + G(2,3)*V3_im) ] ; 
            
con2 = [Q2 == V2_re*( -B(2,1)*E_re(1) - B(2,2)*V2_re - B(2,3)*V3_re - G(2,1)*E_im(1) - G(2,2)*V2_im - G(2,3)*V3_im ) + ...
                V2_im*(G(2,1)*E_re(1) + G(2,2)*V2_re + G(2,3)*V3_re - B(2,1)*E_im(1) - B(2,2)*V2_im - B(2,3)*V3_im) ] ; 

con3 = [ V3_re*(G(3,2)*V2_re + G(3,3)*V3_re - B(3,2)*V2_im - B(3,3)*V3_im) + ...
         V3_im*(B(3,2)*V2_re + B(3,3)*V3_re + G(3,2)*V2_im + G(3,3)*V3_im) == -V4_re * ( G(4,4)*V4_re + G(4,5)*V5_re)   ] ;  
 
con6 = [V4_re == -G(4,5)/(2*G(4,4))* V5_re + sqrt((G(4,5)* V5_re)^2 - 4 * G(5,5) * (V3_re*(G(3,2)*V2_re + G(3,3)*V3_re - B(3,2)*V2_im - B(3,3)*V3_im) + ...
                                              V3_im*(B(3,2)*V2_re + B(3,3)*V3_re + G(3,2)*V2_im + G(3,3)*V3_im)))/(2*G(4,4))];
             
con4 = [ Q3 == V3_re*( -B(3,2)*V2_re - B(3,3)*V3_re - G(3,2)*V2_im - G(3,3)*V3_im) + ...
                 V3_im*(G(3,2)*V2_re + G(3,3)*V3_re - B(3,2)*V2_im - B(3,3)*V3_im)];

con5 = [ P5 == V5_re * ( G(5,4)*V4_re + G(5,5)*V5_re)   ] ;  


% alpha = (G(4,5)* E_star(5))^2 - 0.5*4 * G(5,5) * real(S_star(3));
% 
% P2 == real(E_star(2) * conj(Y(2,1)*E_star(1) + Y(2,2)*E_star(2) + Y(2,3)*E_star(3)))
% Q2 == imag(E_star(2) * conj(Y(2,1)*E_star(1) + Y(2,2)*E_star(2) + Y(2,3)*E_star(3)))
% Q3 == imag(E_star(3) * conj(Y(3,1)*E_star(1) + Y(3,2)*E_star(2) + Y(3,3)*E_star(3)))
% P5 == 2*real(E_star(5) * conj(Y(5,4)*E_star(4) + Y(5,5)*E_star(5) ))
% V4 == -Y(4,5)/(2*Y(4,4))* E_star(5) + sqrt(alpha)/(2*Y(4,4))


con = [con1, con2, con3, con4, con5]%, 0.8<=V2_re<= 1.2,0.8<=V3_re<= 1.2,0.8<=V5_re<= 1.2]
% 
% x0.V2_re = 1; 
% x0.V2_im = 0;
% x0.V3_re = 1; 
% x0.V3_im = 0; 
% x0.V5_re = 1; 
% 
% prob.Constraints.con1 = con1; 
% prob.Constraints.con2 = con2;
% prob.Constraints.con3 = con3; 
% prob.Constraints.con4 = con4; 
% prob.Constraints.con5 = con5; 

cost = 0;
% cost = sum(abs(1-E_star));
%% Solve the problem
assign(V2_re,1)
assign(V3_re,1)
assign(V5_re,1)
assign(V2_im,0)
assign(V3_im,0)


options =  sdpsettings('verbose', 1, 'usex0',1);
optimize(con, cost, options);


E_opt = [E_star(1); complex(double(V2_re),double(V2_im)); complex(double(V3_re),double(V3_im)); V4_re; double(V5_re) ]

E_opt.*conj(Y*E_opt) 













%% partial derivative
G = real(Y);
B = imag(Y);
E_re = real(E_star);
E_im = imag(E_star);
syms V2_re V2_im V3_re V3_im V5_re

P2 = V2_re*( G(2,1)*E_re(1) + G(2,2)*V2_re + G(2,3)*V3_re - B(2,1)*E_im(1) - B(2,2)*V2_im - B(2,3)*V3_im ) + ...
                V2_im*(B(2,1)*E_re(1) + B(2,2)*V2_re + B(2,3)*V3_re + G(2,1)*E_im(1) + G(2,2)*V2_im + G(2,3)*V3_im) ; 
            
Q2 = V2_re*( -B(2,1)*E_re(1) - B(2,2)*V2_re - B(2,3)*V3_re - G(2,1)*E_im(1) - G(2,2)*V2_im - G(2,3)*V3_im ) + ...
                V2_im*(G(2,1)*E_re(1) + G(2,2)*V2_re + G(2,3)*V3_re - B(2,1)*E_im(1) - B(2,2)*V2_im - B(2,3)*V3_im)  ; 

V4 = -G(4,5)/(2*G(4,4))* V5_re + sqrt((G(4,5)* V5_re)^2 - 4 * G(5,5) * (V3_re*(G(3,2)*V2_re + G(3,3)*V3_re - B(3,2)*V2_im - B(3,3)*V3_im) + ...
                                              V3_im*(B(3,2)*V2_re + B(3,3)*V3_re + G(3,2)*V2_im + G(3,3)*V3_im)))/(2*G(4,4));
             
Q3 = V3_re*( -B(3,2)*V2_re - B(3,3)*V3_re - G(3,2)*V2_im - G(3,3)*V3_im) + ...
                 V3_im*(G(3,2)*V2_re + G(3,3)*V3_re - B(3,2)*V2_im - B(3,3)*V3_im);

P5 = V5_re * ( G(5,4)*E_re(4) + G(5,5)*V5_re)   ;  


J_s = [diff(P2,V2_re) diff(P2,V3_re) diff(P2,V5_re) diff(P2,V2_im) diff(P2,V3_im);
      diff(Q2,V2_re) diff(Q2,V3_re) diff(Q2,V5_re) diff(Q2,V2_im) diff(Q2,V3_im)
      diff(V4,V2_re) diff(V4,V3_re) diff(V4,V5_re) diff(V4,V2_im) diff(V4,V3_im)
      diff(Q3,V2_re) diff(Q3,V3_re) diff(Q3,V5_re) diff(Q3,V2_im) diff(Q3,V3_im)
      diff(P5,V2_re) diff(P5,V3_re) diff(P5,V5_re) diff(P5,V2_im) diff(P5,V3_im)]
  
J_sub = subs(J_s, [V2_re  V3_re  V5_re V2_im V3_im], [real(E(2)) real(E(3)) real(E(5)) imag(E(2)) imag(E(2))] )
double(J_sub)


con3 = [ V3_re*(G(3,2)*V2_re + G(3,3)*V3_re - B(3,2)*V2_im - B(3,3)*V3_im) + ...
         V3_im*(B(3,2)*V2_re + B(3,3)*V3_re + G(3,2)*V2_im + G(3,3)*V3_im) == -V4_re * ( G(4,4)*V4_re + G(4,5)*V5_re)   ] ;  
 
diff(con3,V5_re)


G3 = real(Y3);
B3 = imag(Y3);
E_re = real(E_star3);
E_im = imag(E_star3);
syms V2a_re V2a_im V2b_re V2b_im V2c_re V2c_im V3a_re V3a_im V3b_re V3b_im V3c_re V3c_im V5_re


I0 = real( 1/3*(1*(Y3(7,7)*(V3a_re + 1i*V3a_im)+Y3(4,4)*(V2a_re + 1i*V2a_im)) +1*    (Y3(8,8)*(V3b_re + 1i*V3b_im)+Y3(5,5)*(V2b_re + 1i*V2b_im)) + 1*     (Y3(9,9)*(V3c_re + 1i*V3c_im)+Y3(6,6)*(V2c_re + 1i*V2c_im))));

P0 = real(1/3*(1*(V3a_re + 1i*V3a_im)+1*(V3b_re + 1i*V3b_im)+1*(V3c_re + 1i*V3c_im))*       1/3*conj(1*(Y3(7,7)*(V3a_re + 1i*V3a_im)+Y3(4,4)*(V2a_re + 1i*V2a_im)) +1*    (Y3(8,8)*(V3b_re + 1i*V3b_im)+Y3(5,5)*(V2b_re + 1i*V2b_im)) + 1*     (Y3(9,9)*(V3c_re + 1i*V3c_im)+Y3(6,6)*(V2c_re + 1i*V2c_im))));
Pp = real(1/3*(1*(V3a_re + 1i*V3a_im)+alp*(V3b_re + 1i*V3b_im)+alp^2*(V3c_re + 1i*V3c_im))* 1/3*conj(1*(Y3(7,7)*(V3a_re + 1i*V3a_im)+Y3(7,4)*(V2a_re + 1i*V2a_im)) +alp*  (Y3(8,8)*(V3b_re + 1i*V3b_im)+Y3(8,5)*(V2b_re + 1i*V2b_im)) + alp^2* (Y3(9,9)*(V3c_re + 1i*V3c_im)+Y3(9,6)*(V2c_re + 1i*V2c_im))));
Pn = real(1/3*(1*(V3a_re + 1i*V3a_im)+alp^2*(V3b_re + 1i*V3b_im)+alp*(V3c_re + 1i*V3c_im))* 1/3*conj(1*(Y3(7,7)*(V3a_re + 1i*V3a_im)+Y3(7,4)*(V2a_re + 1i*V2a_im)) +alp^2*(Y3(8,8)*(V3b_re + 1i*V3b_im)+Y3(8,5)*(V2b_re + 1i*V2b_im)) + alp*   (Y3(9,9)*(V3c_re + 1i*V3c_im)+Y3(9,6)*(V2c_re + 1i*V2c_im))));

V40 = -G3(10,11)/(2*G3(10,10))* V5_re + sqrt((G3(10,11)* V5_re)^2 - 4 * G3(11,11) * Pp)/(2*G3(10,10));
                                          
                                          
J3 = [diff(P0,V2a_re) diff(P0,V2b_re) diff(P0,V2c_re) diff(P0,V3a_re) diff(P0,V3b_re) diff(P0,V3c_re) diff(P0,V5_re) diff(P0,V2a_im) diff(P0,V2b_im) diff(P0,V2c_im) diff(P0,V3a_im) diff(P0,V3b_im) diff(P0,V3c_im) ;
      diff(V40,V2a_re) diff(V40,V2b_re) diff(V40,V2c_re) diff(V40,V3a_re) diff(V40,V3b_re) diff(V40,V3c_re) diff(V40,V5_re) diff(V40,V2a_im) diff(V40,V2b_im) diff(V40,V2c_im) diff(V40,V3a_im) diff(V40,V3b_im) diff(V40,V3c_im);
      diff(Pn,V2a_re) diff(Pn,V2b_re) diff(Pn,V2c_re) diff(Pn,V3a_re) diff(Pn,V3b_re) diff(Pn,V3c_re) diff(Pn,V5_re) diff(Pn,V2a_im) diff(Pn,V2b_im) diff(Pn,V2c_im) diff(Pn,V3a_im) diff(Pn,V3b_im) diff(Pn,V3c_im) ]
     
J_sub3 = subs(J3, [V2a_re V2a_im V2b_re V2b_im V2c_re V2c_im V3a_re V3a_im V3b_re V3b_im V3c_re V3c_im V5_re], [real(E_star3(4)) imag(E_star3(4))  real(E_star3(5)) imag(E_star3(5)) real(E_star3(6)) imag(E_star3(6)) real(E_star3(7)) imag(E_star3(7)) real(E_star3(8)) imag(E_star3(8)) real(E_star3(9)) imag(E_star3(9)) real(E_star3(11))] )
J3_solved = double(J_sub3)    
    
double(subs(diff(P0,Va_re), [V2a_re V2a_im V2b_re V2b_im V2c_re V2c_im V3a_re V3a_im V3b_re V3b_im V3c_re V3c_im V5_re], [real(E_star3(4)) imag(E_star3(4))  real(E_star3(5)) imag(E_star3(5)) real(E_star3(6)) imag(E_star3(6)) real(E_star3(7)) imag(E_star3(7)) real(E_star3(8)) imag(E_star3(8)) real(E_star3(9)) imag(E_star3(9)) real(E_star3(11))] ))

Ji = [diff(I0,V2a_re) diff(I0,V2b_re) diff(I0,V2c_re) diff(I0,V3a_re) diff(I0,V3b_re) diff(I0,V3c_re) diff(I0,V5_re) diff(I0,V2a_im) diff(I0,V2b_im) diff(I0,V2c_im) diff(I0,V3a_im) diff(I0,V3b_im) diff(I0,V3c_im) ]
J_subi = subs(Ji, [V2a_re V2a_im V2b_re V2b_im V2c_re V2c_im V3a_re V3a_im V3b_re V3b_im V3c_re V3c_im V5_re], [real(E_star3(4)) imag(E_star3(4))  real(E_star3(5)) imag(E_star3(5)) real(E_star3(6)) imag(E_star3(6)) real(E_star3(7)) imag(E_star3(7)) real(E_star3(8)) imag(E_star3(8)) real(E_star3(9)) imag(E_star3(9)) real(E_star3(11))] )
Ji_solved = double(J_subi)    


Jabc = [diff(V40,V2a_re) + 1i*diff(V40,V2a_im); diff(V40,V2b_re) + 1i*diff(V40,V2b_im); diff(V40,V2c_re) + 1i*diff(V40,V2c_im)] 
Jabc_solved = double(subs(Jabc, [V2a_re V2a_im V2b_re V2b_im V2c_re V2c_im V3a_re V3a_im V3b_re V3b_im V3c_re V3c_im V5_re], [real(E_star3(4)) imag(E_star3(4))  real(E_star3(5)) imag(E_star3(5)) real(E_star3(6)) imag(E_star3(6)) real(E_star3(7)) imag(E_star3(7)) real(E_star3(8)) imag(E_star3(8)) real(E_star3(9)) imag(E_star3(9)) real(E_star3(11))] ))  

% symetrical components
syms V20_re V20_im V2p_re V2p_im V2n_re V2n_im V30_re V30_im V3p_re V3p_im V3n_re V3n_im V5_re

E_star3_sym = blkdiag(A,A,A,1,1) * E_star3;

P0 = real((V30_re + 1i*V30_im)* conj((Y3(7,7)*(V30_re + 1i*V30_im)+Y3(4,4)*(V20_re + 1i*V20_im))));
Pp = real((V3p_re + 1i*V3p_im)* conj((Y3(8,8)*(V3p_re + 1i*V3p_im)+Y3(8,5)*(V2p_re + 1i*V2p_im))));
Pn = real((V3n_re + 1i*V3n_im)* conj((Y3(9,9)*(V3n_re + 1i*V3n_im)+Y3(9,6)*(V2n_re + 1i*V2n_im))));


double(subs(P0, [V20_re V20_im V2p_re V2p_im V2n_re V2n_im V30_re V30_im V3p_re V3p_im V3n_re V3n_im V5_re], [real(E_star3_sym(4)) imag(E_star3_sym(4))  real(E_star3_sym(5)) imag(E_star3_sym(5)) real(E_star3_sym(6)) imag(E_star3_sym(6)) real(E_star3_sym(7)) imag(E_star3_sym(7)) real(E_star3_sym(8)) imag(E_star3_sym(8)) real(E_star3_sym(9)) imag(E_star3_sym(9)) real(E_star3_sym(11))] ))


V40 = -G3(10,11)/(2*G3(10,10))* V5_re + sqrt((G3(10,11)* V5_re)^2 - 4 * G3(11,11) * Pp)/(2*G3(10,10));
         
J0pn = [diff(V40,V20_re) + 1i*diff(V40,V20_im); diff(V40,V2p_re) + 1i*diff(V40,V2p_im); diff(V40,V2n_re) + 1i*diff(V40,V2n_im)] 
J0pn_solved = 1/3*double(subs(J0pn, [V20_re V20_im V2p_re V2p_im V2n_re V2n_im V30_re V30_im V3p_re V3p_im V3n_re V3n_im V5_re], [real(E_star3_sym(4)) imag(E_star3_sym(4))  real(E_star3_sym(5)) imag(E_star3_sym(5)) real(E_star3_sym(6)) imag(E_star3_sym(6)) real(E_star3_sym(7)) imag(E_star3_sym(7)) real(E_star3_sym(8)) imag(E_star3_sym(8)) real(E_star3_sym(9)) imag(E_star3_sym(9)) real(E_star3_sym(11))] ))

inv(A) * J0pn_solved


J0pn = [diff(P0,V20_re) + 1i*diff(P0,V20_im); diff(P0,V2p_re) + 1i*diff(P0,V2p_im); diff(P0,V2n_re) + 1i*diff(P0,V2n_im)] 

double(subs(J0pn, [V20_re V20_im V2p_re V2p_im V2n_re V2n_im V30_re V30_im V3p_re V3p_im V3n_re V3n_im V5_re], [real(E_star3_sym(4)) imag(E_star3_sym(4))  real(E_star3_sym(5)) imag(E_star3_sym(5)) real(E_star3_sym(6)) imag(E_star3_sym(6)) real(E_star3_sym(7)) imag(E_star3_sym(7)) real(E_star3_sym(8)) imag(E_star3_sym(8)) real(E_star3_sym(9)) imag(E_star3_sym(9)) real(E_star3_sym(11))] ))


%%



opts = optimoptions('intlinprog','Heuristics','none');
S = solve(con,[V2_re V2_im V3_re V3_im V5_re])


[sol2,fval2,exitstatus2,output2] = solve(con,[V2_re V2_im V3_re V3_im V5_re])% ,x0,'options',opts)

time = tic;
optimize(con, cost, options);
TT = toc(time);


E_opt = [E_star(1); complex(double(V2_re),double(V2_im)); complex(double(V3_re),double(V3_im)); V4_re; double(V5_re) ]




%% visualisation

V5 = [0.8:0.01:1.2]';
Pac = [-2.56:0.1:4];
% V4 = @(V5,Pac) (0.5*(V5+sqrt(V5.^2 + 0.25*Pac)));
figure
% plot3(V5,Pac,V4)

 V4p = (0.5*(V5+sqrt(V5.^2 + 0.25*Pac)));
 V4n = (0.5*(V5-sqrt(V5.^2 + 0.25*Pac)));
 
 [Xm,Ym] = ndgrid(V5,Pac);
 hold on
 surf(Xm,Ym,V4p)
 surf(Xm,Ym,V4n)
 surf([0.8,1.2],[-2.5,4],0.8*ones(2),'FaceColor',[1 0 0], 'FaceAlpha',0.9)
  surf([0.8,1.2],[-2.5,4],1.2*ones(2),'FaceColor',[1 0 0], 'FaceAlpha',0.9)
 xlabel('V5')
 ylabel('Pac')
 zlabel('V4')
 
 V1_re = 1
 V1_im = 0
con1 = @(V2_re,V2_im,V4_re) V2_re*( G(2,1)*V1_re + G(2,2)*V2_re - B(2,1)*V1_im - B(2,2)*V2_im) + ...
         V2_im*(B(2,1)*V1_re + B(2,2)*V2_re + G(2,1)*V1_im + G(2,2)*V2_im) - V3_re * ( G(3,3)*V3_re + G(3,2)*V2_re);  

con2 = @(V2_re,V2_im,V4_re) -V2_re*( B(2,1)*V1_re + B(2,2)*V2_re + G(2,1)*V1_im + G(2,2)*V2_im) + ...
          V2_im*(G(2,1)*V1_re + G(2,1)*V2_re - B(2,1)*V1_im - B(2,2)*V2_im) - Q3;

con3 = @(V2_re,V2_im,V4_re) V4_re * ( G(4,3)*V3_re + G(4,4)*V4_re) - P4;  

V2_re = [0.8:0.01:1.2]';
V2_im = [-0.5:0.1:0.5];
V5_re = 1;
% V4 = @(V5,Pac) (0.5*(V5+sqrt(V5.^2 + 0.25*Pac)));
figure
% plot3(V5,Pac,V4)


 
 [Xm,Ym] = ndgrid(V2_re,V2_im);
 hold on
 surf(Xm,Ym,con2(V2_re,V2_im,V4_re))


Fac = diag(Eac)*conj(Yac)*diag(conj(Eac));
Fdc = diag(Edc)*conj(Ydc)*diag(conj(Edc));
F = blkdiag(Fac,Fdc);


idx_slack = 1;
idx1ph.slack = idx_slack;
idx1ph.pqac = [2]'; %[2:26];
idx1ph.pqdc = [5]'; %[2:26];
idx1ph.vscac = [3]';%[];idx_qv;
idx1ph.vscdc = [4]';%[];idx_qv;
idx1ph.qv = []';
idx3ph.slack = idx_slack;
idx3ph.pqac = idx1ph.pqac;
idx3ph.pqdc = idx1ph.pqdc;
idx3ph.vscac = idx1ph.vscac;
idx3ph.vscdc = idx1ph.vscdc;
idx3ph.qv = idx1ph.qv;
idxCtrl = [1:5];
vdep.alpha = repmat([1 0 0],5,1);
vdep.beta = repmat([1 0 0],5,1);
vdep.lambda = repmat([0 0 0],5,1);
vdep.omega = repmat([0 0 0],5,1);
nph = 1;

% A11 = sum(real(Fac(idx1ph.pqac,:)))./abs(Eac(idx1ph.pqac)) +
% real(Fac(idx1ph.pqac,idx1ph.pqac))./abs(Eac(idx1ph.pqac)); A12 =
% -(sum(imag(Fac(idx1ph.pqac,:))) - imag(Fac(idx1ph.pqac,idx1ph.pqac)));
% A13 = real(Fac(idx1ph.pqac,3))./abs(Eac(3)); A14 =
% -(-imag(Fac(idx1ph.pqac,3))); A15 = 0;
% 
% A21 = sum(imag(Fac(idx1ph.pqac,:)))./abs(Eac(idx1ph.pqac)) +
% imag(Fac(idx1ph.pqac,idx1ph.pqac))./abs(Eac(idx1ph.pqac)); A22 =
% (sum(real(Fac(idx1ph.pqac,:))) - real(Fac(idx1ph.pqac,idx1ph.pqac))); A23
% = imag(Fac(idx1ph.pqac,3))./abs(Eac(3)); A24 = -real(Fac(idx1ph.pqac,3));
% A25 = 0;
% 
% A31 = real(Fac(3,idx1ph.pqac))./abs(Eac(idx1ph.pqac)); A32 =
% -(-imag(Fac(3,idx1ph.pqac))); A33 = sum(real(Fac(3,:)))./abs(Eac(3)) +
% real(Fac(3,3))./abs(Eac(3)); A34 = -(sum(imag(Fac(3,:))) -
% imag(Fac(3,3))); A35 = --real(Fdc(1,2))./abs(Edc(2));
% 
% A41 = imag(Fac(3,idx1ph.pqac))./abs(Eac(idx1ph.pqac)); A42 =
% -real(Fac(3,idx1ph.pqac)); A43 = sum(imag(Fac(3,:)))./abs(Eac(3)) +
% imag(Fac(3,3))./abs(Eac(3)); A44 = (sum(real(Fac(3,:))) -
% real(Fac(3,3))); A45 = 0;
% 
% A51 = 0; A52 = 0; A53 = 0; A54 = 0; A55 = -real(Fdc(2,1))./abs(Edc(2));

alpha = (conj(Yac(3,2))*conj(Eac(2)) + conj(Yac(3,3))*conj(Eac(3)))*conj(Zf(1));
beta = (Yac(3,2)*Eac(2) + Yac(3,3)*Eac(3))*conj(Zf(1));
a = real(alpha*Yac(3,2)*Eac(2) + beta*conj(Yac(3,2))*conj(Eac(2)));
b = real(alpha*Yac(3,3)*Eac(3) + beta*conj(Yac(3,3))*conj(Eac(3)));
c = -imag(alpha*Yac(3,2)*Eac(2) + beta*conj(Yac(3,2))*conj(Eac(2)));
d = -imag(alpha*Yac(3,3)*Eac(3) + beta*conj(Yac(3,3))*conj(Eac(3)));
a = 0;b = 0;c=0;d=0;

A11 = sum(real(F(idx1ph.pqac,:)))./abs(E(idx1ph.pqac)) + real(F(idx1ph.pqac,idx1ph.pqac))./abs(E(idx1ph.pqac));
A12 = -(sum(imag(F(idx1ph.pqac,:))) - imag(F(idx1ph.pqac,idx1ph.pqac)));
A13 = real(F(idx1ph.pqac,idx1ph.vscac))./abs(E(idx1ph.vscac));
A14 = -(-imag(F(idx1ph.pqac,idx1ph.vscac)));
A15 = 0;


A21 = sum(imag(F(idx1ph.pqac,:)))./abs(E(idx1ph.pqac)) + imag(F(idx1ph.pqac,idx1ph.pqac))./abs(E(idx1ph.pqac));
A22 = (sum(real(F(idx1ph.pqac,:))) - real(F(idx1ph.pqac,idx1ph.pqac)));
A23 = imag(F(idx1ph.pqac,idx1ph.vscac))./abs(E(idx1ph.vscac));
A24 = -real(F(idx1ph.pqac,idx1ph.vscac));
A25 = 0;

A31 = real(F(idx1ph.vscac,idx1ph.pqac))./abs(E(idx1ph.pqac));
A32 = -(-imag(F(idx1ph.vscac,idx1ph.pqac)));
A33 = sum(real(F(idx1ph.vscac,:)))./abs(E(idx1ph.vscac)) + real(F(idx1ph.vscac,idx1ph.vscac))./abs(E(idx1ph.vscac));
A34 = -(sum(imag(F(idx1ph.vscac,:))) - imag(F(idx1ph.vscac,idx1ph.vscac)));
A35 = --real(F(idx1ph.vscdc,idx1ph.pqdc))./abs(E(idx1ph.pqdc));

A41 = imag(F(idx1ph.vscac,idx1ph.pqac))./abs(E(idx1ph.pqac));
A42 = -real(F(idx1ph.vscac,idx1ph.pqac));
A43 = sum(imag(F(idx1ph.vscac,:)))./abs(E(idx1ph.vscac)) + imag(F(idx1ph.vscac,idx1ph.vscac))./abs(E(idx1ph.vscac));
A44 = (sum(real(F(idx1ph.vscac,:))) - real(F(idx1ph.vscac,idx1ph.vscac)));
A45 = 0;

A51 = 0;
A52 = 0;
A53 = 0;
A54 = 0;
A55 = -real(F(idx1ph.pqdc,idx1ph.vscdc))./abs(E(idx1ph.pqdc));

A = [ A11 A12 A13 A14 A15; ...
      A21 A22 A23 A24 A25; ...
      A31 A32 A33 A34 A35; ...
      A41 A42 A43 A44 A45; ...
      A51 A52 A53 A54 A55];
  

  
  
% A11
tmp_1 = real(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx3ph.pqac,idx3ph.pqac);
tmp_3 = sum(real(F),2)./abs(E); tmp_3 = tmp_3(idx3ph.pqac);
% tmp_4 = zeros(length(idx3ph.pqac),1); tmp_4(Fl(:,2)) = real(alpha(1));
A11 = tmp_1 + diag(tmp_3);
% A12
tmp_1 = imag(F); tmp_1 = tmp_1(idx3ph.pqac,idx3ph.pqac);
tmp_2 = -sum(imag(F),2); tmp_2 = tmp_2(idx3ph.pqac);
A12 = tmp_1 + diag(tmp_2);
% A13
tmp_1 = real(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx3ph.pqac,idx3ph.vscac);
A13 = tmp_1;

tmp_1 = imag(F); tmp_1 = tmp_1(idx3ph.pqac,idx3ph.vscac);
A14 = tmp_1;

A15 = zeros(size(idx3ph.pqac,1), size(idx3ph.pqdc,1));


% A21
tmp_1 = imag(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx3ph.pqac,idx3ph.pqac);
tmp_3 = sum(imag(F),2)./abs(E); tmp_3 = tmp_3(idx3ph.pqac);
A21 = tmp_1 + diag(tmp_3);
% A22
tmp_1 = -real(F); tmp_1 = tmp_1(idx3ph.pqac,idx3ph.pqac);
tmp_2 = sum(real(F),2); tmp_2 = tmp_2(idx3ph.pqac);
A22 = tmp_1 + diag(tmp_2);
% A23
tmp_1 = imag(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx3ph.pqac,idx3ph.vscac);
A23 = tmp_1;

tmp_1 = -real(F); tmp_1 = tmp_1(idx3ph.pqac,idx1ph.vscac);
A24 = tmp_1;

A25 = zeros(size(idx3ph.pqac,1), size(idx3ph.pqdc,1));


% A31
tmp_1 = real(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx3ph.vscac,idx3ph.pqac);
tmp_3 = sum(real(F),2)./abs(E); tmp_3 = tmp_3(idx3ph.pqac);
% tmp_4 = zeros(length(idx3ph.pqac),1); tmp_4(Fl(:,2)) = real(alpha(1));
A31 = tmp_1 + diag(tmp_3 + a);
% A12
tmp_1 = imag(F); tmp_1 = tmp_1(idx3ph.vscac,idx3ph.pqac);
tmp_2 = -sum(imag(F),2); tmp_2 = tmp_2(idx3ph.pqac);
A32 = tmp_1 + diag(tmp_2 + c);
% A33
tmp_1 = real(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx3ph.vscac,idx3ph.vscac);
tmp_3 = sum(real(F),2)./abs(E); tmp_3 = tmp_3(idx3ph.vscac);
% tmp_4 = real(alpha)*0;
A33 = tmp_1 + diag(tmp_3 + b);
% A34
tmp_1 = imag(F); tmp_1 = tmp_1(idx3ph.vscac,idx3ph.vscac);
tmp_2 = -sum(imag(F),2); tmp_2 = tmp_2(idx3ph.vscac);
% tmp_4 = -imag(beta);
A34 = tmp_1 + diag(tmp_2 + d);
% A35
tmp_1 = real(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx3ph.pqdc,idx3ph.pqdc);
tmp_3 = sum(real(F),2)./abs(E); tmp_3 = tmp_3(idx3ph.vscdc);
A35 = -tmp_1 + -diag(tmp_3); %Added - - 


% A41
tmp_1 = imag(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx3ph.vscac,idx3ph.pqac);
tmp_3 = sum(imag(F),2)./abs(E); tmp_3 = tmp_3(idx3ph.pqac);
A41 = tmp_1 + diag(tmp_3);
% A42
tmp_1 = -real(F); tmp_1 = tmp_1(idx3ph.vscac,idx3ph.pqac);
tmp_2 = sum(real(F),2); tmp_2 = tmp_2(idx3ph.pqac);
A42 = tmp_1 + diag(tmp_2);
% A43
tmp_1 = imag(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx3ph.vscac,idx3ph.vscac);
tmp_3 = sum(imag(F),2)./abs(E); tmp_3 = tmp_3(idx3ph.vscac);
A43 = tmp_1 + diag(tmp_3);
% A44
tmp_1 = -real(Fac); tmp_1 = tmp_1(idx3ph.vscac,idx3ph.vscac);
tmp_2 = sum(real(Fac),2); tmp_2 = tmp_2(idx3ph.vscac);
A44 = tmp_1 + diag(tmp_2);
% A45
A45 = zeros(size(idx3ph.vscac,1), size(idx3ph.pqdc,1));
% A51
A51 = zeros(size(idx3ph.pqdc,1), size(idx3ph.pqac,1));
A52 = zeros(size(idx3ph.pqdc,1), size(idx3ph.pqac,1));
A53 = zeros(size(idx3ph.pqdc,1), size(idx3ph.vscac,1));
A54 = zeros(size(idx3ph.pqdc,1), size(idx3ph.vscdc,1));
% A55
tmp_1 = real(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx3ph.pqdc,idx3ph.pqdc);
tmp_3 = sum(real(F),2)./abs(E); tmp_3 = tmp_3(idx3ph.pqdc);
A55 = tmp_1 + diag(tmp_3); 




A0 = [ A11 A12 A13 A14 A15; ...
      A21 A22 A23 A24 A25; ...
      A31 A32 A33 A34 A35; ...
      A41 A42 A43 A44 A45; ...
      A51 A52 A53 A54 A55];


% M = conj(M_LF(z1,:))';
% Vacg = V_complex_LF(:,z1);
% Vac = V_complex_LF(4:6,z1) + If_complex_LF(4:6,z1).*complex(Rf(4:6),Xf(4:6));%.*[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]'; 
% Iac = If_complex_LF(:,z1);%.*[1 1 1 1]';
% Pac = Nodal_P(z1,:);
% Qac = Nodal_Q(z1,:);
% Vdc = reshape(repmat(Vdc_LF(:,z1)',[3,1]),6,1);
% Idc = -reshape(repmat(Idc_flow(:,z1)/2,[3,1])',6,1); 
% Pdc = Pdc_inj(z1,:);
% 
% M2 = conj(M_LF(z2,:))';
% Vacg2 = V_complex_LF(:,z2);
% Vac2 = V_complex_LF(4:6,z2) + If_complex_LF(4:6,z2).*complex(Rf(4:6),Xf(4:6));%.*[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]'; 
% Iac2 = If_complex_LF(:,z2);%.*[1 1 1 1]';
% Pac2 = Nodal_P(z2,:);
% Qac2 = Nodal_Q(z2,:);
% Vdc2 = reshape(repmat(Vdc_LF(:,z2)',[3,1]),6,1);
% Idc2 = -reshape(repmat(Idc_flow(:,z2)/2,[3,1])',6,1); 
% Pdc2 = Pdc_inj(z2,:);
% 
% 
% A0 = A



M = conj(M_LF(z1,:))';
Vacg = V_complex_LF(:,z1);
% Vac = V_complex_LF(4:6,z1) + If_complex_LF(4:6,z1).*complex(Rf(4:6),Xf(4:6));%.*[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]'; 
Iac = If_complex_LF(:,z1);%.*[1 1 1 1]';
Pac = 3*Nodal_P(z1,:);
Qac = 3*Nodal_Q(z1,:);
Vdc = reshape(repmat(Vdc_LF(:,z1)',[3,1]),6,1);
Idc = -reshape(repmat(Idc_flow(:,z1)/2,[3,1])',6,1); 
Pdc = Pdc_inj(z1,:);

M2 = conj(M_LF(z2,:))';
Vacg2 = V_complex_LF(:,z2);
% Vac2 = V_complex_LF(4:6,z2) + If_complex_LF(4:6,z2).*complex(Rf(4:6),Xf(4:6));%.*[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]'; 
Iac2 = If_complex_LF(:,z2);%.*[1 1 1 1]';
Pac2 = 3*Nodal_P(z2,:);
Qac2 = 3*Nodal_Q(z2,:);
Vdc2 = reshape(repmat(Vdc_LF(:,z2)',[3,1]),6,1);
Idc2 = -reshape(repmat(Idc_flow(:,z2)/2,[3,1])',6,1); 
Pdc2 = Pdc_inj(z2,:);


disp(mode)
switch mode
case 'P5'
    u1 = 0;
    u2 = 0;
    u3 = 0;
    u4 = 0;
    u5 = 1;
    
    
    
    x = linsolve(A0,[u1;u2;u3;u4;u5]);
    
    disp('analytical')
    r1 = Eac(2).*( (1./abs(Eac(2))).*x(1) + 1i*x(2) );
    r2 = Eac(3).*( (1./abs(Eac(3))).*x(3) + 1i*x(4) );
	r3 = x(5);

    disp('EMTP')
    c1 = (Eac(2)-Eac2(2))/(Pdc(2)-Pdc2(2));
    c2 = (Eac(3)-Eac2(3))/(Pdc(2)-Pdc2(2));
    c3 = (Edc(2)-Edc2(2))/(Pdc(2)-Pdc2(2));

    results = table(round([r1; r2; r3],4),round([c1; c2; c3],4),'VariableNames',{'analytical','numerical'},'RowNames',{'dE2/dX','dE3/dX','dE5/dX'})
    
     
case 'E4'
    
    u1 = 0;
    u2 = 0;
    
    tmp_1 = real(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx3ph.vscdc,idx3ph.vscdc);
    tmp_3 = sum(real(F),2)./abs(E); tmp_3 = tmp_3(idx3ph.vscdc);
    u3 = -tmp_1 + -diag(tmp_3); 


%     u3 = -(sum(real(F(idx3ph.vscdc,:)))./abs(E(idx3ph.vscdc)) + real(F(idx3ph.vscdc,idx3ph.vscdc))./abs(E(idx3ph.vscdc)));
    u4 = 0;
    
    tmp_1 = real(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx3ph.pqdc,idx3ph.vscdc);
    u5 = -tmp_1;

%     u5 = -real(F(idx3ph.pqdc,idx3ph.vscdc))./abs(E(idx3ph.vscdc));
    
    x = linsolve(A0,[u1;u2;u3;u4;u5]);
    
    disp('analytical')
    r1 = Eac(2).*( (1./abs(Eac(2))).*x(1) + 1i*x(2) );
    r2 = Eac(3).*( (1./abs(Eac(3))).*x(3) + 1i*x(4) );
	r3 = x(5);

    disp('EMTP')
    c1 = (Eac(2)-Eac2(2))/(Edc(1)-Edc2(1));
    c2 = (Eac(3)-Eac2(3))/(Edc(1)-Edc2(1));
    c3 = (Edc(2)-Edc2(2))/(Edc(1)-Edc2(1));

    results = table(round([r1; r2; r3],4),round([c1; c2; c3],4),'VariableNames',{'analytical','numerical'},'RowNames',{'dE2/dX','dE3/dX','dE5/dX'})

    
case 'Q3'
    u1 = 0;
    u2 = 0;
    u3 = 0;
    u4 = 1;
    u5 = 0;
    
    x = linsolve(A0,[u1;u2;u3;u4;u5]);
    
    disp('analytical')
    r1 = Eac(2).*( (1./abs(Eac(2))).*x(1) + 1i*x(2) );
    r2 = Eac(3).*( (1./abs(Eac(3))).*x(3) + 1i*x(4) );
	r3 = x(5);

    disp('EMTP')
    c1 = (Eac(2)-Eac2(2))/(Qac(1)-Qac2(1));
    c2 = (Eac(3)-Eac2(3))/(Qac(1)-Qac2(1));
    c3 = (Edc(2)-Edc2(2))/(Qac(1)-Qac2(1));

    results = table(round([r1; r2; r3],4),round([c1; c2; c3],4),'VariableNames',{'analytical','numerical'},'RowNames',{'dE2/dX','dE3/dX','dE5/dX'})

case 'P2'
    u1 = 0;
    u2 = 1;
    u3 = 0;
    u4 = 0;
    u5 = 0;
    
    x = linsolve(A0,[u1;u2;u3;u4;u5]);
    
    disp('analytical')
    r1 = Eac(2).*( (1./abs(Eac(2))).*x(1) + 1i*x(2) );
    r2 = Eac(3).*( (1./abs(Eac(3))).*x(3) + 1i*x(4) );
	r3 = x(5);

    disp('EMTP')
    c1 = (Eac(2)-Eac2(2))/(Pac(4)-Pac2(4));
    c2 = (Eac(3)-Eac2(3))/(Pac(4)-Pac2(4));
    c3 = (Edc(2)-Edc2(2))/(Pac(4)-Pac2(4));

    results = table(round([r1; r2; r3],4),round([c1; c2; c3],4),'VariableNames',{'analytical','numerical'},'RowNames',{'dE2/dX','dE3/dX','dE5/dX'})

    
otherwise
    disp('NOPE')
end


idx_slack = 1;
idx1ph.slack = idx_slack;
idx1ph.pqac = [2]'; %[2:26];
idx1ph.pqdc = [5]'; %[2:26];
idx1ph.vscac = [3]';%[];idx_qv;
idx1ph.vscdc = [4]';%[];idx_qv;
idx1ph.qv = []';
idx3ph.slack = idx_slack;
idx3ph.pqac = idx1ph.pqac;
idx3ph.pqdc = idx1ph.pqdc;
idx3ph.vscac = idx1ph.vscac;
idx3ph.vscdc = idx1ph.vscdc;
idx3ph.qv = idx1ph.qv;
idxCtrl = [1:5];
vdep.alpha = repmat([1 0 0],5,1);
vdep.beta = repmat([1 0 0],5,1);
vdep.lambda = repmat([0 0 0],5,1);
vdep.omega = repmat([0 0 0],5,1);
nph = 1;

% [Ydc, YYL, YL, YT, YYT, I_b, Ampacities, y_ih, y_i, A, linedata]  = Ymatrix([1,2,0.1,	1E-09,	1E-09,1,100,0,0],A_b,V_b,ZIN);
% [Yac, YYL, YL, YT, YYT, I_b, Ampacities, y_ih, y_i, A, linedata]  = Ymatrix([1,2,0.00816,0.003581416,0.0000030159,1,100,0,5;2,3,0.00816,0.003581416,0.0000030159,1,100,0,5],A_b,V_b,ZIN);


% S0ac = [0.995+1i*0.8;1+1i*0.8];
% S0dc = [0.995+1i*0.8;1+1i*0.8];
% Eac = Vacg([1,4]);
% Edc = Vdc([1,4]);

S_star = [transpose(complex(Nodal_P(z1,1:3:end), Nodal_Q(z1,1:3:end))); transpose(complex(Pdc_inj(z1,:),0))];
E_star = [Vacg(1:3:end); Vdc(1:3:end)];

S0ac = S_star(1:3);
S0dc = S_star(4:5);
Eac = E_star(1:3);
Edc = E_star(4:5);

Fl = [2,3]; 
% [K, Time] = SC_Voltage_V3(YY_augmented3,S_star,E_star,idx1ph,idx3ph,idxCtrl,nph,vdep);
Eac;
Edc;
filter = 0; unblanced_3ph = 0;
% [K, Time] = SC_Voltage_V5(Yac,S0ac,Eac,Ydc,S0dc,Edc,idx1ph,idx3ph,idxCtrl,nph,vdep,Zf,Fl,unblanced_3ph,filter);
% [K, Time] = SC_Voltage_V5_backup(Yac,S0ac,Eac,Ydc,S0dc,Edc,idx1ph,idx3ph,idxCtrl,nph,vdep,Zf,Fl,unblanced_3ph,filter);
[K, Time] = SC_Voltage_V5(Yac,S0ac,Eac,Ydc,S0dc,Edc,idx1ph,idx3ph,idxCtrl,nph,vdep,Zf,Fl,unblanced_3ph,filter);

n_nodes = 5;
J_PR = zeros(n_nodes);
J_PX = zeros(n_nodes);
J_QR = zeros(n_nodes);
J_QX = zeros(n_nodes);
J_VR = zeros(n_nodes);
J_VX = zeros(n_nodes);
for k = 1:size(K,1)
    if( sum( K{k,1} == idx1ph.pqac ))
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
        J_QR(:,k) = real(K{k,2}{2,1});
        J_QX(:,k) = imag(K{k,2}{2,1});
    elseif( sum( K{k,1} == idx1ph.vscac ))
        J_QR(:,k) = real(K{k,2}{1,1});
        J_QX(:,k) = imag(K{k,2}{1,1});
    elseif( sum( K{k,1} == idx1ph.vscdc ))
        J_VR(:,k) = real(K{k,2}{1,1});
        J_VX(:,k) = imag(K{k,2}{1,1});
        J_QR(:,k) = real(K{k,2}{2,1});
        J_QX(:,k) = imag(K{k,2}{2,1});
    elseif( sum( K{k,1} == idx1ph.pqdc ) )
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
    end
end

J_PR(:,setdiff(1:n_nodes,[idx1ph.pqac; idx1ph.pqdc])) = [];
J_PX(:,setdiff(1:n_nodes,[idx1ph.pqac; idx1ph.pqdc])) = [];
J_QR(:,setdiff(1:n_nodes,[idx1ph.pqac; idx1ph.vscac])) = [];
J_QX(:,setdiff(1:n_nodes,[idx1ph.pqac; idx1ph.vscac])) = [];
J_VR(:,setdiff(1:n_nodes,[idx1ph.vscdc])) = [];
J_VX(:,setdiff(1:n_nodes,[idx1ph.vscdc])) = [];

J_PR(idx1ph.slack,:) = [];
J_PX(idx1ph.slack,:) = [];
J_QR(idx1ph.slack,:) = [];
J_QX(idx1ph.slack,:) = [];
J_VR(idx1ph.slack,:) = [];
J_VX(idx1ph.slack,:) = [];



% E_star = [Vacg(1:3:end); Vdc(1:3:end)];
% E_star2 = [Vacg2(1:3:end); Vdc2(1:3:end)];

dP = [Pac2(1:3:end)-Pac(1:3:end),Pdc2-Pdc]';
dQ = [Qac2(1:3:end)-Qac(1:3:end),0,0]';

dV2 = abs(E_star2).^2-abs(E_star).^2; 
dV = abs(E_star2)-abs(E_star); 

dVr_scs = J_PR * dP(sort([idx1ph.pqac;idx1ph.pqdc])) + J_QR * dQ(sort([idx1ph.pqac;idx1ph.vscac])) + J_VR * dV(sort([idx1ph.vscdc]));
dVi_scs = J_PX * dP(sort([idx1ph.pqac;idx1ph.pqdc])) + J_QX * dQ(sort([idx1ph.pqac;idx1ph.vscac])) + J_VX * dV(sort([idx1ph.vscdc]));

% dVr_scs = J_PR * dP + J_QR * dQ + J_VR * dV;
% dVi_scs = J_PX * dP + J_QX * dQ + J_VX * dV;


%Ground truth
dVr_true = real(E_star2(2:end) - E_star(2:end));
dVi_true = imag(E_star2(2:end) - E_star(2:end));


figure

subplot(2,1,1)
title(strcat('delta V for X = ',mode))
hold on
% scatter(1:25,dVr_jac,'fill')
% scatter(1:25,dVr_sc,'fill')
scatter(1:n_nodes-1,dVr_scs,'fill')
scatter(1:n_nodes-1,dVr_true,'fill')
xticks([1 2 3 4])
xticklabels({'B02','B03','B04','B05'})
ylabel('real part')
legend({'analytical sens coef','ground truth'}) %'analytic senc coef','inverse jacob',
hold off

subplot(2,1,2)
hold on
% scatter(1:25,dVi_jac,'fill')
% scatter(1:25,dVi_sc,'fill')
scatter(1:n_nodes-1,dVi_scs,'fill')
scatter(1:n_nodes-1,dVi_true,'fill')
ylabel('imaginary part')
xticks([1 2 3 4])
xticklabels({'B02','B03','B04','B05'})
hold off
% end

% figure
% subplot(3,1,1); plot(Vdc_LF(2,(10:end))); ylabel('Voltage E3')
% subplot(3,1,2); plot(Nodal_Q((10:end),4)); ylabel('Reactive Power Q2')
% subplot(3,1,3); plot( Pdc_inj((10:end),2)); ylabel('Active Power P4'); xlabel('time')

