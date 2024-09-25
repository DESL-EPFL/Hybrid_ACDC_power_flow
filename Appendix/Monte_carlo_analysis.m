% script for linear power system state estimation
% Monte Carlo analysis

% Willem lambrichts
% General and Unified Model of the Power Flow Problem in Multiterminal AC/DC Networks
% IEEE Transactions on Power Systems 
% 10.1109/TPWRS.2024.3378926

close all
clear all

set(groot,'defaultFigureVisible','on')

set(0, 'DefaultTextInterpreter', 'Latex')
set(0, 'DefaultLegendInterpreter', 'Latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'Latex')


addpath(genpath(pwd))

%% initialize parameters
nb_phases = 1;
[Grid_para,idx,Filter_para,simulation_para] = initialize(nb_phases);


%% Resource limits [complex(min) complex(max)]

% AC - PQ nodes
S3 = [complex(0,-30000), complex(30000,30000)]/Grid_para.A_b;
S5 = [complex(-25000,-25000), complex(25000,25000)]/Grid_para.A_b;
S9 = [complex(-110000,-110000), complex(120000,110000)]/Grid_para.A_b;
S11 = [complex(0,0), complex(30000,0)]/Grid_para.A_b;
S13 = [complex(-5000,0), complex(15000,0)]/Grid_para.A_b;
S14 = [complex(-30000,-30000), complex(30000,30000)]/Grid_para.A_b;

% IC nodes
S15 = [complex(-45000,-45000), complex(45000,45000)]/Grid_para.A_b;
S16 = [complex(-45000,-45000), complex(45000,45000)]/Grid_para.A_b;
S17 = [complex(-45000,-45000), complex(45000,45000)]/Grid_para.A_b;
S18 = [complex(-45000,-45000), complex(45000,45000)]/Grid_para.A_b;

V19 = [760 840]/Grid_para.Vdc_b;
V20 = [760 840]/Grid_para.Vdc_b;
V21 = [760 840]/Grid_para.Vdc_b;
V22 = [760 840]/Grid_para.Vdc_b;

% DC - P nodes
S23 = [complex(-10000,0), complex(15000,0)]/Grid_para.A_b;
S24 = [complex(0,0), complex(5000,0)]/Grid_para.A_b;
S25 = [complex(0,0), complex(5000,0)]/Grid_para.A_b;
S26 = [complex(0,0), complex(5000,0)]/Grid_para.A_b;


%% Sampling

%plant a seed
rng(2,"twister");

%number of scenarios
n = 100;

random_data = rand(32,n);

S_star = [ 
    zeros(Grid_para.n_ph,n);
    zeros(Grid_para.n_ph,n);
    repmat(get_random_profile(random_data([1:2],:),S3),Grid_para.n_ph,1);
    zeros(Grid_para.n_ph,n);
    repmat(get_random_profile(random_data([3:4],:),S5),Grid_para.n_ph,1);
    zeros(Grid_para.n_ph,n);
    zeros(Grid_para.n_ph,n);
    zeros(Grid_para.n_ph,n);
    repmat(get_random_profile(random_data([5:6],:),S9),Grid_para.n_ph,1);
    zeros(Grid_para.n_ph,n);
    repmat(get_random_profile(random_data([7:8],:),S11),Grid_para.n_ph,1);
    zeros(Grid_para.n_ph,n);
    repmat(get_random_profile(random_data([9:10],:),S13),Grid_para.n_ph,1);
    repmat(get_random_profile(random_data([11:12],:),S14),Grid_para.n_ph,1);
    repmat(get_random_profile(random_data([13:14],:),S15),Grid_para.n_ph,1);
    repmat(get_random_profile(random_data([15:16],:),S16),Grid_para.n_ph,1);
    repmat(get_random_profile(random_data([17:18],:),S17),Grid_para.n_ph,1);
    repmat(get_random_profile(random_data([19:20],:),S18),Grid_para.n_ph,1);
    zeros(1,n);
    zeros(1,n);
    zeros(1,n);
    zeros(1,n);
    get_random_profile(random_data([21:22],:),S23);
    get_random_profile(random_data([23:24],:),S24);
    get_random_profile(random_data([25:26],:),S25);
    get_random_profile(random_data([27:28],:),S26);
];

E_star = [
    ones(Grid_para.n_ph,n).*[1]; % exp(-1i*pi*2/3);  exp(1i*pi*2/3) ];
    zeros(Grid_para.n_ph,n);
    zeros(Grid_para.n_ph,n);
    zeros(Grid_para.n_ph,n);
    zeros(Grid_para.n_ph,n);
    zeros(Grid_para.n_ph,n);
    zeros(Grid_para.n_ph,n);
    zeros(Grid_para.n_ph,n);
    zeros(Grid_para.n_ph,n);
    zeros(Grid_para.n_ph,n);
    zeros(Grid_para.n_ph,n);
    zeros(Grid_para.n_ph,n);
    zeros(Grid_para.n_ph,n);
    zeros(Grid_para.n_ph,n);
    zeros(Grid_para.n_ph,n);
    zeros(Grid_para.n_ph,n);
    zeros(Grid_para.n_ph,n);
    zeros(Grid_para.n_ph,n);
    get_random_profile(random_data([21:22],:),V19);
    get_random_profile(random_data([23:24],:),V20);
    get_random_profile(random_data([25:26],:),V21);
    get_random_profile(random_data([27:28],:),V22);
    zeros(1,n);
    zeros(1,n);
    zeros(1,n);
    zeros(1,n);
];




for i = 1:n
%% Solve the PF

t = tic;
[E,J,n_iter] = NR_rectangularACDC_1ph_general(Grid_para,Filter_para,S_star(:,i),E_star(:,i),idx,simulation_para);
time_all(i) = toc(t);

S = E.*conj(Grid_para.YY*E);

E_all(:,i) = E;
S_all(:,i) = S;
n_iter_all(i) = n_iter;

end




%% get color scale
set(groot,'defaultFigureVisible','on')
colorgrd = 'blue_down';
m_cmap = colorgrad(8,colorgrd);
m_cmap2 = [0 0.4470 0.7410 ; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560;0.3010 0.7450 0.9330];

%% save data
folder = './Figures';

%% plot # iter time
f1 = figure('Renderer', 'painters', 'Position', [10 10 1000 300]);

c = cdfplot(n_iter_all);
set( c, 'Linewidth',2.5, 'Color', m_cmap(3,:));

set(0, 'DefaultTextInterpreter', 'Latex')
set(0, 'DefaultLegendInterpreter', 'Latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'Latex')
set(gca,'FontSize',26) 
xlabel('Nb. iterations')
ylabel('Instances')
xlim([0,10])
title(' ')



%% plot cpu time
f2 = figure('Renderer', 'painters', 'Position', [10 10 1000 300]);%

c = cdfplot(time_all);
set( c, 'Linewidth',2.5, 'Color', m_cmap(3,:));

set(0, 'DefaultTextInterpreter', 'Latex')
set(0, 'DefaultLegendInterpreter', 'Latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'Latex')
set(gca,'FontSize',26) 
xlabel('CPU time [sec]')
ylabel('Instances')
% xlim([0,0.01])
title(' ')



%% plot voltage magnitudes
clear gca
Lac = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18];
Ldc = [19,20,21,22,23,24,25,26];

E = abs(E_all)';

Lac3 = [];
for i = 1:length(Lac)
    Lac3 =  [ Lac3; convertCharsToStrings(['$' num2str(Lac(i)) '$'])];
end

Ldc3 = [];
for i = 1:length(Ldc)
    Ldc3 =  [ Ldc3; convertCharsToStrings(['$' num2str(Ldc(i)) '$'])];
end

n_L = [Lac3;Ldc3];

n_name = {'$|$','$V_{AC}$','$|$','$V_{DC}$','$|$'}; 
n_index = [0.5, length(Lac3)/2, length(Lac3)+0.5, length(Lac3)+length(Ldc3)/2, length(Lac3)+length(Ldc3)+0.5];
n_line = [0.75 length(Lac3)+0.25; length(Lac3)+0.75 length(Lac3)+length(Ldc3)+0.25];

f3 = figure('Renderer', 'painters', 'Position', [10 10 1000 500]);
set(gca,'FontSize',30) 
% ylim([min(ylim(gca)),max(ylim(gca))])
% ylim([-1.2e-3,1.2e-3])
    pos = get(gca, 'Position');
    pos(1) = 0.087;
    pos(2) = 0.25;
    pos(3) = 0.9;
    pos(4) = 0.65;
    set(gca, 'Position', pos)
    
aboxplot2(E, 'labels',n_L,'colorgrad','blue_up','Interpreter','Latex');
lowerLevel = min(ylim(gca))-range(ylim(gca))*0.2; 
lowerLabels = n_name; 
text(n_index, repmat(lowerLevel,1,numel(lowerLabels)), lowerLabels,...
    'VerticalAlignment','Top','HorizontalAlignment','Center','Interpreter','Latex', 'FontWeight','bold','FontSize',21)

% gca.Clipping = 'off';
hold on
% yticks([ -1.5e-3 -1e-3 -0.5e-3 0 0.5e-3 1e-3 1.5e-3])



for i = 1:length(n_line)
%     line(n_line(i,:),[lowerLevel,lowerLevel],'Color','k')
    [figx, figy] = axxy2figxy(gca, n_line(i,:), repmat(lowerLevel,1,numel(n_line(i,:))));
    annotation('line',figx,figy)
end
% title('Distribution of the estimation error') 
ylabel('Voltage mag. [p.u.]')


folder = './Figures';
saveas(f3,[folder filesep() 'jpg' filesep() 'Fig7-Monte_Carlo'],'jpg');
saveas(f3,[folder filesep() 'eps' filesep() 'Fig7-Monte_Carlo'],'epsc');




