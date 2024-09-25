set(groot,'defaultFigureVisible','on')

set(0, 'DefaultTextInterpreter', 'Latex')
set(0, 'DefaultLegendInterpreter', 'Latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'Latex')
set(gca,'FontSize',15) 
addpath('/Users/willem/Documents/tools/aboxplot')
%% ESTIMATION

path = '/Users/willem/Documents/phd/State_estimation/matlab_Full_SE_3phase/bad data rejection/validation_estimation_hybrid4.mat';
data =load(path);
clear gca
E_actual = data.data{1};
P_posteriori_diag = data.data{2};
E_expected = randn(size(E_actual)).*sqrt(mean(P_posteriori_diag,2));%sqrt(P_posteriori_diag(:,end));

E = cell(1,2);
E{1} = E_actual';
E{2} = E_expected';

L1 = [1,2,3,4,5,7,8,9,10,11,13,14,15,16,17,18];
L3 = [];

for i = 1:length(L1)
    L3 =  [ L3; convertCharsToStrings(['$' num2str(L1(i)) '_a$']);
                    convertCharsToStrings(['$' num2str(L1(i)) '_b$']);
                    convertCharsToStrings(['$' num2str(L1(i)) '_c$'])];
end


LDC = [convertCharsToStrings(['$' num2str(23) ' \ $']);convertCharsToStrings(['$' num2str(24) ' \ $']);convertCharsToStrings(['$' num2str(25) ' \ $']);convertCharsToStrings(['$' num2str(26) ' \ $'])];
n_L = [L3;L3;LDC];

m_name = {'$|$','Nodal Voltage','$|$','Current injections','$|$','Current flow','$|$','Voltage DC','$|$','Current DC','$|$'}; 
m_index = [0.5, 10.5, 20.5, 31.5, 41.5, 48, 53.5, 58, 62.5, 64, 66];
m_line = [0.75 20.25;20.75 41.25; 41.75 53.25; 53.75 61.25; 61.75 66.75];

n_name = {'$|$','$Nodal \ Voltage_{real}$','$|$','$Nodal \ Voltage_{imag}$','$|$','$I_{DC}$','$|$'}; 
n_index = [0.5, 24, 48.5, 72, 96.5, 98, 100];
n_line = [0.75 48.25;48.75 96.25; 96.75 99.75];

f = figure('Renderer', 'painters', 'Position', [10 10 1500 500])

aboxplot2(E, 'labels',n_L,'colorgrad','blue_down','Interpreter','Latex');
l = legend('Actual','Expected','Interpreter','Latex'); 
% title('Distribution of PECE diag','Interpreter','Latex', 'FontWeight','bold')

lowerLevel = min(ylim(gca))-range(ylim(gca))*.06; 
lowerLabels = n_name; 
text(n_index, repmat(lowerLevel,1,numel(lowerLabels)), lowerLabels,...
    'VerticalAlignment','Top','HorizontalAlignment','Center','Interpreter','Latex', 'FontWeight','bold')
ylim([min(ylim(gca)),max(ylim(gca))])
% gca.Clipping = 'off';
hold on

    pos = get(gca, 'Position');
    pos(1) = 0.055;
    pos(3) = 0.9;
    set(gca, 'Position', pos)

for i = 1:length(n_line)
%     line(n_line(i,:),[lowerLevel,lowerLevel],'Color','k')
    [figx, figy] = axxy2figxy(gca, n_line(i,:), repmat(lowerLevel,1,numel(n_line(i,:))));
    annotation('line',figx,figy)
end
title('Distribution of the estimation error') 

folder = '/Users/willem/Documents/phd/papers'
% saveas(f,[folder filesep() 'estimationerrorvalidation'],'epsc');
% saveas(f,[folder filesep() 'estimationerrorvalidation'],'pdf');


%% PREDICTION

path = '/Users/willem/Documents/phd/State_estimation/matlab_Full_SE_3phase/bad data rejection/validation_prediction_hybrid4.mat';
data =load(path);
clear gca
E_actual = data.data{1};
P_priori_diag = data.data{2};
E_expected = randn(size(E_actual)).*sqrt(mean(P_priori_diag,2));%sqrt(P_priori_diag(:,end));

E = cell(1,2);
E{1} = E_actual';
E{2} = E_expected';

L1 = [1,2,3,4,5,7,8,9,10,11,13,14,15,16,17,18];
L3 = [];

for i = 1:length(L1)
    L3 =  [ L3; convertCharsToStrings(['$' num2str(L1(i)) '_a$']);
                    convertCharsToStrings(['$' num2str(L1(i)) '_b$']);
                    convertCharsToStrings(['$' num2str(L1(i)) '_c$'])];
end


LDC = [convertCharsToStrings(['$' num2str(23) ' \ $']);convertCharsToStrings(['$' num2str(24) ' \ $']);convertCharsToStrings(['$' num2str(25) ' \ $']);convertCharsToStrings(['$' num2str(26) ' \ $'])];
n_L = [L3;L3;LDC];

m_name = {'Nodal Voltage','|','Current injections','|','Current flow','|','Voltage DC','|','Current DC'}; 
m_index = [21, 41.5, 63, 83.5, 96, 107.5, 112, 115.5, 118];

n_name = {'$|$','$Nodal \ Voltage_{real}$','$|$','$Nodal \ Voltage_{imag}$','$|$','$I_{DC}$','$|$'}; 
n_index = [0.5, 24, 48.5, 72, 96.5, 98, 100];
n_line = [0.75 48.25;48.75 96.25; 96.75 99.75];

f2 = figure('Renderer', 'painters', 'Position', [10 10 1500 500])

aboxplot2(E, 'labels',n_L,'colorgrad','blue_down','Interpreter','Latex');
l = legend('Actual','Expected','Interpreter','Latex'); 
% title('Distribution of PECE diag','Interpreter','Latex', 'FontWeight','bold')

lowerLevel = min(ylim(gca))-range(ylim(gca))*.06; 
lowerLabels = n_name; 
text(n_index, repmat(lowerLevel,1,numel(lowerLabels)), lowerLabels,...
    'VerticalAlignment','Top','HorizontalAlignment','Center','Interpreter','Latex', 'FontWeight','bold')
ylim([min(ylim(gca)),max(ylim(gca))])
% gca.Clipping = 'off';
hold on

    pos = get(gca, 'Position');
    pos(1) = 0.055;
    pos(3) = 0.9;
    set(gca, 'Position', pos)

for i = 1:length(n_line)
%     line(n_line(i,:),[lowerLevel,lowerLevel],'Color','k')
    [figx, figy] = axxy2figxy(gca, n_line(i,:), repmat(lowerLevel,1,numel(n_line(i,:))));
    annotation('line',figx,figy)
end
title('Distribution of the prediction error') 

folder = '/Users/willem/Documents/phd/papers';
% saveas(f2,[folder filesep() 'predictionerrorvalidation'],'epsc');
% saveas(f2,[folder filesep() 'predictionerrorvalidation'],'pdf');


%% ESTIMATION reduced
path = '/Users/willem/Documents/phd/State_estimation/matlab_Full_SE_3phase/bad data rejection/validation_estimation_hybrid4.mat';
data =load(path);
clear gca
E_actual = data.data{1};
P_posteriori_diag = data.data{2};
E_expected = randn(size(E_actual)).*sqrt(mean(P_posteriori_diag,2));%sqrt(P_posteriori_diag(:,end));

L1 = [1,2,3,4,5,7,8,9,10,11,13,14,15,16,17,18];
L1 = [2,6,10,14,16,18];


keep = kron_reduce_PMU(polyphase_indices(L1,3),polyphase_indices([6,12],3));
    
ind_keep = [keep(:)',keep(:)'+48,97:100]

E = cell(1,2);
E{1} = E_actual(ind_keep,:)';
E{2} = E_expected(ind_keep,:)';

L3 = [];
for i = 1:length(L1)
    L3 =  [ L3; convertCharsToStrings(['$' num2str(L1(i)) '_a$']);
                    convertCharsToStrings(['$' num2str(L1(i)) '_b$']);
                    convertCharsToStrings(['$' num2str(L1(i)) '_c$'])];
end
LDC = [convertCharsToStrings(['$' num2str(23) ' \ $']);convertCharsToStrings(['$' num2str(24) ' \ $']);convertCharsToStrings(['$' num2str(25) ' \ $']);convertCharsToStrings(['$' num2str(26) ' \ $'])];
n_L = [L3;L3;LDC]; 

m_name = {'$|$','Nodal Voltage','$|$','Current injections','$|$','Current flow','$|$','Voltage DC','$|$','Current DC','$|$'}; 
m_index = [0.5, 10.5, 20.5, 31.5, 41.5, 48, 53.5, 58, 62.5, 64, 66];
m_line = [0.75 20.25;20.75 41.25; 41.75 53.25; 53.75 61.25; 61.75 66.75];

n_name = {'$|$','$Nodal \ Voltage_{real}$','$|$','$Nodal \ Voltage_{imag}$','$|$','$I_{DC}$','$|$'}; 
n_index = [0.5, length(L1)*3/2, length(L1)*3+0.5, length(L1)*3*3/2, length(L1)*3*2+0.5, length(L1)*3*2+2, length(L1)*3*2+4];
n_line = [0.75 length(L1)*3+0.25;length(L1)*3+0.75 length(L1)*3*2+0.25; length(L1)*3*2+0.75 length(L1)*3*2+4-0.25];

f = figure('Renderer', 'painters', 'Position', [10 10 1100 275])
set(gca,'FontSize',15) 
% ylim([min(ylim(gca)),max(ylim(gca))])
ylim([-1.3e-3,1.3e-3]) %ylim([-2.0e-3,2.0e-3])
    pos = get(gca, 'Position');
    pos(1) = 0.055;
    pos(2) = 0.18;
    pos(3) = 0.9;
    pos(4) = 0.750;
    set(gca, 'Position', pos)
aboxplot2(E, 'labels',n_L,'colorgrad','blue_down','Interpreter','Latex');
l = legend('Actual','Expected','Interpreter','Latex'); 
% title('Distribution of PECE diag','Interpreter','Latex', 'FontWeight','bold')
lowerLevel = min(ylim(gca))-range(ylim(gca))*0.135; 
lowerLabels = n_name; 

% gca.Clipping = 'off';
hold on


for i = 1:length(n_line)
%     line(n_line(i,:),[lowerLevel,lowerLevel],'Color','k')
    [figx, figy] = axxy2figxy(gca, n_line(i,:), repmat(lowerLevel,1,numel(n_line(i,:))));
    annotation('line',figx,figy)
end
text(n_index, repmat(lowerLevel,1,numel(lowerLabels)), lowerLabels,...
    'VerticalAlignment','Top','HorizontalAlignment','Center','Interpreter','Latex', 'FontWeight','bold','FontSize',15)
ylim([min(ylim(gca)),max(ylim(gca))])
% title('Distribution of the estimation error') 
ylabel('Estimation error')
% ylim([-0.0022,0.0022])
% folder = '/Users/willem/Documents/phd/papers'
% saveas(f,[folder filesep() 'estimationerrorvalidation2'],'epsc');
% saveas(f,[folder filesep() 'estimationerrorvalidation'],'pdf');

%% PREDICTION reduced
path = '/Users/willem/Documents/phd/State_estimation/matlab_Full_SE_3phase/bad data rejection/validation_prediction_hybrid4.mat';
data =load(path);
clear gca
E_actual = data.data{1};
P_posteriori_diag = data.data{2};
E_expected = randn(size(E_actual)).*sqrt(mean(P_posteriori_diag,2));%sqrt(P_posteriori_diag(:,end));

L1 = [1,2,3,4,5,7,8,9,10,11,13,14,15,16,17,18];
L1 = [2,6,10,14,16,18];



keep = kron_reduce_PMU(polyphase_indices(L1,3),polyphase_indices([6,12],3));
    
ind_keep = [keep(:)',keep(:)'+48,97:100]

E = cell(1,2);
E{1} = E_actual(ind_keep,:)';
E{2} = E_expected(ind_keep,:)';

L3 = [];
for i = 1:length(L1)
    L3 =  [ L3; convertCharsToStrings(['$' num2str(L1(i)) '_a$']);
                    convertCharsToStrings(['$' num2str(L1(i)) '_b$']);
                    convertCharsToStrings(['$' num2str(L1(i)) '_c$'])];
end
LDC = [convertCharsToStrings(['$' num2str(23) ' \ $']);convertCharsToStrings(['$' num2str(24) ' \ $']);convertCharsToStrings(['$' num2str(25) ' \ $']);convertCharsToStrings(['$' num2str(26) ' \ $'])];
n_L = [L3;L3;LDC]; 

m_name = {'$|$','Nodal Voltage','$|$','Current injections','$|$','Current flow','$|$','Voltage DC','$|$','Current DC','$|$'}; 
m_index = [0.5, 10.5, 20.5, 31.5, 41.5, 48, 53.5, 58, 62.5, 64, 66];
m_line = [0.75 20.25;20.75 41.25; 41.75 53.25; 53.75 61.25; 61.75 66.75];

n_name = {'$|$','$Nodal \ Voltage_{real}$','$|$','$Nodal \ Voltage_{imag}$','$|$','$I_{DC}$','$|$'}; 
n_index = [0.5, length(L1)*3/2, length(L1)*3+0.5, length(L1)*3*3/2, length(L1)*3*2+0.5, length(L1)*3*2+2, length(L1)*3*2+4];
n_line = [0.75 length(L1)*3+0.25;length(L1)*3+0.75 length(L1)*3*2+0.25; length(L1)*3*2+0.75 length(L1)*3*2+4-0.25];

f = figure('Renderer', 'painters', 'Position', [10 10 1100 275])
set(gca,'FontSize',15) 
% ylim([min(ylim(gca)),max(ylim(gca))])
ylim([-1.2e-3,1.2e-3])
    pos = get(gca, 'Position');
    pos(1) = 0.055;
    pos(2) = 0.18;
    pos(3) = 0.9;
    pos(4) = 0.750;
    set(gca, 'Position', pos)
    
aboxplot2(E, 'labels',n_L,'colorgrad','blue_down','Interpreter','Latex');
l = legend('Actual','Expected','Interpreter','Latex'); 
% title('Distribution of PECE diag','Interpreter','Latex', 'FontWeight','bold')
% ylim([-0.0022,0.0022])
lowerLevel = min(ylim(gca))-range(ylim(gca))*0.135; 
lowerLabels = n_name; 
text(n_index, repmat(lowerLevel,1,numel(lowerLabels)), lowerLabels,...
    'VerticalAlignment','Top','HorizontalAlignment','Center','Interpreter','Latex', 'FontWeight','bold','FontSize',15)

% gca.Clipping = 'off';
hold on
% yticks([ -1.5e-3 -1e-3 -0.5e-3 0 0.5e-3 1e-3 1.5e-3])



for i = 1:length(n_line)
%     line(n_line(i,:),[lowerLevel,lowerLevel],'Color','k')
    [figx, figy] = axxy2figxy(gca, n_line(i,:), repmat(lowerLevel,1,numel(n_line(i,:))));
    annotation('line',figx,figy)
end
% title('Distribution of the estimation error') 
ylabel('Prediction error')

% folder = '/Users/willem/Documents/phd/papers'
% saveas(f,[folder filesep() 'predictionerrorvalidation2'],'epsc');
% saveas(f,[folder filesep() 'predictionerrorvalidation'],'pdf');