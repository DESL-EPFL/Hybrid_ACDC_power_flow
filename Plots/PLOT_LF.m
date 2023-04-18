%% Set
clear all

set(groot,'defaultFigureVisible','on')
set(0, 'DefaultTextInterpreter', 'Latex')
set(0, 'DefaultLegendInterpreter', 'Latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'Latex')

%% Import paths

colorgrd = 'blue_down';
m_cmap1 = colorgrad(3,colorgrd);

colorgrd = 'red_down';
m_cmap2 = colorgrad(3,colorgrd);

%% Get the data



D_balanced = load('data_balanced.mat')
D_unbalanced = load('data_unbalanced.mat')

E_real_b = real(E_delta);
E_imag_b = imag(E_delta);

E_real_u = real(D_unbalanced.data.deltaE);
E_imag_u = imag(D_unbalanced.data.deltaE);

nac = 1:54;
ndc = 55:62;

%% Plot

folder = './Plots/figures';
file_name = fullfile(folder, 'LoadFlow_Error'); 

f1 = figure('Renderer', 'painters', 'Position', [10 10 850 550])
tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'Compact'); 
clf;
set(gca,'FontSize',35)

nexttile
hold on
    set(gca,'FontSize',25)
    histogram(E_real_b(nac),5,'FaceAlpha',0.4,'FaceColor',[0.1 0.1 0.1],'LineWidth',1)
    histogram(E_imag_b(nac),5,'FaceAlpha',0.6,'FaceColor',m_cmap1(2,:),'LineWidth',1)
    legend({'Real','Imaginary'},'Location','NorthWest','FontSize',19); %
    ylabel('Occurence');
    xlabel('$ \Delta V$ [p.u.]');
    
hold off

nexttile
hold on
    set(gca,'FontSize',25)
    histogram(E_real_b(ndc),5,'FaceAlpha',0.6,'FaceColor',m_cmap1(2,:),'LineWidth',1)
    legend({'DC'},'Location','NorthWest','FontSize',19); %
    xlabel('$ \Delta V$ [p.u.]');
    
hold off


nexttile
hold on
    set(gca,'FontSize',25)
    histogram(E_real_u(nac),4,'FaceAlpha',0.4,'FaceColor',[0.1 0.1 0.1],'LineWidth',1)
    histogram(E_imag_u(nac),5,'FaceAlpha',0.6,'FaceColor',m_cmap1(2,:),'LineWidth',1)
    legend({'Real','Imaginary'},'Location','NorthWest','FontSize',19); %
    xlabel('$ \Delta V$ [p.u.]');
    ylabel('Occurence');
hold off

nexttile
hold on
    set(gca,'FontSize',25)
    histogram(E_real_u(ndc),4,'FaceAlpha',0.6,'FaceColor',m_cmap1(2,:),'LineWidth',1)
    legend({'DC'},'Location','NorthWest','FontSize',19); %
    xlabel('$ \Delta V$ [p.u.]');
hold off

% exportgraphics(f1,'/Users/willem/Documents/phd/Optimal_power_flow/figures/LoadFlow_Error.jpg','BackgroundColor','none')
% saveas(f1,[file_name],'depsc');
saveas(f1,[file_name],'jpg');

%% black 'n white

% 
% f2 = figure('Renderer', 'painters', 'Position', [10 10 800 600])
% clf;
% set(gca,'FontSize',35)
% 
% subplot(2,2,1)
% hold on
%     set(gca,'FontSize',25)
%     histogram(E_real_b(nac),5,'FaceAlpha',0.4,'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'LineWidth',1)
%     histogram(E_imag_b(nac),5,'FaceAlpha',0.1,'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'LineWidth',1)
%     legend({'Real','Imaginary'},'Location','NorthWest','FontSize',20); %
%     ylabel('Occurence');
%     xlabel('$ \Delta V$ [p.u.]');
% hold off
% 
% subplot(2,2,2)
% hold on
%     set(gca,'FontSize',25)
%     histogram(E_real_b(ndc),5,'FaceAlpha',0.6,'FaceColor',m_cmap1(2,:))
%     legend({'DC'},'Location','NorthWest','FontSize',20); %
%     xlabel('$ \Delta V$ [p.u.]');
% hold off
% 
% subplot(2,2,3)
% hold on
%     set(gca,'FontSize',25)
%     histogram(E_real_u(nac),4,'FaceAlpha',0.4,'FaceColor',[0.1 0.1 0.1])
%     histogram(E_imag_u(nac),5,'FaceAlpha',0.6,'FaceColor',m_cmap1(2,:))
%     legend({'Real','Imaginary'},'Location','NorthWest','FontSize',20); %
%     xlabel('$ \Delta V$ [p.u.]');
%     ylabel('Occurence');
% hold off
% 
% subplot(2,2,4)
% hold on
%     set(gca,'FontSize',25)
%     histogram(E_real_u(ndc),4,'FaceAlpha',0.6,'FaceColor',m_cmap1(2,:))
%     legend({'DC'},'Location','NorthWest','FontSize',20); %
%     xlabel('$ \Delta V$ [p.u.]');
% hold off
