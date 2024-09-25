function f = make_plot_Fig5(E_delta_balanced,E_delta_unbalanced,Grid_para)

    set(groot,'defaultFigureVisible','on')
    set(0, 'DefaultTextInterpreter', 'Latex')
    set(0, 'DefaultLegendInterpreter', 'Latex')
    set(0, 'DefaultAxesTickLabelInterpreter', 'Latex')

    colorgrd = 'blue_down';
    m_cmap1 = colorgrad(3,colorgrd);

    % Plot

    idxac = 1:Grid_para.n_ac * Grid_para.n_ph;
    idxdc = 1+Grid_para.n_ac * Grid_para.n_ph: Grid_para.n_ac * Grid_para.n_ph + Grid_para.n_dc;



    f = figure('Renderer', 'painters', 'Position', [10 10 850 550]);
    tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'Compact'); 
    clf;
    set(gca,'FontSize',35)

    nexttile
    hold on
        set(gca,'FontSize',25)
        histogram(real(E_delta_balanced(idxac)),5,'FaceAlpha',0.4,'FaceColor',[0.1 0.1 0.1],'LineWidth',1)
        histogram(imag(E_delta_balanced(idxac)),5,'FaceAlpha',0.6,'FaceColor',m_cmap1(2,:),'LineWidth',1)
        legend({'Real','Imaginary'},'Location','NorthWest','FontSize',19); %
        ylabel('Occurence');
        xlabel('$ \Delta V$ [p.u.]');

    hold off

    nexttile
    hold on
        set(gca,'FontSize',25)
        histogram(real(E_delta_balanced(idxdc)),5,'FaceAlpha',0.6,'FaceColor',m_cmap1(2,:),'LineWidth',1)
        legend({'DC'},'Location','NorthWest','FontSize',19); %
        xlabel('$ \Delta V$ [p.u.]');

    hold off


    nexttile
    hold on
        set(gca,'FontSize',25)
        histogram(real(E_delta_unbalanced(idxac)),4,'FaceAlpha',0.4,'FaceColor',[0.1 0.1 0.1],'LineWidth',1)
        histogram(imag(E_delta_unbalanced(idxac)),5,'FaceAlpha',0.6,'FaceColor',m_cmap1(2,:),'LineWidth',1)
        legend({'Real','Imaginary'},'Location','NorthWest','FontSize',19); %
        xlabel('$ \Delta V$ [p.u.]');
        ylabel('Occurence');
    hold off

    nexttile
    hold on
        set(gca,'FontSize',25)
        histogram(real(E_delta_unbalanced(idxdc)),4,'FaceAlpha',0.6,'FaceColor',m_cmap1(2,:),'LineWidth',1)
        legend({'DC'},'Location','NorthWest','FontSize',19); %
        xlabel('$ \Delta V$ [p.u.]');
    hold off

end
