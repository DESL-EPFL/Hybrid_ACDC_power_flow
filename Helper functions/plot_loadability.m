
function f = plot_loadability(S_all,E_all,SAVE,B,complex)

%%
set(groot,'defaultFigureVisible','on')
colorgrd = 'blue_up';
m_cmap = colorgrad(8,colorgrd);
m_cmap2 = [0 0.4470 0.7410 ; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560;0.3010 0.7450 0.9330];


set(groot,'defaultFigureVisible','on')

set(0, 'DefaultTextInterpreter', 'Latex')
set(0, 'DefaultLegendInterpreter', 'Latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'Latex')



f = figure('Renderer', 'painters', 'Position', [10 10 600 400]);

if complex == 'real'
    [ka,ma]  = max(-real(S_all(B,:)));
    breakdown_max = [ka,abs(E_all(B,ma))];
    [ki,mi]  = min(-real(S_all(B,:)));
    breakdown_min = [ki,abs(E_all(B,mi))];
elseif complex == 'imag'
    [ka,ma]  = max(-imag(S_all(B,:)));
    breakdown_max = [ka,abs(E_all(B,ma))];
    [ki,mi]  = min(-imag(S_all(B,:)));
    breakdown_min = [ki,abs(E_all(B,mi))];
end  

% scatter(imag(S_all(55,:)),abs(E_all(55,:)),80,'filled')
hold on
if complex == 'real'
    
    scatter(-real(S_all(B,2:end-1)),abs(E_all(B,2:end-1)),80,'filled','MarkerEdgeColor',m_cmap(6,:),'MarkerFaceColor',m_cmap(6,:));
    s = scatter(breakdown_max(1),breakdown_max(2),380,'x',"LineWidth",4,'MarkerEdgeColor','r');


elseif complex == 'imag'
    scatter(-imag(S_all(B,2:end-1)),abs(E_all(B,2:end-1)),80,'filled','MarkerEdgeColor',m_cmap(6,:),'MarkerFaceColor',m_cmap(6,:));
    s = scatter(breakdown_max(1),breakdown_max(2),380,'x',"LineWidth",4,'MarkerEdgeColor','r');
    scatter(breakdown_min(1),breakdown_min(2),380,'x',"LineWidth",4,'MarkerEdgeColor','r');
    
end

% xlim([-2,7])
% ylim([0.7,1.6])

set(gca,'FontSize',34) 

if complex == 'real'
    xlabel(['Active power at B'+string(B)+' [p.u.]'],'FontSize',34)
    ylabel(['V mag. at B'+string(B)+' [p.u.]'],'FontSize',34)
    
    if B == 16
        xlabel(['Active power at IC2 [p.u.]'],'FontSize',34)
        ylabel(['V mag. at IC2 [p.u.]'],'FontSize',34)
    end
    
    xlim([0,ka*1.2])
    [h,icons] = legend(s,'unsolveable PF','FontSize',32);
    icons = findobj(icons,'Type','patch');
    icons = findobj(icons,'Marker','none','-xor');
    set(icons,'MarkerSize',20,'LineWidth',4);

    
    
elseif complex == 'imag'
    xlabel(['Reactive power at B'+string(B)+' [p.u.]'],'FontSize',34)
    ylabel(['V mag. at B'+string(B)+' [p.u.]'],'FontSize',34)

    if B == 16
        xlabel(['Reactive power at IC2 [p.u.]'],'FontSize',34)
        ylabel(['V mag. at IC2 [p.u.]'],'FontSize',34)
    end   
    
    xlim([ki*1.1,ka*1.1])
    
    [h,icons] = legend(s,'unsolveable PF','FontSize',34 ,'Location','SouthWest');
    icons = findobj(icons,'Type','patch');
    icons = findobj(icons,'Marker','none','-xor');
    set(icons,'MarkerSize',20,'LineWidth',4);


end

title(' ')


pos = get(gca, 'Position');
    pos(1) = 0.16;
    pos(2) = 0.23;
    pos(3) = 0.83;
    
    if complex == 'real'
       pos(4) = 0.72; 
    elseif complex == 'imag'
        pos(4) = 0.74;
    end
    set(gca, 'Position', pos)
    



%%
if SAVE
    
    if complex == 'real'
        name = ['Voltage_collaps_B'+string(B)+'_p'];
    elseif complex == 'imag'
        name = ['Voltage_collaps_B'+string(B)+'_q'];
    end
    
    folder = '/Users/willem/Documents/phd/papers/TPS - PF/PLOTS/';
    saveas(f,strcat(folder, name),'epsc');
end

end




