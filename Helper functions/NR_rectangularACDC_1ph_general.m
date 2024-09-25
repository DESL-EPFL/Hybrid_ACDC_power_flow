function [E,J,n_iter,time] = NR_rectangularACDC_1ph_general(Grid_para,Filter_para,S_star,E_star,idx,simulation_para)
% INPUT
% - Grid_para   
% - Filter_para
% - S_star      given complex power (active and/or reactive power)
% - E_star      given complex voltage phasors
% - E_0         initial voltages (phasors)
% - idx  	
% - tol         tolerance for convergence criterion
% - n_max       maximum number of iterations
% - pos_ac      array containing the location of the AC side of the AFE and the neighbouring AC node
% - pos_dc      array containing the location of the DC side of the AFE and the neighbouring DC node

% OUTPUT
% - E           solution voltages (phasors)
% - J           Jacobian at the solution
% - n_iter      number of iterations

E_0 = simulation_para.E_0;
tol = simulation_para.tol;
n_max = simulation_para.n_max;

t = tic ;
% Initialization
E_0([idx.slack;idx.vdc;idx.vscdc_vq]) = E_star([idx.slack;idx.vdc;idx.vscdc_vq]);
E_re = real(E_0);
E_im = imag(E_0);
dx = 1;

for k=1:n_max
    
    n_iter = k;
    
    % Compute nodal voltages/currents/powers
    E = complex(E_re,E_im);
    I = complex(Grid_para.G,Grid_para.B)*E;
    S = E.*conj(I);
    
    %% Mismatch calculation
    dF = Mismatch_1ph(E,S,E_star,S_star,Grid_para,Filter_para,idx);
    
    %% Convergence check
    
    if(max(abs(dF))<tol && max(abs(dx)) < tol)
         disp('NR algorithm has converged to a solution!');
        break;
    elseif(k==n_max)
        disp('NR algorithm reached the maximum number of iterations!');
    end
    
    %%  Jacobian construction
    
    % Extract real/imaginary part
    E_re = real(E);
    E_im = imag(E);
        
    n_nodes = Grid_para.n_nodes;
    % Initialization of the Jacobian submatrices
    J_PR = zeros(n_nodes,n_nodes); % derivative: P* versus E_re
    J_PX = zeros(n_nodes,n_nodes); % derivative: P* versus E_im
    J_QR = zeros(n_nodes,n_nodes); % derivative: Q* versus E_re
    J_QX = zeros(n_nodes,n_nodes); % derivative: Q* versus E_im
    J_ER = zeros(n_nodes,n_nodes); % derivative: E^2 versus E_re
    J_EX = zeros(n_nodes,n_nodes); % derivative: E^2 versus E_im
        
    J_AFErR = zeros(n_nodes,n_nodes); % derivative: Pdc* versus E_re
    J_AFErX = zeros(n_nodes,n_nodes); % derivative: Pdc* versus E_im
    J_AFEiR = zeros(n_nodes,n_nodes); % derivative: Edc* versus E_re
    J_AFEiX = zeros(n_nodes,n_nodes); % derivative: Edc* versus E_im
    
    J_PpR = zeros(n_nodes,n_nodes); % derivative: P+* versus E_re
    J_PpX = zeros(n_nodes,n_nodes); % derivative: P+* versus E_im
    J_QpR = zeros(n_nodes,n_nodes); % derivative: Q+* versus E_re
    J_QpX = zeros(n_nodes,n_nodes); % derivative: Q+* versus E_im
    
    J_EdiorR = zeros(n_nodes,n_nodes); % derivative: E0r versus E_re
    J_EdiorX = zeros(n_nodes,n_nodes); % derivative: E0r versus E_im
    J_EdioiR = zeros(n_nodes,n_nodes); % derivative: E0i versus E_re
    J_EdioiX = zeros(n_nodes,n_nodes); % derivative: E0i versus E_im
    
    
    % PQ, PV, Pdc     
    [J_PR, J_PX, J_QR, J_QX, J_ER, J_EX] = Jacobian_Powers_phase(E_re,E_im,Grid_para, J_PR, J_PX, J_QR, J_QX, J_ER, J_EX);

    % E0, E+, E-
%     [J_EdiorR, J_EdiorX, J_EdioiR, J_EdioiX] = Jacobian_Voltage_symmetric(Grid_para, Filter_para, J_EdiorR, J_EdiorX, J_EdioiR, J_EdioiX);
    
    % P+,Q+
%     [J_PpR, J_PpX, J_QpR, J_QpX] = Jacobian_Powers_symmetric(E_re, E_im, Grid_para, J_PpR, J_PpX, J_QpR, J_QpX);

    % Interfacing converters/ Active Front Ends   
    [J_AFEiR, J_AFEiX, J_AFErR, J_AFErX] = Jacobian_Converters_1ph(E_re,E_im, E_star, Grid_para, Filter_para, J_EdiorR, J_EdiorX, J_EdioiR, J_EdioiX, J_PR, J_PX, J_QR, J_QX, J_AFEiR, J_AFEiX, J_AFErR, J_AFErX);

  
%     % Remove extra rows (i.e., unnecessary equations)
%     J_PR(sort([idx.slack;idx.vdc;idx.vscac_pq;idx.vscac_vq ;idx.vscdc_pq;idx.vscdc_vq]),:)=[];
%     J_PX(sort([idx.slack;idx.vdc;idx.vscac_pq;idx.vscac_vq ;idx.vscdc_pq;idx.vscdc_vq]),:)=[];
%     
%     J_QR(sort([idx.slack;idx.pvac;idx.pdc;idx.vdc;idx.vscac_pq;idx.vscac_vq ;idx.vscdc_pq;idx.vscdc_vq]),:)=[];
%     J_QX(sort([idx.slack;idx.pvac;idx.pdc;idx.vdc;idx.vscac_pq;idx.vscac_vq ;idx.vscdc_pq;idx.vscdc_vq]),:)=[];
%     
%     J_ER(sort([idx.slack;idx.pqac;idx.pdc;idx.vdc;idx.vscac_pq;idx.vscac_vq ;idx.vscdc_pq;idx.vscdc_vq]),:)=[];
%     J_EX(sort([idx.slack;idx.pqac;idx.pdc;idx.vdc;idx.vscac_pq;idx.vscac_vq ;idx.vscdc_pq;idx.vscdc_vq]),:)=[];
%     
%     J_AFErR(sort([idx.slack;idx.pqac;idx.pvac;idx.pdc;idx.vdc;idx.vscdc_pq;idx.vscac_vq]),:)=[];%dc
%     J_AFErX(sort([idx.slack;idx.pqac;idx.pvac;idx.pdc;idx.vdc;idx.vscdc_pq;idx.vscac_vq]),:)=[];%dc
%     J_AFEiR(sort([idx.slack;idx.pqac;idx.pvac;idx.pdc;idx.vdc;idx.vscdc_pq;idx.vscdc_vq]),:)=[];%ac
%     J_AFEiX(sort([idx.slack;idx.pqac;idx.pvac;idx.pdc;idx.vdc;idx.vscdc_pq;idx.vscdc_vq]),:)=[];%ac

    % Remove extra Columns (i.e., unnecessary unknown variables)
    % Combination

    J_PR(:,sort([idx.slack;idx.vdc;idx.vscdc_vq]))=[];
    J_PX(:,sort([idx.slack;idx.vdc;idx.pdc;idx.vscdc_pq;idx.vscdc_vq]))=[];
    
    J_QR(:,sort([idx.slack;idx.vdc;idx.vscdc_vq]))=[];
    J_QX(:,sort([idx.slack;idx.vdc;idx.pdc;idx.vscdc_pq;idx.vscdc_vq]))=[];

    J_ER(:,sort([idx.slack;idx.vdc;idx.vscdc_vq]))=[];
    J_EX(:,sort([idx.slack;idx.vdc;idx.pdc;idx.vscdc_pq;idx.vscdc_vq]))=[];
    
    J_AFErR(:,sort([idx.slack;idx.vdc;idx.vscdc_vq]))=[];
    J_AFErX(:,sort([idx.slack;idx.vdc;idx.pdc;idx.vscdc_pq;idx.vscdc_vq]))=[];
    J_AFEiR(:,sort([idx.slack;idx.vdc;idx.vscdc_vq]))=[];
    J_AFEiX(:,sort([idx.slack;idx.vdc;idx.pdc;idx.vscdc_pq;idx.vscdc_vq]))=[];

    
     
     J =[J_PR(idx.pqac,:),J_PX(idx.pqac,:); %PQac real
         J_QR(idx.pqac,:),J_QX(idx.pqac,:); %PQac imag
         J_PR(idx.pvac,:),J_PX(idx.pvac,:); %PVac real
         J_ER(idx.pvac,:),J_EX(idx.pvac,:); %PVac imag
         J_AFErR(idx.vscac_pq,:),J_AFErX(idx.vscac_pq,:);
         J_AFEiR(idx.vscac_pq,:),J_AFEiX(idx.vscac_pq,:); % VSC AC P-Q
         J_AFErR(idx.vscac_vq,:),J_AFErX(idx.vscac_vq,:); % VSC DC Vdc-Q
         J_AFEiR(idx.vscac_vq,:),J_AFEiX(idx.vscac_vq,:); % VSC AC Vdc-Q
         J_PR(idx.vscdc_pq,:),J_PX(idx.vscdc_pq,:); %Pdc
         J_PR(idx.pdc,:),J_PX(idx.pdc,:)];
         %J_ER(idx.vscdc_vq,:),J_EX(idx.vscdc_vq,:)]; %Pdc
     
     
    %% Solution update
    
    % Solve
    dx = J \ dF;
    
    % Reconstruct the solution
    dE_re = zeros(length(E_re),1);
    dE_re(sort([idx.pqac;idx.pvac;idx.pdc;idx.vscac_pq;idx.vscac_vq;idx.vscdc_pq]),1)...
        = dx(1:length(idx.pqac) + length(idx.pvac) + length(idx.pdc) +...
                 length(idx.vscac_pq) + length(idx.vscac_vq) + length(idx.vscdc_pq));
    
    dE_im = zeros(length(E_im),1);
    dE_im(sort([idx.pqac;idx.pvac;idx.vscac_vq;idx.vscac_pq]),1)...
        = dx( length(idx.pqac) + length(idx.pvac) + length(idx.pdc) +...
                 length(idx.vscac_pq) + length(idx.vscac_vq) + length(idx.vscdc_pq) + 1: end);
    
    % Update
    E_re = E_re + dE_re;
    E_im = E_im + dE_im;
end

time = toc(t);
E = complex(E_re,E_im);
end


