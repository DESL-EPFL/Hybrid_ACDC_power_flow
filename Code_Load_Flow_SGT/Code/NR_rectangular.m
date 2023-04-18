function [E,J,n_iter] = NR_rectangular(Y,S_star,E_0,idx_slack,idx_pq,idx_pv,tol,n_max)
%
% INPUT
% - Y           nodal admittance matrix
% - S_star      given complex power (active and/or reactive power)
% - E_0         initial voltages (phasors)
% - idx_slack   index of the slack bus
% - idx_pq      indices of the PQ buses
% - idx_pv      indices of the PV buses
% - tol         tolerance for convergence criterion
% - n_max       maximum number of iterations
%
% OUTPUT
% - E           solution voltages (phasors)
% - J           Jacobian at the solution
% - n_iter      number of iterations

n_nodes = length(E_0);
G = real(Y);
B = imag(Y);
    
% Initialization
E_re = real(E_0);
E_im = imag(E_0);

for k=1:n_max
    n_iter = k;
    
    % Compute nodal voltages/currents/powers
    E = complex(E_re,E_im);
    I = Y*E;
    S = E.*conj(I);
    
    %% Mismatch calculation
    
    % Compute the mismatches for the entire network.
    dS = S_star-S;
    dP = real(dS);
    dQ = imag(dS);
    dV2 = 1-abs(E).^2; % assumption: V^{*}=1
    
    % Keep only the relevant mismatches.
    dP(idx_slack) = [];
    dQ(sort([idx_pv;idx_slack])) = [];
    dV2(sort([idx_pq;idx_slack])) = [];
    
    dF = [dP;dQ;dV2]; % mismatch of the power flow equations
    
    %% Convergence check
    
    if(max(abs(dF))<tol)
        %disp('NR algorithm has converged to a solution!');
        break;
    elseif(k==n_max)
        disp('NR algorithm reached the maximum number of iterations!');
    end
    
    %%  Jacobian construction
    
    % For the sake of simplicity, the blocks of J are constructed
    % for the whole network (i.e., with size n_nodes x n_nodes).
    % The unnecessary rows/columns are removed subsequently
    
    % Extract real/imaginary part
    E_re = real(E);
    E_im = imag(E);
    
    % Initialization
    J_PR = zeros(n_nodes,n_nodes); % derivative: P versus E_re
    J_PX = zeros(n_nodes,n_nodes); % derivative: P versus E_im
    J_QR = zeros(n_nodes,n_nodes); % derivative: Q versus E_re
    J_QX = zeros(n_nodes,n_nodes); % derivative: Q versus E_im
    J_ER = zeros(n_nodes,n_nodes); % derivative: E^2 versus E_re
    J_EX = zeros(n_nodes,n_nodes); % derivative: E^2 versus E_im
    
    % Construction
    for i=1:n_nodes
        % *****************************************
        % ! write your own code here !
        % *****************************************
              
        % Diagonal elements (terms outside the sum)
        J_PR(i,i) = 2*G(i,i)*E_re(i);
        J_PX(i,i) = 2*G(i,i)*E_im(i);
        J_QR(i,i) = -2*B(i,i)*E_re(i);
        J_QX(i,i) = -2*B(i,i)*E_im(i);
        J_ER(i,i) = 2*E_re(i);
        J_EX(i,i) = 2*E_im(i);
        
        for j=1:n_nodes
            if(j~=i)
                % *****************************************
                % ! write your own code here !
                % *****************************************
                
                % Diagonal elements (terms inside the sum)
                J_PR(i,i) = J_PR(i,i) + G(i,j)*E_re(j) - B(i,j)*E_im(j);
                J_PX(i,i) = J_PX(i,i) + B(i,j)*E_re(j) + G(i,j)*E_im(j);
                J_QR(i,i) = J_QR(i,i) - B(i,j)*E_re(j) - G(i,j)*E_im(j);
                J_QX(i,i) = J_QX(i,i) + G(i,j)*E_re(j) - B(i,j)*E_im(j);
                
                % Offdiagonal elements
                J_PR(i,j) = G(i,j)*E_re(i)+B(i,j)*E_im(i);
                J_PX(i,j) = -B(i,j)*E_re(i)+G(i,j)*E_im(i);
                J_QR(i,j) = -B(i,j)*E_re(i)+G(i,j)*E_im(i);
                J_QX(i,j) = -G(i,j)*E_re(i)-B(i,j)*E_im(i);
            end
        end
    end
    
    % Remove extra rows (i.e., unnecessary equations)
    % slack bus: P & Q & E^2, PV buses: Q, PQ buses: E^2
    
    J_PR(idx_slack,:) = [];
    J_PX(idx_slack,:) = [];
    
    J_QR([idx_pv;idx_slack],:) = [];
    J_QX([idx_pv;idx_slack],:) = [];
    
    J_ER([idx_pq;idx_slack],:)=[];
    J_EX([idx_pq;idx_slack],:)=[];
    
    % Remove extra columns (i.e., variables)
    % slack bus: E_re & E_im
    
    J_PR(:,idx_slack)=[];
    J_QR(:,idx_slack)=[];
    J_ER(:,idx_slack)=[];
    
    J_PX(:,idx_slack)=[];
    J_QX(:,idx_slack)=[];
    J_EX(:,idx_slack)=[];
    
    % Combination
    J = [J_PR,J_PX;J_QR,J_QX;J_ER,J_EX];
    
    %% Solution update
    
    % Solve
    dx = J \ dF;
    
    % Reconstruct the solution
    
    dE_re = zeros(length(E_re),1);
    dE_re(sort([idx_pq;idx_pv]),1) = dx(1:(length(idx_pq)+length(idx_pv)));
    
    dE_im = zeros(length(E_im),1);
    dE_im(sort([idx_pq;idx_pv]),1) = dx((length(idx_pq)+length(idx_pv)+1):end);
    
    % Update
    E_re = E_re + dE_re;
    E_im = E_im + dE_im;
end

E = complex(E_re,E_im);

end