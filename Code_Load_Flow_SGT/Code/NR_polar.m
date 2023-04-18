function [E,J,n_iter] = NR_polar(Y,S_star,E_0,idx_slack,idx_pq,idx_pv,tol,n_max)
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

Y_abs = abs(Y);
Y_arg = angle(Y);
n_nodes = length(E_0);

% Initialization
E_abs = abs(E_0);
E_arg = angle(E_0);

for k=1:n_max
    n_iter = k;
    
    % Compute nodal voltages/currents/power
    E = complex(E_abs.*cos(E_arg),E_abs.*sin(E_arg));
    I = Y*E;
    S = E.*conj(I);
    
    %% Mismatch calculation
    
    % Compute the mismatches for the entire network.
    dS = S_star-S;
    dP = real(dS);
    dQ = imag(dS);
    
    % Keep only the relevant mismatches.
    dP(idx_slack) = [];
    dQ(sort([idx_pv;idx_slack])) = [];
    
    dF = [dP;dQ]; % mismatch of the power flow equations
    
    %% Convergence check
    
    if(max(abs(dF))<tol)
        %disp('NR algorithm has converged to a solution!');
        break;
    elseif(k==n_max)
        disp('NR algorithm reached the maximum number of iterations!');
    end
    
    %% Jacobian construction
    
    % For the sake of simplicity, the blocks of J are constructed
    % for the whole network (i.e., with size n_nodes x n_nodes).
    % The unnecessary rows/columns are removed subsequently
    
    % Extract magnitude/angle
    E_abs = abs(E);
    E_arg = angle(E);
    
    % Initialization
    J_PE = zeros(n_nodes,n_nodes); % derivative: P versus E_abs
    J_PT = zeros(n_nodes,n_nodes); % derivative: P versus E_arg (theta)
    J_QE = zeros(n_nodes,n_nodes); % derivative: Q versus E_abs
    J_QT = zeros(n_nodes,n_nodes); % derivative: Q versus E_arg (theta)
    
    % Construction
    for i=1:n_nodes
        % *****************************************
        % ! write your own code here !
        % *****************************************
        
        % Diagonal elements (terms outside the sum)
        J_PE(i,i) = 2*Y_abs(i,i)*E_abs(i)*cos(Y_arg(i,i));
        J_QE(i,i) = -2*Y_abs(i,i)*E_abs(i)*sin(Y_arg(i,i));
        
        for j=1:n_nodes
            if i ~= j
                % *****************************************
                % ! write your own code here !
                % *****************************************
                
                % Diagonal elements (terms inside the sum)
                J_PE(i,i) = J_PE(i,i) + Y_abs(i,j)*E_abs(j)*cos(E_arg(i)-E_arg(j)-Y_arg(i,j));
                J_QE(i,i) = J_QE(i,i) + Y_abs(i,j)*E_abs(j)*sin(E_arg(i)-E_arg(j)-Y_arg(i,j));
                J_PT(i,i) = J_PT(i,i) - E_abs(i)*Y_abs(i,j)*E_abs(j)*sin(E_arg(i)-E_arg(j)-Y_arg(i,j));
                J_QT(i,i) = J_QT(i,i) + E_abs(i)*Y_abs(i,j)*E_abs(j)*cos(E_arg(i)-E_arg(j)-Y_arg(i,j));

                % Offdiagonal elements
                J_PE(i,j) = Y_abs(i,j)*E_abs(i)*cos(E_arg(i)-E_arg(j)-Y_arg(i,j));
                J_QE(i,j) = Y_abs(i,j)*E_abs(i)*sin(E_arg(i)-E_arg(j)-Y_arg(i,j));
                J_PT(i,j) = Y_abs(i,j)*E_abs(i)*E_abs(j)*sin(E_arg(i)-E_arg(j)-Y_arg(i,j));
                J_QT(i,j) = -Y_abs(i,j)*E_abs(i)*E_abs(j)*cos(E_arg(i)-E_arg(j)-Y_arg(i,j));
            end
        end
    end
    
    % Remove extra rows (i.e., unnecessary equations)
    % slack bus: P & Q, PV buses: Q
    
    J_PE(idx_slack,:) = [];
    J_PT(idx_slack,:) = [];
    
    J_QE([idx_pv;idx_slack],:) = [];
    J_QT([idx_pv;idx_slack],:) = [];
    
    % Remove extra columns (i.e., variables)
    % slack bus: E_abs & E_arg, PV nodes: E_abs
    
    J_PE(:,[idx_pv,idx_slack]) = [];
    J_QE(:,[idx_pv,idx_slack]) = [];
    
    J_PT(:,idx_slack) = [];
    J_QT(:,idx_slack) = [];
    
    % Combination
    J = [J_PE,J_PT;J_QE,J_QT];
    
    %% Solution update
    
    % Solve
    dx = J \ dF;
    
    % Reconstruct the solution
    
    dE_abs = zeros(length(E_abs),1);
    dE_abs(idx_pq,1) = dx(1:length(idx_pq));
    
    dE_arg = zeros(length(E_arg),1);
    dE_arg(sort([idx_pq;idx_pv]),1) = dx((length(idx_pq)+1):end);
    
    % Update
    E_abs = E_abs + dE_abs;
    E_arg = E_arg + dE_arg;
end

E = E_abs .* exp(1i*E_arg);

end