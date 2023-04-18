% main script for load flow simulation

clear all;
close all;

% Set Interpreter
set(0,'defaultTextInterpreter','latex');

disp(' ');
disp('********************************************************');
disp('* This script solves the load flow problem *');
disp('********************************************************');
disp(' ');

%%  Step 1: load Y matrix & base values
%   The nodal admittance matrix should saved in a folder 'Data_LF'
%   located in the same directory as this script.

disp('STEP 1: loading nodal admittance matrix and base values ...');
disp(' ');

% ! put the name of the file containing the Y matrix !
Y = importdata ('./Data_LF/Y.mat');

% ! put the names of the files containing the base values !
A_b = importdata ('./Data_LF/Base_Power.mat');
V_b = importdata ('./Data_LF/Base_Voltage.mat');

n_nodes = size(Y,1); % number of nodes

disp(['The network consists of ',num2str(n_nodes), ' nodes.']);
disp(['The base power is ' num2str(A_b/1e6), 'MVA.']);
disp(['The base voltage is ' num2str(V_b/1e3) 'kV.']);
disp(' ');

%%  Step 2: load power profiles of loads and generators
% The nodal power profiles should be also be saved 'Data_LF',
% as matrices of size n_nodes x n_timesteps.

disp('STEP 2: Loading power profiles.')
disp(' ');

% load
P_absorb = importdata('./Data_LF/P_daily_load_curve.mat');
Q_absorb = importdata('./Data_LF/Q_daily_load_curve.mat');

P_inject = importdata('./Data_LF/P_daily_gen_curve.mat');
Q_inject = importdata('./Data_LF/Q_daily_gen_curve.mat');

% sanity check
if (~all([size(P_absorb,1),size(Q_absorb,1),size(P_inject,1),size(Q_inject,1)]==n_nodes))
    error('The P/Q profiles are not compatible in terms of number of nodes.');
elseif (length(unique([size(P_absorb,2),size(Q_absorb,2),size(P_inject,2),size(Q_inject,2)]))~=1)
    error('The P/Q profiles are not compatible in terms of number of timesteps.');
end

% Complex power in p.u. (net)
S_star = complex(P_inject-P_absorb,Q_inject-Q_absorb) / A_b;
S_star = S_star(:,3); % extract a specific timestep

%%  Step 3: Simulation parameters
%%% Set the complex load values in p.u. and initialize the bus voltages
%%% Also define tolerance and max number of allowed iterations for the NR

disp('STEP 3: Configuring the load flow simulation.');
disp(' ');

% ! Configure the bus types !
idx_slack = 1; % index of the slack bus
idx_pq = (2:n_nodes)'; % indices of the PQ buses
idx_pv = []; % indices of the PV buses

% ! Customize the load flow simulation !
coordinate_type = 'polar'; % either 'rectangular' or 'polar'.
start_type = 'previous'; % either 'flat', 'previous'!

% ! Enter Newton-Raphson algorithm parameters !
n_max = 100; % maximum number of iterations
tol = 1e-6; % tolerance for the convergence criterion

%% Step 4: Newton-Raphson algorithm

disp('STEP 4: running Newton-Raphson algorithm.');
disp(' ');

% ! Configure here !

% index of the bus whose net power is to be scaled
idx_scale = 3;

% scaling factor (vector) to be applied to the said bus
%lambda = [linspace(1,500,10),linspace(500,580,80)];
lambda = 1:580;
n_steps = length(lambda);

% scale the complex power of the bus idx_scale
S_scaled = repmat(S_star,[1,n_steps]);
S_scaled(idx_scale,:) = lambda .* S_scaled(idx_scale,:);

% initialize
E = zeros(n_nodes,n_steps);
n_iter = zeros(1,n_steps);
t_exec = zeros(1,n_steps);

for k = 1:n_steps
    %% Initialization   
    switch(start_type)
        case 'flat'
            % always use a flat start
            E_0 = ones(n_nodes,1);
        case 'previous'
            % use the result of a previous timestep (if available)
            if(k<=1)
                E_0 = ones(n_nodes,1);
            else
                E_0 = E(:,k-1);
            end
        otherwise
            error('unknown start type');
    end
    
    %% Call Newton-Raphson function
    
    % The functions NR_rectangular/NR_polar(.) internally make use of
    % rectangular/polar coordinates for the calculation.
    % However, they accept and return complex values (i.e., phasors).

    t_start = tic;
    
    switch(coordinate_type)
        case 'rectangular'
            %E_0 = complex(E_re_0,E_im_0);
            [E(:,k),J,n_iter(k)] = NR_rectangular(Y,S_scaled(:,k),E_0,idx_slack,idx_pq,idx_pv,tol,n_max);
        case 'polar'
            %E_0 = E_abs_0 .* exp(1i*E_arg_0);
            [E(:,k),J,n_iter(k)] = NR_polar(Y,S_scaled(:,k),E_0,idx_slack,idx_pq,idx_pv,tol,n_max);
    end
    
    t_exec(k) = toc(t_start);
end

%% Step 5: Visualization

% Results

E_abs = abs(E);
E_arg = rad2deg(angle(E));
P = real(S_scaled);
Q = imag(S_scaled);
idx_problem = find(n_iter == n_max,1) - 1;

disp([' ']);
disp(['Results']);
disp(['Bus #' num2str(idx_scale) ' has been scaled with a factor lambda = 1:' num2str(lambda(end))]);
disp(['Starting at lambda=' num2str(lambda(idx_problem+1)) ', the NR algorithm stopped converging']);

% Evolution of Voltage Magnitude / Phase Angle

figure(1);
clf;

subplot(2,2,1);
plot(lambda(1:idx_problem),abs(E(:,1:idx_problem))); 
title('Voltage Magnitude Evolution ($n_{iter}<n_{max}$)');
xlabel('$\lambda$ (unitless)');
ylabel('$|V|$ (p.u.)');
grid on;

subplot(2,2,2)
plot(lambda(idx_problem:end),abs(E(:,idx_problem:end)));
title('Voltage Magnitude Evolution ($n_{iter} \geq n_{max}$)');
xlabel('$\lambda$ (unitless)');
ylabel('$|V|$ (p.u.)');
grid on;

subplot(2,2,3)
plot(lambda(1:idx_problem),E_arg(:,1:idx_problem));
title('Phase Angle Evolution ($n_{iter}<n_{max}$)');
xlabel('$\lambda$ (unitless)');
ylabel('Angle (deg)');
grid on;

subplot(2,2,4)
plot(lambda(idx_problem:end),E_arg(:,idx_problem:end));
title('Phase Angle Evolution ($n_{iter} \geq n_{max}$)');
xlabel('$\lambda$ (unitless)');
ylabel('Angle (deg)');
grid on;

subplot(2,2,1);
legend({'Bus 1','Bus 2','Bus 3','Bus 4','Bus 5','Bus 6'},'Location','SouthWest');

% Performance

figure(2);
clf;

subplot(1,2,1)
plot(lambda,n_iter);
title('Number of Iterations');
xlabel('$\lambda$ (unitless)');
ylabel('No. Iterations');
grid on;

subplot(1,2,2)
plot(lambda,t_exec);
title('Execution Time');
xlabel('$\lambda$ (unitless)');
ylabel('Time (sec)');
grid on;

% Nose Curves

figure(3);
clf;

subplot(2,2,1);
plot(abs(P(idx_scale,1:idx_problem))', E_abs(idx_scale,1:idx_problem)'); % ...
title(['PV Nose Curve of Bus ' int2str(idx_scale) ' ($n_{iter}<n_{max}$)']);
xlabel('Load Power -P (p.u.)');
ylabel('Voltage Magnitude $|V|$ (p.u.)');
grid on;

subplot(2,2,2);
plot(abs(P(idx_scale,idx_problem:end))', E_abs(idx_scale,idx_problem:end)'); % ...
title(['PV Nose Curve of Bus ' int2str(idx_scale) ' ($n_{iter} \geq n_{max}$)']);
xlabel('Load Power -P (p.u.)');
ylabel('Voltage Magnitude $|V|$ (p.u.)');
grid on;

subplot(2,2,3);
plot(abs(Q(idx_scale,1:idx_problem))', E_abs(idx_scale,1:idx_problem)'); % ...
title(['QV Nose Curve of Bus ' int2str(idx_scale) ' ($n_{iter}<n_{max}$)']);
xlabel('Load Power -Q (p.u.)');
ylabel('Voltage Magnitude $|V|$ (p.u.)');
grid on;

subplot(2,2,4);
plot(abs(Q(idx_scale,idx_problem:end))', E_abs(idx_scale,idx_problem:end)'); % ...
title(['QV Nose Curve of Bus ' int2str(idx_scale) ' ($n_{iter} \geq n_{max}$)']);
xlabel('Load Power -Q (p.u.)');
ylabel('Voltage Magnitude $|V|$ (p.u.)');
grid on;