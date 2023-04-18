% main script for load flow simulation

clear all;
%close all;

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

% ! Define the profile type !
profile_type = 'daily'; % either 'single' or 'daily'

% load
switch(profile_type)
    case 'single'
        P_absorb = importdata('./Data_LF/P_load_single_value.mat');
        Q_absorb = importdata('./Data_LF/Q_load_single_value.mat');
        
        P_inject = importdata('./Data_LF/P_gen_single_value.mat');
        Q_inject = importdata('./Data_LF/Q_gen_single_value.mat');
    case 'daily'
        P_absorb = importdata('./Data_LF/P_daily_load_curve.mat');
        Q_absorb = importdata('./Data_LF/Q_daily_load_curve.mat');
        
        P_inject = importdata('./Data_LF/P_daily_gen_curve.mat');
        Q_inject = importdata('./Data_LF/Q_daily_gen_curve.mat');
    otherwise
        error('unknown profile type');
end

% sanity check
if (~all([size(P_absorb,1),size(Q_absorb,1),size(P_inject,1),size(Q_inject,1)]==n_nodes))
    error('The P/Q profiles are not compatible in terms of number of nodes.');
elseif (length(unique([size(P_absorb,2),size(Q_absorb,2),size(P_inject,2),size(Q_inject,2)]))~=1)
    error('The P/Q profiles are not compatible in terms of number of timesteps.');
end

n_timesteps = size(P_absorb,2);

disp(['The number of timesteps is ' num2str(n_timesteps) '.']);
disp(' ');

% Complex power in p.u. (net)
S_star = complex(P_inject-P_absorb,Q_inject-Q_absorb) / A_b;

%%  Step 3: Simulation parameters
%%% Set the complex load values in p.u. and initialize the bus voltages
%%% Also define tolerance and max number of allowed iterations for the NR

disp('STEP 3: Configuring the load flow simulation.');
disp(' ');

% ! Configure the bus types !

% index of the slack bus
idx_slack = 1;

% indices of the PQ buses
idx_pq = (2:n_nodes)';

% indices of the PV buses
idx_pv = [];

% ! Customize the load flow simulation !

% either 'rectangular' or 'polar'.
coordinate_type = 'rectangular';

% if coordinate_type='daily' -> 'flat' or 'previous',
% if coordinate_type='single' -> 'flat' or 'bad'
start_type = 'flat';

% ! Enter Newton-Raphson algorithm parameters !

% maximum number of iterations
n_max = 100;

% tolerance for the convergence criterion
tol = 1e-6;

%% Step 4: Newton-Raphson algorithm

disp('STEP 4: running Newton-Raphson algorithm.');
disp(' ');

% initialize
E = zeros(n_nodes,n_timesteps);
n_iter = zeros(1,n_timesteps);
t_exec = zeros(1,n_timesteps);

for k = 1:n_timesteps
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
        case 'bad'
            % bad start (only for the case 'single')
            E_0 = ones(n_nodes,1);
            E_0(3) = exp(-1i*pi/4);
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
            [E(:,k),J,n_iter(k)] = NR_rectangular(Y,S_star(:,k),E_0,idx_slack,idx_pq,idx_pv,tol,n_max);
        case 'polar'
            [E(:,k),J,n_iter(k)] = NR_polar(Y,S_star(:,k),E_0,idx_slack,idx_pq,idx_pv,tol,n_max);
    end
    
    t_exec(k) = toc(t_start);
end

%% Step 5: Visualization

% Power profiles at a specific timestep

figure(1);
clf;

axes('FontSize',18);
hold on;
plot((1:n_nodes),P_absorb(:,1),'--ro','LineWidth',2,'MarkerSize',6);
plot((1:n_nodes),Q_absorb(:,1),'--bo','LineWidth',2,'MarkerSize',6);
plot((1:n_nodes),P_inject(:,1),'--go','LineWidth',2,'MarkerSize',6);
plot((1:n_nodes),Q_inject(:,1),'ko','LineWidth',2,'MarkerSize',6);
hold off;

legend('P_{load}','Q_{load}','P_{generation}','Q_{generation}','Location','SouthEast');
title('Nodal Power Profiles at a Specific Timestep');
xlabel('Bus Index');
ylabel('Power');

set(gca,'XTick',1:n_nodes);
axis([1,n_nodes,-max(max(abs(P_absorb(:,1))),max(abs(P_inject(:,1)))),Inf]);

switch(profile_type)
    case 'single'
        % Voltage magnitude profile at a specific timestep
        
        figure(2);
        clf;
        
        axes('FontSize',18);
        plot((1:n_nodes),abs(E),'--bs',...
            'LineWidth',2,...
            'MarkerSize',6,...
            'MarkerEdgeColor','b',...
            'MarkerFaceColor',[0.5,0.5,0.5]);
        title('Voltage Magnitude Profile at a Specific Timestep');
        xlabel('Bus Index');
        ylabel('Voltage Magnitude (p.u.)');
        set(gca,'XTick',1:n_nodes);
        
        % Phase angle profile at a specific timestep
        
        figure(3);
        clf;
        
        axes('FontSize',18);
        plot((1:n_nodes),rad2deg(angle(E)),'--cs',...
            'LineWidth',2,...
            'MarkerSize',6,...
            'MarkerEdgeColor','c',...
            'MarkerFaceColor',[0.5,0.5,0.5]);
        title('Phase Angle Profile at a Specific Timestep');
        xlabel('Bus Index');
        ylabel('Phase Angle (deg)');
        set(gca,'XTick',1:n_nodes);
        
    case 'daily'
        % Profiles over a day
        
        figure(2);
        clf;
        
        subplot(2,2,1);
        plot(real(S_star)');
        title('Active Power (p.u.)');
        grid on;
        
        subplot(2,2,2);
        plot(imag(S_star)');
        title('Reactive Power (p.u.)');
        grid on;
        
        subplot(2,2,3);
        plot(abs(E)');
        title('Voltage Magnitude (p.u.)');
        grid on;
        
        subplot(2,2,4);
        plot(rad2deg(angle(E))');
        title('Phase Angle (deg)');
        grid on;
        
        legend({'Bus 1','Bus 2','Bus 3','Bus 4','Bus 5','Bus 6'});
        
        % Performance
        
        figure(3);
        clf;
        
        subplot(1,2,1);
        plot(n_iter');
        title('Number of Iterations');
        grid on;
        
        subplot(1,2,2);
        plot(t_exec');
        title('Execution Time (s)');
        grid on;
end