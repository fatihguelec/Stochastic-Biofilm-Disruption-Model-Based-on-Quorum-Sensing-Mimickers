% Stochastic Biofilm Disruption 1-D Model based on Quorum Sensing Mimickers
% This model is created to model the effect of rosmarinic acid (QS
% mimicker) in QS and biofilm disruption given in the paper 'Corral-Lugo, Andrés, et al. "Rosmarinic acid is a 
% homoserine lactone mimic produced by plants that activates a bacterial quorum-sensing regulator." Science Signaling 9.409 (2016)
close all
clear 
clc


tic
%% Definitions
% S1 is the downregulation state
% S2 is the upregulation state
% S3 is the EPS disruption state
% S4 is the biofilm disruption state
% S5 is the no disruption state
% A = Autoinducer AHL molecules
% B = Bacteria
% E = Extracellular Polymeric Substance (EPS)
% M = Quorum sensing mimicker/biofilm disrupter
% S = Nutrient substrate
% ASSUMPTION: Each unit defined below by U represents a group consisting of particles

%% Parameters
% Simulation
MC = 1000; % Number of Monte Carlo loops for Stochastic Tau Leap Simulation
N_av = 6.02214076 * 10^23; % Avogadro constant (mol^-1)
U = 10^-9*N_av; % Number of particles in one unit 
% U = 1;
V = 2*10^-3; % Volume of the domain (l) (Corral-Lugo et al., 2016)
% V = 1*10^-1; % Volume of the domain (l)
delta_t = 0.1; % step size (h)
t_s = 24; % Total simulation time (h)
t_vec = 0:delta_t:t_s; % Time vector (h)

% Initial conditions
init_num = 10; % Initial unit values 
B = init_num; % initial number of B (units)
A = init_num; % initial number of A (units)
E = init_num; % initial number of E (units)
S = init_num*100; % initial number of S (units)
M = 0; % initial number of M (units)
C = 0; % initial number of C (units)

% All states
r_sigma = 1.55*10^-6*N_av*V/U; % Degradation rate of A (h^-1) (Henkel et al., 2013)
mu_max = 0.29; % Maximum specific growth rate of B (h^-1 or g l^-1 h^-1) (Beyenal et al., 2003)
K_M = 26.9*10^-3; % Monod constant (g/l) (Beyenal et al., 2003)
Y_BS = 0.628; % Yield coefficient (-) (Beyenal et al., 2003)
mass_PA = 384*10^-15; % Mass of a P. aeruginosa bacterium (g) (Ghanberi et al.,2016)


C_glu = 0.005; % Glucose concentration of glucose (g/l) (Beyenal et al., 2003 - S_g in Table 2)
k_g = mu_max*(1+Y_BS)/Y_BS;
k_c = k_g/((1+Y_BS)*K_M);
r_g = k_g; % Bacterial growth rate constant (h^-1)
r_c = k_c*C_glu;  % Rate constant of complex formation (h^-1) 
r_dm = 0.001*r_sigma; % Degradation rate of M (h^-1)

% State S1
r_a_1 = 3.8*10^-6*N_av*V/U; % Production rate of A (h^-1) at S1 (Henkel et al., 2013)
mass_eps = 1510*1.66054e-24; % Molecular mass of EPS (g) (Banerjee et al., 2020)
r_e_1 = 0.84/24; % Production rate of E at S1 (Frederick et al., 2011)

% State S2
r_a_2 = 10.9*10^-6*N_av*V/U; % Production rate of A (h^-1) at S2 (Henkel et al., 2013)
r_e_2 = 8.4/24; % Production rate of E at S2 (Frederick et al., 2011)

% State S3
r_e_d = r_e_2; % Disruption rate of EPS (h^-1)

% State S4
r_d = 10*r_g; % Disruption rate of bacteria (h^-1)

% States S1 and S2
r_m = 40*r_a_2; % Production of M (h^-1)

% State decision thresholds
gamma_QS = 50000*10^-9*N_av/U; % Quorum sensing threshold (units l^-1) 5-200 nM in the literature
gamma_DE = 2*10^-3*N_av/U; % EPS disruption threshold (units l^-1) around 2 - 5mM in (Corral-Lugo et al.,2016)
gamma_DB = 7.8*10^-3*N_av/U; % Biofilm disruption threshold (units l^-1) 7.8 - 15.6 mM in (Corral-Lugo et al.,2016)

% Generation of Stoichiometric change matrix
% Columns represent species: A B E M S C
% Rows represent reactions R_1 ... R_7
N_reac = 9; % number of reactions
N_spec = 6; % number of species - A B E M S C
nu = zeros(N_reac, N_spec); % Stoichiometric change matrix (state)
nu(1, 1) = 1; % Reaction 1 - Production of A
nu(2, 3) = 1; % Reaction 2 - Production of E
nu(3, 3) = -1; % Reaction 3 - Disruption of E 
nu(4, 2) = -1; % Reaction 4 - Disruption of B
nu(5, 4) = 1; % Reaction 5 - Production of M
nu(6, 1) = -1; % Reaction 6 - Degradation of A
nu(7, [2 5 6]) = [-1 -1 1]; % Reaction 7 - Production of C with degradation of B and S
nu(8, [2 6]) = [(1+Y_BS) -1]; % Reaction 8 - Production of B with degradation of C
nu(9, 4) = -1; % Reaction 9 - Degradation of M

%Initialize MC variables
C_A_t = zeros(length(t_vec),MC); 
C_B_t = zeros(length(t_vec),MC);
C_E_t = zeros(length(t_vec),MC);
C_S_t = zeros(length(t_vec),MC);
C_M_t = zeros(length(t_vec),MC);
C_C_t = zeros(length(t_vec),MC);
state = zeros(MC,length(t_vec)); % State matrix 

% parpool('local',5) %starts the parallel pool 
% parfor i_mc = 1:MC
for i_mc = 1:MC
    %Initialize variables
    N_A = zeros(1,length(t_vec)); N_A(1) = A; % Number of A
    N_B = zeros(1,length(t_vec)); N_B(1) = B; % Number of B
    N_E = zeros(1,length(t_vec)); N_E(1) = E; % Number of E
    N_M = zeros(1,length(t_vec)); N_M(1) = M; % Number of M
    N_S = zeros(1,length(t_vec)); N_S(1) = S; % Number of S
    N_C = zeros(1,length(t_vec)); N_C(1) = C; % Number of C
    C_A = N_A(1)/V; % concentration of A
    C_M = N_M(1)/V; % concentration of M

    a = zeros(1, N_reac); % Propensity functions matrix
    X = [A B E M S C]; % State vector holding the number of molecules for each species (A B E M S C)    
    state(i_mc,1) = 1; % Assume that initial state is S1

    i = 1;
    t = 0;
    while t <= t_s
        C_A = X(1)/V; % concentration of A
        C_M = X(4)/V; % concentration of M        
        % Decision of the states and updates of propensities depending on state
        % a(1) for Production of A
        % a(2) for Production of E
        % a(3) for Disruption of E
        % a(4) for Disruption of B 
        if C_M < gamma_DE && C_A + C_M < gamma_QS
            state_temp = 1; % State S1
            a(1:4) = [r_a_1*X(2) r_e_1*X(2) 0 0];
        elseif C_M < gamma_DE && C_A + C_M >= gamma_QS
            state_temp = 2; % State S2
            a(1:4) = [r_a_2*X(2) r_e_2*X(2) 0 0];
        elseif C_M >= gamma_DE && C_M < gamma_DB
            state_temp = 3; % State S3
            a(1:4) = [0 0 r_e_d*X(3) 0];            
        else
            state_temp = 4; % State S4
            a(1:4) = [0 0 r_e_d*X(3) r_d*X(2)];
        end

        % Propensity Functions independent of the state
        a(5) = r_m; % Reaction 5 - Production of M
        a(6) = r_sigma*X(1); % Reaction 6 - Degradation of A
        a(7) = r_c*X(2)*X(5); % Reaction 7 - Production of C
        a(8) = r_g*X(6); % Reaction 8 - Production of B
        a(9) = r_dm*X(4); % Reaction 9 - Degradation of M
        
        % Gillespie Direct Algorithm
        [j, tau] = gillespie_direct(a, N_reac);
        X = X + nu(j,:); % Reaction takes place
        X = round(X); % Round X to make all the numbers integer
        t = t + tau; % Update the time
   
        % Sample the number of each species at each time step
            if t >= i*delta_t
                if t >= (i+1)*delta_t % Control if tau is greater than delta_t
                    i_temp = ceil(t/(delta_t));
                    N_A(i+1:i_temp) = X(1);
                    N_B(i+1:i_temp) = X(2);
                    N_E(i+1:i_temp) = X(3);
                    N_M(i+1:i_temp) = X(4);
                    N_S(i+1:i_temp) = X(5);
                    N_C(i+1:i_temp) = X(6);
                    state(i_mc,i+1:i_temp) = state_temp; 
                    i = i_temp;
                else
                    i = i + 1;
                    N_A(i) = X(1);
                    N_B(i) = X(2);
                    N_E(i) = X(3);
                    N_M(i) = X(4);
                    N_S(i) = X(5);
                    N_C(i) = X(6);
                    state(i_mc,i) = state_temp;
                end
            end
    end
        C_A_t(:,i_mc) = N_A(1:length(t_vec))./V; % concentration of A (# l^-1)
        C_B_t(:,i_mc) = N_B(1:length(t_vec))./V; % concentration of B (# l^-1)
        C_E_t(:,i_mc) = N_E(1:length(t_vec))./V; % concentration of E (# l^-1)
        C_M_t(:,i_mc) = N_M(1:length(t_vec))./V; % concentration of M (# l^-1)
        C_S_t(:,i_mc) = N_S(1:length(t_vec))./V; % concentration of S (# l^-1)
        C_C_t(:,i_mc) = N_C(1:length(t_vec))./V; % concentration of C (# l^-1)
%     toc
end

%% Results
t_vec2 = 0:t_s;
t_vec_b = 1:delta_t:t_s+1;


% Bacteria concentration
h_f = figure;
h1 = boxplot(C_B_t(1:10:end,:)',t_vec2,'Whisker', Inf); 
xlabel('Time (h)'); ylabel('Bacterial concentration (CFU/l)');
grid on; 
% ax = gca;
% ax.YAxis.Scale ="log";
yyaxis right
h_s = scatter(t_vec_b, state,10,'filled','m'); ylim([0 5]);
ylabel('States'); set(gca,'ycolor','m');
legend([h1(1,1),h_s], {'Bacteria', 'States'});

% EPS  Concentration
h_f = figure;
h2 = boxplot(C_E_t(1:10:end,:)',t_vec2,'Whisker', Inf); 
xlabel('Time (h)'); ylabel('EPS concentration (particles/l)');
% legend([h1(1,1),h2(1,1)], {'Bacteria', 'EPS'});
grid on;
yyaxis right
h_s = scatter(t_vec_b, state,10,'filled','m'); ylim([0 5]);
ylabel('States'); set(gca,'ycolor','m');
legend([h2(1,1),h_s], {'EPS', 'States'});

% Concentration ratios of the biofilm
C_E_t_mean = mean(C_E_t,2); 
C_B_t_mean = mean(C_B_t,2); 
h_f = figure;
C_bio_t_mean = C_E_t_mean + C_B_t_mean;
C_bio_t_mean_n = C_bio_t_mean./ max(C_bio_t_mean);
C_E_t_mean_n = C_E_t_mean./ max(C_E_t_mean);
C_B_t_mean_n = C_B_t_mean./ max(C_B_t_mean);
plot(t_vec, C_B_t_mean_n, 'b--', 'LineWidth',1.25); hold on;
plot(t_vec, C_E_t_mean_n, 'Color', [0.8500 0.3250 0.0980], 'LineStyle', '-.', 'LineWidth',1.25);
plot(t_vec, C_bio_t_mean_n, 'Color', [0.1 0.6 0.3], 'LineStyle', '-', 'LineWidth',1.25);
xlabel('Time (h)'); ylabel('Normalized concentration (particles/l)');
ylim([-0.1 1.4]);
yyaxis right
scatter(t_vec,state,10,'filled','m'); 
ylabel('States'); set(gca,'ycolor','m');
grid on; legend('Bacteria (Stochastic-mean)','EPS (Stochastic-mean)', 'Biofilm (Stochastic-mean)','States');
xlim([-0.01 24.25]); 
ylim([0 6.5]);

set(h_f,'Units','Inches');
pos = get(h_f,'Position');
set(h_f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(h_f,sprintf('Plot_Biofilm2.pdf'),'-dpdf','-r0') %save as pdf 
toc

% save(['data_MC_', num2str(MC), '.mat']); %save all workspace variables
