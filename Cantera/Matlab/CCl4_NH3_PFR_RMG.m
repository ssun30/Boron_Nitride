%% CCl4_NH3_PFR_RMG - A Gas Phase Reaction Simulation of the C-N-H-Cl system.
% 
% A plug-flow simulation of the C-N-H-Cl system using gas-phase mechanism
% data generated by RMG. The first part prints out standard enthalpy of reaction for 
% all reactions in the kinetics model. The second part creates concentration vs. 
% temperature plots for CCl4 and NH3 decomposition over 1 s residence time. 
% The third part models the evolution of concentrations of selected species 
% over a period of 0 to 300 ms. 
% 
% Requires: cantera >= 2.6.0
% Keywords: kinetics, thermochemistry.

help BN_PFR_Allendorf1995
LoadCantera;
clear all
close all

% Initial Parameters

c1 = 2.390e-7; % J/kmol to kcal/mol

% Create the gas phase object
fname = '/home/ssun30/Work/Boron_Nitride/YAML_Files/CCl4-NH3_012523.yaml';
phasename = 'gas';

g = Solution(fname, phasename);
g.TP = {300, oneatm};

% Print all the reactions and their dH

for i = 1:g.nReactions
    fprintf('%50s  %13.5g\n', ... 
            g.reactionEqn(i), g.dH_standard(i)*c1);
end

%% Plug-flow prediction of decomposition of BCl3 and NH3 vs. Temperature

t0 = cputime;

T_array1 = linspace(1000, 2000, 21);
P0 = 760 * 133.32;
X0 = 'CCl4(1):0.1,H2(3):0.9';

X_CCL4 = [];
X_CCL3 = [];
X_CCL2 = [];
X_CH4 = [];
X_HCL = [];
X_NH3 = [];

for i = 1:length(T_array1)
    fprintf('Solving for temperature = %d\n', T_array1(i));
    g.TPX = {T_array1(i), P0, X0};
    output = Plug_Flow(g, {'CCL4(1)', 'CCL3(363)', 'CCL2(362)', ...
                           'NH3(2)', 'HCL(293)'}, 0.1);
    X_CCL4 = [X_CCL4, output(3, end)];
    X_CCL3 = [X_CCL3, output(4, end)];
    X_CCL2 = [X_CCL2, output(5, end)];
    X_NH3 = [X_NH3, output(6, end)];
    X_HCL = [X_HCL, output(7, end)];
end

disp(['CPU time = ' num2str(cputime - t0)]);

%% Plot the results
figure(1)
hold on
plot(T_array1, X_CCL4, 'k');
plot(T_array1, X_CCL3, 'r');
plot(T_array1, X_CCL2, 'g');
plot(T_array1, X_NH3, 'bo')
plot(T_array1, X_HCL, 'b');
legend('CCl4', 'CCl3', 'CCl2', 'CH4', 'HCl');
xlabel('Temperature (K)');
ylabel('Mole Fraction');
axis square

% Plug-flow preduction of gas-phase concentrations vs. time.

t0 = cputime;

T_array3 = [1125, 1325];
P0 = 2.0 * 133.32;
X0 = 'CCl4(1):0.1,NH3(2):0.9';
speciesList = {'CCl4(1)', 'NH3(2)', 'HCL(293)'}; 

figure(3)

for i = 1:length(T_array3)
    fprintf('Solving for temperature = %d\n', T_array3(i));
    g.TPX = {T_array3(i), P0, X0};
    output = Plug_Flow(g, speciesList, 0.1);
    tim = output(1, 1:end);
    XX_CCL4 = output(3, 1:end);
    XX_NH3 = output(4, 1:end);
    XX_HCL = output(5, 1:end);
    subplot(1, length(T_array3), i)
    title(sprintf('Temperature = %d K', T_array3(i)))
    hold on
    plot(tim, XX_CCL4, 'r');
    plot(tim, XX_NH3, 'g');
    plot(tim, XX_HCL, 'b');
    legend('CCl4', 'NH3', 'HCl');
    xlabel('Time (s)');
    ylabel('Mole Fraction');
    axis square
end

disp(['CPU time = ' num2str(cputime - t0)]);


