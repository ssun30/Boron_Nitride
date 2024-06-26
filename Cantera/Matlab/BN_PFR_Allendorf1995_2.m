%% BN_PFR_Allendorf1995_2 - A Gas Phase Reaction Simulation of the B-N-H-Cl system.
% 
% A full plug-flow reactor simulation of the B-N-H-Cl system using gas-phase mechanism
% data from M. Allendorf et al. 1995 and Thermochemistry data from M. Allendorf
% et al. 1997. The first part prints out standard enthalpy of reaction for 
% all reactions in the kinetics model. The second part creates concentration vs. 
% temperature plots for BCl3 and NH3 decomposition over 1 s residence time. 
% The third part models the evolution of concentrations of selected species 
% over a period of 0 to 300 ms. 
% 
% Requires: cantera >= 2.6.0
% Keywords: kinetics, thermochemistry, boron nitride. 

help BN_PFR_Allendorf1995_2
LoadCantera;
clear all
close all

% Initial Parameters

c1 = 2.390e-7; % J/kmol to kcal/mol
dx = 0.001; % Step size 

% Create the gas phase object

fname = '/home/ssun30/Work/Boron_Nitride/YAML_Files/BNHCL_Allendorf_NASA7.yaml';
phasename = 'gas';

g = Solution(fname, phasename);

g.TP = {298.15, oneatm};

% Print all the reactions and their dH

for i = 1:g.nReactions
    fprintf('%50s  %13.5g\n', ... 
            g.reactionEqn(i), g.dH_standard(i)*c1);
end

% Plug-flow prediction of decomposition of BCl3 and NH3 vs. Temperature

% t0 = cputime;
% 
% T_array1 = linspace(1000, 2000, 21);
% P0 = 760 * 133.32;
% X0 = 'BCl3:0.05,H2:0.95';
% 
% X_BCL3 = [];
% X_BCL2 = [];
% X_BCL = [];
% X_HCL = [];
% X_NH3 = [];
% X_N2 = [];
% X_H2 = [];
% 
% for i = 1:length(T_array1)
%     fprintf('Solving for temperature = %d\n', T_array1(i));
%     g.TPX = {T_array1(i), P0, X0};
%     output = Plug_Flow_Reactor(g, {'BCl3', 'BCl2', 'BCl', 'HCl'}, 1.0, dx);
%     X_BCL3 = [X_BCL3, output(3, end)];
%     X_BCL2 = [X_BCL2, output(4, end)];
%     X_BCL = [X_BCL, output(5, end)];
%     X_HCL = [X_HCL, output(6, end)];
% end
% 
% T_array2 = linspace(800, 2000, 25);
% X0 = zeros(1, g.nSpecies);
% X0(g.speciesIndex('NH3')) = 1.00;
% 
% X_NH3 = [];
% 
% for i = 1:length(T_array2)
%     fprintf('Solving for temperature = %d\n', T_array2(i));
%     g.TPX = {T_array2(i), P0, X0};
%     output = Plug_Flow_Reactor(g, {'NH3', 'N2', 'H2'}, 1.0, dx);
%     X_NH3 = [X_NH3, output(3, end)];
%     X_N2 = [X_N2, output(4, end)];
%     X_H2 = [X_H2, output(5, end)];
% end
% 
% disp(['CPU time = ' num2str(cputime - t0)]);
% 
% %% Plot the results
% figure(1)
% subplot(1, 2, 1);
% hold on
% plot(T_array1, X_BCL3, 'k');
% plot(T_array1, X_BCL2, 'r');
% plot(T_array1, X_BCL, 'g');
% plot(T_array1, X_HCL, 'b');
% legend('BCl3', 'BCl2', 'BCl', 'HCl');
% xlabel('Temperature (K)');
% ylabel('Mole Fraction');
% axis square
% 
% subplot(1, 2, 2);
% hold on
% plot(T_array2, X_NH3, 'k');
% plot(T_array2, X_N2, 'r');
% plot(T_array2, X_H2, 'g');
% legend('NH3', 'N2', 'H2');
% xlabel('Temperature (K)');
% ylabel('Mole Fraction');
% axis square

%% Plug-flow preduction of gas-phase concentrations vs. time.

t0 = cputime;

T_array3 = [1125];
P0 = 2.0 * 133.32;
X0 = 'BCl3:0.4,NH3:0.6';
speciesList = {'BCl3', 'NH3', 'Cl2BNH2', 'ClB(NH2)2', 'B(NH2)3', 'ClBNH', 'HCl'}; 

figure(2)

for i = 1:length(T_array3)
    fprintf('Solving for temperature = %d\n', T_array3(i));
    g.TPX = {T_array3(i), P0, X0};
    output = Plug_Flow_Reactor(g, speciesList, 0.3, dx);
    tim = output(1, 1:end);
    XX_BCL3 = output(3, 1:end);
    XX_NH3 = output(4, 1:end);
    XX_ADCB = output(5, 1:end);
    XX_DACB = output(6, 1:end);
    XX_TAB = output(7, 1:end);
    XX_ClBNH = output(8, 1:end);
    XX_HCL = output(9, 1:end);
    subplot(1, length(T_array3), i)
    title(sprintf('Temperature = %d K', T_array3(i)))
    hold on
    plot(tim, XX_BCL3, 'k');
    plot(tim, XX_NH3, '--k');
    plot(tim, XX_HCL, '--b');
    plot(tim, XX_ADCB, 'r');
    plot(tim, XX_DACB, '--r')
    plot(tim, XX_TAB, 'g');
    plot(tim, XX_ClBNH, '--g');
    legend('BCl3', 'NH3', 'HCl', 'ADCB', 'DACB', 'TAB', 'ClBNH');
    xlabel('Time (s)');
    ylabel('Mole Fraction');
    axis square
end

figure(3)
plot(tim, output(2, 1:end));

disp(['CPU time = ' num2str(cputime - t0)]);
