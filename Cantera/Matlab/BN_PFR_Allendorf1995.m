%% BN_PFR_Allendorf1995 - A Gas Phase Reaction Simulation of the B-N-H-Cl system.
% 
% A plug-flow simulation of the B-N-H-Cl system using gas-phase mechanism
% data from M. Allendorf et al. 1995 and Thermochemistry data from M. Allendorf
% et al. 1997. The first part prints out standard enthalpy of reaction for 
% all reactions in the kinetics model. The second part creates concentration vs. 
% temperature plots for BCl3 and NH3 decomposition over 1 s residence time. 
% The third part models the evolution of concentrations of selected species 
% over a period of 0 to 300 ms. 
% 
% Requires: cantera >= 2.6.0
% Keywords: kinetics, thermochemistry, boron nitride. 

help BN_PFR_Allendorf1995
LoadCantera;
clear all
close all

% Initial Parameters

c1 = 2.390e-7; % J/kmol to kcal/mol

% Create the gas phase object

fname1 = '/home/ssun30/Work/Boron_Nitride/YAML_Files/BNHCL_Allendorf_NASA7.yaml';
fname2 = '/home/ssun30/Work/Boron_Nitride/YAML_Files/BNHCL_Allendorf_ConstCp.yaml';
fname3 = '/home/ssun30/Work/Boron_Nitride/YAML_Files/BNHCL_Allendorf_PiecewiseGibbs.yaml';
phasename = 'gas';

g = Solution(fname1, phasename);
g2 = Solution(fname2, phasename);
g3 = Solution(fname3, phasename);

g.TP = {300, oneatm};
g2.TP = {300, oneatm};
g3.TP = {300, oneatm};

% Print all the reactions and their dH

% for i = 1:g.nReactions
%     fprintf('%50s  %13.5g  %13.5g\n', ... 
%             g.reactionEqn(i), g.dH_standard(i)*c1, g2.dH_standard(i)*c1);
% end
% 
% T_array0 = [300, 600, 1000, 1500, 2000, 2500];
% GG1 = [];
% GG2 = [];
% GG3 = [];
% 
% prompt = 'Please select a species:\n';
% Species = input(prompt); 
% 
% while ~ismember(Species, g.speciesNames)
%     prompt = 'The species is not in the list, please enter again:\n';
%     Species = input(prompt);
% end
% 
% fprintf('Species is %s\n', Species);
% 
% for i = 1:length(T_array0)
%     g.TPX = {T_array0(i), oneatm, [Species, ':1.0']};
%     g2.TPX = {T_array0(i), oneatm, [Species, ':1.0']};
%     g3.TPX = {T_array0(i), oneatm, [Species, ':1.0']};
%     fprintf('Temperature is %.1d\n', T_array0(i))
%     fprintf('Enthalpy is %7.5g  %7.5g\n  %7.5g\n', g.H*c1, g2.H*c1, g3.H*c1)
%     fprintf('Entropy is %7.5g  %7.5g\n  %7.5g\n', g.S*c1*1000, g2.S*c1*1000, g3.S*c1*1000)
%     fprintf('Gibbs Free Energy is %7.5g  %7.5g\n  %7.5g\n', g.G*c1, g2.G*c1, g3.G*c1)
%     GG1 = [GG1, g.G*c1];
%     GG2 = [GG2, g2.G*c1];
%     GG3 = [GG3, g3.G*c1];
% end
% 
% figure(1);
% hold on
% plot(T_array0, GG1, '*r');
% plot(T_array0, GG2, '-g');
% plot(T_array0, GG3, 'ob');
% legend('NASA7', 'ConstCp', 'PiecewiseGibbs');
% xlabel('Temperature (K)');
% ylabel('Gibbs Free Energy(kcal/mol)');
% title(Species);
% axis square

%% Plug-flow prediction of decomposition of BCl3 and NH3 vs. Temperature

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
%     output = Plug_Flow(g, {'BCl3', 'BCl2', 'BCl', 'HCl'}, 1.0);
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
%     output = Plug_Flow(g, {'NH3', 'N2', 'H2'}, 1.0);
%     X_NH3 = [X_NH3, output(3, end)];
%     X_N2 = [X_N2, output(4, end)];
%     X_H2 = [X_H2, output(5, end)];
%     HH = g.enthalpies_RT.*1.9872.*g.T;
% end
% 
% disp(['CPU time = ' num2str(cputime - t0)]);
% 
% %% Plot the results
% figure(2)
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

% Plug-flow preduction of gas-phase concentrations vs. time.

t0 = cputime;

T_array3 = [1125, 1325];
P0 = 2.0 * 133.32;
X0 = 'BCl3:0.4,NH3:0.6';
speciesList = {'BCl3', 'NH3', 'Cl2BNH2', 'ClB(NH2)2', 'B(NH2)3', 'ClBNH', 'HCl'}; 

figure(3)

for i = 1:length(T_array3)
    fprintf('Solving for temperature = %d\n', T_array3(i));
    g.TPX = {T_array3(i), P0, X0};
    output = Plug_Flow(g, speciesList, 0.3);
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

disp(['CPU time = ' num2str(cputime - t0)]);
