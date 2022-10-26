%% BN_Gas - A Gas Phase Reaction Simulation of the B-N-H-Cl system.
% 
% See definition in BN.yaml
% 
% Requires: cantera >= 2.6.0
% Keywords: kinetics

help BN_gas
LoadCantera;
clear all
close all

% Initial Parameters

c1 = 2.390e-7; % J/kmol to kcal/mol

% Create the gas phase object

g = Solution('BNHCL-Allendorf1995.yaml', 'gas');

% Print all the reactions and their dH

for i = 1:g.nReactions
    fprintf('%50s  %13.5g\n', ... 
            g.reactionEqn(i), g.dH_standard(i)*c1);
end

%% Plug-flow prediction of decomposition of BCl3 and NH3 vs. Temperature

t0 = cputime;

T_array1 = linspace(1000, 2000, 21);
P0 = 760 * 133.32;
X0 = 'BCL3:0.05,H2:0.95';

X_BCL3 = [];
X_BCL2 = [];
X_BCL = [];
X_HCL = [];
X_NH3 = [];
X_N2 = [];
X_H2 = [];

for i = 1:length(T_array1)
    fprintf('Solving for temperature = %d\n', T_array1(i));
    g.TPX = {T_array1(i), P0, X0};
    output = Plug_Flow(g, {'BCL3', 'BCL2', 'BCL', 'HCL'}, 1.0);
    X_BCL3 = [X_BCL3, output(2, end)];
    X_BCL2 = [X_BCL2, output(3, end)];
    X_BCL = [X_BCL, output(4, end)];
    X_HCL = [X_HCL, output(5, end)];
end

T_array2 = linspace(800, 2000, 25);
X0 = zeros(1, g.nSpecies);
X0(g.speciesIndex('NH3')) = 1.00;

X_NH3 = [];

for i = 1:length(T_array2)
    fprintf('Solving for temperature = %d\n', T_array2(i));
    g.TPX = {T_array2(i), P0, X0};
    output = Plug_Flow(g, {'NH3', 'N2', 'H2'}, 1.0);
    X_NH3 = [X_NH3, output(2, end)];
    X_N2 = [X_N2, output(3, end)];
    X_H2 = [X_H2, output(4, end)];
end

disp(['CPU time = ' num2str(cputime - t0)]);

%% Plot the results
figure(1)
subplot(1, 2, 1);
hold on
plot(T_array1, X_BCL3, 'k');
plot(T_array1, X_BCL2, 'r');
plot(T_array1, X_BCL, 'g');
plot(T_array1, X_HCL, 'b');
legend('BCl3', 'BCl2', 'BCl', 'HCl');
xlabel('Temperature (K)');
ylabel('Mole Fraction');
axis square

subplot(1, 2, 2);
hold on
plot(T_array2, X_NH3, 'k');
plot(T_array2, X_N2, 'r');
plot(T_array2, X_H2, 'g');
legend('NH3', 'N2', 'H2');
xlabel('Temperature (K)');
ylabel('Mole Fraction');
axis square

%% Plug-flow preduction of gas-phase concentrations vs. time.

t0 = cputime;

T_array3 = [1125, 1325];
P0 = 2.0 * 133.32;
X0 = 'BCL3:0.4,NH3:0.6';
speciesList = {'BCL3', 'NH3', 'Cl2BNH2', 'ClB(NH2)2', 'B(NH2)3', 'ClBNH', 'HCl'}; 

figure(2)

for i = 1:length(T_array3)
    fprintf('Solving for temperature = %d\n', T_array3(i));
    g.TPX = {T_array3(i), P0, X0};
    output = Plug_Flow(g, speciesList, 0.3);
    tim = output(1, 1:end);
    XX_BCL3 = output(2, 1:end);
    XX_NH3 = output(3, 1:end);
    XX_ADCB = output(4, 1:end);
    XX_DACB = output(5, 1:end);
    XX_TAB = output(6, 1:end);
    XX_ClBNH = output(7, 1:end);
    XX_HCL = output(8, 1:end);
    subplot(1, length(T_array3), i)
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
