function output = Plug_Flow_Reactor(gas_calc, species_name, residence_time, dx)
    % Simulation of an ideal plug flow reactor for the BNHCL system.
    % Currently for gas-phase reactions only.
    %
    % :param gas_calc:
    %    Object of class 'Solution'.
    % :param species_name:
    %    String or array of string of species names.
    % :param residence_time:
    %    Residence time of the PFR.
    % :param dx:
    %    Length of each slice of the PFR. 
    % :return:
    %    Mole fraction of species at the exist of the PFR. 

    %% Initial Setup. Reactor parameters are taken from Allendorf 1998. 

    % Inlet Area, in m^2
    A_in = 0.064 * 0.064 * pi / 4;
    % Volumetric flow rate into the reactor, in m^3
    VFR = 2000e-6; % 2000 sccm
    % Mass flow velocity into the reactor, in m/s
    v = VFR / A_in;
    % Mass flow rate into the reactor, in kg/s.
    mdot_calc = A_in * v * gas_calc.D;
    % Create a vector of all the PFR slices. 
    x_calc = [0, dx];

    nsp = gas_calc.nSpecies;

    % These parameters are carried over from the Cantera PFR_Solver example. 
    % Set to zero for cylindrical PFR. 
    k = 0;
    dAdx = 0;

    % Initialize arrays for T, Y, and rho at each location:
    T_calc = [gas_calc.T];
    Y_calc = [gas_calc.Y];
    rho_calc = [gas_calc.D];
    XX_calc = [gas_calc.X(gas_calc.speciesIndex(species_name))];

    t = [0, 0.001];
    i = 2;

    while t(i) <= residence_time

        % Solver location indicator
        % fprintf('Solving reactor %d\n', i)

    %--------------------------------------------------------------------------
    %------The values of variables at previous location are given as initial---
    %------values to the current iteration and the limits of the current-------
    %--------------reactor and the gas entering it are being set---------------
    %--------------------------------------------------------------------------
        inlet_soln(1) = rho_calc(i-1);
        inlet_soln(2) = T_calc(i-1);
        inlet_soln(3: nsp+2) = Y_calc(i-1, :);
        limits = [x_calc(i-1), x_calc(i)];
        gas_calc.TDY = {T_calc(i-1), rho_calc(i-1), Y_calc(i-1, :)};
        options = odeset('RelTol', 1.e-10, 'AbsTol', 1e-10,...
                        'InitialStep', 1e-8, 'NonNegative', 1);
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
        % These values are passed onto the ode15s solver
        [~,y] = ode15s(@PFR_solver, limits, inlet_soln, options, ...
                        gas_calc, mdot_calc, A_in, dAdx, k);

        T_calc = [T_calc, y(end, 2)];
        rho_calc = [rho_calc, y(end, 1)];
        Y_calc = [Y_calc; y(end, 3:nsp+2)];
        XX_calc = [XX_calc; gas_calc.X(gas_calc.speciesIndex(species_name))];
        x_calc = [x_calc, dx * i];
        vx = mdot_calc ./ (A_in * rho_calc(end));
        t = [t, t(i) + dx / vx];
        i = i + 1;
    end

    output = zeros(2 + length(species_name), i - 1);

    output(1, 1:end) = t(1:end-1);
    output(2, 1:end) = T_calc(1:end);
    output(3:end, 1:end) = XX_calc';

    fprintf('Reactor length is %.3f\n', (i - 1) * dx);

end
