function output = Plug_Flow(gas_calc, species_name, total_time)
    % Simulation of a plug flow for the BNHCL system.
    % Currently for gas-phase reactions only.
    %
    % :param gas_calc:
    %    Object of class 'Solution'.
    % :param species_name:
    %    String or array of string of species names.
    % :param total_time:
    %    Total amount of time for simulation. 
    % :return:
    %    Mole fraction of species at the end of 1 second. 
   
    % create a reactor, and insert the gas
    r = IdealGasReactor(gas_calc);

    % create a reactor network and insert the reactor
    network = ReactorNet({r});

    dt = 1.0e-3;
    nSteps = total_time/dt;
    tim = zeros(1, nSteps);
    xx = zeros(length(species_name), nSteps);
    output = zeros(length(species_name) + 1, nSteps);
    t = 0;
    dt = 1.0e-3;
    for n = 1:nSteps
      t = t + dt;
      network.advance(t);
      tim(n) = network.time;
      xx(1:end, n) = gas_calc.moleFraction(species_name);
    end

    output(1, 1:end) = tim(1:end);
    output(2:end, 1:end) = xx;
end
