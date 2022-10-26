function output = Plug_Flow_1(gas_calc, species_name)
    % Simulation of a plug flow for the BNHCL system.
    % Currently for gas-phase reactions only.
    %
    % :param gas_calc:
    %    Object of class 'Solution'.
    % :param species_name:
    %    String or array of string of species names.
    % :return:
    %    Mole fraction of species at the end of 1 second. 
   
    % create a reactor, and insert the gas
    r = IdealGasReactor(gas_calc);

    % create a reactor network and insert the reactor
    network = ReactorNet({r});

    nSteps = 1000;
    tim = zeros(nSteps);
    temp = zeros(nSteps);
    xx = zeros(nSteps, length(species_name));
    t = 0;
    dt = 1.0e-3;
    for n = 1:nSteps
      t = t + dt;
      network.advance(t);
      tim(n) = network.time;
      temp(n) = r.T;
      xx(n, 1:end) = gas_calc.moleFraction(species_name);
    end

    output = xx(1000, 1:end);
end
