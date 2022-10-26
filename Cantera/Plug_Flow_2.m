function output = Plug_Flow_2(gas_calc, species_name)
    % Simulation of a plug flow for the BNHCL system. Currently for 
    % gas-phase reactions only. 
   
    % create a reactor, and insert the gas
    r = IdealGasReactor(gas_calc);

    % create a reactor network and insert the reactor
    network = ReactorNet({r});

    nSteps = 100;
    tim = zeros(nSteps);
    temp = zeros(nSteps);
    xx = zeros(nSteps, length(species_name));
    t = 0;
    dt = 1.0e-3;
    for n = 1:1000
      t = t + dt;
      network.advance(t);
      tim(n) = network.time;
      temp(n) = r.T;
      xx(n, 1:end) = gas_calc.moleFraction(species_name);
    end

    output = xx(1000, 1:end);
end
