
    % Simulation of an intergrate-and-fire neuron


    % Computation and cognition undergrad - ex2
    % See word document, for instructions

    clear; close all; clc;

    %% Declare simulation parameters (we will work in units of Volts, Siemens, Farads and Sec)
    gL = 1e-4; % conductance (S/cm^2), 1/R
    C = 1e-6; % capacitance (F/cm^2)
    tau = C/gL; % sec - time constant, equal to RC
    TH_in_mV = 30; % mV
    TH = TH_in_mV/1000; % V
    tau_R = 0.004; % sec - Refractoriness

    dt = 0.0001; % time step for numerical integration
    t_final = 1; % sec, duration of numerical simulation
    n = round(t_final/dt); % number of iterations (steps)
    V = NaN(1,n+1); % initialize array for the voltage trace
    V(1) = 0; % initial condition for the voltage - This is V(0)
    %In matlab we arrays are 1-indexed (V(1) is the first element)
    %This is why we have 10001 columns in I

    t = (0:n)*dt; % time vector (n+1 components)

    %% load external current
    %TODO: After you implement 'the plotting section', load the currents
    %(one by one) and plot the applied current and the neuron's response

    %load I_const; % loading I
    %load I_sin;
    %load I_exp;
    load I_step;

    %% Numerical integration
    %TODO: add a comment to each line that explains what it does
    % HINT: Explain it to me like I know matlab pretty well, but I don't have a
    %   clue in neuroscience or differntial equations.
    
    RP_flag = 0; % a binary varient that defines whether there is a refrectorial period. 
    
    for idx = 1:n
        VV = V(idx) + (-V(idx)/tau + I(idx)/C)*dt; % dynamical equation of an RC circuit
        
%if the voltage (VV) is larger then the threshold (TH) the neuron
%wiil fire. Afterwards it will enter a refrectorial period.
        if VV >= TH 
            VV = 0; 
            RP_flag = 1; 
            t_RP_start = t(idx); %records the start time of the RP
        end
        if RP_flag
            VV = 0; %while RP the neuron's voltage is 0.
            if (t(idx) - t_RP_start >= tau_R) % after the neuron stays in RP for tau_R(the time that takes for the neuron to be able again to produce an action potential.) the RP will end.
                RP_flag = 0;
            end
        end

        V(idx+1) = VV; 
    end

    %% The plotting section
    figure('Color','w');
    
    subplot(2,1,1)
    
    y = I*(10^6);
    plot(t,y)
    set(gca,'FontSize',16)
    ylabel('Current [\muA]')
    
    subplot(2,1,2)
    y = V*10^3;
    plot(t,y)
    yline(TH_in_mV,['r','--'])
    hold on

    ylim([0 32]);
    set(gca,'FontSize',16)
    xlabel('Time [S]');
    ylabel('Relative Voltage [mV]')
