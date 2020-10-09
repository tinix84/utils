%% Evolutionary curve fitting using PSO
%
% Bartlomiej UFNALSKI
% Institute of Control and Industrial Electronics
% Warsaw University of Technology
% bartlomiej.ufnalskix@xee.pw.edu.pl (please remove both exes)
% rev. 05.10.2014 for Matlab Central
%
% Obviously, it is nothing new. You can use Matlab's fminsearch or Curve
% Fitting Toolbox. There are also many alternatives such as EzyFit for Matlab, Scilab's optimization tools,
% Octave's optimization tools, etc. As long as your current tool uses a gradient-based approach,
% its success rate strongly depends on starting point in the case of non-convex problems.
% It is then your job to select this point. I found this quite
% challenging when trying to identify Foster-type representation of thermal
% transient impedance of transistors, diodes and heat sinks. So I have
% switched to PSO.
%
% This particular example assumes the 3rd order dynamics. In most cases
% 3rd--5th order system is sufficient to reliably model thermal aspects of your power electronic
% converter. It's rather easy to modify this script to perform any curve
% fitting task you want. And don't be affraid of even 200-dimensional
% problems (see e.g. my PDPSRC or PDMSRC on Matlab Central).
%
%% Initialization
clear all
%
% Your digitized plot goes here. You can use e.g. Plot Digitizer or many
% other digitizers available on Matlab Central.
% Here the trench IGBT module SKM145GB066D from Semikron will serve as an
% example. Note that the thermal chain parameters are given explicitely for
% this particular module but this is not always the case. In fact, it's usually not
% the case. More on thermal modelling:
% http://www.plexim.com/download/documentation .
x_igbt_zth = [1.02e-05, 1.46e-05, 1.93e-05, 2.52e-05, 3.14e-05, 3.92e-05,...
    5.03e-05, 6.71e-05, 8.89e-05, 0.0001185, 0.000159193, 0.000215443, 0.00029157,...
    0.000403431, 0.000541969, 0.000722727, 0.000970911, 0.00132372, 0.00184512, 0.002553,...
    0.00350648, 0.00514681, 0.00733474, 0.0112534, 0.01664, 0.0225198, 0.0332994, 0.0444054,...
    0.0632823, 0.0831518, 0.104528, 0.135335, 0.171386, 0.206112, 0.272833, 0.361154,...
    0.492388, 0.65661, 0.831518, 0.992647];
y_igbt_zth = [0.000438281, 0.000600275, 0.000779386, 0.00100595, 0.00123818, 0.00152402,...
    0.00189824, 0.0024793, 0.0031811, 0.00408155, 0.0052994, 0.00675922, 0.00877605,...
    0.0111273, 0.0141926, 0.0176776, 0.0216298, 0.026623, 0.0333577, 0.0410584,...
    0.0487691, 0.0603848, 0.0747672, 0.0936806, 0.117378, 0.138597, 0.168579, 0.195541,...
    0.228164, 0.252384, 0.271013, 0.28419, 0.292749, 0.299781, 0.30336, 0.30336,...
    0.30336, 0.305166, 0.30336, 0.30336];
rng('shuffle');
asur = 0; % 1 --> the asynchronous update rule; 0 --> the synchronous update rule
evaporation_constant = 1; % Static optimization problem (SOP)
swarm_size = 50;
init_position = 3; % note that parameters are coded using exponential function;
% the initial area covered by the swarm is then [0.001, 1000];
% this should work for most thermal models of
% semiconductor switches and heat sinks in power electronic systems. If the
% optimizer fails to identify your device, you should try setting
% init_position = 5 .
init_speed = 1;
diversity_limit = 0; % Static optimization problem (SOP)
velocity_clamping = 2;
num_of_iter = 10000; % don't worry -- it will take seconds to complete the task at hand
% Constricted PSO
correction_factor = 2.05;
Kappa=2/abs(2-2*correction_factor-sqrt((2*correction_factor)^2-8*correction_factor));
% The problem
num_of_D = 6; % 6 parameters of the Foster-type RC ladder network to be identified
% More declarations
swarm = zeros(swarm_size,4,num_of_D);
swarm_diversity = zeros(1,num_of_D);
swarm_dir = zeros(1,num_of_D);
%% Initial swarm position and velocity
for index=1:swarm_size,
    swarm(index, 1, :) = init_position*(rand(num_of_D,1)-0.5)*2;
    swarm(index, 2, :) = init_speed*(rand(num_of_D,1)-0.5)*2;
end
swarm(:, 4, 1) = 10e10;          % initial best value --> some value out of practical limits
%% The swarm is turned on
for iter = 1 : num_of_iter,
    for n = 1 : swarm_size,
        r1 = 10^swarm(n, 1, 1); % resistances
        r2 = 10^swarm(n, 1, 2);
        r3 = 10^swarm(n, 1, 3);
        c1 = 10^swarm(n, 1, 4); % capacitances
        c2 = 10^swarm(n, 1, 5);
        c3 = 10^swarm(n, 1, 6);
        %% In this example the 3rd order Foster chain is used.
        % Particles are rated here using sum of squared errors.
        % Obviously you can define your own cost function.
        fitness = 0;
        for i = 1:numel(x_igbt_zth),
            t = x_igbt_zth(i);
            fitness = fitness + (r1*(1-exp(-t/c1/r1))+r2*(1-exp(-t/c2/r2))+r3*(1-exp(-t/c3/r3))-y_igbt_zth(i))^2;
            % note that in the loglog scale (see final fig) similar fitting errors
            % seem to be bigger at the beginning of the thermal transient
            % impedance curve; be carefull when interpreting them -- don't
            % be fooled by your linear mind ;)
            % You can switch to time weighted cont function such as
            %fitness = fitness + 1/t * (r1*(1-exp(-t/c1/r1))+r2*(1-exp(-t/c2/r2))+r3*(1-exp(-t/c3/r3))-y_igbt_zth(i))^2;
            % to make loglog graph look nicer, but the question is what are
            % your real objectives. Remember to choose your cost functions
            % wisely. Optimization is all about selecting the cost function
            % that best reflects your needs.
        end
        % pbest evaporation
        if fitness < swarm(n, 4, 1)*evaporation_constant,
            for d=1:num_of_D,
                swarm(n, 3, d) = swarm(n, 1, d);
            end
            swarm(n, 4, 1) = fitness;
        else
            swarm(n, 4, 1)=swarm(n, 4, 1)*evaporation_constant;
        end
        if asur == 1, % gbest selection in ASUR and the update
            [~, gbest] = min(swarm(:, 4, 1));
            % diversity control
            for d=1:num_of_D,
                swarm_diversity(d) = 0.5*(max(swarm(:, 1, d))-min(swarm(:, 1, d)));
            end
            for d=1:num_of_D,
                if swarm_diversity(d) < diversity_limit,
                    swarm_dir(d) = -1;
                else
                    swarm_dir(d) = 1;
                end
            end
            % Updating velocity vectors and positions
            for d = 1 : num_of_D,
                rand1=rand;
                rand2=rand;
                % Position
                swarm(n, 2, d) = Kappa*(swarm(n, 2, d)...
                    + swarm_dir(d)*correction_factor*rand1*(swarm(n, 3, d) - swarm(n, 1, d))...
                    + swarm_dir(d)*correction_factor*rand2*(swarm(gbest, 3, d) - swarm(n, 1, d)));
                % Velocity
                swarm(n, 2, d)= min(max(-velocity_clamping,swarm(n, 2, d)),velocity_clamping);
                swarm(n, 1, d) = swarm(n, 1, d) + swarm(n, 2, d);
            end
        end % ASUR
    end
    if asur == 0, % gbest selection in SUR
        [~, gbest] = min(swarm(:, 4, 1));
    end
    if asur == 0, % SUR
        % diversity control
        for d=1:num_of_D,
            swarm_diversity(d) = 0.5*(max(swarm(:, 1, d))-min(swarm(:, 1, d)));
        end
        for d=1:num_of_D,
            if swarm_diversity(d) < diversity_limit,
                swarm_dir(d) = -1;
            else
                swarm_dir(d) = 1;
            end
        end
        % Updating velocity vectors and positions
        for n = 1 : swarm_size,
            for d = 1 : num_of_D,
                rand1=rand;
                rand2=rand;
                % Position
                swarm(n, 2, d) = Kappa*(swarm(n, 2, d)...
                    + swarm_dir(d)*correction_factor*rand1*(swarm(n, 3, d) - swarm(n, 1, d))...
                    + swarm_dir(d)*correction_factor*rand2*(swarm(gbest, 3, d) - swarm(n, 1, d)));
                % Velocity
                swarm(n, 2, d)= min(max(-velocity_clamping,swarm(n, 2, d)),velocity_clamping);
                swarm(n, 1, d) = swarm(n, 1, d) + swarm(n, 2, d);
            end
        end
    end % SUR
end
r1 = 10^swarm(gbest, 3, 1); % resistances
r2 = 10^swarm(gbest, 3, 2);
r3 = 10^swarm(gbest, 3, 3);
c1 = 10^swarm(gbest, 3, 4); % capacitances
c2 = 10^swarm(gbest, 3, 5);
c3 = 10^swarm(gbest, 3, 6);
figure(1)
plot(x_igbt_zth,y_igbt_zth,'o','Linewidth',2); grid; hold on;
plot(x_igbt_zth,r1*(1-exp(-x_igbt_zth/c1/r1))+r2*(1-exp(-x_igbt_zth/c2/r2))+r3*(1-exp(-x_igbt_zth/c3/r3)),'r','Linewidth',2); hold off;
legend('digitized','identified','Location','best');
xlabel('t [s]','Fontsize',12); ylabel('Z_{th(j-c)} [K/W]','Fontsize',12);
title('Transient thermal impedance of IGBT');
figure(2)
loglog(x_igbt_zth,y_igbt_zth,'o','Linewidth',2); grid; hold on;
loglog(x_igbt_zth,r1*(1-exp(-x_igbt_zth/c1/r1))+r2*(1-exp(-x_igbt_zth/c2/r2))+r3*(1-exp(-x_igbt_zth/c3/r3)),'r','Linewidth',2); hold off;
legend('digitized','identified','Location','best');
xlabel('t [s]','Fontsize',12); ylabel('Z_{th(j-c)} [K/W]','Fontsize',12);
title('Transient thermal impedance of IGBT','Fontsize',12);
print(gcf,'-djpeg90','-r100','Zth_jc_IGBT');
% That's all folks! Good luck with your thermal models in power electronics.