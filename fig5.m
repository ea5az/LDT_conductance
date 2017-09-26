%%%
% Visualisation of the approximated optimal voltage trajectory for varying values of T. and fixed value of RE.
%
close all; clear all;

Veq = 0.5; % balanced excitation and inhibition produces this equilibrium potential
V0 = 0; % After spike reset potential
Vth = 0.9; % Spike generation threshold


% One feasible pair of stimulation
RE = 2;
RI = (-Veq - RE*(Veq - 1))/(Veq + 1);

% Derived quantities
gE0 = RE; gI0 = RI;
g0 = gE0 + gI0;
E0 = (1./g0).*(gE0 - gI0);

sigE = sqrt(RE/2);
sigI = sqrt(RI/2);

% Simulate up to time 3
TT = linspace(0,3,20);

% Approximated optimal coefficient
c1 = (Vth*exp(g0.*TT) - V0 - E0.*(exp(g0.*TT) - 1))./((2./(g0+1)).*(exp(TT.*(g0 + 1)) - 1) - (2./(g0 - 1)).*(exp(TT.*(g0 - 1)) - 1));
% Plotting
fig = figure();
set(fig,'defaultAxesColorOrder',[[0 0 0];[0 0 0]]);
hold on;
for ii = 1:length(c1)
    % Closed form solution for voltage trajectory
    Taug = linspace(0,TT(ii));
    cc = c1(ii);
    vt = exp(-g0.*Taug).*( ...
        ((2*cc)/(g0 + 1))*(exp((g0+1).*Taug)-1)  ...
        - ((2*cc)/(g0 - 1))*(exp((g0-1).*Taug)-1)  ...
        +E0*(exp(g0*Taug) - 1) + V0 ...
        );
    % Conductance process
    get = cc.*(exp(Taug)-exp(-Taug));
    yyaxis left
    p1=plot(Taug , vt , 'Marker','none' ,'LineStyle' , '-' , 'Color' ,  [0.7 0.7 0.7] , 'LineWidth' , 2);
    yyaxis right
    p3=plot(Taug , get , '—m');
end
p2=plot(TT,c1,'—r');
ylim([0,2])
xlabel('Spike time T')
ylabel('conductance process g_E(t)')
yyaxis left
ylabel('Membrane potential v(t)')
legend([p1,p3,p2],{'v(t)','g_E(t)','c_1(t)'},'Location','SouthEast')
set(gca,'FontSize',20)
