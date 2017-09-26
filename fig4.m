%%%
% Visualisation of the action functional surface with equality end constraints for varying parameter values of T
%
close all; clear all;

Veq = 0.5; % balanced excitation and inhibition produces this equilibrium potential
V0 = 0; % After spike reset potential
Vth = 0.9; % Spike generation threshold

cspace = linspace(-2,2,200); % influence of excitatory / inhibitory processes

% One feasible pair of stimulation
RE = 2;
RI = (-Veq - RE*(Veq - 1))/(Veq + 1);

% Derived quantities
gE0 = RE; gI0 = RI;
g0 = gE0 + gI0;
E0 = (1./g0).*(gE0 - gI0);

sigE = sqrt(RE/2);
sigI = sqrt(RI/2);

% Try these values of T
TTs = [0.4, 0.6 , 0.8 , 1 , 1.2 , 1.4];%, 1.6];
% Plotting
figure();
for kk = 1:length(TTs)
    % Compute approximate minimizing value of c1
    TT = TTs(kk);
    tt = linspace(0,TT,200);
    dt = tt(2)-tt(1);
    c1 = (Vth*exp(g0.*TT) - V0 - E0.*(exp(g0.*TT) - 1))./((2./(g0+1)).*(exp(TT.*(g0 + 1)) - 1) - (2./(g0 - 1)).*(exp(TT.*(g0 - 1)) - 1));
    % Compute via Euler integration the end point of the solution process and value of action functional
    [XX,YY] = meshgrid(cspace,cspace);
    ZZ = zeros(size(XX));
    VV = zeros(size(XX));
    for ii = 1:length(cspace)
        for jj = 1:length(cspace)
                % Value of action functional
                ZZ(ii,jj) =  XX(ii,jj).^2.*(exp(2*TT) - 1)./(2*sigE^2) + YY(ii,jj).^2.*(exp(2*TT)-1)./(2*sigI.^2);
                
                % Euler integration
                gEt = XX(ii,jj)*(exp(tt) - exp(-tt));
                gIt = YY(ii,jj)*(exp(tt) - exp(-tt));

                mut = (g0 + gEt + gIt);
                nut = (g0*E0 + gEt - gIt);
                
                gat = cumsum(mut).*dt;
                vt = exp(-gat).*cumsum(nut.*exp(gat)).*dt + V0.*exp(-gat);
                VV(ii,jj) = vt(end);
        end
    end
    % determine minimum from euler simulation
    % select only values that satisfie end constraint & minimize
    constsat = abs(VV - 0.9) < 0.01; 
    ZZsat = ZZ(constsat);
    minval = min(ZZsat); 
    [rmin,cmin] = find(ZZ.*constsat == minval);
    % Plotting
    subplot(2,3,kk)
    set(gca,'FontSize',18)
    hold on; grid on;
    % action functional contour plot
    contour(XX , YY , ZZ , linspace(0,14,15).^2 )
    % end point contour plot
    [~,cconst] = contour(XX , YY , VV, [0.9 , 0.9] , 'LineWidth',3 , 'Color' , [0.7,0.7,0.7]);
    % constraint satisfaction contour plot
    [~,cact] = contour(XX , YY , VV, linspace(0.1,1.9,10) ,  'ShowText','on');
    % theoretical minimum
    ctheo = scatter(c1,-c1,50,'r','filled');
    % empirical minimum
    cempi = scatter(cspace(cmin),cspace(rmin),50, 'g','^','filled','MarkerFaceAlpha',.8);
    if kk == 1
        ylabel('c_3')
    end
    if kk == 5
        xlabel('c_1')
    end
    if kk == 3
        legend([cconst,cact,ctheo,cempi],{'constraint contour','action functional contour','theory c','simulation c'})
    end
    title(sprintf('T = %1.1f',TT))
    
end
