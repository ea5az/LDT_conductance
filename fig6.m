%%%
% Visualisation of the large deviation theory cdf approximation
%
close all; clear all;

Veq = 0.5; % balanced excitation and inhibition produces this equilibrium potential
V0 = 0; % After spike reset potential
Vth = 0.9; % Spike generation threshold

TT = linspace(0,3);

dt = TT(2)-TT(1);

REs = 10.^linspace(log10(1.1),log10(100));

surfmat = []; surfmatpdf = []; EFPT = [];

for ii = 1:length(REs)
    RE = REs(ii);
    RI = (-Veq - RE*(Veq - 1))/(Veq + 1);
    
    %plot(RE , RI)
    
    gE0 = RE; gI0 = RI;
    g0 = gE0 + gI0;
    E0 = (1./g0).*(gE0 - gI0);
    
    sigE = sqrt(RE/2);
    sigI = sqrt(RI/2);
   
    
    c1 = (Vth*exp(g0.*TT) - V0 - E0.*(exp(g0.*TT) - 1))./ ...
      ((2./(g0+1)).*(exp(TT.*(g0 + 1)) - 1) - (2./(g0 - 1)).*(exp(TT.*(g0 - 1)) - 1));
    
    J = c1.^2.*(exp(2*TT) - 1)./(2*sigE^2) + c1.^2.*(exp(2*TT)-1)./(2*sigI.^2);
    
    surfmat = [surfmat ; exp(-(J-min(J)))];
    
    surfmatpdf = [surfmatpdf ; diff(exp(-(J-min(J))))];
    pdfApp = diff(exp(-(J-min(J))));
    pdfApp(isnan(pdfApp)) = 0;
    
    EFPT = [EFPT dot(pdfApp,TT(1:end-1))];
end
figure();
subplot(1,2,1)
hold on;
[X , Y] = meshgrid(log10(REs) , TT);
[~,c1] = contour(X , Y , transp(surfmat), linspace(0,0.9,10),  'ShowText','on','LineWidth',2);
c2 = plot(log10(REs),EFPT ,'—r','LineWidth',3);
ylabel('First Passage Time T')
legend([c1,c2],{‚Contour Line of FPT CDF','EFPT'},'Location','SouthEast')
xlabel('Excitatory Firing Rate R_E')
xtick = xticks();
xticklabels(round(10.^xtick,2,'significant'));
ylim([0,2.5])
set(gca,'FontSize',20)
axes('Position',[.2 .8 .15 .15])
box on; grid on;
plot(log10(REs) , 1./EFPT , 'r' ,'LineWidth',3)
xlabel('R_E')
ylabel('1/EFPT')
subplot(1,2,2)
hold on;
[X2 , Y2] = meshgrid(log10(REs) , TT(1:end-1));
[~,c1] = contour(X2 , Y2 , transp(surfmatpdf), round(1000*linspace(0.001,0.04,9))/1000, 'ShowText','on','LineWidth',2);
c2 = plot(log10(REs),EFPT,'—r','LineWidth',3);
legend([c1,c2],{'Contour Line of FPT PDF','EFPT'},'Location','SouthEast')
xlabel('Excitatory Firing Rate R_E')
xtick = xticks();
xticklabels(round(10.^xtick,2,'significant'));
ylim([0,2.5])
set(gca,'FontSize',20)
