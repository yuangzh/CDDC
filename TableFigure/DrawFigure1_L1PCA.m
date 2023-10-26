% clear all
% clc;clear all;close all;
% % delete('*.eps')
load('demo_cpu1_L1PCA');
addpath('util');

POSs = {'NorthEast','NorthEast','NorthEast'};
   
max_ys = [150;150;200];
   
for idata = 1:3
One = result{idata};

pcolor = loadcolor;
fig = figure('color','w');
% One.his1(1:1) = []; One.his2(1:1) = []; One.his3(1:1) = []; One.his4(1:1) = []; One.his5(1:1) = []; 
% One.ts1(1:1) = []; One.ts2(1:1) = []; One.ts3(1:1) = []; One.ts4(1:1) = []; One.ts5(1:1) = [];

CCC = min([One.his1(:);One.his2(:);One.his3(:);One.his4(:);One.his5(:)])  ;
One.his1=One.his1 - CCC; One.his2=One.his2 - CCC; One.his3=One.his3 - CCC; One.his4=One.his4 - CCC;One.his5=One.his5 - CCC;

myplot = @plot;
myplot(One.ts1,One.his1,'--ks','LineWidth',5,'MarkerSize',3,'color', pcolor.red); hold on;
myplot(One.ts2,One.his2,'-.ms','LineWidth',5,'MarkerSize',3,'color',  pcolor.blue); hold on;
myplot(One.ts3,One.his3,':k*','LineWidth',4,'MarkerSize',3,'color', pcolor.purple); hold on;
myplot(One.ts4,One.his4,'-b+','LineWidth',4,'MarkerSize',3,'color', pcolor.yellow); hold on;
myplot(One.ts5,One.his5,'-ro','LineWidth',4,'MarkerSize',3,'color', pcolor.green  ); hold on;

hleg=legend('MSCR','PDCA','SubGrad','CD-SCA','CD-SNCA') ;
set(hleg,'FontSize',20,'FontWeight','normal');
set(hleg,'Fontname','times new Roman');
set(hleg,'Location',POSs{idata});
set(gca,'Fontsize', 19);
xlabel('Time (seconds)','FontSize',20)
ylabel('$F(\mathbf{x})-F_{\min}$','FontSize',20,'interpreter','latex')

% set(gca,'xtick',[5 10 20 30 40 50 60 70 80 90]);
% set(gca,'XTickLabel',{'64','128','256','512','1024'})
% set(gcf,'position', [300 100 600 500]);

hiss = [One.his1(:);One.his2(:);One.his3(:);One.his4(:);One.his5(:)];
tss = [One.ts1(:);One.ts2(:);One.ts3(:);One.ts4(:);One.ts5(:)];
max_y = max_ys(idata);  %max(hiss)
min_y = 0; %0.99*min(hiss)
axis([min(tss) max(tss) min_y max_y])
grid on;
fprintf('\n');
set(gcf,'paperpositionmode','auto')
 print(sprintf('%s_%d.eps',mfilename,idata),'-depsc2','-loose');
%  print(sprintf('%s_%d.png',mfilename,idata),'-dpng');
% print(sprintf('%s_%d.pdf',mfilename,idata),'-dpdf', '-r0');


 


end




