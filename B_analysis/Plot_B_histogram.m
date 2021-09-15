%Print B histograms
B = readtable('B_6.txt');
B = B{:,:};
B_nvar = readtable('B_6_nvar.txt');
B_nvar = B_nvar{:,:};
hist_1 = subplot(1,2,1);
hist_2 = subplot(1,2,2);
subplot(hist_2);
h = histogram(B_nvar,100);
title('Without variance reduction');
subplot(hist_1);
histogram(B,'BinWidth',h.BinWidth);
title('With variance reduction');
linkaxes([hist_1 hist_2],'xy')
saveas(gcf,'B_6.png')