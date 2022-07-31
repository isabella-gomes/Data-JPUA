k = -20:5:20;

%potmn20VLC30aux = load('curvaPot.mat');
%potmn20VLC30 = sum(potmn20VLC30aux.sumexp)/10;
pot0VLC30aux = load('curva030.mat');
pot0VLC30 = pot0VLC30aux.sumTotal;
pot20VLC30aux = load('curva2030.mat');
pot20VLC30=pot20VLC30aux.sumTotal;
pot0VLC40aux = load('curva040.mat');
pot0VLC40 = pot0VLC40aux.sumTotal;
pot20VLC40aux = load('curva2040.mat');
pot20VLC40=pot20VLC40aux.sumTotal;

pot0VLC30auxCSI = load('curvaPot0VLC30.mat');
pot0VLC30CSI = pot0VLC30auxCSI.sumexp;
pot20VLC30auxCSI = load('curvaPot20VLC30.mat');
pot20VLC30CSI=pot20VLC30auxCSI.sumexp;
pot0VLC40auxCSI = load('curvaPot0VLC40.mat');
pot0VLC40CSI = pot0VLC40auxCSI.sumexp;
pot20VLC40auxCSI = load('curvaPot20VLC40.mat');
pot20VLC40CSI=pot20VLC40auxCSI.sumexp;

figure
hold on

bCSI = plot(k,pot0VLC30CSI/100,'Color',[0.4940 0.1840 0.5560]);
cCSI = plot(k,pot20VLC30CSI/100,'Color',[0.6350 0.0780 0.1840]);
eCSI = plot(k,pot0VLC40CSI/100,'-o','Color',[0.4940 0.1840 0.5560]);
fCSI = plot(k,pot20VLC40CSI/100,'-o','Color',[0.6350 0.0780 0.1840]);
b1CSI = plot(k,pot0VLC30CSI/100,'o','Color',[0.4940 0.1840 0.5560]);
c1CSI = plot(k,pot20VLC30CSI/100,'o','Color',[0.6350 0.0780 0.1840]);

c1CSIa = plot(k,-1*ones(1,length(k)),'ko');

b = plot(k,pot0VLC30/100,'Color',[0.4940 0.1840 0.5560]);
c = plot(k,pot20VLC30/100,'Color',[0.6350 0.0780 0.1840]);
e = plot(k,pot0VLC40/100,'-s','Color',[0.4940 0.1840 0.5560]);
f = plot(k,pot20VLC40/100,'-s','Color',[0.6350 0.0780 0.1840]);
b1 = plot(k,pot0VLC30/100,'s','Color',[0.4940 0.1840 0.5560]);
c1 = plot(k,pot20VLC30/100,'s','Color',[0.6350 0.0780 0.1840]);

c1a = plot(k,-2*ones(1,length(k)),'ks');

h = [c1CSIa,c1a,b,c];
 legend(h,'JPUA Perfect CSI','JPUA','$P_s^{\mathrm{RF}}$ = 0 dBm','$P_s^{\mathrm{RF}}$ = 20 dBm')