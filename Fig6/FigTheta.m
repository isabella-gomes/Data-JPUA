phi12 = 10:10:120;

%potmn20VLC30aux = load('curvaPot.mat');
%potmn20VLC30 = sum(potmn20VLC30aux.sumexp)/10;
pot0VLC30aux = load('curvaRF0VLC30.mat');
pot0VLC30 = pot0VLC30aux.sumtotal;
pot20VLC30aux = load('curvaRF20VLC30.mat');
pot20VLC30=pot20VLC30aux.sumtotal;
pot0VLC40aux = load('curvaRF0VLC40.mat');
pot0VLC40 = pot0VLC40aux.sumtotal;
pot20VLC40aux = load('curvaRF20VLC40.mat');
pot20VLC40=pot20VLC40aux.sumtotal;

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
%bCSI = plot(phi12,pot0VLC30CSI/100,'Color',[0.4940 0.1840 0.5560]);
cCSI = plot(phi12,pot20VLC30CSI/100,'Color',[0.6350 0.0780 0.1840]);
%eCSI = plot(phi12,pot0VLC40CSI/100,'-o','Color',[0.4940 0.1840 0.5560]);
fCSI = plot(phi12,pot20VLC40CSI/100,'-o','Color',[0.6350 0.0780 0.1840]);
%b1CSI = plot(phi12,pot0VLC30CSI/100,'o','Color',[0.4940 0.1840 0.5560]);
c1CSI = plot(phi12,pot20VLC30CSI/100,'o','Color',[0.6350 0.0780 0.1840]);

c1CSIa = plot(phi12,-1*ones(1,length(phi12)),'ko');

% a = plot(k,potmn20VLC30/100,'Color',[0 0.4470 0.7410]);
%b = plot(phi12,pot0VLC30/100,'Color',[0.4940 0.1840 0.5560]);
c = plot(phi12,pot20VLC30/100,'Color',[0.6350 0.0780 0.1840]);
%e = plot(phi12,pot0VLC40/100,'-s','Color',[0.4940 0.1840 0.5560]);
f = plot(phi12,pot20VLC40/100,'Color',[0.6350 0.0780 0.1840]);
%b1 = plot(phi12,pot0VLC30/100,'s','Color',[0.4940 0.1840 0.5560]);
f1 = plot(phi12,pot20VLC40/100,'s','Color',[0.6350 0.0780 0.1840]);
c1 = plot(phi12,pot20VLC30/100,'s','Color',[0.6350 0.0780 0.1840]);
c1a = plot(phi12,-2*ones(1,length(phi12)),'ks');

h = [c1CSIa,c1a,c,f];
 legend(h,'JPUA Perfect CSI','JPUA','$p_n^{\mathrm{VLC}}$ = 30 dBm','$p_n^{\mathrm{RF}}$ = 40 dBm')