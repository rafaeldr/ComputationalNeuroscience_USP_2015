% Primeira lista de Exercícios - Introdução à Neurociência Computacional
% Aluno: Rafael Delalibera Rodrigues
% Programa de Pós-Graduação em Computação Aplicada [FFCLRP-USP Ribeirão Preto] - Mestrado

% Exercício 4: Convolução

figure;
f = rand(100,1);
s1 = subplot(6,3,[1 2]);
p1 = plot(f);
title(s1,'Function');
set(s1,'xlim',[0 120]);
% BOXCAR
w1 = rectwin(10);
c1 = conv(f,w1);
sc1 = subplot(6,3,[4 5]);
pc1 = plot(c1);
sw1 = subplot(6,3,6);
pw1 = plot(w1);
%title(sc1,'Convolution');
title(sw1,'Window Type: Boxcar');
% TaylorWin
w2 = taylorwin(10);
c2 = conv(f,w2);
sc2 = subplot(6,3,[7 8]);
pc2 = plot(c2);
sw2 = subplot(6,3,9);
pw2 = plot(w2);
%title(sc2,'Convolution');
title(sw2,'Window Type: TaylorWin');
% Hann
w3 = hann(10);
c3 = conv(f,w3);
sc3 = subplot(6,3,[10 11]);
pc3 = plot(c3);
sw3 = subplot(6,3,12);
pw3 = plot(w3);
%title(sc3,'Convolution');
title(sw3,'Window Type: Hann');
% Blackman
w4 = blackman(10);
c4 = conv(f,w4);
sc4 = subplot(6,3,[13 14]);
pc4 = plot(c4);
sw4 = subplot(6,3,15);
pw4 = plot(w4);
%title(sc4,'Convolution');
title(sw4,'Window Type: Blackman');
% Bartlett
w5 = bartlett(10);
c5 = conv(f,w5);
sc5 = subplot(6,3,[16 17]);
pc5 = plot(c5);
sw5 = subplot(6,3,18);
pw5 = plot(w5);
%title(sc5,'Convolution');
title(sw5,'Window Type: Bartlett');


% Exercício 6: Convolução Matricial - Núcleo 10 Pontos

figure;
f = rand(100,1);
n = length(f);
w = [10 9 8 7 6 5 4 3 2 1];
wt = transpose(w);
y = convmtx(wt,n);
c = y*f;
s1 = subplot(2,1,1);
p1 = plot(f);
s2 = subplot(2,1,2);
p2 = plot(c);
title(s1,'Function');
title(s2,'Convolution');
set(s1,'xlim',[0 120]);


% Exercício 8: MATLAB - Projeto 1

load('projcomput1.mat')
whos

% 8.a) Spikes
sum(disparos)

% 8.b) Duração do experimento (em segundos)
max(tempo)/1000

% 8.c) Taxa de disparos (em Hz)
sum(disparos)/(max(tempo)/1000)

% 8.d) Spikes durante primeira metade do experimento
sum(disparos(1:length(disparos)/2))

% 8.e) Spikes durante primeira metade do experimento
sum(disparos(1:length(disparos)/2))/ ((max(tempo(1:length(tempo)/2)))/1000)

% 8.f) Máximo do estímulo
max(estimulo)

% 8.g) Mínimo do estímulo
min(estimulo)

% 8.h) Tempo do centésimo spike (em segundos)
tempo(max(find(disparos,100)))/1000

% 8.i) Valor médio do estímulo
sum(estimulo)/length(estimulo)
mean(estimulo)

% 8.j) Desvio padrão do estímulo
std(estimulo)

% 8.k) Histograma
figure;
data = diff(tempo(find(disparos)));
media = mean(data);
desv = std(data);
CV = desv/media;
xbins = 0.5:0.5:max(data);
hist(data,xbins);
dim = [.6 .8 .2 .1];
str = strcat('CV =',{' '},num2str(CV,2));
annotation('textbox',dim,'String',str);

% 8.l) Estímulos X Trem de Disparos
figure;
x = 1:1000;
y1 = estimulo(1:1000);
y2 = disparos(1:1000);
[n, e, s] = plotyy(x,y1,x,y2,'plot','stem');
set(n,'xlim',[0 1000]);
set(n(1),'ylim',[-25 50]);
set(n(2),'ylim',[0 4]);
ylabel(n(1),'Estímulo','color','red');
ylabel(n(2),'Spike','color','green');
set(e,'LineWidth', 2);
set(e,'color', 'red');
set(s, 'color', 'green');
set(s, 'marker','none');

% 8.m) Convolução Boxcar
figure;
w = ones(101,1)/101;
d = double(disparos); %vetor não pode ser logical
c = conv(d,w);
prob = c(51:end-50);
ajuste = (sum(disparos)/(max(tempo)/1000))/mean(prob); %Ajuste para Hz
taxa = prob*ajuste;
[p, t, e] = plotyy(tempo,taxa,tempo,estimulo,'plot','plot');
ylabel(p(1),'Taxa de Disparos (Hz)');
ylabel(p(2),'Estimulo','color','red');
set(p(2),'ycolor', 'red');
xlabel(p(1),'Tempo (ms)');
set(e, 'color', 'red');

% 8.n) Modelo Linear
w = ones(101,1)/101;
d = double(disparos); %vetor não pode ser logical
c = conv(d,w);
prob = c(51:end-50);
ajuste = (sum(disparos)/(max(tempo)/1000))/mean(prob); %Ajuste para Hz
taxa = prob*ajuste;
p = polyfit(taxa,estimulo,1);
estimulof = repmat(p,5000,1);
estimulof(:,1) = estimulo * p(:,1);
%Gráfico 1 (versus Tempo)
figure;
[p, t, ef] = plotyy(tempo,taxa,tempo,estimulof,'plot','plot');
ylabel(p(1),'Taxa de Disparos (Hz)');
ylabel(p(2),'Estimulo Fit','color','red');
set(p(2),'ycolor', 'red');
xlabel(p(1),'Tempo (ms)');
set(ef, 'color', 'red');
set(p(2),'ylim',[-20 30]);
%Gráfico 2 (versus Estimulo)
figure;
[p2, t2, ef2] = plotyy(estimulo,taxa,estimulo,estimulof,'plot','plot');
ylabel(p2(1),'Taxa de Disparos (Hz)');
ylabel(p2(2),'Estimulo Fit','color','red');
set(p2(2),'ycolor', 'red');
xlabel(p2(1),'Estimulo');
set(ef2, 'color', 'red');
set(p2(2),'ylim',[-20 30]);
set(t2,'linestyle','none','marker','o');
set(ef2,'LineWidth', 3);
R101 = corrcoef(taxa,estimulo);

% 8.o) Diferentes Núcleos
% NÚCLEO 1001
w = ones(1001,1)/1001;
d = double(disparos); %vetor não pode ser logical
c = conv(d,w);
prob = c(501:end-500);
ajuste = (sum(disparos)/(max(tempo)/1000))/mean(prob); %Ajuste para Hz
taxa = prob*ajuste;
p = polyfit(taxa,estimulo,1);
estimulof = repmat(p,5000,1);
estimulof(:,1) = estimulo * p(:,1);
%Gráfico 1 (versus Tempo)
figure;
[p, t, ef] = plotyy(tempo,taxa,tempo,estimulof,'plot','plot');
ylabel(p(1),'Taxa de Disparos (Hz)');
ylabel(p(2),'Estimulo Fit','color','red');
set(p(2),'ycolor', 'red');
xlabel(p(1),'Tempo (ms)');
set(ef, 'color', 'red');
%set(p(2),'ylim',[-20 30]);
%Gráfico 2 (versus Estimulo)
figure;
[p2, t2, ef2] = plotyy(estimulo,taxa,estimulo,estimulof,'plot','plot');
ylabel(p2(1),'Taxa de Disparos (Hz)');
ylabel(p2(2),'Estimulo Fit','color','red');
set(p2(2),'ycolor', 'red');
xlabel(p2(1),'Estimulo');
set(ef2, 'color', 'red');
%set(p2(2),'ylim',[-20 30]);
set(t2,'linestyle','none','marker','o');
set(ef2,'LineWidth', 3);
R1001 = corrcoef(taxa,estimulo);
% ###################
% NÚCLEO 11
w = ones(11,1)/11;
d = double(disparos); %vetor não pode ser logical
c = conv(d,w);
prob = c(6:end-5);
ajuste = (sum(disparos)/(max(tempo)/1000))/mean(prob); %Ajuste para Hz
taxa = prob*ajuste;
p = polyfit(taxa,estimulo,1);
estimulof = repmat(p,5000,1);
estimulof(:,1) = estimulo * p(:,1);
%Gráfico 1 (versus Tempo)
figure;
[p, t, ef] = plotyy(tempo,taxa,tempo,estimulof,'plot','plot');
ylabel(p(1),'Taxa de Disparos (Hz)');
ylabel(p(2),'Estimulo Fit','color','red');
set(p(2),'ycolor', 'red');
xlabel(p(1),'Tempo (ms)');
set(ef, 'color', 'red');
%Gráfico 2 (versus Estimulo)
figure;
[p2, t2, ef2] = plotyy(estimulo,taxa,estimulo,estimulof,'plot','plot');
ylabel(p2(1),'Taxa de Disparos (Hz)');
ylabel(p2(2),'Estimulo Fit','color','red');
set(p2(2),'ycolor', 'red');
xlabel(p2(1),'Estimulo');
set(ef2, 'color', 'red');
set(t2,'linestyle','none','marker','o');
set(ef2,'LineWidth', 3);
R11 = corrcoef(taxa,estimulo);

% 8.p) Autocorrelação (Estímulo)
figure;
ac = xcorr(estimulo, 1000);
plot(ac);

% 8.q) Autocorrelação (Trem de Disparos)
figure;
ac = xcorr(disparos, 1000);
plot(ac);

% 8.r) Correlação Cruzada (Estímulo x Trem de Disparos)
figure;
disparos = double(disparos);
ac = xcorr(estimulo, disparos, 1000);
p = plot(ac);

% Exercício 9: MATLAB - Projeto Dados Reais

load('projcomput1_Dipteralab.mat')
whos

% 9.a) Spikes
sum(disparos)

% 9.b) Duração do experimento (em segundos)
max(tempo)

% 9.c) Taxa de disparos (em Hz)
sum(disparos)/max(tempo)

% 9.d) Spikes durante primeira metade do experimento
sum(disparos(1:length(disparos)/2))

% 9.e) Spikes durante primeira metade do experimento
sum(disparos(1:length(disparos)/2))/max(tempo(1:length(tempo)/2))

% 9.f) Máximo do estímulo
max(estimulo)

% 9.g) Mínimo do estímulo
min(estimulo)

% 9.h) Tempo do centésimo spike (em segundos)
tempo(max(find(disparos,100)))

% 9.i) Valor médio do estímulo
sum(estimulo)/length(estimulo)
mean(estimulo)

% 9.j) Desvio padrão do estímulo
std(estimulo)

% 9.k) Histograma
figure;
data = diff(tempo(find(disparos)))*1000; % para ms
media = mean(data);
desv = std(data);
CV = desv/media;
xbins = 0.5:0.5:max(data);
hist(data,xbins);
dim = [.6 .8 .2 .1];
str = strcat('CV =',{' '},num2str(CV,2));
annotation('textbox',dim,'String',str);
tabulate(data); % para auxiliar a analise

% 9.l) Estímulos X Trem de Disparos
figure;
x = 1:1000;
y1 = estimulo(1:1000);
y2 = disparos(1:1000);
[n, e, s] = plotyy(x,y1,x,y2,'plot','stem');
set(n,'xlim',[0 1000]);
set(n(2),'ylim',[0 4]);
ylabel(n(1),'Estímulo','color','red');
ylabel(n(2),'Spike','color','green');
set(e,'LineWidth', 2);
set(e,'color', 'red');
set(s, 'color', 'green');
set(s, 'marker','none');

% 9.m) Convolução Boxcar
figure;
w = ones(101,1)/101;
d = double(disparos); %vetor não pode ser logical
c = conv(d,w);
prob = c(51:end-50);
ajuste = (sum(disparos)/max(tempo))/mean(prob); %Ajuste para Hz
taxa = prob*ajuste;
[p, e, t] = plotyy(tempo,estimulo,tempo,taxa,'plot','plot');
ylabel(p(1),'Estimulo');
ylabel(p(2),'Taxa de Disparos (Hz)','color','red');
set(p(2),'ycolor', 'red');
xlabel(p(1),'Tempo (s)');
set(t, 'color', 'red');
set(t,'LineWidth', 2);

% 9.n) Modelo Linear
w = ones(101,1)/101;
d = double(disparos); %vetor não pode ser logical
c = conv(d,w);
prob = c(51:end-50);
ajuste = (sum(disparos)/max(tempo))/mean(prob); %Ajuste para Hz
taxa = prob*ajuste;
p = polyfit(taxa,estimulo,1);
estimulof = repmat(p,5000,1);
estimulof(:,1) = estimulo * p(:,1);
%Gráfico 1 (versus Tempo)
figure;
[p, ef, t] = plotyy(tempo,estimulof,tempo,taxa,'plot','plot');
ylabel(p(2),'Taxa de Disparos (Hz)','color','red');
ylabel(p(1),'Estimulo Fit');
set(p(2),'ycolor', 'red');
xlabel(p(1),'Tempo (s)');
set(t, 'color', 'red');
%set(p(1),'ylim',[-20 30]);
%Gráfico 2 (versus Estimulo)
figure;
[p2, t2, ef2] = plotyy(estimulo,taxa,estimulo,estimulof,'plot','plot');
ylabel(p2(1),'Taxa de Disparos (Hz)');
ylabel(p2(2),'Estimulo Fit','color','red');
set(p2(2),'ycolor', 'red');
xlabel(p2(1),'Estimulo');
set(ef2, 'color', 'red');
%set(p2(2),'ylim',[-20 30]);
set(t2,'linestyle','none','marker','o');
set(ef2,'LineWidth', 3);
R101 = corrcoef(taxa,estimulo);

% 9.o) Diferentes Núcleos
% NÚCLEO 1001
w = ones(1001,1)/1001;
d = double(disparos); %vetor não pode ser logical
c = conv(d,w);
prob = c(501:end-500);
ajuste = (sum(disparos)/max(tempo))/mean(prob); %Ajuste para Hz
taxa = prob*ajuste;
p = polyfit(taxa,estimulo,1);
estimulof = repmat(p,5000,1);
estimulof(:,1) = estimulo * p(:,1);
%Gráfico 1 (versus Tempo)
figure;
[p, ef, t] = plotyy(tempo,estimulof,tempo,taxa,'plot','plot');
ylabel(p(2),'Taxa de Disparos (Hz)','color','red');
ylabel(p(1),'Estimulo Fit');
set(p(2),'ycolor', 'red');
xlabel(p(1),'Tempo (s)');
set(t, 'color', 'red');
%Gráfico 2 (versus Estimulo)
figure;
[p2, t2, ef2] = plotyy(estimulo,taxa,estimulo,estimulof,'plot','plot');
ylabel(p2(1),'Taxa de Disparos (Hz)');
ylabel(p2(2),'Estimulo Fit','color','red');
set(p2(2),'ycolor', 'red');
xlabel(p2(1),'Estimulo');
set(ef2, 'color', 'red');
set(t2,'linestyle','none','marker','o');
set(ef2,'LineWidth', 3);
R1001 = corrcoef(taxa,estimulo);
% ###################
% NÚCLEO 11
w = ones(11,1)/11;
d = double(disparos); %vetor não pode ser logical
c = conv(d,w);
prob = c(6:end-5);
ajuste = (sum(disparos)/max(tempo))/mean(prob); %Ajuste para Hz
taxa = prob*ajuste;
p = polyfit(taxa,estimulo,1);
estimulof = repmat(p,5000,1);
estimulof(:,1) = estimulo * p(:,1);
%Gráfico 1 (versus Tempo)
figure;
[p, ef, t] = plotyy(tempo,estimulof,tempo,taxa,'plot','plot');
ylabel(p(2),'Taxa de Disparos (Hz)','color','red');
ylabel(p(1),'Estimulo Fit','color','blue');
set(p(2),'ycolor', 'red');
xlabel(p(1),'Tempo (s)');
set(t, 'color', 'red');
%Gráfico 2 (versus Estimulo)
figure;
[p2, t2, ef2] = plotyy(estimulo,taxa,estimulo,estimulof,'plot','plot');
ylabel(p2(1),'Taxa de Disparos (Hz)');
ylabel(p2(2),'Estimulo Fit','color','red');
set(p2(2),'ycolor', 'red');
xlabel(p2(1),'Estimulo');
set(ef2, 'color', 'red');
set(t2,'linestyle','none','marker','o');
set(ef2,'LineWidth', 3);
R11 = corrcoef(taxa,estimulo);

% 9.p) Autocorrelação (Estímulo)
figure;
ac = xcorr(estimulo, 1000);
plot(ac);

% 9.q) Autocorrelação (Trem de Disparos)
figure;
ac = xcorr(disparos, 1000);
plot(ac);

% 9.r) Correlação Cruzada (Estímulo x Trem de Disparos)
figure;
disparos = double(disparos);
ac = xcorr(estimulo, disparos, 1000);
p = plot(ac);
