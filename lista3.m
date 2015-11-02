% Terceira lista de Exercícios - Introdução à Neurociência Computacional
% Aluno: Rafael Delalibera Rodrigues
% Programa de Pós-Graduação em Computação Aplicada [FFCLRP-USP Ribeirão Preto] - Mestrado

% Exercício 1: Modelo de Hodgkin-Huxley
% 1.b) Implementação: Solução do Sistema de Equações

% Método de Euler
[v,t,n,m,h,g_Na,g_k] = HH_Euler(0,0,100,0.025,20,60,1);
% Euler: V x t
figure;
s1 = subplot(2,2,1);
p1 = plot(t,v);
title(s1,'V x t');
% Euler: n x t
s2 = subplot(2,2,2);
p2 = plot(t,n);
title(s2,'n x t');
% Euler: m x t
s3 = subplot(2,2,3);
p3 = plot(t,m);
title(s3,'m x t');
% Euler: h x t
s4 = subplot(2,2,4);
p4 = plot(t,h);
title(s4,'h x t');
set(gcf,'name','Método de Euler','numbertitle','off');


% Método de Range-Kutta de Quarta Ordem
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,100,0.05,20,60,1);
% Range-Kutta: V x t
figure;
s1 = subplot(2,2,1);
p1 = plot(t,v);
title(s1,'V x t');
% Range-Kutta: n x t
s2 = subplot(2,2,2);
p2 = plot(t,n);
title(s2,'n x t');
% Range-Kutta: m x t
s3 = subplot(2,2,3);
p3 = plot(t,m);
title(s3,'m x t');
% Range-Kutta: h x t
s4 = subplot(2,2,4);
p4 = plot(t,h);
title(s4,'h x t');
set(gcf,'name','Método de Range-Kutta','numbertitle','off');


% Exercício 1: Modelo de Hodgkin-Huxley
% 1.c) Simulação: Diferentes Densidades de Corrente

% Método de Range-Kutta de Quarta Ordem
% Densidade de Corrente Despolarizante = 1
figure;
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,105,0.05,5,105,1);
s1 = subplot(2,2,1);
p1 = plot(t,v);
title(s1,'V x t (J = 1)');
set(s1,'xlim',[0 105]);
% Densidade de Corrente Despolarizante = 5
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,105,0.05,5,105,5);
s2 = subplot(2,2,2);
p2 = plot(t,v);
title(s2,'V x t (J = 5)');
set(s2,'xlim',[0 105]);
% Densidade de Corrente Despolarizante = 10
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,105,0.05,5,105,10);
s3 = subplot(2,2,3);
p3 = plot(t,v);
title(s3,'V x t (J = 10)');
set(s3,'xlim',[0 105]);
% Densidade de Corrente Despolarizante = 50
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,105,0.05,5,105,50);
s4 = subplot(2,2,4);
p4 = plot(t,v);
title(s4,'V x t (J = 50)');
set(s4,'xlim',[0 105]);
set(gcf,'name','Range-Kutta (Diferentes Densidades de Corrente)','numbertitle','off');


% Exercício 1: Modelo de Hodgkin-Huxley
% 1.d) Análise: Valores Limiares de Densidade de Corrente

% Comportamento 1 para 2 [J=2,4]

% Encontrar Valor de Transição de 1 para 2
% Queremos o valor onde há um SALTO
delta = 10;
j_ant = 0;
j_pos = 0;
j_mark = 0;
for j=1:0.1:5
	[v,tv,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,105,0.05,5,105,j)
	j_pos = max(v);
	if (abs(j_pos-j_ant)>delta)
		 j_mark = j;
    end
	j_ant = j_pos;
end

% Gráficos
% Anterior
figure;
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,105,0.05,5,105,j_mark-0.1);
s1 = subplot(2,2,1);
p1 = plot(t,v);
title(s1,strcat('V x t (J = ',num2str(j_mark-0.1),')'));
set(s1,'xlim',[0 105]);
% Valor de Transição
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,105,0.05,5,105,j_mark);
s2 = subplot(2,2,2);
p2 = plot(t,v);
title(s2,strcat('V x t (J = ',num2str(j_mark),') VALOR TRANSIÇÃO'));
set(s2,'xlim',[0 105]);
% Valor Posterior
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,105,0.05,5,105,j_mark+0.1);
s3 = subplot(2,2,3);
p3 = plot(t,v);
title(s3,strcat('V x t (J = ',num2str(j_mark+0.1),')'));
set(s3,'xlim',[0 105]);
% Valor Posterior ao posterior
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,105,0.05,5,105,j_mark+0.2);
s4 = subplot(2,2,4);
p4 = plot(t,v);
title(s4,strcat('V x t (J = ',num2str(j_mark+0.2),')'));
set(s4,'xlim',[0 105]);
set(gcf,'name','Comportamento 1 para 2','numbertitle','off');


% Comportamento 2 para 3 [J=6,5]

% Encontrar Valor de Transição de 2 para 3
% Queremos o valor onde há a presença de mais de um spike
cont_spike = 0;
j_mark = 0;
flag = true;
for j=2.3:0.1:10
	[v,tv,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,105,0.05,5,105,j)
	
	for i = 1:length(v)
		if(v(i)>90)
			if(flag == true)
				cont_spike = cont_spike + 1;
				flag = false;
			end
		else
			flag = true;
		end
	end
	
	if(cont_spike>2) % Comportamento único ignorado
		j_mark = j;
		break;
	else
		cont_spike = 0;
	end;
end

% Gráficos
% Anterior
figure;
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,105,0.05,5,105,j_mark-0.1);
s1 = subplot(2,2,1);
p1 = plot(t,v);
title(s1,strcat('V x t (J = ',num2str(j_mark-0.1),')'));
set(s1,'xlim',[0 105]);
% Valor de Transição
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,105,0.05,5,105,j_mark);
s2 = subplot(2,2,2);
p2 = plot(t,v);
title(s2,strcat('V x t (J = ',num2str(j_mark),') VALOR TRANSIÇÃO'));
set(s2,'xlim',[0 105]);
% Valor Posterior
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,105,0.05,5,105,j_mark+0.1);
s3 = subplot(2,2,3);
p3 = plot(t,v);
title(s3,strcat('V x t (J = ',num2str(j_mark+0.1),')'));
set(s3,'xlim',[0 105]);
% Valor Posterior ao posterior
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,105,0.05,5,105,j_mark+0.2);
s4 = subplot(2,2,4);
p4 = plot(t,v);
title(s4,strcat('V x t (J = ',num2str(j_mark+0.2),')'));
set(s4,'xlim',[0 105]);
set(gcf,'name','Comportamento 2 para 3','numbertitle','off');


% Comportamento 3 para 4 [Maior que ~158]

% Encontrar Valor de Transição de 3 para 4
% Estratégia: Encontrar picos e vales da função v x t
% Verificar as diferenças e analisar os resultados
% OBS: Primeiro pico/vale é ignorado, pois não há corrente injetada
fim = 1000
vet_correntes = [6.4:0.5:200];

%valores após 1000 ms
vet_diferencas = zeros(1,length(vet_correntes)); % entre pico e vale
vet_picos = zeros(1,length(vet_correntes));
vet_vales = zeros(1,length(vet_correntes));

i = 1;
for j=6.4:0.5:200
%for j=150:0.1:200
	
	[v,tv,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,fim,0.05,5,fim,j)
	[picos,locs] = findpeaks(v);
	[vales,locs] = findpeaks(-v);
	vales = -vales; %ajuste dos valores
	if(length(picos)>length(vales)) % ajuste das dimensões
		diferenca = picos(1:length(vales))-vales;
	else
		diferenca = picos-vales(1:length(picos));
	end
	
	vet_picos(i) = picos(length(picos));
	vet_vales(i) = vales(length(vales));
	vet_diferencas(i) = diferenca(length(diferenca));
	i = i + 1;
end

%Diferença Pico Vale
figure;
s1 = subplot(2,2,1);
p1 = plot(vet_correntes,vet_diferencas);
xlabel('Corrente Injetada');
ylabel('Amplitude Sinal');
title(s1,'Amplitude x J (Após 1000ms)');
%Picos
s1 = subplot(2,2,2);
p1 = plot(vet_correntes,vet_picos);
xlabel('Corrente Injetada');
ylabel('Pico Sinal');
title(s1,'Pico x J (Após 1000ms)');
%Vales
s1 = subplot(2,2,3);
p1 = plot(vet_correntes,vet_vales);
xlabel('Corrente Injetada');
ylabel('Vale Sinal');
title(s1,'Vale x J (Após 1000ms)');
%Picos&Vales
s1 = subplot(2,2,4);
p1 = plot(vet_correntes,vet_picos,vet_correntes,vet_vales);
xlabel('Corrente Injetada');
ylabel('Sinal');
title(s1,'Vale x J (Após 1000ms)');
legend('Picos','Vales');


% Gráficos TRANSIÇÃO
j_mark = 157.3;
% Anterior
fim = 2000;
ymin = 0;
ymax = 40;
figure;
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,fim,0.05,5,fim,j_mark-10);
s1 = subplot(2,2,1);
p1 = plot(t,v);
title(s1,strcat('V x t (J = ',num2str(j_mark-10),')'));
set(s1,'xlim',[0 fim]);
set(s1,'ylim',[ymin ymax]);
% Valor de Transição
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,fim,0.05,5,fim,j_mark);
s2 = subplot(2,2,2);
p2 = plot(t,v);
title(s2,strcat('V x t (J = ',num2str(j_mark),') VALOR TRANSIÇÃO'));
set(s2,'xlim',[0 fim]);
set(s2,'ylim',[ymin ymax]);
% Valor Posterior
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,fim,0.05,5,fim,j_mark+10);
s3 = subplot(2,2,3);
p3 = plot(t,v);
title(s3,strcat('V x t (J = ',num2str(j_mark+10),')'));
set(s3,'xlim',[0 fim]);
set(s3,'ylim',[ymin ymax]);
% Valor Posterior ao posterior
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,fim,0.05,5,fim,j_mark+20);
s4 = subplot(2,2,4);
p4 = plot(t,v);
title(s4,strcat('V x t (J = ',num2str(j_mark+20),')'));
set(s4,'xlim',[0 fim]);
set(s4,'ylim',[ymin ymax]);
set(gcf,'name','Comportamento 1 para 2','numbertitle','off');


% Exercício 1: Modelo de Hodgkin-Huxley
% 1.e) Análise: Frequência de Disparos VS Corrente Injetada

fim = 1000;
vet_correntes = [6.4:0.5:200];
vet_frequencias = zeros(1,length(vet_correntes));

i = 1;
for j=6.4:0.5:200
	[v,tv,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,fim,0.05,5,fim,j)
	[picos,locs] = findpeaks(v);
	frequencia = (length(picos)*1.0)/(fim*0.001); % freq em Hz
	
	vet_frequencias(i) = frequencia;
	i = i + 1;
end

figure;
plot(vet_correntes, vet_frequencias);
xlabel('Densidade de Corrente');
ylabel('Frequência (Hz)');
xlim([0 160]);

%Frequencia Mínima
fim =1000;
[v,tv,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,fim,0.05,5,fim,6.4)
[picos,locs] = findpeaks(v);
frequencia_minima = (length(picos)*1.0)/(fim*0.001); % freq em Hz

%Frequencia Máxima
fim =1000;
[v,tv,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,fim,0.05,5,fim,157)
[picos,locs] = findpeaks(v);
frequencia_maxima = (length(picos)*1.0)/(fim*0.001); % freq em Hz


% Exercício 1: Modelo de Hodgkin-Huxley
% 1.f) Análise: Correntes Hiperpolarizantes

% Encontrar valor de transição
j_mark = 0;
for j=-0.5:-0.1:-10
	[v,tv,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,50,0.05,5,30,j)
	if(max(v)>90)
		j_mark = j;
		break;
	end
end

% Gráficos
% Anterior
figure;
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,50,0.05,5,30,j_mark-0.1);
s1 = subplot(2,2,1);
p1 = plot(t,v);
title(s1,strcat('V x t (J = ',num2str(j_mark-0.1),')'));
set(s1,'xlim',[0 50]);
% Valor de Transição
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,50,0.05,5,30,j_mark);
s2 = subplot(2,2,2);
p2 = plot(t,v);
title(s2,strcat('V x t (J = ',num2str(j_mark),') VALOR TRANSIÇÃO'));
set(s2,'xlim',[0 50]);
% Valor Posterior
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,50,0.05,5,30,j_mark+0.1);
s3 = subplot(2,2,3);
p3 = plot(t,v);
title(s3,strcat('V x t (J = ',num2str(j_mark+0.1),')'));
set(s3,'xlim',[0 50]);
% Valor Posterior ao posterior
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,50,0.05,5,30,j_mark+0.2);
s4 = subplot(2,2,4);
p4 = plot(t,v);
title(s4,strcat('V x t (J = ',num2str(j_mark+0.2),')'));
set(s4,'xlim',[0 50]);


figure;
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,50,0.05,5,30,1);
% V x t
s1 = subplot(2,2,1);
p1 = plot(t,v);
set(s1,'xlim',[0 50]);
title(s1,'V x t');
% n x t
s2 = subplot(2,2,2);
p2 = plot(t,n);
set(s2,'xlim',[0 50]);
title(s2,'n x t');
% m x t
s3 = subplot(2,2,3);
p3 = plot(t,m);
set(s3,'xlim',[0 50]);
title(s3,'m x t');
% h x t
s4 = subplot(2,2,4);
p4 = plot(t,h);
set(s4,'xlim',[0 50]);
title(s4,'h x t');
set(gcf,'name','Comportamento 1','numbertitle','off');

figure;
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,50,0.05,5,30,5);
% V x t
s1 = subplot(2,2,1);
p1 = plot(t,v);
set(s1,'xlim',[0 50]);
title(s1,'V x t');
% n x t
s2 = subplot(2,2,2);
p2 = plot(t,n);
set(s2,'xlim',[0 50]);
title(s2,'n x t');
% m x t
s3 = subplot(2,2,3);
p3 = plot(t,m);
set(s3,'xlim',[0 50]);
title(s3,'m x t');
% h x t
s4 = subplot(2,2,4);
p4 = plot(t,h);
set(s4,'xlim',[0 50]);
title(s4,'h x t');
set(gcf,'name','Comportamento 2','numbertitle','off');


% Exercício 1: Modelo de Hodgkin-Huxley
% 1.g) Análise: Período Refratário

% Amplitude de pulso de corrente de 1 ms para elicitar 1 spike

% Queremos o valor onde há um SALTO
delta = 10;
j_ant = 0;
j_pos = 0;
j_mark = 0;
for j=1:0.1:10
	[v,tv,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,105,0.05,5,6,j)
	j_pos = max(v);
	if (abs(j_pos-j_ant)>delta)
		 j_mark = j;
    end
	j_ant = j_pos;
end

% J_MARK = 7.8
% Gráficos
% Anterior
figure;
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,105,0.05,5,6,j_mark-0.1);
s1 = subplot(2,2,1);
p1 = plot(t,v);
title(s1,strcat('V x t (J = ',num2str(j_mark-0.1),')'));
set(s1,'xlim',[0 105]);
% Valor de Transição
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,105,0.05,5,6,j_mark);
s2 = subplot(2,2,2);
p2 = plot(t,v);
title(s2,strcat('V x t (J = ',num2str(j_mark),') VALOR TRANSIÇÃO'));
set(s2,'xlim',[0 105]);
% Valor Posterior
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,105,0.05,5,6,j_mark+0.1);
s3 = subplot(2,2,3);
p3 = plot(t,v);
title(s3,strcat('V x t (J = ',num2str(j_mark+0.1),')'));
set(s3,'xlim',[0 105]);
% Valor Posterior ao posterior
[v,t,n,m,h,g_Na,g_k] = HH_RKUTA4(0,0,105,0.05,5,6,j_mark+0.2);
s4 = subplot(2,2,4);
p4 = plot(t,v);
title(s4,strcat('V x t (J = ',num2str(j_mark+0.2),')'));
set(s4,'xlim',[0 105]);


% Latência para gerar mais 1 spike (com pulso de mesma amplitude):

%Encontrar Latência
l_mark = 0;
for L=1:0.1:30
	[v,tv,n,m,h,g_Na,g_k,J] = HH_RKUTA4_AL(0,0,105,0.05,5,6,7.8,7.8,L);
	[p] = findpeaks(v);
	counter = 0;
	for i=1:length(p)
		if(p(i)>90)
			counter = counter + 1;
		end
		if(counter>=2)
			l_mark = L;
			break;
		end
	end
	
	if(l_mark ~= 0)
		break;
	end
end
figure;
plot(tv,v);
xlim([0 105]);
% L = 19.2 ms


% Análise parâmetros A e L:

l_mark = 0;
L = 20;
vet_latencias = zeros(1,length(8:1:50));
vet_correntes = [8:1:50];

for A=8:1:70 %Variação Amplitude de Corrente (acima da mínima para 1 spike)

	for L=6:0.1:20  %Variação tempo de latência (terminando abaixo do ref. obtido acima)
		[v,tv,n,m,h,g_Na,g_k,J] = HH_RKUTA4_AL(0,0,105,0.05,5,6,7.8,A,L);
		[p] = findpeaks(v);
		counter = 0;

		for i=1:length(p)
			if(p(i)>90)
				counter = counter + 1;
			end
			if(counter>=2)
				l_mark = L;
				break;
			end
		end
		
		if(l_mark ~= 0)
			vet_latencias(A-7) = l_mark; %AJUSTADO
			l_mark = 0;
			break;
		end
	end
end

figure;
plot(vet_latencias(1:length(vet_correntes)),vet_correntes,'.');
hold on;  % mantem a figura
ha = area([0 min(vet_latencias)], [max(vet_correntes) max(vet_correntes)],'FaceColor',[.5 .5 .5]);
plot([-5 35],[7.8 7.8],':red');
xlabel('Período Refratário (ms)');
ylabel('Amplitude de Corrente necessária para produzir segundo spike');
xlim([-5 35]);


%anotações
dim = [.25 .15 .2 .2];
str = {'absolute', 'refractory'};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
dim = [.46 .15 .2 .2];
str = {'relative', 'refractory'};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS (Must be defined at the end of the script) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exercício 1: Modelo de Hodgkin-Huxley           [FUNCTIONS]
% 1.a) Subrotinas: Taxas de Abertura e Fechamento [FUNCTIONS]

%%%%%%%%%%%%%
%% ALPHA N %%
%%%%%%%%%%%%%
function tx=alpha_n(V)
    s = 0.01*(10-V);
    i = exp(1-0.1*V) - 1.0;
    %Indeterminação
    if (s==0 && i==0)
        s= -0.01;
        i= -0.1*exp(1-0.1*V);
    end
    tx=s/i;
end

%%%%%%%%%%%%%
%% ALPHA M %%
%%%%%%%%%%%%%
function tx=alpha_m(V)
    s = 0.1*(25-V);
    i = exp(2.5-0.1*V)-1.0;
    %Indeterminacao
    if ( (s==0) && (i==0))
        s= -0.1;
        i= -0.1*exp(2.5-0.1*V);
    end
    tx=s/i;
end

%%%%%%%%%%%%%
%% ALPHA H %%
%%%%%%%%%%%%%
function tx=alpha_h(V)
    tx=0.07*exp(-V/20.0);
end

%%%%%%%%%%%%
%% BETA N %%
%%%%%%%%%%%%
function tx=beta_n(V)
    tx=0.125*exp(-V/80.0);
end

%%%%%%%%%%%%
%% BETA M %%
%%%%%%%%%%%%
function tx=beta_m(V)
    tx=4.0*exp(-V/18.0);
end


%%%%%%%%%%%%
%% BETA H %%
%%%%%%%%%%%%
function tx=beta_h(V)
    s = 1.0;
    i = exp(3.0-0.1*V) + 1;
    tx=s/i;
end


% Exercício 1: Modelo de Hodgkin-Huxley              [FUNCTION]
% 1.b) Implementação: Solução do Sistema de Equações [FUNCTION]

% Método de Euler

function [v,tv,n,m,h,g_Na,g_k] = HH_Euler(v0,ti,tf,Dt,ti_J,tf_J,Jin)
    
	%parametros eletricos da membrana
    Ena  = 115;
    Ek   = -12;
    Evaz = 10.16; 
    GNa  = 120;
    GK   = 36;
    GL   = 0.3; 
    cm   = 1;
    
	%Dimensão do vetores (conforme granularidade temporal)
    length = ceil((tf-ti)/Dt);
    
    %Vetores das Condutâncias 
    g_Na = zeros(length,1);
    g_k  = zeros(length,1);
    g_l  = GL;
    
    %Vetor de Correte Injetada
    J=zeros(length);
    for i=1:length
       if ((i*Dt>ti_J) && (i*Dt<tf_J))
            J(i)=Jin;
       else
            J(i)=0;
       end
    end
    
    %Vetores de Probabilidades (Gates)
    n  = zeros(length,1);
    m  = zeros(length,1);
    h  = zeros(length,1);
   
    %Vetor Voltagens e Tempo
    v    = zeros(length,1);
    tv   = zeros(length,1);
    v(1) = v0+(-v0);
    tv(1)= Dt;
    
    %Probabilidade de Abertura dos Canais (para célular em repouso)
    n(1) = alpha_n(v(1))/(alpha_n(v(1))+beta_n(v(1)));
    m(1) = alpha_m(v(1))/(alpha_m(v(1))+beta_m(v(1)));
    h(1) = alpha_h(v(1))/(alpha_h(v(1))+beta_h(v(1)));
    
    for t=2:length
        % Cálculo: Condutâncias
        g_Na(t-1) = GNa*(m(t-1)^3)*h(t-1);
        g_k(t-1)  = GK*n(t-1)^4;
        % Cálculo: Voltagem por Euler
        v(t)= v(t-1) + (Dt/cm)*(-g_Na(t-1)*(v(t-1)-Ena)-g_k(t-1)*(v(t-1)-Ek)-g_l*(v(t-1)-Evaz)+J(t-1));
        % Cálculo: Probabilidades
        n(t) = n(t-1) + Dt*((1.0-n(t- 1))*alpha_n(v(t-1))-n(t-1)*beta_n(v(t-1))); 
        m(t) = m(t-1) + Dt*((1.0-m(t-1))*alpha_m(v(t-1))-m(t-1)*beta_m(v(t-1)));
        h(t) = h(t-1) + Dt*((1.0-h(t-1))*alpha_h(v(t-1))-h(t-1)*beta_h(v(t-1)));
        %Vetor de tempo
        tv(t)=t*Dt;
    end
    % Arrumando o potencial de repouso (lembrar aqui que as equacoes acima sao definidas em termos dos potencias de repouso)
    v = v+v0;
end

% Método de Range-Kutta de Quarta Ordem

function [v,tv,n,m,h,g_Na,g_k] = HH_RKUTA4(v0,ti,tf,Dt,ti_J,tf_J,Jin)
    
	%parametros eletricos da membrana
    Ena  = 115;
    Ek   = -12;
    Evaz = 10.16; 
    GNa  = 120;
    GK   = 36;
    GL   = 0.3; 
    cm   = 1;
    
    %Dimensão do vetores (conforme granularidade temporal)
    length = ceil((tf-ti)/Dt);
    
    %Vetores das Condutâncias 
    g_Na = zeros(length,1);
    g_k  = zeros(length,1);
    g_l  = GL;
    
    %Vetor de Correte Injetada
    J=zeros(length);
    for i=1:length
       if ((i*Dt>ti_J) && (i*Dt<tf_J))
            J(i)=Jin;
       else
            J(i)=0;
       end
    end
    
    %Vetores de Probabilidades (Gates)
    n  = zeros(length,1);
    m  = zeros(length,1);
    h  = zeros(length,1);
   
    %Vetor Voltagens e Tempo
    v    = zeros(length,1);
    tv   = zeros(length,1);
    v(1) = v0+(-v0);
    tv(1)= Dt;
    
	%Probabilidade de Abertura dos Canais (para célular em repouso)
    n(1) = alpha_n(v(1))/(alpha_n(v(1))+beta_n(v(1)));
    m(1) = alpha_m(v(1))/(alpha_m(v(1))+beta_m(v(1)));
    h(1) = alpha_h(v(1))/(alpha_h(v(1))+beta_h(v(1)));
    
    for t=2:length
        % Cálculo: Condutâncias
        g_Na(t-1) = GNa*(m(t-1)^3)*h(t-1);
        g_k(t-1)  = GK*n(t-1)^4;
        % Cálculo: Voltagem por Runge Kutta 
        f1 = (-g_Na(t-1)*(v(t-1)-Ena)-g_k(t-1)*(v(t-1)-Ek)-g_l*(v(t-1)-Evaz)+J(t-1));
        f2 = (-g_Na(t-1)*(v(t-1)+(Dt/2)*f1-Ena)-g_k(t-1)*(v(t-1)+(Dt/2)*f1-Ek)-g_l*(v(t-1)+(Dt/2)*f1-Evaz)+J(t-1));
        f3 = (-g_Na(t-1)*(v(t-1)+(Dt/2)*f2-Ena)-g_k(t-1)*(v(t-1)+(Dt/2)*f2-Ek)-g_l*(v(t-1)+(Dt/2)*f2-Evaz)+J(t-1));
        f4 = (-g_Na(t-1)*(v(t-1)+(Dt)*f3-Ena)-g_k(t-1)*(v(t-1)+(Dt)*f3-Ek)-g_l*(v(t-1)+(Dt)*f3-Evaz)+J(t-1));
        v(t)= v(t-1) + (Dt/(6*cm))*(f1+2*f2+2*f3+f4);
        % Cálculo: Probabilidade n por Runge-Kutta
        n1=      ((1.0-n(t-1))*alpha_n(v(t-1))-n(t-1)*beta_n(v(t-1)));
        n2=      ((1.0-n(t-1)-(Dt/2)*n1)*alpha_n(v(t-1))-(n(t-1)+(Dt/2)*n1)*beta_n(v(t-1)));
        n3=      ((1.0-n(t-1)-(Dt/2)*n2)*alpha_n(v(t-1))-(n(t-1)+(Dt/2)*n2)*beta_n(v(t-1)));
        n4=      ((1.0-n(t-1)-(Dt)*n3)*alpha_n(v(t-1))-(n(t-1)+(Dt)*n3)*beta_n(v(t-1)));
        n(t) = n(t-1) + (Dt/6)*(n1+2*n2+2*n3+n4);
        % Cálculo: Probabilidade m por Runge-Kutta
        m1=       ((1.0-m(t-1))*alpha_m(v(t-1))-m(t-1)*beta_m(v(t-1)));
        m2=       ((1.0-m(t-1)-(Dt/2)*m1)*alpha_m(v(t-1))-(m(t-1)+(Dt/2)*m1)*beta_m(v(t-1)));
        m3=       ((1.0-m(t-1)-(Dt/2)*m2)*alpha_m(v(t-1))-(m(t-1)+(Dt/2)*m2)*beta_m(v(t-1)));
        m4=       ((1.0-m(t-1)-Dt*m3)*alpha_m(v(t-1))-(m(t-1)+Dt*m3)*beta_m(v(t-1)));
        m(t) = m(t-1) + (Dt/6)*(m1+2*m2+2*m3+m4);
        % Cálculo: Probabilidade n por Runge-Kutta
        h1=       ((1.0-h(t-1))*alpha_h(v(t-1))-(h(t-1))*beta_h(v(t-1)));
        h2=       ((1.0-h(t-1)-(Dt/2)*h1)*alpha_h(v(t-1))-(h(t-1)+(Dt/2)*h1)*beta_h(v(t-1)));
        h3=       ((1.0-h(t-1)-(Dt/2)*h2)*alpha_h(v(t-1))-(h(t-1)+(Dt/2)*h2)*beta_h(v(t-1)));
        h4=       ((1.0-h(t-1)-Dt*h3)*alpha_h(v(t-1))-(h(t-1)+Dt*h3)*beta_h(v(t-1)));
        h(t) = h(t-1) + (Dt/6)*(h1+2*h2+2*h3+h4);
        %Vetor de tempo
        tv(t)=t*Dt;
    end
    % Arrumando o potencial de repouso (lembrar aqui que as equacoes acima sao definidas em termos dos potencias de repouso)
    v = v+v0;
end


% Latência para gerar mais 1 spike (com pulso de mesma amplitude): [FUNCTION]

function [v,tv,n,m,h,g_Na,g_k,J] = HH_RKUTA4_AL(v0,ti,tf,Dt,ti_J,tf_J,Jin,A,L)
    
	%parametros eletricos da membrana
    Ena  = 115;
    Ek   = -12;
    Evaz = 10.16; 
    GNa  = 120;
    GK   = 36;
    GL   = 0.3; 
    cm   = 1;
    
    %Dimensão do vetores (conforme granularidade temporal)
    length = ceil((tf-ti)/Dt);
    
    %Vetores das Condutâncias 
    g_Na = zeros(length,1);
    g_k  = zeros(length,1);
    g_l  = GL;
    
	%Vetor de Correte Injetada
    J=zeros(length);
    for i=1:length
       if ((i*Dt>ti_J) && (i*Dt<tf_J))
            J(i)=Jin;
       else
            J(i)=0;
       end
    end
	
	%Ajuste do Vetor de Corrente Injetada
	%Para parâmetros A e L
	for i=1:length
       if ((i*Dt>tf_J+L)&&(i*Dt<tf_J+L+1))
            J(i)=A;
       end
    end
    
    %Vetores de Probabilidades (Gates)
    n  = zeros(length,1);
    m  = zeros(length,1);
    h  = zeros(length,1);
   
    %Vetor Voltagens e Tempo
    v    = zeros(length,1);
    tv   = zeros(length,1);
    v(1) = v0+(-v0);
    tv(1)= Dt;
    
	%Probabilidade de Abertura dos Canais (para célular em repouso)
    n(1) = alpha_n(v(1))/(alpha_n(v(1))+beta_n(v(1)));
    m(1) = alpha_m(v(1))/(alpha_m(v(1))+beta_m(v(1)));
    h(1) = alpha_h(v(1))/(alpha_h(v(1))+beta_h(v(1)));
    
    for t=2:length
        % Cálculo: Condutâncias
        g_Na(t-1) = GNa*(m(t-1)^3)*h(t-1);
        g_k(t-1)  = GK*n(t-1)^4;
        % Cálculo: Voltagem por Runge Kutta 
        f1 = (-g_Na(t-1)*(v(t-1)-Ena)-g_k(t-1)*(v(t-1)-Ek)-g_l*(v(t-1)-Evaz)+J(t-1));
        f2 = (-g_Na(t-1)*(v(t-1)+(Dt/2)*f1-Ena)-g_k(t-1)*(v(t-1)+(Dt/2)*f1-Ek)-g_l*(v(t-1)+(Dt/2)*f1-Evaz)+J(t-1));
        f3 = (-g_Na(t-1)*(v(t-1)+(Dt/2)*f2-Ena)-g_k(t-1)*(v(t-1)+(Dt/2)*f2-Ek)-g_l*(v(t-1)+(Dt/2)*f2-Evaz)+J(t-1));
        f4 = (-g_Na(t-1)*(v(t-1)+(Dt)*f3-Ena)-g_k(t-1)*(v(t-1)+(Dt)*f3-Ek)-g_l*(v(t-1)+(Dt)*f3-Evaz)+J(t-1));
        v(t)= v(t-1) + (Dt/(6*cm))*(f1+2*f2+2*f3+f4);
        % Cálculo: Probabilidade n por Runge-Kutta
        n1=      ((1.0-n(t-1))*alpha_n(v(t-1))-n(t-1)*beta_n(v(t-1)));
        n2=      ((1.0-n(t-1)-(Dt/2)*n1)*alpha_n(v(t-1))-(n(t-1)+(Dt/2)*n1)*beta_n(v(t-1)));
        n3=      ((1.0-n(t-1)-(Dt/2)*n2)*alpha_n(v(t-1))-(n(t-1)+(Dt/2)*n2)*beta_n(v(t-1)));
        n4=      ((1.0-n(t-1)-(Dt)*n3)*alpha_n(v(t-1))-(n(t-1)+(Dt)*n3)*beta_n(v(t-1)));
        n(t) = n(t-1) + (Dt/6)*(n1+2*n2+2*n3+n4);
        % Cálculo: Probabilidade m por Runge-Kutta
        m1=       ((1.0-m(t-1))*alpha_m(v(t-1))-m(t-1)*beta_m(v(t-1)));
        m2=       ((1.0-m(t-1)-(Dt/2)*m1)*alpha_m(v(t-1))-(m(t-1)+(Dt/2)*m1)*beta_m(v(t-1)));
        m3=       ((1.0-m(t-1)-(Dt/2)*m2)*alpha_m(v(t-1))-(m(t-1)+(Dt/2)*m2)*beta_m(v(t-1)));
        m4=       ((1.0-m(t-1)-Dt*m3)*alpha_m(v(t-1))-(m(t-1)+Dt*m3)*beta_m(v(t-1)));
        m(t) = m(t-1) + (Dt/6)*(m1+2*m2+2*m3+m4);
        % Cálculo: Probabilidade n por Runge-Kutta
        h1=       ((1.0-h(t-1))*alpha_h(v(t-1))-(h(t-1))*beta_h(v(t-1)));
        h2=       ((1.0-h(t-1)-(Dt/2)*h1)*alpha_h(v(t-1))-(h(t-1)+(Dt/2)*h1)*beta_h(v(t-1)));
        h3=       ((1.0-h(t-1)-(Dt/2)*h2)*alpha_h(v(t-1))-(h(t-1)+(Dt/2)*h2)*beta_h(v(t-1)));
        h4=       ((1.0-h(t-1)-Dt*h3)*alpha_h(v(t-1))-(h(t-1)+Dt*h3)*beta_h(v(t-1)));
        h(t) = h(t-1) + (Dt/6)*(h1+2*h2+2*h3+h4);
        %Vetor de tempo
        tv(t)=t*Dt;
    end
    % Arrumando o potencial de repouso (lembrar aqui que as equacoes acima sao definidas em termos dos potencias de repouso)
    v = v+v0;
end