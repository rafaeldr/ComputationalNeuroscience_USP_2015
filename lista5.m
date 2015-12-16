% Quinta lista de Exercícios - Introdução à Neurociência Computacional
% Aluno: Rafael Delalibera Rodrigues
% Programa de Pós-Graduação em Computação Aplicada [FFCLRP-USP Ribeirão Preto] - Mestrado

% Exercício 1: Teoria do Cabo - Atenuação de Voltagem

% Reprodução do Gráfico da Aula 7 - Pág. 14

% Reprodução da Figura da aula 7 (pág. 14):
V0 = 10; %Qualquer
figure;

% Cabo Finito (Extremidade Selada) L=1
L = 1;
Z = [0:0.1:L];
cfs = (V0*cosh(L-Z))/cosh(L);
VV0_cfs = cfs/V0;
plot(Z,VV0_cfs);
hold on;

% Cabo Finito (Extremidade Aberta) L=1
L = 1;
Z = [0:0.1:L];
cfa = (V0*sinh(L-Z))/sinh(L);
VV0_cfa = cfa/V0;
plot(Z,VV0_cfa);

% Cabo Finito (Extremidade Selada) L=0,5
L = 0.5;
Z = [0:0.1:L];
cfs = (V0*cosh(L-Z))/cosh(L);
VV0_cfs = cfs/V0;
plot(Z,VV0_cfs);
hold on;

% Cabo Finito (Extremidade Aberta) L=0,5
L = 0.5;
Z = [0:0.1:L];
cfa = (V0*sinh(L-Z))/sinh(L);
VV0_cfa = cfa/V0;
plot(Z,VV0_cfa);

% Cabo Finito (Extremidade Selada) L=2
L = 2;
Z = [0:0.1:L];
cfs = (V0*cosh(L-Z))/cosh(L);
VV0_cfs = cfs/V0;
plot(Z,VV0_cfs);
hold on;

% Cabo Finito (Extremidade Aberta) L=2
L = 2;
Z = [0:0.1:L];
cfa = (V0*sinh(L-Z))/sinh(L);
VV0_cfa = cfa/V0;
plot(Z,VV0_cfa);


% Cabo Semi-Infinito
Z = [0:0.1:3]; % Usamos 3, pois é o valor limite usado no gráfico da aula
csi = V0*exp(-Z);  % z/lambda = Z
VV0_csi = csi/V0;
plot(Z,VV0_csi);
legend('Extremidade Selada (L = 1)','Extremidade Aberta (L = 1)','Extremidade Selada (L = 0,5)','Extremidade Aberta (L = 0,5)','Extremidade Selada (L = 2)','Extremidade Aberta (L = 2)','Cabo Infinito');
xlabel('Z=z/\lambda');
ylabel('V/V_0');

hold off;


% Exercício 3: Simulação com 10 Compartimentos
% 3.b) Simulação: Único Compartimento

% 3.b.i) Corrente para incremento em 10mV

Dt = 0.01;
for i=0.001:0.00001:0.01
	[v,t,I] = COMP1_RKUTA4(0,0,100,Dt,20,80,i);

    if (v(60/Dt)>=10)
		break;
    end
end
figure;
plot(t,v);
xlabel('t (ms)');
ylabel('V (mV)');
legend(sprintf('Corrente Injetada: %.4f \\muA',i));
title('Compartimento Único');

% 3.b.ii) Resistência de Entrada: Rin

Rin = v(60/Dt)/i


% Exercício 3: Simulação com 10 Compartimentos
% 3.c) Simulação: Dez Compartimento

% 3.c.i) Corrente para incremento em 10mV (Compartimento 1)

Dt = 0.01;
for i=0.00001:0.00001:0.01
	[v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,t,I] = COMP10_RKUTA4(0,0,100,Dt,5,100,i);

    if (v1(80/Dt)>=10)
		break;
    end
end
figure;
plot(t,v1);
xlabel('t (ms)');
ylabel('V (mV)');
legend(sprintf('Corrente Injetada: %.4f \\muA',i));
title('Compartimento 1');

% 3.c.ii) Resistência de Entrada: Rin

Rin = v1(80/Dt)/i

% 3.c.iii) Gráficos do Modelo
% Gráfico A - Variação Temporal das Voltagens dos 10 Compartimentos

Dt = 0.01;
i = 1e-04; %uA
[v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,t,I] = COMP10_RKUTA4(0,0,100,Dt,0,100,i);
figure;
plot(t,v1);
hold on;
plot(t,v2);
plot(t,v3);
plot(t,v4);
plot(t,v5);
plot(t,v6);
plot(t,v7);
plot(t,v8);
plot(t,v9);
plot(t,v10);
xlabel('t (ms)');
ylabel('V (mV)');
legend('Compartimento 1','Compartimento 2','Compartimento 3','Compartimento 4','Compartimento 5','Compartimento 6','Compartimento 7','Compartimento 8','Compartimento 9','Compartimento 10');
title('Cabo de 10 Compartimentos');
axis([0 100 0 15]);
dim = [.2 .8 .2 .1];
str = sprintf('Corrente Injetada: %.4f \\muA',i);
a = annotation('textbox',dim,'String',str,'FitBoxToText','on');
a.Color = 'red';

% Gráfico B - Variação do Potencial Estacionário ao Longo de z

Dt = 0.01;
i = 1e-04; %uA
[v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,t,I] = COMP10_RKUTA4(0,0,100,Dt,0,100,i);
V0 = 10;
L = 1;
Z = [0:0.01:L];
cfs = (V0*cosh(L-Z))/cosh(L);
figure;
plot(Z,cfs);
hold on;

simZ = [0.05:0.1:0.95];
L=0.1;
Z=L/2;
V0 = [v1(80/Dt) v2(80/Dt) v3(80/Dt) v4(80/Dt) v5(80/Dt) v6(80/Dt) v7(80/Dt) v8(80/Dt) v9(80/Dt) v10(80/Dt)];
f = (V0*cosh(L-Z))/cosh(L);
plot(simZ,f,'+');
axis([0 1 6 11]);
legend('Solução Analítica', 'Simulação');
title('Decaimento da voltagem estacionária com a distância eletrotônica');
xlabel('Distância Eletrotônica (Z)');
ylabel('Voltagem estacionária (mV)');


% Exercício 4: Experimento de Rall

tau = 10; %ms
dt = tau/4;
ti_Ia = 0;
tf_Ia = ti_Ia + dt;
ti_Ib = tf_Ia;
tf_Ib = ti_Ib + dt;
ti_Ic = tf_Ib;
tf_Ic = ti_Ic + dt;
ti_Id = tf_Ic;
tf_Id = ti_Id + dt;
i = 1e-04; %uA
[v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,t,Ia,Ib,Ic,Id] = RALL_RKUTA4(0,0,2.5*tau,0.001,ti_Ia,tf_Ia,i,ti_Ib,tf_Ib,i,ti_Ic,tf_Ic,i,ti_Id,tf_Id,i);
figure;
plot(t/tau,v1);
hold on;

ti_Id = 0;
tf_Id = ti_Id + dt;
ti_Ic = tf_Id;
tf_Ic = ti_Ic + dt;
ti_Ib = tf_Ic;
tf_Ib = ti_Ib + dt;
ti_Ia = tf_Ib;
tf_Ia = ti_Ia + dt;
[v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,t,Ia,Ib,Ic,Id] = RALL_RKUTA4(0,0,2.5*tau,0.001,ti_Ia,tf_Ia,i,ti_Ib,tf_Ib,i,ti_Ic,tf_Ic,i,ti_Id,tf_Id,i);
plot(t/tau,v1);
legend('A-B-C-D','D-C-B-A');
%axis([0 2.5 6 11]);
title('Efeito do padrão espaço-temporal do estímulo');
xlabel('t/\tau');
ylabel('Voltagem no Soma (mV) [Compartimento 1]');
y=10+(t/tau)*0;
plot(t/tau,y,'black:');
x = [0.3 0.3];
y = [0.8 0.7];
annotation('textarrow',x,y,'String','Possível Limiar de Disparo (Threshold)');


% Exercício 5: Modelo de 10 Compartimentos Ativos
% 5.a) Modelo no MATLAB

Dt = 0.01;
i = 1e-03; %uA
[v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,t,I] = COMP10ATIV_RKUTA4(0,0,8,Dt,0.5,8,i);
figure;
plot(t,v1);
hold on;
plot(t,v2);
plot(t,v3);
plot(t,v4);
plot(t,v5);
plot(t,v6);
plot(t,v7);
plot(t,v8);
plot(t,v9);
plot(t,v10);
xlabel('t (ms)');
ylabel('V (mV)');
legend('Compartimento 1','Compartimento 2','Compartimento 3','Compartimento 4','Compartimento 5','Compartimento 6','Compartimento 7','Compartimento 8','Compartimento 9','Compartimento 10');
title('Propagação do Potencial de Ação');
dim = [.2 .8 .2 .1];
str = sprintf('Corrente Injetada: %.4f \\muA',i);
a = annotation('textbox',dim,'String',str,'FitBoxToText','on');
a.Color = 'red';


% Exercício 5: Modelo de 10 Compartimentos Ativos
% 5.b) Velocidade de Propagação do Potencial de Ação

%Encontrar os picos do potencial de ação (ponto de referência)
[pico1,loc1] = findpeaks(v1);
[pico2,loc2] = findpeaks(v2);
[pico3,loc3] = findpeaks(v3);
[pico4,loc4] = findpeaks(v4);
[pico5,loc5] = findpeaks(v5);
[pico6,loc6] = findpeaks(v6);
[pico7,loc7] = findpeaks(v7);
[pico8,loc8] = findpeaks(v8);
[pico9,loc9] = findpeaks(v9);
[pico10,loc10] = findpeaks(v10);
%Diferenças de tempo para exibição do potencial em cada compartimento
d1 = loc2-loc1; % Extremidade
d2 = loc3-loc2;
d3 = loc4-loc3;
d4 = loc5-loc4;
d5 = loc6-loc5;
d6 = loc7-loc6;
d7 = loc8-loc7;
d8 = loc9-loc8;
d9 = loc10-loc9; % Extremidade
% As extremidades apresentaram valores pouco menores, porém os compartimentos internos apresentaram o mesmo valor de dt.
% d2 = d3 = d4 = d5 = d6 = d7 = d8
dt = (d5/length(v1))*8;
l = 100; % Comprimento do compartimento (em micrometros)
vprop = l/dt


% Exercício 5: Modelo de 10 Compartimentos Ativos
% 5.c) Análise: Velocidade de Propagação do Potencial de Ação

l = 100; % Comprimento do compartimento (em micrometros)
Dt = 0.01;
i = 1e-03; %uA
% a original = 2 um
% a reduzido 20 vezes (0,1 um)
[v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,t,I] = COMP10ATIVv2_RKUTA4(0,0,8,Dt,0.5,8,i,0.05);
[pico5,loc5] = findpeaks(v5);
[pico6,loc6] = findpeaks(v6);
d5 = loc6-loc5;
dt = (d5/length(v1))*8;
vprop_reduzido = l/dt;
% a aumentado 20 vezes (40 um)
[v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,t,I] = COMP10ATIVv2_RKUTA4(0,0,8,Dt,0.5,8,i,20);
[pico5,loc5] = findpeaks(v5);
[pico6,loc6] = findpeaks(v6);
d5 = loc6-loc5;
dt = (d5/length(v1))*8;
vprop_aumentado = l/dt;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS (Must be defined at the end of the script) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exercício 3: Simulação com 10 Compartimentos [FUNCTION]
% 3.b) Simulação: Único Compartimento          [FUNCTION]

function [v,tv,I] = COMP1_RKUTA4(v0,ti,tf,Dt,ti_I,tf_I,Iin)
    
	%parametros do modelo
    G  = 1.25e-4;
    Gv = 1.25e-6; 
    Cm = 1.25e-5;
    
    %Dimensão do vetores (conforme granularidade temporal)
    len = ceil((tf-ti)/Dt);
    
    %Vetor de Correte Injetada
    I=zeros(len,1);
    for i=1:len
        % Decisao quanto a entrada de corrente
       if ((i*Dt>ti_I) && (i*Dt<tf_I))
            I(i)=Iin;
       else
            I(i)=0;
       end
    end
   
    %Vetor Voltagens e Tempo
    v    = zeros(len,1);
    tv   = zeros(len,1);
    v(1) = v0+(-v0);
    tv(1)= Dt;
    
    for t=2:len
        
        f1 = -2*G*v(t-1)-Gv*v(t-1)+I(t-1);
        f2 = -2*G*(v(t-1)+(Dt/2)*f1)-Gv*(v(t-1)+(Dt/2)*f1)+I(t-1);
        f3 = -2*G*(v(t-1)+(Dt/2)*f2)-Gv*(v(t-1)+(Dt/2)*f2)+I(t-1);
        f4 = -2*G*(v(t-1)+Dt*f3)-Gv*(v(t-1)+Dt*f3)+I(t-1);
        v(t)= v(t-1) + (Dt/(6*Cm))*(f1+2*f2+2*f3+f4);
        
        tv(t)=t*Dt;
    end
end


% Exercício 3: Simulação com 10 Compartimentos [FUNCTION]
% 3.c) Simulação: Dez Compartimento            [FUNCTION]

function [v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,tv,I] = COMP10_RKUTA4(v0,ti,tf,Dt,ti_I,tf_I,Iin)
    
	%parametros do modelo
    G  = 1.25e-4;
    Gv = 1.25e-6;
    Cm = 1.25e-5;
    Ev = 0;
    
    %Dimensão do vetores (conforme granularidade temporal)
    len = ceil((tf-ti)/Dt);
    
    %Vetor de Correte Injetada
    I=zeros(len,1);
    for i=1:len
        % Decisao quanto a entrada de corrente
       if ((i*Dt>ti_I) && (i*Dt<tf_I))
            I(i)=Iin;
       else
            I(i)=0;
       end
    end
   
    %Vetor Voltagens e Tempo
    v1   = zeros(len,1);
    v2   = zeros(len,1);
    v3   = zeros(len,1);
    v4   = zeros(len,1);
    v5   = zeros(len,1);
    v6   = zeros(len,1);
    v7   = zeros(len,1);
    v8   = zeros(len,1);
    v9   = zeros(len,1);
    v10  = zeros(len,1);
    v1(1) = v0;
    v2(1) = 0;
    v3(1) = 0;
    v4(1) = 0;
    v5(1) = 0;
    v6(1) = 0;
    v7(1) = 0;
    v8(1) = 0;
    v9(1) = 0;
    v10(1)= 0;
    
    tv   = zeros(len,1);   
    tv(1)= Dt;
    
    for t=2:len
        
        %Compartimento 1
        fc1_1 = -G*(v1(t-1)-v2(t-1))-Gv*(v1(t-1)-Ev)+I(t-1);
        fc1_2 = -G*(v1(t-1)+(Dt/2)*fc1_1-v2(t-1))-Gv*(v1(t-1)+(Dt/2)*fc1_1-Ev)+I(t-1);
        fc1_3 = -G*(v1(t-1)+(Dt/2)*fc1_2-v2(t-1))-Gv*(v1(t-1)+(Dt/2)*fc1_2-Ev)+I(t-1);
        fc1_4 = -G*(v1(t-1)+Dt*fc1_3-v2(t-1))-Gv*(v1(t-1)+Dt*fc1_3-Ev)+I(t-1);
        v1(t) = v1(t-1) + (Dt/(6*Cm))*(fc1_1+2*fc1_2+2*fc1_3+fc1_4);
        
        %Compartimento 2
        fc2_1 = G*(v1(t-1)-v2(t-1))-G*(v2(t-1)-v3(t-1))-Gv*(v2(t-1)-Ev);
        fc2_2 = G*(v1(t-1)-(v2(t-1)+(Dt/2)*fc2_1))-G*(v2(t-1)+(Dt/2)*fc2_1-v3(t-1))-Gv*(v2(t-1)+(Dt/2)*fc2_1-Ev);
        fc2_3 = G*(v1(t-1)-(v2(t-1)+(Dt/2)*fc2_2))-G*(v2(t-1)+(Dt/2)*fc2_2-v3(t-1))-Gv*(v2(t-1)+(Dt/2)*fc2_2-Ev);
        fc2_4 = G*(v1(t-1)-(v2(t-1)+Dt*fc2_3))-G*(v2(t-1)+Dt*fc2_3-v3(t-1))-Gv*(v2(t-1)+Dt*fc2_3-Ev);
        v2(t) = v2(t-1) + (Dt/(6*Cm))*(fc2_1+2*fc2_2+2*fc2_3+fc2_4);
        
        %Compartimento 3
        fc3_1 = G*(v2(t-1)-v3(t-1))-G*(v3(t-1)-v4(t-1))-Gv*(v3(t-1)-Ev);
        fc3_2 = G*(v2(t-1)-(v3(t-1)+(Dt/2)*fc3_1))-G*(v3(t-1)+(Dt/2)*fc3_1-v4(t-1))-Gv*(v3(t-1)+(Dt/2)*fc3_1-Ev);
        fc3_3 = G*(v2(t-1)-(v3(t-1)+(Dt/2)*fc3_2))-G*(v3(t-1)+(Dt/2)*fc3_2-v4(t-1))-Gv*(v3(t-1)+(Dt/2)*fc3_2-Ev);
        fc3_4 = G*(v2(t-1)-(v3(t-1)+Dt*fc3_3))-G*(v3(t-1)+Dt*fc3_3-v4(t-1))-Gv*(v3(t-1)+Dt*fc3_3-Ev);
        v3(t) = v3(t-1) + (Dt/(6*Cm))*(fc3_1+2*fc3_2+2*fc3_3+fc3_4);
        
        %Compartimento 4
        fc4_1 = G*(v3(t-1)-v4(t-1))-G*(v4(t-1)-v5(t-1))-Gv*(v4(t-1)-Ev);
        fc4_2 = G*(v3(t-1)-(v4(t-1)+(Dt/2)*fc4_1))-G*(v4(t-1)+(Dt/2)*fc4_1-v5(t-1))-Gv*(v4(t-1)+(Dt/2)*fc4_1-Ev);
        fc4_3 = G*(v3(t-1)-(v4(t-1)+(Dt/2)*fc4_2))-G*(v4(t-1)+(Dt/2)*fc4_2-v5(t-1))-Gv*(v4(t-1)+(Dt/2)*fc4_2-Ev);
        fc4_4 = G*(v3(t-1)-(v4(t-1)+Dt*fc4_3))-G*(v4(t-1)+Dt*fc4_3-v5(t-1))-Gv*(v4(t-1)+Dt*fc4_3-Ev);
        v4(t) = v4(t-1) + (Dt/(6*Cm))*(fc4_1+2*fc4_2+2*fc4_3+fc4_4);
        
        %Compartimento 5
        fc5_1 = G*(v4(t-1)-v5(t-1))-G*(v5(t-1)-v6(t-1))-Gv*(v5(t-1)-Ev);
        fc5_2 = G*(v4(t-1)-(v5(t-1)+(Dt/2)*fc5_1))-G*(v5(t-1)+(Dt/2)*fc5_1-v6(t-1))-Gv*(v5(t-1)+(Dt/2)*fc5_1-Ev);
        fc5_3 = G*(v4(t-1)-(v5(t-1)+(Dt/2)*fc5_2))-G*(v5(t-1)+(Dt/2)*fc5_2-v6(t-1))-Gv*(v5(t-1)+(Dt/2)*fc5_2-Ev);
        fc5_4 = G*(v4(t-1)-(v5(t-1)+Dt*fc5_3))-G*(v5(t-1)+Dt*fc5_3-v6(t-1))-Gv*(v5(t-1)+Dt*fc5_3-Ev);
        v5(t) = v5(t-1) + (Dt/(6*Cm))*(fc5_1+2*fc5_2+2*fc5_3+fc5_4);
        
        %Compartimento 6
        fc6_1 = G*(v5(t-1)-v6(t-1))-G*(v6(t-1)-v7(t-1))-Gv*(v6(t-1)-Ev);
        fc6_2 = G*(v5(t-1)-(v6(t-1)+(Dt/2)*fc6_1))-G*(v6(t-1)+(Dt/2)*fc6_1-v7(t-1))-Gv*(v6(t-1)+(Dt/2)*fc6_1-Ev);
        fc6_3 = G*(v5(t-1)-(v6(t-1)+(Dt/2)*fc6_2))-G*(v6(t-1)+(Dt/2)*fc6_2-v7(t-1))-Gv*(v6(t-1)+(Dt/2)*fc6_2-Ev);
        fc6_4 = G*(v5(t-1)-(v6(t-1)+Dt*fc6_3))-G*(v6(t-1)+Dt*fc6_3-v7(t-1))-Gv*(v6(t-1)+Dt*fc6_3-Ev);
        v6(t) = v6(t-1) + (Dt/(6*Cm))*(fc6_1+2*fc6_2+2*fc6_3+fc6_4);
        
        %Compartimento 7
        fc7_1 = G*(v6(t-1)-v7(t-1))-G*(v7(t-1)-v8(t-1))-Gv*(v7(t-1)-Ev);
        fc7_2 = G*(v6(t-1)-(v7(t-1)+(Dt/2)*fc7_1))-G*(v7(t-1)+(Dt/2)*fc7_1-v8(t-1))-Gv*(v7(t-1)+(Dt/2)*fc7_1-Ev);
        fc7_3 = G*(v6(t-1)-(v7(t-1)+(Dt/2)*fc7_2))-G*(v7(t-1)+(Dt/2)*fc7_2-v8(t-1))-Gv*(v7(t-1)+(Dt/2)*fc7_2-Ev);
        fc7_4 = G*(v6(t-1)-(v7(t-1)+Dt*fc7_3))-G*(v7(t-1)+Dt*fc7_3-v8(t-1))-Gv*(v7(t-1)+Dt*fc7_3-Ev);
        v7(t) = v7(t-1) + (Dt/(6*Cm))*(fc7_1+2*fc7_2+2*fc7_3+fc7_4);
        
        %Compartimento 8
        fc8_1 = G*(v7(t-1)-v8(t-1))-G*(v8(t-1)-v9(t-1))-Gv*(v8(t-1)-Ev);
        fc8_2 = G*(v7(t-1)-(v8(t-1)+(Dt/2)*fc8_1))-G*(v8(t-1)+(Dt/2)*fc8_1-v9(t-1))-Gv*(v8(t-1)+(Dt/2)*fc8_1-Ev);
        fc8_3 = G*(v7(t-1)-(v8(t-1)+(Dt/2)*fc8_2))-G*(v8(t-1)+(Dt/2)*fc8_2-v9(t-1))-Gv*(v8(t-1)+(Dt/2)*fc8_2-Ev);
        fc8_4 = G*(v7(t-1)-(v8(t-1)+Dt*fc8_3))-G*(v8(t-1)+Dt*fc8_3-v9(t-1))-Gv*(v8(t-1)+Dt*fc8_3-Ev);
        v8(t) = v8(t-1) + (Dt/(6*Cm))*(fc8_1+2*fc8_2+2*fc8_3+fc8_4);
        
        %Compartimento 9
        fc9_1 = G*(v8(t-1)-v9(t-1))-G*(v9(t-1)-v10(t-1))-Gv*(v9(t-1)-Ev);
        fc9_2 = G*(v8(t-1)-(v9(t-1)+(Dt/2)*fc9_1))-G*(v9(t-1)+(Dt/2)*fc9_1-v10(t-1))-Gv*(v9(t-1)+(Dt/2)*fc9_1-Ev);
        fc9_3 = G*(v8(t-1)-(v9(t-1)+(Dt/2)*fc9_2))-G*(v9(t-1)+(Dt/2)*fc9_2-v10(t-1))-Gv*(v9(t-1)+(Dt/2)*fc9_2-Ev);
        fc9_4 = G*(v8(t-1)-(v9(t-1)+Dt*fc9_3))-G*(v9(t-1)+Dt*fc9_3-v10(t-1))-Gv*(v9(t-1)+Dt*fc9_3-Ev);
        v9(t) = v9(t-1) + (Dt/(6*Cm))*(fc9_1+2*fc9_2+2*fc9_3+fc9_4);
        
        %Compartimento 10
        fc10_1 = G*(v9(t-1)-v10(t-1))-Gv*(v10(t-1)-Ev);
        fc10_2 = G*(v9(t-1)-(v10(t-1)+(Dt/2)*fc10_1))-Gv*(v10(t-1)+(Dt/2)*fc10_1-Ev);
        fc10_3 = G*(v9(t-1)-(v10(t-1)+(Dt/2)*fc10_2))-Gv*(v10(t-1)+(Dt/2)*fc10_2-Ev);
        fc10_4 = G*(v9(t-1)-(v10(t-1)+Dt*fc10_3))-Gv*(v10(t-1)+Dt*fc10_3-Ev);
        v10(t) = v10(t-1) + (Dt/(6*Cm))*(fc10_1+2*fc10_2+2*fc10_3+fc10_4);
        
        tv(t)=t*Dt;
    end
end


% Exercício 4: Experimento de Rall [FUNCTION]

function [v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,tv,Ia,Ib,Ic,Id] = RALL_RKUTA4(v0,ti,tf,Dt,ti_Ia,tf_Ia,Iain,ti_Ib,tf_Ib,Ibin,ti_Ic,tf_Ic,Icin,ti_Id,tf_Id,Idin)
    
	%parametros do modelo
    G  = 1.25e-4;
    Gv = 1.25e-6;
    Cm = 1.25e-5;
    Ev = 0;
    
    %Dimensão do vetores (conforme granularidade temporal)
    len = ceil((tf-ti)/Dt);
    
    %Vetor de Correte Injetada
    Ia=zeros(len,1);
    for i=1:len
        % Decisao quanto a entrada de corrente
       if ((i*Dt>ti_Ia) && (i*Dt<tf_Ia))
            Ia(i)=Iain;
       else
            Ia(i)=0;
       end
    end
    
    Ib=zeros(len,1);
    for i=1:len
        % Decisao quanto a entrada de corrente
       if ((i*Dt>ti_Ib) && (i*Dt<tf_Ib))
            Ib(i)=Ibin;
       else
            Ib(i)=0;
       end
    end
    
    Ic=zeros(len,1);
    for i=1:len
        % Decisao quanto a entrada de corrente
       if ((i*Dt>ti_Ic) && (i*Dt<tf_Ic))
            Ic(i)=Icin;
       else
            Ic(i)=0;
       end
    end
    
    Id=zeros(len,1);
    for i=1:len
        % Decisao quanto a entrada de corrente
       if ((i*Dt>ti_Id) && (i*Dt<tf_Id))
            Id(i)=Idin;
       else
            Id(i)=0;
       end
    end
   
    %Vetor Voltagens e Tempo
    v1   = zeros(len,1);
    v2   = zeros(len,1);
    v3   = zeros(len,1);
    v4   = zeros(len,1);
    v5   = zeros(len,1);
    v6   = zeros(len,1);
    v7   = zeros(len,1);
    v8   = zeros(len,1);
    v9   = zeros(len,1);
    v10  = zeros(len,1);
    v1(1) = v0;
    v2(1) = 0;
    v3(1) = 0;
    v4(1) = 0;
    v5(1) = 0;
    v6(1) = 0;
    v7(1) = 0;
    v8(1) = 0;
    v9(1) = 0;
    v10(1)= 0;
    
    tv   = zeros(len,1);   
    tv(1)= Dt;
    
    for t=2:len
        
        %Compartimento 1
        fc1_1 = -G*(v1(t-1)-v2(t-1))-Gv*(v1(t-1)-Ev);
        fc1_2 = -G*(v1(t-1)+(Dt/2)*fc1_1-v2(t-1))-Gv*(v1(t-1)+(Dt/2)*fc1_1-Ev);
        fc1_3 = -G*(v1(t-1)+(Dt/2)*fc1_2-v2(t-1))-Gv*(v1(t-1)+(Dt/2)*fc1_2-Ev);
        fc1_4 = -G*(v1(t-1)+Dt*fc1_3-v2(t-1))-Gv*(v1(t-1)+Dt*fc1_3-Ev);
        v1(t) = v1(t-1) + (Dt/(6*Cm))*(fc1_1+2*fc1_2+2*fc1_3+fc1_4);
        
        %Compartimento 2
        fc2_1 = G*(v1(t-1)-v2(t-1))-G*(v2(t-1)-v3(t-1))-Gv*(v2(t-1)-Ev)+Ia(t-1);
        fc2_2 = G*(v1(t-1)-(v2(t-1)+(Dt/2)*fc2_1))-G*(v2(t-1)+(Dt/2)*fc2_1-v3(t-1))-Gv*(v2(t-1)+(Dt/2)*fc2_1-Ev)+Ia(t-1);
        fc2_3 = G*(v1(t-1)-(v2(t-1)+(Dt/2)*fc2_2))-G*(v2(t-1)+(Dt/2)*fc2_2-v3(t-1))-Gv*(v2(t-1)+(Dt/2)*fc2_2-Ev)+Ia(t-1);
        fc2_4 = G*(v1(t-1)-(v2(t-1)+Dt*fc2_3))-G*(v2(t-1)+Dt*fc2_3-v3(t-1))-Gv*(v2(t-1)+Dt*fc2_3-Ev)+Ia(t-1);
        v2(t) = v2(t-1) + (Dt/(6*Cm))*(fc2_1+2*fc2_2+2*fc2_3+fc2_4);
        
        %Compartimento 3
        fc3_1 = G*(v2(t-1)-v3(t-1))-G*(v3(t-1)-v4(t-1))-Gv*(v3(t-1)-Ev)+Ia(t-1);
        fc3_2 = G*(v2(t-1)-(v3(t-1)+(Dt/2)*fc3_1))-G*(v3(t-1)+(Dt/2)*fc3_1-v4(t-1))-Gv*(v3(t-1)+(Dt/2)*fc3_1-Ev)+Ia(t-1);
        fc3_3 = G*(v2(t-1)-(v3(t-1)+(Dt/2)*fc3_2))-G*(v3(t-1)+(Dt/2)*fc3_2-v4(t-1))-Gv*(v3(t-1)+(Dt/2)*fc3_2-Ev)+Ia(t-1);
        fc3_4 = G*(v2(t-1)-(v3(t-1)+Dt*fc3_3))-G*(v3(t-1)+Dt*fc3_3-v4(t-1))-Gv*(v3(t-1)+Dt*fc3_3-Ev)+Ia(t-1);
        v3(t) = v3(t-1) + (Dt/(6*Cm))*(fc3_1+2*fc3_2+2*fc3_3+fc3_4);
        
        %Compartimento 4
        fc4_1 = G*(v3(t-1)-v4(t-1))-G*(v4(t-1)-v5(t-1))-Gv*(v4(t-1)-Ev)+Ib(t-1);
        fc4_2 = G*(v3(t-1)-(v4(t-1)+(Dt/2)*fc4_1))-G*(v4(t-1)+(Dt/2)*fc4_1-v5(t-1))-Gv*(v4(t-1)+(Dt/2)*fc4_1-Ev)+Ib(t-1);
        fc4_3 = G*(v3(t-1)-(v4(t-1)+(Dt/2)*fc4_2))-G*(v4(t-1)+(Dt/2)*fc4_2-v5(t-1))-Gv*(v4(t-1)+(Dt/2)*fc4_2-Ev)+Ib(t-1);
        fc4_4 = G*(v3(t-1)-(v4(t-1)+Dt*fc4_3))-G*(v4(t-1)+Dt*fc4_3-v5(t-1))-Gv*(v4(t-1)+Dt*fc4_3-Ev)+Ib(t-1);
        v4(t) = v4(t-1) + (Dt/(6*Cm))*(fc4_1+2*fc4_2+2*fc4_3+fc4_4);
        
        %Compartimento 5
        fc5_1 = G*(v4(t-1)-v5(t-1))-G*(v5(t-1)-v6(t-1))-Gv*(v5(t-1)-Ev)+Ib(t-1);
        fc5_2 = G*(v4(t-1)-(v5(t-1)+(Dt/2)*fc5_1))-G*(v5(t-1)+(Dt/2)*fc5_1-v6(t-1))-Gv*(v5(t-1)+(Dt/2)*fc5_1-Ev)+Ib(t-1);
        fc5_3 = G*(v4(t-1)-(v5(t-1)+(Dt/2)*fc5_2))-G*(v5(t-1)+(Dt/2)*fc5_2-v6(t-1))-Gv*(v5(t-1)+(Dt/2)*fc5_2-Ev)+Ib(t-1);
        fc5_4 = G*(v4(t-1)-(v5(t-1)+Dt*fc5_3))-G*(v5(t-1)+Dt*fc5_3-v6(t-1))-Gv*(v5(t-1)+Dt*fc5_3-Ev)+Ib(t-1);
        v5(t) = v5(t-1) + (Dt/(6*Cm))*(fc5_1+2*fc5_2+2*fc5_3+fc5_4);
        
        %Compartimento 6
        fc6_1 = G*(v5(t-1)-v6(t-1))-G*(v6(t-1)-v7(t-1))-Gv*(v6(t-1)-Ev)+Ic(t-1);
        fc6_2 = G*(v5(t-1)-(v6(t-1)+(Dt/2)*fc6_1))-G*(v6(t-1)+(Dt/2)*fc6_1-v7(t-1))-Gv*(v6(t-1)+(Dt/2)*fc6_1-Ev)+Ic(t-1);
        fc6_3 = G*(v5(t-1)-(v6(t-1)+(Dt/2)*fc6_2))-G*(v6(t-1)+(Dt/2)*fc6_2-v7(t-1))-Gv*(v6(t-1)+(Dt/2)*fc6_2-Ev)+Ic(t-1);
        fc6_4 = G*(v5(t-1)-(v6(t-1)+Dt*fc6_3))-G*(v6(t-1)+Dt*fc6_3-v7(t-1))-Gv*(v6(t-1)+Dt*fc6_3-Ev)+Ic(t-1);
        v6(t) = v6(t-1) + (Dt/(6*Cm))*(fc6_1+2*fc6_2+2*fc6_3+fc6_4);
        
        %Compartimento 7
        fc7_1 = G*(v6(t-1)-v7(t-1))-G*(v7(t-1)-v8(t-1))-Gv*(v7(t-1)-Ev)+Ic(t-1);
        fc7_2 = G*(v6(t-1)-(v7(t-1)+(Dt/2)*fc7_1))-G*(v7(t-1)+(Dt/2)*fc7_1-v8(t-1))-Gv*(v7(t-1)+(Dt/2)*fc7_1-Ev)+Ic(t-1);
        fc7_3 = G*(v6(t-1)-(v7(t-1)+(Dt/2)*fc7_2))-G*(v7(t-1)+(Dt/2)*fc7_2-v8(t-1))-Gv*(v7(t-1)+(Dt/2)*fc7_2-Ev)+Ic(t-1);
        fc7_4 = G*(v6(t-1)-(v7(t-1)+Dt*fc7_3))-G*(v7(t-1)+Dt*fc7_3-v8(t-1))-Gv*(v7(t-1)+Dt*fc7_3-Ev)+Ic(t-1);
        v7(t) = v7(t-1) + (Dt/(6*Cm))*(fc7_1+2*fc7_2+2*fc7_3+fc7_4);
        
        %Compartimento 8
        fc8_1 = G*(v7(t-1)-v8(t-1))-G*(v8(t-1)-v9(t-1))-Gv*(v8(t-1)-Ev)+Id(t-1);
        fc8_2 = G*(v7(t-1)-(v8(t-1)+(Dt/2)*fc8_1))-G*(v8(t-1)+(Dt/2)*fc8_1-v9(t-1))-Gv*(v8(t-1)+(Dt/2)*fc8_1-Ev)+Id(t-1);
        fc8_3 = G*(v7(t-1)-(v8(t-1)+(Dt/2)*fc8_2))-G*(v8(t-1)+(Dt/2)*fc8_2-v9(t-1))-Gv*(v8(t-1)+(Dt/2)*fc8_2-Ev)+Id(t-1);
        fc8_4 = G*(v7(t-1)-(v8(t-1)+Dt*fc8_3))-G*(v8(t-1)+Dt*fc8_3-v9(t-1))-Gv*(v8(t-1)+Dt*fc8_3-Ev)+Id(t-1);
        v8(t) = v8(t-1) + (Dt/(6*Cm))*(fc8_1+2*fc8_2+2*fc8_3+fc8_4);
        
        %Compartimento 9
        fc9_1 = G*(v8(t-1)-v9(t-1))-G*(v9(t-1)-v10(t-1))-Gv*(v9(t-1)-Ev)+Id(t-1);
        fc9_2 = G*(v8(t-1)-(v9(t-1)+(Dt/2)*fc9_1))-G*(v9(t-1)+(Dt/2)*fc9_1-v10(t-1))-Gv*(v9(t-1)+(Dt/2)*fc9_1-Ev)+Id(t-1);
        fc9_3 = G*(v8(t-1)-(v9(t-1)+(Dt/2)*fc9_2))-G*(v9(t-1)+(Dt/2)*fc9_2-v10(t-1))-Gv*(v9(t-1)+(Dt/2)*fc9_2-Ev)+Id(t-1);
        fc9_4 = G*(v8(t-1)-(v9(t-1)+Dt*fc9_3))-G*(v9(t-1)+Dt*fc9_3-v10(t-1))-Gv*(v9(t-1)+Dt*fc9_3-Ev)+Id(t-1);
        v9(t) = v9(t-1) + (Dt/(6*Cm))*(fc9_1+2*fc9_2+2*fc9_3+fc9_4);
        
        %Compartimento 10
        fc10_1 = G*(v9(t-1)-v10(t-1))-Gv*(v10(t-1)-Ev);
        fc10_2 = G*(v9(t-1)-(v10(t-1)+(Dt/2)*fc10_1))-Gv*(v10(t-1)+(Dt/2)*fc10_1-Ev);
        fc10_3 = G*(v9(t-1)-(v10(t-1)+(Dt/2)*fc10_2))-Gv*(v10(t-1)+(Dt/2)*fc10_2-Ev);
        fc10_4 = G*(v9(t-1)-(v10(t-1)+Dt*fc10_3))-Gv*(v10(t-1)+Dt*fc10_3-Ev);
        v10(t) = v10(t-1) + (Dt/(6*Cm))*(fc10_1+2*fc10_2+2*fc10_3+fc10_4);
        
        tv(t)=t*Dt;
    end
end


% Exercício 5: Modelo de 10 Compartimentos Ativos [FUNCTION]
% 5.a) Modelo no MATLAB                           [FUNCTION]

function [v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,tv,I] = COMP10ATIV_RKUTA4(v0,ti,tf,Dt,ti_I,tf_I,Iin)
    
	%parametros do modelo
    A  = 1.25e-5;
    G  = 1.25e-4;
    Cm = 1.25e-5;    
    
    %parametros do modelo ativo
    Ena  = 115;
    Ek   = -12;
    Evaz = 10.163; 
    gNa  = 120;
    gK   = 36;
    gV   = 0.3;
    
    %Dimensão do vetores (conforme granularidade temporal)
    len = ceil((tf-ti)/Dt);
    
    %Vetores das Condutâncias 
    g_Na1 = zeros(len,1);
    g_k1  = zeros(len,1);
    g_Na2 = zeros(len,1);
    g_k2  = zeros(len,1);
    g_Na3 = zeros(len,1);
    g_k3  = zeros(len,1);
    g_Na4 = zeros(len,1);
    g_k4  = zeros(len,1);
    g_Na5 = zeros(len,1);
    g_k5  = zeros(len,1);
    g_Na6 = zeros(len,1);
    g_k6  = zeros(len,1);
    g_Na7 = zeros(len,1);
    g_k7  = zeros(len,1);
    g_Na8 = zeros(len,1);
    g_k8  = zeros(len,1);
    g_Na9 = zeros(len,1);
    g_k9  = zeros(len,1);
    g_Na10= zeros(len,1);
    g_k10 = zeros(len,1);
    g_l  = gV;
    
    %Vetor de Correte Injetada
    I=zeros(len,1);
    for i=1:len
        % Decisao quanto a entrada de corrente
       if ((i*Dt>ti_I) && (i*Dt<tf_I))
            I(i)=Iin;
       else
            I(i)=0;
       end
    end
   
    %Vetores de Probabilidades (Gates)
    n1  = zeros(len,1);
    m1  = zeros(len,1);
    h1  = zeros(len,1);
    n2  = zeros(len,1);
    m2  = zeros(len,1);
    h2  = zeros(len,1);
    n3  = zeros(len,1);
    m3  = zeros(len,1);
    h3  = zeros(len,1);
    n4  = zeros(len,1);
    m4  = zeros(len,1);
    h4  = zeros(len,1);
    n5  = zeros(len,1);
    m5  = zeros(len,1);
    h5  = zeros(len,1);
    n6  = zeros(len,1);
    m6  = zeros(len,1);
    h6  = zeros(len,1);
    n7  = zeros(len,1);
    m7  = zeros(len,1);
    h7  = zeros(len,1);
    n8  = zeros(len,1);
    m8  = zeros(len,1);
    h8  = zeros(len,1);
    n9  = zeros(len,1);
    m9  = zeros(len,1);
    h9  = zeros(len,1);
    n10 = zeros(len,1);
    m10 = zeros(len,1);
    h10 = zeros(len,1);
    
    %Vetor Voltagens e Tempo
    v1   = zeros(len,1);
    v2   = zeros(len,1);
    v3   = zeros(len,1);
    v4   = zeros(len,1);
    v5   = zeros(len,1);
    v6   = zeros(len,1);
    v7   = zeros(len,1);
    v8   = zeros(len,1);
    v9   = zeros(len,1);
    v10  = zeros(len,1);
    v1(1) = v0;
    v2(1) = 0;
    v3(1) = 0;
    v4(1) = 0;
    v5(1) = 0;
    v6(1) = 0;
    v7(1) = 0;
    v8(1) = 0;
    v9(1) = 0;
    v10(1)= 0;
    
    tv   = zeros(len,1);   
    tv(1)= Dt;
    
    %Probabilidade de Abertura dos Canais (para célula em repouso)
    n1(1) = alpha_n(v1(1))/(alpha_n(v1(1))+beta_n(v1(1)));
    m1(1) = alpha_m(v1(1))/(alpha_m(v1(1))+beta_m(v1(1)));
    h1(1) = alpha_h(v1(1))/(alpha_h(v1(1))+beta_h(v1(1)));
    n2(1) = alpha_n(v2(1))/(alpha_n(v2(1))+beta_n(v2(1)));
    m2(1) = alpha_m(v2(1))/(alpha_m(v2(1))+beta_m(v2(1)));
    h2(1) = alpha_h(v2(1))/(alpha_h(v2(1))+beta_h(v2(1)));
    n3(1) = alpha_n(v3(1))/(alpha_n(v3(1))+beta_n(v3(1)));
    m3(1) = alpha_m(v3(1))/(alpha_m(v3(1))+beta_m(v3(1)));
    h3(1) = alpha_h(v3(1))/(alpha_h(v3(1))+beta_h(v3(1)));
    n4(1) = alpha_n(v4(1))/(alpha_n(v4(1))+beta_n(v4(1)));
    m4(1) = alpha_m(v4(1))/(alpha_m(v4(1))+beta_m(v4(1)));
    h4(1) = alpha_h(v4(1))/(alpha_h(v4(1))+beta_h(v4(1)));
    n5(1) = alpha_n(v5(1))/(alpha_n(v5(1))+beta_n(v5(1)));
    m5(1) = alpha_m(v5(1))/(alpha_m(v5(1))+beta_m(v5(1)));
    h5(1) = alpha_h(v5(1))/(alpha_h(v5(1))+beta_h(v5(1)));
    n6(1) = alpha_n(v6(1))/(alpha_n(v6(1))+beta_n(v6(1)));
    m6(1) = alpha_m(v6(1))/(alpha_m(v6(1))+beta_m(v6(1)));
    h6(1) = alpha_h(v6(1))/(alpha_h(v6(1))+beta_h(v6(1)));
    n7(1) = alpha_n(v7(1))/(alpha_n(v7(1))+beta_n(v7(1)));
    m7(1) = alpha_m(v7(1))/(alpha_m(v7(1))+beta_m(v7(1)));
    h7(1) = alpha_h(v7(1))/(alpha_h(v7(1))+beta_h(v7(1)));
    n8(1) = alpha_n(v8(1))/(alpha_n(v8(1))+beta_n(v8(1)));
    m8(1) = alpha_m(v8(1))/(alpha_m(v8(1))+beta_m(v8(1)));
    h8(1) = alpha_h(v8(1))/(alpha_h(v8(1))+beta_h(v8(1)));
    n9(1) = alpha_n(v9(1))/(alpha_n(v9(1))+beta_n(v9(1)));
    m9(1) = alpha_m(v9(1))/(alpha_m(v9(1))+beta_m(v9(1)));
    h9(1) = alpha_h(v9(1))/(alpha_h(v9(1))+beta_h(v9(1)));
    n10(1) = alpha_n(v10(1))/(alpha_n(v10(1))+beta_n(v10(1)));
    m10(1) = alpha_m(v10(1))/(alpha_m(v10(1))+beta_m(v10(1)));
    h10(1) = alpha_h(v10(1))/(alpha_h(v10(1))+beta_h(v10(1)));
    
    for t=2:len
        % Cálculo: Condutâncias
        g_Na1(t-1) = gNa*(m1(t-1)^3)*h1(t-1);
        g_k1(t-1)  = gK*n1(t-1)^4;
        g_Na2(t-1) = gNa*(m2(t-1)^3)*h2(t-1);
        g_k2(t-1)  = gK*n2(t-1)^4;
        g_Na3(t-1) = gNa*(m3(t-1)^3)*h3(t-1);
        g_k3(t-1)  = gK*n3(t-1)^4;
        g_Na4(t-1) = gNa*(m4(t-1)^3)*h4(t-1);
        g_k4(t-1)  = gK*n4(t-1)^4;
        g_Na5(t-1) = gNa*(m5(t-1)^3)*h5(t-1);
        g_k5(t-1)  = gK*n5(t-1)^4;
        g_Na6(t-1) = gNa*(m6(t-1)^3)*h6(t-1);
        g_k6(t-1)  = gK*n6(t-1)^4;
        g_Na7(t-1) = gNa*(m7(t-1)^3)*h7(t-1);
        g_k7(t-1)  = gK*n7(t-1)^4;
        g_Na8(t-1) = gNa*(m8(t-1)^3)*h8(t-1);
        g_k8(t-1)  = gK*n8(t-1)^4;
        g_Na9(t-1) = gNa*(m9(t-1)^3)*h9(t-1);
        g_k9(t-1)  = gK*n9(t-1)^4;
        g_Na10(t-1) = gNa*(m10(t-1)^3)*h10(t-1);
        g_k10(t-1)  = gK*n10(t-1)^4;
        
        %% Equações Diferenciais para Voltagem
        %Compartimento 1
        g_Na = g_Na1(t-1);
        g_k  = g_k1(t-1);
        v = v1(t-1);
        vP = v2(t-1); % Posterior
        fc1_1 = -A*g_Na*(v-Ena)-A*g_k*(v-Ek)-A*g_l*(v-Evaz)-G*(v-vP)+I(t-1);
        fc1_2 = -A*g_Na*(v+(Dt/2)*fc1_1-Ena)-A*g_k*(v+(Dt/2)*fc1_1-Ek)-A*g_l*(v+(Dt/2)*fc1_1-Evaz) -G*(v+(Dt/2)*fc1_1-vP)+I(t-1);
        fc1_3 = -A*g_Na*(v+(Dt/2)*fc1_2-Ena)-A*g_k*(v+(Dt/2)*fc1_2-Ek)-A*g_l*(v+(Dt/2)*fc1_2-Evaz) -G*(v+(Dt/2)*fc1_2-vP)+I(t-1);
        fc1_4 = -A*g_Na*(v+Dt*fc1_3-Ena)-A*g_k*(v+Dt*fc1_3-Ek)-A*g_l*(v+Dt*fc1_3-Evaz) -G*(v+Dt*fc1_3-vP)+I(t-1);
        v1(t) = v1(t-1) + (Dt/(6*Cm))*(fc1_1+2*fc1_2+2*fc1_3+fc1_4);
        
        %Compartimento 2
        g_Na = g_Na2(t-1);
        g_k  = g_k2(t-1);
        v = v2(t-1);
        vA = v1(t-1); % Anterior
        vP = v3(t-1); % Posterior
        fc2_1 = -A*g_Na*(v-Ena)-A*g_k*(v-Ek)-A*g_l*(v-Evaz) +G*(vA-v)-G*(v-vP);
        fc2_2 = -A*g_Na*(v+(Dt/2)*fc2_1-Ena)-A*g_k*(v+(Dt/2)*fc2_1-Ek)-A*g_l*(v+(Dt/2)*fc2_1-Evaz) +G*(vA-(v+(Dt/2)*fc2_1))-G*(v+(Dt/2)*fc2_1-vP);
        fc2_3 = -A*g_Na*(v+(Dt/2)*fc2_2-Ena)-A*g_k*(v+(Dt/2)*fc2_2-Ek)-A*g_l*(v+(Dt/2)*fc2_2-Evaz) +G*(vA-(v+(Dt/2)*fc2_2))-G*(v+(Dt/2)*fc2_2-vP);
        fc2_4 = -A*g_Na*(v+Dt*fc2_3-Ena)-A*g_k*(v+Dt*fc2_3-Ek)-A*g_l*(v+Dt*fc2_3-Evaz) +G*(vA-(v+Dt*fc2_3))-G*(v+Dt*fc2_3-vP);
        v2(t) = v2(t-1) + (Dt/(6*Cm))*(fc2_1+2*fc2_2+2*fc2_3+fc2_4);
        
        %Compartimento 3
        g_Na = g_Na3(t-1);
        g_k  = g_k3(t-1);
        v = v3(t-1);
        vA = v2(t-1); % Anterior
        vP = v4(t-1); % Posterior
        fc3_1 = -A*g_Na*(v-Ena)-A*g_k*(v-Ek)-A*g_l*(v-Evaz) +G*(vA-v)-G*(v-vP);
        fc3_2 = -A*g_Na*(v+(Dt/2)*fc3_1-Ena)-A*g_k*(v+(Dt/2)*fc3_1-Ek)-A*g_l*(v+(Dt/2)*fc3_1-Evaz) +G*(vA-(v+(Dt/2)*fc3_1))-G*(v+(Dt/2)*fc3_1-vP);
        fc3_3 = -A*g_Na*(v+(Dt/2)*fc3_2-Ena)-A*g_k*(v+(Dt/2)*fc3_2-Ek)-A*g_l*(v+(Dt/2)*fc3_2-Evaz) +G*(vA-(v+(Dt/2)*fc3_2))-G*(v+(Dt/2)*fc3_2-vP);
        fc3_4 = -A*g_Na*(v+Dt*fc3_3-Ena)-A*g_k*(v+Dt*fc3_3-Ek)-A*g_l*(v+Dt*fc3_3-Evaz) +G*(vA-(v+Dt*fc3_3))-G*(v+Dt*fc3_3-vP);
        v3(t) = v3(t-1) + (Dt/(6*Cm))*(fc3_1+2*fc3_2+2*fc3_3+fc3_4);
        
        %Compartimento 4
        g_Na = g_Na4(t-1);
        g_k  = g_k4(t-1);
        v = v4(t-1);
        vA = v3(t-1); % Anterior
        vP = v5(t-1); % Posterior
        fc4_1 = -A*g_Na*(v-Ena)-A*g_k*(v-Ek)-A*g_l*(v-Evaz) +G*(vA-v)-G*(v-vP);
        fc4_2 = -A*g_Na*(v+(Dt/2)*fc4_1-Ena)-A*g_k*(v+(Dt/2)*fc4_1-Ek)-A*g_l*(v+(Dt/2)*fc4_1-Evaz) +G*(vA-(v+(Dt/2)*fc4_1))-G*(v+(Dt/2)*fc4_1-vP);
        fc4_3 = -A*g_Na*(v+(Dt/2)*fc4_2-Ena)-A*g_k*(v+(Dt/2)*fc4_2-Ek)-A*g_l*(v+(Dt/2)*fc4_2-Evaz) +G*(vA-(v+(Dt/2)*fc4_2))-G*(v+(Dt/2)*fc4_2-vP);
        fc4_4 = -A*g_Na*(v+Dt*fc4_3-Ena)-A*g_k*(v+Dt*fc4_3-Ek)-A*g_l*(v+Dt*fc4_3-Evaz) +G*(vA-(v+Dt*fc4_3))-G*(v+Dt*fc4_3-vP);
        v4(t) = v4(t-1) + (Dt/(6*Cm))*(fc4_1+2*fc4_2+2*fc4_3+fc4_4);
        
        %Compartimento 5
        g_Na = g_Na5(t-1);
        g_k  = g_k5(t-1);
        v = v5(t-1);
        vA = v4(t-1); % Anterior
        vP = v6(t-1); % Posterior
        fc5_1 = -A*g_Na*(v-Ena)-A*g_k*(v-Ek)-A*g_l*(v-Evaz) +G*(vA-v)-G*(v-vP);
        fc5_2 = -A*g_Na*(v+(Dt/2)*fc5_1-Ena)-A*g_k*(v+(Dt/2)*fc5_1-Ek)-A*g_l*(v+(Dt/2)*fc5_1-Evaz) +G*(vA-(v+(Dt/2)*fc5_1))-G*(v+(Dt/2)*fc5_1-vP);
        fc5_3 = -A*g_Na*(v+(Dt/2)*fc5_2-Ena)-A*g_k*(v+(Dt/2)*fc5_2-Ek)-A*g_l*(v+(Dt/2)*fc5_2-Evaz) +G*(vA-(v+(Dt/2)*fc5_2))-G*(v+(Dt/2)*fc5_2-vP);
        fc5_4 = -A*g_Na*(v+Dt*fc5_3-Ena)-A*g_k*(v+Dt*fc5_3-Ek)-A*g_l*(v+Dt*fc5_3-Evaz) +G*(vA-(v+Dt*fc5_3))-G*(v+Dt*fc5_3-vP);
        v5(t) = v5(t-1) + (Dt/(6*Cm))*(fc5_1+2*fc5_2+2*fc5_3+fc5_4);
        
        %Compartimento 6
        g_Na = g_Na6(t-1);
        g_k  = g_k6(t-1);
        v = v6(t-1);
        vA = v5(t-1); % Anterior
        vP = v7(t-1); % Posterior
        fc6_1 = -A*g_Na*(v-Ena)-A*g_k*(v-Ek)-A*g_l*(v-Evaz) +G*(vA-v)-G*(v-vP);
        fc6_2 = -A*g_Na*(v+(Dt/2)*fc6_1-Ena)-A*g_k*(v+(Dt/2)*fc6_1-Ek)-A*g_l*(v+(Dt/2)*fc6_1-Evaz) +G*(vA-(v+(Dt/2)*fc6_1))-G*(v+(Dt/2)*fc6_1-vP);
        fc6_3 = -A*g_Na*(v+(Dt/2)*fc6_2-Ena)-A*g_k*(v+(Dt/2)*fc6_2-Ek)-A*g_l*(v+(Dt/2)*fc6_2-Evaz) +G*(vA-(v+(Dt/2)*fc6_2))-G*(v+(Dt/2)*fc6_2-vP);
        fc6_4 = -A*g_Na*(v+Dt*fc6_3-Ena)-A*g_k*(v+Dt*fc6_3-Ek)-A*g_l*(v+Dt*fc6_3-Evaz) +G*(vA-(v+Dt*fc6_3))-G*(v+Dt*fc6_3-vP);
        v6(t) = v6(t-1) + (Dt/(6*Cm))*(fc6_1+2*fc6_2+2*fc6_3+fc6_4);
        
        %Compartimento 7
        g_Na = g_Na7(t-1);
        g_k  = g_k7(t-1);
        v = v7(t-1);
        vA = v6(t-1); % Anterior
        vP = v8(t-1); % Posterior
        fc7_1 = -A*g_Na*(v-Ena)-A*g_k*(v-Ek)-A*g_l*(v-Evaz) +G*(vA-v)-G*(v-vP);
        fc7_2 = -A*g_Na*(v+(Dt/2)*fc7_1-Ena)-A*g_k*(v+(Dt/2)*fc7_1-Ek)-A*g_l*(v+(Dt/2)*fc7_1-Evaz) +G*(vA-(v+(Dt/2)*fc7_1))-G*(v+(Dt/2)*fc7_1-vP);
        fc7_3 = -A*g_Na*(v+(Dt/2)*fc7_2-Ena)-A*g_k*(v+(Dt/2)*fc7_2-Ek)-A*g_l*(v+(Dt/2)*fc7_2-Evaz) +G*(vA-(v+(Dt/2)*fc7_2))-G*(v+(Dt/2)*fc7_2-vP);
        fc7_4 = -A*g_Na*(v+Dt*fc7_3-Ena)-A*g_k*(v+Dt*fc7_3-Ek)-A*g_l*(v+Dt*fc7_3-Evaz) +G*(vA-(v+Dt*fc7_3))-G*(v+Dt*fc7_3-vP);
        v7(t) = v7(t-1) + (Dt/(6*Cm))*(fc7_1+2*fc7_2+2*fc7_3+fc7_4);
        
        %Compartimento 8
        g_Na = g_Na8(t-1);
        g_k  = g_k8(t-1);
        v = v8(t-1);
        vA = v7(t-1); % Anterior
        vP = v9(t-1); % Posterior
        fc8_1 = -A*g_Na*(v-Ena)-A*g_k*(v-Ek)-A*g_l*(v-Evaz) +G*(vA-v)-G*(v-vP);
        fc8_2 = -A*g_Na*(v+(Dt/2)*fc8_1-Ena)-A*g_k*(v+(Dt/2)*fc8_1-Ek)-A*g_l*(v+(Dt/2)*fc8_1-Evaz) +G*(vA-(v+(Dt/2)*fc8_1))-G*(v+(Dt/2)*fc8_1-vP);
        fc8_3 = -A*g_Na*(v+(Dt/2)*fc8_2-Ena)-A*g_k*(v+(Dt/2)*fc8_2-Ek)-A*g_l*(v+(Dt/2)*fc8_2-Evaz) +G*(vA-(v+(Dt/2)*fc8_2))-G*(v+(Dt/2)*fc8_2-vP);
        fc8_4 = -A*g_Na*(v+Dt*fc8_3-Ena)-A*g_k*(v+Dt*fc8_3-Ek)-A*g_l*(v+Dt*fc8_3-Evaz) +G*(vA-(v+Dt*fc8_3))-G*(v+Dt*fc8_3-vP);
        v8(t) = v8(t-1) + (Dt/(6*Cm))*(fc8_1+2*fc8_2+2*fc8_3+fc8_4);
        
        %Compartimento 9
        g_Na = g_Na9(t-1);
        g_k  = g_k9(t-1);
        v = v9(t-1);
        vA = v8(t-1); % Anterior
        vP = v10(t-1); % Posterior
        fc9_1 = -A*g_Na*(v-Ena)-A*g_k*(v-Ek)-A*g_l*(v-Evaz) +G*(vA-v)-G*(v-vP);
        fc9_2 = -A*g_Na*(v+(Dt/2)*fc9_1-Ena)-A*g_k*(v+(Dt/2)*fc9_1-Ek)-A*g_l*(v+(Dt/2)*fc9_1-Evaz) +G*(vA-(v+(Dt/2)*fc9_1))-G*(v+(Dt/2)*fc9_1-vP);
        fc9_3 = -A*g_Na*(v+(Dt/2)*fc9_2-Ena)-A*g_k*(v+(Dt/2)*fc9_2-Ek)-A*g_l*(v+(Dt/2)*fc9_2-Evaz) +G*(vA-(v+(Dt/2)*fc9_2))-G*(v+(Dt/2)*fc9_2-vP);
        fc9_4 = -A*g_Na*(v+Dt*fc9_3-Ena)-A*g_k*(v+Dt*fc9_3-Ek)-A*g_l*(v+Dt*fc9_3-Evaz) +G*(vA-(v+Dt*fc9_3))-G*(v+Dt*fc9_3-vP);
        v9(t) = v9(t-1) + (Dt/(6*Cm))*(fc9_1+2*fc9_2+2*fc9_3+fc9_4);
        
        %Compartimento 10
        g_Na = g_Na10(t-1);
        g_k  = g_k10(t-1);
        v = v10(t-1);
        vA = v9(t-1); % Anterior
        fc10_1 = -A*g_Na*(v-Ena)-A*g_k*(v-Ek)-A*g_l*(v-Evaz)+G*(vA-v);
        fc10_2 = -A*g_Na*(v+(Dt/2)*fc10_1-Ena)-A*g_k*(v+(Dt/2)*fc10_1-Ek)-A*g_l*(v+(Dt/2)*fc10_1-Evaz) +G*(vA-(v+(Dt/2)*fc10_1));
        fc10_3 = -A*g_Na*(v+(Dt/2)*fc10_2-Ena)-A*g_k*(v+(Dt/2)*fc10_2-Ek)-A*g_l*(v+(Dt/2)*fc10_2-Evaz) +G*(vA-(v+(Dt/2)*fc10_2));
        fc10_4 = -A*g_Na*(v+Dt*fc10_3-Ena)-A*g_k*(v+Dt*fc10_3-Ek)-A*g_l*(v+Dt*fc10_3-Evaz) +G*(vA-(v+Dt*fc10_3));
        v10(t) = v10(t-1) + (Dt/(6*Cm))*(fc10_1+2*fc10_2+2*fc10_3+fc10_4);

        
        %% Cálculo: Probabilidade n por Runge-Kutta
        %Compartimento 1
        v = v1(t-1);
        n = n1(t-1);
        n1_1=      ((1.0-n)*alpha_n(v)-n*beta_n(v));
        n1_2=      ((1.0-n-(Dt/2)*n1_1)*alpha_n(v)-(n+(Dt/2)*n1_1)*beta_n(v));
        n1_3=      ((1.0-n-(Dt/2)*n1_2)*alpha_n(v)-(n+(Dt/2)*n1_2)*beta_n(v));
        n1_4=      ((1.0-n-(Dt)*n1_3)*alpha_n(v)-(n+(Dt)*n1_3)*beta_n(v));
        n1(t) = n1(t-1) + (Dt/6)*(n1_1+2*n1_2+2*n1_3+n1_4);
        %Compartimento 2
        v = v2(t-1);
        n = n2(t-1);
        n2_1=      ((1.0-n)*alpha_n(v)-n*beta_n(v));
        n2_2=      ((1.0-n-(Dt/2)*n2_1)*alpha_n(v)-(n+(Dt/2)*n2_1)*beta_n(v));
        n2_3=      ((1.0-n-(Dt/2)*n2_2)*alpha_n(v)-(n+(Dt/2)*n2_2)*beta_n(v));
        n2_4=      ((1.0-n-(Dt)*n2_3)*alpha_n(v)-(n+(Dt)*n2_3)*beta_n(v));
        n2(t) = n2(t-1) + (Dt/6)*(n2_1+2*n2_2+2*n2_3+n2_4);        
        %Compartimento 3
        v = v3(t-1);
        n = n3(t-1);
        n3_1=      ((1.0-n)*alpha_n(v)-n*beta_n(v));
        n3_2=      ((1.0-n-(Dt/2)*n3_1)*alpha_n(v)-(n+(Dt/2)*n3_1)*beta_n(v));
        n3_3=      ((1.0-n-(Dt/2)*n3_2)*alpha_n(v)-(n+(Dt/2)*n3_2)*beta_n(v));
        n3_4=      ((1.0-n-(Dt)*n3_3)*alpha_n(v)-(n+(Dt)*n3_3)*beta_n(v));
        n3(t) = n3(t-1) + (Dt/6)*(n3_1+2*n3_2+2*n3_3+n3_4);
        %Compartimento 4
        v = v4(t-1);
        n = n4(t-1);
        n4_1=      ((1.0-n)*alpha_n(v)-n*beta_n(v));
        n4_2=      ((1.0-n-(Dt/2)*n4_1)*alpha_n(v)-(n+(Dt/2)*n4_1)*beta_n(v));
        n4_3=      ((1.0-n-(Dt/2)*n4_2)*alpha_n(v)-(n+(Dt/2)*n4_2)*beta_n(v));
        n4_4=      ((1.0-n-(Dt)*n4_3)*alpha_n(v)-(n+(Dt)*n4_3)*beta_n(v));
        n4(t) = n4(t-1) + (Dt/6)*(n4_1+2*n4_2+2*n4_3+n4_4);
        %Compartimento 5
        v = v5(t-1);
        n = n5(t-1);
        n5_1=      ((1.0-n)*alpha_n(v)-n*beta_n(v));
        n5_2=      ((1.0-n-(Dt/2)*n5_1)*alpha_n(v)-(n+(Dt/2)*n5_1)*beta_n(v));
        n5_3=      ((1.0-n-(Dt/2)*n5_2)*alpha_n(v)-(n+(Dt/2)*n5_2)*beta_n(v));
        n5_4=      ((1.0-n-(Dt)*n5_3)*alpha_n(v)-(n+(Dt)*n5_3)*beta_n(v));
        n5(t) = n5(t-1) + (Dt/6)*(n5_1+2*n5_2+2*n5_3+n5_4);
        %Compartimento 6
        v = v6(t-1);
        n = n6(t-1);
        n6_1=      ((1.0-n)*alpha_n(v)-n*beta_n(v));
        n6_2=      ((1.0-n-(Dt/2)*n6_1)*alpha_n(v)-(n+(Dt/2)*n6_1)*beta_n(v));
        n6_3=      ((1.0-n-(Dt/2)*n6_2)*alpha_n(v)-(n+(Dt/2)*n6_2)*beta_n(v));
        n6_4=      ((1.0-n-(Dt)*n6_3)*alpha_n(v)-(n+(Dt)*n6_3)*beta_n(v));
        n6(t) = n6(t-1) + (Dt/6)*(n6_1+2*n6_2+2*n6_3+n6_4);       
        %Compartimento 7
        v = v7(t-1);
        n = n7(t-1);
        n7_1=      ((1.0-n)*alpha_n(v)-n*beta_n(v));
        n7_2=      ((1.0-n-(Dt/2)*n7_1)*alpha_n(v)-(n+(Dt/2)*n7_1)*beta_n(v));
        n7_3=      ((1.0-n-(Dt/2)*n7_2)*alpha_n(v)-(n+(Dt/2)*n7_2)*beta_n(v));
        n7_4=      ((1.0-n-(Dt)*n7_3)*alpha_n(v)-(n+(Dt)*n7_3)*beta_n(v));
        n7(t) = n7(t-1) + (Dt/6)*(n7_1+2*n7_2+2*n7_3+n7_4);       
        %Compartimento 8
        v = v8(t-1);
        n = n8(t-1);
        n8_1=      ((1.0-n)*alpha_n(v)-n*beta_n(v));
        n8_2=      ((1.0-n-(Dt/2)*n8_1)*alpha_n(v)-(n+(Dt/2)*n8_1)*beta_n(v));
        n8_3=      ((1.0-n-(Dt/2)*n8_2)*alpha_n(v)-(n+(Dt/2)*n8_2)*beta_n(v));
        n8_4=      ((1.0-n-(Dt)*n8_3)*alpha_n(v)-(n+(Dt)*n8_3)*beta_n(v));
        n8(t) = n8(t-1) + (Dt/6)*(n8_1+2*n8_2+2*n8_3+n8_4); 
        %Compartimento 9
        v = v9(t-1);
        n = n9(t-1);
        n9_1=      ((1.0-n)*alpha_n(v)-n*beta_n(v));
        n9_2=      ((1.0-n-(Dt/2)*n9_1)*alpha_n(v)-(n+(Dt/2)*n9_1)*beta_n(v));
        n9_3=      ((1.0-n-(Dt/2)*n9_2)*alpha_n(v)-(n+(Dt/2)*n9_2)*beta_n(v));
        n9_4=      ((1.0-n-(Dt)*n9_3)*alpha_n(v)-(n+(Dt)*n9_3)*beta_n(v));
        n9(t) = n9(t-1) + (Dt/6)*(n9_1+2*n9_2+2*n9_3+n9_4); 
        %Compartimento 10
        v = v10(t-1);
        n = n10(t-1);
        n10_1=      ((1.0-n)*alpha_n(v)-n*beta_n(v));
        n10_2=      ((1.0-n-(Dt/2)*n10_1)*alpha_n(v)-(n+(Dt/2)*n10_1)*beta_n(v));
        n10_3=      ((1.0-n-(Dt/2)*n10_2)*alpha_n(v)-(n+(Dt/2)*n10_2)*beta_n(v));
        n10_4=      ((1.0-n-(Dt)*n10_3)*alpha_n(v)-(n+(Dt)*n10_3)*beta_n(v));
        n10(t) = n10(t-1) + (Dt/6)*(n10_1+2*n10_2+2*n10_3+n10_4);        
        
        
        %% Cálculo: Probabilidade m por Runge-Kutta
        %Compartimento 1
        v = v1(t-1);
        m = m1(t-1);
        m1_1=       ((1.0-m)*alpha_m(v)-m*beta_m(v));
        m1_2=       ((1.0-m-(Dt/2)*m1_1)*alpha_m(v)-(m+(Dt/2)*m1_1)*beta_m(v));
        m1_3=       ((1.0-m-(Dt/2)*m1_2)*alpha_m(v)-(m+(Dt/2)*m1_2)*beta_m(v));
        m1_4=       ((1.0-m-Dt*m1_3)*alpha_m(v)-(m+Dt*m1_3)*beta_m(v));
        m1(t) = m1(t-1) + (Dt/6)*(m1_1+2*m1_2+2*m1_3+m1_4); 
        %Compartimento 2
        v = v2(t-1);
        m = m2(t-1);
        m2_1=       ((1.0-m)*alpha_m(v)-m*beta_m(v));
        m2_2=       ((1.0-m-(Dt/2)*m2_1)*alpha_m(v)-(m+(Dt/2)*m2_1)*beta_m(v));
        m2_3=       ((1.0-m-(Dt/2)*m2_2)*alpha_m(v)-(m+(Dt/2)*m2_2)*beta_m(v));
        m2_4=       ((1.0-m-Dt*m2_3)*alpha_m(v)-(m+Dt*m2_3)*beta_m(v));
        m2(t) = m2(t-1) + (Dt/6)*(m2_1+2*m2_2+2*m2_3+m2_4); 
        %Compartimento 3
        v = v3(t-1);
        m = m3(t-1);
        m3_1=       ((1.0-m)*alpha_m(v)-m*beta_m(v));
        m3_2=       ((1.0-m-(Dt/2)*m3_1)*alpha_m(v)-(m+(Dt/2)*m3_1)*beta_m(v));
        m3_3=       ((1.0-m-(Dt/2)*m3_2)*alpha_m(v)-(m+(Dt/2)*m3_2)*beta_m(v));
        m3_4=       ((1.0-m-Dt*m3_3)*alpha_m(v)-(m+Dt*m3_3)*beta_m(v));
        m3(t) = m3(t-1) + (Dt/6)*(m3_1+2*m3_2+2*m3_3+m3_4); 
        %Compartimento 4
        v = v4(t-1);
        m = m4(t-1);
        m4_1=       ((1.0-m)*alpha_m(v)-m*beta_m(v));
        m4_2=       ((1.0-m-(Dt/2)*m4_1)*alpha_m(v)-(m+(Dt/2)*m4_1)*beta_m(v));
        m4_3=       ((1.0-m-(Dt/2)*m4_2)*alpha_m(v)-(m+(Dt/2)*m4_2)*beta_m(v));
        m4_4=       ((1.0-m-Dt*m4_3)*alpha_m(v)-(m+Dt*m4_3)*beta_m(v));
        m4(t) = m4(t-1) + (Dt/6)*(m4_1+2*m4_2+2*m4_3+m4_4); 
        %Compartimento 5
        v = v5(t-1);
        m = m5(t-1);
        m5_1=       ((1.0-m)*alpha_m(v)-m*beta_m(v));
        m5_2=       ((1.0-m-(Dt/2)*m5_1)*alpha_m(v)-(m+(Dt/2)*m5_1)*beta_m(v));
        m5_3=       ((1.0-m-(Dt/2)*m5_2)*alpha_m(v)-(m+(Dt/2)*m5_2)*beta_m(v));
        m5_4=       ((1.0-m-Dt*m5_3)*alpha_m(v)-(m+Dt*m5_3)*beta_m(v));
        m5(t) = m5(t-1) + (Dt/6)*(m5_1+2*m5_2+2*m5_3+m5_4); 
        %Compartimento 6
        v = v6(t-1);
        m = m6(t-1);
        m6_1=       ((1.0-m)*alpha_m(v)-m*beta_m(v));
        m6_2=       ((1.0-m-(Dt/2)*m6_1)*alpha_m(v)-(m+(Dt/2)*m6_1)*beta_m(v));
        m6_3=       ((1.0-m-(Dt/2)*m6_2)*alpha_m(v)-(m+(Dt/2)*m6_2)*beta_m(v));
        m6_4=       ((1.0-m-Dt*m6_3)*alpha_m(v)-(m+Dt*m6_3)*beta_m(v));
        m6(t) = m6(t-1) + (Dt/6)*(m6_1+2*m6_2+2*m6_3+m6_4);
        %Compartimento 7
        v = v7(t-1);
        m = m7(t-1);
        m7_1=       ((1.0-m)*alpha_m(v)-m*beta_m(v));
        m7_2=       ((1.0-m-(Dt/2)*m7_1)*alpha_m(v)-(m+(Dt/2)*m7_1)*beta_m(v));
        m7_3=       ((1.0-m-(Dt/2)*m7_2)*alpha_m(v)-(m+(Dt/2)*m7_2)*beta_m(v));
        m7_4=       ((1.0-m-Dt*m7_3)*alpha_m(v)-(m+Dt*m7_3)*beta_m(v));
        m7(t) = m7(t-1) + (Dt/6)*(m7_1+2*m7_2+2*m7_3+m7_4);
        %Compartimento 8
        v = v8(t-1);
        m = m8(t-1);
        m8_1=       ((1.0-m)*alpha_m(v)-m*beta_m(v));
        m8_2=       ((1.0-m-(Dt/2)*m8_1)*alpha_m(v)-(m+(Dt/2)*m8_1)*beta_m(v));
        m8_3=       ((1.0-m-(Dt/2)*m8_2)*alpha_m(v)-(m+(Dt/2)*m8_2)*beta_m(v));
        m8_4=       ((1.0-m-Dt*m8_3)*alpha_m(v)-(m+Dt*m8_3)*beta_m(v));
        m8(t) = m8(t-1) + (Dt/6)*(m8_1+2*m8_2+2*m8_3+m8_4);     
        %Compartimento 9
        v = v9(t-1);
        m = m9(t-1);
        m9_1=       ((1.0-m)*alpha_m(v)-m*beta_m(v));
        m9_2=       ((1.0-m-(Dt/2)*m9_1)*alpha_m(v)-(m+(Dt/2)*m9_1)*beta_m(v));
        m9_3=       ((1.0-m-(Dt/2)*m9_2)*alpha_m(v)-(m+(Dt/2)*m9_2)*beta_m(v));
        m9_4=       ((1.0-m-Dt*m9_3)*alpha_m(v)-(m+Dt*m9_3)*beta_m(v));
        m9(t) = m9(t-1) + (Dt/6)*(m9_1+2*m9_2+2*m9_3+m9_4);  
        %Compartimento 10
        v = v10(t-1);
        m = m10(t-1);
        m10_1=       ((1.0-m)*alpha_m(v)-m*beta_m(v));
        m10_2=       ((1.0-m-(Dt/2)*m10_1)*alpha_m(v)-(m+(Dt/2)*m10_1)*beta_m(v));
        m10_3=       ((1.0-m-(Dt/2)*m10_2)*alpha_m(v)-(m+(Dt/2)*m10_2)*beta_m(v));
        m10_4=       ((1.0-m-Dt*m10_3)*alpha_m(v)-(m+Dt*m10_3)*beta_m(v));
        m10(t) = m10(t-1) + (Dt/6)*(m10_1+2*m10_2+2*m10_3+m10_4);  
                
        %% Cálculo: Probabilidade n por Runge-Kutta
        %Compartimento 1
        v = v1(t-1);
        h = h1(t-1);
        h1_1=       ((1.0-h)*alpha_h(v)-(h)*beta_h(v));
        h1_2=       ((1.0-h-(Dt/2)*h1_1)*alpha_h(v)-(h+(Dt/2)*h1_1)*beta_h(v));
        h1_3=       ((1.0-h-(Dt/2)*h1_2)*alpha_h(v)-(h+(Dt/2)*h1_2)*beta_h(v));
        h1_4=       ((1.0-h-Dt*h1_3)*alpha_h(v)-(h+Dt*h1_3)*beta_h(v));
        h1(t) = h1(t-1) + (Dt/6)*(h1_1+2*h1_2+2*h1_3+h1_4);
        %Compartimento 2
        v = v2(t-1);
        h = h2(t-1);
        h2_1=       ((1.0-h)*alpha_h(v)-(h)*beta_h(v));
        h2_2=       ((1.0-h-(Dt/2)*h2_1)*alpha_h(v)-(h+(Dt/2)*h2_1)*beta_h(v));
        h2_3=       ((1.0-h-(Dt/2)*h2_2)*alpha_h(v)-(h+(Dt/2)*h2_2)*beta_h(v));
        h2_4=       ((1.0-h-Dt*h2_3)*alpha_h(v)-(h+Dt*h2_3)*beta_h(v));
        h2(t) = h2(t-1) + (Dt/6)*(h2_1+2*h2_2+2*h2_3+h2_4);
        %Compartimento 3
        v = v3(t-1);
        h = h3(t-1);
        h3_1=       ((1.0-h)*alpha_h(v)-(h)*beta_h(v));
        h3_2=       ((1.0-h-(Dt/2)*h3_1)*alpha_h(v)-(h+(Dt/2)*h3_1)*beta_h(v));
        h3_3=       ((1.0-h-(Dt/2)*h3_2)*alpha_h(v)-(h+(Dt/2)*h3_2)*beta_h(v));
        h3_4=       ((1.0-h-Dt*h3_3)*alpha_h(v)-(h+Dt*h3_3)*beta_h(v));
        h3(t) = h3(t-1) + (Dt/6)*(h3_1+2*h3_2+2*h3_3+h3_4);
        %Compartimento 4
        v = v4(t-1);
        h = h4(t-1);
        h4_1=       ((1.0-h)*alpha_h(v)-(h)*beta_h(v));
        h4_2=       ((1.0-h-(Dt/2)*h4_1)*alpha_h(v)-(h+(Dt/2)*h4_1)*beta_h(v));
        h4_3=       ((1.0-h-(Dt/2)*h4_2)*alpha_h(v)-(h+(Dt/2)*h4_2)*beta_h(v));
        h4_4=       ((1.0-h-Dt*h4_3)*alpha_h(v)-(h+Dt*h4_3)*beta_h(v));
        h4(t) = h4(t-1) + (Dt/6)*(h4_1+2*h4_2+2*h4_3+h4_4);
        %Compartimento 5
        v = v5(t-1);
        h = h5(t-1);
        h5_1=       ((1.0-h)*alpha_h(v)-(h)*beta_h(v));
        h5_2=       ((1.0-h-(Dt/2)*h5_1)*alpha_h(v)-(h+(Dt/2)*h5_1)*beta_h(v));
        h5_3=       ((1.0-h-(Dt/2)*h5_2)*alpha_h(v)-(h+(Dt/2)*h5_2)*beta_h(v));
        h5_4=       ((1.0-h-Dt*h5_3)*alpha_h(v)-(h+Dt*h5_3)*beta_h(v));
        h5(t) = h5(t-1) + (Dt/6)*(h5_1+2*h5_2+2*h5_3+h5_4);
        %Compartimento 6
        v = v6(t-1);
        h = h6(t-1);
        h6_1=       ((1.0-h)*alpha_h(v)-(h)*beta_h(v));
        h6_2=       ((1.0-h-(Dt/2)*h6_1)*alpha_h(v)-(h+(Dt/2)*h6_1)*beta_h(v));
        h6_3=       ((1.0-h-(Dt/2)*h6_2)*alpha_h(v)-(h+(Dt/2)*h6_2)*beta_h(v));
        h6_4=       ((1.0-h-Dt*h6_3)*alpha_h(v)-(h+Dt*h6_3)*beta_h(v));
        h6(t) = h6(t-1) + (Dt/6)*(h6_1+2*h6_2+2*h6_3+h6_4);
        %Compartimento 7
        v = v7(t-1);
        h = h7(t-1);
        h7_1=       ((1.0-h)*alpha_h(v)-(h)*beta_h(v));
        h7_2=       ((1.0-h-(Dt/2)*h7_1)*alpha_h(v)-(h+(Dt/2)*h7_1)*beta_h(v));
        h7_3=       ((1.0-h-(Dt/2)*h7_2)*alpha_h(v)-(h+(Dt/2)*h7_2)*beta_h(v));
        h7_4=       ((1.0-h-Dt*h7_3)*alpha_h(v)-(h+Dt*h7_3)*beta_h(v));
        h7(t) = h7(t-1) + (Dt/6)*(h7_1+2*h7_2+2*h7_3+h7_4);
        %Compartimento 8
        v = v8(t-1);
        h = h8(t-1);
        h8_1=       ((1.0-h)*alpha_h(v)-(h)*beta_h(v));
        h8_2=       ((1.0-h-(Dt/2)*h8_1)*alpha_h(v)-(h+(Dt/2)*h8_1)*beta_h(v));
        h8_3=       ((1.0-h-(Dt/2)*h8_2)*alpha_h(v)-(h+(Dt/2)*h8_2)*beta_h(v));
        h8_4=       ((1.0-h-Dt*h8_3)*alpha_h(v)-(h+Dt*h8_3)*beta_h(v));
        h8(t) = h8(t-1) + (Dt/6)*(h8_1+2*h8_2+2*h8_3+h8_4);
        %Compartimento 9
        v = v9(t-1);
        h = h9(t-1);
        h9_1=       ((1.0-h)*alpha_h(v)-(h)*beta_h(v));
        h9_2=       ((1.0-h-(Dt/2)*h9_1)*alpha_h(v)-(h+(Dt/2)*h9_1)*beta_h(v));
        h9_3=       ((1.0-h-(Dt/2)*h9_2)*alpha_h(v)-(h+(Dt/2)*h9_2)*beta_h(v));
        h9_4=       ((1.0-h-Dt*h9_3)*alpha_h(v)-(h+Dt*h9_3)*beta_h(v));
        h9(t) = h9(t-1) + (Dt/6)*(h9_1+2*h9_2+2*h9_3+h9_4);
        %Compartimento 10
        v = v10(t-1);
        h = h10(t-1);
        h10_1=       ((1.0-h)*alpha_h(v)-(h)*beta_h(v));
        h10_2=       ((1.0-h-(Dt/2)*h10_1)*alpha_h(v)-(h+(Dt/2)*h10_1)*beta_h(v));
        h10_3=       ((1.0-h-(Dt/2)*h10_2)*alpha_h(v)-(h+(Dt/2)*h10_2)*beta_h(v));
        h10_4=       ((1.0-h-Dt*h10_3)*alpha_h(v)-(h+Dt*h10_3)*beta_h(v));
        h10(t) = h10(t-1) + (Dt/6)*(h10_1+2*h10_2+2*h10_3+h10_4);      

        
        tv(t)=t*Dt;
    end
end


% Exercício 5: Modelo de 10 Compartimentos Ativos             [FUNCTION]
% 5.c) Análise: Velocidade de Propagação do Potencial de Ação [FUNCTION]

% Função praticamente igual à COMP10ATIV_RKUTA4; foram alterados apenas os parâmetros do modelo (A,G,Cm)
function [v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,tv,I] = COMP10ATIVv2_RKUTA4(v0,ti,tf,Dt,ti_I,tf_I,Iin,a_alt)

	%parametros do modelo
    %a_alt = 1; %Proporção de alteração no parâmetro a.
    A  = a_alt*1.25e-5;
    G  = a_alt*1.25e-4;
    Cm = a_alt*1.25e-5;   
    
    %parametros do modelo ativo
    Ena  = 115;
    Ek   = -12;
    Evaz = 10.163; 
    gNa  = 120;
    gK   = 36;
    gV   = 0.3;
    
    %Dimensão do vetores (conforme granularidade temporal)
    len = ceil((tf-ti)/Dt);
    
    %Vetores das Condutâncias 
    g_Na1 = zeros(len,1);
    g_k1  = zeros(len,1);
    g_Na2 = zeros(len,1);
    g_k2  = zeros(len,1);
    g_Na3 = zeros(len,1);
    g_k3  = zeros(len,1);
    g_Na4 = zeros(len,1);
    g_k4  = zeros(len,1);
    g_Na5 = zeros(len,1);
    g_k5  = zeros(len,1);
    g_Na6 = zeros(len,1);
    g_k6  = zeros(len,1);
    g_Na7 = zeros(len,1);
    g_k7  = zeros(len,1);
    g_Na8 = zeros(len,1);
    g_k8  = zeros(len,1);
    g_Na9 = zeros(len,1);
    g_k9  = zeros(len,1);
    g_Na10= zeros(len,1);
    g_k10 = zeros(len,1);
    g_l  = gV;
    
    %Vetor de Correte Injetada
    I=zeros(len,1);
    for i=1:len
        % Decisao quanto a entrada de corrente
       if ((i*Dt>ti_I) && (i*Dt<tf_I))
            I(i)=Iin;
       else
            I(i)=0;
       end
    end
   
    %Vetores de Probabilidades (Gates)
    n1  = zeros(len,1);
    m1  = zeros(len,1);
    h1  = zeros(len,1);
    n2  = zeros(len,1);
    m2  = zeros(len,1);
    h2  = zeros(len,1);
    n3  = zeros(len,1);
    m3  = zeros(len,1);
    h3  = zeros(len,1);
    n4  = zeros(len,1);
    m4  = zeros(len,1);
    h4  = zeros(len,1);
    n5  = zeros(len,1);
    m5  = zeros(len,1);
    h5  = zeros(len,1);
    n6  = zeros(len,1);
    m6  = zeros(len,1);
    h6  = zeros(len,1);
    n7  = zeros(len,1);
    m7  = zeros(len,1);
    h7  = zeros(len,1);
    n8  = zeros(len,1);
    m8  = zeros(len,1);
    h8  = zeros(len,1);
    n9  = zeros(len,1);
    m9  = zeros(len,1);
    h9  = zeros(len,1);
    n10 = zeros(len,1);
    m10 = zeros(len,1);
    h10 = zeros(len,1);
    
    %Vetor Voltagens e Tempo
    v1   = zeros(len,1);
    v2   = zeros(len,1);
    v3   = zeros(len,1);
    v4   = zeros(len,1);
    v5   = zeros(len,1);
    v6   = zeros(len,1);
    v7   = zeros(len,1);
    v8   = zeros(len,1);
    v9   = zeros(len,1);
    v10  = zeros(len,1);
    v1(1) = v0;
    v2(1) = 0;
    v3(1) = 0;
    v4(1) = 0;
    v5(1) = 0;
    v6(1) = 0;
    v7(1) = 0;
    v8(1) = 0;
    v9(1) = 0;
    v10(1)= 0;
    
    tv   = zeros(len,1);   
    tv(1)= Dt;
    
    %Probabilidade de Abertura dos Canais (para célula em repouso)
    n1(1) = alpha_n(v1(1))/(alpha_n(v1(1))+beta_n(v1(1)));
    m1(1) = alpha_m(v1(1))/(alpha_m(v1(1))+beta_m(v1(1)));
    h1(1) = alpha_h(v1(1))/(alpha_h(v1(1))+beta_h(v1(1)));
    n2(1) = alpha_n(v2(1))/(alpha_n(v2(1))+beta_n(v2(1)));
    m2(1) = alpha_m(v2(1))/(alpha_m(v2(1))+beta_m(v2(1)));
    h2(1) = alpha_h(v2(1))/(alpha_h(v2(1))+beta_h(v2(1)));
    n3(1) = alpha_n(v3(1))/(alpha_n(v3(1))+beta_n(v3(1)));
    m3(1) = alpha_m(v3(1))/(alpha_m(v3(1))+beta_m(v3(1)));
    h3(1) = alpha_h(v3(1))/(alpha_h(v3(1))+beta_h(v3(1)));
    n4(1) = alpha_n(v4(1))/(alpha_n(v4(1))+beta_n(v4(1)));
    m4(1) = alpha_m(v4(1))/(alpha_m(v4(1))+beta_m(v4(1)));
    h4(1) = alpha_h(v4(1))/(alpha_h(v4(1))+beta_h(v4(1)));
    n5(1) = alpha_n(v5(1))/(alpha_n(v5(1))+beta_n(v5(1)));
    m5(1) = alpha_m(v5(1))/(alpha_m(v5(1))+beta_m(v5(1)));
    h5(1) = alpha_h(v5(1))/(alpha_h(v5(1))+beta_h(v5(1)));
    n6(1) = alpha_n(v6(1))/(alpha_n(v6(1))+beta_n(v6(1)));
    m6(1) = alpha_m(v6(1))/(alpha_m(v6(1))+beta_m(v6(1)));
    h6(1) = alpha_h(v6(1))/(alpha_h(v6(1))+beta_h(v6(1)));
    n7(1) = alpha_n(v7(1))/(alpha_n(v7(1))+beta_n(v7(1)));
    m7(1) = alpha_m(v7(1))/(alpha_m(v7(1))+beta_m(v7(1)));
    h7(1) = alpha_h(v7(1))/(alpha_h(v7(1))+beta_h(v7(1)));
    n8(1) = alpha_n(v8(1))/(alpha_n(v8(1))+beta_n(v8(1)));
    m8(1) = alpha_m(v8(1))/(alpha_m(v8(1))+beta_m(v8(1)));
    h8(1) = alpha_h(v8(1))/(alpha_h(v8(1))+beta_h(v8(1)));
    n9(1) = alpha_n(v9(1))/(alpha_n(v9(1))+beta_n(v9(1)));
    m9(1) = alpha_m(v9(1))/(alpha_m(v9(1))+beta_m(v9(1)));
    h9(1) = alpha_h(v9(1))/(alpha_h(v9(1))+beta_h(v9(1)));
    n10(1) = alpha_n(v10(1))/(alpha_n(v10(1))+beta_n(v10(1)));
    m10(1) = alpha_m(v10(1))/(alpha_m(v10(1))+beta_m(v10(1)));
    h10(1) = alpha_h(v10(1))/(alpha_h(v10(1))+beta_h(v10(1)));
    
    for t=2:len
        % Cálculo: Condutâncias
        g_Na1(t-1) = gNa*(m1(t-1)^3)*h1(t-1);
        g_k1(t-1)  = gK*n1(t-1)^4;
        g_Na2(t-1) = gNa*(m2(t-1)^3)*h2(t-1);
        g_k2(t-1)  = gK*n2(t-1)^4;
        g_Na3(t-1) = gNa*(m3(t-1)^3)*h3(t-1);
        g_k3(t-1)  = gK*n3(t-1)^4;
        g_Na4(t-1) = gNa*(m4(t-1)^3)*h4(t-1);
        g_k4(t-1)  = gK*n4(t-1)^4;
        g_Na5(t-1) = gNa*(m5(t-1)^3)*h5(t-1);
        g_k5(t-1)  = gK*n5(t-1)^4;
        g_Na6(t-1) = gNa*(m6(t-1)^3)*h6(t-1);
        g_k6(t-1)  = gK*n6(t-1)^4;
        g_Na7(t-1) = gNa*(m7(t-1)^3)*h7(t-1);
        g_k7(t-1)  = gK*n7(t-1)^4;
        g_Na8(t-1) = gNa*(m8(t-1)^3)*h8(t-1);
        g_k8(t-1)  = gK*n8(t-1)^4;
        g_Na9(t-1) = gNa*(m9(t-1)^3)*h9(t-1);
        g_k9(t-1)  = gK*n9(t-1)^4;
        g_Na10(t-1) = gNa*(m10(t-1)^3)*h10(t-1);
        g_k10(t-1)  = gK*n10(t-1)^4;
        
        %% Equações Diferenciais para Voltagem
        %Compartimento 1
        g_Na = g_Na1(t-1);
        g_k  = g_k1(t-1);
        v = v1(t-1);
        vP = v2(t-1); % Posterior
        fc1_1 = -A*g_Na*(v-Ena)-A*g_k*(v-Ek)-A*g_l*(v-Evaz)-G*(v-vP)+I(t-1);
        fc1_2 = -A*g_Na*(v+(Dt/2)*fc1_1-Ena)-A*g_k*(v+(Dt/2)*fc1_1-Ek)-A*g_l*(v+(Dt/2)*fc1_1-Evaz) -G*(v+(Dt/2)*fc1_1-vP)+I(t-1);
        fc1_3 = -A*g_Na*(v+(Dt/2)*fc1_2-Ena)-A*g_k*(v+(Dt/2)*fc1_2-Ek)-A*g_l*(v+(Dt/2)*fc1_2-Evaz) -G*(v+(Dt/2)*fc1_2-vP)+I(t-1);
        fc1_4 = -A*g_Na*(v+Dt*fc1_3-Ena)-A*g_k*(v+Dt*fc1_3-Ek)-A*g_l*(v+Dt*fc1_3-Evaz) -G*(v+Dt*fc1_3-vP)+I(t-1);
        v1(t) = v1(t-1) + (Dt/(6*Cm))*(fc1_1+2*fc1_2+2*fc1_3+fc1_4);
        
        %Compartimento 2
        g_Na = g_Na2(t-1);
        g_k  = g_k2(t-1);
        v = v2(t-1);
        vA = v1(t-1); % Anterior
        vP = v3(t-1); % Posterior
        fc2_1 = -A*g_Na*(v-Ena)-A*g_k*(v-Ek)-A*g_l*(v-Evaz) +G*(vA-v)-G*(v-vP);
        fc2_2 = -A*g_Na*(v+(Dt/2)*fc2_1-Ena)-A*g_k*(v+(Dt/2)*fc2_1-Ek)-A*g_l*(v+(Dt/2)*fc2_1-Evaz) +G*(vA-(v+(Dt/2)*fc2_1))-G*(v+(Dt/2)*fc2_1-vP);
        fc2_3 = -A*g_Na*(v+(Dt/2)*fc2_2-Ena)-A*g_k*(v+(Dt/2)*fc2_2-Ek)-A*g_l*(v+(Dt/2)*fc2_2-Evaz) +G*(vA-(v+(Dt/2)*fc2_2))-G*(v+(Dt/2)*fc2_2-vP);
        fc2_4 = -A*g_Na*(v+Dt*fc2_3-Ena)-A*g_k*(v+Dt*fc2_3-Ek)-A*g_l*(v+Dt*fc2_3-Evaz) +G*(vA-(v+Dt*fc2_3))-G*(v+Dt*fc2_3-vP);
        v2(t) = v2(t-1) + (Dt/(6*Cm))*(fc2_1+2*fc2_2+2*fc2_3+fc2_4);
        
        %Compartimento 3
        g_Na = g_Na3(t-1);
        g_k  = g_k3(t-1);
        v = v3(t-1);
        vA = v2(t-1); % Anterior
        vP = v4(t-1); % Posterior
        fc3_1 = -A*g_Na*(v-Ena)-A*g_k*(v-Ek)-A*g_l*(v-Evaz) +G*(vA-v)-G*(v-vP);
        fc3_2 = -A*g_Na*(v+(Dt/2)*fc3_1-Ena)-A*g_k*(v+(Dt/2)*fc3_1-Ek)-A*g_l*(v+(Dt/2)*fc3_1-Evaz) +G*(vA-(v+(Dt/2)*fc3_1))-G*(v+(Dt/2)*fc3_1-vP);
        fc3_3 = -A*g_Na*(v+(Dt/2)*fc3_2-Ena)-A*g_k*(v+(Dt/2)*fc3_2-Ek)-A*g_l*(v+(Dt/2)*fc3_2-Evaz) +G*(vA-(v+(Dt/2)*fc3_2))-G*(v+(Dt/2)*fc3_2-vP);
        fc3_4 = -A*g_Na*(v+Dt*fc3_3-Ena)-A*g_k*(v+Dt*fc3_3-Ek)-A*g_l*(v+Dt*fc3_3-Evaz) +G*(vA-(v+Dt*fc3_3))-G*(v+Dt*fc3_3-vP);
        v3(t) = v3(t-1) + (Dt/(6*Cm))*(fc3_1+2*fc3_2+2*fc3_3+fc3_4);
        
        %Compartimento 4
        g_Na = g_Na4(t-1);
        g_k  = g_k4(t-1);
        v = v4(t-1);
        vA = v3(t-1); % Anterior
        vP = v5(t-1); % Posterior
        fc4_1 = -A*g_Na*(v-Ena)-A*g_k*(v-Ek)-A*g_l*(v-Evaz) +G*(vA-v)-G*(v-vP);
        fc4_2 = -A*g_Na*(v+(Dt/2)*fc4_1-Ena)-A*g_k*(v+(Dt/2)*fc4_1-Ek)-A*g_l*(v+(Dt/2)*fc4_1-Evaz) +G*(vA-(v+(Dt/2)*fc4_1))-G*(v+(Dt/2)*fc4_1-vP);
        fc4_3 = -A*g_Na*(v+(Dt/2)*fc4_2-Ena)-A*g_k*(v+(Dt/2)*fc4_2-Ek)-A*g_l*(v+(Dt/2)*fc4_2-Evaz) +G*(vA-(v+(Dt/2)*fc4_2))-G*(v+(Dt/2)*fc4_2-vP);
        fc4_4 = -A*g_Na*(v+Dt*fc4_3-Ena)-A*g_k*(v+Dt*fc4_3-Ek)-A*g_l*(v+Dt*fc4_3-Evaz) +G*(vA-(v+Dt*fc4_3))-G*(v+Dt*fc4_3-vP);
        v4(t) = v4(t-1) + (Dt/(6*Cm))*(fc4_1+2*fc4_2+2*fc4_3+fc4_4);
        
        %Compartimento 5
        g_Na = g_Na5(t-1);
        g_k  = g_k5(t-1);
        v = v5(t-1);
        vA = v4(t-1); % Anterior
        vP = v6(t-1); % Posterior
        fc5_1 = -A*g_Na*(v-Ena)-A*g_k*(v-Ek)-A*g_l*(v-Evaz) +G*(vA-v)-G*(v-vP);
        fc5_2 = -A*g_Na*(v+(Dt/2)*fc5_1-Ena)-A*g_k*(v+(Dt/2)*fc5_1-Ek)-A*g_l*(v+(Dt/2)*fc5_1-Evaz) +G*(vA-(v+(Dt/2)*fc5_1))-G*(v+(Dt/2)*fc5_1-vP);
        fc5_3 = -A*g_Na*(v+(Dt/2)*fc5_2-Ena)-A*g_k*(v+(Dt/2)*fc5_2-Ek)-A*g_l*(v+(Dt/2)*fc5_2-Evaz) +G*(vA-(v+(Dt/2)*fc5_2))-G*(v+(Dt/2)*fc5_2-vP);
        fc5_4 = -A*g_Na*(v+Dt*fc5_3-Ena)-A*g_k*(v+Dt*fc5_3-Ek)-A*g_l*(v+Dt*fc5_3-Evaz) +G*(vA-(v+Dt*fc5_3))-G*(v+Dt*fc5_3-vP);
        v5(t) = v5(t-1) + (Dt/(6*Cm))*(fc5_1+2*fc5_2+2*fc5_3+fc5_4);
        
        %Compartimento 6
        g_Na = g_Na6(t-1);
        g_k  = g_k6(t-1);
        v = v6(t-1);
        vA = v5(t-1); % Anterior
        vP = v7(t-1); % Posterior
        fc6_1 = -A*g_Na*(v-Ena)-A*g_k*(v-Ek)-A*g_l*(v-Evaz) +G*(vA-v)-G*(v-vP);
        fc6_2 = -A*g_Na*(v+(Dt/2)*fc6_1-Ena)-A*g_k*(v+(Dt/2)*fc6_1-Ek)-A*g_l*(v+(Dt/2)*fc6_1-Evaz) +G*(vA-(v+(Dt/2)*fc6_1))-G*(v+(Dt/2)*fc6_1-vP);
        fc6_3 = -A*g_Na*(v+(Dt/2)*fc6_2-Ena)-A*g_k*(v+(Dt/2)*fc6_2-Ek)-A*g_l*(v+(Dt/2)*fc6_2-Evaz) +G*(vA-(v+(Dt/2)*fc6_2))-G*(v+(Dt/2)*fc6_2-vP);
        fc6_4 = -A*g_Na*(v+Dt*fc6_3-Ena)-A*g_k*(v+Dt*fc6_3-Ek)-A*g_l*(v+Dt*fc6_3-Evaz) +G*(vA-(v+Dt*fc6_3))-G*(v+Dt*fc6_3-vP);
        v6(t) = v6(t-1) + (Dt/(6*Cm))*(fc6_1+2*fc6_2+2*fc6_3+fc6_4);
        
        %Compartimento 7
        g_Na = g_Na7(t-1);
        g_k  = g_k7(t-1);
        v = v7(t-1);
        vA = v6(t-1); % Anterior
        vP = v8(t-1); % Posterior
        fc7_1 = -A*g_Na*(v-Ena)-A*g_k*(v-Ek)-A*g_l*(v-Evaz) +G*(vA-v)-G*(v-vP);
        fc7_2 = -A*g_Na*(v+(Dt/2)*fc7_1-Ena)-A*g_k*(v+(Dt/2)*fc7_1-Ek)-A*g_l*(v+(Dt/2)*fc7_1-Evaz) +G*(vA-(v+(Dt/2)*fc7_1))-G*(v+(Dt/2)*fc7_1-vP);
        fc7_3 = -A*g_Na*(v+(Dt/2)*fc7_2-Ena)-A*g_k*(v+(Dt/2)*fc7_2-Ek)-A*g_l*(v+(Dt/2)*fc7_2-Evaz) +G*(vA-(v+(Dt/2)*fc7_2))-G*(v+(Dt/2)*fc7_2-vP);
        fc7_4 = -A*g_Na*(v+Dt*fc7_3-Ena)-A*g_k*(v+Dt*fc7_3-Ek)-A*g_l*(v+Dt*fc7_3-Evaz) +G*(vA-(v+Dt*fc7_3))-G*(v+Dt*fc7_3-vP);
        v7(t) = v7(t-1) + (Dt/(6*Cm))*(fc7_1+2*fc7_2+2*fc7_3+fc7_4);
        
        %Compartimento 8
        g_Na = g_Na8(t-1);
        g_k  = g_k8(t-1);
        v = v8(t-1);
        vA = v7(t-1); % Anterior
        vP = v9(t-1); % Posterior
        fc8_1 = -A*g_Na*(v-Ena)-A*g_k*(v-Ek)-A*g_l*(v-Evaz) +G*(vA-v)-G*(v-vP);
        fc8_2 = -A*g_Na*(v+(Dt/2)*fc8_1-Ena)-A*g_k*(v+(Dt/2)*fc8_1-Ek)-A*g_l*(v+(Dt/2)*fc8_1-Evaz) +G*(vA-(v+(Dt/2)*fc8_1))-G*(v+(Dt/2)*fc8_1-vP);
        fc8_3 = -A*g_Na*(v+(Dt/2)*fc8_2-Ena)-A*g_k*(v+(Dt/2)*fc8_2-Ek)-A*g_l*(v+(Dt/2)*fc8_2-Evaz) +G*(vA-(v+(Dt/2)*fc8_2))-G*(v+(Dt/2)*fc8_2-vP);
        fc8_4 = -A*g_Na*(v+Dt*fc8_3-Ena)-A*g_k*(v+Dt*fc8_3-Ek)-A*g_l*(v+Dt*fc8_3-Evaz) +G*(vA-(v+Dt*fc8_3))-G*(v+Dt*fc8_3-vP);
        v8(t) = v8(t-1) + (Dt/(6*Cm))*(fc8_1+2*fc8_2+2*fc8_3+fc8_4);
        
        %Compartimento 9
        g_Na = g_Na9(t-1);
        g_k  = g_k9(t-1);
        v = v9(t-1);
        vA = v8(t-1); % Anterior
        vP = v10(t-1); % Posterior
        fc9_1 = -A*g_Na*(v-Ena)-A*g_k*(v-Ek)-A*g_l*(v-Evaz) +G*(vA-v)-G*(v-vP);
        fc9_2 = -A*g_Na*(v+(Dt/2)*fc9_1-Ena)-A*g_k*(v+(Dt/2)*fc9_1-Ek)-A*g_l*(v+(Dt/2)*fc9_1-Evaz) +G*(vA-(v+(Dt/2)*fc9_1))-G*(v+(Dt/2)*fc9_1-vP);
        fc9_3 = -A*g_Na*(v+(Dt/2)*fc9_2-Ena)-A*g_k*(v+(Dt/2)*fc9_2-Ek)-A*g_l*(v+(Dt/2)*fc9_2-Evaz) +G*(vA-(v+(Dt/2)*fc9_2))-G*(v+(Dt/2)*fc9_2-vP);
        fc9_4 = -A*g_Na*(v+Dt*fc9_3-Ena)-A*g_k*(v+Dt*fc9_3-Ek)-A*g_l*(v+Dt*fc9_3-Evaz) +G*(vA-(v+Dt*fc9_3))-G*(v+Dt*fc9_3-vP);
        v9(t) = v9(t-1) + (Dt/(6*Cm))*(fc9_1+2*fc9_2+2*fc9_3+fc9_4);
        
        %Compartimento 10
        g_Na = g_Na10(t-1);
        g_k  = g_k10(t-1);
        v = v10(t-1);
        vA = v9(t-1); % Anterior
        fc10_1 = -A*g_Na*(v-Ena)-A*g_k*(v-Ek)-A*g_l*(v-Evaz)+G*(vA-v);
        fc10_2 = -A*g_Na*(v+(Dt/2)*fc10_1-Ena)-A*g_k*(v+(Dt/2)*fc10_1-Ek)-A*g_l*(v+(Dt/2)*fc10_1-Evaz) +G*(vA-(v+(Dt/2)*fc10_1));
        fc10_3 = -A*g_Na*(v+(Dt/2)*fc10_2-Ena)-A*g_k*(v+(Dt/2)*fc10_2-Ek)-A*g_l*(v+(Dt/2)*fc10_2-Evaz) +G*(vA-(v+(Dt/2)*fc10_2));
        fc10_4 = -A*g_Na*(v+Dt*fc10_3-Ena)-A*g_k*(v+Dt*fc10_3-Ek)-A*g_l*(v+Dt*fc10_3-Evaz) +G*(vA-(v+Dt*fc10_3));
        v10(t) = v10(t-1) + (Dt/(6*Cm))*(fc10_1+2*fc10_2+2*fc10_3+fc10_4);

        
        %% Cálculo: Probabilidade n por Runge-Kutta
        %Compartimento 1
        v = v1(t-1);
        n = n1(t-1);
        n1_1=      ((1.0-n)*alpha_n(v)-n*beta_n(v));
        n1_2=      ((1.0-n-(Dt/2)*n1_1)*alpha_n(v)-(n+(Dt/2)*n1_1)*beta_n(v));
        n1_3=      ((1.0-n-(Dt/2)*n1_2)*alpha_n(v)-(n+(Dt/2)*n1_2)*beta_n(v));
        n1_4=      ((1.0-n-(Dt)*n1_3)*alpha_n(v)-(n+(Dt)*n1_3)*beta_n(v));
        n1(t) = n1(t-1) + (Dt/6)*(n1_1+2*n1_2+2*n1_3+n1_4);
        %Compartimento 2
        v = v2(t-1);
        n = n2(t-1);
        n2_1=      ((1.0-n)*alpha_n(v)-n*beta_n(v));
        n2_2=      ((1.0-n-(Dt/2)*n2_1)*alpha_n(v)-(n+(Dt/2)*n2_1)*beta_n(v));
        n2_3=      ((1.0-n-(Dt/2)*n2_2)*alpha_n(v)-(n+(Dt/2)*n2_2)*beta_n(v));
        n2_4=      ((1.0-n-(Dt)*n2_3)*alpha_n(v)-(n+(Dt)*n2_3)*beta_n(v));
        n2(t) = n2(t-1) + (Dt/6)*(n2_1+2*n2_2+2*n2_3+n2_4);        
        %Compartimento 3
        v = v3(t-1);
        n = n3(t-1);
        n3_1=      ((1.0-n)*alpha_n(v)-n*beta_n(v));
        n3_2=      ((1.0-n-(Dt/2)*n3_1)*alpha_n(v)-(n+(Dt/2)*n3_1)*beta_n(v));
        n3_3=      ((1.0-n-(Dt/2)*n3_2)*alpha_n(v)-(n+(Dt/2)*n3_2)*beta_n(v));
        n3_4=      ((1.0-n-(Dt)*n3_3)*alpha_n(v)-(n+(Dt)*n3_3)*beta_n(v));
        n3(t) = n3(t-1) + (Dt/6)*(n3_1+2*n3_2+2*n3_3+n3_4);
        %Compartimento 4
        v = v4(t-1);
        n = n4(t-1);
        n4_1=      ((1.0-n)*alpha_n(v)-n*beta_n(v));
        n4_2=      ((1.0-n-(Dt/2)*n4_1)*alpha_n(v)-(n+(Dt/2)*n4_1)*beta_n(v));
        n4_3=      ((1.0-n-(Dt/2)*n4_2)*alpha_n(v)-(n+(Dt/2)*n4_2)*beta_n(v));
        n4_4=      ((1.0-n-(Dt)*n4_3)*alpha_n(v)-(n+(Dt)*n4_3)*beta_n(v));
        n4(t) = n4(t-1) + (Dt/6)*(n4_1+2*n4_2+2*n4_3+n4_4);
        %Compartimento 5
        v = v5(t-1);
        n = n5(t-1);
        n5_1=      ((1.0-n)*alpha_n(v)-n*beta_n(v));
        n5_2=      ((1.0-n-(Dt/2)*n5_1)*alpha_n(v)-(n+(Dt/2)*n5_1)*beta_n(v));
        n5_3=      ((1.0-n-(Dt/2)*n5_2)*alpha_n(v)-(n+(Dt/2)*n5_2)*beta_n(v));
        n5_4=      ((1.0-n-(Dt)*n5_3)*alpha_n(v)-(n+(Dt)*n5_3)*beta_n(v));
        n5(t) = n5(t-1) + (Dt/6)*(n5_1+2*n5_2+2*n5_3+n5_4);
        %Compartimento 6
        v = v6(t-1);
        n = n6(t-1);
        n6_1=      ((1.0-n)*alpha_n(v)-n*beta_n(v));
        n6_2=      ((1.0-n-(Dt/2)*n6_1)*alpha_n(v)-(n+(Dt/2)*n6_1)*beta_n(v));
        n6_3=      ((1.0-n-(Dt/2)*n6_2)*alpha_n(v)-(n+(Dt/2)*n6_2)*beta_n(v));
        n6_4=      ((1.0-n-(Dt)*n6_3)*alpha_n(v)-(n+(Dt)*n6_3)*beta_n(v));
        n6(t) = n6(t-1) + (Dt/6)*(n6_1+2*n6_2+2*n6_3+n6_4);       
        %Compartimento 7
        v = v7(t-1);
        n = n7(t-1);
        n7_1=      ((1.0-n)*alpha_n(v)-n*beta_n(v));
        n7_2=      ((1.0-n-(Dt/2)*n7_1)*alpha_n(v)-(n+(Dt/2)*n7_1)*beta_n(v));
        n7_3=      ((1.0-n-(Dt/2)*n7_2)*alpha_n(v)-(n+(Dt/2)*n7_2)*beta_n(v));
        n7_4=      ((1.0-n-(Dt)*n7_3)*alpha_n(v)-(n+(Dt)*n7_3)*beta_n(v));
        n7(t) = n7(t-1) + (Dt/6)*(n7_1+2*n7_2+2*n7_3+n7_4);       
        %Compartimento 8
        v = v8(t-1);
        n = n8(t-1);
        n8_1=      ((1.0-n)*alpha_n(v)-n*beta_n(v));
        n8_2=      ((1.0-n-(Dt/2)*n8_1)*alpha_n(v)-(n+(Dt/2)*n8_1)*beta_n(v));
        n8_3=      ((1.0-n-(Dt/2)*n8_2)*alpha_n(v)-(n+(Dt/2)*n8_2)*beta_n(v));
        n8_4=      ((1.0-n-(Dt)*n8_3)*alpha_n(v)-(n+(Dt)*n8_3)*beta_n(v));
        n8(t) = n8(t-1) + (Dt/6)*(n8_1+2*n8_2+2*n8_3+n8_4); 
        %Compartimento 9
        v = v9(t-1);
        n = n9(t-1);
        n9_1=      ((1.0-n)*alpha_n(v)-n*beta_n(v));
        n9_2=      ((1.0-n-(Dt/2)*n9_1)*alpha_n(v)-(n+(Dt/2)*n9_1)*beta_n(v));
        n9_3=      ((1.0-n-(Dt/2)*n9_2)*alpha_n(v)-(n+(Dt/2)*n9_2)*beta_n(v));
        n9_4=      ((1.0-n-(Dt)*n9_3)*alpha_n(v)-(n+(Dt)*n9_3)*beta_n(v));
        n9(t) = n9(t-1) + (Dt/6)*(n9_1+2*n9_2+2*n9_3+n9_4); 
        %Compartimento 10
        v = v10(t-1);
        n = n10(t-1);
        n10_1=      ((1.0-n)*alpha_n(v)-n*beta_n(v));
        n10_2=      ((1.0-n-(Dt/2)*n10_1)*alpha_n(v)-(n+(Dt/2)*n10_1)*beta_n(v));
        n10_3=      ((1.0-n-(Dt/2)*n10_2)*alpha_n(v)-(n+(Dt/2)*n10_2)*beta_n(v));
        n10_4=      ((1.0-n-(Dt)*n10_3)*alpha_n(v)-(n+(Dt)*n10_3)*beta_n(v));
        n10(t) = n10(t-1) + (Dt/6)*(n10_1+2*n10_2+2*n10_3+n10_4);        
        
        
        %% Cálculo: Probabilidade m por Runge-Kutta
        %Compartimento 1
        v = v1(t-1);
        m = m1(t-1);
        m1_1=       ((1.0-m)*alpha_m(v)-m*beta_m(v));
        m1_2=       ((1.0-m-(Dt/2)*m1_1)*alpha_m(v)-(m+(Dt/2)*m1_1)*beta_m(v));
        m1_3=       ((1.0-m-(Dt/2)*m1_2)*alpha_m(v)-(m+(Dt/2)*m1_2)*beta_m(v));
        m1_4=       ((1.0-m-Dt*m1_3)*alpha_m(v)-(m+Dt*m1_3)*beta_m(v));
        m1(t) = m1(t-1) + (Dt/6)*(m1_1+2*m1_2+2*m1_3+m1_4); 
        %Compartimento 2
        v = v2(t-1);
        m = m2(t-1);
        m2_1=       ((1.0-m)*alpha_m(v)-m*beta_m(v));
        m2_2=       ((1.0-m-(Dt/2)*m2_1)*alpha_m(v)-(m+(Dt/2)*m2_1)*beta_m(v));
        m2_3=       ((1.0-m-(Dt/2)*m2_2)*alpha_m(v)-(m+(Dt/2)*m2_2)*beta_m(v));
        m2_4=       ((1.0-m-Dt*m2_3)*alpha_m(v)-(m+Dt*m2_3)*beta_m(v));
        m2(t) = m2(t-1) + (Dt/6)*(m2_1+2*m2_2+2*m2_3+m2_4); 
        %Compartimento 3
        v = v3(t-1);
        m = m3(t-1);
        m3_1=       ((1.0-m)*alpha_m(v)-m*beta_m(v));
        m3_2=       ((1.0-m-(Dt/2)*m3_1)*alpha_m(v)-(m+(Dt/2)*m3_1)*beta_m(v));
        m3_3=       ((1.0-m-(Dt/2)*m3_2)*alpha_m(v)-(m+(Dt/2)*m3_2)*beta_m(v));
        m3_4=       ((1.0-m-Dt*m3_3)*alpha_m(v)-(m+Dt*m3_3)*beta_m(v));
        m3(t) = m3(t-1) + (Dt/6)*(m3_1+2*m3_2+2*m3_3+m3_4); 
        %Compartimento 4
        v = v4(t-1);
        m = m4(t-1);
        m4_1=       ((1.0-m)*alpha_m(v)-m*beta_m(v));
        m4_2=       ((1.0-m-(Dt/2)*m4_1)*alpha_m(v)-(m+(Dt/2)*m4_1)*beta_m(v));
        m4_3=       ((1.0-m-(Dt/2)*m4_2)*alpha_m(v)-(m+(Dt/2)*m4_2)*beta_m(v));
        m4_4=       ((1.0-m-Dt*m4_3)*alpha_m(v)-(m+Dt*m4_3)*beta_m(v));
        m4(t) = m4(t-1) + (Dt/6)*(m4_1+2*m4_2+2*m4_3+m4_4); 
        %Compartimento 5
        v = v5(t-1);
        m = m5(t-1);
        m5_1=       ((1.0-m)*alpha_m(v)-m*beta_m(v));
        m5_2=       ((1.0-m-(Dt/2)*m5_1)*alpha_m(v)-(m+(Dt/2)*m5_1)*beta_m(v));
        m5_3=       ((1.0-m-(Dt/2)*m5_2)*alpha_m(v)-(m+(Dt/2)*m5_2)*beta_m(v));
        m5_4=       ((1.0-m-Dt*m5_3)*alpha_m(v)-(m+Dt*m5_3)*beta_m(v));
        m5(t) = m5(t-1) + (Dt/6)*(m5_1+2*m5_2+2*m5_3+m5_4); 
        %Compartimento 6
        v = v6(t-1);
        m = m6(t-1);
        m6_1=       ((1.0-m)*alpha_m(v)-m*beta_m(v));
        m6_2=       ((1.0-m-(Dt/2)*m6_1)*alpha_m(v)-(m+(Dt/2)*m6_1)*beta_m(v));
        m6_3=       ((1.0-m-(Dt/2)*m6_2)*alpha_m(v)-(m+(Dt/2)*m6_2)*beta_m(v));
        m6_4=       ((1.0-m-Dt*m6_3)*alpha_m(v)-(m+Dt*m6_3)*beta_m(v));
        m6(t) = m6(t-1) + (Dt/6)*(m6_1+2*m6_2+2*m6_3+m6_4);
        %Compartimento 7
        v = v7(t-1);
        m = m7(t-1);
        m7_1=       ((1.0-m)*alpha_m(v)-m*beta_m(v));
        m7_2=       ((1.0-m-(Dt/2)*m7_1)*alpha_m(v)-(m+(Dt/2)*m7_1)*beta_m(v));
        m7_3=       ((1.0-m-(Dt/2)*m7_2)*alpha_m(v)-(m+(Dt/2)*m7_2)*beta_m(v));
        m7_4=       ((1.0-m-Dt*m7_3)*alpha_m(v)-(m+Dt*m7_3)*beta_m(v));
        m7(t) = m7(t-1) + (Dt/6)*(m7_1+2*m7_2+2*m7_3+m7_4);
        %Compartimento 8
        v = v8(t-1);
        m = m8(t-1);
        m8_1=       ((1.0-m)*alpha_m(v)-m*beta_m(v));
        m8_2=       ((1.0-m-(Dt/2)*m8_1)*alpha_m(v)-(m+(Dt/2)*m8_1)*beta_m(v));
        m8_3=       ((1.0-m-(Dt/2)*m8_2)*alpha_m(v)-(m+(Dt/2)*m8_2)*beta_m(v));
        m8_4=       ((1.0-m-Dt*m8_3)*alpha_m(v)-(m+Dt*m8_3)*beta_m(v));
        m8(t) = m8(t-1) + (Dt/6)*(m8_1+2*m8_2+2*m8_3+m8_4);     
        %Compartimento 9
        v = v9(t-1);
        m = m9(t-1);
        m9_1=       ((1.0-m)*alpha_m(v)-m*beta_m(v));
        m9_2=       ((1.0-m-(Dt/2)*m9_1)*alpha_m(v)-(m+(Dt/2)*m9_1)*beta_m(v));
        m9_3=       ((1.0-m-(Dt/2)*m9_2)*alpha_m(v)-(m+(Dt/2)*m9_2)*beta_m(v));
        m9_4=       ((1.0-m-Dt*m9_3)*alpha_m(v)-(m+Dt*m9_3)*beta_m(v));
        m9(t) = m9(t-1) + (Dt/6)*(m9_1+2*m9_2+2*m9_3+m9_4);  
        %Compartimento 10
        v = v10(t-1);
        m = m10(t-1);
        m10_1=       ((1.0-m)*alpha_m(v)-m*beta_m(v));
        m10_2=       ((1.0-m-(Dt/2)*m10_1)*alpha_m(v)-(m+(Dt/2)*m10_1)*beta_m(v));
        m10_3=       ((1.0-m-(Dt/2)*m10_2)*alpha_m(v)-(m+(Dt/2)*m10_2)*beta_m(v));
        m10_4=       ((1.0-m-Dt*m10_3)*alpha_m(v)-(m+Dt*m10_3)*beta_m(v));
        m10(t) = m10(t-1) + (Dt/6)*(m10_1+2*m10_2+2*m10_3+m10_4);  
                
        %% Cálculo: Probabilidade n por Runge-Kutta
        %Compartimento 1
        v = v1(t-1);
        h = h1(t-1);
        h1_1=       ((1.0-h)*alpha_h(v)-(h)*beta_h(v));
        h1_2=       ((1.0-h-(Dt/2)*h1_1)*alpha_h(v)-(h+(Dt/2)*h1_1)*beta_h(v));
        h1_3=       ((1.0-h-(Dt/2)*h1_2)*alpha_h(v)-(h+(Dt/2)*h1_2)*beta_h(v));
        h1_4=       ((1.0-h-Dt*h1_3)*alpha_h(v)-(h+Dt*h1_3)*beta_h(v));
        h1(t) = h1(t-1) + (Dt/6)*(h1_1+2*h1_2+2*h1_3+h1_4);
        %Compartimento 2
        v = v2(t-1);
        h = h2(t-1);
        h2_1=       ((1.0-h)*alpha_h(v)-(h)*beta_h(v));
        h2_2=       ((1.0-h-(Dt/2)*h2_1)*alpha_h(v)-(h+(Dt/2)*h2_1)*beta_h(v));
        h2_3=       ((1.0-h-(Dt/2)*h2_2)*alpha_h(v)-(h+(Dt/2)*h2_2)*beta_h(v));
        h2_4=       ((1.0-h-Dt*h2_3)*alpha_h(v)-(h+Dt*h2_3)*beta_h(v));
        h2(t) = h2(t-1) + (Dt/6)*(h2_1+2*h2_2+2*h2_3+h2_4);
        %Compartimento 3
        v = v3(t-1);
        h = h3(t-1);
        h3_1=       ((1.0-h)*alpha_h(v)-(h)*beta_h(v));
        h3_2=       ((1.0-h-(Dt/2)*h3_1)*alpha_h(v)-(h+(Dt/2)*h3_1)*beta_h(v));
        h3_3=       ((1.0-h-(Dt/2)*h3_2)*alpha_h(v)-(h+(Dt/2)*h3_2)*beta_h(v));
        h3_4=       ((1.0-h-Dt*h3_3)*alpha_h(v)-(h+Dt*h3_3)*beta_h(v));
        h3(t) = h3(t-1) + (Dt/6)*(h3_1+2*h3_2+2*h3_3+h3_4);
        %Compartimento 4
        v = v4(t-1);
        h = h4(t-1);
        h4_1=       ((1.0-h)*alpha_h(v)-(h)*beta_h(v));
        h4_2=       ((1.0-h-(Dt/2)*h4_1)*alpha_h(v)-(h+(Dt/2)*h4_1)*beta_h(v));
        h4_3=       ((1.0-h-(Dt/2)*h4_2)*alpha_h(v)-(h+(Dt/2)*h4_2)*beta_h(v));
        h4_4=       ((1.0-h-Dt*h4_3)*alpha_h(v)-(h+Dt*h4_3)*beta_h(v));
        h4(t) = h4(t-1) + (Dt/6)*(h4_1+2*h4_2+2*h4_3+h4_4);
        %Compartimento 5
        v = v5(t-1);
        h = h5(t-1);
        h5_1=       ((1.0-h)*alpha_h(v)-(h)*beta_h(v));
        h5_2=       ((1.0-h-(Dt/2)*h5_1)*alpha_h(v)-(h+(Dt/2)*h5_1)*beta_h(v));
        h5_3=       ((1.0-h-(Dt/2)*h5_2)*alpha_h(v)-(h+(Dt/2)*h5_2)*beta_h(v));
        h5_4=       ((1.0-h-Dt*h5_3)*alpha_h(v)-(h+Dt*h5_3)*beta_h(v));
        h5(t) = h5(t-1) + (Dt/6)*(h5_1+2*h5_2+2*h5_3+h5_4);
        %Compartimento 6
        v = v6(t-1);
        h = h6(t-1);
        h6_1=       ((1.0-h)*alpha_h(v)-(h)*beta_h(v));
        h6_2=       ((1.0-h-(Dt/2)*h6_1)*alpha_h(v)-(h+(Dt/2)*h6_1)*beta_h(v));
        h6_3=       ((1.0-h-(Dt/2)*h6_2)*alpha_h(v)-(h+(Dt/2)*h6_2)*beta_h(v));
        h6_4=       ((1.0-h-Dt*h6_3)*alpha_h(v)-(h+Dt*h6_3)*beta_h(v));
        h6(t) = h6(t-1) + (Dt/6)*(h6_1+2*h6_2+2*h6_3+h6_4);
        %Compartimento 7
        v = v7(t-1);
        h = h7(t-1);
        h7_1=       ((1.0-h)*alpha_h(v)-(h)*beta_h(v));
        h7_2=       ((1.0-h-(Dt/2)*h7_1)*alpha_h(v)-(h+(Dt/2)*h7_1)*beta_h(v));
        h7_3=       ((1.0-h-(Dt/2)*h7_2)*alpha_h(v)-(h+(Dt/2)*h7_2)*beta_h(v));
        h7_4=       ((1.0-h-Dt*h7_3)*alpha_h(v)-(h+Dt*h7_3)*beta_h(v));
        h7(t) = h7(t-1) + (Dt/6)*(h7_1+2*h7_2+2*h7_3+h7_4);
        %Compartimento 8
        v = v8(t-1);
        h = h8(t-1);
        h8_1=       ((1.0-h)*alpha_h(v)-(h)*beta_h(v));
        h8_2=       ((1.0-h-(Dt/2)*h8_1)*alpha_h(v)-(h+(Dt/2)*h8_1)*beta_h(v));
        h8_3=       ((1.0-h-(Dt/2)*h8_2)*alpha_h(v)-(h+(Dt/2)*h8_2)*beta_h(v));
        h8_4=       ((1.0-h-Dt*h8_3)*alpha_h(v)-(h+Dt*h8_3)*beta_h(v));
        h8(t) = h8(t-1) + (Dt/6)*(h8_1+2*h8_2+2*h8_3+h8_4);
        %Compartimento 9
        v = v9(t-1);
        h = h9(t-1);
        h9_1=       ((1.0-h)*alpha_h(v)-(h)*beta_h(v));
        h9_2=       ((1.0-h-(Dt/2)*h9_1)*alpha_h(v)-(h+(Dt/2)*h9_1)*beta_h(v));
        h9_3=       ((1.0-h-(Dt/2)*h9_2)*alpha_h(v)-(h+(Dt/2)*h9_2)*beta_h(v));
        h9_4=       ((1.0-h-Dt*h9_3)*alpha_h(v)-(h+Dt*h9_3)*beta_h(v));
        h9(t) = h9(t-1) + (Dt/6)*(h9_1+2*h9_2+2*h9_3+h9_4);
        %Compartimento 10
        v = v10(t-1);
        h = h10(t-1);
        h10_1=       ((1.0-h)*alpha_h(v)-(h)*beta_h(v));
        h10_2=       ((1.0-h-(Dt/2)*h10_1)*alpha_h(v)-(h+(Dt/2)*h10_1)*beta_h(v));
        h10_3=       ((1.0-h-(Dt/2)*h10_2)*alpha_h(v)-(h+(Dt/2)*h10_2)*beta_h(v));
        h10_4=       ((1.0-h-Dt*h10_3)*alpha_h(v)-(h+Dt*h10_3)*beta_h(v));
        h10(t) = h10(t-1) + (Dt/6)*(h10_1+2*h10_2+2*h10_3+h10_4);      

        
        tv(t)=t*Dt;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS Required from Previous Lits %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Funcao para o calculo de alpha_n(V) 
function aln=alpha_n(V)
    % Parte superior e inferior da equacao para alpha_n
    s = 0.01*(10-V);
    i = exp(1-0.1*V) - 1.0;
    %verificar a condicao para nao haver indeterminacao
    if ((s==0) && (i==0))
        s= -0.01;
        i= -0.1*exp(1-0.1*V);
    end
    aln=s/i;
end

%% Funcao para o calculo de alpha_m(V) 
function aln=alpha_m(V)
    % Parte superior e inferior da equacao para alpha_n
    s = 0.1*(25-V);
    i = exp(2.5-0.1*V)-1.0;
    %verificar a condicao para nao haver indeterminacao
    if ( (s==0) && (i==0))
        s= -0.1;
        i= -0.1*exp(2.5-0.1*V);
    end
    aln=s/i;
end

%% Funcao para o calculo de alpha_h(V) 
function aln=alpha_h(V)
    % Para esse canal nao ha a necessidade de verificacoes
    aln=0.07*exp(-V/20.0);
end

%% Funcao para o calculo de beta_n(V) 
function aln=beta_n(V)
    % Para esse canal nao ha a necessidade de verificacoes
    aln=0.125*exp(-V/80.0);
end

%% Funcao para o calculo de beta_m(V) 
function aln=beta_m(V)
    % Para esse canal nao ha a necessidade de verificacoes
    aln=4.0*exp(-V/18.0);
end

%% Funcao para o calculo de beta_h(V) 
function aln=beta_h(V)
    % Nao ha problema para o calculo do valor de beta_h
    s = 1.0;
    i = exp(3.0-0.1*V) + 1;
    aln=s/i;
end
