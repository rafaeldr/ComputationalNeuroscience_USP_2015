% Quarta lista de Exercícios - Introdução à Neurociência Computacional
% Aluno: Rafael Delalibera Rodrigues
% Programa de Pós-Graduação em Computação Aplicada [FFCLRP-USP Ribeirão Preto] - Mestrado

% Exercício 4: Hodgkin-Huxley
% 4.a) Análise da Variável de Ativação do Sódio

% Assinatura da Função: [v,tv,n,m,h,ninf,minf,hinf,ntau,mtau,htau,g_Na,g_k,J] = HH_RKUTA4_inf(v0,ti,tf,Dt,ti_J,tf_J,Jin)

% Potencial de Repouso = -70mV
% Tempo do experimento: 0 -> 200ms
% Corrente injetada: 10uA/cm2; 10->160ms
[v,t,n,m,h,ninf,minf,hinf,ntau,mtau,htau,g_Na,g_k] = HH_RKUTA4_inf(-70,0,200,0.05,10,160,10);

% V x t
figure(1);
plot(t,v);
title('V x t');
xlabel('Tempo (ms)');
ylabel('mV');

% m(t) e minf(t)
figure(2);
hold on;
plot(t,m,'blue');
plot(t,minf,'-red');
title('m x t VS minf x t');
legend('m(t)','minf(t)');
xlabel('Tempo (ms)');
ylabel('Ativação');
hold off;


% Estudo adicional das constantes temporáis das variáveis do modelo HH
figure(3);
s1 = subplot(2,2,[1 2]);
hold on;
% Tau N: Constante de Tempo da Variável de Ativação do Potássio
p1 = plot(v,ntau,'blue');
% Tau M: Constante de Tempo da Variável de Ativação do Sódio
p1 = plot(v,mtau,'black');
% Tau H: Constante de Tempo da Variável de Inativação do Sódio
p1 = plot(v,htau,'red');
legend('n','m','h');
xlabel('V (mV)');
ylabel('Tau (ms)');
set(s1,'xlim',[-80 40]);
hold off;

% n(t) e ninf(t)
s2 = subplot(2,2,3);
hold on;
plot(t,n,'blue');
plot(t,ninf,'-red');
title('n x t VS ninf x t');
legend('n(t)','ninf(t)');
xlabel('Tempo (ms)');
ylabel('Ativação');
hold off;

% h(t) e hinf(t)
s3 = subplot(2,2,4);
hold on;
plot(t,h,'blue');
plot(t,hinf,'-red');
title('h x t VS hinf x t');
legend('h(t)','hinf(t)');
xlabel('Tempo (ms)');
ylabel('Ativação');
hold off;


% Exercício 4: Hodgkin-Huxley
% 4.b) Gráfico h X n

% Gráfico h X n
figure(4);
plot(n,h);
xlabel('n');
ylabel('h');
title('h X n');


% Exercício 5: Modelo de Connor-Stevens
% 5.a) A-Current

% Figura 6.1 (Gráfico B)
[v,t,n,m,h,a,b,ninf,minf,hinf,ainf,binf,ntau,mtau,htau,atau,btau,g_Na,g_k,J] = CS_RKUTA4(-68,0,100,0.01,0,100,10,0,0);
% V x t
figure(5);
plot(t,v);
title('V x t');
xlabel('t (ms)');
ylabel('V (mV)');
legend(sprintf('%.2f uA/cm²',10));

% Figura 6.1 (Gráfico D)
[v,t,n,m,h,a,b,ninf,minf,hinf,ainf,binf,ntau,mtau,htau,atau,btau,g_Na,g_k,J] = CS_RKUTA4(0,0,200,0.01,0,200,-50,1,10);
% V x t
figure(6);
plot(t,v);
title('V x t');
xlabel('Tempo (ms)');
ylabel('mV');
legend(sprintf('%.2f uA/cm²   ',-50,10));


% Exercício 5: Modelo de Connor-Stevens
% 5.b) Corrente Transiente de Cálcio

[v,t,n,m,h,a,b,ninf,minf,hinf,ainf,binf,ntau,mtau,htau,atau,btau,g_Na,g_k,J] = CS_CAT_RKUTA4(0,0,350,0.01,0,350,-50,1,8);
% V x t
figure(7);
plot(t,v);
title('V x t');
xlabel('Tempo (ms)');
ylabel('mV');


% Exercício 6: Modelo de Traub-Miles
% 6.a) Corrente Adaptativa z

[v,t,n,m,h,g_Na,g_k,J] = HH_RKUTA4(0,0,400,0.05,50,250,5);
figure(8);
s1 = subplot(2,2,1);
plot(t,v);
title('[Hodgkin-Huxley Model] V x t');
xlabel('Tempo (ms)');
ylabel('V(mV)');

s2 = subplot(2,2,3);
plot(t,J);
title('[Hodgkin-Huxley Model] J x t');
xlabel('Tempo (ms)');
ylabel('J(uA/cm²)');
set(s2,'ylim',[-10 40]);
clear;
[v,t,n,m,h,z,g_Na,g_k,J] = TM_RKUTA4(0,0,400,0.01,50,250,5);
s3 = subplot(2,2,2);
plot(t,v);
title('[Traub-Miles Model] V x t');
xlabel('Tempo (ms)');
ylabel('V(mV)');

s4 = subplot(2,2,4);
plot(t,J);
title('[Traub-Miles Model] J x t');
xlabel('Tempo (ms)');
ylabel('J(uA/cm²)');
set(s4,'ylim',[-10 40]);
clear;

[v,t,n,m,h,z,g_Na,g_k,J] = TM_RKUTA4(0,0,400,0.01,50,250,5);
figure(9);
plot(t,z);
title('[Traub-Miles Model] z x t');
xlabel('Tempo (ms)');
ylabel('z');


% Exercício 6: Modelo de Traub-Miles
% 6.b) Corrente Adaptativa: Bi-estabilidade

[v,t,n,m,h,z,g_Na,g_k,J] = TM_RKUTA4_v2(0,0,1000,0.01,150,200,5,350,650,-10);
figure(10);
s1 = subplot(2,1,1);
plot(t,v);
title('[Traub-Miles Model] V x t');
xlabel('Tempo (ms)');
ylabel('V(mV)');

s2 = subplot(2,1,2);
plot(t,J);
title('[Traub-Miles Model] J x t');
xlabel('Tempo (ms)');
ylabel('J(uS/cm²)');
set(s2,'ylim',[-20 20]);


% Exercício 6: Modelo de Traub-Miles
% 6.c) Corrente de Potássio Dependente de Cálcio

[v,t,n,m,h,z,g_Na,g_k,J] = TM_RKUTA4(0,0,400,0.01,50,250,5);
figure(11);
s1 = subplot(2,2,1);
plot(t,v);
title('[Traub-Miles Model] V x t');
xlabel('Tempo (ms)');
ylabel('V(mV)');

s2 = subplot(2,2,3);
plot(t,J);
title('[Traub-Miles Model] J x t');
xlabel('Tempo (ms)');
ylabel('J(uA/cm²)');
set(s2,'ylim',[-10 40]);
clear;
[v,t,n,m,h,ca,mca,minf,g_Na,g_k,J] = TM_RKUTA4_v3(0,0,400,0.01,50,250,5);
s3 = subplot(2,2,2);
plot(t,v);
title('[Corrente de Potássio Dependente de Cálcio] V x t');
xlabel('Tempo (ms)');
ylabel('V(mV)');

s4 = subplot(2,2,4);
plot(t,J);
title('[Corrente de Potássio Dependente de Cálcio] J x t');
xlabel('Tempo (ms)');
ylabel('J(uA/cm²)');
set(s4,'ylim',[-10 40]);
clear;

[v,t,n,m,h,ca,mca,minf,g_Na,g_k,J] = TM_RKUTA4_v3(0,0,400,0.01,50,250,5);
figure(12);
plot(t,mca);
title('[Corrente de Potássio Dependente de Cálcio] m(Ca) x t');
xlabel('Tempo (ms)');
ylabel('m(Ca)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS (Must be defined at the end of the script) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exercício 4: Hodgkin-Huxley                    [FUNCTION]
% 4.a) Análise da Variável de Ativação do Sódio  [FUNCTION]

function [v,tv,n,m,h,ninf,minf,hinf,ntau,mtau,htau,g_Na,g_k,J] = HH_RKUTA4_inf(v0,ti,tf,Dt,ti_J,tf_J,Jin)
    
	%parametros eletricos da membrana
    Ena  = 115;
    Ek   = -12;
    Evaz = 10.16; 
    GNa  = 120;
    GK   = 36;
    GL   = 0.3; 
    cm   = 1;
    
    %Dimensão do vetores (conforme granularidade temporal)
    len = ceil((tf-ti)/Dt);
    
    %Vetores das Condutâncias 
    g_Na = zeros(len,1);
    g_k  = zeros(len,1);
    g_l  = GL;
    
    %Vetor de Correte Injetada
    J=zeros(len);
    for i=1:len
        % Decisao quanto a entrada de corrente
       if ((i*Dt>ti_J) && (i*Dt<tf_J))
            J(i)=Jin;
       else
            J(i)=0;
       end
    end
    
    %Vetores de Probabilidades (Gates)
    n  = zeros(len,1);
    m  = zeros(len,1);
    h  = zeros(len,1);
    
    %Inicialização de Vetores
    ninf = zeros(len,1);
    minf = zeros(len,1);
    hinf = zeros(len,1);
    ntau = zeros(len,1);
    mtau = zeros(len,1);
    htau = zeros(len,1);
    
    %Vetor Voltagens e Tempo
    v    = zeros(len,1);
    tv   = zeros(len,1);
    v(1) = v0+(-v0);
    tv(1)= Dt;
    
	%Probabilidade de Abertura dos Canais (para célula em repouso)
    n(1) = alpha_n(v(1))/(alpha_n(v(1))+beta_n(v(1)));
    m(1) = alpha_m(v(1))/(alpha_m(v(1))+beta_m(v(1)));
    h(1) = alpha_h(v(1))/(alpha_h(v(1))+beta_h(v(1)));
    
    %Taxas de Mudanças (para célula em repouso)
    ninf(1) = alpha_n(v(1))/(alpha_n(v(1))+beta_n(v(1)));
    minf(1) = alpha_m(v(1))/(alpha_m(v(1))+beta_m(v(1)));
    hinf(1) = alpha_h(v(1))/(alpha_h(v(1))+beta_h(v(1)));
    
    %Constantes de tempo (para célula em repouso)
    ntau(1) = 1.0/(alpha_n(v(1))+beta_n(v(1)));
    mtau(1) = 1.0/(alpha_m(v(1))+beta_m(v(1)));
    htau(1) = 1.0/(alpha_h(v(1))+beta_h(v(1)));
    
    for t=2:len
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
        
        % Cálculo: Taxas de mudança: ninf, minf, hinf
        ninf(t) = alpha_n(v(t))/(alpha_n(v(t))+beta_n(v(t)));
        minf(t) = alpha_m(v(t))/(alpha_m(v(t))+beta_m(v(t)));
        hinf(t) = alpha_h(v(t))/(alpha_h(v(t))+beta_h(v(t)));
        
        % Cálculo: Constantes de Tempo
        ntau(t) = 1.0/(alpha_n(v(t))+beta_n(v(t)));
        mtau(t) = 1.0/(alpha_m(v(t))+beta_m(v(t)));
        htau(t) = 1.0/(alpha_h(v(t))+beta_h(v(t)));       
        
        %Vetor de tempo
        tv(t)=t*Dt;
    end
    % Arrumando o potencial de repouso (lembrar aqui que as equacoes acima sao definidas em termos dos potencias de repouso)
    v = v+v0;
end


% Exercício 5: Modelo de Connor-Stevens [FUNCTION]
% 5.a) A-Current                        [FUNCTION]

%% Funcao para o calculo de alpha_m(V) 
function aln=CSalpha_m(V)
    % Parte superior e inferior da equacao para alpha_n
    s = 0.38*(V+29.7);
    i = 1.0-exp(-0.1*(V+29.7));
    %verificar a condicao para nao haver indeterminacao
    if ( (s==0) && (i==0))
        s = 0.38;
        i = 0.1*exp(-0.1*(V+29.7));
    end
    aln=s/i;
end

%% Funcao para o calculo de beta_m(V) 
function aln=CSbeta_m(V)
    % Para esse canal nao ha a necessidade de verificacoes
    aln=15.2*exp(-0.0556*(V+54.7));
end

%% Funcao para o calculo de alpha_n(V) 
function aln=CSalpha_n(V)
    % Parte superior e inferior da equacao para alpha_n
    s = 0.02*(V+45.7);
    i = 1.0-exp(-0.1*(V+45.7));
    %verificar a condicao para nao haver indeterminacao
    if ((s==0) && (i==0))
        s = 0.02;
        i = 0.1*exp(-0.1*(V+45.7));
    end
    aln=s/i;
end

%% Funcao para o calculo de beta_n(V) 
function aln=CSbeta_n(V)
    % Para esse canal nao ha a necessidade de verificacoes
    aln=0.25*exp(-0.0125*(V+55.7));
end

%% Funcao para o calculo de alpha_h(V) 
function aln=CSalpha_h(V)
    % Para esse canal nao ha a necessidade de verificacoes
    aln=0.266*exp(-0.05*(V+48));
end

%% Funcao para o calculo de beta_h(V) 
function aln=CSbeta_h(V)
    % Nao ha problema para o calculo do valor de beta_h
    s = 3.8;
    i = 1.0+exp(-0.1*(V+18.0));
    aln=s/i;
end

function [v,tv,n,m,h,a,b,ninf,minf,hinf,ainf,binf,ntau,mtau,htau,atau,btau,g_Na,g_k,J] = CS_RKUTA4(v0,ti,tf,Dt,ti_J,tf_J,Jin,J2,Jin2)
    
	%parametros eletricos da membrana
    Ena  = 55;
    Ek   = -72;
    Evaz = -17; %El
    Ea   = -75;
    GNa  = 120;
    GK   = 20;
    GL   = 0.3; 
    GA   = 47.7;
    cm   = 1;
    
    %Dimensão do vetores (conforme granularidade temporal)
    len = ceil((tf-ti)/Dt);
    
    %Vetores das Condutâncias 
    g_Na = zeros(len,1);
    g_k  = zeros(len,1);
    g_a  = zeros(len,1);
    g_l  = GL;
    
    
    %Vetor de Correte Injetada
    J=zeros(len,1);
    for i=1:len
        % Decisao quanto a entrada de corrente
       if ((i*Dt>ti_J) && (i*Dt<tf_J))
           if ((J2 == 1) && (i*Dt>50))
               J(i) = Jin2;
           else
               J(i)=Jin;
           end
       else
            J(i)=0;
       end
    end
    
    %Vetores de Probabilidades (Gates)
    n  = zeros(len,1);
    m  = zeros(len,1);
    h  = zeros(len,1);
    a  = zeros(len,1); 
    b  = zeros(len,1);
    
    %Inicialização de Vetores
    ninf = zeros(len,1);
    minf = zeros(len,1);
    hinf = zeros(len,1);
    ainf = zeros(len,1);
    binf = zeros(len,1);
    ntau = zeros(len,1);
    mtau = zeros(len,1);
    htau = zeros(len,1);
    atau = zeros(len,1);
    btau = zeros(len,1);
    
    %Vetor Voltagens e Tempo
    v    = zeros(len,1);
    tv   = zeros(len,1);
    v(1) = v0+(-v0);
    tv(1)= Dt;
    
	%Probabilidade de Abertura dos Canais (para célula em repouso)
    n(1) = CSalpha_n(v(1))/(CSalpha_n(v(1))+CSbeta_n(v(1)));
    m(1) = CSalpha_m(v(1))/(CSalpha_m(v(1))+CSbeta_m(v(1)));
    h(1) = CSalpha_h(v(1))/(CSalpha_h(v(1))+CSbeta_h(v(1)));
    
    %Taxas de Mudanças (para célula em repouso)
    ninf(1) = CSalpha_n(v(1))/(CSalpha_n(v(1))+CSbeta_n(v(1)));
    minf(1) = CSalpha_m(v(1))/(CSalpha_m(v(1))+CSbeta_m(v(1)));
    hinf(1) = CSalpha_h(v(1))/(CSalpha_h(v(1))+CSbeta_h(v(1)));
    
    %Constantes de tempo (para célula em repouso)
    ntau(1) = 1.0/(CSalpha_n(v(1))+CSbeta_n(v(1)));
    mtau(1) = 1.0/(CSalpha_m(v(1))+CSbeta_m(v(1)));
    htau(1) = 1.0/(CSalpha_h(v(1))+CSbeta_h(v(1)));
    
    %A-Current
    ainf(1) = ((0.0761*exp(0.0314*(v(1)+94.22)))/(1.0+exp(0.0346*(v(1)+1.17))))^(1.0/3);
    binf(1) = (1.0/(1.0+exp(0.0688*(v(1)+53.3))))^4;
    a(1) = ainf(1);
    b(1) = binf(1);
    atau(1) = 0.3632+1.158/(1.0+exp(0.0497*(v(1)+55.96)));
    btau(1) = 1.24+2.678/(1.0+exp(0.0624*(v(1)+50)));
    
    for t=2:len    
                
        % Cálculo: Condutâncias
        g_Na(t-1) = GNa*(m(t-1)^3)*h(t-1);
        g_k(t-1)  = GK*n(t-1)^4;
        g_a(t-1)  = GA*(a(t-1)^3)*b(t-1);  % A-current
        % Cálculo: Voltagem por Runge Kutta 
        f1 = (-g_Na(t-1)*(v(t-1)-Ena)-g_k(t-1)*(v(t-1)-Ek)-g_l*(v(t-1)-Evaz) -g_a(t-1)*(v(t-1)-Ea) +J(t-1));
        f2 = (-g_Na(t-1)*(v(t-1)+(Dt/2)*f1-Ena)-g_k(t-1)*(v(t-1)+(Dt/2)*f1-Ek)-g_l*(v(t-1)+(Dt/2)*f1-Evaz) -g_a(t-1)*(v(t-1)+(Dt/2)*f1-Ea) +J(t-1));
        f3 = (-g_Na(t-1)*(v(t-1)+(Dt/2)*f2-Ena)-g_k(t-1)*(v(t-1)+(Dt/2)*f2-Ek)-g_l*(v(t-1)+(Dt/2)*f2-Evaz) -g_a(t-1)*(v(t-1)+(Dt/2)*f2-Ea) +J(t-1));
        f4 = (-g_Na(t-1)*(v(t-1)+(Dt)*f3-Ena)-g_k(t-1)*(v(t-1)+(Dt)*f3-Ek)-g_l*(v(t-1)+(Dt)*f3-Evaz) -g_a(t-1)*(v(t-1)+(Dt)*f3-Ea) +J(t-1));       
        v(t)= v(t-1) + (Dt/(6*cm))*(f1+2*f2+2*f3+f4);
        % Cálculo: Probabilidade n por Runge-Kutta
        n1=      ((1.0-n(t-1))*CSalpha_n(v(t-1))-n(t-1)*CSbeta_n(v(t-1)));
        n2=      ((1.0-n(t-1)-(Dt/2)*n1)*CSalpha_n(v(t-1))-(n(t-1)+(Dt/2)*n1)*CSbeta_n(v(t-1)));
        n3=      ((1.0-n(t-1)-(Dt/2)*n2)*CSalpha_n(v(t-1))-(n(t-1)+(Dt/2)*n2)*CSbeta_n(v(t-1)));
        n4=      ((1.0-n(t-1)-(Dt)*n3)*CSalpha_n(v(t-1))-(n(t-1)+(Dt)*n3)*CSbeta_n(v(t-1)));
        n(t) = n(t-1) + (Dt/6)*(n1+2*n2+2*n3+n4);
        % Cálculo: Probabilidade m por Runge-Kutta
        m1=       ((1.0-m(t-1))*CSalpha_m(v(t-1))-m(t-1)*CSbeta_m(v(t-1)));
        m2=       ((1.0-m(t-1)-(Dt/2)*m1)*CSalpha_m(v(t-1))-(m(t-1)+(Dt/2)*m1)*CSbeta_m(v(t-1)));
        m3=       ((1.0-m(t-1)-(Dt/2)*m2)*CSalpha_m(v(t-1))-(m(t-1)+(Dt/2)*m2)*CSbeta_m(v(t-1)));
        m4=       ((1.0-m(t-1)-Dt*m3)*CSalpha_m(v(t-1))-(m(t-1)+Dt*m3)*CSbeta_m(v(t-1)));
        m(t) = m(t-1) + (Dt/6)*(m1+2*m2+2*m3+m4);
        % Cálculo: Probabilidade n por Runge-Kutta
        h1=       ((1.0-h(t-1))*CSalpha_h(v(t-1))-(h(t-1))*CSbeta_h(v(t-1)));
        h2=       ((1.0-h(t-1)-(Dt/2)*h1)*CSalpha_h(v(t-1))-(h(t-1)+(Dt/2)*h1)*CSbeta_h(v(t-1)));
        h3=       ((1.0-h(t-1)-(Dt/2)*h2)*CSalpha_h(v(t-1))-(h(t-1)+(Dt/2)*h2)*CSbeta_h(v(t-1)));
        h4=       ((1.0-h(t-1)-Dt*h3)*CSalpha_h(v(t-1))-(h(t-1)+Dt*h3)*CSbeta_h(v(t-1)));
        h(t) = h(t-1) + (Dt/6)*(h1+2*h2+2*h3+h4);
        
        
        %A-Current
        ainf(t) = ((0.0761*exp(0.0314*(v(t)+94.22)))/(1.0+exp(0.0346*(v(t)+1.17))))^(1.0/3);
        binf(t) = (1.0/(1.0+exp(0.0688*(v(t)+53.3))))^4;
        atau(t) = 0.3632+1.158/(1.0+exp(0.0497*(v(t)+55.96)));
        btau(t) = 1.24+2.678/(1.0+exp(0.0624*(v(t)+50)));
        a(t) = a(t-1) + (ainf(t)-a(t-1))*Dt/atau(t);
        b(t) = b(t-1) + (binf(t)-b(t-1))*Dt/btau(t);
        
        
        % Cálculo: Taxas de mudança: ninf, minf, hinf
        ninf(t) = CSalpha_n(v(t))/(CSalpha_n(v(t))+CSbeta_n(v(t)));
        minf(t) = CSalpha_m(v(t))/(CSalpha_m(v(t))+CSbeta_m(v(t)));
        hinf(t) = CSalpha_h(v(t))/(CSalpha_h(v(t))+CSbeta_h(v(t)));
        
        % Cálculo: Constantes de Tempo
        ntau(t) = 1.0/(CSalpha_n(v(t))+CSbeta_n(v(t)));
        mtau(t) = 1.0/(CSalpha_m(v(t))+CSbeta_m(v(t)));
        htau(t) = 1.0/(CSalpha_h(v(t))+CSbeta_h(v(t)));       
        
        %Vetor de tempo
        tv(t)=t*Dt;
    end
    % Arrumando o potencial de repouso (lembrar aqui que as equacoes acima sao definidas em termos dos potencias de repouso)
    v = v+v0;
end


% Exercício 5: Modelo de Connor-Stevens [FUNCTION]
% 5.b) Corrente Transiente de Cálcio    [FUNCTION]

function [v,tv,n,m,h,a,b,ninf,minf,hinf,ainf,binf,ntau,mtau,htau,atau,btau,g_Na,g_k,J] = CS_CAT_RKUTA4(v0,ti,tf,Dt,ti_J,tf_J,Jin,J2,Jin2)
    
	%parametros eletricos da membrana
    Ena  = 55;
    Ek   = -72;
    Evaz = -17; %El
    Ea   = -75;
    Eca  = 120;
    GCA  = 1.3;
    GNa  = 120;
    GK   = 20;
    GL   = 0.3; 
    GA   = 47.7;
    cm   = 1;
    
    %Dimensão do vetores (conforme granularidade temporal)
    len = ceil((tf-ti)/Dt);
    
    %Vetores das Condutâncias 
    g_Na = zeros(len,1);
    g_k  = zeros(len,1);
    g_a  = zeros(len,1);
    g_ca = zeros(len,1);
    g_l  = GL;
    
    
    %Vetor de Correte Injetada
    J=zeros(len,1);
    for i=1:len
        % Decisao quanto a entrada de corrente
       if ((i*Dt>ti_J) && (i*Dt<tf_J))
           if ((J2 == 1) && (i*Dt>50))
               J(i) = Jin2;
           else
               J(i)=Jin;
           end
       else
            J(i)=0;
       end
    end
    
    %Vetores de Probabilidades (Gates)
    n  = zeros(len,1);
    m  = zeros(len,1);
    h  = zeros(len,1);
    a  = zeros(len,1); 
    b  = zeros(len,1);
    m_ca = zeros(len,1);
    h_ca = zeros(len,1);
    
    %Inicialização de Vetores
    ninf = zeros(len,1);
    minf = zeros(len,1);
    hinf = zeros(len,1);
    ainf = zeros(len,1);
    binf = zeros(len,1);
    ntau = zeros(len,1);
    mtau = zeros(len,1);
    htau = zeros(len,1);
    atau = zeros(len,1);
    btau = zeros(len,1);
    m_cainf = zeros(len,1);
    h_cainf = zeros(len,1);
    m_catau = zeros(len,1);
    h_catau = zeros(len,1);
    
    %Vetor Voltagens e Tempo
    v    = zeros(len,1);
    tv   = zeros(len,1);
    v(1) = v0+(-v0);
    tv(1)= Dt;
    
	%Probabilidade de Abertura dos Canais (para célula em repouso)
    n(1) = CSalpha_n(v(1))/(CSalpha_n(v(1))+CSbeta_n(v(1)));
    m(1) = CSalpha_m(v(1))/(CSalpha_m(v(1))+CSbeta_m(v(1)));
    h(1) = CSalpha_h(v(1))/(CSalpha_h(v(1))+CSbeta_h(v(1)));
    
    %Taxas de Mudanças (para célula em repouso)
    ninf(1) = CSalpha_n(v(1))/(CSalpha_n(v(1))+CSbeta_n(v(1)));
    minf(1) = CSalpha_m(v(1))/(CSalpha_m(v(1))+CSbeta_m(v(1)));
    hinf(1) = CSalpha_h(v(1))/(CSalpha_h(v(1))+CSbeta_h(v(1)));
    
    %Constantes de tempo (para célula em repouso)
    ntau(1) = 1.0/(CSalpha_n(v(1))+CSbeta_n(v(1)));
    mtau(1) = 1.0/(CSalpha_m(v(1))+CSbeta_m(v(1)));
    htau(1) = 1.0/(CSalpha_h(v(1))+CSbeta_h(v(1)));
    
    %A-Current
    ainf(1) = ((0.0761*exp(0.0314*(v(1)+94.22)))/(1.0+exp(0.0346*(v(1)+1.17))))^(1.0/3);
    binf(1) = (1.0/(1.0+exp(0.0688*(v(1)+53.3))))^4;
    a(1) = ainf(1);
    b(1) = binf(1);
    atau(1) = 0.3632+1.158/(1.0+exp(0.0497*(v(1)+55.96)));
    btau(1) = 1.24+2.678/(1.0+exp(0.0624*(v(1)+50)));
    
    %CaT Current
    m_cainf(1) = 1.0/(1.0+exp(-(v(1)+57)/6.2));
    h_cainf(1) = 1.0/(1.0+exp((v(1)+81)/4));
    m_ca(1) = m_cainf(1);
    h_ca(1) = h_cainf(1);
    m_catau(1) = 0.612+(exp(-(v(1)+132)/16.7) + exp((v(1)+16.8)/18.2))^(-1);
    if (v(1)<-80)
        h_catau(1) = exp((v(1)+467)/66.6);
    else
        h_catau(1) = 28+exp(-(v(1)+22)/10.5);
    end
    
    for t=2:len    
        
        % Cálculo: Condutâncias
        g_Na(t-1) = GNa*(m(t-1)^3)*h(t-1);
        g_k(t-1)  = GK*n(t-1)^4;
        g_a(t-1)  = GA*(a(t-1)^3)*b(t-1);  % A-current
        g_ca(t-1) = GCA*(m_ca(t-1)^2)*h_ca(t-1);
        % Cálculo: Voltagem por Runge Kutta 
        f1 = (-g_Na(t-1)*(v(t-1)-Ena)-g_k(t-1)*(v(t-1)-Ek)-g_l*(v(t-1)-Evaz) -g_a(t-1)*(v(t-1)-Ea) -g_ca(t-1)*(v(t-1)-Eca) +J(t-1));
        f2 = (-g_Na(t-1)*(v(t-1)+(Dt/2)*f1-Ena)-g_k(t-1)*(v(t-1)+(Dt/2)*f1-Ek)-g_l*(v(t-1)+(Dt/2)*f1-Evaz) -g_a(t-1)*(v(t-1)+(Dt/2)*f1-Ea) -g_ca(t-1)*(v(t-1)+(Dt/2)*f1-Eca) +J(t-1));
        f3 = (-g_Na(t-1)*(v(t-1)+(Dt/2)*f2-Ena)-g_k(t-1)*(v(t-1)+(Dt/2)*f2-Ek)-g_l*(v(t-1)+(Dt/2)*f2-Evaz) -g_a(t-1)*(v(t-1)+(Dt/2)*f2-Ea) -g_ca(t-1)*(v(t-1)+(Dt/2)*f2-Eca) +J(t-1));
        f4 = (-g_Na(t-1)*(v(t-1)+(Dt)*f3-Ena)-g_k(t-1)*(v(t-1)+(Dt)*f3-Ek)-g_l*(v(t-1)+(Dt)*f3-Evaz) -g_a(t-1)*(v(t-1)+(Dt)*f3-Ea) -g_ca(t-1)*(v(t-1)+(Dt)*f3-Eca) +J(t-1));       
        v(t)= v(t-1) + (Dt/(6*cm))*(f1+2*f2+2*f3+f4);
        % Cálculo: Probabilidade n por Runge-Kutta
        n1=      ((1.0-n(t-1))*CSalpha_n(v(t-1))-n(t-1)*CSbeta_n(v(t-1)));
        n2=      ((1.0-n(t-1)-(Dt/2)*n1)*CSalpha_n(v(t-1))-(n(t-1)+(Dt/2)*n1)*CSbeta_n(v(t-1)));
        n3=      ((1.0-n(t-1)-(Dt/2)*n2)*CSalpha_n(v(t-1))-(n(t-1)+(Dt/2)*n2)*CSbeta_n(v(t-1)));
        n4=      ((1.0-n(t-1)-(Dt)*n3)*CSalpha_n(v(t-1))-(n(t-1)+(Dt)*n3)*CSbeta_n(v(t-1)));
        n(t) = n(t-1) + (Dt/6)*(n1+2*n2+2*n3+n4);
        % Cálculo: Probabilidade m por Runge-Kutta
        m1=       ((1.0-m(t-1))*CSalpha_m(v(t-1))-m(t-1)*CSbeta_m(v(t-1)));
        m2=       ((1.0-m(t-1)-(Dt/2)*m1)*CSalpha_m(v(t-1))-(m(t-1)+(Dt/2)*m1)*CSbeta_m(v(t-1)));
        m3=       ((1.0-m(t-1)-(Dt/2)*m2)*CSalpha_m(v(t-1))-(m(t-1)+(Dt/2)*m2)*CSbeta_m(v(t-1)));
        m4=       ((1.0-m(t-1)-Dt*m3)*CSalpha_m(v(t-1))-(m(t-1)+Dt*m3)*CSbeta_m(v(t-1)));
        m(t) = m(t-1) + (Dt/6)*(m1+2*m2+2*m3+m4);
        % Cálculo: Probabilidade n por Runge-Kutta
        h1=       ((1.0-h(t-1))*CSalpha_h(v(t-1))-(h(t-1))*CSbeta_h(v(t-1)));
        h2=       ((1.0-h(t-1)-(Dt/2)*h1)*CSalpha_h(v(t-1))-(h(t-1)+(Dt/2)*h1)*CSbeta_h(v(t-1)));
        h3=       ((1.0-h(t-1)-(Dt/2)*h2)*CSalpha_h(v(t-1))-(h(t-1)+(Dt/2)*h2)*CSbeta_h(v(t-1)));
        h4=       ((1.0-h(t-1)-Dt*h3)*CSalpha_h(v(t-1))-(h(t-1)+Dt*h3)*CSbeta_h(v(t-1)));
        h(t) = h(t-1) + (Dt/6)*(h1+2*h2+2*h3+h4);
        
        
        %A-Current
        ainf(t) = ((0.0761*exp(0.0314*(v(t)+94.22)))/(1.0+exp(0.0346*(v(t)+1.17))))^(1.0/3);
        binf(t) = (1.0/(1.0+exp(0.0688*(v(t)+53.3))))^4;
        atau(t) = 0.3632+1.158/(1.0+exp(0.0497*(v(t)+55.96)));
        btau(t) = 1.24+2.678/(1.0+exp(0.0624*(v(t)+50)));
        a(t) = a(t-1) + (ainf(t)-a(t-1))*Dt/atau(t);
        b(t) = b(t-1) + (binf(t)-b(t-1))*Dt/btau(t);
        
        %CaT Current
        m_cainf(t) = 1.0/(1.0+exp(-(v(t)+57)/6.2));
        h_cainf(t) = 1.0/(1.0+exp((v(t)+81)/4));
        m_catau(t) = 0.612+(exp(-(v(t)+132)/16.7) + exp((v(t)+16.8)/18.2))^(-1);
        if (v(t)<-80)
            h_catau(t) = exp((v(t)+467)/66.6);
        else
            h_catau(t) = 28+exp(-(v(t)+22)/10.5);
        end
        m_ca(t) = m_ca(t-1) + (m_cainf(t)-m_ca(t-1))*Dt/m_catau(t);
        h_ca(t) = h_ca(t-1) + (h_cainf(t)-h_ca(t-1))*Dt/h_catau(t);
        
        
        % Cálculo: Taxas de mudança: ninf, minf, hinf
        ninf(t) = CSalpha_n(v(t))/(CSalpha_n(v(t))+CSbeta_n(v(t)));
        minf(t) = CSalpha_m(v(t))/(CSalpha_m(v(t))+CSbeta_m(v(t)));
        hinf(t) = CSalpha_h(v(t))/(CSalpha_h(v(t))+CSbeta_h(v(t)));
        
        % Cálculo: Constantes de Tempo
        ntau(t) = 1.0/(CSalpha_n(v(t))+CSbeta_n(v(t)));
        mtau(t) = 1.0/(CSalpha_m(v(t))+CSbeta_m(v(t)));
        htau(t) = 1.0/(CSalpha_h(v(t))+CSbeta_h(v(t)));       
        
        %Vetor de tempo
        tv(t)=t*Dt;
    end
    % Arrumando o potencial de repouso (lembrar aqui que as equacoes acima sao definidas em termos dos potencias de repouso)
    v = v+v0;
end


% Exercício 6: Modelo de Traub-Miles [FUNCTION]
% 6.a) Corrente Adaptativa z         [FUNCTION]

%% Funcao para o calculo de alpha_m(V) 
function aln=TMalpha_m(V)
    % Parte superior e inferior da equacao para alpha_n
    s = 0.32*(54.0+V);
    i = 1-exp(-(V+54.0)/4.0);
    %verificar a condicao para nao haver indeterminacao
    if ( (s==0) && (i==0))      
        s = 0.32;
        i = (1.0/4.0)*exp(-(V+54.0)/4.0);
    end
    aln=s/i;
end

%% Funcao para o calculo de beta_m(V) 
function aln=TMbeta_m(V)
    % Parte superior e inferior da equacao para alpha_n
    s = 0.28*(V+27.0);
    i = exp((V+27.0)/5.0)-1.0;
    %verificar a condicao para nao haver indeterminacao
    if ( (s==0) && (i==0))              
        s = 0.28;
        i = (1.0/5.0)*exp((V+27.0)/5.0);
    end
    aln=s/i;
end

%% Funcao para o calculo de alpha_n(V) 
function aln=TMalpha_n(V)
    % Parte superior e inferior da equacao para alpha_n
    s = 0.032*(V+52.0);
    i = 1-exp(-(V+52.0)/5.0);
    %verificar a condicao para nao haver indeterminacao
    if ((s==0) && (i==0))        
        s = 0.032;
        i = (1.0/5.0)*exp(-(V+52.0)/5.0);
    end
    aln=s/i;
end

%% Funcao para o calculo de beta_n(V) 
function aln=TMbeta_n(V)
    % Para esse canal nao ha a necessidade de verificacoes
    aln = 0.5*exp(-(57.0+V)/40.0);
end

%% Funcao para o calculo de alpha_h(V) 
function aln=TMalpha_h(V)
    % Para esse canal nao ha a necessidade de verificacoes
    aln = 0.128*exp(-(50.0+V)/18.0);
end

%% Funcao para o calculo de beta_h(V) 
function aln=TMbeta_h(V)
    % Nao ha problema para o calculo do valor de beta_h
    s = 4.0;
    i = 1+exp(-(V+27.0)/5.0);
    aln=s/i;
end

function [v,tv,n,m,h,z,g_Na,g_k,J] = TM_RKUTA4(v0,ti,tf,Dt,ti_J,tf_J,Jin)
    
	%parametros eletricos da membrana
    Ena  = 50;
    Ek   = -100;
    Evaz = -67; 
    GNa  = 100;
    GK   = 80;
    GL   = 0.1; 
    G    = 5;  %Adaptive Current
    cm   = 1;
    
    %Dimensão do vetores (conforme granularidade temporal)
    len = ceil((tf-ti)/Dt);
    
    %Vetores das Condutâncias 
    g_Na = zeros(len,1);
    g_k  = zeros(len,1);
    g_z  = zeros(len,1);
    g_l  = GL;
    
    %Vetor de Correte Injetada
    J=zeros(len,1);
    for i=1:len
        % Decisao quanto a entrada de corrente
       if ((i*Dt>ti_J) && (i*Dt<tf_J))
            J(i)=Jin;
       else
            J(i)=0;
       end
    end
    
    %Vetores de Probabilidades (Gates)
    n  = zeros(len,1);
    m  = zeros(len,1);
    h  = zeros(len,1);
    z  = zeros(len,1);
    
    %Vetor Voltagens e Tempo
    v    = zeros(len,1);
    tv   = zeros(len,1);
    v(1) = v0+(-v0);
    tv(1)= Dt;
    
	%Probabilidade de Abertura dos Canais (para célula em repouso)
    n(1) = TMalpha_n(v(1))/(TMalpha_n(v(1))+TMbeta_n(v(1)));
    m(1) = TMalpha_m(v(1))/(TMalpha_m(v(1))+TMbeta_m(v(1)));
    h(1) = TMalpha_h(v(1))/(TMalpha_h(v(1))+TMbeta_h(v(1)));
    z(1) = (0.01/(1+exp(-(v(1)+20)/5)))-0.01*(0);  % or 0?
    
    for t=2:len
        % Cálculo: Condutâncias
        g_Na(t-1) = GNa*(m(t-1)^3)*h(t-1);
        g_k(t-1)  = GK*n(t-1)^4;
        g_z(t-1)  = G*z(t-1);
        % Cálculo: Voltagem por Runge Kutta 
        f1 = (-g_Na(t-1)*(v(t-1)-Ena)-g_k(t-1)*(v(t-1)-Ek)-g_l*(v(t-1)-Evaz) -g_z(t-1)*(v(t-1)-Ek) +J(t-1));
        f2 = (-g_Na(t-1)*(v(t-1)+(Dt/2)*f1-Ena)-g_k(t-1)*(v(t-1)+(Dt/2)*f1-Ek)-g_l*(v(t-1)+(Dt/2)*f1-Evaz) -g_z(t-1)*(v(t-1)+(Dt/2)*f1-Ek) +J(t-1));
        f3 = (-g_Na(t-1)*(v(t-1)+(Dt/2)*f2-Ena)-g_k(t-1)*(v(t-1)+(Dt/2)*f2-Ek)-g_l*(v(t-1)+(Dt/2)*f2-Evaz) -g_z(t-1)*(v(t-1)+(Dt/2)*f2-Ek) +J(t-1));
        f4 = (-g_Na(t-1)*(v(t-1)+(Dt)*f3-Ena)-g_k(t-1)*(v(t-1)+(Dt)*f3-Ek)-g_l*(v(t-1)+(Dt)*f3-Evaz) -g_z(t-1)*(v(t-1)+(Dt)*f3-Ek) +J(t-1));
        v(t)= v(t-1) + (Dt/(6*cm))*(f1+2*f2+2*f3+f4);
        % Cálculo: Probabilidade n por Runge-Kutta
        n1=      ((1.0-n(t-1))*TMalpha_n(v(t-1))-n(t-1)*TMbeta_n(v(t-1)));
        n2=      ((1.0-n(t-1)-(Dt/2)*n1)*TMalpha_n(v(t-1))-(n(t-1)+(Dt/2)*n1)*TMbeta_n(v(t-1)));
        n3=      ((1.0-n(t-1)-(Dt/2)*n2)*TMalpha_n(v(t-1))-(n(t-1)+(Dt/2)*n2)*TMbeta_n(v(t-1)));
        n4=      ((1.0-n(t-1)-(Dt)*n3)*TMalpha_n(v(t-1))-(n(t-1)+(Dt)*n3)*TMbeta_n(v(t-1)));
        n(t) = n(t-1) + (Dt/6)*(n1+2*n2+2*n3+n4);
        % Cálculo: Probabilidade m por Runge-Kutta
        m1=       ((1.0-m(t-1))*TMalpha_m(v(t-1))-m(t-1)*TMbeta_m(v(t-1)));
        m2=       ((1.0-m(t-1)-(Dt/2)*m1)*TMalpha_m(v(t-1))-(m(t-1)+(Dt/2)*m1)*TMbeta_m(v(t-1)));
        m3=       ((1.0-m(t-1)-(Dt/2)*m2)*TMalpha_m(v(t-1))-(m(t-1)+(Dt/2)*m2)*TMbeta_m(v(t-1)));
        m4=       ((1.0-m(t-1)-Dt*m3)*TMalpha_m(v(t-1))-(m(t-1)+Dt*m3)*TMbeta_m(v(t-1)));
        m(t) = m(t-1) + (Dt/6)*(m1+2*m2+2*m3+m4);
        % Cálculo: Probabilidade n por Runge-Kutta
        h1=       ((1.0-h(t-1))*TMalpha_h(v(t-1))-(h(t-1))*TMbeta_h(v(t-1)));
        h2=       ((1.0-h(t-1)-(Dt/2)*h1)*TMalpha_h(v(t-1))-(h(t-1)+(Dt/2)*h1)*TMbeta_h(v(t-1)));
        h3=       ((1.0-h(t-1)-(Dt/2)*h2)*TMalpha_h(v(t-1))-(h(t-1)+(Dt/2)*h2)*TMbeta_h(v(t-1)));
        h4=       ((1.0-h(t-1)-Dt*h3)*TMalpha_h(v(t-1))-(h(t-1)+Dt*h3)*TMbeta_h(v(t-1)));
        h(t) = h(t-1) + (Dt/6)*(h1+2*h2+2*h3+h4);
        % Cálculo: Probabilidade z por Runge-Kutta
        z1=       (0.01/(1+exp(-(v(t-1)+20)/5)))-0.01*(z(t-1));
        z2=       (0.01/(1+exp(-(v(t-1)+20)/5)))-0.01*(z(t-1)+(Dt/2)*z1);
        z3=       (0.01/(1+exp(-(v(t-1)+20)/5)))-0.01*(z(t-1)+(Dt/2)*z2);
        z4=       (0.01/(1+exp(-(v(t-1)+20)/5)))-0.01*(z(t-1)+Dt*z3);
        z(t) = z(t-1) + (Dt/6)*(z1+2*z2+2*z3+z4);
        %Vetor de tempo
        tv(t)=t*Dt;
    end
    % Arrumando o potencial de repouso (lembrar aqui que as equacoes acima sao definidas em termos dos potencias de repouso)
    v = v+v0;
end


% Exercício 6: Modelo de Traub-Miles        [FUNCTION]
% 6.b) Corrente Adaptativa: Bi-estabilidade [FUNCTION]

function [v,tv,n,m,h,z,g_Na,g_k,J] = TM_RKUTA4_v2(v0,ti,tf,Dt,ti_J1,tf_J1,Jin1,ti_J2,tf_J2,Jin2)
    
	%parametros eletricos da membrana
    Ena  = 50;
    Ek   = -100;
    Evaz = -67; 
    GNa  = 100;
    GK   = 80;
    GL   = 0.1; 
    G    = 1;  %Adaptive Current
    cm   = 1;
    
    %Dimensão do vetores (conforme granularidade temporal)
    len = ceil((tf-ti)/Dt);
    
    %Vetores das Condutâncias 
    g_Na = zeros(len,1);
    g_k  = zeros(len,1);
    g_z  = zeros(len,1);
    g_l  = GL;
    
    %Vetor de Correte Injetada
    J=zeros(len,1);
    for i=1:len
        % Decisao quanto a entrada de corrente
       if ((i*Dt>ti_J1) && (i*Dt<tf_J1))
            J(i)=Jin1;
       elseif ((i*Dt>ti_J2) && (i*Dt<tf_J2))    
            J(i) = Jin2;
       else
            J(i)=0;
       end
    end
    
    %Vetores de Probabilidades (Gates)
    n  = zeros(len,1);
    m  = zeros(len,1);
    h  = zeros(len,1);
    z  = zeros(len,1);
    
    %Vetor Voltagens e Tempo
    v    = zeros(len,1);
    tv   = zeros(len,1);
    v(1) = v0+(-v0);
    tv(1)= Dt;
    
	%Probabilidade de Abertura dos Canais (para célula em repouso)
    n(1) = TMalpha_n(v(1))/(TMalpha_n(v(1))+TMbeta_n(v(1)));
    m(1) = TMalpha_m(v(1))/(TMalpha_m(v(1))+TMbeta_m(v(1)));
    h(1) = TMalpha_h(v(1))/(TMalpha_h(v(1))+TMbeta_h(v(1)));
    z(1) = 0;%(0.01/(1+exp(-(v(1)+20)/5)))-0.01*(0);  % or 0?
    
    for t=2:len
        % Cálculo: Condutâncias
        g_Na(t-1) = GNa*(m(t-1)^3)*h(t-1);
        g_k(t-1)  = GK*n(t-1)^4;
        g_z(t-1)  = G*z(t-1);
        % Cálculo: Voltagem por Runge Kutta 
        f1 = (-g_Na(t-1)*(v(t-1)-Ena)-g_k(t-1)*(v(t-1)-Ek)-g_l*(v(t-1)-Evaz) -g_z(t-1)*(v(t-1)-0) +J(t-1));
        f2 = (-g_Na(t-1)*(v(t-1)+(Dt/2)*f1-Ena)-g_k(t-1)*(v(t-1)+(Dt/2)*f1-Ek)-g_l*(v(t-1)+(Dt/2)*f1-Evaz) -g_z(t-1)*(v(t-1)+(Dt/2)*f1-0) +J(t-1));
        f3 = (-g_Na(t-1)*(v(t-1)+(Dt/2)*f2-Ena)-g_k(t-1)*(v(t-1)+(Dt/2)*f2-Ek)-g_l*(v(t-1)+(Dt/2)*f2-Evaz) -g_z(t-1)*(v(t-1)+(Dt/2)*f2-0) +J(t-1));
        f4 = (-g_Na(t-1)*(v(t-1)+(Dt)*f3-Ena)-g_k(t-1)*(v(t-1)+(Dt)*f3-Ek)-g_l*(v(t-1)+(Dt)*f3-Evaz) -g_z(t-1)*(v(t-1)+(Dt)*f3-0) +J(t-1));
        v(t)= v(t-1) + (Dt/(6*cm))*(f1+2*f2+2*f3+f4);
        % Cálculo: Probabilidade n por Runge-Kutta
        n1=      ((1.0-n(t-1))*TMalpha_n(v(t-1))-n(t-1)*TMbeta_n(v(t-1)));
        n2=      ((1.0-n(t-1)-(Dt/2)*n1)*TMalpha_n(v(t-1))-(n(t-1)+(Dt/2)*n1)*TMbeta_n(v(t-1)));
        n3=      ((1.0-n(t-1)-(Dt/2)*n2)*TMalpha_n(v(t-1))-(n(t-1)+(Dt/2)*n2)*TMbeta_n(v(t-1)));
        n4=      ((1.0-n(t-1)-(Dt)*n3)*TMalpha_n(v(t-1))-(n(t-1)+(Dt)*n3)*TMbeta_n(v(t-1)));
        n(t) = n(t-1) + (Dt/6)*(n1+2*n2+2*n3+n4);
        % Cálculo: Probabilidade m por Runge-Kutta
        m1=       ((1.0-m(t-1))*TMalpha_m(v(t-1))-m(t-1)*TMbeta_m(v(t-1)));
        m2=       ((1.0-m(t-1)-(Dt/2)*m1)*TMalpha_m(v(t-1))-(m(t-1)+(Dt/2)*m1)*TMbeta_m(v(t-1)));
        m3=       ((1.0-m(t-1)-(Dt/2)*m2)*TMalpha_m(v(t-1))-(m(t-1)+(Dt/2)*m2)*TMbeta_m(v(t-1)));
        m4=       ((1.0-m(t-1)-Dt*m3)*TMalpha_m(v(t-1))-(m(t-1)+Dt*m3)*TMbeta_m(v(t-1)));
        m(t) = m(t-1) + (Dt/6)*(m1+2*m2+2*m3+m4);
        % Cálculo: Probabilidade n por Runge-Kutta
        h1=       ((1.0-h(t-1))*TMalpha_h(v(t-1))-(h(t-1))*TMbeta_h(v(t-1)));
        h2=       ((1.0-h(t-1)-(Dt/2)*h1)*TMalpha_h(v(t-1))-(h(t-1)+(Dt/2)*h1)*TMbeta_h(v(t-1)));
        h3=       ((1.0-h(t-1)-(Dt/2)*h2)*TMalpha_h(v(t-1))-(h(t-1)+(Dt/2)*h2)*TMbeta_h(v(t-1)));
        h4=       ((1.0-h(t-1)-Dt*h3)*TMalpha_h(v(t-1))-(h(t-1)+Dt*h3)*TMbeta_h(v(t-1)));
        h(t) = h(t-1) + (Dt/6)*(h1+2*h2+2*h3+h4);
        % Cálculo: Probabilidade z por Runge-Kutta
        z1=       (0.01/(1+exp(-(v(t-1)+20)/5)))-0.01*(z(t-1));
        z2=       (0.01/(1+exp(-(v(t-1)+20)/5)))-0.01*(z(t-1)+(Dt/2)*z1);
        z3=       (0.01/(1+exp(-(v(t-1)+20)/5)))-0.01*(z(t-1)+(Dt/2)*z2);
        z4=       (0.01/(1+exp(-(v(t-1)+20)/5)))-0.01*(z(t-1)+Dt*z3);
        z(t) = z(t-1) + (Dt/6)*(z1+2*z2+2*z3+z4);
        %Vetor de tempo
        tv(t)=t*Dt;
    end
    % Arrumando o potencial de repouso (lembrar aqui que as equacoes acima sao definidas em termos dos potencias de repouso)
    v = v+v0;
end


% Exercício 6: Modelo de Traub-Miles             [FUNCTION]
% 6.c) Corrente de Potássio Dependente de Cálcio [FUNCTION]

function [v,tv,n,m,h,ca,mca,minf,g_Na,g_k,J] = TM_RKUTA4_v3(v0,ti,tf,Dt,ti_J,tf_J,Jin)
    
	%parametros eletricos da membrana
    Ena  = 50;
    Ek   = -100;
    Evaz = -67; 
    GNa  = 100;
    GK   = 80;
    GL   = 0.1; 
    GAHP = 5;  %Adaptive Current
    GCA  = 1;
    Eca  = 120;
    cm   = 1;
    
    %Dimensão do vetores (conforme granularidade temporal)
    len = ceil((tf-ti)/Dt);
    
    %Vetores das Condutâncias 
    g_Na  = zeros(len,1);
    g_k   = zeros(len,1);
    g_ahp = zeros(len,1);
    g_ca  = zeros(len,1);
    g_l   = GL;
    
    %Vetor de Correte Injetada
    J=zeros(len,1);
    for i=1:len
        % Decisao quanto a entrada de corrente
       if ((i*Dt>ti_J) && (i*Dt<tf_J))
            J(i)=Jin;
       else
            J(i)=0;
       end
    end
    
    %Vetores de Probabilidades (Gates)
    n  = zeros(len,1);
    m  = zeros(len,1);
    h  = zeros(len,1);
    ca = zeros(len,1);
    mca = zeros(len,1);
    minf = zeros(len,1);
    
    %Vetor Voltagens e Tempo
    v    = zeros(len,1);
    tv   = zeros(len,1);
    v(1) = v0+(-v0);
    tv(1)= Dt;
    
	%Probabilidade de Abertura dos Canais (para célula em repouso)
    n(1) = TMalpha_n(v(1))/(TMalpha_n(v(1))+TMbeta_n(v(1)));
    m(1) = TMalpha_m(v(1))/(TMalpha_m(v(1))+TMbeta_m(v(1)));
    h(1) = TMalpha_h(v(1))/(TMalpha_h(v(1))+TMbeta_h(v(1)));
    ca(1)= 0;
    mca(1) = 0;
    minf(1) = 1/(1+exp(-(v(1)+25)/5));
    
    for t=2:len
        % Cálculo: Condutâncias
        g_Na(t-1) = GNa*(m(t-1)^3)*h(t-1);
        g_k(t-1)  = GK*n(t-1)^4;
        g_ca(t-1) = GCA*minf(t-1);
        g_ahp(t-1)= GAHP*(ca(t-1)/30);
        mca(t-1) = (ca(t-1)/30);
        % Cálculo: Voltagem por Runge Kutta 
        f1 = (-g_Na(t-1)*(v(t-1)-Ena)-g_k(t-1)*(v(t-1)-Ek)-g_l*(v(t-1)-Evaz) -g_ca(t-1)*(v(t-1)-Eca) -g_ahp(t-1)*(v(t-1)-Ek) +J(t-1));
        f2 = (-g_Na(t-1)*(v(t-1)+(Dt/2)*f1-Ena)-g_k(t-1)*(v(t-1)+(Dt/2)*f1-Ek)-g_l*(v(t-1)+(Dt/2)*f1-Evaz) -g_ca(t-1)*(v(t-1)+(Dt/2)*f1-Eca) -g_ahp(t-1)*(v(t-1)+(Dt/2)*f1-Ek) +J(t-1));
        f3 = (-g_Na(t-1)*(v(t-1)+(Dt/2)*f2-Ena)-g_k(t-1)*(v(t-1)+(Dt/2)*f2-Ek)-g_l*(v(t-1)+(Dt/2)*f2-Evaz) -g_ca(t-1)*(v(t-1)+(Dt/2)*f2-Eca) -g_ahp(t-1)*(v(t-1)+(Dt/2)*f2-Ek) +J(t-1));
        f4 = (-g_Na(t-1)*(v(t-1)+(Dt)*f3-Ena)-g_k(t-1)*(v(t-1)+(Dt)*f3-Ek)-g_l*(v(t-1)+(Dt)*f3-Evaz) -g_ca(t-1)*(v(t-1)+(Dt)*f3-Eca) -g_ahp(t-1)*(v(t-1)+(Dt)*f3-Ek) +J(t-1));
        v(t)= v(t-1) + (Dt/(6*cm))*(f1+2*f2+2*f3+f4);
        % Cálculo: Probabilidade n por Runge-Kutta
        n1=      ((1.0-n(t-1))*TMalpha_n(v(t-1))-n(t-1)*TMbeta_n(v(t-1)));
        n2=      ((1.0-n(t-1)-(Dt/2)*n1)*TMalpha_n(v(t-1))-(n(t-1)+(Dt/2)*n1)*TMbeta_n(v(t-1)));
        n3=      ((1.0-n(t-1)-(Dt/2)*n2)*TMalpha_n(v(t-1))-(n(t-1)+(Dt/2)*n2)*TMbeta_n(v(t-1)));
        n4=      ((1.0-n(t-1)-(Dt)*n3)*TMalpha_n(v(t-1))-(n(t-1)+(Dt)*n3)*TMbeta_n(v(t-1)));
        n(t) = n(t-1) + (Dt/6)*(n1+2*n2+2*n3+n4);
        % Cálculo: Probabilidade m por Runge-Kutta
        m1=       ((1.0-m(t-1))*TMalpha_m(v(t-1))-m(t-1)*TMbeta_m(v(t-1)));
        m2=       ((1.0-m(t-1)-(Dt/2)*m1)*TMalpha_m(v(t-1))-(m(t-1)+(Dt/2)*m1)*TMbeta_m(v(t-1)));
        m3=       ((1.0-m(t-1)-(Dt/2)*m2)*TMalpha_m(v(t-1))-(m(t-1)+(Dt/2)*m2)*TMbeta_m(v(t-1)));
        m4=       ((1.0-m(t-1)-Dt*m3)*TMalpha_m(v(t-1))-(m(t-1)+Dt*m3)*TMbeta_m(v(t-1)));
        m(t) = m(t-1) + (Dt/6)*(m1+2*m2+2*m3+m4);
        % Cálculo: Probabilidade n por Runge-Kutta
        h1=       ((1.0-h(t-1))*TMalpha_h(v(t-1))-(h(t-1))*TMbeta_h(v(t-1)));
        h2=       ((1.0-h(t-1)-(Dt/2)*h1)*TMalpha_h(v(t-1))-(h(t-1)+(Dt/2)*h1)*TMbeta_h(v(t-1)));
        h3=       ((1.0-h(t-1)-(Dt/2)*h2)*TMalpha_h(v(t-1))-(h(t-1)+(Dt/2)*h2)*TMbeta_h(v(t-1)));
        h4=       ((1.0-h(t-1)-Dt*h3)*TMalpha_h(v(t-1))-(h(t-1)+Dt*h3)*TMbeta_h(v(t-1)));
        h(t) = h(t-1) + (Dt/6)*(h1+2*h2+2*h3+h4);
        
        % Cálculo: Probabilidade ca por Runge-Kutta
        ica = g_ca(t-1)*minf(t-1)*(v(t-1)-Eca);
        ca1=      (-0.002*ica)-0.0125*(ca(t-1));
        ca2=      (-0.002*ica)-0.0125*(ca(t-1)+(Dt/2)*ca1);
        ca3=      (-0.002*ica)-0.0125*(ca(t-1)+(Dt/2)*ca2);
        ca4=      (-0.002*ica)-0.0125*(ca(t-1)+Dt*ca3);
        ca(t) = ca(t-1) + (Dt/6)*(ca1+2*ca2+2*ca3+ca4);
                
        minf(t) = 1/(1+exp(-(v(t-1)+25)/5));
        
        %Vetor de tempo
        tv(t)=t*Dt;
    end
    % Arrumando o potencial de repouso (lembrar aqui que as equacoes acima sao definidas em termos dos potencias de repouso)
    v = v+v0;
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

function [v,tv,n,m,h,g_Na,g_k,J] = HH_RKUTA4(v0,ti,tf,Dt,ti_J,tf_J,Jin)
    
	%parametros eletricos da membrana
    Ena  = 115;
    Ek   = -12;
    Evaz = 10.16; 
    GNa  = 120;
    GK   = 36;
    GL   = 0.3; 
    cm   = 1;
    
    %Dimensão do vetores (conforme granularidade temporal)
    len = ceil((tf-ti)/Dt);
    
    %Vetores das Condutâncias 
    g_Na = zeros(len,1);
    g_k  = zeros(len,1);
    g_l  = GL;
    
    %Vetor de Correte Injetada
    J=zeros(len,1);
    for i=1:len
        % Decisao quanto a entrada de corrente
       if ((i*Dt>ti_J) && (i*Dt<tf_J))
            J(i)=Jin;
       else
            J(i)=0;
       end
    end
    
    %Vetores de Probabilidades (Gates)
    n  = zeros(len,1);
    m  = zeros(len,1);
    h  = zeros(len,1);
   
    %Vetor Voltagens e Tempo
    v    = zeros(len,1);
    tv   = zeros(len,1);
    v(1) = v0+(-v0);
    tv(1)= Dt;
    
	%Probabilidade de Abertura dos Canais (para célular em repouso)
    n(1) = alpha_n(v(1))/(alpha_n(v(1))+beta_n(v(1)));
    m(1) = alpha_m(v(1))/(alpha_m(v(1))+beta_m(v(1)));
    h(1) = alpha_h(v(1))/(alpha_h(v(1))+beta_h(v(1)));
    
    for t=2:len
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
