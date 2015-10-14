% Segunda lista de Exercícios - Introdução à Neurociência Computacional
% Aluno: Rafael Delalibera Rodrigues
% Programa de Pós-Graduação em Computação Aplicada [FFCLRP-USP Ribeirão Preto] - Mestrado

% Exercício 3: Equação da Membrana (Solução Numérica – Euler)

euler(@funct,0,1,0,80)


% Exercício 4: Equação da Membrana (Solução Numérica – Euler) Diferentes condutâncias

[out1,t] = euler_v2(@funct_v2,0,1,0,80,200);
[out2,t] = euler_v2(@funct_v2,0,1,0,80,100);
[out3,t] = euler_v2(@funct_v2,0,1,0,80,50);
[out4,t] = euler_v2(@funct_v2,0,1,0,80,20);
[out5,t] = euler_v2(@funct_v2,0,1,0,80,10);
figure;
plot(t,out1,t,out2,t,out3,t,out4,t,out5);
legend('g=200','g=100','g=50','g=20','g=10');


% Exercício 5: Estudo Euler - Precisão

% Solução Analítica
tic;
ta=0:0.01:20; %intervalo de tempo
za = -(-1)*exp(-ta/2);
exect0 = toc;
 
%Euler
[out1, t1, exect1] = euler_v3(@funct_v3,1,0.02,0,20);
[out2, t2, exect2] = euler_v3(@funct_v3,1,0.2,0,20);
[out3, t3, exect3] = euler_v3(@funct_v3,1,1,0,20);
[out4, t4, exect4] = euler_v3(@funct_v3,1,2,0,20);
[out5, t5, exect5] = euler_v3(@funct_v3,1,6,0,20);
figure;
plot(ta,za,t1,out1,t2,out2,t3,out3,t4,out4,t5,out5);
xlabel('t')
ylabel('V')
legend('exata','step=0.02','step=0.2','step=1','step=2','step=6');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS (Must be defined at the end of the script) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exercício 3: Equação da Membrana (Solução Numérica – Euler) [FUNCTIONS]

function euler(func,y0,dt,t0,tf)

    t=t0:dt:tf; %intervalo de tempo
 
    y(1)=y0;
 
    %euler
    for i=1:length(t)-1
        y(i+1) = y(i)+dt*(feval(func,t(i),y(i)));
    end
 
    plot(t,y)
    xlabel('t')
    ylabel('V')
 
end

function z = funct(t,y)
 
    j=0;
    if t>=5 && t<=50
        j = 2;
    else
        j = 0;
    end
 
    z = (-y+j*10)/10;
 
end


% Exercício 4: Equação da Membrana (Solução Numérica – Euler) Diferentes condutâncias [FUNCTIONS]

function [out,t] = euler_v2(func,y0,dt,t0,tf,g)
    t=t0:dt:tf; %intervalo de tempo
    y(1)=y0;
    %euler
    for i=1:length(t)-1
        y(i+1) = y(i)+dt*(feval(func,t(i),y(i),g));
    end
    out = y;
end

function z = funct_v2(t,y,g)
    j=0;
    if t>=5 && t<=50
        j = 2;
    else
        j = 0;
    end
    z = (-y+j*g^(-1)*1E3)/(g^(-1)*1E3);
end


% Exercício 5: Estudo Euler - Precisão [FUNCTIONS]

function [out,t, exec_t] = euler_v3(func,y0,dt,t0,tf)
    tic;
    t=t0:dt:tf; %intervalo de tempo
    y(1)=y0;
    %euler
    for i=1:length(t)-1
        y(i+1) = y(i)+dt*(feval(func,t(i),y(i)));
    end
    out = y;
    exec_t = toc;
end

function z = funct_v3(t,y)
    z = -y/2;  
end