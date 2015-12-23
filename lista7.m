% S�tima lista de Exerc�cios - Introdu��o � Neuroci�ncia Computacional
% Aluno: Rafael Delalibera Rodrigues
% Programa de P�s-Gradua��o em Computa��o Aplicada [FFCLRP-USP Ribeir�o Preto] - Mestrado

% Exerc�cio 1: Modelo LIF

Vl = -50;
Iinj = [0.1:0.1:2.0];
freq_v = zeros(length(Iinj),1);
disp = 0;

for i=1:length(Iinj)

	[v,t,I] = LIF_RKUTA4(-60,0,350,0.01,50,250,Iinj(i));
	
	temp = v>Vl;
	disp = sum(temp);
	freq_v(i) = disp/(250-50)*1000; % segundos
	disp = 0;
end

figure(1);
plot(Iinj, freq_v);
xlabel('I (nA)');
ylabel('f (disparos/s)');
title('Curva FI do modelo LIF');


% Exerc�cio 2: Modelo LIF + Termo de Ru�do
% 2.a) Disparos do Neur�nio com Ru�do

figure(2);
s1 = subplot(2,1,1);
Vl = -50;
Dt = 0.01;
desv_pad = 6;
[v,t,I] = LIFR_RKUTA4(-60,0,5000,Dt,0,5000,0.25,desv_pad);
disp = v>Vl;
t_new = t(disp==1);
disp = disp(disp>=1);
stem(t_new,disp,'Marker','none');
xlim([2000 3000]);
ylim([0 1.5]);
set(gca,'ytick',[]);
xlabel('Tempo (ms)');
ylabel('Spikes');
title('Disparos do Modelo LIF COM Ru�do');
freq1 = length(t_new(t_new>= 2000 & t_new<=3000))/(3000-2000)*1000; % em segundos
dim = [.75 .8 .2 .1];
str = sprintf('Frequ�ncia: %i Disparos/s',freq1);
a = annotation('textbox',dim,'String',str,'FitBoxToText','on');
a.Color = 'red';

s2 = subplot(2,1,2); 
[v,t,I] = LIF_RKUTA4(-60,0,5000,Dt,0,5000,0.25);
disp = v>Vl;
t_new = t(disp==1);
disp = disp(disp>=1);
stem(t_new,disp,'Marker','none');
xlim([2000 3000]);
ylim([0 1.5]);
set(gca,'ytick',[]);
xlabel('Tempo (ms)');
ylabel('Spikes');
title('Disparos do Modelo LIF {\color{red}SEM} Ru�do');
freq2 = length(t_new(t_new>= 2000 & t_new<=3000))/(3000-2000)*1000; % em segundos
dim = [.75 .3 .2 .1];
str = sprintf('Frequ�ncia: %i Disparos/s',freq2);
b = annotation('textbox',dim,'String',str,'FitBoxToText','on');
b.Color = 'red';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS (Must be defined at the end of the script) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exerc�cio 1: Modelo LIF [FUNCTION]

function [v,tv,I] = LIF_RKUTA4(v0,ti,tf,Dt,ti_I,tf_I,Iin)
    
	% Par�metros do Modelo
    tau     = 20;   %ms
    tau_ref = 5;    %ms
    Vr      = -60;  %mV
    R       = 100;  %M*Ohm
    Vl      = -50;  %mV
    
    % Flag do per�odo refrat�rio
    is_ref  = false;
    % Per�odo refrat�rio ajustado
    tau_ref_Dt = tau_ref/Dt;
    
    %Dimens�o do vetores (conforme granularidade temporal)
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
    v(1) = v0;
    tv(1)= Dt;
    
    for t=2:len

        % Calcula caso n�o esteja no per�odo refrat�rio
        if (v(t-1)<=Vl) && (is_ref == false)
            % C�lculo: Voltagem por Runge Kutta 
            f1 = -v(t-1)+Vr+R*I(t-1);
            f2 = -(v(t-1)+(Dt/2)*f1)+Vr+R*I(t-1);
            f3 = -(v(t-1)+(Dt/2)*f2)+Vr+R*I(t-1);
            f4 = -(v(t-1)+Dt*f3)+Vr+R*I(t-1);
            v(t)= v(t-1) + (Dt/(6*tau))*(f1+2*f2+2*f3+f4);
        end
        
        % Disparo!
        % Prepara in�cio do per�odo refrat�rio
        if (v(t-1)>Vl)
            is_ref = true;
            tau_ref_Dt = tau_ref/Dt;
        end
        
        % Per�odo Refrat�rio
        if (is_ref==true)
            v(t) = Vr;
            
            if (tau_ref_Dt>0)
                % Contados do TAU ref (decremento)
                tau_ref_Dt = tau_ref_Dt - 1;
            else
                % Sinaliza o fim do per�odo refrat�rio
                is_ref = false;
            end
        end
        
        %Vetor de tempo
        tv(t)=t*Dt;
    end
end


% Exerc�cio 2: Modelo LIF + Termo de Ru�do [FUNCTION]
% 2.a) Disparos do Neur�nio com Ru�do      [FUNCTION]

function [v,tv,I] = LIFR_RKUTA4(v0,ti,tf,Dt,ti_I,tf_I,Iin,desv_pad)
    
	% Par�metros do Modelo
    tau     = 20;   %ms
    tau_ref = 5;    %ms
    Vr      = -60;  %mV
    R       = 100;  %M*Ohm
    Vl      = -50;  %mV
    
    % Flag do per�odo refrat�rio
    is_ref  = false;
    % Per�odo refrat�rio ajustado
    tau_ref_Dt = tau_ref/Dt;
    
    %Dimens�o do vetores (conforme granularidade temporal)
    len = ceil((tf-ti)/Dt);
    
    %Vetor de Correte Injetada
    I=zeros(len,1);
    for i=1:len
        % Decisao quanto a entrada de corrente
       if ((i*Dt>ti_I) && (i*Dt<tf_I))
            I(i)=Iin + (rand*desv_pad);
       else
            I(i)=0;
       end
    end
   
    %Vetor Voltagens e Tempo
    v    = zeros(len,1);
    tv   = zeros(len,1);
    v(1) = v0;
    tv(1)= Dt;
    
    for t=2:len

        % Calcula caso n�o esteja no per�odo refrat�rio
        if (v(t-1)<=Vl) && (is_ref == false)
            % C�lculo: Voltagem por Runge Kutta 
            f1 = -v(t-1)+Vr+R*I(t-1);
            f2 = -(v(t-1)+(Dt/2)*f1)+Vr+R*I(t-1);
            f3 = -(v(t-1)+(Dt/2)*f2)+Vr+R*I(t-1);
            f4 = -(v(t-1)+Dt*f3)+Vr+R*I(t-1);
            v(t)= v(t-1) + (Dt/(6*tau))*(f1+2*f2+2*f3+f4);
        end
        
        % Disparo!
        % Prepara in�cio do per�odo refrat�rio
        if (v(t-1)>Vl)
            is_ref = true;
            tau_ref_Dt = tau_ref/Dt;
        end
        
        % Per�odo Refrat�rio
        if (is_ref==true)
            v(t) = Vr;
            
            if (tau_ref_Dt>0)
                % Contados do TAU ref (decremento)
                tau_ref_Dt = tau_ref_Dt - 1;
            else
                % Sinaliza o fim do per�odo refrat�rio
                is_ref = false;
            end
        end
        
        %Vetor de tempo
        tv(t)=t*Dt;
    end
end