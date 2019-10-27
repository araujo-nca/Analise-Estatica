%%% FLUXO DE CARGA NAO LINEARIZADO - NEWTON-RAPHSON %%%

%%

clear
close all
format longEng
clc

%%  leitura do arquivo de exemplo

% Sistema_14_barra_2;
% Sistema_24_barras;
% Sistema_33_barras;
% Sistema_4_barras_Monticelli;
% SISTEMA 107 BARRAS;

%%  teste para implementacao

% ****   Sistema teste Monticelli  ****

% ------------------------------ DADOS DE BARRA -----------------------------------

%  No TB G (V)     (Ang)     (Pg)    (Qg)     (Qn)    (Qm)    (Pl)   (Ql)    bshbar
DBAR = [
    1 3 1 1.0         .0       .00   .00     -999.9  9999.9   0.00   0.00    0.0
    2 1 1 1.0         .0       .00   .07     -999.9  9999.9   0.30   0.00    0.0
    ];

% TB = 1: carga ; 3: referencia
% Tipos de barra: 1 - carga (PQ), 2 - geracao (PV), 3 - referencia (V-theta)

% ------------------------------ DADOS DE LINHA -----------------------------------

DLIN = [
    %FROM  TO   %R(pu)  %X(pu)   %Bsh     %TAP     %PHI                                   CH
    1      2     0.20    1.00    0.02     1.00     .000     .000    .0     900     .0     7
    ];

PB = 1;

%%  declaracao de variaveis

dados_barra = DBAR; % matriz de entrada com informacoes das barras do sistema
dados_linha = DLIN; % matriz de entrada com informacoes das linhas do sistema
S_base = PB;    % potencia base do sistema

n_barras = size(dados_barra, 1);    % numero de barras do sistema

erro_admitido = 0.003;  % erro admitido para fim das iteracoes

% apenas para entendimento da matriz
% bsh_bus = dados_barra(:, 12);
% r_km = dados_linha(:, 3);
% x_km = dados_linha(:, 4);
% bsh_km = dados_linha(:, 5);
% tap_km = dados_linha(:, 6);
% def_km = dados_linha(:, 7);

%%  matriz admitancia nodal

% pre-alocacao de matriz
Ybus = zeros(n_barras);

y_km = 1 ./ (dados_linha(:, 3) + j*dados_linha(:, 4));

% criacao matriz admitancia nodal
for k = 1:1:size(dados_linha, 1)
    Ybus(dados_linha(k, 1), dados_linha(k, 1)) = Ybus(dados_linha(k, 1), dados_linha(k, 1)) + dados_linha(k, 6)^2 * y_km(k) + j*dados_barra(dados_linha(k, 1), 12) + j*dados_linha(k, 5);
    Ybus(dados_linha(k, 1), dados_linha(k, 2)) = dados_linha(k, 6) * (- y_km(k));
    Ybus(dados_linha(k, 2), dados_linha(k, 1)) = Ybus(dados_linha(k, 1), dados_linha(k, 2));
    Ybus(dados_linha(k, 2), dados_linha(k, 2)) = Ybus(dados_linha(k, 2), dados_linha(k, 2)) + y_km(k) + j*dados_barra(dados_linha(k, 1), 12) + j*dados_linha(k, 5);
end

Gbus = real(Ybus);  % matriz de condutancias
Bbus = imag(Ybus);  % matriz de susceptancias

%%  matriz S_barra

S_geracao = dados_barra(:, 6) + j*dados_barra(:, 7);
S_consumo = dados_barra(:, 10) + j*dados_barra(:, 11);

S_barra = (S_geracao - S_consumo) / S_base; % matriz potencia aparente em pu

P_barra = real(S_barra);
Q_barra = imag(S_barra);

% incluir trafo defasador (nao existente em todos os exemplos)

%%  chutes iniciais (V e theta)

V_inicial = dados_barra(:, 4);
theta_inicial = dados_barra(:, 5);

for k = 1:1:n_barras
    if dados_barra(k, 2) == 1
        V_inicial(k) = 1;
        theta_inicial(k) = 0;
    elseif dados_barra(k, 2) == 2
        theta_inicial(k) = 0;
    end
end

% matrizes V e theta do sistema com os chutes iniciais aplicados, de acordo com os tipos das barras
V_calc(:,:,1) = V_inicial;
theta_calc(:,:,1) = theta_inicial;

%%  subsistema 1 (calculo de V e theta)

% pre-alocacao das matrizes utilizadas no processo de iteracao
[Vtheta_calc, I_calc, S_calc, P_calc, Q_calc, delta_P, delta_Q, delta_theta, delta_V] = deal(zeros(n_barras, 1));
[delta_PQ, delta_thetaV] = deal([zeros(n_barras, 1); zeros(n_barras, 1)]);
[H, N, M, L] = deal(zeros(n_barras));
[Jacob, Jacob_bn] = deal(zeros(2 * n_barras));

delta_P_esp = P_barra;
delta_Q_esp = Q_barra;

for k = 1:1:n_barras
    if dados_barra(k, 2) == 2
        delta_Q_esp(k) = 0;        
    elseif dados_barra(k, 2) == 3
        delta_P_esp(k) = 0;
        delta_Q_esp(k) = 0;
    end
end

% inicializacao do contador
i = 1;

for k = 1:1:n_barras
    Vtheta_calc(k,:,i) = V_calc(k,:,i)*exp(j*theta_calc(k,:,i));
end

% calculo da matriz de correntes da primeira iteracao
I_calc(:,:,i) = Ybus * Vtheta_calc(:,:,i);
% calculo da matriz potencia da primeira iteracao
S_calc(:,:,i) = Vtheta_calc(:,:,i) .* conj(I_calc(:,:,i));

P_calc(:,:,i) = real(S_calc(:,:,i));    % matriz potencia ativa da primeira iteracao
Q_calc(:,:,i) = imag(S_calc(:,:,i));    % matriz potencia reativa da primeira iteracao

delta_P(:,:,i) = delta_P_esp - P_calc(:,:,i);   % matriz delta P
delta_Q(:,:,i) = delta_Q_esp - Q_calc(:,:,i);   % matriz delta Q
delta_PQ(:,:,i) = [delta_P(:,:,i); delta_Q(:,:,i)];

for k = 1:1:n_barras
    if dados_barra(k, 2) == 2
        delta_PQ(k + n_barras,:,i) = 0;
    elseif dados_barra(k, 2) == 3
        delta_PQ(k,:,i) = 0;
        delta_PQ(k + n_barras,:,i) = 0;
    end
end

% inicio do processo de iteracoes
while max(abs(delta_PQ(:,:,i))) >= erro_admitido
    
    % criacao/atualizacao da metriz jacobiana
    for m = 1:1:n_barras
        for n = 1:1:n_barras
            
            % caso os indices sejam iguais
            if m == n
                H(m, n) = -(V_calc(m,:,i)^2) * Bbus(m, m) - Q_calc(m,:,i);
                N(m, n) = ( P_calc(m,:,i) + (V_calc(m,:,i)^2) * Gbus(m, m) ) / V_calc(m,:,i);
                M(m, n) = -(V_calc(m,:,i)^2) * Gbus(m, m) + P_calc(m,:,i);
                L(m, n) = ( Q_calc(m,:,i) - (V_calc(m,:,i)^2) * Bbus(m, m) ) / V_calc(m,:,i);
            else
                % caso os indices sejam diferentes
                H(m, n) = V_calc(m,:,i) * V_calc(n,:,i) * ( Gbus(m, n) * sin(theta_calc(m,:,i)-theta_calc(n,:,i)) - Bbus(m, n) * cos(theta_calc(m,:,i)-theta_calc(n,:,i)) );
                N(m, n) = V_calc(m,:,i) * ( Gbus(m, n) * cos(theta_calc(m,:,i)-theta_calc(n,:,i)) + Bbus(m, n) * sin(theta_calc(m,:,i)-theta_calc(n,:,i)) );
                M(m, n) = -V_calc(m,:,i) * V_calc(n,:,i) * ( Gbus(m, n) * cos(theta_calc(m,:,i)-theta_calc(n,:,i)) + Bbus(m, n) * sin(theta_calc(m,:,i)-theta_calc(n,:,i)) );
                L(m, n) = V_calc(m,:,i) * ( Gbus(m, n) * sin(theta_calc(m,:,i)-theta_calc(n,:,i)) - Bbus(m, n) * cos(theta_calc(m,:,i)-theta_calc(n,:,i)) );
            end
            
        end
    end
    
    Jacob(:,:,i) = [H N; M L];
    
    % aplicacao do big number de acordo de com o tipo das barras
    for k = 1:1:n_barras
        if dados_barra(k, 2) == 2
            L(k, k) = 10e12;
        elseif dados_barra(k, 2) == 3
            H(k, k) = 10e12;
            L(k, k) = 10e12;
        end
    end
    
    % matriz jacobiana com o big number aplicado
    Jacob_bn(:,:,i) = [H N; M L];
    
    % calculo dos componentes delta theta e delta V
    delta_thetaV(:,:,i) = Jacob_bn(:,:,i) \ delta_PQ(:,:,i);
    delta_theta(:,:,i) = delta_thetaV(1:size(delta_thetaV)/2,:,i);
    delta_V(:,:,i) = delta_thetaV(size(delta_thetaV)/2 + 1:size(delta_thetaV),:,i);
    
    % calculo dos componentes finais delta P e delta Q para analise de uma possivel nova iteracao
    theta_calc(:,:,i + 1) = theta_calc(:,:,i) + delta_theta(:,:,i);
    V_calc(:,:,i + 1) = V_calc(:,:,i) + delta_V(:,:,i);
    
    for k = 1:1:n_barras
        Vtheta_calc(k,:,i + 1) = V_calc(k,:,i + 1)*exp(j*theta_calc(k,:,i + 1));
    end
    
    I_calc(:,:,i + 1) = Ybus * Vtheta_calc(:,:,i + 1);
    S_calc(:,:,i + 1) = Vtheta_calc(:,:,i + 1) .* conj(I_calc(:,:,i + 1));
    
    P_calc(:,:,i + 1) = real(S_calc(:,:,i + 1));
    Q_calc(:,:,i + 1) = imag(S_calc(:,:,i + 1));
    
    delta_P(:,:,i + 1) = delta_P_esp - P_calc(:,:,i + 1);
    delta_Q(:,:,i + 1) = delta_Q_esp - Q_calc(:,:,i + 1);
    delta_PQ(:,:,i + 1) = [delta_P(:,:,i + 1); delta_Q(:,:,i + 1)];
    
    for k = 1:1:n_barras
        if dados_barra(k, 2) == 2
            delta_PQ(k + n_barras,:,i + 1) = 0;
        elseif dados_barra(k, 2) == 3
            delta_PQ(k,:,i + 1) = 0;
            delta_PQ(k + n_barras,:,i + 1) = 0;
        end
    end
    
    % acrescimo de uma unidade ao contador
    i = i + 1;
end


