%%% FLUXO DE CARGA NAO LINEARIZADO - NEWTON-RAPHSON - METODO ATRAVES DO CALCULO DAS CORRENTES %%%

%%

clear
close all
format long
clc

%%  leitura do arquivo de exemplo

% Sistema_4_barras_Monticelli;
% Sistema_14_barra_2;
% Sistema_24_barras;
% Sistema_24_barras_naocorrigido;
% Sistema_33_barras;
Sistema_107_barras;

%%  declaracao de variaveis

dados_barra = DBAR; % matriz de entrada com informacoes das barras do sistema
dados_linha = DLIN; % matriz de entrada com informacoes das linhas do sistema
S_base = PB;    % potencia base do sistema

n_barras = size(dados_barra, 1);    % numero de barras do sistema

erro_admitido = 1e-6;  % erro admitido para fim das iteracoes

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
    Ybus(dados_linha(k, 1), dados_linha(k, 1)) = Ybus(dados_linha(k, 1), dados_linha(k, 1)) + dados_linha(k, 6)^2 * y_km(k) + j*dados_linha(k, 5);
    Ybus(dados_linha(k, 1), dados_linha(k, 2)) = Ybus(dados_linha(k, 1), dados_linha(k, 2)) + dados_linha(k, 6) * exp(- j * dados_linha(k, 7)) * (- y_km(k));
    Ybus(dados_linha(k, 2), dados_linha(k, 1)) = Ybus(dados_linha(k, 2), dados_linha(k, 1)) + dados_linha(k, 6) * exp(j * dados_linha(k, 7)) * (- y_km(k));
    Ybus(dados_linha(k, 2), dados_linha(k, 2)) = Ybus(dados_linha(k, 2), dados_linha(k, 2)) + y_km(k) + j*dados_linha(k, 5);
end

for k = 1:1:n_barras
    Ybus(k,k) = Ybus(k,k) + j*dados_barra(k, 12);
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

V_inicial = dados_barra(:, 4);  % valores de tensao indicados na matriz de entrada
theta_inicial = dados_barra(:, 5);  % valores de abertura angular indicados na matriz de entrada

% aplicacao dos chutes iniciais nas barras em que e necessario o calculo de V e theta
for k = 1:1:n_barras
    if dados_barra(k, 2) == 1   % barra PQ
        V_inicial(k) = 1;
        theta_inicial(k) = 0;
    elseif dados_barra(k, 2) == 2   % barra PV
        theta_inicial(k) = 0;
    end
end

% matrizes V e theta do sistema com os chutes iniciais aplicados, de acordo com os tipos das barras
V_calc(:,:,1) = V_inicial;
theta_calc(:,:,1) = theta_inicial;

%%  subsistema 1 (calculo de V e theta)

% pre-alocacao das matrizes utilizadas no processo de iteracao
[Vtheta_calc, I_calc, S_calc, P_calc, Q_calc, delta_P, delta_Q, delta_theta, delta_V, graph_delta_P, graph_delta_Q] = deal(zeros(n_barras, 1));
[delta_PQ, delta_thetaV] = deal([zeros(n_barras, 1); zeros(n_barras, 1)]);
[H, N, M, L] = deal(zeros(n_barras));
[Jacob, Jacob_bn] = deal(zeros(2 * n_barras));

P_esp = P_barra;
Q_esp = Q_barra;

% aplicacao de zeros nos indices das barras que nao necessitam de calculo a fim de evitar problemas
for k = 1:1:n_barras
    if dados_barra(k, 2) == 2   % barra PV
        Q_esp(k) = 0;        
    elseif dados_barra(k, 2) == 3   % barra V-theta
        P_esp(k) = 0;
        Q_esp(k) = 0;
    end
end

% inicializacao do contador
i = 1;

% juncao de V e theta (coordenada polar) para calculo das correntes e potencias
for k = 1:1:n_barras
    Vtheta_calc(k,:,i) = V_calc(k,:,i) * exp(j*theta_calc(k,:,i));
end

% calculo da matriz de correntes
I_calc(:,:,i) = Ybus * Vtheta_calc(:,:,i);
% calculo da matriz potencia aparente
S_calc(:,:,i) = Vtheta_calc(:,:,i) .* conj(I_calc(:,:,i));

P_calc(:,:,i) = real(S_calc(:,:,i));    % matriz potencia ativa
Q_calc(:,:,i) = imag(S_calc(:,:,i));    % matriz potencia reativa

delta_P(:,:,i) = P_esp - P_calc(:,:,i);   % matriz delta P
delta_Q(:,:,i) = Q_esp - Q_calc(:,:,i);   % matriz delta Q
delta_PQ(:,:,i) = [delta_P(:,:,i); delta_Q(:,:,i)];

% aplicacao de zeros nos indices das barras que nao necessitam de calculo a fim de evitar problemas
for k = 1:1:n_barras
    if dados_barra(k, 2) == 2   % barra PV
        delta_PQ(k + n_barras,:,i) = 0;
    elseif dados_barra(k, 2) == 3   % barra V-theta
        delta_PQ(k,:,i) = 0;
        delta_PQ(k + n_barras,:,i) = 0;
    end
end

% inicio do processo de iteracoes
while max(abs(delta_PQ(:,:,i))) >= erro_admitido
    
    % criacao/atualizacao da matriz jacobiana
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
        if dados_barra(k, 2) == 2   % barra PV
            L(k, k) = 10e12;
        elseif dados_barra(k, 2) == 3   % barra V-theta
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
    
    % calculo dos componentes theta e V
    theta_calc(:,:,i + 1) = theta_calc(:,:,i) + delta_theta(:,:,i);
    V_calc(:,:,i + 1) = V_calc(:,:,i) + delta_V(:,:,i);
    
    % juncao de V e theta (coordenada polar) para calculo das correntes e potencias
    for k = 1:1:n_barras
        Vtheta_calc(k,:,i + 1) = V_calc(k,:,i + 1) * exp(j*theta_calc(k,:,i + 1));
    end
    
    % novos calculos de correntes e potencias aparentes
    I_calc(:,:,i + 1) = Ybus * Vtheta_calc(:,:,i + 1);
    S_calc(:,:,i + 1) = Vtheta_calc(:,:,i + 1) .* conj(I_calc(:,:,i + 1));
    
    P_calc(:,:,i + 1) = real(S_calc(:,:,i + 1));
    Q_calc(:,:,i + 1) = imag(S_calc(:,:,i + 1));
    
    % calculo final de delta P e delta Q para analise de erro
    delta_P(:,:,i + 1) = P_esp - P_calc(:,:,i + 1);
    delta_Q(:,:,i + 1) = Q_esp - Q_calc(:,:,i + 1);
    delta_PQ(:,:,i + 1) = [delta_P(:,:,i + 1); delta_Q(:,:,i + 1)];
    
    % preparacao para possivel nova iteracao
    for k = 1:1:n_barras
        if dados_barra(k, 2) == 2   % barra PV
            delta_PQ(k + n_barras,:,i + 1) = 0;
        elseif dados_barra(k, 2) == 3   % barra V-theta
            delta_PQ(k,:,i + 1) = 0;
            delta_PQ(k + n_barras,:,i + 1) = 0;
        end
    end
    
    % acrescimo de uma unidade ao contador
    i = i + 1;
end

% numero final de iteracoes
iteracoes = i - 1;

for k = 1:1:i
    graph_delta_P(:,k) = delta_PQ(1:size(delta_PQ)/2,:,k);
    graph_delta_Q(:,k) = delta_PQ(size(delta_PQ)/2+1:size(delta_PQ),:,k);
end

%%  subsistema 2 (calculo de P e Q + fluxo de potencia e perdas)

% pre-alocacao de matrizes
Vtheta_calc_final  = zeros(n_barras, 1);
[P_km, Q_km, perdas_P, perdas_Q] = deal(zeros(n_barras));

% juncao de V e theta (coordenada polar) para calculo das correntes e potencias
for k = 1:1:n_barras
    Vtheta_calc_final(k) = V_calc(k,:,i) * exp(j*theta_calc(k,:,i));
end

% novos calculos de correntes e potencias aparentes
I_calc_final = Ybus * Vtheta_calc_final;
S_calc_final = Vtheta_calc_final .* conj(I_calc_final);

P_calc_final = real(S_calc_final);
Q_calc_final = imag(S_calc_final);

% calculo dos fluxos de potencia entre as barras que possuem conexao
for k = 1:1:size(dados_linha, 1)
    P_km(dados_linha(k, 1), dados_linha(k, 2)) = ( dados_linha(k, 6) * V_calc(dados_linha(k, 1),:,i) )^2 * Gbus(dados_linha(k, 1), dados_linha(k, 2)) ...
        - ( dados_linha(k, 6) * V_calc(dados_linha(k, 1),:,i) ) * V_calc(dados_linha(k, 2),:,i) * Gbus(dados_linha(k, 1), dados_linha(k, 2)) * cos(theta_calc(dados_linha(k, 1),:,i) - theta_calc(dados_linha(k, 2),:,i) + dados_linha(k, 7)) ...
        - ( dados_linha(k, 6) * V_calc(dados_linha(k, 1),:,i) ) * V_calc(dados_linha(k, 2),:,i) * Bbus(dados_linha(k, 1), dados_linha(k, 2)) * sin(theta_calc(dados_linha(k, 1),:,i) - theta_calc(dados_linha(k, 2),:,i) + dados_linha(k, 7));
    P_km(dados_linha(k, 2), dados_linha(k, 1)) = ( dados_linha(k, 6) * V_calc(dados_linha(k, 2),:,i) )^2 * Gbus(dados_linha(k, 2), dados_linha(k, 1)) ...
        - ( dados_linha(k, 6) * V_calc(dados_linha(k, 2),:,i) ) * V_calc(dados_linha(k, 1),:,i) * Gbus(dados_linha(k, 2), dados_linha(k, 1)) * cos(theta_calc(dados_linha(k, 2),:,i) - theta_calc(dados_linha(k, 1),:,i) + dados_linha(k, 7)) ...
        - ( dados_linha(k, 6) * V_calc(dados_linha(k, 2),:,i) ) * V_calc(dados_linha(k, 1),:,i) * Bbus(dados_linha(k, 2), dados_linha(k, 1)) * sin(theta_calc(dados_linha(k, 2),:,i) - theta_calc(dados_linha(k, 1),:,i) + dados_linha(k, 7));
    Q_km(dados_linha(k, 1), dados_linha(k, 2)) = -( dados_linha(k, 6) * V_calc(dados_linha(k, 1),:,i) )^2 * ( Bbus(dados_linha(k, 1), dados_linha(k, 2)) + dados_linha(k, 5) ) ...
        + ( dados_linha(k, 6) * V_calc(dados_linha(k, 1),:,i) ) * V_calc(dados_linha(k, 2),:,i) * Bbus(dados_linha(k, 1), dados_linha(k, 2)) * cos(theta_calc(dados_linha(k, 1),:,i) - theta_calc(dados_linha(k, 2),:,i) + dados_linha(k, 7)) ...
        - ( dados_linha(k, 6) * V_calc(dados_linha(k, 1),:,i) ) * V_calc(dados_linha(k, 2),:,i) * Gbus(dados_linha(k, 1), dados_linha(k, 2)) * sin(theta_calc(dados_linha(k, 1),:,i) - theta_calc(dados_linha(k, 2),:,i) + dados_linha(k, 7));
    Q_km(dados_linha(k, 2), dados_linha(k, 1)) = -( dados_linha(k, 6) * V_calc(dados_linha(k, 2),:,i) )^2 * ( Bbus(dados_linha(k, 2), dados_linha(k, 1)) + dados_linha(k, 5) ) ...
        + ( dados_linha(k, 6) * V_calc(dados_linha(k, 2),:,i) ) * V_calc(dados_linha(k, 1),:,i) * Bbus(dados_linha(k, 2), dados_linha(k, 1)) * cos(theta_calc(dados_linha(k, 2),:,i) - theta_calc(dados_linha(k, 1),:,i) + dados_linha(k, 7)) ...
        - ( dados_linha(k, 6) * V_calc(dados_linha(k, 2),:,i) ) * V_calc(dados_linha(k, 1),:,i) * Gbus(dados_linha(k, 2), dados_linha(k, 1)) * sin(theta_calc(dados_linha(k, 2),:,i) - theta_calc(dados_linha(k, 1),:,i) + dados_linha(k, 7));
    
    perdas_P(dados_linha(k, 1), dados_linha(k, 2)) = P_km(dados_linha(k, 1), dados_linha(k, 2)) + P_km(dados_linha(k, 2), dados_linha(k, 1)); 
    perdas_Q(dados_linha(k, 1), dados_linha(k, 2)) = Q_km(dados_linha(k, 1), dados_linha(k, 2)) + Q_km(dados_linha(k, 2), dados_linha(k, 1));
end

P_km = sparse(P_km);
Q_km = sparse(Q_km);
perdas_P = sparse(perdas_P);
perdas_Q = sparse(perdas_Q);

%%  impressao de resultados

figure
subplot(2,1,1)
bar(graph_delta_P)
title('Evolucao \DeltaP / iteracao')
xlabel('nº da barra')
ylabel('Residuo \DeltaP (pu)')
grid
subplot(2,1,2)
bar(graph_delta_Q)
title('Evolucao \DeltaQ / iteracao')
xlabel('nº da barra')
ylabel('Residuo \DeltaQ (pu)')
grid

format short

fprintf('---------- RESULTADOS ----------\n\n\n')


fprintf('*** SUBSISTEMA 1 ***\n\n')

disp('Abertura angular theta (°):')
display(rad2deg(theta_calc(:,:,i)))
figure
bar(rad2deg(theta_calc(:,:,i)))
xlabel('nº da barra')
ylabel('Theta (°)')
grid

disp('Modulo da tensao (pu):')
display(V_calc(:,:,i))
figure
bar(V_calc(:,:,i))
xlabel('nº da barra')
ylabel('Tensao (pu)')
grid


fprintf('\n*** SUBSISTEMA 2 ***\n\n')

disp('Potencia ativa (MW):')
display(P_calc_final * S_base)

disp('Potencia reativa (MVAr):')
display(Q_calc_final * S_base)

disp('Perdas de potencia ativa (MW):')
display(perdas_P * S_base)

disp('Perdas de potencia reativa (MVAr):')
display(perdas_Q * S_base)
