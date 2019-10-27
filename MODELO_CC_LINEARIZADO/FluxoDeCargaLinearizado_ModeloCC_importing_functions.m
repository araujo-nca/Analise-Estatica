%%% FLUXO DE CARGA LINEARIZADO - MODELO CC %%%

%%

clear
close all
format longEng
clc

%% leitura do arquivo de exemplo

Sistema_2_barras;
% Sistema_3_barras;
% Sistema_3_barras_M;
% Sistema_14_barras_Winter;

%% declaracao de variaveis

dados_barra = DBAR; 
dados_linha = DLIN;
P_base = PB;

n_barras = size(dados_barra, 1);
P_geracao = dados_barra(:, 6);
P_consumo = dados_barra(:, 10);

% apenas para entendimento da matriz
% r_km = dados_linha(:, 3);
% x_km = dados_linha(:, 4);
% tap_km = dados_linha(:, 6);
% def_km = dados_linha(:, 7);

%% criacao da matriz de reatancia e aplicacao do 'big number' --- funcao matriz_reatancia

[B, B_bignumber] = matriz_reatancia(n_barras, dados_barra, dados_linha);

%% matriz P_barra em pu e sua alteracao caso haja algum trafo defasador no sistema --- funcao alteracao_Pbarra_trafo_defasador

P_barra_pu = (P_geracao - P_consumo) / P_base;

P_barra_pu = alteracao_Pbarra_trafo_defasador(P_barra_pu, dados_linha);

%% primeira resolucao do sistema sem a consideracao das perdas

theta = B_bignumber \ P_barra_pu;

%% fluxo de potencia sem perdas --- funcao calculo_fluxo_potencia

fluxo_potencia = calculo_fluxo_potencia(n_barras, dados_linha, theta);

%% calculo de g_km, perdas e o novo P_barra com as perdas --- funcao perdas_e_novo_Pbarra

[g_km, perdas, total_perdas, P_barra_pu_novo] = perdas_e_novo_Pbarra(P_barra_pu, dados_linha, theta);

%% segunda resolucao do sistema com a consideracao das perdas

theta_novo = B_bignumber \ P_barra_pu_novo;

%% fluxo de potencia com perdas --- funcao calculo_fluxo_potencia

fluxo_potencia_perdas = calculo_fluxo_potencia(n_barras, dados_linha, theta_novo);

%% impressao de resultados

format short

fprintf('---------- RESULTADOS ----------\n\n')

disp('Matriz de reatancia do sistema:')
display(B)

disp('Injecao de potencia nas barras (pu):')
display(P_barra_pu)

disp('Angulos theta sem a inclusao das perdas (rad):')
theta_deg = rad2deg(theta);
display(theta_deg)

disp('Fluxo de potencia entre as barras sem as perdas (MW):')
fluxo_potencia_MW  = fluxo_potencia * P_base;
display(fluxo_potencia_MW)

disp('Perdas (MW) em cada LT e o total de perdas (MW):')
perdas_MW = perdas * P_base;
total_perdas_MW = total_perdas * P_base;
display(perdas_MW)
display(total_perdas_MW)

disp('Injecao de potencia nas barras apos a inclusao das perdas (pu):')
display(P_barra_pu_novo)

disp('Novos angulos theta apos a inclusao das perdas (rad):')
theta_novo_deg = rad2deg(theta_novo);
display(theta_novo_deg)

disp('Fluxo de potencia entre as barras com as perdas (MW):')
fluxo_potencia_perdas_MW  = fluxo_potencia_perdas * P_base;
display(fluxo_potencia_perdas_MW)
