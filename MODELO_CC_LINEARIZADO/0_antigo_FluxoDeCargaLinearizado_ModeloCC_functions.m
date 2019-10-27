%%% FLUXO DE CARGA LINEARIZADO - MODELO CC

clear
close all
format longEng
clc

% No TB G (V)     (Ang)      (Pg)    (Qg)  (Qn)    (Qm)     (Pl)   (Ql)    (A(Cf))  bshbar 
dados_barra = [                                                                           
   1 3 1 1.0         .0       232.4   .0  -999.9  9999.9     0.0    0.0      0.0      0
   2 1 1 1.0         .0       40.0    .0  -999.9  9999.9     21.7   1.6      0.0      0
   3 1 1 1.0         .0       .0      .0  -999.9  9999.9     94.2   1.5      1.1      0
   4 1 1 1.0         .0       .0      .0  -999.9  9999.9     47.8   0.8      1.2      0
   5 1 1 1.0         .0       .0      .0  -999.9  9999.9     7.6    1.2      0.0      0
   6 1 1 1.0         .0       .0      .0  -999.9  9999.9     11.2   2.7      0.0      0
   7 1 1 1.0         .0       .0      .0  -999.9  9999.9     0.0    3.0      1.2      0
   8 1 1 1.0         .0       .0      .0  -999.9  9999.9     0.0    0.9      0.0      0
   9 1 1 1.0         .0       .0      .0  -999.9  9999.9     29.5   0.1      0.6      0.19
  10 1 1 1.0         .0       .0      .0  -999.9  9999.9     9.0    2.0      3.7      0
  11 1 1 1.0         .0       .0      .0  -999.9  9999.9     3.5    0.9      0.0      0
  12 1 1 1.0         .0       .0      .0  -999.9  9999.9     6.1    0.7      1.8      0
  13 1 1 1.0         .0       .0      .0  -999.9  9999.9     13.5   0.9      0.0      0
  14 1 1 1.0         .0       .0      .0  -999.9  9999.9     14.9   1.0      1.8      0
];
% TB = 1: carga ; 3: referencia


% FROM TO     %R(pu)   %X(pu)    %Bsh(pu)     %Tap    %Def(graus)                              CH
dados_linha = [
   1    2    0.01938  0.05917     0.0264      1.000     .000     .000    .0     900     .0     7
   1    5    0.05403  0.22304     0.0246      1.000     .000     .000    .0     900     .0     7
   2    3    0.04699  0.19797     0.0219      1.000     .000     .000    .0     900     .0     7
   2    4    0.05811  0.17632     0.0170      1.000     .000     .000    .0     900     .0     7
   2    5    0.05695  0.17388     0.0173      1.000     .000     .000    .0     900     .0     7
   3    4    0.06701  0.17103     0.0064      1.000     .000     .000    .0     900     .0     7
   4    5    0.01335  0.04211     0.0000      1.000     .000     .000    .0     900     .0     7
   4    7    0.00000  0.20912     0.0000      1.000     .000     .000    .0     900     .0     7
   4    9    0.00000  0.55618     0.0000      1.000     .000     .000    .0     900     .0     7
   5    6    0.00000  0.25202     0.0000      1.000     .000     .000    .0     900     .0     7
   6   11    0.09498  0.19890     0.0000      1.000     .000     .000    .0     900     .0     7
   6   12    0.12291  0.25581     0.0000      1.000     .000     .000    .0     900     .0     7
   6   13    0.06615  0.13027     0.0000      1.000     .000     .000    .0     900     .0     7
   7    8    0.00000  0.17615     0.0000      1.000     .000     .000    .0     900     .0     7
   7    9    0.00000  0.11001     0.0000      1.000     .000     .000    .0     900     .0     7
   9   10    0.03181  0.08450     0.0000      1.000     .000     .000    .0     900     .0     7
   9   14    0.12711  0.27038     0.0000      1.000     .000     .000    .0     900     .0     7
   10  11    0.08205  0.19207     0.0000      1.000     .000     .000    .0     900     .0     7
   12  13    0.22092  0.19988     0.0000      1.000     .000     .000    .0     900     .0     7
   13  14    0.17093  0.34802     0.0000      1.000     .000     .000    .0     900     .0     7
 ];

n_barras = size(dados_barra, 1);
P_geracao = dados_barra(:, 6);
P_consumo = dados_barra(:, 10);
r_km = dados_linha(:, 3);
x_km = dados_linha(:, 4);
tap_km = dados_linha(:, 6);
def_km = dados_linha(:, 7);

B = matriz_reatancia(n_barras, dados_barra, dados_linha);

P_barra = P_geracao - P_consumo;

P_barra = alteracao_Pbarra_trafo_defasador(P_barra, dados_linha);

theta = B \ P_barra;

[g_km, perdas, P_barra_novo] = perdas_e_novo_Pbarra(P_barra, dados_linha, theta);

total_perdas = sum(sum(perdas));

theta_novo = B \ P_barra_novo;

fluxo_potencia = calculo_fluxo_potencia(n_barras, dados_linha, theta_novo);

%%%%% functions %%%%%

function B = matriz_reatancia(n_barras, dados_barra, dados_linha)
%matriz_reatancia - calcula a matriz de reatancias e aplica a tecnica do 'big number'
    
b = zeros(n_barras);

% agrupa as reatancias entre linhas de acordo com os indices (i, j)
for k = 1:1:size(dados_linha, 1)
    b(dados_linha(k, 1), dados_linha(k, 2)) = dados_linha(k, 6) * (dados_linha(k, 4))^(-1);
    b(dados_linha(k, 2), dados_linha(k, 1)) = b(dados_linha(k, 1), dados_linha(k, 2));
end

% criacao da matriz de reatancias
B = -b;

% somatorio das reatancias em indices iguais (k, k) e inclusao do 'big number' na barra de referencia
for k = 1:1:n_barras
    B(k, k) = sum(b(k, :));
    
    if dados_barra(k, 2) == 3
        B(k, k) = 10e9;
    end
end

end

function P_barra = alteracao_Pbarra_trafo_defasador(P_barra, dados_linha)
%alteracao_Pbarra_trafo_defasador - alteracao de P_barra caso haja trafo defasador

for k = 1:1:size(dados_linha, 1)
    P_barra(dados_linha(k, 1)) = P_barra(dados_linha(k, 1)) - dados_linha(k, 7) / dados_linha(k, 4);
    P_barra(dados_linha(k, 2)) = P_barra(dados_linha(k, 2)) + dados_linha(k, 7) / dados_linha(k, 4);
end

end

function [g_km, perdas, P_barra_novo] = perdas_e_novo_Pbarra(P_barra, dados_linha, theta)
%perdas_e_novo_Pbarra - calcula as perdas e a nova matriz P_barra

P_barra_novo = P_barra;

for k = 1:1:size(dados_linha, 1)
    if dados_linha(k, 3) == 0
        g_km(k) = 0;
    else
        g_km(k) = 1 / dados_linha(k, 3);
    end
    
    perdas(dados_linha(k, 1), dados_linha(k, 2)) = g_km(k) * (theta(dados_linha(k, 1)) - theta(dados_linha(k, 2)))^2;
    P_barra_novo(dados_linha(k, 1)) = P_barra_novo(dados_linha(k, 1)) - perdas(dados_linha(k, 1), dados_linha(k, 2)) / 2;
    P_barra_novo(dados_linha(k, 2)) = P_barra_novo(dados_linha(k, 2)) - perdas(dados_linha(k, 1), dados_linha(k, 2)) / 2;
end

end

function fluxo_potencia = calculo_fluxo_potencia(n_barras, dados_linha, theta_novo)
%calculo_fluxo_potencia - calculo dos fluxos de potencia entre barras

fluxo_potencia = zeros(n_barras);

for k = 1:1:size(dados_linha, 1)
    fluxo_potencia(dados_linha(k, 1), dados_linha(k, 2)) = dados_linha(k, 6) * (theta_novo(dados_linha(k, 1)) - theta_novo(dados_linha(k, 2)) + dados_linha(k, 7)) / dados_linha(k, 4);
    fluxo_potencia(dados_linha(k, 2), dados_linha(k, 1)) =  - fluxo_potencia(dados_linha(k, 1), dados_linha(k, 2));
end
    
end