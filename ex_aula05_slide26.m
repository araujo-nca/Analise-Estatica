%%% FLUXO DE CARGA LINEARIZADO - MODELO CC

clear
close all
format longEng
clc

% No TB G (V)     (Ang)      (Pg)    (Qg)  (Qn)    (Qm)     (Pl)   (Ql)    (A(Cf))  bshbar 
dados_barra = [                                                                           
   1 3 1 1.0         .0       0.4   .0  -999.9  9999.9     0.0    0.0      0.0      0
   2 1 1 1.0         .0       0.4   .0  -999.9  9999.9     0.0    1.6      0.0      0
   3 1 1 1.0         .0       .0    .0  -999.9  9999.9     0.8    1.5      1.1      0
];
% TB = 1: carga ; 3: referencia


% FROM TO  %R(pu)   %X(pu)  %Bsh(pu)     %Tap    %Def(graus)                              CH
dados_linha = [
   1    2    0       2       0.0264      1.200     .000     .000    .0     900     .0     7
   1    3    0       1       0.0246      1.000     .000     .000    .0     900     .0     7
   2    3    0       1/2     0.0219      1.000     .000     .000    .0     900     .0     7
 ];

n_barras = size(dados_barra, 1);
x_km = dados_linha(:, 4);

b = zeros(n_barras);

for k = 1:1:size(dados_linha, 1)
    b(dados_linha(k, 1), dados_linha(k, 2)) = dados_linha(k, 6) * (dados_linha(k, 4))^(-1);
    b(dados_linha(k, 2), dados_linha(k, 1)) = b(dados_linha(k, 1), dados_linha(k, 2));
end

B = -b;

for k = 1:1:n_barras
    B(k, k) = sum(b(k, :));
    
    if dados_barra(k, 2) == 3
        B(k, k) = 10e9;
    end
end

P_barra = dados_barra(:, 6) - dados_barra(:, 10);
P_barra_novo = P_barra;

theta = B \ P_barra;

for k = 1:1:size(dados_linha, 1)
    if dados_linha(k, 3) == 0
        g(k) = 0;
    else
        g(k) = 1 / dados_linha(k, 3);
    end
    
    perdas(dados_linha(k, 1), dados_linha(k, 2)) = g(k) * (theta(dados_linha(k, 1)) - theta(dados_linha(k, 2)))^2;
    P_barra_novo(dados_linha(k, 1)) = P_barra_novo(dados_linha(k, 1)) - perdas(dados_linha(k, 1), dados_linha(k, 2)) / 2;
    P_barra_novo(dados_linha(k, 2)) = P_barra_novo(dados_linha(k, 2)) - perdas(dados_linha(k, 1), dados_linha(k, 2)) / 2;
end

total_perdas = sum(sum(perdas));

theta_novo = B \ P_barra_novo;

fluxo_potencia = zeros(n_barras);

for k = 1:1:size(dados_linha, 1)
    fluxo_potencia(dados_linha(k, 1), dados_linha(k, 2)) = dados_linha(k, 6) * (theta_novo(dados_linha(k, 1)) - theta_novo(dados_linha(k, 2))) / dados_linha(k, 4);
    fluxo_potencia(dados_linha(k, 2), dados_linha(k, 1)) =  - fluxo_potencia(dados_linha(k, 1), dados_linha(k, 2));
end