function [B, B_bignumber] = matriz_reatancia(n_barras, dados_barra, dados_linha)
%matriz_reatancia - calcula a matriz de reatancias e aplica a tecnica do 'big number'
    
b = zeros(n_barras);

% agrupa as reatancias entre linhas de acordo com os indices (i, j)
for k = 1:1:size(dados_linha, 1)
    b(dados_linha(k, 1), dados_linha(k, 2)) = dados_linha(k, 6) / dados_linha(k, 4);
    b(dados_linha(k, 2), dados_linha(k, 1)) = b(dados_linha(k, 1), dados_linha(k, 2));
end

% criacao da matriz de reatancias
B = -b;
B_bignumber = B;

% somatorio das reatancias em indices iguais (k, k) e inclusao do 'big number' na barra de referencia
for k = 1:1:n_barras
    B(k, k) = sum(b(k, :));
    B_bignumber(k, k) = B(k, k);
    
    if dados_barra(k, 2) == 3
        B_bignumber(k, k) = 10e9;
    end
end

end