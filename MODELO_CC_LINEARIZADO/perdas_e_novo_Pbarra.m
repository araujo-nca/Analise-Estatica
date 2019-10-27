function [g_km, perdas, total_perdas, P_barra_novo] = perdas_e_novo_Pbarra(P_barra, dados_linha, theta)
%perdas_e_novo_Pbarra - calcula as perdas e a nova matriz P_barra

P_barra_novo = P_barra;

for k = 1:1:size(dados_linha, 1)
    g_km(k) = dados_linha(k, 3) / ( (dados_linha(k, 3))^2 + (dados_linha(k, 4))^2 );
    
    perdas(dados_linha(k, 1), dados_linha(k, 2)) = g_km(k) * ( theta(dados_linha(k, 1)) - theta(dados_linha(k, 2)) )^2;

    P_barra_novo(dados_linha(k, 1)) = P_barra_novo(dados_linha(k, 1)) - perdas(dados_linha(k, 1), dados_linha(k, 2)) / 2;
    P_barra_novo(dados_linha(k, 2)) = P_barra_novo(dados_linha(k, 2)) - perdas(dados_linha(k, 1), dados_linha(k, 2)) / 2;
end

total_perdas = sum(sum(perdas));
perdas = sparse(perdas);

end