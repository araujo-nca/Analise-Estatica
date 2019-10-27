function fluxo_potencia = calculo_fluxo_potencia(n_barras, dados_linha, theta_novo)
%calculo_fluxo_potencia - calculo dos fluxos de potencia entre barras

for k = 1:1:size(dados_linha, 1)
    fluxo_potencia(dados_linha(k, 1), dados_linha(k, 2)) = dados_linha(k, 6) * (theta_novo(dados_linha(k, 1)) - theta_novo(dados_linha(k, 2)) + deg2rad(dados_linha(k, 7))) / dados_linha(k, 4);
end

fluxo_potencia = sparse(fluxo_potencia);

end