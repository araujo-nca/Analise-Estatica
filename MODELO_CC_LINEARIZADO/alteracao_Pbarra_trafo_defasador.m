function P_barra = alteracao_Pbarra_trafo_defasador(P_barra, dados_linha)
%alteracao_Pbarra_trafo_defasador - alteracao de P_barra caso haja trafo defasador

% aplicacao do trafo defasador em ambas as barras cujo o trafo afeta (representa geracao em uma e carga na outra de acordo com o sinal)
for k = 1:1:size(dados_linha, 1)
    P_barra(dados_linha(k, 1)) = P_barra(dados_linha(k, 1)) - deg2rad(dados_linha(k, 7)) / dados_linha(k, 4);
    P_barra(dados_linha(k, 2)) = P_barra(dados_linha(k, 2)) + deg2rad(dados_linha(k, 7)) / dados_linha(k, 4);
end

end