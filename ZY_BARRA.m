%   Z E Y BARRA DE UM SISTEMA DE TRANSMISSAO

clear
clc

a = input('Quantas barras possui o sistema?   ');
disp(' ')

y = inf.*ones(a);

for m = 1:1:a
    for n = 1:1:a
        if (y(m,n) == inf)
            y(m,n) = input(['Valor da admitancia y[',num2str(m),',',num2str(n),']: ']);
            y(n,m) = y(m,n);
        end
    end
end

Ybarra = -y;

for m = 1:1:a
    Ybarra(m,m) = sum(y(m,:));
end

Zbarra = inv(Ybarra);

display(y)
display(Ybarra)
display(Zbarra)