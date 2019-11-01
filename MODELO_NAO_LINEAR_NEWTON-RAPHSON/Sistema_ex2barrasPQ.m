% ****   Sistema exemplo aula 2 barras PQ  ****

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
