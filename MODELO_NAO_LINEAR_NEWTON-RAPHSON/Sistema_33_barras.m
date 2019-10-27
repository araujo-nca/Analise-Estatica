% ****   Sistema teste de 33 BARRAS - IEEE   ****                             

% ------------------------------ DADOS DE BARRA -----------------------------------

%  No TB G (V)     (Ang)     (Pg)    (Qg)   (Qn)    (Qm)      (Pl)    (Ql)    (A(Cf)
DBAR = [
   1 1 1 0.9927      .0       .0      .0  -999.9  9999.9     100.0    60.0      0.0     0
   2 1 1 0.9574      .0       .0      .0  -999.9  9999.9     90.0     40.0      0.0     0
   3 1 1 0.9374      .0       .0      .0  -999.9  9999.9     120.0    80.0      0.0     0
   4 1 1 0.9176      .0       .0      .0  -999.9  9999.9     60.0     30.0      0.0     0
   5 1 1 0.8707      .0       .0      .0  -999.9  9999.9     60.0     20.0      0.0     0
   6 1 1 0.8641      .0       .0      .0  -999.9  9999.9     200.0    100.0     0.0     0
   7 1 1 0.8550      .0       .0      .0  -999.9  9999.9     200.0    100.0     0.0     0
   8 1 1 0.8432      .0       .0      .0  -999.9  9999.9     60.0     20.0      0.0     0
   9 1 1 0.8324      .0       .0      .0  -999.9  9999.9     60.0     20.0      0.0     0 
  10 1 1 0.8308      .0       .0      .0  -999.9  9999.9     45.0     30.0      0.0     0
  11 1 1 0.8280      .0       .0      .0  -999.9  9999.9     60.0     35.0      0.0     0
  12 1 1 0.8167      .0       .0      .0  -999.9  9999.9     60.0     35.0      0.0     0
  13 1 1 0.8125      .0       .0      .0  -999.9  9999.9     120.0    80.0      0.0     0 
  14 1 1 0.8099      .0       .0      .0  -999.9  9999.9     60.0     10.0      0.0     0
  15 1 1 0.8074      .0       .0      .0  -999.9  9999.9     60.0     20.0      0.0     0
  16 1 1 0.8037      .0       .0      .0  -999.9  9999.9     60.0     20.0      0.0     0
  17 1 1 0.8026      .0       .0      .0  -999.9  9999.9     90.0     40.0      0.0     0
  18 1 1 0.9916      .0       .0      .0  -999.9  9999.9     90.0     40.0      0.0     0
  19 1 1 0.9845      .0       .0      .0  -999.9  9999.9     90.0     40.0      0.0     0
  20 1 1 0.9831      .0       .0      .0  -999.9  9999.9     90.0     40.0      0.0     0
  21 1 1 0.9818      .0       .0      .0  -999.9  9999.9     90.0     40.0      0.0     0
  22 1 1 0.9504      .0       .0      .0  -999.9  9999.9     90.0     50.0      0.0     0
  23 1 1 0.9373      .0       .0      .0  -999.9  9999.9     420.0    200.0     0.0     0
  24 1 1 0.9309      .0       .0      .0  -999.9  9999.9     420.0    200.0     0.0     0
  25 1 1 0.8643      .0       .0      .0  -999.9  9999.9     60.0     25.0      0.0     0
  26 1 1 0.8557      .0       .0      .0  -999.9  9999.9     60.0     25.0      0.0     0
  27 1 1 0.8201      .0       .0      .0  -999.9  9999.9     60.0     20.0      0.0     0
  28 1 1 0.7945      .0       .0      .0  -999.9  9999.9     120.0    70.0      0.0     0
  29 1 1 0.7816      .0       .0      .0  -999.9  9999.9     200.0    600.0     0.0     0
  30 1 1 0.7739      .0       .0      .0  -999.9  9999.9     150.0    70.0      0.0     0
  31 1 1 0.7723      .0       .0      .0  -999.9  9999.9     210.0    100.0     0.0     0
  32 1 1 0.7717      .0       .0      .0  -999.9  9999.9     60.0     40.0      0.0     0
  33 3 1 1.0         .0       .0      .0  -999.9  9999.9     0.0      0.0       0.0     0	];
        
% kW -> MW
DBAR(:,10) = DBAR(:,10) .* 1e-3;
DBAR(:,11) = DBAR(:,11) .* 1e-3;
DBAR(:,12) = DBAR(:,12) .* 1e-3;

% Tipos de barra: 1 - carga (PQ), 2 - geracao (PV), 3 - referencia (V-theta)

PB = 100;

% ------------------------------ DADOS DE LINHA -----------------------------------

 DLIN = [
  33    1    0.0922   0.0470      0.0      .000     .000     .000    .0     900     .0     7
   1    2    0.4930   0.2511      0.0      .000     .000     .000    .0     900     .0     7
   2    3    0.3660   0.1864      0.0      .000     .000     .000    .0     900     .0     7
   3    4    0.3811   0.1941      0.0      .000     .000     .000    .0     900     .0     7
   4    5    0.8190   0.7070      0.0      .000     .000     .000    .0     900     .0     7
   5    6    0.1872   0.6188      0.0      .000     .000     .000    .0     900     .0     7
   6    7    0.7114   0.2351      0.0      .000     .000     .000    .0     900     .0     7
   7    8    1.0300   0.7400      0.0      .000     .000     .000    .0     900     .0     7
   8    9    1.0440   0.7400      0.0      .000     .000     .000    .0     900     .0     7
   9    10   0.1966   0.0650      0.0      .000     .000     .000    .0     900     .0     7
   10   11   0.3744   0.1238      0.0      .000     .000     .000    .0     900     .0     7
   11   12   1.4680   1.1550      0.0      .000     .000     .000    .0     900     .0     7
   12   13   0.5416   0.7129      0.0      .000     .000     .000    .0     900     .0     7
   13   14   0.5910   0.5260      0.0      .000     .000     .000    .0     900     .0     7
   14   15   0.7463   0.5450      0.0      .000     .000     .000    .0     900     .0     7
   15   16   1.2890   1.7210      0.0      .000     .000     .000    .0     900     .0     7
   16   17   0.7320   0.5740      0.0      .000     .000     .000    .0     900     .0     7
    1   18   0.1640   0.1565      0.0      .000     .000     .000    .0     900     .0     7
   18   19   1.5042   1.3554      0.0      .000     .000     .000    .0     900     .0     7
   19   20   0.4095   0.4784      0.0      .000     .000     .000    .0     900     .0     7
   20   21   0.7089   0.9373      0.0      .000     .000     .000    .0     900     .0     7
   2    22   0.4512   0.3083      0.0      .000     .000     .000    .0     900     .0     7
   22   23   0.8980   0.7091      0.0      .000     .000     .000    .0     900     .0     7
   23   24   0.8960   0.7011      0.0      .000     .000     .000    .0     900     .0     7
   5    25   0.2030   0.1034      0.0      .000     .000     .000    .0     900     .0     7
   25   26   0.2842   0.1447      0.0      .000     .000     .000    .0     900     .0     7
   26   27   1.0590   0.9337      0.0      .000     .000     .000    .0     900     .0     7
   27   28   0.8042   0.7006      0.0      .000     .000     .000    .0     900     .0     7
   28   29   0.5075   0.2585      0.0      .000     .000     .000    .0     900     .0     7
   29   30   0.9744   0.9630      0.0      .000     .000     .000    .0     900     .0     7
   30   31   0.3105   0.3619      0.0      .000     .000     .000    .0     900     .0     7
   31   32   0.3410   0.5302      0.0      .000     .000     .000    .0     900     .0     7	];


