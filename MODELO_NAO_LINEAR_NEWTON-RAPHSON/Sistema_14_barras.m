% ****   Sistema teste de 14 BARRAS   ****                             

% ------------------------------ DADOS DE BARRA -----------------------------------

%  No TB G (V)     (Ang)     (Pg)    (Qg)   (Qn)  (Qm)      (Pl)    (Ql)    bshbar  
 DBAR = [                                                                           
   1 3 1 1.060      0.0      232.4  -16.9  -999.9  9999.9    0.0    0.0      0.0      
   2 2 1 1.045     -0.0855   40.0    42.4  -999.9  9999.9    21.7   12.7     0.0      
   3 2 1 1.010     -0.2094   .0      23.4  -999.9  9999.9    94.2   19.0     0.0      
   4 1 1 1.019     -0.1745   .0      .0    -999.9  9999.9    47.8   -3.9     0.0      
   5 1 1 1.020     -0.1518   .0      .0    -999.9  9999.9    7.6    1.6      0.0      
   6 2 1 1.070     -0.2443   .0      12.2  -999.9  9999.9    11.2   7.5      0.0      
   7 1 1 1.062     -0.2269   .0      .0    -999.9  9999.9    0.0    0.0      0.0      
   8 2 1 1.090     -0.2269   .0      17.4  -999.9  9999.9    0.0    0.0      0.0      
   9 1 1 1.056     -0.2443   .0      .0    -999.9  9999.9    29.5   16.6     0.19      
  10 1 1 1.051     -0.2618   .0      .0    -999.9  9999.9    9.0    5.8      0.0      
  11 1 1 1.057     -0.2443   .0      .0    -999.9  9999.9    3.5    1.8      0.0      
  12 1 1 1.055     -0.2618   .0      .0    -999.9  9999.9    6.1    1.6      0.0      
  13 1 1 1.050     -0.2618   .0      .0    -999.9  9999.9    13.5   5.8      0.0      
  14 1 1 1.036     -0.2793   .0      .0    -999.9  9999.9    14.9   5.0      0.0      ];

% Tipos de barra: 1 - carga (PQ), 2 - geracao (PV), 3 - referencia (V-theta)
 
PB = 100;


% ------------------------------ DADOS DE LINHA -----------------------------------

 DLIN = [
 %FROM  TO    %R(pu)   %X(pu)     %Bsh(pu)    %Tap      %PHI                                  CH
   1    2    0.01938  0.05917     0.0264      1.000     .000     .000    .0     900     .0     7
   1    5    0.05403  0.22304     0.0246      1.000     .000     .000    .0     900     .0     7
   2    3    0.04699  0.19797     0.0219      1.000     .000     .000    .0     900     .0     7
   2    4    0.05811  0.17632     0.0170      1.000     .000     .000    .0     900     .0     7
   2    5    0.05695  0.17388     0.0173      1.000     .000     .000    .0     900     .0     7
   3    4    0.06701  0.17103     0.0064      1.000     .000     .000    .0     900     .0     7
   4    5    0.01335  0.04211     0.0000      1.000     .000     .000    .0     900     .0     7
   4    7    0.00000  0.20912     0.0000      1.022     .000     .000    .0     900     .0     7
   4    9    0.00000  0.55618     0.0000      1.032     .000     .000    .0     900     .0     7
   5    6    0.00000  0.25202     0.0000      1.073     .000     .000    .0     900     .0     7
   6   11    0.09498  0.19890     0.0000      1.000     .000     .000    .0     900     .0     7
   6   12    0.12291  0.25581     0.0000      1.000     .000     .000    .0     900     .0     7
   6   13    0.06615  0.13027     0.0000      1.000     .000     .000    .0     900     .0     7
   7    8    0.00000  0.17615     0.0000      1.000     .000     .000    .0     900     .0     7
   7    9    0.00000  0.11001     0.0000      1.000     .000     .000    .0     900     .0     7
   9   10    0.03181  0.08450     0.0000      1.000     .000     .000    .0     900     .0     7
   9   14    0.12711  0.27038     0.0000      1.000     .000     .000    .0     900     .0     7
   10  11    0.08205  0.19207     0.0000      1.000     .000     .000    .0     900     .0     7
   12  13    0.22092  0.19988     0.0000      1.000     .000     .000    .0     900     .0     7
   13  14    0.17093  0.34802     0.0000      1.000     .000     .000    .0     900     .0     7	];

