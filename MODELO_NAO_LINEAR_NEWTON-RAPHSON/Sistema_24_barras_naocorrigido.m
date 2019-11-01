% ****   Sistema teste de 24 BARRAS   ****                             

% ------------------------------ DADOS DE BARRA -----------------------------------

%  No TB G (V)  (Ang(graus))  (Pg)    (Qg)   (Qn)  (Qm)     (Pl)   (Ql)    bshbar  
 DBAR = [                                                                           
   1 2 1 1.040      0.0      172.0   34.8   -999.9  9999.9  108.0  22.0     0.0      
   2 2 1 1.040     -4.98     172.0   20.0   -999.9  9999.9  97.0   19.7     0.0      
   3 1 1 0.983     -12.72      .0      .0   -999.9  9999.9  180.0  36.5     0.0      
   4 1 1 1.005     -10.33      .0      .0   -999.9  9999.9  74.0   15.0     0.0      
   5 1 1 1.021     -8.78       .0      .0   -999.9  9999.9  71.0   14.5     0.0      
   6 1 1 1.013     -14.22      .0      .0   -999.9  9999.9  136.0  27.8     0.0      
   7 2 1 0.990     -13.37    240.0   19.0   -999.9  9999.9  125.0  25.5     0.0      
   8 1 1 1.090     -13.36      .0      .0   -999.9  9999.9  171.0  34.7     0.0      
   9 1 1 1.056     -14.94      .0      .0   -999.9  9999.9  175.0  35.3     0.0      
  10 1 1 1.051     -15.10      .0      .0   -999.9  9999.9  195.0  39.4     0.0      
  11 1 1 1.057     -14.79      .0      .0   -999.9  9999.9  0.0    0.00     0.0      
  12 1 1 1.055     -15.07      .0      .0   -999.9  9999.9  0.0    0.00     0.0      
  13 3 1 1.040      0.0      190.0     .0   -999.9  9999.9  265.0  0.00     0.0      
  14 2 1 0.995     -16.04      .0    15.0   -999.9  9999.9  194.0  39.4     0.0
  15 2 1 1.006     -13.36    215.0   13.10  -999.9  9999.9  317.0  64.2     0.0      
  16 2 1 1.010     -14.94    155.0   06.60  -999.9  9999.9  100.0  20.3     0.0      
  17 1 1 1.022     -15.10      .0      .0   -999.9  9999.9  0.0    0.00     0.0      
  18 2 1 1.025     -14.79    400.0   39.70  -999.9  9999.9  333.0  67.7     0.0      
  19 2 1 1.055     -15.07      .0      .0   -999.9  9999.9  181.0  37.0     0.0      
  20 1 1 1.050     -15.16      .0      .0   -999.9  9999.9  128.0  26.0     0.0      
  21 2 1 1.030      0.0      400.0   61.30  -999.9  9999.9  0.0    0.00     0.0      
  22 2 1 1.050     -15.07    300.0   18.50  -999.9  9999.9  0.0    0.00     0.0      
  23 2 1 1.050     -15.16    660.0   106.3  -999.9  9999.9  0.0    0.00     0.0      
  24 1 1 1.000     -16.04      .0      .0   -999.9  9999.9  0.0    0.00     0.0	];

% Tipos de barra: 1 - carga (PQ), 2 - geracao (PV), 3 - referencia (V-theta)
 
 PB = 100;

% ------------------------------ DADOS DE LINHA -----------------------------------

 DLIN = [
 %FROM  TO    %R(pu)   %X(pu)     %Bsh(pu)    %Tap      %PHI                                  CH
   1    2    0.00260   0.01390    0.46110      1.000     .000     .000    .0     900     .0     7
   1    3    0.05460   0.21120    0.05720      1.000     .000     .000    .0     900     .0     7
   1    5    0.02180   0.08450    0.02290      1.000     .000     .000    .0     900     .0     7
   2    4    0.03280   0.12670    0.03430      1.000     .000     .000    .0     900     .0     7
   2    6    0.04970   0.19200    0.05200      1.000     .000     .000    .0     900     .0     7
   3    9    0.03080   0.11900    0.03220      1.000     .000     .000    .0     900     .0     7
   3    24   0.00230   0.08390    0.00000      1.000     .000     .000    .0     900     .0     7
   4    9    0.02680   0.10370    0.02810      1.000     .000     .000    .0     900     .0     7
   5    10   0.02280   0.08830    0.02390      1.000     .000     .000    .0     900     .0     7
   6    10   0.01390   0.06050    2.45900      1.000     .000     .000    .0     900     .0     7
   7    8    0.01590   0.06140    0.01660      1.000     .000     .000    .0     900     .0     7
   8    9    0.04270   0.16510    0.04470      1.000     .000     .000    .0     900     .0     7
   8   10    0.04270   0.16510    0.04470      1.000     .000     .000    .0     900     .0     7
   9   11    0.00230   0.08390    0.00000      1.000     .000     .000    .0     900     .0     7
   9   12    0.00230   0.08390    0.00000      1.000     .000     .000    .0     900     .0     7
   10  11    0.00230   0.08390    0.00000      1.000     .000     .000    .0     900     .0     7
   10  12    0.00230   0.08390    0.00000      1.000     .000     .000    .0     900     .0     7
   11  13    0.00610   0.04760    0.09990      1.000     .000     .000    .0     900     .0     7
   11  14    0.00540   0.04180    0.08790      1.000     .000     .000    .0     900     .0     7
   12  13    0.00610   0.04760    0.09990      1.000     .000     .000    .0     900     .0     7
   12  23    0.01240   0.09660    0.20300      1.000     .000     .000    .0     900     .0     7
   13  23    0.01100   0.08650    0.18200      1.000     .000     .000    .0     900     .0     7
   14  16    0.00500   0.03890    0.08180      1.000     .000     .000    .0     900     .0     7
   15  16    0.00220   0.01730    0.03640      1.000     .000     .000    .0     900     .0     7
   15  21    0.00630   0.04900    0.10300      1.000     .000     .000    .0     900     .0     7
   15  21    0.00630   0.04900    0.10300      1.000     .000     .000    .0     900     .0     7
   15  24    0.00670   0.05190    0.10900      1.000     .000     .000    .0     900     .0     7
   16  17    0.00330   0.02590    0.05450      1.000     .000     .000    .0     900     .0     7
   16  19    0.00300   0.02310    0.04850      1.000     .000     .000    .0     900     .0     7
   17  18    0.00180   0.01400    0.03030      1.000     .000     .000    .0     900     .0     7
   17  22    0.01350   0.10530    0.22100      1.000     .000     .000    .0     900     .0     7
   18  21    0.00330   0.02590    0.05450      1.000     .000     .000    .0     900     .0     7
   18  21    0.00330   0.02590    0.05450      1.000     .000     .000    .0     900     .0     7
   19  20    0.00510   0.03960    0.08300      1.000     .000     .000    .0     900     .0     7
   19  20    0.00510   0.03960    0.08300      1.000     .000     .000    .0     900     .0     7
   20  23    0.00280   0.02160    0.04550      1.000     .000     .000    .0     900     .0     7
   20  23    0.00280   0.02160    0.04550      1.000     .000     .000    .0     900     .0     7
   21  22    0.00870   0.06780    0.14240      1.000     .000     .000    .0     900     .0     7	];