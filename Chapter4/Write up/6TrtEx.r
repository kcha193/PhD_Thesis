> design.summary.RBD(design.df)
Design parameters:
Trt: 6 Ani: 24 Cag: 4 Tag: 4 Run: 12 
Animal design:
      [,1] [,2] [,3] [,4]
 [1,] "1B" "2A" "2C" "1D"
 [2,] "2A" "1B" "1D" "2C"
 [3,] "2D" "1C" "2F" "1E"
 [4,] "1C" "2D" "1E" "2F"
 [5,] "1F" "2E" "1A" "2B"
 [6,] "2E" "1F" "2B" "1A"
 [7,] "3C" "4B" "3A" "4F"
 [8,] "4B" "3C" "4F" "3A"
 [9,] "3D" "4A" "3E" "4C"
[10,] "4A" "3D" "4C" "3E"
[11,] "3F" "4E" "3B" "4D"
[12,] "4E" "3F" "4D" "3B"
Animal efficiency:
$nCan
[1] 16

$can.eff
 [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

Treatment design:
      [,1] [,2] [,3] [,4]
 [1,] "b"  "a"  "c"  "d" 
 [2,] "a"  "b"  "d"  "c" 
 [3,] "d"  "c"  "f"  "e" 
 [4,] "c"  "d"  "e"  "f" 
 [5,] "f"  "e"  "a"  "b" 
 [6,] "e"  "f"  "b"  "a" 
 [7,] "c"  "b"  "a"  "f" 
 [8,] "b"  "c"  "f"  "a" 
 [9,] "d"  "a"  "e"  "c" 
[10,] "a"  "d"  "c"  "e" 
[11,] "f"  "e"  "b"  "d" 
[12,] "e"  "f"  "d"  "b" 
Treatment incidence matrix:
   Run
Trt 1 2 3 4 5 6 7 8 9 10 11 12
  a 1 1 0 0 1 1 1 1 1  1  0  0
  b 1 1 0 0 1 1 1 1 0  0  1  1
  c 1 1 1 1 0 0 1 1 1  1  0  0
  d 1 1 1 1 0 0 0 0 1  1  1  1
  e 0 0 1 1 1 1 0 0 1  1  1  1
  f 0 0 1 1 1 1 1 1 0  0  1  1
Treatment concurrence matrix:
   Trt
Trt a b c d e f
  a 8 6 6 4 4 4
  b 6 8 4 4 4 6
  c 6 4 8 6 4 4
  d 4 4 6 8 6 4
  e 4 4 4 6 8 6
  f 4 6 4 4 6 8
Treatment efficiency:
$trace
[1] 36

$nCan
[1] 5

$can.eff
[1] 1.0000 0.9375 0.9375 0.8125 0.8125

$ave.eff
[1] 0.8936755

$e.vec
           [,1]       [,2]        [,3]       [,4]       [,5]       [,6]
[1,]  0.4082483  0.5748886 -0.05325797 -0.5678682  0.1042065 -0.4082483
[2,] -0.4082483 -0.3335671 -0.47123916 -0.1936886  0.5438916 -0.4082483
[3,] -0.4082483 -0.2413216  0.52449713 -0.3741796 -0.4396850 -0.4082483
[4,]  0.4082483 -0.3335671 -0.47123916  0.1936886 -0.5438916 -0.4082483
[5,] -0.4082483  0.5748886 -0.05325797  0.5678682 -0.1042065 -0.4082483
[6,]  0.4082483 -0.2413216  0.52449713  0.3741796  0.4396850 -0.4082483

Phase 1 theoretical ANOVA:
$ANOVA
                DF e Cag:Ani Cag
Between Cag     3  1 2       12 
Between Cag:Ani                 
   Trt          5  1 2       0  
   Residual     15 1 2       0  
Within          24 1 0       0  

$EF
                Trt eff.Trt
Between Cag                
Between Cag:Ani            
   Trt          8   1      
Within                     

Phase 2 theoretical ANOVA:
$ANOVA
                   DF e Cag:Ani Cag Run
Between Run                            
   Between Cag     1  1 2       12  4  
   Between Cag:Ani                     
      Trt          4  1 2       0   4  
   Residual        6  1 0       0   4  
Within                                 
   Between Cag     2  1 2       12  0  
   Between Cag:Ani                     
      Tag          1  1 2       0   0  
      Trt          5  1 2       0   0  
      Residual     10 1 2       0   0  
   Residual                            
      Tag          2  1 0       0   0  
      Residual     16 1 0       0   0  

$EF
                   Tag Trt       eff.Tag eff.Trt 
Between Run                                      
   Between Cag                                   
   Between Cag:Ani                               
      Trt              3/4               3/32    
    Residual                                     
Within                                           
   Between Cag                                   
   Between Cag:Ani                               
      Tag          12            1               
      Trt              7800/1091         975/1091
    Residual                                     
      Tag          12            1               


Phase 1 theoretical ANOVA:
\begin{table}[ht]
\centering
 \caption{Theoretical ANOVA table}
 \begin{tabular}[t]{lrll} 
 \toprule 
 \multicolumn{1}{l}{\textbf{Source of Variation}} & \multicolumn{1}{l}{\textbf{DF}} & \multicolumn{1}{l}{\textbf{EMS}}& \multicolumn{1}{l}{$\bm{E_{\gamma}}$}\\ 
 \midrule 
 Between Cag & $3$ & $\sigma^2+2\sigma_{A}^2+12\sigma_{C}^2$ &\\ 
 Between Cag:Ani &  &  &\\ 
 \quad Trt & $5$ & $\sigma^2+2\sigma_{A}^2+8\theta_{\gamma}$ &$1$\\ 
 \quad Residual & $15$ & $\sigma^2+2\sigma_{A}^2$ &\\ 
 Within & $24$ & $\sigma^2$ &\\ 
 \bottomrule 
 \end{tabular} 
 \label{tab:} 
\end{table} 
 NULL
Phase 2 theoretical ANOVA:
\begin{table}[ht]
\centering
 \caption{Theoretical ANOVA table}
 \begin{tabular}[t]{lrlll} 
 \toprule 
 \multicolumn{1}{l}{\textbf{Source of Variation}} & \multicolumn{1}{l}{\textbf{DF}} & \multicolumn{1}{l}{\textbf{EMS}}& \multicolumn{1}{l}{$\bm{E_{\gamma}}$}&\multicolumn{1}{l}{$\bm{E_{\tau}}$}\\ 
 \midrule 
 Between Run &  &  & & \\ 
 \quad Between Cag & $1$ & $\sigma^2+2\sigma_{A}^2+12\sigma_{C}^2+4\sigma_{R}^2$ & & \\ 
 \quad Between Cag:Ani &  &  & & \\ 
 \quad \quad Trt & $4$ & $\sigma^2+2\sigma_{A}^2+4\sigma_{R}^2+3/4\theta_{\tau}$ & & $3/32$\\ 
 \quad Residual & $6$ & $\sigma^2+4\sigma_{R}^2$ & & \\ 
 Within &  &  & & \\ 
 \quad Between Cag & $2$ & $\sigma^2+2\sigma_{A}^2+12\sigma_{C}^2$ & & \\ 
 \quad Between Cag:Ani &  &  & & \\ 
 \quad \quad Tag & $1$ & $\sigma^2+2\sigma_{A}^2+12\theta_{\gamma}$ &$1$ & \\ 
 \quad \quad Trt & $5$ & $\sigma^2+2\sigma_{A}^2+7800/1091\theta_{\tau}$ & & $975/1091$\\ 
 \quad \quad Residual & $10$ & $\sigma^2+2\sigma_{A}^2$ & & \\ 
 \quad Residual &  &  & & \\ 
 \quad \quad Tag & $2$ & $\sigma^2+12\theta_{\gamma}$ &$1$ & \\ 
 \quad \quad Residual & $16$ & $\sigma^2$ & & \\ 
 \bottomrule 
 \end{tabular} 
 \label{tab:} 
\end{table} 
 
Design parameters:
Trt: 6 Ani: 24 Cag: 4 Tag: 4 Run: 12 
Animal design:
      [,1] [,2] [,3] [,4]
 [1,] "1A" "1D" "3E" "3F"
 [2,] "1D" "1A" "3F" "3E"
 [3,] "1F" "1E" "3B" "3C"
 [4,] "1E" "1F" "3C" "3B"
 [5,] "1B" "1C" "3D" "3A"
 [6,] "1C" "1B" "3A" "3D"
 [7,] "2D" "2E" "4F" "4B"
 [8,] "2E" "2D" "4B" "4F"
 [9,] "2A" "2B" "4C" "4E"
[10,] "2B" "2A" "4E" "4C"
[11,] "2C" "2F" "4A" "4D"
[12,] "2F" "2C" "4D" "4A"
Animal efficiency:
$nCan
[1] 16

$can.eff
 [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

Treatment design:
      [,1] [,2] [,3] [,4]
 [1,] "a"  "d"  "e"  "f" 
 [2,] "d"  "a"  "f"  "e" 
 [3,] "f"  "e"  "b"  "c" 
 [4,] "e"  "f"  "c"  "b" 
 [5,] "b"  "c"  "d"  "a" 
 [6,] "c"  "b"  "a"  "d" 
 [7,] "d"  "e"  "f"  "b" 
 [8,] "e"  "d"  "b"  "f" 
 [9,] "a"  "b"  "c"  "e" 
[10,] "b"  "a"  "e"  "c" 
[11,] "c"  "f"  "a"  "d" 
[12,] "f"  "c"  "d"  "a" 
Treatment incidence matrix:
   Run
Trt 1 2 3 4 5 6 7 8 9 10 11 12
  a 1 1 0 0 1 1 0 0 1  1  1  1
  b 0 0 1 1 1 1 1 1 1  1  0  0
  c 0 0 1 1 1 1 0 0 1  1  1  1
  d 1 1 0 0 1 1 1 1 0  0  1  1
  e 1 1 1 1 0 0 1 1 1  1  0  0
  f 1 1 1 1 0 0 1 1 0  0  1  1
Treatment concurrence matrix:
   Trt
Trt a b c d e f
  a 8 4 6 6 4 4
  b 4 8 6 4 6 4
  c 6 6 8 4 4 4
  d 6 4 4 8 4 6
  e 4 6 4 4 8 6
  f 4 4 4 6 6 8
Treatment efficiency:
$trace
[1] 36

$nCan
[1] 5

$can.eff
[1] 1.0000 0.9375 0.9375 0.8125 0.8125

$ave.eff
[1] 0.8936755

$e.vec
           [,1]       [,2]       [,3]        [,4]       [,5]       [,6]
[1,]  0.4082483 -0.1321076 -0.5620328 -0.45717138  0.3526013 -0.4082483
[2,]  0.4082483 -0.4206809  0.3954250 -0.07677598 -0.5722227 -0.4082483
[3,] -0.4082483  0.5527885  0.1666078 -0.53394736 -0.2196214 -0.4082483
[4,] -0.4082483 -0.4206809  0.3954250  0.07677598  0.5722227 -0.4082483
[5,] -0.4082483 -0.1321076 -0.5620328  0.45717138 -0.3526013 -0.4082483
[6,]  0.4082483  0.5527885  0.1666078  0.53394736  0.2196214 -0.4082483

Phase 1 theoretical ANOVA:
$ANOVA
                DF e Cag:Ani Cag
Between Cag     3  1 2       12 
Between Cag:Ani                 
   Trt          5  1 2       0  
   Residual     15 1 2       0  
Within          24 1 0       0  

$EF
                Trt eff.Trt
Between Cag                
Between Cag:Ani            
   Trt          8   1      
Within                     

Phase 2 theoretical ANOVA:
$ANOVA
                   DF e Cag:Ani Cag Run
Between Run                            
   Between Cag     1  1 2       12  4  
   Between Cag:Ani                     
      Trt          4  1 2       0   4  
   Residual        6  1 0       0   4  
Within                                 
   Between Cag                         
      Tag          1  1 2       12  0  
      Residual     1  1 2       12  0  
   Between Cag:Ani                     
      Trt          5  1 2       0   0  
      Residual     11 1 2       0   0  
   Residual                            
      Tag          2  1 0       0   0  
      Residual     16 1 0       0   0  

$EF
                   Tag Trt       eff.Tag eff.Trt 
Between Run                                      
   Between Cag                                   
   Between Cag:Ani                               
      Trt              3/4               3/32    
    Residual                                     
Within                                           
   Between Cag                                   
      Tag          12            1               
   Between Cag:Ani                               
      Trt              7800/1091         975/1091
    Residual                                     
      Tag          12            1         
      
      
\begin{table}[ht]
\centering
 \caption{Theoretical ANOVA table}
 \begin{tabular}[t]{lrlll} 
 \toprule 
 \multicolumn{1}{l}{\textbf{Source of Variation}} & \multicolumn{1}{l}{\textbf{DF}} & \multicolumn{1}{l}{\textbf{EMS}}& \multicolumn{1}{l}{$\bm{E_{\gamma}}$}&\multicolumn{1}{l}{$\bm{E_{\tau}}$}\\ 
 \midrule 
 Between Run &  &  & & \\ 
 \quad Between Cag & $1$ & $\sigma^2+2\sigma_{A}^2+12\sigma_{C}^2+4\sigma_{R}^2$ & & \\ 
 \quad Between Cag:Ani &  &  & & \\ 
 \quad \quad Trt & $4$ & $\sigma^2+2\sigma_{A}^2+4\sigma_{R}^2+3/4\theta_{\tau}$ & & $3/32$\\ 
 \quad Residual & $6$ & $\sigma^2+4\sigma_{R}^2$ & & \\ 
 Within &  &  & & \\ 
 \quad Between Cag &  &  & & \\ 
 \quad \quad Tag & $1$ & $\sigma^2+2\sigma_{A}^2+12\sigma_{C}^2+12\theta_{\gamma}$ &$1$ & \\ 
 \quad \quad Residual & $1$ & $\sigma^2+2\sigma_{A}^2+12\sigma_{C}^2$ & & \\ 
 \quad Between Cag:Ani &  &  & & \\ 
 \quad \quad Trt & $5$ & $\sigma^2+2\sigma_{A}^2+7800/1091\theta_{\tau}$ & & $975/1091$\\ 
 \quad \quad Residual & $11$ & $\sigma^2+2\sigma_{A}^2$ & & \\ 
 \quad Residual &  &  & & \\ 
 \quad \quad Tag & $2$ & $\sigma^2+12\theta_{\gamma}$ &$1$ & \\ 
 \quad \quad Residual & $16$ & $\sigma^2$ & & \\ 
 \bottomrule 
 \end{tabular} 
 \label{tab:} 
\end{table} 
      