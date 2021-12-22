#Special case
1  & a & c & b & d \\  
2  & c & a & d & b \\  
3  & b & d & a & c \\  
4  & d & b & c & a \\  
5  & a & c & d & b \\  
6  & c & a & b & d \\  
7  & b & d & a & c \\  
8  & d & b & c & a \\             



> design.summary.RBD(design.df)
Design parameters:
Trt: 4 Ani: 16 Cag: 4 Tag: 4 Run: 8
Animal design:
     [,1] [,2] [,3] [,4]
[1,] "2A" "2C" "3B" "3D"
[2,] "2C" "2A" "3D" "3B"
[3,] "1B" "1D" "4A" "4C"
[4,] "1D" "1B" "4C" "4A"
[5,] "1A" "1C" "4D" "4B"
[6,] "1C" "1A" "4B" "4D"
[7,] "2B" "2D" "3A" "3C"
[8,] "2D" "2B" "3C" "3A"
Animal efficiency:
$nCan
[1] 10

$can.eff
 [1] 1 1 1 1 1 1 1 1 1 1

Treatment design:
     [,1] [,2] [,3] [,4]
[1,] "a"  "c"  "b"  "d"
[2,] "c"  "a"  "d"  "b"
[3,] "b"  "d"  "a"  "c"
[4,] "d"  "b"  "c"  "a"
[5,] "a"  "c"  "d"  "b"
[6,] "c"  "a"  "b"  "d"
[7,] "b"  "d"  "a"  "c"
[8,] "d"  "b"  "c"  "a"
Treatment incidence matrix:
   Run
Trt 1 2 3 4 5 6 7 8
  a 1 1 1 1 1 1 1 1
  b 1 1 1 1 1 1 1 1
  c 1 1 1 1 1 1 1 1
  d 1 1 1 1 1 1 1 1
Treatment concurrence matrix:
   Trt
Trt a b c d
  a 8 8 8 8
  b 8 8 8 8
  c 8 8 8 8
  d 8 8 8 8
Treatment efficiency:
$trace
[1] 24

$nCan
[1] 3

$can.eff
[1] 1 1 1

$ave.eff
[1] 1

$e.vec
            [,1]       [,2]       [,3] [,4]
[1,] -0.55043070  0.0000000  0.6686001  0.5
[2,] -0.26225647 -0.5773503 -0.5898205  0.5
[3,]  0.02032716  0.7886751 -0.3571811  0.5
[4,]  0.79236001 -0.2113249  0.2784015  0.5

Phase 1 theoretical ANOVA:
$ANOVA
                DF e Cag:Ani Cag
Between Cag     3  1 2       8
Between Cag:Ani
   Trt          3  1 2       0
   Residual     9  1 2       0
Within          16 1 0       0

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
   Between Cag     1  1 2       8   4
   Between Cag:Ani 2  1 2       0   4
   Residual        4  1 0       0   4
Within
   Between Cag
      Tag          1  1 2       8   0
      Residual     1  1 2       8   0
   Between Cag:Ani
      Trt          3  1 2       0   0
      Residual     7  1 2       0   0
   Residual
      Tag          2  1 0       0   0
      Residual     10 1 0       0   0

$EF
                   Tag Trt eff.Tag eff.Trt
Between Run
   Between Cag
   Between Cag:Ani
    Residual
Within
   Between Cag
      Tag          8       1
   Between Cag:Ani
      Trt              8           1
    Residual
      Tag          8       1

> design.summary.RBD(des)
Design parameters:
Trt: 4 Ani: 16 Cag: 4 Tag: 4 Run: 8
Animal design:
     [,1] [,2] [,3] [,4]
[1,] "1D" "1B" "1A" "1C"
[2,] "1B" "1D" "1C" "1A"
[3,] "2A" "2D" "2C" "2B"
[4,] "2D" "2A" "2B" "2C"
[5,] "3B" "3C" "3D" "3A"
[6,] "3C" "3B" "3A" "3D"
[7,] "4C" "4A" "4D" "4B"
[8,] "4A" "4C" "4B" "4D"

1  & d & b & a & c \\  
2  & b & d & c & a \\  
3  & a & d & c & b \\  
4  & d & a & b & c \\  
5  & b & c & d & a \\  
6  & c & b & a & d \\  
7  & c & a & d & b \\  
8  & a & c & b & d \\             


Animal efficiency:
$nCan
[1] 12

$can.eff
 [1] 1 1 1 1 1 1 1 1 1 1 1 1

Treatment design:
     [,1] [,2] [,3] [,4]
[1,] "d"  "b"  "a"  "c"
[2,] "b"  "d"  "c"  "a"
[3,] "a"  "d"  "c"  "b"
[4,] "d"  "a"  "b"  "c"
[5,] "b"  "c"  "d"  "a"
[6,] "c"  "b"  "a"  "d"
[7,] "c"  "a"  "d"  "b"
[8,] "a"  "c"  "b"  "d"
Treatment incidence matrix:
   Run
Trt 1 2 3 4 5 6 7 8
  a 1 1 1 1 1 1 1 1
  b 1 1 1 1 1 1 1 1
  c 1 1 1 1 1 1 1 1
  d 1 1 1 1 1 1 1 1
Treatment concurrence matrix:
   Trt
Trt a b c d
  a 8 8 8 8
  b 8 8 8 8
  c 8 8 8 8
  d 8 8 8 8
Treatment efficiency:
$trace
[1] 24

$nCan
[1] 3

$can.eff
[1] 1 1 1

$ave.eff
[1] 1

$e.vec
           [,1]        [,2]       [,3] [,4]
[1,] -0.4538160  0.71804685  0.1687003 -0.5
[2,]  0.7628165  0.03802978  0.4082458 -0.5
[3,] -0.4412243 -0.69199415  0.2765235 -0.5
[4,]  0.1322237 -0.06408248 -0.8534696 -0.5

Phase 1 theoretical ANOVA:
$ANOVA
                DF e Cag:Ani Cag
Between Cag     3  1 2       8
Between Cag:Ani
   Trt          3  1 2       0
   Residual     9  1 2       0
Within          16 1 0       0

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
   Between Cag     3  1 2       8   4
   Residual        4  1 0       0   4
Within
   Between Cag:Ani
      Tag          1  1 2       0   0
      Trt          3  1 2       0   0
      Residual     8  1 2       0   0
   Residual
      Tag          2  1 0       0   0
      Residual     10 1 0       0   0

$EF
                   Tag Trt eff.Tag eff.Trt
Between Run
   Between Cag
    Residual
Within
   Between Cag:Ani
      Tag          8       1
      Trt              8           1
    Residual
      Tag          8       1

Phase 1 theoretical ANOVA:
\begin{table}[ht]
\centering
 \caption{Theoretical ANOVA table}
 \begin{tabular}[t]{lrll}
 \toprule
 \multicolumn{1}{l}{\textbf{Source of Variation}} & \multicolumn{1}{l}{\textbf{DF}} & \multicolumn{1}{l}{\textbf{EMS}}& \multicolumn{1}{l}{$\bm{E_{\gamma}}$}\\
 \midrule
 Between Cag & $3$ & $\sigma^2+2\sigma_{A}^2+8\sigma_{C}^2$ &\\
 Between Cag:Ani &  &  &\\
 \quad Trt & $3$ & $\sigma^2+2\sigma_{A}^2+8\theta_{\gamma}$ &$1$\\
 \quad Residual & $9$ & $\sigma^2+2\sigma_{A}^2$ &\\
 Within & $16$ & $\sigma^2$ &\\
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
 \quad Between Cag & $3$ & $\sigma^2+2\sigma_{A}^2+8\sigma_{C}^2+4\sigma_{R}^2$ & & \\
 \quad Residual & $4$ & $\sigma^2+4\sigma_{R}^2$ & & \\
 Within &  &  & & \\
 \quad Between Cag:Ani &  &  & & \\
 \quad \quad Tag & $1$ & $\sigma^2+2\sigma_{A}^2+8\theta_{\gamma}$ &$1$ & \\
 \quad \quad Trt & $3$ & $\sigma^2+2\sigma_{A}^2+8\theta_{\tau}$ & & $1$\\
 \quad \quad Residual & $8$ & $\sigma^2+2\sigma_{A}^2$ & & \\
 \quad Residual &  &  & & \\
 \quad \quad Tag & $2$ & $\sigma^2+8\theta_{\gamma}$ &$1$ & \\
 \quad \quad Residual & $10$ & $\sigma^2$ & & \\
 \bottomrule
 \end{tabular}
 \label{tab:}
\end{table}


\begin{table}[ht]
\centering
 \caption{Theoretical ANOVA table}
 \begin{tabular}[t]{lrlll}
 \toprule
 \multicolumn{1}{l}{\textbf{Source of Variation}} & \multicolumn{1}{l}{\textbf{DF}} & \multicolumn{1}{l}{\textbf{EMS}}& \multicolumn{1}{l}{$\bm{E_{\gamma}}$}&\multicolumn{1}{l}{$\bm{E_{\tau}}$}\\
 \midrule
 Between Run &  &  & & \\
 \quad Between Cag & $1$ & $\sigma^2+2\sigma_{A}^2+8\sigma_{C}^2+4\sigma_{R}^2$ & & \\
 \quad Between Cag:Ani & $2$ & $\sigma^2+2\sigma_{A}^2+4\sigma_{R}^2$ & & \\
 \quad Residual & $4$ & $\sigma^2+4\sigma_{R}^2$ & & \\
 Within &  &  & & \\
 \quad Between Cag &  &  & & \\
 \quad \quad Tag & $1$ & $\sigma^2+2\sigma_{A}^2+8\sigma_{C}^2+8\theta_{\gamma}$ &$1$ & \\
 \quad \quad Residual & $1$ & $\sigma^2+2\sigma_{A}^2+8\sigma_{C}^2$ & & \\
 \quad Between Cag:Ani &  &  & & \\
 \quad \quad Trt & $3$ & $\sigma^2+2\sigma_{A}^2+8\theta_{\tau}$ & & $1$\\
 \quad \quad Residual & $7$ & $\sigma^2+2\sigma_{A}^2$ & & \\
 \quad Residual &  &  & & \\
 \quad \quad Tag & $2$ & $\sigma^2+8\theta_{\gamma}$ &$1$ & \\
 \quad \quad Residual & $10$ & $\sigma^2$ & & \\
 \bottomrule
 \end{tabular}
 \label{tab:}
\end{table}
