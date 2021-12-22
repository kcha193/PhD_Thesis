

Design parameters:
Trt: 4 Ani: 8 Cag: 2 Tag: 4 Run: 4
Animal design:
     [,1] [,2] [,3] [,4]
[1,] "1B" "1A" "1D" "1C"
[2,] "1A" "1B" "1C" "1D"
[3,] "2C" "2D" "2A" "2B"
[4,] "2D" "2C" "2B" "2A"
Animal efficiency:
$nCan
[1] 6

$can.eff
[1] 1 1 1 1 1 1

Treatment design:
     [,1] [,2] [,3] [,4]
[1,] "b"  "a"  "d"  "c"
[2,] "a"  "b"  "c"  "d"
[3,] "c"  "d"  "a"  "b"
[4,] "d"  "c"  "b"  "a"
Treatment incidence matrix:
   Run
Trt 1 2 3 4
  a 1 1 1 1
  b 1 1 1 1
  c 1 1 1 1
  d 1 1 1 1
Treatment concurrence matrix:
   Trt
Trt a b c d
  a 4 4 4 4
  b 4 4 4 4
  c 4 4 4 4
  d 4 4 4 4
Treatment efficiency:
$trace
[1] 12

$nCan
[1] 3

$can.eff
[1] 1 1 1

$ave.eff
[1] 1

$e.vec
           [,1]       [,2]        [,3] [,4]
[1,] -0.1638689  0.2747319  0.80477908 -0.5
[2,]  0.2026489 -0.8417897  0.01798546 -0.5
[3,] -0.7017859  0.1174941 -0.49365134 -0.5
[4,]  0.6630060  0.4495638 -0.32911321 -0.5

Phase 1 theoretical ANOVA:
$ANOVA
                DF e Cag:Ani Cag
Between Cag     1  1 2       8
Between Cag:Ani
   Trt          3  1 2       0
   Residual     3  1 2       0
Within Cag.Ani  8  1 0       0

$EF
                Trt eff.Trt
Between Cag
Between Cag:Ani
   Trt          4   1
Within Cag.Ani

Phase 2 theoretical ANOVA:
$ANOVA
                   DF e Cag:Ani Cag Run
Between Run
   Between Cag     1  1 2       8   4
   Within Cag.Ani  2  1 0       0   4
Within Run
   Between Cag:Ani
      Tag          1  1 2       0   0
      Trt          3  1 2       0   0
      Residual     2  1 2       0   0
   Within Cag.Ani
      Tag          2  1 0       0   0
      Residual     4  1 0       0   0

$EF
                   Tag Trt eff.Tag eff.Trt
Between Run
   Between Cag
   Within Cag.Ani
Within Run
   Between Cag:Ani
      Tag          4       1
      Trt              4           1
   Within Cag.Ani
      Tag          4       1



 NULL
Phase 2 theoretical ANOVA:
\begin{table}[ht]
\centering
 \caption{Theoretical ANOVA table}
 \begin{tabular}[t]{lrlll}
 \toprule
 \multicolumn{1}{l}{\textbf{Source of Variation}} & \multicolumn{1}{l}{\textbf{DF}} & \multicolumn{1}{l}{\textbf{EMS}}& \multicolumn{1}{l}{$\bm{E_{\gamma}}$}&\multicolumn{1}{l}{$\bm{E_{\tau}}$}\\
 \midrule
 Between Run &  &  & & \\ \hline
 \quad Between Cag & $1$ & $\sigma^2+2\sigma_{CA}^2+8\sigma_{C}^2+4\sigma_{R}^2$ & & \\ \hline
 \quad Within Cag.Ani & $2$ & $\sigma^2+4\sigma_{R}^2$ & & \\ \hline
 Within Run &  &  & & \\ \hline
 \quad Between Cag:Ani &  &  & & \\
 \quad \quad Tag & $1$ & $\sigma^2+2\sigma_{CA}^2+4\theta_{\gamma}$ &$1$ & \\
 \quad \quad Trt & $3$ & $\sigma^2+2\sigma_{CA}^2+4\theta_{\tau}$ & & $1$\\
 \quad \quad Residual & $2$ & $\sigma^2+2\sigma_{CA}^2$ & & \\ \hline
 \quad Within Cag.Ani &  &  & & \\
 \quad \quad Tag & $2$ & $\sigma^2+4\theta_{\gamma}$ &$1$ & \\
 \quad \quad Residual & $4$ & $\sigma^2$ & & \\
 \bottomrule
 \end{tabular}
 \label{tab:}
\end{table}

1  & a & b & d & d \\
2  & b & a & c & c \\
3  & c & d & b & b \\
4  & d & c & a & a \\




Design parameters:
Trt: 4 Ani: 8 Cag: 2 Tag: 4 Run: 4
Animal design:
     [,1] [,2] [,3] [,4]
[1,] "1A" "1B" "2D" "2C"
[2,] "1B" "1A" "2C" "2D"
[3,] "1C" "1D" "2B" "2A"
[4,] "1D" "1C" "2A" "2B"
Animal efficiency:
$nCan
[1] 5

$can.eff
[1] 1 1 1 1 1

Treatment design:
     [,1] [,2] [,3] [,4]
[1,] "a"  "b"  "d"  "c"
[2,] "b"  "a"  "c"  "d"
[3,] "c"  "d"  "b"  "a"
[4,] "d"  "c"  "a"  "b"
Treatment incidence matrix:
   Run
Trt 1 2 3 4
  a 1 1 1 1
  b 1 1 1 1
  c 1 1 1 1
  d 1 1 1 1
Treatment concurrence matrix:
   Trt
Trt a b c d
  a 4 4 4 4
  b 4 4 4 4
  c 4 4 4 4
  d 4 4 4 4
Treatment efficiency:
$trace
[1] 12

$nCan
[1] 3

$can.eff
[1] 1 1 1

$ave.eff
[1] 1

$e.vec
           [,1]        [,2]       [,3] [,4]
[1,] -0.2076070  0.40285260  0.7379764 -0.5
[2,]  0.2449835 -0.81231211  0.1735860 -0.5
[3,] -0.6880850 -0.01209997 -0.5257306 -0.5
[4,]  0.6507084  0.42155949 -0.3858317 -0.5

Phase 1 theoretical ANOVA:
$ANOVA
                DF e Cag:Ani Cag
Between Cag     1  1 2       8
Between Cag:Ani
   Trt          3  1 2       0
   Residual     3  1 2       0
Within Cag.Ani  8  1 0       0

$EF
                Trt eff.Trt
Between Cag
Between Cag:Ani
   Trt          4   1
Within Cag.Ani

Phase 2 theoretical ANOVA:
$ANOVA
                   DF e Cag:Ani Cag Run
Between Run
   Between Cag:Ani 1  1 2       0   4
   Within Cag.Ani  2  1 0       0   4
Within Run
   Between Cag
      Tag          1  1 2       8   0
   Between Cag:Ani
      Trt          3  1 2       0   0
      Residual     2  1 2       0   0
   Within Cag.Ani
      Tag          2  1 0       0   0
      Residual     4  1 0       0   0

$EF
                   Tag Trt eff.Tag eff.Trt
Between Run
   Between Cag:Ani
   Within Cag.Ani
Within Run
   Between Cag
      Tag          4       1
   Between Cag:Ani
      Trt              4           1
   Within Cag.Ani
      Tag          4       1


\begin{table}[ht]
\centering
 \caption{Theoretical ANOVA table}
 \begin{tabular}[t]{lrll}
 \toprule
 \multicolumn{1}{l}{\textbf{Source of Variation}} & \multicolumn{1}{l}{\textbf{DF}} & \multicolumn{1}{l}{\textbf{EMS}}& \multicolumn{1}{l}{$\bm{E_{\gamma}}$}\\
 \midrule
 Between Cag & $1$ & $\sigma^2+2\sigma_{CA}^2+8\sigma_{C}^2$ &\\ \hline
 Between Cag:Ani &  &  &\\
 \quad Trt & $3$ & $\sigma^2+2\sigma_{CA}^2+4\theta_{\gamma}$ &$1$\\
 \quad Residual & $3$ & $\sigma^2+2\sigma_{CA}^2$ &\\ \hline
 Within Cag.Ani & $8$ & $\sigma^2$ &\\
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
 Between Run &  &  & & \\ \hline
 \quad Between Cag:Ani & $1$ & $\sigma^2+2\sigma_{CA}^2+4\sigma_{R}^2$ & & \\ \hline
 \quad Within Cag.Ani & $2$ & $\sigma^2+4\sigma_{R}^2$ & & \\ \hline
 Within Run &  &  & & \\ \hline
 \quad Between Cag &  &  & & \\
 \quad \quad Tag & $1$ & $\sigma^2+2\sigma_{CA}^2+8\sigma_{C}^2+4\theta_{\gamma}$ &$1$ & \\ \hline
 \quad Between Cag:Ani &  &  & & \\
 \quad \quad Trt & $3$ & $\sigma^2+2\sigma_{CA}^2+4\theta_{\tau}$ & & $1$\\
 \quad \quad Residual & $2$ & $\sigma^2+2\sigma_{CA}^2$ & & \\ \hline
 \quad Within Cag.Ani &  &  & & \\
 \quad \quad Tag & $2$ & $\sigma^2+4\theta_{\gamma}$ &$1$ & \\
 \quad \quad Residual & $4$ & $\sigma^2$ & & \\
 \bottomrule
 \end{tabular}
 \label{tab:}
\end{table}




