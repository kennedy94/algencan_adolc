# algencan_adolc
Algencan 3.1.1 using Automatic Differenciation library ADOL-C 2.7.3

In this file we present the complete instalation of Algencan 3.1.1 and ADOL-C 2.7.3 and how to use those libraries together in Linux Ubuntu.

Algencan 3.1.1 is a nonlinear programming solver and ADOL-C 2.7.3 is used to evaluate sparses Hessians in Algencan methods.

# INSTALLATION

Before installing ADOL-C, we first need to install two libraries: 1. Boost 1.75.0; 2. Colpack 1.0.10.

## Boost 1.75.0 INSTALL
\begin{enumerate}
    \item Download boost\_1\_75\_0.tar.bz2 \url{https://www.boost.org/users/download/};
    \item Extract the file in the desired directory \emph{tar --bzip2 -xf /path/to/boost\_1\_75\_0.tar.bz2};
    \item Open boost directory: \emph{cd path/to/boost\_1\_75\_0};
    \item Run \emph{bootstrap.sh}: \emph{./bootstrap.sh --prefix=path/to/installation/prefix};
    \item Run \emph{b2} to generate external libraries \textbackslash lib\textbackslash: \emph{./b2 install} 
\end{enumerate}

\subsection{Instalação ColPack}

\begin{enumerate}
    \item Baixar o arquivo no link \url{https://github.com/CSCsw/ColPack/releases}
    \item Extrair o arquivo;
    \item Executar na pasta do Colpack os seguintes comandos
    \begin{description}
        \item \emph{autoreconf -vif}
        \item \emph{./configure --prefix=/path/to/install/}
        \item \emph{make}
        \item \emph{make install}
    \end{description}
    
\end{enumerate}

\subsection{Instalação ADOL-C}

Com a instalação bem sucedida das duas bibliotecas anteriores podemos, então, instalar o ADOL-C.

\begin{enumerate}
    \item Baixar o código no link \url{https://github.com/coin-or/ADOL-C};
    \item Extrair o arquivo e executar os seguintes comandos na pasta de instalação:
    \begin{description}
        \item \emph{autoreconf -fi}
        \item \emph{./configure --enable-sparse --with-boost=CAMINHOPARABOOST --with-colpack=CAMINHOPARACOLPACK}
        \item \emph{make}
        \item \emph{make install}
    \end{description}
\end{enumerate}

Se a instalação for bem sucedida as bibliotecas referentes aos drivers do ADOL-C serão geradas na pasta \emph{/home/USER/adolc\_base} 

\subsection{Instalação ALGENCAN (Versão 3.1.1 usada como referência)}

\begin{itemize}
    \item Baixar o arquivo referente ao ALGENCAN no link \url{https://www.ime.usp.br/~egbirgin/tango/codes.php};
    \item Extrair arquivo e executar o comando \emph{make} no diretório do ALGENCAN;
    \item Se a instalação for bem sucedida a biblioteca libalgencan.a será gerada na pasta ALGENCAN/lib/
\end{itemize}

\section{Linkagem e compilação}

\subsection{Linkagem via terminal}

Para compilar um código main.cpp que usa a biblioteca ADOL-C, basta utilizar o seguinte comando:

\emph{g++ -w -I/home/USER/adolc\_base/include -o main main.cpp -Wl,--rpath -Wl,/home/USER/adolc\_base/lib64 -L/home/USER/adolc\_base/lib64 -ladolc}

Para compilar um código main.c que usa o solver ALGENCAN basta utilizar o seguinte comando:

\emph{gcc -O3 main.c -L\$ALGENCAN/lib -lalgencan -lgfortran -lm -o algencan}

\subsection{Linkagem das bibliotecas com Codeblocks (Versão 20.03 usada como referência)}

\begin{enumerate}
    \item Criar projeto de C/C++ Console Application.
    \item Em \emph{Build Options$>>$Linker Settings}, adicione as bibliotecas:
    \begin{description}
        \item adolc.so
        \item libColPack.so
        \item libalgencan.a
        \item gfortran
    \end{description}
    \item Em \emph{Build Options$>>$Search Directories}, adicione os caminhos:
    \begin{description}
        \item /home/USER/adolc\_base/include;
        \item /home/USER/adolc\_base/lib64;
    \end{description}
    \noindent em \emph{Compiler, Linker,} e \emph{Resource compiler}
    
\end{enumerate}
