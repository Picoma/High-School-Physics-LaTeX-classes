\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{HSP-Experiment}[2021/06/17]

\LoadClass[a4paper]{article}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Packages %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defaults
\RequirePackage{amsmath}
\RequirePackage{amsfonts}
\RequirePackage{amssymb}
\RequirePackage{graphicx}

\RequirePackage[left=2.00cm, right=2.00cm, top=2.00cm, bottom=2.00cm]{geometry}
\RequirePackage[colorlinks=true,linkcolor=black]{hyperref}

\RequirePackage[dvipsnames,table,xcdraw]{xcolor}
\definecolor{cPrim}{RGB}{15, 115, 145}

% TikZ + custom environment for frames
\RequirePackage{tikz}

\RequirePackage[framemethod=TikZ]{mdframed}
\mdfdefinestyle{theoremstyle}{%
	linecolor=cPrim,linewidth=1pt,%
	roundcorner=2pt,
	frametitlebackgroundcolor=cPrim!50}
\mdtheorem[style=theoremstyle,
%		   startinnercode={\hspace{\parindent}}
		  ]{doc}{\hspace{-\parindent}Document}

\newenvironment{ToRemember}[1]
{\begin{mdframed}[
				frametitle={\colorbox{cPrim!10}{\space Remember : #1 \space}},
				frametitleaboveskip=-\ht\strutbox,
				backgroundcolor=cPrim!10,
				linecolor=cPrim,
				startinnercode={\hspace{\parindent}}]}
{\end{mdframed}}

% Fonts
\renewcommand{\familydefault}{\sfdefault}
\RequirePackage{sfmath}

% Custom sections styling
\usepackage{titlesec}
% section format
%\renewcommand{\thesection}{\normalsize\Roman{section}}
\titleformat{\section}%
{\sffamily\bfseries\Large}% format applied to label+text
{Partie \thesection : }% label
{0em}% horizontal separation between label and title body
{}% before the title body
[]% after the title body

% Fancy Header
\newcommand{\@class}{} % Variable @classe
\newcommand{\class}[1]{% Modifies @classe, declared above
	\renewcommand{\@class}{#1}
}

\RequirePackage{lastpage}
\RequirePackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{}
\lhead{\@class}
\rhead{Travaux pratiques}
\lfoot{\@date}
\rfoot{\thepage/\pageref{LastPage}}

% Multicols
\RequirePackage{multicol}

% SI Units
\RequirePackage{siunitx}
\sisetup{
	output-decimal-marker = {,},
	inter-unit-product = { \text{$\,\cdot\,$}},
	separate-uncertainty=true
}
\newcommand{\SIGray}[2]{\SI[unit-color=Gray]{#1}{#2}}

% Better enumerate (with custom numbering)
\usepackage{enumitem}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\DeclareOption{code}{% import listings
	\AtBeginDocument{
		\RequirePackage{listings}
		\definecolor{codegreen}{rgb}{0,0.6,0}
		\definecolor{codegray}{rgb}{0.5,0.5,0.5}
		\definecolor{codepurple}{rgb}{0.58,0,0.82}
		\definecolor{backcolour}{rgb}{0.95,0.95,0.92}
		
		\lstdefinestyle{mystyle}{
			belowcaptionskip=1\baselineskip,
			breaklines=true,
			frame=L,
			xleftmargin=\parindent,
			language=C,
			showstringspaces=false,
			basicstyle=\footnotesize\ttfamily\normalsize,
			keywordstyle=\bfseries\color{green!40!black},
			commentstyle=\itshape\color{purple!40!black},
			identifierstyle=\color{blue},
			stringstyle=\color{orange},
		}
		\lstset{style=mystyle}
	}
}

\ProcessOptions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Commands %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\renewcommand{\maketitle}{
	\begin{center}
		{\huge \@title}\\
		\rule{5cm}{0.5pt}
	\end{center}
}