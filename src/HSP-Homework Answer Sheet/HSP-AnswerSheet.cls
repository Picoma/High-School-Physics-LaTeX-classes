 
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{HSP-AnswerSheet}[2021/07/27]

\LoadClass[12pt,a4paper,landscape]{article}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Paquets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defaults
\RequirePackage{amsmath}
\RequirePackage{amsfonts}
\RequirePackage{amssymb}
\RequirePackage{graphicx}

\RequirePackage[margin=2cm]{geometry}
\RequirePackage[colorlinks=true,linkcolor=black]{hyperref}

\RequirePackage[dvipsnames,table,xcdraw]{xcolor}
\definecolor{cPrim}{RGB}{15, 115, 145}

% Fonts
\RequirePackage[no-math]{fontspec}
\setmainfont{TeX Gyre Heros}
\setsansfont{TeX Gyre Heros}
\RequirePackage[helvet]{sfmath}

% Sections and title style
\usepackage{titlesec}
\titleformat*{\subsection}%
{\bfseries\large\color{White!0!cPrim}}% format applied to label+text

% Multicols
\RequirePackage{multicol}

% SI Units
\RequirePackage{cancel}
\RequirePackage{siunitx}
\sisetup{
	output-decimal-marker = {,},
	inter-unit-product = {\cdot},
	separate-uncertainty=true,
	math-ohm=\text{\textohm}
}
\newcommand{\SIgray}[2]{\SI[unit-color=Gray]{#1}{#2}}

% Enumitem (pretty lists, better labeling)
\usepackage{enumitem}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

\DeclareOption{OpenDyslexic}{% Opendyslexic font
	\AtBeginDocument{
		\setmainfont{OpenDyslexicAlta}
	}
}

\DeclareOption{doublespacing}{% 1.5 spacing
	\AtBeginDocument{
		\usepackage{setspace}
		\doublespacing
	}
}

\ProcessOptions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Commandes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\renewcommand{\maketitle}{
	\begin{center}
		{\huge Homework Answer Sheet}\\[0.5\baselineskip]
		\large\@title
	\end{center}
	
	\setlength{\columnseprule}{1pt}
	\setlength{\columnsep}{1cm}
}