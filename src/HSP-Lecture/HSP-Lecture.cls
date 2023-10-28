\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{HSP-Lecture}[2021/07/27]

\LoadClass[12pt,a4paper]{article}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Paquets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defauts
\RequirePackage{amsmath}
\RequirePackage{amsfonts}
\RequirePackage{amssymb}
\RequirePackage{graphicx}

\RequirePackage[margin=2cm]{geometry}
\RequirePackage[colorlinks=true,linkcolor=black]{hyperref}

\RequirePackage[dvipsnames,table,xcdraw]{xcolor}
\definecolor{cPrim}{RGB}{15, 115, 145}

% Polices
\RequirePackage{newtxsf}
\RequirePackage{mathspec}
\setmainfont{TeX Gyre Heros}
\setsansfont{TeX Gyre Heros}
\setmathfont{TeX Gyre Heros}
\RequirePackage[italic]{mathastext}


% Enumitem, pour de meilleures listes
\RequirePackage{enumitem}

% Style de titre de sections
\usepackage{titlesec}
% section format
\renewcommand{\thesection}{\normalsize\Roman{section}}
\titleformat{\section}%
{\sffamily\bfseries\Large\color{White!0!cPrim}}% format applied to label+text
{\llap{\colorbox{White!0!cPrim}{\parbox{1.25cm}{\hfill\color{white}\thesection}}}}% label
{1em}% horizontal separation between label and title body
{}% before the title body
[]% after the title body

% subsection format
%\renewcommand{\thesubsection}{\thesection.\arabic{section}}
\titleformat{\subsection}%
{\sffamily\bfseries\large\color{White!37!cPrim}}% format applied to label+text
{\llap{\colorbox{White!37!cPrim}{\parbox{1.25cm}{\hfill\color{white}\thesubsection}}}}% label
{1em}% horizontal separation between label and title body
{}% before the title body
[]% after the title body

% Conditionnels et "environ" pour sélection conditionnelle du contenu
\RequirePackage{environ}
\RequirePackage{etoolbox}
\newtoggle{versionProf}

% TikZ + Environnement "coursBase" (Cadre de cours, utilisé avec conditionnels
\RequirePackage{tikz}

\RequirePackage[framemethod=TikZ]{mdframed}
\newmdenv[
	skipabove=7pt,
	skipbelow=7pt,
	rightline=false,
	leftline=true,
	topline=false,
	bottomline=false,
	linecolor=cPrim,
	backgroundcolor=White!90!cPrim,
	innerleftmargin=5pt,
	innerrightmargin=5pt,
	innertopmargin=5pt,
	leftmargin=0cm,
	rightmargin=0cm,
	linewidth=2pt,
	innerbottommargin=5pt,
%	startinnercode={\hspace{\parindent}},
	]
	{coursBase}

% Multicols
\RequirePackage{multicol}

% SI Units
\RequirePackage{cancel}
\RequirePackage{siunitx}
\sisetup{
	output-decimal-marker = {,},
	inter-unit-product = { \text{$\,\cdot\,$}}
}
\newcommand{\SIgray}[2]{\SI[unit-color=Gray]{#1}{#2}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

\DeclareOption{tikz-cross}{% import listings
	\usetikzlibrary{shapes.misc}
	
	\tikzset{cross/.style={cross out, draw=black, minimum size=2*(#1-\pgflinewidth), inner sep=0pt, outer sep=0pt},
		%default radius will be 1pt. 
		cross/.default={1pt}}
}

\DeclareOption{OpenDyslexic}{% Opendyslexic font, 1.5 spacing
	\AtBeginDocument{
		\setmainfont{OpenDyslexicAlta}
	}
}

\DeclareOption{doublespacing}{% Opendyslexic font, 1.5 spacing
	\AtBeginDocument{
		\usepackage{setspace}
		\doublespacing
	}
}

\ProcessOptions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Commandes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Style de titre principal
\newcommand{\@chapNumber}{} % Variable @chapNumber

\newcommand{\chapNumber}[1]{% Modifie @chapNumber, déclaré en début de document
	\renewcommand{\@chapNumber}{#1}
}

\renewcommand{\maketitle}{
	\noindent Chapter \@chapNumber\\
	\begin{center}
		{\huge \@title}\\[\baselineskip]
		\rule{5cm}{0.5pt}
	\end{center}
}

% Environnement remarque et attention, pour détails graphiques
\newenvironment{remarque}{\par\small\vspace{0pt} % Vertical white space above the remark and smaller font size
\begin{list}{}{
		\leftmargin=35pt % Indentation on the left
		\rightmargin=25pt}\item\ignorespaces % Indentation on the right
	\makebox[-2.5pt]{\begin{tikzpicture}[overlay]
		\node[draw=cPrim!60,line width=1pt,circle,fill=cPrim!25,font=\sffamily\bfseries,inner sep=2pt,outer sep=0pt] at (-15pt,4pt){\textcolor{cPrim}{R}};\end{tikzpicture}} % Orange R in a circle
	\advance\baselineskip -1pt}{\end{list}\vskip5pt} % Tighter line spacing and white space after remark

\newenvironment{attention}{\par\small % smaller font size
	\begin{list}{}{
			\leftmargin=35pt % Indentation on the left
			\rightmargin=25pt}\item\ignorespaces % Indentation on the right
		\makebox[-2.5pt]{\begin{tikzpicture}[overlay, baseline=(current bounding box)]
			\node[draw=cPrim!85,line width=1pt,circle,fill=cPrim!25,font=\sffamily\bfseries,inner sep=2pt,outer sep=0pt] at (-15pt,4pt){\textcolor{cPrim}{\hspace{-2pt}\textbf{!}}};\end{tikzpicture}} % Orange R in a circle
		\advance\baselineskip -1pt}{\end{list}} % Tighter line spacing

% Cours, affiché si versionProf=true
\NewEnviron{cours}[1]{
	\begin{coursBase}
		\iftoggle{versionProf}
		{\BODY}
		{\vspace{#1}}
	\end{coursBase}
}