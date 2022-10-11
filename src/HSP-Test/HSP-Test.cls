\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{HSP-Test}[2021/06/17]

\LoadClass[12pt,a4paper]{article}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Packages %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defaults
\RequirePackage{amsmath}
\RequirePackage{amsfonts}
\RequirePackage{amssymb}
\RequirePackage{graphicx}

\RequirePackage[left=2.00cm, right=2.00cm, top=2.00cm, bottom=2.00cm]{geometry}
\savegeometry{normal}
\RequirePackage[colorlinks=true,linkcolor=cPrim,urlcolor=cPrim]{hyperref}

\RequirePackage[dvipsnames,table,xcdraw]{xcolor}
\definecolor{cPrim}{RGB}{15, 115, 145}

% Booleans to enable switching behaviour :
\usepackage{etoolbox}

% Custom environments, controlled by etoolbox's switches
\usepackage{environ}

% Fonts
\RequirePackage{mathptmx}
\RequirePackage{roboto}
\RequirePackage{fontspec}
\setmainfont{TeX Gyre Termes}
\setsansfont{Roboto}
%\renewcommand{\familydefault}{\sfdefault}
%\RequirePackage{sfmath}

% TikZ + custom environment for documents : "doc"
\RequirePackage{tikz}

\RequirePackage[framemethod=TikZ]{mdframed}
\mdfdefinestyle{theoremstyle}{%
	linecolor=cPrim,linewidth=1pt,%
	roundcorner=2pt,
	frametitlebackgroundcolor=cPrim!50,
	startinnercode=\sffamily\small}
\mdtheorem[style=theoremstyle]{doc}{\sffamily\small Document}

% Custom sections styling
\usepackage{titlesec}

%\renewcommand{\thesection}{\normalsize\Roman{section}}
\titleformat{\section}%
{\sffamily\bfseries\Large\color{cPrim}}% Format applied to label+text
{Exercise \thesection : }% Exercise label
{0em}% Distance between label and title
{}% Text before the title
[]% Text after the title

\renewcommand{\thesubsection}{\arabic{subsection}}
\titleformat{\subsection}%
{\sffamily\bfseries\large\color{cPrim}}% Format applied to label+text
{\thesubsection. }% Exercise label
{0em}% Distance between label and title
{}% Text before the title
[]% Text after the title

\titleformat{\paragraph}[runin]%
{\sffamily\bfseries\small}% Format applied to label+text
{\theparagraph}% Exercise label
{0em}% Distance between label and title
{}% Text before the title
[]% Text after the title

% Fancy Header
\newcommand{\@class}{} % Variable @classe
\newcommand{\class}[1]{% Modifies @classe, declared above
	\renewcommand{\@class}{#1}
}
\newcommand{\@version}{} % Variable @version
\newcommand{\version}[1]{% Modifies @version, declared above
\renewcommand{\@version}{#1}
}

\newtoggle{withVersions}

\RequirePackage{lastpage}
\RequirePackage{fancyhdr}
\pagestyle{fancy}
	\fancyhf{}
	\lhead{\sffamily\@class}
	\rhead{\sffamily\@title}
	\lfoot{\sffamily\@date}
	\rfoot{\sffamily\thepage/\pageref{LastPage}}
	\cfoot{\sffamily\iftoggle{withVersions}{Sujet \fbox{\@version}{}}}

% Multicols
\RequirePackage{multicol}

% Better enumerate (with custom numbering)
\usepackage{enumitem}

%%%% Automatic calculations with units
% Setting up the formatting with siunitx
\RequirePackage{siunitx}
\sisetup{
	output-decimal-marker = {,},
	inter-unit-product = {\,\cdot\,},
	separate-uncertainty=true,
%	math-ohm=\text{\textohm}
}
\newcommand{\SIgray}[2]{\SI[unit-color=Gray]{#1}{#2}}

% PythonTeX : setup & wrappers for conviniently and concisely
\RequirePackage{pythontex}
\begin{pycode}
# see physics.py for more details
from physics import Q
\end{pycode}

\newcommand{\Py}[2][]{\pyc{print((#2).latex(addition='#1'))}}
\newcommand{\PyGray}[2][]{\pyc{print((#2).latex(gray=True, addition='#1'))}}

% If quantities cancellation is needed
\RequirePackage{cancel}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%\ProcessOptions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Commands %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\@duration}{} % Variable @classe
\newcommand{\duration}[1]{% Modifies @classe, declared above
	\renewcommand{\@duration}{#1}
}

\newtoggle{answerSheet}
\geometry{left=2.00cm, right=2.00cm, top=2.00cm, bottom=2.00cm, marginpar=1.50cm, includemp=true}
\AtBeginDocument{
	\iftoggle{answerSheet}{}{\loadgeometry{normal}}
}

\NewEnviron{answer}{
	% Affichage de la réponse
	\iftoggle{answerSheet}{
		\begin{mdframed}[frametitle={\colorbox{white}{\sffamily\small\space Answer\space}},
					 frametitleaboveskip=-\ht\strutbox]
			\sffamily\small\BODY
		\end{mdframed}
	}
	{}
}

% Points in margin
\RequirePackage{marginnote}
\newcommand{\pts}[2][/]{
	% Affichage, dans la marge, du nombre de points
	\iftoggle{answerSheet}{\marginnote{\hfill#1#2}}{}
}

\renewcommand{\maketitle}[1]{
	% Don't display student's reminders if it's the answer sheet
	\iftoggle{answerSheet}
	{\begin{center}\sffamily
			{\huge \@title{} -- Answer sheet}\\[0.5\baselineskip]
			\iftoggle{withVersions}{\large \emph{Sujet \fbox{\@version}}}{}
		\end{center}
		
		\begin{center}
			\rule{10cm}{0.5pt}
	\end{center}}
	% 
	{\begin{center}\sffamily
		{\huge \@title}\\
		\small Duration : \@duration
	\end{center}
	
	#1
		
	\begin{center}
		\rule{10cm}{0.5pt}
	\end{center}}
}
