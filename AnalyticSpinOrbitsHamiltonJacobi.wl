(* ::Package:: *)

(* ::Section:: *)
(*Begin package*)


BeginPackage["AnalyticSpinOrbitsHamiltonJacobi`"];


(* ::Text:: *)
(*If you make use of this package, please acknowledge "Piovano" arXiv:2510.09597 (https://arxiv.org/abs/2510.09597 ), and "Piovano" (https://arxiv.org/abs/2603.04682 )*)


(* ::Text:: *)
(*IMPORTANT NOTE: the following package include all the functions to compute the analytic spin-corrections to the orbits (trajectories and velocities), constants of motion and frequencies for periodic, homoclinic, and plunging orbits in the case of near equatorial motion.  The package provide the contributions to both the parallel and orthogonal components of the secondary spin. *)
(*For periodic orbits, the analytic spin-corrections are available in the fixed turning points (or "FT"), fixed constants of motion (or "FC"), and fixed eccentricity (or "FE") spin-gauges*)
(*Homoclinic and plunging orbits are only available for the FE spin-gauge since other spin-gauge introduce discontinuities in parts of the phase-space.*)


(* ::Subsection::Closed:: *)
(*Periodic and homoclinic orbits*)


\[Delta]orbelements::usage = "\[Delta]orbelements[a,p,e,x,spingauge] calculates the linear corrections to the orbital elements (p,e,x) for fixed frequency (spingauge=FF) and fixed eccentricity (spingauge=FE)";


(* ::Subsubsection::Closed:: *)
(*Fixed turning points*)


KerrNearEqSpinOrbitCorrFTPer::usage = "KerrNearEqSpinOrbitCorrFTPer[a, p, e, x] calculates the linear corrections to periodic orbits in the fixed turning points spin-gauge. The initial radius at \[Lambda]=0 is the geodesic periastron.";


(* ::Subsubsection::Closed:: *)
(*Fixed constants of motion*)


KerrNearEqSpinOrbitCorrFCPer::usage = "KerrNearEqSpinOrbitCorrFCPer[a, p, e, x] calculates the linear corrections to periodic orbits in the fixed constants of motion spin-gauge. The initial radius at \[Lambda]=0 is the geodesic periastron.";


(* ::Subsubsection::Closed:: *)
(*Fixed frequency*)


KerrNearEqSpinOrbitCorrFFPer::usage = "KerrNearEqSpinOrbitCorrFFPer[a, p, e, x] calculates the linear corrections to periodic orbits in the fixed frequency spin-gauge. The initial radius at \[Lambda]=0 is the geodesic apoastron.";


(* ::Subsubsection::Closed:: *)
(*Fixed eccentricity*)


KerrNearEqSpinOrbitCorrFEPer::usage = "KerrNearEqSpinOrbitCorrFEPer[a, p, e, x] calculates the linear corrections to periodic orbits in the fixed eccentricity spin-gauge. The initial radius at \[Lambda]=0 is the geodesic apoastron.";


KerrNearEqSpinOrbitCorrFEHom::usage = "KerrNearEqSpinOrbitCorrFEHom[a, p, e, x] calculates the linear corrections to homoclinic orbits.";


(* ::Subsubsection::Closed:: *)
(*R gauge*)


KerrNearEqSpinOrbitCorrRPer::usage = "KerrNearEqSpinOrbitCorrRPer[a, p, e, x] calculates the linear corrections to periodic orbits in the fixed eccentricity spin-gauge. The initial radius at \[Lambda]=0 is the geodesic apoastron.";


KerrNearEqSpinOrbitCorrRHom::usage = "KerrNearEqSpinOrbitCorrRHom[a, p, e, x] calculates the linear corrections to homoclinic orbits.";


(* ::Subsection::Closed:: *)
(*Bound plunges*)


KerrNearEqSpinOrbitCorr\[Delta]IBCO::usage = "KerrNearEqSpinOrbitCorr\[Delta]IBCO[a, p, x] provides the linear shifts to constants of motion and IBCO radius for parabolic orbits.";


KerrNearEqSpinOrbitCorrISCOPlunge::usage = "KerrNearEqSpinOrbitCorrISCOPlunge[a, x] calculates the linear corrections to ISCO plunges.";


KerrEqISCOPlungePartialExp::usage = "KerrEqISCOPlungePartialExp[a, S, x] calculates the partially expanded solutions for ISCO plunges.";


(* ::Subsubsection::Closed:: *)
(*Fixed eccentricity*)


KerrNearEqSpinOrbitCorrFECritPlunge::usage = "KerrNearEqSpinOrbitCorrFECritPlunge[a, p, e, x] calculates the linear corrections to critical plunges in the fixed eccentricity gauge.";


KerrNearEqSpinOrbitCorrFEPlunge::usage = "KerrNearEqSpinOrbitCorrFEPlunge[a, p, e, x] calculates the linear corrections to generic plunges in the fixed eccentricity gauge.";


KerrNearEqSpinOrbitCorrFEPlungeBoundOrbit::usage = "KerrNearEqSpinOrbitCorrFEBoundOrbit[a, p, e, x] calculates the linear corrections to plunges related to bound orbits (PBOs) in the fixed eccentricity gauge.";


KerrNearEqSpinOrbitCorrFEPlungeBoundOrbitDoubleRoot::usage = "KerrNearEqSpinOrbitCorrFEDoubleRoot[a, p, x] calculates the linear corrections to plunges related to bound orbits (PBOs) with a double root in the fixed eccentricity gauge.";


(* ::Subsubsection::Closed:: *)
(*Fixed eccentricity - partial expansion*)


KerrEqFECritplungePartialExp::usage = "KerrEqFECriticalplungePartialExp[a, S, p, e, x] calculates the partially expanded solutions for critical plunges in the fixed eccentricity gauge.";


KerrEqFEPlungeBoundOrbitPartialExp::usage = "KerrEqFEPlungeBoundOrbitPartialExp[a, S, p, e, x] calculates the partially expanded solutions fo plunges related to bound orbits (PBOs) in the fixed eccentricity gauge.";


KerrEqFEPlungePartialExp::usage = "KerrEqFEPlungePartialExp[a, S, p, e, x] calculates the partially expanded solutions fo plunges related to bound orbits (PBOs) in the fixed eccentricity gauge.";


(* ::Subsection::Closed:: *)
(*Begin package*)


Begin["`Private`"];


(* ::Section::Closed:: *)
(*Geodesic quantities*)


(* ::Subsubsection::Closed:: *)
(*Coordinate-time geodesic velocity - general expression*)


Vtrgfun[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];

	EEg((r^2+a^2)^2/(r^2-2r+a^2))-(2a*r)/(r^2-2r+a^2)Lzg-EEg*a^2
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal geodesic velocity - general expression*)


V\[Phi]rgfun[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];

	a/(r^2-2r+a^2)(EEg(r^2+a^2)-a*Lzg)-a*EEg+Lzg
]


(* ::Subsection::Closed:: *)
(*ISCO radius*)


(* ::Text:: *)
(*See "Rotating Black Holes: Locally Non-rotating Frames, Energy Extraction, and Scalar Synchrotron Radiation" , Bardeen, J. M., Press, W. H., & Teukolsky, S. A., Astrophysical Journal, Vol. 178, pp. 347-370 (1972)*)


ISCOradius[a_,x_]:=Module[{Z1,Z2},
	Z1=1+(1-a^2)^(1/3)((1+a)^(1/3)+(1-a)^(1/3));
	Z2=Sqrt[(3a^2+Z1^2)];

	3+Z2-RealSign[x]Sqrt[(3-Z1)(3+Z1+2Z2)]
]


(* ::Subsection::Closed:: *)
(*IBCO radius*)


IBCOradius[a_,x_]:=2-RealSign[x]a+2Sqrt[1-RealSign[x]a]


(* ::Subsection::Closed:: *)
(*Photon radius*)


Photonradius[a_,x_]:=4Cos[1/3ArcCos[-RealSign[x]*a]]^2


(* ::Subsection::Closed:: *)
(*Constants of motion*)


EEgfun[a_,p_,e_,xg_]:=Sqrt[(p(-a^2(1-e^2)^2(-5+e^2+3p)+p(4e^4+(-3+p)(p-2)^2-e^2(-8+p^2)))-2(1-e^2)^2*Sqrt[a^2*p(a^2(1+e)^2+p(-2-2e+p))(a^2(1-e)^2+p(-2+2e+p))]RealSign[xg])/(p^2(-4a^2(1-e^2)^2+(3+e^2-p)^2*p))]


Lzgfun[a_,p_,e_,xg_]:=Module[{r1g,r2g,EEg,Lzgzaux,ddinc,ffinc,gginc,hhinc},
	r1g=p/(1-e);
	r2g=p/(1+e);

	EEg=EEgfun[a,p,e,xg];

	ddinc[r_]:=(r^2-2r+a^2)r^2;
	ffinc[r_]:=r^4+a^2(r+2)r;
	gginc[r_]:=2 a r;
	hhinc[r_]:=r(r-2);
	Lzgzaux=Sqrt[(gginc[r2g]^2+hhinc[r2g]ffinc[r2g])EEg^2-hhinc[r2g]ddinc[r2g]];

	(-gginc[r2g] EEg+RealSign[xg]*Lzgzaux)/hhinc[r2g]
]


(* ::Subsubsection::Closed:: *)
(*Equatorial orbits  - complex conjugate roots*)


EEgfuncc[p_,e_]:=Module[{r1g,\[Rho]rg},
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);

	Sqrt[(-2+r1g+2\[Rho]rg)/(r1g+2\[Rho]rg)]
]


Lzgfuncc[a_,p_,e_,x_]:=Module[{r1g,\[Rho]rg,\[CapitalDelta]r1g},
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[CapitalDelta]r1g=a^2+(-2+r1g)r1g;

	2RealSign[x]*Sqrt[((-2r1g+r1g^2)r1g*\[Rho]rg-2RealSign[x]Sqrt[a^2*r1g*\[CapitalDelta]r1g*\[Rho]rg(-2+r1g+2 \[Rho]rg)]+a^2(-2+r1g+2\[Rho]rg+r1g*\[Rho]rg))/((-2+r1g)^2(r1g+2\[Rho]rg))]
]


\[Rho]igfun[a_,p_,e_,x_]:=Module[{r1g,\[Rho]rg,\[CapitalDelta]r1g},
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[CapitalDelta]r1g=a^2+(-2+r1g)r1g;

	-Sqrt[(1/(-2+r1g)^2(-((-2+r1g)(r1g(-4+\[Rho]rg)-2\[Rho]rg)\[Rho]rg)-4RealSign[x]Sqrt[a^2*r1g*\[CapitalDelta]r1g*\[Rho]rg(-2+r1g+2 \[Rho]rg)]+a^2((-2+r1g)r1g+2(2+r1g)\[Rho]rg)))]
]


dEEdpfuncc[p_,e_]:=(1-e^2)/((3-e)p^2*EEgfuncc[p,e])


dLzdpfuncc[a_,p_,e_,x_]:=Module[{Lzg,r1g,\[Rho]rg,\[CapitalDelta]r1g},
	Lzg=Lzgfuncc[a,p,e,x];
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[CapitalDelta]r1g=a^2+(-2+r1g)r1g;

	2/Lzg((-p^3(-2+2e+p)(-4+4e+p)+a^2(-1+e)^2(4(-1+e)^2 (1+e)+6(-1+e^2)p-4(-2+e)p^2+p^3))/((-3+e)p^2(-2+2e+p)^3)+(2(-1+e)^3 (1+e)(-2+2e+3p)RealSign[x]Sqrt[a^2*r1g*\[CapitalDelta]r1g*\[Rho]rg(-2+r1g+2\[Rho]rg)])/((-3+e)p^2(-2+2e+p)^3)- RealSign[x]/((-3+e)p (-2+2e+p)^2 Sqrt[a^2*r1g*\[CapitalDelta]r1g*\[Rho]rg(-2+r1g+2\[Rho]rg)]) (a^2*p(a^2 (-1+e)^2 (-4+4e^2+9p-3e*p)+p(12(-1+e)^2(1+e)+32(-1+e)p-5(-3+e)p^2)))/((-1+e)(1+e)))
]


(* ::Subsection::Closed:: *)
(*Derivatives constants of motion*)


(* ::Subsubsection::Closed:: *)
(*Derivatives of the energy*)


dEgdr1gfun[a_,p_,e_,xg_]:=Module[{r1g,r2g,EEg,faux,dfauxdr1g,daux,ddauxdr1g,gaux,dgauxdr1g},
	r1g=p/(1-e);
	r2g=p/(1+e);

	EEg=EEgfun[a,p,e,xg];

	faux=4r1g*r2g((-2+r1g)(-2+r2g)(r1g+r2g)(r1g^2(-2+r2g)+r1g(-2+r2g)r2g-2r2g^2)+a^2(r1g^2(4-6r2g)-6r1g(-2+r2g)r2g+4r2g^2));
	gaux=(r1g+r2g)RealSign[xg]Sqrt[(a^2r1g^3r2g^3(a^2+(-2+r1g)r1g)(a^2+(-2+r2g)r2g))/(r1g+r2g)];
	daux=-8a^2*r1g^2*r2g^2(r1g+r2g)+r1g*r2g (r1g^2(-2+r2g)+r1g(-2+r2g)r2g-2r2g^2)^2;


	dfauxdr1g=4r2g(2a^2(r1g^2(6-9r2g)-6r1g(-2+r2g)r2g+2r2g^2)+(-2+r2g)(5r1g^4(-2+r2g)+8r1g^3(-2+r2g)(-1+r2g)-8r1g(-2+r2g)r2g^2+4r2g^3+3r1g^2*r2g(8+(-8+r2g)r2g)));
	ddauxdr1g=r2g (5r1g^4(-2+r2g)^2+8r1g^3(-2+r2g)^2*r2g+4r2g^4-8r1g*r2g^2(2a^2+(-2+r2g)r2g)+3r1g^2*r2g(-8a^2+(-6+r2g)(-2+r2g)r2g));
	
	dgauxdr1g=(r2g^(3/2)(r1g(6r1g^2+5r1g(-2+r2g)-8r2g)+a^2(4r1g+3r2g)))/(2(a^2+(-2+r1g)r1g)) Sqrt[(a^2*r1g(a^2+(-2+r1g)r1g)(a^2+(-2+r2g)r2g))/(r1g+r2g)]RealSign[xg];

	1/2 EEg(-(ddauxdr1g/daux)+(dfauxdr1g-16Sqrt[2]*dgauxdr1g)/(faux-16Sqrt[2]*gaux))
]


dEgdr2gfun[a_,p_,e_,xg_]:=Module[{r1g,r2g,EEg,faux,dfauxdr2g,daux,ddauxdr2g,gaux,dgauxdr2g},
	r1g=p/(1-e);
	r2g=p/(1+e);

	EEg=EEgfun[a,p,e,xg];

	faux=4r1g*r2g((-2+r1g)(-2+r2g)(r1g+r2g)(r1g^2(-2+r2g)+r1g(-2+r2g)r2g-2r2g^2)+a^2(r1g^2(4-6r2g)-6r1g(-2+r2g)r2g+4r2g^2));
	gaux=(r1g+r2g)RealSign[xg]Sqrt[(a^2r1g^3r2g^3(a^2+(-2+r1g)r1g)(a^2+(-2+r2g)r2g))/(r1g+r2g)];
	daux=-8a^2*r1g^2*r2g^2(r1g+r2g)+r1g*r2g (r1g^2(-2+r2g)+r1g(-2+r2g)r2g-2r2g^2)^2;


	dfauxdr2g=4r1g(2a^2(r1g^2(2-6r2g)+3r1g(4-3r2g)r2g+6r2g^2)+(-2+r1g)(4r1g^3-8(-2+r1g)r1g^2*r2g+3r1g(8+(-8+r1g)r1g)r2g^2+8(-2+r1g)(-1+r1g)r2g^3+5(-2+r1g)r2g^4));
	ddauxdr2g=r1g(4r1g^4-8r1g^2(2a^2+(-2+r1g)r1g)r2g+3r1g(-8a^2+(-6+r1g)(-2+r1g)r1g)r2g^2+8(-2+r1g)^2*r1g*r2g^3+5(-2+r1g)^2*r2g^4);

	dgauxdr2g=(r1g^(3/2)(a^2(3r1g+4r2g)+r2g(-8r1g+5(-2+r1g)r2g+6r2g^2)))/(2(a^2+(-2+r2g)r2g)) Sqrt[(a^2(a^2+(-2+r1g)r1g)r2g(a^2+(-2+r2g)r2g))/(r1g+r2g)]RealSign[xg];

	1/2 EEg(-(ddauxdr2g/daux)+(dfauxdr2g-16Sqrt[2]*dgauxdr2g)/(faux-16Sqrt[2]*gaux))
]


dEEdpfun[a_,p_,e_,xg_]:=Module[{EEg,auxsqr,auxden},
	EEg=EEgfun[a,p,e,xg];
	auxsqr=Sqrt[p^11(a^4(1-e^2)^2+(-4e^2+(-2+p)^2)p^2+2a^2*p(-2+p+e^2(2+p)))];
	auxden=(-4a^2(1-e^2)^2+(3+e^2-p)^2*p);
	
	-(((-14a^2(1-e^2)^2+(12+4e^2-5p)(3+e^2-p)p)((2+2e-p)(3+e^2-p)p(-2+2e+p)-a^2(1-e^2)^2(-5+e^2+3p)))/(EEg*p^2auxden^2))+(2(1-e^2)^2(-14a^2(1-e^2)^2+(12+4e^2-5p)(3+e^2-p)p)a*auxsqr*RealSign[xg])/(EEg*p^8*auxden^2)-(a (1-e^2)^2*p^3(11a^4(1-e^2)^2+48a^2(-1+e^2)p+26(2+a^2+(-2+a^2)e^2)p^2-56p^3+15p^4)RealSign[xg])/(2EEg*auxden*auxsqr)+(-3a^2(1-e^2)^2(-10+2e^2+7p)+p(28(-3+2e^2+e^4)+128p-9(7+e^2)p^2+10p^3))/(2EEg*p^2*auxden)
]


dEEdefun[a_,p_,e_,x_]:=Module[{aux},
	aux=p(a^2 (1+e)^2+p (-2-2 e+p)) (a^2 (-1+e)^2+p (-2+2 e+p));
	(-4 e (-4 a^2 (-1+e^2)+(3+e^2-p) p) (p (-a^2 (-1+e^2)^2 (-5+e^2+3 p)+p (4 e^4+(-3+p) (-2+p)^2-e^2 (-8+p^2)))-2 a Sqrt[aux] (-1+e^2)^2 RealSign[x])+2 e (-4 a^2 (-1+e^2)^2+(3+e^2-p)^2 p) (p (-a^2 (-1+e^2) (-11+3 e^2+6 p)+p (8+8 e^2-p^2))-4 a Sqrt[aux] (-1+e^2) RealSign[x]-(2 a (-1+e^2)^2 p (a^4 (-1+e^2)-2 p^2+a^2 p (2+p)) RealSign[x])/Sqrt[aux]))/(2 p (-4 a^2 (-1+e^2)^2+(3+e^2-p)^2 p)^2 Sqrt[(p (-a^2 (-1+e^2)^2 (-5+e^2+3 p)+p (4 e^4+(-3+p) (-2+p)^2-e^2 (-8+p^2)))-2 a Sqrt[aux] (-1+e^2)^2 RealSign[x])/(-4 a^2 (-1+e^2)^2+(3+e^2-p)^2 p)])
]


(* ::Subsubsection::Closed:: *)
(*Derivatives of the energy - complex conjugate roots*)


dEEdpfuncc[p_,e_]:=(1-e^2)/((3-e)p^2*EEgfuncc[p,e])


dEEdefuncc[p_,e_]:=(-1+6e-e^2)/((e-3)^2*p*EEgfuncc[p,e])


(* ::Subsubsection::Closed:: *)
(*Derivatives of the angular momentum*)


dLzgdr1gfun[a_,p_,e_,xg_]:=Module[{r1g,r2g,EEg,Lzg,dEgdr1g ,faux,dfauxdr1g,daux,ddauxdr1g,gaux,dgauxdr1g},
	r1g=p/(1-e);
	r2g=p/(1+e);

	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	dEgdr1g=dEgdr1gfun[a,p,e,xg];

	(dEgdr1g(-2a*Lzg+EEg*r1g^3+a^2*EEg(2+r1g)))/(2a*EEg+Lzg(-2+r1g))+((r1g-r2g)(2+2(-1+EEg^2)r1g+(-1+EEg^2)r2g))/(4a*EEg+2Lzg(-2+r1g))
]


dLzgdr2gfun[a_,p_,e_,xg_]:=Module[{r1g,r2g,EEg,Lzg,dEgdr2g,faux,dfauxdr2g,daux,ddauxdr2g,gaux,dgauxdr2g},
	r1g=p/(1-e);
	r2g=p/(1+e);

	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	dEgdr2g=dEgdr2gfun[a,p,e,xg];

	(dEgdr2g(-2a*Lzg+EEg*r2g^3+a^2*EEg(2+r2g)))/(2a*EEg+Lzg(-2+r2g))-((r1g-r2g)(2+(-1+EEg^2)r1g+2(-1+EEg^2)r2g))/(4a*EEg+2Lzg(-2+r2g))
]


dLzdpfun[a_,p_,e_,xg_]:=Module[{r2g,EEg,dEEdp,Lzgzaux,ddinc,ffinc,gginc,hhinc},
	r2g=p/(1+e);
	EEg=EEgfun[a,p,e,xg];
	dEEdp=dEEdpfun[a,p,e,xg];

	ddinc[r_]:=(r^2-2r+a^2)r^2;
	ffinc[r_]:=r^4+a^2(r+2)r;
	gginc[r_]:=2 a r;
	hhinc[r_]:=r(r-2);
	Lzgzaux=Sqrt[(gginc[r2g]^2+hhinc[r2g]ffinc[r2g])EEg^2-hhinc[r2g]ddinc[r2g]];

	(-EEg*gginc[r2g]+Lzgzaux*RealSign[xg])/((1+e)hhinc[r2g]^2)(2(1+e-p))/(1+e)-(gginc[r2g] dEEdp)/hhinc[r2g]-EEg/hhinc[r2g](2a)/(1+e)+RealSign[xg]/hhinc[r2g]((hhinc[r2g]*EEg*2ffinc[r2g]*dEEdp)/(2Lzgzaux)-hhinc[r2g]/(2Lzgzaux)(2p(a^2(1+e)^2+p(-3-3e+2p)))/(1+e)^4+(hhinc[r2g]EEg^2)/(2Lzgzaux)((4p^3+2a^2(1+e)^2(1+e+p))/(1+e)^4)+(2 EEg^2*gginc[r2g])/(2Lzgzaux) (2a)/(1+e)+(2gginc[r2g]^2EEg*dEEdp)/(2Lzgzaux)-(-ddinc[r2g]+EEg^2*ffinc[r2g])/(2Lzgzaux)(2(1+e-p))/(1+e)^2)
]


dLzdefun[a_,p_,e_,x_]:=Module[{r2g,EEg,dEEde,Lzgzaux,ddinc,ffinc,gginc,hhinc,dggincde,dhhincde,dLzgzauxde},
	r2g=p/(1+e);
	EEg=EEgfun[a,p,e,x];
	dEEde=dEEdefun[a,p,e,x];

	ddinc[r_]:=(r^2-2r+a^2)r^2;
	ffinc[r_]:=r^4+a^2(r+2)r;
	gginc[r_]:=2 a r;
	hhinc[r_]:=r(r-2);
	Lzgzaux=Sqrt[(gginc[r2g]^2+hhinc[r2g]ffinc[r2g])EEg^2-hhinc[r2g]ddinc[r2g]];

	dLzgzauxde=(p^(3/2) (-a^2 (1+e)^2 (3+3 e-2 p)+dEEde (1+e) EEg p (a^2 (1+e)^2+p (-2-2 e+p))-EEg^2 p (2 a^2 (1+e)^2+p (-5-5 e+3 p))+p (8+8 e^2-10 p+3 p^2-2 e (-8+5 p))))/((1+e)^4 Sqrt[(2+2 e-p+EEg^2 p) (a^2 (1+e)^2+p (-2-2 e+p))]);

	(-2 a dEEde (1+e)^2 (-2+r2g) r2g-2 a (2 p+(1+e) (-2+r2g) r2g) EEg+(dLzgzauxde (1+e)^2 (-2+r2g)+(-2+p+e (-2+r2g)+r2g) Lzgzaux) RealSign[x])/((1+e)p (-2+r2g)^2)
]


(* ::Subsubsection::Closed:: *)
(*Derivatives of the energy - complex conjugate roots*)


dLzdpfuncc[a_,p_,e_,x_]:=Module[{Lzg,r1g,\[Rho]rg,\[CapitalDelta]r1g},
	Lzg=Lzgfuncc[a,p,e,x];
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[CapitalDelta]r1g=a^2+(-2+r1g)r1g;

	2/Lzg((-p^3(-2+2e+p)(-4+4e+p)+a^2(-1+e)^2(4(-1+e)^2 (1+e)+6(-1+e^2)p-4(-2+e)p^2+p^3))/((-3+e)p^2(-2+2e+p)^3)+(2(-1+e)^3 (1+e)(-2+2e+3p)RealSign[x]Sqrt[a^2*r1g*\[CapitalDelta]r1g*\[Rho]rg(-2+r1g+2\[Rho]rg)])/((-3+e)p^2(-2+2e+p)^3)- RealSign[x]/((-3+e)p (-2+2e+p)^2 Sqrt[a^2*r1g*\[CapitalDelta]r1g*\[Rho]rg(-2+r1g+2\[Rho]rg)]) (a^2*p(a^2 (-1+e)^2 (-4+4e^2+9p-3e*p)+p(12(-1+e)^2(1+e)+32(-1+e)p-5(-3+e)p^2)))/((-1+e)(1+e)))
]


dLzdefuncc[a_,p_,e_,x_]:=Module[{Lzg,r1g,\[Rho]rg,\[CapitalDelta]r1g, auxvar},
	Lzg=Lzgfuncc[a,p,e,x];
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[CapitalDelta]r1g=a^2+(-2+r1g)r1g;
  auxvar=Sqrt[(a^2 p^2 (a^2+(-2+r1g) r1g) (-2+r1g+2 \[Rho]rg))/((1-e) (1+e))];
	2/Lzg 1/((-2+r1g)^3 (r1g+2 \[Rho]rg)^2) (1/(-1+e^2)^3 (1-6 e+e^2) p (-2+r1g) (-p^2 (-2+r1g) r1g+a^2 (2-p^2-r1g-2 \[Rho]rg+e^2 (-2+r1g+2 \[Rho]rg))-2 auxvar (-1+e^2) RealSign[x])+1/((-1+e)^3 (1+e)) 2 p (r1g+2 \[Rho]rg) (p^2 (-2+r1g) r1g+a^2 (-2+p^2+r1g+2 \[Rho]rg-e^2 (-2+r1g+2 \[Rho]rg))+2 auxvar (-1+e^2) RealSign[x])+(-2+r1g) (r1g+2 \[Rho]rg) ((2 p^3 (-1+e+p))/((-1+e)^4 (1+e))+((-2+r1g) r1g^3)/(1+e)+((-2+r1g) r1g \[Rho]rg^2)/(-1+e)+a^2 ((1/(-1+e)^2-2/(1+e)^2) p+r1g^2/(1+e)+\[Rho]rg^2/(-1+e))-1/ auxvar a^2 (((1-6 e+e^2) p^3 (a^2+(-2+r1g) r1g))/(-1+e^2)^3+(p^3 (-2+p-e (-2+r1g)+r1g) (-2+r1g+2 \[Rho]rg))/((-1+e)^4 (1+e))+(r1g^2 (a^2+(-2+r1g) r1g) (-2+r1g+2 \[Rho]rg))/(1+e)+((a^2+(-2+r1g) r1g) \[Rho]rg^2 (-2+r1g+2 \[Rho]rg))/(-1+e)) RealSign[x]))
]


(* ::Subsection::Closed:: *)
(*Geodesic frequencies*)


(* ::Subsubsection::Closed:: *)
(*Mino-time radial and polar geodesic frequencies*)


\[CapitalUpsilon]rgfun[a_,p_,e_,xg_]:=Module[{r1g,r2g,EEg,Lzg,r3g,krg},
	r1g=p/(1-e);
	r2g=p/(1+e);

	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];

	r3g=2/(1-EEg^2)-(r1g+r2g);
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);

	\[Pi]/(2EllipticK[krg]) Sqrt[(1-EEg^2)(r1g-r3g)r2g]
]


\[CapitalUpsilon]zgfun[a_,p_,e_,xg_]:=Module[{r1g,r2g,EEg,Lzg},
	r1g=p/(1-e);
	r2g=p/(1+e);

	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];

	Sqrt[a^2(1-EEg^2)+Lzg^2]
]


(* ::Subsubsection::Closed:: *)
(*Mino-time coordinate time geodesic frequency*)


\[CapitalUpsilon]tgrfun[a_,p_,e_,xg_]:=Module[{EEg,Lzg,r1g,r2g,r3g,rp,rm,krg,ellK,hr,hp,hm},
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	ellK=EllipticK[krg];
	hr=(r1g-r2g)/(r1g-r3g);
	hp=hr (r3g-rp)/(r2g-rp);
	hm=hr (r3g-rm)/(r2g-rm);

	(4+a^2)EEg+EEg(1/2 (4+r1g+r2g+r3g)r3g- (r1g*r2g)/2+(r1g-r3g)r2g EllipticE[krg]/(2ellK)+1/2 (4+2/(1-EEg^2))(r2g-r3g) EllipticPi[hr,krg]/ellK)+(2EEg)/(rp-rm) (((4-a Lzg/EEg)rp-2a^2)/(r3g-rp) (1-(r2g-r3g)/(r2g-rp) EllipticPi[hp,krg]/ellK))-(2EEg)/(rp-rm) (((4-a Lzg/EEg)rm-2a^2)/(r3g-rm) (1-(r2g-r3g)/(r2g-rm) EllipticPi[hm,krg]/ellK))
]


\[CapitalUpsilon]tgzfun[a_,p_,e_,xg_]:=-a^2EEgfun[a,p,e,xg]


\[CapitalUpsilon]tgfunLim[a_,p_,e_,xg_]:=Module[{EEg,Lzg,r1g,r2g,rp,rm,krg,ellK,hr,hp,hm},
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];

	r1g=p/(1-e);
	r2g=p/(1+e);
	4EEg+2EEg*r2g+EEg*r2g^2-(2(2a^2*EEg-4EEg*r2g+a*Lzg*r2g))/(a^2-2r2g+r2g^2)
]


(* ::Subsubsection::Closed:: *)
(*Mino-time azimuthal geodesic frequency*)


\[CapitalUpsilon]\[Phi]grfun[a_,p_,e_,xg_]:=Module[{EEg,Lzg,r1g,r2g,r3g,rp,rm,krg,ellK,hr,hp,hm},
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	ellK=EllipticK[krg];
	hr=(r1g-r2g)/(r1g-r3g);
	hp=hr (r3g-rp)/(r2g-rp);
	hm=hr (r3g-rm)/(r2g-rm);

	a/(rp-rm)((2EEg*rp-a*Lzg)/(r3g-rp)(1-(r2g-r3g)/(r2g-rp)EllipticPi[hp,krg]/ellK))-a/(rp-rm)((2EEg*rm-a*Lzg)/(r3g-rm)(1-(r2g-r3g)/(r2g-rm)EllipticPi[hm,krg]/ellK))
]


\[CapitalUpsilon]\[Phi]gzfun[a_,p_,e_,xg_]:=Lzgfun[a,p,e,xg]


\[CapitalUpsilon]\[Phi]gfunLim[a_,p_,e_,xg_]:=Module[{EEg,Lzg,r1g,r2g},
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];

	r1g=p/(1-e);
	r2g=p/(1+e);

	-((a(a*Lzg-2EEg*r2g))/(a^2-2r2g+r2g^2))+Lzg
]


(* ::Subsubsection::Closed:: *)
(*Mino-time spin precession frequency*)


\[CapitalUpsilon]pfun[a_,p_,e_,xg_]:=Module[{r1g,r2g,EEg,Lzg,r3g,krg,ellK,Lzred,sgn,\[Psi]freq},
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	ellK=EllipticK[krg];
	Lzred=Lzg-a*EEg;
	sgn=RealSign[Lzg-a*EEg];
	\[Psi]freq=1/Sqrt[(1-EEg^2)r2g(r1g-r3g)] ((2r3g^2(a*sgn+EEg*Abs[Lzred]))/(r3g^2+Lzred^2) ellK+I*Abs[Lzred](-a*sgn-EEg*Abs[Lzred])(r2g-r3g)( EllipticPi[((r1g-r2g)(r3g-I*Abs[Lzred]))/((r1g-r3g)(r2g-I*Abs[Lzred])),krg]/((r2g-I*Abs[Lzred])(r3g-I*Abs[Lzred]))- EllipticPi[((r1g-r2g)(r3g+I*Abs[Lzred]))/((r1g-r3g)(r2g+I*Abs[Lzred])),krg]/((r2g+I*Abs[Lzred])(r3g+I*Abs[Lzred]))));

	\[CapitalUpsilon]rgfun[a,p,e,xg]/(2\[Pi]) 2Re[\[Psi]freq]
]


\[CapitalUpsilon]pfunLim[a_,p_,e_,xg_]:=Module[{r2g,EEg,Lzg,Lzred},
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	r2g=p/(1+e);
	Lzred=Lzg-a*EEg;
	
	Lzred(a+EEg*Lzred)r2g^2/(Abs[Lzred](r2g^2+Lzred^2))
]


(* ::Subsubsection::Closed:: *)
(*BL radial frequency*)


\[CapitalOmega]rgfun[a_,p_,e_,xg_]:=\[CapitalUpsilon]rgfun[a,p,e,xg]/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])


(* ::Subsubsection::Closed:: *)
(*BL azimuthal geodesic frequency*)


\[CapitalOmega]\[Phi]gfun[a_,p_,e_,xg_]:=(\[CapitalUpsilon]\[Phi]grfun[a,p,e,xg]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,xg])/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])


(* ::Subsubsection::Closed:: *)
(*BL azimuthal BL time geodesic frequency - plunge orbits*)


\[CapitalOmega]\[Phi]gfunISCOplunge[r_,a_,x_]:=Module[{r1g,EEg,Lzg},
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,x];
	Lzg=Lzgfun[a,r1g,0,x];

	(2a*EEg+Lzg(-2+r))/(-2a*Lzg+EEg*r^3+a^2*EEg(2+r))
]


\[CapitalOmega]\[Phi]gfun[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];

	(2a*EEg+Lzg(-2+r))/(-2a*Lzg+EEg*r^3+a^2*EEg(2+r))
]


\[CapitalOmega]\[Phi]gfunPlunge[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg},
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];

	(2a*EEg+Lzg(-2+r))/(-2a*Lzg+EEg*r^3+a^2*EEg(2+r))
]


(* ::Subsubsection::Closed:: *)
(*BL spin precession frequency - plunge orbits*)


\[CapitalOmega]pfunISCOplunge[r_,a_,x_]:=Module[{r1g,EEg,Lzg,Lzred},
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,x];
	Lzg=Lzgfun[a,r1g,0,x];
	Lzred=Lzg-a*EEg;
	
	Lzred(a+EEg*Lzred)r^2/(Abs[Lzred](r^2+Lzred^2))*(a^2+(-2+r)r)/(r(-2a*Lzg+EEg*r^3+a^2*EEg(2+r)))
]


\[CapitalOmega]pfun[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,Lzred},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	Lzred=Lzg-a*EEg;
	
	Lzred(a+EEg*Lzred)r^2/(Abs[Lzred](r^2+Lzred^2))*(a^2+(-2+r)r)/(r(-2a*Lzg+EEg*r^3+a^2*EEg(2+r)))
]


\[CapitalOmega]pfunPlunge[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,Lzred},
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	Lzred=Lzg-a*EEg;
	
	Lzred(a+EEg*Lzred)r^2/(Abs[Lzred](r^2+Lzred^2))*(a^2+(-2+r)r)/(r(-2a*Lzg+EEg*r^3+a^2*EEg(2+r)))
]


(* ::Subsection::Closed:: *)
(*Derivatives geodesic frequencies wrt to the constants of motion*)


(* ::Subsubsection::Closed:: *)
(*Mino-time radial frequency*)


d\[CapitalUpsilon]rgdEgfun[a_,p_,e_,xg_]:=Module[{r1g,r2g,EEg,r3g,krg,ellK,ellE},
	r1g=p/(1-e);
	r2g=p/(1+e);
	EEg=EEgfun[a,p,e,xg];

	r3g=2/(1-EEg^2)-(r1g+r2g);
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	ellK=EllipticK[krg];
	ellE=EllipticE[krg];

	(EEg*\[Pi] Sqrt[(1-EEg^2)r2g(r1g-r3g)](-2ellE*r2g+ellK(r2g-r3g)(2-(1-EEg^2)r3g)))/(2(1-EEg^2)^2*ellK^2*(r2g-r3g)r3g)
]


(* ::Subsubsection::Closed:: *)
(*Mino-time coordinate-time frequency*)


d\[CapitalUpsilon]tgdEgfun[a_,p_,e_,xg_]:=Module[{r1g,r2g,rp,rm,EEg,Lzg,r3g,krg,hr,hp,hm,ellK,ellE,ellPirp,ellPirm,ellPi},
	r1g=p/(1-e);
	r2g=p/(1+e);
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];

	r3g=2/(1-EEg^2)-(r1g+r2g);
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	hr=(r1g-r2g)/(r1g-r3g);
	hp=hr (r3g-rp)/(r2g-rp);
	hm=hr (r3g-rm)/(r2g-rm);
	ellK=EllipticK[krg];
	ellPirp=EllipticPi[hp,krg];
	ellPirm=EllipticPi[hm,krg];
	ellPi=EllipticPi[hr,krg];
	ellK=EllipticK[krg];
	ellE=EllipticE[krg];

	(EEg^2r2g^2(-r1g+r3g))/((1-EEg^2)^2*ellK^2(r2g-r3g)r3g) ellE^2+(r2g(r1g-r3g)(-2EEg^2(-2+r3g)+r3g+EEg^4*r3g))/(2 (-1+EEg^2)^2* ellK*r3g) ellE+1/((-1+EEg^2)^3*ellK^2*r3g) (6EEg^2(ellE*r2g+ellK(r2g-r3g)(-1+r3g))+2EEg^6*ellK(r2g-r3g)r3g+3ellK*r3g(-r2g+r3g)+EEg^4(-4ellE*r2g-ellK(r2g-r3g)(-4+5r3g)))ellPi+4 /((-1+EEg^2)^2 ellK^2 r3g (r2g-rp) (r3g-rp)^2 (-rm+rp)) (a EEg Lzg rp (ellK (-r2g+r3g) rp+ellE r2g (-r3g+rp))+a^2 (ellK (r2g-r3g) r3g (r3g-rp)+EEg^4 ellK (r2g-r3g) r3g (r3g-rp)-2 EEg^2 (ellE r2g (r3g-rp)+ellK (r2g-r3g) (r3g^2+rp-r3g rp)))-2 rp (ellK (r2g-r3g) r3g (r3g-rp)+EEg^4 ellK (r2g-r3g) r3g (r3g-rp)-2 EEg^2 (ellE r2g (r3g-rp)+ellK (r2g-r3g) (r3g^2+rp-r3g rp))))ellPirp+4/((-1+EEg^2)^2 ellK^2 r3g (r2g-rm) (r3g-rm)^2 (rm-rp)) (a EEg Lzg rm (ellK (-r2g+r3g) rm+ellE r2g (-r3g+rm))+a^2 (ellK (r2g-r3g) r3g (r3g-rm)+EEg^4 ellK (r2g-r3g) r3g (r3g-rm)-2 EEg^2 (ellE r2g (r3g-rm)+ellK (r2g-r3g) (r3g^2+rm-r3g rm)))-2 rm (ellK (r2g-r3g) r3g (r3g-rm)+EEg^4 ellK (r2g-r3g) r3g (r3g-rm)-2 EEg^2 (ellE r2g (r3g-rm)+ellK (r2g-r3g) (r3g^2+rm-r3g rm))))ellPirm+4-(r1g r2g)/2+((-3+2 EEg^2) r3g)/(-1+EEg^2)+(EEg^2 (r1g (r2g-EEg^2 r2g)+2 r3g (-3-r3g+EEg^2 (2+r3g))))/((-1+EEg^2)^3 r3g)+(2 a Lzg rm)/(EEg (r3g-rm) (rm-rp))+(4 EEg (2 a^2 EEg-4 EEg rm+a Lzg rm))/((-1+EEg^2)^2 (r3g-rm)^2 (rm-rp))-(2 (2 a^2 EEg-4 EEg rm+a Lzg rm))/(EEg (r3g-rm) (rm-rp))+(2a*Lzg*rp)/(EEg(r3g-rp)(-rm+rp))+(4 EEg (2 a^2 EEg-4 EEg rp+a Lzg rp))/((-1+EEg^2)^2 (r3g-rp)^2 (-rm+rp))+(2 (-2 a^2+4 rp-(a Lzg rp)/EEg))/((r3g-rp) (-rm+rp))
]


d\[CapitalUpsilon]tgdLzgfun[a_,p_,e_,xg_]:=Module[{r1g,r2g,rp,rm,EEg,r3g,krg,hr,hp,hm,ellK,ellPirp,ellPirm},
	r1g=p/(1-e);
	r2g=p/(1+e);
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,xg];

	r3g=2/(1-EEg^2)-(r1g+r2g);
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	hr=(r1g-r2g)/(r1g-r3g);
	hp=hr (r3g-rp)/(r2g-rp);
	hm=hr (r3g-rm)/(r2g-rm);
	ellK=EllipticK[krg];
	ellPirp=EllipticPi[hp,krg];
	ellPirm=EllipticPi[hm,krg];

	(2a )/(rp-rm) (1/(r3g-rm) (1+(ellPirm(-r2g+r3g))/(ellK(r2g-rm))) rm-1/(r3g-rp) (1+(ellPirp(-r2g+r3g))/(ellK(r2g-rp))) rp)
]


(* ::Subsubsection::Closed:: *)
(*Mino-time azimuthal frequency*)


d\[CapitalUpsilon]\[Phi]gdEgfun[a_,p_,e_,xg_]:=Module[{r1g,r2g,rp,rm,EEg,Lzg,r3g,krg,hr,hp,hm,ellK,ellE,ellPirp,ellPirm,ellPi},
	r1g=p/(1-e);
	r2g=p/(1+e);
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];

	r3g=2/(1-EEg^2)-(r1g+r2g);
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	hr=(r1g-r2g)/(r1g-r3g);
	hp=hr (r3g-rp)/(r2g-rp);
	hm=hr (r3g-rm)/(r2g-rm);
	ellK=EllipticK[krg];
	ellPirp=EllipticPi[hp,krg];
	ellPirm=EllipticPi[hm,krg];
	ellPi=EllipticPi[hr,krg];
	ellK=EllipticK[krg];
	ellE=EllipticE[krg];

	((2a(a EEg Lzg (ellK (-r2g+r3g) rp+ellE r2g (-r3g+rp))+rp (-ellK (r2g-r3g) r3g (r3g-rp)-EEg^4 ellK (r2g-r3g) r3g (r3g-rp)+2 EEg^2 (ellE r2g (r3g-rp)+ellK (r2g-r3g) (r3g^2+rp-r3g rp)))))/((-1+EEg^2)^2 ellK^2 r3g (r2g-rp) (r3g-rp)^2 (-rm+rp)))ellPirp+((2 a (a EEg Lzg (ellK (-r2g+r3g) rm+ellE r2g (-r3g+rm))+rm (-ellK (r2g-r3g) r3g (r3g-rm)-EEg^4 ellK (r2g-r3g) r3g (r3g-rm)+2 EEg^2 (ellE r2g (r3g-rm)+ellK (r2g-r3g) (r3g^2+rm-r3g rm)))))/((-1+EEg^2)^2 ellK^2 r3g (r2g-rm) (r3g-rm)^2 (rm-rp)))ellPirm+(2 a (r3g (r3g-rm) (r3g-rp)+EEg^4 r3g (r3g-rm) (r3g-rp)+a EEg Lzg (2 r3g-rm-rp)+2 EEg^2 (r3g^2 (-1-r3g+rm)+(r3g^2+rm-r3g rm) rp)))/((-1+EEg^2)^2 (r3g-rm)^2 (r3g-rp)^2)
]


d\[CapitalUpsilon]\[Phi]gdLzgfun[a_,p_,e_,xg_]:=Module[{r1g,r2g,rp,rm,EEg,r3g,krg,hr,hp,hm,ellK,ellPirp,ellPirm},
	r1g=p/(1-e);
	r2g=p/(1+e);
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,xg];

	r3g=2/(1-EEg^2)-(r1g+r2g);
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	hr=(r1g-r2g)/(r1g-r3g);
	hp=hr (r3g-rp)/(r2g-rp);
	hm=hr (r3g-rm)/(r2g-rm);
	ellK=EllipticK[krg];
	ellPirp=EllipticPi[hp,krg];
	ellPirm=EllipticPi[hm,krg];
	
	1+a^2/((r3g-rm)(-rm+rp)) (1-(ellPirm(r2g-r3g))/(ellK(r2g-rm)))-a^2/((r3g-rp)(-rm+rp)) (1-(ellPirp(r2g-r3g))/(ellK(r2g-rp)))
]


(* ::Subsubsection::Closed:: *)
(*BL time radial frequency*)


d\[CapitalOmega]rgdEgfun[a_,p_,e_,xg_]:=d\[CapitalUpsilon]rgdEgfun[a,p,e,xg]/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])-(\[CapitalOmega]rgfun[a,p,e,xg] d\[CapitalUpsilon]tgdEgfun[a,p,e,xg])/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])


d\[CapitalOmega]rgdLzgfun[a_,p_,e_,xg_]:=0-(\[CapitalOmega]rgfun[a,p,e,xg] d\[CapitalUpsilon]tgdLzgfun[a,p,e,xg])/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])


(* ::Subsubsection::Closed:: *)
(*BL time azimuthal frequency*)


d\[CapitalOmega]\[Phi]gdEgfun[a_,p_,e_,xg_]:=d\[CapitalUpsilon]\[Phi]gdEgfun[a,p,e,xg]/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])-(\[CapitalOmega]\[Phi]gfun[a,p,e,xg] d\[CapitalUpsilon]tgdEgfun[a,p,e,xg])/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])


d\[CapitalOmega]\[Phi]gdLzgfun[a_,p_,e_,xg_]:=d\[CapitalUpsilon]\[Phi]gdLzgfun[a,p,e,xg]/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])-(\[CapitalOmega]\[Phi]gfun[a,p,e,xg] d\[CapitalUpsilon]tgdLzgfun[a,p,e,xg])/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])


(* ::Subsection::Closed:: *)
(*Derivatives geodesic frequencies wrt to geodesics radial roots*)


(* ::Subsubsection::Closed:: *)
(*Mino-time radial frequency*)


d\[CapitalUpsilon]rgdr1gfun[a_,p_,e_,xg_]:=Module[{r1g,r2g,EEg,r3g,krg,ellK,ellE},
	r1g=p/(1-e);
	r2g=p/(1+e);
	EEg=EEgfun[a,p,e,xg];

	r3g=2/(1-EEg^2)-(r1g+r2g);
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	ellK=EllipticK[krg];
	ellE=EllipticE[krg];
	
	d\[CapitalUpsilon]rgdEgfun[a,p,e,xg]dEgdr1gfun[a,p,e,xg]+(\[Pi] Sqrt[(1-EEg^2)r2g(r1g-r3g)](ellK*r1g(r2g-r3g)(-r1g+r2g+r3g)+ellE*r2g(r1g^2-r1g*r2g+r3g(-r2g+r3g))))/(4ellK^2*r1g(r1g-r2g)(r2g-r3g)r3g)
]


d\[CapitalUpsilon]rgdr2gfun[a_,p_,e_,xg_]:=Module[{r1g,r2g,EEg,r3g,krg,ellK,ellE},
	r1g=p/(1-e);
	r2g=p/(1+e);
	EEg=EEgfun[a,p,e,xg];

	r3g=2/(1-EEg^2)-(r1g+r2g);
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	ellK=EllipticK[krg];
	ellE=EllipticE[krg];
	
	d\[CapitalUpsilon]rgdEgfun[a,p,e,xg]dEgdr2gfun[a,p,e,xg]+(\[Pi] Sqrt[(1-EEg^2)r2g(r1g-r3g)](-ellK(r2g-r3g)(r1g-r2g+r3g)+ellE(-r2g^2-r3g^2+r1g(r2g+r3g))))/(4ellK^2*(r1g-r2g)(r2g-r3g)r3g)
]


(* ::Subsubsection::Closed:: *)
(*Mino-time coordinate-time frequency*)


d\[CapitalUpsilon]tgdr1gfun[a_,p_,e_,xg_]:=Module[{r1g,r2g,rp,rm,EEg,Lzg,r3g,krg,hr,hp,hm,ellK,ellE,ellPirp,ellPirm,ellPi},
	r1g=p/(1-e);
	r2g=p/(1+e);
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];

	r3g=2/(1-EEg^2)-(r1g+r2g);
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	hr=(r1g-r2g)/(r1g-r3g);
	hp=hr (r3g-rp)/(r2g-rp);
	hm=hr (r3g-rm)/(r2g-rm);
	ellK=EllipticK[krg];
	ellPirp=EllipticPi[hp,krg];
	ellPirm=EllipticPi[hm,krg];
	ellPi=EllipticPi[hr,krg];
	ellK=EllipticK[krg];
	ellE=EllipticE[krg];
d\[CapitalUpsilon]tgdEgfun[a,p,e,xg]dEgdr1gfun[a,p,e,xg]+d\[CapitalUpsilon]tgdLzgfun[a,p,e,xg]dLzgdr1gfun[a,p,e,xg]+(EEg r2g^2 (r1g-r3g) (r1g^2-r1g r2g+r3g (-r2g+r3g)))/(4 ellK^2 r1g (r1g-r2g) (r2g-r3g) r3g) ellE^2+(-((r2g (r1g-r3g) (a^2 EEg (-1+EEg^2) r1g (a^2+(-2+r1g) r1g) (r1g-r2g)+(2 a^3 (-1+EEg^2) Lzg+a^4 EEg (3+r1g-EEg^2 (2+r1g))+2 EEg (-1+EEg^2) (-2+r1g) r1g^2 (-r1g+r2g)+a^2 EEg r1g (r1g (3+r1g-EEg^2 (2+r1g))+2 EEg^2 r2g-2 (1+r2g))) r3g+(-2 a (-1+EEg^2) Lzg r1g+a^2 EEg (-2+(-1+EEg^2) r1g (2+r1g-r2g))+EEg r1g (4-2 r1g-r1g^3+EEg^2 r1g^3-(-1+EEg^2) (-2+r1g) r1g r2g)) r3g^2-EEg (a^2+(-2+r1g) r1g) (-3-r1g+EEg^2 (2+r1g)) r3g^3))/(2 (-1+EEg^2) ellK r1g (a^2+(-2+r1g) r1g) (r1g-r2g) r3g (a^2+(-2+r3g) r3g))))ellE+(EEg (-3+2 EEg^2) (ellK r1g (r2g-r3g) (-r1g+r2g+r3g)+ellE r2g (r1g^2-r1g r2g+r3g (-r2g+r3g))))/(2 (-1+EEg^2) ellK^2 r1g (r1g-r2g) r3g) ellPi+(((2 a^2 EEg-4 EEg rp+a Lzg rp) (ellE r2g (r1g^2-r1g r2g+r3g (-r2g+r3g)) (r1g-rp) (r3g-rp)+ellK r1g (r2g-r3g) ((r1g-r3g) (r1g+r3g-rp) rp+r2g (r3g^2-r3g rp+rp (-r1g+rp)))))/(ellK^2 r1g (r1g-r2g) r3g (r1g-rp) (r2g-rp) (r3g-rp)^2 (-rm+rp)))ellPirp+(((2 a^2 EEg-4 EEg rm+a Lzg rm) (ellE r2g (r1g^2-r1g r2g+r3g (-r2g+r3g)) (r1g-rm) (r3g-rm)+ellK r1g (r2g-r3g) ((r1g-r3g) (r1g+r3g-rm) rm+r2g (r3g^2-r3g rm+rm (-r1g+rm)))))/(ellK^2 r1g (r1g-r2g) r3g (r1g-rm) (r2g-rm) (r3g-rm)^2 (rm-rp)))ellPirm+1/4 (-8 EEg+(4 EEg)/(-1+EEg^2)-2 EEg r2g+(2 EEg (-3+2 EEg^2) (r1g-2 r2g+r3g))/((-1+EEg^2) (r1g-r2g))+(EEg r2g (r1g^2-r1g r2g+r3g (-r2g+r3g)))/((r1g-r2g) r3g)-(8 (2 a^2 EEg-4 EEg rm+a Lzg rm))/((r3g-rm)^2 (rm-rp))-(4 (2 a^2 EEg-4 EEg rm+a Lzg rm) (-r1g^2+r2g (r3g-2 rm)+r1g (r2g+rm)+r3g (-r3g+rm)))/((r1g-r2g) (r1g-rm) (r3g-rm)^2 (rm-rp))-(8 (2 a^2 EEg-4 EEg rp+a Lzg rp))/((r3g-rp)^2 (-rm+rp))-(4 (2 a^2 EEg-4 EEg rp+a Lzg rp) (-r1g^2+r2g (r3g-2 rp)+r1g (r2g+rp)+r3g (-r3g+rp)))/((r1g-r2g) (r1g-rp) (r3g-rp)^2 (-rm+rp)))
]


d\[CapitalUpsilon]tgdr2gfun[a_,p_,e_,xg_]:=Module[{r1g,r2g,rp,rm,EEg,Lzg,r3g,krg,hr,hp,hm,ellK,ellE,ellPirp,ellPirm,ellPi},
	r1g=p/(1-e);
	r2g=p/(1+e);
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];

	r3g=2/(1-EEg^2)-(r1g+r2g);
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	hr=(r1g-r2g)/(r1g-r3g);
	hp=hr (r3g-rp)/(r2g-rp);
	hm=hr (r3g-rm)/(r2g-rm);
	ellK=EllipticK[krg];
	ellPirp=EllipticPi[hp,krg];
	ellPirm=EllipticPi[hm,krg];
	ellPi=EllipticPi[hr,krg];
	ellK=EllipticK[krg];
	ellE=EllipticE[krg];
	d\[CapitalUpsilon]tgdEgfun[a,p,e,xg]dEgdr2gfun[a,p,e,xg]+d\[CapitalUpsilon]tgdLzgfun[a,p,e,xg]dLzgdr2gfun[a,p,e,xg]-(EEg r2g (r1g-r3g) (r2g^2+r3g^2-r1g (r2g+r3g)))/(4 ellK^2 (r1g-r2g) (r2g-r3g) r3g) ellE^2+(-(((r1g-r3g) (a^2 EEg (-1+EEg^2) (r1g-r2g) r2g (a^2+(-2+r2g) r2g)+(-2 a^3 (-1+EEg^2) Lzg-2 EEg (-1+EEg^2) (r1g-r2g) (-2+r2g) r2g^2+a^4 EEg (-3-r2g+EEg^2 (2+r2g))+a^2 EEg r2g (2-2 (-1+EEg^2) r1g+r2g (-3-r2g+EEg^2 (2+r2g)))) r3g+(2 a (-1+EEg^2) Lzg r2g+a^2 EEg (2+(-1+EEg^2) (-2+r1g-r2g) r2g)+EEg r2g (-4+r2g (2+(-1+EEg^2) (r1g (-2+r2g)-r2g^2)))) r3g^2+EEg (a^2+(-2+r2g) r2g) (-3-r2g+EEg^2 (2+r2g)) r3g^3))/(2 (-1+EEg^2) ellK (r1g-r2g) (a^2+(-2+r2g) r2g) r3g (a^2+(-2+r3g) r3g))))ellE+(EEg (-3+2 EEg^2) (-ellK (r2g-r3g) (r1g-r2g+r3g)+ellE (-r2g^2-r3g^2+r1g (r2g+r3g))))/(2 (-1+EEg^2) ellK^2 (r1g-r2g) r3g) ellPi+(((2 a^2 EEg-4 EEg rp+a Lzg rp) (ellE (-r2g^2-r3g^2+r1g (r2g+r3g)) (r2g-rp) (r3g-rp)+ellK (r2g-r3g) (-((r2g-r3g) (r2g+r3g-rp) rp)-r1g (r3g^2-r3g rp+rp (-r2g+rp)))))/(ellK^2 (r1g-r2g) r3g (r2g-rp)^2 (r3g-rp)^2 (-rm+rp)))ellPirp+(((2 a^2 EEg-4 EEg rm+a Lzg rm) (ellE (-r2g^2-r3g^2+r1g (r2g+r3g)) (r2g-rm) (r3g-rm)+ellK (r2g-r3g) (-((r2g-r3g) (r2g+r3g-rm) rm)-r1g (r3g^2-r3g rm+rm (-r2g+rm)))))/(ellK^2 (r1g-r2g) r3g (r2g-rm)^2 (r3g-rm)^2 (rm-rp)))ellPirm+1/4 (-8 EEg+(4 EEg)/(-1+EEg^2)-2 EEg r1g+(2 EEg (-3+2 EEg^2) (2 r1g-r2g-r3g))/((-1+EEg^2) (r1g-r2g))+(EEg r1g (-r2g^2-r3g^2+r1g (r2g+r3g)))/((r1g-r2g) r3g)-(8 (2 a^2 EEg-4 EEg rm+a Lzg rm))/((r3g-rm)^2 (rm-rp))+(4 (2 a^2 EEg-4 EEg rm+a Lzg rm) (-r2g^2+r1g (r2g+r3g-2 rm)+r2g rm+r3g (-r3g+rm)))/((r1g-r2g) (r2g-rm) (r3g-rm)^2 (rm-rp))-(8 (2 a^2 EEg-4 EEg rp+a Lzg rp))/((r3g-rp)^2 (-rm+rp))+(4 (2 a^2 EEg-4 EEg rp+a Lzg rp) (-r2g^2+r1g (r2g+r3g-2 rp)+r2g rp+r3g (-r3g+rp)))/((r1g-r2g) (r2g-rp) (r3g-rp)^2 (-rm+rp)))
]


(* ::Subsubsection::Closed:: *)
(*Mino-time azimuthal frequency*)


d\[CapitalUpsilon]\[Phi]gdr1gfun[a_,p_,e_,xg_]:=Module[{r1g,r2g,rp,rm,EEg,Lzg,r3g,krg,hr,hp,hm,ellK,ellE,ellPirp,ellPirm,ellPi},
	r1g=p/(1-e);
	r2g=p/(1+e);
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];

	r3g=2/(1-EEg^2)-(r1g+r2g);
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	hr=(r1g-r2g)/(r1g-r3g);
	hp=hr (r3g-rp)/(r2g-rp);
	hm=hr (r3g-rm)/(r2g-rm);
	ellK=EllipticK[krg];
	ellPirp=EllipticPi[hp,krg];
	ellPirm=EllipticPi[hm,krg];
	ellPi=EllipticPi[hr,krg];
	ellK=EllipticK[krg];
	ellE=EllipticE[krg];

	(d\[CapitalUpsilon]\[Phi]gdEgfun[a,p,e,xg]dEgdr1gfun[a,p,e,xg]+d\[CapitalUpsilon]\[Phi]gdLzgfun[a,p,e,xg]dLzgdr1gfun[a,p,e,xg])+((a (a Lzg-2 EEg rp) (ellE r2g (r1g^2-r1g r2g+r3g (-r2g+r3g)) (r1g-rp) (r3g-rp)+ellK r1g (r2g-r3g) ((r1g-r3g) (r1g+r3g-rp) rp+r2g (r3g^2-r3g rp+rp (-r1g+rp)))))/(2 ellK^2 r1g (r1g-r2g) r3g (r1g-rp) (r2g-rp) (r3g-rp)^2 (-rm+rp)))ellPirp+((a (a Lzg-2 EEg rm) (ellE r2g (r1g^2-r1g r2g+r3g (-r2g+r3g)) (r1g-rm) (r3g-rm)+ellK r1g (r2g-r3g) ((r1g-r3g) (r1g+r3g-rm) rm+r2g (r3g^2-r3g rm+rm (-r1g+rm)))))/(2 ellK^2 r1g (r1g-r2g) r3g (r1g-rm) (r2g-rm) (r3g-rm)^2 (rm-rp)))ellPirm+(a r2g (r1g-r3g) (a Lzg (-r1g-r3g+rm+rp)+2 EEg (r1g r3g-rm rp)))/(2 ellK r1g (r1g-r2g) (r1g-rm) (-r3g+rm) (r1g-rp) (r3g-rp)) ellE-1/(2 (r1g-r2g) (r1g-rm) (r3g-rm)^2 (r1g-rp) (r3g-rp)^2) a (r1g-r3g) (a Lzg (-((r2g-r3g) (r3g-rm)^2)+r1g^2 (2 r3g-rm-rp)+(2 r3g-rm) (r2g-r3g+rm) rp-(r2g-r3g+rm) rp^2-r1g (2 r3g-rm-rp) (r2g-r3g+rm+rp))+2 EEg (r1g (r2g-r3g+rm+rp) (r3g^2-rm rp)+r1g^2 (-r3g^2+rm rp)+rm rp ((r3g-rm) (r3g-rp)+r2g (-2 r3g+rm+rp))))
]


d\[CapitalUpsilon]\[Phi]gdr2gfun[a_,p_,e_,xg_]:=Module[{r1g,r2g,rp,rm,EEg,Lzg,r3g,krg,hr,hp,hm,ellK,ellE,ellPirp,ellPirm,ellPi},
	r1g=p/(1-e);
	r2g=p/(1+e);
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];

	r3g=2/(1-EEg^2)-(r1g+r2g);
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	hr=(r1g-r2g)/(r1g-r3g);
	hp=hr (r3g-rp)/(r2g-rp);
	hm=hr (r3g-rm)/(r2g-rm);
	ellK=EllipticK[krg];
	ellPirp=EllipticPi[hp,krg];
	ellPirm=EllipticPi[hm,krg];
	ellPi=EllipticPi[hr,krg];
	ellK=EllipticK[krg];
	ellE=EllipticE[krg];
	(d\[CapitalUpsilon]\[Phi]gdEgfun[a,p,e,xg]dEgdr2gfun[a,p,e,xg]+d\[CapitalUpsilon]\[Phi]gdLzgfun[a,p,e,xg]dLzgdr2gfun[a,p,e,xg])+((a (a Lzg-2 EEg rp) (ellE (-r2g^2-r3g^2+r1g (r2g+r3g)) (r2g-rp) (r3g-rp)+ellK (r2g-r3g) (-((r2g-r3g) (r2g+r3g-rp) rp)-r1g (r3g^2-r3g rp+rp (-r2g+rp)))))/(2 ellK^2 (r1g-r2g) r3g (r2g-rp)^2 (r3g-rp)^2 (-rm+rp)))ellPirp+((a (a Lzg-2 EEg rm) (ellE (-r2g^2-r3g^2+r1g (r2g+r3g)) (r2g-rm) (r3g-rm)+ellK (r2g-r3g) (-((r2g-r3g) (r2g+r3g-rm) rm)-r1g (r3g^2-r3g rm+rm (-r2g+rm)))))/(2 ellK^2 (r1g-r2g) r3g (r2g-rm)^2 (r3g-rm)^2 (rm-rp)))ellPirm+(a (r1g-r3g) (a Lzg (r2g+r3g-rm-rp)+EEg (-2 r2g r3g+2 rm rp)))/(2 ellK (r1g-r2g) (r2g-rm) (-r3g+rm) (r2g-rp) (r3g-rp)) ellE+1/(2 (r1g-r2g) (r2g-rm) (r3g-rm)^2 (r2g-rp) (r3g-rp)^2) a (r2g-r3g) (2 EEg (r2g r3g^2 (r1g-r2g-r3g+rm)+(r2g r3g^2+(r2g^2+r2g r3g+r3g^2-r1g (r2g+2 r3g)) rm-(-r1g+r2g+r3g) rm^2) rp+rm (r1g-r2g-r3g+rm) rp^2)+a Lzg ((r3g-rm) (r3g-rp) (r3g-rm-rp)+r2g^2 (2 r3g-rm-rp)+r2g (r3g-rm-rp) (2 r3g-rm-rp)-r1g (r3g^2+rm^2+r2g (2 r3g-rm-rp)+rm rp+rp^2-2 r3g (rm+rp))))
]


(* ::Subsubsection::Closed:: *)
(*BL time radial frequency*)


d\[CapitalOmega]rgdr1gfun[a_,p_,e_,xg_]:=d\[CapitalUpsilon]rgdr1gfun[a,p,e,xg]/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])-(\[CapitalOmega]rgfun[a,p,e,xg] d\[CapitalUpsilon]tgdr1gfun[a,p,e,xg])/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])


d\[CapitalOmega]rgdr2gfun[a_,p_,e_,xg_]:=d\[CapitalUpsilon]rgdr2gfun[a,p,e,xg]/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])-(\[CapitalOmega]rgfun[a,p,e,xg] d\[CapitalUpsilon]tgdr2gfun[a,p,e,xg])/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])


(* ::Subsubsection::Closed:: *)
(*BL time azimuthal frequency*)


d\[CapitalOmega]\[Phi]gdr1gfun[a_,p_,e_,xg_]:=d\[CapitalUpsilon]\[Phi]gdr1gfun[a,p,e,xg]/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])-(\[CapitalOmega]\[Phi]gfun[a,p,e,xg] d\[CapitalUpsilon]tgdr1gfun[a,p,e,xg])/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])


d\[CapitalOmega]\[Phi]gdr2gfun[a_,p_,e_,xg_]:=d\[CapitalUpsilon]\[Phi]gdr2gfun[a,p,e,xg]/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])-(\[CapitalOmega]\[Phi]gfun[a,p,e,xg] d\[CapitalUpsilon]tgdr2gfun[a,p,e,xg])/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])


(* ::Subsection::Closed:: *)
(*Derivatives geodesic frequencies wrt to (p,e)*)


(* ::Subsubsection::Closed:: *)
(*Mino-time radial frequency*)


d\[CapitalUpsilon]rgdpfun[a_,p_,e_,xg_]:=d\[CapitalUpsilon]rgdr1gfun[a,p,e,xg] 1/(1-e)+d\[CapitalUpsilon]rgdr2gfun[a,p,e,xg] 1/(1+e)


d\[CapitalUpsilon]rgdefun[a_,p_,e_,xg_]:=p(d\[CapitalUpsilon]rgdr1gfun[a,p,e,xg] 1/(1-e)^2-d\[CapitalUpsilon]rgdr2gfun[a,p,e,xg] 1/(1+e)^2)


(* ::Subsubsection::Closed:: *)
(*Mino-time coordinate-time frequency*)


d\[CapitalUpsilon]tgdpfun[a_,p_,e_,xg_]:=d\[CapitalUpsilon]tgdr1gfun[a,p,e,xg] 1/(1-e)+d\[CapitalUpsilon]tgdr2gfun[a,p,e,xg] 1/(1+e)


d\[CapitalUpsilon]tgdefun[a_,p_,e_,xg_]:=p(d\[CapitalUpsilon]tgdr1gfun[a,p,e,xg] 1/(1-e)^2-d\[CapitalUpsilon]tgdr2gfun[a,p,e,xg] 1/(1+e)^2)


(* ::Subsubsection::Closed:: *)
(*Mino-time azimuthal frequency*)


d\[CapitalUpsilon]\[Phi]gdpfun[a_,p_,e_,xg_]:=d\[CapitalUpsilon]\[Phi]gdr1gfun[a,p,e,xg] 1/(1-e)+d\[CapitalUpsilon]\[Phi]gdr2gfun[a,p,e,xg] 1/(1+e)


d\[CapitalUpsilon]\[Phi]gdefun[a_,p_,e_,xg_]:=p(d\[CapitalUpsilon]\[Phi]gdr1gfun[a,p,e,xg] 1/(1-e)^2-d\[CapitalUpsilon]\[Phi]gdr2gfun[a,p,e,xg] 1/(1+e)^2)


(* ::Subsubsection::Closed:: *)
(*BL time radial frequency*)


d\[CapitalOmega]rgdpfun[a_,p_,e_,xg_]:=d\[CapitalUpsilon]rgdpfun[a,p,e,xg]/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])-(\[CapitalOmega]rgfun[a,p,e,xg] d\[CapitalUpsilon]tgdpfun[a,p,e,xg])/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])


d\[CapitalOmega]rgdefun[a_,p_,e_,xg_]:=d\[CapitalUpsilon]rgdefun[a,p,e,xg]/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])-(\[CapitalOmega]rgfun[a,p,e,xg] d\[CapitalUpsilon]tgdefun[a,p,e,xg])/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])


(* ::Subsubsection::Closed:: *)
(*BL time azimuthal frequency*)


d\[CapitalOmega]\[Phi]gdpfun[a_,p_,e_,xg_]:=d\[CapitalUpsilon]\[Phi]gdpfun[a,p,e,xg]/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])-(\[CapitalOmega]\[Phi]gfun[a,p,e,xg] d\[CapitalUpsilon]tgdpfun[a,p,e,xg])/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])


d\[CapitalOmega]\[Phi]gdefun[a_,p_,e_,xg_]:=d\[CapitalUpsilon]\[Phi]gdefun[a,p,e,xg]/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])-(\[CapitalOmega]\[Phi]gfun[a,p,e,xg] d\[CapitalUpsilon]tgdefun[a,p,e,xg])/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])


(* ::Subsection::Closed:: *)
(*Geodesics periodic trajectory*)


(* ::Subsubsection::Closed:: *)
(*Geodesic radial trajectory*)


rgICr2gfun[wr_,a_,p_,e_,xg_]:=Module[{EEg,r1g,r2g,r3g,krg,ellK,jSN},
	EEg=EEgfun[a,p,e,xg];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	ellK=EllipticK[krg];
	jSN=JacobiSN[1/\[Pi]*ellK*wr,krg];

	(r3g(r1g-r2g)jSN^2-r2g(r1g-r3g))/((r1g-r2g)jSN^2-(r1g-r3g))
]


rgICr1gfun[wr_,a_,p_,e_,x_]:=Module[{EEg,r1g,r2g,r3g,krg,ellK,jSN},
	EEg=EEgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	ellK=EllipticK[krg];
	jSN=JacobiSN[1/\[Pi]*ellK*wr,krg];

	r1g*r2g/(r2g+(r1g-r2g)*jSN^2)
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory*)


tgfunref[a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1g,r2g,r3g,krg,\[Gamma]r,ellK,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over,\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Gamma]r=1-r1g/r2g;
	ellK=2EllipticK[krg];
	ellE=2EllipticE[krg];
	ellPi=2EllipticPi[\[Gamma]r,krg];

	\[ScriptCapitalI]=(2ellK)/Sqrt[r2g(r1g-r3g)];
	\[ScriptCapitalI]rover=2/Sqrt[r2g(r1g-r3g)]r1g*ellPi;
	\[ScriptCapitalI]r2over=-Sqrt[r2g(r1g-r3g)](-ellE+r1g*ellK/(r1g-r3g)-r1g(r1g+r2g+r3g)ellPi/(r2g(r1g-r3g)));
	\[ScriptCapitalI]p=2/Sqrt[(r1g-r3g)r2g] 1/rp (r1g/(r1g-rp) 2EllipticPi[(rp(r1g-r2g))/((r1g-rp)r2g),krg]- ellK);
	\[ScriptCapitalI]mreg=2/Sqrt[(r1g-r3g)r2g] (r1g/(r1g-rm) 2EllipticPi[(rm(r1g-r2g))/((r1g-rm)r2g),krg]- ellK);

	1/(2\[Pi]) 1/Sqrt[1-EEg^2] (4EEg*\[ScriptCapitalI]+2EEg*\[ScriptCapitalI]rover+EEg*\[ScriptCapitalI]r2over-4EEg rp/(rp-rm) (rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)-2(-4EEg+a Lzg) 1/(rp-rm) (rp \[ScriptCapitalI]p-\[ScriptCapitalI]mreg))
]


(* ::Text:: *)
(*Purely oscillatory part of the geodesic coordinate time trajectory*)


tgICr2gfun[wr_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1g,r2g,r3g,krg,krghold,\[Gamma]r,\[Phi],ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over,\[ScriptCapitalI]p,\[ScriptCapitalI]m},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Gamma]r=(r1g-r2g)/(r1g-r3g);
	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result. Having \[Phi]=\[Pi]/2 with MachinePrecision cause some problems with the integrals \[ScriptCapitalI]p and \[ScriptCapitalI]m*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	ellPi=EllipticPi[\[Gamma]r,\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[r2g(r1g-r3g)];
	\[ScriptCapitalI]rover=2/Sqrt[r2g(r1g-r3g)]((r2g-r3g)ellPi+r3g ellF);
	\[ScriptCapitalI]r2over=2r2g^2/Sqrt[r2g(r1g-r3g)]((r3g/r2g)^2*ellF+2 r3g/r2g (1-r3g/r2g)ellPi+(1-r3g/r2g)^2 1/(2(\[Gamma]r-1)(krg-\[Gamma]r)) (\[Gamma]r*ellE+(krg-\[Gamma]r)ellF+(2\[Gamma]r*krg+2\[Gamma]r-\[Gamma]r^2-3krg)ellPi-(\[Gamma]r^2Sin[\[Phi]]Cos[\[Phi]]Sqrt[1-krg Sin[\[Phi]]^2])/(1-\[Gamma]r*Sin[\[Phi]]^2)));

	\[ScriptCapitalI]p=2/(Sqrt[r2g(r1g-r3g)](r3g-rp))(ellF-(r2g-r3g)/(r2g-rp)EllipticPi[(r3g-rp)(r1g-r2g)/((r1g-r3g)(r2g-rp)),\[Phi],krg]);
	\[ScriptCapitalI]m=2/(Sqrt[r2g(r1g-r3g)](r3g-rm))(ellF-(r2g-r3g)/(r2g-rm)EllipticPi[(r3g-rm)(r1g-r2g)/((r1g-r3g)(r2g-rm)),\[Phi],krg]);

	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	1/Sqrt[1-EEg^2](4EEg*\[ScriptCapitalI]+2EEg*\[ScriptCapitalI]rover+EEg*\[ScriptCapitalI]r2over-4a^2*EEg 1/(rp-rm)(\[ScriptCapitalI]p-\[ScriptCapitalI]m)-2(-4EEg+a*Lzg)1/(rp-rm)(rp*\[ScriptCapitalI]p-rm*\[ScriptCapitalI]m))-wr*tgfunref[a,p,e,x]
]


tgICr1gfun[wr_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1g,r2g,r3g,krg,krghold,\[Gamma]r,\[Phi],ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over,\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Gamma]r=1-r1g/r2g;
	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result. Having \[Phi]=\[Pi]/2 with MachinePrecision cause some problems with the integrals \[ScriptCapitalI]p and \[ScriptCapitalI]m*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	ellPi=EllipticPi[\[Gamma]r,\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[r2g(r1g-r3g)];
	\[ScriptCapitalI]rover=2/Sqrt[r2g(r1g-r3g)]r1g*ellPi;
	\[ScriptCapitalI]r2over=-Sqrt[r2g(r1g-r3g)](-ellE+r1g*ellF/(r1g-r3g)-r1g(r1g+r2g+r3g)ellPi/(r2g(r1g-r3g))+((r1g-r2g)Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(-r2g+(-r1g+r2g)Sin[\[Phi]]^2));
	\[ScriptCapitalI]p=2/Sqrt[(r1g-r3g)r2g] 1/rp (r1g/(r1g-rp) EllipticPi[(rp(r1g-r2g))/((r1g-rp)r2g),\[Phi],krg]- ellF);
	\[ScriptCapitalI]mreg=2/Sqrt[(r1g-r3g)r2g] (r1g/(r1g-rm) EllipticPi[(rm(r1g-r2g))/((r1g-rm)r2g),\[Phi],krg]- ellF);
	
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	1/Sqrt[1-EEg^2](4EEg*\[ScriptCapitalI]+2EEg*\[ScriptCapitalI]rover+EEg*\[ScriptCapitalI]r2over-4rp*EEg 1/(rp-rm)(rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)-2(-4EEg+a*Lzg)1/(rp-rm)(rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg))-wr*tgfunref[a,p,e,x]
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal trajectory*)


\[Phi]gfunref[a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1g,r2g,r3g,krg,\[Gamma]r,ellK,\[ScriptCapitalI],\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Gamma]r=1-r1g/r2g;
	ellK=2EllipticK[krg];
	\[ScriptCapitalI]=(2 ellK)/Sqrt[r2g (r1g-r3g)];
	\[ScriptCapitalI]p=2/Sqrt[(r1g-r3g)r2g] 1/rp (r1g/(r1g-rp) 2EllipticPi[(rp(r1g-r2g))/((r1g-rp)r2g),krg]- ellK);
	\[ScriptCapitalI]mreg=2/Sqrt[(r1g-r3g)r2g] (r1g/(r1g-rm) 2EllipticPi[(rm(r1g-r2g))/((r1g-rm)r2g),krg]- ellK);

	1/(2\[Pi]) 1/Sqrt[1-EEg^2] (Lzg*\[ScriptCapitalI]-Lzg rp/(rp-rm) (rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)+2a*EEg 1/(rp-rm) (rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg))
]


(* ::Text:: *)
(*Purely oscillatory part of the geodesic azimuthal trajectory*)


\[Phi]gICr2gfun[wr_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1g,r2g,r3g,krg,krghold,\[Phi],ellF,\[ScriptCapitalI],\[ScriptCapitalI]p,\[ScriptCapitalI]m},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result. Having \[Phi]=\[Pi]/2 with MachinePrecision cause some problems with the integrals \[ScriptCapitalI]p and \[ScriptCapitalI]m*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	ellF=EllipticF[\[Phi],krg];

	\[ScriptCapitalI]=2ellF/Sqrt[r2g(r1g-r3g)];
	\[ScriptCapitalI]p=2/(Sqrt[r2g(r1g-r3g)](r3g-rp))(ellF-(r2g-r3g)/(r2g-rp)EllipticPi[(r3g-rp)(r1g-r2g)/((r1g-r3g)(r2g-rp)),\[Phi],krg]);
	\[ScriptCapitalI]m=2/(Sqrt[r2g(r1g-r3g)](r3g-rm))(ellF-(r2g-r3g)/(r2g-rm)EllipticPi[(r3g-rm)(r1g-r2g)/((r1g-r3g)(r2g-rm)),\[Phi],krg]);

	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	1/Sqrt[1-EEg^2](Lzg*\[ScriptCapitalI]-a^2*Lzg 1/(rp-rm)(\[ScriptCapitalI]p-\[ScriptCapitalI]m)+2a*EEg 1/(rp-rm)(rp*\[ScriptCapitalI]p-rm*\[ScriptCapitalI]m))-wr*\[Phi]gfunref[a,p,e,x]
]


\[Phi]gICr1gfun[wr_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1g,r2g,r3g,krg,krghold,\[Gamma]r,\[Phi],ellF,\[ScriptCapitalI],\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Gamma]r=1-r1g/r2g;
	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result. Having \[Phi]=\[Pi]/2 with MachinePrecision cause some problems with the integrals \[ScriptCapitalI]p and \[ScriptCapitalI]m*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	ellF=EllipticF[\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[r2g(r1g-r3g)];
	\[ScriptCapitalI]p=2/Sqrt[(r1g-r3g)r2g] 1/rp (r1g/(r1g-rp) EllipticPi[(rp(r1g-r2g))/((r1g-rp)r2g),\[Phi],krg]- ellF);
	\[ScriptCapitalI]mreg=2/Sqrt[(r1g-r3g)r2g] (r1g/(r1g-rm) EllipticPi[(rm(r1g-r2g))/((r1g-rm)r2g),\[Phi],krg]- ellF);

	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	1/Sqrt[1-EEg^2](Lzg*\[ScriptCapitalI]-rp*Lzg 1/(rp-rm)(rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)+2a*EEg 1/(rp-rm)(rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg))-wr*\[Phi]gfunref[a,p,e,x]
]


(* ::Subsubsection::Closed:: *)
(*Fourier coefficients of the geodesic coordinate time and azimuthal trajectory*)


tgcoeffICr2gfun[nmax_,a_,p_,e_,x_]:=Module[{stepsr,wrlist,tglist,ExpniTable},
	stepsr=4*nmax;
	wrlist=Table[i,{i,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	tglist=tgICr2gfun[wrlist,a,p,e,x];

	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e}]],{n,-nmax,nmax},{i,1,stepsr}];
	Chop[(ExpniTable . tglist)/stepsr,10^(-16)]
]


\[Phi]gcoeffICr2gfun[nmax_,a_,p_,e_,x_]:=Module[{stepsr,wrlist,\[Phi]glist,ExpniTable},
	stepsr=4*nmax;
	wrlist=Table[i,{i,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	\[Phi]glist=\[Phi]gICr2gfun[wrlist,a,p,e,x];

	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e}]],{n,-nmax,nmax},{i,1,stepsr}];
	Chop[(ExpniTable . \[Phi]glist)/stepsr,10^(-16)]
]


tgcoeffICr1gfun[nmax_,a_,p_,e_,x_]:=Module[{stepsr,wrlist,tglist,ExpniTable},
	stepsr=4*nmax;
	wrlist=Table[i,{i,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	tglist=tgICr1gfun[wrlist,a,p,e,x];

	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e}]],{n,-nmax,nmax},{i,1,stepsr}];
	Chop[(ExpniTable . tglist)/stepsr,10^(-16)]
]


\[Phi]gcoeffICr1gfun[nmax_,a_,p_,e_,x_]:=Module[{stepsr,wrlist,\[Phi]glist,ExpniTable},
	stepsr=4*nmax;
	wrlist=Table[i,{i,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	\[Phi]glist=\[Phi]gICr1gfun[wrlist,a,p,e,x];

	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e}]],{n,-nmax,nmax},{i,1,stepsr}];
	Chop[(ExpniTable . \[Phi]glist)/stepsr,10^(-16)]
]


(* ::Subsubsection::Closed:: *)
(*Expansion t_g and \[Phi]_g as Fourier series. Fourier expansion from analytic solutions*)


\[CapitalDelta]tgFourier[wr_,coeff_]:=Module[{dim},
	dim=(Length[coeff]-1)/2;
	Re[Sum[2Sin[n*wr]Im[coeff][[n+dim+1]],{n,1,dim}]]
];


\[CapitalDelta]\[Phi]gFourier[wr_,coeff_]:=Module[{dim},
	dim=(Length[coeff]-1)/2;
	Re[Sum[2Sin[n*wr]Im[coeff][[n+dim+1]],{n,1,dim}]]
];


(* ::Subsubsection::Closed:: *)
(*Expansion t_g and \[Phi]_g as Fourier series. Fourier expansion from integration velocity*)


\[CapitalDelta]trgIntVel[wr_,\[CapitalUpsilon]rg_,coeff_]:=Module[{dim},
	dim=(Length[coeff]-1)/2;
	Re[Sum[2Sin[n wr]coeff[[n+dim+1]]/(n \[CapitalUpsilon]rg),{n,1,dim}]]
];


\[CapitalDelta]\[Phi]rgIntVel[wr_,\[CapitalUpsilon]rg_,coeff_]:=Module[{dim},
	dim=(Length[coeff]-1)/2;
	Re[Sum[2Sin[n wr]coeff[[n+dim+1]]/(n \[CapitalUpsilon]rg),{n,1,dim}]]
];


(* ::Subsection::Closed:: *)
(*Geodesics periodic velocity*)


(* ::Subsubsection::Closed:: *)
(*Radial geodesic velocity*)


drgd\[Lambda]fun[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,r1g,r2g,r3g,krg,ellK,jSN,jCN},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	ellK=EllipticK[krg];
	jSN=JacobiSN[1/\[Pi]*ellK*wr,krg];
	jCN=JacobiCN[1/\[Pi]*ellK*wr,krg];

	Sqrt[1-EEg^2]((r1g-r2g)(r1g-r3g)(r2g-r3g)jCN*jSN Sqrt[r2g(r1g-r3g)+(-r1g+r2g)r3g*jSN^2])/(-r1g+r3g+(r1g-r2g)jSN^2)^2
]


drgd\[Lambda]funOfr[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,r1g,r2g,r3g,krg,ellK,jSN,jCN},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	Sqrt[1-EEg^2]Sqrt[(r1g-r)(r-r2g)(r-r3g)r]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate-time geodesic velocity*)


VtrgICr2gfun[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,rg},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	rg=rgICr2gfun[wr,a,p,e,x];

	EEg((rg^2+a^2)^2/(rg^2-2rg+a^2))-(2a*rg)/(rg^2-2rg+a^2)Lzg-EEg*a^2
]


VtrgICr1gfun[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,rg},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	rg=rgICr1gfun[wr,a,p,e,x];

	EEg((rg^2+a^2)^2/(rg^2-2rg+a^2))-(2a*rg)/(rg^2-2rg+a^2)Lzg-EEg*a^2
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal geodesic velocity*)


V\[Phi]rgICr2gfun[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,rg},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	rg=rgICr2gfun[wr,a,p,e,x];

	a/(rg^2-2rg+a^2)(EEg(rg^2+a^2)-a*Lzg)-a*EEg+Lzg
]


V\[Phi]rgICr1gfun[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,rg},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	rg=rgICr1gfun[wr,a,p,e,x];

	a/(rg^2-2rg+a^2)(EEg(rg^2+a^2)-a*Lzg)-a*EEg+Lzg
]


(* ::Subsection::Closed:: *)
(*Geodesics homoclinic*)


(* ::Subsubsection::Closed:: *)
(*Geodesic radial trajectory*)


rgfunHom[\[Lambda]_,a_,p_,e_,xg_]:=Module[{r1g,r2g,EEg,aux},
	EEg=EEgfun[a,p,e,xg];
	r1g=p/(1-e);
	r2g=p/(1+e);
	aux=1/2*Sqrt[(1-EEg^2)r2g(r1g-r2g)];
	
	r1g*r2g/(r2g+(r1g-r2g)Tanh[aux*\[Lambda]]^2)
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory*)


tgfunHom[r_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1g,r2g,\[ScriptCapitalI],\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over,\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	If[e==0,
		0,
		rp=1+Sqrt[1-a^2];
		rm=1-Sqrt[1-a^2];
		EEg=EEgfun[a,p,e,x];
		Lzg=Lzgfun[a,p,e,x];
		r1g=p/(1-e);
		r2g=p/(1+e);
		\[ScriptCapitalI]=-1/(2Sqrt[r2g(r1g-r2g)])Log[(Sqrt[r(r1g-r2g)]-Sqrt[(-r+r1g)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(-r+r1g)r2g])^2];
		\[ScriptCapitalI]rover=2 ArcCos[Sqrt[r/r1g]]+r2g*\[ScriptCapitalI];
		\[ScriptCapitalI]r2over=Sqrt[r(-r+r1g)]+r1g ArcCos[Sqrt[r/r1g]]+r2g*\[ScriptCapitalI]rover;
		\[ScriptCapitalI]p=\[ScriptCapitalI]/(r2g-rp)+1/(2Sqrt[r1g-rp](r2g-rp)Sqrt[rp])Log[(Sqrt[r(r1g-rp)]-Sqrt[(-r+r1g)rp])^2/(Sqrt[r(r1g-rp)]+Sqrt[(-r+r1g)rp])^2];
		\[ScriptCapitalI]mreg=rm*\[ScriptCapitalI]/(r2g-rm)+Sqrt[rm]/(2Sqrt[r1g-rm](r2g-rm))Log[(Sqrt[r (r1g-rm)]-Sqrt[(-r+r1g) rm])^2/(Sqrt[r(r1g-rm)]+Sqrt[(-r+r1g) rm])^2];

		1/Sqrt[1-EEg^2](4EEg*\[ScriptCapitalI]+2EEg*\[ScriptCapitalI]rover+EEg*\[ScriptCapitalI]r2over-4EEg rp/(rp-rm) (rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)-2(-4EEg+a*Lzg) 1/(rp-rm) (rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg))
	]
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal trajectory*)


\[Phi]gfunHom[r_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1g,r2g,\[ScriptCapitalI],\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	If[e==0,
		0,
		rp=1+Sqrt[1-a^2];
		rm=1-Sqrt[1-a^2];
		EEg=EEgfun[a,p,e,x];
		Lzg=Lzgfun[a,p,e,x];
		r1g=p/(1-e);
		r2g=p/(1+e);
		\[ScriptCapitalI]=-1/(2Sqrt[r2g(r1g-r2g)])Log[(Sqrt[r(r1g-r2g)]-Sqrt[(-r+r1g)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(-r+r1g)r2g])^2];
		\[ScriptCapitalI]p=\[ScriptCapitalI]/(r2g-rp)+1/(2Sqrt[r1g-rp](r2g-rp)Sqrt[rp])Log[(Sqrt[r(r1g-rp)]-Sqrt[(-r+r1g)rp])^2/(Sqrt[r(r1g-rp)]+Sqrt[(-r+r1g)rp])^2];
		\[ScriptCapitalI]mreg=rm*\[ScriptCapitalI]/(r2g-rm)+Sqrt[rm]/(2Sqrt[r1g-rm](r2g-rm))Log[(Sqrt[r(r1g-rm)]-Sqrt[(-r+r1g)rm])^2/(Sqrt[r(r1g-rm)]+Sqrt[(-r+r1g)rm])^2];

		1/Sqrt[1-EEg^2] (Lzg*\[ScriptCapitalI]-Lzg rp/(rp-rm) (rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)+2a*EEg 1/(rp-rm)(rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg))
	]
]


(* ::Subsubsection::Closed:: *)
(*Geodesic radial velocity*)


drgd\[Lambda]funHom[r_,a_,p_,e_,x_]:=Module[{EEg,r1g,r2g},
	EEg=EEgfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);

	Sqrt[1-EEg^2](r-r2g)Sqrt[r(r1g-r)]
]


(* ::Subsection::Closed:: *)
(*Geodesics ISCO plunge*)


(* ::Subsubsection::Closed:: *)
(*Geodesic radial trajectory*)


rgfunISCOplunge[\[Lambda]_,a_,x_]:=Module[{r1g,EEg ,Lzg},
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,x];
	Lzg=Lzgfun[a,r1g,0,x];

	r1g-(4r1g)/(4+(1-EEg^2)r1g^2\[Lambda]^2)
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory*)


tgfunISCOplunge[r_,a_,x_]:=Module[{r1g,rp,rm,EEg ,Lzg,\[ScriptCapitalI],\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over,\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	r1g=ISCOradius[a,x];
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,x];
	Lzg=Lzgfun[a,r1g,0,x];

	\[ScriptCapitalI]=2Sqrt[r]/(r1g*Sqrt[r1g-r]);
	\[ScriptCapitalI]rover=2ArcSin[Sqrt[1-r/r1g]]+r1g*\[ScriptCapitalI];
	\[ScriptCapitalI]r2over=Sqrt[r(-r+r1g)]+r1g*ArcSin[Sqrt[1-r/r1g]]+r1g*\[ScriptCapitalI]rover;
	\[ScriptCapitalI]p=1/(r1g-rp)\[ScriptCapitalI]+1/(2(r1g-rp)^(3/2) Sqrt[rp])Log[(Sqrt[r(r1g-rp)]-Sqrt[(r1g-r)rp])^2/(Sqrt[r(r1g-rp)]+Sqrt[(r1g-r) rp])^2];
	\[ScriptCapitalI]mreg=rm/(r1g-rm)\[ScriptCapitalI]+Sqrt[rm]/(2(r1g-rm)^(3/2))Log[(Sqrt[r(r1g-rm)]-Sqrt[(r1g-r)rm])^2/(Sqrt[r(r1g-rm)]+Sqrt[(r1g-r)rm])^2];

	-1/Sqrt[1-EEg^2](4EEg*\[ScriptCapitalI]+2EEg*\[ScriptCapitalI]rover+EEg*\[ScriptCapitalI]r2over-4EEg*rp/(rp-rm)(rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)-2(-4EEg+a*Lzg)1/(rp-rm)(rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg))
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal trajectory*)


\[Phi]gfunISCOplunge[r_,a_,x_]:=Module[{r1g,rp,rm,EEg ,Lzg,\[ScriptCapitalI],\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	r1g=ISCOradius[a,x];
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,x];
	Lzg=Lzgfun[a,r1g,0,x];

	\[ScriptCapitalI]=2Sqrt[r]/(r1g*Sqrt[r1g-r]);
	\[ScriptCapitalI]p=1/(r1g-rp)\[ScriptCapitalI]+1/(2(r1g-rp)^(3/2) Sqrt[rp])Log[(Sqrt[r(r1g-rp)]-Sqrt[(r1g-r)rp])^2/(Sqrt[r(r1g-rp)]+Sqrt[(r1g-r) rp])^2];
	\[ScriptCapitalI]mreg=rm/(r1g-rm)\[ScriptCapitalI]+Sqrt[rm]/(2(r1g-rm)^(3/2))Log[(Sqrt[r(r1g-rm)]-Sqrt[(r1g-r)rm])^2/(Sqrt[r(r1g-rm)]+Sqrt[(r1g-r)rm])^2];

	-1/Sqrt[1-EEg^2](Lzg*\[ScriptCapitalI]-Lzg*rp/(rp-rm)(rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)+2a*EEg 1/(rp-rm)(rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg))
]


(* ::Subsubsection::Closed:: *)
(*Geodesic radial velocity*)


drgd\[Lambda]funISCOplunge[r_,a_,x_]:=Module[{EEg,r1g},
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,x];

	(r-r1g)Sqrt[1-EEg^2]Sqrt[r(r1g-r)]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate-time geodesic velocity*)


VtrgfunISCOplunge[r_,a_,x_]:=Module[{r1g,EEg,Lzg},
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,x];
	Lzg=Lzgfun[a,r1g,0,x];

	EEg((r^2+a^2)^2/(r^2-2r+a^2))-(2a*r)/(r^2-2r+a^2)Lzg-EEg*a^2
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal geodesic velocity*)


V\[Phi]rgfunISCOplunge[r_,a_,x_]:=Module[{r1g,EEg,Lzg},
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,x];
	Lzg=Lzgfun[a,r1g,0,x];

	a/(r^2-2r+a^2)(EEg(r^2+a^2)-a*Lzg)-a*EEg+Lzg
]


(* ::Subsection::Closed:: *)
(*Geodesics Critical plunge*)


(* ::Subsubsection::Closed:: *)
(*Geodesic radial trajectory*)


rgfunCritplunge[\[Lambda]_,a_,p_,e_,x_]:=Module[{r1g,r2g,EEg,aux},
	EEg=EEgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	aux=1/2*Sqrt[(1-EEg^2)r2g(r1g-r2g)];

	(r1g(r1g-r2g)Tanh[aux*\[Lambda]]^2)/(r2g+(r1g-r2g)Tanh[aux*\[Lambda]]^2)
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory*)


tgfunCritplunge[r_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1g,r2g,\[ScriptCapitalI],\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over,\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	\[ScriptCapitalI]=-1/(2Sqrt[(r1g-r2g)r2g]) Log[(Sqrt[r(r1g-r2g)]-Sqrt[(r1g-r)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(r1g-r)r2g])^2];
	\[ScriptCapitalI]rover=2ArcCos[Sqrt[r/r1g]]+r2g*\[ScriptCapitalI];
	\[ScriptCapitalI]r2over=Sqrt[r(r1g-r)]+r1g*ArcCos[Sqrt[r/r1g]]+r2g*\[ScriptCapitalI]rover;
	\[ScriptCapitalI]p=\[ScriptCapitalI]/(r2g-rp)+1/(2Sqrt[r1g-rp](r2g-rp)Sqrt[rp])Log[(Sqrt[r(r1g-rp)]-Sqrt[(r1g-r)rp])^2/(Sqrt[r(r1g-rp)]+Sqrt[(r1g-r)rp])^2];
	\[ScriptCapitalI]mreg=(rm*\[ScriptCapitalI])/(r2g-rm)+Sqrt[rm]/(2Sqrt[r1g-rm](r2g-rm))Log[(Sqrt[r(r1g-rm)]-Sqrt[(r1g-r)rm])^2/(Sqrt[r(r1g-rm)]+Sqrt[(r1g-r)rm])^2];

	-1/Sqrt[1-EEg^2](4EEg*\[ScriptCapitalI]+2EEg*\[ScriptCapitalI]rover+EEg*\[ScriptCapitalI]r2over-4EEg rp/(rp-rm)(rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)-2(-4EEg+a*Lzg) 1/(rp-rm)(rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg))
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal trajectory*)


\[Phi]gfunCritplunge[r_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1g,r2g,\[ScriptCapitalI],\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	\[ScriptCapitalI]=-1/(2Sqrt[(r1g-r2g)r2g])Log[(Sqrt[r(r1g-r2g)]-Sqrt[(r1g-r)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(r1g-r)r2g])^2];
	\[ScriptCapitalI]p=\[ScriptCapitalI]/(r2g-rp)+1/(2Sqrt[r1g-rp](r2g-rp)Sqrt[rp]) Log[(Sqrt[r(r1g-rp)]-Sqrt[(r1g-r)rp])^2/(Sqrt[r(r1g-rp)]+Sqrt[(r1g-r)rp])^2];
	\[ScriptCapitalI]mreg=(rm*\[ScriptCapitalI])/(r2g-rm)+Sqrt[rm]/(2Sqrt[r1g-rm](r2g-rm))Log[(Sqrt[r(r1g-rm)]-Sqrt[(r1g-r)rm])^2/(Sqrt[r(r1g-rm)]+Sqrt[(r1g-r)rm])^2];

	-1/Sqrt[1-EEg^2] (Lzg*\[ScriptCapitalI]-Lzg*rp/(rp-rm)(rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)+2a*EEg 1/(rp-rm)(rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg))
]


(* ::Subsection::Closed:: *)
(*Geodesics Generic plunge*)


(* ::Subsubsection::Closed:: *)
(*Geodesic radial trajectory*)


rgfunPlunge[\[Lambda]_,a_,p_,e_,x_]:=Module[{EEg,r1g,\[Rho]rg,\[Rho]ig,Ar1g,B0,krg},
	EEg=EEgfuncc[p,e];
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];

	Ar1g=Sqrt[\[Rho]ig^2+(r1g-\[Rho]rg)^2];
	B0=Sqrt[\[Rho]ig^2+\[Rho]rg^2];
	krg=(-(Ar1g-B0)^2+r1g^2)/(4Ar1g*B0);

	(B0*r1g(1+JacobiCN[Sqrt[Ar1g*B0]Sqrt[1-EEg^2]\[Lambda],krg]))/(Ar1g+B0-(Ar1g-B0)JacobiCN[Sqrt[Ar1g*B0]Sqrt[1-EEg^2]\[Lambda],krg])
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory*)


tgfunPlunge[r_,a_,p_,e_,x_]:=Module[{r1g,\[Rho]rg,\[Rho]ig,rp,rm,EEg,Lzg,Ar1g,B0,krg,\[Alpha]p,\[Alpha]m,\[Gamma]r,\[Phi],ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over,\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];

	Ar1g=Sqrt[\[Rho]ig^2+(r1g-\[Rho]rg)^2];
	B0=Sqrt[\[Rho]ig^2+\[Rho]rg^2];
	krg=(-(Ar1g-B0)^2+r1g^2)/(4Ar1g*B0);
	\[Alpha]p=(-B0*r1g-Ar1g*rp+B0*rp)/(B0*r1g-Ar1g*rp-B0*rp);
	\[Alpha]m=(-B0*r1g-Ar1g*rm+B0*rm)/(B0*r1g-Ar1g*rm-B0*rm);
	\[Gamma]r=-((Ar1g-B0)^2/(4Ar1g*B0));
	\[Phi]=ArcCos[(-Ar1g*r+B0(r1g-r))/(Ar1g*r+B0(r1g-r))];
	
	ellF=EllipticF[\[Pi]-\[Phi],krg];
	ellE=EllipticE[\[Pi]-\[Phi],krg];
	ellPi=EllipticPi[\[Gamma]r,\[Pi]-\[Phi],krg];

	\[ScriptCapitalI]=ellF/Sqrt[Ar1g*B0];
	\[ScriptCapitalI]rover=ArcTan[(r1g*Sin[\[Phi]])/(2Sqrt[Ar1g*B0] Sqrt[1-krg*Sin[\[Phi]]^2])]-(B0*r1g)/(Ar1g-B0)*\[ScriptCapitalI]+((Ar1g+B0)r1g*ellPi)/(2(Ar1g-B0)Sqrt[Ar1g*B0]);
	\[ScriptCapitalI]r2over=((Ar1g-B0)Sqrt[r(r1g-r)(\[Rho]ig^2+(r-\[Rho]rg)^2)])/(Ar1g*r+B0(r1g-r))+Sqrt[Ar1g*B0]*ellE+1/2(r1g+2\[Rho]rg)\[ScriptCapitalI]rover-(B0*r1g(r1g-2\[Rho]rg))/(2(Ar1g-B0))*\[ScriptCapitalI];
	\[ScriptCapitalI]p=(4r1g)/(B0^2(r1g-rp)^2-Ar1g^2*rp^2) (Ar1g*B0)/Sqrt[2(Ar1g*B0+B0^2-r1g*\[Rho]rg)] (EllipticPi[\[Alpha]p^2,krg/(-1+krg)]+EllipticPi[\[Alpha]p^2,(\[Pi]/2-\[Phi]),krg/(-1+krg)])-(Ar1g-B0)/(B0(r1g-rp)+Ar1g*rp) \[ScriptCapitalI]-1/(4Sqrt[(r1g-rp)rp(\[Rho]ig^2+(rp-\[Rho]rg)^2)]) Log[(Sqrt[-\[Alpha]p^2+krg(-1+\[Alpha]p^2)]Sin[\[Phi]]-Sqrt[1-\[Alpha]p^2] Sqrt[1-krg*Sin[\[Phi]]^2])^2/(Sqrt[-\[Alpha]p^2+krg(-1+\[Alpha]p^2)]Sin[\[Phi]]+Sqrt[1-\[Alpha]p^2] Sqrt[1-krg*Sin[\[Phi]]^2])^2];
	\[ScriptCapitalI]mreg=(4r1g*rm)/(B0^2(r1g-rm)^2-Ar1g^2*rm^2) (Ar1g*B0)/Sqrt[2(Ar1g*B0+B0^2-r1g*\[Rho]rg)] (EllipticPi[\[Alpha]m^2,krg/(-1+krg)]+EllipticPi[\[Alpha]m^2,(\[Pi]/2-\[Phi]),krg/(-1+krg)])-(rm(Ar1g-B0) )/(B0(r1g-rm)+Ar1g*rm) \[ScriptCapitalI]-Sqrt[rm]/(4Sqrt[(r1g-rm)(\[Rho]ig^2+(rm-\[Rho]rg)^2)]) Log[(Sqrt[-\[Alpha]m^2+krg (-1+\[Alpha]m^2)]Sin[\[Phi]]-Sqrt[1-\[Alpha]m^2] Sqrt[1-krg*Sin[\[Phi]]^2])^2/(Sqrt[-\[Alpha]m^2+krg (-1+\[Alpha]m^2)]Sin[\[Phi]]+Sqrt[1-\[Alpha]m^2] Sqrt[1-krg*Sin[\[Phi]]^2])^2];
	
	1/Sqrt[1-EEg^2](4EEg*\[ScriptCapitalI]+2EEg*\[ScriptCapitalI]rover+EEg*\[ScriptCapitalI]r2over-4EEg*rp/(rp-rm)(rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)-2(-4EEg+a*Lzg)/(rp-rm)(rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg))
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal trajectory*)


\[Phi]gfunPlunge[r_,a_,p_,e_,x_]:=Module[{r1g,\[Rho]rg,\[Rho]ig,rp,rm,EEg,Lzg,Ar1g,B0,krg,\[Alpha]p,\[Alpha]m,\[Phi],\[ScriptCapitalI],\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];

	Ar1g=Sqrt[\[Rho]ig^2+(r1g-\[Rho]rg)^2];
	B0=Sqrt[\[Rho]ig^2+\[Rho]rg^2];
	krg=(-(Ar1g-B0)^2+r1g^2)/(4Ar1g*B0);
	\[Alpha]p=(-B0*r1g-Ar1g*rp+B0*rp)/(B0*r1g-Ar1g*rp-B0*rp);
	\[Alpha]m=(-B0*r1g-Ar1g*rm+B0*rm)/(B0*r1g-Ar1g*rm-B0*rm);
	\[Phi]=ArcCos[(-Ar1g*r+B0(r1g-r))/(Ar1g*r+B0(r1g-r))];

	\[ScriptCapitalI]=EllipticF[\[Pi]-\[Phi],krg]/Sqrt[Ar1g*B0];
	\[ScriptCapitalI]p=(4r1g)/(B0^2(r1g-rp)^2-Ar1g^2*rp^2) (Ar1g*B0)/Sqrt[2(Ar1g*B0+B0^2-r1g*\[Rho]rg)] (EllipticPi[\[Alpha]p^2,krg/(-1+krg)]+EllipticPi[\[Alpha]p^2,(\[Pi]/2-\[Phi]),krg/(-1+krg)])-(Ar1g-B0)/(B0(r1g-rp)+Ar1g*rp) \[ScriptCapitalI]-1/(4Sqrt[(r1g-rp)rp(\[Rho]ig^2+(rp-\[Rho]rg)^2)]) Log[(Sqrt[-\[Alpha]p^2+krg(-1+\[Alpha]p^2)]Sin[\[Phi]]-Sqrt[1-\[Alpha]p^2] Sqrt[1-krg*Sin[\[Phi]]^2])^2/(Sqrt[-\[Alpha]p^2+krg(-1+\[Alpha]p^2)]Sin[\[Phi]]+Sqrt[1-\[Alpha]p^2] Sqrt[1-krg*Sin[\[Phi]]^2])^2];
	\[ScriptCapitalI]mreg=(4r1g*rm)/(B0^2(r1g-rm)^2-Ar1g^2*rm^2) (Ar1g*B0)/Sqrt[2(Ar1g*B0+B0^2-r1g*\[Rho]rg)] (EllipticPi[\[Alpha]m^2,krg/(-1+krg)]+EllipticPi[\[Alpha]m^2,(\[Pi]/2-\[Phi]),krg/(-1+krg)])-(rm(Ar1g-B0) )/(B0(r1g-rm)+Ar1g*rm) \[ScriptCapitalI]-Sqrt[rm]/(4Sqrt[(r1g-rm)(\[Rho]ig^2+(rm-\[Rho]rg)^2)]) Log[(Sqrt[-\[Alpha]m^2+krg (-1+\[Alpha]m^2)]Sin[\[Phi]]-Sqrt[1-\[Alpha]m^2] Sqrt[1-krg*Sin[\[Phi]]^2])^2/(Sqrt[-\[Alpha]m^2+krg (-1+\[Alpha]m^2)]Sin[\[Phi]]+Sqrt[1-\[Alpha]m^2] Sqrt[1-krg*Sin[\[Phi]]^2])^2];

	1/Sqrt[1-EEg^2](Lzg*\[ScriptCapitalI]-Lzg*rp/(rp-rm)(rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)+2a*EEg/(rp-rm)(rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg))
]


(* ::Subsubsection::Closed:: *)
(*Geodesic radial velocity*)


drgd\[Lambda]funPlunge[r_,a_,p_,e_,x_]:=Module[{EEg,r1g,\[Rho]rg,\[Rho]ig},
	EEg=EEgfuncc[p,e];
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];

	-Sqrt[1-EEg^2]Sqrt[(r-\[Rho]rg-I*\[Rho]ig)(r-\[Rho]rg+I*\[Rho]ig)]Sqrt[(r1g-r)r]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate-time geodesic velocity*)


VtrgfunPlunge[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg},
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];

	EEg((r^2+a^2)^2/(r^2-2r+a^2))-(2a*r)/(r^2-2r+a^2)Lzg-EEg*a^2
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal geodesic velocity*)


V\[Phi]rgfunPlunge[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg},
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];

	a/(r^2-2r+a^2)(EEg(r^2+a^2)-a*Lzg)-a*EEg+Lzg
]


(* ::Subsection::Closed:: *)
(*Geodesics Plunge related to bound orbits*)


(* ::Subsubsection::Closed:: *)
(*Geodesic radial trajectory*)


rgfunPlungeBoundOrbit[\[Lambda]_,a_,p_,e_,x_]:=Module[{EEg,r1g,r2g,r3g,krg,Yrg},
	EEg=EEgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	Yrg=1/2 Sqrt[(1-EEg^2)r2g(r1g-r3g)];

	(r2g*r3g*JacobiCN[Yrg*\[Lambda],krg]^2)/(r2g-r3g*JacobiSN[Yrg*\[Lambda],krg]^2)
]


rgfunPlungeBoundOrbitDoubleRoot[\[Lambda]_,a_,p_,x_]:=Module[{r1g,r3g,EEg,Yrg},
	EEg=EEgfun[a,p,0,x];
	r1g=p;
	r3g=2/(1-EEg^2)-2p;
	Yrg=1/2 Sqrt[(1-EEg^2)r1g(r1g-r3g)];

	(r1g*r3g*Cos[Yrg*\[Lambda]]^2)/(r1g-r3g*Sin[Yrg*\[Lambda]]^2)
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory*)


tgfunPlungeBoundOrbit[r_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1g,r2g,r3g,krg,\[Gamma]r,\[Phi],ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over,\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Gamma]r=r3g/r2g;
	\[Phi]=ArcSin[Sqrt[(r2g(r3g-r))/(r3g(r2g-r))]];
	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	ellPi=EllipticPi[\[Gamma]r,\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[r2g(r1g-r3g)];
	\[ScriptCapitalI]rover=2 /Sqrt[r2g(r1g-r3g)] (r2g*ellF+(-r2g+r3g)ellPi);
	\[ScriptCapitalI]r2over=-r3g Sqrt[((r1g-r)(r2g-r))/(r2g(r2g-r3g))]Sin[\[Phi]]Cos[\[Phi]]+Sqrt[r2g(r1g-r3g)]ellE-((r2g-r3g)(r1g+r2g+r3g))/Sqrt[r2g(r1g-r3g)] ellPi+(r2g(r2g+r3g))/Sqrt[r2g(r1g-r3g)] ellF;
	\[ScriptCapitalI]p=2/(Sqrt[r2g (r1g-r3g)](r2g-rp)) (ellF+(r2g-r3g)/(r3g-rp) EllipticPi[(r3g (-r2g+rp))/(r2g (-r3g+rp)),\[Phi],krg]);
	\[ScriptCapitalI]mreg=(2rm)/(Sqrt[r2g (r1g-r3g)](r2g-rm)) (ellF+(r2g-r3g) /(r3g-rm) EllipticPi[(r3g (-r2g+rm))/(r2g (-r3g+rm)),\[Phi],krg]);

	-(1/Sqrt[1-EEg^2])(4EEg*\[ScriptCapitalI]+2EEg*\[ScriptCapitalI]rover+EEg*\[ScriptCapitalI]r2over-4EEg rp/(rp-rm) (rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)-2(-4EEg+a*Lzg) 1/(rp-rm) (rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg))
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory - double root*)


tgfunPlungeBoundOrbitDoubleRoot[r_,a_,p_,x_]:=Module[{rp,rm,EEg,Lzg,r1g,r3g,krg,\[Phi],\[ScriptCapitalI],\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over,\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,0,x];
	Lzg=Lzgfun[a,p,0,x];
	r1g=p;
	r3g=2/(1-EEg^2)-2p;

	krg=0;
	\[Phi]=ArcSin[Sqrt[(r1g(r3g-r))/(r3g(r1g-r))]];

	\[ScriptCapitalI]=(2\[Phi])/Sqrt[r1g(r1g-r3g)];
	\[ScriptCapitalI]rover=(2r1g*\[Phi])/Sqrt[r1g(r1g-r3g)]-2 1/2 I*Log[(Sqrt[r]-I*Sqrt[r3g-r])/(Sqrt[r]+I*Sqrt[r3g-r])];
	\[ScriptCapitalI]r2over=-Sqrt[(r3g-r)r]+(2r1g^2*\[Phi])/Sqrt[r1g(r1g-r3g)]-(r1g(2r1g+r3g))/Sqrt[r1g(r1g-r3g)] 1/2 I*Log[(Sqrt[r]-I*Sqrt[r3g-r])/(Sqrt[r]+I*Sqrt[r3g-r])]Sqrt[1-r3g/r1g];

	\[ScriptCapitalI]p=2/(Sqrt[r1g(r1g-r3g)](r1g-rp)) (\[Phi]-1/4 r1g/rp Sqrt[((r1g-r3g)rp)/(r1g(r3g-rp))]Log[((Sqrt[r1g(r3g-rp)]-Sqrt[rp(r1g-r3g)]*Tan[\[Phi]])/(Sqrt[r1g(r3g-rp)]+Sqrt[rp(r1g-r3g)]*Tan[\[Phi]]))^2]);
	\[ScriptCapitalI]mreg=2/(Sqrt[r1g(r1g-r3g)](r1g-rm)) (rm*\[Phi]-1/4 r1g Sqrt[((r1g-r3g)rm)/(r1g(r3g-rm))]Log[((Sqrt[r1g(r3g-rm)]-Sqrt[rm(r1g-r3g)]*Tan[\[Phi]])/(Sqrt[r1g(r3g-rm)]+Sqrt[rm(r1g-r3g)]*Tan[\[Phi]]))^2]);

	-(1/Sqrt[1-EEg^2])(4EEg*\[ScriptCapitalI]+2EEg*\[ScriptCapitalI]rover+EEg*\[ScriptCapitalI]r2over-4EEg rp/(rp-rm) (rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)-2(-4EEg+a*Lzg) 1/(rp-rm) (rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg))
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal trajectory*)


\[Phi]gfunPlungeBoundOrbit[r_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1g,r2g,r3g,krg,\[Phi],ellF,\[ScriptCapitalI],\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Phi]=ArcSin[Sqrt[(r2g(r3g-r))/(r3g(r2g-r))]];
	ellF=EllipticF[\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[r2g(r1g-r3g)];
	\[ScriptCapitalI]p=2/(Sqrt[r2g (r1g-r3g)](r2g-rp)) (ellF+(r2g-r3g)/(r3g-rp) EllipticPi[(r3g(-r2g+rp))/(r2g(-r3g+rp)),\[Phi],krg]);
	\[ScriptCapitalI]mreg=(2rm)/(Sqrt[r2g (r1g-r3g)](r2g-rm)) (ellF+(r2g-r3g) /(r3g-rm) EllipticPi[(r3g(-r2g+rm))/(r2g(-r3g+rm)),\[Phi],krg]);

	-(1/Sqrt[1-EEg^2])(Lzg*\[ScriptCapitalI]-Lzg rp/(rp-rm) (rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)+2a*EEg 1/(rp-rm) (rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg))
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal trajectory - double root*)


\[Phi]gfunPlungeBoundOrbitDoubleRoot[r_,a_,p_,x_]:=Module[{rp,rm,EEg,Lzg,r1g,r3g,krg,\[Phi],\[ScriptCapitalI],\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,0,x];
	Lzg=Lzgfun[a,p,0,x];
	r1g=p;
	r3g=2/(1-EEg^2)-2p;

	krg=0;
	\[Phi]=ArcSin[Sqrt[(r1g(r3g-r))/(r3g(r1g-r))]];
	\[ScriptCapitalI]=(2\[Phi])/Sqrt[r1g(r1g-r3g)];
	\[ScriptCapitalI]p=2/(Sqrt[r1g(r1g-r3g)](r1g-rp)) (\[Phi]-1/4 r1g/rp Sqrt[((r1g-r3g)rp)/(r1g(r3g-rp))]Log[((Sqrt[r1g(r3g-rp)]-Sqrt[rp(r1g-r3g)]*Tan[\[Phi]])/(Sqrt[r1g(r3g-rp)]+Sqrt[rp(r1g-r3g)]*Tan[\[Phi]]))^2]);
	\[ScriptCapitalI]mreg=2/(Sqrt[r1g(r1g-r3g)](r1g-rm)) (rm*\[Phi]-1/4 r1g Sqrt[((r1g-r3g)rm)/(r1g(r3g-rm))]Log[((Sqrt[r1g(r3g-rm)]-Sqrt[rm(r1g-r3g)]*Tan[\[Phi]])/(Sqrt[r1g(r3g-rm)]+Sqrt[rm(r1g-r3g)]*Tan[\[Phi]]))^2]);

	-(1/Sqrt[1-EEg^2])(Lzg*\[ScriptCapitalI]-Lzg rp/(rp-rm) (rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)+2a*EEg 1/(rp-rm) (rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg))
]


(* ::Subsubsection::Closed:: *)
(*Geodesic radial velocity*)


drgd\[Lambda]funPlungeBoundOrbit[r_,a_,p_,e_,x_]:=Module[{EEg,r1g,r2g,r3g},
	EEg=EEgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	-Sqrt[1-EEg^2]Sqrt[(r-r2g)(r-r3g)]Sqrt[(r1g-r)r]
]


(* ::Subsection::Closed:: *)
(*Scattering orbits - single roots*)


(* ::Subsubsection:: *)
(*Geodesic radial trajectory*)


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory*)


tgfunScattering[r_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1gp,r2g,r3g,krg,\[Phi],ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r,\[ScriptCapitalI]r2,\[ScriptCapitalI]p,\[ScriptCapitalI]m},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	r1gp=p/(e-1);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(-r1gp+r2g);

	krg=(r3g(r2g+r1gp))/(r2g(r3g+r1gp));

	\[Phi]=ArcSin[Sqrt[((r3g+r1gp)(r-r2g))/((r2g+r1gp)(r-r3g))]];
	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	ellPi=EllipticPi[(r1gp+r2g)/(r1gp+r3g),\[Phi],krg];

	\[ScriptCapitalI]=(2 ellF)/Sqrt[r2g (r1gp+r3g)];
	\[ScriptCapitalI]r= 2/Sqrt[r2g(r1gp+r3g)] ((r2g-r3g)ellPi+r3g ellF);
	\[ScriptCapitalI]r2= (r1gp (r2g-r3g)+r3g (r2g+r3g))/Sqrt[r2g(r1gp+r3g)] EllipticF[\[Phi],krg]-Sqrt[r2g(r1gp+r3g)]EllipticE[\[Phi],krg]+((r2g-r3g)(-r1gp+r2g+r3g))/Sqrt[r2g(r1gp+r3g)] ellPi+((r1gp+r2g)Sqrt[r2g(r1gp+r3g)]Sin[\[Phi]]Cos[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(r1gp+r3g-(r1gp+r2g)Sin[\[Phi]]^2);

	\[ScriptCapitalI]p=1/(r3g-rp) 2/Sqrt[r2g(r1gp+r3g)] (ellF-(r2g-r3g)/(r2g-rp) EllipticPi[((r3g-rp)(r1gp+r2g))/((r1gp+r3g)(r2g-rp)),\[Phi],krg]);
	\[ScriptCapitalI]m=1/(r3g-rm) 2/Sqrt[r2g(r1gp+r3g)] (ellF-(r2g-r3g)/(r2g-rm) EllipticPi[((r3g-rm)(r1gp+r2g))/((r1gp+r3g)(r2g-rm)),\[Phi],krg]);

	1/Sqrt[EEg^2-1] (4EEg*\[ScriptCapitalI]+2EEg*\[ScriptCapitalI]r+EEg*\[ScriptCapitalI]r2-4a^2*EEg 1/(rp-rm) (\[ScriptCapitalI]p-\[ScriptCapitalI]m)-2(-4EEg+a*Lzg) 1/(rp-rm) (rp*\[ScriptCapitalI]p-rm*\[ScriptCapitalI]m))
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal trajectory*)


\[Phi]gfunScattering[r_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1gp,r2g,r3g,krg,\[Phi],ellF,\[ScriptCapitalI],\[ScriptCapitalI]p,\[ScriptCapitalI]m},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	r1gp=p/(e-1);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(-r1gp+r2g);

	krg=(r3g(r2g+r1gp))/(r2g(r3g+r1gp));
	\[Phi]=ArcSin[Sqrt[((r3g+r1gp)(r-r2g))/((r2g+r1gp)(r-r3g))]];
	ellF=EllipticF[\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[r2g(r1gp+r3g)];
	\[ScriptCapitalI]p=1/(r3g-rp) 2/Sqrt[r2g(r1gp+r3g)] (ellF-(r2g-r3g)/(r2g-rp) EllipticPi[((r3g-rp)(r1gp+r2g))/((r1gp+r3g)(r2g-rp)),\[Phi],krg]);
	\[ScriptCapitalI]m=1/(r3g-rm) 2/Sqrt[r2g(r1gp+r3g)] (ellF-(r2g-r3g)/(r2g-rm) EllipticPi[((r3g-rm)(r1gp+r2g))/((r1gp+r3g)(r2g-rm)),\[Phi],krg]);

	1/Sqrt[EEg^2-1] (Lzg*\[ScriptCapitalI]-a^2*Lzg 1/(rp-rm) (\[ScriptCapitalI]p-\[ScriptCapitalI]m)+2a*EEg 1/(rp-rm) (rp*\[ScriptCapitalI]p-rm*\[ScriptCapitalI]m))
]


(* ::Subsubsection:: *)
(*Geodesic radial velocity*)


(* ::Subsection::Closed:: *)
(*Scattering orbits - double root*)


(* ::Subsubsection:: *)
(*Geodesic radial trajectory*)


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory*)


tgfunCapture[r_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1gp,r2g,\[ScriptCapitalI],\[ScriptCapitalI]r,\[ScriptCapitalI]r2,\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	r1gp=p/(e-1);
	r2g=p/(1+e);
	\[ScriptCapitalI]=-(1/(2Sqrt[r2g] Sqrt[r1gp+r2g]))Log[(Sqrt[r(r1gp+r2g)]-Sqrt[(r+r1gp)r2g])^2/(Sqrt[r(r1gp+r2g)]+Sqrt[(r+r1gp)r2g])^2];
	\[ScriptCapitalI]r=1/2 Log[((Sqrt[r]-Sqrt[r+r1gp])/(Sqrt[r]+Sqrt[r+r1gp]))^2]+r2g*\[ScriptCapitalI];
	\[ScriptCapitalI]r2=-Sqrt[r(r+r1gp)]-1/4 (r1gp-2r2g)Log[((Sqrt[r]-Sqrt[r+r1gp])/(Sqrt[r]+Sqrt[r+r1gp]))^2]+r2g^2*\[ScriptCapitalI];
	\[ScriptCapitalI]p=\[ScriptCapitalI]/ (r2g-rp)+1/(2Sqrt[r1gp+rp](r2g-rp)Sqrt[rp]) Log[(Sqrt[r(r1gp+rp)]-Sqrt[(r+r1gp) rp])^2/(Sqrt[(r+r1gp)rp]+Sqrt[r(r1gp+rp)])^2];
	\[ScriptCapitalI]mreg=(rm*\[ScriptCapitalI])/ (r2g-rm)+Sqrt[rm]/(2Sqrt[r1gp+rm](r2g-rm) ) Log[(Sqrt[r(r1gp+rm)]-Sqrt[(r+r1gp)rm])^2/(Sqrt[(r+r1gp) rm]+Sqrt[r(r1gp+rm)])^2];

	1/Sqrt[EEg^2-1] (4EEg*\[ScriptCapitalI]+2EEg*\[ScriptCapitalI]r+EEg*\[ScriptCapitalI]r2-4*EEg rp/(rp-rm) (rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)-2(-4EEg+a*Lzg) 1/(rp-rm) (rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg))
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal trajectory*)


\[Phi]gfunCapture[r_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1gp,r2g,\[ScriptCapitalI],\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	r1gp=p/(e-1);
	r2g=p/(1+e);
	\[ScriptCapitalI]=-(1/(2Sqrt[r2g] Sqrt[r1gp+r2g]))Log[(Sqrt[r(r1gp+r2g)]-Sqrt[(r+r1gp)r2g])^2/(Sqrt[r(r1gp+r2g)]+Sqrt[(r+r1gp)r2g])^2];
	\[ScriptCapitalI]p=\[ScriptCapitalI]/ (r2g-rp)+1/(2Sqrt[r1gp+rp](r2g-rp)Sqrt[rp]) Log[(Sqrt[r(r1gp+rp)]-Sqrt[(r+r1gp) rp])^2/(Sqrt[(r+r1gp)rp]+Sqrt[r(r1gp+rp)])^2];
	\[ScriptCapitalI]mreg=(rm*\[ScriptCapitalI])/ (r2g-rm)+Sqrt[rm]/(2Sqrt[r1gp+rm](r2g-rm) ) Log[(Sqrt[r(r1gp+rm)]-Sqrt[(r+r1gp)rm])^2/(Sqrt[(r+r1gp) rm]+Sqrt[r(r1gp+rm)])^2];

	1/Sqrt[EEg^2-1] (Lzg*\[ScriptCapitalI]-Lzg rp/(rp-rm) (rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)+2a*EEg 1/(rp-rm) (rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg))
]


(* ::Subsubsection:: *)
(*Geodesic radial velocity*)


(* ::Subsection::Closed:: *)
(*Plunge related to scattering orbits*)


(* ::Subsubsection:: *)
(*Geodesic radial trajectory*)


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory*)


tgfunPlungeScatteringOrbitDoubleRoot[r_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1gp,r2g,\[ScriptCapitalI],\[ScriptCapitalI]r,\[ScriptCapitalI]r2,\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	r1gp=p/(e-1);
	r2g=p/(1+e);
	\[ScriptCapitalI]=-(1/(2Sqrt[r2g] Sqrt[r1gp+r2g]))Log[(Sqrt[r(r1gp+r2g)]-Sqrt[(r+r1gp)r2g])^2/(Sqrt[r(r1gp+r2g)]+Sqrt[(r+r1gp)r2g])^2];
	\[ScriptCapitalI]r=1/2 Log[((Sqrt[r]-Sqrt[r+r1gp])/(Sqrt[r]+Sqrt[r+r1gp]))^2]+r2g*\[ScriptCapitalI];
	\[ScriptCapitalI]r2=-Sqrt[r(r+r1gp)]-1/4 (r1gp-2r2g)Log[((Sqrt[r]-Sqrt[r+r1gp])/(Sqrt[r]+Sqrt[r+r1gp]))^2]+r2g^2*\[ScriptCapitalI];
	\[ScriptCapitalI]p=\[ScriptCapitalI]/ (r2g-rp)+1/(2Sqrt[r1gp+rp](r2g-rp)Sqrt[rp]) Log[(Sqrt[r(r1gp+rp)]-Sqrt[(r+r1gp) rp])^2/(Sqrt[(r+r1gp)rp]+Sqrt[r(r1gp+rp)])^2];
	\[ScriptCapitalI]mreg=(rm*\[ScriptCapitalI])/ (r2g-rm)+Sqrt[rm]/(2Sqrt[r1gp+rm](r2g-rm) ) Log[(Sqrt[r(r1gp+rm)]-Sqrt[(r+r1gp)rm])^2/(Sqrt[(r+r1gp) rm]+Sqrt[r(r1gp+rm)])^2];

	-(1/Sqrt[EEg^2-1])(4EEg*\[ScriptCapitalI]+2EEg*\[ScriptCapitalI]r+EEg*\[ScriptCapitalI]r2-4*EEg rp/(rp-rm) (rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)-2(-4EEg+a*Lzg) 1/(rp-rm) (rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg))
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal trajectory*)


\[Phi]gfunPlungeScatteringOrbitDoubleRoot[r_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1gp,r2g,\[ScriptCapitalI],\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	r1gp=p/(e-1);
	r2g=p/(1+e);
	\[ScriptCapitalI]=-(1/(2Sqrt[r2g] Sqrt[r1gp+r2g]))Log[(Sqrt[r(r1gp+r2g)]-Sqrt[(r+r1gp)r2g])^2/(Sqrt[r(r1gp+r2g)]+Sqrt[(r+r1gp)r2g])^2];
	\[ScriptCapitalI]p=\[ScriptCapitalI]/ (r2g-rp)+1/(2Sqrt[r1gp+rp](r2g-rp)Sqrt[rp]) Log[(Sqrt[r(r1gp+rp)]-Sqrt[(r+r1gp) rp])^2/(Sqrt[(r+r1gp)rp]+Sqrt[r(r1gp+rp)])^2];
	\[ScriptCapitalI]mreg=(rm*\[ScriptCapitalI])/ (r2g-rm)+Sqrt[rm]/(2Sqrt[r1gp+rm](r2g-rm) ) Log[(Sqrt[r(r1gp+rm)]-Sqrt[(r+r1gp)rm])^2/(Sqrt[(r+r1gp) rm]+Sqrt[r(r1gp+rm)])^2];

	-(1/Sqrt[EEg^2-1])(Lzg*\[ScriptCapitalI]-Lzg rp/(rp-rm) (rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)+2a*EEg 1/(rp-rm) (rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg))
]


(* ::Subsubsection:: *)
(*Geodesic radial velocity*)


(* ::Subsection:: *)
(*Unbound plunges*)


(* ::Chapter::Closed:: *)
(*Back-end: spin corrections to constants of motion, radial potential roots, shift frequencies*)


(* ::Section::Closed:: *)
(*Shifts to the constants of motion*)


(* ::Subsection::Closed:: *)
(*Bound orbits - fixed turning points*)


\[Delta]EEfunFT[a_,p_,e_,x_]:=Module[{EEg,Lzg,r1g,r2g,Lzrd,\[ScriptCapitalD]},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	Lzrd=Lzg-a*EEg;
	
	\[ScriptCapitalD]=r1g^2r2g^2(-2Lzrd*EEg (r1g+r2g)^2+2EEg*Lzrd*r1g*r2g+2a*Lzrd^2+EEg*Lzg*r1g*r2g(r1g+r2g));

	Lzrd/\[ScriptCapitalD](-2a*Lzrd^2(r1g+r2g)-Lzg*EEg*r1g^2r2g^2-a*Lzg*Lzrd*r1g*r2g+a*Lzrd*Lzg (r1g+r2g)^2)
]


\[Delta]LzfunFT[a_,p_,e_,x_]:=Module[{EEg,Lzg,r1g,r2g,Lzrd,\[ScriptCapitalD]},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	Lzrd=Lzg-a*EEg;
	
	\[ScriptCapitalD]=r1g^2r2g^2(-2Lzrd*EEg (r1g+r2g)^2+2EEg*Lzrd*r1g*r2g+2a*Lzrd^2+EEg*Lzg*r1g*r2g(r1g+r2g));

	EEg/\[ScriptCapitalD]*r1g*r2g*Lzrd(a^3*Lzrd-a(a*EEg-3Lzrd)r1g*r2g-3EEg*r1g^2*r2g^2)+(a*Lzrd^2)/\[ScriptCapitalD](-2a*Lzrd(r1g+r2g)+a^2*EEg(r1g^2+r2g^2)+EEg(r1g^2+r2g^2)r1g*r2g+EEg(r1g^4+r2g^4))+EEg/\[ScriptCapitalD] r1g*r2g(EEg*Lzg*r1g^2*r2g^2(r1g+r2g)-3EEg*Lzrd*r1g*r2g(r1g^2+r2g^2))
]


(* ::Subsection::Closed:: *)
(*Generic plunging orbits - fixed eccentricity*)


\[Delta]EEfuncc[a_,p_,e_,x_]:=Module[{EEg,Lzg},
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];

	((-1+e^2)^2(a^3(-1+EEg^2)-a*Lzg^2))/(2(-3+e)^2*EEg (-a*EEg+Lzg)^2*p^2)
]


\[Delta]Lzfuncc[a_,p_,e_,x_]:=Module[{r1g,EEg,Lzg,\[Delta]EE},
	r1g=p/(1-e);
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];

	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x];

	(a (-a*EEg+Lzg)^2+3EEg(a*EEg-Lzg)r1g^2+EEg*Lzg*r1g^3)/((2a*EEg+Lzg(-2+r1g))r1g^2)+((-2a*Lzg+EEg*r1g^3+a^2*EEg(2+r1g))\[Delta]EE)/(2a*EEg+Lzg(-2+r1g))
]


(* ::Subsection::Closed:: *)
(*Shift ISCO -  fixed eccentricity*)


\[Delta]EEISCOfun[a_,x_]:=Module[{r1g,EEg,Lzg,y},
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,x];
	Lzg=Lzgfun[a,r1g,0,x];
	y=Lzg-a*EEg;

	-((y(4a*y^2-3a*y*Lzg*r1g+EEg*Lzg*r1g^3))/(2a*y^2*r1g^3-6EEg*y*r1g^5+2EEg*Lzg*r1g^6))
]


\[Delta]LzISCOfun[a_,x_]:=Module[{r1g,EEg,Lzg,y},
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,x];
	Lzg=Lzgfun[a,r1g,0,x];
	y=Lzg-a*EEg;

	(-4a^2y^3+3a^3*EEg*y^2*r1g-a*EEg(8a*EEg-7Lzg)y*r1g^3-9EEg^2*y*r1g^5+2EEg^2*Lzg*r1g^6)/(2a*y^2*r1g^3-6EEg*y*r1g^5+2EEg*Lzg*r1g^6)
]


(* ::Section::Closed:: *)
(*Shifts to the radial potential roots*)


\[Delta]r41fun[a_,p_,e_,x_]:=Module[{EEg,Lzg,Lzred},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	Lzred=Lzg-a*EEg;

	-a(a^2(1-EEg^2)+Lzg^2)/(4Lzred^2)
]


\[Delta]r41funcc[a_,p_,e_,x_]:=Module[{EEg,Lzg,Lzred},
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	Lzred=Lzg-a*EEg;

	-a(a^2(1-EEg^2)+Lzg^2)/(4Lzred^2)
]


\[Delta]r41funreg[a_,p_,e_,x_]:=Module[{EEg,Lzg,Lzred},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	Lzred=Lzg-a*EEg;

	-(a^2(1-EEg^2)+Lzg^2)/(4Lzred^2)
]


(* ::Subsection::Closed:: *)
(*Fixed turning points*)


\[Delta]r3funFT[a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,Lzred},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	Lzred=Lzg-a*EEg;

	a (a^2(1-EEg^2)+Lzg^2)/(2Lzred^2)+(4EEg*\[Delta]EE)/(1-EEg^2)^2
]


(* ::Subsection::Closed:: *)
(*Fixed constants of motion*)


\[Delta]r1funFC[a_,p_,e_,x_]:=Module[{r1g,r2g,EEg,Lzg,Lzred},
	r1g=p/(1-e);
	r2g=p/(1+e);
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	Lzred=Lzg-a*EEg;

	-((a*Lzred^2-3EEg*Lzred*r1g^2+EEg*Lzg*r1g^3)/(r1g(3Lzred^2-2a^2(1-EEg^2)r1g-2Lzg^2r1g+5r1g^2-3(1-EEg^2)r1g^3)))
]


\[Delta]r2funFC[a_,p_,e_,x_]:=Module[{r1g,r2g,EEg,Lzg,\[Delta]EE,Lzred},
	r1g=p/(1-e);
	r2g=p/(1+e);
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	Lzred=Lzg-a*EEg;

	-((a*Lzred^2-3EEg*Lzred*r2g^2+EEg*Lzg*r2g^3)/(r2g(3Lzred^2-2a^2(1-EEg^2)r2g-2Lzg^2*r2g+5r2g^2-3(1-EEg^2)r2g^3)))
]


\[Delta]r3funFC[a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,Lzred},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	Lzred=Lzg-a*EEg;

	a (a^2(1-EEg^2)+Lzg^2)/(2Lzred^2)-\[Delta]r1funFC[a,p,e,x]-\[Delta]r2funFC[a,p,e,x]
]


(* ::Subsection::Closed:: *)
(*Fixed frequency *)


\[Delta]turnpointsfunFF[a_,p_,e_,x_]:=Module[{\[Delta]\[CapitalOmega]rFT,\[Delta]\[CapitalOmega]\[Phi]FT,jacob},
	\[Delta]\[CapitalOmega]rFT=(\[Delta]\[CapitalUpsilon]rfunFT[a,p,e,x]-\[CapitalOmega]rgfun[a,p,e,x] \[Delta]\[CapitalUpsilon]tfunFT[a,p,e,x])/(\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x]);
	\[Delta]\[CapitalOmega]\[Phi]FT=(\[Delta]\[CapitalUpsilon]\[Phi]funFT[a,p,e,x]-\[CapitalOmega]\[Phi]gfun[a,p,e,x]\[Delta]\[CapitalUpsilon]tfunFT[a,p,e,x])/(\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x]);
	 
	jacob=-{{d\[CapitalOmega]rgdr1gfun[a,p,e,x],d\[CapitalOmega]rgdr2gfun[a,p,e,x]},{d\[CapitalOmega]\[Phi]gdr1gfun[a,p,e,x],d\[CapitalOmega]\[Phi]gdr2gfun[a,p,e,x]}};

	Inverse[jacob] . {\[Delta]\[CapitalOmega]rFT,\[Delta]\[CapitalOmega]\[Phi]FT}
]


\[Delta]r3funFF[a_,p_,e_,x_,{\[Delta]r1_,\[Delta]r2_}]:=Module[{EEg,Lzg,\[Delta]EE,Lzred},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	Lzred=Lzg-a*EEg;
	(*{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];*)
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;

	a (a^2(1-EEg^2)+Lzg^2)/(2Lzred^2)+(4EEg*\[Delta]EE)/(1-EEg^2)^2-(\[Delta]r1+\[Delta]r2)
]


(* ::Subsection::Closed:: *)
(*Fixed eccentricity*)


\[Delta]r1funFE[a_,p_,e_,x_]:=\[Delta]pfunFE[a,p,e,x]/(1-e)


\[Delta]r2funFE[a_,p_,e_,x_]:=\[Delta]pfunFE[a,p,e,x]/(1+e)


\[Delta]rISCOfun[a_,x_]:=Module[{r1g,EEg,Lzg,y},
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,x];
	Lzg=Lzgfun[a,r1g,0,x];
	y=Lzg-a*EEg;

	1/6((a^3(1-EEg^2)+a*Lzg^2)/y^2-(4EEg*y(4a*y^2-3a*y*Lzg*r1g+EEg*Lzg*r1g^3))/((1-EEg^2)^2*r1g^3(a*y^2-3EEg*y*r1g^2+EEg*Lzg*r1g^3)))
]


(* ::Subsection::Closed:: *)
(*Regular gauge*)


(* ::Subsubsection::Closed:: *)
(*Periodic orbits and plunge related to periodic orbits*)


\[Delta]r1funRPer[a_,p_,e_,x_]:=\[Delta]pfunRPer[a,p,e,x]/(1-e)+p*\[Delta]efunRPer[a,p,e,x]/(1-e)^2


\[Delta]r2funRPer[a_,p_,e_,x_]:=\[Delta]pfunRPer[a,p,e,x]/(1+e)-p*\[Delta]efunRPer[a,p,e,x]/(1+e)^2


(* ::Subsubsection::Closed:: *)
(*Homoclinic orbits and critical plunge *)


\[Delta]r1funRHom[a_,p_,e_,x_]:=\[Delta]pfunRHom[a,p,e,x]/(1-e)+p*\[Delta]efunRHom[a,p,e,x]/(1-e)^2


\[Delta]r2funRHom[a_,p_,e_,x_]:=\[Delta]pfunRHom[a,p,e,x]/(1+e)-p*\[Delta]efunRHom[a,p,e,x]/(1+e)^2


(* ::Subsection::Closed:: *)
(*Shift semilatus rectum and eccentricity - Fixed eccentricity*)


\[Delta]pfunFE[a_,p_,e_,x_]:=Module[{EEg,dEdp,\[Delta]EE,\[Delta]\[Rho]r4},
	EEg=EEgfun[a,p,e,x];
	dEdp=dEEdpfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];

	(2(1-e^2)(2*EEg*\[Delta]EE-(1-EEg^2)^2\[Delta]\[Rho]r4))/(-4dEdp*EEg+4*dEdp*e^2*EEg+3(1-EEg^2)^2-e (1-EEg^2)^2)
]


(* ::Text:: *)
(*For generic plunge*)


\[Delta]pfunFEcc[a_,p_,e_,x_]:=Module[{r1g,\[Rho]rg,EEg,Lzg,Lzred,dEdp,dLzdp,\[Delta]\[Rho]r4,num,den},
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);

	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	Lzred=Lzg-a*EEg;

	dEdp=dEEdpfuncc[p,e];
	dLzdp=dLzdpfuncc[a,p,e,x];
	\[Delta]\[Rho]r4=\[Delta]r41funcc[a,p,e,x];

	num=(1-e)^2/(-3+e)^2(8\[Delta]\[Rho]r4)/\[Rho]rg^2(-2+((1-EEg^2)(2Lzred^2+r1g^2))/r1g+(2(1-EEg^2)^2*Lzred(EEg*r1g^2-a*Lzred))/(EEg(Lzg*r1g-2Lzred)))-(8(1-EEg^2)^2*Lzred^2\[Delta]\[Rho]r4)/r1g^2-(4a (1-EEg^2)^2)/r1g+(8a (1-EEg^2)^2Lzred^3)/((Lzg*r1g-2Lzred)r1g^3)-(4(1-EEg^2)^2*Lzred(EEg Lzg))/(Lzg*r1g-2Lzred)-(1-EEg^2)^2*2\[Delta]\[Rho]r4(2-(1-EEg^2)r1g);

	den=-((8dLzdp (1-EEg^2)^2*Lzred)/r1g)+4dEdp(2EEg+(2(1-EEg^2)(a(1-EEg^2)Lzred-EEg*Lzred^2))/r1g-EEg(1-EEg^2)r1g)+((1-EEg^2)^2(-4Lzred^2+r1g^2(2+(-1+EEg^2)r1g)))/((-1+e)r1g^2);

	num/den
]


(* ::Subsection::Closed:: *)
(*Shift semilatus rectum and eccentricity - Regular gauge*)


(* ::Subsubsection::Closed:: *)
(*Periodic orbits and plunge related to bound orbits*)


\[Delta]pfunRPer[a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,rp,EEg,dEdp,dEde,\[Delta]EE,\[Delta]r41,krg,ellK,ellEc,\[Phi],ellFsqr,ellEsqr},
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	rp=1+Sqrt[1-a^2];
	
	EEg=EEgfun[a,p,e,x];
	dEdp=dEEdpfun[a,p,e,x];
	dEde=dEEdefun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));
	ellK=EllipticK[krg];
	ellEc=EllipticE[krg];
	\[Phi]=ArcSin[Sqrt[(r2g (r3g-rp))/(r3g (r2g-rp))]];
	ellFsqr=Sqrt[(rp(-rp+r1g)(-rp+r2g)(-rp+r3g))]/Sqrt[r2g (r1g-r3g)]EllipticF[\[Phi],krg];
	ellEsqr=Sqrt[(rp(-rp+r1g)(-rp+r2g)(-rp+r3g))]/Sqrt[r2g (r1g-r3g)]EllipticE[\[Phi],krg];
	
	-((a (-1+e^2) (4 dEde (-1+e^2)^2 EEg+(1+(-6+e) e) (-1+EEg^2)^2 p) (r1g-r2g) (r1g-r3g) (ellK r1g^2 r2g r3g^2+(2 ellEc ellFsqr (r1g-r3g)+ellK r1g r3g+2 ellEsqr ellK (-r1g+r3g)) (r2g r3g+r1g (r2g+r3g)) rp-ellK (2 r2g r3g^2+r1g^2 (2 r2g+r3g)+r1g r3g (3 r2g+r3g)) rp^2+ellK (2 r2g r3g+r1g (2 r2g+r3g)) rp^3))/(3 r1g r2g r3g rp ((-3+e-4 dEdp (-1+e^2) EEg-2 (-3+e) EEg^2+(-3+e) EEg^4) p ((-1+e)^2 ellEc ellFsqr r1g r2g (r1g-r3g)-(-1+e)^2 ellEsqr ellK r1g r2g (r1g-r3g)-(-1+e)^2 ellEc ellFsqr r1g (r1g-r3g)^2+(-1+e)^2 ellEsqr ellK r1g (r1g-r3g)^2+(1+e)^2 ellEc ellFsqr r2g (r1g-r3g) r3g-(1+e)^2 ellEsqr ellK r2g (r1g-r3g) r3g+(-1+e)^2 ellK r1g (r1g-r2g) (r1g-rp) rp+(1+e)^2 ellK (r1g-r2g) r3g (r3g-rp) rp)+(4 dEde (-1+e^2)^2 EEg+(1+(-6+e) e) (-1+EEg^2)^2 p) ((ellEc ellFsqr-ellEsqr ellK) (r1g-r3g) ((-1+e) r1g^2-(1+e) r2g r3g-(-1+e) r1g (r2g+r3g))-ellK (r1g-r2g) ((-1+e) r1g^2+(1+e) r3g^2) rp+ellK (r1g-r2g) ((-1+e) r1g+r3g+e r3g) rp^2))))+(4 (-1+e^2) EEg p ((-1+e)^2 ellEc ellFsqr r1g r2g (r1g-r3g)-(-1+e)^2 ellEsqr ellK r1g r2g (r1g-r3g)-(-1+e)^2 ellEc ellFsqr r1g (r1g-r3g)^2+(-1+e)^2 ellEsqr ellK r1g (r1g-r3g)^2+(1+e)^2 ellEc ellFsqr r2g (r1g-r3g) r3g-(1+e)^2 ellEsqr ellK r2g (r1g-r3g) r3g+(-1+e)^2 ellK r1g (r1g-r2g) (r1g-rp) rp+(1+e)^2 ellK (r1g-r2g) r3g (r3g-rp) rp) \[Delta]EE)/((-3+e-4 dEdp (-1+e^2) EEg-2 (-3+e) EEg^2+(-3+e) EEg^4) p ((-1+e)^2 ellEc ellFsqr r1g r2g (r1g-r3g)-(-1+e)^2 ellEsqr ellK r1g r2g (r1g-r3g)-(-1+e)^2 ellEc ellFsqr r1g (r1g-r3g)^2+(-1+e)^2 ellEsqr ellK r1g (r1g-r3g)^2+(1+e)^2 ellEc ellFsqr r2g (r1g-r3g) r3g-(1+e)^2 ellEsqr ellK r2g (r1g-r3g) r3g+(-1+e)^2 ellK r1g (r1g-r2g) (r1g-rp) rp+(1+e)^2 ellK (r1g-r2g) r3g (r3g-rp) rp)+(4 dEde (-1+e^2)^2 EEg+(1+(-6+e) e) (-1+EEg^2)^2 p) ((ellEc ellFsqr-ellEsqr ellK) (r1g-r3g) ((-1+e) r1g^2-(1+e) r2g r3g-(-1+e) r1g (r2g+r3g))-ellK (r1g-r2g) ((-1+e) r1g^2+(1+e) r3g^2) rp+ellK (r1g-r2g) ((-1+e) r1g+r3g+e r3g) rp^2))-(2 (-1+e^2) (r1g-r2g) (r1g-r3g) (-((4 dEde (-1+e^2)^2 EEg+(1+(-6+e) e) (-1+EEg^2)^2 p) (ellEc ellFsqr (r1g-r3g)+ellEsqr ellK (-r1g+r3g)+ellK (r1g-rp) (r3g-rp)))+1/((r1g-r2g) (r1g-r3g)) (-1+EEg^2)^2 p ((-1+e)^2 ellEc ellFsqr r1g r2g (r1g-r3g)-(-1+e)^2 ellEsqr ellK r1g r2g (r1g-r3g)-(-1+e)^2 ellEc ellFsqr r1g (r1g-r3g)^2+(-1+e)^2 ellEsqr ellK r1g (r1g-r3g)^2+(1+e)^2 ellEc ellFsqr r2g (r1g-r3g) r3g-(1+e)^2 ellEsqr ellK r2g (r1g-r3g) r3g+(-1+e)^2 ellK r1g (r1g-r2g) (r1g-rp) rp+(1+e)^2 ellK (r1g-r2g) r3g (r3g-rp) rp)) \[Delta]r41)/((-3+e-4 dEdp (-1+e^2) EEg-2 (-3+e) EEg^2+(-3+e) EEg^4) p ((-1+e)^2 ellEc ellFsqr r1g r2g (r1g-r3g)-(-1+e)^2 ellEsqr ellK r1g r2g (r1g-r3g)-(-1+e)^2 ellEc ellFsqr r1g (r1g-r3g)^2+(-1+e)^2 ellEsqr ellK r1g (r1g-r3g)^2+(1+e)^2 ellEc ellFsqr r2g (r1g-r3g) r3g-(1+e)^2 ellEsqr ellK r2g (r1g-r3g) r3g+(-1+e)^2 ellK r1g (r1g-r2g) (r1g-rp) rp+(1+e)^2 ellK (r1g-r2g) r3g (r3g-rp) rp)+(4 dEde (-1+e^2)^2 EEg+(1+(-6+e) e) (-1+EEg^2)^2 p) ((ellEc ellFsqr-ellEsqr ellK) (r1g-r3g) ((-1+e) r1g^2-(1+e) r2g r3g-(-1+e) r1g (r2g+r3g))-ellK (r1g-r2g) ((-1+e) r1g^2+(1+e) r3g^2) rp+ellK (r1g-r2g) ((-1+e) r1g+r3g+e r3g) rp^2))
]


\[Delta]efunRPer[a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,rp,EEg,dEdp,dEde,\[Delta]EE,\[Delta]r41,krg,ellK,ellEc,\[Phi],ellFsqr,ellEsqr},
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	rp=1+Sqrt[1-a^2];
	
	EEg=EEgfun[a,p,e,x];
	dEdp=dEEdpfun[a,p,e,x];
	dEde=dEEdefun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));
	ellK=EllipticK[krg];
	ellEc=EllipticE[krg];
	\[Phi]=ArcSin[Sqrt[(r2g (r3g-rp))/(r3g (r2g-rp))]];
	ellFsqr=Sqrt[(rp(-rp+r1g)(-rp+r2g)(-rp+r3g))]/Sqrt[r2g (r1g-r3g)]EllipticF[\[Phi],krg];
	ellEsqr=Sqrt[(rp(-rp+r1g)(-rp+r2g)(-rp+r3g))]/Sqrt[r2g (r1g-r3g)]EllipticE[\[Phi],krg];
	
	-((a (-1+e^2) (-4 dEdp EEg+4 dEdp e^2 EEg+3 (-1+EEg^2)^2-e (-1+EEg^2)^2) (r1g-r2g) (r1g-r3g) (ellK r1g^2 r2g r3g^2+(2 ellEc ellFsqr (r1g-r3g)+ellK r1g r3g+2 ellEsqr ellK (-r1g+r3g)) (r2g r3g+r1g (r2g+r3g)) rp-ellK (2 r2g r3g^2+r1g^2 (2 r2g+r3g)+r1g r3g (3 r2g+r3g)) rp^2+ellK (2 r2g r3g+r1g (2 r2g+r3g)) rp^3))/(6 r1g r2g r3g rp (-2 dEde (-1+e^2) EEg ((ellEc ellFsqr-ellEsqr ellK) (r1g-r3g) ((-1+e) r1g^2-(1+e) r2g r3g-(-1+e) r1g (r2g+r3g))-ellK (r1g-r2g) ((-1+e) r1g^2+(1+e) r3g^2) rp+ellK (r1g-r2g) ((-1+e) r1g+r3g+e r3g) rp^2)+p (-ellEc ellFsqr (r1g-r3g) ((2 dEdp (-1+e)^2 EEg-(-1+EEg^2)^2) r1g (r1g-r2g)+(-2 dEdp (-1+e)^2 EEg+(-1+EEg^2)^2) r1g r3g-2 (dEdp (1+e)^2 EEg-(-1+EEg^2)^2) r2g r3g)+ellK (ellEsqr (r1g-r3g) ((2 dEdp (-1+e)^2 EEg-(-1+EEg^2)^2) r1g (r1g-r2g)+(-2 dEdp (-1+e)^2 EEg+(-1+EEg^2)^2) r1g r3g-2 (dEdp (1+e)^2 EEg-(-1+EEg^2)^2) r2g r3g)+(r1g-r2g) rp ((2 dEdp (-1+e)^2 EEg-(-1+EEg^2)^2) r1g^2+2 (dEdp (1+e)^2 EEg-(-1+EEg^2)^2) r3g (r3g-rp)+(-2 dEdp (-1+e)^2 EEg+(-1+EEg^2)^2) r1g rp))))))-(2 (-1+e^2) EEg (-((-1+e) ellEc ellFsqr r1g (r1g-r2g) (r1g-r3g))+(-1+e) ellEsqr ellK r1g (r1g-r2g) (r1g-r3g)+ellEc ellFsqr ((-1+e) r1g+r2g+e r2g) (r1g-r3g) r3g-ellEsqr ellK ((-1+e) r1g+r2g+e r2g) (r1g-r3g) r3g+ellK (r1g-r2g) rp ((-1+e) r1g^2+(1+e) r3g (r3g-rp)+r1g (rp-e rp))) \[Delta]EE)/(-2 dEde (-1+e^2) EEg ((ellEc ellFsqr-ellEsqr ellK) (r1g-r3g) ((-1+e) r1g^2-(1+e) r2g r3g-(-1+e) r1g (r2g+r3g))-ellK (r1g-r2g) ((-1+e) r1g^2+(1+e) r3g^2) rp+ellK (r1g-r2g) ((-1+e) r1g+r3g+e r3g) rp^2)+p (-ellEc ellFsqr (r1g-r3g) ((2 dEdp (-1+e)^2 EEg-(-1+EEg^2)^2) r1g (r1g-r2g)+(-2 dEdp (-1+e)^2 EEg+(-1+EEg^2)^2) r1g r3g-2 (dEdp (1+e)^2 EEg-(-1+EEg^2)^2) r2g r3g)+ellK (ellEsqr (r1g-r3g) ((2 dEdp (-1+e)^2 EEg-(-1+EEg^2)^2) r1g (r1g-r2g)+(-2 dEdp (-1+e)^2 EEg+(-1+EEg^2)^2) r1g r3g-2 (dEdp (1+e)^2 EEg-(-1+EEg^2)^2) r2g r3g)+(r1g-r2g) rp ((2 dEdp (-1+e)^2 EEg-(-1+EEg^2)^2) r1g^2+2 (dEdp (1+e)^2 EEg-(-1+EEg^2)^2) r3g (r3g-rp)+(-2 dEdp (-1+e)^2 EEg+(-1+EEg^2)^2) r1g rp))))-((-1+e^2) (-4 dEdp (-1+e^2) EEg ellEc ellFsqr (r1g-r2g) (r1g-r3g)^2+4 dEdp (-1+e^2) EEg ellEsqr ellK (r1g-r2g) (r1g-r3g)^2+2 ellEc ellFsqr (r1g-r3g) ((-2+e) r1g^2-2 r2g r3g-(-2+e) r1g (r2g+r3g))+4 dEdp (-1+e^2) EEg ellK (r1g-r2g) (r1g-r3g) (r1g-rp) (-r3g+rp)-ellK (2 ellEsqr (r1g-r3g) ((-2+e) r1g^2-2 r2g r3g-(-2+e) r1g (r2g+r3g))+(r1g-r2g) (-((-3+e) r1g (r1g-r3g) r3g)+2 ((-2+e) r1g^2+2 r3g^2) rp-2 ((-2+e) r1g+2 r3g) rp^2))-EEg^4 (-2 ellEc ellFsqr (r1g-r3g) ((-2+e) r1g^2-2 r2g r3g-(-2+e) r1g (r2g+r3g))+ellK (2 ellEsqr (r1g-r3g) ((-2+e) r1g^2-2 r2g r3g-(-2+e) r1g (r2g+r3g))+(r1g-r2g) (-((-3+e) r1g (r1g-r3g) r3g)+2 ((-2+e) r1g^2+2 r3g^2) rp-2 ((-2+e) r1g+2 r3g) rp^2)))-2 EEg^2 (2 ellEc ellFsqr (r1g-r3g) ((-2+e) r1g^2-2 r2g r3g-(-2+e) r1g (r2g+r3g))+ellK (-2 ellEsqr (r1g-r3g) ((-2+e) r1g^2-2 r2g r3g-(-2+e) r1g (r2g+r3g))+(r1g-r2g) ((-3+e) r1g (r1g-r3g) r3g-2 ((-2+e) r1g^2+2 r3g^2) rp+2 ((-2+e) r1g+2 r3g) rp^2)))) \[Delta]r41)/(-2 dEde (-1+e^2) EEg ((ellEc ellFsqr-ellEsqr ellK) (r1g-r3g) ((-1+e) r1g^2-(1+e) r2g r3g-(-1+e) r1g (r2g+r3g))-ellK (r1g-r2g) ((-1+e) r1g^2+(1+e) r3g^2) rp+ellK (r1g-r2g) ((-1+e) r1g+r3g+e r3g) rp^2)+p (-ellEc ellFsqr (r1g-r3g) ((2 dEdp (-1+e)^2 EEg-(-1+EEg^2)^2) r1g (r1g-r2g)+(-2 dEdp (-1+e)^2 EEg+(-1+EEg^2)^2) r1g r3g-2 (dEdp (1+e)^2 EEg-(-1+EEg^2)^2) r2g r3g)+ellK (ellEsqr (r1g-r3g) ((2 dEdp (-1+e)^2 EEg-(-1+EEg^2)^2) r1g (r1g-r2g)+(-2 dEdp (-1+e)^2 EEg+(-1+EEg^2)^2) r1g r3g-2 (dEdp (1+e)^2 EEg-(-1+EEg^2)^2) r2g r3g)+(r1g-r2g) rp ((2 dEdp (-1+e)^2 EEg-(-1+EEg^2)^2) r1g^2+2 (dEdp (1+e)^2 EEg-(-1+EEg^2)^2) r3g (r3g-rp)+(-2 dEdp (-1+e)^2 EEg+(-1+EEg^2)^2) r1g rp))))
]	


(* ::Subsubsection::Closed:: *)
(*Homoclinic orbits and critical plunge *)


\[Delta]pfunRHom[a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,rp,EEg,dEdp,dEde,\[Delta]EE,\[Delta]r41,krg,ellK,ellEc,\[Phi],ellFsqr,ellEsqr},
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g)/2;
	rp=1+Sqrt[1-a^2];
	
	EEg=EEgfun[a,p,e,x];
	dEdp=dEEdpfun[a,p,e,x];
	dEde=dEEdefun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));
	ellK=EllipticK[krg];
	ellEc=EllipticE[krg];
	\[Phi]=ArcSin[Sqrt[(r2g (r3g-rp))/(r3g (r2g-rp))]];
	ellFsqr=Sqrt[(rp(-rp+r1g)(-rp+r2g)(-rp+r3g))]EllipticF[\[Phi],krg];
	ellEsqr=Sqrt[(rp(-rp+r1g)(-rp+r2g)(-rp+r3g))]EllipticE[\[Phi],krg];
	
	(a (4 dEde (-1+e^2)^2 EEg+(1-6 e+e^2) (-1+EEg^2)^2 p) (r1g-r2g) (r1g-rp) (r2g-rp) (2 r2g rp+r1g (r2g+3 rp)))/(6 r1g r2g rp^2 (p ((2 dEdp (-1+e)^2 EEg-(-1+EEg^2)^2) r1g^2+2 (dEdp (1+e)^2 EEg-(-1+EEg^2)^2) r2g (r2g-rp)+(-2 dEdp (-1+e)^2 EEg+(-1+EEg^2)^2) r1g rp)+2 dEde (-1+e^2) EEg ((-1+e) r1g^2+(1+e) r2g (r2g-rp)+r1g (rp-e rp))))-(2 EEg p ((-1+e)^2 r1g^2+(1+e)^2 r2g (r2g-rp)-(-1+e)^2 r1g rp) \[Delta]EE)/(p ((2 dEdp (-1+e)^2 EEg-(-1+EEg^2)^2) r1g^2+2 (dEdp (1+e)^2 EEg-(-1+EEg^2)^2) r2g (r2g-rp)+(-2 dEdp (-1+e)^2 EEg+(-1+EEg^2)^2) r1g rp)+2 dEde (-1+e^2) EEg ((-1+e) r1g^2+(1+e) r2g (r2g-rp)+r1g (rp-e rp)))+((4 dEde (-1+e^2)^2 EEg (r1g-r2g) (r1g-rp) (r2g-rp)+(-1+EEg^2)^2 p (-((1-6 e+e^2) r1g r2g^2)+2 (1-4 e+e^2) r1g rp^2+8 e r2g rp (-r2g+rp)+r1g^2 ((1-6 e+e^2) r2g-2 (1-4 e+e^2) rp))) \[Delta]r41)/(rp (p ((-2 dEdp (-1+e)^2 EEg+(-1+EEg^2)^2) r1g^2-2 (dEdp (1+e)^2 EEg-(-1+EEg^2)^2) r2g (r2g-rp)+(2 dEdp (-1+e)^2 EEg-(-1+EEg^2)^2) r1g rp)-2 dEde (-1+e^2) EEg ((-1+e) r1g^2+(1+e) r2g (r2g-rp)+r1g (rp-e rp))))
]


\[Delta]efunRHom[a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,rp,EEg,dEdp,dEde,\[Delta]EE,\[Delta]r41,krg,ellK,ellEc,\[Phi],ellFsqr,ellEsqr},
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g)/2;
	rp=1+Sqrt[1-a^2];
	
	EEg=EEgfun[a,p,e,x];
	dEdp=dEEdpfun[a,p,e,x];
	dEde=dEEdefun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));
	ellK=EllipticK[krg];
	ellEc=EllipticE[krg];
	\[Phi]=ArcSin[Sqrt[(r2g (r3g-rp))/(r3g (r2g-rp))]];
	ellFsqr=Sqrt[(rp(-rp+r1g)(-rp+r2g)(-rp+r3g))]EllipticF[\[Phi],krg];
	ellEsqr=Sqrt[(rp(-rp+r1g)(-rp+r2g)(-rp+r3g))]EllipticE[\[Phi],krg];
	
	-((a (-1+e^2) (-4 dEdp EEg+4 dEdp e^2 EEg+3 (-1+EEg^2)^2-e (-1+EEg^2)^2) (r1g-r2g) (r1g-rp) (r2g-rp) (2 r2g rp+r1g (r2g+3 rp)))/(6 r1g r2g rp^2 (p ((2 dEdp (-1+e)^2 EEg-(-1+EEg^2)^2) r1g^2+2 (dEdp (1+e)^2 EEg-(-1+EEg^2)^2) r2g (r2g-rp)+(-2 dEdp (-1+e)^2 EEg+(-1+EEg^2)^2) r1g rp)+2 dEde (-1+e^2) EEg ((-1+e) r1g^2+(1+e) r2g (r2g-rp)+r1g (rp-e rp)))))-(2 (-1+e^2) EEg ((-1+e) r1g^2+(1+e) r2g (r2g-rp)+r1g (rp-e rp)) \[Delta]EE)/(p ((2 dEdp (-1+e)^2 EEg-(-1+EEg^2)^2) r1g^2+2 (dEdp (1+e)^2 EEg-(-1+EEg^2)^2) r2g (r2g-rp)+(-2 dEdp (-1+e)^2 EEg+(-1+EEg^2)^2) r1g rp)+2 dEde (-1+e^2) EEg ((-1+e) r1g^2+(1+e) r2g (r2g-rp)+r1g (rp-e rp)))-((-1+e^2) (4 (dEdp (-1+e^2) EEg+(-1+EEg^2)^2) r2g (r2g-rp) rp+r1g^2 ((3-e+4 dEdp (-1+e^2) EEg+2 (-3+e) EEg^2-(-3+e) EEg^4) r2g+2 (-2+e-2 dEdp (-1+e^2) EEg-2 (-2+e) EEg^2+(-2+e) EEg^4) rp)+r1g ((-3+e-4 dEdp (-1+e^2) EEg-2 (-3+e) EEg^2+(-3+e) EEg^4) r2g^2+2 (2-e+2 dEdp (-1+e^2) EEg+2 (-2+e) EEg^2-(-2+e) EEg^4) rp^2)) \[Delta]r41)/(rp (p ((-2 dEdp (-1+e)^2 EEg+(-1+EEg^2)^2) r1g^2-2 (dEdp (1+e)^2 EEg-(-1+EEg^2)^2) r2g (r2g-rp)+(2 dEdp (-1+e)^2 EEg-(-1+EEg^2)^2) r1g rp)-2 dEde (-1+e^2) EEg ((-1+e) r1g^2+(1+e) r2g (r2g-rp)+r1g (rp-e rp))))
]


(* ::Subsubsection::Closed:: *)
(*Generic plunge*)


\[Delta]pfunRPlunge[a_,p_,e_,x_]:=Module[{r1g,\[Rho]rg,\[Rho]ig,Ar1g,B0,rp,EEg,Lzg,dEdp,dEde,dLzdp,dLzde,\[Delta]EE,\[Delta]r41,krg,\[Phi],ellE},
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];
	Ar1g=Sqrt[\[Rho]ig^2+(r1g-\[Rho]rg)^2];
	B0=Sqrt[\[Rho]ig^2+\[Rho]rg^2];
	rp=1+Sqrt[1-a^2];
	
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	dEdp=dEEdpfuncc[p,e];
	dEde=dEEdefuncc[p,e];
	dLzdp=dLzdpfuncc[a,p,e,x];
	dLzde=dLzdefuncc[a,p,e,x];
	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x];
	\[Delta]r41=\[Delta]r41funcc[a,p,e,x];
		
	krg=(-(Ar1g-B0)^2+r1g^2)/(4Ar1g*B0);
	\[Phi]=ArcCos[(-Ar1g*rp+B0(r1g-rp))/(Ar1g*rp+B0(r1g-rp))];
	ellE=EllipticE[\[Pi]-\[Phi],krg];
	
	-((\[Delta]r41 (-(((-((4 a dEde (a EEg-Lzg))/(-1+EEg^2))+(4 dLzde (a EEg-Lzg))/(-1+EEg^2)+(4 dEde EEg (-a EEg+Lzg)^2)/(-1+EEg^2)^2+(2 (-a EEg+Lzg)^2 p)/((-1+e)^2 (-1+EEg^2) r1g)+(2 dEde EEg r1g (2+(-1+EEg^2) r1g))/(-1+EEg^2)^3+(p r1g (2+(-1+EEg^2) r1g))/((-1+e)^2 (2-2 EEg^2))) (r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2])/(B0 r1g^2 (B0 (r1g-rp)+Ar1g rp) \[Rho]ig))+(Ar1g^2 ellE (-((4 a dEde (a EEg-Lzg))/(-1+EEg^2))+(4 dLzde (a EEg-Lzg))/(-1+EEg^2)+(4 dEde EEg (-a EEg+Lzg)^2)/(-1+EEg^2)^2+(2 (-a EEg+Lzg)^2 p)/((-1+e)^2 (-1+EEg^2) r1g)+(2 dEde EEg r1g (2+(-1+EEg^2) r1g))/(-1+EEg^2)^3+(p r1g (2+(-1+EEg^2) r1g))/((-1+e)^2 (2-2 EEg^2))) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) r1g^2 \[Rho]ig)))/(1/(4 (-1+EEg^2)^3 r1g \[Rho]ig) (-8 a dEdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dLzdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dEdp EEg (-1+EEg^2) (-a EEg+Lzg)^2-(4 (-1+EEg^2)^2 (-a EEg+Lzg)^2)/((-1+e) r1g)+4 dEdp EEg r1g (2+(-1+EEg^2) r1g)+((-1+EEg^2)^2 r1g (2+(-1+EEg^2) r1g))/(-1+e)) (p/(-1+e)^2-(p (r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (Ar1g B0+rp (r1g+rp-2 \[Rho]rg)))/(Ar1g (-1+e)^2 r1g (B0 (r1g-rp)+Ar1g rp) Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2])+(B0^2 ellE p Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (-1+e)^2 r1g)-(ellE p (r1g-2 \[Rho]rg) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (1+e)^2)+(p r1g (r1g-rp) rp Sqrt[rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2] (2 Ar1g B0 (rp-\[Rho]rg)-(r1g-2 \[Rho]rg) (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)))/(Ar1g B0 (1+e)^2 Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (B0^2 (r1g-2 rp)-rp^2 (r1g-2 \[Rho]rg))))-1/(2 r1g \[Rho]ig) (-((4 a dEde (a EEg-Lzg))/(-1+EEg^2))+(4 dLzde (a EEg-Lzg))/(-1+EEg^2)+(4 dEde EEg (-a EEg+Lzg)^2)/(-1+EEg^2)^2+(2 (-a EEg+Lzg)^2 p)/((-1+e)^2 (-1+EEg^2) r1g)+(2 dEde EEg r1g (2+(-1+EEg^2) r1g))/(-1+EEg^2)^3+(p r1g (2+(-1+EEg^2) r1g))/((-1+e)^2 (2-2 EEg^2))) (1/(1-e)+((r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (Ar1g B0+rp (r1g+rp-2 \[Rho]rg)))/(Ar1g (-1+e) r1g (B0 (r1g-rp)+Ar1g rp) Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2])+(B0^2 ellE Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (r1g-e r1g))+(ellE (r1g-2 \[Rho]rg) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (1+e))+(r1g (r1g-rp) rp Sqrt[rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2] (2 Ar1g B0 (-rp+\[Rho]rg)+(r1g-2 \[Rho]rg) (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)))/(Ar1g B0 (1+e) Sqrt[r1g (B0^2 (r1g-rp)+rp (Ar1g^2+r1g (-r1g+rp)))] (B0^2 (r1g-2 rp)-rp^2 (r1g-2 \[Rho]rg))))))-(-((a (-((4 a dEde (a EEg-Lzg))/(-1+EEg^2))+(4 dLzde (a EEg-Lzg))/(-1+EEg^2)+(4 dEde EEg (-a EEg+Lzg)^2)/(-1+EEg^2)^2+(2 (-a EEg+Lzg)^2 p)/((-1+e)^2 (-1+EEg^2) r1g)+(2 dEde EEg r1g (2+(-1+EEg^2) r1g))/(-1+EEg^2)^3+(p r1g (2+(-1+EEg^2) r1g))/((-1+e)^2 (2-2 EEg^2))) (-r1g+rp))/(6 r1g^2 rp \[Rho]ig))+(a (-((4 a dEde (a EEg-Lzg))/(-1+EEg^2))+(4 dLzde (a EEg-Lzg))/(-1+EEg^2)+(4 dEde EEg (-a EEg+Lzg)^2)/(-1+EEg^2)^2+(2 (-a EEg+Lzg)^2 p)/((-1+e)^2 (-1+EEg^2) r1g)+(2 dEde EEg r1g (2+(-1+EEg^2) r1g))/(-1+EEg^2)^3+(p r1g (2+(-1+EEg^2) r1g))/((-1+e)^2 (2-2 EEg^2))) (r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2] (B0^2+r1g \[Rho]rg))/(3 B0^3 r1g^3 (B0 (r1g-rp)+Ar1g rp) \[Rho]ig)-(a Sqrt[Ar1g B0] ellE (-((4 a dEde (a EEg-Lzg))/(-1+EEg^2))+(4 dLzde (a EEg-Lzg))/(-1+EEg^2)+(4 dEde EEg (-a EEg+Lzg)^2)/(-1+EEg^2)^2+(2 (-a EEg+Lzg)^2 p)/((-1+e)^2 (-1+EEg^2) r1g)+(2 dEde EEg r1g (2+(-1+EEg^2) r1g))/(-1+EEg^2)^3+(p r1g (2+(-1+EEg^2) r1g))/((-1+e)^2 (2-2 EEg^2))) (B0^2+2 r1g \[Rho]rg) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/(3 B0^4 r1g^3 \[Rho]ig)-(a (-((4 a dEde (a EEg-Lzg))/(-1+EEg^2))+(4 dLzde (a EEg-Lzg))/(-1+EEg^2)+(4 dEde EEg (-a EEg+Lzg)^2)/(-1+EEg^2)^2+(2 (-a EEg+Lzg)^2 p)/((-1+e)^2 (-1+EEg^2) r1g)+(2 dEde EEg r1g (2+(-1+EEg^2) r1g))/(-1+EEg^2)^3+(p r1g (2+(-1+EEg^2) r1g))/((-1+e)^2 (2-2 EEg^2))) (r1g-rp) rp Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2] (-((4 Ar1g krg (Ar1g rp+B0 (-r1g+rp)) \[Rho]rg)/(r1g (B0 (r1g-rp)+Ar1g rp)))+(-\[Rho]ig^2+\[Rho]rg^2)/B0))/(6 B0^3 r1g Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] \[Rho]ig)+1/(2 r1g \[Rho]ig) (-((6 EEg (a EEg-Lzg))/(-1+EEg^2))-(2 a^2 (-1+e^2)^2 (a^2 (-1+EEg^2)-Lzg^2))/((-3+e)^2 EEg (-1+EEg^2) (a EEg-Lzg) p^2)+(2 (-1+e^2)^2 (a^3 (-1+EEg^2)-a Lzg^2))/((-3+e)^2 (-1+EEg^2)^2 p^2)+a^3/r1g-a r1g+a (2/(-1+EEg^2)+r1g)+(a Lzg^2)/(r1g-EEg^2 r1g)+(a (a^2 (-1+EEg^2)-Lzg^2) r1g (2+(-1+EEg^2) r1g))/(2 (-1+EEg^2) (-a EEg+Lzg)^2)-((a^3 (-1+EEg^2)-a Lzg^2) r1g (2+(-1+EEg^2) r1g))/(4 (-1+EEg^2) (-a EEg+Lzg)^2)+((-1+e^2)^2 (a^3 (-1+EEg^2)-a Lzg^2) r1g (2+(-1+EEg^2) r1g))/((-3+e)^2 (-1+EEg^2)^3 (-a EEg+Lzg)^2 p^2)+(4 (a EEg-Lzg) (a^3 EEg^2+EEg Lzg (-3+r1g) r1g^2+((-1+e^2)^2 (a^3 (-1+EEg^2)-a Lzg^2) r1g^5)/(2 (-3+e)^2 (-a EEg+Lzg)^2 p^2)+a (Lzg^2+3 EEg^2 r1g^2+((-1+e^2)^2 Lzg (-a^3 (-1+EEg^2)+a Lzg^2) r1g^2)/((-3+e)^2 EEg (-a EEg+Lzg)^2 p^2))+1/2 a^2 (-4 EEg Lzg+((-1+e^2)^2 (a^3 (-1+EEg^2)-a Lzg^2) r1g^2 (2+r1g))/((-3+e)^2 (-a EEg+Lzg)^2 p^2))))/((-1+EEg^2) (2 a EEg+Lzg (-2+r1g)) r1g^2)) (p/(-1+e)^2-(p (r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (Ar1g B0+rp (r1g+rp-2 \[Rho]rg)))/(Ar1g (-1+e)^2 r1g (B0 (r1g-rp)+Ar1g rp) Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2])+(B0^2 ellE p Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (-1+e)^2 r1g)-(ellE p (r1g-2 \[Rho]rg) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (1+e)^2)+(p r1g (r1g-rp) rp Sqrt[rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2] (2 Ar1g B0 (rp-\[Rho]rg)-(r1g-2 \[Rho]rg) (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)))/(Ar1g B0 (1+e)^2 Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (B0^2 (r1g-2 rp)-rp^2 (r1g-2 \[Rho]rg)))))/(1/(4 (-1+EEg^2)^3 r1g \[Rho]ig) (-8 a dEdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dLzdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dEdp EEg (-1+EEg^2) (-a EEg+Lzg)^2-(4 (-1+EEg^2)^2 (-a EEg+Lzg)^2)/((-1+e) r1g)+4 dEdp EEg r1g (2+(-1+EEg^2) r1g)+((-1+EEg^2)^2 r1g (2+(-1+EEg^2) r1g))/(-1+e)) (p/(-1+e)^2-(p (r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (Ar1g B0+rp (r1g+rp-2 \[Rho]rg)))/(Ar1g (-1+e)^2 r1g (B0 (r1g-rp)+Ar1g rp) Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2])+(B0^2 ellE p Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (-1+e)^2 r1g)-(ellE p (r1g-2 \[Rho]rg) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (1+e)^2)+(p r1g (r1g-rp) rp Sqrt[rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2] (2 Ar1g B0 (rp-\[Rho]rg)-(r1g-2 \[Rho]rg) (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)))/(Ar1g B0 (1+e)^2 Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (B0^2 (r1g-2 rp)-rp^2 (r1g-2 \[Rho]rg))))-1/(2 r1g \[Rho]ig) (-((4 a dEde (a EEg-Lzg))/(-1+EEg^2))+(4 dLzde (a EEg-Lzg))/(-1+EEg^2)+(4 dEde EEg (-a EEg+Lzg)^2)/(-1+EEg^2)^2+(2 (-a EEg+Lzg)^2 p)/((-1+e)^2 (-1+EEg^2) r1g)+(2 dEde EEg r1g (2+(-1+EEg^2) r1g))/(-1+EEg^2)^3+(p r1g (2+(-1+EEg^2) r1g))/((-1+e)^2 (2-2 EEg^2))) (1/(1-e)+((r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (Ar1g B0+rp (r1g+rp-2 \[Rho]rg)))/(Ar1g (-1+e) r1g (B0 (r1g-rp)+Ar1g rp) Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2])+(B0^2 ellE Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (r1g-e r1g))+(ellE (r1g-2 \[Rho]rg) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (1+e))+(r1g (r1g-rp) rp Sqrt[rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2] (2 Ar1g B0 (-rp+\[Rho]rg)+(r1g-2 \[Rho]rg) (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)))/(Ar1g B0 (1+e) Sqrt[r1g (B0^2 (r1g-rp)+rp (Ar1g^2+r1g (-r1g+rp)))] (B0^2 (r1g-2 rp)-rp^2 (r1g-2 \[Rho]rg)))))
]	


\[Delta]efunRPlunge[a_,p_,e_,x_]:=Module[{r1g,\[Rho]rg,\[Rho]ig,Ar1g,B0,rp,EEg,Lzg,dEdp,dEde,dLzdp,dLzde,\[Delta]EE,\[Delta]r41,krg,\[Phi],ellE},
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];
	Ar1g=Sqrt[\[Rho]ig^2+(r1g-\[Rho]rg)^2];
	B0=Sqrt[\[Rho]ig^2+\[Rho]rg^2];
	rp=1+Sqrt[1-a^2];
	
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	dEdp=dEEdpfuncc[p,e];
	dEde=dEEdefuncc[p,e];
	dLzdp=dLzdpfuncc[a,p,e,x];
	dLzde=dLzdefuncc[a,p,e,x];
	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x];
	\[Delta]r41=\[Delta]r41funcc[a,p,e,x];
		
	krg=(-(Ar1g-B0)^2+r1g^2)/(4Ar1g*B0);
	\[Phi]=ArcCos[(-Ar1g*rp+B0(r1g-rp))/(Ar1g*rp+B0(r1g-rp))];
	ellE=EllipticE[\[Pi]-\[Phi],krg];
	
	-((-((6 EEg (a EEg-Lzg))/(-1+EEg^2))-(2 a^2 (-1+e^2)^2 (a^2 (-1+EEg^2)-Lzg^2))/((-3+e)^2 EEg (-1+EEg^2) (a EEg-Lzg) p^2)+(2 (-1+e^2)^2 (a^3 (-1+EEg^2)-a Lzg^2))/((-3+e)^2 (-1+EEg^2)^2 p^2)+a^3/r1g-a r1g+a (2/(-1+EEg^2)+r1g)+(a Lzg^2)/(r1g-EEg^2 r1g)+(a (a^2 (-1+EEg^2)-Lzg^2) r1g (2+(-1+EEg^2) r1g))/(2 (-1+EEg^2) (-a EEg+Lzg)^2)-((a^3 (-1+EEg^2)-a Lzg^2) r1g (2+(-1+EEg^2) r1g))/(4 (-1+EEg^2) (-a EEg+Lzg)^2)+((-1+e^2)^2 (a^3 (-1+EEg^2)-a Lzg^2) r1g (2+(-1+EEg^2) r1g))/((-3+e)^2 (-1+EEg^2)^3 (-a EEg+Lzg)^2 p^2)+(4 (a EEg-Lzg) (a^3 EEg^2+EEg Lzg (-3+r1g) r1g^2+((-1+e^2)^2 (a^3 (-1+EEg^2)-a Lzg^2) r1g^5)/(2 (-3+e)^2 (-a EEg+Lzg)^2 p^2)+a (Lzg^2+3 EEg^2 r1g^2+((-1+e^2)^2 Lzg (-a^3 (-1+EEg^2)+a Lzg^2) r1g^2)/((-3+e)^2 EEg (-a EEg+Lzg)^2 p^2))+1/2 a^2 (-4 EEg Lzg+((-1+e^2)^2 (a^3 (-1+EEg^2)-a Lzg^2) r1g^2 (2+r1g))/((-3+e)^2 (-a EEg+Lzg)^2 p^2))))/((-1+EEg^2) (2 a EEg+Lzg (-2+r1g)) r1g^2))/(-((4 a dEde (a EEg-Lzg))/(-1+EEg^2))+(4 dLzde (a EEg-Lzg))/(-1+EEg^2)+(4 dEde EEg (-a EEg+Lzg)^2)/(-1+EEg^2)^2+(2 (-a EEg+Lzg)^2 p)/((-1+e)^2 (-1+EEg^2) r1g)+(2 dEde EEg r1g (2+(-1+EEg^2) r1g))/(-1+EEg^2)^3+(p r1g (2+(-1+EEg^2) r1g))/((-1+e)^2 (2-2 EEg^2))))-(a (-8 a dEdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dLzdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dEdp EEg (-1+EEg^2) (-a EEg+Lzg)^2-(4 (-1+EEg^2)^2 (-a EEg+Lzg)^2)/((-1+e) r1g)+4 dEdp EEg r1g (2+(-1+EEg^2) r1g)+((-1+EEg^2)^2 r1g (2+(-1+EEg^2) r1g))/(-1+e)) (-r1g+rp))/(12 (-1+EEg^2)^3 r1g^2 rp \[Rho]ig (1/(4 (-1+EEg^2)^3 r1g \[Rho]ig) (-8 a dEdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dLzdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dEdp EEg (-1+EEg^2) (-a EEg+Lzg)^2-(4 (-1+EEg^2)^2 (-a EEg+Lzg)^2)/((-1+e) r1g)+4 dEdp EEg r1g (2+(-1+EEg^2) r1g)+((-1+EEg^2)^2 r1g (2+(-1+EEg^2) r1g))/(-1+e)) (p/(-1+e)^2-(p (r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (Ar1g B0+rp (r1g+rp-2 \[Rho]rg)))/(Ar1g (-1+e)^2 r1g (B0 (r1g-rp)+Ar1g rp) Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2])+(B0^2 ellE p Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (-1+e)^2 r1g)-(ellE p (r1g-2 \[Rho]rg) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (1+e)^2)+(p r1g (r1g-rp) rp Sqrt[rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2] (2 Ar1g B0 (rp-\[Rho]rg)-(r1g-2 \[Rho]rg) (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)))/(Ar1g B0 (1+e)^2 Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (B0^2 (r1g-2 rp)-rp^2 (r1g-2 \[Rho]rg))))-1/(2 r1g \[Rho]ig) (-((4 a dEde (a EEg-Lzg))/(-1+EEg^2))+(4 dLzde (a EEg-Lzg))/(-1+EEg^2)+(4 dEde EEg (-a EEg+Lzg)^2)/(-1+EEg^2)^2+(2 (-a EEg+Lzg)^2 p)/((-1+e)^2 (-1+EEg^2) r1g)+(2 dEde EEg r1g (2+(-1+EEg^2) r1g))/(-1+EEg^2)^3+(p r1g (2+(-1+EEg^2) r1g))/((-1+e)^2 (2-2 EEg^2))) (1/(1-e)+((r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (Ar1g B0+rp (r1g+rp-2 \[Rho]rg)))/(Ar1g (-1+e) r1g (B0 (r1g-rp)+Ar1g rp) Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2])+(B0^2 ellE Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (r1g-e r1g))+(ellE (r1g-2 \[Rho]rg) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (1+e))+(r1g (r1g-rp) rp Sqrt[rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2] (2 Ar1g B0 (-rp+\[Rho]rg)+(r1g-2 \[Rho]rg) (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)))/(Ar1g B0 (1+e) Sqrt[r1g (B0^2 (r1g-rp)+rp (Ar1g^2+r1g (-r1g+rp)))] (B0^2 (r1g-2 rp)-rp^2 (r1g-2 \[Rho]rg))))))+(a (-8 a dEdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dLzdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dEdp EEg (-1+EEg^2) (-a EEg+Lzg)^2-(4 (-1+EEg^2)^2 (-a EEg+Lzg)^2)/((-1+e) r1g)+4 dEdp EEg r1g (2+(-1+EEg^2) r1g)+((-1+EEg^2)^2 r1g (2+(-1+EEg^2) r1g))/(-1+e)) (r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2] (B0^2+r1g \[Rho]rg))/(6 B0^3 (-1+EEg^2)^3 r1g^3 (B0 (r1g-rp)+Ar1g rp) \[Rho]ig (1/(4 (-1+EEg^2)^3 r1g \[Rho]ig) (-8 a dEdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dLzdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dEdp EEg (-1+EEg^2) (-a EEg+Lzg)^2-(4 (-1+EEg^2)^2 (-a EEg+Lzg)^2)/((-1+e) r1g)+4 dEdp EEg r1g (2+(-1+EEg^2) r1g)+((-1+EEg^2)^2 r1g (2+(-1+EEg^2) r1g))/(-1+e)) (p/(-1+e)^2-(p (r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (Ar1g B0+rp (r1g+rp-2 \[Rho]rg)))/(Ar1g (-1+e)^2 r1g (B0 (r1g-rp)+Ar1g rp) Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2])+(B0^2 ellE p Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (-1+e)^2 r1g)-(ellE p (r1g-2 \[Rho]rg) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (1+e)^2)+(p r1g (r1g-rp) rp Sqrt[rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2] (2 Ar1g B0 (rp-\[Rho]rg)-(r1g-2 \[Rho]rg) (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)))/(Ar1g B0 (1+e)^2 Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (B0^2 (r1g-2 rp)-rp^2 (r1g-2 \[Rho]rg))))-1/(2 r1g \[Rho]ig) (-((4 a dEde (a EEg-Lzg))/(-1+EEg^2))+(4 dLzde (a EEg-Lzg))/(-1+EEg^2)+(4 dEde EEg (-a EEg+Lzg)^2)/(-1+EEg^2)^2+(2 (-a EEg+Lzg)^2 p)/((-1+e)^2 (-1+EEg^2) r1g)+(2 dEde EEg r1g (2+(-1+EEg^2) r1g))/(-1+EEg^2)^3+(p r1g (2+(-1+EEg^2) r1g))/((-1+e)^2 (2-2 EEg^2))) (1/(1-e)+((r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (Ar1g B0+rp (r1g+rp-2 \[Rho]rg)))/(Ar1g (-1+e) r1g (B0 (r1g-rp)+Ar1g rp) Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2])+(B0^2 ellE Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (r1g-e r1g))+(ellE (r1g-2 \[Rho]rg) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (1+e))+(r1g (r1g-rp) rp Sqrt[rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2] (2 Ar1g B0 (-rp+\[Rho]rg)+(r1g-2 \[Rho]rg) (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)))/(Ar1g B0 (1+e) Sqrt[r1g (B0^2 (r1g-rp)+rp (Ar1g^2+r1g (-r1g+rp)))] (B0^2 (r1g-2 rp)-rp^2 (r1g-2 \[Rho]rg))))))-(a Sqrt[Ar1g B0] ellE (-8 a dEdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dLzdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dEdp EEg (-1+EEg^2) (-a EEg+Lzg)^2-(4 (-1+EEg^2)^2 (-a EEg+Lzg)^2)/((-1+e) r1g)+4 dEdp EEg r1g (2+(-1+EEg^2) r1g)+((-1+EEg^2)^2 r1g (2+(-1+EEg^2) r1g))/(-1+e)) (B0^2+2 r1g \[Rho]rg) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/(6 B0^4 (-1+EEg^2)^3 r1g^3 \[Rho]ig (1/(4 (-1+EEg^2)^3 r1g \[Rho]ig) (-8 a dEdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dLzdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dEdp EEg (-1+EEg^2) (-a EEg+Lzg)^2-(4 (-1+EEg^2)^2 (-a EEg+Lzg)^2)/((-1+e) r1g)+4 dEdp EEg r1g (2+(-1+EEg^2) r1g)+((-1+EEg^2)^2 r1g (2+(-1+EEg^2) r1g))/(-1+e)) (p/(-1+e)^2-(p (r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (Ar1g B0+rp (r1g+rp-2 \[Rho]rg)))/(Ar1g (-1+e)^2 r1g (B0 (r1g-rp)+Ar1g rp) Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2])+(B0^2 ellE p Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (-1+e)^2 r1g)-(ellE p (r1g-2 \[Rho]rg) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (1+e)^2)+(p r1g (r1g-rp) rp Sqrt[rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2] (2 Ar1g B0 (rp-\[Rho]rg)-(r1g-2 \[Rho]rg) (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)))/(Ar1g B0 (1+e)^2 Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (B0^2 (r1g-2 rp)-rp^2 (r1g-2 \[Rho]rg))))-1/(2 r1g \[Rho]ig) (-((4 a dEde (a EEg-Lzg))/(-1+EEg^2))+(4 dLzde (a EEg-Lzg))/(-1+EEg^2)+(4 dEde EEg (-a EEg+Lzg)^2)/(-1+EEg^2)^2+(2 (-a EEg+Lzg)^2 p)/((-1+e)^2 (-1+EEg^2) r1g)+(2 dEde EEg r1g (2+(-1+EEg^2) r1g))/(-1+EEg^2)^3+(p r1g (2+(-1+EEg^2) r1g))/((-1+e)^2 (2-2 EEg^2))) (1/(1-e)+((r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (Ar1g B0+rp (r1g+rp-2 \[Rho]rg)))/(Ar1g (-1+e) r1g (B0 (r1g-rp)+Ar1g rp) Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2])+(B0^2 ellE Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (r1g-e r1g))+(ellE (r1g-2 \[Rho]rg) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (1+e))+(r1g (r1g-rp) rp Sqrt[rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2] (2 Ar1g B0 (-rp+\[Rho]rg)+(r1g-2 \[Rho]rg) (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)))/(Ar1g B0 (1+e) Sqrt[r1g (B0^2 (r1g-rp)+rp (Ar1g^2+r1g (-r1g+rp)))] (B0^2 (r1g-2 rp)-rp^2 (r1g-2 \[Rho]rg))))))-(a (-8 a dEdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dLzdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dEdp EEg (-1+EEg^2) (-a EEg+Lzg)^2-(4 (-1+EEg^2)^2 (-a EEg+Lzg)^2)/((-1+e) r1g)+4 dEdp EEg r1g (2+(-1+EEg^2) r1g)+((-1+EEg^2)^2 r1g (2+(-1+EEg^2) r1g))/(-1+e)) (r1g-rp) rp Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2] (-((4 Ar1g krg (Ar1g rp+B0 (-r1g+rp)) \[Rho]rg)/(r1g (B0 (r1g-rp)+Ar1g rp)))+(-\[Rho]ig^2+\[Rho]rg^2)/B0))/(12 B0^3 (-1+EEg^2)^3 r1g Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] \[Rho]ig (1/(4 (-1+EEg^2)^3 r1g \[Rho]ig) (-8 a dEdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dLzdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dEdp EEg (-1+EEg^2) (-a EEg+Lzg)^2-(4 (-1+EEg^2)^2 (-a EEg+Lzg)^2)/((-1+e) r1g)+4 dEdp EEg r1g (2+(-1+EEg^2) r1g)+((-1+EEg^2)^2 r1g (2+(-1+EEg^2) r1g))/(-1+e)) (p/(-1+e)^2-(p (r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (Ar1g B0+rp (r1g+rp-2 \[Rho]rg)))/(Ar1g (-1+e)^2 r1g (B0 (r1g-rp)+Ar1g rp) Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2])+(B0^2 ellE p Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (-1+e)^2 r1g)-(ellE p (r1g-2 \[Rho]rg) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (1+e)^2)+(p r1g (r1g-rp) rp Sqrt[rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2] (2 Ar1g B0 (rp-\[Rho]rg)-(r1g-2 \[Rho]rg) (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)))/(Ar1g B0 (1+e)^2 Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (B0^2 (r1g-2 rp)-rp^2 (r1g-2 \[Rho]rg))))-1/(2 r1g \[Rho]ig) (-((4 a dEde (a EEg-Lzg))/(-1+EEg^2))+(4 dLzde (a EEg-Lzg))/(-1+EEg^2)+(4 dEde EEg (-a EEg+Lzg)^2)/(-1+EEg^2)^2+(2 (-a EEg+Lzg)^2 p)/((-1+e)^2 (-1+EEg^2) r1g)+(2 dEde EEg r1g (2+(-1+EEg^2) r1g))/(-1+EEg^2)^3+(p r1g (2+(-1+EEg^2) r1g))/((-1+e)^2 (2-2 EEg^2))) (1/(1-e)+((r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (Ar1g B0+rp (r1g+rp-2 \[Rho]rg)))/(Ar1g (-1+e) r1g (B0 (r1g-rp)+Ar1g rp) Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2])+(B0^2 ellE Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (r1g-e r1g))+(ellE (r1g-2 \[Rho]rg) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (1+e))+(r1g (r1g-rp) rp Sqrt[rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2] (2 Ar1g B0 (-rp+\[Rho]rg)+(r1g-2 \[Rho]rg) (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)))/(Ar1g B0 (1+e) Sqrt[r1g (B0^2 (r1g-rp)+rp (Ar1g^2+r1g (-r1g+rp)))] (B0^2 (r1g-2 rp)-rp^2 (r1g-2 \[Rho]rg))))))+((-8 a dEdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dLzdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dEdp EEg (-1+EEg^2) (-a EEg+Lzg)^2-(4 (-1+EEg^2)^2 (-a EEg+Lzg)^2)/((-1+e) r1g)+4 dEdp EEg r1g (2+(-1+EEg^2) r1g)+((-1+EEg^2)^2 r1g (2+(-1+EEg^2) r1g))/(-1+e)) (-((6 EEg (a EEg-Lzg))/(-1+EEg^2))-(2 a^2 (-1+e^2)^2 (a^2 (-1+EEg^2)-Lzg^2))/((-3+e)^2 EEg (-1+EEg^2) (a EEg-Lzg) p^2)+(2 (-1+e^2)^2 (a^3 (-1+EEg^2)-a Lzg^2))/((-3+e)^2 (-1+EEg^2)^2 p^2)+a^3/r1g-a r1g+a (2/(-1+EEg^2)+r1g)+(a Lzg^2)/(r1g-EEg^2 r1g)+(a (a^2 (-1+EEg^2)-Lzg^2) r1g (2+(-1+EEg^2) r1g))/(2 (-1+EEg^2) (-a EEg+Lzg)^2)-((a^3 (-1+EEg^2)-a Lzg^2) r1g (2+(-1+EEg^2) r1g))/(4 (-1+EEg^2) (-a EEg+Lzg)^2)+((-1+e^2)^2 (a^3 (-1+EEg^2)-a Lzg^2) r1g (2+(-1+EEg^2) r1g))/((-3+e)^2 (-1+EEg^2)^3 (-a EEg+Lzg)^2 p^2)+(4 (a EEg-Lzg) (a^3 EEg^2+EEg Lzg (-3+r1g) r1g^2+((-1+e^2)^2 (a^3 (-1+EEg^2)-a Lzg^2) r1g^5)/(2 (-3+e)^2 (-a EEg+Lzg)^2 p^2)+a (Lzg^2+3 EEg^2 r1g^2+((-1+e^2)^2 Lzg (-a^3 (-1+EEg^2)+a Lzg^2) r1g^2)/((-3+e)^2 EEg (-a EEg+Lzg)^2 p^2))+1/2 a^2 (-4 EEg Lzg+((-1+e^2)^2 (a^3 (-1+EEg^2)-a Lzg^2) r1g^2 (2+r1g))/((-3+e)^2 (-a EEg+Lzg)^2 p^2))))/((-1+EEg^2) (2 a EEg+Lzg (-2+r1g)) r1g^2)) (p/(-1+e)^2-(p (r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (Ar1g B0+rp (r1g+rp-2 \[Rho]rg)))/(Ar1g (-1+e)^2 r1g (B0 (r1g-rp)+Ar1g rp) Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2])+(B0^2 ellE p Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (-1+e)^2 r1g)-(ellE p (r1g-2 \[Rho]rg) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (1+e)^2)+(p r1g (r1g-rp) rp Sqrt[rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2] (2 Ar1g B0 (rp-\[Rho]rg)-(r1g-2 \[Rho]rg) (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)))/(Ar1g B0 (1+e)^2 Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (B0^2 (r1g-2 rp)-rp^2 (r1g-2 \[Rho]rg)))))/(4 (-1+EEg^2)^3 r1g (-((4 a dEde (a EEg-Lzg))/(-1+EEg^2))+(4 dLzde (a EEg-Lzg))/(-1+EEg^2)+(4 dEde EEg (-a EEg+Lzg)^2)/(-1+EEg^2)^2+(2 (-a EEg+Lzg)^2 p)/((-1+e)^2 (-1+EEg^2) r1g)+(2 dEde EEg r1g (2+(-1+EEg^2) r1g))/(-1+EEg^2)^3+(p r1g (2+(-1+EEg^2) r1g))/((-1+e)^2 (2-2 EEg^2))) \[Rho]ig (1/(4 (-1+EEg^2)^3 r1g \[Rho]ig) (-8 a dEdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dLzdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dEdp EEg (-1+EEg^2) (-a EEg+Lzg)^2-(4 (-1+EEg^2)^2 (-a EEg+Lzg)^2)/((-1+e) r1g)+4 dEdp EEg r1g (2+(-1+EEg^2) r1g)+((-1+EEg^2)^2 r1g (2+(-1+EEg^2) r1g))/(-1+e)) (p/(-1+e)^2-(p (r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (Ar1g B0+rp (r1g+rp-2 \[Rho]rg)))/(Ar1g (-1+e)^2 r1g (B0 (r1g-rp)+Ar1g rp) Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2])+(B0^2 ellE p Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (-1+e)^2 r1g)-(ellE p (r1g-2 \[Rho]rg) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (1+e)^2)+(p r1g (r1g-rp) rp Sqrt[rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2] (2 Ar1g B0 (rp-\[Rho]rg)-(r1g-2 \[Rho]rg) (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)))/(Ar1g B0 (1+e)^2 Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (B0^2 (r1g-2 rp)-rp^2 (r1g-2 \[Rho]rg))))-1/(2 r1g \[Rho]ig) (-((4 a dEde (a EEg-Lzg))/(-1+EEg^2))+(4 dLzde (a EEg-Lzg))/(-1+EEg^2)+(4 dEde EEg (-a EEg+Lzg)^2)/(-1+EEg^2)^2+(2 (-a EEg+Lzg)^2 p)/((-1+e)^2 (-1+EEg^2) r1g)+(2 dEde EEg r1g (2+(-1+EEg^2) r1g))/(-1+EEg^2)^3+(p r1g (2+(-1+EEg^2) r1g))/((-1+e)^2 (2-2 EEg^2))) (1/(1-e)+((r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (Ar1g B0+rp (r1g+rp-2 \[Rho]rg)))/(Ar1g (-1+e) r1g (B0 (r1g-rp)+Ar1g rp) Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2])+(B0^2 ellE Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (r1g-e r1g))+(ellE (r1g-2 \[Rho]rg) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (1+e))+(r1g (r1g-rp) rp Sqrt[rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2] (2 Ar1g B0 (-rp+\[Rho]rg)+(r1g-2 \[Rho]rg) (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)))/(Ar1g B0 (1+e) Sqrt[r1g (B0^2 (r1g-rp)+rp (Ar1g^2+r1g (-r1g+rp)))] (B0^2 (r1g-2 rp)-rp^2 (r1g-2 \[Rho]rg))))))+\[Delta]r41 (-(((-8 a dEdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dLzdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dEdp EEg (-1+EEg^2) (-a EEg+Lzg)^2-(4 (-1+EEg^2)^2 (-a EEg+Lzg)^2)/((-1+e) r1g)+4 dEdp EEg r1g (2+(-1+EEg^2) r1g)+((-1+EEg^2)^2 r1g (2+(-1+EEg^2) r1g))/(-1+e)) (r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2])/(2 B0 (-1+EEg^2)^3 r1g^2 (B0 (r1g-rp)+Ar1g rp) \[Rho]ig (1/(4 (-1+EEg^2)^3 r1g \[Rho]ig) (-8 a dEdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dLzdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dEdp EEg (-1+EEg^2) (-a EEg+Lzg)^2-(4 (-1+EEg^2)^2 (-a EEg+Lzg)^2)/((-1+e) r1g)+4 dEdp EEg r1g (2+(-1+EEg^2) r1g)+((-1+EEg^2)^2 r1g (2+(-1+EEg^2) r1g))/(-1+e)) (p/(-1+e)^2-(p (r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (Ar1g B0+rp (r1g+rp-2 \[Rho]rg)))/(Ar1g (-1+e)^2 r1g (B0 (r1g-rp)+Ar1g rp) Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2])+(B0^2 ellE p Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (-1+e)^2 r1g)-(ellE p (r1g-2 \[Rho]rg) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (1+e)^2)+(p r1g (r1g-rp) rp Sqrt[rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2] (2 Ar1g B0 (rp-\[Rho]rg)-(r1g-2 \[Rho]rg) (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)))/(Ar1g B0 (1+e)^2 Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (B0^2 (r1g-2 rp)-rp^2 (r1g-2 \[Rho]rg))))-1/(2 r1g \[Rho]ig) (-((4 a dEde (a EEg-Lzg))/(-1+EEg^2))+(4 dLzde (a EEg-Lzg))/(-1+EEg^2)+(4 dEde EEg (-a EEg+Lzg)^2)/(-1+EEg^2)^2+(2 (-a EEg+Lzg)^2 p)/((-1+e)^2 (-1+EEg^2) r1g)+(2 dEde EEg r1g (2+(-1+EEg^2) r1g))/(-1+EEg^2)^3+(p r1g (2+(-1+EEg^2) r1g))/((-1+e)^2 (2-2 EEg^2))) (1/(1-e)+((r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (Ar1g B0+rp (r1g+rp-2 \[Rho]rg)))/(Ar1g (-1+e) r1g (B0 (r1g-rp)+Ar1g rp) Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2])+(B0^2 ellE Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (r1g-e r1g))+(ellE (r1g-2 \[Rho]rg) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (1+e))+(r1g (r1g-rp) rp Sqrt[rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2] (2 Ar1g B0 (-rp+\[Rho]rg)+(r1g-2 \[Rho]rg) (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)))/(Ar1g B0 (1+e) Sqrt[r1g (B0^2 (r1g-rp)+rp (Ar1g^2+r1g (-r1g+rp)))] (B0^2 (r1g-2 rp)-rp^2 (r1g-2 \[Rho]rg)))))))+(Ar1g^2 ellE (-8 a dEdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dLzdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dEdp EEg (-1+EEg^2) (-a EEg+Lzg)^2-(4 (-1+EEg^2)^2 (-a EEg+Lzg)^2)/((-1+e) r1g)+4 dEdp EEg r1g (2+(-1+EEg^2) r1g)+((-1+EEg^2)^2 r1g (2+(-1+EEg^2) r1g))/(-1+e)) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/(2 (Ar1g B0)^(3/2) (-1+EEg^2)^3 r1g^2 \[Rho]ig (1/(4 (-1+EEg^2)^3 r1g \[Rho]ig) (-8 a dEdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dLzdp (-1+EEg^2)^2 (a EEg-Lzg)+8 dEdp EEg (-1+EEg^2) (-a EEg+Lzg)^2-(4 (-1+EEg^2)^2 (-a EEg+Lzg)^2)/((-1+e) r1g)+4 dEdp EEg r1g (2+(-1+EEg^2) r1g)+((-1+EEg^2)^2 r1g (2+(-1+EEg^2) r1g))/(-1+e)) (p/(-1+e)^2-(p (r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (Ar1g B0+rp (r1g+rp-2 \[Rho]rg)))/(Ar1g (-1+e)^2 r1g (B0 (r1g-rp)+Ar1g rp) Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2])+(B0^2 ellE p Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (-1+e)^2 r1g)-(ellE p (r1g-2 \[Rho]rg) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (1+e)^2)+(p r1g (r1g-rp) rp Sqrt[rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2] (2 Ar1g B0 (rp-\[Rho]rg)-(r1g-2 \[Rho]rg) (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)))/(Ar1g B0 (1+e)^2 Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (B0^2 (r1g-2 rp)-rp^2 (r1g-2 \[Rho]rg))))-1/(2 r1g \[Rho]ig) (-((4 a dEde (a EEg-Lzg))/(-1+EEg^2))+(4 dLzde (a EEg-Lzg))/(-1+EEg^2)+(4 dEde EEg (-a EEg+Lzg)^2)/(-1+EEg^2)^2+(2 (-a EEg+Lzg)^2 p)/((-1+e)^2 (-1+EEg^2) r1g)+(2 dEde EEg r1g (2+(-1+EEg^2) r1g))/(-1+EEg^2)^3+(p r1g (2+(-1+EEg^2) r1g))/((-1+e)^2 (2-2 EEg^2))) (1/(1-e)+((r1g-rp) Sqrt[r1g (Ar1g^2 rp+(-r1g+rp) (-B0^2+r1g rp))] (Ar1g B0+rp (r1g+rp-2 \[Rho]rg)))/(Ar1g (-1+e) r1g (B0 (r1g-rp)+Ar1g rp) Sqrt[\[Rho]ig^2+(rp-\[Rho]rg)^2])+(B0^2 ellE Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (r1g-e r1g))+(ellE (r1g-2 \[Rho]rg) Sqrt[(r1g-rp) rp (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)])/((Ar1g B0)^(3/2) (1+e))+(r1g (r1g-rp) rp Sqrt[rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2] (2 Ar1g B0 (-rp+\[Rho]rg)+(r1g-2 \[Rho]rg) (rp^2+\[Rho]ig^2-2 rp \[Rho]rg+\[Rho]rg^2)))/(Ar1g B0 (1+e) Sqrt[r1g (B0^2 (r1g-rp)+rp (Ar1g^2+r1g (-r1g+rp)))] (B0^2 (r1g-2 rp)-rp^2 (r1g-2 \[Rho]rg)))))))
]


(* ::Subsection::Closed:: *)
(*Function to compute shifts orbital elements based on the spin gauge*)


\[Delta]orbelements[a_,p_,e_,x_,spingauge_?StringQ]:=Module[{\[Delta]r1,\[Delta]r2},
	Which[
		spingauge=="FC",
		{\[Delta]r1,\[Delta]r2}={\[Delta]r1funFC[a,p,e,x],\[Delta]r2funFC[a,p,e,x]};
		{1/2((1-e)^2*\[Delta]r1+(1+e)^2*\[Delta]r2),((1-e^2)((1-e)\[Delta]r1-\[Delta]r2-e*\[Delta]r2))/(2p)}
		,
		spingauge=="FF",
		If[e==0,
			{1/p(a-Sqrt[p]),0}
			,
			{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];
			{1/2((1-e)^2*\[Delta]r1+(1+e)^2*\[Delta]r2),((1-e^2)((1-e)\[Delta]r1-\[Delta]r2-e*\[Delta]r2))/(2p)}
		]
		,
		spingauge=="FE",
		{\[Delta]pfunFE[a,p,e,x],0}
		,
		spingauge=="RPer",
		{\[Delta]pfunRPer[a,p,e,x],\[Delta]efunRPer[a,p,e,x]}
		,
		spingauge=="RHom",
		{\[Delta]pfunRHom[a,p,e,x],\[Delta]efunRHom[a,p,e,x]}
		,
		spingauge=="RPlunge",
		{\[Delta]pfunRPlunge[a,p,e,x],\[Delta]efunRPlunge[a,p,e,x]}
	]
]


(* ::Subsubsection::Closed:: *)
(*Shifts constants of motion and IBCO radius for parabolic orbits*)


KerrNearEqSpinOrbitCorr\[Delta]IBCO[a_, x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,r,\[CapitalDelta],\[Delta]rIBCO,Lzred},
	(* geodesic constants of motion *)
	EEg=1;
	Lzg=RealSign[x]2(1+Sqrt[1-RealSign[x]a]);
	r=2-RealSign[x]a+2Sqrt[1-RealSign[x]a];
	\[CapitalDelta]=r^2-2r+a^2;
	Lzred=Lzg-a;
	(* spin correction constants of motion *)
	\[Delta]EE=0;
	\[Delta]Lz=(1-RealSign[x]a+Sqrt[1-RealSign[x]a])/(2+2Sqrt[1-a]-RealSign[x]a);
	\[Delta]rIBCO=-(1/2)((Lzg (r-2)^3 Sqrt[\[CapitalDelta]](2a*Lzred^4-2Lzred^3 r^2-2Lzred*a*Lzg*r^2+a*Lzg^2*r^3))/(Lzred^2(2a+Lzg(-2+r))r(-r Sqrt[\[CapitalDelta]](10a^2-4(4+\[CapitalDelta])+r(8+\[CapitalDelta]))+Sqrt[2]a Sqrt[r](2a^2(-1+r)+(4+r)\[CapitalDelta])RealSign[x])));
			
	<|
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->RealSign[x]*\[Delta]EE,
	 "Jzs"->RealSign[x]*\[Delta]Lz,
	 "rIBCO"->r,
	 "\[Delta]rIBCO"->RealSign[x]*\[Delta]rIBCO
	|>
]


(* ::Section::Closed:: *)
(*Shift frequencies*)


(* ::Subsection::Closed:: *)
(*Shift radial frequency*)


(* ::Subsubsection::Closed:: *)
(*Fixed turning points*)


\[Delta]\[CapitalUpsilon]rfunFT[a_,p_,e_,x_]:=Module[{\[CapitalUpsilon]rg,EEg,Lzg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,krg,ellK,ellE},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	\[Delta]r3=\[Delta]r3funFT[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));

	ellK=EllipticK[krg];
	ellE=EllipticE[krg];

	-2\[CapitalUpsilon]rg^2/(2\[Pi])/Sqrt[1-EEg^2]2/Sqrt[(r1g-r3g)r2g]((EEg*\[Delta]EE)/(1-EEg^2)ellK+\[Delta]r3/2(1/r3g(r2g/(r2g-r3g)ellE-ellK))+\[Delta]r41((r1g*ellK-(r1g-r3g)ellE)/(r1g*r3g ))-a/2(-(1/(3r1g*r2g*r3g))((r1g+r2g+r3g)ellK-2(r2g*r3g+r1g(r2g+r3g))((-((r1g-r3g)ellE)+r1g*ellK)/(r1g*r3g)))))
]


(* ::Subsubsection::Closed:: *)
(*Fixed constants of motion*)


\[Delta]\[CapitalUpsilon]rfunFC[a_,p_,e_,x_]:=Module[{\[CapitalUpsilon]rg,EEg,Lzg,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,krg,ellK,ellE},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	EEg=EEgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	\[Delta]r1=\[Delta]r1funFC[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFC[a,p,e,x];
	\[Delta]r3=\[Delta]r3funFC[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));

	ellK=EllipticK[krg];
	ellE=EllipticE[krg];

	-2 \[CapitalUpsilon]rg^2/(2\[Pi]) 1/Sqrt[1-EEg^2] 2/Sqrt[(r1g-r3g)r2g]((r2g*ellE-r1g*ellK)/(2r1g(r1g-r2g))\[Delta]r1+1/2(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellK/(r1g-r2g))\[Delta]r2+1/2 \[Delta]r3(1/r3g (r2g/(r2g-r3g) ellE-ellK))+\[Delta]r41((r1g*ellK-(r1g-r3g)ellE)/(r1g*r3g ))-1/2 a(-(1/(3r1g*r2g*r3g))((r1g+r2g+r3g)ellK-2(r2g*r3g+r1g(r2g+r3g))( (-((r1g-r3g)ellE)+r1g*ellK)/(r1g*r3g)))))
]


(* ::Subsubsection::Closed:: *)
(*Fixed frequency*)


\[Delta]\[CapitalUpsilon]rfunFF[a_,p_,e_,x_,{\[Delta]r1_,\[Delta]r2_}]:=Module[{\[CapitalUpsilon]rg,EEg,Lzg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,Yint,krg,ellK,ellE},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	EEg=EEgfun[a,p,e,x];
	
	(*{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];*)
	\[Delta]r3=\[Delta]r3funFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];	
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	
	Yint=Sqrt[(r1g-r3g)r2g];
	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));

	ellK=EllipticK[krg];
	ellE=EllipticE[krg];

	-2\[CapitalUpsilon]rg^2/(2\[Pi])1/Sqrt[1-EEg^2] 2/Yint((EEg*\[Delta]EE)/(1-EEg^2)ellK+(r2g*ellE-r1g*ellK)/(2r1g(r1g-r2g))\[Delta]r1+1/2(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellK/(r1g-r2g))\[Delta]r2+1/2*\[Delta]r3(1/r3g(r2g/(r2g-r3g)ellE-ellK))+\[Delta]r41((r1g*ellK-(r1g-r3g)ellE)/(r1g*r3g))-1/2*a(-(1/(3r1g*r2g*r3g))((r1g+r2g+r3g)ellK-2(r2g*r3g+r1g(r2g+r3g))((-((r1g-r3g)ellE)+r1g*ellK)/(r1g*r3g)))))
]


(* ::Subsubsection::Closed:: *)
(*Fixed eccentricity*)


\[Delta]\[CapitalUpsilon]rfunFE[a_,p_,e_,x_]:=Module[{\[CapitalUpsilon]rg,EEg,Lzg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,krg,ellK,ellE},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));

	ellK=EllipticK[krg];
	ellE=EllipticE[krg];

	-2\[CapitalUpsilon]rg^2/(2\[Pi])1/Sqrt[1-EEg^2] 2/Sqrt[(r1g-r3g)r2g]((EEg \[Delta]EE)/(1-EEg^2)ellK+(r2g*ellE-r1g*ellK)/(2r1g(r1g-r2g))\[Delta]r1+1/2(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellK/(r1g-r2g))\[Delta]r2+1/2 \[Delta]r3(1/r3g(r2g/(r2g-r3g)ellE-ellK))+\[Delta]r41((r1g ellK-(r1g-r3g)ellE)/(r1g*r3g))-1/2*a(-(1/(3r1g*r2g*r3g))((r1g+r2g+r3g)ellK-2(r2g*r3g+r1g(r2g+r3g))((-((r1g-r3g)ellE)+r1g*ellK)/(r1g*r3g)))))
]


\[Delta]\[CapitalUpsilon]rover\[CapitalUpsilon]rgfunLimFE[a_,p_,e_,x_]:=Module[{\[CapitalUpsilon]rg,EEg,Lzg,\[Delta]EE,r1g,r2g,\[Delta]r1,\[Delta]r2,\[Delta]r41},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	If[e==0,
		1/2(a/r1g^2-(2EEg*\[Delta]EE)/(1-EEg^2)+2\[Delta]r1/r1g-(2\[Delta]r41)/r1g),
		1/2(a/r2g^2-(2EEg*\[Delta]EE)/(1-EEg^2)+\[Delta]r1/(r1g-r2g)-\[Delta]r2/(r1g-r2g)+\[Delta]r2/r2g-(2\[Delta]r41)/r2g)
	]
]


(* ::Subsubsection::Closed:: *)
(*R gauge*)


\[Delta]\[CapitalUpsilon]rfunR[a_,p_,e_,x_]:=Module[{\[CapitalUpsilon]rg,EEg,Lzg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,krg,ellK,ellE},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	\[Delta]r1=\[Delta]r1funRPer[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRPer[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));

	ellK=EllipticK[krg];
	ellE=EllipticE[krg];

	-2\[CapitalUpsilon]rg^2/(2\[Pi])1/Sqrt[1-EEg^2] 2/Sqrt[(r1g-r3g)r2g]((EEg \[Delta]EE)/(1-EEg^2)ellK+(r2g*ellE-r1g*ellK)/(2r1g(r1g-r2g))\[Delta]r1+1/2(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellK/(r1g-r2g))\[Delta]r2+1/2 \[Delta]r3(1/r3g(r2g/(r2g-r3g)ellE-ellK))+\[Delta]r41((r1g ellK-(r1g-r3g)ellE)/(r1g*r3g))-1/2*a(-(1/(3r1g*r2g*r3g))((r1g+r2g+r3g)ellK-2(r2g*r3g+r1g(r2g+r3g))((-((r1g-r3g)ellE)+r1g*ellK)/(r1g*r3g)))))
]


\[Delta]\[CapitalUpsilon]rover\[CapitalUpsilon]rgfunLimR[a_,p_,e_,x_]:=Module[{\[CapitalUpsilon]rg,EEg,Lzg,\[Delta]EE,r1g,r2g,\[Delta]r1,\[Delta]r2,\[Delta]r41},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	
	\[Delta]r1=\[Delta]r1funRHom[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRHom[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	If[e==0,
		1/2(a/r1g^2-(2EEg*\[Delta]EE)/(1-EEg^2)+2\[Delta]r1/r1g-(2\[Delta]r41)/r1g),
		1/2(a/r2g^2-(2EEg*\[Delta]EE)/(1-EEg^2)+\[Delta]r1/(r1g-r2g)-\[Delta]r2/(r1g-r2g)+\[Delta]r2/r2g-(2\[Delta]r41)/r2g)
	]
]


(* ::Subsection::Closed:: *)
(*Shift "coordinate time frequency"*)


(* ::Subsubsection::Closed:: *)
(*Fixed turning points*)


\[Delta]\[CapitalUpsilon]tfunFT[a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[CapitalUpsilon]rg,\[CapitalUpsilon]tg,r1g,r2g,r3g,\[Delta]r3,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,krg,\[Gamma]r,ellK,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r3=\[Delta]r3funFT[a,p,e,x];
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];
	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));
	\[Gamma]r=(r1g-r2g)/(r1g-r3g);

	ellK=EllipticK[krg];
	ellE=EllipticE[krg];
	ellPi=EllipticPi[\[Gamma]r,krg];

	\[ScriptCapitalI]=(2ellK)/(Sqrt[1-EEg^2]Sqrt[r2g(r1g-r3g)]);
	\[ScriptCapitalI]r3g=2/(Sqrt[1-EEg^2]Sqrt[r2g(r1g-r3g)]*r3g)((r2g*ellE)/(r2g-r3g)-ellK);

	\[ScriptCapitalI]rover=(2(r3g*ellK+(r2g-r3g)ellPi))/(Sqrt[1-EEg^2]*Sqrt[r2g(r1g-r3g)]);
	\[ScriptCapitalI]r2over=(2r2g^2)/(Sqrt[(1-EEg^2)]Sqrt[(r1g-r3g)r2g])((r3g/r2g)^2*ellK+2r3g/r2g(1-r3g/r2g)ellPi+(1-r3g/r2g)^2 1/(2(\[Gamma]r-1)(krg-\[Gamma]r))(\[Gamma]r*ellE+(krg-\[Gamma]r)ellK+(2\[Gamma]r*krg+2\[Gamma]r-\[Gamma]r^2-3krg)ellPi));
	(2\[CapitalUpsilon]rg)/(2\[Pi])((4\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+(r3g*\[Delta]r3)/(2(r3g-rm)(r3g-rp))(-2a*Lzg+EEg*r3g^3+a^2*EEg(2+r3g))\[ScriptCapitalI]r3g+EEg/2((2+r3g)\[Delta]r3-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+(EEg*\[Delta]r3)/2+EEg*\[Delta]\[Rho]r4)\[ScriptCapitalI]rover)+\[CapitalUpsilon]tg/\[CapitalUpsilon]rg*\[Delta]\[CapitalUpsilon]rfunFT[a,p,e,x]
]


(* ::Subsubsection::Closed:: *)
(*Fixed constants of motion*)


\[Delta]\[CapitalUpsilon]tfunFC[a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[CapitalUpsilon]rg,\[CapitalUpsilon]tg,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,krg,\[Gamma]r,ellK,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFC[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFC[a,p,e,x];
	\[Delta]r3=\[Delta]r3funFC[a,p,e,x];
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];

	krg=((r1g-r2g) r3g)/(r2g(r1g-r3g));
	\[Gamma]r=1-r1g/r2g;

	ellK=EllipticK[krg];
	ellE=EllipticE[krg];
	ellPi=EllipticPi[\[Gamma]r,krg];
	
	\[ScriptCapitalI]=(2 ellK)/(Sqrt[1-EEg^2] Sqrt[r2g (r1g-r3g)]);
	\[ScriptCapitalI]r1g=2 /(Sqrt[1-EEg^2] Sqrt[r2g (r1g-r3g)] r1g) ((r2g*ellE-r1g*ellK)/(r1g-r2g) );
	\[ScriptCapitalI]r2g=2 /(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]) (-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellK/(r1g-r2g));
	\[ScriptCapitalI]r3g=2 /(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]r3g) ((r2g*ellE)/(r2g-r3g)-ellK);
	\[ScriptCapitalI]rover=(2r1g*ellPi)/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]);

	(2\[CapitalUpsilon]rg)/(2\[Pi]) (1/2 EEg((2+r1g)\[Delta]r1+(2+r2g)\[Delta]r2+(2+r3g)\[Delta]r3-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+(1/2 EEg(\[Delta]r1+\[Delta]r2+\[Delta]r3+2\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+(r1g(-2a*Lzg+EEg*r1g^3+a^2*EEg(2+r1g)))/(2(r1g-rm)(r1g-rp)) \[ScriptCapitalI]r1g*\[Delta]r1+(r2g(-2a*Lzg+EEg*r2g^3+a^2*EEg(2+r2g)))/(2(r2g-rm)(r2g-rp)) \[ScriptCapitalI]r2g*\[Delta]r2+(r3g(-2a*Lzg+EEg*r3g^3+a^2*EEg(2+r3g)))/(2(r3g-rm)(r3g-rp)) \[ScriptCapitalI]r3g*\[Delta]r3)+\[CapitalUpsilon]tg/\[CapitalUpsilon]rg \[Delta]\[CapitalUpsilon]rfunFC[a,p,e,x]
]


(* ::Subsubsection::Closed:: *)
(*Fixed frequency*)


\[Delta]\[CapitalUpsilon]tfunFF[a_,p_,e_,x_,{\[Delta]r1_,\[Delta]r2_}]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[CapitalUpsilon]rg,\[CapitalUpsilon]tg,r1g,r2g,r3g,\[Delta]r3,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,krg,\[Gamma]r,ellK,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	(*{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];*)
	\[Delta]r3=\[Delta]r3funFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}];
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;
	
	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));
	\[Gamma]r=1-r1g/r2g;

	ellK=EllipticK[krg];
	ellE=EllipticE[krg];
	ellPi=EllipticPi[\[Gamma]r,krg];

	\[ScriptCapitalI]=(2ellK)/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]);
	\[ScriptCapitalI]r1g=2/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]r1g) ((r2g *ellE-r1g*ellK)/(r1g-r2g) );
	\[ScriptCapitalI]r2g=2/(Sqrt[1-EEg^2] Sqrt[r2g (r1g-r3g)]) (-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellK/(r1g-r2g));
	\[ScriptCapitalI]r3g=2/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]r3g) ((r2g*ellE)/(r2g-r3g)-ellK);
	\[ScriptCapitalI]rover=(2r1g*ellPi)/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]);
	\[ScriptCapitalI]r2over=(r2g((r1g-r3g)ellE-r1g*ellK))/(Sqrt[1-EEg^2]*Sqrt[r2g(r1g-r3g)])+(r1g(r1g+r2g+r3g)ellPi)/(Sqrt[1-EEg^2]*Sqrt[r2g(r1g-r3g)]);

	(2\[CapitalUpsilon]rg)/(2\[Pi])((4\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+1/2 EEg(\[Delta]r1+\[Delta]r2+\[Delta]r3+2\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+1/2 EEg((2+r1g)\[Delta]r1+(2+r2g)\[Delta]r2+(2+r3g)\[Delta]r3-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+(r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(2(r1g-rm)(r1g-rp)) \[ScriptCapitalI]r1g*\[Delta]r1+(r2g(-2a*Lzg+EEg(2a^2+a^2*r2g+r2g^3)))/(2(r2g-rm)(r2g-rp))\[ScriptCapitalI]r2g*\[Delta]r2+(r3g(-2a*Lzg+EEg(2a^2+a^2*r3g+r3g^3)))/(2(r3g-rm)(r3g-rp))\[ScriptCapitalI]r3g*\[Delta]r3)+\[CapitalUpsilon]tg/\[CapitalUpsilon]rg*\[Delta]\[CapitalUpsilon]rfunFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}]
]


(* ::Subsubsection::Closed:: *)
(*Fixed eccentricity*)


\[Delta]\[CapitalUpsilon]tfunFE[a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[CapitalUpsilon]rg,\[CapitalUpsilon]tg,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,krg,\[Gamma]r,ellK,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];
	
	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));
	\[Gamma]r=1-r1g/r2g;

	ellK=EllipticK[krg];
	ellE=EllipticE[krg];
	ellPi=EllipticPi[\[Gamma]r,krg];

	\[ScriptCapitalI]=(2ellK)/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]);
	\[ScriptCapitalI]r1g=2 /(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]r1g) ((r2g *ellE-r1g*ellK)/(r1g-r2g) );
	\[ScriptCapitalI]r2g=2 /(Sqrt[1-EEg^2] Sqrt[r2g (r1g-r3g)]) (-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellK/(r1g-r2g));
	\[ScriptCapitalI]r3g=2 /(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]r3g) ((r2g*ellE)/(r2g-r3g)-ellK);
	\[ScriptCapitalI]rover=(2r1g*ellPi)/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]);
	\[ScriptCapitalI]r2over=(r2g((r1g-r3g)ellE-r1g*ellK))/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)])+(r1g(r1g+r2g+r3g)ellPi)/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]);

	(2\[CapitalUpsilon]rg)/(2\[Pi])((4\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+1/2 EEg(\[Delta]r1+\[Delta]r2+\[Delta]r3+2\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+1/2 EEg((2+r1g)\[Delta]r1+(2+r2g)\[Delta]r2+(2+r3g)\[Delta]r3-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+(r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(2(r1g-rm)(r1g-rp)) \[ScriptCapitalI]r1g*\[Delta]r1+(r2g(-2a*Lzg+EEg(2a^2+a^2*r2g+r2g^3)))/(2(r2g-rm)(r2g-rp))\[ScriptCapitalI]r2g*\[Delta]r2+(r3g(-2a*Lzg+EEg(2a^2+a^2*r3g+r3g^3)))/(2(r3g-rm)(r3g-rp))\[ScriptCapitalI]r3g*\[Delta]r3)+\[CapitalUpsilon]tg/\[CapitalUpsilon]rg*\[Delta]\[CapitalUpsilon]rfunFE[a,p,e,x]
]


\[Delta]\[CapitalUpsilon]tfunLimFE[a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[CapitalUpsilon]rg,\[CapitalUpsilon]tg,r1g,r2g,\[Delta]r1,\[Delta]r2,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];
	
	If[e==0,
		((4\[Delta]EE)/(1-EEg^2)+1/2 EEg(3(2+r1g)\[Delta]r1-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)+((2\[Delta]EE)/(1-EEg^2)+1/2 EEg(3\[Delta]r1+2\[Delta]\[Rho]r4))r2g+ \[Delta]EE/(1-EEg^2) r1g^2-((4a^3*EEg-4a*EEg*r1g+Lzg (-2+r1g)^2*r1g+a^2*Lzg(-4+3r1g))\[Delta]r1)/(2(a^2+(-2+r1g)r1g)^2))+\[CapitalUpsilon]tgfunLim[a,p,e,x] \[Delta]\[CapitalUpsilon]rover\[CapitalUpsilon]rgfunLimFE[a,p,e,x],
		((4\[Delta]EE)/(1-EEg^2)+1/2 EEg((2+r1g)\[Delta]r1+(2+r2g)\[Delta]r2+(2+r2g)\[Delta]r2-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)+((2\[Delta]EE)/(1-EEg^2)+1/2 EEg(\[Delta]r1+\[Delta]r2+\[Delta]r2+2\[Delta]\[Rho]r4))r2g+ \[Delta]EE/(1-EEg^2) 1/2 (-r1g*r2g+(r1g+r2g+r2g)r2g)+(r1g(-2a*Lzg+EEg*r1g^3+a^2*EEg(2+r1g)))/(2(r1g-rm)(r1g-rp)) (-(1/(r1g-r2g)))\[Delta]r1+(r2g(-2a*Lzg+EEg*r2g^3+a^2*EEg(2+r2g)))/(2(r2g-rm)(r2g-rp)) 1/(r1g-r2g) \[Delta]r2+(r2g(-2a*Lzg+EEg*r2g^3+a^2*EEg(2+r2g)))/(2(r2g-rm)(r2g-rp)) (-(1/r2g))\[Delta]r2)+\[CapitalUpsilon]tgfunLim[a,p,e,x] \[Delta]\[CapitalUpsilon]rover\[CapitalUpsilon]rgfunLimFE[a,p,e,x]
	]
]


(* ::Subsubsection::Closed:: *)
(*R gauge*)


\[Delta]\[CapitalUpsilon]tfunR[a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[CapitalUpsilon]rg,\[CapitalUpsilon]tg,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,krg,\[Gamma]r,ellK,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funRPer[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRPer[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];
	
	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));
	\[Gamma]r=1-r1g/r2g;

	ellK=EllipticK[krg];
	ellE=EllipticE[krg];
	ellPi=EllipticPi[\[Gamma]r,krg];

	\[ScriptCapitalI]=(2ellK)/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]);
	\[ScriptCapitalI]r1g=2/(Sqrt[1-EEg^2]*Sqrt[r2g(r1g-r3g)]r1g)((r2g *ellE-r1g*ellK)/(r1g-r2g));
	\[ScriptCapitalI]r2g=2/(Sqrt[1-EEg^2]*Sqrt[r2g (r1g-r3g)])(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellK/(r1g-r2g));
	\[ScriptCapitalI]r3g=2/(Sqrt[1-EEg^2]*Sqrt[r2g(r1g-r3g)]r3g)((r2g*ellE)/(r2g-r3g)-ellK);
	\[ScriptCapitalI]rover=(2r1g*ellPi)/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]);
	\[ScriptCapitalI]r2over=(r2g((r1g-r3g)ellE-r1g*ellK))/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)])+(r1g(r1g+r2g+r3g)ellPi)/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]);

	(2\[CapitalUpsilon]rg)/(2\[Pi])((4\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+1/2 EEg(\[Delta]r1+\[Delta]r2+\[Delta]r3+2\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+1/2 EEg((2+r1g)\[Delta]r1+(2+r2g)\[Delta]r2+(2+r3g)\[Delta]r3-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+(r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(2(r1g-rm)(r1g-rp)) \[ScriptCapitalI]r1g*\[Delta]r1+(r2g(-2a*Lzg+EEg(2a^2+a^2*r2g+r2g^3)))/(2(r2g-rm)(r2g-rp))\[ScriptCapitalI]r2g*\[Delta]r2+(r3g(-2a*Lzg+EEg(2a^2+a^2*r3g+r3g^3)))/(2(r3g-rm)(r3g-rp))\[ScriptCapitalI]r3g*\[Delta]r3)+\[CapitalUpsilon]tg/\[CapitalUpsilon]rg*\[Delta]\[CapitalUpsilon]rfunR[a,p,e,x]
]


\[Delta]\[CapitalUpsilon]tfunLimR[a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[CapitalUpsilon]rg,\[CapitalUpsilon]tg,r1g,r2g,\[Delta]r1,\[Delta]r2,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	\[Delta]r1=\[Delta]r1funRHom[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRHom[a,p,e,x];
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];
	
	If[e==0,
		((4\[Delta]EE)/(1-EEg^2)+1/2 EEg(3(2+r1g)\[Delta]r1-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)+((2\[Delta]EE)/(1-EEg^2)+1/2 EEg(3\[Delta]r1+2\[Delta]\[Rho]r4))r2g+ \[Delta]EE/(1-EEg^2) r1g^2-((4a^3*EEg-4a*EEg*r1g+Lzg (-2+r1g)^2*r1g+a^2*Lzg(-4+3r1g))\[Delta]r1)/(2(a^2+(-2+r1g)r1g)^2))+\[CapitalUpsilon]tgfunLim[a,p,e,x] \[Delta]\[CapitalUpsilon]rover\[CapitalUpsilon]rgfunLimR[a,p,e,x],
		((4\[Delta]EE)/(1-EEg^2)+1/2 EEg((2+r1g)\[Delta]r1+(2+r2g)\[Delta]r2+(2+r2g)\[Delta]r2-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)+((2\[Delta]EE)/(1-EEg^2)+1/2 EEg(\[Delta]r1+\[Delta]r2+\[Delta]r2+2\[Delta]\[Rho]r4))r2g+ \[Delta]EE/(1-EEg^2) 1/2 (-r1g*r2g+(r1g+r2g+r2g)r2g)+(r1g(-2a*Lzg+EEg*r1g^3+a^2*EEg(2+r1g)))/(2(r1g-rm)(r1g-rp)) (-(1/(r1g-r2g)))\[Delta]r1+(r2g(-2a*Lzg+EEg*r2g^3+a^2*EEg(2+r2g)))/(2(r2g-rm)(r2g-rp)) 1/(r1g-r2g) \[Delta]r2+(r2g(-2a*Lzg+EEg*r2g^3+a^2*EEg(2+r2g)))/(2(r2g-rm)(r2g-rp)) (-(1/r2g))\[Delta]r2)+\[CapitalUpsilon]tgfunLim[a,p,e,x] \[Delta]\[CapitalUpsilon]rover\[CapitalUpsilon]rgfunLimR[a,p,e,x]
	]
]


(* ::Subsection::Closed:: *)
(*Shift azimuthal frequency*)


(* ::Subsubsection::Closed:: *)
(*Fixed turning points*)


\[Delta]\[CapitalUpsilon]\[Phi]funFT[a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g,r1g,r2g,r3g,\[Delta]r3,krg,ellK,ellE,\[ScriptCapitalI],\[ScriptCapitalI]r3g},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r3=\[Delta]r3funFT[a,p,e,x];

	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));

	ellK=EllipticK[krg];
	ellE=EllipticE[krg];

	\[ScriptCapitalI]=(2ellK)/(Sqrt[1-EEg^2]Sqrt[r2g(r1g-r3g)]);
	\[ScriptCapitalI]r3g=2/(Sqrt[1-EEg^2]Sqrt[r2g(r1g-r3g)]r3g)((r2g*ellE)/(r2g-r3g)-ellK);

	(2\[CapitalUpsilon]rg)/(2\[Pi])(-EEg(1-(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+((2a*EEg+Lzg(-2+r3g))r3g*\[Delta]r3)/(2(r3g-rm)(r3g-rp))\[ScriptCapitalI]r3g)+\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]rg*\[Delta]\[CapitalUpsilon]rfunFT[a,p,e,x]
]


(* ::Subsubsection::Closed:: *)
(*Fixed constants of motion*)


\[Delta]\[CapitalUpsilon]\[Phi]funFC[a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,krg,ellK,ellE,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFC[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFC[a,p,e,x];
	\[Delta]r3=\[Delta]r3funFC[a,p,e,x];

	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));

	ellK=EllipticK[krg];
	ellE=EllipticE[krg];

	\[ScriptCapitalI]=(2ellK)/(Sqrt[1-EEg^2] Sqrt[r2g (r1g-r3g)]);
	\[ScriptCapitalI]r1g=2/(Sqrt[1-EEg^2] Sqrt[r2g (r1g-r3g)] r1g)((r2g*ellE-r1g*ellK)/(r1g-r2g));
	\[ScriptCapitalI]r2g=2/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)])(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellK/(r1g-r2g));
	\[ScriptCapitalI]r3g=2/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]r3g)((r2g*ellE)/(r2g-r3g)-ellK);

	(2\[CapitalUpsilon]rg)/(2\[Pi])(-EEg*\[ScriptCapitalI]+((2a*EEg+Lzg(-2+r1g))r1g*\[ScriptCapitalI]r1g*\[Delta]r1)/(2(r1g-rm)(r1g-rp))+((2a*EEg+Lzg(-2+r2g))r2g*\[ScriptCapitalI]r2g*\[Delta]r2)/(2(r2g-rm)(r2g-rp))+((2a*EEg+Lzg(-2+r3g))r3g*\[ScriptCapitalI]r3g*\[Delta]r3)/(2(r3g-rm)(r3g-rp)))+\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]rg \[Delta]\[CapitalUpsilon]rfunFC[a,p,e,x]
]


(* ::Subsubsection::Closed:: *)
(*Fixed frequency*)


\[Delta]\[CapitalUpsilon]\[Phi]funFF[a_,p_,e_,x_,{\[Delta]r1_,\[Delta]r2_}]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g,r1g,r2g,r3g,\[Delta]r3,krg,ellK,ellE,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	(*{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];*)
	\[Delta]r3=\[Delta]r3funFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzgdr1gfun[a,p,e,x]\[Delta]r1+dLzgdr2gfun[a,p,e,x]\[Delta]r2;
	
	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));

	ellK=EllipticK[krg];
	ellE=EllipticE[krg];

	\[ScriptCapitalI]=(2ellK)/(Sqrt[1-EEg^2]*Sqrt[r2g(r1g-r3g)]);
	\[ScriptCapitalI]r1g=2/(Sqrt[1-EEg^2]*Sqrt[r2g(r1g-r3g)]r1g)((r2g *ellE-r1g*ellK)/(r1g-r2g));
	\[ScriptCapitalI]r2g=2/(Sqrt[1-EEg^2]*Sqrt[r2g (r1g-r3g)])(-(((r1g-r3g)*ellE)/((r1g-r2g)(r2g-r3g)))+ellK/(r1g-r2g));
	\[ScriptCapitalI]r3g=2/(Sqrt[1-EEg^2]*Sqrt[r2g (r1g-r3g)]r3g)((r2g*ellE)/(r2g-r3g)-ellK);

	(2\[CapitalUpsilon]rg)/(2\[Pi])(EEg(-1+(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+((2a*EEg+Lzg(-2+r1g))r1g*\[ScriptCapitalI]r1g*\[Delta]r1)/(2(r1g-rm)(r1g-rp))+((2a*EEg+Lzg(-2+r2g))r2g*\[ScriptCapitalI]r2g*\[Delta]r2)/(2(r2g-rm)(r2g-rp))+((2a*EEg+Lzg(-2+r3g))r3g*\[ScriptCapitalI]r3g*\[Delta]r3)/(2(r3g-rm)(r3g-rp)))+\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]rg*\[Delta]\[CapitalUpsilon]rfunFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}]
]


(* ::Subsubsection::Closed:: *)
(*Fixed eccentricity*)


\[Delta]\[CapitalUpsilon]\[Phi]funFE[a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,krg,ellK,ellE,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;

	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));

	ellK=EllipticK[krg];
	ellE=EllipticE[krg];

	\[ScriptCapitalI]=(2ellK)/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]);
	\[ScriptCapitalI]r1g=2 /(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]r1g) ((r2g *ellE-r1g*ellK)/(r1g-r2g) );
	\[ScriptCapitalI]r2g=2 /(Sqrt[1-EEg^2] Sqrt[r2g (r1g-r3g)]) (-(((r1g-r3g)*ellE)/((r1g-r2g)(r2g-r3g)))+ellK/(r1g-r2g));
	\[ScriptCapitalI]r3g=2 /(Sqrt[1-EEg^2] Sqrt[r2g (r1g-r3g)]r3g) ((r2g*ellE)/(r2g-r3g)-ellK);

	(2\[CapitalUpsilon]rg)/(2\[Pi])(EEg(-1+(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+((2a*EEg+Lzg(-2+r1g))r1g*\[ScriptCapitalI]r1g*\[Delta]r1)/(2(r1g-rm)(r1g-rp))+((2a*EEg+Lzg(-2+r2g))r2g*\[ScriptCapitalI]r2g*\[Delta]r2)/(2(r2g-rm)(r2g-rp))+((2a*EEg+Lzg(-2+r3g))r3g*\[ScriptCapitalI]r3g*\[Delta]r3)/(2(r3g-rm)(r3g-rp)))+\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]rg \[Delta]\[CapitalUpsilon]rfunFE[a,p,e,x]
]


\[Delta]\[CapitalUpsilon]\[Phi]funLimFE[a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,\[Delta]r1,\[Delta]r2,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];
	
	If[e==0,
		(EEg(-1+(Lzg*\[Delta]EE)/(1-EEg^2))+\[Delta]Lz-((4a^3*EEg-4a*EEg*r1g+Lzg (-2+r1g)^2*r1g+a^2*Lzg(-4+3 r1g))\[Delta]r1)/(2(a^2+(-2+r1g)r1g)^2))+\[CapitalUpsilon]\[Phi]gfunLim[a,p,e,x] \[Delta]\[CapitalUpsilon]rover\[CapitalUpsilon]rgfunLimFE[a,p,e,x],
		(EEg(-1+(Lzg*\[Delta]EE)/(1-EEg^2))+\[Delta]Lz+((2a*EEg+Lzg(-2+r1g))r1g*\[Delta]r1)/(2(r1g-rm)(r1g-rp)) (1/(-r1g+r2g))+((2a*EEg+Lzg(-2+r2g))r2g*\[Delta]r2)/(2(r2g-rm)(r2g-rp)) 1/(r1g-r2g)+((2a*EEg+Lzg(-2+r2g))r2g*\[Delta]r2)/(2(r2g-rm)(r2g-rp)) (-(1/r2g)))+\[CapitalUpsilon]\[Phi]gfunLim[a,p,e,x] \[Delta]\[CapitalUpsilon]rover\[CapitalUpsilon]rgfunLimFE[a,p,e,x]
	]
]


(* ::Subsubsection::Closed:: *)
(*R gauge*)


\[Delta]\[CapitalUpsilon]\[Phi]funR[a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,krg,ellK,ellE,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dLzdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funRPer[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRPer[a,p,e,x];
	\[Delta]r3=\[Delta]r2;

	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));

	ellK=EllipticK[krg];
	ellE=EllipticE[krg];

	\[ScriptCapitalI]=(2ellK)/(Sqrt[1-EEg^2]*Sqrt[r2g(r1g-r3g)]);
	\[ScriptCapitalI]r1g=2/(Sqrt[1-EEg^2]*Sqrt[r2g(r1g-r3g)]r1g)((r2g *ellE-r1g*ellK)/(r1g-r2g));
	\[ScriptCapitalI]r2g=2/(Sqrt[1-EEg^2]*Sqrt[r2g (r1g-r3g)])(-(((r1g-r3g)*ellE)/((r1g-r2g)(r2g-r3g)))+ellK/(r1g-r2g));
	\[ScriptCapitalI]r3g=2/(Sqrt[1-EEg^2]*Sqrt[r2g (r1g-r3g)]r3g)((r2g*ellE)/(r2g-r3g)-ellK);

	(2\[CapitalUpsilon]rg)/(2\[Pi])(EEg(-1+(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+((2a*EEg+Lzg(-2+r1g))r1g*\[ScriptCapitalI]r1g*\[Delta]r1)/(2(r1g-rm)(r1g-rp))+((2a*EEg+Lzg(-2+r2g))r2g*\[ScriptCapitalI]r2g*\[Delta]r2)/(2(r2g-rm)(r2g-rp))+((2a*EEg+Lzg(-2+r3g))r3g*\[ScriptCapitalI]r3g*\[Delta]r3)/(2(r3g-rm)(r3g-rp)))+\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]rg \[Delta]\[CapitalUpsilon]rfunR[a,p,e,x]
]


\[Delta]\[CapitalUpsilon]\[Phi]funLimR[a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,\[Delta]r1,\[Delta]r2,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dLzdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	\[Delta]r1=\[Delta]r1funRHom[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRHom[a,p,e,x];
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];
	
	If[e==0,
		(EEg(-1+(Lzg*\[Delta]EE)/(1-EEg^2))+\[Delta]Lz-((4a^3*EEg-4a*EEg*r1g+Lzg (-2+r1g)^2*r1g+a^2*Lzg(-4+3 r1g))\[Delta]r1)/(2(a^2+(-2+r1g)r1g)^2))+\[CapitalUpsilon]\[Phi]gfunLim[a,p,e,x] \[Delta]\[CapitalUpsilon]rover\[CapitalUpsilon]rgfunLimR[a,p,e,x],
		(EEg(-1+(Lzg*\[Delta]EE)/(1-EEg^2))+\[Delta]Lz+((2a*EEg+Lzg(-2+r1g))r1g*\[Delta]r1)/(2(r1g-rm)(r1g-rp)) (1/(-r1g+r2g))+((2a*EEg+Lzg(-2+r2g))r2g*\[Delta]r2)/(2(r2g-rm)(r2g-rp)) 1/(r1g-r2g)+((2a*EEg+Lzg(-2+r2g))r2g*\[Delta]r2)/(2(r2g-rm)(r2g-rp)) (-(1/r2g)))+\[CapitalUpsilon]\[Phi]gfunLim[a,p,e,x] \[Delta]\[CapitalUpsilon]rover\[CapitalUpsilon]rgfunLimR[a,p,e,x]
	]
]


(* ::Subsection::Closed:: *)
(*Shift azimuthal BL frequency bound plunges*)


\[Delta]\[CapitalOmega]\[Phi]funISCOplunge[r_,a_,x_]:=Module[{r1g,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]r,\[CapitalDelta]r,Lzred},
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,x];
	Lzg=Lzgfun[a,r1g,0,x];
	\[Delta]EE=\[Delta]EEISCOfun[a,x];
	\[Delta]Lz=\[Delta]LzISCOfun[a,x];
	\[Delta]r=\[Delta]rfunISCOplungePar[r,a,x];
	\[CapitalDelta]r=a^2+(-2+r)r ;
	Lzred=Lzg-a*EEg;

	(\[CapitalDelta]r(Lzred^2-Lzg*r^3\[Delta]EE+EEg*r^3(-EEg+\[Delta]Lz))-2r(a*Lzred^2-3EEg*Lzred*r^2+EEg*Lzg*r^3)\[Delta]r)/(r (r*EEg(a^2+r^2)-2a*Lzred)^2)
]


(* ::Subsubsection::Closed:: *)
(*Fixed eccentricity*)


\[Delta]\[CapitalOmega]\[Phi]funFECritplunge[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]r,\[CapitalDelta]r,Lzred},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]r=\[Delta]rfunFECritplungePar[r,a,p,e,x];
	\[CapitalDelta]r=a^2+(-2+r)r ;
	Lzred=Lzg-a*EEg;

	(\[CapitalDelta]r(Lzred^2-Lzg*r^3\[Delta]EE+EEg*r^3(-EEg+\[Delta]Lz))-2r(a*Lzred^2-3EEg*Lzred*r^2+EEg*Lzg*r^3)\[Delta]r)/(r (r*EEg(a^2+r^2)-2a*Lzred)^2)
]


\[Delta]\[CapitalOmega]\[Phi]funFEPlunge[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]EE,\[Delta]Lz,\[Delta]r,\[CapitalDelta]r,Lzred},
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	\[Delta]p=\[Delta]pfunFEcc[a,p,e,x];
	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x]+dEEdpfuncc[p,e]\[Delta]p;
	\[Delta]Lz=\[Delta]Lzfuncc[a,p,e,x]+dLzdpfuncc[a,p,e,x]\[Delta]p;
	\[Delta]r=\[Delta]rfunFEPlungePar[r,a,p,e,x];
	\[CapitalDelta]r=a^2+(-2+r)r ;
	Lzred=Lzg-a*EEg;

	(\[CapitalDelta]r(Lzred^2-Lzg*r^3\[Delta]EE+EEg*r^3(-EEg+\[Delta]Lz))-2r(a*Lzred^2-3EEg*Lzred*r^2+EEg*Lzg*r^3)\[Delta]r)/(r (r*EEg(a^2+r^2)-2a*Lzred)^2)
]


\[Delta]\[CapitalOmega]\[Phi]funFEPlungeBoundOrbit[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]EE,\[Delta]Lz,\[Delta]r,\[CapitalDelta]r,Lzred},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]p=\[Delta]pfunFEcc[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]p;
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]p;
	\[Delta]r=\[Delta]rfunFEPlungeBoundOrbitPar[r,a,p,e,x];
	\[CapitalDelta]r=a^2+(-2+r)r ;
	Lzred=Lzg-a*EEg;

	(\[CapitalDelta]r(Lzred^2-Lzg*r^3\[Delta]EE+EEg*r^3(-EEg+\[Delta]Lz))-2r(a*Lzred^2-3EEg*Lzred*r^2+EEg*Lzg*r^3)\[Delta]r)/(r (r*EEg(a^2+r^2)-2a*Lzred)^2)
]


\[Delta]\[CapitalOmega]\[Phi]funFEPlungeBoundOrbitDoubleRoot[r_,a_,p_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]EE,\[Delta]Lz,\[Delta]r,\[CapitalDelta]r,Lzred},
	EEg=EEgfun[a,p,0,x];
	Lzg=Lzgfun[a,p,0,x];
	\[Delta]p=\[Delta]pfunFEcc[a,p,0,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,0,x]+dEEdpfun[a,p,0,x]\[Delta]p;
	\[Delta]Lz=\[Delta]LzfunFT[a,p,0,x]+dLzdpfun[a,p,0,x]\[Delta]p;
	\[Delta]r=\[Delta]rfunFEPlungeBoundOrbitDoubleRootPar[r,a,p,x];
	\[CapitalDelta]r=a^2+(-2+r)r ;
	Lzred=Lzg-a*EEg;

	(\[CapitalDelta]r(Lzred^2-Lzg*r^3\[Delta]EE+EEg*r^3(-EEg+\[Delta]Lz))-2r(a*Lzred^2-3EEg*Lzred*r^2+EEg*Lzg*r^3)\[Delta]r)/(r (r*EEg(a^2+r^2)-2a*Lzred)^2)
]


(* ::Subsubsection::Closed:: *)
(*R gauge*)


\[Delta]\[CapitalOmega]\[Phi]funRCritplunge[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]e,\[Delta]EE,\[Delta]Lz,\[Delta]r,\[CapitalDelta]r,Lzred},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]p=\[Delta]pfunRHom[a,p,e,x];
	\[Delta]e=\[Delta]efunRHom[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]p+dEEdefun[a,p,e,x]\[Delta]e;
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]p+dLzdefun[a,p,e,x]\[Delta]e;
	\[Delta]r=\[Delta]rfunRCritplungePar[r,a,p,e,x];
	\[CapitalDelta]r=a^2+(-2+r)r ;
	Lzred=Lzg-a*EEg;

	(\[CapitalDelta]r(Lzred^2-Lzg*r^3\[Delta]EE+EEg*r^3(-EEg+\[Delta]Lz))-2r(a*Lzred^2-3EEg*Lzred*r^2+EEg*Lzg*r^3)\[Delta]r)/(r (r*EEg(a^2+r^2)-2a*Lzred)^2)
]


\[Delta]\[CapitalOmega]\[Phi]funRPlunge[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]e,\[Delta]EE,\[Delta]Lz,\[Delta]r,\[CapitalDelta]r,Lzred},
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	\[Delta]p=\[Delta]pfunRPlunge[a,p,e,x];
	\[Delta]e=\[Delta]efunRPlunge[a,p,e,x];
	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x]+dEEdpfuncc[p,e]\[Delta]p+dEEdefuncc[p,e]\[Delta]e;
	\[Delta]Lz=\[Delta]Lzfuncc[a,p,e,x]+dLzdpfuncc[a,p,e,x]\[Delta]p+dLzdefuncc[a,p,e,x]\[Delta]e;
	\[Delta]r=\[Delta]rfunRPlungePar[r,a,p,e,x];
	\[CapitalDelta]r=a^2+(-2+r)r ;
	Lzred=Lzg-a*EEg;

	(\[CapitalDelta]r(Lzred^2-Lzg*r^3\[Delta]EE+EEg*r^3(-EEg+\[Delta]Lz))-2r(a*Lzred^2-3EEg*Lzred*r^2+EEg*Lzg*r^3)\[Delta]r)/(r (r*EEg(a^2+r^2)-2a*Lzred)^2)
]


\[Delta]\[CapitalOmega]\[Phi]funRPlungeBoundOrbit[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]e,\[Delta]EE,\[Delta]Lz,\[Delta]r,\[CapitalDelta]r,Lzred},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]p=\[Delta]pfunRPer[a,p,e,x];
	\[Delta]e=\[Delta]efunRPer[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]p+dEEdefun[a,p,e,x]\[Delta]e;
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]p+dLzdefun[a,p,e,x]\[Delta]e;
	\[Delta]r=\[Delta]rfunRPlungeBoundOrbitPar[r,a,p,e,x];
	\[CapitalDelta]r=a^2+(-2+r)r ;
	Lzred=Lzg-a*EEg;

	(\[CapitalDelta]r(Lzred^2-Lzg*r^3\[Delta]EE+EEg*r^3(-EEg+\[Delta]Lz))-2r(a*Lzred^2-3EEg*Lzred*r^2+EEg*Lzg*r^3)\[Delta]r)/(r (r*EEg(a^2+r^2)-2a*Lzred)^2)
]


(* ::Chapter::Closed:: *)
(*Back-end: spin corrections to the orbit *)


(* ::Section::Closed:: *)
(*Spin corrections to the orbit - parallel component of the spin, fixed turning points*)


(* ::Subsection::Closed:: *)
(*Bound trajectory*)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunFTPerPar[wr_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,EEg,\[Delta]EE,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunFT[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r3=\[Delta]r3funFT[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Gamma]rg=(r1g-r2g)/(r1g-r3g);
	Yint=Sqrt[(r1g-r3g)r2g];

	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result.*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	jSN=Sin[\[Phi]];
	jCN=Cos[\[Phi]];

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];

	-((jSN*jCN(r1g-r2g)(r1g-r3g)(r2g-r3g)Sqrt[r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g])/((1-jSN^2) r1g+jSN^2 r2g-r3g)^2)(\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg (2ellF)/Yint+(EEg \[Delta]EE)/(1-EEg^2)(2ellF)/Yint+\[Delta]r3/2 1/r3g  2/Yint((r2g ellE)/(r2g-r3g)-ellF)+2((jSN*jCN(r1g-r2g))/(r1g*r2g Sqrt[r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g])-((r1g-r3g)ellE-r1g*ellF)/(r1g*r3g*Yint))\[Delta]r41-a/2 (-(2/(3r1g*r2g*r3g))((jSN*jCN(r1g-r2g)(r1g-r3g)(r2g-r3g)Sqrt[(r1g*r2g-jSN^2*r1g*r3g-r2g*r3g*jCN^2)])/(r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g)^2+((r1g+r2g+r3g)ellF)/Sqrt[r2g(r1g-r3g)]-2(r2g*r3g+r1g(r2g+r3g))((jSN*jCN(r1g-r2g))/(r1g*r2g Sqrt[r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g])+(-((r1g-r3g)ellE)+r1g*ellF)/(r1g*r3g Sqrt[r2g(r1g-r3g)] )))))
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunFTPerParRes[wr_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,EEg,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunFT[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r3=\[Delta]r3funFT[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Gamma]rg=(r1g-r2g)/(r1g-r3g);
	Yint=Sqrt[(r1g-r3g)r2g];
	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result. Having \[Phi]=\[Pi]/2 with MachinePrecision cause some problems with the integrals \[ScriptCapitalI]p and \[ScriptCapitalI]m*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	jSN=Sin[\[Phi]];
	jCN=Cos[\[Phi]];

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];

	-(\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg(2ellF)/Yint+(EEg \[Delta]EE)/(1-EEg^2)(2ellF)/Yint+\[Delta]r3/2 1/r3g  2/Yint((r2g ellE)/(r2g-r3g)-ellF)+2((jSN*jCN(r1g-r2g))/(r1g*r2g Sqrt[r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g])-((r1g-r3g)ellE-r1g*ellF)/(r1g*r3g*Yint))\[Delta]r41-a/2 (-(2/(3r1g*r2g*r3g))((jSN*jCN(r1g-r2g)(r1g-r3g)(r2g-r3g)Sqrt[(r1g*r2g-jSN^2*r1g*r3g-r2g*r3g*jCN^2)])/(r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g)^2+((r1g+r2g+r3g)ellF)/Sqrt[r2g(r1g-r3g)]-2(r2g*r3g+r1g(r2g+r3g))((jSN*jCN(r1g-r2g))/(r1g*r2g Sqrt[r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g])+(-((r1g-r3g)ellE)+r1g*ellF)/(r1g*r3g Sqrt[r2g(r1g-r3g)] )))))
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory - purely oscillatory part*)


\[Delta]tfunFTPerPar[wr_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[CapitalUpsilon]rg,\[CapitalUpsilon]tg,r1g,r2g,r3g,\[Delta]r3,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,\[Phi],rg,dtgd\[Lambda]fun,krghold,krg,\[Gamma]r,ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r3=\[Delta]r3funFT[a,p,e,x];
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];

	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));
	\[Gamma]r=(r1g-r2g)/(r1g-r3g);

	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result.*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	rg=rgICr2gfun[wr,a,p,e,x];

	dtgd\[Lambda]fun=(rg(-2a*Lzg+EEg(2a^2+a^2*rg+rg^3)))/(a^2-2rg+rg^2);

	ellF=EllipticF[\[Phi],krg];
	ellE= EllipticE[\[Phi],krg];
	ellPi= EllipticPi[\[Gamma]r,\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/(Sqrt[1-EEg^2]Sqrt[r2g(r1g-r3g)]);
	\[ScriptCapitalI]r3g=2/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]*r3g)((r2g ellE)/(r2g-r3g)-ellF);
	\[ScriptCapitalI]rover=(2(r3g*ellF+(r2g-r3g)ellPi))/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]);
	\[ScriptCapitalI]r2over=1/ Sqrt[(1-EEg^2)]((2r2g^2)/Sqrt[(r1g-r3g)r2g]((r3g/r2g)^2*ellF+2r3g/r2g(1-r3g/r2g)ellPi+(1-r3g/r2g)^2*1/(2(\[Gamma]r-1)(krg-\[Gamma]r))(\[Gamma]r*ellE+(krg-\[Gamma]r)ellF+(2\[Gamma]r krg+2\[Gamma]r-\[Gamma]r^2-3krg)ellPi-(\[Gamma]r^2Sin[\[Phi]]Cos[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(1-\[Gamma]r*Sin[\[Phi]]^2))));

	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	(4\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+(r3g*\[Delta]r3)/(2(r3g-rm)(r3g-rp))(-2a*Lzg+EEg*r3g^3+a^2*EEg(2+r3g))*\[ScriptCapitalI]r3g+EEg/2((2+r3g)\[Delta]r3-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+(EEg*\[Delta]r3)/2+EEg*\[Delta]\[Rho]r4)\[ScriptCapitalI]rover-1/\[CapitalUpsilon]rg(\[Delta]\[CapitalUpsilon]tfunFT[a,p,e,x]-\[CapitalUpsilon]tg/\[CapitalUpsilon]rg \[Delta]\[CapitalUpsilon]rfunFT[a,p,e,x])wr+dtgd\[Lambda]fun/Sqrt[1-EEg^2]\[Delta]rfunFTPerParRes[wr,a,p,e,x]
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal  trajectory - purely oscillatory part*)


\[Delta]\[Phi]funFTPerPar[wr_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g,r1g,r2g,r3g,\[Delta]r3,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,\[Phi],rg,d\[Phi]gd\[Lambda]fun,krghold,krg,\[Gamma]r,ellF,ellE,\[ScriptCapitalI],\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r3=\[Delta]r3funFT[a,p,e,x];
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];

	krg=((r1g-r2g) r3g)/(r2g (r1g-r3g));
	\[Gamma]r=(r1g-r2g)/(r1g-r3g);

	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result.*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	rg=rgICr2gfun[wr,a,p,e,x];

	d\[Phi]gd\[Lambda]fun=((2a*EEg+Lzg(rg-2))rg)/(a^2-2 rg+rg^2);

	ellF=EllipticF[\[Phi],krg];
	ellE= EllipticE[\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]);
	\[ScriptCapitalI]r3g=2 /(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]*r3g) ((r2g*ellE)/(r2g-r3g)-ellF);
	
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	-EEg(1-(Lzg \[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+((2a*EEg+Lzg(r3g-2))r3g*\[Delta]r3)/(2(r3g-rm)(r3g-rp)) \[ScriptCapitalI]r3g-1/\[CapitalUpsilon]rg(\[Delta]\[CapitalUpsilon]\[Phi]funFT[a,p,e,x]-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]rg*\[Delta]\[CapitalUpsilon]rfunFT[a,p,e,x])wr+d\[Phi]gd\[Lambda]fun/Sqrt[1-EEg^2]*\[Delta]rfunFTPerParRes[wr,a,p,e,x]
]


(* ::Subsection::Closed:: *)
(*Bound velocity*)


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunFTPerPar[wr_,a_,p_,e_,x_]:=Module[{EEg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,dRgdr,rg},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	\[Delta]r3=\[Delta]r3funFT[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	
	dRgdr=Sqrt[(1-EEg^2)](-4rg^3+r1g*r2g*r3g+3rg^2(r1g+r2g+r3g)-2rg(r2g*r3g+r1g(r2g+r3g)));
	rg=rgICr2gfun[wr,a,p,e,x];

	-Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(1-EEg^2)rg(r1g-rg)(rg-r2g)(rg-r3g)](\[Delta]r3/(2(rg-r3g))+\[Delta]r41/rg-a/(2rg^2)+(EEg*\[Delta]EE)/(1-EEg^2))+1/2dRgdr*\[Delta]rfunFTPerParRes[wr,a,p,e,x]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity*)


\[Delta]vtfunFTPerPar[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,r3g,rg,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x];

	rg=rgICr2gfun[wr,a,p,e,x];
	\[CapitalDelta]=a^2-2rg+rg^2;
	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+rg^2)+a*EEg(a^2+3rg^2))/(rg*\[CapitalDelta])+(rg(2a^2+a^2*rg+rg^3)\[Delta]EE)/\[CapitalDelta]-(2a*rg*\[Delta]Lz)/\[CapitalDelta];
	ddtgd\[Lambda]drfun=(EEg*rg(a^2+3rg^2))/\[CapitalDelta]-((rg^2-a^2)(-2a*Lzg+EEg*rg^3+a^2*EEg(2+rg)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunFTPerPar[wr,a,p,e,x]+d\[Delta]td\[Lambda]fun
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


\[Delta]v\[Phi]funFTPerPar[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,rg,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x];

	rg=rgICr2gfun[wr,a,p,e,x];
	\[CapitalDelta]=a^2-2rg+rg^2;
	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+rg)))/(rg*\[CapitalDelta])+(2a*rg*\[Delta]EE)/\[CapitalDelta]+((rg-2)rg*\[Delta]Lz)/\[CapitalDelta];
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(rg-1)-EEg*rg^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunFTPerPar[wr,a,p,e,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Section::Closed:: *)
(*Spin corrections to the orbit - parallel component of the spin, fixed constants of motion*)


(* ::Subsection::Closed:: *)
(*Bound trajectory*)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunFCPerPar[wr_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,EEg,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE,rg},
	EEg=EEgfun[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunFC[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFC[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFC[a,p,e,x];
	\[Delta]r3=\[Delta]r3funFC[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Gamma]rg=(r1g-r2g)/(r1g-r3g);
	Yint=Sqrt[(r1g-r3g)r2g];

	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result.*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	jSN=Sin[\[Phi]];
	jCN=Cos[\[Phi]];

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	
	rg=rgICr2gfun[wr,a,p,e,x];
	
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	-((jSN*jCN(r1g-r2g)(r1g-r3g)(r2g-r3g)Sqrt[r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g])/((1-jSN^2)r1g+jSN^2 r2g-r3g)^2)(\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg(2ellF)/Yint+(r2g*ellE-r1g*ellF)/(r1g(r1g-r2g)Yint)\[Delta]r1+(-(r1g-r3g)ellE/((r1g-r2g)(r2g-r3g))+ellF/(r1g-r2g))/Yint*\[Delta]r2+\[Delta]r3/2*1/r3g  2/Yint((r2g ellE)/(r2g-r3g)-ellF)+2((jSN*jCN(r1g-r2g))/(r1g*r2g Sqrt[r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g])-((r1g-r3g)ellE-r1g*ellF)/(r1g*r3g*Yint))\[Delta]r41-a/2 (-(2/(3r1g*r2g*r3g))((jSN*jCN(r1g-r2g)(r1g-r3g)(r2g-r3g)Sqrt[(r1g*r2g-jSN^2*r1g*r3g-r2g*r3g*jCN^2)])/(r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g)^2+((r1g+r2g+r3g)ellF)/Sqrt[r2g(r1g-r3g)]-2(r2g*r3g+r1g(r2g+r3g))((jSN*jCN(r1g-r2g))/(r1g*r2g Sqrt[r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g])+(-((r1g-r3g)ellE)+r1g*ellF)/(r1g*r3g Sqrt[r2g(r1g-r3g)])))))+(rg(r2g(-r2g+rg)\[Delta]r1+r1g(r1g-rg)\[Delta]r2))/(r1g(r1g-r2g)r2g)
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunFCPerParRes[wr_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,EEg,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE,rg},
	EEg=EEgfun[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunFC[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFC[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFC[a,p,e,x];
	\[Delta]r3=\[Delta]r3funFC[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Gamma]rg=(r1g-r2g)/(r1g-r3g);
	Yint=Sqrt[(r1g-r3g)r2g];
	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result. Having \[Phi]=\[Pi]/2 with MachinePrecision cause some problems with the integrals \[ScriptCapitalI]p and \[ScriptCapitalI]m*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	jSN=Sin[\[Phi]];
	jCN=Cos[\[Phi]];

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	
	rg=rgICr2gfun[wr,a,p,e,x];
	
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	-(\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg(2ellF)/Yint+(r2g*ellE-r1g*ellF)/(r1g(r1g-r2g)Yint)\[Delta]r1+(-(r1g-r3g)ellE/((r1g-r2g)(r2g-r3g))+ellF/(r1g-r2g))/Yint*\[Delta]r2+\[Delta]r3/2 1/r3g  2/Yint((r2g*ellE)/(r2g-r3g)-ellF)+2((jSN*jCN(r1g-r2g))/(r1g*r2g*Sqrt[r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g])-((r1g-r3g)ellE-r1g*ellF)/(r1g*r3g*Yint))\[Delta]r41-a/2(-(2/(3r1g*r2g*r3g))((jSN*jCN(r1g-r2g)(r1g-r3g)(r2g-r3g)Sqrt[(r1g*r2g-jSN^2*r1g*r3g-r2g*r3g*jCN^2)])/(r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g)^2+((r1g+r2g+r3g)ellF)/Sqrt[r2g(r1g-r3g)]-2(r2g*r3g+r1g(r2g+r3g))((jSN*jCN(r1g-r2g))/(r1g*r2g Sqrt[r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g])+(-((r1g-r3g)ellE)+r1g*ellF)/(r1g*r3g Sqrt[r2g(r1g-r3g)]))))+Sign[\[Pi]-Mod[wr,2 \[Pi]]](\[Delta]r1/(r1g(r1g-r2g)) (-((Sqrt[rg] Sqrt[rg-r2g])/(Sqrt[-rg+r1g] Sqrt[rg-r3g])))- \[Delta]r2/((r1g-r2g)r2g) ((Sqrt[rg] Sqrt[-rg+r1g] )/(Sqrt[rg-r2g] Sqrt[rg-r3g]))))
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory - purely oscillatory part*)


\[Delta]tfunFCPerPar[wr_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[CapitalUpsilon]rg,\[CapitalUpsilon]tg,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,\[Phi],rg,dtgd\[Lambda]fun,krghold,krg,\[Gamma]r,ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFC[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFC[a,p,e,x];
	\[Delta]r3=\[Delta]r3funFC[a,p,e,x];
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];

	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));
	\[Gamma]r=(r1g-r2g)/(r1g-r3g);

	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result.*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	rg=rgICr2gfun[wr,a,p,e,x];

	dtgd\[Lambda]fun=(rg(-2a*Lzg+EEg(2a^2+a^2*rg+rg^3)))/(a^2-2rg+rg^2);

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	ellPi=EllipticPi[\[Gamma]r,\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/(Sqrt[1-EEg^2] Sqrt[r2g (r1g-r3g)]);
	\[ScriptCapitalI]r1g=2 /(Sqrt[1-EEg^2] Sqrt[r2g (r1g-r3g)] r1g) ((r2g*ellE-r1g*ellF)/(r1g-r2g) );
	\[ScriptCapitalI]r2g=2 /(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]) (-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g));
	\[ScriptCapitalI]r3g=2 /(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]r3g) ((r2g*ellE)/(r2g-r3g)-ellF);

	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	1/2 EEg((2+r1g)\[Delta]r1+(2+r2g)\[Delta]r2+(2+r3g)\[Delta]r3-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+(r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(2(r1g-rm)(r1g-rp)) \[ScriptCapitalI]r1g*\[Delta]r1+(r2g(-2a*Lzg+EEg(2a^2+a^2*r2g+r2g^3)))/(2(r2g-rm)(r2g-rp)) \[ScriptCapitalI]r2g*\[Delta]r2+(r3g(-2a*Lzg+EEg(2a^2+a^2*r3g+r3g^3)) )/(2(r3g-rm)(r3g-rp)) \[ScriptCapitalI]r3g*\[Delta]r3-1/\[CapitalUpsilon]rg (\[Delta]\[CapitalUpsilon]tfunFC[a,p,e,x]-\[CapitalUpsilon]tg/\[CapitalUpsilon]rg \[Delta]\[CapitalUpsilon]rfunFC[a,p,e,x])wr+Sign[\[Pi]-Mod[wr,2 \[Pi]]] (r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(Sqrt[1-EEg^2](r1g-rm)(r1g-rp)) (-((Sqrt[rg] Sqrt[rg-r2g])/(Sqrt[-rg+r1g] Sqrt[rg-r3g])) 1/(r1g(r1g-r2g)))\[Delta]r1+Sign[\[Pi]-Mod[wr,2 \[Pi]]] (r2g(-2a*Lzg+EEg(2a^2+a^2*r2g+r2g^3)))/(Sqrt[1-EEg^2](r2g-rm)(r2g-rp)) (-((Sqrt[rg] Sqrt[-rg+r1g] )/(Sqrt[rg-r2g] Sqrt[rg-r3g])) 1/((r1g-r2g)r2g))\[Delta]r2+dtgd\[Lambda]fun/Sqrt[1-EEg^2] \[Delta]rfunFCPerParRes[wr,a,p,e,x]
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal  trajectory - purely oscillatory part*)


\[Delta]\[Phi]funFCPerPar[wr_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Phi],rg,d\[Phi]gd\[Lambda]fun,krghold,krg,\[Gamma]r,ellF,ellE,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFC[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFC[a,p,e,x];
	\[Delta]r3=\[Delta]r3funFC[a,p,e,x];

	krg=((r1g-r2g) r3g)/(r2g (r1g-r3g));
	\[Gamma]r=(r1g-r2g)/(r1g-r3g);

	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result.*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	rg=rgICr2gfun[wr,a,p,e,x];

	d\[Phi]gd\[Lambda]fun=((2a*EEg+Lzg(rg-2))rg)/(a^2-2rg+rg^2);

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/(Sqrt[1-EEg^2] Sqrt[r2g (r1g-r3g)]);
	\[ScriptCapitalI]r1g=2 /(Sqrt[1-EEg^2] Sqrt[r2g (r1g-r3g)] r1g) ((r2g*ellE-r1g*ellF)/(r1g-r2g) );
	\[ScriptCapitalI]r2g=2 /(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]) (-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g));
	\[ScriptCapitalI]r3g=2 /(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]r3g) ((r2g*ellE)/(r2g-r3g)-ellF);
	
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	-EEg*\[ScriptCapitalI]+((2a*EEg+Lzg(-2+r1g))r1g*\[Delta]r1)/(2(r1g-rm)(r1g-rp)) \[ScriptCapitalI]r1g+((2a*EEg+Lzg(-2+r2g))r2g*\[Delta]r2)/(2(r2g-rm)(r2g-rp)) \[ScriptCapitalI]r2g+((2a*EEg+Lzg(-2+r3g))r3g*\[Delta]r3)/(2(r3g-rm)(r3g-rp)) \[ScriptCapitalI]r3g-1/\[CapitalUpsilon]rg (\[Delta]\[CapitalUpsilon]\[Phi]funFC[a,p,e,x]-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]rg \[Delta]\[CapitalUpsilon]rfunFC[a,p,e,x])wr+Sign[\[Pi]-Mod[wr,2 \[Pi]]] ((2a*EEg+Lzg(-2+r1g))r1g*\[Delta]r1)/(Sqrt[1-EEg^2](r1g-rm)(r1g-rp)) (-((Sqrt[rg] Sqrt[rg-r2g])/(Sqrt[-rg+r1g] Sqrt[rg-r3g])) 1/(r1g(r1g-r2g))) +Sign[\[Pi]-Mod[wr,2 \[Pi]]] ((2a*EEg+Lzg(-2+r2g))r2g*\[Delta]r2)/(Sqrt[1-EEg^2](r2g-rm)(r2g-rp)) (-((Sqrt[rg] Sqrt[-rg+r1g] )/(Sqrt[rg-r2g] Sqrt[rg-r3g])) 1/((r1g-r2g)r2g))+d\[Phi]gd\[Lambda]fun/Sqrt[1-EEg^2] \[Delta]rfunFCPerParRes[wr,a,p,e,x]
]


(* ::Subsection::Closed:: *)
(*Bound velocity*)


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunFCPerPar[wr_,a_,p_,e_,x_]:=Module[{EEg,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,dRgdr,rg},
	EEg=EEgfun[a,p,e,x];
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFC[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFC[a,p,e,x];
	\[Delta]r3=\[Delta]r3funFC[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	dRgdr=Sqrt[(1-EEg^2)](-4rg^3+r1g*r2g*r3g+3rg^2(r1g+r2g+r3g)-2rg(r2g*r3g+r1g(r2g+r3g)));
	rg=rgICr2gfun[wr,a,p,e,x];

	-Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(1-EEg^2)rg(r1g-rg)(rg-r2g)(rg-r3g)](\[Delta]r1/(2(rg-r1g))+\[Delta]r2/(2(rg-r2g))+\[Delta]r3/(2(rg-r3g))+\[Delta]r41/rg-a/(2rg^2))+1/2dRgdr*\[Delta]rfunFCPerParRes[wr,a,p,e,x]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity*)


\[Delta]vtfunFCPerPar[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,rg,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];

	rg=rgICr2gfun[wr,a,p,e,x];	
	\[CapitalDelta]=a^2-2rg+rg^2;
	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+rg^2)+a*EEg(a^2+3rg^2))/(rg*\[CapitalDelta]);
	ddtgd\[Lambda]drfun=(EEg*rg(a^2+3rg^2))/\[CapitalDelta]-((rg^2-a^2)(-2a*Lzg+EEg*rg^3+a^2*EEg(2+rg)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunFCPerPar[wr,a,p,e,x]+d\[Delta]td\[Lambda]fun
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


\[Delta]v\[Phi]funFCPerPar[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,rg,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];

	rg=rgICr2gfun[wr,a,p,e,x];	
	\[CapitalDelta]=a^2-2rg+rg^2;
	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+rg)))/(rg*\[CapitalDelta]);
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(rg-1)-EEg*rg^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunFCPerPar[wr,a,p,e,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Section::Closed:: *)
(*Spin corrections to the orbit - parallel component of the spin, fixed frequency*)


(* ::Subsection::Closed:: *)
(*Bound trajectory*)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunFFPerPar[wr_,a_,p_,e_,x_,{\[Delta]r1_,\[Delta]r2_}]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,EEg,\[Delta]EE,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE,rg},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	(*{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];*)
	\[Delta]r3=\[Delta]r3funFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;
	
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Gamma]rg=(r1g-r2g)/(r1g-r3g);
	Yint=Sqrt[(r1g-r3g)r2g];

	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result.*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	jSN=Sin[\[Phi]];
	jCN=Cos[\[Phi]];

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	
	rg=rgICr1gfun[wr,a,p,e,x];
	
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	((jSN*jCN*r1g(r1g-r2g)r2g*Sqrt[r1g(r2g-jSN^2*r3g)-(1-jSN^2)r2g*r3g])/(jSN^2(r1g-r2g)+r2g)^2)(\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg*(2ellF)/Sqrt[(r1g-r3g)r2g]+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Sqrt[(r1g-r3g)r2g]+((r2g*ellE-r1g*ellF)\[Delta]r1)/(Sqrt[(r1g-r3g)r2g]r1g(r1g-r2g))+1/Sqrt[(r1g-r3g)r2g](-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g))\[Delta]r2+1/Sqrt[(r1g-r3g)r2g] 1/r3g(((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]))\[Delta]r3+2((-(r1g-r3g)ellE+r1g*ellF)/(r1g Sqrt[r2g(r1g-r3g)]r3g))\[Delta]r41-a/2 (-((4ellE Sqrt[r2g(r1g-r3g)](r2g*r3g+r1g(r2g+r3g)))/(3r1g^2r2g^2r3g^2))+(2ellF((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g Sqrt[r2g (r1g-r3g)]r3g^2)-(2(-r1g+r2g)(r1g-r3g)^2Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3r1g^2(r2g(r1g-r3g))^(3/2)r3g)))+(((rg-r2g)(rg-r3g))/((r1g-r2g)(r1g-r3g))\[Delta]r1+((r1g-rg)(rg-r3g))/((r1g-r2g)(r2g-r3g))\[Delta]r2)
]


\[Delta]rfunOfrFFPerPar[r_,a_,p_,e_,x_,{\[Delta]r1_,\[Delta]r2_}]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,EEg,\[Delta]EE,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],krg,ellF,ellE},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	(*{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];*)
	\[Delta]r3=\[Delta]r3funFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Phi]=ArcSin[Sqrt[(r2g(r1g-r))/((r1g-r2g)r)]];
	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	
	Sqrt[(r1g-r)(r-r2g)r(r-r3g)](\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg*(2ellF)/Sqrt[(r1g-r3g)r2g]+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Sqrt[(r1g-r3g)r2g]+(r2g*ellE-r1g*ellF)/(Sqrt[(r1g-r3g)r2g]r1g(r1g-r2g) ) \[Delta]r1+1/Sqrt[(r1g-r3g)r2g] (-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g))\[Delta]r2+1/Sqrt[(r1g-r3g)r2g] 1/r3g (((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]))\[Delta]r3+ 2((-(r1g-r3g)ellE+r1g*ellF)/(r1g*r3g Sqrt[r2g(r1g-r3g)] ))\[Delta]r41-a/2 (-((4ellE*Sqrt[r2g(r1g-r3g)](r2g*r3g+r1g(r2g+r3g)))/(3r1g^2r2g^2r3g^2))+(2ellF((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g*r3g^2 Sqrt[r2g (r1g-r3g)] )-(2(-r1g+r2g)(r1g-r3g)Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3r1g^2r2g*r3g Sqrt[r2g(r1g-r3g)])))+(((r-r2g)(r-r3g))/((r1g-r2g)(r1g-r3g) ) \[Delta]r1+((r1g-r)(r-r3g))/((r1g-r2g)(r2g-r3g))\[Delta]r2)
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunFFPerParRes[wr_,a_,p_,e_,x_,{\[Delta]r1_,\[Delta]r2_}]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,EEg,\[Delta]EE,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE,rg},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	(*{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];*)
	\[Delta]r3=\[Delta]r3funFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Gamma]rg=(r1g-r2g)/(r1g-r3g);
	Yint=Sqrt[(r1g-r3g)r2g];
	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result. Having \[Phi]=\[Pi]/2 with MachinePrecision cause some problems with the integrals \[ScriptCapitalI]p and \[ScriptCapitalI]m*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	jSN=Sin[\[Phi]];
	jCN=Cos[\[Phi]];

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	
	rg=rgICr1gfun[wr,a,p,e,x];
	
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg(2ellF)/Sqrt[(r1g-r3g)r2g]+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Sqrt[(r1g-r3g)r2g]+((r2g*ellE-r1g*ellF)\[Delta]r1)/(Sqrt[(r1g-r3g)r2g]r1g(r1g-r2g) )+1/Sqrt[(r1g-r3g)r2g](-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g))\[Delta]r2+1/Sqrt[(r1g-r3g)r2g] 1/r3g (((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]))\[Delta]r3+2((-(r1g-r3g)ellE+r1g*ellF)/(r1g*r3g*Sqrt[r2g(r1g-r3g)]))\[Delta]r41-a/2(-((4ellE Sqrt[r2g(r1g-r3g)](r2g*r3g+r1g(r2g+r3g)))/(3r1g^2r2g^2r3g^2))+(2ellF((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g*r3g^2 Sqrt[r2g(r1g-r3g)])-(2(-r1g+r2g)(r1g-r3g)Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3r1g^2r2g*r3g Sqrt[r2g(r1g-r3g)]))+Sign[\[Pi]-Mod[wr,2\[Pi]]](Sqrt[(-r2g+rg)(-r3g+rg)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-rg)rg])\[Delta]r1+Sqrt[(r1g-rg)(-r3g+rg)]/((r1g-r2g)(r2g-r3g)Sqrt[rg(-r2g+rg)])\[Delta]r2)
]


\[Delta]rfunOfrFFPerParRes[r_,a_,p_,e_,x_,{\[Delta]r1_,\[Delta]r2_}]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,EEg,\[Delta]EE,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],krg,\[Gamma]rg,ellF,ellE},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	(*{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];*)
	\[Delta]r3=\[Delta]r3funFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;

	
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Phi]=ArcSin[Sqrt[(r2g (r1g-r))/((r1g-r2g)r)]];
	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	
	(\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg*(2ellF)/Sqrt[(r1g-r3g)r2g]+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Sqrt[(r1g-r3g)r2g]+(r2g*ellE-r1g*ellF)/(Sqrt[(r1g-r3g)r2g]r1g(r1g-r2g))\[Delta]r1+1/Sqrt[(r1g-r3g)r2g](-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g))\[Delta]r2+1/Sqrt[(r1g-r3g)r2g] 1/r3g (((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]))\[Delta]r3+ 2((-(r1g-r3g)ellE+r1g*ellF)/(r1g*r3g*Sqrt[r2g(r1g-r3g)] ))\[Delta]r41-a/2(-((4ellE*Sqrt[r2g(r1g-r3g)](r2g*r3g+r1g(r2g+r3g)))/(3r1g^2r2g^2r3g^2))+(2ellF((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g*r3g^2*Sqrt[r2g (r1g-r3g)])-(2(-r1g+r2g)(r1g-r3g)Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3r1g^2r2g*r3g Sqrt[r2g(r1g-r3g)])))+(Sqrt[(-r2g+r)(-r3g+r)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-r)r])\[Delta]r1+Sqrt[(r1g-r)(-r3g+r)]/((r1g-r2g)(r2g-r3g)Sqrt[r(-r2g+r)])\[Delta]r2)
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory - purely oscillatory part*)


\[Delta]tfunFFPerPar[wr_,a_,p_,e_,x_,{\[Delta]r1_,\[Delta]r2_}]:=Module[{EEg,Lzg,\[Delta]EE,\[CapitalUpsilon]rg,\[CapitalUpsilon]tg,r1g,r2g,r3g,\[Delta]r3,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,\[Phi],rg,dtgd\[Lambda]fun,dtgd\[Lambda]funr1g,dtgd\[Lambda]funr2g,dtgd\[Lambda]funr3g,krghold,krg,\[Gamma]r,ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	(*{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];*)
	\[Delta]r3=\[Delta]r3funFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}];
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;

	krg=((r1g-r2g)r3g)/(r2g (r1g-r3g));
	\[Gamma]r=1-r1g/r2g;

	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result.*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	rg=rgICr1gfun[wr,a,p,e,x];

	dtgd\[Lambda]fun=(rg(-2a*Lzg+EEg(2a^2+a^2*rg+rg^3)))/(a^2-2rg+rg^2);
	dtgd\[Lambda]funr1g=(r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(a^2-2r1g+r1g^2);
	dtgd\[Lambda]funr2g=(r2g(-2a*Lzg+EEg(2a^2+a^2*r2g+r2g^3)))/(a^2-2r2g+r2g^2);
	dtgd\[Lambda]funr3g=(r3g(-2a*Lzg+EEg(2a^2+a^2*r3g+r3g^3)))/(a^2-2r3g+r3g^2);

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	ellPi=EllipticPi[\[Gamma]r,\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[r2g(r1g-r3g)];
	\[ScriptCapitalI]r1g=2/(Sqrt[r2g(r1g-r3g)]r1g)((r2g*ellE-r1g*ellF)/(r1g-r2g));
	\[ScriptCapitalI]r2g=2/(Sqrt[r2g (r1g-r3g)])(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g));
	\[ScriptCapitalI]r3g=2/(Sqrt[r2g(r1g-r3g)]r3g)((r2g*ellE)/(r2g-r3g)-ellF+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]));
	\[ScriptCapitalI]rover=(2r1g)/(Sqrt[(r1g-r3g)r2g])ellPi;
	\[ScriptCapitalI]r2over=(r2g((r1g-r3g)ellE-r1g*ellF))/(Sqrt[r2g(r1g-r3g)])+(r1g(r1g+r2g+r3g)ellPi)/(Sqrt[r2g(r1g-r3g)])+((r1g-r2g)Sqrt[r2g(r1g-r3g)]Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/((r2g+(r1g-r2g)Sin[\[Phi]]^2));

	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	1/Sqrt[1-EEg^2]((4\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+1/2 EEg(\[Delta]r1+\[Delta]r2+\[Delta]r3+2\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+1/2 EEg((2+r1g)\[Delta]r1+(2+r2g)\[Delta]r2+(2+r3g)\[Delta]r3-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+1/2dtgd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+1/2dtgd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2+1/2dtgd\[Lambda]funr3g*\[ScriptCapitalI]r3g*\[Delta]r3+Sign[\[Pi]-Mod[wr,2\[Pi]]]dtgd\[Lambda]funr1g*(Sqrt[(-r2g+rg)(-r3g+rg)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-rg)rg]))\[Delta]r1+Sign[\[Pi]-Mod[wr,2\[Pi]]]dtgd\[Lambda]funr2g*(Sqrt[(r1g-rg)(-r3g+rg)]/((r1g-r2g)(r2g-r3g)Sqrt[rg(-r2g+rg)]))\[Delta]r2-dtgd\[Lambda]fun*\[Delta]rfunFFPerParRes[wr,a,p,e,x,{\[Delta]r1,\[Delta]r2}])-1/\[CapitalUpsilon]rg(\[Delta]\[CapitalUpsilon]tfunFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}]-\[CapitalUpsilon]tg/\[CapitalUpsilon]rg*\[Delta]\[CapitalUpsilon]rfunFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}])wr
]	


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory as function of r*)


\[Delta]tfunOfrFFPerPar[r_,a_,p_,e_,x_,{\[Delta]r1_,\[Delta]r2_}]:=Module[{EEg,Lzg,\[Delta]EE,\[CapitalUpsilon]rg,\[CapitalUpsilon]tg,r1g,r2g,r3g,\[Delta]r3,\[Phi],\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,dtgd\[Lambda]fun,dtgd\[Lambda]funr1g,dtgd\[Lambda]funr2g,dtgd\[Lambda]funr3g,krg,\[Gamma]r,ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	(*{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];*)
	\[Delta]r3=\[Delta]r3funFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}];
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;

	krg=((r1g-r2g)r3g)/(r2g (r1g-r3g));
	\[Gamma]r=1-r1g/r2g;
	\[Phi]=ArcSin[Sqrt[(r2g(r1g-r))/((r1g-r2g)r)]];

	dtgd\[Lambda]fun=(r(-2a*Lzg+EEg(2a^2+a^2*r+r^3)))/(a^2-2r+r^2);
	dtgd\[Lambda]funr1g=(r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(a^2-2r1g+r1g^2);
	dtgd\[Lambda]funr2g=(r2g(-2a*Lzg+EEg(2a^2+a^2*r2g+r2g^3)))/(a^2-2r2g+r2g^2);
	dtgd\[Lambda]funr3g=(r3g(-2a*Lzg+EEg(2a^2+a^2*r3g+r3g^3)))/(a^2-2r3g+r3g^2);

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	ellPi=EllipticPi[\[Gamma]r,\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[r2g(r1g-r3g)];
	\[ScriptCapitalI]r1g=2/(Sqrt[r2g(r1g-r3g)]r1g)((r2g*ellE-r1g*ellF)/(r1g-r2g));
	\[ScriptCapitalI]r2g=2/(Sqrt[r2g (r1g-r3g)])(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g));
	\[ScriptCapitalI]r3g=2/(Sqrt[r2g(r1g-r3g)]r3g)((r2g*ellE)/(r2g-r3g)-ellF+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]));
	\[ScriptCapitalI]rover=(2r1g)/(Sqrt[(r1g-r3g)r2g])ellPi;
	\[ScriptCapitalI]r2over=(r2g((r1g-r3g)ellE-r1g*ellF))/(Sqrt[r2g(r1g-r3g)])+(r1g(r1g+r2g+r3g)ellPi)/(Sqrt[r2g(r1g-r3g)])+((r1g-r2g)Sqrt[r2g(r1g-r3g)]Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/((r2g+(r1g-r2g)Sin[\[Phi]]^2));

	1/Sqrt[1-EEg^2]((4\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+1/2 EEg(\[Delta]r1+\[Delta]r2+\[Delta]r3+2\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+1/2EEg((2+r1g)\[Delta]r1+(2+r2g)\[Delta]r2+(2+r3g)\[Delta]r3-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+1/2dtgd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+1/2dtgd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2+1/2dtgd\[Lambda]funr3g*\[ScriptCapitalI]r3g*\[Delta]r3+dtgd\[Lambda]funr1g*(Sqrt[(-r2g+r)(-r3g+r)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-r)r]))\[Delta]r1+dtgd\[Lambda]funr2g*(Sqrt[(r1g-r)(-r3g+r)]/((r1g-r2g)(r2g-r3g)Sqrt[r(-r2g+r)]))\[Delta]r2-dtgd\[Lambda]fun*\[Delta]rfunOfrFFPerParRes[r,a,p,e,x,{\[Delta]r1,\[Delta]r2}])
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal  trajectory - purely oscillatory part*)


\[Delta]\[Phi]funFFPerPar[wr_,a_,p_,e_,x_,{\[Delta]r1_,\[Delta]r2_}]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g,r1g,r2g,r3g,\[Delta]r3,\[Phi],rg,d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,d\[Phi]gd\[Lambda]funr2g,d\[Phi]gd\[Lambda]funr3g,krghold,krg,\[Gamma]r,ellF,ellE,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	(*{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];*)
	\[Delta]r3=\[Delta]r3funFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzgdr1gfun[a,p,e,x]\[Delta]r1+dLzgdr2gfun[a,p,e,x]\[Delta]r2;

	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));
	
	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result.*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	rg=rgICr1gfun[wr,a,p,e,x];

	d\[Phi]gd\[Lambda]fun=((2a*EEg+Lzg(rg-2))rg)/(a^2-2rg+rg^2);
	d\[Phi]gd\[Lambda]funr1g=((2a*EEg+Lzg(-2+r1g))r1g)/(a^2-2r1g+r1g^2);
	d\[Phi]gd\[Lambda]funr2g=((2a*EEg+Lzg(-2+r2g))r2g)/(a^2-2r2g+r2g^2);
	d\[Phi]gd\[Lambda]funr3g=((2a*EEg+Lzg(-2+r3g))r3g)/(a^2-2r3g+r3g^2);

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[r2g(r1g-r3g)];
	\[ScriptCapitalI]r1g=2/(Sqrt[r2g(r1g-r3g)]r1g)((r2g*ellE-r1g*ellF)/(r1g-r2g));
	\[ScriptCapitalI]r2g=2/(Sqrt[r2g(r1g-r3g)])(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g));
	\[ScriptCapitalI]r3g=2/(Sqrt[r2g(r1g-r3g)]r3g)(((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]));
	
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	1/Sqrt[1-EEg^2](-EEg(1-(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+1/2d\[Phi]gd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+1/2d\[Phi]gd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2+1/2d\[Phi]gd\[Lambda]funr3g*\[ScriptCapitalI]r3g*\[Delta]r3+Sign[\[Pi]-Mod[wr,2\[Pi]]]d\[Phi]gd\[Lambda]funr1g*\[Delta]r1*(Sqrt[(-r2g+rg)(-r3g+rg)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-rg)rg]))+Sign[\[Pi]-Mod[wr,2\[Pi]]]d\[Phi]gd\[Lambda]funr2g*\[Delta]r2*(Sqrt[(r1g-rg)(-r3g+rg)]/((r1g-r2g)(r2g-r3g)Sqrt[rg(-r2g+rg)]))-d\[Phi]gd\[Lambda]fun*\[Delta]rfunFFPerParRes[wr,a,p,e,x,{\[Delta]r1,\[Delta]r2}])-1/\[CapitalUpsilon]rg(\[Delta]\[CapitalUpsilon]\[Phi]funFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}]-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]rg \[Delta]\[CapitalUpsilon]rfunFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}])wr
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal  trajectory as function of r*)


\[Delta]\[Phi]funOfrFFPerPar[r_,a_,p_,e_,x_,{\[Delta]r1_,\[Delta]r2_}]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g,r1g,r2g,r3g,\[Delta]r3,\[Phi],d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,d\[Phi]gd\[Lambda]funr2g,d\[Phi]gd\[Lambda]funr3g,krg,\[Gamma]r,ellF,ellE,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	(*{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];*)
	\[Delta]r3=\[Delta]r3funFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzgdr1gfun[a,p,e,x]\[Delta]r1+dLzgdr2gfun[a,p,e,x]\[Delta]r2;

	krg=((r1g-r2g)r3g)/(r2g (r1g-r3g));
	\[Phi]=ArcSin[Sqrt[(r2g(r1g-r))/((r1g-r2g)r)]];

	d\[Phi]gd\[Lambda]fun=((2a*EEg+Lzg(r-2))r)/(a^2-2r+r^2);
	d\[Phi]gd\[Lambda]funr1g=((2a*EEg+Lzg(-2+r1g))r1g)/(a^2-2r1g+r1g^2);
	d\[Phi]gd\[Lambda]funr2g=((2a*EEg+Lzg(-2+r2g))r2g)/(a^2-2r2g+r2g^2);
	d\[Phi]gd\[Lambda]funr3g=((2a*EEg+Lzg(-2+r3g))r3g)/(a^2-2r3g+r3g^2);

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[r2g(r1g-r3g)];
	\[ScriptCapitalI]r1g=2/(Sqrt[r2g(r1g-r3g)]r1g)((r2g *ellE-r1g*ellF)/(r1g-r2g));
	\[ScriptCapitalI]r2g=2/(Sqrt[r2g(r1g-r3g)])(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g));
	\[ScriptCapitalI]r3g=2/(Sqrt[r2g(r1g-r3g)]r3g)(((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]));

	1/Sqrt[1-EEg^2](-EEg(1-(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+1/2d\[Phi]gd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+1/2d\[Phi]gd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2+1/2d\[Phi]gd\[Lambda]funr3g*\[ScriptCapitalI]r3g*\[Delta]r3+d\[Phi]gd\[Lambda]funr1g*\[Delta]r1(Sqrt[(-r2g+r)(-r3g+r)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-r)r]))+d\[Phi]gd\[Lambda]funr2g*\[Delta]r2(Sqrt[(r1g-r)(-r3g+r)]/((r1g-r2g)(r2g-r3g)Sqrt[r(-r2g+r)]))-d\[Phi]gd\[Lambda]fun*\[Delta]rfunOfrFFPerParRes[r,a,p,e,x,{\[Delta]r1,\[Delta]r2}])
]


(* ::Subsection::Closed:: *)
(*Bound velocity*)


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunFFPerPar[wr_,a_,p_,e_,x_,{\[Delta]r1_,\[Delta]r2_}]:=Module[{EEg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,dRgdr,rg},	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	(*{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];*)
	\[Delta]r3=\[Delta]r3funFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;

	dRgdr=Sqrt[(1-EEg^2)](-4rg^3+r1g*r2g*r3g+3rg^2(r1g+r2g+r3g)-2rg(r2g*r3g+r1g(r2g+r3g)));
	rg=rgICr1gfun[wr,a,p,e,x];

	-Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(1-EEg^2)rg(r1g-rg)(rg-r2g)(rg-r3g)](\[Delta]r1/(2(rg-r1g))+\[Delta]r2/(2(rg-r2g))+\[Delta]r3/(2(rg-r3g))+\[Delta]r41/rg-a/(2rg^2)+(EEg*\[Delta]EE)/(1-EEg^2))+1/2dRgdr*\[Delta]rfunFFPerParRes[wr,a,p,e,x,{\[Delta]r1,\[Delta]r2}]
]


\[Delta]vrfunOfrFFPerPar[r_,a_,p_,e_,x_,{\[Delta]r1_,\[Delta]r2_}]:=Module[{EEg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,dRgdr},
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	(*{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];*)
	\[Delta]r3=\[Delta]r3funFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;
	
	dRgdr=Sqrt[(1-EEg^2)](-4r^3+r1g*r2g*r3g+3r^2(r1g+r2g+r3g)-2r(r2g*r3g+r1g(r2g+r3g)));
	
	-Sqrt[(1-EEg^2)r(r1g-r)(r-r2g)(r-r3g)](\[Delta]r1/(2(r-r1g))+\[Delta]r2/(2(r-r2g))+\[Delta]r3/(2(r-r3g))+\[Delta]r41/r-a/(2r^2)+(EEg \[Delta]EE)/(1-EEg^2))+1/2dRgdr*\[Delta]rfunOfrFFPerParRes[r,a,p,e,x,{\[Delta]r1,\[Delta]r2}]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity*)


\[Delta]vtfunFFPerPar[wr_,a_,p_,e_,x_,{\[Delta]r1_,\[Delta]r2_}]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,rg,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	(*{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];*)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzgdr1gfun[a,p,e,x]\[Delta]r1+dLzgdr2gfun[a,p,e,x]\[Delta]r2;

	rg=rgICr1gfun[wr,a,p,e,x];	
	\[CapitalDelta]=a^2-2rg+rg^2;
	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+rg^2)+a*EEg(a^2+3rg^2))/(rg*\[CapitalDelta])+(rg(2a^2+a^2*rg+rg^3)\[Delta]EE)/\[CapitalDelta]-(2a*rg*\[Delta]Lz)/\[CapitalDelta];	
	ddtgd\[Lambda]drfun=(EEg*rg(a^2+3rg^2))/\[CapitalDelta]-((rg^2-a^2)(-2a*Lzg+EEg*rg^3+a^2*EEg(2+rg)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunFFPerPar[wr,a,p,e,x,{\[Delta]r1,\[Delta]r2}]+d\[Delta]td\[Lambda]fun
]


\[Delta]vtfunOfrFFPerPar[r_,a_,p_,e_,x_,{\[Delta]r1_,\[Delta]r2_}]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	(*{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];*)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzgdr1gfun[a,p,e,x]\[Delta]r1+dLzgdr2gfun[a,p,e,x]\[Delta]r2;

	\[CapitalDelta]=a^2-2r+r^2;
	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+r^2)+a*EEg(a^2+3r^2))/(r*\[CapitalDelta])+(r(2a^2+a^2*r+r^3)\[Delta]EE)/\[CapitalDelta]-(2a*r*\[Delta]Lz)/\[CapitalDelta];	
	ddtgd\[Lambda]drfun=(EEg*r(a^2+3r^2))/\[CapitalDelta]-((r^2-a^2)(-2a*Lzg+EEg*r^3+a^2*EEg(2+r)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunOfrFFPerPar[r,a,p,e,x,{\[Delta]r1,\[Delta]r2}]+d\[Delta]td\[Lambda]fun
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


\[Delta]v\[Phi]funFFPerPar[wr_,a_,p_,e_,x_,{\[Delta]r1_,\[Delta]r2_}]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,rg,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	(*{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];*)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzgdr1gfun[a,p,e,x]\[Delta]r1+dLzgdr2gfun[a,p,e,x]\[Delta]r2;

	rg=rgICr1gfun[wr,a,p,e,x];
	\[CapitalDelta]=a^2-2rg+rg^2;
	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+rg)))/(rg*\[CapitalDelta])+(2a*rg*\[Delta]EE)/\[CapitalDelta]+((rg-2)rg*\[Delta]Lz)/\[CapitalDelta];	
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(rg-1)-EEg*rg^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunFFPerPar[wr,a,p,e,x,{\[Delta]r1,\[Delta]r2}]+d\[Delta]\[Phi]d\[Lambda]fun
]


\[Delta]v\[Phi]funOfrFFPerPar[r_,a_,p_,e_,x_,{\[Delta]r1_,\[Delta]r2_}]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	(*{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];*)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzgdr1gfun[a,p,e,x]\[Delta]r1+dLzgdr2gfun[a,p,e,x]\[Delta]r2;

	\[CapitalDelta]=a^2-2r+r^2;
	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+r)))/(r*\[CapitalDelta])+(2a*r*\[Delta]EE)/\[CapitalDelta]+((r-2)r*\[Delta]Lz)/\[CapitalDelta];	
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(r-1)-EEg*r^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunOfrFFPerPar[r,a,p,e,x,{\[Delta]r1,\[Delta]r2}]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Section::Closed:: *)
(*Spin corrections to the orbit - parallel component of the spin, fixed eccentricity*)


(* ::Text:: *)
(*The ISCO plunge is universal to any spin gauge. I implemented here for convenience*)


(* ::Subsection::Closed:: *)
(*Bound trajectory*)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunFEPerPar[wr_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,\[Delta]r1,\[Delta]r2,EEg,\[Delta]EE,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE,rg},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunFE[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Gamma]rg=(r1g-r2g)/(r1g-r3g);
	Yint=Sqrt[(r1g-r3g)r2g];

	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result.*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	jSN=Sin[\[Phi]];
	jCN=Cos[\[Phi]];

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	
	rg=rgICr1gfun[wr,a,p,e,x];
	
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	((jSN*jCN*r1g(r1g-r2g)r2g*Sqrt[r1g(r2g-jSN^2*r3g)-(1-jSN^2)r2g*r3g])/(jSN^2(r1g-r2g)+r2g)^2)(\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg*(2ellF)/Sqrt[(r1g-r3g)r2g]+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Sqrt[(r1g-r3g)r2g]+((r2g*ellE-r1g*ellF)\[Delta]r1)/(Sqrt[(r1g-r3g)r2g]r1g(r1g-r2g))+1/Sqrt[(r1g-r3g)r2g](-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g))\[Delta]r2+1/Sqrt[(r1g-r3g)r2g] 1/r3g(((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]))\[Delta]r3+2((-(r1g-r3g)ellE+r1g*ellF)/(r1g Sqrt[r2g(r1g-r3g)]r3g))\[Delta]r41-a/2 (-((4ellE Sqrt[r2g(r1g-r3g)](r2g*r3g+r1g(r2g+r3g)))/(3r1g^2r2g^2r3g^2))+(2ellF((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g Sqrt[r2g (r1g-r3g)]r3g^2)-(2(-r1g+r2g)(r1g-r3g)^2Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3r1g^2(r2g(r1g-r3g))^(3/2)r3g)))+(((rg-r2g)(rg-r3g))/((r1g-r2g)(r1g-r3g))\[Delta]r1+((r1g-rg)(rg-r3g))/((r1g-r2g)(r2g-r3g))\[Delta]r2)
]


\[Delta]rfunOfrFEPerPar[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,EEg,\[Delta]EE,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],krg,ellF,ellE},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunFE[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Phi]=ArcSin[Sqrt[(r2g(r1g-r))/((r1g-r2g)r)]];
	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	
	Sqrt[(r1g-r)(r-r2g)r(r-r3g)](\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg*(2ellF)/Sqrt[(r1g-r3g)r2g]+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Sqrt[(r1g-r3g)r2g]+(r2g*ellE-r1g*ellF)/(Sqrt[(r1g-r3g)r2g]r1g(r1g-r2g) ) \[Delta]r1+1/Sqrt[(r1g-r3g)r2g] (-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g))\[Delta]r2+1/Sqrt[(r1g-r3g)r2g] 1/r3g (((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]))\[Delta]r3+ 2((-(r1g-r3g)ellE+r1g*ellF)/(r1g*r3g Sqrt[r2g(r1g-r3g)] ))\[Delta]r41-a/2 (-((4ellE*Sqrt[r2g(r1g-r3g)](r2g*r3g+r1g(r2g+r3g)))/(3r1g^2r2g^2r3g^2))+(2ellF((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g*r3g^2 Sqrt[r2g (r1g-r3g)] )-(2(-r1g+r2g)(r1g-r3g)Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3r1g^2r2g*r3g Sqrt[r2g(r1g-r3g)])))+(((r-r2g)(r-r3g))/((r1g-r2g)(r1g-r3g) ) \[Delta]r1+((r1g-r)(r-r3g))/((r1g-r2g)(r2g-r3g))\[Delta]r2)
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunFEPerParRes[wr_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,\[Delta]r1,\[Delta]r2,EEg,\[Delta]EE,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE,rg},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunFE[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Gamma]rg=(r1g-r2g)/(r1g-r3g);
	Yint=Sqrt[(r1g-r3g)r2g];
	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result. Having \[Phi]=\[Pi]/2 with MachinePrecision cause some problems with the integrals \[ScriptCapitalI]p and \[ScriptCapitalI]m*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	jSN=Sin[\[Phi]];
	jCN=Cos[\[Phi]];

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	
	rg=rgICr1gfun[wr,a,p,e,x];
	
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg(2ellF)/Sqrt[(r1g-r3g)r2g]+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Sqrt[(r1g-r3g)r2g]+((r2g*ellE-r1g*ellF)\[Delta]r1)/(Sqrt[(r1g-r3g)r2g]r1g(r1g-r2g) )+1/Sqrt[(r1g-r3g)r2g](-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g))\[Delta]r2+1/Sqrt[(r1g-r3g)r2g] 1/r3g (((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]))\[Delta]r3+2((-(r1g-r3g)ellE+r1g*ellF)/(r1g*r3g*Sqrt[r2g(r1g-r3g)]))\[Delta]r41-a/2(-((4ellE Sqrt[r2g(r1g-r3g)](r2g*r3g+r1g(r2g+r3g)))/(3r1g^2r2g^2r3g^2))+(2ellF((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g*r3g^2 Sqrt[r2g(r1g-r3g)])-(2(-r1g+r2g)(r1g-r3g)Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3r1g^2r2g*r3g Sqrt[r2g(r1g-r3g)]))+Sign[\[Pi]-Mod[wr,2\[Pi]]](Sqrt[(-r2g+rg)(-r3g+rg)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-rg)rg])\[Delta]r1+Sqrt[(r1g-rg)(-r3g+rg)]/((r1g-r2g)(r2g-r3g)Sqrt[rg(-r2g+rg)])\[Delta]r2)
]


\[Delta]rfunOfrFEPerParRes[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,EEg,\[Delta]EE,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],krg,\[Gamma]rg,ellF,ellE},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunFE[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Phi]=ArcSin[Sqrt[(r2g (r1g-r))/((r1g-r2g)r)]];
	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	
	(\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg*(2ellF)/Sqrt[(r1g-r3g)r2g]+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Sqrt[(r1g-r3g)r2g]+(r2g*ellE-r1g*ellF)/(Sqrt[(r1g-r3g)r2g]r1g(r1g-r2g))\[Delta]r1+1/Sqrt[(r1g-r3g)r2g](-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g))\[Delta]r2+1/Sqrt[(r1g-r3g)r2g] 1/r3g (((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]))\[Delta]r3+ 2((-(r1g-r3g)ellE+r1g*ellF)/(r1g*r3g*Sqrt[r2g(r1g-r3g)] ))\[Delta]r41-a/2(-((4ellE*Sqrt[r2g(r1g-r3g)](r2g*r3g+r1g(r2g+r3g)))/(3r1g^2r2g^2r3g^2))+(2ellF((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g*r3g^2*Sqrt[r2g (r1g-r3g)])-(2(-r1g+r2g)(r1g-r3g)Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3r1g^2r2g*r3g Sqrt[r2g(r1g-r3g)])))+(Sqrt[(-r2g+r)(-r3g+r)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-r)r])\[Delta]r1+Sqrt[(r1g-r)(-r3g+r)]/((r1g-r2g)(r2g-r3g)Sqrt[r(-r2g+r)])\[Delta]r2)
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory - purely oscillatory part*)


\[Delta]tfunFEPerPar[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[CapitalUpsilon]rg,\[CapitalUpsilon]tg,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,\[Phi],rg,dtgd\[Lambda]fun,dtgd\[Lambda]funr1g,dtgd\[Lambda]funr2g,dtgd\[Lambda]funr3g,krghold,krg,\[Gamma]r,ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];

	krg=((r1g-r2g)r3g)/(r2g (r1g-r3g));
	\[Gamma]r=1-r1g/r2g;

	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result.*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	rg=rgICr1gfun[wr,a,p,e,x];

	dtgd\[Lambda]fun=(rg(-2a*Lzg+EEg(2a^2+a^2*rg+rg^3)))/(a^2-2rg+rg^2);
	dtgd\[Lambda]funr1g=(r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(a^2-2r1g+r1g^2);
	dtgd\[Lambda]funr2g=(r2g(-2a*Lzg+EEg(2a^2+a^2*r2g+r2g^3)))/(a^2-2r2g+r2g^2);
	dtgd\[Lambda]funr3g=(r3g(-2a*Lzg+EEg(2a^2+a^2*r3g+r3g^3)))/(a^2-2r3g+r3g^2);

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	ellPi=EllipticPi[\[Gamma]r,\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[r2g(r1g-r3g)];
	\[ScriptCapitalI]r1g=2/(Sqrt[r2g(r1g-r3g)]r1g)((r2g*ellE-r1g*ellF)/(r1g-r2g));
	\[ScriptCapitalI]r2g=2/(Sqrt[r2g (r1g-r3g)])(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g));
	\[ScriptCapitalI]r3g=2/(Sqrt[r2g(r1g-r3g)]r3g)((r2g*ellE)/(r2g-r3g)-ellF+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]));
	\[ScriptCapitalI]rover=(2r1g)/(Sqrt[(r1g-r3g)r2g])ellPi;
	\[ScriptCapitalI]r2over=(r2g((r1g-r3g)ellE-r1g*ellF))/(Sqrt[r2g(r1g-r3g)])+(r1g(r1g+r2g+r3g)ellPi)/(Sqrt[r2g(r1g-r3g)])+((r1g-r2g)Sqrt[r2g(r1g-r3g)]Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/((r2g+(r1g-r2g)Sin[\[Phi]]^2));

	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	1/Sqrt[1-EEg^2]((4\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+1/2 EEg(\[Delta]r1+\[Delta]r2+\[Delta]r3+2\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+1/2 EEg((2+r1g)\[Delta]r1+(2+r2g)\[Delta]r2+(2+r3g)\[Delta]r3-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+1/2dtgd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+1/2dtgd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2+1/2dtgd\[Lambda]funr3g*\[ScriptCapitalI]r3g*\[Delta]r3+Sign[\[Pi]-Mod[wr,2\[Pi]]]dtgd\[Lambda]funr1g*(Sqrt[(-r2g+rg)(-r3g+rg)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-rg)rg]))\[Delta]r1+Sign[\[Pi]-Mod[wr,2\[Pi]]]dtgd\[Lambda]funr2g*(Sqrt[(r1g-rg)(-r3g+rg)]/((r1g-r2g)(r2g-r3g)Sqrt[rg(-r2g+rg)]))\[Delta]r2-dtgd\[Lambda]fun*\[Delta]rfunFEPerParRes[wr,a,p,e,x])-1/\[CapitalUpsilon]rg(\[Delta]\[CapitalUpsilon]tfunFE[a,p,e,x]-\[CapitalUpsilon]tg/\[CapitalUpsilon]rg*\[Delta]\[CapitalUpsilon]rfunFE[a,p,e,x])wr
]	


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory as function of r*)


\[Delta]tfunOfrFEPerPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[CapitalUpsilon]rg,\[CapitalUpsilon]tg,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Phi],\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,dtgd\[Lambda]fun,dtgd\[Lambda]funr1g,dtgd\[Lambda]funr2g,dtgd\[Lambda]funr3g,krg,\[Gamma]r,ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];

	krg=((r1g-r2g)r3g)/(r2g (r1g-r3g));
	\[Gamma]r=1-r1g/r2g;
	\[Phi]=ArcSin[Sqrt[(r2g(r1g-r))/((r1g-r2g)r)]];

	dtgd\[Lambda]fun=(r(-2a*Lzg+EEg(2a^2+a^2*r+r^3)))/(a^2-2r+r^2);
	dtgd\[Lambda]funr1g=(r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(a^2-2r1g+r1g^2);
	dtgd\[Lambda]funr2g=(r2g(-2a*Lzg+EEg(2a^2+a^2*r2g+r2g^3)))/(a^2-2r2g+r2g^2);
	dtgd\[Lambda]funr3g=(r3g(-2a*Lzg+EEg(2a^2+a^2*r3g+r3g^3)))/(a^2-2r3g+r3g^2);

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	ellPi=EllipticPi[\[Gamma]r,\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[r2g(r1g-r3g)];
	\[ScriptCapitalI]r1g=2/(Sqrt[r2g(r1g-r3g)]r1g)((r2g*ellE-r1g*ellF)/(r1g-r2g));
	\[ScriptCapitalI]r2g=2/(Sqrt[r2g (r1g-r3g)])(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g));
	\[ScriptCapitalI]r3g=2/(Sqrt[r2g(r1g-r3g)]r3g)((r2g*ellE)/(r2g-r3g)-ellF+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]));
	\[ScriptCapitalI]rover=(2r1g)/(Sqrt[(r1g-r3g)r2g])ellPi;
	\[ScriptCapitalI]r2over=(r2g((r1g-r3g)ellE-r1g*ellF))/(Sqrt[r2g(r1g-r3g)])+(r1g(r1g+r2g+r3g)ellPi)/(Sqrt[r2g(r1g-r3g)])+((r1g-r2g)Sqrt[r2g(r1g-r3g)]Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/((r2g+(r1g-r2g)Sin[\[Phi]]^2));

	1/Sqrt[1-EEg^2]((4\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+1/2 EEg(\[Delta]r1+\[Delta]r2+\[Delta]r3+2\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+1/2EEg((2+r1g)\[Delta]r1+(2+r2g)\[Delta]r2+(2+r3g)\[Delta]r3-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+1/2dtgd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+1/2dtgd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2+1/2dtgd\[Lambda]funr3g*\[ScriptCapitalI]r3g*\[Delta]r3+dtgd\[Lambda]funr1g*(Sqrt[(-r2g+r)(-r3g+r)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-r)r]))\[Delta]r1+dtgd\[Lambda]funr2g*(Sqrt[(r1g-r)(-r3g+r)]/((r1g-r2g)(r2g-r3g)Sqrt[r(-r2g+r)]))\[Delta]r2-dtgd\[Lambda]fun*\[Delta]rfunOfrFEPerParRes[r,a,p,e,x])
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal  trajectory - purely oscillatory part*)


\[Delta]\[Phi]funFEPerPar[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Phi],rg,d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,d\[Phi]gd\[Lambda]funr2g,d\[Phi]gd\[Lambda]funr3g,krghold,krg,\[Gamma]r,ellF,ellE,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;

	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));
	
	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result.*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	rg=rgICr1gfun[wr,a,p,e,x];

	d\[Phi]gd\[Lambda]fun=((2a*EEg+Lzg(rg-2))rg)/(a^2-2rg+rg^2);
	d\[Phi]gd\[Lambda]funr1g=((2a*EEg+Lzg(-2+r1g))r1g)/(a^2-2r1g+r1g^2);
	d\[Phi]gd\[Lambda]funr2g=((2a*EEg+Lzg(-2+r2g))r2g)/(a^2-2r2g+r2g^2);
	d\[Phi]gd\[Lambda]funr3g=((2a*EEg+Lzg(-2+r3g))r3g)/(a^2-2r3g+r3g^2);

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[r2g(r1g-r3g)];
	\[ScriptCapitalI]r1g=2/(Sqrt[r2g(r1g-r3g)]r1g)((r2g*ellE-r1g*ellF)/(r1g-r2g));
	\[ScriptCapitalI]r2g=2/(Sqrt[r2g(r1g-r3g)])(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g));
	\[ScriptCapitalI]r3g=2/(Sqrt[r2g(r1g-r3g)]r3g)(((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]));
	
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	1/Sqrt[1-EEg^2](-EEg(1-(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+1/2d\[Phi]gd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+1/2d\[Phi]gd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2+1/2d\[Phi]gd\[Lambda]funr3g*\[ScriptCapitalI]r3g*\[Delta]r3+Sign[\[Pi]-Mod[wr,2\[Pi]]]d\[Phi]gd\[Lambda]funr1g*\[Delta]r1*(Sqrt[(-r2g+rg)(-r3g+rg)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-rg)rg]))+Sign[\[Pi]-Mod[wr,2\[Pi]]]d\[Phi]gd\[Lambda]funr2g*\[Delta]r2*(Sqrt[(r1g-rg)(-r3g+rg)]/((r1g-r2g)(r2g-r3g)Sqrt[rg(-r2g+rg)]))-d\[Phi]gd\[Lambda]fun*\[Delta]rfunFEPerParRes[wr,a,p,e,x])-1/\[CapitalUpsilon]rg(\[Delta]\[CapitalUpsilon]\[Phi]funFE[a,p,e,x]-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]rg \[Delta]\[CapitalUpsilon]rfunFE[a,p,e,x])wr
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal  trajectory as function of r*)


\[Delta]\[Phi]funOfrFEPerPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Phi],d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,d\[Phi]gd\[Lambda]funr2g,d\[Phi]gd\[Lambda]funr3g,krg,\[Gamma]r,ellF,ellE,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;

	krg=((r1g-r2g)r3g)/(r2g (r1g-r3g));
	\[Phi]=ArcSin[Sqrt[(r2g(r1g-r))/((r1g-r2g)r)]];

	d\[Phi]gd\[Lambda]fun=((2a*EEg+Lzg(r-2))r)/(a^2-2r+r^2);
	d\[Phi]gd\[Lambda]funr1g=((2a*EEg+Lzg(-2+r1g))r1g)/(a^2-2r1g+r1g^2);
	d\[Phi]gd\[Lambda]funr2g=((2a*EEg+Lzg(-2+r2g))r2g)/(a^2-2r2g+r2g^2);
	d\[Phi]gd\[Lambda]funr3g=((2a*EEg+Lzg(-2+r3g))r3g)/(a^2-2r3g+r3g^2);

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[r2g(r1g-r3g)];
	\[ScriptCapitalI]r1g=2/(Sqrt[r2g(r1g-r3g)]r1g)((r2g *ellE-r1g*ellF)/(r1g-r2g));
	\[ScriptCapitalI]r2g=2/(Sqrt[r2g(r1g-r3g)])(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g));
	\[ScriptCapitalI]r3g=2/(Sqrt[r2g(r1g-r3g)]r3g)(((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]));

	1/Sqrt[1-EEg^2](-EEg(1-(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+1/2d\[Phi]gd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+1/2d\[Phi]gd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2+1/2d\[Phi]gd\[Lambda]funr3g*\[ScriptCapitalI]r3g*\[Delta]r3+d\[Phi]gd\[Lambda]funr1g*\[Delta]r1(Sqrt[(-r2g+r)(-r3g+r)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-r)r]))+d\[Phi]gd\[Lambda]funr2g*\[Delta]r2(Sqrt[(r1g-r)(-r3g+r)]/((r1g-r2g)(r2g-r3g)Sqrt[r(-r2g+r)]))-d\[Phi]gd\[Lambda]fun*\[Delta]rfunOfrFEPerParRes[r,a,p,e,x])
]


(* ::Subsection::Closed:: *)
(*Bound velocity*)


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunFEPerPar[wr_,a_,p_,e_,x_]:=Module[{EEg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,dRgdr,rg},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	dRgdr=Sqrt[(1-EEg^2)](-4rg^3+r1g*r2g*r3g+3rg^2(r1g+r2g+r3g)-2rg(r2g*r3g+r1g(r2g+r3g)));
	rg=rgICr1gfun[wr,a,p,e,x];

	-Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(1-EEg^2)rg(r1g-rg)(rg-r2g)(rg-r3g)](\[Delta]r1/(2(rg-r1g))+\[Delta]r2/(2(rg-r2g))+\[Delta]r3/(2(rg-r3g))+\[Delta]r41/rg-a/(2rg^2)+(EEg*\[Delta]EE)/(1-EEg^2))+1/2dRgdr*\[Delta]rfunFEPerParRes[wr,a,p,e,x]
]


\[Delta]vrfunOfrFEPerPar[r_,a_,p_,e_,x_]:=Module[{EEg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,dRgdr},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	dRgdr=Sqrt[(1-EEg^2)](-4r^3+r1g*r2g*r3g+3r^2(r1g+r2g+r3g)-2r(r2g*r3g+r1g(r2g+r3g)));
	
	-Sqrt[(1-EEg^2)r(r1g-r)(r-r2g)(r-r3g)](\[Delta]r1/(2(r-r1g))+\[Delta]r2/(2(r-r2g))+\[Delta]r3/(2(r-r3g))+\[Delta]r41/r-a/(2r^2)+(EEg \[Delta]EE)/(1-EEg^2))+1/2dRgdr*\[Delta]rfunOfrFEPerParRes[r,a,p,e,x]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity*)


\[Delta]vtfunFEPerPar[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,rg,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];

	rg=rgICr1gfun[wr,a,p,e,x];
	\[CapitalDelta]=a^2-2rg+rg^2;
	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+rg^2)+a*EEg(a^2+3rg^2))/(rg*\[CapitalDelta])+(rg(2a^2+a^2*rg+rg^3)\[Delta]EE)/\[CapitalDelta]-(2a*rg*\[Delta]Lz)/\[CapitalDelta];
	ddtgd\[Lambda]drfun=(EEg*rg(a^2+3rg^2))/\[CapitalDelta]-((rg^2-a^2)(-2a*Lzg+EEg*rg^3+a^2*EEg(2+rg)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunFEPerPar[wr,a,p,e,x]+d\[Delta]td\[Lambda]fun
]


\[Delta]vtfunOfrFEPerPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];

	\[CapitalDelta]=a^2-2r+r^2;
	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+r^2)+a*EEg(a^2+3r^2))/(r*\[CapitalDelta])+(r(2a^2+a^2*r+r^3)\[Delta]EE)/\[CapitalDelta]-(2a*r*\[Delta]Lz)/\[CapitalDelta];
	ddtgd\[Lambda]drfun=(EEg*r(a^2+3r^2))/\[CapitalDelta]-((r^2-a^2)(-2a*Lzg+EEg*r^3+a^2*EEg(2+r)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunOfrFEPerPar[r,a,p,e,x]+d\[Delta]td\[Lambda]fun
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


\[Delta]v\[Phi]funFEPerPar[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,rg,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];

	rg=rgICr1gfun[wr,a,p,e,x];
	\[CapitalDelta]=a^2-2rg+rg^2;
	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+rg)))/(rg*\[CapitalDelta])+(2a*rg*\[Delta]EE)/\[CapitalDelta]+((rg-2)rg*\[Delta]Lz)/\[CapitalDelta];	
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(rg-1)-EEg*rg^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunFEPerPar[wr,a,p,e,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


\[Delta]v\[Phi]funOfrFEPerPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];

	\[CapitalDelta]=a^2-2r+r^2;
	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+r)))/(r*\[CapitalDelta])+(2a*r*\[Delta]EE)/\[CapitalDelta]+((r-2)r*\[Delta]Lz)/\[CapitalDelta];
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(r-1)-EEg*r^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunOfrFEPerPar[r,a,p,e,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Subsection::Closed:: *)
(*Homoclinic orbits*)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunFEHomPar[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,\[Delta]r3,\[Delta]r41,\[Delta]r1,\[Delta]r2,EEg,\[Delta]EE,\[Delta]Lz,\[ScriptCapitalI],\[ScriptCapitalI]r,\[ScriptCapitalI]r2},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	If[e==0,
		\[Delta]r1
		,
		\[ScriptCapitalI]=-(1/(2Sqrt[(r1g-r2g)r2g]))Log[(Sqrt[r(r1g-r2g)]-Sqrt[(r1g-r)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(r1g-r)r2g])^2];
		\[ScriptCapitalI]r=-2Sqrt[r1g-r]/(r1g*r2g*Sqrt[r])+\[ScriptCapitalI]/r2g;
		\[ScriptCapitalI]r2=-((2Sqrt[r(-r+r1g)])/(3r*r1g*r2g))(2/r1g+3/r2g+1/r)+\[ScriptCapitalI]/r2g^2;

		(r-r2g)Sqrt[(r1g-r)r](\[Delta]\[CapitalUpsilon]rover\[CapitalUpsilon]rgfunLimFE[a,p,e,x]\[ScriptCapitalI]+(EEg*\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+(r2g*Sqrt[-r+r1g])/(r1g (r1g-r2g)^2*Sqrt[r])\[Delta]r1-\[ScriptCapitalI]/(2(r1g-r2g))\[Delta]r1-(\[ScriptCapitalI](r1g-2r2g)\[Delta]r2)/(2(r1g-r2g)r2g)+1/((r1g-r2g)r2g)Sqrt[-1+r1g/r]\[Delta]r2+\[Delta]r41*\[ScriptCapitalI]r-a/2*\[ScriptCapitalI]r2)+(r-r2g)((r-r2g)/(r1g-r2g)^2*\[Delta]r1)+Sqrt[(r1g-r)](Sqrt[r1g-r]/(r1g-r2g)\[Delta]r2)
	]
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunFEHomParRes[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,\[Delta]r3,\[Delta]r41,\[Delta]r1,\[Delta]r2,EEg,\[Delta]EE,\[Delta]Lz,\[ScriptCapitalI],\[ScriptCapitalI]r,\[ScriptCapitalI]r2},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	\[ScriptCapitalI]=-(1/(2Sqrt[(r1g-r2g)r2g]))Log[(Sqrt[r(r1g-r2g)]-Sqrt[(r1g-r)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(r1g-r)r2g])^2];
	\[ScriptCapitalI]r=-2Sqrt[r1g-r]/(r1g*r2g*Sqrt[r])+\[ScriptCapitalI]/r2g;
	\[ScriptCapitalI]r2=-((2Sqrt[r(-r+r1g)])/(3r*r1g*r2g))(2/r1g+3/r2g+1/r)+\[ScriptCapitalI]/r2g^2;

	(\[Delta]\[CapitalUpsilon]rover\[CapitalUpsilon]rgfunLimFE[a,p,e,x]\[ScriptCapitalI]+(EEg*\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+(r2g*Sqrt[-r+r1g])/(r1g (r1g-r2g)^2*Sqrt[r])\[Delta]r1-\[ScriptCapitalI]/(2(r1g-r2g))\[Delta]r1-(\[ScriptCapitalI](r1g-2r2g)\[Delta]r2)/(2(r1g-r2g)r2g)+1/((r1g-r2g)r2g)Sqrt[-1+r1g/r]\[Delta]r2+\[Delta]r41*\[ScriptCapitalI]r-a/2*\[ScriptCapitalI]r2)+((r-r2g)/(Sqrt[(r1g-r)r](r1g-r2g)^2)\[Delta]r1)+(Sqrt[(r1g-r)r]/((r1g-r2g)(r-r2g)r)\[Delta]r2)
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory *)


\[Delta]tfunFEHomPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,dtgd\[Lambda]fun,dtgd\[Lambda]funr1g,dtgd\[Lambda]funr2g,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];
	
	If[e==0,
		0,
		dtgd\[Lambda]fun=(r(-2a*Lzg+EEg(2a^2+a^2*r+r^3)))/(a^2-2r+r^2);
		dtgd\[Lambda]funr1g=(r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(a^2-2r1g+r1g^2);
		dtgd\[Lambda]funr2g=(r2g(-2a*Lzg+EEg(2a^2+a^2*r2g+r2g^3)))/(a^2-2r2g+r2g^2);
		
		\[ScriptCapitalI]=-(1/(2Sqrt[(r1g-r2g)r2g]))Log[(Sqrt[r(r1g-r2g)]-Sqrt[(r1g-r)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(r1g-r)r2g])^2];
		\[ScriptCapitalI]rover=2ArcCos[Sqrt[r/r1g]]+r2g*\[ScriptCapitalI];
		\[ScriptCapitalI]r2over=Sqrt[r(-r+r1g)]+(r1g+2r2g)ArcCos[Sqrt[r/r1g]]+r2g^2*\[ScriptCapitalI];
		\[ScriptCapitalI]r1g=2/(r1g(r1g-r2g))Sqrt[r/(r1g-r)]-\[ScriptCapitalI]/(r1g-r2g);
		\[ScriptCapitalI]r2g=(Sqrt[r(r1g-r)])/((r-r2g)r2g(r1g-r2g))-((r1g-2r2g)\[ScriptCapitalI])/(2(r1g-r2g)r2g);
		
		1/Sqrt[1-EEg^2]((4*\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+1/2EEg(\[Delta]r1+\[Delta]r2+\[Delta]r2+2\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+1/2EEg((2+r1g)\[Delta]r1+(2+r2g)\[Delta]r2+(2+r2g)\[Delta]r2-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+(1/2)dtgd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+dtgd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2-dtgd\[Lambda]fun*\[Delta]rfunFEHomParRes[r,a,p,e,x])
	]
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal trajectory*)


\[Delta]\[Phi]funFEHomPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,d\[Phi]gd\[Lambda]funr2g,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	
	If[e==0,
		0,
		d\[Phi]gd\[Lambda]fun=((2a*EEg+Lzg(-2+r))r)/(a^2-2r+r^2);
		d\[Phi]gd\[Lambda]funr1g=((2a*EEg+Lzg(-2+r1g))r1g)/(a^2-2r1g+r1g^2);
		d\[Phi]gd\[Lambda]funr2g=((2a*EEg+Lzg(-2+r2g))r2g)/(a^2-2r2g+r2g^2);
		
		\[ScriptCapitalI]=-(1/(2Sqrt[(r1g-r2g)r2g]))Log[(Sqrt[r(r1g-r2g)]-Sqrt[(r1g-r)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(r1g-r)r2g])^2];
		\[ScriptCapitalI]r1g=2/(r1g(r1g-r2g)) Sqrt[r/(r1g-r)]-\[ScriptCapitalI]/(r1g-r2g);
		\[ScriptCapitalI]r2g=(Sqrt[r(r1g-r)])/((r-r2g)r2g(r1g-r2g))-((r1g-2r2g)\[ScriptCapitalI])/(2(r1g-r2g)r2g);

		1/Sqrt[1-EEg^2](-EEg(1-(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+(1/2)d\[Phi]gd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+d\[Phi]gd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2-d\[Phi]gd\[Lambda]fun*\[Delta]rfunFEHomParRes[r,a,p,e,x])
	]
]


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunFEHomPar[r_,a_,p_,e_,x_]:=Module[{EEg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,\[Phi],dRgdr},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	dRgdr=Sqrt[(1-EEg^2)](r2g-r)(r1g(r2g-3r)-2(r2g-2r)r);
	
	If[e==0,
		0,
		-(r-r2g)Sqrt[(1-EEg^2)r(r1g-r)](\[Delta]r1/(2(r-r1g))+\[Delta]r2/(r-r2g)+\[Delta]r41/r-a/(2r^2)+(EEg*\[Delta]EE)/(1-EEg^2))+1/2dRgdr*\[Delta]rfunFEHomParRes[r,a,p,e,x]
	]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity*)


\[Delta]vtfunFEHomPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[CapitalDelta]=a^2-2r+r^2;

	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+r^2)+a*EEg(a^2+3r^2))/(r*\[CapitalDelta])+(r(2a^2+a^2*r+r^3)\[Delta]EE)/\[CapitalDelta]-(2a*r*\[Delta]Lz)/\[CapitalDelta];
	ddtgd\[Lambda]drfun=(EEg*r(a^2+3r^2))/\[CapitalDelta]-((r^2-a^2)(-2a*Lzg+EEg*r^3+a^2*EEg(2+r)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunFEHomPar[r,a,p,e,x]+d\[Delta]td\[Lambda]fun
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


\[Delta]v\[Phi]funFEHomPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[CapitalDelta]=a^2-2r+r^2;

	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+r)))/(r*\[CapitalDelta])+(2a*r*\[Delta]EE)/\[CapitalDelta]+((r-2)r*\[Delta]Lz)/\[CapitalDelta];
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(r-1)-EEg*r^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunFEHomPar[r,a,p,e,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Subsection::Closed:: *)
(*ISCO plunge*)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunISCOplungePar[r_,a_,x_]:=Module[{r1g,\[Delta]r41,EEg,\[Delta]EE,\[Delta]rISCO},
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,1];
	\[Delta]rISCO=\[Delta]rISCOfun[a,x];
	\[Delta]r41=\[Delta]r41fun[a,r1g,0,x];

	-((a (r1g-r)^2(r1g+5r))/(3r1g^3*r))+r/r1g*\[Delta]rISCO+(2(r1g-r)^2*\[Delta]r41)/r1g^2
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunISCOplungeParRes[r_,a_,x_]:=Module[{r1g,\[Delta]r41,EEg,\[Delta]rISCO},
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,1];
	\[Delta]rISCO=\[Delta]rISCOfun[a,x];
	\[Delta]r41=\[Delta]r41fun[a,r1g,0,x];

	(a Sqrt[r1g-r](r1g+5r))/(3r1g^3*r^(3/2))-(Sqrt[r]\[Delta]rISCO)/((r1g-r)^(3/2)r1g)-(2Sqrt[r1g-r]\[Delta]r41)/(r1g^2*Sqrt[r])
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory *)


\[Delta]tfunISCOplungePar[r_,a_,x_]:=Module[{r1g,rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]rISCO,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,dtgd\[Lambda]fun,dtgd\[Lambda]funr1g,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	r1g=ISCOradius[a,x];
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,r1g,0,1];
	Lzg=Lzgfun[a,r1g,0,1];
	\[Delta]EE=\[Delta]EEISCOfun[a,x];
	\[Delta]Lz=\[Delta]LzISCOfun[a,x];
	\[Delta]rISCO=\[Delta]rISCOfun[a,x];
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,r1g,0,x];
	\[Delta]\[Rho]i4=Sqrt[a];

	dtgd\[Lambda]fun=(r(-2a*Lzg+EEg(2a^2+a^2*r+r^3)))/(a^2-2r+r^2);
	dtgd\[Lambda]funr1g=(r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(a^2-2r1g+r1g^2);

	\[ScriptCapitalI]=2Sqrt[r]/(r1g*Sqrt[r1g-r]);
	\[ScriptCapitalI]r1g=(2Sqrt[r](2r-3r1g))/(3Sqrt[(r1g-r)^3]r1g^2);
	\[ScriptCapitalI]rover=(2Sqrt[r])/Sqrt[r1g-r]+2ArcCos[Sqrt[r/r1g]];
	\[ScriptCapitalI]r2over=(Sqrt[r(r1g-r)]+r1g*ArcCos[Sqrt[r/r1g]])+r1g*\[ScriptCapitalI]rover;

	-1/Sqrt[(1-EEg^2)]((4\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+1/2EEg(3(2+r1g)\[Delta]rISCO-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+EEg((3\[Delta]rISCO)/2+\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+3/2dtgd\[Lambda]funr1g*\[Delta]rISCO*\[ScriptCapitalI]r1g-dtgd\[Lambda]fun*\[Delta]rfunISCOplungeParRes[r,a,x])
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal trajectory*)


\[Delta]\[Phi]funISCOplungePar[r_,a_,x_]:=Module[{r1g,rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]rISCO,d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,\[ScriptCapitalI],\[ScriptCapitalI]r1g},
	r1g=ISCOradius[a,x];
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,r1g,0,1];
	Lzg=Lzgfun[a,r1g,0,1];
	\[Delta]EE=\[Delta]EEISCOfun[a,x];
	\[Delta]Lz=\[Delta]LzISCOfun[a,x];
	\[Delta]rISCO=\[Delta]rISCOfun[a,x];

	d\[Phi]gd\[Lambda]fun=((2a*EEg+Lzg(-2+r))r)/(a^2-2r+r^2);
	d\[Phi]gd\[Lambda]funr1g=((2a*EEg+Lzg(-2+r1g))r1g)/(a^2-2r1g+r1g^2);
	
	\[ScriptCapitalI]=2Sqrt[r]/(r1g*Sqrt[r1g-r]);
	\[ScriptCapitalI]r1g=(2Sqrt[r](2r-3r1g))/(3Sqrt[(-r+r1g)^3]r1g^2);

	-1/Sqrt[(1-EEg^2)](-EEg(1-(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+3/2d\[Phi]gd\[Lambda]funr1g*\[Delta]rISCO*\[ScriptCapitalI]r1g-d\[Phi]gd\[Lambda]fun*\[Delta]rfunISCOplungeParRes[r,a,x])
]


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunISCOplungePar[r_,a_,x_]:=Module[{r1g,\[Delta]r41,EEg,\[Delta]EE,\[Delta]rISCO,dRgdr},
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,1];
	\[Delta]EE=\[Delta]EEISCOfun[a,x];
	\[Delta]rISCO=\[Delta]rISCOfun[a,x];
	\[Delta]r41=\[Delta]r41fun[a,r1g,0,x];
	dRgdr=Sqrt[(1-EEg^2)](r1g-4r)(r1g-r)^2;

	Sqrt[(1-EEg^2)r (r1g-r)^3]((3\[Delta]rISCO)/(2(r-r1g))+\[Delta]r41/r-a/(2r^2)+(EEg*\[Delta]EE)/(1-EEg^2))+1/2dRgdr*\[Delta]rfunISCOplungeParRes[r,a,x]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity*)


\[Delta]vtfunISCOplungePar[r_,a_,x_]:=Module[{r1g,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,1];
	Lzg=Lzgfun[a,r1g,0,1];
	\[Delta]EE=\[Delta]EEISCOfun[a,x];
	\[Delta]Lz=\[Delta]LzISCOfun[a,x];

	\[CapitalDelta]=a^2-2r+r^2;

	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+r^2)+a*EEg(a^2+3r^2))/(r*\[CapitalDelta])+(r(2a^2+a^2*r+r^3)\[Delta]EE)/\[CapitalDelta]-(2a*r*\[Delta]Lz)/\[CapitalDelta];
	ddtgd\[Lambda]drfun=(EEg*r(a^2+3r^2))/\[CapitalDelta]-((r^2-a^2)(-2a*Lzg+EEg*r^3+a^2*EEg(2+r)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunISCOplungePar[r,a,x]+d\[Delta]td\[Lambda]fun
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


\[Delta]v\[Phi]funISCOplungePar[r_,a_,x_]:=Module[{r1g,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,1];
	Lzg=Lzgfun[a,r1g,0,1];
	\[Delta]EE=\[Delta]EEISCOfun[a,x];
	\[Delta]Lz=\[Delta]LzISCOfun[a,x];

	\[CapitalDelta]=a^2-2r+r^2;

	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+r)))/(r*\[CapitalDelta])+(2a*r*\[Delta]EE)/\[CapitalDelta]+((r-2)r*\[Delta]Lz)/\[CapitalDelta];
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(r-1)-EEg*r^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunISCOplungePar[r,a,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Subsection::Closed:: *)
(*Critical plunge *)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunFECritplungePar[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,\[Delta]r41,\[Delta]r1,\[Delta]r2,EEg},
	EEg=EEgfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	(r-r2g)Sqrt[r(-r+r1g)](-((2Sqrt[r1g-r])/(r1g*r2g*Sqrt[r]))\[Delta]r41-a/2(-((2Sqrt[r(-r+r1g)])/(3r*r1g*r2g))(2/r1g+3/r2g+1/r)))+(r(r-r2g)\[Delta]r1)/(r1g(r1g-r2g))+(r(-r+r1g)\[Delta]r2)/((r1g-r2g)r2g)
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunFECritplungeParRes[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,\[Delta]r41,\[Delta]r1,\[Delta]r2,EEg,Lzg},
	EEg=EEgfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	 -((2Sqrt[r1g-r])/(r1g*r2g*Sqrt[r]))\[Delta]r41-a/2 (-((2Sqrt[r(-r+r1g)])/(3r*r1g*r2g))(2/r1g+3/r2g+1/r))+ \[Delta]r1/(r1g(r1g-r2g))Sqrt[r/(-r+r1g)]+(Sqrt[r(-r+r1g)]\[Delta]r2)/((r-r2g)r2g(r1g-r2g))
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory *)


\[Delta]tfunFECritplungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,dtgd\[Lambda]fun,dtgd\[Lambda]funr1g,dtgd\[Lambda]funr2g,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];
	
	dtgd\[Lambda]fun=(r(-2a*Lzg+EEg(2a^2+a^2*r+r^3)))/(a^2-2r+r^2);
	dtgd\[Lambda]funr1g=(r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(a^2-2r1g+r1g^2);
	dtgd\[Lambda]funr2g=(r2g(-2a*Lzg+EEg(2a^2+a^2*r2g+r2g^3)))/(a^2-2r2g+r2g^2);
	\[ScriptCapitalI]=-1/(2Sqrt[(r1g-r2g)r2g])Log[(Sqrt[r(r1g-r2g)]-Sqrt[(r1g-r)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(r1g-r)r2g])^2];
	\[ScriptCapitalI]r1g=2/(r1g(r1g-r2g))Sqrt[r/(r1g-r)]-\[ScriptCapitalI]/(r1g-r2g);
	\[ScriptCapitalI]r2g=(Sqrt[r(r1g-r)])/((r-r2g)r2g(r1g-r2g))-((r1g-2r2g)\[ScriptCapitalI])/(2(r1g-r2g)r2g);
	\[ScriptCapitalI]rover=2ArcCos[Sqrt[r/r1g]]+r2g*\[ScriptCapitalI];
	\[ScriptCapitalI]r2over=Sqrt[r(r1g-r)]+(r1g+2r2g)ArcCos[Sqrt[r/r1g]]+r2g^2\[ScriptCapitalI];

	-1/Sqrt[1-EEg^2]((4\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+1/2*EEg(\[Delta]r1+2\[Delta]r2+2\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+1/2*EEg((2+r1g)\[Delta]r1+2(2+r2g)\[Delta]r2-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+(1/2)dtgd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+dtgd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2-dtgd\[Lambda]fun*\[Delta]rfunFECritplungeParRes[r,a,p,e,x])
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal  trajectory*)


\[Delta]\[Phi]funFECritplungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,\[Delta]r1,\[Delta]r2,d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,d\[Phi]gd\[Lambda]funr2g,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	
	d\[Phi]gd\[Lambda]fun=((2a*EEg+Lzg(-2+r))r)/(a^2-2r+r^2);
	d\[Phi]gd\[Lambda]funr1g=((2a*EEg+Lzg(-2+r1g))r1g)/(a^2-2r1g+r1g^2);
	d\[Phi]gd\[Lambda]funr2g=((2a*EEg+Lzg(-2+r2g))r2g)/(a^2-2r2g+r2g^2);
	\[ScriptCapitalI]=-1/(2Sqrt[(r1g-r2g)r2g])Log[(Sqrt[r(r1g-r2g)]-Sqrt[(r1g-r)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(r1g-r)r2g])^2];
	\[ScriptCapitalI]r1g=2/(r1g(r1g-r2g)) Sqrt[r/(r1g-r)]-\[ScriptCapitalI]/(r1g-r2g);
	\[ScriptCapitalI]r2g=(Sqrt[r(r1g-r)])/((r-r2g)r2g(r1g-r2g))-((r1g-2r2g)\[ScriptCapitalI])/(2(r1g-r2g)r2g);

	-1/Sqrt[1-EEg^2](-EEg(1-(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+(1/2)d\[Phi]gd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+d\[Phi]gd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2-d\[Phi]gd\[Lambda]fun*\[Delta]rfunFECritplungeParRes[r,a,p,e,x])
]


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunFECritplungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,r1g,r2g,\[Delta]r1,\[Delta]r2,\[Delta]r41,dRgdr},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	dRgdr=Sqrt[(1-EEg^2)](r2g-r)(r1g(r2g-3r)-2(r2g-2r)r);
	
	-(r-r2g)Sqrt[(1-EEg^2)r(r1g-r)](\[Delta]r1/(2(r-r1g))+\[Delta]r2/(r-r2g)+\[Delta]r41/r-a/(2r^2)+(EEg*\[Delta]EE)/(1-EEg^2))+1/2dRgdr*\[Delta]rfunFECritplungeParRes[r,a,p,e,x]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity*)


\[Delta]vtfunFECritplungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];

	\[CapitalDelta]=a^2-2r+r^2;
	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+r^2)+a*EEg(a^2+3r^2))/(r*\[CapitalDelta])+(r(2a^2+a^2*r+r^3)\[Delta]EE)/\[CapitalDelta]-(2a*r*\[Delta]Lz)/\[CapitalDelta];
	ddtgd\[Lambda]drfun=(EEg*r(a^2+3r^2))/\[CapitalDelta]-((r^2-a^2)(-2a*Lzg+EEg*r^3+a^2*EEg(2+r)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunFECritplungePar[r,a,p,e,x]+d\[Delta]td\[Lambda]fun
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


\[Delta]v\[Phi]funFECritplungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];

	\[CapitalDelta]=a^2-2r+r^2;
	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+r)))/(r*\[CapitalDelta])+(2a*r*\[Delta]EE)/\[CapitalDelta]+((r-2)r*\[Delta]Lz)/\[CapitalDelta];
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(r-1)-EEg*r^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunFECritplungePar[r,a,p,e,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Subsection::Closed:: *)
(*Generic plunge*)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunFEPlungePar[r_,a_,p_,e_,x_]:=Module[{EEg,\[Delta]p,\[Delta]r41,r1g,\[Rho]rg,\[Rho]ig,Ar1g,B0,\[Delta]r1,\[Delta]\[Rho]r,\[Phi],krg,ellF,ellE},
	EEg=EEgfuncc[p,e];
	\[Delta]p=\[Delta]pfunFEcc[a,p,e,x];

	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];

	Ar1g=Sqrt[\[Rho]ig^2+(r1g-\[Rho]rg)^2];
	B0=Sqrt[\[Rho]ig^2+\[Rho]rg^2];
	krg=(-(Ar1g-B0)^2+r1g^2)/(4Ar1g*B0);
	\[Delta]r1=\[Delta]p/(1-e);
	\[Delta]\[Rho]r=\[Delta]p/(1+e);
	\[Delta]r41=\[Delta]r41funcc[a,p,e,x];

	\[Phi]=ArcCos[(-Ar1g*r+B0(r1g-r))/(Ar1g*r+B0(r1g-r))];
	ellE=EllipticE[\[Pi]-\[Phi],krg];
	
	Sqrt[(r1g-r)(\[Rho]ig^2+(r-\[Rho]rg)^2)r](1/Sqrt[Ar1g*B0])(\[Delta]r1/2((2B0*ellE)/(Ar1g*r1g)-((Ar1g*r+B0(r1g-r))(Ar1g*B0+r(r+r1g-2\[Rho]rg))Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(Ar1g*r1g*r(\[Rho]ig^2+(r-\[Rho]rg)^2)))+\[Delta]\[Rho]r((r1g-2\[Rho]rg)/(Ar1g*B0)ellE+(r1g(-2Ar1g*B0*r+(\[Rho]ig^2+(r-\[Rho]rg)^2)(r1g-2\[Rho]rg)+2Ar1g*B0*\[Rho]rg)Sin[\[Phi]])/(2Ar1g*B0(B0^2(-2r+r1g)-r^2(r1g-2\[Rho]rg))Sqrt[1-krg*Sin[\[Phi]]^2]))+\[Delta]r41(-((2Ar1g*ellE)/(B0*r1g))+((Ar1g*r+B0(-r+r1g))Sin[\[Phi]]Sqrt[1-krg Sin[\[Phi]]^2])/(B0*r*r1g))-a/2((2Sqrt[Ar1g*B0]Sqrt[r1g-r])/(3r1g*r Sqrt[r(\[Rho]ig^2+(r-\[Rho]rg)^2)])-(4Ar1g(B0^2+2r1g*\[Rho]rg)ellE)/(3B0^3*r1g^2)+(2(Ar1g*r+B0(-r+r1g))(B0^2+r1g*\[Rho]rg)Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3B0^3*r1g^2*r)-Sin[\[Phi]]/(3B0^3 Sqrt[1-krg*Sin[\[Phi]]^2]) ((-\[Rho]ig^2+\[Rho]rg^2)/B0+(4Ar1g*krg*\[Rho]rg*Cos[\[Phi]])/r1g)))+\[Delta]r1
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunFEPlungeParRes[r_,a_,p_,e_,x_]:=Module[{EEg,\[Delta]p,\[Delta]r41,r1g,\[Rho]rg,\[Rho]ig,Ar1g,B0,\[Delta]r1,\[Delta]\[Rho]r,\[Phi],krg,ellF,ellE},
	EEg=EEgfuncc[p,e];
	\[Delta]p=\[Delta]pfunFEcc[a,p,e,x];

	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];

	Ar1g=Sqrt[\[Rho]ig^2+(r1g-\[Rho]rg)^2];
	B0=Sqrt[\[Rho]ig^2+\[Rho]rg^2];
	krg=(-(Ar1g-B0)^2+r1g^2)/(4Ar1g*B0);
	\[Delta]r1=\[Delta]p/(1-e);
	\[Delta]\[Rho]r=\[Delta]p/(1+e);
	\[Delta]r41=\[Delta]r41funcc[a,p,e,x];

	\[Phi]=ArcCos[(-Ar1g*r+B0(r1g-r))/(Ar1g*r+B0(r1g-r))];
	ellE=EllipticE[\[Pi]-\[Phi],krg];

	(1/Sqrt[Ar1g*B0])(\[Delta]r1/2((2B0*ellE)/(Ar1g*r1g)-((Ar1g*r+B0(r1g-r))(Ar1g*B0+r(r+r1g-2\[Rho]rg))Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(Ar1g*r1g*r(\[Rho]ig^2+(r-\[Rho]rg)^2)))+\[Delta]\[Rho]r((r1g-2\[Rho]rg)/(Ar1g*B0)ellE+(r1g(-2Ar1g*B0*r+(\[Rho]ig^2+(r-\[Rho]rg)^2)(r1g-2\[Rho]rg)+2Ar1g*B0*\[Rho]rg)Sin[\[Phi]])/(2Ar1g*B0(B0^2(-2r+r1g)-r^2(r1g-2\[Rho]rg))Sqrt[1-krg*Sin[\[Phi]]^2]))+\[Delta]r41(-((2Ar1g*ellE)/(B0*r1g))+((Ar1g*r+B0(-r+r1g))Sin[\[Phi]]Sqrt[1-krg Sin[\[Phi]]^2])/(B0*r*r1g))-a/2((2Sqrt[Ar1g*B0]Sqrt[r1g-r])/(3r1g*r Sqrt[r(\[Rho]ig^2+(r-\[Rho]rg)^2)])-(4Ar1g(B0^2+2r1g*\[Rho]rg)ellE)/(3B0^3*r1g^2)+(2(Ar1g*r+B0(-r+r1g))(B0^2+r1g*\[Rho]rg)Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3B0^3*r1g^2*r)-Sin[\[Phi]]/(3B0^3 Sqrt[1-krg*Sin[\[Phi]]^2]) ((-\[Rho]ig^2+\[Rho]rg^2)/B0+(4Ar1g*krg*\[Rho]rg*Cos[\[Phi]])/r1g)))+\[Delta]r1/Sqrt[(r1g-r)(\[Rho]ig^2+(r-\[Rho]rg)^2)r]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory*)


\[Delta]tfunFEPlungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]EE,\[Delta]Lz,r1g,r2g,r2gCC,\[Rho]rg,\[Rho]ig,\[Delta]r1,\[Delta]\[Rho]r,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,dtgd\[Lambda]fun,dtgd\[Lambda]funr1g,dtgd\[Lambda]funr2g,dtgd\[Lambda]funr2gCC,Ar1g,B0,krg,\[Gamma]r,\[Phi],ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over,\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r2gCC},
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	\[Delta]p=\[Delta]pfunFEcc[a,p,e,x];
	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x]+dEEdpfuncc[p,e]\[Delta]p;
	\[Delta]Lz=\[Delta]Lzfuncc[a,p,e,x]+dLzdpfuncc[a,p,e,x]\[Delta]p;

	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];
	r2g=\[Rho]rg+I \[Rho]ig;
	r2gCC=\[Rho]rg-I \[Rho]ig;

	Ar1g=Sqrt[\[Rho]ig^2+(r1g-\[Rho]rg)^2];
	B0=Sqrt[\[Rho]ig^2+\[Rho]rg^2];
	krg=(-(Ar1g-B0)^2+r1g^2)/(4Ar1g*B0);
	\[Delta]r1=\[Delta]p/(1-e);
	\[Delta]\[Rho]r=\[Delta]p/(1+e);
	\[Delta]\[Rho]r4=\[Delta]r41funcc[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];

	\[Gamma]r=-((Ar1g-B0)^2/(4Ar1g*B0));
	\[Phi]=ArcCos[(-Ar1g*r+B0(r1g-r))/(Ar1g*r+B0(r1g-r))];

	dtgd\[Lambda]fun=(r(-2a*Lzg+EEg(2a^2+a^2*r+r^3)))/(a^2-2r+r^2);
	dtgd\[Lambda]funr1g=(r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(a^2-2r1g+r1g^2);
	dtgd\[Lambda]funr2g=(r2g(-2a*Lzg+EEg*r2g^3+a^2*EEg(r2g+2)))/(a^2+(r2g-2)r2g);
	dtgd\[Lambda]funr2gCC=(r2gCC(-2a*Lzg+EEg*r2gCC^3+a^2 EEg(r2gCC+2)))/(a^2+(r2gCC-2)r2gCC);
	ellF=EllipticF[\[Pi]-\[Phi],krg];
	ellE=EllipticE[\[Pi]-\[Phi],krg];
	ellPi=EllipticPi[\[Gamma]r,\[Pi]-\[Phi],krg];

	\[ScriptCapitalI]=ellF/(Sqrt[Ar1g*B0]);
	\[ScriptCapitalI]r1g=(2B0*ellE)/(Sqrt[Ar1g*B0]Ar1g*r1g)-(1/r1g+B0/(Ar1g*r1g))\[ScriptCapitalI]-((Ar1g*r+B0(r1g-r))(Ar1g*B0+r(r+r1g-2\[Rho]rg))Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(Sqrt[Ar1g*B0]*Ar1g*r1g*r(\[Rho]ig^2+(r-\[Rho]rg)^2));
	\[ScriptCapitalI]rover=ArcTan[(r1g*Sin[\[Phi]])/(2Sqrt[Ar1g*B0] Sqrt[1-krg*Sin[\[Phi]]^2])]-(B0*r1g)/(Ar1g-B0)*\[ScriptCapitalI]+((Ar1g+B0)r1g*ellPi)/(2(Ar1g-B0)Sqrt[Ar1g*B0]);
	\[ScriptCapitalI]r2over=((Ar1g-B0)Sqrt[r(r1g-r)(\[Rho]ig^2+(r-\[Rho]rg)^2)])/(Ar1g*r+B0(r1g-r))+Sqrt[Ar1g*B0]*ellE+1/2(r1g+2\[Rho]rg)\[ScriptCapitalI]rover-(B0*r1g(r1g-2\[Rho]rg))/(2(Ar1g-B0))*\[ScriptCapitalI];
	\[ScriptCapitalI]r2gCC=(2Sqrt[Ar1g*B0]r1g*ellE)/(B0^2(r1g-r2gCC)^2-Ar1g^2*r2gCC^2)-((Ar1g-B0)ellF)/(Sqrt[Ar1g*B0](B0(r1g-r2gCC)+Ar1g*r2gCC))+(r1g(r-r2g)Sin[\[Phi]])/(Sqrt[Ar1g*B0](-Ar1g*B0(r-r2gCC)+(r1g-r2gCC)r2gCC(r-r2g))Sqrt[1-krg*Sin[\[Phi]]^2]);
	\[ScriptCapitalI]r2g=(2Sqrt[Ar1g*B0]r1g*ellE)/(B0^2(r1g-r2g)^2-Ar1g^2*r2g^2)-((Ar1g-B0)ellF)/(Sqrt[Ar1g*B0](B0(r1g-r2g)+Ar1g*r2g))+(r1g(r-r2gCC)Sin[\[Phi]])/(Sqrt[Ar1g*B0](-Ar1g*B0(r-r2g)+(r1g-r2g)r2g(r-r2gCC))Sqrt[1-krg*Sin[\[Phi]]^2]);

	1/Sqrt[1-EEg^2](((4\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+1/2dtgd\[Lambda]funr1g*\[Delta]r1*\[ScriptCapitalI]r1g+1/2dtgd\[Lambda]funr2gCC*\[ScriptCapitalI]r2gCC*\[Delta]\[Rho]r+1/2dtgd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]\[Rho]r+1/2EEg((2+r1g)\[Delta]r1+2(2+\[Rho]rg)\[Delta]\[Rho]r-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+\[ScriptCapitalI]rover((2\[Delta]EE)/(1-EEg^2)+(EEg*\[Delta]r1)/2+EEg*\[Delta]\[Rho]r+EEg*\[Delta]\[Rho]r4)+dtgd\[Lambda]funr1g*\[Delta]r1(1/Sqrt[r(r1g-r)(r-I*\[Rho]ig-\[Rho]rg)(r+I*\[Rho]ig-\[Rho]rg)])-dtgd\[Lambda]fun*\[Delta]rfunFEPlungeParRes[r,a,p,e,x])
]	


(* ::Subsubsection::Closed:: *)
(*Azimuthal trajectory*)


\[Delta]\[Phi]funFEPlungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]EE,\[Delta]Lz,r1g,r2g,r2gCC,\[Rho]rg,\[Rho]ig,\[Delta]r1,\[Delta]\[Rho]r,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,d\[Phi]gd\[Lambda]funr2g,d\[Phi]gd\[Lambda]funr2gCC,Ar1g,B0,krg,\[Gamma]r,\[Phi],ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r2gCC},
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	\[Delta]p=\[Delta]pfunFEcc[a,p,e,x];
	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x]+dEEdpfuncc[p,e]\[Delta]p;
	\[Delta]Lz=\[Delta]Lzfuncc[a,p,e,x]+dLzdpfuncc[a,p,e,x]\[Delta]p;

	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];
	r2g=\[Rho]rg+I \[Rho]ig;
	r2gCC=\[Rho]rg-I \[Rho]ig;

	Ar1g=Sqrt[\[Rho]ig^2+(r1g-\[Rho]rg)^2];
	B0=Sqrt[\[Rho]ig^2+\[Rho]rg^2];
	krg=(-(Ar1g-B0)^2+r1g^2)/(4Ar1g*B0);
	\[Delta]r1=\[Delta]p/(1-e);
	\[Delta]\[Rho]r=\[Delta]p/(1+e);

	\[Gamma]r=-((Ar1g-B0)^2/(4Ar1g*B0));
	\[Phi]=ArcCos[(-Ar1g*r+B0(r1g-r))/(Ar1g*r+B0(r1g-r))];

	d\[Phi]gd\[Lambda]fun=((2a*EEg+Lzg(-2+r))r)/(a^2-2r+r^2);
	d\[Phi]gd\[Lambda]funr1g=((2a*EEg+Lzg(-2+r1g))r1g)/(a^2-2r1g+r1g^2);
	d\[Phi]gd\[Lambda]funr2g=((2a*EEg+Lzg(-2+r2g))r2g)/(a^2-2r2g+r2g^2);
	d\[Phi]gd\[Lambda]funr2gCC=((2a*EEg+Lzg(-2+r2gCC))r2gCC)/(a^2-2r2gCC+r2gCC^2);
	ellF=EllipticF[\[Pi]-\[Phi],krg];
	ellE=EllipticE[\[Pi]-\[Phi],krg];
	ellPi=EllipticPi[\[Gamma]r,\[Pi]-\[Phi],krg];

	\[ScriptCapitalI]=ellF/(Sqrt[Ar1g*B0]);
	\[ScriptCapitalI]r1g=(2B0*ellE)/(Sqrt[Ar1g*B0]Ar1g*r1g)-(1/r1g+B0/(Ar1g*r1g))\[ScriptCapitalI]-((Ar1g*r+B0(r1g-r))(Ar1g*B0+r(r+r1g-2\[Rho]rg))Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(Sqrt[Ar1g*B0]*Ar1g*r1g*r(\[Rho]ig^2+(r-\[Rho]rg)^2));
	\[ScriptCapitalI]r2gCC=(2Sqrt[Ar1g*B0]r1g*ellE)/(B0^2(r1g-r2gCC)^2-Ar1g^2*r2gCC^2)-((Ar1g-B0)ellF)/(Sqrt[Ar1g*B0](B0(r1g-r2gCC)+Ar1g*r2gCC))+(r1g(r-r2g)Sin[\[Phi]])/(Sqrt[Ar1g*B0](-Ar1g*B0(r-r2gCC)+(r1g-r2gCC)r2gCC(r-r2g))Sqrt[1-krg*Sin[\[Phi]]^2]);
	\[ScriptCapitalI]r2g=(2Sqrt[Ar1g*B0]r1g*ellE)/(B0^2(r1g-r2g)^2-Ar1g^2*r2g^2)-((Ar1g-B0)ellF)/(Sqrt[Ar1g*B0](B0(r1g-r2g)+Ar1g*r2g))+(r1g(r-r2gCC)Sin[\[Phi]])/(Sqrt[Ar1g*B0](-Ar1g*B0(r-r2g)+(r1g-r2g)r2g(r-r2gCC))Sqrt[1-krg*Sin[\[Phi]]^2]);
	
	1/Sqrt[1-EEg^2]((-EEg+(EEg*Lzg*\[Delta]EE)/(1-EEg^2)+\[Delta]Lz)\[ScriptCapitalI]+1/2d\[Phi]gd\[Lambda]funr1g*\[Delta]r1*\[ScriptCapitalI]r1g+1/2\[ScriptCapitalI]r2g*d\[Phi]gd\[Lambda]funr2g*\[Delta]\[Rho]r+1/2\[ScriptCapitalI]r2gCC*d\[Phi]gd\[Lambda]funr2gCC*\[Delta]\[Rho]r+d\[Phi]gd\[Lambda]funr1g*\[Delta]r1(1/(Sqrt[r(r1g-r)(r-I*\[Rho]ig-\[Rho]rg)(r+I*\[Rho]ig-\[Rho]rg)]))-d\[Phi]gd\[Lambda]fun*\[Delta]rfunFEPlungeParRes[r,a,p,e,x])
]


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunFEPlungePar[r_,a_,p_,e_,x_]:=Module[{EEg,\[Delta]p,\[Delta]EE,\[Delta]r41,r1g,\[Rho]rg,\[Rho]ig,Ar1g,B0,\[Delta]r1,\[Delta]\[Rho]r,dRgdr},
	EEg=EEgfuncc[p,e];
	\[Delta]p=\[Delta]pfunFEcc[a,p,e,x];
	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x]+dEEdpfuncc[p,e]\[Delta]p;

	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];
	\[Delta]r1=\[Delta]p/(1-e);
	\[Delta]\[Rho]r=\[Delta]p/(1+e);
	\[Delta]r41=\[Delta]r41funcc[a,p,e,x];
	dRgdr=Sqrt[1-EEg^2](-4r^3+3r^2(r1g+2\[Rho]rg)+r1g(\[Rho]ig^2+\[Rho]rg^2)-2r(\[Rho]ig^2+\[Rho]rg(2 r1g+\[Rho]rg)));

	-Sqrt[(1-EEg^2)(r1g-r)(\[Rho]ig^2+(r-\[Rho]rg)^2)r](\[Delta]r1/(2(r-r1g))+(\[Delta]\[Rho]r(r-\[Rho]rg))/(r^2+\[Rho]ig^2-2r*\[Rho]rg+\[Rho]rg^2)+\[Delta]r41/r-a/(2r^2)+(EEg*\[Delta]EE)/(1-EEg^2))+1/2dRgdr*\[Delta]rfunFEPlungeParRes[r,a,p,e,x]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity*)


\[Delta]vtfunFEPlungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	\[Delta]p=\[Delta]pfunFEcc[a,p,e,x];
	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x]+dEEdpfuncc[p,e]\[Delta]p;
	\[Delta]Lz=\[Delta]Lzfuncc[a,p,e,x]+dLzdpfuncc[a,p,e,x]\[Delta]p;
	
	\[CapitalDelta]=a^2-2r+r^2;
	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+r^2)+a*EEg(a^2+3r^2))/(r*\[CapitalDelta])+(r(2a^2+a^2*r+r^3)\[Delta]EE)/\[CapitalDelta]-(2a*r*\[Delta]Lz)/\[CapitalDelta];
	ddtgd\[Lambda]drfun=(EEg*r(a^2+3r^2))/\[CapitalDelta]-((r^2-a^2)(-2a*Lzg+EEg*r^3+a^2*EEg(2+r)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunFEPlungePar[r,a,p,e,x]+d\[Delta]td\[Lambda]fun
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


\[Delta]v\[Phi]funFEPlungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	\[Delta]p=\[Delta]pfunFEcc[a,p,e,x];
	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x]+dEEdpfuncc[p,e]\[Delta]p;
	\[Delta]Lz=\[Delta]Lzfuncc[a,p,e,x]+dLzdpfuncc[a,p,e,x]\[Delta]p;

	\[CapitalDelta]=a^2-2r+r^2;
	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+r)))/(r*\[CapitalDelta])+(2a*r*\[Delta]EE)/\[CapitalDelta]+((r-2)r*\[Delta]Lz)/\[CapitalDelta];
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(r-1)-EEg*r^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunFEPlungePar[r,a,p,e,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Subsection::Closed:: *)
(*Plunge related to bound orbits*)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunFEPlungeBoundOrbitPar[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,EEg,\[Delta]EE,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],krg,ellF,ellE},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunFE[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Phi]=ArcSin[Sqrt[(r2g (r3g-r))/(r3g (r2g-r))]];
	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	
	Sqrt[(r1g-r)(r-r2g)r(r-r3g)](\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg(2ellF)/Sqrt[(r1g-r3g)r2g]+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Sqrt[(r1g-r3g)r2g]+1/2 (2/Sqrt[r2g(r1g-r3g)] 1/(r1g-r2g)((r2g ellE)/r1g-ellF-((r1g-r2g)r3g*Cos[\[Phi]]*Sin[\[Phi]])/(r1g(r1g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2])))\[Delta]r1+1/2(2/((r1g-r2g)Sqrt[r2g(r1g-r3g)]) (-(((r1g-r3g)ellE)/(r2g-r3g))+ellF))\[Delta]r2+1/2 ((2r2g*ellE)/(Sqrt[r2g(r1g-r3g)](r2g-r3g)r3g)-(2 ellF)/(Sqrt[r2g(r1g-r3g)]r3g))\[Delta]r3+(2/(Sqrt[r2g (r1g-r3g)]r3g)(ellF+(-1+r3g/r1g)(ellE-Sqrt[1-krg*Sin[\[Phi]]^2]Tan[\[Phi]])))\[Delta]r41-a/2((2((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g Sqrt[r2g(r1g-r3g)]r3g^2) ellF-2/(3Sqrt[r2g(r1g-r3g)]) (r2g*r3g+r1g(r2g+r3g))/(r1g^2r2g^2r3g^2) r2g(r1g-r3g)(2ellE-Tan[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])+(2r3g)/(3r1g^2(r1g-r3g)Sqrt[r2g(r2g-r3g)]) Sqrt[(r2g-r)/(r1g-r)]Sin[\[Phi]]Cos[\[Phi]]-(2 (r3g^2+r3g*r+r^2-r1g(r3g+r)))/(3(r1g-r3g)r3g^2*r) Sqrt[(r3g-r)/(r(r1g-r)(r2g-r))]))+1/2 ((2(r1g-r)r)/((r1g-r3g)r3g)\[Delta]r3)
]


\[Delta]rfunFEPlungeBoundOrbitDoubleRootPar[r_,a_,p_,x_]:=Module[{r1g,r3g,\[Delta]r1,\[Delta]r41,\[Delta]EE,EEg,\[Phi]},
	EEg=EEgfun[a,p,0,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,0,x]+dEEdpfun[a,p,0,x]\[Delta]pfunFE[a,p,0,x];

	r1g=p;
	r3g=2/(1-EEg^2)-2p;
	\[Delta]r1=\[Delta]r1funFE[a,p,0,x];
	\[Delta]r41=\[Delta]r41fun[a,p,0,x];

	\[Phi]=ArcSin[Sqrt[(r1g (r3g-r))/(r3g (r1g-r))]];
	
	(r1g-r)Sqrt[(r3g-r)r]((-((r3g*Cos[\[Phi]]*Sin[\[Phi]])/(r1g(r1g-r3g))^(3/2)))\[Delta]r1+((2(r1g-r3g)Tan[\[Phi]])/(r1g*r3g Sqrt[r1g(r1g-r3g)]))\[Delta]r41-a/2(-((2(r^2+r*r3g+r3g^2-r1g(r+r3g)))/(3(r(r1g-r3g)r3g^2)(r1g-r))) Sqrt[r3g/r-1]+2/3 1/r1g^2 1/Sqrt[r1g(r1g-r3g)] ((r3g*Cos[\[Phi]]*Sin[\[Phi]])/ (r1g-r3g)+((r1g-r3g)(r1g+2r3g)Tan[\[Phi]])/ r3g^2)))+1/2 ((2(r1g-r)r)/((r1g-r3g)r3g) \[Delta]r1)
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunFEPlungeBoundOrbitParRes[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]KK,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],krg,\[Gamma]rg,ellF,ellE},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunFE[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Phi]=ArcSin[Sqrt[(r2g(r3g-r))/(r3g (r2g-r))]];
	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	
	\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg(2ellF)/Sqrt[(r1g-r3g)r2g]+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Sqrt[(r1g-r3g)r2g]+1/2(2/Sqrt[r2g(r1g-r3g)] 1/(r1g-r2g)((r2g ellE)/r1g-ellF-((r1g-r2g)r3g*Cos[\[Phi]]*Sin[\[Phi]])/(r1g(r1g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2])))\[Delta]r1+1/2(2/((r1g-r2g)Sqrt[r2g(r1g-r3g)])(-(((r1g-r3g)ellE)/(r2g-r3g))+ellF))\[Delta]r2+1/2((2r2g*ellE)/(Sqrt[r2g(r1g-r3g)](r2g-r3g)r3g)-(2 ellF)/(Sqrt[r2g(r1g-r3g)]r3g))\[Delta]r3+(2/(Sqrt[r2g(r1g-r3g)]r3g) (ellF+(-1+r3g/r1g)(ellE-Sqrt[1-krg*Sin[\[Phi]]^2]Tan[\[Phi]])))\[Delta]r41-a/2(((2((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g Sqrt[r2g(r1g-r3g)]r3g^2) ellF-2/(3Sqrt[r2g(r1g-r3g)]) (r2g*r3g+r1g(r2g+r3g))/(r1g^2r2g^2r3g^2) r2g(r1g-r3g)(2ellE-Tan[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])+(2r3g)/(3r1g^2(r1g-r3g)Sqrt[r2g(r2g-r3g)]) Sqrt[(r2g-r)/(r1g-r)]Sin[\[Phi]]Cos[\[Phi]]-(2 (r3g^2+r3g*r+r^2-r1g(r3g+r)))/(3(r1g-r3g)r3g^2*r) Sqrt[(r3g-r)/(r(r1g-r)(r2g-r))]))+1/2(2Sqrt[(r1g-r)r])/((r1g-r3g)r3g Sqrt[(r2g-r)(r3g-r)])\[Delta]r3
]


\[Delta]rfunFEPlungeBoundOrbitDoubleRootParRes[r_,a_,p_,x_]:=Module[{r1g,r3g,\[Delta]r1,\[Delta]r41,\[Delta]EE,EEg,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi]},
	EEg=EEgfun[a,p,0,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,0,x]+dEEdpfun[a,p,0,x]\[Delta]pfunFE[a,p,0,x];

	r1g=p;
	r3g=2/(1-EEg^2)-2p;
	\[Delta]r1=\[Delta]r1funFE[a,p,0,x];
	\[Delta]r41=\[Delta]r41fun[a,p,0,x];

	\[Phi]=ArcSin[Sqrt[(r1g (r3g-r))/(r3g (r1g-r))]];
	
	(-((r3g*Cos[\[Phi]]*Sin[\[Phi]])/(r1g(r1g-r3g))^(3/2)))\[Delta]r1+((2(r1g-r3g)Tan[\[Phi]])/(r1g*r3g Sqrt[r1g(r1g-r3g)]))\[Delta]r41-a/2(-((2(r^2+r*r3g+r3g^2-r1g(r+r3g)))/(3(r(r1g-r3g)r3g^2)(r1g-r))) Sqrt[r3g/r-1]+2/3 1/r1g^2 1/Sqrt[r1g(r1g-r3g)] ((r3g*Cos[\[Phi]]*Sin[\[Phi]])/ (r1g-r3g)+((r1g-r3g)(r1g+2r3g)Tan[\[Phi]])/ r3g^2))+1/2((2Sqrt[(r1g-r)r])/((r1g-r3g)*r3g Sqrt[(r1g-r)(r3g-r)])\[Delta]r1)
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory*)


\[Delta]tfunFEPlungeBoundOrbitPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Phi],\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,dtgd\[Lambda]fun,dtgd\[Lambda]funr1g,dtgd\[Lambda]funr2g,dtgd\[Lambda]funr3g,krg,\[Gamma]r,ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];

	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));
	\[Gamma]r=r3g/r2g;
	\[Phi]=ArcSin[\[Sqrt]((r2g(r3g-r))/(r3g(r2g-r)))];

	dtgd\[Lambda]fun=(r(-2a*Lzg+EEg(2a^2+a^2*r+r^3)))/(a^2-2r+r^2);
	dtgd\[Lambda]funr1g=(r1g(-2a*Lzg+EEg(2a^2+a^2r1g+r1g^3)))/(a^2-2r1g+r1g^2);
	dtgd\[Lambda]funr2g=(r2g(-2a*Lzg+EEg(2a^2+a^2r2g+r2g^3)))/(a^2-2r2g+r2g^2);
	dtgd\[Lambda]funr3g=(r3g(-2a*Lzg+EEg(2a^2+a^2r3g+r3g^3)))/(a^2-2r3g+r3g^2);

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	ellPi=EllipticPi[\[Gamma]r,\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[(r2g(r1g-r3g))];
	\[ScriptCapitalI]r1g=(2/Sqrt[(r2g(r1g-r3g))])(1/(r1g-r2g))((r2g*ellE)/r1g-ellF-((r1g-r2g)r3g*Cos[\[Phi]]*Sin[\[Phi]])/(r1g(r1g-r3g)Sqrt[(1-krg*Sin[\[Phi]]^2)]));
	\[ScriptCapitalI]r2g=(2/Sqrt[(r2g(r1g-r3g))])(1/(r1g-r2g))(-(((r1g-r3g)ellE)/(r2g-r3g))+ellF);
	\[ScriptCapitalI]r3g=(2/Sqrt[(r2g(r1g-r3g))]) 1/r3g ((r2g*ellE)/(r2g-r3g)-ellF);
	\[ScriptCapitalI]rover=2 /Sqrt[r2g(r1g-r3g)] (r2g*ellF+(-r2g+r3g)ellPi);
	\[ScriptCapitalI]r2over=-r3g Sqrt[((r1g-r)(r2g-r))/(r2g(r2g-r3g))]Sin[\[Phi]]Cos[\[Phi]]+Sqrt[r2g(r1g-r3g)]ellE-((r2g-r3g)(r1g+r2g+r3g))/Sqrt[r2g(r1g-r3g)] ellPi+(r2g(r2g+r3g))/Sqrt[r2g(r1g-r3g)] ellF;

	-1/Sqrt[1-EEg^2]((4\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+1/2 EEg(\[Delta]r1+\[Delta]r2+\[Delta]r3+2\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+1/2EEg((2+r1g)\[Delta]r1+(2+r2g)\[Delta]r2+(2+r3g)\[Delta]r3-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+1/2dtgd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+1/2dtgd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2+1/2dtgd\[Lambda]funr3g*\[ScriptCapitalI]r3g*\[Delta]r3+1/2 dtgd\[Lambda]funr3g((2/(r1g-r3g))1/r3g((\[Sqrt]((r1g-r)r))/(\[Sqrt]((r2g-r)(r3g-r)))))\[Delta]r3-dtgd\[Lambda]fun*\[Delta]rfunFEPlungeBoundOrbitParRes[r,a,p,e,x])
]


\[Delta]tfunFEPlungeBoundOrbitDoubleRootPar[r_,a_,p_,x_]:=Module[{EEg,Lzg,\[Delta]EE,r1g,r3g,\[Delta]r1,\[Phi],\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,dtgd\[Lambda]fun,dtgd\[Lambda]funr1g,dtgd\[Lambda]funr2g,dtgd\[Lambda]funr3g,\[ScriptCapitalI],\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	EEg=EEgfun[a,p,0,x];
	Lzg=Lzgfun[a,p,0,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,0,x]+dEEdpfun[a,p,0,x]\[Delta]pfunFE[a,p,0,x];

	r1g=p;
	r3g=2/(1-EEg^2)-2p;
	\[Delta]r1=\[Delta]r1funFE[a,p,0,x];
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,0,x];
	\[Delta]\[Rho]i4=Sqrt[a];

	\[Phi]=ArcSin[Sqrt[(r1g(r3g-r))/(r3g(r1g-r))]];

	dtgd\[Lambda]fun=(r(-2a*Lzg+EEg(2a^2+a^2*r+r^3)))/(a^2-2r+r^2);
	dtgd\[Lambda]funr1g=(r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(a^2-2r1g+r1g^2);
	dtgd\[Lambda]funr3g=(r3g(-2a*Lzg+EEg(2a^2+a^2*r3g+r3g^3)))/(a^2-2r3g+r3g^2);
	\[ScriptCapitalI]=(2\[Phi])/Sqrt[r1g(r1g-r3g)];
	\[ScriptCapitalI]r3g=2 /(Sqrt[r1g(r1g-r3g)]r3g) ((r3g*\[Phi])/(r1g-r3g));
	\[ScriptCapitalI]rover=(2r1g*\[Phi])/Sqrt[r1g(r1g-r3g)]-2 1/2 I*Log[(Sqrt[r]-I*Sqrt[r3g-r])/(Sqrt[r]+I*Sqrt[r3g-r])];
	\[ScriptCapitalI]r2over=-Sqrt[(r3g-r)r]+(2r1g^2*\[Phi])/Sqrt[r1g(r1g-r3g)]-(r1g(2r1g+r3g))/Sqrt[r1g(r1g-r3g)] 1/2 I*Log[(Sqrt[r]-I*Sqrt[r3g-r])/(Sqrt[r]+I*Sqrt[r3g-r])]Sqrt[1-r3g/r1g];

	-1/Sqrt[1-EEg^2]((4\[Delta]EE)/(1-EEg^2) \[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+1/2 EEg(3\[Delta]r1+2\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+\[Delta]EE/(1-EEg^2) \[ScriptCapitalI]r2over+1/2 EEg(2(2+r1g)\[Delta]r1+(2+r3g)\[Delta]r1-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+(dtgd\[Lambda]funr1g(-2r1g+r3g)\[Phi])/(r1g(r1g-r3g))^(3/2) \[Delta]r1+1/2 dtgd\[Lambda]funr3g*\[ScriptCapitalI]r3g*\[Delta]r1+1/2 dtgd\[Lambda]funr3g((2Sqrt[r])/((r1g-r3g)r3g Sqrt[-r+r3g]))\[Delta]r1-dtgd\[Lambda]fun*\[Delta]rfunFEPlungeBoundOrbitDoubleRootParRes[r,a,p,x])
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal trajectory*)


\[Delta]\[Phi]funFEPlungeBoundOrbitPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Phi],d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,d\[Phi]gd\[Lambda]funr2g,d\[Phi]gd\[Lambda]funr3g,krg,\[Gamma]r,ellF,ellE,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;

	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));
	\[Phi]=ArcSin[\[Sqrt]((r2g(r3g-r))/(r3g(r2g-r)))];

	d\[Phi]gd\[Lambda]fun=((2a*EEg+Lzg(-2+r))r)/(a^2-2r+r^2);
	d\[Phi]gd\[Lambda]funr1g=((2a*EEg+Lzg(-2+r1g))r1g)/(a^2-2r1g+r1g^2);
	d\[Phi]gd\[Lambda]funr2g=((2a*EEg+Lzg(-2+r2g))r2g)/(a^2-2r2g+r2g^2);
	d\[Phi]gd\[Lambda]funr3g=((2a*EEg+Lzg(-2+r3g))r3g)/(a^2-2r3g+r3g^2);

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[(r2g(r1g-r3g))];
	\[ScriptCapitalI]r1g=(2/Sqrt[(r2g(r1g-r3g))])(1/(r1g-r2g))((r2g*ellE)/r1g-ellF-((r1g-r2g)r3g*Cos[\[Phi]]*Sin[\[Phi]])/(r1g(r1g-r3g)Sqrt[(1-krg*Sin[\[Phi]]^2)]));
	\[ScriptCapitalI]r2g=(2/Sqrt[(r2g(r1g-r3g))])(1/(r1g-r2g))(-(((r1g-r3g)ellE)/(r2g-r3g))+ellF);
	\[ScriptCapitalI]r3g=(2/Sqrt[(r2g(r1g-r3g))]) 1/r3g ((r2g*ellE)/(r2g-r3g)-ellF);

	-1/Sqrt[1-EEg^2](-EEg(1-(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+1/2d\[Phi]gd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+1/2d\[Phi]gd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2+1/2d\[Phi]gd\[Lambda]funr3g*\[ScriptCapitalI]r3g*\[Delta]r3+1/2d\[Phi]gd\[Lambda]funr3g*\[Delta]r3((2/(r1g-r3g))1/r3g((\[Sqrt]((r1g-r)r))/(\[Sqrt]((r2g-r)(r3g-r)))))-d\[Phi]gd\[Lambda]fun*\[Delta]rfunFEPlungeBoundOrbitParRes[r,a,p,e,x])
]


\[Delta]\[Phi]funFEPlungeBoundOrbitDoubleRootPar[r_,a_,p_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r3g,\[Delta]r1,\[Phi],d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,d\[Phi]gd\[Lambda]funr2g,d\[Phi]gd\[Lambda]funr3g,\[ScriptCapitalI],\[ScriptCapitalI]r3g},
	EEg=EEgfun[a,p,0,x];
	Lzg=Lzgfun[a,p,0,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,0,x]+dEEdpfun[a,p,0,x]\[Delta]pfunFE[a,p,0,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,0,x]+dLzdpfun[a,p,0,x]\[Delta]pfunFE[a,p,0,x];

	r1g=p;
	r3g=2/(1-EEg^2)-2p;
	\[Delta]r1=\[Delta]r1funFE[a,p,0,x];

	\[Phi]=ArcSin[Sqrt[(r1g(r3g-r))/(r3g(r1g-r))]];

	d\[Phi]gd\[Lambda]fun=((2a*EEg+Lzg(-2+r))r)/(a^2-2r+r^2);
	d\[Phi]gd\[Lambda]funr1g=((2a*EEg+Lzg(-2+r1g))r1g)/(a^2-2r1g+r1g^2);
	d\[Phi]gd\[Lambda]funr3g=((2 a EEg+Lzg (-2+r3g))r3g)/(a^2-2r3g+r3g^2);
	
	\[ScriptCapitalI]=(2\[Phi])/Sqrt[r1g(r1g-r3g)];
	\[ScriptCapitalI]r3g=2 /(Sqrt[r1g(r1g-r3g)]r3g) ((r3g*\[Phi])/(r1g-r3g));

	-1/Sqrt[1-EEg^2](-EEg(1-(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+(d\[Phi]gd\[Lambda]funr1g(-2r1g+r3g)\[Phi])/(r1g(r1g-r3g))^(3/2)\[Delta]r1+1/2 d\[Phi]gd\[Lambda]funr3g*\[ScriptCapitalI]r3g*\[Delta]r1+1/2 d\[Phi]gd\[Lambda]funr3g*\[Delta]r1((2Sqrt[r])/((r1g-r3g)r3g Sqrt[-r+r3g]))-d\[Phi]gd\[Lambda]fun*\[Delta]rfunFEPlungeBoundOrbitDoubleRootParRes[r,a,p,x])
]


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunFEPlungeBoundOrbitPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,dRgdr},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	dRgdr=Sqrt[(1-EEg^2)](-4r^3+r1g*r2g*r3g+3r^2(r1g+r2g+r3g)-2r(r2g*r3g+r1g(r2g+r3g)));

	Sqrt[(1-EEg^2)r(r1g-r)(r-r2g)(r-r3g)](\[Delta]r1/(2(r-r1g))+\[Delta]r2/(2(r-r2g))+\[Delta]r3/(2(r-r3g))+\[Delta]r41/r-a/(2r^2)+(EEg*\[Delta]EE)/(1-EEg^2))-1/2dRgdr*\[Delta]rfunFEPlungeBoundOrbitParRes[r,a,p,e,x]
]


\[Delta]vrfunFEPlungeBoundOrbitDoubleRootPar[r_,a_,p_,x_]:=Module[{EEg,\[Delta]EE,r1g,r3g,\[Delta]r1,\[Delta]r41,dRgdr},
	EEg=EEgfun[a,p,0,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,0,x]+dEEdpfun[a,p,0,x]\[Delta]pfunFE[a,p,0,x];
	\[Delta]r1=\[Delta]r1funFE[a,p,0,x];
	\[Delta]r41=\[Delta]r41fun[a,p,0,x];

	r1g=p;
	r3g=2/(1-EEg^2)-2p;
	dRgdr=Sqrt[(1-EEg^2)](-4r^3+r1g^2*r3g+3r^2(2r1g+r3g)-2r(r1g*r3g+r1g(r1g+r3g)));

	(r1g-r)Sqrt[(1-EEg^2)r(r3g-r)](\[Delta]r1/(r-r1g)+\[Delta]r1/(2(r-r3g))+\[Delta]r41/r-a/(2 r^2)+(EEg*\[Delta]EE)/(1-EEg^2))-1/2dRgdr*\[Delta]rfunFEPlungeBoundOrbitDoubleRootParRes[r,a,p,x]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity*)


\[Delta]vtfunFEPlungeBoundOrbitPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]td\[Lambda]fun,d\[Delta]td\[Lambda]fun2,ddtgd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];

	\[CapitalDelta]=a^2-2r+r^2;
	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+r^2)+a*EEg(a^2+3r^2))/(r*\[CapitalDelta])+(r(2a^2+a^2*r+r^3)\[Delta]EE)/\[CapitalDelta]-(2a*r*\[Delta]Lz)/\[CapitalDelta];
	ddtgd\[Lambda]drfun=(EEg*r(a^2+3r^2))/\[CapitalDelta]-((r^2-a^2)(-2a*Lzg+EEg*r^3+a^2*EEg(2+r)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunFEPlungeBoundOrbitPar[r,a,p,e,x]+d\[Delta]td\[Lambda]fun
]


\[Delta]vtfunFEPlungeBoundOrbitDoubleRootPar[r_,a_,p_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	EEg=EEgfun[a,p,0,x];
	Lzg=Lzgfun[a,p,0,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,0,x]+dEEdpfun[a,p,0,x]\[Delta]pfunFE[a,p,0,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,0,x]+dLzdpfun[a,p,0,x]\[Delta]pfunFE[a,p,0,x];

	\[CapitalDelta]=a^2-2r+r^2;
	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+r^2)+a*EEg(a^2+3r^2))/(r*\[CapitalDelta])+(r(2a^2+a^2*r+r^3)\[Delta]EE)/\[CapitalDelta]-(2a*r*\[Delta]Lz)/\[CapitalDelta];
	ddtgd\[Lambda]drfun=(EEg*r(a^2+3r^2))/\[CapitalDelta]-((r^2-a^2)(-2a*Lzg+EEg*r^3+a^2*EEg(2+r)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunFEPlungeBoundOrbitDoubleRootPar[r,a,p,x]+d\[Delta]td\[Lambda]fun
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


\[Delta]v\[Phi]funFEPlungeBoundOrbitPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];

	\[CapitalDelta]=a^2-2r+r^2;
	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+r)))/(r*\[CapitalDelta])+(2a*r*\[Delta]EE)/\[CapitalDelta]+((r-2)r*\[Delta]Lz)/\[CapitalDelta];
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(r-1)-EEg*r^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunFEPlungeBoundOrbitPar[r,a,p,e,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


\[Delta]v\[Phi]funFEPlungeBoundOrbitDoubleRootPar[r_,a_,p_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	EEg=EEgfun[a,p,0,x];
	Lzg=Lzgfun[a,p,0,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,0,x]+dEEdpfun[a,p,0,x]\[Delta]pfunFE[a,p,0,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,0,x]+dLzdpfun[a,p,0,x]\[Delta]pfunFE[a,p,0,x];

	\[CapitalDelta]=a^2-2r+r^2;
	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+r)))/(r*\[CapitalDelta])+(2a*r*\[Delta]EE)/\[CapitalDelta]+((r-2)r*\[Delta]Lz)/\[CapitalDelta];
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(r-1)-EEg*r^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunFEPlungeBoundOrbitDoubleRootPar[r,a,p,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Section::Closed:: *)
(*Spin corrections to the orbit - parallel component of the spin, fixed eccentricity, non-linearized*)


(* ::Text:: *)
(*The ISCO plunge is universal to any spin gauge. I implemented here for convenience*)


(* ::Subsubsection::Closed:: *)
(*Coordinate-time velocity - general expression, no expansion*)


VtrgfunPartialExp[r_,a_,S_,p_,e_,x_]:=Module[{EEg,\[Delta]EE,Lzg,\[Delta]Lz,EE,Lz,\[CapitalDelta]},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	
	EE=EEg+S*RealSign[x]*\[Delta]EE;
	Lz=Lzg+S*RealSign[x]*\[Delta]Lz;
	\[CapitalDelta]=a^2-2r+r^2;
	
	EE (r^2+a^2)^2/\[CapitalDelta]-Lz (2a*r)/\[CapitalDelta]-EE*a^2
]


\[Delta]vtfunParPartialExp[r_,a_,S_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,EE,Lz,\[CapitalDelta],d\[Delta]td\[Lambda]fun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	
	EE=EEg+S*RealSign[x]*\[Delta]EE;
	Lz=Lzg+S*RealSign[x]*\[Delta]Lz;
	\[CapitalDelta]=a^2-2r+r^2;
	
	(-Lz(a^2+r^2)+a*EE(a^2+3r^2))/(r*\[CapitalDelta])
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity - general expression, no expansion*)


V\[Phi]rgfunPartialExp[r_,a_,S_,p_,e_,x_]:=Module[{EEg,\[Delta]EE,Lzg,\[Delta]Lz,EE,Lz,\[CapitalDelta]},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	
	EE=EEg+S*RealSign[x]*\[Delta]EE;
	Lz=Lzg+S*RealSign[x]*\[Delta]Lz;
	\[CapitalDelta]=a^2-2r+r^2;
	
	a/\[CapitalDelta](EE(r^2+a^2)-a*Lz)-a*EE+Lz
]


\[Delta]v\[Phi]funParPartialExp[r_,a_,S_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]EE,\[Delta]Lz,EE,Lz,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	
	EE=EEg+S*RealSign[x]*\[Delta]EE;
	Lz=Lzg+S*RealSign[x]*\[Delta]Lz;
	\[CapitalDelta]=a^2-2r+r^2;
	
	-EE+(a(-Lz+a*EE(1+r)))/(r*\[CapitalDelta])
]


(* ::Subsection::Closed:: *)
(*ISCO plunge*)


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory *)


tfunISCOplungePartialExp[r_,a_,S_,x_]:=Module[{r1g,rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]rISCO,\[Delta]r4,EE,Lz,r1,\[CapitalDelta],tgPartialExp,\[Delta]tPartialExp,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]p,\[ScriptCapitalI]mreg,\[ScriptCapitalI]mreg2,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	r1g=ISCOradius[a,x];
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,r1g,0,x];
	Lzg=Lzgfun[a,r1g,0,x];
	\[Delta]EE=\[Delta]EEISCOfun[a,x];
	\[Delta]Lz=\[Delta]LzISCOfun[a,x];
	\[Delta]rISCO=\[Delta]rISCOfun[a,x];
	\[Delta]r4=\[Delta]r41fun[a,r1g,0,x];
	
	EE=EEg+S*RealSign[x]*\[Delta]EE;
	Lz=Lzg+S*RealSign[x]*\[Delta]Lz;
	r1=r1g+S*RealSign[x]*\[Delta]rISCO;
	\[CapitalDelta]=a^2-2r+r^2;

	\[ScriptCapitalI]=2Sqrt[r]/(r1*Sqrt[r1-r]);
	\[ScriptCapitalI]rover=2ArcSin[Sqrt[1-r/r1]]+r1*\[ScriptCapitalI];
	\[ScriptCapitalI]r2over=Sqrt[r(-r+r1)]+r1*ArcSin[Sqrt[1-r/r1]]+r1*\[ScriptCapitalI]rover;
	\[ScriptCapitalI]p=1/(r1-rp)\[ScriptCapitalI]+1/(2(r1-rp)^(3/2)Sqrt[rp])Log[(Sqrt[r(r1-rp)]-Sqrt[(r1-r)rp])^2/(Sqrt[r(r1-rp)]+Sqrt[(r1-r)rp])^2];
	\[ScriptCapitalI]mreg=rm/(r1-rm)\[ScriptCapitalI]+Sqrt[rm]/(2(r1-rm)^(3/2))Log[(Sqrt[r(r1-rm)]-Sqrt[(r1-r)rm])^2/(Sqrt[r(r1-rm)]+Sqrt[(r1-r)rm])^2];
	(*\[ScriptCapitalI]mreg2=Sqrt[rm]/(r1-rm)\[ScriptCapitalI]+1/(2(r1-rm)^(3/2))Log[(Sqrt[r(r1-rm)]-Sqrt[(r1-r)rm])^2/(Sqrt[r(r1-rm)]+Sqrt[(r1-r)rm])^2];*)
	
	tgPartialExp=-1/Sqrt[1-EE^2](4EE*\[ScriptCapitalI]+2EE*\[ScriptCapitalI]rover+EE*\[ScriptCapitalI]r2over-4EE*rp/(rp-rm)(rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)-2(-4EE+a*Lz)1/(rp-rm)(rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg));
	\[Delta]tPartialExp=-1/Sqrt[(1-EE^2)](1/2EE(4\[Delta]r4-a)\[ScriptCapitalI]+EE*\[Delta]r4*\[ScriptCapitalI]rover-\[ScriptCapitalI]p/( 2(rp-rm)) ((a^3*EE+2Lz*rp+a*EE(-6+rp)rp)-2(-2a*Lz+EE*rp^3+a^2*EE(2+rp))\[Delta]r4)+(\[ScriptCapitalI]mreg/rm)/( 2(rp-rm)) ((a^3*EE+2Lz*rm+a*EE(-6+rm)rm)-2(-2a*Lz+EE*rm^3+a^2*EE(2+rm))\[Delta]r4));
	
	tgPartialExp+S*\[Delta]tPartialExp
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal trajectory*)


\[Phi]funISCOplungePartialExp[r_,a_,S_,x_]:=Module[{r1g,rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]rISCO,\[Delta]r4,EE,Lz,r1,\[CapitalDelta],\[Phi]gPartialExp,\[Delta]\[Phi]PartialExp,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]p,\[ScriptCapitalI]mreg,\[ScriptCapitalI]mreg2},
	r1g=ISCOradius[a,x];
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,r1g,0,x];
	Lzg=Lzgfun[a,r1g,0,x];
	\[Delta]EE=\[Delta]EEISCOfun[a,x];
	\[Delta]Lz=\[Delta]LzISCOfun[a,x];
	\[Delta]rISCO=\[Delta]rISCOfun[a,x];
	\[Delta]r4=\[Delta]r41fun[a,r1g,0,x];
	
	EE=EEg+S*RealSign[x]*\[Delta]EE;
	Lz=Lzg+S*RealSign[x]*\[Delta]Lz;
	r1=r1g+S*RealSign[x]*\[Delta]rISCO;
	\[CapitalDelta]=a^2-2r+r^2;
	
	\[ScriptCapitalI]=2Sqrt[r]/(r1*Sqrt[r1-r]);
	\[ScriptCapitalI]p=1/(r1-rp)\[ScriptCapitalI]+1/(2(r1-rp)^(3/2)Sqrt[rp])Log[(Sqrt[r(r1-rp)]-Sqrt[(r1-r)rp])^2/(Sqrt[r(r1-rp)]+Sqrt[(r1-r)rp])^2];
	\[ScriptCapitalI]mreg=rm/(r1-rm)\[ScriptCapitalI]+Sqrt[rm]/(2(r1-rm)^(3/2))Log[(Sqrt[r(r1-rm)]-Sqrt[(r1-r)rm])^2/(Sqrt[r(r1-rm)]+Sqrt[(r1-r)rm])^2];
	(*\[ScriptCapitalI]mreg2=Sqrt[rm]/(r1-rm)\[ScriptCapitalI]+1/(2(r1-rm)^(3/2))Log[(Sqrt[r(r1-rm)]-Sqrt[(r1-r)rm])^2/(Sqrt[r(r1-rm)]+Sqrt[(r1-r)rm])^2];*)

	\[Phi]gPartialExp=-1/Sqrt[1-EE^2](Lz*\[ScriptCapitalI]-Lz*rp/(rp-rm)(rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)+2a*EE/(rp-rm)(rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg));
	
	\[Delta]\[Phi]PartialExp=-1/Sqrt[(1-EE^2)](-EE*\[ScriptCapitalI]-\[ScriptCapitalI]p/( 2(rp-rm)) (a(Lz-2a*EE)-2(2a*EE+Lz(rp-2))\[Delta]r4)+(\[ScriptCapitalI]mreg/rm)/( 2(rp-rm)) (-a(2a*EE-Lz)-2(2a*EE+Lz(rm-2))\[Delta]r4));
	
	\[Phi]gPartialExp+S*\[Delta]\[Phi]PartialExp
]


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


drd\[Lambda]funISCOplungePartialExp[r_,a_,S_,x_]:=Module[{EEg,\[Delta]EE,r1g,\[Delta]rISCO,r1,EE},
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,x];
	\[Delta]EE=\[Delta]EEISCOfun[a,x];
	\[Delta]rISCO=\[Delta]rISCOfun[a,x];

	EE=EEg+S*RealSign[x]*\[Delta]EE;
	r1=r1g+S*RealSign[x]*\[Delta]rISCO;

	(r-r1)Sqrt[1-EE^2]Sqrt[r(r1-r)]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity - dt/dr*)


dtdrfunISCOplungeOPartialExp[r_,a_,S_,p_,e_,x_]:=1/drd\[Lambda]funISCOplungePartialExp[r,a,S,x](VtrgfunPartialExp[r,a,S,p,e,x]+S*\[Delta]vtfunParPartialExp[r,a,S,p,e,x]+S*(\[Delta]r41fun[a,p,e,x]/r-a/(2r^2))VtrgfunPartialExp[r,a,S,p,e,x])


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity - d\[Phi]/dr*)


d\[Phi]drfunISCOplungePartialExp[r_,a_,S_,p_,e_,x_]:=1/drd\[Lambda]funISCOplungePartialExp[r,a,S,x](V\[Phi]rgfunPartialExp[r,a,S,p,e,x]+S*\[Delta]v\[Phi]funParPartialExp[r,a,S,p,e,x]+S*(\[Delta]r41fun[a,p,e,x]/r-a/(2r^2))V\[Phi]rgfunPartialExp[r,a,S,p,e,x])


(* ::Subsection::Closed:: *)
(*Critical plunge *)


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory *)


tfunCritplungePartialExp[r_,a_,S_,p_,e_,x_]:=Module[{r1g,r2g,rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]r1,\[Delta]r2,\[Delta]r4,EE,Lz,r1,r2,\[CapitalDelta],tgPartialExp,\[Delta]tPartialExp,\[ScriptCapitalI],\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over,\[ScriptCapitalI]r1g,\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	r1g=p/(1-e);
	r2g=p/(1+e);
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r4=\[Delta]r41fun[a,p,e,x];
	
	EE=EEg+S*RealSign[x]*\[Delta]EE;
	Lz=Lzg+S*RealSign[x]*\[Delta]Lz;
	r1=r1g+S*RealSign[x]*\[Delta]r1;
	r2=r2g+S*RealSign[x]*\[Delta]r2;
	\[CapitalDelta]=a^2-2r+r^2;

	\[ScriptCapitalI]=-1/(2Sqrt[(r1-r2)r2])Log[(Sqrt[r(r1-r2)]-Sqrt[(r1-r)r2])^2/(Sqrt[r(r1-r2)]+Sqrt[(r1-r)r2])^2];
	\[ScriptCapitalI]rover=2ArcCos[Sqrt[r/r1]]+r2*\[ScriptCapitalI];
	\[ScriptCapitalI]r2over=Sqrt[r(r1-r)]+r1*ArcCos[Sqrt[r/r1]]+r2*\[ScriptCapitalI]rover;
	\[ScriptCapitalI]p=\[ScriptCapitalI]/(r2-rp)+1/(2Sqrt[r1-rp](r2-rp)Sqrt[rp]) Log[(Sqrt[r(r1-rp)]-Sqrt[(r1-r)rp])^2/(Sqrt[r(r1-rp)]+Sqrt[(r1-r)rp])^2];
	\[ScriptCapitalI]mreg=(rm*\[ScriptCapitalI])/(r2-rm)+Sqrt[rm]/(2Sqrt[r1-rm](r2-rm))Log[(Sqrt[r(r1-rm)]-Sqrt[(r1-r)rm])^2/(Sqrt[r(r1-rm)]+Sqrt[(r1-r)rm])^2];
	
	tgPartialExp=-1/Sqrt[1-EE^2](4EE*\[ScriptCapitalI]+2EE*\[ScriptCapitalI]rover+EE*\[ScriptCapitalI]r2over-4EE*rp/(rp-rm)(rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)-2(-4EE+a*Lz)1/(rp-rm)(rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg));
	\[Delta]tPartialExp=-1/Sqrt[(1-EE^2)](1/2EE(4\[Delta]r4-a)\[ScriptCapitalI]+EE*\[Delta]r4*\[ScriptCapitalI]rover-\[ScriptCapitalI]p/( 2(rp-rm)) ((a^3*EE+2Lz*rp+a*EE(-6+rp)rp)-2(-2a*Lz+EE*rp^3+a^2*EE(2+rp))\[Delta]r4)+(\[ScriptCapitalI]mreg/rm)/( 2(rp-rm)) ((a^3*EE+2Lz*rm+a*EE(-6+rm)rm)-2(-2a*Lz+EE*rm^3+a^2*EE(2+rm))\[Delta]r4));
	
	tgPartialExp+S*\[Delta]tPartialExp
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal  trajectory*)


\[Phi]funCritplungePartialExp[r_,a_,S_,p_,e_,x_]:=Module[{r1g,r2g,rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]r1,\[Delta]r2,\[Delta]r4,EE,Lz,r1,r2,\[CapitalDelta],\[Phi]gPartialExp,\[Delta]\[Phi]PartialExp,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	r1g=p/(1-e);
	r2g=p/(1+e);
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r4=\[Delta]r41fun[a,p,e,x];
	
	EE=EEg+S*RealSign[x]*\[Delta]EE;
	Lz=Lzg+S*RealSign[x]*\[Delta]Lz;
	r1=r1g+S*RealSign[x]*\[Delta]r1;
	r2=r2g+S*RealSign[x]*\[Delta]r2;
	\[CapitalDelta]=a^2-2r+r^2;
	
	\[ScriptCapitalI]=-1/(2Sqrt[(r1-r2)r2])Log[(Sqrt[r(r1-r2)]-Sqrt[(r1-r)r2])^2/(Sqrt[r(r1-r2)]+Sqrt[(r1-r)r2])^2];
	\[ScriptCapitalI]p=\[ScriptCapitalI]/(r2-rp)+1/(2Sqrt[r1-rp](r2-rp)Sqrt[rp]) Log[(Sqrt[r(r1-rp)]-Sqrt[(r1-r)rp])^2/(Sqrt[r(r1-rp)]+Sqrt[(r1-r)rp])^2];
	\[ScriptCapitalI]mreg=(rm*\[ScriptCapitalI])/(r2-rm)+Sqrt[rm]/(2Sqrt[r1-rm](r2-rm))Log[(Sqrt[r(r1-rm)]-Sqrt[(r1-r)rm])^2/(Sqrt[r(r1-rm)]+Sqrt[(r1-r)rm])^2];

	\[Phi]gPartialExp=-1/Sqrt[1-EE^2](Lz*\[ScriptCapitalI]-Lz*rp/(rp-rm)(rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)+2a*EE/(rp-rm)(rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg));
	
	\[Delta]\[Phi]PartialExp=-1/Sqrt[(1-EE^2)](-EE*\[ScriptCapitalI]-\[ScriptCapitalI]p/( 2(rp-rm)) (a(Lz-2a*EE)-2(2a*EE+Lz(rp-2))\[Delta]r4)+(\[ScriptCapitalI]mreg/rm)/( 2(rp-rm)) (-a(2a*EE-Lz)-2(2a*EE+Lz(rm-2))\[Delta]r4));
	
	\[Phi]gPartialExp+S*\[Delta]\[Phi]PartialExp
]


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


drd\[Lambda]funHomPartialExp[r_,a_,S_,p_,e_,x_]:=Module[{EEg,\[Delta]EE,r1g,r2g,\[Delta]r1,\[Delta]r2,EE,r1,r2},
	EEg=EEgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	
	EE=EEg+S*RealSign[x]*\[Delta]EE;
	r1=r1g+S*RealSign[x]*\[Delta]r1;
	r2=r2g+S*RealSign[x]*\[Delta]r2;
	
	Sqrt[1-EE^2](r-r2)Sqrt[r(r1-r)]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity - dt/dr*)


dtdrfunCritplungePartialExp[r_,a_,S_,p_,e_,x_]:=1/drd\[Lambda]funHomPartialExp[r,a,S,p,e,x](VtrgfunPartialExp[r,a,S,p,e,x]+S*\[Delta]vtfunParPartialExp[r,a,S,p,e,x]+S*(\[Delta]r41fun[a,p,e,x]/r-a/(2r^2))VtrgfunPartialExp[r,a,S,p,e,x])


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity - d\[Phi]/dr*)


d\[Phi]drfunCritplungePartialExp[r_,a_,S_,p_,e_,x_]:=1/drd\[Lambda]funHomPartialExp[r,a,S,p,e,x](V\[Phi]rgfunPartialExp[r,a,S,p,e,x]+S*\[Delta]v\[Phi]funParPartialExp[r,a,S,p,e,x]+S*(\[Delta]r41fun[a,p,e,x]/r-a/(2r^2))V\[Phi]rgfunPartialExp[r,a,S,p,e,x])


(* ::Subsection::Closed:: *)
(*Generic plunge*)


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory*)


tfunPlungePartialExp[r_,a_,S_,p_,e_,x_]:=Module[{r1g,\[Rho]rg,\[Rho]ig,rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]r1,\[Delta]\[Rho]r,\[Delta]p,\[Delta]r4,EE,Lz,r1,\[Rho]r,\[Rho]i,\[CapitalDelta],Ar1,B0,kr,\[Alpha]p,\[Alpha]m,\[Gamma]r,\[Phi],ellF,ellE,ellPi,tgPartialExp,\[Delta]tPartialExp,\[ScriptCapitalI],\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over,\[ScriptCapitalI]r1g,\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	\[Delta]p=\[Delta]pfunFEcc[a,p,e,x];
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x]+dEEdpfuncc[p,e]\[Delta]p;
	\[Delta]Lz=\[Delta]Lzfuncc[a,p,e,x]+dLzdpfuncc[a,p,e,x]\[Delta]p;
	\[Delta]r1=\[Delta]p/(1-e);
	\[Delta]\[Rho]r=\[Delta]p/(1+e);
	\[Delta]r4=\[Delta]r41funcc[a,p,e,x];
	
	EE=EEg+S*RealSign[x]*\[Delta]EE;
	Lz=Lzg+S*RealSign[x]*\[Delta]Lz;
	r1=r1g+S*RealSign[x]*\[Delta]r1;
	\[Rho]r=\[Rho]rg+S*RealSign[x]*\[Delta]\[Rho]r;
	\[Rho]i=\[Rho]ig;
	\[CapitalDelta]=a^2-2r+r^2;

	Ar1=Sqrt[\[Rho]i^2+(r1-\[Rho]r)^2];
	B0=Sqrt[\[Rho]i^2+\[Rho]r^2];
	kr=(-(Ar1-B0)^2+r1^2)/(4Ar1*B0);
	\[Alpha]p=(-B0*r1-Ar1*rp+B0*rp)/(B0*r1-Ar1*rp-B0*rp);
	\[Alpha]m=(-B0*r1-Ar1*rm+B0*rm)/(B0*r1-Ar1*rm-B0*rm);
	\[Gamma]r=-((Ar1-B0)^2/(4Ar1*B0));
	\[Phi]=ArcCos[(-Ar1*r+B0(r1-r))/(Ar1*r+B0(r1-r))];
	
	ellF=EllipticF[\[Pi]-\[Phi],kr];
	ellE=EllipticE[\[Pi]-\[Phi],kr];
	ellPi=EllipticPi[\[Gamma]r,\[Pi]-\[Phi],kr];

	\[ScriptCapitalI]=ellF/Sqrt[Ar1*B0];
	\[ScriptCapitalI]rover=ArcTan[(r1*Sin[\[Phi]])/(2Sqrt[Ar1*B0] Sqrt[1-kr*Sin[\[Phi]]^2])]-(B0*r1)/(Ar1-B0)*\[ScriptCapitalI]+((Ar1+B0)r1*ellPi)/(2(Ar1-B0)Sqrt[Ar1*B0]);
	\[ScriptCapitalI]r2over=((Ar1-B0)Sqrt[r(r1-r)(\[Rho]i^2+(r-\[Rho]r)^2)])/(Ar1*r+B0(r1-r))+Sqrt[Ar1*B0]*ellE+1/2(r1+2\[Rho]r)\[ScriptCapitalI]rover-(B0*r1(r1-2\[Rho]r))/(2(Ar1-B0))*\[ScriptCapitalI];
	\[ScriptCapitalI]p=(4r1)/(B0^2(r1-rp)^2-Ar1^2*rp^2)(Ar1*B0)/Sqrt[2(Ar1*B0+B0^2-r1*\[Rho]r)](EllipticPi[\[Alpha]p^2,kr/(-1+kr)]+EllipticPi[\[Alpha]p^2,(\[Pi]/2-\[Phi]),kr/(-1+kr)])-(Ar1-B0)/(B0(r1-rp)+Ar1*rp)\[ScriptCapitalI]-1/(4Sqrt[(r1-rp)rp(\[Rho]i^2+(rp-\[Rho]r)^2)])Log[(Sqrt[-\[Alpha]p^2+kr(-1+\[Alpha]p^2)]Sin[\[Phi]]-Sqrt[1-\[Alpha]p^2] Sqrt[1-kr*Sin[\[Phi]]^2])^2/(Sqrt[-\[Alpha]p^2+kr(-1+\[Alpha]p^2)]Sin[\[Phi]]+Sqrt[1-\[Alpha]p^2]Sqrt[1-kr*Sin[\[Phi]]^2])^2];
	\[ScriptCapitalI]mreg=(4r1*rm)/(B0^2(r1-rm)^2-Ar1^2*rm^2)(Ar1*B0)/Sqrt[2(Ar1*B0+B0^2-r1*\[Rho]r)](EllipticPi[\[Alpha]m^2,kr/(-1+kr)]+EllipticPi[\[Alpha]m^2,(\[Pi]/2-\[Phi]),kr/(-1+kr)])-(rm(Ar1-B0))/(B0(r1-rm)+Ar1*rm)\[ScriptCapitalI]-Sqrt[rm]/(4Sqrt[(r1-rm)(\[Rho]i^2+(rm-\[Rho]r)^2)])Log[(Sqrt[-\[Alpha]m^2+kr(-1+\[Alpha]m^2)]Sin[\[Phi]]-Sqrt[1-\[Alpha]m^2] Sqrt[1-kr*Sin[\[Phi]]^2])^2/(Sqrt[-\[Alpha]m^2+kr(-1+\[Alpha]m^2)]Sin[\[Phi]]+Sqrt[1-\[Alpha]m^2] Sqrt[1-kr*Sin[\[Phi]]^2])^2];
	
	tgPartialExp=1/Sqrt[1-EE^2](4EE*\[ScriptCapitalI]+2EE*\[ScriptCapitalI]rover+EE*\[ScriptCapitalI]r2over-4EE*rp/(rp-rm)(rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)-2(-4EE+a*Lz)1/(rp-rm)(rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg));
	\[Delta]tPartialExp=1/Sqrt[(1-EE^2)](1/2EE(4\[Delta]r4-a)\[ScriptCapitalI]+EE*\[Delta]r4*\[ScriptCapitalI]rover-\[ScriptCapitalI]p/( 2(rp-rm)) ((a^3*EE+2Lz*rp+a*EE(-6+rp)rp)-2(-2a*Lz+EE*rp^3+a^2*EE(2+rp))\[Delta]r4)+(\[ScriptCapitalI]mreg/rm)/( 2(rp-rm)) ((a^3*EE+2Lz*rm+a*EE(-6+rm)rm)-2(-2a*Lz+EE*rm^3+a^2*EE(2+rm))\[Delta]r4));
	
	tgPartialExp+S*\[Delta]tPartialExp
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal trajectory*)


\[Phi]funPlungePartialExp[r_,a_,S_,p_,e_,x_]:=Module[{r1g,\[Rho]rg,\[Rho]ig,rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]r1,\[Delta]\[Rho]r,\[Delta]p,\[Delta]r4,EE,Lz,r1,\[Rho]r,\[Rho]i,\[CapitalDelta],Ar1,B0,kr,\[Alpha]p,\[Alpha]m,\[Gamma]r,\[Phi],ellF,ellE,\[Phi]gPartialExp,\[Delta]\[Phi]PartialExp,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	\[Delta]p=\[Delta]pfunFEcc[a,p,e,x];
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x]+dEEdpfuncc[p,e]\[Delta]p;
	\[Delta]Lz=\[Delta]Lzfuncc[a,p,e,x]+dLzdpfuncc[a,p,e,x]\[Delta]p;
	\[Delta]r1=\[Delta]p/(1-e);
	\[Delta]\[Rho]r=\[Delta]p/(1+e);
	\[Delta]r4=\[Delta]r41funcc[a,p,e,x];
	
	EE=EEg+S*RealSign[x]*\[Delta]EE;
	Lz=Lzg+S*RealSign[x]*\[Delta]Lz;
	r1=r1g+S*RealSign[x]*\[Delta]r1;
	\[Rho]r=\[Rho]rg+S*RealSign[x]*\[Delta]\[Rho]r;
	\[Rho]i=\[Rho]ig;
	\[CapitalDelta]=a^2-2r+r^2;

	Ar1=Sqrt[\[Rho]i^2+(r1-\[Rho]r)^2];
	B0=Sqrt[\[Rho]i^2+\[Rho]r^2];
	kr=(-(Ar1-B0)^2+r1^2)/(4Ar1*B0);
	\[Alpha]p=(-B0*r1-Ar1*rp+B0*rp)/(B0*r1-Ar1*rp-B0*rp);
	\[Alpha]m=(-B0*r1-Ar1*rm+B0*rm)/(B0*r1-Ar1*rm-B0*rm);
	\[Gamma]r=-((Ar1-B0)^2/(4Ar1*B0));
	\[Phi]=ArcCos[(-Ar1*r+B0(r1-r))/(Ar1*r+B0(r1-r))];
	
	ellF=EllipticF[\[Pi]-\[Phi],kr];
	ellE=EllipticE[\[Pi]-\[Phi],kr];

	\[ScriptCapitalI]=ellF/Sqrt[Ar1*B0];
	\[ScriptCapitalI]p=(4r1)/(B0^2(r1-rp)^2-Ar1^2*rp^2)(Ar1*B0)/Sqrt[2(Ar1*B0+B0^2-r1*\[Rho]r)](EllipticPi[\[Alpha]p^2,kr/(-1+kr)]+EllipticPi[\[Alpha]p^2,(\[Pi]/2-\[Phi]),kr/(-1+kr)])-(Ar1-B0)/(B0(r1-rp)+Ar1*rp)\[ScriptCapitalI]-1/(4Sqrt[(r1-rp)rp(\[Rho]i^2+(rp-\[Rho]r)^2)])Log[(Sqrt[-\[Alpha]p^2+kr(-1+\[Alpha]p^2)]Sin[\[Phi]]-Sqrt[1-\[Alpha]p^2] Sqrt[1-kr*Sin[\[Phi]]^2])^2/(Sqrt[-\[Alpha]p^2+kr(-1+\[Alpha]p^2)]Sin[\[Phi]]+Sqrt[1-\[Alpha]p^2]Sqrt[1-kr*Sin[\[Phi]]^2])^2];
	\[ScriptCapitalI]mreg=(4r1*rm)/(B0^2(r1-rm)^2-Ar1^2*rm^2)(Ar1*B0)/Sqrt[2(Ar1*B0+B0^2-r1*\[Rho]r)](EllipticPi[\[Alpha]m^2,kr/(-1+kr)]+EllipticPi[\[Alpha]m^2,(\[Pi]/2-\[Phi]),kr/(-1+kr)])-(rm(Ar1-B0))/(B0(r1-rm)+Ar1*rm)\[ScriptCapitalI]-Sqrt[rm]/(4Sqrt[(r1-rm)(\[Rho]i^2+(rm-\[Rho]r)^2)])Log[(Sqrt[-\[Alpha]m^2+kr(-1+\[Alpha]m^2)]Sin[\[Phi]]-Sqrt[1-\[Alpha]m^2] Sqrt[1-kr*Sin[\[Phi]]^2])^2/(Sqrt[-\[Alpha]m^2+kr(-1+\[Alpha]m^2)]Sin[\[Phi]]+Sqrt[1-\[Alpha]m^2] Sqrt[1-kr*Sin[\[Phi]]^2])^2];

	\[Phi]gPartialExp=1/Sqrt[1-EE^2](Lz*\[ScriptCapitalI]-Lz*rp/(rp-rm)(rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)+2a*EE/(rp-rm)(rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg));
	
	\[Delta]\[Phi]PartialExp=1/Sqrt[(1-EE^2)](-EE*\[ScriptCapitalI]-\[ScriptCapitalI]p/( 2(rp-rm)) (a(Lz-2a*EE)-2(2a*EE+Lz(rp-2))\[Delta]r4)+(\[ScriptCapitalI]mreg/rm)/( 2(rp-rm)) (-a(2a*EE-Lz)-2(2a*EE+Lz(rm-2))\[Delta]r4));
	
	\[Phi]gPartialExp+S*\[Delta]\[Phi]PartialExp
]


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


drd\[Lambda]funPlungePartialExp[r_,a_,S_,p_,e_,x_]:=Module[{EEg,\[Delta]EE,r1g,\[Rho]rg,\[Rho]ig,\[Delta]p,\[Delta]r1,\[Delta]\[Rho]r,EE,r1,\[Rho]r,\[Rho]i},
	EEg=EEgfuncc[p,e];
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];
	\[Delta]p=\[Delta]pfunFEcc[a,p,e,x];
	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x]+dEEdpfuncc[p,e]\[Delta]p;
	\[Delta]r1=\[Delta]p/(1-e);
	\[Delta]\[Rho]r=\[Delta]p/(1+e);
	
	EE=EEg+S*RealSign[x]*\[Delta]EE;
	r1=r1g+S*RealSign[x]*\[Delta]r1;
	\[Rho]r=\[Rho]rg+S*RealSign[x]*\[Delta]\[Rho]r;
	\[Rho]i=\[Rho]ig;
	
	-Sqrt[1-EE^2]Sqrt[(r-\[Rho]r-I*\[Rho]i)(r-\[Rho]r+I*\[Rho]i)]Sqrt[(r1-r)r]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity*)


VtrgfunPlungePartialExp[r_,a_,S_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]EE,\[Delta]Lz,EE,Lz,\[CapitalDelta]},
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	\[Delta]p=\[Delta]pfunFEcc[a,p,e,x];
	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x]+dEEdpfuncc[p,e]\[Delta]p;
	\[Delta]Lz=\[Delta]Lzfuncc[a,p,e,x]+dLzdpfuncc[a,p,e,x]\[Delta]p;
	
	EE=EEg+S*RealSign[x]*\[Delta]EE;
	Lz=Lzg+S*RealSign[x]*\[Delta]Lz;
	\[CapitalDelta]=a^2-2r+r^2;
	
	EE (r^2+a^2)^2/\[CapitalDelta]-Lz (2a*r)/\[CapitalDelta]-EE*a^2
]


\[Delta]vtfunPlungeParPartialExp[r_,a_,S_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]EE,\[Delta]Lz,EE,Lz,\[CapitalDelta]},
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	\[Delta]p=\[Delta]pfunFEcc[a,p,e,x];
	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x]+dEEdpfuncc[p,e]\[Delta]p;
	\[Delta]Lz=\[Delta]Lzfuncc[a,p,e,x]+dLzdpfuncc[a,p,e,x]\[Delta]p;
	
	EE=EEg+S*RealSign[x]*\[Delta]EE;
	Lz=Lzg+S*RealSign[x]*\[Delta]Lz;
	\[CapitalDelta]=a^2-2r+r^2;
	
	(-Lz(a^2+r^2)+a*EE(a^2+3r^2))/(r*\[CapitalDelta])
]


dtdrfunPlungePartialExp[r_,a_,S_,p_,e_,x_]:=1/drd\[Lambda]funPlungePartialExp[r,a,S,p,e,x](VtrgfunPlungePartialExp[r,a,S,p,e,x]+S*\[Delta]vtfunPlungeParPartialExp[r,a,S,p,e,x]+S*(\[Delta]r41funcc[a,p,e,x]/r-a/(2r^2))VtrgfunPlungePartialExp[r,a,S,p,e,x])


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


V\[Phi]rgfunPlungePartialExp[r_,a_,S_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]EE,\[Delta]Lz,EE,Lz,\[CapitalDelta]},
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	\[Delta]p=\[Delta]pfunFEcc[a,p,e,x];
	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x]+dEEdpfuncc[p,e]\[Delta]p;
	\[Delta]Lz=\[Delta]Lzfuncc[a,p,e,x]+dLzdpfuncc[a,p,e,x]\[Delta]p;
	
	EE=EEg+S*RealSign[x]*\[Delta]EE;
	Lz=Lzg+S*RealSign[x]*\[Delta]Lz;
	\[CapitalDelta]=a^2-2r+r^2;

	a/\[CapitalDelta](EE(r^2+a^2)-a*Lz)-a*EE+Lz
]


\[Delta]v\[Phi]funPlungeParPartialExp[r_,a_,S_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]EE,\[Delta]Lz,EE,Lz,\[CapitalDelta]},
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	\[Delta]p=\[Delta]pfunFEcc[a,p,e,x];
	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x]+dEEdpfuncc[p,e]\[Delta]p;
	\[Delta]Lz=\[Delta]Lzfuncc[a,p,e,x]+dLzdpfuncc[a,p,e,x]\[Delta]p;
	
	EE=EEg+S*RealSign[x]*\[Delta]EE;
	Lz=Lzg+S*RealSign[x]*\[Delta]Lz;
	\[CapitalDelta]=a^2-2r+r^2;
	
	-EE+(a(-Lz+a*EE(1+r)))/(r*\[CapitalDelta])
]


d\[Phi]drfunPlungePartialExp[r_,a_,S_,p_,e_,x_]:=1/drd\[Lambda]funPlungePartialExp[r,a,S,p,e,x](V\[Phi]rgfunPlungePartialExp[r,a,S,p,e,x]+S*\[Delta]v\[Phi]funPlungeParPartialExp[r,a,S,p,e,x]+S*(\[Delta]r41funcc[a,p,e,x]/r-a/(2r^2))V\[Phi]rgfunPlungePartialExp[r,a,S,p,e,x])


(* ::Subsection::Closed:: *)
(*Plunge related to bound orbits*)


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory*)


tfunPlungeBoundOrbitPartialExp[r_,a_,S_,p_,e_,x_]:=Module[{r1g,r2g,r3g,rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]r1,\[Delta]r2,\[Delta]r4,EE,Lz,r1,r2,r3,\[CapitalDelta],kr,\[Gamma]r,\[Phi],ellF,ellE,ellPi,tgPartialExp,\[Delta]tPartialExp,\[ScriptCapitalI],\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over,\[ScriptCapitalI]r1g,\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r4=\[Delta]r41fun[a,p,e,x];
	
	EE=EEg+S*RealSign[x]*\[Delta]EE;
	Lz=Lzg+S*RealSign[x]*\[Delta]Lz;
	r1=r1g+S*RealSign[x]*\[Delta]r1;
	r2=r2g+S*RealSign[x]*\[Delta]r2;
	r3=r3g+S*RealSign[x]*\[Delta]r2;
	\[CapitalDelta]=a^2-2r+r^2;
	
	kr=((r1-r2)r3)/(r2(r1-r3));
	\[Gamma]r=r3/r2;
	\[Phi]=ArcSin[\[Sqrt]((r2(r3-r))/(r3(r2-r)))];
	ellF=EllipticF[\[Phi],kr];
	ellE=EllipticE[\[Phi],kr];
	ellPi=EllipticPi[\[Gamma]r,\[Phi],kr];

	\[ScriptCapitalI]=(2ellF)/Sqrt[r2(r1-r3)];
	\[ScriptCapitalI]rover=2/Sqrt[r2(r1-r3)] (r2*ellF+(-r2+r3)ellPi);
	\[ScriptCapitalI]r2over=-r3*Sqrt[((r1-r)(r2-r))/(r2(r2-r3))]Sin[\[Phi]]Cos[\[Phi]]+Sqrt[r2(r1-r3)]ellE-((r2-r3)(r1+r2+r3))/Sqrt[r2(r1-r3)] ellPi+(r2(r2+r3))/Sqrt[r2(r1-r3)] ellF;
	\[ScriptCapitalI]p=2/(Sqrt[r2(r1-r3)](r2-rp))(ellF+(r2-r3)/(r3-rp) EllipticPi[(r3(-r2+rp))/(r2(-r3+rp)),\[Phi],kr]);
	\[ScriptCapitalI]mreg=(2rm)/(Sqrt[r2(r1-r3)](r2-rm))(ellF+(r2-r3)/(r3-rm) EllipticPi[(r3(-r2+rm))/(r2(-r3+rm)),\[Phi],kr]);
	
	tgPartialExp=-1/Sqrt[1-EE^2](4EE*\[ScriptCapitalI]+2EE*\[ScriptCapitalI]rover+EE*\[ScriptCapitalI]r2over-4EE*rp/(rp-rm)(rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)-2(-4EE+a*Lz)1/(rp-rm)(rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg));
	\[Delta]tPartialExp=-1/Sqrt[(1-EE^2)](1/2EE(4\[Delta]r4-a)\[ScriptCapitalI]+EE*\[Delta]r4*\[ScriptCapitalI]rover-\[ScriptCapitalI]p/( 2(rp-rm)) ((a^3*EE+2Lz*rp+a*EE(-6+rp)rp)-2(-2a*Lz+EE*rp^3+a^2*EE(2+rp))\[Delta]r4)+(\[ScriptCapitalI]mreg/rm)/( 2(rp-rm)) ((a^3*EE+2Lz*rm+a*EE(-6+rm)rm)-2(-2a*Lz+EE*rm^3+a^2*EE(2+rm))\[Delta]r4));
	
	tgPartialExp+S*\[Delta]tPartialExp
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal trajectory*)


\[Phi]funPlungeBoundOrbitPartialExp[r_,a_,S_,p_,e_,x_]:=Module[{r1g,r2g,r3g,rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]r1,\[Delta]r2,\[Delta]r4,EE,Lz,r1,r2,r3,\[CapitalDelta],kr,\[Phi],ellF,ellE,\[Phi]gPartialExp,\[Delta]\[Phi]PartialExp,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r4=\[Delta]r41fun[a,p,e,x];
	
	EE=EEg+S*RealSign[x]*\[Delta]EE;
	Lz=Lzg+S*RealSign[x]*\[Delta]Lz;
	r1=r1g+S*RealSign[x]*\[Delta]r1;
	r2=r2g+S*RealSign[x]*\[Delta]r2;
	r3=r3g+S*RealSign[x]*\[Delta]r2;
	\[CapitalDelta]=a^2-2r+r^2;
	
	kr=((r1-r2)r3)/(r2(r1-r3));
	\[Phi]=ArcSin[\[Sqrt]((r2(r3-r))/(r3(r2-r)))];
	ellF=EllipticF[\[Phi],kr];
	ellE=EllipticE[\[Phi],kr];
	
	\[ScriptCapitalI]=(2ellF)/Sqrt[r2(r1-r3)];
	\[ScriptCapitalI]p=2/(Sqrt[r2(r1-r3)](r2-rp))(ellF+(r2-r3)/(r3-rp) EllipticPi[(r3(-r2+rp))/(r2(-r3+rp)),\[Phi],kr]);
	\[ScriptCapitalI]mreg=(2rm)/(Sqrt[r2(r1-r3)](r2-rm))(ellF+(r2-r3)/(r3-rm) EllipticPi[(r3(-r2+rm))/(r2(-r3+rm)),\[Phi],kr]);

	\[Phi]gPartialExp=-1/Sqrt[1-EE^2](Lz*\[ScriptCapitalI]-Lz*rp/(rp-rm)(rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)+2a*EE/(rp-rm)(rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg));
	
	\[Delta]\[Phi]PartialExp=-1/Sqrt[(1-EE^2)](-EE*\[ScriptCapitalI]-\[ScriptCapitalI]p/( 2(rp-rm)) (a(Lz-2a*EE)-2(2a*EE+Lz(rp-2))\[Delta]r4)+(\[ScriptCapitalI]mreg/rm)/( 2(rp-rm)) (-a(2a*EE-Lz)-2(2a*EE+Lz(rm-2))\[Delta]r4));
	
	\[Phi]gPartialExp+S*\[Delta]\[Phi]PartialExp
]


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


drd\[Lambda]funPlungeBoundOrbitPartialExp[r_,a_,S_,p_,e_,x_]:=Module[{EEg,r1g,r2g,r3g,\[Delta]EE,\[Delta]r1,\[Delta]r2,EE,r1,r2,r3},
	EEg=EEgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	
	EE=EEg+S*RealSign[x]*\[Delta]EE;
	r1=r1g+S*RealSign[x]*\[Delta]r1;
	r2=r2g+S*RealSign[x]*\[Delta]r2;
	r3=r3g+S*RealSign[x]*\[Delta]r2;

	Sqrt[1-EE^2]Sqrt[(r-r2)(r-r3)]Sqrt[(r1-r)r]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity - dt/dr*)


dtdrfunPlungeBoundOrbitPartialExp[r_,a_,S_,p_,e_,x_]:=1/drd\[Lambda]funPlungeBoundOrbitPartialExp[r,a,S,p,e,x](VtrgfunPartialExp[r,a,S,p,e,x]+S*\[Delta]vtfunParPartialExp[r,a,S,p,e,x]+S*(\[Delta]r41fun[a,p,e,x]/r-a/(2r^2))VtrgfunPartialExp[r,a,S,p,e,x])


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity - d\[Phi]/dr*)


d\[Phi]drfunPlungeBoundOrbitPartialExp[r_,a_,S_,p_,e_,x_]:=1/drd\[Lambda]funPlungeBoundOrbitPartialExp[r,a,S,p,e,x](V\[Phi]rgfunPartialExp[r,a,S,p,e,x]+S*\[Delta]v\[Phi]funParPartialExp[r,a,S,p,e,x]+S*(\[Delta]r41fun[a,p,e,x]/r-a/(2r^2))V\[Phi]rgfunPartialExp[r,a,S,p,e,x])


(* ::Section::Closed:: *)
(*Spin corrections to the orbit - parallel component of the spin, regular gauge*)


(* ::Subsection::Closed:: *)
(*Bound trajectory*)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunRPerPar[wr_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,\[Delta]r1,\[Delta]r2,EEg,\[Delta]EE,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE,rg},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunR[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funRPer[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRPer[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Gamma]rg=(r1g-r2g)/(r1g-r3g);
	Yint=Sqrt[(r1g-r3g)r2g];

	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result.*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	jSN=Sin[\[Phi]];
	jCN=Cos[\[Phi]];

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	
	rg=rgICr1gfun[wr,a,p,e,x];
	
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	((jSN*jCN*r1g(r1g-r2g)r2g*Sqrt[r1g(r2g-jSN^2*r3g)-(1-jSN^2)r2g*r3g])/(jSN^2(r1g-r2g)+r2g)^2)(\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg*(2ellF)/Sqrt[(r1g-r3g)r2g]+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Sqrt[(r1g-r3g)r2g]+((r2g*ellE-r1g*ellF)\[Delta]r1)/(Sqrt[(r1g-r3g)r2g]r1g(r1g-r2g))+1/Sqrt[(r1g-r3g)r2g](-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g))\[Delta]r2+1/Sqrt[(r1g-r3g)r2g] 1/r3g(((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]))\[Delta]r3+2((-(r1g-r3g)ellE+r1g*ellF)/(r1g Sqrt[r2g(r1g-r3g)]r3g))\[Delta]r41-a/2 (-((4ellE Sqrt[r2g(r1g-r3g)](r2g*r3g+r1g(r2g+r3g)))/(3r1g^2r2g^2r3g^2))+(2ellF((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g Sqrt[r2g (r1g-r3g)]r3g^2)-(2(-r1g+r2g)(r1g-r3g)^2Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3r1g^2(r2g(r1g-r3g))^(3/2)r3g)))+(((rg-r2g)(rg-r3g))/((r1g-r2g)(r1g-r3g))\[Delta]r1+((r1g-rg)(rg-r3g))/((r1g-r2g)(r2g-r3g))\[Delta]r2)
]


\[Delta]rfunOfrRPerPar[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,EEg,\[Delta]EE,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],krg,ellF,ellE},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunR[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funRPer[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRPer[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Phi]=ArcSin[Sqrt[(r2g(r1g-r))/((r1g-r2g)r)]];
	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	
	Sqrt[(r1g-r)(r-r2g)r(r-r3g)](\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg*(2ellF)/Sqrt[(r1g-r3g)r2g]+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Sqrt[(r1g-r3g)r2g]+(r2g*ellE-r1g*ellF)/(Sqrt[(r1g-r3g)r2g]r1g(r1g-r2g) ) \[Delta]r1+1/Sqrt[(r1g-r3g)r2g] (-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g))\[Delta]r2+1/Sqrt[(r1g-r3g)r2g] 1/r3g (((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]))\[Delta]r3+ 2((-(r1g-r3g)ellE+r1g*ellF)/(r1g*r3g Sqrt[r2g(r1g-r3g)] ))\[Delta]r41-a/2 (-((4ellE*Sqrt[r2g(r1g-r3g)](r2g*r3g+r1g(r2g+r3g)))/(3r1g^2r2g^2r3g^2))+(2ellF((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g*r3g^2 Sqrt[r2g (r1g-r3g)] )-(2(-r1g+r2g)(r1g-r3g)Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3r1g^2r2g*r3g Sqrt[r2g(r1g-r3g)])))+(((r-r2g)(r-r3g))/((r1g-r2g)(r1g-r3g) ) \[Delta]r1+((r1g-r)(r-r3g))/((r1g-r2g)(r2g-r3g))\[Delta]r2)
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunRPerParRes[wr_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,\[Delta]r1,\[Delta]r2,EEg,\[Delta]EE,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE,rg},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunR[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funRPer[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRPer[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Gamma]rg=(r1g-r2g)/(r1g-r3g);
	Yint=Sqrt[(r1g-r3g)r2g];
	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result. Having \[Phi]=\[Pi]/2 with MachinePrecision cause some problems with the integrals \[ScriptCapitalI]p and \[ScriptCapitalI]m*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	jSN=Sin[\[Phi]];
	jCN=Cos[\[Phi]];

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	
	rg=rgICr1gfun[wr,a,p,e,x];
	
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg(2ellF)/Sqrt[(r1g-r3g)r2g]+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Sqrt[(r1g-r3g)r2g]+((r2g*ellE-r1g*ellF)\[Delta]r1)/(Sqrt[(r1g-r3g)r2g]r1g(r1g-r2g) )+1/Sqrt[(r1g-r3g)r2g](-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g))\[Delta]r2+1/Sqrt[(r1g-r3g)r2g] 1/r3g (((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]))\[Delta]r3+2((-(r1g-r3g)ellE+r1g*ellF)/(r1g*r3g*Sqrt[r2g(r1g-r3g)]))\[Delta]r41-a/2(-((4ellE Sqrt[r2g(r1g-r3g)](r2g*r3g+r1g(r2g+r3g)))/(3r1g^2r2g^2r3g^2))+(2ellF((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g*r3g^2 Sqrt[r2g(r1g-r3g)])-(2(-r1g+r2g)(r1g-r3g)Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3r1g^2r2g*r3g Sqrt[r2g(r1g-r3g)]))+Sign[\[Pi]-Mod[wr,2\[Pi]]](Sqrt[(-r2g+rg)(-r3g+rg)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-rg)rg])\[Delta]r1+Sqrt[(r1g-rg)(-r3g+rg)]/((r1g-r2g)(r2g-r3g)Sqrt[rg(-r2g+rg)])\[Delta]r2)
]


\[Delta]rfunOfrRPerParRes[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,EEg,\[Delta]EE,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],krg,\[Gamma]rg,ellF,ellE},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunR[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funRPer[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRPer[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Phi]=ArcSin[Sqrt[(r2g (r1g-r))/((r1g-r2g)r)]];
	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	
	(\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg*(2ellF)/Sqrt[(r1g-r3g)r2g]+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Sqrt[(r1g-r3g)r2g]+(r2g*ellE-r1g*ellF)/(Sqrt[(r1g-r3g)r2g]r1g(r1g-r2g))\[Delta]r1+1/Sqrt[(r1g-r3g)r2g](-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g))\[Delta]r2+1/Sqrt[(r1g-r3g)r2g] 1/r3g (((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]))\[Delta]r3+ 2((-(r1g-r3g)ellE+r1g*ellF)/(r1g*r3g*Sqrt[r2g(r1g-r3g)] ))\[Delta]r41-a/2(-((4ellE*Sqrt[r2g(r1g-r3g)](r2g*r3g+r1g(r2g+r3g)))/(3r1g^2r2g^2r3g^2))+(2ellF((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g*r3g^2*Sqrt[r2g (r1g-r3g)])-(2(-r1g+r2g)(r1g-r3g)Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3r1g^2r2g*r3g Sqrt[r2g(r1g-r3g)])))+(Sqrt[(-r2g+r)(-r3g+r)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-r)r])\[Delta]r1+Sqrt[(r1g-r)(-r3g+r)]/((r1g-r2g)(r2g-r3g)Sqrt[r(-r2g+r)])\[Delta]r2)
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory - purely oscillatory part*)


\[Delta]tfunRPerPar[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[CapitalUpsilon]rg,\[CapitalUpsilon]tg,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,\[Phi],rg,dtgd\[Lambda]fun,dtgd\[Lambda]funr1g,dtgd\[Lambda]funr2g,dtgd\[Lambda]funr3g,krghold,krg,\[Gamma]r,ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funRPer[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRPer[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];

	krg=((r1g-r2g)r3g)/(r2g (r1g-r3g));
	\[Gamma]r=1-r1g/r2g;

	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result.*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	rg=rgICr1gfun[wr,a,p,e,x];

	dtgd\[Lambda]fun=(rg(-2a*Lzg+EEg(2a^2+a^2*rg+rg^3)))/(a^2-2rg+rg^2);
	dtgd\[Lambda]funr1g=(r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(a^2-2r1g+r1g^2);
	dtgd\[Lambda]funr2g=(r2g(-2a*Lzg+EEg(2a^2+a^2*r2g+r2g^3)))/(a^2-2r2g+r2g^2);
	dtgd\[Lambda]funr3g=(r3g(-2a*Lzg+EEg(2a^2+a^2*r3g+r3g^3)))/(a^2-2r3g+r3g^2);

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	ellPi=EllipticPi[\[Gamma]r,\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[r2g(r1g-r3g)];
	\[ScriptCapitalI]r1g=2/(Sqrt[r2g(r1g-r3g)]r1g)((r2g*ellE-r1g*ellF)/(r1g-r2g));
	\[ScriptCapitalI]r2g=2/(Sqrt[r2g (r1g-r3g)])(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g));
	\[ScriptCapitalI]r3g=2/(Sqrt[r2g(r1g-r3g)]r3g)((r2g*ellE)/(r2g-r3g)-ellF+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]));
	\[ScriptCapitalI]rover=(2r1g)/(Sqrt[(r1g-r3g)r2g])ellPi;
	\[ScriptCapitalI]r2over=(r2g((r1g-r3g)ellE-r1g*ellF))/(Sqrt[r2g(r1g-r3g)])+(r1g(r1g+r2g+r3g)ellPi)/(Sqrt[r2g(r1g-r3g)])+((r1g-r2g)Sqrt[r2g(r1g-r3g)]Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/((r2g+(r1g-r2g)Sin[\[Phi]]^2));

	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	1/Sqrt[1-EEg^2]((4\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+1/2 EEg(\[Delta]r1+\[Delta]r2+\[Delta]r3+2\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+1/2 EEg((2+r1g)\[Delta]r1+(2+r2g)\[Delta]r2+(2+r3g)\[Delta]r3-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+1/2dtgd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+1/2dtgd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2+1/2dtgd\[Lambda]funr3g*\[ScriptCapitalI]r3g*\[Delta]r3+Sign[\[Pi]-Mod[wr,2\[Pi]]]dtgd\[Lambda]funr1g*(Sqrt[(-r2g+rg)(-r3g+rg)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-rg)rg]))\[Delta]r1+Sign[\[Pi]-Mod[wr,2\[Pi]]]dtgd\[Lambda]funr2g*(Sqrt[(r1g-rg)(-r3g+rg)]/((r1g-r2g)(r2g-r3g)Sqrt[rg(-r2g+rg)]))\[Delta]r2-dtgd\[Lambda]fun*\[Delta]rfunRPerParRes[wr,a,p,e,x])-1/\[CapitalUpsilon]rg(\[Delta]\[CapitalUpsilon]tfunR[a,p,e,x]-\[CapitalUpsilon]tg/\[CapitalUpsilon]rg*\[Delta]\[CapitalUpsilon]rfunR[a,p,e,x])wr
]	


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory as function of r*)


\[Delta]tfunOfrRPerPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[CapitalUpsilon]rg,\[CapitalUpsilon]tg,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Phi],\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,dtgd\[Lambda]fun,dtgd\[Lambda]funr1g,dtgd\[Lambda]funr2g,dtgd\[Lambda]funr3g,krg,\[Gamma]r,ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funRPer[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRPer[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];

	krg=((r1g-r2g)r3g)/(r2g (r1g-r3g));
	\[Gamma]r=1-r1g/r2g;
	\[Phi]=ArcSin[Sqrt[(r2g(r1g-r))/((r1g-r2g)r)]];

	dtgd\[Lambda]fun=(r(-2a*Lzg+EEg(2a^2+a^2*r+r^3)))/(a^2-2r+r^2);
	dtgd\[Lambda]funr1g=(r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(a^2-2r1g+r1g^2);
	dtgd\[Lambda]funr2g=(r2g(-2a*Lzg+EEg(2a^2+a^2*r2g+r2g^3)))/(a^2-2r2g+r2g^2);
	dtgd\[Lambda]funr3g=(r3g(-2a*Lzg+EEg(2a^2+a^2*r3g+r3g^3)))/(a^2-2r3g+r3g^2);

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	ellPi=EllipticPi[\[Gamma]r,\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[r2g(r1g-r3g)];
	\[ScriptCapitalI]r1g=2/(Sqrt[r2g(r1g-r3g)]r1g)((r2g*ellE-r1g*ellF)/(r1g-r2g));
	\[ScriptCapitalI]r2g=2/(Sqrt[r2g (r1g-r3g)])(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g));
	\[ScriptCapitalI]r3g=2/(Sqrt[r2g(r1g-r3g)]r3g)((r2g*ellE)/(r2g-r3g)-ellF+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]));
	\[ScriptCapitalI]rover=(2r1g)/(Sqrt[(r1g-r3g)r2g])ellPi;
	\[ScriptCapitalI]r2over=(r2g((r1g-r3g)ellE-r1g*ellF))/(Sqrt[r2g(r1g-r3g)])+(r1g(r1g+r2g+r3g)ellPi)/(Sqrt[r2g(r1g-r3g)])+((r1g-r2g)Sqrt[r2g(r1g-r3g)]Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/((r2g+(r1g-r2g)Sin[\[Phi]]^2));

	1/Sqrt[1-EEg^2]((4\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+1/2 EEg(\[Delta]r1+\[Delta]r2+\[Delta]r3+2\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+1/2EEg((2+r1g)\[Delta]r1+(2+r2g)\[Delta]r2+(2+r3g)\[Delta]r3-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+1/2dtgd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+1/2dtgd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2+1/2dtgd\[Lambda]funr3g*\[ScriptCapitalI]r3g*\[Delta]r3+dtgd\[Lambda]funr1g*(Sqrt[(-r2g+r)(-r3g+r)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-r)r]))\[Delta]r1+dtgd\[Lambda]funr2g*(Sqrt[(r1g-r)(-r3g+r)]/((r1g-r2g)(r2g-r3g)Sqrt[r(-r2g+r)]))\[Delta]r2-dtgd\[Lambda]fun*\[Delta]rfunOfrRPerParRes[r,a,p,e,x])
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal  trajectory - purely oscillatory part*)


\[Delta]\[Phi]funRPerPar[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Phi],rg,d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,d\[Phi]gd\[Lambda]funr2g,d\[Phi]gd\[Lambda]funr3g,krghold,krg,\[Gamma]r,ellF,ellE,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dLzdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funRPer[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRPer[a,p,e,x];
	\[Delta]r3=\[Delta]r2;

	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));
	
	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result.*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	rg=rgICr1gfun[wr,a,p,e,x];

	d\[Phi]gd\[Lambda]fun=((2a*EEg+Lzg(rg-2))rg)/(a^2-2rg+rg^2);
	d\[Phi]gd\[Lambda]funr1g=((2a*EEg+Lzg(-2+r1g))r1g)/(a^2-2r1g+r1g^2);
	d\[Phi]gd\[Lambda]funr2g=((2a*EEg+Lzg(-2+r2g))r2g)/(a^2-2r2g+r2g^2);
	d\[Phi]gd\[Lambda]funr3g=((2a*EEg+Lzg(-2+r3g))r3g)/(a^2-2r3g+r3g^2);

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[r2g(r1g-r3g)];
	\[ScriptCapitalI]r1g=2/(Sqrt[r2g(r1g-r3g)]r1g)((r2g*ellE-r1g*ellF)/(r1g-r2g));
	\[ScriptCapitalI]r2g=2/(Sqrt[r2g(r1g-r3g)])(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g));
	\[ScriptCapitalI]r3g=2/(Sqrt[r2g(r1g-r3g)]r3g)(((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]));
	
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	1/Sqrt[1-EEg^2](-EEg(1-(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+1/2d\[Phi]gd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+1/2d\[Phi]gd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2+1/2d\[Phi]gd\[Lambda]funr3g*\[ScriptCapitalI]r3g*\[Delta]r3+Sign[\[Pi]-Mod[wr,2\[Pi]]]d\[Phi]gd\[Lambda]funr1g*\[Delta]r1*(Sqrt[(-r2g+rg)(-r3g+rg)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-rg)rg]))+Sign[\[Pi]-Mod[wr,2\[Pi]]]d\[Phi]gd\[Lambda]funr2g*\[Delta]r2*(Sqrt[(r1g-rg)(-r3g+rg)]/((r1g-r2g)(r2g-r3g)Sqrt[rg(-r2g+rg)]))-d\[Phi]gd\[Lambda]fun*\[Delta]rfunRPerParRes[wr,a,p,e,x])-1/\[CapitalUpsilon]rg(\[Delta]\[CapitalUpsilon]\[Phi]funR[a,p,e,x]-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]rg \[Delta]\[CapitalUpsilon]rfunR[a,p,e,x])wr
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal  trajectory as function of r*)


\[Delta]\[Phi]funOfrRPerPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Phi],d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,d\[Phi]gd\[Lambda]funr2g,d\[Phi]gd\[Lambda]funr3g,krg,\[Gamma]r,ellF,ellE,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dLzdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funRPer[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRPer[a,p,e,x];
	\[Delta]r3=\[Delta]r2;

	krg=((r1g-r2g)r3g)/(r2g (r1g-r3g));
	\[Phi]=ArcSin[Sqrt[(r2g(r1g-r))/((r1g-r2g)r)]];

	d\[Phi]gd\[Lambda]fun=((2a*EEg+Lzg(r-2))r)/(a^2-2r+r^2);
	d\[Phi]gd\[Lambda]funr1g=((2a*EEg+Lzg(-2+r1g))r1g)/(a^2-2r1g+r1g^2);
	d\[Phi]gd\[Lambda]funr2g=((2a*EEg+Lzg(-2+r2g))r2g)/(a^2-2r2g+r2g^2);
	d\[Phi]gd\[Lambda]funr3g=((2a*EEg+Lzg(-2+r3g))r3g)/(a^2-2r3g+r3g^2);

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[r2g(r1g-r3g)];
	\[ScriptCapitalI]r1g=2/(Sqrt[r2g(r1g-r3g)]r1g)((r2g *ellE-r1g*ellF)/(r1g-r2g));
	\[ScriptCapitalI]r2g=2/(Sqrt[r2g(r1g-r3g)])(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g));
	\[ScriptCapitalI]r3g=2/(Sqrt[r2g(r1g-r3g)]r3g)(((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]));

	1/Sqrt[1-EEg^2](-EEg(1-(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+1/2d\[Phi]gd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+1/2d\[Phi]gd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2+1/2d\[Phi]gd\[Lambda]funr3g*\[ScriptCapitalI]r3g*\[Delta]r3+d\[Phi]gd\[Lambda]funr1g*\[Delta]r1(Sqrt[(-r2g+r)(-r3g+r)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-r)r]))+d\[Phi]gd\[Lambda]funr2g*\[Delta]r2(Sqrt[(r1g-r)(-r3g+r)]/((r1g-r2g)(r2g-r3g)Sqrt[r(-r2g+r)]))-d\[Phi]gd\[Lambda]fun*\[Delta]rfunOfrRPerParRes[r,a,p,e,x])
]


(* ::Subsection::Closed:: *)
(*Bound velocity*)


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunRPerPar[wr_,a_,p_,e_,x_]:=Module[{EEg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,dRgdr,rg},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funRPer[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRPer[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	dRgdr=Sqrt[(1-EEg^2)](-4rg^3+r1g*r2g*r3g+3rg^2(r1g+r2g+r3g)-2rg(r2g*r3g+r1g(r2g+r3g)));
	rg=rgICr1gfun[wr,a,p,e,x];

	-Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(1-EEg^2)rg(r1g-rg)(rg-r2g)(rg-r3g)](\[Delta]r1/(2(rg-r1g))+\[Delta]r2/(2(rg-r2g))+\[Delta]r3/(2(rg-r3g))+\[Delta]r41/rg-a/(2rg^2)+(EEg*\[Delta]EE)/(1-EEg^2))+1/2dRgdr*\[Delta]rfunRPerParRes[wr,a,p,e,x]
]


\[Delta]vrfunOfrRPerPar[r_,a_,p_,e_,x_]:=Module[{EEg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,dRgdr},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funRPer[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRPer[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	dRgdr=Sqrt[(1-EEg^2)](-4r^3+r1g*r2g*r3g+3r^2(r1g+r2g+r3g)-2r(r2g*r3g+r1g(r2g+r3g)));
	
	-Sqrt[(1-EEg^2)r(r1g-r)(r-r2g)(r-r3g)](\[Delta]r1/(2(r-r1g))+\[Delta]r2/(2(r-r2g))+\[Delta]r3/(2(r-r3g))+\[Delta]r41/r-a/(2r^2)+(EEg \[Delta]EE)/(1-EEg^2))+1/2dRgdr*\[Delta]rfunOfrRPerParRes[r,a,p,e,x]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity*)


\[Delta]vtfunRPerPar[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,rg,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dLzdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];

	rg=rgICr1gfun[wr,a,p,e,x];
	\[CapitalDelta]=a^2-2rg+rg^2;
	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+rg^2)+a*EEg(a^2+3rg^2))/(rg*\[CapitalDelta])+(rg(2a^2+a^2*rg+rg^3)\[Delta]EE)/\[CapitalDelta]-(2a*rg*\[Delta]Lz)/\[CapitalDelta];
	ddtgd\[Lambda]drfun=(EEg*rg(a^2+3rg^2))/\[CapitalDelta]-((rg^2-a^2)(-2a*Lzg+EEg*rg^3+a^2*EEg(2+rg)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunRPerPar[wr,a,p,e,x]+d\[Delta]td\[Lambda]fun
]


\[Delta]vtfunOfrRPerPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dLzdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];

	\[CapitalDelta]=a^2-2r+r^2;
	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+r^2)+a*EEg(a^2+3r^2))/(r*\[CapitalDelta])+(r(2a^2+a^2*r+r^3)\[Delta]EE)/\[CapitalDelta]-(2a*r*\[Delta]Lz)/\[CapitalDelta];
	ddtgd\[Lambda]drfun=(EEg*r(a^2+3r^2))/\[CapitalDelta]-((r^2-a^2)(-2a*Lzg+EEg*r^3+a^2*EEg(2+r)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunOfrRPerPar[r,a,p,e,x]+d\[Delta]td\[Lambda]fun
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


\[Delta]v\[Phi]funRPerPar[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,rg,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dLzdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];

	rg=rgICr1gfun[wr,a,p,e,x];
	\[CapitalDelta]=a^2-2rg+rg^2;
	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+rg)))/(rg*\[CapitalDelta])+(2a*rg*\[Delta]EE)/\[CapitalDelta]+((rg-2)rg*\[Delta]Lz)/\[CapitalDelta];
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(rg-1)-EEg*rg^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunRPerPar[wr,a,p,e,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


\[Delta]v\[Phi]funOfrRPerPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dLzdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];

	\[CapitalDelta]=a^2-2r+r^2;
	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+r)))/(r*\[CapitalDelta])+(2a*r*\[Delta]EE)/\[CapitalDelta]+((r-2)r*\[Delta]Lz)/\[CapitalDelta];
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(r-1)-EEg*r^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunOfrRPerPar[r,a,p,e,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Subsection::Closed:: *)
(*Homoclinic orbits*)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunRHomPar[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,\[Delta]r3,\[Delta]r41,\[Delta]r1,\[Delta]r2,EEg,\[Delta]EE,\[Delta]Lz,\[ScriptCapitalI],\[ScriptCapitalI]r,\[ScriptCapitalI]r2},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	\[Delta]r1=\[Delta]r1funRHom[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRHom[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	If[e==0,
		\[Delta]r1
		,
		\[ScriptCapitalI]=-(1/(2Sqrt[(r1g-r2g)r2g]))Log[(Sqrt[r(r1g-r2g)]-Sqrt[(r1g-r)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(r1g-r)r2g])^2];
		\[ScriptCapitalI]r=-2Sqrt[r1g-r]/(r1g*r2g*Sqrt[r])+\[ScriptCapitalI]/r2g;
		\[ScriptCapitalI]r2=-((2Sqrt[r(-r+r1g)])/(3r*r1g*r2g))(2/r1g+3/r2g+1/r)+\[ScriptCapitalI]/r2g^2;

		(r-r2g)Sqrt[(r1g-r)r](\[Delta]\[CapitalUpsilon]rover\[CapitalUpsilon]rgfunLimR[a,p,e,x]\[ScriptCapitalI]+(EEg*\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+(r2g*Sqrt[-r+r1g])/(r1g (r1g-r2g)^2*Sqrt[r])\[Delta]r1-\[ScriptCapitalI]/(2(r1g-r2g))\[Delta]r1-(\[ScriptCapitalI](r1g-2r2g)\[Delta]r2)/(2(r1g-r2g)r2g)+1/((r1g-r2g)r2g)Sqrt[-1+r1g/r]\[Delta]r2+\[Delta]r41*\[ScriptCapitalI]r-a/2*\[ScriptCapitalI]r2)+(r-r2g)((r-r2g)/(r1g-r2g)^2*\[Delta]r1)+Sqrt[(r1g-r)](Sqrt[r1g-r]/(r1g-r2g)\[Delta]r2)
	]
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunRHomParRes[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,\[Delta]r3,\[Delta]r41,\[Delta]r1,\[Delta]r2,EEg,\[Delta]EE,\[Delta]Lz,\[ScriptCapitalI],\[ScriptCapitalI]r,\[ScriptCapitalI]r2},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	\[Delta]r1=\[Delta]r1funRHom[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRHom[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	\[ScriptCapitalI]=-(1/(2Sqrt[(r1g-r2g)r2g]))Log[(Sqrt[r(r1g-r2g)]-Sqrt[(r1g-r)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(r1g-r)r2g])^2];
	\[ScriptCapitalI]r=-2Sqrt[r1g-r]/(r1g*r2g*Sqrt[r])+\[ScriptCapitalI]/r2g;
	\[ScriptCapitalI]r2=-((2Sqrt[r(-r+r1g)])/(3r*r1g*r2g))(2/r1g+3/r2g+1/r)+\[ScriptCapitalI]/r2g^2;

	(\[Delta]\[CapitalUpsilon]rover\[CapitalUpsilon]rgfunLimR[a,p,e,x]\[ScriptCapitalI]+(EEg*\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+(r2g*Sqrt[-r+r1g])/(r1g (r1g-r2g)^2*Sqrt[r])\[Delta]r1-\[ScriptCapitalI]/(2(r1g-r2g))\[Delta]r1-(\[ScriptCapitalI](r1g-2r2g)\[Delta]r2)/(2(r1g-r2g)r2g)+1/((r1g-r2g)r2g)Sqrt[-1+r1g/r]\[Delta]r2+\[Delta]r41*\[ScriptCapitalI]r-a/2*\[ScriptCapitalI]r2)+((r-r2g)/(Sqrt[(r1g-r)r](r1g-r2g)^2)\[Delta]r1)+(Sqrt[(r1g-r)r]/((r1g-r2g)(r-r2g)r)\[Delta]r2)
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory *)


\[Delta]tfunRHomPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,dtgd\[Lambda]fun,dtgd\[Lambda]funr1g,dtgd\[Lambda]funr2g,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funRHom[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRHom[a,p,e,x];
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];
	
	If[e==0,
		0,
		dtgd\[Lambda]fun=(r(-2a*Lzg+EEg(2a^2+a^2*r+r^3)))/(a^2-2r+r^2);
		dtgd\[Lambda]funr1g=(r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(a^2-2r1g+r1g^2);
		dtgd\[Lambda]funr2g=(r2g(-2a*Lzg+EEg(2a^2+a^2*r2g+r2g^3)))/(a^2-2r2g+r2g^2);
		
		\[ScriptCapitalI]=-(1/(2Sqrt[(r1g-r2g)r2g]))Log[(Sqrt[r(r1g-r2g)]-Sqrt[(r1g-r)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(r1g-r)r2g])^2];
		\[ScriptCapitalI]rover=2ArcCos[Sqrt[r/r1g]]+r2g*\[ScriptCapitalI];
		\[ScriptCapitalI]r2over=Sqrt[r(-r+r1g)]+(r1g+2r2g)ArcCos[Sqrt[r/r1g]]+r2g^2*\[ScriptCapitalI];
		\[ScriptCapitalI]r1g=2/(r1g(r1g-r2g)) Sqrt[r/(r1g-r)]-\[ScriptCapitalI]/(r1g-r2g);
		\[ScriptCapitalI]r2g=(Sqrt[r(r1g-r)])/((r-r2g)r2g(r1g-r2g))-((r1g-2r2g)\[ScriptCapitalI])/(2(r1g-r2g)r2g);
		
		1/Sqrt[1-EEg^2]((4*\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+1/2EEg(\[Delta]r1+\[Delta]r2+\[Delta]r2+2\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+1/2EEg((2+r1g)\[Delta]r1+(2+r2g)\[Delta]r2+(2+r2g)\[Delta]r2-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+(1/2)dtgd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+dtgd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2-dtgd\[Lambda]fun*\[Delta]rfunRHomParRes[r,a,p,e,x])
	]
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal trajectory*)


\[Delta]\[Phi]funRHomPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,d\[Phi]gd\[Lambda]funr2g,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dLzdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funRHom[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRHom[a,p,e,x];
	
	If[e==0,
		0,
		d\[Phi]gd\[Lambda]fun=((2a*EEg+Lzg(-2+r))r)/(a^2-2r+r^2);
		d\[Phi]gd\[Lambda]funr1g=((2a*EEg+Lzg(-2+r1g))r1g)/(a^2-2r1g+r1g^2);
		d\[Phi]gd\[Lambda]funr2g=((2a*EEg+Lzg(-2+r2g))r2g)/(a^2-2r2g+r2g^2);
		
		\[ScriptCapitalI]=-(1/(2Sqrt[(r1g-r2g)r2g]))Log[(Sqrt[r(r1g-r2g)]-Sqrt[(r1g-r)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(r1g-r)r2g])^2];
		\[ScriptCapitalI]r1g=2/(r1g(r1g-r2g))Sqrt[r/(r1g-r)]-\[ScriptCapitalI]/(r1g-r2g);
		\[ScriptCapitalI]r2g=(Sqrt[r(r1g-r)])/((r-r2g)r2g(r1g-r2g))-((r1g-2r2g)\[ScriptCapitalI])/(2(r1g-r2g)r2g);

		1/Sqrt[1-EEg^2](-EEg(1-(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+(1/2)d\[Phi]gd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+d\[Phi]gd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2-d\[Phi]gd\[Lambda]fun*\[Delta]rfunRHomParRes[r,a,p,e,x])
	]
]


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunRHomPar[r_,a_,p_,e_,x_]:=Module[{EEg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,\[Phi],dRgdr},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funRHom[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRHom[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	dRgdr=Sqrt[(1-EEg^2)](r2g-r)(r1g(r2g-3r)-2(r2g-2r)r);
	
	If[e==0,
		0,
		-(r-r2g)Sqrt[(1-EEg^2)r(r1g-r)](\[Delta]r1/(2(r-r1g))+\[Delta]r2/(r-r2g)+\[Delta]r41/r-a/(2r^2)+(EEg*\[Delta]EE)/(1-EEg^2))+1/2dRgdr*\[Delta]rfunRHomParRes[r,a,p,e,x]
	]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity*)


\[Delta]vtfunRHomPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dLzdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];
	
	\[CapitalDelta]=a^2-2r+r^2;
	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+r^2)+a*EEg(a^2+3r^2))/(r*\[CapitalDelta])+(r(2a^2+a^2*r+r^3)\[Delta]EE)/\[CapitalDelta]-(2a*r*\[Delta]Lz)/\[CapitalDelta];
	ddtgd\[Lambda]drfun=(EEg*r(a^2+3r^2))/\[CapitalDelta]-((r^2-a^2)(-2a*Lzg+EEg*r^3+a^2*EEg(2+r)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunRHomPar[r,a,p,e,x]+d\[Delta]td\[Lambda]fun
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


\[Delta]v\[Phi]funRHomPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dLzdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];
	
	\[CapitalDelta]=a^2-2r+r^2;
	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+r)))/(r*\[CapitalDelta])+(2a*r*\[Delta]EE)/\[CapitalDelta]+((r-2)r*\[Delta]Lz)/\[CapitalDelta];
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(r-1)-EEg*r^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunRHomPar[r,a,p,e,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Subsection::Closed:: *)
(*Critical plunge *)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunRCritplungePar[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,\[Delta]r41,\[Delta]r1,\[Delta]r2,EEg},
	EEg=EEgfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	\[Delta]r1=\[Delta]r1funRHom[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRHom[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	(r-r2g)Sqrt[r(-r+r1g)](-((2Sqrt[r1g-r])/(r1g*r2g*Sqrt[r]))\[Delta]r41-a/2(-((2Sqrt[r(-r+r1g)])/(3r*r1g*r2g))(2/r1g+3/r2g+1/r)))+(r(r-r2g)\[Delta]r1)/(r1g(r1g-r2g))+(r(-r+r1g)\[Delta]r2)/((r1g-r2g)r2g)
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunRCritplungeParRes[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,\[Delta]r41,\[Delta]r1,\[Delta]r2,EEg,Lzg},
	EEg=EEgfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	\[Delta]r1=\[Delta]r1funRHom[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRHom[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	 -((2Sqrt[r1g-r])/(r1g*r2g*Sqrt[r]))\[Delta]r41-a/2(-((2Sqrt[r(-r+r1g)])/(3r*r1g*r2g))(2/r1g+3/r2g+1/r))+\[Delta]r1/(r1g(r1g-r2g))Sqrt[r/(-r+r1g)]+(Sqrt[r(-r+r1g)]\[Delta]r2)/((r-r2g)r2g(r1g-r2g))
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory *)


\[Delta]tfunRCritplungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,dtgd\[Lambda]fun,dtgd\[Lambda]funr1g,dtgd\[Lambda]funr2g,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	\[Delta]r1=\[Delta]r1funRHom[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRHom[a,p,e,x];
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];
	
	dtgd\[Lambda]fun=(r(-2a*Lzg+EEg(2a^2+a^2*r+r^3)))/(a^2-2r+r^2);
	dtgd\[Lambda]funr1g=(r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(a^2-2r1g+r1g^2);
	dtgd\[Lambda]funr2g=(r2g(-2a*Lzg+EEg(2a^2+a^2*r2g+r2g^3)))/(a^2-2r2g+r2g^2);
	
	\[ScriptCapitalI]=-1/(2Sqrt[(r1g-r2g)r2g])Log[(Sqrt[r(r1g-r2g)]-Sqrt[(r1g-r)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(r1g-r)r2g])^2];
	\[ScriptCapitalI]r1g=2/(r1g(r1g-r2g))Sqrt[r/(r1g-r)]-\[ScriptCapitalI]/(r1g-r2g);
	\[ScriptCapitalI]r2g=(Sqrt[r(r1g-r)])/((r-r2g)r2g(r1g-r2g))-((r1g-2r2g)\[ScriptCapitalI])/(2(r1g-r2g)r2g);
	\[ScriptCapitalI]rover=2ArcCos[Sqrt[r/r1g]]+r2g*\[ScriptCapitalI];
	\[ScriptCapitalI]r2over=Sqrt[r(r1g-r)]+(r1g+2r2g)ArcCos[Sqrt[r/r1g]]+r2g^2\[ScriptCapitalI];

	-1/Sqrt[1-EEg^2]((4\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+1/2 EEg(\[Delta]r1+2\[Delta]r2+2\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+1/2*EEg((2+r1g)\[Delta]r1+2(2+r2g)\[Delta]r2-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+(1/2)dtgd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+dtgd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2-dtgd\[Lambda]fun*\[Delta]rfunRCritplungeParRes[r,a,p,e,x])
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal trajectory*)


\[Delta]\[Phi]funRCritplungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,\[Delta]r1,\[Delta]r2,d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,d\[Phi]gd\[Lambda]funr2g,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dLzdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	\[Delta]r1=\[Delta]r1funRHom[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRHom[a,p,e,x];
	
	d\[Phi]gd\[Lambda]fun=((2a*EEg+Lzg(-2+r))r)/(a^2-2r+r^2);
	d\[Phi]gd\[Lambda]funr1g=((2a*EEg+Lzg(-2+r1g))r1g)/(a^2-2r1g+r1g^2);
	d\[Phi]gd\[Lambda]funr2g=((2a*EEg+Lzg(-2+r2g))r2g)/(a^2-2r2g+r2g^2);
	
	\[ScriptCapitalI]=-1/(2Sqrt[(r1g-r2g)r2g])Log[(Sqrt[r(r1g-r2g)]-Sqrt[(r1g-r)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(r1g-r)r2g])^2];
	\[ScriptCapitalI]r1g=2/(r1g(r1g-r2g)) Sqrt[r/(r1g-r)]-\[ScriptCapitalI]/(r1g-r2g);
	\[ScriptCapitalI]r2g=(Sqrt[r(r1g-r)])/((r-r2g)r2g(r1g-r2g))-((r1g-2r2g)\[ScriptCapitalI])/(2(r1g-r2g)r2g);

	-1/Sqrt[1-EEg^2](-EEg(1-(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+(1/2)d\[Phi]gd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+d\[Phi]gd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2-d\[Phi]gd\[Lambda]fun*\[Delta]rfunRCritplungeParRes[r,a,p,e,x])
]


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunRCritplungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,r1g,r2g,\[Delta]r1,\[Delta]r2,\[Delta]r41,dRgdr},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];
	\[Delta]r1=\[Delta]r1funRHom[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRHom[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	dRgdr=Sqrt[(1-EEg^2)](r2g-r)(r1g(r2g-3r)-2(r2g-2r)r);
	
	-(r-r2g)Sqrt[(1-EEg^2)r(r1g-r)](\[Delta]r1/(2(r-r1g))+\[Delta]r2/(r-r2g)+\[Delta]r41/r-a/(2r^2)+(EEg*\[Delta]EE)/(1-EEg^2))+1/2dRgdr*\[Delta]rfunRCritplungeParRes[r,a,p,e,x]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity*)


\[Delta]vtfunRCritplungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dLzdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];

	\[CapitalDelta]=a^2-2r+r^2;
	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+r^2)+a*EEg(a^2+3r^2))/(r*\[CapitalDelta])+(r(2a^2+a^2*r+r^3)\[Delta]EE)/\[CapitalDelta]-(2a*r*\[Delta]Lz)/\[CapitalDelta];
	ddtgd\[Lambda]drfun=(EEg*r(a^2+3r^2))/\[CapitalDelta]-((r^2-a^2)(-2a*Lzg+EEg*r^3+a^2*EEg(2+r)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunRCritplungePar[r,a,p,e,x]+d\[Delta]td\[Lambda]fun
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


\[Delta]v\[Phi]funRCritplungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dLzdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];

	\[CapitalDelta]=a^2-2r+r^2;
	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+r)))/(r*\[CapitalDelta])+(2a*r*\[Delta]EE)/\[CapitalDelta]+((r-2)r*\[Delta]Lz)/\[CapitalDelta];
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(r-1)-EEg*r^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunRCritplungePar[r,a,p,e,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Subsection::Closed:: *)
(*Generic plunge*)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunRPlungePar[r_,a_,p_,e_,x_]:=Module[{EEg,\[Delta]p,\[Delta]e,\[Delta]r41,r1g,\[Rho]rg,\[Rho]ig,Ar1g,B0,\[Delta]r1,\[Delta]\[Rho]r,\[Phi],krg,ellF,ellE},
	EEg=EEgfuncc[p,e];
	\[Delta]p=\[Delta]pfunRPlunge[a,p,e,x];
	\[Delta]e=\[Delta]efunRPlunge[a,p,e,x];

	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];

	Ar1g=Sqrt[\[Rho]ig^2+(r1g-\[Rho]rg)^2];
	B0=Sqrt[\[Rho]ig^2+\[Rho]rg^2];
	krg=(-(Ar1g-B0)^2+r1g^2)/(4Ar1g*B0);
	\[Delta]r1=\[Delta]p/(1-e)+p*\[Delta]e/(1-e)^2;
	\[Delta]\[Rho]r=\[Delta]p/(1+e)-p*\[Delta]e/(1+e)^2;
	\[Delta]r41=\[Delta]r41funcc[a,p,e,x];

	\[Phi]=ArcCos[(-Ar1g*r+B0(r1g-r))/(Ar1g*r+B0(r1g-r))];
	ellE=EllipticE[\[Pi]-\[Phi],krg];
	
	Sqrt[(r1g-r)(\[Rho]ig^2+(r-\[Rho]rg)^2)r](1/Sqrt[Ar1g*B0])(\[Delta]r1/2((2B0*ellE)/(Ar1g*r1g)-((Ar1g*r+B0(r1g-r))(Ar1g*B0+r(r+r1g-2\[Rho]rg))Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(Ar1g*r1g*r(\[Rho]ig^2+(r-\[Rho]rg)^2)))+\[Delta]\[Rho]r((r1g-2\[Rho]rg)/(Ar1g*B0)ellE+(r1g(-2Ar1g*B0*r+(\[Rho]ig^2+(r-\[Rho]rg)^2)(r1g-2\[Rho]rg)+2Ar1g*B0*\[Rho]rg)Sin[\[Phi]])/(2Ar1g*B0(B0^2(-2r+r1g)-r^2(r1g-2\[Rho]rg))Sqrt[1-krg*Sin[\[Phi]]^2]))+\[Delta]r41(-((2Ar1g*ellE)/(B0*r1g))+((Ar1g*r+B0(-r+r1g))Sin[\[Phi]]Sqrt[1-krg Sin[\[Phi]]^2])/(B0*r*r1g))-a/2((2Sqrt[Ar1g*B0]Sqrt[r1g-r])/(3r1g*r Sqrt[r(\[Rho]ig^2+(r-\[Rho]rg)^2)])-(4Ar1g(B0^2+2r1g*\[Rho]rg)ellE)/(3B0^3*r1g^2)+(2(Ar1g*r+B0(-r+r1g))(B0^2+r1g*\[Rho]rg)Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3B0^3*r1g^2*r)-Sin[\[Phi]]/(3B0^3 Sqrt[1-krg*Sin[\[Phi]]^2]) ((-\[Rho]ig^2+\[Rho]rg^2)/B0+(4Ar1g*krg*\[Rho]rg*Cos[\[Phi]])/r1g)))+\[Delta]r1
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunRPlungeParRes[r_,a_,p_,e_,x_]:=Module[{EEg,\[Delta]p,\[Delta]e,\[Delta]r41,r1g,\[Rho]rg,\[Rho]ig,Ar1g,B0,\[Delta]r1,\[Delta]\[Rho]r,\[Phi],krg,ellF,ellE},
	EEg=EEgfuncc[p,e];
	\[Delta]p=\[Delta]pfunRPlunge[a,p,e,x];
	\[Delta]e=\[Delta]efunRPlunge[a,p,e,x];

	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];

	Ar1g=Sqrt[\[Rho]ig^2+(r1g-\[Rho]rg)^2];
	B0=Sqrt[\[Rho]ig^2+\[Rho]rg^2];
	krg=(-(Ar1g-B0)^2+r1g^2)/(4Ar1g*B0);
	\[Delta]r1=\[Delta]p/(1-e)+p*\[Delta]e/(1-e)^2;
	\[Delta]\[Rho]r=\[Delta]p/(1+e)-p*\[Delta]e/(1+e)^2;
	\[Delta]r41=\[Delta]r41funcc[a,p,e,x];

	\[Phi]=ArcCos[(-Ar1g*r+B0(r1g-r))/(Ar1g*r+B0(r1g-r))];
	ellE=EllipticE[\[Pi]-\[Phi],krg];

	(1/Sqrt[Ar1g*B0])(\[Delta]r1/2((2B0*ellE)/(Ar1g*r1g)-((Ar1g*r+B0(r1g-r))(Ar1g*B0+r(r+r1g-2\[Rho]rg))Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(Ar1g*r1g*r(\[Rho]ig^2+(r-\[Rho]rg)^2)))+\[Delta]\[Rho]r((r1g-2\[Rho]rg)/(Ar1g*B0)ellE+(r1g(-2Ar1g*B0*r+(\[Rho]ig^2+(r-\[Rho]rg)^2)(r1g-2\[Rho]rg)+2Ar1g*B0*\[Rho]rg)Sin[\[Phi]])/(2Ar1g*B0(B0^2(-2r+r1g)-r^2(r1g-2\[Rho]rg))Sqrt[1-krg*Sin[\[Phi]]^2]))+\[Delta]r41(-((2Ar1g*ellE)/(B0*r1g))+((Ar1g*r+B0(-r+r1g))Sin[\[Phi]]Sqrt[1-krg Sin[\[Phi]]^2])/(B0*r*r1g))-a/2((2Sqrt[Ar1g*B0]Sqrt[r1g-r])/(3r1g*r Sqrt[r(\[Rho]ig^2+(r-\[Rho]rg)^2)])-(4Ar1g(B0^2+2r1g*\[Rho]rg)ellE)/(3B0^3*r1g^2)+(2(Ar1g*r+B0(-r+r1g))(B0^2+r1g*\[Rho]rg)Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3B0^3*r1g^2*r)-Sin[\[Phi]]/(3B0^3 Sqrt[1-krg*Sin[\[Phi]]^2]) ((-\[Rho]ig^2+\[Rho]rg^2)/B0+(4Ar1g*krg*\[Rho]rg*Cos[\[Phi]])/r1g)))+\[Delta]r1/Sqrt[(r1g-r)(\[Rho]ig^2+(r-\[Rho]rg)^2)r]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory*)


\[Delta]tfunRPlungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]e,\[Delta]EE,\[Delta]Lz,r1g,r2g,r2gCC,\[Rho]rg,\[Rho]ig,\[Delta]r1,\[Delta]\[Rho]r,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,dtgd\[Lambda]fun,dtgd\[Lambda]funr1g,dtgd\[Lambda]funr2g,dtgd\[Lambda]funr2gCC,Ar1g,B0,krg,\[Gamma]r,\[Phi],ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over,\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r2gCC},
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	\[Delta]p=\[Delta]pfunRPlunge[a,p,e,x];
	\[Delta]e=\[Delta]efunRPlunge[a,p,e,x];
	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x]+dEEdpfuncc[p,e]\[Delta]p+dEEdefuncc[p,e]\[Delta]e;
	\[Delta]Lz=\[Delta]Lzfuncc[a,p,e,x]+dLzdpfuncc[a,p,e,x]\[Delta]p+dLzdefuncc[a,p,e,x]\[Delta]e;

	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];
	r2g=\[Rho]rg+I \[Rho]ig;
	r2gCC=\[Rho]rg-I \[Rho]ig;

	Ar1g=Sqrt[\[Rho]ig^2+(r1g-\[Rho]rg)^2];
	B0=Sqrt[\[Rho]ig^2+\[Rho]rg^2];
	krg=(-(Ar1g-B0)^2+r1g^2)/(4Ar1g*B0);
	\[Delta]r1=\[Delta]p/(1-e)+p*\[Delta]e/(1-e)^2;
	\[Delta]\[Rho]r=\[Delta]p/(1+e)-p*\[Delta]e/(1+e)^2;
	\[Delta]\[Rho]r4=\[Delta]r41funcc[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];

	\[Gamma]r=-((Ar1g-B0)^2/(4Ar1g*B0));
	\[Phi]=ArcCos[(-Ar1g*r+B0(r1g-r))/(Ar1g*r+B0(r1g-r))];

	dtgd\[Lambda]fun=(r(-2a*Lzg+EEg(2a^2+a^2*r+r^3)))/(a^2-2r+r^2);
	dtgd\[Lambda]funr1g=(r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(a^2-2r1g+r1g^2);
	dtgd\[Lambda]funr2g=(r2g(-2a*Lzg+EEg*r2g^3+a^2*EEg(r2g+2)))/(a^2+(r2g-2)r2g);
	dtgd\[Lambda]funr2gCC=(r2gCC(-2a*Lzg+EEg*r2gCC^3+a^2 EEg(r2gCC+2)))/(a^2+(r2gCC-2)r2gCC);
	ellF=EllipticF[\[Pi]-\[Phi],krg];
	ellE=EllipticE[\[Pi]-\[Phi],krg];
	ellPi=EllipticPi[\[Gamma]r,\[Pi]-\[Phi],krg];

	\[ScriptCapitalI]=ellF/(Sqrt[Ar1g*B0]);
	\[ScriptCapitalI]r1g=(2B0*ellE)/(Sqrt[Ar1g*B0]Ar1g*r1g)-(1/r1g+B0/(Ar1g*r1g))\[ScriptCapitalI]-((Ar1g*r+B0(r1g-r))(Ar1g*B0+r(r+r1g-2\[Rho]rg))Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(Sqrt[Ar1g*B0]*Ar1g*r1g*r(\[Rho]ig^2+(r-\[Rho]rg)^2));
	\[ScriptCapitalI]rover=ArcTan[(r1g*Sin[\[Phi]])/(2Sqrt[Ar1g*B0] Sqrt[1-krg*Sin[\[Phi]]^2])]-(B0*r1g)/(Ar1g-B0)*\[ScriptCapitalI]+((Ar1g+B0)r1g*ellPi)/(2(Ar1g-B0)Sqrt[Ar1g*B0]);
	\[ScriptCapitalI]r2over=((Ar1g-B0)Sqrt[r(r1g-r)(\[Rho]ig^2+(r-\[Rho]rg)^2)])/(Ar1g*r+B0(r1g-r))+Sqrt[Ar1g*B0]*ellE+1/2(r1g+2\[Rho]rg)\[ScriptCapitalI]rover-(B0*r1g(r1g-2\[Rho]rg))/(2(Ar1g-B0))*\[ScriptCapitalI];
	\[ScriptCapitalI]r2gCC=(2Sqrt[Ar1g*B0]r1g*ellE)/(B0^2(r1g-r2gCC)^2-Ar1g^2*r2gCC^2)-((Ar1g-B0)ellF)/(Sqrt[Ar1g*B0](B0(r1g-r2gCC)+Ar1g*r2gCC))+(r1g(r-r2g)Sin[\[Phi]])/(Sqrt[Ar1g*B0](-Ar1g*B0(r-r2gCC)+(r1g-r2gCC)r2gCC(r-r2g))Sqrt[1-krg*Sin[\[Phi]]^2]);
	\[ScriptCapitalI]r2g=(2Sqrt[Ar1g*B0]r1g*ellE)/(B0^2(r1g-r2g)^2-Ar1g^2*r2g^2)-((Ar1g-B0)ellF)/(Sqrt[Ar1g*B0](B0(r1g-r2g)+Ar1g*r2g))+(r1g(r-r2gCC)Sin[\[Phi]])/(Sqrt[Ar1g*B0](-Ar1g*B0(r-r2g)+(r1g-r2g)r2g(r-r2gCC))Sqrt[1-krg*Sin[\[Phi]]^2]);

	1/Sqrt[1-EEg^2](((4\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+1/2dtgd\[Lambda]funr1g*\[Delta]r1*\[ScriptCapitalI]r1g+1/2dtgd\[Lambda]funr2gCC*\[ScriptCapitalI]r2gCC*\[Delta]\[Rho]r+1/2dtgd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]\[Rho]r+1/2EEg((2+r1g)\[Delta]r1+2(2+\[Rho]rg)\[Delta]\[Rho]r-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+\[ScriptCapitalI]rover((2\[Delta]EE)/(1-EEg^2)+(EEg*\[Delta]r1)/2+EEg*\[Delta]\[Rho]r+EEg*\[Delta]\[Rho]r4)+dtgd\[Lambda]funr1g*\[Delta]r1(1/Sqrt[r(r1g-r)(r-I*\[Rho]ig-\[Rho]rg)(r+I*\[Rho]ig-\[Rho]rg)])-dtgd\[Lambda]fun*\[Delta]rfunRPlungeParRes[r,a,p,e,x])
]	


(* ::Subsubsection::Closed:: *)
(*Azimuthal trajectory*)


\[Delta]\[Phi]funRPlungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]e,\[Delta]EE,\[Delta]Lz,r1g,r2g,r2gCC,\[Rho]rg,\[Rho]ig,\[Delta]r1,\[Delta]\[Rho]r,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,d\[Phi]gd\[Lambda]funr2g,d\[Phi]gd\[Lambda]funr2gCC,Ar1g,B0,krg,\[Gamma]r,\[Phi],ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r2gCC},
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	\[Delta]p=\[Delta]pfunRPlunge[a,p,e,x];
	\[Delta]e=\[Delta]efunRPlunge[a,p,e,x];
	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x]+dEEdpfuncc[p,e]\[Delta]p+dEEdefuncc[p,e]\[Delta]e;
	\[Delta]Lz=\[Delta]Lzfuncc[a,p,e,x]+dLzdpfuncc[a,p,e,x]\[Delta]p+dLzdefuncc[a,p,e,x]\[Delta]e;

	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];
	r2g=\[Rho]rg+I \[Rho]ig;
	r2gCC=\[Rho]rg-I \[Rho]ig;

	Ar1g=Sqrt[\[Rho]ig^2+(r1g-\[Rho]rg)^2];
	B0=Sqrt[\[Rho]ig^2+\[Rho]rg^2];
	krg=(-(Ar1g-B0)^2+r1g^2)/(4Ar1g*B0);
	\[Delta]r1=\[Delta]p/(1-e)+p*\[Delta]e/(1-e)^2;
	\[Delta]\[Rho]r=\[Delta]p/(1+e)-p*\[Delta]e/(1+e)^2;

	\[Gamma]r=-((Ar1g-B0)^2/(4Ar1g*B0));
	\[Phi]=ArcCos[(-Ar1g*r+B0(r1g-r))/(Ar1g*r+B0(r1g-r))];

	d\[Phi]gd\[Lambda]fun=((2a*EEg+Lzg(-2+r))r)/(a^2-2r+r^2);
	d\[Phi]gd\[Lambda]funr1g=((2a*EEg+Lzg(-2+r1g))r1g)/(a^2-2r1g+r1g^2);
	d\[Phi]gd\[Lambda]funr2g=((2a*EEg+Lzg(-2+r2g))r2g)/(a^2-2r2g+r2g^2);
	d\[Phi]gd\[Lambda]funr2gCC=((2a*EEg+Lzg(-2+r2gCC))r2gCC)/(a^2-2r2gCC+r2gCC^2);
	ellF=EllipticF[\[Pi]-\[Phi],krg];
	ellE=EllipticE[\[Pi]-\[Phi],krg];
	ellPi=EllipticPi[\[Gamma]r,\[Pi]-\[Phi],krg];

	\[ScriptCapitalI]=ellF/(Sqrt[Ar1g*B0]);
	\[ScriptCapitalI]r1g=(2B0*ellE)/(Sqrt[Ar1g*B0]Ar1g*r1g)-(1/r1g+B0/(Ar1g*r1g))\[ScriptCapitalI]-((Ar1g*r+B0(r1g-r))(Ar1g*B0+r(r+r1g-2\[Rho]rg))Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(Sqrt[Ar1g*B0]*Ar1g*r1g*r(\[Rho]ig^2+(r-\[Rho]rg)^2));
	\[ScriptCapitalI]r2gCC=(2Sqrt[Ar1g*B0]r1g*ellE)/(B0^2(r1g-r2gCC)^2-Ar1g^2*r2gCC^2)-((Ar1g-B0)ellF)/(Sqrt[Ar1g*B0](B0(r1g-r2gCC)+Ar1g*r2gCC))+(r1g(r-r2g)Sin[\[Phi]])/(Sqrt[Ar1g*B0](-Ar1g*B0(r-r2gCC)+(r1g-r2gCC)r2gCC(r-r2g))Sqrt[1-krg*Sin[\[Phi]]^2]);
	\[ScriptCapitalI]r2g=(2Sqrt[Ar1g*B0]r1g*ellE)/(B0^2(r1g-r2g)^2-Ar1g^2*r2g^2)-((Ar1g-B0)ellF)/(Sqrt[Ar1g*B0](B0(r1g-r2g)+Ar1g*r2g))+(r1g(r-r2gCC)Sin[\[Phi]])/(Sqrt[Ar1g*B0](-Ar1g*B0(r-r2g)+(r1g-r2g)r2g(r-r2gCC))Sqrt[1-krg*Sin[\[Phi]]^2]);
	
	1/Sqrt[1-EEg^2]((-EEg+(EEg*Lzg*\[Delta]EE)/(1-EEg^2)+\[Delta]Lz)\[ScriptCapitalI]+1/2d\[Phi]gd\[Lambda]funr1g*\[Delta]r1*\[ScriptCapitalI]r1g+1/2\[ScriptCapitalI]r2g*d\[Phi]gd\[Lambda]funr2g*\[Delta]\[Rho]r+1/2\[ScriptCapitalI]r2gCC*d\[Phi]gd\[Lambda]funr2gCC*\[Delta]\[Rho]r+d\[Phi]gd\[Lambda]funr1g*\[Delta]r1(1/(Sqrt[r(r1g-r)(r-I*\[Rho]ig-\[Rho]rg)(r+I*\[Rho]ig-\[Rho]rg)]))-d\[Phi]gd\[Lambda]fun*\[Delta]rfunRPlungeParRes[r,a,p,e,x])
]


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunRPlungePar[r_,a_,p_,e_,x_]:=Module[{EEg,\[Delta]p,\[Delta]e,\[Delta]EE,\[Delta]r41,r1g,\[Rho]rg,\[Rho]ig,Ar1g,B0,\[Delta]r1,\[Delta]\[Rho]r,dRgdr},
	EEg=EEgfuncc[p,e];
	\[Delta]p=\[Delta]pfunRPlunge[a,p,e,x];
	\[Delta]e=\[Delta]efunRPlunge[a,p,e,x];
	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x]+dEEdpfuncc[p,e]\[Delta]p+dEEdefuncc[p,e]\[Delta]e;

	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];
	\[Delta]r1=\[Delta]p/(1-e)+p*\[Delta]e/(1-e)^2;
	\[Delta]\[Rho]r=\[Delta]p/(1+e)-p*\[Delta]e/(1+e)^2;
	\[Delta]r41=\[Delta]r41funcc[a,p,e,x];
	dRgdr=Sqrt[1-EEg^2](-4r^3+3r^2(r1g+2\[Rho]rg)+r1g(\[Rho]ig^2+\[Rho]rg^2)-2r(\[Rho]ig^2+\[Rho]rg(2 r1g+\[Rho]rg)));

	-Sqrt[(1-EEg^2)(r1g-r)(\[Rho]ig^2+(r-\[Rho]rg)^2)r](\[Delta]r1/(2(r-r1g))+(\[Delta]\[Rho]r(r-\[Rho]rg))/(r^2+\[Rho]ig^2-2r*\[Rho]rg+\[Rho]rg^2)+\[Delta]r41/r-a/(2r^2)+(EEg*\[Delta]EE)/(1-EEg^2))+1/2dRgdr*\[Delta]rfunRPlungeParRes[r,a,p,e,x]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity*)


\[Delta]vtfunRPlungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]e,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	\[Delta]p=\[Delta]pfunRPlunge[a,p,e,x];
	\[Delta]e=\[Delta]efunRPlunge[a,p,e,x];
	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x]+dEEdpfuncc[p,e]\[Delta]p+dEEdefuncc[p,e]\[Delta]e;
	\[Delta]Lz=\[Delta]Lzfuncc[a,p,e,x]+dLzdpfuncc[a,p,e,x]\[Delta]p+dLzdefuncc[a,p,e,x]\[Delta]e;
	
	\[CapitalDelta]=a^2-2r+r^2;	
	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+r^2)+a*EEg(a^2+3r^2))/(r*\[CapitalDelta])+(r(2a^2+a^2*r+r^3)\[Delta]EE)/\[CapitalDelta]-(2a*r*\[Delta]Lz)/\[CapitalDelta];
	ddtgd\[Lambda]drfun=(EEg*r(a^2+3r^2))/\[CapitalDelta]-((r^2-a^2)(-2a*Lzg+EEg*r^3+a^2*EEg(2+r)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunRPlungePar[r,a,p,e,x]+d\[Delta]td\[Lambda]fun
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


\[Delta]v\[Phi]funRPlungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]e,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	\[Delta]p=\[Delta]pfunRPlunge[a,p,e,x];
	\[Delta]e=\[Delta]efunRPlunge[a,p,e,x];
	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x]+dEEdpfuncc[p,e]\[Delta]p+dEEdefuncc[p,e]\[Delta]e;
	\[Delta]Lz=\[Delta]Lzfuncc[a,p,e,x]+dLzdpfuncc[a,p,e,x]\[Delta]p+dLzdefuncc[a,p,e,x]\[Delta]e;

	\[CapitalDelta]=a^2-2r+r^2;
	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+r)))/(r*\[CapitalDelta])+(2a*r*\[Delta]EE)/\[CapitalDelta]+((r-2)r*\[Delta]Lz)/\[CapitalDelta];
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(r-1)-EEg*r^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunRPlungePar[r,a,p,e,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Subsection::Closed:: *)
(*Plunge related to bound orbits*)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunRPlungeBoundOrbitPar[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,EEg,\[Delta]EE,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],krg,ellF,ellE},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunR[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funRPer[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRPer[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Phi]=ArcSin[Sqrt[(r2g (r3g-r))/(r3g (r2g-r))]];
	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	
	Sqrt[(r1g-r)(r-r2g)r(r-r3g)](\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg(2ellF)/Sqrt[(r1g-r3g)r2g]+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Sqrt[(r1g-r3g)r2g]+1/2 (2/Sqrt[r2g(r1g-r3g)] 1/(r1g-r2g)((r2g ellE)/r1g-ellF-((r1g-r2g)r3g*Cos[\[Phi]]*Sin[\[Phi]])/(r1g(r1g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2])))\[Delta]r1+1/2(2/((r1g-r2g)Sqrt[r2g(r1g-r3g)]) (-(((r1g-r3g)ellE)/(r2g-r3g))+ellF))\[Delta]r2+1/2 ((2r2g*ellE)/(Sqrt[r2g(r1g-r3g)](r2g-r3g)r3g)-(2 ellF)/(Sqrt[r2g(r1g-r3g)]r3g))\[Delta]r3+(2/(Sqrt[r2g (r1g-r3g)]r3g)(ellF+(-1+r3g/r1g)(ellE-Sqrt[1-krg*Sin[\[Phi]]^2]Tan[\[Phi]])))\[Delta]r41-a/2((2((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g Sqrt[r2g(r1g-r3g)]r3g^2) ellF-2/(3Sqrt[r2g(r1g-r3g)]) (r2g*r3g+r1g(r2g+r3g))/(r1g^2r2g^2r3g^2) r2g(r1g-r3g)(2ellE-Tan[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])+(2r3g)/(3r1g^2(r1g-r3g)Sqrt[r2g(r2g-r3g)]) Sqrt[(r2g-r)/(r1g-r)]Sin[\[Phi]]Cos[\[Phi]]-(2 (r3g^2+r3g*r+r^2-r1g(r3g+r)))/(3(r1g-r3g)r3g^2*r) Sqrt[(r3g-r)/(r(r1g-r)(r2g-r))]))+1/2((2(r1g-r)r)/((r1g-r3g)r3g)\[Delta]r3)
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunRPlungeBoundOrbitParRes[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]KK,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],krg,\[Gamma]rg,ellF,ellE},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunR[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funRPer[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRPer[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Phi]=ArcSin[Sqrt[(r2g(r3g-r))/(r3g (r2g-r))]];
	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	
	\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg(2ellF)/Sqrt[(r1g-r3g)r2g]+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Sqrt[(r1g-r3g)r2g]+1/2(2/Sqrt[r2g(r1g-r3g)] 1/(r1g-r2g)((r2g ellE)/r1g-ellF-((r1g-r2g)r3g*Cos[\[Phi]]*Sin[\[Phi]])/(r1g(r1g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2])))\[Delta]r1+1/2(2/((r1g-r2g)Sqrt[r2g(r1g-r3g)])(-(((r1g-r3g)ellE)/(r2g-r3g))+ellF))\[Delta]r2+1/2((2r2g*ellE)/(Sqrt[r2g(r1g-r3g)](r2g-r3g)r3g)-(2 ellF)/(Sqrt[r2g(r1g-r3g)]r3g))\[Delta]r3+(2/(Sqrt[r2g(r1g-r3g)]r3g) (ellF+(-1+r3g/r1g)(ellE-Sqrt[1-krg*Sin[\[Phi]]^2]Tan[\[Phi]])))\[Delta]r41-a/2(((2((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g Sqrt[r2g(r1g-r3g)]r3g^2) ellF-2/(3Sqrt[r2g(r1g-r3g)]) (r2g*r3g+r1g(r2g+r3g))/(r1g^2r2g^2r3g^2) r2g(r1g-r3g)(2ellE-Tan[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])+(2r3g)/(3r1g^2(r1g-r3g)Sqrt[r2g(r2g-r3g)]) Sqrt[(r2g-r)/(r1g-r)]Sin[\[Phi]]Cos[\[Phi]]-(2 (r3g^2+r3g*r+r^2-r1g(r3g+r)))/(3(r1g-r3g)r3g^2*r) Sqrt[(r3g-r)/(r(r1g-r)(r2g-r))]))+1/2(2Sqrt[(r1g-r)r])/((r1g-r3g)r3g Sqrt[(r2g-r)(r3g-r)])\[Delta]r3
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory*)


\[Delta]tfunRPlungeBoundOrbitPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Phi],\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,dtgd\[Lambda]fun,dtgd\[Lambda]funr1g,dtgd\[Lambda]funr2g,dtgd\[Lambda]funr3g,krg,\[Gamma]r,ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funRPer[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRPer[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];
	\[Delta]\[Rho]i4=Sqrt[a];

	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));
	\[Gamma]r=r3g/r2g;
	\[Phi]=ArcSin[\[Sqrt]((r2g(r3g-r))/(r3g(r2g-r)))];

	dtgd\[Lambda]fun=(r(-2a*Lzg+EEg(2a^2+a^2*r+r^3)))/(a^2-2r+r^2);
	dtgd\[Lambda]funr1g=(r1g(-2a*Lzg+EEg(2a^2+a^2r1g+r1g^3)))/(a^2-2r1g+r1g^2);
	dtgd\[Lambda]funr2g=(r2g(-2a*Lzg+EEg(2a^2+a^2r2g+r2g^3)))/(a^2-2r2g+r2g^2);
	dtgd\[Lambda]funr3g=(r3g(-2a*Lzg+EEg(2a^2+a^2r3g+r3g^3)))/(a^2-2r3g+r3g^2);

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	ellPi=EllipticPi[\[Gamma]r,\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[(r2g(r1g-r3g))];
	\[ScriptCapitalI]r1g=(2/Sqrt[(r2g(r1g-r3g))])(1/(r1g-r2g))((r2g*ellE)/r1g-ellF-((r1g-r2g)r3g*Cos[\[Phi]]*Sin[\[Phi]])/(r1g(r1g-r3g)Sqrt[(1-krg*Sin[\[Phi]]^2)]));
	\[ScriptCapitalI]r2g=(2/Sqrt[(r2g(r1g-r3g))])(1/(r1g-r2g))(-(((r1g-r3g)ellE)/(r2g-r3g))+ellF);
	\[ScriptCapitalI]r3g=(2/Sqrt[(r2g(r1g-r3g))]) 1/r3g ((r2g*ellE)/(r2g-r3g)-ellF);
	\[ScriptCapitalI]rover=2/Sqrt[r2g(r1g-r3g)] (r2g*ellF+(-r2g+r3g)ellPi);
	\[ScriptCapitalI]r2over=-r3g Sqrt[((r1g-r)(r2g-r))/(r2g(r2g-r3g))]Sin[\[Phi]]Cos[\[Phi]]+Sqrt[r2g(r1g-r3g)]ellE-((r2g-r3g)(r1g+r2g+r3g))/Sqrt[r2g(r1g-r3g)] ellPi+(r2g(r2g+r3g))/Sqrt[r2g(r1g-r3g)] ellF;

	-1/Sqrt[1-EEg^2]((4\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+1/2 EEg(\[Delta]r1+\[Delta]r2+\[Delta]r3+2\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+1/2 EEg((2+r1g)\[Delta]r1+(2+r2g)\[Delta]r2+(2+r3g)\[Delta]r3-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+1/2dtgd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+1/2dtgd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2+1/2dtgd\[Lambda]funr3g*\[ScriptCapitalI]r3g*\[Delta]r3+1/2dtgd\[Lambda]funr3g((2/(r1g-r3g))1/r3g((\[Sqrt]((r1g-r)r))/(\[Sqrt]((r2g-r)(r3g-r)))))\[Delta]r3-dtgd\[Lambda]fun*\[Delta]rfunRPlungeBoundOrbitParRes[r,a,p,e,x])
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal trajectory*)


\[Delta]\[Phi]funRPlungeBoundOrbitPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Phi],d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,d\[Phi]gd\[Lambda]funr2g,d\[Phi]gd\[Lambda]funr3g,krg,ellF,ellE,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dLzdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	(*\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];*)

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funRPer[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRPer[a,p,e,x];
	(*\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];*)
	\[Delta]r3=\[Delta]r2;

	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));
	\[Phi]=ArcSin[\[Sqrt]((r2g(r3g-r))/(r3g(r2g-r)))];

	d\[Phi]gd\[Lambda]fun=((2a*EEg+Lzg(r-2))r)/(a^2-2r+r^2);
	d\[Phi]gd\[Lambda]funr1g=((2a*EEg+Lzg(-2+r1g))r1g)/(a^2-2r1g+r1g^2);
	d\[Phi]gd\[Lambda]funr2g=((2a*EEg+Lzg(-2+r2g))r2g)/(a^2-2r2g+r2g^2);
	d\[Phi]gd\[Lambda]funr3g=((2a*EEg+Lzg(-2+r3g))r3g)/(a^2-2r3g+r3g^2);

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[(r2g(r1g-r3g))];
	\[ScriptCapitalI]r1g=(2/Sqrt[(r2g(r1g-r3g))])(1/(r1g-r2g))((r2g*ellE)/r1g-ellF-((r1g-r2g)r3g*Cos[\[Phi]]*Sin[\[Phi]])/(r1g(r1g-r3g)Sqrt[(1-krg*Sin[\[Phi]]^2)]));
	\[ScriptCapitalI]r2g=(2/Sqrt[(r2g(r1g-r3g))])(1/(r1g-r2g))(-(((r1g-r3g)ellE)/(r2g-r3g))+ellF);
	\[ScriptCapitalI]r3g=(2/Sqrt[(r2g(r1g-r3g))]) 1/r3g ((r2g*ellE)/(r2g-r3g)-ellF);

	-1/Sqrt[1-EEg^2](-EEg(1-(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+1/2d\[Phi]gd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+1/2d\[Phi]gd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2+1/2d\[Phi]gd\[Lambda]funr3g*\[ScriptCapitalI]r3g*\[Delta]r3+1/2d\[Phi]gd\[Lambda]funr3g((2/(r1g-r3g))1/r3g((\[Sqrt]((r1g-r)r))/(\[Sqrt]((r2g-r)(r3g-r)))))\[Delta]r3-d\[Phi]gd\[Lambda]fun*\[Delta]rfunRPlungeBoundOrbitParRes[r,a,p,e,x])
]


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunRPlungeBoundOrbitPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,dRgdr},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[Delta]r1=\[Delta]r1funRPer[a,p,e,x];
	\[Delta]r2=\[Delta]r2funRPer[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	dRgdr=Sqrt[(1-EEg^2)](-4r^3+r1g*r2g*r3g+3r^2(r1g+r2g+r3g)-2r(r2g*r3g+r1g(r2g+r3g)));

	Sqrt[(1-EEg^2)r(r1g-r)(r-r2g)(r-r3g)](\[Delta]r1/(2(r-r1g))+\[Delta]r2/(2(r-r2g))+\[Delta]r3/(2(r-r3g))+\[Delta]r41/r-a/(2r^2)+(EEg*\[Delta]EE)/(1-EEg^2))-1/2dRgdr*\[Delta]rfunRPlungeBoundOrbitParRes[r,a,p,e,x]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity*)


\[Delta]vtfunRPlungeBoundOrbitPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]td\[Lambda]fun,d\[Delta]td\[Lambda]fun2,ddtgd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dLzdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];

	\[CapitalDelta]=a^2-2r+r^2;	
	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+r^2)+a*EEg(a^2+3r^2))/(r*\[CapitalDelta])+(r(2a^2+a^2*r+r^3)\[Delta]EE)/\[CapitalDelta]-(2a*r*\[Delta]Lz)/\[CapitalDelta];
	ddtgd\[Lambda]drfun=(EEg*r(a^2+3r^2))/\[CapitalDelta]-((r^2-a^2)(-2a*Lzg+EEg*r^3+a^2*EEg(2+r)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunRPlungeBoundOrbitPar[r,a,p,e,x]+d\[Delta]td\[Lambda]fun
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


\[Delta]v\[Phi]funRPlungeBoundOrbitPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dLzdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];

	\[CapitalDelta]=a^2-2r+r^2;
	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+r)))/(r*\[CapitalDelta])+(2a*r*\[Delta]EE)/\[CapitalDelta]+((r-2)r*\[Delta]Lz)/\[CapitalDelta];
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(r-1)-EEg*r^2))/\[CapitalDelta]^2;
	
	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunRPlungeBoundOrbitPar[r,a,p,e,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Section::Closed:: *)
(*Spin corrections to the orbit - orthogonal component of the spin *)


(* ::Subsection::Closed:: *)
(*Spin precession phase*)


(* ::Subsubsection::Closed:: *)
(*Bound orbits*)


\[Psi]pref[a_,p_,e_,x_]:=Module[{EEg,Lzg,r1g,r2g,r3g,\[Phi],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun,krg,Lzred,\[Gamma]r,\[Psi]freq},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));
	\[Gamma]r=1-r1g/r2g;
	Lzred=Lzg-a*EEg;
	\[Psi]freq=(2r1g)/Sqrt[(1-EEg^2)r2g(r1g-r3g)](a*Lzred+EEg Lzred^2)/Abs[Lzred] (EllipticPi[-\[Gamma]r*I*Abs[Lzred]/((r1g-I*Abs[Lzred])),krg]/(r1g-I*Abs[Lzred])+EllipticPi[\[Gamma]r*I*Abs[Lzred]/((r1g+I*Abs[Lzred])),krg]/(r1g+I*Abs[Lzred]));
	1/(2\[Pi])*Re[\[Psi]freq]
]


(* ::Text:: *)
(*Purely oscillatory part of the precession phase*)


\[Psi]pICr2g[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,\[Phi],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun,krghold,krg,Lzred,ellF,\[Psi]freq},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	KKg=Lzg-a*EEg;
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));

	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result.*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	
	Lzred=Lzg-a*EEg;
	ellF=EllipticF[\[Phi],krg];
	\[Psi]freq=RealSign[Lzred](2(a+EEg*Lzred)ellF)/(Sqrt[1-EEg^2]*Sqrt[r2g(r1g-r3g)])-(I(-a*Lzred-EEg*Lzred^2))/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)])(ellF/(r3g-I*Abs[Lzred])-ellF/(r3g+I*Abs[Lzred])+((-r2g+r3g)EllipticPi[((r1g-r2g)(r3g-I*Abs[Lzred]))/((r1g-r3g)(r2g-I*Abs[Lzred])),\[Phi],krg])/((r3g-I*Abs[Lzred])(r2g-I*Abs[Lzred]))-((-r2g+r3g)EllipticPi[((r1g-r2g)(r3g+I*Abs[Lzred]))/((r1g-r3g)(r2g+I*Abs[Lzred])),\[Phi],krg])/((r3g+I*Abs[Lzred])(r2g+I*Abs[Lzred])));
	
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	Re[\[Psi]freq]-\[Psi]pref[a,p,e,x]wr
]


\[Psi]pICr1g[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,r1g,r2g,r3g,\[Phi],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun,krghold,krg,\[Gamma]r,Lzred,sgn,\[Psi]freq},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));
	\[Gamma]r=1-r1g/r2g;

	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result.*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;

	Lzred=Lzg-a*EEg;
	sgn=RealSign[Lzg-a*EEg];
	\[Psi]freq=r1g/Sqrt[(1-EEg^2)r2g(r1g-r3g)](a*Lzred+EEg*Lzred^2)/Abs[Lzred](EllipticPi[-\[Gamma]r*I*Abs[Lzred]/((r1g-I*Abs[Lzred])),\[Phi],krg]/(r1g-I*Abs[Lzred])+EllipticPi[\[Gamma]r*I*Abs[Lzred]/((r1g+I*Abs[Lzred])),\[Phi],krg]/(r1g+I*Abs[Lzred]));
	
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobi`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobi`Private`*"],True]]]];
	Re[\[Psi]freq]-\[Psi]pref[a,p,e,x]wr
]


\[Psi]pOfrICr1g[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,r1g,r2g,r3g,\[Phi],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun,krg,\[Gamma]r,Lzred,sgn,\[Psi]freq},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));
	\[Gamma]r=1-r1g/r2g;
	
	\[Phi]=ArcSin[Sqrt[(r2g (r1g-r))/((r1g-r2g)r)]];

	Lzred=Lzg-a*EEg;
	sgn=RealSign[Lzg-a*EEg];
	\[Psi]freq=r1g/Sqrt[(1-EEg^2)r2g(r1g-r3g)](a*Lzred+EEg*Lzred^2)/Abs[Lzred](EllipticPi[-\[Gamma]r*I*Abs[Lzred]/((r1g-I*Abs[Lzred])),\[Phi],krg]/(r1g-I*Abs[Lzred])+EllipticPi[\[Gamma]r*I*Abs[Lzred]/((r1g+I*Abs[Lzred])),\[Phi],krg]/(r1g+I*Abs[Lzred]));
	
	Re[\[Psi]freq]
]


(* ::Subsubsection::Closed:: *)
(*Homoclinic orbits*)


\[Psi]pHom[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,r1g,r2g,\[ScriptCapitalI],Lzred,sgn,\[Psi]freq},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);

	Lzred=Lzg-a*EEg;
	\[ScriptCapitalI]=-(1/(2Sqrt[(r1g-r2g)r2g]))Log[(Sqrt[r(r1g-r2g)]-Sqrt[(-r+r1g)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(-r+r1g)r2g])^2];
	\[Psi]freq=r2g^2/(r2g^2+Abs[Lzred]^2)\[ScriptCapitalI]-Exp[I*\[Pi]/4]/2 Sqrt[Abs[Lzred]](Log[(Sqrt[r(r1g+I*Abs[Lzred])]-Sqrt[I(r-r1g)] Sqrt[Abs[Lzred]])/(Sqrt[r(r1g+I*Abs[Lzred])]+Sqrt[I(r-r1g)] Sqrt[Abs[Lzred]])]/(Sqrt[r1g+I*Abs[Lzred]](-I*r2g+Abs[Lzred]))-Log[(Sqrt[r(r1g-I*Abs[Lzred])]-Sqrt[I*(-r+r1g)] Sqrt[Abs[Lzred]])/(Sqrt[r(r1g-I*Abs[Lzred])]+Sqrt[I*(-r+r1g)] Sqrt[Abs[Lzred]])]/(Sqrt[r1g-I*Abs[Lzred]](r2g-I*Abs[Lzred])));
	
	RealSign[Lzred]/Sqrt[(1-EEg^2)](a+EEg*Lzred)Re[\[Psi]freq]
]


(* ::Subsubsection::Closed:: *)
(*ISCO plunge*)


\[Psi]pISCOplunge[r_,a_,x_]:=Module[{EEg,Lzg,r1g,\[ScriptCapitalI],Lzred,sgn,\[Psi]freq},
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,x];
	Lzg=Lzgfun[a,r1g,0,x];
	Lzred=Lzg-a*EEg;

	\[ScriptCapitalI]=2Sqrt[r]/(r1g*Sqrt[r1g-r]);
	\[Psi]freq=r1g^2/(r1g^2+Abs[Lzred]^2)\[ScriptCapitalI]-Exp[I*\[Pi]/4]/2 Sqrt[Abs[Lzred]](Log[(Sqrt[r(r1g+I*Abs[Lzred])]-Sqrt[I(r-r1g)] Sqrt[Abs[Lzred]])/(Sqrt[r(r1g+I*Abs[Lzred])]+Sqrt[I(r-r1g)] Sqrt[Abs[Lzred]])]/(Sqrt[r1g+I*Abs[Lzred]](-I*r1g+Abs[Lzred]))-Log[(Sqrt[r(r1g-I*Abs[Lzred])]-Sqrt[I*(-r+r1g)] Sqrt[Abs[Lzred]])/(Sqrt[r(r1g-I*Abs[Lzred])]+Sqrt[I*(-r+r1g)] Sqrt[Abs[Lzred]])]/(Sqrt[r1g-I*Abs[Lzred]](r1g-I*Abs[Lzred])));
	
	-RealSign[Lzred]/Sqrt[(1-EEg^2)](a+EEg*Lzred)Re[\[Psi]freq]
]


(* ::Subsubsection::Closed:: *)
(*Critical plunge*)


\[Psi]pCritplunge[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,r1g,r2g,\[ScriptCapitalI],Lzred,sgn,\[Psi]freq},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	Lzred=Lzg-a*EEg;
	r1g=p/(1-e);
	r2g=p/(1+e);

	\[ScriptCapitalI]=-1/(2Sqrt[(r1g-r2g)r2g])Log[(Sqrt[r(r1g-r2g)]-Sqrt[(r1g-r)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(r1g-r)r2g])^2];
	\[Psi]freq=r2g^2/(r2g^2+Abs[Lzred]^2)\[ScriptCapitalI]-Exp[I*\[Pi]/4]/2 Sqrt[Abs[Lzred]](Log[(Sqrt[r(r1g+I*Abs[Lzred])]-Sqrt[I(r-r1g)] Sqrt[Abs[Lzred]])/(Sqrt[r(r1g+I*Abs[Lzred])]+Sqrt[I(r-r1g)] Sqrt[Abs[Lzred]])]/(Sqrt[r1g+I*Abs[Lzred]](-I*r2g+Abs[Lzred]))-Log[(Sqrt[r(r1g-I*Abs[Lzred])]-Sqrt[I*(-r+r1g)] Sqrt[Abs[Lzred]])/(Sqrt[r(r1g-I*Abs[Lzred])]+Sqrt[I*(-r+r1g)] Sqrt[Abs[Lzred]])]/(Sqrt[r1g-I*Abs[Lzred]](r2g-I*Abs[Lzred])));
	
	-RealSign[Lzred]/Sqrt[1-EEg^2](a+EEg*Lzred)Re[\[Psi]freq]
]


(* ::Subsubsection::Closed:: *)
(*Generic plunge*)


\[Psi]pPlunge[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,r1g,\[Rho]rg,\[Rho]ig,Ar1g,B0,krg,mkrg,\[Alpha]pLred,\[Alpha]mLred,\[Phi],ellPi,\[ScriptCapitalI],Lzred,\[Psi]freq},
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	Lzred=Lzg-a*EEg;
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];

	Ar1g=Sqrt[\[Rho]ig^2+(r1g-\[Rho]rg)^2];
	B0=Sqrt[\[Rho]ig^2+\[Rho]rg^2];
	krg=(-(Ar1g-B0)^2+r1g^2)/(4Ar1g*B0);
	mkrg=krg/(-1+krg);
	\[Alpha]pLred=(-I*B0*r1g+(Ar1g-B0)Abs[Lzred])/(I*B0*r1g+(Ar1g+B0)Abs[Lzred]);
	\[Alpha]mLred=(I*B0*r1g+(Ar1g-B0)Abs[Lzred])/(-I*B0*r1g+(Ar1g+B0)Abs[Lzred]);
	\[Phi]=ArcCos[(-Ar1g*r+B0(r1g-r))/(Ar1g*r+B0(r1g-r))];

	\[ScriptCapitalI]=EllipticF[\[Pi]-\[Phi],krg]/Sqrt[Ar1g*B0];
	\[Psi]freq=RealSign[Lzred] (B0^2(a+EEg*Lzred)r1g^2)/((Ar1g-B0)^2*Lzred^2+B0^2*r1g^2) \[ScriptCapitalI]+(-a*Lzred-EEg*Lzred^2)/(2I) (4Ar1g*B0*r1g)/Sqrt[(Ar1g+B0-r1g)(Ar1g+B0+r1g)] (-(EllipticPi[\[Alpha]mLred^2,mkrg]/(Ar1g^2*Lzred^2+B0^2(r1g+I*Abs[Lzred])^2))+EllipticPi[\[Alpha]pLred^2,mkrg]/(Ar1g^2*Lzred^2+B0^2(r1g-I*Abs[Lzred])^2)-EllipticPi[\[Alpha]mLred^2,\[Pi]/2-\[Phi],mkrg]/(Ar1g^2*Lzred^2+B0^2(r1g+I*Abs[Lzred])^2)+EllipticPi[\[Alpha]pLred^2,\[Pi]/2-\[Phi],mkrg]/(Ar1g^2*Lzred^2+B0^2(r1g-I*Abs[Lzred])^2))+(-a*Lzred-EEg*Lzred^2)/(4I*Abs[Lzred]) 1/ Sqrt[(-Lzred^2+\[Rho]ig^2+\[Rho]rg(2r1g+\[Rho]rg))] (Log[-(Sqrt[(-B0^2+2I \[Rho]rg Abs[Lzred]+Lzred^2)] Sqrt[r(r1g-r)]-Sqrt[Abs[Lzred]] Sqrt[-(I r1g+Abs[Lzred])] Sqrt[\[Rho]ig^2+(r-\[Rho]rg)^2])/(Sqrt[(-B0^2+2I \[Rho]rg Abs[Lzred]+Lzred^2)] Sqrt[r(r1g-r)]+Sqrt[Abs[Lzred]] Sqrt[-(I r1g+Abs[Lzred])] Sqrt[\[Rho]ig^2+(r-\[Rho]rg)^2])]-Log[-(Sqrt[(-B0^2-2I \[Rho]rg Abs[Lzred]+Lzred^2)] Sqrt[r(r1g-r)]-Sqrt[Abs[Lzred]] Sqrt[(I r1g- Abs[Lzred])] Sqrt[\[Rho]ig^2+(r-\[Rho]rg)^2])/(Sqrt[(-B0^2-2I \[Rho]rg Abs[Lzred]+Lzred^2)] Sqrt[r(r1g-r)]+Sqrt[Abs[Lzred]] Sqrt[I r1g- Abs[Lzred]] Sqrt[\[Rho]ig^2+(r-\[Rho]rg)^2])]);
	
	1/Sqrt[(1-EEg^2)]Re[\[Psi]freq]
]


(* ::Subsubsection::Closed:: *)
(*Plunge related to bound orbits*)


\[Psi]pPlungeBoundOrbit[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,r1g,r2g,r3g,\[Phi],\[Gamma]p,\[Gamma]m,krg,Lzred,\[Psi]freq},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	Lzred=Lzg-a*EEg;
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	Lzred=Lzg-a*EEg;

	\[Phi]=ArcSin[Sqrt[(r2g(r3g-r))/(r3g(r2g-r))]];
	\[Gamma]p=(r3g(r2g+I*Abs[Lzred]))/(r2g(r3g+I*Abs[Lzred]));
	\[Gamma]m=(r3g(r2g-I*Abs[Lzred]))/(r2g(r3g-I*Abs[Lzred]));
	\[Psi]freq=RealSign[Lzred] (2(a+EEg*Lzred)r2g^2*EllipticF[\[Phi],krg])/((Lzred^2+r2g^2)Sqrt[r2g(r1g-r3g)])+I*Lzred(a+EEg*Lzred) (r2g-r3g)/Sqrt[r2g(r1g-r3g)] (EllipticPi[\[Gamma]m,\[Phi],krg]/((r2g-I*Abs[Lzred])(r3g-I*Abs[Lzred]))-EllipticPi[\[Gamma]p,\[Phi],krg]/((r2g+I*Abs[Lzred])(r3g+I*Abs[Lzred])));
	-(1/Sqrt[1-EEg^2])Re[\[Psi]freq]
]


\[Psi]pPlungeBoundOrbitDoubleRoot[r_,a_,p_,x_]:=Module[{EEg,Lzg,r1g,r3g,\[Phi],\[Gamma]p,\[Gamma]m,Lzred,\[Psi]freq},
	EEg=EEgfun[a,p,0,x];
	Lzg=Lzgfun[a,p,0,x];
	r1g=p;
	r3g=2/(1-EEg^2)-2p;

	Lzred=Lzg-a*EEg;

	\[Phi]=ArcSin[Sqrt[(r1g(r3g-r))/(r3g(r1g-r))]];
	\[Gamma]p=(r3g(r1g+I*Abs[Lzred]))/(r1g(r3g+I*Abs[Lzred]));
	\[Gamma]m=(r3g(r1g-I*Abs[Lzred]))/(r1g(r3g-I*Abs[Lzred]));

	\[Psi]freq=RealSign[Lzred] (2(a+EEg*Lzred)r1g^2*\[Phi])/((Lzred^2+r1g^2)Sqrt[r1g(r1g-r3g)])+I*Lzred(a+EEg*Lzred) (r1g-r3g)/Sqrt[r1g(r1g-r3g)] (1/((r1g-I*Abs[Lzred])(r3g-I*Abs[Lzred])) ArcTanh[Sqrt[\[Gamma]m-1]Tan[\[Phi]]]/Sqrt[\[Gamma]m-1]-1/((r1g+I*Abs[Lzred])(r3g+I*Abs[Lzred])) ArcTanh[Sqrt[\[Gamma]p-1]Tan[\[Phi]]]/Sqrt[\[Gamma]p-1]);

	-(1/Sqrt[1-EEg^2])Re[\[Psi]freq]
]


(* ::Subsection::Closed:: *)
(*Polar trajectory*)


(* ::Subsubsection::Closed:: *)
(*Bound orbits*)


\[Delta]zfunICr2gPer[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,Lzred,\[CapitalUpsilon]rg,rg,\[Psi]p,\[CapitalUpsilon]p},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	Lzred=Lzg-a EEg;
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];

	rg=rgICr2gfun[wr,a,p,e,x];

	\[Psi]p=\[Psi]pICr2g[wr,a,p,e,x];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,x];

	-Cos[\[Psi]p+\[CapitalUpsilon]p/\[CapitalUpsilon]rg*wr]Sqrt[rg^2+Lzred^2]/(rg*Lzred)
]


\[Delta]zfunICr1gPer[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,Lzred,\[CapitalUpsilon]rg,r1g,r2g,r3g,\[Phi],rg,krg,\[Psi]p,\[CapitalUpsilon]p},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	Lzred=Lzg-a EEg;
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];

	rg=rgICr1gfun[wr,a,p,e,x];

	\[Psi]p=\[Psi]pICr1g[wr,a,p,e,x];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,x];

	-Cos[\[Psi]p+\[CapitalUpsilon]p/\[CapitalUpsilon]rg*wr]Sqrt[rg^2+Lzred^2]/(rg*Lzred)
]


\[Delta]zfunOfrICr1gPer[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,Lzred,\[CapitalUpsilon]rg,r1g,r2g,r3g,\[Psi]p,\[CapitalUpsilon]p},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	Lzred=Lzg-a EEg;
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];

	\[Psi]p=\[Psi]pOfrICr1g[r,a,p,e,x];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,x];

	-Cos[\[Psi]p]Sqrt[r^2+Lzred^2]/(r*Lzred)
]


(* ::Subsubsection::Closed:: *)
(*Homoclinic orbits*)


\[Delta]zfunHom[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,Lzred,\[Psi]p},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	Lzred=Lzg-a*EEg;

	\[Psi]p=\[Psi]pHom[r,a,p,e,x];

	-Cos[\[Psi]p]Sqrt[r^2+Lzred^2]/(r*Lzred)
]


(* ::Subsubsection::Closed:: *)
(*ISCO plunge*)


\[Delta]zfunISCOplunge[r_,a_,x_]:=Module[{EEg,Lzg,Lzred,r1g,\[Psi]p},
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,x];
	Lzg=Lzgfun[a,r1g,0,x];
	Lzred=Lzg-a*EEg;

	\[Psi]p=\[Psi]pISCOplunge[r,a,x];

	-Cos[\[Psi]p]Sqrt[r^2+Lzred^2]/(r*Lzred)
]


(* ::Subsubsection::Closed:: *)
(*Critical plunge*)


\[Delta]zfunCritplunge[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,Lzred,\[Psi]p},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	Lzred=Lzg-a*EEg;

	\[Psi]p=\[Psi]pCritplunge[r,a,p,e,x];

	-Cos[\[Psi]p]Sqrt[r^2+Lzred^2]/(r*Lzred)
]


(* ::Subsubsection::Closed:: *)
(*Generic plunge*)


\[Delta]zfunPlunge[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,Lzred,\[Psi]p},
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];
	Lzred=Lzg-a*EEg;

	\[Psi]p=\[Psi]pPlunge[r,a,p,e,x];

	-Cos[\[Psi]p]Sqrt[r^2+Lzred^2]/(r*Lzred)
]


(* ::Subsubsection::Closed:: *)
(*Plunge related to bound orbits*)


\[Delta]zfunPlungeBoundOrbit[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,Lzred,\[Psi]p},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	Lzred=Lzg-a*EEg;

	\[Psi]p=\[Psi]pPlungeBoundOrbit[r,a,p,e,x];

	-Cos[\[Psi]p]Sqrt[r^2+Lzred^2]/(r*Lzred)
]


\[Delta]zfunPlungeBoundOrbitDoubleRoot[r_,a_,p_,x_]:=Module[{EEg,Lzg,Lzred,\[Psi]p},
	EEg=EEgfun[a,p,0,x];
	Lzg=Lzgfun[a,p,0,x];
	Lzred=Lzg-a*EEg;

	\[Psi]p=\[Psi]pPlungeBoundOrbitDoubleRoot[r,a,p,x];

	-Cos[\[Psi]p]Sqrt[r^2+Lzred^2]/(r*Lzred)
]


(* ::Chapter::Closed:: *)
(*Front-end*)


(* ::Section::Closed:: *)
(*Near equatorial orbits - fixed turning points*)


(* ::Subsection::Closed:: *)
(*Analytic solutions*)


(* ::Subsubsection::Closed:: *)
(*Bound solutions*)


KerrNearEqSpinOrbitCorrFTPer[a_, p_, e_, x_]:=Module[{EEg,Lzg,\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g,\[CapitalUpsilon]p,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]zs,\[CapitalUpsilon]\[Phi]s},
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	
	(* geodesic frequencies *)
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,x];
	
	(* spin correction constants of motion *)
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x];
	  
	(* spin correction radial and polar frequencies *)
	\[CapitalUpsilon]rs=\[Delta]\[CapitalUpsilon]rfunFT[a,p,e,x];
	\[CapitalUpsilon]ts=\[Delta]\[CapitalUpsilon]tfunFT[a,p,e,x];
	\[CapitalUpsilon]\[Phi]s=\[Delta]\[CapitalUpsilon]\[Phi]funFT[a,p,e,x];
	
	(*The RealSign serves to match the sign convetion of generic orbits solutions in the limit z->0*)
	<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->RealSign[x]*\[Delta]EE,
	 "Jzs"->RealSign[x]*\[Delta]Lz,
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 "MinoFrequenciesCorrection"->RealSign[x]*{\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]\[Phi]s},
	 "BLFrequenciesCorrection"->RealSign[x]*{\[CapitalUpsilon]rs/\[CapitalUpsilon]tg-\[CapitalUpsilon]rg/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts,\[CapitalUpsilon]\[Phi]s/\[CapitalUpsilon]tg-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts},
	 "MinoPrecessionFrequency"->\[CapitalUpsilon]p,
	 "BLPrecessionFrequency"->\[CapitalUpsilon]p/\[CapitalUpsilon]tg,
	  (*Keys purely oscillatory part geodesic coordinate time and azimuthal trajectory*)
	 "\[CapitalDelta]trg"->Function[{wr},tgICr2gfun[wr,a,p,e,x]],
	 "\[CapitalDelta]\[Phi]rg"->Function[{wr},\[Phi]gICr2gfun[wr,a,p,e,x]],
	 (*Keys geodesic radial trajectory*)
	 "rg"->Function[{wr},RealSign[x]*rgICr2gfun[wr,a,p,e,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{wr},VtrgICr2gfun[wr,a,p,e,x]],
	 "vrg"->Function[{wr},drgd\[Lambda]fun[wr,a,p,e,x]],
	 "v\[Phi]g"->Function[{wr},V\[Phi]rgICr2gfun[wr,a,p,e,x]],
	 (*Keys corrections trajectory*)
	 "\[CapitalDelta]\[Delta]tpar"->Function[{wr},RealSign[x]*\[Delta]tfunFTPerPar[wr,a,p,e,x]],
	 "\[Delta]rpar"->Function[{wr},RealSign[x]*\[Delta]rfunFTPerPar[wr,a,p,e,x]],
	 "\[Psi]p"->Function[{wr},\[Psi]pICr2g[wr,a,p,e,x]],
	 "\[Delta]zort"->Function[{wr},\[Delta]zfunICr2gPer[wr,a,p,e,x]],
	 "\[CapitalDelta]\[Delta]\[Phi]par"->Function[{wr},RealSign[x]*\[Delta]\[Phi]funFTPerPar[wr,a,p,e,x]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{wr},RealSign[x]*\[Delta]vtfunFTPerPar[wr,a,p,e,x]],
	 "\[Delta]vrpar"->Function[{wr},RealSign[x]*\[Delta]vrfunFTPerPar[wr,a,p,e,x]],
	 "\[Delta]v\[Phi]par"->Function[{wr},RealSign[x]*\[Delta]v\[Phi]funFTPerPar[wr,a,p,e,x]]
	|>
]


(* ::Section::Closed:: *)
(*Near equatorial orbits - fixed constants of motion*)


(* ::Subsection::Closed:: *)
(*Analytic solutions*)


(* ::Subsubsection::Closed:: *)
(*Bound solutions*)


KerrNearEqSpinOrbitCorrFCPer[a_, p_, e_, x_]:=Module[{EEg,Lzg,\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g,\[CapitalUpsilon]p,\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]zs,\[CapitalUpsilon]\[Phi]s,\[Delta]r1,\[Delta]r2},
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	
	(* geodesic frequencies *)
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,x];
	   
	(* spin corrections radial and polar frequencies *)
	\[CapitalUpsilon]rs=\[Delta]\[CapitalUpsilon]rfunFC[a,p,e,x];
	\[CapitalUpsilon]ts=\[Delta]\[CapitalUpsilon]tfunFC[a,p,e,x];
	\[CapitalUpsilon]\[Phi]s=\[Delta]\[CapitalUpsilon]\[Phi]funFC[a,p,e,x];
	
	(* spin corrections turning points *)
	\[Delta]r1=\[Delta]r1funFC[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFC[a,p,e,x];
	
	(*The RealSign serves to match the sign convetion of generic orbits solutions in the limit z->0*)		
	<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->0,
	 "Jzs"->0,
	 (*Shift orbital elements*)
	 "\[Delta]p"->RealSign[x]*(1/2((1-e)^2*\[Delta]r1+(1+e)^2*\[Delta]r2)),
	 "\[Delta]e"->RealSign[x]*((1-e^2)((1-e)\[Delta]r1-\[Delta]r2-e*\[Delta]r2))/(2p),
	 (*Shift frequencies*)
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 "MinoFrequenciesCorrection"->RealSign[x]*{\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]\[Phi]s},
	 "BLFrequenciesCorrection"->{RealSign[x]*\[CapitalUpsilon]rs/\[CapitalUpsilon]tg-\[CapitalUpsilon]rg/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts,\[CapitalUpsilon]\[Phi]s/\[CapitalUpsilon]tg-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts},
	 "MinoPrecessionFrequency"->\[CapitalUpsilon]p,
	 "BLPrecessionFrequency"->\[CapitalUpsilon]p/\[CapitalUpsilon]tg,
	  (*Keys purely oscillatory part geodesic coordinate time and azimuthal trajectory*)
	 "\[CapitalDelta]trg"->Function[{wr},tgICr2gfun[wr,a,p,e,x]],
	 "\[CapitalDelta]\[Phi]rg"->Function[{wr},\[Phi]gICr2gfun[wr,a,p,e,x]],
	 (*Keys geodesic radial trajectory*)
	 "rg"->Function[{wr},rgICr2gfun[wr,a,p,e,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{wr},VtrgICr2gfun[wr,a,p,e,x]],
	 "vrg"->Function[{wr},drgd\[Lambda]fun[wr,a,p,e,x]],
	 "v\[Phi]g"->Function[{wr},V\[Phi]rgICr2gfun[wr,a,p,e,x]],
	 (*Keys corrections trajectory*)
	 "\[CapitalDelta]\[Delta]tpar"->Function[{wr},RealSign[x]*\[Delta]tfunFCPerPar[wr,a,p,e,x]],
	 "\[Delta]rpar"->Function[{wr},RealSign[x]*\[Delta]rfunFCPerPar[wr,a,p,e,x]],
	 "\[Psi]p"->Function[{wr},\[Psi]pICr2g[wr,a,p,e,x]],
	 "\[Delta]zort"->Function[{wr},\[Delta]zfunICr2gPer[wr,a,p,e,x]],
	 "\[CapitalDelta]\[Delta]\[Phi]par"->Function[{wr},\[Delta]\[Phi]funFCPerPar[wr,a,p,e,x]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{wr},RealSign[x]*\[Delta]vtfunFCPerPar[wr,a,p,e,x]],
	 "\[Delta]vrpar"->Function[{wr},RealSign[x]*\[Delta]vrfunFCPerPar[wr,a,p,e,x]],
	 "\[Delta]v\[Phi]par"->Function[{wr},RealSign[x]*\[Delta]v\[Phi]funFCPerPar[wr,a,p,e,x]]
	|>
]


(* ::Section::Closed:: *)
(*Near equatorial orbits - fixed frequency*)


(* ::Subsection::Closed:: *)
(*Analytic solutions*)


(* ::Subsubsection::Closed:: *)
(*Bound solutions*)


KerrNearEqSpinOrbitCorrFFPer[a_, p_, e_, x_]:=Module[{EEg,Lzg,\[Delta]r1,\[Delta]r2,\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g,\[CapitalUpsilon]p,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]zs,\[CapitalUpsilon]\[Phi]s},
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];	
	(* geodesic frequencies *)
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,x];	
	(* spin correction constants of motion *)
	{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzgdr1gfun[a,p,e,x]\[Delta]r1+dLzgdr2gfun[a,p,e,x]\[Delta]r2;	   
	(* spin correction radial and polar frequencies *)
	\[CapitalUpsilon]rs=\[Delta]\[CapitalUpsilon]rfunFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}];
	\[CapitalUpsilon]ts=\[Delta]\[CapitalUpsilon]tfunFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}];
	\[CapitalUpsilon]\[Phi]s=\[Delta]\[CapitalUpsilon]\[Phi]funFF[a,p,e,x,{\[Delta]r1,\[Delta]r2}];
	
	(*The RealSign serves to match the sign convetion of generic orbits solutions in the limit z->0*)				
	<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->RealSign[x]*\[Delta]EE,
	 "Jzs"->RealSign[x]*\[Delta]Lz,
	 (*Shift orbital elements*)
	 "\[Delta]p"->RealSign[x]*(1/2((1-e)^2*\[Delta]r1+(1+e)^2*\[Delta]r2)),
	 "\[Delta]e"->RealSign[x]*((1-e^2)((1-e)\[Delta]r1-\[Delta]r2-e*\[Delta]r2))/(2p),
	 (*Shift frequencies*)
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 "MinoFrequenciesCorrection"->RealSign[x]*{\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]\[Phi]s},
	 "BLFrequenciesCorrection"->RealSign[x]*{\[CapitalUpsilon]rs/\[CapitalUpsilon]tg-\[CapitalUpsilon]rg/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts,\[CapitalUpsilon]\[Phi]s/\[CapitalUpsilon]tg-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts},
	 "MinoPrecessionFrequency"->\[CapitalUpsilon]p,
	 "BLPrecessionFrequency"->\[CapitalUpsilon]p/\[CapitalUpsilon]tg,
	 (*Keys purely oscillatory part geodesic coordinate time and azimuthal trajectory*)
	 "\[CapitalDelta]trg"->Function[{wr},tgICr1gfun[wr,a,p,e,x]],
	 "\[CapitalDelta]\[Phi]rg"->Function[{wr},\[Phi]gICr1gfun[wr,a,p,e,x]],
	 (*Keys geodesic radial trajectory*)
	 "rg"->Function[{wr},rgICr1gfun[wr,a,p,e,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{wr},VtrgICr1gfun[wr,a,p,e,x]],
	 "vrg"->Function[{wr},drgd\[Lambda]fun[wr,a,p,e,x]],
	 "vrgOfr"->Function[{r},drgd\[Lambda]funOfr[r,a,p,e,x]],
	 "v\[Phi]g"->Function[{wr},V\[Phi]rgICr1gfun[wr,a,p,e,x]],
	 (*Keys corrections trajectory*)
	 "\[CapitalDelta]\[Delta]tpar"->Function[{wr},RealSign[x]*\[Delta]tfunFFPerPar[wr,a,p,e,x,{\[Delta]r1,\[Delta]r2}]],
	 "\[Delta]tOfrpar"->Function[{r},RealSign[x]*\[Delta]tfunOfrFFPerPar[r,a,p,e,x,{\[Delta]r1,\[Delta]r2}]],
	 "\[Delta]rpar"->Function[{wr},RealSign[x]*\[Delta]rfunFFPerPar[wr,a,p,e,x,{\[Delta]r1,\[Delta]r2}]],
	 "\[Delta]rOfrpar"->Function[{r},RealSign[x]*\[Delta]rfunOfrFFPerPar[r,a,p,e,x,{\[Delta]r1,\[Delta]r2}]],
	 "\[Psi]p"->Function[{wr},RealSign[x]*\[Psi]pICr1g[wr,a,p,e,x]],
	 "\[Psi]pOfr"->Function[{r},RealSign[x]*\[Psi]pOfrICr1g[r,a,p,e,x]],
	 "\[Delta]zort"->Function[{wr},\[Delta]zfunICr1gPer[wr,a,p,e,x]],
	 "\[Delta]zOfrort"->Function[{r},\[Delta]zfunOfrICr1gPer[r,a,p,e,x]],
	 "\[CapitalDelta]\[Delta]\[Phi]par"->Function[{wr},RealSign[x]*\[Delta]\[Phi]funFFPerPar[wr,a,p,e,x,{\[Delta]r1,\[Delta]r2}]],
	 "\[Delta]\[Phi]Ofrpar"->Function[{r},RealSign[x]*\[Delta]\[Phi]funOfrFFPerPar[r,a,p,e,x,{\[Delta]r1,\[Delta]r2}]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{wr},RealSign[x]*\[Delta]vtfunFFPerPar[wr,a,p,e,x,{\[Delta]r1,\[Delta]r2}]],
	 "\[Delta]vtOfrpar"->Function[{r},RealSign[x]*\[Delta]vtfunOfrFFPerPar[r,a,p,e,x,{\[Delta]r1,\[Delta]r2}]],
	 "\[Delta]vrpar"->Function[{wr},RealSign[x]*\[Delta]vrfunFFPerPar[wr,a,p,e,x,{\[Delta]r1,\[Delta]r2}]],
	 "\[Delta]vrOfrpar"->Function[{r},RealSign[x]*\[Delta]vrfunOfrFFPerPar[r,a,p,e,x,{\[Delta]r1,\[Delta]r2}]],
	 "\[Delta]v\[Phi]par"->Function[{wr},RealSign[x]*\[Delta]v\[Phi]funFFPerPar[wr,a,p,e,x,{\[Delta]r1,\[Delta]r2}]],
	 "\[Delta]v\[Phi]Ofrpar"->Function[{r},RealSign[x]*\[Delta]v\[Phi]funOfrFFPerPar[r,a,p,e,x,{\[Delta]r1,\[Delta]r2}]]
	|>
]


(* ::Section::Closed:: *)
(*Near equatorial orbits - fixed eccentricity*)


(* ::Subsection::Closed:: *)
(*Analytic solutions*)


(* ::Subsubsection::Closed:: *)
(*Bound solutions*)


KerrNearEqSpinOrbitCorrFEPer[a_, p_, e_, x_]:=Module[{EEg,Lzg,\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g,\[CapitalUpsilon]p,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]zs,\[CapitalUpsilon]\[Phi]s},
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];	
	(* geodesic frequencies *)
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,x];	
	(* spin correction constants of motion *)
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];	  
	(* spin correction radial and polar frequencies *)
	\[CapitalUpsilon]rs=\[Delta]\[CapitalUpsilon]rfunFE[a,p,e,x];
	\[CapitalUpsilon]ts=\[Delta]\[CapitalUpsilon]tfunFE[a,p,e,x];
	\[CapitalUpsilon]\[Phi]s=\[Delta]\[CapitalUpsilon]\[Phi]funFE[a,p,e,x];
	
	(*The RealSign serves to match the sign convetion of generic orbits solutions in the limit z->0*)				
	<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->RealSign[x]*\[Delta]EE,
	 "Jzs"->RealSign[x]*\[Delta]Lz,
	 (*Shift orbital elements*)
	 "\[Delta]p"->RealSign[x]*\[Delta]pfunFE[a,p,e,x],
	 "\[Delta]e"->0,
	 (*Shift frequencies*)
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 "MinoFrequenciesCorrection"->RealSign[x]*{\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]\[Phi]s},
	 "BLFrequenciesCorrection"->RealSign[x]*{\[CapitalUpsilon]rs/\[CapitalUpsilon]tg-\[CapitalUpsilon]rg/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts,\[CapitalUpsilon]\[Phi]s/\[CapitalUpsilon]tg-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts},
	 "MinoPrecessionFrequency"->\[CapitalUpsilon]p,
	 "BLPrecessionFrequency"->\[CapitalUpsilon]p/\[CapitalUpsilon]tg,
	 (*Keys purely oscillatory part geodesic coordinate time and azimuthal trajectory*)
	 "\[CapitalDelta]trg"->Function[{wr},tgICr1gfun[wr,a,p,e,x]],
	 "\[CapitalDelta]\[Phi]rg"->Function[{wr},\[Phi]gICr1gfun[wr,a,p,e,x]],
	 (*Keys geodesic radial trajectory*)
	 "rg"->Function[{wr},rgICr1gfun[wr,a,p,e,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{wr},VtrgICr1gfun[wr,a,p,e,x]],
	 "vrg"->Function[{wr},drgd\[Lambda]fun[wr,a,p,e,x]],
	 "vrgOfr"->Function[{r},drgd\[Lambda]funOfr[r,a,p,e,x]],
	 "v\[Phi]g"->Function[{wr},V\[Phi]rgICr1gfun[wr,a,p,e,x]],
	 (*Keys corrections trajectory*)
	 "\[CapitalDelta]\[Delta]tpar"->Function[{wr},RealSign[x]*\[Delta]tfunFEPerPar[wr,a,p,e,x]],
	 "\[Delta]tOfrpar"->Function[{r},RealSign[x]*\[Delta]tfunOfrFEPerPar[r,a,p,e,x]],
	 "\[Delta]rpar"->Function[{wr},RealSign[x]*\[Delta]rfunFEPerPar[wr,a,p,e,x]],
	 "\[Delta]rOfrpar"->Function[{r},RealSign[x]*\[Delta]rfunOfrFEPerPar[r,a,p,e,x]],
	 "\[Psi]p"->Function[{wr},RealSign[x]*\[Psi]pICr1g[wr,a,p,e,x]],
	 "\[Psi]pOfr"->Function[{r},RealSign[x]*\[Psi]pOfrICr1g[r,a,p,e,x]],
	 "\[Delta]zort"->Function[{wr},\[Delta]zfunICr1gPer[wr,a,p,e,x]],
	 "\[Delta]zOfrort"->Function[{r},\[Delta]zfunOfrICr1gPer[r,a,p,e,x]],
	 "\[CapitalDelta]\[Delta]\[Phi]par"->Function[{wr},RealSign[x]*\[Delta]\[Phi]funFEPerPar[wr,a,p,e,x]],
	 "\[Delta]\[Phi]Ofrpar"->Function[{r},RealSign[x]*\[Delta]\[Phi]funOfrFEPerPar[r,a,p,e,x]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{wr},RealSign[x]*\[Delta]vtfunFEPerPar[wr,a,p,e,x]],
	 "\[Delta]vtOfrpar"->Function[{r},RealSign[x]*\[Delta]vtfunOfrFEPerPar[r,a,p,e,x]],
	 "\[Delta]vrpar"->Function[{wr},RealSign[x]*\[Delta]vrfunFEPerPar[wr,a,p,e,x]],
	 "\[Delta]vrOfrpar"->Function[{r},RealSign[x]*\[Delta]vrfunOfrFEPerPar[r,a,p,e,x]],
	 "\[Delta]v\[Phi]par"->Function[{wr},RealSign[x]*\[Delta]v\[Phi]funFEPerPar[wr,a,p,e,x]],
	 "\[Delta]v\[Phi]Ofrpar"->Function[{r},RealSign[x]*\[Delta]v\[Phi]funOfrFEPerPar[r,a,p,e,x]]
	|>
]


(* ::Subsubsection::Closed:: *)
(*Homoclinic solutions*)


KerrNearEqSpinOrbitCorrFEHom[a_, p_, e_, x_]:=Module[{EEg,Lzg,\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g,\[CapitalUpsilon]p,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]zs,\[CapitalUpsilon]\[Phi]s},
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];	
	(* geodesic frequencies *)
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgfunLim[a,p,e,x];
	\[CapitalUpsilon]rg=0;
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]gfunLim[a,p,e,x];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfunLim[a,p,e,x];	
	(* spin correction constants of motion *)
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];	  
	(* spin correction radial and polar frequencies *)
	\[CapitalUpsilon]ts=\[Delta]\[CapitalUpsilon]tfunLimFE[a,p,e,x];
	\[CapitalUpsilon]\[Phi]s=\[Delta]\[CapitalUpsilon]\[Phi]funLimFE[a,p,e,x];
	
	(*The RealSign serves to match the sign convetion of generic orbits solutions in the limit z->0*)
	<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->RealSign[x]*\[Delta]EE,
	 "Jzs"->RealSign[x]*\[Delta]Lz,
	 (*Shift separatrix*)
	 "\[Delta]p"->RealSign[x]*\[Delta]pfunFE[a,p,e,x],
	 "\[Delta]e"->0,
	 (*Shift frequencies*)
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 "MinoFrequenciesCorrection"->RealSign[x]*{\[CapitalUpsilon]ts,\[CapitalUpsilon]\[Phi]s},
	 "BLFrequenciesCorrection"->RealSign[x]*{\[CapitalUpsilon]\[Phi]s/\[CapitalUpsilon]tg-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts},
	 "MinoPrecessionFrequency"->\[CapitalUpsilon]p,
	 "BLPrecessionFrequency"->\[CapitalUpsilon]p/\[CapitalUpsilon]tg,
	 "\[Delta]\[CapitalUpsilon]rhat"->\[Delta]\[CapitalUpsilon]rover\[CapitalUpsilon]rgfunLimFE[a,p,e,x],
	 (*Keys geodesic trajectories*)
	 "trg"->Function[{r},tgfunHom[r,a,p,e,x]],
	 "\[Phi]rg"->Function[{r},\[Phi]gfunHom[r,a,p,e,x]],
	 "rg"->Function[{\[Lambda]},rgfunHom[\[Lambda],a,p,e,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{r},Vtrgfun[r,a,p,e,x]],
	 "vrg"->Function[{r},drgd\[Lambda]funHom[r,a,p,e,x]],
	 "v\[Phi]g"->Function[{r},V\[Phi]rgfun[r,a,p,e,x]],
	 (*Keys corrections trajectory*)
	 "\[Delta]tpar"->Function[{r},RealSign[x]*\[Delta]tfunFEHomPar[r,a,p,e,x]],
	 "\[Delta]rpar"->Function[{r},RealSign[x]*\[Delta]rfunFEHomPar[r,a,p,e,x]],
	 "\[Psi]p"->Function[{r},\[Psi]pHom[r,a,p,e,x]],
	 "\[Delta]zort"->Function[{r},\[Delta]zfunHom[r,a,p,e,x]],
	 "\[Delta]\[Phi]par"->Function[{r},RealSign[x]*\[Delta]\[Phi]funFEHomPar[r,a,p,e,x]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{r},RealSign[x]*\[Delta]vtfunFEHomPar[r,a,p,e,x]],
	 "\[Delta]vrpar"->Function[{r},RealSign[x]*\[Delta]vrfunFEHomPar[r,a,p,e,x]],
	 "\[Delta]v\[Phi]par"->Function[{r},RealSign[x]*\[Delta]v\[Phi]funFEHomPar[r,a,p,e,x]]
	|>
]


(* ::Subsubsection::Closed:: *)
(*ISCO plunge solutions*)


KerrNearEqSpinOrbitCorrISCOPlunge[a_, x_]:=Module[{r1g,EEg,Lzg,\[Delta]EE,\[Delta]Lz},
	r1g=ISCOradius[a,x];
	(* geodesic constants of motion *)
	EEg=EEgfun[a,r1g,0,x];
	Lzg=Lzgfun[a,r1g,0,x];	
	(* spin correction constants of motion *)
	\[Delta]EE=\[Delta]EEISCOfun[a,x];
	\[Delta]Lz=\[Delta]LzISCOfun[a,x];
	
	(*The RealSign serves to match the sign convetion of generic orbits solutions in the limit z->0*)    
	<|
	 "OrbitalElements"->{a,r1g,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->RealSign[x]*\[Delta]EE,
	 "Jzs"->RealSign[x]*\[Delta]Lz,
	 (*Shift separatrix*)
	 "\[Delta]p"->RealSign[x]*\[Delta]pfunFE[a,r1g,0,x],
	 "\[Delta]e"->0,
	 (*Shift frequencies*)
	 "BLFrequenciesGeo"->Function[{r},\[CapitalOmega]\[Phi]gfunISCOplunge[r,a,x]],
	 "BLFrequenciesCorrection"->Function[{r},\[Delta]\[CapitalOmega]\[Phi]funISCOplunge[r,a,x]],
	 "BLPrecessionFrequency"->Function[{r},\[CapitalOmega]pfunISCOplunge[r,a,x]],
	 "\[Delta]\[CapitalUpsilon]rhat"->RealSign[x]*\[Delta]\[CapitalUpsilon]rover\[CapitalUpsilon]rgfunLimFE[a,r1g,0,x],
	 (*Keys geodesic trajectories*)
	 "trg"->Function[{r},tgfunISCOplunge[r,a,x]],
	 "\[Phi]rg"->Function[{r},\[Phi]gfunISCOplunge[r,a,x]],
	 "rg"->Function[{\[Lambda]},rgfunISCOplunge[\[Lambda],a,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{r},VtrgfunISCOplunge[r,a,x]],
	 "vrg"->Function[{r},drgd\[Lambda]funISCOplunge[r,a,x]],
	 "v\[Phi]g"->Function[{r},V\[Phi]rgfunISCOplunge[r,a,x]],
	 (*Keys corrections trajectories*)
	 "\[Delta]tpar"->Function[{r},RealSign[x]*\[Delta]tfunISCOplungePar[r,a,x]],
	 "\[Delta]rpar"->Function[{r},RealSign[x]*\[Delta]rfunISCOplungePar[r,a,x]],
	 "\[Psi]p"->Function[{r},\[Psi]pISCOplunge[r,a,x]],
	 "\[Delta]zort"->Function[{r},\[Delta]zfunISCOplunge[r,a,x]],
	 "\[Delta]\[Phi]par"->Function[{r},RealSign[x]*\[Delta]\[Phi]funISCOplungePar[r,a,x]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{r},RealSign[x]*\[Delta]vtfunISCOplungePar[r,a,x]],
	 "\[Delta]vrpar"->Function[{r},RealSign[x]*\[Delta]vrfunISCOplungePar[r,a,x]],
	 "\[Delta]v\[Phi]par"->Function[{r},RealSign[x]*\[Delta]v\[Phi]funISCOplungePar[r,a,x]],
	 (*Shift ISCO*)
	 "\[Delta]ISCO"->RealSign[x]*\[Delta]rISCOfun[a,x]
	|>
]


(* ::Subsubsection::Closed:: *)
(*Critical plunge solutions*)


KerrNearEqSpinOrbitCorrFECritPlunge[a_, p_, e_, x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz},
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];	
	(* spin correction constants of motion *)
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	
	(*The RealSign serves to match the sign convetion of generic orbits solutions in the limit z->0*)		
	<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->RealSign[x]*\[Delta]EE,
	 "Jzs"->RealSign[x]*\[Delta]Lz,
	 (*Shift separatrix*)
	 "\[Delta]p"->RealSign[x]*\[Delta]pfunFE[a,p,e,x],
	 "\[Delta]e"->0,
	 (*Shift frequencies*)
	 "BLFrequenciesGeo"->Function[{r},\[CapitalOmega]\[Phi]gfun[r,a,p,e,x]],
	 "BLFrequenciesCorrection"->Function[{r},RealSign[x]*\[Delta]\[CapitalOmega]\[Phi]funFECritplunge[r,a,p,e,x]],
	 "BLPrecessionFrequency"->Function[{r},\[CapitalOmega]pfun[r,a,p,e,x]],
	 "\[Delta]\[CapitalUpsilon]rhat"->RealSign[x]*\[Delta]\[CapitalUpsilon]rover\[CapitalUpsilon]rgfunLimFE[a,p,e,x],
	 (*Keys purely oscillatory part geodesic coordinate time and azimuthal trajectory*)
	 "trg"->Function[{r},tgfunCritplunge[r,a,p,e,x]],
	 "\[Phi]rg"->Function[{r},\[Phi]gfunCritplunge[r,a,p,e,x]],
	 (*Keys geodesic radial trajectory*)
	 "rg"->Function[{\[Lambda]},rgfunCritplunge[\[Lambda],a,p,e,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{r},Vtrgfun[r,a,p,e,x]],
	 "vrg"->Function[{r},drgd\[Lambda]funHom[r,a,p,e,x]],
	 "v\[Phi]g"->Function[{r},V\[Phi]rgfun[r,a,p,e,x]],
	 (*Keys corrections trajectory*)
	 "\[Delta]tpar"->Function[{r},RealSign[x]*\[Delta]tfunFECritplungePar[r,a,p,e,x]],
	 "\[Delta]rpar"->Function[{r},RealSign[x]*\[Delta]rfunFECritplungePar[r,a,p,e,x]],
	 "\[Psi]p"->Function[{r},\[Psi]pCritplunge[r,a,p,e,x]],
	 "\[Delta]zort"->Function[{r},\[Delta]zfunCritplunge[r,a,p,e,x]],
	 "\[Delta]\[Phi]par"->Function[{r},RealSign[x]*\[Delta]\[Phi]funFECritplungePar[r,a,p,e,x]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{r},RealSign[x]*\[Delta]vtfunFECritplungePar[r,a,p,e,x]],
	 "\[Delta]vrpar"->Function[{r},RealSign[x]*\[Delta]vrfunFECritplungePar[r,a,p,e,x]],
	 "\[Delta]v\[Phi]par"->Function[{r},RealSign[x]*\[Delta]v\[Phi]funFECritplungePar[r,a,p,e,x]]
	|>
]


(* ::Subsubsection::Closed:: *)
(*Generic plunge solutions*)


KerrNearEqSpinOrbitCorrFEPlunge[a_, p_, e_, x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,\[Rho]rg,\[Rho]ig,\[Delta]p,\[Delta]r1,\[Delta]\[Rho]r,Ar1g,B0},
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];
	Ar1g=Sqrt[\[Rho]ig^2+(r1g-\[Rho]rg)^2];
	B0=Sqrt[\[Rho]ig^2+\[Rho]rg^2];
	
	\[Delta]p=\[Delta]pfunFEcc[a,p,e,x];
	\[Delta]r1=\[Delta]p/(1-e);
	\[Delta]\[Rho]r=\[Delta]p/(1+e);
	
	(* geodesic constants of motion *)
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];	
	(* spin correction constants of motion *)
	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x]+dEEdpfuncc[p,e]\[Delta]pfunFEcc[a,p,e,x];
	\[Delta]Lz=\[Delta]Lzfuncc[a,p,e,x]+dLzdpfuncc[a,p,e,x]\[Delta]pfunFEcc[a,p,e,x];
	
	(*The RealSign serves to match the sign convetion of generic orbits solutions in the limit z->0*)		
	<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->RealSign[x]*\[Delta]EE,
	 "Jzs"->RealSign[x]*\[Delta]Lz,
	 (*Shift orbital elements*)
	 "\[Delta]p"->RealSign[x]*\[Delta]pfunFEcc[a,p,e,x],
	 "\[Delta]e"->0,
	 (*Shift frequencies*)
	 "BLFrequenciesGeo"->Function[{r},\[CapitalOmega]\[Phi]gfunPlunge[r,a,p,e,x]],
	 "BLFrequenciesCorrection"->Function[{r},RealSign[x]*\[Delta]\[CapitalOmega]\[Phi]funFEPlunge[r,a,p,e,x]],
	 "BLPrecessionFrequency"->Function[{r},\[CapitalOmega]pfunPlunge[r,a,p,e,x]],
	 "\[Delta]\[CapitalUpsilon]rhat"->RealSign[x]*1/2 ((1+B0/Ar1g) \[Delta]r1/r1g-(2(Ar1g+B0)\[Delta]r41funcc[a,p,e,x])/(B0*r1g)+(\[Delta]\[Rho]r(r1g-2\[Rho]rg))/(Ar1g*B0)+(a(2B0^2(Ar1g+B0)-B0*r1g^2+2(2Ar1g+B0)r1g*\[Rho]rg))/(3B0^3*r1g^2)-(2EEg*\[Delta]EE)/(1-EEg^2)),
	 (*Keys geodesic trajectories*)
	 "trg"->Function[{r},tgfunPlunge[r,a,p,e,x]],
	 "\[Phi]rg"->Function[{r},\[Phi]gfunPlunge[r,a,p,e,x]],
	 "rg"->Function[{\[Lambda]},rgfunPlunge[\[Lambda],a,p,e,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{r},VtrgfunPlunge[r,a,p,e,x]],
	 "vrg"->Function[{r},drgd\[Lambda]funPlunge[r,a,p,e,x]],
	 "v\[Phi]g"->Function[{r},V\[Phi]rgfunPlunge[r,a,p,e,x]],
	 (*Keys corrections trajectory*)
	 "\[Delta]tpar"->Function[{r},RealSign[x]*\[Delta]tfunFEPlungePar[r,a,p,e,x]],
	 "\[Delta]rpar"->Function[{r},RealSign[x]*\[Delta]rfunFEPlungePar[r,a,p,e,x]],
	 "\[Psi]p"->Function[r,\[Psi]pPlunge[r,a,p,e,x]],
	 "\[Delta]zort"->Function[{r},\[Delta]zfunPlunge[r,a,p,e,x]],
	 "\[Delta]\[Phi]par"->Function[{r},RealSign[x]*\[Delta]\[Phi]funFEPlungePar[r,a,p,e,x]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{r},RealSign[x]*\[Delta]vtfunFEPlungePar[r,a,p,e,x]],
	 "\[Delta]vrpar"->Function[{r},RealSign[x]*\[Delta]vrfunFEPlungePar[r,a,p,e,x]],
	 "\[Delta]v\[Phi]par"->Function[{r},RealSign[x]*\[Delta]v\[Phi]funFEPlungePar[r,a,p,e,x]]
	|>
]


(* ::Subsubsection::Closed:: *)
(*Plunge related to bound orbits*)


KerrNearEqSpinOrbitCorrFEPlungeBoundOrbit[a_, p_, e_, x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz},
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	(* spin correction constants of motion *)
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunFE[a,p,e,x];
	
	(*The RealSign serves to match the sign convetion of generic orbits solutions in the limit z->0*)				
	<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->RealSign[x]*\[Delta]EE,
	 "Jzs"->RealSign[x]*\[Delta]Lz,
	 (*Shift separatrix*)
	 "\[Delta]p"->RealSign[x]*\[Delta]pfunFE[a,p,e,x],
	 "\[Delta]e"->0,
	 (*Shift frequencies*)
	 "BLFrequenciesGeo"->Function[{r},\[CapitalOmega]\[Phi]gfun[r,a,p,e,x]],
	 "BLFrequenciesCorrection"->Function[{r},RealSign[x]*\[Delta]\[CapitalOmega]\[Phi]funFEPlungeBoundOrbit[r,a,p,e,x]],
	 "BLPrecessionFrequency"->Function[{r},\[CapitalOmega]pfun[r,a,p,e,x]],
	 "\[Delta]\[CapitalUpsilon]rhat"->RealSign[x]*\[Delta]\[CapitalUpsilon]rfunFE[a,p,e,x]/\[CapitalUpsilon]rgfun[a,p,e,x],
	 (*Keys geodesic trajectories*)
	 "trg"->Function[{r},tgfunPlungeBoundOrbit[r,a,p,e,x]],
	 "\[Phi]rg"->Function[{r},\[Phi]gfunPlungeBoundOrbit[r,a,p,e,x]], 
	 "rg"->Function[{\[Lambda]},rgfunPlungeBoundOrbit[\[Lambda],a,p,e,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{r},Vtrgfun[r,a,p,e,x]],
	 "vrg"->Function[{r},drgd\[Lambda]funPlungeBoundOrbit[r,a,p,e,x]],
	 "v\[Phi]g"->Function[{r},V\[Phi]rgfun[r,a,p,e,x]],
	 (*Keys corrections trajectory*)
	 "\[Delta]tpar"->Function[{r},RealSign[x]*\[Delta]tfunFEPlungeBoundOrbitPar[r,a,p,e,x]],
	 "\[Delta]rpar"->Function[{r},RealSign[x]*\[Delta]rfunFEPlungeBoundOrbitPar[r,a,p,e,x]],
	 "\[Psi]p"->Function[{r},\[Psi]pPlungeBoundOrbit[r,a,p,e,x]],
	 "\[Delta]zort"->Function[{r},\[Delta]zfunPlungeBoundOrbit[r,a,p,e,x]],
	 "\[Delta]\[Phi]par"->Function[{r},RealSign[x]*\[Delta]\[Phi]funFEPlungeBoundOrbitPar[r,a,p,e,x]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{r},RealSign[x]*\[Delta]vtfunFEPlungeBoundOrbitPar[r,a,p,e,x]],
	 "\[Delta]vrpar"->Function[{r},RealSign[x]*\[Delta]vrfunFEPlungeBoundOrbitPar[r,a,p,e,x]],
	 "\[Delta]v\[Phi]par"->Function[{r},RealSign[x]*\[Delta]v\[Phi]funFEPlungeBoundOrbitPar[r,a,p,e,x]]
	|>
]


(* ::Subsubsection::Closed:: *)
(*Plunge related to bound orbits - double root*)


KerrNearEqSpinOrbitCorrFEPlungeBoundOrbitDoubleRoot[a_, p_, x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz},
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,0,x];
	Lzg=Lzgfun[a,p,0,x];	
	(* spin correction constants of motion *)
	\[Delta]EE=\[Delta]EEfunFT[a,p,0,x]+dEEdpfun[a,p,0,x]\[Delta]pfunFE[a,p,0,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,0,x]+dLzdpfun[a,p,0,x]\[Delta]pfunFE[a,p,0,x];
	
	(*The RealSign serves to match the sign convetion of generic orbits solutions in the limit z->0*)		
	<|
	 "OrbitalElements"->{a,p,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->RealSign[x]*\[Delta]EE,
	 "Jzs"->RealSign[x]*\[Delta]Lz,
	 (*Shift separatrix*)
	 "\[Delta]p"->RealSign[x]*\[Delta]pfunFE[a,p,0,x],
	 "\[Delta]e"->0,
	 (*Shift frequencies*)
	 "BLFrequenciesGeo"->Function[{r},\[CapitalOmega]\[Phi]gfun[r,a,p,0,x]],
	 "BLFrequenciesCorrection"->Function[{r},RealSign[x]*\[Delta]\[CapitalOmega]\[Phi]funFEPlungeBoundOrbitDoubleRoot[r,a,p,x]],
	 "BLPrecessionFrequency"->Function[{r},\[CapitalOmega]pfun[r,a,p,0,x]],
	 (*Keys geodesic trajectories*)
	 "trg"->Function[{r},tgfunPlungeBoundOrbitDoubleRoot[r,a,p,x]],
	 "\[Phi]rg"->Function[{r},\[Phi]gfunPlungeBoundOrbitDoubleRoot[r,a,p,x]],
	 "rg"->Function[{\[Lambda]},rgfunPlungeBoundOrbitDoubleRoot[\[Lambda],a,p,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{r},Vtrgfun[r,a,0,p,x]],
	 "vrg"->Function[{r},drgd\[Lambda]funPlungeBoundOrbit[r,a,0,p,x]],
	 "v\[Phi]g"->Function[{r},V\[Phi]rgfun[r,a,0,p,x]],
	 (*Keys corrections trajectory*)
	 "\[Delta]tpar"->Function[{r},RealSign[x]*\[Delta]tfunFEPlungeBoundOrbitDoubleRootPar[r,a,p,x]],
	 "\[Delta]rpar"->Function[{r},RealSign[x]*\[Delta]rfunFEPlungeBoundOrbitDoubleRootPar[r,a,p,x]],
	 "\[Psi]p"->Function[{r},\[Psi]pPlungeBoundOrbitDoubleRoot[r,a,p,x]],
	 "\[Delta]zort"->Function[{r},\[Delta]zfunPlungeBoundOrbitDoubleRoot[r,a,p,x]],
	 "\[Delta]\[Phi]par"->Function[{r},RealSign[x]*\[Delta]\[Phi]funFEPlungeBoundOrbitDoubleRootPar[r,a,p,x]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{r},RealSign[x]*\[Delta]vtfunFEPlungeBoundOrbitDoubleRootPar[r,a,p,x]],
	 "\[Delta]vrpar"->Function[{r},RealSign[x]*\[Delta]vrfunFEPlungeBoundOrbitDoubleRootPar[r,a,p,x]],
	 "\[Delta]v\[Phi]par"->Function[{r},RealSign[x]*\[Delta]v\[Phi]funFEPlungeBoundOrbitDoubleRootPar[r,a,p,x]]
	|>
]


(* ::Section::Closed:: *)
(*Near equatorial orbits - fixed eccentricity, non-linearized solutions*)


(* ::Subsection::Closed:: *)
(*Analytic solutions*)


(* ::Subsubsection::Closed:: *)
(*ISCO plunge solutions*)


KerrEqISCOPlungePartialExp[a_, S_, x_]:=Module[{r1g,\[Delta]rISCO,r1},
	r1g=ISCOradius[a,x];
	\[Delta]rISCO=\[Delta]rISCOfun[a,x];
	r1=r1g+S*RealSign[x]*\[Delta]rISCO;
		
	(*The RealSign serves to match the sign convetion of generic orbits solutions in the limit z->0*)    
	<|
	 (*Keys trajectories*)
	 "tr"->Function[{r},tfunISCOplungePartialExp[r,a,S,x]],
	 "\[Phi]r"->Function[{r},\[Phi]funISCOplungePartialExp[r,a,S,x]],
	 (*Keys velocities*)
	 "dtdr"->Function[{r},dtdrfunISCOplungeOPartialExp[r,a,S,r1g,0,x]],
	 "vr"->Function[{r},drd\[Lambda]funISCOplungePartialExp[r,a,S,x]],
	 "d\[Phi]dr"->Function[{r},d\[Phi]drfunISCOplungePartialExp[r,a,S,r1g,0,x]],
	 "rISCO"->r1
	|>
]


(* ::Subsubsection::Closed:: *)
(*Critical plunge solutions*)


KerrEqFECritplungePartialExp[a_, S_, p_, e_, x_]:=Module[{\[Delta]p,r1g,r2g,r1,r2},
	r1g=p/(1-e);
	r2g=p/(1+e);
	\[Delta]p=\[Delta]pfunFE[a,p,e,x];
	r1=r1g+S*RealSign[x]*\[Delta]p/(1-e);
	r2=r2g+S*RealSign[x]*\[Delta]p/(1+e);
	
	(*The RealSign serves to match the sign convetion of generic orbits solutions in the limit z->0*)    
	<|
	 (*Keys trajectories*)
	 "tr"->Function[{r},tfunCritplungePartialExp[r,a,S,p,e,x]],
	 "\[Phi]r"->Function[{r},\[Phi]funCritplungePartialExp[r,a,S,p,e,x]],
	 (*Keys velocities*)
	 "dtdr"->Function[{r},dtdrfunCritplungePartialExp[r,a,S,p,e,x]],
	 "vr"->Function[{r},drd\[Lambda]funHomPartialExp[r,a,S,p,e,x]],
	 "d\[Phi]dr"->Function[{r},d\[Phi]drfunCritplungePartialExp[r,a,S,p,e,x]],
	 "r1"->r1,
	 "r2"->r2
	|>
]


(* ::Subsubsection::Closed:: *)
(*Generic plunge solutions*)


KerrEqFEPlungePartialExp[a_, S_, p_, e_, x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,\[Delta]p,r1},
	r1g=p/(1-e);
	\[Delta]p=\[Delta]pfunFEcc[a,p,e,x];
	r1=r1g+S*RealSign[x]*\[Delta]p/(1-e);
		
	(*The RealSign serves to match the sign convetion of generic orbits solutions in the limit z->0*)    
	<|
	 (*Keys trajectories*)
	 "tr"->Function[{r},tfunPlungePartialExp[r,a,S,p,e,x]],
	 "\[Phi]r"->Function[{r},\[Phi]funPlungePartialExp[r,a,S,p,e,x]],
	 (*Keys velocities*)
	 "dtdr"->Function[{r},dtdrfunPlungePartialExp[r,a,S,p,e,x]],
	 "vr"->Function[{r},drd\[Lambda]funPlungePartialExp[r,a,S,p,e,x]],
	 "d\[Phi]dr"->Function[{r},d\[Phi]drfunPlungePartialExp[r,a,S,p,e,x]],
	 "r1"->r1
	|>
]


(* ::Subsubsection::Closed:: *)
(*Plunge related to bound orbits*)


KerrEqFEPlungeBoundOrbitPartialExp[a_, S_, p_, e_,x_]:=Module[{EEg,r1g,r2g,\[Delta]p,r1,r2,r3},
	r1g=p/(1-e);
	r2g=p/(1+e);
	
	EEg=EEgfun[a,p,e,x];
	\[Delta]p=\[Delta]pfunFE[a,p,e,x];
	r1=r1g+S*RealSign[x]*\[Delta]p/(1-e);
	r2=r2g+S*RealSign[x]*\[Delta]p/(1+e);
	r3=2/(1-EEg^2)-(r1g+r2g)+S*RealSign[x]*\[Delta]p/(1+e);
	
	(*The RealSign serves to match the sign convetion of generic orbits solutions in the limit z->0*)    
	<|
	 (*Keys trajectories*)
	 "tr"->Function[{r},tfunPlungeBoundOrbitPartialExp[r,a,S,p,e,x]],
	 "\[Phi]r"->Function[{r},\[Phi]funPlungeBoundOrbitPartialExp[r,a,S,p,e,x]],
	 (*Keys velocities*)
	 "dtdr"->Function[{r},dtdrfunPlungeBoundOrbitPartialExp[r,a,S,p,e,x]],
	 "vr"->Function[{r},drd\[Lambda]funPlungeBoundOrbitPartialExp[r,a,S,p,e,x]],
	 "d\[Phi]dr"->Function[{r},d\[Phi]drfunPlungeBoundOrbitPartialExp[r,a,S,p,e,x]],
	 "r1"->r1,
	 "r2"->r2,
	 "r3"->r3
	|>
]


(* ::Section::Closed:: *)
(*Near equatorial orbits - regular gauge*)


(* ::Subsection::Closed:: *)
(*Analytic solutions*)


(* ::Subsubsection::Closed:: *)
(*Bound solutions*)


KerrNearEqSpinOrbitCorrRPer[a_, p_, e_, x_]:=Module[{EEg,Lzg,\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g,\[CapitalUpsilon]p,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]zs,\[CapitalUpsilon]\[Phi]s},
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];	
	(* geodesic frequencies *)
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,x];	
	(* spin correction constants of motion *)
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dLzdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];	  
	(* spin correction radial and polar frequencies *)
	\[CapitalUpsilon]rs=\[Delta]\[CapitalUpsilon]rfunR[a,p,e,x];
	\[CapitalUpsilon]ts=\[Delta]\[CapitalUpsilon]tfunR[a,p,e,x];
	\[CapitalUpsilon]\[Phi]s=\[Delta]\[CapitalUpsilon]\[Phi]funR[a,p,e,x];
	
	(*The RealSign serves to match the sign convetion of generic orbits solutions in the limit z->0*)				
	<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->RealSign[x]*\[Delta]EE,
	 "Jzs"->RealSign[x]*\[Delta]Lz,
	 (*Shift orbital elements*)
	 "\[Delta]p"->RealSign[x]*\[Delta]pfunRPer[a,p,e,x],
	 "\[Delta]e"->RealSign[x]*\[Delta]efunRPer[a,p,e,x],
	 (*Shift frequencies*)
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 "MinoFrequenciesCorrection"->RealSign[x]*{\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]\[Phi]s},
	 "BLFrequenciesCorrection"->RealSign[x]*{\[CapitalUpsilon]rs/\[CapitalUpsilon]tg-\[CapitalUpsilon]rg/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts,\[CapitalUpsilon]\[Phi]s/\[CapitalUpsilon]tg-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts},
	 "MinoPrecessionFrequency"->\[CapitalUpsilon]p,
	 "BLPrecessionFrequency"->\[CapitalUpsilon]p/\[CapitalUpsilon]tg,
	 (*Keys purely oscillatory part geodesic coordinate time and azimuthal trajectory*)
	 "\[CapitalDelta]trg"->Function[{wr},tgICr1gfun[wr,a,p,e,x]],
	 "\[CapitalDelta]\[Phi]rg"->Function[{wr},\[Phi]gICr1gfun[wr,a,p,e,x]],
	 (*Keys geodesic radial trajectory*)
	 "rg"->Function[{wr},rgICr1gfun[wr,a,p,e,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{wr},VtrgICr1gfun[wr,a,p,e,x]],
	 "vrg"->Function[{wr},drgd\[Lambda]fun[wr,a,p,e,x]],
	 "vrgOfr"->Function[{r},drgd\[Lambda]funOfr[r,a,p,e,x]],
	 "v\[Phi]g"->Function[{wr},V\[Phi]rgICr1gfun[wr,a,p,e,x]],
	 (*Keys corrections trajectory*)
	 "\[CapitalDelta]\[Delta]tpar"->Function[{wr},RealSign[x]*\[Delta]tfunRPerPar[wr,a,p,e,x]],
	 "\[Delta]tOfrpar"->Function[{r},RealSign[x]*\[Delta]tfunOfrRPerPar[r,a,p,e,x]],
	 "\[Delta]rpar"->Function[{wr},RealSign[x]*\[Delta]rfunRPerPar[wr,a,p,e,x]],
	 "\[Delta]rOfrpar"->Function[{r},RealSign[x]*\[Delta]rfunOfrRPerPar[r,a,p,e,x]],
	 "\[Psi]p"->Function[{wr},RealSign[x]*\[Psi]pICr1g[wr,a,p,e,x]],
	 "\[Psi]pOfr"->Function[{r},RealSign[x]*\[Psi]pOfrICr1g[r,a,p,e,x]],
	 "\[Delta]zort"->Function[{wr},\[Delta]zfunICr1gPer[wr,a,p,e,x]],
	 "\[Delta]zOfrort"->Function[{r},\[Delta]zfunOfrICr1gPer[r,a,p,e,x]],
	 "\[CapitalDelta]\[Delta]\[Phi]par"->Function[{wr},RealSign[x]*\[Delta]\[Phi]funRPerPar[wr,a,p,e,x]],
	 "\[Delta]\[Phi]Ofrpar"->Function[{r},RealSign[x]*\[Delta]\[Phi]funOfrRPerPar[r,a,p,e,x]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{wr},RealSign[x]*\[Delta]vtfunRPerPar[wr,a,p,e,x]],
	 "\[Delta]vtOfrpar"->Function[{r},RealSign[x]*\[Delta]vtfunOfrRPerPar[r,a,p,e,x]],
	 "\[Delta]vrpar"->Function[{wr},RealSign[x]*\[Delta]vrfunRPerPar[wr,a,p,e,x]],
	 "\[Delta]vrOfrpar"->Function[{r},RealSign[x]*\[Delta]vrfunOfrRPerPar[r,a,p,e,x]],
	 "\[Delta]v\[Phi]par"->Function[{wr},RealSign[x]*\[Delta]v\[Phi]funRPerPar[wr,a,p,e,x]],
	 "\[Delta]v\[Phi]Ofrpar"->Function[{r},RealSign[x]*\[Delta]v\[Phi]funOfrRPerPar[r,a,p,e,x]]
	|>
]


(* ::Subsubsection::Closed:: *)
(*Homoclinic solutions*)


KerrNearEqSpinOrbitCorrRHom[a_, p_, e_, x_]:=Module[{EEg,Lzg,\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g,\[CapitalUpsilon]p,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]zs,\[CapitalUpsilon]\[Phi]s},
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];	
	(* geodesic frequencies *)
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgfunLim[a,p,e,x];
	\[CapitalUpsilon]rg=0;
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]gfunLim[a,p,e,x];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfunLim[a,p,e,x];	
	(* spin correction constants of motion *)
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dLzdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];	  
	(* spin correction radial and polar frequencies *)
	\[CapitalUpsilon]ts=\[Delta]\[CapitalUpsilon]tfunLimR[a,p,e,x];
	\[CapitalUpsilon]\[Phi]s=\[Delta]\[CapitalUpsilon]\[Phi]funLimR[a,p,e,x];
	
	(*The RealSign serves to match the sign convetion of generic orbits solutions in the limit z->0*)
	<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->RealSign[x]*\[Delta]EE,
	 "Jzs"->RealSign[x]*\[Delta]Lz,
	 (*Shift separatrix*)
	 "\[Delta]p"->RealSign[x]*\[Delta]pfunRHom[a,p,e,x],
	 "\[Delta]e"->RealSign[x]*\[Delta]efunRHom[a,p,e,x],
	 (*Shift frequencies*)
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 "MinoFrequenciesCorrection"->RealSign[x]*{\[CapitalUpsilon]ts,\[CapitalUpsilon]\[Phi]s},
	 "BLFrequenciesCorrection"->RealSign[x]*{\[CapitalUpsilon]\[Phi]s/\[CapitalUpsilon]tg-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts},
	 "MinoPrecessionFrequency"->\[CapitalUpsilon]p,
	 "BLPrecessionFrequency"->\[CapitalUpsilon]p/\[CapitalUpsilon]tg,
	 (*Keys geodesic trajectories*)
	 "trg"->Function[{r},tgfunHom[r,a,p,e,x]],
	 "\[Phi]rg"->Function[{r},\[Phi]gfunHom[r,a,p,e,x]],
	 "rg"->Function[{\[Lambda]},rgfunHom[\[Lambda],a,p,e,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{r},Vtrgfun[r,a,p,e,x]],
	 "vrg"->Function[{r},drgd\[Lambda]funHom[r,a,p,e,x]],
	 "v\[Phi]g"->Function[{r},V\[Phi]rgfun[r,a,p,e,x]],
	 (*Keys corrections trajectory*)
	 "\[Delta]tpar"->Function[{r},RealSign[x]*\[Delta]tfunRHomPar[r,a,p,e,x]],
	 "\[Delta]rpar"->Function[{r},RealSign[x]*\[Delta]rfunRHomPar[r,a,p,e,x]],
	 "\[Psi]p"->Function[{r},\[Psi]pHom[r,a,p,e,x]],
	 "\[Delta]zort"->Function[{r},\[Delta]zfunHom[r,a,p,e,x]],
	 "\[Delta]\[Phi]par"->Function[{r},RealSign[x]*\[Delta]\[Phi]funRHomPar[r,a,p,e,x]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{r},RealSign[x]*\[Delta]vtfunRHomPar[r,a,p,e,x]],
	 "\[Delta]vrpar"->Function[{r},RealSign[x]*\[Delta]vrfunRHomPar[r,a,p,e,x]],
	 "\[Delta]v\[Phi]par"->Function[{r},RealSign[x]*\[Delta]v\[Phi]funRHomPar[r,a,p,e,x]]
	|>
]


(* ::Subsubsection::Closed:: *)
(*Critical plunge solutions*)


KerrNearEqSpinOrbitCorrRCritPlunge[a_, p_, e_, x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz},
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];	
	(* spin correction constants of motion *)
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunRHom[a,p,e,x]+dLzdefun[a,p,e,x]\[Delta]efunRHom[a,p,e,x];
	
	(*The RealSign serves to match the sign convetion of generic orbits solutions in the limit z->0*)		
	<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->RealSign[x]*\[Delta]EE,
	 "Jzs"->RealSign[x]*\[Delta]Lz,
	 (*Shift separatrix*)
	 "\[Delta]p"->RealSign[x]*\[Delta]pfunRHom[a,p,e,x],
	 "\[Delta]e"->RealSign[x]*\[Delta]efunRHom[a,p,e,x],
	 (*Shift frequencies*)
	 "BLFrequenciesGeo"->Function[{r},\[CapitalOmega]\[Phi]gfun[r,a,p,e,x]],
	 "BLFrequenciesCorrection"->Function[{r},RealSign[x]*\[Delta]\[CapitalOmega]\[Phi]funRCritplunge[r,a,p,e,x]],
	 "BLPrecessionFrequency"->Function[{r},\[CapitalOmega]pfun[r,a,p,e,x]],
	 "\[Delta]\[CapitalUpsilon]rhat"->RealSign[x]*\[Delta]\[CapitalUpsilon]rover\[CapitalUpsilon]rgfunLimR[a,p,e,x],
	 (*Keys geodesic trajectories*)
	 "trg"->Function[{r},tgfunCritplunge[r,a,p,e,x]],
	 "\[Phi]rg"->Function[{r},\[Phi]gfunCritplunge[r,a,p,e,x]],
	 "rg"->Function[{\[Lambda]},rgfunCritplunge[\[Lambda],a,p,e,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{r},Vtrgfun[r,a,p,e,x]],
	 "vrg"->Function[{r},drgd\[Lambda]funHom[r,a,p,e,x]],
	 "v\[Phi]g"->Function[{r},V\[Phi]rgfun[r,a,p,e,x]],
	 (*Keys corrections trajectory*)
	 "\[Delta]tpar"->Function[{r},RealSign[x]*\[Delta]tfunRCritplungePar[r,a,p,e,x]],
	 "\[Delta]rpar"->Function[{r},RealSign[x]*\[Delta]rfunRCritplungePar[r,a,p,e,x]],
	 "\[Psi]p"->Function[{r},\[Psi]pCritplunge[r,a,p,e,x]],
	 "\[Delta]zort"->Function[{r},\[Delta]zfunCritplunge[r,a,p,e,x]],
	 "\[Delta]\[Phi]par"->Function[{r},RealSign[x]*\[Delta]\[Phi]funRCritplungePar[r,a,p,e,x]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{r},RealSign[x]*\[Delta]vtfunRCritplungePar[r,a,p,e,x]],
	 "\[Delta]vrpar"->Function[{r},RealSign[x]*\[Delta]vrfunRCritplungePar[r,a,p,e,x]],
	 "\[Delta]v\[Phi]par"->Function[{r},RealSign[x]*\[Delta]v\[Phi]funRCritplungePar[r,a,p,e,x]]
	|>
]


(* ::Subsubsection::Closed:: *)
(*Generic plunge solutions*)


KerrNearEqSpinOrbitCorrRPlunge[a_, p_, e_, x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,\[Rho]rg,\[Rho]ig,\[Delta]p,\[Delta]e,\[Delta]r1,\[Delta]\[Rho]r,Ar1g,B0},
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];
	Ar1g=Sqrt[\[Rho]ig^2+(r1g-\[Rho]rg)^2];
	B0=Sqrt[\[Rho]ig^2+\[Rho]rg^2];
	
	\[Delta]p=\[Delta]pfunRPlunge[a,p,e,x];
	\[Delta]e=\[Delta]efunRPlunge[a,p,e,x];
	\[Delta]r1=\[Delta]p/(1-e)+p*\[Delta]e/(1-e)^2;
	\[Delta]\[Rho]r=\[Delta]p/(1+e)-p*\[Delta]e/(1+e)^2;
	
	(* geodesic constants of motion *)
	EEg=EEgfuncc[p,e];
	Lzg=Lzgfuncc[a,p,e,x];	
	(* spin correction constants of motion *)
	\[Delta]EE=\[Delta]EEfuncc[a,p,e,x]+dEEdpfuncc[p,e]\[Delta]p+dEEdefuncc[p,e]\[Delta]e;
	\[Delta]Lz=\[Delta]Lzfuncc[a,p,e,x]+dLzdpfuncc[a,p,e,x]\[Delta]p+dLzdefuncc[a,p,e,x]\[Delta]e;
	
	(*The RealSign serves to match the sign convetion of generic orbits solutions in the limit z->0*)		
	<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->RealSign[x]*\[Delta]EE,
	 "Jzs"->RealSign[x]*\[Delta]Lz,
	 (*Shift orbital elements*)
	 "\[Delta]p"->RealSign[x]*\[Delta]p,
	 "\[Delta]e"->RealSign[x]*\[Delta]e,
	 (*Shift frequencies*)
	 "BLFrequenciesGeo"->Function[{r},\[CapitalOmega]\[Phi]gfunPlunge[r,a,p,e,x]],
	 "BLFrequenciesCorrection"->Function[{r},RealSign[x]*\[Delta]\[CapitalOmega]\[Phi]funRPlunge[r,a,p,e,x]],
	 "BLPrecessionFrequency"->Function[{r},\[CapitalOmega]pfunPlunge[r,a,p,e,x]],
	 "\[Delta]\[CapitalUpsilon]rhat"->RealSign[x]*1/2 ((1+B0/Ar1g) \[Delta]r1/r1g-(2(Ar1g+B0)\[Delta]r41funcc[a,p,e,x])/(B0*r1g)+(\[Delta]\[Rho]r(r1g-2\[Rho]rg))/(Ar1g*B0)+(a(2B0^2(Ar1g+B0)-B0*r1g^2+2(2Ar1g+B0)r1g*\[Rho]rg))/(3B0^3*r1g^2)-(2EEg*\[Delta]EE)/(1-EEg^2)),
	 (*Keys geodesic trajectories*)
	 "trg"->Function[{r},tgfunPlunge[r,a,p,e,x]],
	 "\[Phi]rg"->Function[{r},\[Phi]gfunPlunge[r,a,p,e,x]],
	 "rg"->Function[{\[Lambda]},rgfunPlunge[\[Lambda],a,p,e,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{r},VtrgfunPlunge[r,a,p,e,x]],
	 "vrg"->Function[{r},drgd\[Lambda]funPlunge[r,a,p,e,x]],
	 "v\[Phi]g"->Function[{r},V\[Phi]rgfunPlunge[r,a,p,e,x]],
	 (*Keys corrections trajectory*)
	 "\[Delta]tpar"->Function[{r},RealSign[x]*\[Delta]tfunRPlungePar[r,a,p,e,x]],
	 "\[Delta]rpar"->Function[{r},RealSign[x]*\[Delta]rfunRPlungePar[r,a,p,e,x]],
	 "\[Psi]p"->Function[r,\[Psi]pPlunge[r,a,p,e,x]],
	 "\[Delta]zort"->Function[{r},\[Delta]zfunPlunge[r,a,p,e,x]],
	 "\[Delta]\[Phi]par"->Function[{r},RealSign[x]*\[Delta]\[Phi]funRPlungePar[r,a,p,e,x]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{r},RealSign[x]*\[Delta]vtfunRPlungePar[r,a,p,e,x]],
	 "\[Delta]vrpar"->Function[{r},RealSign[x]*\[Delta]vrfunRPlungePar[r,a,p,e,x]],
	 "\[Delta]v\[Phi]par"->Function[{r},RealSign[x]*\[Delta]v\[Phi]funRPlungePar[r,a,p,e,x]]
	|>
]


(* ::Subsubsection::Closed:: *)
(*Plunge related to bound orbits*)


KerrNearEqSpinOrbitCorrRPlungeBoundOrbit[a_, p_, e_, x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz},
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	(* spin correction constants of motion *)
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dEEdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfunRPer[a,p,e,x]+dLzdefun[a,p,e,x]\[Delta]efunRPer[a,p,e,x];
	
	(*The RealSign serves to match the sign convetion of generic orbits solutions in the limit z->0*)				
	<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->RealSign[x]*\[Delta]EE,
	 "Jzs"->RealSign[x]*\[Delta]Lz,
	 (*Shift separatrix*)
	 "\[Delta]p"->RealSign[x]*\[Delta]pfunRPer[a,p,e,x],
	 "\[Delta]e"->RealSign[x]*\[Delta]efunRPer[a,p,e,x],
	 (*Shift frequencies*)
	 "BLFrequenciesGeo"->Function[{r},\[CapitalOmega]\[Phi]gfun[r,a,p,e,x]],
	 "BLFrequenciesCorrection"->Function[{r},RealSign[x]*\[Delta]\[CapitalOmega]\[Phi]funRPlungeBoundOrbit[r,a,p,e,x]],
	 "BLPrecessionFrequency"->Function[{r},\[CapitalOmega]pfun[r,a,p,e,x]],
	 "\[Delta]\[CapitalUpsilon]rhat"->RealSign[x]*\[Delta]\[CapitalUpsilon]rfunR[a,p,e,x]/\[CapitalUpsilon]rgfun[a,p,e,x],
	 (*Keys geodesic trajectories*)
	 "trg"->Function[{r},tgfunPlungeBoundOrbit[r,a,p,e,x]],
	 "\[Phi]rg"->Function[{r},\[Phi]gfunPlungeBoundOrbit[r,a,p,e,x]],
	 "rg"->Function[{\[Lambda]},rgfunPlungeBoundOrbit[\[Lambda],a,p,e,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{r},Vtrgfun[r,a,p,e,x]],
	 "vrg"->Function[{r},drgd\[Lambda]funPlungeBoundOrbit[r,a,p,e,x]],
	 "v\[Phi]g"->Function[{r},V\[Phi]rgfun[r,a,p,e,x]],
	 (*Keys corrections trajectory*)
	 "\[Delta]tpar"->Function[{r},RealSign[x]*\[Delta]tfunRPlungeBoundOrbitPar[r,a,p,e,x]],
	 "\[Delta]rpar"->Function[{r},RealSign[x]*\[Delta]rfunRPlungeBoundOrbitPar[r,a,p,e,x]],
	 "\[Psi]p"->Function[{r},\[Psi]pPlungeBoundOrbit[r,a,p,e,x]],
	 "\[Delta]zort"->Function[{r},\[Delta]zfunPlungeBoundOrbit[r,a,p,e,x]],
	 "\[Delta]\[Phi]par"->Function[{r},RealSign[x]*\[Delta]\[Phi]funRPlungeBoundOrbitPar[r,a,p,e,x]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{r},RealSign[x]*\[Delta]vtfunRPlungeBoundOrbitPar[r,a,p,e,x]],
	 "\[Delta]vrpar"->Function[{r},RealSign[x]*\[Delta]vrfunRPlungeBoundOrbitPar[r,a,p,e,x]],
	 "\[Delta]v\[Phi]par"->Function[{r},RealSign[x]*\[Delta]v\[Phi]funRPlungeBoundOrbitPar[r,a,p,e,x]]
	|>
]


(* ::Section::Closed:: *)
(*End package*)


End[];


EndPackage[];
