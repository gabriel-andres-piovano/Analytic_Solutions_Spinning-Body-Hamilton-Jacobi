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


KerrNearEqSpinOrbitCorrFTPer::usage = "KerrNearEqSpinOrbitCorrFTPer[a, p, e, x] calculates the linear corrections to periodic orbits in the fixed turning points spin-gauge. The initial radius at \[Lambda]=0 is the geodesic periastron.";


KerrNearEqSpinOrbitCorrFTPerFourier::usage = "KerrNearEqSpinOrbitCorrFTPerFourier[a, p, e, x, nmax] calculates the linear corrections to periodic orbits in the fixed turning points spin-gauge. The coordinate time and azimuthal trajectories are expanded in Fourier series. The initial radius at \[Lambda]=0 is the geodesic periastron.";


KerrNearEqSpinOrbitCorrFCPer::usage = "KerrNearEqSpinOrbitCorrFCPer[a, p, e, x] calculates the linear corrections to periodic orbits in the fixed constants of motion spin-gauge. The initial radius at \[Lambda]=0 is the geodesic periastron.";


KerrNearEqSpinOrbitCorrFCPerFourier::usage = "KerrNearEqSpinOrbitCorrFCPerFourier[a, p, e, x, nmax]  calculates the linear corrections to periodic orbits in the fixed constants of motion spin-gauge. The coordinate time and azimuthal trajectories are expanded in Fourier series. The initial radius at \[Lambda]=0 is the geodesic periastron.";


KerrNearEqSpinOrbitCorrFEPer::usage = "KerrNearEqSpinOrbitCorrFEPer[a, p, e, x]  calculates the linear corrections to periodic orbits in the fixed eccentricity spin-gauge. The initial radius at \[Lambda]=0 is the geodesic apoastron.";


KerrNearEqSpinOrbitCorrFEPerFourier::usage = "KerrNearEqSpinOrbitCorrFEPerFourier[a, p, e, x, nmax]  calculates the linear corrections to periodic orbits in the fixed eccentricity spin-gauge. The coordinate time and azimuthal trajectories are expanded in Fourier series. The initial radius at \[Lambda]=0 is the geodesic periastron.";


KerrNearEqSpinOrbitCorrFEHom::usage = "KerrNearEqSpinOrbitCorrFEHom[a, p, e, x] calculates the linear corrections to homoclinic orbits.";


KerrNearEqSpinOrbitCorrFEISCOPlunge::usage = "KerrNearEqSpinOrbitCorrFEISCOPlunge[a, p, x] calculates the linear corrections to homoclinic orbits.";


KerrNearEqSpinOrbitCorrFECritPlunge::usage = "KerrNearEqSpinOrbitCorrFECritPlunge[a, p, e, x] calculates the linear corrections to homoclinic orbits.";


KerrNearEqSpinOrbitCorrFEPlunge::usage = "KerrNearEqSpinOrbitCorrFEPlunge[a, p, e, x] calculates the linear corrections to homoclinic orbits.";


KerrNearEqSpinOrbitCorrFEPlungeBoundOrbit::usage = "KerrNearEqSpinOrbitCorrFEBoundOrbit[a, p, e, x] calculates the linear corrections to homoclinic orbits.";


KerrNearEqSpinOrbitCorrFEPlungeBoundOrbitDoubleRoot::usage = "KerrNearEqSpinOrbitCorrFEDoubleRoot[a, p, x] calculates the linear corrections to homoclinic orbits.";


KerrNearEqSpinOrbitCorrFE\[Delta]IBCO::usage = "KerrNearEqSpinOrbitCorrFEHom[a, p, x] provides the linear shifts to constants of motion and IBCO radius for parabolic orbits.";


Begin["`Private`"];


(*ADD ME: add algorithms to automatic select radial modes based on the sought precision for Fourier series expansions.*)


(*IMPORTANT: still to check for memory leakages on Fourier series expansion*)


(* ::Subsection::Closed:: *)
(*Building Fourier series*)


(* ::Text:: *)
(*Direct Fourier expansion of the trajectories. *)


FourierTra[wr_,coeff_]:=Module[{dim},
	dim=(Length[coeff]-1)/2;
	Re[Sum[2Sin[n*wr]Im[coeff][[n+dim+1]],{n,1,dim}]]
];


(* ::Text:: *)
(*Direct Fourier expansion of velocities*)


FourierVel[wr_,coeff_]:=Module[{dim},
	dim=(Length[coeff]-1)/2;
	Re[coeff][[dim+1]]+Re[Sum[2Cos[n*wr]Re[coeff][[n+dim+1]],{n,1,dim}]]
];


(* ::Subsubsection::Closed:: *)
(*Building Fourier series for spin corrections coordinate time and azimuthal orbits*)


(* ::Text:: *)
(*Integration Fourier expansion  of the velocities*)


\[CapitalDelta]\[Delta]IntvelPerPar[wr_,\[CapitalUpsilon]rg_,\[Delta]\[CapitalUpsilon]r_,\[Delta]coeff_,coeffVgr_]:=Module[{dimr},
	dimr=(Length[\[Delta]coeff]-1)/2;
	
	Re[Sum[2Sin[n*wr]\[Delta]coeff[[n+dimr+1]]/(n*\[CapitalUpsilon]rg),{n,1,dimr}]-\[Delta]\[CapitalUpsilon]r*Sum[2Sin[n*wr]coeffVgr[[n+dimr+1]]/(n*\[CapitalUpsilon]rg^2),{n,1,dimr}]]
]


(* ::Section:: *)
(*Geodesic quantities*)


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


dEEdpfun[a_,p_,e_,xg_]:=Module[{EEg,auxsqr,auxden},
	EEg=EEgfun[a,p,e,xg];
	auxsqr=Sqrt[p^11(a^4(1-e^2)^2+(-4e^2+(-2+p)^2)p^2+2a^2*p(-2+p+e^2(2+p)))];
	auxden=(-4a^2(1-e^2)^2+(3+e^2-p)^2*p);
	
	-(((-14a^2(1-e^2)^2+(12+4e^2-5p)(3+e^2-p)p)((2+2e-p)(3+e^2-p)p(-2+2e+p)-a^2(1-e^2)^2(-5+e^2+3p)))/(EEg*p^2auxden^2))+(2(1-e^2)^2(-14a^2(1-e^2)^2+(12+4e^2-5p)(3+e^2-p)p)a*auxsqr*RealSign[xg])/(EEg*p^8*auxden^2)-(a (1-e^2)^2*p^3(11a^4(1-e^2)^2+48a^2(-1+e^2)p+26(2+a^2+(-2+a^2)e^2)p^2-56p^3+15p^4)RealSign[xg])/(2EEg*auxden*auxsqr)+(-3a^2(1-e^2)^2(-10+2e^2+7p)+p(28(-3+2e^2+e^4)+128p-9(7+e^2)p^2+10p^3))/(2EEg*p^2*auxden)
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

	(-EEg*gginc[r2g]+Lzgzaux*RealSign[xg])/((1+e)hhinc[r2g]^2) (2(1+e-p))/(1+e)-(gginc[r2g] dEEdp)/hhinc[r2g]-EEg/hhinc[r2g] (2a)/(1+e)+RealSign[xg]/hhinc[r2g] ((hhinc[r2g]*EEg*2ffinc[r2g]*dEEdp)/(2Lzgzaux)-hhinc[r2g] /(2Lzgzaux) (2p(a^2(1+e)^2+p(-3-3e+2p)))/(1+e)^4+(hhinc[r2g]EEg^2)/(2Lzgzaux) ((4p^3+2a^2(1+e)^2(1+e+p))/(1+e)^4)+(2 EEg^2*gginc[r2g])/(2Lzgzaux) (2a)/(1+e)+(2gginc[r2g]^2EEg*dEEdp)/(2Lzgzaux)-(-ddinc[r2g]+EEg^2*ffinc[r2g]) /(2Lzgzaux) (2(1+e-p))/(1+e)^2)
]


(* ::Subsubsection::Closed:: *)
(*Equatorial orbits  - complex conjugate roots*)


EEgfunCC[p_,e_]:=Module[{r1g,\[Rho]rg},
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);

	Sqrt[(-2+r1g+2\[Rho]rg)/(r1g+2\[Rho]rg)]
]


LzgfunCC[a_,p_,e_,x_]:=Module[{r1g,\[Rho]rg,\[CapitalDelta]r1g},
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[CapitalDelta]r1g=a^2+(-2+r1g)r1g;

	2 RealSign[x] Sqrt[((-2r1g+r1g^2)r1g*\[Rho]rg-2RealSign[x]Sqrt[a^2*r1g*\[CapitalDelta]r1g*\[Rho]rg(-2+r1g+2 \[Rho]rg)]+a^2(-2+r1g+2\[Rho]rg+r1g*\[Rho]rg))/((-2+r1g)^2(r1g+2\[Rho]rg))]
]


\[Rho]igfun[a_,p_,e_,x_]:=Module[{r1g,\[Rho]rg,\[CapitalDelta]r1g},
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[CapitalDelta]r1g=a^2+(-2+r1g)r1g;

	-Sqrt[(1/(-2+r1g)^2 (-((-2+r1g)(r1g(-4+\[Rho]rg)-2\[Rho]rg)\[Rho]rg)-4RealSign[x]Sqrt[a^2*r1g*\[CapitalDelta]r1g*\[Rho]rg(-2+r1g+2 \[Rho]rg)]+a^2((-2+r1g)r1g+2(2+r1g)\[Rho]rg)))]
]


dEEdpfunCC[p_,e_]:=(1-e^2)/((3-e)p^2*EEgfunCC[p,e])


dLzdpfunCC[a_,p_,e_,x_]:=Module[{Lzg,r1g,\[Rho]rg,\[CapitalDelta]r1g},
	Lzg=LzgfunCC[a,p,e,x];
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[CapitalDelta]r1g=a^2+(-2+r1g)r1g;

	2/Lzg ((-p^3(-2+2e+p)(-4+4e+p)+a^2 (-1+e)^2 (4(-1+e)^2 (1+e)+6(-1+e^2)p-4(-2+e)p^2+p^3)) /((-3+e)p^2(-2+2e+p)^3)+(2(-1+e)^3 (1+e)(-2+2e+3p)RealSign[x]Sqrt[a^2*r1g*\[CapitalDelta]r1g*\[Rho]rg(-2+r1g+2\[Rho]rg)])/((-3+e)p^2(-2+2e+p)^3)- RealSign[x]/((-3+e)p (-2+2e+p)^2 Sqrt[a^2*r1g*\[CapitalDelta]r1g*\[Rho]rg(-2+r1g+2\[Rho]rg)]) (a^2*p(a^2 (-1+e)^2 (-4+4e^2+9p-3e*p)+p(12(-1+e)^2(1+e)+32(-1+e)p-5(-3+e)p^2)))/((-1+e)(1+e)))
]


(* ::Subsection::Closed:: *)
(*Geodesic frequencies*)


(* ::Subsubsection::Closed:: *)
(*Radial and polar geodesic frequencies*)


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
(*Coordinate time geodesic frequency*)


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
(*Azimuthal geodesic frequency*)


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
	hp=hr ( r3g-rp)/(r2g-rp);
	hm=hr ( r3g-rm)/(r2g-rm);

	a/(rp-rm) ((2EEg*rp-a*Lzg)/(r3g-rp) (1-(r2g-r3g)/(r2g-rp) EllipticPi[hp,krg]/ellK))-a/(rp-rm) ((2EEg*rm-a*Lzg)/(r3g-rm) (1-(r2g-r3g)/(r2g-rm) EllipticPi[hm,krg]/ellK))
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
(*Spin precession frequency*)


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
(*Azimuthal BL time geodesic frequency - plunge orbits*)


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
	EEg=EEgfunCC[p,e];
	Lzg=LzgfunCC[a,p,e,x];

	(2a*EEg+Lzg(-2+r))/(-2a*Lzg+EEg*r^3+a^2*EEg(2+r))
]


(* ::Subsubsection::Closed:: *)
(*Spin precession frequency - plunge orbits*)


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
	EEg=EEgfunCC[p,e];
	Lzg=LzgfunCC[a,p,e,x];
	Lzred=Lzg-a*EEg;
	
	Lzred(a+EEg*Lzred)r^2/(Abs[Lzred](r^2+Lzred^2))*(a^2+(-2+r)r)/(r(-2a*Lzg+EEg*r^3+a^2*EEg(2+r)))
]


(* ::Subsection:: *)
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


(* ::Subsubsection:: *)
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


(* ::Subsubsection:: *)
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


(* ::Subsubsection::Closed:: *)
(*Fourier coefficients of the geodesic coordinate time and azimuthal velocities*)


VtrgcoeffICr2g[nmax_,a_,p_,e_,x_]:=Module[{stepsr,wrlist,Vtrglist,ExpniTable},
	stepsr=4*nmax;
	wrlist=Table[i,{i,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	Vtrglist=VtrgICr2gfun[wrlist,a,p,e,x];
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e}]],{n,-nmax,nmax},{i,1,stepsr}];
	Chop[(ExpniTable . Vtrglist)/stepsr,10^(-16)]
]


V\[Phi]rgcoeffICr2g[nmax_,a_,p_,e_,x_]:=Module[{stepsr,wrlist,V\[Phi]rglist,ExpniTable},
	stepsr=4*nmax;
	wrlist=Table[i,{i,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	V\[Phi]rglist=V\[Phi]rgICr2gfun[wrlist,a,p,e,x];
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e}]],{n,-nmax,nmax},{i,1,stepsr}];
	Chop[(ExpniTable . V\[Phi]rglist)/stepsr,10^(-16)]
]


VtrgcoeffICr1g[nmax_,a_,p_,e_,x_]:=Module[{stepsr,wrlist,Vtrglist,ExpniTable},
	stepsr=4*nmax;
	wrlist=Table[i,{i,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	Vtrglist=VtrgICr1gfun[wrlist,a,p,e,x];
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e}]],{n,-nmax,nmax},{i,1,stepsr}];
	Chop[(ExpniTable . Vtrglist)/stepsr,10^(-16)]
]


V\[Phi]rgcoeffICr1g[nmax_,a_,p_,e_,x_]:=Module[{stepsr,wrlist,V\[Phi]rglist,ExpniTable},
	stepsr=4*nmax;
	wrlist=Table[i,{i,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	V\[Phi]rglist=V\[Phi]rgICr1gfun[wrlist,a,p,e,x];
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e}]],{n,-nmax,nmax},{i,1,stepsr}];
	Chop[(ExpniTable . V\[Phi]rglist)/stepsr,10^(-16)]
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


(* ::Subsubsection::Closed:: *)
(*Coordinate-time geodesic velocity*)


VtrgfunHom[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];

	EEg((r^2+a^2)^2/(r^2-2r+a^2))-(2a*r)/(r^2-2r+a^2)Lzg-EEg*a^2
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal geodesic velocity*)


V\[Phi]rgfunHom[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];

	a/(r^2-2r+a^2)(EEg(r^2+a^2)-a*Lzg)-a*EEg+Lzg
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


(* ::Subsubsection::Closed:: *)
(*Geodesic radial velocity*)


drgd\[Lambda]funCritplunge[r_,a_,p_,e_,x_]:=Module[{EEg,r1g,r2g},
	EEg=EEgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);

	Sqrt[1-EEg^2](r-r2g)Sqrt[r(r1g-r)]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate-time geodesic velocity*)


VtrgfunCritplunge[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];

	EEg((r^2+a^2)^2/(r^2-2r+a^2))-(2a*r)/(r^2-2r+a^2)Lzg-EEg*a^2
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal geodesic velocity*)


V\[Phi]rgfunCritplunge[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];

	a/(r^2-2r+a^2)(EEg(r^2+a^2)-a*Lzg)-a*EEg+Lzg
]


(* ::Subsection::Closed:: *)
(*Geodesics Generic plunge*)


(* ::Subsubsection::Closed:: *)
(*Geodesic radial trajectory*)


rgfunPlunge[\[Lambda]_,a_,p_,e_,x_]:=Module[{EEg,r1g,\[Rho]rg,\[Rho]ig,Ar1g,B0,krg},
	EEg=EEgfunCC[p,e];
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
	EEg=EEgfunCC[p,e];
	Lzg=LzgfunCC[a,p,e,x];
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
	EEg=EEgfunCC[p,e];
	Lzg=LzgfunCC[a,p,e,x];
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
	EEg=EEgfunCC[p,e];
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];

	Sqrt[1-EEg^2]Sqrt[(r-\[Rho]rg-I*\[Rho]ig)(r-\[Rho]rg+I*\[Rho]ig)]Sqrt[(r1g-r)r]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate-time geodesic velocity*)


VtrgfunPlunge[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg},
	EEg=EEgfunCC[p,e];
	Lzg=LzgfunCC[a,p,e,x];

	EEg((r^2+a^2)^2/(r^2-2r+a^2))-(2a*r)/(r^2-2r+a^2)Lzg-EEg*a^2
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal geodesic velocity*)


V\[Phi]rgfunPlunge[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg},
	EEg=EEgfunCC[p,e];
	Lzg=LzgfunCC[a,p,e,x];

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


(* ::Subsubsection::Closed:: *)
(*Coordinate-time geodesic velocity*)


VtrgfunPlungeBoundOrbit[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];

	EEg((r^2+a^2)^2/(r^2-2r+a^2))-(2a*r)/(r^2-2r+a^2)Lzg-EEg*a^2
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal geodesic velocity*)


V\[Phi]rgfunPlungeBoundOrbit[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];

	a/(r^2-2r+a^2)(EEg(r^2+a^2)-a*Lzg)-a*EEg+Lzg
]


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


\[Delta]EEfunCC[a_,p_,e_,x_]:=Module[{EEg,Lzg},
	EEg=EEgfunCC[p,e];
	Lzg=LzgfunCC[a,p,e,x];

	((-1+e^2)^2(a^3(-1+EEg^2)-a*Lzg^2))/(2(-3+e)^2*EEg (-a*EEg+Lzg)^2*p^2)
]


\[Delta]LzfunCC[a_,p_,e_,x_]:=Module[{r1g,EEg,Lzg,\[Delta]EE},
	r1g=p/(1-e);
	EEg=EEgfunCC[p,e];
	Lzg=LzgfunCC[a,p,e,x];

	\[Delta]EE=\[Delta]EEfunCC[a,p,e,x];

	(a (-a*EEg+Lzg)^2+3EEg(a*EEg-Lzg)r1g^2+EEg*Lzg*r1g^3)/((2a*EEg+Lzg(-2+r1g))r1g^2)+((-2a*Lzg+EEg*r1g^3+a^2*EEg(2+r1g))\[Delta]EE)/(2a*EEg+Lzg(-2+r1g))
]


(* ::Subsection::Closed:: *)
(*Shift ISCO- fixed eccentricity*)


\[Delta]EEISCOfunFE[a_,x_]:=Module[{r1g,EEg,Lzg,y},
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,x];
	Lzg=Lzgfun[a,r1g,0,x];
	y=Lzg-a*EEg;

	-((y(4a*y^2-3a*y*Lzg*r1g+EEg*Lzg*r1g^3))/(2a*y^2*r1g^3-6EEg*y*r1g^5+2EEg*Lzg*r1g^6))
]


\[Delta]LzISCOfunFE[a_,x_]:=Module[{r1g,EEg,Lzg,y},
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,x];
	Lzg=Lzgfun[a,r1g,0,x];
	y=Lzg-a*EEg;

	(-4a^2y^3+3a^3*EEg*y^2*r1g-a*EEg(8a*EEg-7Lzg)y*r1g^3-9EEg^2*y*r1g^5+2EEg^2*Lzg*r1g^6)/(2a*y^2*r1g^3-6EEg*y*r1g^5+2EEg*Lzg*r1g^6)
]


(* ::Section::Closed:: *)
(*Shifts to the radial potential roots*)


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


\[Delta]r3funFT[a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,Lzred},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	Lzred=Lzg-a*EEg;

	a (a^2(1-EEg^2)+Lzg^2)/(2Lzred^2)+(4EEg*\[Delta]EE)/(1-EEg^2)^2
]


\[Delta]r41fun[a_,p_,e_,x_]:=Module[{EEg,Lzg,Lzred},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	Lzred=Lzg-a*EEg;

	-a(a^2(1-EEg^2)+Lzg^2)/(4Lzred^2)
]


\[Delta]r41funCC[a_,p_,e_,x_]:=Module[{EEg,Lzg,Lzred},
	EEg=EEgfunCC[p,e];
	Lzg=LzgfunCC[a,p,e,x];
	Lzred=Lzg-a*EEg;

	-a(a^2(1-EEg^2)+Lzg^2)/(4Lzred^2)
]


\[Delta]r1funFE[a_,p_,e_,x_]:=\[Delta]pfun[a,p,e,x]/(1-e)


\[Delta]r2funFE[a_,p_,e_,x_]:=\[Delta]pfun[a,p,e,x]/(1+e)


\[Delta]rISCOfunFE[a_,x_]:=Module[{r1g,EEg,Lzg,y},
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,x];
	Lzg=Lzgfun[a,r1g,0,x];
	y=Lzg-a*EEg;

	1/6 ((a^3(1-EEg^2)+a*Lzg^2)/y^2-(4EEg*y(4a*y^2-3a*y*Lzg*r1g+EEg*Lzg*r1g^3))/((1-EEg^2)^2*r1g^3(a*y^2-3EEg*y*r1g^2+EEg*Lzg*r1g^3)))
]


(* ::Subsection::Closed:: *)
(*Shift semilatus rectum*)


\[Delta]pfun[a_,p_,e_,x_]:=Module[{EEg,dEdp,\[Delta]EE,\[Delta]\[Rho]r4},
	EEg=EEgfun[a,p,e,x];
	dEdp=dEEdpfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];

	(2(1-e^2)(2*EEg*\[Delta]EE-(1-EEg^2)^2\[Delta]\[Rho]r4))/(-4dEdp*EEg+4*dEdp*e^2*EEg+3(1-EEg^2)^2-e (1-EEg^2)^2)
]


(* ::Text:: *)
(*For generic plunge*)


\[Delta]pfunCC[a_,p_,e_,x_]:=Module[{r1g,\[Rho]rg,EEg,Lzg,Lzred,\[Delta]EE,dEdp,dLzdp,\[Delta]\[Rho]r4,num,den},
	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);

	EEg=EEgfunCC[p,e];
	Lzg=LzgfunCC[a,p,e,x];
	Lzred=Lzg-a*EEg;

	\[Delta]EE=\[Delta]EEfunCC[a,p,e,x];
	dEdp=dEEdpfunCC[p,e];
	dLzdp=dLzdpfunCC[a,p,e,x];
	\[Delta]\[Rho]r4=\[Delta]r41funCC[a,p,e,x];

	num=(1-e)^2/(-3+e)^2 (8\[Delta]\[Rho]r4)/\[Rho]rg^2 (-2+((1-EEg^2)(2Lzred^2+r1g^2))/r1g+(2(1-EEg^2)^2*Lzred(EEg*r1g^2-a*Lzred))/( EEg(Lzg*r1g-2Lzred)))-(8(1-EEg^2)^2*Lzred^2\[Delta]\[Rho]r4)/r1g^2-(4a (1-EEg^2)^2)/r1g+(8a (1-EEg^2)^2Lzred^3)/((Lzg*r1g-2Lzred)r1g^3)-(4(1-EEg^2)^2*Lzred(EEg Lzg))/(Lzg*r1g-2Lzred)-(1-EEg^2)^2*2\[Delta]\[Rho]r4(2-(1-EEg^2)r1g);

	den=-((8dLzdp (1-EEg^2)^2*Lzred)/r1g)+4dEdp(2EEg+(2(1-EEg^2)(a(1-EEg^2)Lzred-EEg*Lzred^2))/r1g-EEg(1-EEg^2)r1g)+((1-EEg^2)^2(-4Lzred^2+r1g^2(2+(-1+EEg^2)r1g)))/((-1+e)r1g^2);

	num/den
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

	-2 \[CapitalUpsilon]rg^2/(2\[Pi]) 1/Sqrt[1-EEg^2] 2/Sqrt[(r1g-r3g)r2g] ((r2g*ellE-r1g*ellK)/(2r1g(r1g-r2g) ) \[Delta]r1+1/2 (-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellK/(r1g-r2g))\[Delta]r2+1/2 \[Delta]r3(1/r3g (r2g/(r2g-r3g) ellE-ellK))+\[Delta]r41((r1g*ellK-(r1g-r3g)ellE)/(r1g*r3g ))-1/2 a(-(1/(3r1g*r2g*r3g))((r1g+r2g+r3g)ellK-2(r2g*r3g+r1g(r2g+r3g))( (-((r1g-r3g)ellE)+r1g*ellK)/(r1g*r3g)))))
]


(* ::Subsubsection::Closed:: *)
(*Fixed eccentricity*)


\[Delta]\[CapitalUpsilon]rfunFE[a_,p_,e_,x_]:=Module[{\[CapitalUpsilon]rg,EEg,Lzg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,krg,ellK,ellE},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
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
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	If[e==0,
		1/2 (a/r1g^2-(2EEg*\[Delta]EE)/(1-EEg^2)+\[Delta]r1/r1g-(2\[Delta]r41)/r1g),
		1/2 (a/r2g^2-(2EEg*\[Delta]EE)/(1-EEg^2)+\[Delta]r1/(r1g-r2g)-\[Delta]r2/(r1g-r2g)+\[Delta]r2/r2g-(2\[Delta]r41)/r2g)
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
(*Fixed eccentricity*)


\[Delta]\[CapitalUpsilon]tfunFE[a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[CapitalUpsilon]rg,\[CapitalUpsilon]tg,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,krg,\[Gamma]r,ellK,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
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
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

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
	\[ScriptCapitalI]r1g=2 /(Sqrt[1-EEg^2] Sqrt[r2g (r1g-r3g)] r1g) ((r2g*ellE-r1g*ellK)/(r1g-r2g) );
	\[ScriptCapitalI]r2g=2 /(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]) (-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellK/(r1g-r2g));
	\[ScriptCapitalI]r3g=2 /(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]r3g) ((r2g*ellE)/(r2g-r3g)-ellK);

	(2\[CapitalUpsilon]rg)/(2\[Pi]) (-EEg*\[ScriptCapitalI]+((2a*EEg+Lzg(-2+r1g))r1g*\[ScriptCapitalI]r1g*\[Delta]r1)/(2(r1g-rm)(r1g-rp))+((2a*EEg+Lzg(-2+r2g))r2g*\[ScriptCapitalI]r2g*\[Delta]r2)/(2(r2g-rm)(r2g-rp))+((2a*EEg+Lzg(-2+r3g))r3g*\[ScriptCapitalI]r3g*\[Delta]r3)/(2(r3g-rm)(r3g-rp)))+\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]rg \[Delta]\[CapitalUpsilon]rfunFC[a,p,e,x]
]


(* ::Subsubsection::Closed:: *)
(*Fixed eccentricity*)


\[Delta]\[CapitalUpsilon]\[Phi]funFE[a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,krg,ellK,ellE,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
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
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

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


(* ::Subsection::Closed:: *)
(*Shift azimuthal BL frequency*)


\[Delta]\[CapitalOmega]\[Phi]funISCOplunge[r_,a_,x_]:=Module[{r1g,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]r,\[CapitalDelta]r,Lzred},
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,x];
	Lzg=Lzgfun[a,r1g,0,x];
	\[Delta]EE=\[Delta]EEISCOfunFE[a,x];
	\[Delta]Lz=\[Delta]LzISCOfunFE[a,x];
	\[Delta]r=\[Delta]rfunISCOplungePar[r,a,x];
	\[CapitalDelta]r=a^2+(-2+r)r ;
	Lzred=Lzg-a*EEg;

	(\[CapitalDelta]r(Lzred^2-Lzg*r^3\[Delta]EE+EEg*r^3(-EEg+\[Delta]Lz))-2r(a*Lzred^2-3EEg*Lzred*r^2+EEg*Lzg*r^3)\[Delta]r)/(r (r*EEg(a^2+r^2)-2a*Lzred)^2)
]


\[Delta]\[CapitalOmega]\[Phi]funCritplunge[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]r,\[CapitalDelta]r,Lzred},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]r=\[Delta]rfunCritplungePar[r,a,p,e,x];
	\[CapitalDelta]r=a^2+(-2+r)r ;
	Lzred=Lzg-a*EEg;

	(\[CapitalDelta]r(Lzred^2-Lzg*r^3\[Delta]EE+EEg*r^3(-EEg+\[Delta]Lz))-2r(a*Lzred^2-3EEg*Lzred*r^2+EEg*Lzg*r^3)\[Delta]r)/(r (r*EEg(a^2+r^2)-2a*Lzred)^2)
]


\[Delta]\[CapitalOmega]\[Phi]funPlunge[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]EE,\[Delta]Lz,\[Delta]r,\[CapitalDelta]r,Lzred},
	EEg=EEgfunCC[p,e];
	Lzg=LzgfunCC[a,p,e,x];
	\[Delta]p=\[Delta]pfunCC[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunCC[a,p,e,x]+dEEdpfunCC[p,e]\[Delta]p;
	\[Delta]Lz=\[Delta]LzfunCC[a,p,e,x]+dLzdpfunCC[a,p,e,x]\[Delta]p;
	\[Delta]r=\[Delta]rfunPlungePar[r,a,p,e,x];
	\[CapitalDelta]r=a^2+(-2+r)r ;
	Lzred=Lzg-a*EEg;

	(\[CapitalDelta]r(Lzred^2-Lzg*r^3\[Delta]EE+EEg*r^3(-EEg+\[Delta]Lz))-2r(a*Lzred^2-3EEg*Lzred*r^2+EEg*Lzg*r^3)\[Delta]r)/(r (r*EEg(a^2+r^2)-2a*Lzred)^2)
]


\[Delta]\[CapitalOmega]\[Phi]funPlungeBoundOrbit[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]EE,\[Delta]Lz,\[Delta]r,\[CapitalDelta]r,Lzred},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]p=\[Delta]pfunCC[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]p;
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]p;
	\[Delta]r=\[Delta]rfunPlungeBoundOrbitPar[r,a,p,e,x];
	\[CapitalDelta]r=a^2+(-2+r)r ;
	Lzred=Lzg-a*EEg;

	(\[CapitalDelta]r(Lzred^2-Lzg*r^3\[Delta]EE+EEg*r^3(-EEg+\[Delta]Lz))-2r(a*Lzred^2-3EEg*Lzred*r^2+EEg*Lzg*r^3)\[Delta]r)/(r (r*EEg(a^2+r^2)-2a*Lzred)^2)
]


\[Delta]\[CapitalOmega]\[Phi]funPlungeBoundOrbitDoubleRoot[r_,a_,p_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]EE,\[Delta]Lz,\[Delta]r,\[CapitalDelta]r,Lzred},
	EEg=EEgfun[a,p,0,x];
	Lzg=Lzgfun[a,p,0,x];
	\[Delta]p=\[Delta]pfunCC[a,p,0,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,0,x]+dEEdpfun[a,p,0,x]\[Delta]p;
	\[Delta]Lz=\[Delta]LzfunFT[a,p,0,x]+dLzdpfun[a,p,0,x]\[Delta]p;
	\[Delta]r=\[Delta]rfunPlungeBoundOrbitDoubleRootPar[r,a,p,x];
	\[CapitalDelta]r=a^2+(-2+r)r ;
	Lzred=Lzg-a*EEg;

	(\[CapitalDelta]r(Lzred^2-Lzg*r^3\[Delta]EE+EEg*r^3(-EEg+\[Delta]Lz))-2r(a*Lzred^2-3EEg*Lzred*r^2+EEg*Lzg*r^3)\[Delta]r)/(r (r*EEg(a^2+r^2)-2a*Lzred)^2)
]


(* ::Section:: *)
(*Spin corrections to the orbit - parallel component of the spin, fixed turning points*)


(* ::Subsection:: *)
(*Bound trajectory*)


(* ::Subsubsection:: *)
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


(* ::Subsubsection:: *)
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


(* ::Subsubsection:: *)
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
(*Fourier coefficients of the spin-correction to the coordinate time trajectory*)


\[Delta]tcoefffunFTPerPar[nmax_,a_,p_,e_,x_]:=Module[{stepsr,wrlist,\[Delta]tlist,ExpniTable},
	stepsr=4*nmax;
	wrlist=Table[i,{i,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	\[Delta]tlist=\[Delta]tfunFTPerPar[wrlist,a,p,e,x];

	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e}]],{n,-nmax,nmax},{i,1,stepsr}];
	Chop[(ExpniTable . \[Delta]tlist)/stepsr,10^(-16)]
]


(* ::Subsubsection::Closed:: *)
(*Fourier coefficients of the spin-correction to the coordinate time velocity*)


d\[Delta]td\[Lambda]coefffunFTPerPar[nmax_,a_,p_,e_,x_]:=Module[{stepsr,wrlist,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]vttlist,ExpniTable},
	stepsr=4*nmax;
	wrlist=Table[i,{i,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	\[Delta]vttlist=\[Delta]vtfunFTPerPar[wrlist,a,p,e,x];

	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e}]],{n,-nmax,nmax},{i,1,stepsr}];
	Chop[(ExpniTable . \[Delta]vttlist)/stepsr,10^(-16)]
]


(* ::Subsubsection:: *)
(*Azimuthal correction to the trajectory - purely oscillatory part*)


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


(* ::Subsubsection::Closed:: *)
(*Fourier coefficients of the spin-correction to the azimuthal trajectory*)


\[Delta]\[Phi]coefffunFTPerPar[nmax_,a_,p_,e_,x_]:=Module[{stepsr,wrlist,\[Delta]\[Phi]list,ExpniTable},
	stepsr=4*nmax;
	wrlist=Table[i,{i,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	\[Delta]\[Phi]list=\[Delta]\[Phi]funFTPerPar[wrlist,a,p,e,x];

	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e}]],{n,-nmax,nmax},{i,1,stepsr}];
	Chop[(ExpniTable . \[Delta]\[Phi]list)/stepsr,10^(-16)]
]


(* ::Subsubsection::Closed:: *)
(*Fourier coefficients of the spin-correction to the azimuthal velocity*)


d\[Delta]\[Phi]d\[Lambda]coefffunFTPerPar[nmax_,a_,p_,e_,x_]:=Module[{stepsr,wrlist,\[Delta]v\[Phi]tlist,ExpniTable},
	stepsr=4*nmax;
	wrlist=Table[i,{i,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	\[Delta]v\[Phi]tlist=\[Delta]v\[Phi]funFTPerPar[wrlist,a,p,e,x];

	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e}]],{n,-nmax,nmax},{i,1,stepsr}];
	Chop[(ExpniTable . \[Delta]v\[Phi]tlist)/stepsr,10^(-16)]
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


\[Delta]vtfunFTPerPar[wr_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,r3g,rg,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	rg=rgICr2gfun[wr,a,p,e,x];
	\[CapitalDelta]=a^2-2rg+rg^2;

	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+rg^2)+a*EEg(a^2+3rg^2))/(rg*\[CapitalDelta])+(rg(2a^2+a^2*rg+rg^3)\[Delta]EE)/\[CapitalDelta]-(2a*rg*\[Delta]Lz)/\[CapitalDelta];
	ddtgd\[Lambda]drfun=(EEg*rg(a^2+3rg^2))/\[CapitalDelta]-((rg^2-a^2)(-2a*Lzg+EEg*rg^3+a^2*EEg(2+rg)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunFTPerPar[wr,a,p,e,x]+d\[Delta]td\[Lambda]fun
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


\[Delta]v\[Phi]funFTPerPar[wr_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,r3g,rg,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	rg=rgICr2gfun[wr,a,p,e,x];
	\[CapitalDelta]=a^2-2rg+rg^2;

	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+rg)))/(rg*\[CapitalDelta])+(2a*rg*\[Delta]EE)/\[CapitalDelta]+((rg-2)rg*\[Delta]Lz)/\[CapitalDelta];
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(rg-1)-EEg*rg^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunFTPerPar[wr,a,p,e,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Section:: *)
(*Spin corrections to the orbit - parallel component of the spin, fixed constants of motion*)


(* ::Subsection:: *)
(*Bound trajectory*)


(* ::Subsubsection:: *)
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


(* ::Subsubsection:: *)
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


(* ::Subsubsection:: *)
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
(*Fourier coefficients of the spin-correction to the coordinate time trajectory*)


\[Delta]tcoefffunFCPerPar[nmax_,a_,p_,e_,x_]:=Module[{stepsr,wrlist,\[Delta]tlist,ExpniTable},
	stepsr=4*nmax;
	wrlist=Table[i,{i,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	\[Delta]tlist=\[Delta]tfunFCPerPar[wrlist,a,p,e,x];

	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e}]],{n,-nmax,nmax},{i,1,stepsr}];
	Chop[(ExpniTable . \[Delta]tlist)/stepsr,10^(-16)]
]


(* ::Subsubsection::Closed:: *)
(*Fourier coefficients of the spin-correction to the coordinate time velocity*)


d\[Delta]td\[Lambda]coefffunFCPerPar[nmax_,a_,p_,e_,x_]:=Module[{stepsr,wrlist,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]vttlist,ExpniTable},
	stepsr=4*nmax;
	wrlist=Table[i,{i,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	\[Delta]vttlist=\[Delta]vtfunFCPerPar[wrlist,a,p,e,x];

	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e}]],{n,-nmax,nmax},{i,1,stepsr}];
	Chop[(ExpniTable . \[Delta]vttlist)/stepsr,10^(-16)]
]


(* ::Subsubsection:: *)
(*Azimuthal correction to the trajectory - purely oscillatory part*)


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


(* ::Subsubsection::Closed:: *)
(*Fourier coefficients of the spin-correction to the azimuthal trajectory*)


\[Delta]\[Phi]coefffunFCPerPar[nmax_,a_,p_,e_,x_]:=Module[{stepsr,wrlist,\[Delta]\[Phi]list,ExpniTable},
	stepsr=4*nmax;
	wrlist=Table[i,{i,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	\[Delta]\[Phi]list=\[Delta]\[Phi]funFCPerPar[wrlist,a,p,e,x];

	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e}]],{n,-nmax,nmax},{i,1,stepsr}];
	Chop[(ExpniTable . \[Delta]\[Phi]list)/stepsr,10^(-16)]
]


(* ::Subsubsection::Closed:: *)
(*Fourier coefficients of the spin-correction to the azimuthal velocity*)


d\[Delta]\[Phi]d\[Lambda]coefffunFCPerPar[nmax_,a_,p_,e_,x_]:=Module[{stepsr,wrlist,\[Delta]v\[Phi]tlist,ExpniTable},
	stepsr=4*nmax;
	wrlist=Table[i,{i,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	\[Delta]v\[Phi]tlist=\[Delta]v\[Phi]funFCPerPar[wrlist,a,p,e,x];

	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e}]],{n,-nmax,nmax},{i,1,stepsr}];
	Chop[(ExpniTable . \[Delta]v\[Phi]tlist)/stepsr,10^(-16)]
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


\[Delta]vtfunFCPerPar[wr_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1g,r2g,r3g,rg,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	rg=rgICr2gfun[wr,a,p,e,x];
	
	\[CapitalDelta]=a^2-2rg+rg^2;

	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+rg^2)+a*EEg(a^2+3rg^2))/(rg*\[CapitalDelta]);
	ddtgd\[Lambda]drfun=(EEg*rg(a^2+3rg^2))/\[CapitalDelta]-((rg^2-a^2)(-2a*Lzg+EEg*rg^3+a^2*EEg(2+rg)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunFCPerPar[wr,a,p,e,x]+d\[Delta]td\[Lambda]fun
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


\[Delta]v\[Phi]funFCPerPar[wr_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1g,r2g,r3g,rg,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	rg=rgICr2gfun[wr,a,p,e,x];
	
	\[CapitalDelta]=a^2-2rg+rg^2;

	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+rg)))/(rg*\[CapitalDelta]);
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(rg-1)-EEg*rg^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunFCPerPar[wr,a,p,e,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Section:: *)
(*Spin corrections to the orbit - parallel component of the spin, fixed eccentricity*)


(* ::Subsection:: *)
(*Bound trajectory*)


(* ::Subsubsection:: *)
(*Radial trajectory*)


\[Delta]rfunFEPerPar[wr_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,\[Delta]r1,\[Delta]r2,EEg,\[Delta]EE,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE,rg},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
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
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
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


(* ::Subsubsection:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunFEPerParRes[wr_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,\[Delta]r1,\[Delta]r2,EEg,\[Delta]EE,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE,rg},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
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
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
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


(* ::Subsubsection:: *)
(*Coordinate time trajectory - purely oscillatory part*)


\[Delta]tfunFEPerPar[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[CapitalUpsilon]rg,\[CapitalUpsilon]tg,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,\[Phi],rg,dtgd\[Lambda]fun,dtgd\[Lambda]funr1g,dtgd\[Lambda]funr2g,dtgd\[Lambda]funr3g,krghold,krg,\[Gamma]r,ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
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
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
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
(*Fourier coefficients of the spin-correction to the coordinate time trajectory*)


\[Delta]tcoefffunFEPerPar[nmax_,a_,p_,e_,x_]:=Module[{stepsr,wrlist,\[Delta]tlist,ExpniTable},
	stepsr=4*nmax;
	wrlist=Table[i,{i,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	\[Delta]tlist=\[Delta]tfunFEPerPar[wrlist,a,p,e,x];

	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e}]],{n,-nmax,nmax},{i,1,stepsr}];
	Chop[(ExpniTable . \[Delta]tlist)/stepsr,10^(-16)]
]


(* ::Subsubsection::Closed:: *)
(*Fourier coefficients of the spin-correction to the coordinate time velocity*)


d\[Delta]td\[Lambda]coefffunFEPerPar[nmax_,a_,p_,e_,x_]:=Module[{stepsr,wrlist,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]vttlist,ExpniTable},
	stepsr=4*nmax;
	wrlist=Table[i,{i,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	\[Delta]vttlist=\[Delta]vtfunFEPerPar[wrlist,a,p,e,x];

	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e}]],{n,-nmax,nmax},{i,1,stepsr}];
	Chop[(ExpniTable . \[Delta]vttlist)/stepsr,10^(-16)]
]


(* ::Subsubsection:: *)
(*Azimuthal correction to the trajectory - purely oscillatory part*)


\[Delta]\[Phi]funFEPerPar[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Phi],rg,d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,d\[Phi]gd\[Lambda]funr2g,d\[Phi]gd\[Lambda]funr3g,krghold,krg,\[Gamma]r,ellF,ellE,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
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
(*Azimuthal correction to the trajectory as function of r*)


\[Delta]\[Phi]funOfrFEPerPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Phi],d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,d\[Phi]gd\[Lambda]funr2g,d\[Phi]gd\[Lambda]funr3g,krg,\[Gamma]r,ellF,ellE,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
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


(* ::Subsubsection::Closed:: *)
(*Fourier coefficients of the spin-correction to the azimuthal trajectory*)


\[Delta]\[Phi]coefffunFEPerPar[nmax_,a_,p_,e_,x_]:=Module[{stepsr,wrlist,\[Delta]\[Phi]list,ExpniTable},
	stepsr=4*nmax;
	wrlist=Table[i,{i,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	\[Delta]\[Phi]list=\[Delta]\[Phi]funFEPerPar[wrlist,a,p,e,x];

	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e}]],{n,-nmax,nmax},{i,1,stepsr}];
	Chop[(ExpniTable . \[Delta]\[Phi]list)/stepsr,10^(-16)]
]


(* ::Subsubsection::Closed:: *)
(*Fourier coefficients of the spin-correction to the azimuthal velocity*)


d\[Delta]\[Phi]d\[Lambda]coefffunFEPerPar[nmax_,a_,p_,e_,x_]:=Module[{stepsr,wrlist,\[Delta]v\[Phi]tlist,ExpniTable},
	stepsr=4*nmax;
	wrlist=Table[i,{i,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	\[Delta]v\[Phi]tlist=\[Delta]v\[Phi]funFEPerPar[wrlist,a,p,e,x];

	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e}]],{n,-nmax,nmax},{i,1,stepsr}];
	Chop[(ExpniTable . \[Delta]v\[Phi]tlist)/stepsr,10^(-16)]
]


(* ::Subsection::Closed:: *)
(*Bound velocity*)


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunFEPerPar[wr_,a_,p_,e_,x_]:=Module[{EEg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,dRgdr,rg},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	
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
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	
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


\[Delta]vtfunFEPerPar[wr_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,r3g,rg,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	rg=rgICr1gfun[wr,a,p,e,x];
	
	\[CapitalDelta]=a^2-2rg+rg^2;

	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+rg^2)+a*EEg(a^2+3rg^2))/(rg*\[CapitalDelta])+(rg(2a^2+a^2*rg+rg^3)\[Delta]EE)/\[CapitalDelta]-(2a*rg*\[Delta]Lz)/\[CapitalDelta];
	
	ddtgd\[Lambda]drfun=(EEg*rg(a^2+3rg^2))/\[CapitalDelta]-((rg^2-a^2)(-2a*Lzg+EEg*rg^3+a^2*EEg(2+rg)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunFEPerPar[wr,a,p,e,x]+d\[Delta]td\[Lambda]fun
]


\[Delta]vtfunOfrFEPerPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

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
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

	rg=rgICr1gfun[wr,a,p,e,x];
	
	\[CapitalDelta]=a^2-2rg+rg^2;

	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+rg)))/(rg*\[CapitalDelta])+(2a*rg*\[Delta]EE)/\[CapitalDelta]+((rg-2)rg*\[Delta]Lz)/\[CapitalDelta];
	
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(rg-1)-EEg*rg^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunFEPerPar[wr,a,p,e,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


\[Delta]v\[Phi]funOfrFEPerPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

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
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

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
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

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
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

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
		\[ScriptCapitalI]r1g=2/(r1g(r1g-r2g)) Sqrt[r/(r1g-r)]-\[ScriptCapitalI]/(r1g-r2g);
		\[ScriptCapitalI]r2g=(2Sqrt[r(r1g-r)])/((r-r2g)r2g(r1g-r2g))-((r1g-2r2g)\[ScriptCapitalI])/((r1g-r2g)r2g);
		
		1/Sqrt[1-EEg^2]((4*\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+1/2EEg(\[Delta]r1+\[Delta]r2+\[Delta]r2+2\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+1/2EEg((2+r1g)\[Delta]r1+(2+r2g)\[Delta]r2+(2+r2g)\[Delta]r2-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+1/2dtgd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+1/2dtgd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2-dtgd\[Lambda]fun*\[Delta]rfunFEHomParRes[r,a,p,e,x])
	]
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal correction to the trajectory*)


\[Delta]\[Phi]funFEHomPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,d\[Phi]gd\[Lambda]funr2g,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

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
		\[ScriptCapitalI]r2g=(2Sqrt[r(r1g-r)])/((r-r2g)r2g(r1g-r2g))-((r1g-2r2g)\[ScriptCapitalI])/((r1g-r2g)r2g);

		1/Sqrt[1-EEg^2](-EEg(1-(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+(1/2)d\[Phi]gd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+1/2d\[Phi]gd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2-d\[Phi]gd\[Lambda]fun*\[Delta]rfunFEHomParRes[r,a,p,e,x])
	]
]


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunFEHomPar[r_,a_,p_,e_,x_]:=Module[{EEg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,\[Phi],dRgdr},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	
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
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
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
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
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
	\[Delta]rISCO=\[Delta]rISCOfunFE[a,x];
	\[Delta]r41=\[Delta]r41fun[a,r1g,0,x];

	-((a (r1g-r)^2(r1g+5r))/(3r1g^3*r))+r/r1g*\[Delta]rISCO+(2(r1g-r)^2*\[Delta]r41)/r1g^2
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunISCOplungeParRes[r_,a_,x_]:=Module[{r1g,\[Delta]r41,EEg,\[Delta]rISCO},
	r1g=ISCOradius[a,x];
	EEg=EEgfun[a,r1g,0,1];
	\[Delta]rISCO=\[Delta]rISCOfunFE[a,x];
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
	\[Delta]EE=\[Delta]EEISCOfunFE[a,x];
	\[Delta]Lz=\[Delta]LzISCOfunFE[a,x];
	\[Delta]rISCO=\[Delta]rISCOfunFE[a,x];
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
(*Azimuthal correction to the trajectory*)


\[Delta]\[Phi]funISCOplungePar[r_,a_,x_]:=Module[{r1g,rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]rISCO,d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,\[ScriptCapitalI],\[ScriptCapitalI]r1g},
	r1g=ISCOradius[a,x];
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,r1g,0,1];
	Lzg=Lzgfun[a,r1g,0,1];
	\[Delta]EE=\[Delta]EEISCOfunFE[a,x];
	\[Delta]Lz=\[Delta]LzISCOfunFE[a,x];
	\[Delta]rISCO=\[Delta]rISCOfunFE[a,x];

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
	\[Delta]EE=\[Delta]EEISCOfunFE[a,x];
	\[Delta]rISCO=\[Delta]rISCOfunFE[a,x];
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
	\[Delta]EE=\[Delta]EEISCOfunFE[a,x];
	\[Delta]Lz=\[Delta]LzISCOfunFE[a,x];

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
	\[Delta]EE=\[Delta]EEISCOfunFE[a,x];
	\[Delta]Lz=\[Delta]LzISCOfunFE[a,x];

	\[CapitalDelta]=a^2-2r+r^2;

	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+r)))/(r*\[CapitalDelta])+(2a*r*\[Delta]EE)/\[CapitalDelta]+((r-2)r*\[Delta]Lz)/\[CapitalDelta];

	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(r-1)-EEg*r^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunISCOplungePar[r,a,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Subsection::Closed:: *)
(*Critical plunge *)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunCritplungePar[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,\[Delta]r41,\[Delta]r1,\[Delta]r2,EEg},
	EEg=EEgfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	(r-r2g)Sqrt[r(-r+r1g)](-((2Sqrt[r1g-r])/(r1g*r2g*Sqrt[r]))\[Delta]r41-a/2 (-((2Sqrt[r(-r+r1g)])/(3r*r1g*r2g))(2/r1g+3/r2g+1/r)))+(r(r-r2g)\[Delta]r1)/(r1g(r1g-r2g))+(r(-r+r1g)\[Delta]r2)/((r1g-r2g)r2g)
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunCritplungeParRes[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,\[Delta]r41,\[Delta]r1,\[Delta]r2,EEg,Lzg},
	EEg=EEgfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	 -((2Sqrt[r1g-r])/(r1g*r2g*Sqrt[r]))\[Delta]r41-a/2 (-((2Sqrt[r(-r+r1g)])/(3r*r1g*r2g))(2/r1g+3/r2g+1/r))+ \[Delta]r1/(r1g(r1g-r2g)) Sqrt[r/(-r+r1g)]+(Sqrt[r(-r+r1g)]\[Delta]r2)/((r-r2g)r2g(r1g-r2g))
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory *)


\[Delta]tfunCritplungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,dtgd\[Lambda]fun,dtgd\[Lambda]funr1g,dtgd\[Lambda]funr2g,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

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
	\[ScriptCapitalI]r1g=2/(r1g(r1g-r2g)) Sqrt[r/(r1g-r)]-\[ScriptCapitalI]/(r1g-r2g);
	\[ScriptCapitalI]r2g=(2Sqrt[r(r1g-r)])/((r-r2g)r2g(r1g-r2g))-((r1g-2r2g)\[ScriptCapitalI])/((r1g-r2g)r2g);
	\[ScriptCapitalI]rover=2ArcCos[Sqrt[r/r1g]]+r2g*\[ScriptCapitalI];
	\[ScriptCapitalI]r2over=Sqrt[r(r1g-r)]+(r1g+2r2g)ArcCos[Sqrt[r/r1g]]+r2g^2\[ScriptCapitalI];

	-1/Sqrt[1-EEg^2]((4\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+1/2 EEg(\[Delta]r1+2\[Delta]r2+2\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+1/2*EEg((2+r1g)\[Delta]r1+2(2+r2g)\[Delta]r2-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+1/2dtgd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+1/2dtgd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2-dtgd\[Lambda]fun*\[Delta]rfunCritplungeParRes[r,a,p,e,x])
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal correction to the trajectory*)


\[Delta]\[Phi]funCritplungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,\[Delta]r1,\[Delta]r2,d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,d\[Phi]gd\[Lambda]funr2g,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	
	d\[Phi]gd\[Lambda]fun=((2a*EEg+Lzg(-2+r))r)/(a^2-2r+r^2);
	d\[Phi]gd\[Lambda]funr1g=((2a*EEg+Lzg(-2+r1g))r1g)/(a^2-2r1g+r1g^2);
	d\[Phi]gd\[Lambda]funr2g=((2a*EEg+Lzg(-2+r2g))r2g)/(a^2-2r2g+r2g^2);
	\[ScriptCapitalI]=-1/(2Sqrt[(r1g-r2g)r2g])Log[(Sqrt[r(r1g-r2g)]-Sqrt[(r1g-r)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(r1g-r)r2g])^2];
	\[ScriptCapitalI]r1g=2/(r1g(r1g-r2g)) Sqrt[r/(r1g-r)]-\[ScriptCapitalI]/(r1g-r2g);
	\[ScriptCapitalI]r2g=(2Sqrt[r(r1g-r)])/((r-r2g)r2g(r1g-r2g))-((r1g-2r2g)\[ScriptCapitalI])/((r1g-r2g)r2g);

	-1/Sqrt[1-EEg^2](-EEg(1-(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz \[ScriptCapitalI]+1/2d\[Phi]gd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+1/2d\[Phi]gd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2-d\[Phi]gd\[Lambda]fun*\[Delta]rfunCritplungeParRes[r,a,p,e,x])
]


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunCritplungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,r1g,r2g,\[Delta]r1,\[Delta]r2,\[Delta]r41,dRgdr},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	dRgdr=Sqrt[(1-EEg^2)](r2g-r)(r1g(r2g-3r)-2(r2g-2r)r);
	
	-(r-r2g)Sqrt[(1-EEg^2)r(r1g-r)](\[Delta]r1/(2(r-r1g))+\[Delta]r2/(r-r2g)+\[Delta]r41/r-a/(2r^2)+(EEg*\[Delta]EE)/(1-EEg^2))+1/2dRgdr*\[Delta]rfunCritplungeParRes[r,a,p,e,x]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity*)


\[Delta]vtfunCritplungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

	\[CapitalDelta]=a^2-2r+r^2;

	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+r^2)+a*EEg(a^2+3r^2))/(r*\[CapitalDelta])+(r(2a^2+a^2*r+r^3)\[Delta]EE)/\[CapitalDelta]-(2a*r*\[Delta]Lz)/\[CapitalDelta];

	ddtgd\[Lambda]drfun=(EEg*r(a^2+3r^2))/\[CapitalDelta]-((r^2-a^2)(-2a*Lzg+EEg*r^3+a^2*EEg(2+r)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunCritplungePar[r,a,p,e,x]+d\[Delta]td\[Lambda]fun
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


\[Delta]v\[Phi]funCritplungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

	\[CapitalDelta]=a^2-2r+r^2;

	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+r)))/(r*\[CapitalDelta])+(2a*r*\[Delta]EE)/\[CapitalDelta]+((r-2)r*\[Delta]Lz)/\[CapitalDelta];

	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(r-1)-EEg*r^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunCritplungePar[r,a,p,e,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Subsection::Closed:: *)
(*Generic plunge*)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunPlungePar[r_,a_,p_,e_,x_]:=Module[{EEg,\[Delta]p,\[Delta]r41,r1g,\[Rho]rg,\[Rho]ig,Ar1g,B0,\[Delta]r1,\[Delta]\[Rho]r,\[Phi],krg,ellF,ellE},
	EEg=EEgfunCC[p,e];
	\[Delta]p=\[Delta]pfunCC[a,p,e,x];

	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];

	Ar1g=Sqrt[\[Rho]ig^2+(r1g-\[Rho]rg)^2];
	B0=Sqrt[\[Rho]ig^2+\[Rho]rg^2];
	krg=(-(Ar1g-B0)^2+r1g^2)/(4Ar1g*B0);
	\[Delta]r1=\[Delta]p/(1-e);
	\[Delta]\[Rho]r=\[Delta]p/(1+e);
	\[Delta]r41=\[Delta]r41funCC[a,p,e,x];

	\[Phi]=ArcCos[(-Ar1g*r+B0(r1g-r))/(Ar1g*r+B0(r1g-r))];
	ellE=EllipticE[\[Pi]-\[Phi],krg];
	
	Sqrt[(r1g-r)(\[Rho]ig^2+(r-\[Rho]rg)^2)r](1/Sqrt[Ar1g*B0])(\[Delta]r1/2((2B0*ellE)/(Ar1g*r1g)-((Ar1g*r+B0(r1g-r))(Ar1g*B0+r(r+r1g-2\[Rho]rg))Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(Ar1g*r1g*r(\[Rho]ig^2+(r-\[Rho]rg)^2)))+\[Delta]\[Rho]r((r1g-2\[Rho]rg)/(Ar1g*B0)ellE+(r1g(-2Ar1g*B0*r+(\[Rho]ig^2+(r-\[Rho]rg)^2)(r1g-2\[Rho]rg)+2Ar1g*B0*\[Rho]rg)Sin[\[Phi]])/(2Ar1g*B0(B0^2(-2r+r1g)-r^2(r1g-2\[Rho]rg))Sqrt[1-krg*Sin[\[Phi]]^2]))+\[Delta]r41(-((2Ar1g*ellE)/(B0*r1g))+((Ar1g*r+B0(-r+r1g))Sin[\[Phi]]Sqrt[1-krg Sin[\[Phi]]^2])/(B0*r*r1g))-a/2((2Sqrt[Ar1g*B0]Sqrt[r1g-r])/(3r1g*r Sqrt[r(\[Rho]ig^2+(r-\[Rho]rg)^2)])-(4Ar1g(B0^2+2r1g*\[Rho]rg)ellE)/(3B0^3*r1g^2)+(2(Ar1g*r+B0(-r+r1g))(B0^2+r1g*\[Rho]rg)Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3B0^3*r1g^2*r)-Sin[\[Phi]]/(3B0^3 Sqrt[1-krg*Sin[\[Phi]]^2]) ((-\[Rho]ig^2+\[Rho]rg^2)/B0+(4Ar1g*krg*\[Rho]rg*Cos[\[Phi]])/r1g)))+\[Delta]r1
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunPlungeParRes[r_,a_,p_,e_,x_]:=Module[{EEg,\[Delta]p,\[Delta]r41,r1g,\[Rho]rg,\[Rho]ig,Ar1g,B0,\[Delta]r1,\[Delta]\[Rho]r,\[Phi],krg,ellF,ellE},
	EEg=EEgfunCC[p,e];
	\[Delta]p=\[Delta]pfunCC[a,p,e,x];

	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];

	Ar1g=Sqrt[\[Rho]ig^2+(r1g-\[Rho]rg)^2];
	B0=Sqrt[\[Rho]ig^2+\[Rho]rg^2];
	krg=(-(Ar1g-B0)^2+r1g^2)/(4Ar1g*B0);
	\[Delta]r1=\[Delta]p/(1-e);
	\[Delta]\[Rho]r=\[Delta]p/(1+e);
	\[Delta]r41=\[Delta]r41funCC[a,p,e,x];

	\[Phi]=ArcCos[(-Ar1g*r+B0(r1g-r))/(Ar1g*r+B0(r1g-r))];
	ellE=EllipticE[\[Pi]-\[Phi],krg];

	(1/Sqrt[Ar1g*B0])(\[Delta]r1/2((2B0*ellE)/(Ar1g*r1g)-((Ar1g*r+B0(r1g-r))(Ar1g*B0+r(r+r1g-2\[Rho]rg))Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(Ar1g*r1g*r(\[Rho]ig^2+(r-\[Rho]rg)^2)))+\[Delta]\[Rho]r((r1g-2\[Rho]rg)/(Ar1g*B0)ellE+(r1g(-2Ar1g*B0*r+(\[Rho]ig^2+(r-\[Rho]rg)^2)(r1g-2\[Rho]rg)+2Ar1g*B0*\[Rho]rg)Sin[\[Phi]])/(2Ar1g*B0(B0^2(-2r+r1g)-r^2(r1g-2\[Rho]rg))Sqrt[1-krg*Sin[\[Phi]]^2]))+\[Delta]r41(-((2Ar1g*ellE)/(B0*r1g))+((Ar1g*r+B0(-r+r1g))Sin[\[Phi]]Sqrt[1-krg Sin[\[Phi]]^2])/(B0*r*r1g))-a/2((2Sqrt[Ar1g*B0]Sqrt[r1g-r])/(3r1g*r Sqrt[r(\[Rho]ig^2+(r-\[Rho]rg)^2)])-(4Ar1g(B0^2+2r1g*\[Rho]rg)ellE)/(3B0^3*r1g^2)+(2(Ar1g*r+B0(-r+r1g))(B0^2+r1g*\[Rho]rg)Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3B0^3*r1g^2*r)-Sin[\[Phi]]/(3B0^3 Sqrt[1-krg*Sin[\[Phi]]^2]) ((-\[Rho]ig^2+\[Rho]rg^2)/B0+(4Ar1g*krg*\[Rho]rg*Cos[\[Phi]])/r1g)))+\[Delta]r1/Sqrt[(r1g-r)(\[Rho]ig^2+(r-\[Rho]rg)^2)r]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory*)


\[Delta]tfunPlungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]EE,\[Delta]Lz,r1g,r2g,r2gCC,\[Rho]rg,\[Rho]ig,\[Delta]r1,\[Delta]\[Rho]r,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,dtgd\[Lambda]fun,dtgd\[Lambda]funr1g,dtgd\[Lambda]funr2g,dtgd\[Lambda]funr2gCC,Ar1g,B0,krg,\[Gamma]r,\[Phi],ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over,\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r2gCC},
	EEg=EEgfunCC[p,e];
	Lzg=LzgfunCC[a,p,e,x];
	\[Delta]p=\[Delta]pfunCC[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunCC[a,p,e,x]+dEEdpfunCC[p,e]\[Delta]p;
	\[Delta]Lz=\[Delta]LzfunCC[a,p,e,x]+dLzdpfunCC[a,p,e,x]\[Delta]p;

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
	\[Delta]\[Rho]r4=\[Delta]r41funCC[a,p,e,x];
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

	1/Sqrt[1-EEg^2](((4\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+1/2dtgd\[Lambda]funr1g*\[Delta]r1*\[ScriptCapitalI]r1g+1/2dtgd\[Lambda]funr2gCC*\[ScriptCapitalI]r2gCC*\[Delta]\[Rho]r+1/2dtgd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]\[Rho]r+1/2EEg((2+r1g)\[Delta]r1+2(2+\[Rho]rg)\[Delta]\[Rho]r-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+\[ScriptCapitalI]rover((2\[Delta]EE)/(1-EEg^2)+(EEg*\[Delta]r1)/2+EEg*\[Delta]\[Rho]r+EEg*\[Delta]\[Rho]r4)+dtgd\[Lambda]funr1g*\[Delta]r1(1/Sqrt[r(r1g-r)(r-I*\[Rho]ig-\[Rho]rg)(r+I*\[Rho]ig-\[Rho]rg)])-dtgd\[Lambda]fun*\[Delta]rfunPlungeParRes[r,a,p,e,x])
]	


(* ::Subsubsection::Closed:: *)
(*Azimuthal trajectory*)


\[Delta]\[Phi]funPlungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]EE,\[Delta]Lz,r1g,r2g,r2gCC,\[Rho]rg,\[Rho]ig,\[Delta]r1,\[Delta]\[Rho]r,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,d\[Phi]gd\[Lambda]funr2g,d\[Phi]gd\[Lambda]funr2gCC,Ar1g,B0,krg,\[Gamma]r,\[Phi],ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r2gCC},
	EEg=EEgfunCC[p,e];
	Lzg=LzgfunCC[a,p,e,x];
	\[Delta]p=\[Delta]pfunCC[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunCC[a,p,e,x]+dEEdpfunCC[p,e]\[Delta]p;
	\[Delta]Lz=\[Delta]LzfunCC[a,p,e,x]+dLzdpfunCC[a,p,e,x]\[Delta]p;

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
	
	1/Sqrt[1-EEg^2]((-EEg+(EEg*Lzg*\[Delta]EE)/(1-EEg^2)+\[Delta]Lz)\[ScriptCapitalI]+1/2d\[Phi]gd\[Lambda]funr1g*\[Delta]r1*\[ScriptCapitalI]r1g+1/2\[ScriptCapitalI]r2g*d\[Phi]gd\[Lambda]funr2g*\[Delta]\[Rho]r+1/2\[ScriptCapitalI]r2gCC*d\[Phi]gd\[Lambda]funr2gCC*\[Delta]\[Rho]r+d\[Phi]gd\[Lambda]funr1g*\[Delta]r1(1/(Sqrt[r(r1g-r)(r-I*\[Rho]ig-\[Rho]rg)(r+I*\[Rho]ig-\[Rho]rg)]))-d\[Phi]gd\[Lambda]fun*\[Delta]rfunPlungeParRes[r,a,p,e,x])
]


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunPlungePar[r_,a_,p_,e_,x_]:=Module[{EEg,\[Delta]p,\[Delta]EE,\[Delta]r41,r1g,\[Rho]rg,\[Rho]ig,Ar1g,B0,\[Delta]r1,\[Delta]\[Rho]r,dRgdr},
	EEg=EEgfunCC[p,e];
	\[Delta]p=\[Delta]pfunCC[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunCC[a,p,e,x]+dEEdpfunCC[p,e]\[Delta]p;

	r1g=p/(1-e);
	\[Rho]rg=p/(1+e);
	\[Rho]ig=\[Rho]igfun[a,p,e,x];
	\[Delta]r1=\[Delta]p/(1-e);
	\[Delta]\[Rho]r=\[Delta]p/(1+e);
	\[Delta]r41=\[Delta]r41funCC[a,p,e,x];
	dRgdr=Sqrt[1-EEg^2](-4r^3+3r^2(r1g+2\[Rho]rg)+r1g(\[Rho]ig^2+\[Rho]rg^2)-2r(\[Rho]ig^2+\[Rho]rg(2 r1g+\[Rho]rg)));

	-Sqrt[(1-EEg^2)(r1g-r)(\[Rho]ig^2+(r-\[Rho]rg)^2)r](\[Delta]r1/(2(r-r1g))+(\[Delta]\[Rho]r(r-\[Rho]rg))/(r^2+\[Rho]ig^2-2r*\[Rho]rg+\[Rho]rg^2)+\[Delta]r41/r-a/(2r^2)+(EEg*\[Delta]EE)/(1-EEg^2))+1/2dRgdr*\[Delta]rfunPlungeParRes[r,a,p,e,x]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity*)


\[Delta]vtfunPlungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	EEg=EEgfunCC[p,e];
	Lzg=LzgfunCC[a,p,e,x];
	\[Delta]p=\[Delta]pfunCC[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunCC[a,p,e,x]+dEEdpfunCC[p,e]\[Delta]p;
	\[Delta]Lz=\[Delta]LzfunCC[a,p,e,x]+dLzdpfunCC[a,p,e,x]\[Delta]p;
	
	\[CapitalDelta]=a^2-2r+r^2;
	
	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+r^2)+a*EEg(a^2+3r^2))/(r*\[CapitalDelta])+(r(2a^2+a^2*r+r^3)\[Delta]EE)/\[CapitalDelta]-(2a*r*\[Delta]Lz)/\[CapitalDelta];
	
	ddtgd\[Lambda]drfun=(EEg*r(a^2+3r^2))/\[CapitalDelta]-((r^2-a^2)(-2a*Lzg+EEg*r^3+a^2*EEg(2+r)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunPlungePar[r,a,p,e,x]+d\[Delta]td\[Lambda]fun
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


\[Delta]v\[Phi]funPlungePar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]p,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	EEg=EEgfunCC[p,e];
	Lzg=LzgfunCC[a,p,e,x];
	\[Delta]p=\[Delta]pfunCC[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunCC[a,p,e,x]+dEEdpfunCC[p,e]\[Delta]p;
	\[Delta]Lz=\[Delta]LzfunCC[a,p,e,x]+dLzdpfunCC[a,p,e,x]\[Delta]p;

	\[CapitalDelta]=a^2-2r+r^2;

	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+r)))/(r*\[CapitalDelta])+(2a*r*\[Delta]EE)/\[CapitalDelta]+((r-2)r*\[Delta]Lz)/\[CapitalDelta];

	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(r-1)-EEg*r^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunPlungePar[r,a,p,e,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Subsection::Closed:: *)
(*Plunge related to bound orbits*)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunPlungeBoundOrbitPar[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,EEg,\[Delta]EE,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],krg,ellF,ellE},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
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


\[Delta]rfunPlungeBoundOrbitDoubleRootPar[r_,a_,p_,x_]:=Module[{r1g,r3g,\[Delta]r1,\[Delta]r41,\[Delta]EE,EEg,\[Phi]},
	EEg=EEgfun[a,p,0,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,0,x]+dEEdpfun[a,p,0,x]\[Delta]pfun[a,p,0,x];

	r1g=p;
	r3g=2/(1-EEg^2)-2p;
	\[Delta]r1=\[Delta]r1funFE[a,p,0,x];
	\[Delta]r41=\[Delta]r41fun[a,p,0,x];

	\[Phi]=ArcSin[Sqrt[(r1g (r3g-r))/(r3g (r1g-r))]];
	
	(r1g-r)Sqrt[(r3g-r)r]((-((r3g*Cos[\[Phi]]*Sin[\[Phi]])/(r1g(r1g-r3g))^(3/2)))\[Delta]r1+((2(r1g-r3g)Tan[\[Phi]])/(r1g*r3g Sqrt[r1g(r1g-r3g)]))\[Delta]r41-a/2(-((2(r^2+r*r3g+r3g^2-r1g(r+r3g)))/(3(r(r1g-r3g)r3g^2)(r1g-r))) Sqrt[r3g/r-1]+2/3 1/r1g^2 1/Sqrt[r1g(r1g-r3g)] ((r3g*Cos[\[Phi]]*Sin[\[Phi]])/ (r1g-r3g)+((r1g-r3g)(r1g+2r3g)Tan[\[Phi]])/ r3g^2)))+1/2 ((2(r1g-r)r)/((r1g-r3g)r3g) \[Delta]r1)
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunPlungeBoundOrbitParRes[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]KK,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],krg,\[Gamma]rg,ellF,ellE},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
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


\[Delta]rfunPlungeBoundOrbitDoubleRootParRes[r_,a_,p_,x_]:=Module[{r1g,r3g,\[Delta]r1,\[Delta]r41,\[Delta]EE,EEg,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi]},
	EEg=EEgfun[a,p,0,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,0,x]+dEEdpfun[a,p,0,x]\[Delta]pfun[a,p,0,x];

	r1g=p;
	r3g=2/(1-EEg^2)-2p;
	\[Delta]r1=\[Delta]r1funFE[a,p,0,x];
	\[Delta]r41=\[Delta]r41fun[a,p,0,x];

	\[Phi]=ArcSin[Sqrt[(r1g (r3g-r))/(r3g (r1g-r))]];
	
	(-((r3g*Cos[\[Phi]]*Sin[\[Phi]])/(r1g(r1g-r3g))^(3/2)))\[Delta]r1+((2(r1g-r3g)Tan[\[Phi]])/(r1g*r3g Sqrt[r1g(r1g-r3g)]))\[Delta]r41-a/2(-((2(r^2+r*r3g+r3g^2-r1g(r+r3g)))/(3(r(r1g-r3g)r3g^2)(r1g-r))) Sqrt[r3g/r-1]+2/3 1/r1g^2 1/Sqrt[r1g(r1g-r3g)] ((r3g*Cos[\[Phi]]*Sin[\[Phi]])/ (r1g-r3g)+((r1g-r3g)(r1g+2r3g)Tan[\[Phi]])/ r3g^2))+1/2((2Sqrt[(r1g-r)r])/((r1g-r3g)*r3g Sqrt[(r1g-r)(r3g-r)])\[Delta]r1)
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory*)


\[Delta]tfunPlungeBoundOrbitPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Phi],\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,dtgd\[Lambda]fun,dtgd\[Lambda]funr1g,dtgd\[Lambda]funr2g,dtgd\[Lambda]funr3g,krg,\[Gamma]r,ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

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
	ellE= EllipticE[\[Phi],krg];
	ellPi=EllipticPi[\[Gamma]r,\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[(r2g(r1g-r3g))];
	\[ScriptCapitalI]r1g=(2/Sqrt[(r2g(r1g-r3g))])(1/(r1g-r2g))((r2g*ellE)/r1g-ellF-((r1g-r2g)r3g*Cos[\[Phi]]*Sin[\[Phi]])/(r1g(r1g-r3g)Sqrt[(1-krg*Sin[\[Phi]]^2)]));
	\[ScriptCapitalI]r2g=(2/Sqrt[(r2g(r1g-r3g))])(1/(r1g-r2g))(-(((r1g-r3g)ellE)/(r2g-r3g))+ellF);
	\[ScriptCapitalI]r3g=(2/Sqrt[(r2g(r1g-r3g))]) 1/r3g ((r2g*ellE)/(r2g-r3g)-ellF);
	\[ScriptCapitalI]rover=2 /Sqrt[r2g(r1g-r3g)] (r2g*ellF+(-r2g+r3g)ellPi);
	\[ScriptCapitalI]r2over=-r3g Sqrt[((r1g-r)(r2g-r))/(r2g(r2g-r3g))]Sin[\[Phi]]Cos[\[Phi]]+Sqrt[r2g(r1g-r3g)]ellE-((r2g-r3g)(r1g+r2g+r3g))/Sqrt[r2g(r1g-r3g)] ellPi+(r2g(r2g+r3g))/Sqrt[r2g(r1g-r3g)] ellF;

	-1/Sqrt[1-EEg^2]((4\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+1/2 EEg(\[Delta]r1+\[Delta]r2+\[Delta]r3+2\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+\[Delta]EE/(1-EEg^2)\[ScriptCapitalI]r2over+1/2 EEg((2+r1g)\[Delta]r1+(2+r2g)\[Delta]r2+(2+r3g)\[Delta]r3-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+1/2dtgd\[Lambda]funr1g*\[ScriptCapitalI]r1g*\[Delta]r1+1/2dtgd\[Lambda]funr2g*\[ScriptCapitalI]r2g*\[Delta]r2+1/2dtgd\[Lambda]funr3g*\[ScriptCapitalI]r3g*\[Delta]r3+1/2 dtgd\[Lambda]funr3g((2/(r1g-r3g))1/r3g((\[Sqrt]((r1g-r)r))/(\[Sqrt]((r2g-r)(r3g-r)))))\[Delta]r3-dtgd\[Lambda]fun*\[Delta]rfunPlungeBoundOrbitParRes[r,a,p,e,x])
]


\[Delta]tfunPlungeBoundOrbitDoubleRootPar[r_,a_,p_,x_]:=Module[{EEg,Lzg,\[Delta]EE,r1g,r3g,\[Delta]r1,\[Phi],\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,dtgd\[Lambda]fun,dtgd\[Lambda]funr1g,dtgd\[Lambda]funr2g,dtgd\[Lambda]funr3g,\[ScriptCapitalI],\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
	EEg=EEgfun[a,p,0,x];
	Lzg=Lzgfun[a,p,0,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,0,x]+dEEdpfun[a,p,0,x]\[Delta]pfun[a,p,0,x];

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

	-(1/Sqrt[1-EEg^2])((4\[Delta]EE)/(1-EEg^2) \[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+1/2 EEg(3\[Delta]r1+2\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+\[Delta]EE/(1-EEg^2) \[ScriptCapitalI]r2over+1/2 EEg(2(2+r1g)\[Delta]r1+(2+r3g)\[Delta]r1-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+(dtgd\[Lambda]funr1g(-2r1g+r3g)\[Phi])/(r1g(r1g-r3g))^(3/2) \[Delta]r1+1/2 dtgd\[Lambda]funr3g*\[ScriptCapitalI]r3g*\[Delta]r1+1/2 dtgd\[Lambda]funr3g((2Sqrt[r])/((r1g-r3g)r3g Sqrt[-r+r3g]))\[Delta]r1-dtgd\[Lambda]fun*\[Delta]rfunPlungeBoundOrbitDoubleRootParRes[r,a,p,x])
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal correction to the trajectory*)


\[Delta]\[Phi]funPlungeBoundOrbitPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Phi],d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,d\[Phi]gd\[Lambda]funr2g,d\[Phi]gd\[Lambda]funr3g,krg,\[Gamma]r,ellF,ellE,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

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
	ellE= EllipticE[\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[(r2g(r1g-r3g))];
	\[ScriptCapitalI]r1g=(2/Sqrt[(r2g(r1g-r3g))])(1/(r1g-r2g))((r2g*ellE)/r1g-ellF-((r1g-r2g)r3g*Cos[\[Phi]]*Sin[\[Phi]])/(r1g(r1g-r3g)Sqrt[(1-krg*Sin[\[Phi]]^2)]));
	\[ScriptCapitalI]r2g=(2/Sqrt[(r2g(r1g-r3g))])(1/(r1g-r2g))(-(((r1g-r3g)ellE)/(r2g-r3g))+ellF);
	\[ScriptCapitalI]r3g=(2/Sqrt[(r2g(r1g-r3g))]) 1/r3g ((r2g*ellE)/(r2g-r3g)-ellF);

	-1/Sqrt[1-EEg^2](-EEg(1-(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+1/2 d\[Phi]gd\[Lambda]funr1g \[ScriptCapitalI]r1g+1/2 d\[Phi]gd\[Lambda]funr2g \[ScriptCapitalI]r2g+1/2 d\[Phi]gd\[Lambda]funr3g \[ScriptCapitalI]r3g+1/2 d\[Phi]gd\[Lambda]funr3g*\[Delta]r3((2/(r1g-r3g))1/r3g((\[Sqrt]((r1g-r)r))/(\[Sqrt]((r2g-r)(r3g-r)))))-d\[Phi]gd\[Lambda]fun*\[Delta]rfunPlungeBoundOrbitParRes[r,a,p,e,x])
]


\[Delta]\[Phi]funPlungeBoundOrbitDoubleRootPar[r_,a_,p_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r3g,\[Delta]r1,\[Phi],d\[Phi]gd\[Lambda]fun,d\[Phi]gd\[Lambda]funr1g,d\[Phi]gd\[Lambda]funr2g,d\[Phi]gd\[Lambda]funr3g,\[ScriptCapitalI],\[ScriptCapitalI]r3g},
	EEg=EEgfun[a,p,0,x];
	Lzg=Lzgfun[a,p,0,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,0,x]+dEEdpfun[a,p,0,x]\[Delta]pfun[a,p,0,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,0,x]+dLzdpfun[a,p,0,x]\[Delta]pfun[a,p,0,x];

	r1g=p;
	r3g=2/(1-EEg^2)-2p;
	\[Delta]r1=\[Delta]r1funFE[a,p,0,x];

	\[Phi]=ArcSin[Sqrt[(r1g(r3g-r))/(r3g(r1g-r))]];

	d\[Phi]gd\[Lambda]fun=((2a*EEg+Lzg(-2+r))r)/(a^2-2r+r^2);
	d\[Phi]gd\[Lambda]funr1g=((2a*EEg+Lzg(-2+r1g))r1g)/(a^2-2r1g+r1g^2);
	d\[Phi]gd\[Lambda]funr3g=((2 a EEg+Lzg (-2+r3g))r3g)/(a^2-2r3g+r3g^2);
	
	\[ScriptCapitalI]=(2\[Phi])/Sqrt[r1g(r1g-r3g)];
	\[ScriptCapitalI]r3g=2 /(Sqrt[r1g(r1g-r3g)]r3g) ((r3g*\[Phi])/(r1g-r3g));

	-(1/Sqrt[1-EEg^2])(-EEg(1-(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+(d\[Phi]gd\[Lambda]funr1g(-2r1g+r3g)\[Phi])/(r1g(r1g-r3g))^(3/2) \[Delta]r1+1/2 d\[Phi]gd\[Lambda]funr3g*\[ScriptCapitalI]r3g*\[Delta]r1+1/2 d\[Phi]gd\[Lambda]funr3g*\[Delta]r1((2Sqrt[r])/((r1g-r3g)r3g Sqrt[-r+r3g]))-d\[Phi]gd\[Lambda]fun*\[Delta]rfunPlungeBoundOrbitDoubleRootParRes[r,a,p,x])
]


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunPlungeBoundOrbitPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,dRgdr},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	dRgdr=Sqrt[(1-EEg^2)](-4r^3+r1g*r2g*r3g+3r^2(r1g+r2g+r3g)-2r(r2g*r3g+r1g(r2g+r3g)));

	Sqrt[(1-EEg^2)r(r1g-r)(r-r2g)(r-r3g)](\[Delta]r1/(2(r-r1g))+\[Delta]r2/(2(r-r2g))+\[Delta]r3/(2(r-r3g))+\[Delta]r41/r-a/(2r^2)+(EEg*\[Delta]EE)/(1-EEg^2))-1/2dRgdr*\[Delta]rfunPlungeBoundOrbitParRes[r,a,p,e,x]
]


\[Delta]vrfunPlungeBoundOrbitDoubleRootPar[r_,a_,p_,x_]:=Module[{EEg,\[Delta]EE,r1g,r3g,\[Delta]r1,\[Delta]r41,dRgdr},
	EEg=EEgfun[a,p,0,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,0,x]+dEEdpfun[a,p,0,x]\[Delta]pfun[a,p,0,x];
	\[Delta]r1=\[Delta]r1funFE[a,p,0,x];
	\[Delta]r41=\[Delta]r41fun[a,p,0,x];

	r1g=p;
	r3g=2/(1-EEg^2)-2p;
	dRgdr=Sqrt[(1-EEg^2)](-4r^3+r1g^2*r3g+3r^2(2r1g+r3g)-2r(r1g*r3g+r1g(r1g+r3g)));

	(r1g-r)Sqrt[(1-EEg^2)r(r3g-r)](\[Delta]r1/(r-r1g)+\[Delta]r1/(2(r-r3g))+\[Delta]r41/r-a/(2 r^2)+(EEg*\[Delta]EE)/(1-EEg^2))-1/2dRgdr*\[Delta]rfunPlungeBoundOrbitDoubleRootParRes[r,a,p,x]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity*)


\[Delta]vtfunPlungeBoundOrbitPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]td\[Lambda]fun,d\[Delta]td\[Lambda]fun2,ddtgd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

	\[CapitalDelta]=a^2-2r+r^2;
	
	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+r^2)+a*EEg(a^2+3r^2))/(r*\[CapitalDelta])+(r(2a^2+a^2*r+r^3)\[Delta]EE)/\[CapitalDelta]-(2a*r*\[Delta]Lz)/\[CapitalDelta];
	
	ddtgd\[Lambda]drfun=(EEg*r(a^2+3r^2))/\[CapitalDelta]-((r^2-a^2)(-2a*Lzg+EEg*r^3+a^2*EEg(2+r)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunPlungeBoundOrbitPar[r,a,p,e,x]+d\[Delta]td\[Lambda]fun
]


\[Delta]vtfunPlungeBoundOrbitDoubleRootPar[r_,a_,p_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	EEg=EEgfun[a,p,0,x];
	Lzg=Lzgfun[a,p,0,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,0,x]+dEEdpfun[a,p,0,x]\[Delta]pfun[a,p,0,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,0,x]+dLzdpfun[a,p,0,x]\[Delta]pfun[a,p,0,x];

	\[CapitalDelta]=a^2-2r+r^2;
	
	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+r^2)+a*EEg(a^2+3r^2))/(r*\[CapitalDelta])+(r(2a^2+a^2*r+r^3)\[Delta]EE)/\[CapitalDelta]-(2a*r*\[Delta]Lz)/\[CapitalDelta];
	
	ddtgd\[Lambda]drfun=(EEg*r(a^2+3r^2))/\[CapitalDelta]-((r^2-a^2)(-2a*Lzg+EEg*r^3+a^2*EEg(2+r)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunPlungeBoundOrbitDoubleRootPar[r,a,p,x]+d\[Delta]td\[Lambda]fun
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


\[Delta]v\[Phi]funPlungeBoundOrbitPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

	\[CapitalDelta]=a^2-2r+r^2;
	
	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+r)))/(r*\[CapitalDelta])+(2a*r*\[Delta]EE)/\[CapitalDelta]+((r-2)r*\[Delta]Lz)/\[CapitalDelta];
	
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(r-1)-EEg*r^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunPlungeBoundOrbitPar[r,a,p,e,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


\[Delta]v\[Phi]funPlungeBoundOrbitDoubleRootPar[r_,a_,p_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	EEg=EEgfun[a,p,0,x];
	Lzg=Lzgfun[a,p,0,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,0,x]+dEEdpfun[a,p,0,x]\[Delta]pfun[a,p,0,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,0,x]+dLzdpfun[a,p,0,x]\[Delta]pfun[a,p,0,x];

	\[CapitalDelta]=a^2-2r+r^2;
	
	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+r)))/(r*\[CapitalDelta])+(2a*r*\[Delta]EE)/\[CapitalDelta]+((r-2)r*\[Delta]Lz)/\[CapitalDelta];
	
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(r-1)-EEg*r^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunPlungeBoundOrbitDoubleRootPar[r,a,p,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Section:: *)
(*Spin corrections to the orbit - orthogonal component of the spin *)


(* ::Subsection:: *)
(*Spin precession phase*)


(* ::Subsubsection:: *)
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
	1/(2\[Pi]) Re[\[Psi]freq]
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
	EEg=EEgfunCC[p,e];
	Lzg=LzgfunCC[a,p,e,x];
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
	EEg=EEgfunCC[p,e];
	Lzg=LzgfunCC[a,p,e,x];
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
	
<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->\[Delta]EE,
	 "Jzs"->\[Delta]Lz,
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 "MinoFrequenciesCorrection"->{\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]\[Phi]s},
	 "BLFrequenciesCorrection"->{\[CapitalUpsilon]rs/\[CapitalUpsilon]tg-\[CapitalUpsilon]rg/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts,\[CapitalUpsilon]\[Phi]s/\[CapitalUpsilon]tg-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts},
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
	 "\[CapitalDelta]\[Delta]tpar"->Function[{wr},\[Delta]tfunFTPerPar[wr,a,p,e,x]],
	 "\[Delta]rpar"->Function[{wr},\[Delta]rfunFTPerPar[wr,a,p,e,x]],
	 "\[Psi]p"->Function[{wr},\[Psi]pICr2g[wr,a,p,e,x]],
	 "\[Delta]zort"->Function[{wr},\[Delta]zfunICr2gPer[wr,a,p,e,x]],
	 "\[CapitalDelta]\[Delta]\[Phi]par"->Function[{wr},\[Delta]\[Phi]funFTPerPar[wr,a,p,e,x]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{wr},\[Delta]vtfunFTPerPar[wr,a,p,e,x]],
	 "\[Delta]vrpar"->Function[{wr},\[Delta]vrfunFTPerPar[wr,a,p,e,x]],
	 "\[Delta]v\[Phi]par"->Function[{wr},\[Delta]v\[Phi]funFTPerPar[wr,a,p,e,x]]
	|>
]


(* ::Subsection::Closed:: *)
(*Fourier series expansion*)


KerrNearEqSpinOrbitCorrFTPerFourier[a_, p_, e_, x_, nmax_?(IntegerQ[#] && # > 0&)]:=Module[{EEg,Lzg,\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g,\[CapitalUpsilon]p,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]zs,\[CapitalUpsilon]\[Phi]s,stepsr,\[CapitalDelta]tspar,\[CapitalDelta]\[Phi]spar,\[Delta]rpar,\[Psi]phase,\[Delta]zort,
	ExpniTable,wrlist,dtrgd\[Lambda]coeff,d\[Phi]rgd\[Lambda]coeff,dtspard\[Lambda]coeff,d\[Phi]spard\[Lambda]coeff,dtrgd\[Lambda],d\[Phi]rgd\[Lambda],drgd\[Lambda],rg,\[CapitalDelta]trg,\[CapitalDelta]\[Phi]rg,dtspard\[Lambda],drspard\[Lambda],dzspard\[Lambda],d\[Phi]spard\[Lambda]},
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
	
	(* steps for numerical integration *)
	stepsr=4*nmax;
	(* matrices of discrete Fourier transform *)
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,x}]],{n,-nmax,nmax},{i,1,stepsr}];
	
	wrlist=Table[wr,{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
  
	(* Fourier coeffients geodesic functions *)
	dtrgd\[Lambda]coeff=VtrgcoeffICr2g[nmax,a,p,e,x];
	d\[Phi]rgd\[Lambda]coeff=V\[Phi]rgcoeffICr2g[nmax,a,p,e,x];
	
	dtrgd\[Lambda][wr_]:=FourierVel[wr,dtspard\[Lambda]coeff];
	d\[Phi]rgd\[Lambda][wr_]:=FourierVel[wr,d\[Phi]spard\[Lambda]coeff];
	drgd\[Lambda][wr_]:=drgd\[Lambda]fun[wr,a,p,e,x];

	\[CapitalDelta]trg[wr_]:=\[CapitalDelta]trgIntVel[wr,\[CapitalUpsilon]rg,dtrgd\[Lambda]coeff];
	\[CapitalDelta]\[Phi]rg[wr_]:=\[CapitalDelta]\[Phi]rgIntVel[wr,\[CapitalUpsilon]rg,d\[Phi]rgd\[Lambda]coeff];
	rg[wr_]:=rgICr2gfun[wr,a,p,e,x];
	
	(* spin correction radial and polar frequencies *)
	\[CapitalUpsilon]rs=\[Delta]\[CapitalUpsilon]rfunFT[a,p,e,x];
	
	Print["Calculating Fourier coefficients of \!\(\*SubscriptBox[\(dt\), \(s\)]\)/d\[Lambda]"];
	dtspard\[Lambda]coeff=d\[Delta]td\[Lambda]coefffunFTPerPar[nmax,a,p,e,x];
	\[CapitalUpsilon]ts=dtspard\[Lambda]coeff[[nmax+1]];
	Print["Calculating Fourier coefficients of \!\(\*SubscriptBox[\(d\[Phi]\), \(s\)]\)/d\[Lambda]"];
	d\[Phi]spard\[Lambda]coeff=d\[Delta]\[Phi]d\[Lambda]coefffunFTPerPar[nmax,a,p,e,x];
	\[CapitalUpsilon]\[Phi]s=d\[Phi]spard\[Lambda]coeff[[nmax+1]];
	
	\[CapitalDelta]tspar[wr_]:=\[CapitalDelta]\[Delta]IntvelPerPar[wr,\[CapitalUpsilon]rg,\[CapitalUpsilon]rs,dtspard\[Lambda]coeff,dtrgd\[Lambda]coeff];
	\[CapitalDelta]\[Phi]spar[wr_]:=\[CapitalDelta]\[Delta]IntvelPerPar[wr,\[CapitalUpsilon]rg,\[CapitalUpsilon]rs,d\[Phi]spard\[Lambda]coeff,d\[Phi]rgd\[Lambda]coeff];
	\[Delta]rpar[wr_]:=\[Delta]rfunFTPerPar[wr,a,p,e,x];
	\[Psi]phase[wr_]:=\[Psi]pICr2g[wr,a,p,e,x];
	\[Delta]zort[wr_]:=\[Delta]zfunICr2gPer[wr,a,p,e,x];
	      
	dtspard\[Lambda][wr_]:=FourierVel[wr,dtspard\[Lambda]coeff];
	d\[Phi]spard\[Lambda][wr_]:=FourierVel[wr,d\[Phi]spard\[Lambda]coeff];
	drspard\[Lambda][wr_]:=\[Delta]vrfunFTPerPar[wr,a,p,e,x];
			
	<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->\[Delta]EE,
	 "Jzs"->\[Delta]Lz,
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 "MinoFrequenciesCorrection"->{\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]\[Phi]s},
	 "BLFrequenciesCorrection"->{\[CapitalUpsilon]rs/\[CapitalUpsilon]tg-\[CapitalUpsilon]rg/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts,\[CapitalUpsilon]\[Phi]s/\[CapitalUpsilon]tg-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts},
	 "MinoPrecessionFrequency"->\[CapitalUpsilon]p,
	 "BLPrecessionFrequency"->\[CapitalUpsilon]p/\[CapitalUpsilon]tg,
	 (*Keys purely oscillatory part geodesic coordinate time and azimuthal trajectory*)
	 "\[CapitalDelta]trg"->Function[{wr},\[CapitalDelta]trg[wr]],
	 "\[CapitalDelta]\[Phi]rg"->Function[{wr},\[CapitalDelta]\[Phi]rg[wr]],
	 (*Keys geodesic radial trajectory*)
	 "rg"->Function[{wr},rg[wr]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{wr},dtrgd\[Lambda][wr]],
	 "vrg"->Function[{wr},drgd\[Lambda][wr]],
	 "v\[Phi]g"->Function[{wr},d\[Phi]rgd\[Lambda][wr]],
	 (*Keys corrections trajectory*)
	 "\[CapitalDelta]\[Delta]tpar"->Function[{wr},\[CapitalDelta]tspar[wr]],
	 "\[Delta]rpar"->Function[{wr},\[Delta]rpar[wr]],
	 "\[Psi]p"->Function[wr,\[Psi]phase[wr]],
	 "\[Delta]zort"->Function[{wr},\[Delta]zort[wr]],
	 "\[CapitalDelta]\[Delta]\[Phi]par"->Function[{wr},\[CapitalDelta]\[Phi]spar[wr]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{wr},dtspard\[Lambda][wr]],
	 "\[Delta]vrpar"->Function[{wr},drspard\[Lambda][wr]],
	 "\[Delta]v\[Phi]par"->Function[{wr},d\[Phi]spard\[Lambda][wr]]
	|>
]


(* ::Section::Closed:: *)
(*Near equatorial orbits - fixed constants of motion*)


(* ::Subsection::Closed:: *)
(*Analytic solutions*)


(* ::Subsubsection::Closed:: *)
(*Bound solutions*)


KerrNearEqSpinOrbitCorrFCPer[a_, p_, e_, x_]:=Module[{EEg,Lzg,\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g,\[CapitalUpsilon]p,\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]zs,\[CapitalUpsilon]\[Phi]s},
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	
	(* geodesic frequencies *)
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,x];
	   
	(* spin correction radial and polar frequencies *)
	\[CapitalUpsilon]rs=\[Delta]\[CapitalUpsilon]rfunFC[a,p,e,x];
	\[CapitalUpsilon]ts=\[Delta]\[CapitalUpsilon]tfunFC[a,p,e,x];
	\[CapitalUpsilon]\[Phi]s=\[Delta]\[CapitalUpsilon]\[Phi]funFC[a,p,e,x];
			
<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 "MinoFrequenciesCorrection"->{\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]\[Phi]s},
	 "BLFrequenciesCorrection"->{\[CapitalUpsilon]rs/\[CapitalUpsilon]tg-\[CapitalUpsilon]rg/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts,\[CapitalUpsilon]\[Phi]s/\[CapitalUpsilon]tg-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts},
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
	 "\[CapitalDelta]\[Delta]tpar"->Function[{wr},\[Delta]tfunFCPerPar[wr,a,p,e,x]],
	 "\[Delta]rpar"->Function[{wr},\[Delta]rfunFCPerPar[wr,a,p,e,x]],
	 "\[Psi]p"->Function[{wr},\[Psi]pICr2g[wr,a,p,e,x]],
	 "\[Delta]zort"->Function[{wr},\[Delta]zfunICr2gPer[wr,a,p,e,x]],
	 "\[CapitalDelta]\[Delta]\[Phi]par"->Function[{wr},\[Delta]\[Phi]funFCPerPar[wr,a,p,e,x]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{wr},\[Delta]vtfunFCPerPar[wr,a,p,e,x]],
	 "\[Delta]vrpar"->Function[{wr},\[Delta]vrfunFCPerPar[wr,a,p,e,x]],
	 "\[Delta]v\[Phi]par"->Function[{wr},\[Delta]v\[Phi]funFCPerPar[wr,a,p,e,x]]
	|>
]


(* ::Subsection::Closed:: *)
(*Fourier series expansion*)


KerrNearEqSpinOrbitCorrFCPerFourier[a_, p_, e_, x_, nmax_?(IntegerQ[#] && # > 0&)]:=Module[{EEg,Lzg,\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g,\[CapitalUpsilon]p,\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]zs,\[CapitalUpsilon]\[Phi]s,stepsr,\[CapitalDelta]tspar,\[CapitalDelta]\[Phi]spar,\[Delta]rpar,\[Psi]phase,\[Delta]zort,
	ExpniTable,wrlist,dtrgd\[Lambda]coeff,d\[Phi]rgd\[Lambda]coeff,dtspard\[Lambda]coeff,d\[Phi]spard\[Lambda]coeff,dtrgd\[Lambda],d\[Phi]rgd\[Lambda],drgd\[Lambda],rg,\[CapitalDelta]trg,\[CapitalDelta]\[Phi]rg,dtspard\[Lambda],drspard\[Lambda],dzspard\[Lambda],d\[Phi]spard\[Lambda]},
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	
	(* geodesic frequencies *)
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,x];
	
	(* steps for numerical integration *)
	stepsr=4*nmax;
	(* matrices of discrete Fourier transform *)
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,x}]],{n,-nmax,nmax},{i,1,stepsr}];
	
	wrlist=Table[wr,{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
  
	(* Fourier coeffients geodesic functions *)
	dtrgd\[Lambda]coeff=VtrgcoeffICr2g[nmax,a,p,e,x];
	d\[Phi]rgd\[Lambda]coeff=V\[Phi]rgcoeffICr2g[nmax,a,p,e,x];
	
	dtrgd\[Lambda][wr_]:=FourierVel[wr,dtspard\[Lambda]coeff];
	d\[Phi]rgd\[Lambda][wr_]:=FourierVel[wr,d\[Phi]spard\[Lambda]coeff];
	drgd\[Lambda][wr_]:=drgd\[Lambda]fun[wr,a,p,e,x];

	\[CapitalDelta]trg[wr_]:=\[CapitalDelta]trgIntVel[wr,\[CapitalUpsilon]rg,dtrgd\[Lambda]coeff];
	\[CapitalDelta]\[Phi]rg[wr_]:=\[CapitalDelta]\[Phi]rgIntVel[wr,\[CapitalUpsilon]rg,d\[Phi]rgd\[Lambda]coeff];
	rg[wr_]:=rgICr2gfun[wr,a,p,e,x];
	
	(* spin correction radial and polar frequencies *)
	\[CapitalUpsilon]rs=\[Delta]\[CapitalUpsilon]rfunFC[a,p,e,x];
	
	Print["Calculating Fourier coefficients of \!\(\*SubscriptBox[\(dt\), \(s\)]\)/d\[Lambda]"];
	dtspard\[Lambda]coeff=d\[Delta]td\[Lambda]coefffunFCPerPar[nmax,a,p,e,x];
	\[CapitalUpsilon]ts=dtspard\[Lambda]coeff[[nmax+1]];
	Print["Calculating Fourier coefficients of \!\(\*SubscriptBox[\(d\[Phi]\), \(s\)]\)/d\[Lambda]"];
	d\[Phi]spard\[Lambda]coeff=d\[Delta]\[Phi]d\[Lambda]coefffunFCPerPar[nmax,a,p,e,x];
	\[CapitalUpsilon]\[Phi]s=d\[Phi]spard\[Lambda]coeff[[nmax+1]];
	
	\[CapitalDelta]tspar[wr_]:=\[CapitalDelta]\[Delta]IntvelPerPar[wr,\[CapitalUpsilon]rg,\[CapitalUpsilon]rs,dtspard\[Lambda]coeff,dtrgd\[Lambda]coeff];
	\[CapitalDelta]\[Phi]spar[wr_]:=\[CapitalDelta]\[Delta]IntvelPerPar[wr,\[CapitalUpsilon]rg,\[CapitalUpsilon]rs,d\[Phi]spard\[Lambda]coeff,d\[Phi]rgd\[Lambda]coeff];
	\[Delta]rpar[wr_]:=\[Delta]rfunFCPerPar[wr,a,p,e,x];
	\[Psi]phase[wr_]:=\[Psi]pICr2g[wr,a,p,e,x];
	\[Delta]zort[wr_]:=\[Delta]zfunICr2gPer[wr,a,p,e,x];
	  
	dtspard\[Lambda][wr_]:=FourierVel[wr,dtspard\[Lambda]coeff];
	d\[Phi]spard\[Lambda][wr_]:=FourierVel[wr,d\[Phi]spard\[Lambda]coeff];
	drspard\[Lambda][wr_]:=\[Delta]vrfunFCPerPar[wr,a,p,e,x];
			
	<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 "MinoFrequenciesCorrection"->{\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]\[Phi]s},
	 "BLFrequenciesCorrection"->{\[CapitalUpsilon]rs/\[CapitalUpsilon]tg-\[CapitalUpsilon]rg/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts,\[CapitalUpsilon]\[Phi]s/\[CapitalUpsilon]tg-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts},
	 "MinoPrecessionFrequency"->\[CapitalUpsilon]p,
	 "BLPrecessionFrequency"->\[CapitalUpsilon]p/\[CapitalUpsilon]tg,
	 (*Keys purely oscillatory part geodesic coordinate time and azimuthal trajectory*)
	 "\[CapitalDelta]trg"->Function[{wr},\[CapitalDelta]trg[wr]],
	 "\[CapitalDelta]\[Phi]rg"->Function[{wr},\[CapitalDelta]\[Phi]rg[wr]],
	 (*Keys geodesic radial trajectory*)
	 "rg"->Function[{wr},rg[wr]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{wr},dtrgd\[Lambda][wr]],
	 "vrg"->Function[{wr},drgd\[Lambda][wr]],
	 "v\[Phi]g"->Function[{wr},d\[Phi]rgd\[Lambda][wr]],
	 (*Keys corrections trajectory*)
	 "\[CapitalDelta]\[Delta]tpar"->Function[{wr},\[CapitalDelta]tspar[wr]],
	 "\[Delta]rpar"->Function[{wr},\[Delta]rpar[wr]],
	 "\[Psi]p"->Function[wr,\[Psi]phase[wr]],
	 "\[Delta]zort"->Function[{wr},\[Delta]zort[wr]],
	 "\[CapitalDelta]\[Delta]\[Phi]par"->Function[{wr},\[CapitalDelta]\[Phi]spar[wr]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{wr},dtspard\[Lambda][wr]],
	 "\[Delta]vrpar"->Function[{wr},drspard\[Lambda][wr]],
	 "\[Delta]v\[Phi]par"->Function[{wr},d\[Phi]spard\[Lambda][wr]]
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
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	  
	(* spin correction radial and polar frequencies *)
	\[CapitalUpsilon]rs=\[Delta]\[CapitalUpsilon]rfunFE[a,p,e,x];
	\[CapitalUpsilon]ts=\[Delta]\[CapitalUpsilon]tfunFE[a,p,e,x];
	\[CapitalUpsilon]\[Phi]s=\[Delta]\[CapitalUpsilon]\[Phi]funFE[a,p,e,x];
				
<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->\[Delta]EE,
	 "Jzs"->\[Delta]Lz,
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 "MinoFrequenciesCorrection"->{\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]\[Phi]s},
	 "BLFrequenciesCorrection"->{\[CapitalUpsilon]rs/\[CapitalUpsilon]tg-\[CapitalUpsilon]rg/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts,\[CapitalUpsilon]\[Phi]s/\[CapitalUpsilon]tg-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts},
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
	 "\[CapitalDelta]\[Delta]tpar"->Function[{wr},\[Delta]tfunFEPerPar[wr,a,p,e,x]],
	 "\[Delta]tOfrpar"->Function[{r},\[Delta]tfunOfrFEPerPar[r,a,p,e,x]],
	 "\[Delta]rpar"->Function[{wr},\[Delta]rfunFEPerPar[wr,a,p,e,x]],
	 "\[Delta]rOfrpar"->Function[{r},\[Delta]rfunOfrFEPerPar[r,a,p,e,x]],
	 "\[Psi]p"->Function[{wr},\[Psi]pICr1g[wr,a,p,e,x]],
	 "\[Psi]pOfr"->Function[{r},\[Psi]pOfrICr1g[r,a,p,e,x]],
	 "\[Delta]zort"->Function[{wr},\[Delta]zfunICr1gPer[wr,a,p,e,x]],
	 "\[Delta]zOfrort"->Function[{r},\[Delta]zfunOfrICr1gPer[r,a,p,e,x]],
	 "\[CapitalDelta]\[Delta]\[Phi]par"->Function[{wr},\[Delta]\[Phi]funFEPerPar[wr,a,p,e,x]],
	 "\[Delta]\[Phi]Ofrpar"->Function[{r},\[Delta]\[Phi]funOfrFEPerPar[r,a,p,e,x]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{wr},\[Delta]vtfunFEPerPar[wr,a,p,e,x]],
	 "\[Delta]vtOfrpar"->Function[{r},\[Delta]vtfunOfrFEPerPar[r,a,p,e,x]],
	 "\[Delta]vrpar"->Function[{wr},\[Delta]vrfunFEPerPar[wr,a,p,e,x]],
	 "\[Delta]vrOfrpar"->Function[{r},\[Delta]vrfunOfrFEPerPar[r,a,p,e,x]],
	 "\[Delta]v\[Phi]par"->Function[{wr},\[Delta]v\[Phi]funFEPerPar[wr,a,p,e,x]],
	 "\[Delta]v\[Phi]Ofrpar"->Function[{r},\[Delta]v\[Phi]funOfrFEPerPar[r,a,p,e,x]]
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
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	  
	(* spin correction radial and polar frequencies *)
	\[CapitalUpsilon]ts=\[Delta]\[CapitalUpsilon]tfunLimFE[a,p,e,x];
	\[CapitalUpsilon]\[Phi]s=\[Delta]\[CapitalUpsilon]\[Phi]funLimFE[a,p,e,x];

<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->\[Delta]EE,
	 "Jzs"->\[Delta]Lz,
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 "MinoFrequenciesCorrection"->{\[CapitalUpsilon]ts,\[CapitalUpsilon]\[Phi]s},
	 "BLFrequenciesCorrection"->{\[CapitalUpsilon]\[Phi]s/\[CapitalUpsilon]tg-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts},
	 "MinoPrecessionFrequency"->\[CapitalUpsilon]p,
	 "BLPrecessionFrequency"->\[CapitalUpsilon]p/\[CapitalUpsilon]tg,
	 (*Keys purely oscillatory part geodesic coordinate time and azimuthal trajectory*)
	 "trg"->Function[{r},tgfunHom[r,a,p,e,x]],
	 "\[Phi]rg"->Function[{r},\[Phi]gfunHom[r,a,p,e,x]],
	 (*Keys geodesic radial trajectory*)
	 "rg"->Function[{\[Lambda]},rgfunHom[\[Lambda],a,p,e,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{r},VtrgfunHom[r,a,p,e,x]],
	 "vrg"->Function[{r},drgd\[Lambda]funHom[r,a,p,e,x]],
	 "v\[Phi]g"->Function[{r},V\[Phi]rgfunHom[r,a,p,e,x]],
	 (*Keys corrections trajectory*)
	 "\[Delta]tpar"->Function[{r},\[Delta]tfunFEHomPar[r,a,p,e,x]],
	 "\[Delta]rpar"->Function[{r},\[Delta]rfunFEHomPar[r,a,p,e,x]],
	 "\[Psi]p"->Function[{r},\[Psi]pHom[r,a,p,e,x]],
	 "\[Delta]zort"->Function[{r},\[Delta]zfunHom[r,a,p,e,x]],
	 "\[Delta]\[Phi]par"->Function[{r},\[Delta]\[Phi]funFEHomPar[r,a,p,e,x]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{r},\[Delta]vtfunFEHomPar[r,a,p,e,x]],
	 "\[Delta]vrpar"->Function[{r},\[Delta]vrfunFEHomPar[r,a,p,e,x]],
	 "\[Delta]v\[Phi]par"->Function[{r},\[Delta]v\[Phi]funFEHomPar[r,a,p,e,x]],
	 (*Shift separatrix*)
	 "\[Delta]p"->\[Delta]pfun[a,p,e,x]
	|>
]


(* ::Subsubsection::Closed:: *)
(*ISCO plunge solutions*)


KerrNearEqSpinOrbitCorrFEISCOPlunge[a_, x_]:=Module[{r1g,EEg,Lzg,\[Delta]EE,\[Delta]Lz},
	r1g=ISCOradius[a,x];
	(* geodesic constants of motion *)
	EEg=EEgfun[a,r1g,0,x];
	Lzg=Lzgfun[a,r1g,0,x];
	
	(* spin correction constants of motion *)
	\[Delta]EE=\[Delta]EEISCOfunFE[a,x];
	\[Delta]Lz=\[Delta]LzISCOfunFE[a,x];
	  
<|
	 "OrbitalElements"->{a,r1g,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->\[Delta]EE,
	 "Jzs"->\[Delta]Lz,
	 "BLFrequenciesGeo"->Function[{r},\[CapitalOmega]\[Phi]gfunISCOplunge[r,a,x]],
	 "BLFrequenciesCorrection"->Function[{r},\[Delta]\[CapitalOmega]\[Phi]funISCOplunge[r,a,x]],
	 "BLPrecessionFrequency"->Function[{r},\[CapitalOmega]pfunISCOplunge[r,a,x]],
	 "trg"->Function[{r},tgfunISCOplunge[r,a,x]],
	 "\[Phi]rg"->Function[{r},\[Phi]gfunISCOplunge[r,a,x]],
	 (*Keys geodesic radial trajectory*)
	 "rg"->Function[{\[Lambda]},rgfunISCOplunge[\[Lambda],a,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{r},VtrgfunISCOplunge[r,a,x]],
	 "vrg"->Function[{r},drgd\[Lambda]funISCOplunge[r,a,x]],
	 "v\[Phi]g"->Function[{r},V\[Phi]rgfunISCOplunge[r,a,x]],
	 (*Keys corrections trajectory*)
	 "\[Delta]tpar"->Function[{r},\[Delta]tfunISCOplungePar[r,a,x]],
	 "\[Delta]rpar"->Function[{r},\[Delta]rfunISCOplungePar[r,a,x]],
	 "\[Psi]p"->Function[{r},\[Psi]pISCOplunge[r,a,x]],
	 "\[Delta]zort"->Function[{r},\[Delta]zfunISCOplunge[r,a,x]],
	 "\[Delta]\[Phi]par"->Function[{r},\[Delta]\[Phi]funISCOplungePar[r,a,x]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{r},\[Delta]vtfunISCOplungePar[r,a,x]],
	 "\[Delta]vrpar"->Function[{r},\[Delta]vrfunISCOplungePar[r,a,x]],
	 "\[Delta]v\[Phi]par"->Function[{r},\[Delta]v\[Phi]funISCOplungePar[r,a,x]],
	 (*Shift ISCO*)
	 "\[Delta]ISCO"->\[Delta]rISCOfunFE[a,x]
	|>
]


(* ::Subsubsection::Closed:: *)
(*Critical plunge solutions*)


KerrNearEqSpinOrbitCorrFECritPlunge[a_, p_, e_, x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz},
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	
	(* spin correction constants of motion *)
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
			
<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->\[Delta]EE,
	 "Jzs"->\[Delta]Lz,
	 "BLFrequenciesGeo"->Function[{r},\[CapitalOmega]\[Phi]gfun[r,a,p,e,x]],
	 "BLFrequenciesCorrection"->Function[{r},\[Delta]\[CapitalOmega]\[Phi]funCritplunge[r,a,p,e,x]],
	 "BLPrecessionFrequency"->Function[{r},\[CapitalOmega]pfun[r,a,p,e,x]],
	 (*Keys purely oscillatory part geodesic coordinate time and azimuthal trajectory*)
	 "trg"->Function[{r},tgfunCritplunge[r,a,p,e,x]],
	 "\[Phi]rg"->Function[{r},\[Phi]gfunCritplunge[r,a,p,e,x]],
	 (*Keys geodesic radial trajectory*)
	 "rg"->Function[{\[Lambda]},rgfunCritplunge[\[Lambda],a,p,e,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{r},VtrgfunCritplunge[r,a,p,e,x]],
	 "vrg"->Function[{r},drgd\[Lambda]funCritplunge[r,a,p,e,x]],
	 "v\[Phi]g"->Function[{r},V\[Phi]rgfunCritplunge[r,a,p,e,x]],
	 (*Keys corrections trajectory*)
	 "\[Delta]tpar"->Function[{r},\[Delta]tfunCritplungePar[r,a,p,e,x]],
	 "\[Delta]rpar"->Function[{r},\[Delta]rfunCritplungePar[r,a,p,e,x]],
	 "\[Psi]p"->Function[{r},\[Psi]pCritplunge[r,a,p,e,x]],
	 "\[Delta]zort"->Function[{r},\[Delta]zfunCritplunge[r,a,p,e,x]],
	 "\[Delta]\[Phi]par"->Function[{r},\[Delta]\[Phi]funCritplungePar[r,a,p,e,x]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{r},\[Delta]vtfunCritplungePar[r,a,p,e,x]],
	 "\[Delta]vrpar"->Function[{r},\[Delta]vrfunCritplungePar[r,a,p,e,x]],
	 "\[Delta]v\[Phi]par"->Function[{r},\[Delta]v\[Phi]funCritplungePar[r,a,p,e,x]],
	 (*Shift separatrix*)
	 "\[Delta]p"->\[Delta]pfun[a,p,e,x]
	|>
]


(* ::Subsubsection::Closed:: *)
(*Generic plunge solutions*)


KerrNearEqSpinOrbitCorrFEPlunge[a_, p_, e_, x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz},
	(* geodesic constants of motion *)
	EEg=EEgfunCC[p,e];
	Lzg=LzgfunCC[a,p,e,x];
	
	(* spin correction constants of motion *)
	\[Delta]EE=\[Delta]EEfunCC[a,p,e,x]+dEEdpfunCC[p,e]\[Delta]pfunCC[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunCC[a,p,e,x]+dLzdpfunCC[a,p,e,x]\[Delta]pfunCC[a,p,e,x];
			
<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->\[Delta]EE,
	 "Jzs"->\[Delta]Lz,
	 "BLFrequenciesGeo"->Function[{r},\[CapitalOmega]\[Phi]gfunPlunge[r,a,p,e,x]],
	 "BLFrequenciesCorrection"->Function[{r},\[Delta]\[CapitalOmega]\[Phi]funPlunge[r,a,p,e,x]],
	 "BLPrecessionFrequency"->Function[{r},\[CapitalOmega]pfunPlunge[r,a,p,e,x]],
	 (*Keys purely oscillatory part geodesic coordinate time and azimuthal trajectory*)
	 "trg"->Function[{r},tgfunPlunge[r,a,p,e,x]],
	 "\[Phi]rg"->Function[{r},\[Phi]gfunPlunge[r,a,p,e,x]],
	 (*Keys geodesic radial trajectory*)
	 "rg"->Function[{\[Lambda]},rgfunPlunge[\[Lambda],a,p,e,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{r},VtrgfunPlunge[r,a,p,e,x]],
	 "vrg"->Function[{r},drgd\[Lambda]funPlunge[r,a,p,e,x]],
	 "v\[Phi]g"->Function[{r},V\[Phi]rgfunPlunge[r,a,p,e,x]],
	 (*Keys corrections trajectory*)
	 "\[Delta]tpar"->Function[{r},\[Delta]tfunPlungePar[r,a,p,e,x]],
	 "\[Delta]rpar"->Function[{r},\[Delta]rfunPlungePar[r,a,p,e,x]],
	 "\[Psi]p"->Function[r,\[Psi]pPlunge[r,a,p,e,x]],
	 "\[Delta]zort"->Function[{r},\[Delta]zfunPlunge[r,a,p,e,x]],
	 "\[Delta]\[Phi]par"->Function[{r},\[Delta]\[Phi]funPlungePar[r,a,p,e,x]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{r},\[Delta]vtfunPlungePar[r,a,p,e,x]],
	 "\[Delta]vrpar"->Function[{r},\[Delta]vrfunPlungePar[r,a,p,e,x]],
	 "\[Delta]v\[Phi]par"->Function[{r},\[Delta]v\[Phi]funPlungePar[r,a,p,e,x]],
	 (*Shift separatrix*)
	 "\[Delta]p"->\[Delta]pfunCC[a,p,e,x]
	|>
]


(* ::Subsubsection::Closed:: *)
(*Plunge related to bound orbits*)


KerrNearEqSpinOrbitCorrFEPlungeBoundOrbit[a_, p_, e_, x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz},
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	
	(* spin correction constants of motion *)
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
			
<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->\[Delta]EE,
	 "Jzs"->\[Delta]Lz,
	 "BLFrequenciesGeo"->Function[{r},\[CapitalOmega]\[Phi]gfun[r,a,p,e,x]],
	 "BLFrequenciesCorrection"->Function[{r},\[Delta]\[CapitalOmega]\[Phi]funPlungeBoundOrbit[r,a,p,e,x]],
	 "BLPrecessionFrequency"->Function[{r},\[CapitalOmega]pfun[r,a,p,e,x]],
	 (*Keys purely oscillatory part geodesic coordinate time and azimuthal trajectory*)
	 "trg"->Function[{r},tgfunPlungeBoundOrbit[r,a,p,e,x]],
	 "\[Phi]rg"->Function[{r},\[Phi]gfunPlungeBoundOrbit[r,a,p,e,x]],
	 (*Keys geodesic radial trajectory*)
	 "rg"->Function[{\[Lambda]},rgfunPlungeBoundOrbit[\[Lambda],a,p,e,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{r},VtrgfunPlungeBoundOrbit[r,a,p,e,x]],
	 "vrg"->Function[{r},drgd\[Lambda]funPlungeBoundOrbit[r,a,p,e,x]],
	 "v\[Phi]g"->Function[{r},V\[Phi]rgfunPlungeBoundOrbit[r,a,p,e,x]],
	 (*Keys corrections trajectory*)
	 "\[Delta]tpar"->Function[{r},\[Delta]tfunPlungeBoundOrbitPar[r,a,p,e,x]],
	 "\[Delta]rpar"->Function[{r},\[Delta]rfunPlungeBoundOrbitPar[r,a,p,e,x]],
	 "\[Psi]p"->Function[{r},\[Psi]pPlungeBoundOrbit[r,a,p,e,x]],
	 "\[Delta]zort"->Function[{r},\[Delta]zfunPlungeBoundOrbit[r,a,p,e,x]],
	 "\[Delta]\[Phi]par"->Function[{r},\[Delta]\[Phi]funPlungeBoundOrbitPar[r,a,p,e,x]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{r},\[Delta]vtfunPlungeBoundOrbitPar[r,a,p,e,x]],
	 "\[Delta]vrpar"->Function[{r},\[Delta]vrfunPlungeBoundOrbitPar[r,a,p,e,x]],
	 "\[Delta]v\[Phi]par"->Function[{r},\[Delta]v\[Phi]funPlungeBoundOrbitPar[r,a,p,e,x]],
	 (*Shift separatrix*)
	 "\[Delta]p"->\[Delta]pfun[a,p,e,x]
	|>
]


(* ::Subsubsection::Closed:: *)
(*Plunge related to bound orbits - double root*)


KerrNearEqSpinOrbitCorrFEPlungeBoundOrbitDoubleRoot[a_, p_, x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz},
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,0,x];
	Lzg=Lzgfun[a,p,0,x];
	
	(* spin correction constants of motion *)
	\[Delta]EE=\[Delta]EEfunFT[a,p,0,x]+dEEdpfun[a,p,0,x]\[Delta]pfun[a,p,0,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,0,x]+dLzdpfun[a,p,0,x]\[Delta]pfun[a,p,0,x];
			
<|
	 "OrbitalElements"->{a,p,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->\[Delta]EE,
	 "Jzs"->\[Delta]Lz,
	 "BLFrequenciesGeo"->Function[{r},\[CapitalOmega]\[Phi]gfun[r,a,p,0,x]],
	 "BLFrequenciesCorrection"->Function[{r},\[Delta]\[CapitalOmega]\[Phi]funPlungeBoundOrbitDoubleRoot[r,a,p,x]],
	 "BLPrecessionFrequency"->Function[{r},\[CapitalOmega]pfun[r,a,p,0,x]],
	 (*Keys purely oscillatory part geodesic coordinate time and azimuthal trajectory*)
	 "trg"->Function[{r},tgfunPlungeBoundOrbitDoubleRoot[r,a,p,x]],
	 "\[Phi]rg"->Function[{r},\[Phi]gfunPlungeBoundOrbitDoubleRoot[r,a,p,x]],
	 (*Keys geodesic radial trajectory*)
	 "rg"->Function[{\[Lambda]},rgfunPlungeBoundOrbitDoubleRoot[\[Lambda],a,p,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{r},VtrgfunPlungeBoundOrbit[r,a,0,p,x]],
	 "vrg"->Function[{r},drgd\[Lambda]funPlungeBoundOrbit[r,a,0,p,x]],
	 "v\[Phi]g"->Function[{r},V\[Phi]rgfunPlungeBoundOrbit[r,a,0,p,x]],
	 (*Keys corrections trajectory*)
	 "\[Delta]tpar"->Function[{r},\[Delta]tfunPlungeBoundOrbitDoubleRootPar[r,a,p,x]],
	 "\[Delta]rpar"->Function[{r},\[Delta]rfunPlungeBoundOrbitDoubleRootPar[r,a,p,x]],
	 "\[Psi]p"->Function[{r},\[Psi]pPlungeBoundOrbitDoubleRoot[r,a,p,x]],
	 "\[Delta]zort"->Function[{r},\[Delta]zfunPlungeBoundOrbitDoubleRoot[r,a,p,x]],
	 "\[Delta]\[Phi]par"->Function[{r},\[Delta]\[Phi]funPlungeBoundOrbitDoubleRootPar[r,a,p,x]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{r},\[Delta]vtfunPlungeBoundOrbitDoubleRootPar[r,a,p,x]],
	 "\[Delta]vrpar"->Function[{r},\[Delta]vrfunPlungeBoundOrbitDoubleRootPar[r,a,p,x]],
	 "\[Delta]v\[Phi]par"->Function[{r},\[Delta]v\[Phi]funPlungeBoundOrbitDoubleRootPar[r,a,p,x]],
	 (*Shift separatrix*)
	 "\[Delta]p"->\[Delta]pfun[a,p,0,x]
	|>
]


(* ::Subsubsection::Closed:: *)
(*Shifts constants of motion and IBCO radius for parabolic orbits*)


KerrNearEqSpinOrbitCorrFE\[Delta]IBCO[a_, x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,r,\[CapitalDelta],\[Delta]rIBCO,Lzred},
	(* geodesic constants of motion *)
	EEg=1;
	Lzg=RealSign[x]2(1+Sqrt[1-RealSign[x]a]);
	
	(* spin correction constants of motion *)
	\[Delta]EE=0;
	\[Delta]Lz=(1-RealSign[x]a+Sqrt[1-RealSign[x]a])/(2+2Sqrt[1-a]-RealSign[x]a);
	r=2-RealSign[x]a+2Sqrt[1-RealSign[x]a];
	\[CapitalDelta]=r^2-2r+a^2;
	Lzred=Lzg-a;
	
	\[Delta]rIBCO=-(1/2)((Lzg (r-2)^3 Sqrt[\[CapitalDelta]](2a*Lzred^4-2Lzred^3 r^2-2Lzred*a*Lzg*r^2+a*Lzg^2*r^3))/(Lzred^2(2a+Lzg(-2+r))r(-r Sqrt[\[CapitalDelta]](10a^2-4(4+\[CapitalDelta])+r(8+\[CapitalDelta]))+Sqrt[2]a Sqrt[r](2a^2(-1+r)+(4+r)\[CapitalDelta])RealSign[x])));
			
<|
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->\[Delta]EE,
	 "Jzs"->\[Delta]Lz,
	 "rIBCO"->r,
	 "\[Delta]rIBCO"->\[Delta]rIBCO
	|>
]


(* ::Subsection::Closed:: *)
(*Fourier series expansion*)


KerrNearEqSpinOrbitCorrFEPerFourier[a_, p_, e_, x_, nmax_?(IntegerQ[#] && # > 0&)]:=Module[{EEg,Lzg,\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g,\[CapitalUpsilon]p,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]zs,\[CapitalUpsilon]\[Phi]s,stepsr,\[CapitalDelta]tspar,\[CapitalDelta]\[Phi]spar,\[Delta]rpar,\[Psi]phase,\[Delta]zort,
	ExpniTable,wrlist,dtrgd\[Lambda]coeff,d\[Phi]rgd\[Lambda]coeff,dtspard\[Lambda]coeff,d\[Phi]spard\[Lambda]coeff,dtrgd\[Lambda],d\[Phi]rgd\[Lambda],drgd\[Lambda],rg,\[CapitalDelta]trg,\[CapitalDelta]\[Phi]rg,dtspard\[Lambda],drspard\[Lambda],dzspard\[Lambda],d\[Phi]spard\[Lambda]},
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
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	
	(* steps for numerical integration *)
	stepsr=4*nmax;
	(* matrices of discrete Fourier transform *)
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,x}]],{n,-nmax,nmax},{i,1,stepsr}];
	
	wrlist=Table[wr,{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
  
	(* Fourier coeffients geodesic functions *)
	dtrgd\[Lambda]coeff=VtrgcoeffICr1g[nmax,a,p,e,x];
	d\[Phi]rgd\[Lambda]coeff=V\[Phi]rgcoeffICr1g[nmax,a,p,e,x];
	
	dtrgd\[Lambda][wr_]:=FourierVel[wr,dtspard\[Lambda]coeff];
	d\[Phi]rgd\[Lambda][wr_]:=FourierVel[wr,d\[Phi]spard\[Lambda]coeff];
	drgd\[Lambda][wr_]:=drgd\[Lambda]fun[wr,a,p,e,x];

	\[CapitalDelta]trg[wr_]:=\[CapitalDelta]trgIntVel[wr,\[CapitalUpsilon]rg,dtrgd\[Lambda]coeff];
	\[CapitalDelta]\[Phi]rg[wr_]:=\[CapitalDelta]\[Phi]rgIntVel[wr,\[CapitalUpsilon]rg,d\[Phi]rgd\[Lambda]coeff];
	rg[wr_]:=rgICr1gfun[wr,a,p,e,x];
	
	(* spin correction radial and polar frequencies *)
	\[CapitalUpsilon]rs=\[Delta]\[CapitalUpsilon]rfunFE[a,p,e,x];
	
	Print["Calculating Fourier coefficients of \!\(\*SubscriptBox[\(dt\), \(s\)]\)/d\[Lambda]"];
	dtspard\[Lambda]coeff=d\[Delta]td\[Lambda]coefffunFEPerPar[nmax,a,p,e,x];
	\[CapitalUpsilon]ts=dtspard\[Lambda]coeff[[nmax+1]];
	Print["Calculating Fourier coefficients of \!\(\*SubscriptBox[\(d\[Phi]\), \(s\)]\)/d\[Lambda]"];
	d\[Phi]spard\[Lambda]coeff=d\[Delta]\[Phi]d\[Lambda]coefffunFEPerPar[nmax,a,p,e,x];
	\[CapitalUpsilon]\[Phi]s=d\[Phi]spard\[Lambda]coeff[[nmax+1]];
	
	\[CapitalDelta]tspar[wr_]:=\[CapitalDelta]\[Delta]IntvelPerPar[wr,\[CapitalUpsilon]rg,\[CapitalUpsilon]rs,dtspard\[Lambda]coeff,dtrgd\[Lambda]coeff];
	\[CapitalDelta]\[Phi]spar[wr_]:=\[CapitalDelta]\[Delta]IntvelPerPar[wr,\[CapitalUpsilon]rg,\[CapitalUpsilon]rs,d\[Phi]spard\[Lambda]coeff,d\[Phi]rgd\[Lambda]coeff];
	\[Delta]rpar[wr_]:=\[Delta]rfunFEPerPar[wr,a,p,e,x];
	\[Psi]phase[wr_]:=\[Psi]pICr1g[wr,a,p,e,x];
	\[Delta]zort[wr_]:=\[Delta]zfunICr1gPer[wr,a,p,e,x];
			  
	dtspard\[Lambda][wr_]:=FourierVel[wr,dtspard\[Lambda]coeff];
	d\[Phi]spard\[Lambda][wr_]:=FourierVel[wr,d\[Phi]spard\[Lambda]coeff];
	drspard\[Lambda][wr_]:=\[Delta]vrfunFEPerPar[wr,a,p,e,x];
			
	<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->\[Delta]EE,
	 "Jzs"->\[Delta]Lz,
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 "MinoFrequenciesCorrection"->{\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]\[Phi]s},
	 "BLFrequenciesCorrection"->{\[CapitalUpsilon]rs/\[CapitalUpsilon]tg-\[CapitalUpsilon]rg/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts,\[CapitalUpsilon]\[Phi]s/\[CapitalUpsilon]tg-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts},
	 "MinoPrecessionFrequency"->\[CapitalUpsilon]p,
	 "BLPrecessionFrequency"->\[CapitalUpsilon]p/\[CapitalUpsilon]tg,
	 (*Keys purely oscillatory part geodesic coordinate time and azimuthal trajectory*)
	 "\[CapitalDelta]trg"->Function[{wr},\[CapitalDelta]trg[wr]],
	 "\[CapitalDelta]\[Phi]rg"->Function[{wr},\[CapitalDelta]\[Phi]rg[wr]],
	 (*Keys geodesic radial trajectory*)
	 "rg"->Function[{wr},rg[wr]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{wr},dtrgd\[Lambda][wr]],
	 "vrg"->Function[{wr},drgd\[Lambda][wr]],
	 "v\[Phi]g"->Function[{wr},d\[Phi]rgd\[Lambda][wr]],
	(*Keys corrections trajectory*)
	 "\[CapitalDelta]\[Delta]tpar"->Function[{wr},\[CapitalDelta]tspar[wr]],
	 "\[Delta]rpar"->Function[{wr},\[Delta]rpar[wr]],
	 "\[Psi]p"->Function[wr,\[Psi]phase[wr]],
	 "\[Delta]zort"->Function[{wr},\[Delta]zort[wr]],
	 "\[CapitalDelta]\[Delta]\[Phi]par"->Function[{wr},\[CapitalDelta]\[Phi]spar[wr]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{wr},dtspard\[Lambda][wr]],
	 "\[Delta]vrpar"->Function[{wr},drspard\[Lambda][wr]],
	 "\[Delta]v\[Phi]par"->Function[{wr},d\[Phi]spard\[Lambda][wr]]
	|>
]


(* ::Section::Closed:: *)
(*End package*)


End[];


EndPackage[];
