(* ::Package:: *)

(* ::Section:: *)
(*Begin package*)


BeginPackage["AnalyticSpinningOrbitNearEq`"];


(* ::Text:: *)
(*If you make use of this package, please acknowledge "Piovano" arXiv:2510.09597 (https://arxiv.org/abs/2510.09597 )*)


(* ::Text:: *)
(*IMPORTANT NOTE: the following package include all the functions to compute the analytic spin-corrections to the orbits (trajectories and velocities), constants of motion and frequencies for the fixed turning points (or "FT"), fixed constants of motion (or "FC"), and fixed eccentricity (or "FE") spin-gauges. It includes the contributions to both the parallel and orthogonal components of the secondary spin. *)
(*Homoclinic orbits are For the FE spin-gauge, the*)


(* ::Text:: *)
(*IMPORTANT NOTE 2: the functions for the FT and FC spin-gauges consider as initial radius at \[Lambda]=0 the geodesic periastron, while all functions for the FE spin-gauge the initial *)


KerrNearEqSpinOrbitCorrFTPer::usage = "KerrNearEqSpinOrbitCorrFTPer[a, p, e, x] calculates the linear corrections to periodic orbits in the fixed turning points spin-gauge. The initial radius at \[Lambda]=0 is the geodesic periastron.";


KerrNearEqSpinOrbitCorrFTPerFourier::usage = "KerrNearEqSpinOrbitCorrFTPerFourier[a, p, e, x, nmax] calculates the linear corrections to periodic orbits in the fixed turning points spin-gauge. The coordinate time and azimuthal trajectories are expanded in Fourier series. The initial radius at \[Lambda]=0 is the geodesic periastron.";


KerrNearEqSpinOrbitCorrFCPer::usage = "KerrNearEqSpinOrbitCorrFCPer[a, p, e, x] calculates the linear corrections to periodic orbits in the fixed constants of motion spin-gauge. The initial radius at \[Lambda]=0 is the geodesic periastron.";


KerrNearEqSpinOrbitCorrFCPerFourier::usage = "KerrNearEqSpinOrbitCorrFCPerFourier[a, p, e, x, nmax]  calculates the linear corrections to periodic orbits in the fixed constants of motion spin-gauge. The coordinate time and azimuthal trajectories are expanded in Fourier series. The initial radius at \[Lambda]=0 is the geodesic periastron.";


KerrNearEqSpinOrbitCorrFEPer::usage = "KerrNearEqSpinOrbitCorrFEPer[a, p, e, x]  calculates the linear corrections to periodic orbits in the fixed eccentricity spin-gauge. The initial radius at \[Lambda]=0 is the geodesic apoastron.";


KerrNearEqSpinOrbitCorrFEPerFourier::usage = "KerrNearEqSpinOrbitCorrFEPerFourier[a, p, e, x, nmax]  calculates the linear corrections to periodic orbits in the fixed eccentricity spin-gauge. The coordinate time and azimuthal trajectories are expanded in Fourier series. The initial radius at \[Lambda]=0 is the geodesic periastron.";


KerrNearEqSpinOrbitCorrFEHom::usage = "KerrNearEqSpinOrbitCorrFEHom[a, p, e, x] calculates the linear corrections to homoclinic orbits.";


Begin["`Private`"];


(*ADD ME: add algorithms to automatic select radial modes based on the sought precision for Fourier expansions.*)


(*IMPORTANT: still to check for memory leakages*)


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


(* ::Section::Closed:: *)
(*Geodesic quantities*)


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
(*Azimuthal time geodesic frequency*)


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


\[CapitalUpsilon]pfunLim[a_,p_,e_,xg_]:=Module[{r1g,r2g,EEg,Lzg,krg,ellK,Lzred,\[Gamma]r,\[Psi]freq},
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	r1g=p/(1-e);
	r2g=p/(1+e);

	Lzred=Lzg-a*EEg;
	
 Lzred(a+EEg*Lzred)r2g^2/(Abs[Lzred](r2g^2+Lzred^2))
]


(* ::Subsection::Closed:: *)
(*Geodesics periodic trajectory*)


(* ::Subsubsection::Closed:: *)
(*Geodesic radial trajectory*)


rgICr2gfun[wr_,a_,p_,e_,xg_]:=Module[{r1g,r2g,EEg,Lzg,r3g,krg,ellK,jSN},
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	ellK=EllipticK[krg];
	jSN=JacobiSN[1/\[Pi]*ellK*wr,krg];

	(r3g(r1g-r2g)jSN^2-r2g(r1g-r3g))/((r1g-r2g)jSN^2-(r1g-r3g))
]


rgICr1gfun[wr_,a_,p_,e_,x_]:=Module[{r1g,r2g,EEg,Lzg,KKg,r3g,krg,ellK,jSN},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];

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


tgfunref[a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1g,r2g,r3g,krg,\[Gamma]r,ellK,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r,\[ScriptCapitalI]r2,\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
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
	\[ScriptCapitalI]r=2/Sqrt[r2g(r1g-r3g)]r1g*ellPi;
	\[ScriptCapitalI]r2=-Sqrt[r2g(r1g-r3g)](-ellE+r1g*ellK/(r1g-r3g)-r1g(r1g+r2g+r3g)ellPi/(r2g(r1g-r3g)));
	\[ScriptCapitalI]p=2/Sqrt[(r1g-r3g)r2g] 1/rp (r1g/(r1g-rp) 2EllipticPi[(rp(r1g-r2g))/((r1g-rp)r2g),krg]- ellK);
	\[ScriptCapitalI]mreg=2/Sqrt[(r1g-r3g)r2g] (r1g/(r1g-rm) 2EllipticPi[(rm(r1g-r2g))/((r1g-rm)r2g),krg]- ellK);

	1/(2\[Pi]) 1/Sqrt[1-EEg^2] (4EEg*\[ScriptCapitalI]+2EEg*\[ScriptCapitalI]r+EEg*\[ScriptCapitalI]r2-4EEg rp/(rp-rm) (rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)-2(-4EEg+a Lzg) 1/(rp-rm) (rp \[ScriptCapitalI]p-\[ScriptCapitalI]mreg))
]


(* ::Text:: *)
(*Purely oscillatory part of the geodesic coordinate time trajectory*)


tgICr2gfun[wr_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1g,r2g,r3g,krg,krghold,\[Gamma]r,\[Phi],ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r,\[ScriptCapitalI]r2,\[ScriptCapitalI]p,\[ScriptCapitalI]m},
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
	\[ScriptCapitalI]r=2/Sqrt[r2g(r1g-r3g)]((r2g-r3g)ellPi+r3g ellF);
	\[ScriptCapitalI]r2=2r2g^2/Sqrt[r2g(r1g-r3g)]((r3g/r2g)^2*ellF+2 r3g/r2g (1-r3g/r2g)ellPi+(1-r3g/r2g)^2 1/(2(\[Gamma]r-1)(krg-\[Gamma]r)) (\[Gamma]r*ellE+(krg-\[Gamma]r)ellF+(2\[Gamma]r*krg+2\[Gamma]r-\[Gamma]r^2-3krg)ellPi-(\[Gamma]r^2Sin[\[Phi]]Cos[\[Phi]]Sqrt[1-krg Sin[\[Phi]]^2])/(1-\[Gamma]r*Sin[\[Phi]]^2)));

	\[ScriptCapitalI]p=2/(Sqrt[r2g(r1g-r3g)](r3g-rp))(ellF-(r2g-r3g)/(r2g-rp)EllipticPi[(r3g-rp)(r1g-r2g)/((r1g-r3g)(r2g-rp)),\[Phi],krg]);
	\[ScriptCapitalI]m=2/(Sqrt[r2g(r1g-r3g)](r3g-rm))(ellF-(r2g-r3g)/(r2g-rm)EllipticPi[(r3g-rm)(r1g-r2g)/((r1g-r3g)(r2g-rm)),\[Phi],krg]);


	1/Sqrt[1-EEg^2](4EEg*\[ScriptCapitalI]+2EEg*\[ScriptCapitalI]r+EEg*\[ScriptCapitalI]r2-4a^2*EEg 1/(rp-rm)(\[ScriptCapitalI]p-\[ScriptCapitalI]m)-2(-4EEg+a*Lzg)1/(rp-rm)(rp*\[ScriptCapitalI]p-rm*\[ScriptCapitalI]m))-wr*tgfunref[a,p,e,x]
]


tgICr1gfun[wr_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1g,r2g,r3g,krg,krghold,\[Gamma]r,\[Phi],ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r,\[ScriptCapitalI]r2,\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
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
	\[ScriptCapitalI]r=2/Sqrt[r2g(r1g-r3g)]r1g*ellPi;
	\[ScriptCapitalI]r2=-Sqrt[r2g(r1g-r3g)](-ellE+r1g*ellF/(r1g-r3g)-r1g(r1g+r2g+r3g)ellPi/(r2g(r1g-r3g))+((r1g-r2g)Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(-r2g+(-r1g+r2g)Sin[\[Phi]]^2));
	\[ScriptCapitalI]p=2/Sqrt[(r1g-r3g)r2g] 1/rp (r1g/(r1g-rp) EllipticPi[(rp(r1g-r2g))/((r1g-rp)r2g),\[Phi],krg]- ellF);
	\[ScriptCapitalI]mreg=2/Sqrt[(r1g-r3g)r2g] (r1g/(r1g-rm) EllipticPi[(rm(r1g-r2g))/((r1g-rm)r2g),\[Phi],krg]- ellF);

	1/Sqrt[1-EEg^2](4EEg*\[ScriptCapitalI]+2EEg*\[ScriptCapitalI]r+EEg*\[ScriptCapitalI]r2-4rp*EEg 1/(rp-rm)(rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)-2(-4EEg+a*Lzg)1/(rp-rm)(rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg))-wr*tgfunref[a,p,e,x]
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
(*Geodesics homoclinic trajectory*)


(* ::Subsubsection::Closed:: *)
(*Geodesic radial trajectory*)


rgfunHom[\[Lambda]_,a_,p_,e_,xg_]:=Module[{r1g,r2g,EEg,aux},
	EEg=EEgfun[a,p,e,xg];

	r1g=p/(1-e);
	r2g=p/(1+e);
	aux=1/2 Sqrt[(1-EEg^2)r2g(r1g-r2g)];
	r1g*r2g/(r2g+(r1g-r2g)Tanh[aux*\[Lambda]]^2)
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory*)


tgfunHom[r_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,r1g,r2g,\[ScriptCapitalI],\[ScriptCapitalI]r,\[ScriptCapitalI]r2,\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	If[e==0,
		0,
		rp=1+Sqrt[1-a^2];
		rm=1-Sqrt[1-a^2];
		EEg=EEgfun[a,p,e,x];
		Lzg=Lzgfun[a,p,e,x];
		r1g=p/(1-e);
		r2g=p/(1+e);
		\[ScriptCapitalI]=-1/(2Sqrt[r2g(r1g-r2g)])Log[(Sqrt[r(r1g-r2g)]-Sqrt[(-r+r1g)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(-r+r1g)r2g])^2];
		\[ScriptCapitalI]r=2 ArcCos[Sqrt[r/r1g]]+r2g*\[ScriptCapitalI];
		\[ScriptCapitalI]r2=Sqrt[r(-r+r1g)]+r1g ArcCos[Sqrt[r/r1g]]+r2g*\[ScriptCapitalI]r;
		\[ScriptCapitalI]p=\[ScriptCapitalI]/(r2g-rp)+1/(2Sqrt[r1g-rp](r2g-rp)Sqrt[rp])Log[(Sqrt[r(r1g-rp)]-Sqrt[(-r+r1g)rp])^2/(Sqrt[r(r1g-rp)]+Sqrt[(-r+r1g)rp])^2];
		\[ScriptCapitalI]mreg=\[ScriptCapitalI]/(r2g-rm)+Sqrt[rm]/(2Sqrt[r1g-rm](r2g-rm))Log[(Sqrt[r (r1g-rm)]-Sqrt[(-r+r1g) rm])^2/(Sqrt[r(r1g-rm)]+Sqrt[(-r+r1g) rm])^2];

		1/Sqrt[1-EEg^2] (4EEg*\[ScriptCapitalI]+2EEg*\[ScriptCapitalI]r+EEg*\[ScriptCapitalI]r2-4EEg rp/(rp-rm) (rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)-2(-4EEg+a*Lzg) 1/(rp-rm) (rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg))
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
		\[ScriptCapitalI]mreg=\[ScriptCapitalI]/(r2g-rm)+Sqrt[rm]/(2Sqrt[r1g-rm](r2g-rm))Log[(Sqrt[r(r1g-rm)]-Sqrt[(-r+r1g)rm])^2/(Sqrt[r(r1g-rm)]+Sqrt[(-r+r1g)rm])^2];

		1/Sqrt[1-EEg^2] (Lzg*\[ScriptCapitalI]-Lzg rp/(rp-rm) (rm*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg)+2a*EEg 1/(rp-rm)(rp*\[ScriptCapitalI]p-\[ScriptCapitalI]mreg))
	]
]


(* ::Subsection::Closed:: *)
(*Geodesics homoclinic velocity*)


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


(* ::Section::Closed:: *)
(*Shifts to the constants of motion*)


(* ::Subsection::Closed:: *)
(*Periodic orbits - fixed turning points*)


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


\[Delta]r1funFE[a_,p_,e_,x_]:=\[Delta]pfun[a,p,e,x]/(1-e)


\[Delta]r2funFE[a_,p_,e_,x_]:=\[Delta]pfun[a,p,e,x]/(1+e)


(* ::Subsection::Closed:: *)
(*Shift semilatus rectum*)


\[Delta]pfun[a_,p_,e_,x_]:=Module[{EEg,dEdp,\[Delta]EE,\[Delta]\[Rho]r4},
	EEg=EEgfun[a,p,e,x];
	dEdp= dEEdpfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	\[Delta]\[Rho]r4=\[Delta]r41fun[a,p,e,x];

	(2(1-e^2)(2*EEg*\[Delta]EE-(1-EEg^2)^2\[Delta]\[Rho]r4))/(-4dEdp*EEg+4*dEdp*e^2*EEg+3(1-EEg^2)^2-e (1-EEg^2)^2)
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


(* ::Section::Closed:: *)
(*Spin corrections to the orbit  - parallel component of the spin*)


(* ::Subsubsection::Closed:: *)
(*Building Fourier series for spin corrections coordinate time and azimuthal orbits*)


(* ::Text:: *)
(*Integration Fourier expansion  of the velocities*)


\[CapitalDelta]\[Delta]IntvelPerPar[wr_,\[CapitalUpsilon]rg_,\[Delta]\[CapitalUpsilon]r_,\[Delta]coeff_,coeffVgr_]:=Module[{dimr},
	dimr=(Length[\[Delta]coeff]-1)/2;
	
	Re[Sum[2Sin[n*wr]\[Delta]coeff[[n+dimr+1]]/(n*\[CapitalUpsilon]rg),{n,1,dimr}]-\[Delta]\[CapitalUpsilon]r*Sum[2Sin[n*wr]coeffVgr[[n+dimr+1]]/(n*\[CapitalUpsilon]rg^2),{n,1,dimr}]]
]


(* ::Subsection::Closed:: *)
(*Periodic trajectory - fixed turning points*)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunFTPerPar[wr_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]KK,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
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

	-((jSN*jCN(r1g-r2g)(r1g-r3g)(r2g-r3g)Sqrt[r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g])/((1-jSN^2) r1g+jSN^2 r2g-r3g)^2)(\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg (2ellF)/Yint+(EEg \[Delta]EE)/(1-EEg^2)(2ellF)/Yint+\[Delta]r3/2 1/r3g  2/Yint((r2g ellE)/(r2g-r3g)-ellF)+2((jSN*jCN(r1g-r2g))/(r1g*r2g Sqrt[r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g])-((r1g-r3g)ellE-r1g*ellF)/(r1g*r3g*Yint))\[Delta]r41-a/2 (-(2/(3r1g*r2g*r3g))((jSN*jCN(r1g-r2g)(r1g-r3g)(r2g-r3g)Sqrt[(r1g*r2g-jSN^2*r1g*r3g-r2g*r3g*jCN^2)])/(r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g)^2+((r1g+r2g+r3g)ellF)/Sqrt[r2g(r1g-r3g)]-2(r2g*r3g+r1g(r2g+r3g))((jSN*jCN(r1g-r2g))/(r1g*r2g Sqrt[r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g])+(-((r1g-r3g)ellE)+r1g*ellF)/(r1g*r3g Sqrt[r2g(r1g-r3g)] )))))
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunFTPerParRes[wr_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]KK,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
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


(* ::Subsubsection::Closed:: *)
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
(*Periodic trajectory - fixed constants of motion*)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunFCPerPar[wr_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,EEg,Lzg,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE,rg},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
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
	-((jSN*jCN(r1g-r2g)(r1g-r3g)(r2g-r3g)Sqrt[r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g])/((1-jSN^2)r1g+jSN^2 r2g-r3g)^2)(\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg(2ellF)/Yint+(r2g*ellE-r1g*ellF)/(r1g(r1g-r2g)Yint)\[Delta]r1+(-(r1g-r3g)ellE/((r1g-r2g)(r2g-r3g))+ellF/(r1g-r2g))/Yint*\[Delta]r2+\[Delta]r3/2*1/r3g  2/Yint((r2g ellE)/(r2g-r3g)-ellF)+2((jSN*jCN(r1g-r2g))/(r1g*r2g Sqrt[r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g])-((r1g-r3g)ellE-r1g*ellF)/(r1g*r3g*Yint))\[Delta]r41-a/2 (-(2/(3r1g*r2g*r3g))((jSN*jCN(r1g-r2g)(r1g-r3g)(r2g-r3g)Sqrt[(r1g*r2g-jSN^2*r1g*r3g-r2g*r3g*jCN^2)])/(r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g)^2+((r1g+r2g+r3g)ellF)/Sqrt[r2g(r1g-r3g)]-2(r2g*r3g+r1g(r2g+r3g))((jSN*jCN(r1g-r2g))/(r1g*r2g Sqrt[r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g])+(-((r1g-r3g)ellE)+r1g*ellF)/(r1g*r3g Sqrt[r2g(r1g-r3g)])))))+(rg(r2g(-r2g+rg)\[Delta]r1+r1g(r1g-rg)\[Delta]r2))/(r1g(r1g-r2g)r2g)
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunFCPerParRes[wr_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,EEg,Lzg,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE,rg},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
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


(* ::Subsubsection::Closed:: *)
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
(*Periodic trajectory - fixed eccentricity*)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunFEPerPar[wr_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,\[Delta]r1,\[Delta]r2,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]KK,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE,rg},
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
	\[Gamma]rg=(r1g-r2g)/(r1g-r3g);
	Yint=Sqrt[(r1g-r3g)r2g];

	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result.*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	jSN=Sin[\[Phi]];
	jCN=Cos[\[Phi]];

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	
	rg=rgICr1gfun[wr,a,p,e,x];

	((jSN*jCN*r1g(r1g-r2g)r2g*Sqrt[r1g(r2g-jSN^2*r3g)-(1-jSN^2)r2g*r3g])/(jSN^2(r1g-r2g)+r2g)^2)(\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg*(2ellF)/Sqrt[(r1g-r3g)r2g]+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Sqrt[(r1g-r3g)r2g]+((r2g*ellE-r1g*ellF)\[Delta]r1)/(Sqrt[(r1g-r3g)r2g]r1g(r1g-r2g))+1/Sqrt[(r1g-r3g)r2g](-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g))\[Delta]r2+1/Sqrt[(r1g-r3g)r2g] 1/r3g(((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]))\[Delta]r3+2((-(r1g-r3g)ellE+r1g*ellF)/(r1g Sqrt[r2g(r1g-r3g)]r3g))\[Delta]r41-a/2 (-((4ellE Sqrt[r2g(r1g-r3g)](r2g*r3g+r1g(r2g+r3g)))/(3r1g^2r2g^2r3g^2))+(2ellF((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g Sqrt[r2g (r1g-r3g)]r3g^2)-(2(-r1g+r2g)(r1g-r3g)^2Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3r1g^2(r2g(r1g-r3g))^(3/2)r3g)))+(((rg-r2g)(rg-r3g))/((r1g-r2g)(r1g-r3g))\[Delta]r1+((r1g-rg)(rg-r3g))/((r1g-r2g)(r2g-r3g))\[Delta]r2)
]


\[Delta]rfunOfrFEPerPar[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]KK,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],krghold,krg,ellF,ellE},
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
	\[Phi]=ArcSin[Sqrt[(r2g(r1g-r))/((r1g-r2g)r)]];
	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	Sqrt[(r1g-r)(r-r2g)r(r-r3g)](\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg*(2ellF)/Sqrt[(r1g-r3g)r2g]+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Sqrt[(r1g-r3g)r2g]+(r2g*ellE-r1g*ellF)/(Sqrt[(r1g-r3g)r2g]r1g(r1g-r2g) ) \[Delta]r1+1/Sqrt[(r1g-r3g)r2g] (-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g))\[Delta]r2+1/Sqrt[(r1g-r3g)r2g] 1/r3g (((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]))\[Delta]r3+ 2((-(r1g-r3g)ellE+r1g*ellF)/(r1g*r3g Sqrt[r2g(r1g-r3g)] ))\[Delta]r41-a/2 (-((4ellE*Sqrt[r2g(r1g-r3g)](r2g*r3g+r1g(r2g+r3g)))/(3r1g^2r2g^2r3g^2))+(2ellF((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g*r3g^2 Sqrt[r2g (r1g-r3g)] )-(2(-r1g+r2g)(r1g-r3g)Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3r1g^2r2g*r3g Sqrt[r2g(r1g-r3g)])))+(((r-r2g)(r-r3g))/((r1g-r2g)(r1g-r3g) ) \[Delta]r1+((r1g-r)(r-r3g))/((r1g-r2g)(r2g-r3g) ) \[Delta]r2)
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunFEPerParRes[wr_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,\[Delta]r1,\[Delta]r2,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]KK,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE,rg},
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
	\[Gamma]rg=(r1g-r2g)/(r1g-r3g);
	Yint=Sqrt[(r1g-r3g)r2g];
	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result. Having \[Phi]=\[Pi]/2 with MachinePrecision cause some problems with the integrals \[ScriptCapitalI]p and \[ScriptCapitalI]m*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	jSN=Sin[\[Phi]];
	jCN=Cos[\[Phi]];

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	
	rg=rgICr1gfun[wr,a,p,e,x];
	\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg(2ellF)/Sqrt[(r1g-r3g)r2g]+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Sqrt[(r1g-r3g)r2g]+((r2g*ellE-r1g*ellF)\[Delta]r1)/(Sqrt[(r1g-r3g)r2g]r1g(r1g-r2g) )+1/Sqrt[(r1g-r3g)r2g](-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g))\[Delta]r2+1/Sqrt[(r1g-r3g)r2g] 1/r3g (((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]))\[Delta]r3+2((-(r1g-r3g)ellE+r1g*ellF)/(r1g*r3g*Sqrt[r2g(r1g-r3g)]))\[Delta]r41-a/2(-((4ellE Sqrt[r2g(r1g-r3g)](r2g*r3g+r1g(r2g+r3g)))/(3r1g^2r2g^2r3g^2))+(2ellF((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g*r3g^2 Sqrt[r2g(r1g-r3g)])-(2(-r1g+r2g)(r1g-r3g)Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3r1g^2r2g*r3g Sqrt[r2g(r1g-r3g)]))+Sign[\[Pi]-Mod[wr,2\[Pi]]](Sqrt[(-r2g+rg)(-r3g+rg)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-rg)rg])\[Delta]r1+Sqrt[(r1g-rg)(-r3g+rg)]/((r1g-r2g)(r2g-r3g)Sqrt[rg(-r2g+rg)])\[Delta]r2)
]


\[Delta]rfunOfrFEPerParRes[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[Delta]KK,\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r,\[Phi],krghold,krg,\[Gamma]rg,ellF,ellE},
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
	\[Phi]=ArcSin[Sqrt[(r2g (r1g-r))/((r1g-r2g)r)]];
	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	
	(\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg*(2ellF)/Sqrt[(r1g-r3g)r2g]+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Sqrt[(r1g-r3g)r2g]+(r2g*ellE-r1g*ellF)/(Sqrt[(r1g-r3g)r2g]r1g(r1g-r2g))\[Delta]r1+1/Sqrt[(r1g-r3g)r2g](-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g))\[Delta]r2+1/Sqrt[(r1g-r3g)r2g] 1/r3g (((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]))\[Delta]r3+ 2((-(r1g-r3g)ellE+r1g*ellF)/(r1g*r3g*Sqrt[r2g(r1g-r3g)] ))\[Delta]r41-a/2(-((4ellE*Sqrt[r2g(r1g-r3g)](r2g*r3g+r1g(r2g+r3g)))/(3r1g^2r2g^2r3g^2))+(2ellF((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g*r3g^2*Sqrt[r2g (r1g-r3g)])-(2(-r1g+r2g)(r1g-r3g)Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3r1g^2r2g*r3g Sqrt[r2g(r1g-r3g)])))+(Sqrt[(-r2g+r)(-r3g+r)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-r)r])\[Delta]r1+Sqrt[(r1g-r)(-r3g+r)]/((r1g-r2g)(r2g-r3g)Sqrt[r(-r2g+r)])\[Delta]r2)
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory - purely oscillatory part*)


\[Delta]tfunFEPerPar[wr_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[CapitalUpsilon]rg,\[CapitalUpsilon]tg,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,\[Phi],rg,dtgd\[Lambda]fun,krghold,krg,\[Gamma]r,ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
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

	krg=((r1g-r2g) r3g)/(r2g (r1g-r3g));
	\[Gamma]r=1-r1g/r2g;

	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result.*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;

	rg=rgICr1gfun[wr,a,p,e,x];

	dtgd\[Lambda]fun=(rg(-2a*Lzg+EEg(2a^2+a^2*rg+rg^3)))/(a^2-2rg+rg^2);

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	ellPi=EllipticPi[\[Gamma]r,\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]);
	\[ScriptCapitalI]r1g=2 /(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]r1g) ((r2g *ellE-r1g*ellF)/(r1g-r2g) );
	\[ScriptCapitalI]r2g=2 /(Sqrt[1-EEg^2] Sqrt[r2g (r1g-r3g)]) (-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g));
	\[ScriptCapitalI]r3g=2 /(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]r3g) ((r2g*ellE)/(r2g-r3g)-ellF+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]));
	\[ScriptCapitalI]rover=(2r1g)/(Sqrt[1-EEg^2] Sqrt[(r1g-r3g)r2g]) ellPi;
	\[ScriptCapitalI]r2over=(r2g ((r1g-r3g)ellE-r1g*ellF))/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)])+(r1g(r1g+r2g+r3g)ellPi)/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)])+((r1g-r2g)Sqrt[r2g (r1g-r3g)]Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(Sqrt[1-EEg^2](r2g+(r1g-r2g)Sin[\[Phi]]^2));

	(4\[Delta]EE)/(1-EEg^2) \[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+1/2 EEg(\[Delta]r1+\[Delta]r2+\[Delta]r3+2\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+\[Delta]EE/(1-EEg^2) \[ScriptCapitalI]r2over+1/2 EEg((2+r1g)\[Delta]r1+(2+r2g)\[Delta]r2+(2+r3g)\[Delta]r3-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+(r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(2(r1g-rm)(r1g-rp)) \[ScriptCapitalI]r1g*\[Delta]r1+(r2g(-2a*Lzg+EEg(2a^2+a^2*r2g+r2g^3)))/(2(r2g-rm)(r2g-rp)) \[ScriptCapitalI]r2g*\[Delta]r2+(r3g(-2a*Lzg+EEg(2a^2+a^2*r3g+r3g^3)))/(2(r3g-rm)(r3g-rp)) \[ScriptCapitalI]r3g*\[Delta]r3-1/\[CapitalUpsilon]rg (\[Delta]\[CapitalUpsilon]tfunFE[a,p,e,x]-\[CapitalUpsilon]tg/\[CapitalUpsilon]rg \[Delta]\[CapitalUpsilon]rfunFE[a,p,e,x])wr+Sign[\[Pi]-Mod[wr,2 \[Pi]]] (r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(Sqrt[1-EEg^2](r1g-rm)(r1g-rp)) (Sqrt[(-r2g+rg)(-r3g+rg)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-rg)rg]))\[Delta]r1+Sign[\[Pi]-Mod[wr,2 \[Pi]]] (r2g(-2a*Lzg+EEg(2a^2+a^2*r2g+r2g^3)))/(Sqrt[1-EEg^2](r2g-rm)(r2g-rp)) ( Sqrt[(r1g-rg)(-r3g+rg)]/((r1g-r2g)(r2g-r3g)Sqrt[rg(-r2g+rg)]))\[Delta]r2-dtgd\[Lambda]fun/Sqrt[1-EEg^2] \[Delta]rfunFEPerParRes[wr,a,p,e,x]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory as function of r*)


\[Delta]tfunOfrFEPerPar[r_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[CapitalUpsilon]rg,\[CapitalUpsilon]tg,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Phi],\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,dtgd\[Lambda]fun,krghold,krg,\[Gamma]r,ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
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

	krg=((r1g-r2g) r3g)/(r2g (r1g-r3g));
	\[Gamma]r=1-r1g/r2g;

	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result.*)
	\[Phi]=ArcSin[Sqrt[(r2g (r1g-r))/((r1g-r2g)r)]];

	dtgd\[Lambda]fun=(r(-2a*Lzg+EEg(2a^2+a^2*r+r^3)))/(a^2-2r+r^2);

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	ellPi=EllipticPi[\[Gamma]r,\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]);
	\[ScriptCapitalI]r1g=2 /(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]r1g) ((r2g *ellE-r1g*ellF)/(r1g-r2g) );
	\[ScriptCapitalI]r2g=2 /(Sqrt[1-EEg^2] Sqrt[r2g (r1g-r3g)]) (-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g));
	\[ScriptCapitalI]r3g=2 /(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]r3g) ((r2g*ellE)/(r2g-r3g)-ellF+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]));
	\[ScriptCapitalI]rover=(2r1g)/(Sqrt[1-EEg^2] Sqrt[(r1g-r3g)r2g]) ellPi;
	\[ScriptCapitalI]r2over=(r2g ((r1g-r3g)ellE-r1g*ellF))/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)])+(r1g(r1g+r2g+r3g)ellPi)/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)])+((r1g-r2g)Sqrt[r2g (r1g-r3g)]Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(Sqrt[1-EEg^2](r2g+(r1g-r2g)Sin[\[Phi]]^2));

	(4\[Delta]EE)/(1-EEg^2) \[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+1/2 EEg(\[Delta]r1+\[Delta]r2+\[Delta]r3+2\[Delta]\[Rho]r4))\[ScriptCapitalI]rover+\[Delta]EE/(1-EEg^2) \[ScriptCapitalI]r2over+1/2 EEg((2+r1g)\[Delta]r1+(2+r2g)\[Delta]r2+(2+r3g)\[Delta]r3-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+(r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(2(r1g-rm)(r1g-rp)) \[ScriptCapitalI]r1g*\[Delta]r1+(r2g(-2a*Lzg+EEg(2a^2+a^2*r2g+r2g^3)))/(2(r2g-rm)(r2g-rp)) \[ScriptCapitalI]r2g*\[Delta]r2+(r3g(-2a*Lzg+EEg(2a^2+a^2*r3g+r3g^3)))/(2(r3g-rm)(r3g-rp)) \[ScriptCapitalI]r3g*\[Delta]r3+(r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(Sqrt[1-EEg^2](r1g-rm)(r1g-rp)) (Sqrt[(-r2g+r)(-r3g+r)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-r)r]))\[Delta]r1+(r2g(-2a*Lzg+EEg(2a^2+a^2*r2g+r2g^3)))/(Sqrt[1-EEg^2](r2g-rm)(r2g-rp)) ( Sqrt[(r1g-r)(-r3g+r)]/((r1g-r2g)(r2g-r3g)Sqrt[r(-r2g+r)]))\[Delta]r2-dtgd\[Lambda]fun/Sqrt[1-EEg^2] \[Delta]rfunOfrFEPerParRes[r,a,p,e,x]
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


(* ::Subsubsection::Closed:: *)
(*Azimuthal correction to the trajectory - purely oscillatory part*)


\[Delta]\[Phi]funFEPerPar[wr_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Phi],rg,d\[Phi]gd\[Lambda]fun,krghold,krg,\[Gamma]r,ellF,ellE,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
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

	krg=((r1g-r2g) r3g)/(r2g (r1g-r3g));
	
	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result.*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;

	rg=rgICr1gfun[wr,a,p,e,x];

	d\[Phi]gd\[Lambda]fun=((2a*EEg+Lzg(rg-2))rg)/(a^2-2 rg+rg^2);

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]);
	\[ScriptCapitalI]r1g=2 /(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]r1g) ((r2g *ellE-r1g*ellF)/(r1g-r2g) );
	\[ScriptCapitalI]r2g=2 /(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]) (-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g));
	\[ScriptCapitalI]r3g=2 /(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]r3g) (((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]));

	-EEg(1-(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+((2a*EEg+Lzg(-2+r1g))r1g*\[Delta]r1)/(2(r1g-rm)(r1g-rp)) \[ScriptCapitalI]r1g+((2a*EEg+Lzg(-2+r2g))r2g*\[Delta]r2)/(2(r2g-rm)(r2g-rp)) \[ScriptCapitalI]r2g+((2a*EEg+Lzg(-2+r3g))r3g*\[Delta]r3)/(2(r3g-rm)(r3g-rp)) \[ScriptCapitalI]r3g-1/\[CapitalUpsilon]rg (\[Delta]\[CapitalUpsilon]\[Phi]funFE[a,p,e,x]-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]rg \[Delta]\[CapitalUpsilon]rfunFE[a,p,e,x])wr+Sign[\[Pi]-Mod[wr,2 \[Pi]]] ((2a*EEg+Lzg(-2+r1g))r1g*\[Delta]r1)/(Sqrt[1-EEg^2](r1g-rm)(r1g-rp)) (Sqrt[(-r2g+rg)(-r3g+rg)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-rg)rg])) +Sign[\[Pi]-Mod[wr,2 \[Pi]]] ((2a*EEg+Lzg(-2+r2g))r2g*\[Delta]r2)/(Sqrt[1-EEg^2](r2g-rm)(r2g-rp)) ( Sqrt[(r1g-rg)(-r3g+rg)]/((r1g-r2g)(r2g-r3g)Sqrt[rg(-r2g+rg)]))-d\[Phi]gd\[Lambda]fun/Sqrt[1-EEg^2] \[Delta]rfunFEPerParRes[wr,a,p,e,x]
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal correction to the trajectory as function of r*)


\[Delta]\[Phi]funOfrFEPerPar[r_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Phi],d\[Phi]gd\[Lambda]fun,krghold,krg,\[Gamma]r,ellF,ellE,\[ScriptCapitalI],\[ScriptCapitalI]r1g,\[ScriptCapitalI]r2g,\[ScriptCapitalI]r3g,\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over},
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

	krg=((r1g-r2g) r3g)/(r2g (r1g-r3g));
	
	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result.*)
	\[Phi]=ArcSin[Sqrt[(r2g (r1g-r))/((r1g-r2g)r)]];

	d\[Phi]gd\[Lambda]fun=((2a*EEg+Lzg(r-2))r)/(a^2-2 r+r^2);

	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]);
	\[ScriptCapitalI]r1g=2 /(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]r1g) ((r2g *ellE-r1g*ellF)/(r1g-r2g) );
	\[ScriptCapitalI]r2g=2 /(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]) (-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g));
	\[ScriptCapitalI]r3g=2 /(Sqrt[1-EEg^2] Sqrt[r2g(r1g-r3g)]r3g) (((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]));

	-EEg(1-(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+((2a*EEg+Lzg(-2+r1g))r1g*\[Delta]r1)/(2(r1g-rm)(r1g-rp)) \[ScriptCapitalI]r1g+((2a*EEg+Lzg(-2+r2g))r2g*\[Delta]r2)/(2(r2g-rm)(r2g-rp)) \[ScriptCapitalI]r2g+((2a*EEg+Lzg(-2+r3g))r3g*\[Delta]r3)/(2(r3g-rm)(r3g-rp)) \[ScriptCapitalI]r3g+((2a*EEg+Lzg(-2+r1g))r1g*\[Delta]r1)/(Sqrt[1-EEg^2](r1g-rm)(r1g-rp)) (Sqrt[(-r2g+r)(-r3g+r)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-r)r]))+((2a*EEg+Lzg(-2+r2g))r2g*\[Delta]r2)/(Sqrt[1-EEg^2](r2g-rm)(r2g-rp)) ( Sqrt[(r1g-r)(-r3g+r)]/((r1g-r2g)(r2g-r3g)Sqrt[r(-r2g+r)]))-d\[Phi]gd\[Lambda]fun/Sqrt[1-EEg^2] \[Delta]rfunOfrFEPerParRes[r,a,p,e,x]
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
(*Periodic velocity - fixed turning points*)


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunFTPerPar[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,rg},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	\[Delta]r3=\[Delta]r3funFT[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	rg=rgICr2gfun[wr,a,p,e,x];

	-Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(1-EEg^2)rg(r1g-rg)(rg-r2g)(rg-r3g)](\[Delta]r3/(2(rg-r3g))+\[Delta]r41/rg-a/(2rg^2)+(EEg*\[Delta]EE)/(1-EEg^2))-1/2 Sqrt[1-EEg^2](4rg^3-r1g*r2g*r3g-3rg^2(r1g+r2g+r3g)+2rg(r2g*r3g+r1g(r2g+r3g)))\[Delta]rfunFTPerParRes[wr,a,p,e,x]
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


(* ::Subsection::Closed:: *)
(*Periodic velocity - fixed constants of motion*)


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunFCPerPar[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,rg},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFC[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFC[a,p,e,x];
	\[Delta]r3=\[Delta]r3funFC[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	rg=rgICr2gfun[wr,a,p,e,x];

	-Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(1-EEg^2)rg(r1g-rg)(rg-r2g)(rg-r3g)](\[Delta]r1/(2(rg-r1g))+\[Delta]r2/(2(rg-r2g))+\[Delta]r3/(2(rg-r3g))+\[Delta]r41/rg-a/(2rg^2))-1/2 Sqrt[1-EEg^2](4rg^3-r1g*r2g*r3g-3rg^2(r1g+r2g+r3g)+2rg(r2g*r3g+r1g(r2g+r3g)))\[Delta]rfunFCPerParRes[wr,a,p,e,x]
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


(* ::Subsection::Closed:: *)
(*Periodic velocity - fixed eccentricity*)


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunFEPerPar[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,rg},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	rg=rgICr1gfun[wr,a,p,e,x];

	-Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(1-EEg^2)rg(r1g-rg)(rg-r2g)(rg-r3g)](\[Delta]r1/(2(rg-r1g))+\[Delta]r2/(2(rg-r2g))+\[Delta]r3/(2(rg-r3g))+\[Delta]r41/rg-a/(2rg^2)+(EEg*\[Delta]EE)/(1-EEg^2))-1/2 Sqrt[1-EEg^2](4rg^3-r1g*r2g*r3g-3rg^2(r1g+r2g+r3g)+2rg(r2g*r3g+r1g(r2g+r3g)))\[Delta]rfunFEPerParRes[wr,a,p,e,x]
]


\[Delta]vrfunOfrFEPerPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	-Sqrt[(1-EEg^2)r(r1g-r)(r-r2g)(r-r3g)](\[Delta]r1/(2(r-r1g))+\[Delta]r2/(2(r-r2g))+\[Delta]r3/(2(r-r3g))+\[Delta]r41/r-a/(2r^2)+(EEg \[Delta]EE)/(1-EEg^2))-1/2 Sqrt[1-EEg^2](4r^3-r1g*r2g*r3g-3r^2(r1g+r2g+r3g)+2r(r2g*r3g+r1g(r2g+r3g)))\[Delta]rfunOfrFEPerParRes[r,a,p,e,x]
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


\[Delta]vtfunOfrFEPerPar[r_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,r3g,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	\[CapitalDelta]=a^2-2r+r^2;

	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+r^2)+a*EEg(a^2+3r^2))/(r*\[CapitalDelta])+(r(2a^2+a^2*r+r^3)\[Delta]EE)/\[CapitalDelta]-(2a*r*\[Delta]Lz)/\[CapitalDelta];
	ddtgd\[Lambda]drfun=(EEg*r(a^2+3r^2))/\[CapitalDelta]-((r^2-a^2)(-2a*Lzg+EEg*r^3+a^2*EEg(2+r)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunOfrFEPerPar[r,a,p,e,x]+d\[Delta]td\[Lambda]fun
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


\[Delta]v\[Phi]funFEPerPar[wr_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,r3g,rg,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
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

	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+rg)))/(rg*\[CapitalDelta])+(2a*rg*\[Delta]EE)/\[CapitalDelta]+((rg-2)rg*\[Delta]Lz)/\[CapitalDelta];
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(rg-1)-EEg*rg^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunFEPerPar[wr,a,p,e,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


\[Delta]v\[Phi]funOfrFEPerPar[r_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,r3g,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	\[CapitalDelta]=a^2-2r+r^2;

	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+r)))/(r*\[CapitalDelta])+(2a*r*\[Delta]EE)/\[CapitalDelta]+((r-2)r*\[Delta]Lz)/\[CapitalDelta];
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(r-1)-EEg*r^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunOfrFEPerPar[r,a,p,e,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Subsection::Closed:: *)
(*Homoclinic orbits*)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunFEHomPar[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,\[Delta]r1,\[Delta]r2,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[ScriptCapitalI]},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	If[e==0,
		\[Delta]r1
		,
		\[ScriptCapitalI]=-(1/(2Sqrt[(r1g-r2g)r2g]))Log[(Sqrt[r(r1g-r2g)]-Sqrt[(r1g-r)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(r1g-r)r2g])^2];

		Sqrt[(r1g-r)(r-r2g)^2r](Sqrt[(1-EEg^2)])1/Sqrt[(1-EEg^2)](\[Delta]\[CapitalUpsilon]rover\[CapitalUpsilon]rgfunLimFE[a,p,e,x]\[ScriptCapitalI]+(EEg*\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+(r2g*Sqrt[-r+r1g])/(r1g (r1g-r2g)^2*Sqrt[r])\[Delta]r1-\[ScriptCapitalI]/(2(r1g-r2g))\[Delta]r1-(\[ScriptCapitalI](r1g-2r2g)\[Delta]r2)/(2(r1g-r2g)r2g)+1/((r1g-r2g)r2g)Sqrt[-1+r1g/r]\[Delta]r2+(\[ScriptCapitalI]/r2g-(2Sqrt[-r+r1g])/(Sqrt[r]r1g*r2g))\[Delta]r41-a/2(-((2Sqrt[r(-r+r1g)])/(3r*r1g*r2g))(2/r1g+3/r2g+1/r)+\[ScriptCapitalI]/r2g^2))+Sqrt[(r-r2g)^2]((r-r2g)/(r1g-r2g)^2 *\[Delta]r1)+Sqrt[(r1g-r)](Sqrt[r1g-r]/(r1g-r2g)\[Delta]r2)
	]
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunFEHomParRes[r_,a_,p_,e_,x_]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,\[Delta]r1,\[Delta]r2,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[ScriptCapitalI]},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	\[ScriptCapitalI]=-(1/(2Sqrt[(r1g-r2g)r2g]))Log[(Sqrt[r(r1g-r2g)]-Sqrt[(r1g-r)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(r1g-r)r2g])^2];

	(\[Delta]\[CapitalUpsilon]rover\[CapitalUpsilon]rgfunLimFE[a,p,e,x]\[ScriptCapitalI]+(EEg*\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+(r2g*Sqrt[-r+r1g])/(r1g (r1g-r2g)^2*Sqrt[r])\[Delta]r1-\[ScriptCapitalI]/(2(r1g-r2g))\[Delta]r1-(\[ScriptCapitalI](r1g-2r2g)\[Delta]r2)/(2(r1g-r2g)r2g)+1/((r1g-r2g)r2g)Sqrt[-1+r1g/r]\[Delta]r2+(\[ScriptCapitalI]/r2g-(2Sqrt[-r+r1g])/(Sqrt[r]*r1g*r2g))\[Delta]r41-a/2(-((2Sqrt[r(-r+r1g)])/(3r*r1g*r2g))(2/r1g+3/r2g+1/r)+\[ScriptCapitalI]/r2g^2))+((r-r2g)/(Sqrt[(r1g-r)r](r1g-r2g)^2)\[Delta]r1)+(Sqrt[(r1g-r)r]/((r1g-r2g)(r-r2g)r)\[Delta]r2)
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time trajectory *)


\[Delta]tfunFEHomPar[r_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]\[Rho]r4,\[Delta]\[Rho]i4,dtgd\[Lambda]fun,\[ScriptCapitalI]},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
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
	
	If[e==0,
		0,
		dtgd\[Lambda]fun=(r(-2a*Lzg+EEg(2a^2+a^2*r+r^3)))/(a^2-2r+r^2);
		\[ScriptCapitalI]=-(1/(2Sqrt[(r1g-r2g)r2g]))Log[(Sqrt[r(r1g-r2g)]-Sqrt[(r1g-r)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(r1g-r)r2g])^2];
		1/Sqrt[1-EEg^2]((4*\[Delta]EE)/(1-EEg^2)\[ScriptCapitalI]+((2\[Delta]EE)/(1-EEg^2)+1/2EEg(\[Delta]r1+\[Delta]r2+\[Delta]r2+2\[Delta]\[Rho]r4))(2ArcCos[Sqrt[r/r1g]]+r2g*\[ScriptCapitalI])+\[Delta]EE/(1-EEg^2)(Sqrt[r(-r+r1g)]+(r1g+2r2g)ArcCos[Sqrt[r/r1g]]+r2g^2*\[ScriptCapitalI])+1/2EEg((2+r1g)\[Delta]r1+(2+r2g)\[Delta]r2+(2+r2g)\[Delta]r2-\[Delta]\[Rho]i4^2+4\[Delta]\[Rho]r4)\[ScriptCapitalI]+(r1g(-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/(2(r1g-rm)(r1g-rp))((2r2g*Sqrt[-r+r1g])/(r1g (r1g-r2g)^2)1/Sqrt[r]-\[ScriptCapitalI]/(r1g-r2g))\[Delta]r1+(r2g(-2a*Lzg+EEg(2a^2+a^2*r2g+r2g^3)))/(2(r2g-rm)(r2g-rp))(-(((r1g-2r2g)\[ScriptCapitalI])/((r1g-r2g)r2g))+2/((r1g-r2g)r2g)Sqrt[-1+r1g/r])\[Delta]r2+(r1g (-2a*Lzg+EEg(2a^2+a^2*r1g+r1g^3)))/((r1g-rm)(r1g-rp))((r-r2g)/(Sqrt[(r1g-r)r](r1g-r2g)^2))\[Delta]r1+(r2g(-2a*Lzg+EEg(2a^2+a^2*r2g+r2g^3)))/((r2g-rm)(r2g-rp))(Sqrt[(r1g-r)r]/((r1g-r2g)(r-r2g)r))\[Delta]r2-dtgd\[Lambda]fun*\[Delta]rfunFEHomParRes[r,a,p,e,x])
	]
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal correction to the trajectory*)


\[Delta]\[Phi]funFEHomPar[r_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,d\[Phi]gd\[Lambda]fun,\[ScriptCapitalI]},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
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
	If[e==0,
		0,
		d\[Phi]gd\[Lambda]fun=((2a*EEg+Lzg(-2+r))r)/(a^2-2r+r^2);
		\[ScriptCapitalI]=-(1/(2Sqrt[(r1g-r2g)r2g]))Log[(Sqrt[r(r1g-r2g)]-Sqrt[(r1g-r)r2g])^2/(Sqrt[r(r1g-r2g)]+Sqrt[(r1g-r)r2g])^2];

		1/Sqrt[1-EEg^2](-EEg(1-(Lzg*\[Delta]EE)/(1-EEg^2))\[ScriptCapitalI]+\[Delta]Lz*\[ScriptCapitalI]+((2a*EEg+Lzg(-2+r1g))r1g)/(2(r1g-rm)(r1g-rp))((2r2g*Sqrt[-r+r1g])/(r1g (r1g-r2g)^2) 1/Sqrt[r]-\[ScriptCapitalI]/(r1g-r2g))\[Delta]r1+((2a*EEg+Lzg(-2+r2g))r2g)/(2(r2g-rm)(r2g-rp)) (-(((r1g-2r2g)\[ScriptCapitalI])/((r1g-r2g)r2g))+2/((r1g-r2g) r2g)Sqrt[-1+r1g/r])\[Delta]r2+((2a*EEg+Lzg(-2+r1g))r1g)/((r1g-rm)(r1g-rp))((r-r2g)/(Sqrt[(r1g-r)r](r1g-r2g) ^2))\[Delta]r1+((2a*EEg+Lzg(-2+r2g))r2g)/((r2g-rm)(r2g-rp))(Sqrt[(r1g-r)r]/((r1g-r2g)(r-r2g)r))\[Delta]r2-d\[Phi]gd\[Lambda]fun*\[Delta]rfunFEHomParRes[r,a,p,e,x])
	]
]


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunFEHomPar[r_,a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,\[Phi]},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	\[Delta]r1=\[Delta]r1funFE[a,p,e,x];
	\[Delta]r2=\[Delta]r2funFE[a,p,e,x];
	\[Delta]r3=\[Delta]r2;
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	
	If[e==0,
		0,
		-Sqrt[(1-EEg^2)r(r1g-r)(r-r2g)^2](\[Delta]r1/(2(r-r1g))+\[Delta]r2/(r-r2g)+\[Delta]r41/r-a/(2r^2)+(EEg \[Delta]EE)/(1-EEg^2))-1/2 Sqrt[1-EEg^2](r-r2g)(r1g(r2g-3r)-2(r2g-2r)r)\[Delta]rfunFEHomParRes[r,a,p,e,x]
	]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity*)


\[Delta]vtfunFEHomPar[r_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
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


\[Delta]v\[Phi]funFEHomPar[r_,a_,p_,e_,x_]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

	\[CapitalDelta]=a^2-2r+r^2;

	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+r)))/(r*\[CapitalDelta])+(2a*r*\[Delta]EE)/\[CapitalDelta]+((r-2)r*\[Delta]Lz)/\[CapitalDelta];

	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(r-1)-EEg*r^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunFEHomPar[r,a,p,e,x]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Section::Closed:: *)
(*Spin corrections to the orbit  - orthogonal component of the spin *)


(* ::Subsection::Closed:: *)
(*Spin precession phase*)


(* ::Subsubsection::Closed:: *)
(*Periodic orbits*)


\[Psi]pref[a_,p_,e_,x_]:=Module[{EEg,Lzg,r1g,r2g,r3g,\[Phi],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun,krghold,krg,Lzred,\[Gamma]r,\[Psi]freq},
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
	\[Psi]freq=r2g^2/(r2g^2+Abs[Lzred]^2)\[ScriptCapitalI]+Sqrt[Abs[Lzred]]/4*Abs[1/(Sqrt[r1g+I*Abs[Lzred]](-I*r2g+Abs[Lzred]))Log[(Sqrt[r(r1g+I*Abs[Lzred])]-Sqrt[I*(r-r1g)]*Sqrt[Abs[Lzred]])^2/(Sqrt[r(r1g+I*Abs[Lzred])]+Sqrt[I*(r-r1g)]Sqrt[Abs[Lzred]])^2]-1/(Sqrt[r1g-I*Abs[Lzred]](r2g-I*Abs[Lzred])) Log[(Sqrt[r(r1g-I*Abs[Lzred])]-Sqrt[I(-r+r1g)] Sqrt[Abs[Lzred]])^2/(Sqrt[r(r1g-I*Abs[Lzred])]+Sqrt[I(-r+r1g)] Sqrt[Abs[Lzred]])^2]];
	
	RealSign[Lzred]/Sqrt[(1-EEg^2)](a+EEg*Lzred)Re[\[Psi]freq]
]


(* ::Subsection::Closed:: *)
(*Polar trajectory*)


(* ::Subsubsection::Closed:: *)
(*Periodic orbits*)


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


\[Delta]zfunICr1gPer[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,Lzred,\[CapitalUpsilon]rg,r1g,r2g,r3g,\[Phi],rg,krghold,krg,\[Psi]p,\[CapitalUpsilon]p},
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


(* ::Section::Closed:: *)
(*Near equatorial orbits - fixed turning points*)


(* ::Subsection::Closed:: *)
(*Analytic solutions*)


(* ::Subsubsection::Closed:: *)
(*Periodic solutions*)


KerrNearEqSpinOrbitCorrFTPer[a_, p_, e_, x_]:=Module[{EEg,Lzg,\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g,\[CapitalUpsilon]p,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]zs,\[CapitalUpsilon]\[Phi]s,dtrgd\[Lambda],d\[Phi]rgd\[Lambda],drgd\[Lambda],\[CapitalDelta]tspar,\[CapitalDelta]\[Phi]spar,\[Delta]rpar,
	\[Psi]phase,\[Delta]zort,\[CapitalDelta]trg,\[CapitalDelta]\[Phi]rg,rg,dtspard\[Lambda],drspard\[Lambda],dzspard\[Lambda],d\[Phi]spard\[Lambda]},
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
	
	dtrgd\[Lambda][wr_]:=VtrgICr2gfun[wr,a,p,e,x];
	d\[Phi]rgd\[Lambda][wr_]:=V\[Phi]rgICr2gfun[wr,a,p,e,x];
	drgd\[Lambda][wr_]:=drgd\[Lambda]fun[wr,a,p,e,x];

	\[CapitalDelta]trg[wr_]:=tgICr2gfun[wr,a,p,e,x];
	\[CapitalDelta]\[Phi]rg[wr_]:=\[Phi]gICr2gfun[wr,a,p,e,x];
	rg[wr_]:=rgICr2gfun[wr,a,p,e,x];
	
	\[CapitalDelta]tspar[wr_]:=\[Delta]tfunFTPerPar[wr,a,p,e,x];
	\[CapitalDelta]\[Phi]spar[wr_]:=\[Delta]\[Phi]funFTPerPar[wr,a,p,e,x];
	\[Delta]rpar[wr_]:=\[Delta]rfunFTPerPar[wr,a,p,e,x];
	\[Psi]phase[wr_]:=\[Psi]pICr2g[wr,a,p,e,x];
	\[Delta]zort[wr_]:=\[Delta]zfunICr2gPer[wr,a,p,e,x];
	  
	dtspard\[Lambda][wr_]:=\[Delta]vtfunFTPerPar[wr,a,p,e,x];
	d\[Phi]spard\[Lambda][wr_]:=\[Delta]v\[Phi]funFTPerPar[wr,a,p,e,x];
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
(*Periodic solutions*)


KerrNearEqSpinOrbitCorrFCPer[a_, p_, e_, x_]:=Module[{EEg,Lzg,\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g,\[CapitalUpsilon]p,\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]zs,\[CapitalUpsilon]\[Phi]s,dtrgd\[Lambda],d\[Phi]rgd\[Lambda],drgd\[Lambda],\[CapitalDelta]tspar,\[CapitalDelta]\[Phi]spar,\[Delta]rpar,\[Psi]phase,\[Delta]zort,\[CapitalDelta]trg,\[CapitalDelta]\[Phi]rg,rg,dtspard\[Lambda],drspard\[Lambda],
	dzspard\[Lambda],d\[Phi]spard\[Lambda]},
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
	
	dtrgd\[Lambda][wr_]:=VtrgICr2gfun[wr,a,p,e,x];
	d\[Phi]rgd\[Lambda][wr_]:=V\[Phi]rgICr2gfun[wr,a,p,e,x];
	drgd\[Lambda][wr_]:=drgd\[Lambda]fun[wr,a,p,e,x];

	\[CapitalDelta]trg[wr_]:=tgICr2gfun[wr,a,p,e,x];
	\[CapitalDelta]\[Phi]rg[wr_]:=\[Phi]gICr2gfun[wr,a,p,e,x];
	rg[wr_]:=rgICr2gfun[wr,a,p,e,x];
	
	\[CapitalDelta]tspar[wr_]:=\[Delta]tfunFCPerPar[wr,a,p,e,x];
	\[CapitalDelta]\[Phi]spar[wr_]:=\[Delta]\[Phi]funFCPerPar[wr,a,p,e,x];
	\[Delta]rpar[wr_]:=\[Delta]rfunFCPerPar[wr,a,p,e,x];
	\[Psi]phase[wr_]:=\[Psi]pICr2g[wr,a,p,e,x];
	\[Delta]zort[wr_]:=\[Delta]zfunICr2gPer[wr,a,p,e,x];
	  
	dtspard\[Lambda][wr_]:=\[Delta]vtfunFCPerPar[wr,a,p,e,x];
	d\[Phi]spard\[Lambda][wr_]:=\[Delta]v\[Phi]funFCPerPar[wr,a,p,e,x];
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
(*Periodic solutions*)


KerrNearEqSpinOrbitCorrFEPer[a_, p_, e_, x_]:=Module[{EEg,Lzg,\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g,\[CapitalUpsilon]p,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]zs,\[CapitalUpsilon]\[Phi]s,dtrgd\[Lambda],d\[Phi]rgd\[Lambda],drgd\[Lambda],\[CapitalDelta]tspar,\[Delta]tOfrpar,\[CapitalDelta]\[Phi]spar,\[Delta]\[Phi]Ofrpar,\[Delta]rpar,\[Delta]rOfrpar,\[Psi]phase,\[Psi]phaseOfr,\[Delta]zort,\[Delta]zOfrort,\[CapitalDelta]trg,\[CapitalDelta]\[Phi]rg,rg,dtspard\[Lambda],\[Delta]vtOfrpar,drspard\[Lambda],\[Delta]vrOfrpar,
	dzspard\[Lambda],d\[Phi]spard\[Lambda],\[Delta]v\[Phi]Ofrpar},
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
	
	dtrgd\[Lambda][wr_]:=VtrgICr1gfun[wr,a,p,e,x];
	d\[Phi]rgd\[Lambda][wr_]:=V\[Phi]rgICr1gfun[wr,a,p,e,x];
	drgd\[Lambda][wr_]:=drgd\[Lambda]fun[wr,a,p,e,x];

	\[CapitalDelta]trg[wr_]:=tgICr1gfun[wr,a,p,e,x];
	\[CapitalDelta]\[Phi]rg[wr_]:=\[Phi]gICr1gfun[wr,a,p,e,x];
	rg[wr_]:=rgICr1gfun[wr,a,p,e,x];
	
	\[CapitalDelta]tspar[wr_]:=\[Delta]tfunFEPerPar[wr,a,p,e,x];
	\[Delta]tOfrpar[r_]:=\[Delta]tfunOfrFEPerPar[r,a,p,e,x];
	\[CapitalDelta]\[Phi]spar[wr_]:=\[Delta]\[Phi]funFEPerPar[wr,a,p,e,x];
	\[Delta]\[Phi]Ofrpar[r_]:=\[Delta]\[Phi]funOfrFEPerPar[r,a,p,e,x];
	\[Delta]rpar[wr_]:=\[Delta]rfunFEPerPar[wr,a,p,e,x];
	\[Delta]rOfrpar[r_]:=\[Delta]rfunOfrFEPerPar[r,a,p,e,x];
	\[Psi]phase[wr_]:=\[Psi]pICr1g[wr,a,p,e,x];
	\[Psi]phaseOfr[r_]:=\[Psi]pOfrICr1g[r,a,p,e,x];
	\[Delta]zort[wr_]:=\[Delta]zfunICr1gPer[wr,a,p,e,x];
	\[Delta]zOfrort[wr_]:=\[Delta]zfunOfrICr1gPer[wr,a,p,e,x];
	  
	dtspard\[Lambda][wr_]:=\[Delta]vtfunFEPerPar[wr,a,p,e,x];
	\[Delta]vtOfrpar[r_]:=\[Delta]vtfunOfrFEPerPar[r,a,p,e,x];
	d\[Phi]spard\[Lambda][wr_]:=\[Delta]v\[Phi]funFEPerPar[wr,a,p,e,x];
	\[Delta]v\[Phi]Ofrpar[r_]:=\[Delta]v\[Phi]funOfrFEPerPar[r,a,p,e,x];
	drspard\[Lambda][wr_]:=\[Delta]vrfunFEPerPar[wr,a,p,e,x];
	\[Delta]vrOfrpar[r_]:=\[Delta]vrfunOfrFEPerPar[r,a,p,e,x];
			
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
	 "\[Delta]tOfrpar"->Function[{r},\[Delta]tOfrpar[r]],
	 "\[Delta]rpar"->Function[{wr},\[Delta]rpar[wr]],
	 "\[Delta]rOfrpar"->Function[{r},\[Delta]rOfrpar[r]],
	 "\[Psi]p"->Function[wr,\[Psi]phase[wr]],
	 "\[Psi]pOfr"->Function[r,\[Psi]phaseOfr[r]],
	 "\[Delta]zort"->Function[{wr},\[Delta]zort[wr]],
	 "\[Delta]zOfrort"->Function[{r},\[Delta]zOfrort[r]],
	 "\[CapitalDelta]\[Delta]\[Phi]par"->Function[{wr},\[CapitalDelta]\[Phi]spar[wr]],
	 "\[Delta]\[Phi]Ofrpar"->Function[{r},\[Delta]\[Phi]Ofrpar[r]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{wr},dtspard\[Lambda][wr]],
	 "\[Delta]vtOfrpar"->Function[{r},\[Delta]vtOfrpar[r]],
	 "\[Delta]vrpar"->Function[{wr},drspard\[Lambda][wr]],
	 "\[Delta]vrOfrpar"->Function[{r},\[Delta]vrOfrpar[r]],
	 "\[Delta]v\[Phi]par"->Function[{wr},d\[Phi]spard\[Lambda][wr]],
	 "\[Delta]v\[Phi]Ofrpar"->Function[{r},\[Delta]v\[Phi]Ofrpar[r]]
	|>
]


(* ::Subsubsection::Closed:: *)
(*Homoclinic solutions*)


KerrNearEqSpinOrbitCorrFEHom[a_, p_, e_, x_]:=Module[{EEg,Lzg,\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g,\[CapitalUpsilon]p,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]zs,\[CapitalUpsilon]\[Phi]s,dtrgd\[Lambda],d\[Phi]rgd\[Lambda],drgd\[Lambda],tspar,\[Phi]spar,\[Delta]rpar,\[Psi]phase,\[Delta]zort,trg,\[Phi]rg,rg,dtspard\[Lambda],drspard\[Lambda],dzspard\[Lambda],d\[Phi]spard\[Lambda]},
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
	
	dtrgd\[Lambda][r_]:=VtrgfunHom[r,a,p,e,x];
	d\[Phi]rgd\[Lambda][r_]:=V\[Phi]rgfunHom[r,a,p,e,x];
	drgd\[Lambda][r_]:=drgd\[Lambda]funHom[r,a,p,e,x];

	trg[r_]:=tgfunHom[r,a,p,e,x];
	\[Phi]rg[r_]:=\[Phi]gfunHom[r,a,p,e,x];
	rg[\[Lambda]_]:=rgfunHom[\[Lambda],a,p,e,x];
	
	tspar[r_]:=\[Delta]tfunFEHomPar[r,a,p,e,x];
	\[Phi]spar[r_]:=\[Delta]\[Phi]funFEHomPar[r,a,p,e,x];
	\[Delta]rpar[r_]:=\[Delta]rfunFEHomPar[r,a,p,e,x];
	\[Psi]phase[r_]:=\[Psi]pHom[r,a,p,e,x];
	\[Delta]zort[r_]:=\[Delta]zfunHom[r,a,p,e,x];
	  
	dtspard\[Lambda][r_]:=\[Delta]vtfunFEHomPar[r,a,p,e,x];
	d\[Phi]spard\[Lambda][r_]:=\[Delta]v\[Phi]funFEHomPar[r,a,p,e,x];
	drspard\[Lambda][r_]:=\[Delta]vrfunFEHomPar[r,a,p,e,x];
			
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
	 "trg"->Function[{r},trg[r]],
	 "\[Phi]rg"->Function[{r},\[Phi]rg[r]],
	 (*Keys geodesic radial trajectory*)
	 "rg"->Function[{wr},rg[wr]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{r},dtrgd\[Lambda][r]],
	 "vrg"->Function[{r},drgd\[Lambda][r]],
	 "v\[Phi]g"->Function[{r},d\[Phi]rgd\[Lambda][r]],
	 (*Keys corrections trajectory*)
	 "\[Delta]tpar"->Function[{r},tspar[r]],
	 "\[Delta]rpar"->Function[{r},\[Delta]rpar[r]],
	 "\[Psi]p"->Function[r,\[Psi]phase[r]],
	 "\[Delta]zort"->Function[{r},\[Delta]zort[r]],
	 "\[Delta]\[Phi]par"->Function[{r},\[Phi]spar[r]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{wr},dtspard\[Lambda][wr]],
	 "\[Delta]vrpar"->Function[{wr},drspard\[Lambda][wr]],
	 "\[Delta]v\[Phi]par"->Function[{wr},d\[Phi]spard\[Lambda][wr]],
	 (*Shift separatrix*)
	 "\[Delta]p"->\[Delta]pfun[a,p,e,x]
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
