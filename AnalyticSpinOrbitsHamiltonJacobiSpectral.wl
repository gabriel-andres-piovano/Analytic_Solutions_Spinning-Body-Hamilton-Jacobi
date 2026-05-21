(* ::Package:: *)

(* ::Section:: *)
(*Begin package*)


BeginPackage["AnalyticSpinOrbitsHamiltonJacobiSpectral`"];


(* ::Text:: *)
(*If you make use of this package, please acknowledge "Piovano" arXiv:2510.09597 (https://arxiv.org/abs/2510.09597 )*)


(* ::Text:: *)
(*IMPORTANT NOTE: the following package include all the functions to compute the spin-corrections to the orbits (trajectories and velocities), constants of motion and frequencies for nearly equatorial orbits using spectral methods.*)
(*The package provide the contributions to both the parallel and orthogonal components of the secondary spin (the latter are not implemented with a spectral solver yet). *)
(*The spin-corrections are available in the fixed frequency (or "FF"), fixed turning points (or "FT"), and fixed eccentricity (or "FE") spin-gauges*)


KerrNearEqSpinOrbitCorrFFSpectral::usage = "KerrNearEqSpinOrbitCorrFFPerFourier[a, p, e, x, tolerance] calculates the linear corrections to periodic orbits in the fixed frequency spin-gauge. The coordinate time and azimuthal trajectories are expanded in Fourier series. The initial radius at \[Lambda]=0 is the geodesic periastron.";


KerrNearEqSpinOrbitCorrFESpectral::usage = "KerrNearEqSpinOrbitCorrFEPerFourier[a, p, e, x, tolerance] calculates the linear corrections to periodic orbits in the fixed eccentricity spin-gauge. The coordinate time and azimuthal trajectories are expanded in Fourier series. The initial radius at \[Lambda]=0 is the geodesic periastron.";


KerrNearEqSpinOrbitCorrFTSpectral::usage = "KerrNearEqSpinOrbitCorrFTPerFourier[a, p, e, x, tolerance] calculates the linear corrections to periodic orbits in the fixed turning points spin-gauge. The coordinate time and azimuthal trajectories are expanded in Fourier series. The initial radius at \[Lambda]=0 is the geodesic periastron.";


KerrNearEqFrequencyCorrFF::usage = "KerrNearEqFrequencyCorrFF[a, p, e, x] calculates the analytic geodesic frequencies and their linear spin corrections in the fixed frequency spin-gauge.";


KerrNearEqFrequencyCorrFE::usage = "KerrNearEqFrequencyCorrFE[a, p, e, x] calculates the analytic geodesic frequencies and their linear spin corrections in the fixed eccentricity spin-gauge.";


KerrNearEqFrequencyCorrFT::usage = "KerrNearEqFrequencyCorrFT[a, p, e, x] calculates the analytic geodesic frequencies and their linear spin corrections in the fixed turning point spin-gauge.";


Begin["`Private`"];


(* ::Text:: *)
(*NOTE: still need to find correct Fourier expansion spin-correction radial velocity fixed turning points.*)


(* ::Subsection::Closed:: *)
(*Building Fourier series*)


(* ::Text:: *)
(*Direct Fourier expansion of velocities*)


fourierCos[wr_,growthrate_,coeff_]:=Module[{dim},
	dim=Length[coeff];
	growthrate+Sum[2Cos[n*wr]coeff[[n]],{n,1,dim}]
];


fourierd\[Delta]rd\[Lambda][wr_,\[CapitalUpsilon]rg_,\[Delta]\[CapitalUpsilon]r_,coeff_,\[Delta]coeff_]:=Module[{dim0,dim1},
	dim0=Length[coeff];
	dim1=Length[\[Delta]coeff];
	\[CapitalUpsilon]rg*Sum[2Sin[n*wr]n*\[Delta]coeff[[n]],{n,1,dim1}]+\[Delta]\[CapitalUpsilon]r*Sum[2Sin[n*wr]n*coeff[[n]],{n,1,dim0}]
]


(* ::Text:: *)
(*For geodesic trajectories*)


\[CapitalDelta]integratedFunc[wr_,coeff_]:=Module[{dim},
	dim=Length[coeff];
	Sum[2Sin[n*wr]coeff[[n]],{n,1,dim}]
];


(* ::Text:: *)
(*For spin-corrections to the trajectories*)


\[CapitalDelta]\[Delta]integratedFunc[wr_,\[CapitalUpsilon]rg_,\[Delta]\[CapitalUpsilon]r_,coeff_,\[Delta]coeff_]:=Module[{dim0,dim1},
	dim0=Length[coeff];
	dim1=Length[\[Delta]coeff];
	Sum[2Sin[n*wr]\[Delta]coeff[[n]],{n,1,dim1}]-(\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg)*Sum[2Sin[n*wr]coeff[[n]],{n,1,dim0}]
]


(* ::Subsubsection::Closed:: *)
(*Subroutine for spectral integration*)


CoeffsFourier[func_,\[CapitalUpsilon]rg_,prec_]:=
Module[{N0,eps,nInt,steps,wrlist,sampledFunc,ExpniTable,coeffs,relerr,i,nmax,growthrate,coeffSpec,coeffSpecInt},

	If[prec==MachinePrecision,eps=15;,eps=prec;];
	
	N0=2^4;
	nInt[n_]:=2^n*N0;
	steps[n_]:=4nInt[n];
	wrlist[n_]:=N[Table[wr,{wr,2Pi/(2*steps[n]),2Pi,2Pi/steps[n]}]];
	
	sampledFunc[0]=func[wrlist[0]];
	sampledFunc[n_]:=sampledFunc[n]=func[wrlist[n]];
	
	ExpniTable[0]=Table[N[Exp[2Pi*I*j*(i-1/2)/steps[0]],prec],{j,-nInt[0],nInt[0]},{i,1,steps[0]}];
	ExpniTable[n_]:=ExpniTable[n]=Table[N[Exp[2Pi*I*j*(i-1/2)/steps[n]],prec],{j,-nInt[n],nInt[n]},{i,1,steps[n]}];
	
	coeffs[0]=Re[(ExpniTable[0] . sampledFunc[0])/steps[0]];
	coeffs[n_]:=coeffs[n]=Re[(ExpniTable[n] . sampledFunc[n])/steps[n]];
	
	relerr[0]=Min[1-Abs[coeffs[0][[-1]]/coeffs[0][[nInt[0]+1]]],1-Abs[coeffs[0][[-2]]/coeffs[0][[nInt[0]+1]]]];
	relerr[n_]:=Min[1-Abs[coeffs[n][[-1]]/coeffs[n][[nInt[n]+1]]],1-Abs[coeffs[n][[-2]]/coeffs[n][[nInt[n]+1]]]];
	
	i=1;
	(*Keep evaluating points until the relative error is smaller than then tolerance*)
	While[relerr[i]<10^(-eps)&&i<=10,i++];
	
	nmax=(Length[coeffs[i]]-1)/2;
	growthrate=coeffs[i][[nmax+1]];
	(*With arbitray precision, several terms have zero significant digits, so they can be removed by imposing a larger error tolerance*)
	If[prec==MachinePrecision,
		coeffSpec=Chop[Table[coeffs[i][[n+nmax+1]],{n,1,nmax}],10^(-eps+1)];
		(*Coefficients for spectral integration in Mino-time*)
		coeffSpecInt=Chop[Table[coeffs[i][[n+nmax+1]]/(n*\[CapitalUpsilon]rg),{n,1,nmax}],10^(-eps+1)];
		,
		coeffSpec=Chop[Table[coeffs[i][[n+nmax+1]],{n,1,nmax}],10^(-eps+5)];
		(*Coefficients for spectral integration in Mino-time*)
		coeffSpecInt=Chop[Table[coeffs[i][[n+nmax+1]]/(n*\[CapitalUpsilon]rg),{n,1,nmax}],10^(-eps+5)];
	];
	
	Remove[sampledFunc,wrlist,steps,ExpniTable,nInt,i,coeffs];
	{growthrate,coeffSpec,coeffSpecInt}
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

	(-EEg*gginc[r2g]+Lzgzaux*RealSign[xg])/((1+e)hhinc[r2g]^2) (2(1+e-p))/(1+e)-(gginc[r2g] dEEdp)/hhinc[r2g]-EEg/hhinc[r2g](2a)/(1+e)+RealSign[xg]/hhinc[r2g]((hhinc[r2g]*EEg*2ffinc[r2g]*dEEdp)/(2Lzgzaux)-hhinc[r2g]/(2Lzgzaux) (2p(a^2(1+e)^2+p(-3-3e+2p)))/(1+e)^4+(hhinc[r2g]EEg^2)/(2Lzgzaux) ((4p^3+2a^2(1+e)^2(1+e+p))/(1+e)^4)+(2 EEg^2*gginc[r2g])/(2Lzgzaux)(2a)/(1+e)+(2gginc[r2g]^2EEg*dEEdp)/(2Lzgzaux)-(-ddinc[r2g]+EEg^2*ffinc[r2g])/(2Lzgzaux)(2(1+e-p))/(1+e)^2)
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
(*Geodesics*)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


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


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


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

	Sqrt[1-EEg^2]((r1g-r2g)(r1g-r3g)(r2g-r3g)jCN*jSN*Sqrt[r2g(r1g-r3g)+(-r1g+r2g)r3g*jSN^2])/(-r1g+r3g+(r1g-r2g)jSN^2)^2
]


(* ::Subsubsection::Closed:: *)
(*Coordinate-time velocity*)


VtrgICr1gfun[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,rg},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	rg=rgICr1gfun[wr,a,p,e,x];

	EEg((rg^2+a^2)^2/(rg^2-2rg+a^2))-(2a*rg)/(rg^2-2rg+a^2)Lzg-EEg*a^2
]


VtrgICr2gfun[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,rg},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	rg=rgICr2gfun[wr,a,p,e,x];

	EEg((rg^2+a^2)^2/(rg^2-2rg+a^2))-(2a*rg)/(rg^2-2rg+a^2)Lzg-EEg*a^2
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


V\[Phi]rgICr1gfun[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,rg},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	rg=rgICr1gfun[wr,a,p,e,x];

	a/(rg^2-2rg+a^2)(EEg(rg^2+a^2)-a*Lzg)-a*EEg+Lzg
]


V\[Phi]rgICr2gfun[wr_,a_,p_,e_,x_]:=Module[{EEg,Lzg,rg},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	rg=rgICr2gfun[wr,a,p,e,x];

	a/(rg^2-2rg+a^2)(EEg(rg^2+a^2)-a*Lzg)-a*EEg+Lzg
]


(* ::Subsection::Closed:: *)
(*Mino-time geodesic frequencies*)


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


(* ::Subsection::Closed:: *)
(*BL geodesic frequencies*)


(* ::Subsubsection::Closed:: *)
(*Radial frequency*)


\[CapitalOmega]rgfun[a_,p_,e_,xg_]:=\[CapitalUpsilon]rgfun[a,p,e,xg]/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])


(* ::Subsubsection::Closed:: *)
(*Azimuthal geodesic frequency*)


\[CapitalOmega]\[Phi]gfun[a_,p_,e_,xg_]:=(\[CapitalUpsilon]\[Phi]grfun[a,p,e,xg]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,xg])/(\[CapitalUpsilon]tgrfun[a,p,e,xg]+\[CapitalUpsilon]tgzfun[a,p,e,xg])


(* ::Section::Closed:: *)
(*Shifts to the constants of motion*)


(* ::Subsection::Closed:: *)
(*Fixed turning points*)


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


\[Delta]r41fun[a_,p_,e_,x_]:=Module[{EEg,Lzg,Lzred},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	Lzred=Lzg-a*EEg;

	-a(a^2(1-EEg^2)+Lzg^2)/(4Lzred^2)
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


(* ::Subsection::Closed:: *)
(*Fixed frequency *)


\[Delta]turnpointsfunFF[a_,p_,e_,x_]:=Module[{\[Delta]\[CapitalOmega]rFT,\[Delta]\[CapitalOmega]\[Phi]FT,jacob},
	\[Delta]\[CapitalOmega]rFT=(\[Delta]\[CapitalUpsilon]rfunFT[a,p,e,x]-\[CapitalOmega]rgfun[a,p,e,x] \[Delta]\[CapitalUpsilon]tfunFT[a,p,e,x])/(\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x]);
	\[Delta]\[CapitalOmega]\[Phi]FT=(\[Delta]\[CapitalUpsilon]\[Phi]funFT[a,p,e,x]-\[CapitalOmega]\[Phi]gfun[a,p,e,x]\[Delta]\[CapitalUpsilon]tfunFT[a,p,e,x])/(\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x]);
	 
	jacob=-{{d\[CapitalOmega]rgdr1gfun[a,p,e,x],d\[CapitalOmega]rgdr2gfun[a,p,e,x]},{d\[CapitalOmega]\[Phi]gdr1gfun[a,p,e,x],d\[CapitalOmega]\[Phi]gdr2gfun[a,p,e,x]}};

	Inverse[jacob] . {\[Delta]\[CapitalOmega]rFT,\[Delta]\[CapitalOmega]\[Phi]FT}
]


\[Delta]r3funFF[a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,Lzred,\[Delta]r1,\[Delta]r2},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	Lzred=Lzg-a*EEg;
	{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;

	a (a^2(1-EEg^2)+Lzg^2)/(2Lzred^2)+(4EEg*\[Delta]EE)/(1-EEg^2)^2-(\[Delta]r1+\[Delta]r2)
]


(* ::Subsection::Closed:: *)
(*Fixed eccentricity*)


\[Delta]r1funFE[a_,p_,e_,x_]:=\[Delta]pfun[a,p,e,x]/(1-e)


\[Delta]r2funFE[a_,p_,e_,x_]:=\[Delta]pfun[a,p,e,x]/(1+e)


(* ::Subsection::Closed:: *)
(*Fixed turning points*)


\[Delta]r3funFT[a_,p_,e_,x_]:=Module[{EEg,Lzg,\[Delta]EE,Lzred},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	Lzred=Lzg-a*EEg;

	a (a^2(1-EEg^2)+Lzg^2)/(2Lzred^2)+(4EEg*\[Delta]EE)/(1-EEg^2)^2
]


(* ::Section::Closed:: *)
(*Shift frequencies*)


(* ::Subsection::Closed:: *)
(*Shift radial frequency*)


(* ::Subsubsection::Closed:: *)
(*Fixed frequency*)


\[Delta]\[CapitalUpsilon]rfunFF[a_,p_,e_,x_]:=Module[{\[CapitalUpsilon]rg,EEg,Lzg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,Yint,krg,ellK,ellE},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	EEg=EEgfun[a,p,e,x];
	{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];
	\[Delta]r3=\[Delta]r3funFF[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	
	Yint=Sqrt[(r1g-r3g)r2g];
	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));

	ellK=EllipticK[krg];
	ellE=EllipticE[krg];

	-2\[CapitalUpsilon]rg^2/(2\[Pi])1/Sqrt[1-EEg^2] 2/Yint((EEg \[Delta]EE)/(1-EEg^2)ellK+(r2g*ellE-r1g*ellK)/(2r1g(r1g-r2g))\[Delta]r1+1/2(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellK/(r1g-r2g))\[Delta]r2+1/2*\[Delta]r3(1/r3g(r2g/(r2g-r3g)ellE-ellK))+\[Delta]r41((r1g*ellK-(r1g-r3g)ellE)/(r1g*r3g))-1/2*a(-(1/(3r1g*r2g*r3g))((r1g+r2g+r3g)ellK-2(r2g*r3g+r1g(r2g+r3g))((-((r1g-r3g)ellE)+r1g*ellK)/(r1g*r3g)))))
]


(* ::Subsubsection::Closed:: *)
(*Fixed eccentricity*)


\[Delta]\[CapitalUpsilon]rfunFE[a_,p_,e_,x_]:=Module[{\[CapitalUpsilon]rg,EEg,Lzg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,Yint,krg,ellK,ellE},
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
	Yint=Sqrt[(r1g-r3g)r2g];
	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));

	ellK=EllipticK[krg];
	ellE=EllipticE[krg];

	-2\[CapitalUpsilon]rg^2/(2\[Pi])1/Sqrt[1-EEg^2] 2/Yint((EEg \[Delta]EE)/(1-EEg^2)ellK+(r2g*ellE-r1g*ellK)/(2r1g(r1g-r2g))\[Delta]r1+1/2(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellK/(r1g-r2g))\[Delta]r2+1/2*\[Delta]r3(1/r3g(r2g/(r2g-r3g)ellE-ellK))+\[Delta]r41((r1g*ellK-(r1g-r3g)ellE)/(r1g*r3g))-1/2*a(-(1/(3r1g*r2g*r3g))((r1g+r2g+r3g)ellK-2(r2g*r3g+r1g(r2g+r3g))((-((r1g-r3g)ellE)+r1g*ellK)/(r1g*r3g)))))
]


(* ::Subsubsection::Closed:: *)
(*Fixed turning points*)


\[Delta]\[CapitalUpsilon]rfunFT[a_,p_,e_,x_]:=Module[{\[CapitalUpsilon]rg,EEg,Lzg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,Yint,krg,ellK,ellE},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	\[Delta]r3=\[Delta]r3funFT[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	Yint=Sqrt[(r1g-r3g)r2g];
	krg=((r1g-r2g)r3g)/(r2g(r1g-r3g));

	ellK=EllipticK[krg];
	ellE=EllipticE[krg];

	-2\[CapitalUpsilon]rg^2/(2\[Pi])/Sqrt[1-EEg^2]2/Yint((EEg*\[Delta]EE)/(1-EEg^2)ellK+\[Delta]r3/2(1/r3g(r2g/(r2g-r3g)ellE-ellK))+\[Delta]r41((r1g*ellK-(r1g-r3g)ellE)/(r1g*r3g))-a/2(-(1/(3r1g*r2g*r3g))((r1g+r2g+r3g)ellK-2(r2g*r3g+r1g(r2g+r3g))((-((r1g-r3g)ellE)+r1g*ellK)/(r1g*r3g)))))
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


(* ::Section::Closed:: *)
(*Spin corrections to the orbit - parallel component of the spin, fixed frequency*)


(* ::Subsection::Closed:: *)
(*Trajectory*)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunFFPerPar[wr_,a_,p_,e_,x_,{\[Delta]r1_,\[Delta]r2_},{\[CapitalUpsilon]rg_,\[Delta]\[CapitalUpsilon]r_}]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,EEg,\[Delta]EE,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE,rg},
	EEg=EEgfun[a,p,e,x];
	(*\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];*)
	(*\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunFF[a,p,e,x];*)

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	(*{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];*)
	\[Delta]r3=\[Delta]r3funFF[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
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
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobiSpectral`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobiSpectral`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobiSpectral`Private`*"],True]]]];

	((jSN*jCN*r1g(r1g-r2g)r2g*Sqrt[r1g(r2g-jSN^2*r3g)-(1-jSN^2)r2g*r3g])/(jSN^2(r1g-r2g)+r2g)^2)(\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg*(2ellF)/Yint+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Yint+((r2g*ellE-r1g*ellF)\[Delta]r1)/(Yint*r1g(r1g-r2g))+1/Yint(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g))\[Delta]r2+1/Yint 1/r3g(((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]))\[Delta]r3+2((-(r1g-r3g)ellE+r1g*ellF)/(r1g*Yint*r3g))\[Delta]r41-a/2(-((4ellE*Yint(r2g*r3g+r1g(r2g+r3g)))/(3r1g^2r2g^2r3g^2))+(2ellF((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g*Yint*r3g^2)-(2(-r1g+r2g)(r1g-r3g)^2Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3r1g^2(r2g(r1g-r3g))^(3/2)r3g)))+(((rg-r2g)(rg-r3g))/((r1g-r2g)(r1g-r3g))\[Delta]r1+((r1g-rg)(rg-r3g))/((r1g-r2g)(r2g-r3g))\[Delta]r2)
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunFFPerParRes[wr_,a_,p_,e_,x_,{\[Delta]r1_,\[Delta]r2_},{\[CapitalUpsilon]rg_,\[Delta]\[CapitalUpsilon]r_}]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,EEg,\[Delta]EE,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE,rg},
	EEg=EEgfun[a,p,e,x];
	(*\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];*)
	(*\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunFF[a,p,e,x];*)

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	(*{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];*)
	\[Delta]r3=\[Delta]r3funFF[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
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
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobiSpectral`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobiSpectral`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobiSpectral`Private`*"],True]]]];
	
	\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg(2ellF)/Yint+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Yint+((r2g*ellE-r1g*ellF)\[Delta]r1)/(Yint*r1g(r1g-r2g))+1/Yint(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g))\[Delta]r2+1/Yint 1/r3g(((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]))\[Delta]r3+2((-(r1g-r3g)ellE+r1g*ellF)/(r1g*r3g*Yint))\[Delta]r41-a/2(-((4ellE*Yint(r2g*r3g+r1g(r2g+r3g)))/(3r1g^2r2g^2r3g^2))+(2ellF((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g*r3g^2*Yint)-(2(-r1g+r2g)(r1g-r3g)Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3r1g^2r2g*r3g*Yint))+Sign[\[Pi]-Mod[wr,2\[Pi]]](Sqrt[(-r2g+rg)(-r3g+rg)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-rg)rg])\[Delta]r1+Sqrt[(r1g-rg)(-r3g+rg)]/((r1g-r2g)(r2g-r3g)Sqrt[rg(-r2g+rg)])\[Delta]r2)
]


(* ::Subsection::Closed:: *)
(*Velocity*)


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunFFPerPar[wr_,a_,p_,e_,x_,{\[Delta]r1_,\[Delta]r2_},{\[CapitalUpsilon]rg_,\[Delta]\[CapitalUpsilon]r_}]:=Module[{EEg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,dRgdr,rg},
	EEg=EEgfun[a,p,e,x];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	(*{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];*)
	\[Delta]r3=\[Delta]r3funFF[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;

	dRgdr=Sqrt[(1-EEg^2)](-4rg^3+r1g*r2g*r3g+3rg^2(r1g+r2g+r3g)-2rg(r2g*r3g+r1g(r2g+r3g)));
	rg=rgICr1gfun[wr,a,p,e,x];

	-Sign[\[Pi]-Mod[wr,2\[Pi]]]Sqrt[(1-EEg^2)rg(r1g-rg)(rg-r2g)(rg-r3g)](\[Delta]r1/(2(rg-r1g))+\[Delta]r2/(2(rg-r2g))+\[Delta]r3/(2(rg-r3g))+\[Delta]r41/rg-a/(2rg^2)+(EEg*\[Delta]EE)/(1-EEg^2))+1/2dRgdr*\[Delta]rfunFFPerParRes[wr,a,p,e,x,{\[Delta]r1,\[Delta]r2},{\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r}]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity*)


\[Delta]vtfunFFPerPar[wr_,a_,p_,e_,x_,{\[Delta]r1_,\[Delta]r2_},{\[CapitalUpsilon]rg_,\[Delta]\[CapitalUpsilon]r_}]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,rg,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	(*{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];*)
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzgdr1gfun[a,p,e,x]\[Delta]r1+dLzgdr2gfun[a,p,e,x]\[Delta]r2;

	rg=rgICr1gfun[wr,a,p,e,x];
	\[CapitalDelta]=a^2-2rg+rg^2;
	d\[Delta]td\[Lambda]fun=(-Lzg(a^2+rg^2)+a*EEg(a^2+3rg^2))/(rg*\[CapitalDelta])+(rg(2a^2+a^2*rg+rg^3)\[Delta]EE)/\[CapitalDelta]-(2a*rg*\[Delta]Lz)/\[CapitalDelta];
	ddtgd\[Lambda]drfun=(EEg*rg(a^2+3rg^2))/\[CapitalDelta]-((rg^2-a^2)(-2a*Lzg+EEg*rg^3+a^2*EEg(2+rg)))/\[CapitalDelta]^2;

	ddtgd\[Lambda]drfun*\[Delta]rfunFFPerPar[wr,a,p,e,x,{\[Delta]r1,\[Delta]r2},{\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r}]+d\[Delta]td\[Lambda]fun
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


\[Delta]v\[Phi]funFFPerPar[wr_,a_,p_,e_,x_,{\[Delta]r1_,\[Delta]r2_},{\[CapitalUpsilon]rg_,\[Delta]\[CapitalUpsilon]r_}]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,rg,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	(*{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];*)
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzgdr1gfun[a,p,e,x]\[Delta]r1+dLzgdr2gfun[a,p,e,x]\[Delta]r2;

	rg=rgICr1gfun[wr,a,p,e,x];
	\[CapitalDelta]=a^2-2rg+rg^2;
	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+rg)))/(rg*\[CapitalDelta])+(2a*rg*\[Delta]EE)/\[CapitalDelta]+((rg-2)rg*\[Delta]Lz)/\[CapitalDelta];
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(rg-1)-EEg*rg^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunFFPerPar[wr,a,p,e,x,{\[Delta]r1,\[Delta]r2},{\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r}]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Section::Closed:: *)
(*Spin corrections to the orbit - parallel component of the spin, fixed eccentricity*)


(* ::Subsection::Closed:: *)
(*Trajectory*)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunFEPerPar[wr_,a_,p_,e_,x_,{\[CapitalUpsilon]rg_,\[Delta]\[CapitalUpsilon]r_}]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,\[Delta]r1,\[Delta]r2,EEg,\[Delta]EE,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE,rg},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	(*\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];*)
	(*\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunFE[a,p,e,x];*)

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
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobiSpectral`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobiSpectral`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobiSpectral`Private`*"],True]]]];

	((jSN*jCN*r1g(r1g-r2g)r2g*Sqrt[r1g(r2g-jSN^2*r3g)-(1-jSN^2)r2g*r3g])/(jSN^2(r1g-r2g)+r2g)^2)(\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg*(2ellF)/Yint+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Yint+((r2g*ellE-r1g*ellF)\[Delta]r1)/(Yint*r1g(r1g-r2g))+1/Yint(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g))\[Delta]r2+1/Yint 1/r3g(((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]))\[Delta]r3+2((-(r1g-r3g)ellE+r1g*ellF)/(r1g*Yint*r3g))\[Delta]r41-a/2 (-((4ellE Sqrt[r2g(r1g-r3g)](r2g*r3g+r1g(r2g+r3g)))/(3r1g^2r2g^2r3g^2))+(2ellF((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g Sqrt[r2g (r1g-r3g)]r3g^2)-(2(-r1g+r2g)(r1g-r3g)^2Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3r1g^2(r2g(r1g-r3g))^(3/2)r3g)))+(((rg-r2g)(rg-r3g))/((r1g-r2g)(r1g-r3g))\[Delta]r1+((r1g-rg)(rg-r3g))/((r1g-r2g)(r2g-r3g))\[Delta]r2)
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunFEPerParRes[wr_,a_,p_,e_,x_,{\[CapitalUpsilon]rg_,\[Delta]\[CapitalUpsilon]r_}]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,\[Delta]r1,\[Delta]r2,EEg,\[Delta]EE,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE,rg},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	(*\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];*)
	(*\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunFE[a,p,e,x];*)

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
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobiSpectral`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobiSpectral`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobiSpectral`Private`*"],True]]]];
	
	\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg(2ellF)/Yint+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Yint+((r2g*ellE-r1g*ellF)\[Delta]r1)/(Yint*r1g(r1g-r2g))+1/Yint(-(((r1g-r3g)ellE)/((r1g-r2g)(r2g-r3g)))+ellF/(r1g-r2g))\[Delta]r2+1/Yint 1/r3g (((r2g*ellE)/(r2g-r3g)-ellF)+(r3g*r2g(-r1g+r2g)Cos[\[Phi]]*Sin[\[Phi]])/((r2g(r1g-r3g))(r2g-r3g)Sqrt[1-krg*Sin[\[Phi]]^2]))\[Delta]r3+2((-(r1g-r3g)ellE+r1g*ellF)/(r1g*r3g*Sqrt[r2g(r1g-r3g)]))\[Delta]r41-a/2(-((4ellE*Yint(r2g*r3g+r1g(r2g+r3g)))/(3r1g^2r2g^2r3g^2))+(2ellF((r2g-r3g)r3g+r1g(2r2g+r3g)))/(3r1g*r2g*r3g^2*Yint)-(2(-r1g+r2g)(r1g-r3g)Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(3r1g^2r2g*r3g Sqrt[r2g(r1g-r3g)]))+Sign[\[Pi]-Mod[wr,2\[Pi]]](Sqrt[(-r2g+rg)(-r3g+rg)]/((r1g-r2g)(r1g-r3g)Sqrt[(r1g-rg)rg])\[Delta]r1+Sqrt[(r1g-rg)(-r3g+rg)]/((r1g-r2g)(r2g-r3g)Sqrt[rg(-r2g+rg)])\[Delta]r2)
]


(* ::Subsection::Closed:: *)
(*Velocity*)


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunFEPerPar[wr_,a_,p_,e_,x_,{\[CapitalUpsilon]rg_,\[Delta]\[CapitalUpsilon]r_}]:=Module[{EEg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r1,\[Delta]r2,\[Delta]r3,\[Delta]r41,dRgdr,rg},
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

	-Sign[\[Pi]-Mod[wr,2\[Pi]]]Sqrt[(1-EEg^2)rg(r1g-rg)(rg-r2g)(rg-r3g)](\[Delta]r1/(2(rg-r1g))+\[Delta]r2/(2(rg-r2g))+\[Delta]r3/(2(rg-r3g))+\[Delta]r41/rg-a/(2rg^2)+(EEg*\[Delta]EE)/(1-EEg^2))+1/2dRgdr*\[Delta]rfunFEPerParRes[wr,a,p,e,x,{\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r}]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity*)


\[Delta]vtfunFEPerPar[wr_,a_,p_,e_,x_,{\[CapitalUpsilon]rg_,\[Delta]\[CapitalUpsilon]r_}]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,r3g,rg,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
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

	ddtgd\[Lambda]drfun*\[Delta]rfunFEPerPar[wr,a,p,e,x,{\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r}]+d\[Delta]td\[Lambda]fun
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


\[Delta]v\[Phi]funFEPerPar[wr_,a_,p_,e_,x_,{\[CapitalUpsilon]rg_,\[Delta]\[CapitalUpsilon]r_}]:=Module[{EEg,Lzg,\[Delta]EE,\[Delta]Lz,rg,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];

	rg=rgICr1gfun[wr,a,p,e,x];
	\[CapitalDelta]=a^2-2rg+rg^2;
	d\[Delta]\[Phi]d\[Lambda]fun=-EEg+(a(-Lzg+a*EEg(1+rg)))/(rg*\[CapitalDelta])+(2a*rg*\[Delta]EE)/\[CapitalDelta]+((rg-2)rg*\[Delta]Lz)/\[CapitalDelta];
	dd\[Phi]gd\[Lambda]drfun=(2a(a^2*EEg+a*Lzg(rg-1)-EEg*rg^2))/\[CapitalDelta]^2;

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunFEPerPar[wr,a,p,e,x,{\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r}]+d\[Delta]\[Phi]d\[Lambda]fun
]


(* ::Section::Closed:: *)
(*Spin corrections to the orbit - parallel component of the spin, fixed turning points*)


(* ::Subsection::Closed:: *)
(*Trajectory*)


(* ::Subsubsection::Closed:: *)
(*Radial trajectory*)


\[Delta]rfunFTPerPar[wr_,a_,p_,e_,x_,{\[CapitalUpsilon]rg_,\[Delta]\[CapitalUpsilon]r_}]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,EEg,\[Delta]EE,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	(*\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];*)
	(*\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunFT[a,p,e,x];*)

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
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobiSpectral`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobiSpectral`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobiSpectral`Private`*"],True]]]];

	-((jSN*jCN(r1g-r2g)(r1g-r3g)(r2g-r3g)Sqrt[r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g])/((1-jSN^2)r1g+jSN^2*r2g-r3g)^2)(\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg(2ellF)/Yint+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Yint+\[Delta]r3/2 1/r3g  2/Yint((r2g*ellE)/(r2g-r3g)-ellF)+2((jSN*jCN(r1g-r2g))/(r1g*r2g*Sqrt[r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g])-((r1g-r3g)ellE-r1g*ellF)/(r1g*r3g*Yint))\[Delta]r41-a/2 (-(2/(3r1g*r2g*r3g))((jSN*jCN(r1g-r2g)(r1g-r3g)(r2g-r3g)Sqrt[(r1g*r2g-jSN^2*r1g*r3g-r2g*r3g*jCN^2)])/(r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g)^2+((r1g+r2g+r3g)ellF)/Yint-2(r2g*r3g+r1g(r2g+r3g))((jSN*jCN(r1g-r2g))/(r1g*r2g*Sqrt[r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g])+(-((r1g-r3g)ellE)+r1g*ellF)/(r1g*r3g*Yint)))))
]


(* ::Subsubsection::Closed:: *)
(*Radial trajectory - rescaled version*)


\[Delta]rfunFTPerParRes[wr_,a_,p_,e_,x_,{\[CapitalUpsilon]rg_,\[Delta]\[CapitalUpsilon]r_}]:=Module[{r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,EEg,\[Delta]EE,\[Delta]Lz,\[Phi],jSN,jCN,krghold,krg,\[Gamma]rg,Yint,ellF,ellE},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	(*\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];*)
	(*\[Delta]\[CapitalUpsilon]r=\[Delta]\[CapitalUpsilon]rfunFT[a,p,e,x];*)

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
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobiSpectral`Private`*"],StringMatchQ[#,"AnalyticSpinOrbitsHamiltonJacobiSpectral`Private`krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobiSpectral`Private`*"],True]]]];

	-(\[Delta]\[CapitalUpsilon]r/\[CapitalUpsilon]rg(2ellF)/Yint+(EEg*\[Delta]EE)/(1-EEg^2)(2ellF)/Yint+\[Delta]r3/2 1/r3g  2/Yint((r2g*ellE)/(r2g-r3g)-ellF)+2((jSN*jCN(r1g-r2g))/(r1g*r2g*Sqrt[r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g])-((r1g-r3g)ellE-r1g*ellF)/(r1g*r3g*Yint))\[Delta]r41-a/2(-(2/(3r1g*r2g*r3g))((jSN*jCN(r1g-r2g)(r1g-r3g)(r2g-r3g)Sqrt[(r1g*r2g-jSN^2*r1g*r3g-r2g*r3g*jCN^2)])/(r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g)^2+((r1g+r2g+r3g)ellF)/Sqrt[r2g(r1g-r3g)]-2(r2g*r3g+r1g(r2g+r3g))((jSN*jCN(r1g-r2g))/(r1g*r2g*Sqrt[r1g(r2g-jSN^2*r3g)-jCN^2*r2g*r3g])+(-((r1g-r3g)ellE)+r1g*ellF)/(r1g*r3g*Yint)))))
]


(* ::Subsection::Closed:: *)
(*Velocity*)


(* ::Subsubsection::Closed:: *)
(*Radial velocity*)


\[Delta]vrfunFTPerPar[wr_,a_,p_,e_,x_,{\[CapitalUpsilon]rg_,\[Delta]\[CapitalUpsilon]r_}]:=Module[{EEg,\[Delta]EE,r1g,r2g,r3g,\[Delta]r3,\[Delta]r41,dRgdr,rg},
	EEg=EEgfun[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	\[Delta]r3=\[Delta]r3funFT[a,p,e,x];
	\[Delta]r41=\[Delta]r41fun[a,p,e,x];

	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);
	
	dRgdr=Sqrt[(1-EEg^2)](-4rg^3+r1g*r2g*r3g+3rg^2(r1g+r2g+r3g)-2rg(r2g*r3g+r1g(r2g+r3g)));
	rg=rgICr2gfun[wr,a,p,e,x];

	-Sign[\[Pi]-Mod[wr,2\[Pi]]]Sqrt[(1-EEg^2)rg(r1g-rg)(rg-r2g)(rg-r3g)](\[Delta]r3/(2(rg-r3g))+\[Delta]r41/rg-a/(2rg^2)+(EEg*\[Delta]EE)/(1-EEg^2))+1/2dRgdr*\[Delta]rfunFTPerParRes[wr,a,p,e,x,{\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r}]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity*)


\[Delta]vtfunFTPerPar[wr_,a_,p_,e_,x_,{\[CapitalUpsilon]rg_,\[Delta]\[CapitalUpsilon]r_}]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,r3g,rg,\[CapitalDelta],d\[Delta]td\[Lambda]fun,ddtgd\[Lambda]drfun},
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

	ddtgd\[Lambda]drfun*\[Delta]rfunFTPerPar[wr,a,p,e,x,{\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r}]+d\[Delta]td\[Lambda]fun
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity*)


\[Delta]v\[Phi]funFTPerPar[wr_,a_,p_,e_,x_,{\[CapitalUpsilon]rg_,\[Delta]\[CapitalUpsilon]r_}]:=Module[{rp,rm,EEg,Lzg,\[Delta]EE,\[Delta]Lz,r1g,r2g,r3g,rg,\[CapitalDelta],d\[Delta]\[Phi]d\[Lambda]fun,dd\[Phi]gd\[Lambda]drfun},
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

	dd\[Phi]gd\[Lambda]drfun*\[Delta]rfunFTPerPar[wr,a,p,e,x,{\[CapitalUpsilon]rg,\[Delta]\[CapitalUpsilon]r}]+d\[Delta]\[Phi]d\[Lambda]fun
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
	
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobiSpectral`Private`*"],StringMatchQ[#,"krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobiSpectral`Private`*"],True]]]];
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
	
	Remove[Evaluate[ToExpression[Pick[Names["AnalyticSpinOrbitsHamiltonJacobiSpectral`Private`*"],StringMatchQ[#,"krghold$"~~__]&/@Names["AnalyticSpinOrbitsHamiltonJacobiSpectral`Private`*"],True]]]];
	Re[\[Psi]freq]-\[Psi]pref[a,p,e,x]wr
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


(* ::Section::Closed:: *)
(*Analytical functions for geodesic frequencies and spin-corrections*)


(* ::Subsection::Closed:: *)
(*Fixed frequency*)


KerrNearEqFrequencyCorrFF[a_, p_, e_, x_]:=Module[{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g,\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]\[Phi]s},
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];
	
	\[CapitalUpsilon]rs=\[Delta]\[CapitalUpsilon]rfunFF[a,p,e,x];
					
	<|
	(* "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},*)
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 "MinoFrequenciesCorrection"->\[CapitalUpsilon]rs
	 (*"MinoPrecessionFrequency"->\[CapitalUpsilon]p,
	 "BLPrecessionFrequency"->\[CapitalUpsilon]p/\[CapitalUpsilon]tg,*)
	|>
]


(* ::Subsection::Closed:: *)
(*Fixed eccentricity*)


KerrNearEqFrequencyCorrFE[a_, p_, e_, x_]:=Module[{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g,\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]\[Phi]s},
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];
	
	\[CapitalUpsilon]ts=\[Delta]\[CapitalUpsilon]tfunFE[a,p,e,x];
	\[CapitalUpsilon]rs=\[Delta]\[CapitalUpsilon]rfunFE[a,p,e,x];
	\[CapitalUpsilon]\[Phi]s=\[Delta]\[CapitalUpsilon]\[Phi]funFE[a,p,e,x];
					
	<|
	(* "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},*)
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 "MinoFrequenciesCorrection"->{\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]\[Phi]s},
	 "BLFrequenciesCorrection"->{\[CapitalUpsilon]rs/\[CapitalUpsilon]tg-\[CapitalUpsilon]rg/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts,\[CapitalUpsilon]\[Phi]s/\[CapitalUpsilon]tg-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts}
	 (*"MinoPrecessionFrequency"->\[CapitalUpsilon]p,
	 "BLPrecessionFrequency"->\[CapitalUpsilon]p/\[CapitalUpsilon]tg,*)
	|>
]


(* ::Subsection::Closed:: *)
(*Fixed turning points*)


KerrNearEqFrequencyCorrFT[a_, p_, e_, x_]:=Module[{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g,\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]\[Phi]s},
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];
	
	\[CapitalUpsilon]ts=\[Delta]\[CapitalUpsilon]tfunFT[a,p,e,x];
	\[CapitalUpsilon]rs=\[Delta]\[CapitalUpsilon]rfunFT[a,p,e,x];
	\[CapitalUpsilon]\[Phi]s=\[Delta]\[CapitalUpsilon]\[Phi]funFT[a,p,e,x];
					
	<|
	(* "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},*)
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 "MinoFrequenciesCorrection"->{\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]\[Phi]s},
	 "BLFrequenciesCorrection"->{\[CapitalUpsilon]rs/\[CapitalUpsilon]tg-\[CapitalUpsilon]rg/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts,\[CapitalUpsilon]\[Phi]s/\[CapitalUpsilon]tg-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts}
	 (*"MinoPrecessionFrequency"->\[CapitalUpsilon]p,
	 "BLPrecessionFrequency"->\[CapitalUpsilon]p/\[CapitalUpsilon]tg,*)
	|>
]


(* ::Section::Closed:: *)
(*Near equatorial orbits - fixed frequency*)


(* ::Subsection::Closed:: *)
(*Fourier series expansion*)


KerrNearEqSpinOrbitCorrFFSpectral[a_, p_, e_, x_, prec_]:=
Module[{precsol,EEg,Lzg,\[Delta]r1,\[Delta]r2,\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g,\[CapitalUpsilon]p,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]zs,\[CapitalUpsilon]\[Phi]s,
	\[Psi]phase,\[Delta]zort,dtrgd\[Lambda]coeff,\[CapitalDelta]trgcoeff,\[CapitalDelta]\[Phi]rgcoeff,d\[Phi]rgd\[Lambda]coeff,growthraterg,rgcoeff,dtspard\[Lambda]coeff,\[Delta]rsparcoeff,\[CapitalDelta]tsparcoeff,d\[Phi]spard\[Lambda]coeff,\[CapitalDelta]\[Phi]sparcoeff,growthrate\[Delta]r},
	
	precsol=Min[Precision[{a,p,e,x}],prec];
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	(* spin correction constants of motion *)
	{\[Delta]r1,\[Delta]r2}=\[Delta]turnpointsfunFF[a,p,e,x];
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEgdr1gfun[a,p,e,x]\[Delta]r1+dEgdr2gfun[a,p,e,x]\[Delta]r2;
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzgdr1gfun[a,p,e,x]\[Delta]r1+dLzgdr2gfun[a,p,e,x]\[Delta]r2;
	(* geodesic radial frequency *)
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	(*\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,x];*)
	(* spin correction radial frequency *)
	\[CapitalUpsilon]rs=\[Delta]\[CapitalUpsilon]rfunFF[a,p,e,x];
	(*\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,x];*)
	(*\[Psi]phase[wr_]:=\[Psi]pICr1g[wr,a,p,e,x];
	\[Delta]zort[wr_]:=\[Delta]zfunICr1gPer[wr,a,p,e,x];*)
	
	(* Fourier coefficients geodesic functions *)
	{\[CapitalUpsilon]tg,dtrgd\[Lambda]coeff,\[CapitalDelta]trgcoeff}=CoeffsFourier[VtrgICr1gfun[#,a,p,e,x]&,\[CapitalUpsilon]rg,precsol];
	{\[CapitalUpsilon]\[Phi]g,d\[Phi]rgd\[Lambda]coeff,\[CapitalDelta]\[Phi]rgcoeff}=CoeffsFourier[V\[Phi]rgICr1gfun[#,a,p,e,x]&,\[CapitalUpsilon]rg,precsol];
	{growthraterg,rgcoeff}=CoeffsFourier[rgICr1gfun[#,a,p,e,x]&,\[CapitalUpsilon]rg,precsol][[1;;2]];
	
	(* Fourier coefficients spin-corrections coordinate-time and azimuthal velocities *)
	{\[CapitalUpsilon]ts,dtspard\[Lambda]coeff,\[CapitalDelta]tsparcoeff}=CoeffsFourier[\[Delta]vtfunFFPerPar[#,a,p,e,x,{\[Delta]r1,\[Delta]r2},{\[CapitalUpsilon]rg,\[CapitalUpsilon]rs}]&,\[CapitalUpsilon]rg,precsol];
	{\[CapitalUpsilon]\[Phi]s,d\[Phi]spard\[Lambda]coeff,\[CapitalDelta]\[Phi]sparcoeff}=CoeffsFourier[\[Delta]v\[Phi]funFFPerPar[#,a,p,e,x,{\[Delta]r1,\[Delta]r2},{\[CapitalUpsilon]rg,\[CapitalUpsilon]rs}]&,\[CapitalUpsilon]rg,precsol];
	{growthrate\[Delta]r,\[Delta]rsparcoeff}=CoeffsFourier[\[Delta]rfunFFPerPar[#,a,p,e,x,{\[Delta]r1,\[Delta]r2},{\[CapitalUpsilon]rg,\[CapitalUpsilon]rs}]&,\[CapitalUpsilon]rg,precsol][[1;;2]];
	
	(*The RealSign serves to match the sign convetion of generic orbits solutions in the limit z->0*)		
	<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->RealSign[x]*\[Delta]EE,
	 "Jzs"->RealSign[x]*\[Delta]Lz,
	 "\[Delta]r1"->RealSign[x]*\[Delta]r1,
	 "\[Delta]r2"->RealSign[x]*\[Delta]r2,
	(* "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},*)
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 "MinoFrequenciesCorrection"->RealSign[x]*{\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]\[Phi]s},
	 "BLFrequenciesCorrection"->RealSign[x]*{\[CapitalUpsilon]rs/\[CapitalUpsilon]tg-\[CapitalUpsilon]rg/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts,\[CapitalUpsilon]\[Phi]s/\[CapitalUpsilon]tg-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts},
	 (*"MinoPrecessionFrequency"->\[CapitalUpsilon]p,
	 "BLPrecessionFrequency"->\[CapitalUpsilon]p/\[CapitalUpsilon]tg,*)
	 (*Keys purely oscillatory part geodesic coordinate time and azimuthal trajectory*)
	 "\[CapitalDelta]trg"->Function[{wr},Evaluate[\[CapitalDelta]integratedFunc[wr,\[CapitalDelta]trgcoeff]]],
	 "\[CapitalDelta]\[Phi]rg"->Function[{wr},Evaluate[\[CapitalDelta]integratedFunc[wr,\[CapitalDelta]\[Phi]rgcoeff]]],
	 (*Keys geodesic radial trajectory*)
	 "rg"->Function[{wr},rgICr1gfun[wr,a,p,e,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{wr},Evaluate[fourierCos[wr,\[CapitalUpsilon]tg,dtrgd\[Lambda]coeff]]],
	 "vrg"->Function[{wr},drgd\[Lambda]fun[wr,a,p,e,x]],
	 "v\[Phi]g"->Function[{wr},Evaluate[fourierCos[wr,\[CapitalUpsilon]\[Phi]g,d\[Phi]rgd\[Lambda]coeff]]],
	(*Keys corrections trajectory*)
	 "\[CapitalDelta]\[Delta]tpar"->Function[{wr},Evaluate[RealSign[x]*\[CapitalDelta]\[Delta]integratedFunc[wr,\[CapitalUpsilon]rg,\[CapitalUpsilon]rs,\[CapitalDelta]trgcoeff,\[CapitalDelta]tsparcoeff]]],
	 "\[Delta]rpar"->Function[{wr},Evaluate[RealSign[x]*fourierCos[wr,growthrate\[Delta]r,\[Delta]rsparcoeff]]],
	 (*"\[Psi]p"->Function[wr,\[Psi]phase[wr]],*)
	(* "\[Delta]zort"->Function[{wr},\[Delta]zort[wr]],*)
	 "\[CapitalDelta]\[Delta]\[Phi]par"->Function[{wr},Evaluate[RealSign[x]*\[CapitalDelta]\[Delta]integratedFunc[wr,\[CapitalUpsilon]rg,\[CapitalUpsilon]rs,\[CapitalDelta]\[Phi]rgcoeff,\[CapitalDelta]\[Phi]sparcoeff]]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{wr},Evaluate[RealSign[x]*fourierCos[wr,\[CapitalUpsilon]ts,dtspard\[Lambda]coeff]]],
	 "\[Delta]vrpar"->Function[{wr},Evaluate[RealSign[x]*fourierd\[Delta]rd\[Lambda][wr,\[CapitalUpsilon]rg,\[CapitalUpsilon]rs,rgcoeff,\[Delta]rsparcoeff]]],
	 "\[Delta]v\[Phi]par"->Function[{wr},Evaluate[RealSign[x]*fourierCos[wr,\[CapitalUpsilon]\[Phi]s,d\[Phi]spard\[Lambda]coeff]]]
	|>
]


(* ::Section::Closed:: *)
(*Near equatorial orbits - fixed eccentricity*)


(* ::Subsection::Closed:: *)
(*Fourier series expansion*)


KerrNearEqSpinOrbitCorrFESpectral[a_, p_, e_, x_, prec_]:=
Module[{precsol,EEg,Lzg,\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g,\[CapitalUpsilon]p,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]zs,\[CapitalUpsilon]\[Phi]s,
	\[Psi]phase,\[Delta]zort,dtrgd\[Lambda]coeff,\[CapitalDelta]trgcoeff,\[CapitalDelta]\[Phi]rgcoeff,d\[Phi]rgd\[Lambda]coeff,growthraterg,rgcoeff,dtspard\[Lambda]coeff,\[Delta]rsparcoeff,\[CapitalDelta]tsparcoeff,d\[Phi]spard\[Lambda]coeff,\[CapitalDelta]\[Phi]sparcoeff,growthrate\[Delta]r},
	
	precsol=Min[Precision[{a,p,e,x}],prec];
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	(* spin correction constants of motion *)
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x]+dEEdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x]+dLzdpfun[a,p,e,x]\[Delta]pfun[a,p,e,x];
	(* geodesic radial frequency *)
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	(*\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,x];*)
	(* spin correction radial frequency *)
	\[CapitalUpsilon]rs=\[Delta]\[CapitalUpsilon]rfunFE[a,p,e,x];
	(*\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,x];*)
	(*\[Psi]phase[wr_]:=\[Psi]pICr1g[wr,a,p,e,x];
	\[Delta]zort[wr_]:=\[Delta]zfunICr1gPer[wr,a,p,e,x];*)
	
	(* Fourier coefficients geodesic functions *)
	{\[CapitalUpsilon]tg,dtrgd\[Lambda]coeff,\[CapitalDelta]trgcoeff}=CoeffsFourier[VtrgICr1gfun[#,a,p,e,x]&,\[CapitalUpsilon]rg,precsol];
	{\[CapitalUpsilon]\[Phi]g,d\[Phi]rgd\[Lambda]coeff,\[CapitalDelta]\[Phi]rgcoeff}=CoeffsFourier[V\[Phi]rgICr1gfun[#,a,p,e,x]&,\[CapitalUpsilon]rg,precsol];
	{growthraterg,rgcoeff}=CoeffsFourier[rgICr1gfun[#,a,p,e,x]&,\[CapitalUpsilon]rg,precsol][[1;;2]];
	
	(* Fourier coefficients spin-corrections coordinate-time and azimuthal velocities *)
	{\[CapitalUpsilon]ts,dtspard\[Lambda]coeff,\[CapitalDelta]tsparcoeff}=CoeffsFourier[\[Delta]vtfunFEPerPar[#,a,p,e,x,{\[CapitalUpsilon]rg,\[CapitalUpsilon]rs}]&,\[CapitalUpsilon]rg,precsol];
	{\[CapitalUpsilon]\[Phi]s,d\[Phi]spard\[Lambda]coeff,\[CapitalDelta]\[Phi]sparcoeff}=CoeffsFourier[\[Delta]v\[Phi]funFEPerPar[#,a,p,e,x,{\[CapitalUpsilon]rg,\[CapitalUpsilon]rs}]&,\[CapitalUpsilon]rg,precsol];
	{growthrate\[Delta]r,\[Delta]rsparcoeff}=CoeffsFourier[\[Delta]rfunFEPerPar[#,a,p,e,x,{\[CapitalUpsilon]rg,\[CapitalUpsilon]rs}]&,\[CapitalUpsilon]rg,precsol][[1;;2]];
	
	(*The RealSign serves to match the sign convetion of generic orbits solutions in the limit z->0*)				
	<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->RealSign[x]*\[Delta]EE,
	 "Jzs"->RealSign[x]*\[Delta]Lz,
	 "\[Delta]p"->RealSign[x]*\[Delta]pfun[a,p,e,x],
	(* "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},*)
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 "MinoFrequenciesCorrection"->RealSign[x]*{\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]\[Phi]s},
	 "BLFrequenciesCorrection"->RealSign[x]*{\[CapitalUpsilon]rs/\[CapitalUpsilon]tg-\[CapitalUpsilon]rg/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts,\[CapitalUpsilon]\[Phi]s/\[CapitalUpsilon]tg-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts},
	 (*"MinoPrecessionFrequency"->\[CapitalUpsilon]p,
	 "BLPrecessionFrequency"->\[CapitalUpsilon]p/\[CapitalUpsilon]tg,*)
	 (*Keys purely oscillatory part geodesic coordinate time and azimuthal trajectory*)
	 "\[CapitalDelta]trg"->Function[{wr},Evaluate[\[CapitalDelta]integratedFunc[wr,\[CapitalDelta]trgcoeff]]],
	 "\[CapitalDelta]\[Phi]rg"->Function[{wr},Evaluate[\[CapitalDelta]integratedFunc[wr,\[CapitalDelta]\[Phi]rgcoeff]]],
	 (*Keys geodesic radial trajectory*)
	 "rg"->Function[{wr},rgICr1gfun[wr,a,p,e,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{wr},Evaluate[fourierCos[wr,\[CapitalUpsilon]tg,dtrgd\[Lambda]coeff]]],
	 "vrg"->Function[{wr},drgd\[Lambda]fun[wr,a,p,e,x]],
	 "v\[Phi]g"->Function[{wr},Evaluate[fourierCos[wr,\[CapitalUpsilon]\[Phi]g,d\[Phi]rgd\[Lambda]coeff]]],
	(*Keys corrections trajectory*)
	 "\[CapitalDelta]\[Delta]tpar"->Function[{wr},Evaluate[RealSign[x]*\[CapitalDelta]\[Delta]integratedFunc[wr,\[CapitalUpsilon]rg,\[CapitalUpsilon]rs,\[CapitalDelta]trgcoeff,\[CapitalDelta]tsparcoeff]]],
	 "\[Delta]rpar"->Function[{wr},Evaluate[RealSign[x]*fourierCos[wr,growthrate\[Delta]r,\[Delta]rsparcoeff]]],
	 (*"\[Psi]p"->Function[wr,\[Psi]phase[wr]],*)
	(* "\[Delta]zort"->Function[{wr},\[Delta]zort[wr]],*)
	 "\[CapitalDelta]\[Delta]\[Phi]par"->Function[{wr},Evaluate[RealSign[x]*\[CapitalDelta]\[Delta]integratedFunc[wr,\[CapitalUpsilon]rg,\[CapitalUpsilon]rs,\[CapitalDelta]\[Phi]rgcoeff,\[CapitalDelta]\[Phi]sparcoeff]]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{wr},Evaluate[RealSign[x]*fourierCos[wr,\[CapitalUpsilon]ts,dtspard\[Lambda]coeff]]],
	 "\[Delta]vrpar"->Function[{wr},Evaluate[RealSign[x]*fourierd\[Delta]rd\[Lambda][wr,\[CapitalUpsilon]rg,\[CapitalUpsilon]rs,rgcoeff,\[Delta]rsparcoeff]]],
	 "\[Delta]v\[Phi]par"->Function[{wr},Evaluate[RealSign[x]*fourierCos[wr,\[CapitalUpsilon]\[Phi]s,d\[Phi]spard\[Lambda]coeff]]]
	|>
]


(* ::Section::Closed:: *)
(*Near equatorial orbits - fixed turning points*)


(* ::Subsection::Closed:: *)
(*Fourier series expansion*)


KerrNearEqSpinOrbitCorrFTSpectral[a_, p_, e_, x_, prec_]:=
Module[{precsol,EEg,Lzg,\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g,\[CapitalUpsilon]p,\[Delta]EE,\[Delta]Lz,\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]zs,\[CapitalUpsilon]\[Phi]s,
	\[Psi]phase,\[Delta]zort,dtrgd\[Lambda]coeff,\[CapitalDelta]trgcoeff,\[CapitalDelta]\[Phi]rgcoeff,d\[Phi]rgd\[Lambda]coeff,growthraterg,rgcoeff,dtspard\[Lambda]coeff,\[Delta]rsparcoeff,\[CapitalDelta]tsparcoeff,d\[Phi]spard\[Lambda]coeff,\[CapitalDelta]\[Phi]sparcoeff,growthrate\[Delta]r},
	
	precsol=Min[Precision[{a,p,e,x}],prec];
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	(* spin correction constants of motion *)
	\[Delta]EE=\[Delta]EEfunFT[a,p,e,x];
	\[Delta]Lz=\[Delta]LzfunFT[a,p,e,x];
	(* geodesic radial frequency *)
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	(*\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,x];*)
	(* spin correction radial frequency *)
	\[CapitalUpsilon]rs=\[Delta]\[CapitalUpsilon]rfunFT[a,p,e,x];
	(*\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,x];*)
	(*\[Psi]phase[wr_]:=\[Psi]pICr1g[wr,a,p,e,x];
	\[Delta]zort[wr_]:=\[Delta]zfunICr1gPer[wr,a,p,e,x];*)
	
	(* Fourier coefficients geodesic functions *)
	{\[CapitalUpsilon]tg,dtrgd\[Lambda]coeff,\[CapitalDelta]trgcoeff}=CoeffsFourier[VtrgICr2gfun[#,a,p,e,x]&,\[CapitalUpsilon]rg,precsol];
	{\[CapitalUpsilon]\[Phi]g,d\[Phi]rgd\[Lambda]coeff,\[CapitalDelta]\[Phi]rgcoeff}=CoeffsFourier[V\[Phi]rgICr2gfun[#,a,p,e,x]&,\[CapitalUpsilon]rg,precsol];
	{growthraterg,rgcoeff}=CoeffsFourier[rgICr2gfun[#,a,p,e,x]&,\[CapitalUpsilon]rg,precsol][[1;;2]];
	
	(* Fourier coefficients spin-corrections coordinate-time and azimuthal velocities *)
	{\[CapitalUpsilon]ts,dtspard\[Lambda]coeff,\[CapitalDelta]tsparcoeff}=CoeffsFourier[\[Delta]vtfunFTPerPar[#,a,p,e,x,{\[CapitalUpsilon]rg,\[CapitalUpsilon]rs}]&,\[CapitalUpsilon]rg,precsol];
	{\[CapitalUpsilon]\[Phi]s,d\[Phi]spard\[Lambda]coeff,\[CapitalDelta]\[Phi]sparcoeff}=CoeffsFourier[\[Delta]v\[Phi]funFTPerPar[#,a,p,e,x,{\[CapitalUpsilon]rg,\[CapitalUpsilon]rs}]&,\[CapitalUpsilon]rg,precsol];
	{growthrate\[Delta]r,\[Delta]rsparcoeff}=CoeffsFourier[\[Delta]rfunFTPerPar[#,a,p,e,x,{\[CapitalUpsilon]rg,\[CapitalUpsilon]rs}]&,\[CapitalUpsilon]rg,precsol][[1;;2]];
	
	(*The RealSign serves to match the sign convetion of generic orbits solutions in the limit z->0*)				
	<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Es"->RealSign[x]*\[Delta]EE,
	 "Jzs"->RealSign[x]*\[Delta]Lz,
	(* "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},*)
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 "MinoFrequenciesCorrection"->RealSign[x]*{\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]\[Phi]s},
	 "BLFrequenciesCorrection"->RealSign[x]*{\[CapitalUpsilon]rs/\[CapitalUpsilon]tg-\[CapitalUpsilon]rg/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts,\[CapitalUpsilon]\[Phi]s/\[CapitalUpsilon]tg-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts},
	 (*"MinoPrecessionFrequency"->\[CapitalUpsilon]p,
	 "BLPrecessionFrequency"->\[CapitalUpsilon]p/\[CapitalUpsilon]tg,*)
	 (*Keys purely oscillatory part geodesic coordinate time and azimuthal trajectory*)
	 "\[CapitalDelta]trg"->Function[{wr},Evaluate[\[CapitalDelta]integratedFunc[wr,\[CapitalDelta]trgcoeff]]],
	 "\[CapitalDelta]\[Phi]rg"->Function[{wr},Evaluate[\[CapitalDelta]integratedFunc[wr,\[CapitalDelta]\[Phi]rgcoeff]]],
	 (*Keys geodesic radial trajectory*)
	 "rg"->Function[{wr},rgICr2gfun[wr,a,p,e,x]],
	 (*Keys geodesic velocities*)
	 "vtg"->Function[{wr},Evaluate[fourierCos[wr,\[CapitalUpsilon]tg,dtrgd\[Lambda]coeff]]],
	 "vrg"->Function[{wr},drgd\[Lambda]fun[wr,a,p,e,x]],
	 "v\[Phi]g"->Function[{wr},Evaluate[fourierCos[wr,\[CapitalUpsilon]\[Phi]g,d\[Phi]rgd\[Lambda]coeff]]],
	(*Keys corrections trajectory*)
	 "\[CapitalDelta]\[Delta]tpar"->Function[{wr},Evaluate[RealSign[x]*\[CapitalDelta]\[Delta]integratedFunc[wr,\[CapitalUpsilon]rg,\[CapitalUpsilon]rs,\[CapitalDelta]trgcoeff,\[CapitalDelta]tsparcoeff]]],
	 "\[Delta]rpar"->Function[{wr},Evaluate[RealSign[x]*fourierCos[wr,growthrate\[Delta]r,\[Delta]rsparcoeff]]],
	 (*"\[Psi]p"->Function[wr,\[Psi]phase[wr]],*)
	(* "\[Delta]zort"->Function[{wr},\[Delta]zort[wr]],*)
	 "\[CapitalDelta]\[Delta]\[Phi]par"->Function[{wr},Evaluate[RealSign[x]*\[CapitalDelta]\[Delta]integratedFunc[wr,\[CapitalUpsilon]rg,\[CapitalUpsilon]rs,\[CapitalDelta]\[Phi]rgcoeff,\[CapitalDelta]\[Phi]sparcoeff]]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{wr},Evaluate[RealSign[x]*fourierCos[wr,\[CapitalUpsilon]ts,dtspard\[Lambda]coeff]]],
	 "\[Delta]vrpar"->Function[{wr},Evaluate[RealSign[x]*\[Delta]vrfunFTPerPar[wr,a,p,e,x,{\[CapitalUpsilon]rg,\[CapitalUpsilon]rs}]]],
	 "\[Delta]v\[Phi]par"->Function[{wr},Evaluate[RealSign[x]*fourierCos[wr,\[CapitalUpsilon]\[Phi]s,d\[Phi]spard\[Lambda]coeff]]]
	|>
]


(* ::Section::Closed:: *)
(*End package*)


End[];


EndPackage[];
