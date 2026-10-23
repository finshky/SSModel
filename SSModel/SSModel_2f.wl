(* ::Package:: *)

(* ::Title:: *)
(*2 Fluid Sound Shell Model*)


BeginPackage["SSModel2f`",{"SSModel`"}]


ClearAll["SSModel2f`*"]


(* ::Subsection:: *)
(*Available public functions*)


OmegaGW2f;
OmegaGW2fHighK;
OmegaGW2fLowK;
OmegaGW2fPrefLowK;
{\[CapitalDelta]CEi,\[CapitalDelta]SEi,\[CapitalDelta]mn2f,\[CapitalDelta]2f,\[CapitalDelta]02f};
{F,FArg,\[CapitalUpsilon]2f,\[CapitalUpsilon]2fOver1f};
ComputeOmegaGW2f;


OmegaGW2F;
Integrand;
IntegrandP;


(* ::Section:: *)
(*Initialize*)


(* ::Subsection::Closed:: *)
(*Private*)


Begin["`Private`"];


(* ::Subsection::Closed:: *)
(*Plotting*)


(* ::Input::Initialization:: *)
PlotGWithNLOa[vw_,\[Delta]\[Tau]H_,HR_,\[Gamma]_,opts:OptionsPattern[]]:=Module[{sel=Select[#vw==vw\[And]#\[Delta]\[Tau]H==\[Delta]\[Tau]H\[And]#HR==HR\[And]#\[Gamma]==\[Gamma]&],LO,NLO,in},
	LO=data[sel][-1,"KOmega"]//Normal;
	NLO=dataNLOa[sel][-1,"KOmega"]//Normal;
	in=Interpolation[LO];
	PlotGWData[{LO,Abs[NLO],NLO/.{x_,y_}:>{x,y+in[x]}},
	opts//Evaluate,
	Joined->True,
	ImageSize->Medium,
	PlotLegends->{"LO","|NLO|","LO+NLO"},
	PlotStyle->{Dashing[None],Dashing[None],Dashed},
	PlotLabel->StringTemplate["\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\)=``   \[Delta]\[Tau]H=``   HR=``   \[Gamma]=``"][vw,\[Delta]\[Tau]H,HR,\[Gamma]]]
	]
PlotGWithNLOf[vw_,\[Delta]\[Tau]H_,HR_,\[Gamma]_,opts:OptionsPattern[]]:=Module[{sel=Select[#vw==vw\[And]#\[Delta]\[Tau]H==\[Delta]\[Tau]H\[And]#HR==HR\[And]#\[Gamma]==\[Gamma]&],LO,NLO,in},
	LO=data[sel][-1,"KOmega"]//Normal;
	NLO=dataNLOf[sel][-1,"KOmega"]//Normal;
	in=Interpolation[LO];
	PlotGWData[{LO,NLO/.{x_,y_}:>{x,Abs[y-in[x]]},NLO},
	opts//Evaluate,
	Joined->True,
	ImageSize->Medium,
	PlotLegends->{"LO","|\!\(\*SubscriptBox[\(NLO\), \(f\)]\)|","LO+\!\(\*SubscriptBox[\(NLO\), \(f\)]\)"},
	PlotStyle->{Dashing[None],Dashing[None],Dashed},
	PlotLabel->StringTemplate["\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\)=``   \[Delta]\[Tau]H=``   HR=``   \[Gamma]=``"][vw,\[Delta]\[Tau]H,HR,\[Gamma]]]
	]


(* ::Section:: *)
(*Time integrals (\[CapitalDelta] functions)*)


(* ::Subsection::Closed:: *)
(*LO*)


CosExpIntegral[\[Theta]_,q_,\[Gamma]_]:=\[Integral]Cos[q \[Theta]]/\[Theta] E^(-2\[Gamma](\[Theta]-1)) \[DifferentialD]\[Theta]//Evaluate (* 1/2 \[ExponentialE]^(2 \[Gamma]) (ExpIntegralEi[-\[ImaginaryI] q \[Theta]-2 \[Gamma] \[Theta]]+ExpIntegralEi[\[ImaginaryI] q \[Theta]-2 \[Gamma] \[Theta]]) *)
SinExpIntegral[\[Theta]_,q_,\[Gamma]_]:=\[Integral]Sin[q \[Theta]]/\[Theta] E^(-2\[Gamma](\[Theta]-1)) \[DifferentialD]\[Theta]//Simplify//Evaluate (* 1/2 \[ImaginaryI] \[ExponentialE]^(2 \[Gamma]) (ExpIntegralEi[-\[ImaginaryI] q \[Theta]-2 \[Gamma] \[Theta]]-ExpIntegralEi[\[ImaginaryI] q \[Theta]-2 \[Gamma] \[Theta]]) *)
\[CapitalDelta]CEi[\[Delta]\[Theta]_,q_,\[Gamma]_]:=CosExpIntegral[1+\[Delta]\[Theta],q,\[Gamma]]-CosExpIntegral[1,q,\[Gamma]]//Simplify//Evaluate
\[CapitalDelta]SEi[\[Delta]\[Theta]_,q_,\[Gamma]_]:=SinExpIntegral[1+\[Delta]\[Theta],q,\[Gamma]]-SinExpIntegral[1,q,\[Gamma]]//Simplify//Evaluate

(* original implementation *)
\[CapitalDelta]mnArg[\[Tau]H_,q_,{\[Gamma]_,r_}]:=CosIntegral[q \[Tau]H]+I SinIntegral[q \[Tau]H]+E^(2 \[Gamma]) r ExpIntegralE[1,-I q \[Tau]H+2 \[Gamma] \[Tau]H]^2
\[CapitalDelta]mnHOld[\[Delta]\[Tau]H_,Qmn_,{\[Gamma]_,r_}]:=1/4 1/(1+r)^2 Abs[\[CapitalDelta]mnArg[1+\[Delta]\[Tau]H,Qmn,{\[Gamma],r}]-\[CapitalDelta]mnArg[1,Qmn,{\[Gamma],r}]]^2;
\[CapitalDelta]mn2f[\[Delta]\[Theta]_,qmn_,{\[Gamma]_,r_}]:=1/4 1/(1+r)^2 ((\[CapitalDelta]Ci[\[Delta]\[Theta],qmn]+r \[CapitalDelta]CEi[\[Delta]\[Theta],qmn,\[Gamma]])^2+(\[CapitalDelta]Si[\[Delta]\[Theta],qmn]+r \[CapitalDelta]SEi[\[Delta]\[Theta],qmn,\[Gamma]])^2)
\[CapitalDelta]mn2f[\[Delta]\[Tau]H_,{m_Integer,n_Integer},pmnArgs_List,{\[Gamma]_,r_}]:=\[CapitalDelta]mn2f[\[Delta]\[Tau]H,pmn[{m,n},pmnArgs],{\[Gamma],r}];
\[CapitalDelta]2f[\[Delta]\[Tau]H_,pmnArgs_List,{\[Gamma]_,r_}]:=Sum[\[CapitalDelta]mn2f[\[Delta]\[Tau]H,{m,n},pmnArgs,{\[Gamma],r}],{m,{-1,1}},{n,{-1,1}}]//Evaluate


(* ::Subsection::Closed:: *)
(*NLO*)


(* ::Input::Initialization:: *)
ClearAll[\[CurlyEpsilon]]

(* defines \[CurlyEpsilon]s *)
\[CurlyEpsilon][1,z_,cs_]:=(fInt'[z])^2
\[CurlyEpsilon][2,z_,cs_]:=(cs fInt'[z]lInt[z])^2
\[Zeta]kFun["NLO",i_,K_?NumericQ,cs_?NumericQ,HR_?NumericQ,\[Gamma]_?NumericQ]:=3/4 \[Gamma] HR/cs K NIntegrate[ \[Nu][T]T^6 \[CurlyEpsilon][i,T K/R\[Beta],cs]//Evaluate,{T,0,\[Infinity]},Method->"AdaptiveMonteCarlo"] (* assuming Subscript[\[CapitalGamma], v]/(2Overscript[w, _])=3\[Gamma]/4 *)

orders={"LO","NLO1","NLO2"};
Options[Compute\[Zeta]k]={Verbose->True};
Compute\[Zeta]k[HR_,\[Gamma]_,OptionsPattern[]]:=Module[{out,EkPeaks,print=OptionValue[Verbose]},
ClearAll[\[Zeta]k];
If[print,Echo[HR,"HR\[Rule]"];
Echo[\[Gamma],"\[Gamma]\[Rule]"]];
out["LO"]=Interpolate\[Zeta]k[\[Zeta]kFun,"NormalizePeak"->False];
out["NLO1"]=Interpolate\[Zeta]k[\[Zeta]kFun["NLO",1,#,cs,HR,\[Gamma]]&,"NormalizePeak"->False];
out["NLO2"]=Interpolate\[Zeta]k[\[Zeta]kFun["NLO",2,#,cs,HR,\[Gamma]]&,"NormalizePeak"->False];
EkPeaks=Table[out[o]["EkOverRstar"],{o,{"LO","NLO1","NLO2"}}];
EkOverRstar=EkPeaks//First; (* they're all the same *)
Do[\[Zeta]k[o,K_]:=out[o]["\[Zeta]k"][K]//Evaluate,{o,orders}]
(*\[Zeta]k[K_]:=\[Zeta]k["LO"][K]+\[Zeta]k["NLO"][K]*)
]


(* ::Input::Initialization:: *)
(* NLO amplitude *)
\[CapitalDelta]mnH["NLO1",\[Delta]\[Theta]_,qmn_,{\[Gamma]_,r_}]:=r/(6(1+r)^2) (\[CapitalDelta]CEi[\[Delta]\[Theta],qmn,\[Gamma]]\[CapitalDelta]Si[\[Delta]\[Theta],qmn]-\[CapitalDelta]Ci[\[Delta]\[Theta],qmn]\[CapitalDelta]SEi[\[Delta]\[Theta],qmn,\[Gamma]])//Simplify//Evaluate
\[CapitalDelta]mnH["NLO2",\[Delta]\[Theta]_,qmn_,{\[Gamma]_,r_}]:=r/(6(1+r)^2) (\[CapitalDelta]Ci[\[Delta]\[Theta],qmn]\[CapitalDelta]CEi[\[Delta]\[Theta],qmn,\[Gamma]]+\[CapitalDelta]Si[\[Delta]\[Theta],qmn]\[CapitalDelta]SEi[\[Delta]\[Theta],qmn,\[Gamma]]+r(\[CapitalDelta]CEi[\[Delta]\[Theta],qmn,\[Gamma]]^2+\[CapitalDelta]SEi[\[Delta]\[Theta],qmn,\[Gamma]]^2))//Simplify//Evaluate
\[CapitalDelta]mnH["NLO1",\[Delta]\[Theta]_,{m_Integer,n_Integer},pmnArgs_List,{\[Gamma]_,r_}]:=\[CapitalDelta]mnH["NLO1",\[Delta]\[Theta],pmn[{m,n},pmnArgs],{\[Gamma],r}]//Simplify//Evaluate
\[CapitalDelta]mnH["NLO2",\[Delta]\[Theta]_,{m_Integer,n_Integer},pmnArgs_List,{\[Gamma]_,r_}]:=\[CapitalDelta]mnH["NLO2",\[Delta]\[Theta],pmn[{m,n},pmnArgs],{\[Gamma],r}]//Simplify//Evaluate
\[CapitalDelta]H["NLO1",\[Delta]\[Theta]_,{\[ScriptK]_,q_,\[ScriptQ]_,c_},{\[Gamma]_,r_}]:=Re@Sum[\[CapitalDelta]mnH["NLO1",\[Delta]\[Theta],{m,n},{\[ScriptK],q,\[ScriptQ],c},{\[Gamma],r}],{m,{-1,1}},{n,{-1,1}}]//Simplify//Evaluate
\[CapitalDelta]H["NLO2",\[Delta]\[Theta]_,{\[ScriptK]_,q_,\[ScriptQ]_,c_},{\[Gamma]_,r_}]:=Re@Sum[\[CapitalDelta]mnH["NLO2",\[Delta]\[Theta],{m,n},{\[ScriptK],q,\[ScriptQ],c},{\[Gamma],r}],{m,{-1,1}},{n,{-1,1}}]//Simplify//Evaluate


(* ::Input::Initialization:: *)
(* NLO frequency *)
NLOfpiece0pmnPos[pmn_,pmnx_,Gam_,r_,deltau_]:=(CosIntegral[pmn]-CosIntegral[(1+deltau) pmn])^2+(SinIntegral[pmn]-SinIntegral[(1+deltau) pmn])^2;
NLOfpiece0[pmn_,pmnx_,Gam_,r_,deltau_]:=NLOfpiece0pmnPos[Abs[pmn],pmnx,Gam,r,deltau];
NLOfpiece1[pmn_,pmnx_,Gam_,r_,deltau_]:=1/2 E^(2 Gam-I pmnx) r ((ExpIntegralEi[-I pmn]-ExpIntegralEi[-I (1+deltau) pmn]) (ExpIntegralEi[-2 Gam+I (pmn+pmnx)]-ExpIntegralEi[(1+deltau) (-2 Gam+I (pmn+pmnx))])+E^(2 I pmnx) (ExpIntegralEi[I pmn]-ExpIntegralEi[I (1+deltau) pmn]) (ExpIntegralEi[-2 Gam-I (pmn+pmnx)]-ExpIntegralEi[-((1+deltau) (2 Gam+I (pmn+pmnx)))]));
NLOfpiece3[pmn_,pmnx_,Gam_,r_,deltau_]:=E^(4 Gam) r^2 (ExpIntegralEi[-2 Gam+I (pmn+pmnx)]-ExpIntegralEi[(1+deltau) (-2 Gam+I (pmn+pmnx))]) (ExpIntegralEi[-2 Gam-I (pmn+pmnx)]-ExpIntegralEi[-((1+deltau) (2 Gam+I (pmn+pmnx)))]);
delmn2fNLOf[pmn_,pmnx_,Gam_,r_,deltau_]:=NLOfpiece0[pmn,pmnx,Gam,r,deltau]+2*NLOfpiece1[pmn,pmnx,Gam,r,deltau]+NLOfpiece3[pmn,pmnx,Gam,r,deltau];
cs=1/Sqrt[3.];
pmnFun[m_,n_,k_,p_,ptilde_]:=k+m*cs*p+n*cs*ptilde;
(*let's take the benchmark \bar\omega = Gam^2/(8 cs p)*)
pmnxFun[m_,n_,k_,p_,ptilde_,Gam_]:=m*Gam^2/(8*cs*p)+n*Gam^2/(8*cs*ptilde);
delmn2fNLOfmn[m_,n_,k_,p_,ptilde_,Gam_,r_,deltau_]:=delmn2fNLOf[pmnFun[m,n,k,p,ptilde],pmnxFun[m,n,k,p,ptilde,Gam],Gam,r,deltau];
(*this is the result for the function \[CapitalDelta]^(2f,NLO) from the notes up to a factor of 1/4*1/(1+r)^2*)
\[CapitalDelta]H["NLOf",\[Delta]\[Theta]_,{\[ScriptK]_,q_,\[ScriptQ]_,c_},{\[Gamma]_,r_}]:=Re[1/(4(1+r)^2) Sum[delmn2fNLOfmn[m,n,\[ScriptK],q,\[ScriptQ],\[Gamma],r,\[Delta]\[Theta]],{m,{-1,1}},{n,{-1,1}}]];


(* ::Subsection:: *)
(*\[CapitalDelta]0*)


\[CapitalDelta]02f[\[Delta]\[Tau]H_,{\[Gamma]_,r_}]:=.5/(1+r)^2 (Log[1+\[Delta]\[Tau]H] + r E^(2\[Gamma])(ExpIntegralEi[-2 \[Gamma] (1+\[Delta]\[Tau]H)]-ExpIntegralEi[-2 \[Gamma]]))^2


(* ::Section:: *)
(*Gravitational wave spectrum (\[CapitalOmega]_GW)*)


(* ::Subsection:: *)
(*LO*)


(* ::Subsubsection:: *)
(*\[CapitalOmega]_GW (fixed K)*)


OmegaGW2f::usage="OmegaGW2f[K,{\!\(\*SubscriptBox[\(c\), \(s\)]\),\[Delta]\[Tau]H,HR},{\[Gamma],\!\(\*SubscriptBox[\(r\), \(w\)]\)},kineticSpectrum] gives the GW spectrum at given K\[Congruent]R*k.
OmegaGW2f[K,{\!\(\*SubscriptBox[\(c\), \(s\)]\),\[Delta]\[Tau]H,HR},{\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\),\[Alpha]},{\[Gamma],\!\(\*SubscriptBox[\(r\), \(w\)]\)}] computes/loads the kinetic spectrum for given vw, \[Alpha] and cs parameters. 
";
Options[OmegaGW2f]={AccuracyGoal->\[Infinity],PrecisionGoal->2,MaxRecursion->\[Infinity],"RelativisticDOF"->None};
OmegaGW2f[K_?NumericQ,{cs_?NumericQ,\[Delta]\[Tau]H_?NumericQ,HR_?NumericQ},{\[Gamma]_?NumericQ,r_?NumericQ},kinSpectrum_Association,opts:OptionsPattern[]]:=Module[{
	pref,\[Rho]Frac,\[Zeta]k,EkOverRstar,\[CapitalDelta]2fN,optsPref,optsInt
	},
	optsPref=FilterRules[{opts},Options[OmegaGWPref]];
	optsInt=FilterRules[{opts},Options[NIntegrate]];
	{\[Zeta]k,EkOverRstar}=Lookup[kinSpectrum,{"KSpectrum","EkOverRstar"}];
	\[Rho]Frac=1/(1+r);
	pref=OmegaGWPref[\[CapitalGamma]meanAdiabIndex,EkOverRstar,"DensityFraction"->\[Rho]Frac,optsPref//Evaluate];
	\[CapitalDelta]2fN[p_,\[ScriptP]_]:=Re@\[CapitalDelta]2f[\[Delta]\[Tau]H,{K/HR,p/HR,\[ScriptP]/HR,cs},{\[Gamma],r}]//Simplify//Evaluate;
	pref K^2 
	NIntegrate[E^(2logP)\[Zeta]k[E^logP]
		OmegaGWIntegrandP[\[Zeta]k,\[CapitalDelta]2fN,{K,logP},optsInt//Evaluate],
		{logP,-\[Infinity],\[Infinity]},
		optsInt//Evaluate
		]
	]
OmegaGW2f[K_?NumericQ,{cs_?NumericQ,\[Delta]\[Tau]H_?NumericQ,HR_?NumericQ},{vw_?NumericQ,\[Alpha]_?NumericQ},{\[Gamma]_?NumericQ,r_?NumericQ},opts:OptionsPattern[]]:=Module[{
	pref,kinSpectrum,\[Zeta]k,EkOverRstar,
	optsKin,optsPref,optsInt
	},
	optsKin=FilterRules[{opts},Options[GetKineticSpectrum]]//Evaluate;
	kinSpectrum=GetKineticSpectrum[vw,\[Alpha],cs,optsKin//Evaluate];
	OmegaGW2f[K,{cs,\[Delta]\[Tau]H,HR},{\[Gamma],r},kinSpectrum,opts//Evaluate]
	]
AutoComplete[OmegaGW2f]


z[K_,P_,\[ScriptCapitalP]_]:=(K^2+P^2-\[ScriptCapitalP]^2)/(2 P K)
Integrand[\[Zeta]k_,\[CapitalDelta]N_,{K_?NumericQ,P_?NumericQ,\[ScriptCapitalP]_?NumericQ}]:=(1-z[K,P,\[ScriptCapitalP]]^2)^2/\[ScriptCapitalP]^3 \[Zeta]k[\[ScriptCapitalP]]\[CapitalDelta]N[P,\[ScriptCapitalP]]
IntegrandP::usage="IntegrandP[\[Zeta]k,\[CapitalDelta]N,{K,z}] gives the GW spectrum integrated over P but not z, for a given kinetic spectrum \[Zeta]k[K], and a numerical time integral \[CapitalDelta]N[p,\[ScriptP]].";
Options[IntegrandP]={AccuracyGoal->\[Infinity],PrecisionGoal->2,MaxRecursion->\[Infinity]};
IntegrandP[\[Zeta]k_,\[CapitalDelta]N_,{K_?NumericQ,P_?NumericQ},opts:OptionsPattern[]]:=NIntegrate[
	Integrand[\[Zeta]k,\[CapitalDelta]N,{K,P,\[ScriptCapitalP]}],
	{\[ScriptCapitalP],Abs[K-P],K+P},
	opts//Evaluate,
	AccuracyGoal->\[Infinity],
	PrecisionGoal->2,
	MaxRecursion->\[Infinity]
	]
OmegaGW2F::usage="OmegaGW2F[K,{\!\(\*SubscriptBox[\(c\), \(s\)]\),\[Delta]\[Tau]H,HR},{\[Gamma],\!\(\*SubscriptBox[\(r\), \(w\)]\)},kineticSpectrum] gives the GW spectrum at given K\[Congruent]R*k.
OmegaGW2F[K,{\!\(\*SubscriptBox[\(c\), \(s\)]\),\[Delta]\[Tau]H,HR},{\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\),\[Alpha]},{\[Gamma],\!\(\*SubscriptBox[\(r\), \(w\)]\)}] computes/loads the kinetic spectrum for given vw, \[Alpha] and cs parameters. 
";
Options[OmegaGW2F]={AccuracyGoal->\[Infinity],PrecisionGoal->2,MaxRecursion->\[Infinity],"RelativisticDOF"->None};
OmegaGW2F[K_?NumericQ,{cs_?NumericQ,\[Delta]\[Tau]H_?NumericQ,HR_?NumericQ},{\[Gamma]_?NumericQ,r_?NumericQ},kinSpectrum_Association,opts:OptionsPattern[]]:=Module[{
	pref,\[Rho]Frac,\[Zeta]k,EkOverRstar,\[CapitalDelta]2fN,optsPref,optsInt
	},
	optsPref=FilterRules[{opts},Options[OmegaGWPref]];
	optsInt=FilterRules[{opts},Options[NIntegrate]];
	{\[Zeta]k,EkOverRstar}=Lookup[kinSpectrum,{"KSpectrum","EkOverRstar"}];
	\[Rho]Frac=1/(1+r);
	pref=OmegaGWPref[\[CapitalGamma]meanAdiabIndex,EkOverRstar,"DensityFraction"->\[Rho]Frac,optsPref//Evaluate];
	\[CapitalDelta]2fN[p_,\[ScriptP]_]:=Re@\[CapitalDelta]2f[\[Delta]\[Tau]H,{K/HR,p/HR,\[ScriptP]/HR,cs},{\[Gamma],r}]//Simplify//Evaluate;
	pref K^2 
	NIntegrate[P \[Zeta]k[P]
		IntegrandP[\[Zeta]k,\[CapitalDelta]2fN,{K,P},optsInt//Evaluate],
		{P,0,\[Infinity]},
		optsInt//Evaluate
		]
	]
OmegaGW2F[K_?NumericQ,{cs_?NumericQ,\[Delta]\[Tau]H_?NumericQ,HR_?NumericQ},{vw_?NumericQ,\[Alpha]_?NumericQ},{\[Gamma]_?NumericQ,r_?NumericQ},opts:OptionsPattern[]]:=Module[{
	pref,kinSpectrum,\[Zeta]k,EkOverRstar,
	optsKin,optsPref,optsInt
	},
	optsKin=FilterRules[{opts},Options[GetKineticSpectrum]]//Evaluate;
	kinSpectrum=GetKineticSpectrum[vw,\[Alpha],cs,optsKin//Evaluate];
	OmegaGW2F[K,{cs,\[Delta]\[Tau]H,HR},{\[Gamma],r},kinSpectrum,opts//Evaluate]
	]
AutoComplete[OmegaGW2F]


(* ::Subsubsection:: *)
(*Compute spectrum*)


$GWKeyOrder={"vw","\[Alpha]","cs","SourceDuration","HR","\[Gamma]","rw","R\[Beta]","EkOverRstar","KSpectrumIntegral","KEnergyDensityFraction","\[ScriptCapitalC]","EddyTurnover","KOmega","KTime","Metadata"};
Options[ComputeOmegaGW2f]=Options[ComputeOmegaGW]~Join~Options[OmegaGW2f];
ComputeOmegaGW2f[{cs_,\[Delta]\[Tau]H_,HR_},{vw_,\[Alpha]_},{\[Gamma]_,r_},opts:OptionsPattern[$dependencies[ComputeOmegaGW2f]]]:=Module[{
	K\[CapitalOmega],KTime,KSpectrum,GWSpectrum,\[Delta]\[Tau]Hnl,store,
	optsKin,optsGW,optsComp},
	optsKin=FilterOptions[opts,GetKineticSpectrum];
	optsGW=FilterOptions[opts,OmegaGW];
	optsComp=FilterOptions[opts,ComputeOmegaGW];
	store=OptionValue["Store"]/.True->$GWFile;
	(* extract kinetic spectrum (fixed vw and \[Alpha]) *)
	KSpectrum=GetKineticSpectrum[vw,\[Alpha],cs,optsKin//Evaluate];
	(* compute spectrum *)
	K\[CapitalOmega]=ComputeOmegaGW[Function[K,OmegaGW2f[K,{cs,\[Delta]\[Tau]H,HR},{\[Gamma],r},KSpectrum,optsGW//Evaluate]],
		optsComp//Evaluate];
	{K\[CapitalOmega],KTime}={K\[CapitalOmega][[All,;;2]],K\[CapitalOmega][[All,{1,3}]]};
	(* create spectrum object *)
	GWSpectrum=Select[KSpectrum,NumericQ];
	\[Delta]\[Tau]Hnl=EddyTurnover[HR,GWSpectrum["KEnergyDensityFraction"]];
	AppendTo[GWSpectrum,<|
		"SourceDuration"->\[Delta]\[Tau]H,
		"HR"->HR,
		"\[Gamma]"->\[Gamma],
		"rw"->r,
		"EddyTurnover"->\[Delta]\[Tau]Hnl,
		"KOmega"->K\[CapitalOmega],
		"KTime"->KTime,
		"Metadata"->OptionValue["Metadata"]
		|>];
	If[Equal@@Length/@{$GWKeyOrder,GWSpectrum},GWSpectrum=KeyTake[GWSpectrum,$GWKeyOrder]];
	(* optionally store *)
	If[MatchQ[store,_String],AppendAssociation[store,GWSpectrum]];
	GWSpectrum
	]
$dependencies[ComputeOmegaGW2f]={ComputeOmegaGW2f,OmegaGW2f,OmegaGWPref,GetKineticSpectrum};


PlotGWithNLOa[vw_,\[Delta]\[Tau]H_,HR_,\[Gamma]_,opts:OptionsPattern[]]:=Module[{sel=Select[#vw==vw\[And]#\[Delta]\[Tau]H==\[Delta]\[Tau]H\[And]#HR==HR\[And]#\[Gamma]==\[Gamma]&],LO,NLO,in},
	LO=data[sel][-1,"KOmega"]//Normal;
	NLO=dataNLOa[sel][-1,"KOmega"]//Normal;
	in=Interpolation[LO];
	PlotGWData[{LO,Abs[NLO],NLO/.{x_,y_}:>{x,y+in[x]}},
	opts//Evaluate,
	Joined->True,
	ImageSize->Medium,
	PlotLegends->{"LO","|NLO|","LO+NLO"},
	PlotStyle->{Dashing[None],Dashing[None],Dashed},
	PlotLabel->StringTemplate["\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\)=``   \[Delta]\[Tau]H=``   HR=``   \[Gamma]=``"][vw,\[Delta]\[Tau]H,HR,\[Gamma]]]
	]
PlotGWithNLOf[vw_,\[Delta]\[Tau]H_,HR_,\[Gamma]_,opts:OptionsPattern[]]:=Module[{sel=Select[#vw==vw\[And]#\[Delta]\[Tau]H==\[Delta]\[Tau]H\[And]#HR==HR\[And]#\[Gamma]==\[Gamma]&],LO,NLO,in},
	LO=data[sel][-1,"KOmega"]//Normal;
	NLO=dataNLOf[sel][-1,"KOmega"]//Normal;
	in=Interpolation[LO];
	PlotGWData[{LO,NLO/.{x_,y_}:>{x,Abs[y-in[x]]},NLO},
	opts//Evaluate,
	Joined->True,
	ImageSize->Medium,
	PlotLegends->{"LO","|\!\(\*SubscriptBox[\(NLO\), \(f\)]\)|","LO+\!\(\*SubscriptBox[\(NLO\), \(f\)]\)"},
	PlotStyle->{Dashing[None],Dashing[None],Dashed},
	PlotLabel->StringTemplate["\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\)=``   \[Delta]\[Tau]H=``   HR=``   \[Gamma]=``"][vw,\[Delta]\[Tau]H,HR,\[Gamma]]]
	]


(* ::Subsection::Closed:: *)
(*NLO*)


(* ::Input::Initialization:: *)
(* GW: integrals over momenta *)
\[CapitalOmega]GWIntegrandHR["NLOa",\[CapitalDelta]HN_,{K_?NumericQ,P_?NumericQ,z_?NumericQ}]:=With[{\[ScriptCapitalP]=\[Sqrt](P^2+K^2-2P K z)},
(1-z^2)^2 P^2/\[ScriptCapitalP]^4 (
(\[Zeta]k["LO",P]\[Zeta]k["NLO1",\[ScriptCapitalP]]\[CapitalDelta]HN["NLO1",P,\[ScriptCapitalP]]+\[Zeta]k["LO",\[ScriptCapitalP]]\[Zeta]k["NLO1",P]\[CapitalDelta]HN["NLO1",\[ScriptCapitalP],P])+
(\[Zeta]k["LO",P]\[Zeta]k["NLO2",\[ScriptCapitalP]]+\[Zeta]k["LO",\[ScriptCapitalP]]\[Zeta]k["NLO2",P])\[CapitalDelta]HN["NLO2",P,\[ScriptCapitalP]]
)]
\[CapitalOmega]GWIntegrandHRNorm["NLOa",\[CapitalDelta]HN_,{K_?NumericQ,P_?NumericQ,z_?NumericQ}]:=\[CapitalOmega]GWPref[\[CapitalGamma]meanAdiabIndex,1.,EkOverRstar] K^3 \[CapitalOmega]GWIntegrandHR["NLOa",\[CapitalDelta]HN,{K,P,z}]//Evaluate
(*\[CapitalOmega]GWIntHR[\[Delta]\[Tau]H_?NumericQ,K_?NumericQ,HR_?NumericQ,{\[Gamma]_?NumericQ,r_?NumericQ}]:=\[CapitalOmega]GWPref[\[CapitalGamma]meanAdiabIndex,1.,EkOverRstar]K^3NIntegrate[\[CapitalOmega]GWIntegrandHR[\[Delta]\[Tau]H,{K,P,z,cs},HR,{\[Gamma],r}],{z,-1,1},{P,Pm,PM},Method->"GlobalAdaptive",MaxRecursion->10,AccuracyGoal->\[Infinity],PrecisionGoal->2]*)
\[CapitalOmega]GWIntP["NLOa",\[CapitalDelta]HN_,{K_?NumericQ,z_?NumericQ}]:=NIntegrate[\[CapitalOmega]GWIntegrandHR["NLOa",\[CapitalDelta]HN,{K,P,z}],{P,Pm,PM},
	Method->"GlobalAdaptive",
	MaxRecursion->3,
	AccuracyGoal->\[Infinity],
	PrecisionGoal->2
	]
\[CapitalOmega]GWIntHR["NLOa",\[Delta]\[Tau]H_?NumericQ,K_?NumericQ,cs_?NumericQ,HR_?NumericQ,{\[Gamma]_?NumericQ,r_?NumericQ}]:=Module[{\[CapitalDelta]HN},
\[CapitalDelta]HN["NLO1",p_,\[ScriptP]_]:=\[CapitalDelta]H["NLO1",\[Delta]\[Tau]H,{K/HR,p/HR,\[ScriptP]/HR,cs},{\[Gamma],r}]//Simplify//Evaluate;
\[CapitalDelta]HN["NLO2",p_,\[ScriptP]_]:=\[CapitalDelta]H["NLO2",\[Delta]\[Tau]H,{K/HR,p/HR,\[ScriptP]/HR,cs},{\[Gamma],r}]//Simplify//Evaluate;
(\[CapitalOmega]GWPref[\[CapitalGamma]meanAdiabIndex,1.,EkOverRstar] K^3 
	NIntegrate[\[CapitalOmega]GWIntP["NLOa",\[CapitalDelta]HN,{K,z}],{z,-.999,.999},
	Method->"LocalAdaptive",
	MaxRecursion->1,
	AccuracyGoal->\[Infinity],
	PrecisionGoal->2
	])]


(* ::Input::Initialization:: *)
(* NLO frequency *)
\[CapitalOmega]GWIntegrandHR["NLOf",\[CapitalDelta]HN_,{K_?NumericQ,P_?NumericQ,z_?NumericQ}]:=With[{\[ScriptCapitalP]=\[Sqrt](P^2+K^2-2P K z)},
(1-z^2)^2 P^2/\[ScriptCapitalP]^4 \[Zeta]k[P]\[Zeta]k[\[ScriptCapitalP]]\[CapitalDelta]HN[P,\[ScriptCapitalP]]]
(*\[CapitalOmega]GWIntHR[\[Delta]\[Tau]H_?NumericQ,K_?NumericQ,HR_?NumericQ,{\[Gamma]_?NumericQ,r_?NumericQ}]:=\[CapitalOmega]GWPref[\[CapitalGamma]meanAdiabIndex,1.,EkOverRstar]K^3NIntegrate[\[CapitalOmega]GWIntegrandHR[\[Delta]\[Tau]H,{K,P,z,cs},HR,{\[Gamma],r}],{z,-1,1},{P,Pm,PM},Method->"GlobalAdaptive",MaxRecursion->10,AccuracyGoal->\[Infinity],PrecisionGoal->2]*)
Options[\[CapitalOmega]GWIntP]={MaxRecursion->3};
\[CapitalOmega]GWIntP["NLOf",\[CapitalDelta]HN_,{K_?NumericQ,z_?NumericQ},opts:OptionsPattern[]]:=NIntegrate[\[CapitalOmega]GWIntegrandHR["NLOf",\[CapitalDelta]HN,{K,P,z}],{P,Pm,PM},
	Method->"GlobalAdaptive",
	opts,
	AccuracyGoal->\[Infinity],
	PrecisionGoal->2
	]
Options[\[CapitalOmega]GWIntHR]={MaxRecursion->3};
\[CapitalOmega]GWIntHR["NLOf",\[Delta]\[Tau]H_?NumericQ,K_?NumericQ,cs_?NumericQ,HR_?NumericQ,{\[Gamma]_?NumericQ,r_?NumericQ},opts:OptionsPattern[]]:=Module[{\[CapitalDelta]HN},
\[CapitalDelta]HN[p_,\[ScriptP]_]:=\[CapitalDelta]H["NLOf",\[Delta]\[Tau]H,{K/HR,p/HR,\[ScriptP]/HR,cs},{\[Gamma],r}]//Simplify//Evaluate;
(\[CapitalOmega]GWPref[\[CapitalGamma]meanAdiabIndex,1.,EkOverRstar] K^3 
	NIntegrate[\[CapitalOmega]GWIntP["NLOf",\[CapitalDelta]HN,{K,z},opts],{z,-.999,.999},
	Method->"LocalAdaptive",
	opts,
	AccuracyGoal->\[Infinity],
	PrecisionGoal->2
	])]

dataNLOa=Dataset[<<data/GWSpectra.NLOa.wl];
dataNLOf=Dataset[<<data/GWSpectra.NLOf.wl];


(* ::Subsection:: *)
(*\[Delta]-approximations (HH19)*)


OmegaGW2fHighK::usage="OmegaGW2fHighK[\!\(\*SubscriptBox[\(kR\), \(*\)]\),{\!\(\*SubscriptBox[\(c\), \(s\)]\),\[Delta]\[Tau]H,HR},{\[Gamma],\!\(\*SubscriptBox[\(r\), \(w\)]\)},kineticSpectrum] gives the GW spectrum at given K\[Congruent]R*k.
OmegaGW2fHighK[\!\(\*SubscriptBox[\(kR\), \(*\)]\),{\!\(\*SubscriptBox[\(c\), \(s\)]\),\[Delta]\[Tau]H,HR},{\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\),\[Alpha]},{\[Gamma],\!\(\*SubscriptBox[\(r\), \(w\)]\)}] computes/loads the kinetic spectrum for given vw, \[Alpha] and cs parameters. 
";
Options[OmegaGW2fHighK]={AccuracyGoal->\[Infinity],PrecisionGoal->2,MaxRecursion->\[Infinity],"RelativisticDOF"->None};
OmegaGW2fHighK[K_?NumericQ,{cs_?NumericQ,\[Delta]\[Tau]H_?NumericQ,HR_?NumericQ},{\[Gamma]_?NumericQ,r_?NumericQ},kinSpectrum_Association,opts:OptionsPattern[]]:=Module[{
	\[Rho]Frac
	},
	\[Rho]Frac=1/(1+r);
	\[CapitalUpsilon]2fOver1f[\[Delta]\[Tau]H,{\[Gamma],r}] OmegaGWHighK[K,{cs,\[Delta]\[Tau]H,HR},kinSpectrum,"DensityFraction"->\[Rho]Frac,opts//Evaluate]
	]
OmegaGW2fHighK[K_?NumericQ,{cs_?NumericQ,\[Delta]\[Tau]H_?NumericQ,HR_?NumericQ},{vw_?NumericQ,\[Alpha]_?NumericQ},{\[Gamma]_?NumericQ,r_?NumericQ},opts:OptionsPattern[]]:=Module[{
	kinSpectrum,optsKin
	},
	optsKin=FilterRules[Join[{opts},Options[OmegaGW2fHighK]],Options[GetKineticSpectrum]]//Evaluate;
	kinSpectrum=GetKineticSpectrum[vw,\[Alpha],cs,optsKin//Evaluate];
	OmegaGW2fHighK[K,{cs,\[Delta]\[Tau]H,HR},{\[Gamma],r},kinSpectrum,opts//Evaluate]
	]


Options[OmegaGW2fPrefLowK]=Options[OmegaGWPref];
OmegaGW2fPrefLowK[\[CapitalGamma]meanAdiabIndex_,EkOverRstar_,\[ScriptCapitalC]_,\[Delta]\[Tau]H_,{\[Gamma]_,r_},opts:OptionsPattern[]]:=\[ScriptCapitalC] \[CapitalDelta]02f[\[Delta]\[Tau]H,{\[Gamma],r}] OmegaGWPref[\[CapitalGamma]meanAdiabIndex,EkOverRstar,opts//Evaluate]


OmegaGW2fLowK::usage="OmegaGW2fLowK[k SubStar[R],\[Delta]\[Tau]H,kineticSpectrum,{\[Gamma],r}] gives the GW spectrum in the low-K\[Congruent]R*k limit.
OmegaGW2fLowK[k SubStar[R],{\!\(\*SubscriptBox[\(c\), \(s\)]\),\[Delta]\[Tau]H},{\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\),\[Alpha]},{\[Gamma],r}] computes/loads the kinetic spectrum for given vw, \[Alpha] and cs parameters. 
";
Options[OmegaGW2fLowK]={
	"RelativisticDOF"->None,
	"Verbose"->False
	};
OmegaGW2fLowK[K_,\[Delta]\[Tau]H_?NumericQ,kinSpectrum_Association,{\[Gamma]_?NumericQ,r_?NumericQ},opts:OptionsPattern[]]:=Module[{
	pref,\[Zeta]k,EkOverRstar,\[ScriptCapitalC],optsPref
	},
	optsPref=FilterRules[Join[{opts},Options[OmegaGW2fLowK]],Options[OmegaGWPref]];
	{\[Zeta]k,EkOverRstar,\[ScriptCapitalC]}=Lookup[kinSpectrum,{"KSpectrum","EkOverRstar","\[ScriptCapitalC]"}];
	pref=OmegaGW2fPrefLowK[\[CapitalGamma]meanAdiabIndex,EkOverRstar,\[ScriptCapitalC],\[Delta]\[Tau]H,{\[Gamma],r},optsPref//Evaluate];
	pref K^3
	]
OmegaGW2fLowK[K_,\[Delta]\[Tau]H_?NumericQ,{vw_?NumericQ,\[Alpha]_?NumericQ,cs_?NumericQ},{\[Gamma]_?NumericQ,r_?NumericQ},opts:OptionsPattern[]]:=Module[{
	pref,kinSpectrum,\[Zeta]k,EkOverRstar,
	optsKin,optsPref
	},
	optsKin=FilterRules[Join[{opts},Options[OmegaGW2fLowK]],Options[GetKineticSpectrum]]//Evaluate;
	kinSpectrum=GetKineticSpectrum[vw,\[Alpha],cs,optsKin//Evaluate];
	OmegaGW2fLowK[K,\[Delta]\[Tau]H,kinSpectrum,{\[Gamma],r},opts//Evaluate]
	]
$dependencies[OmegaGW2fLowK]={OmegaGW2fLowK,OmegaGWPref,GetKineticSpectrum};
AutoComplete[OmegaGW2fLowK,$dependencies[OmegaGW2fLowK]];


(* ::Section::Closed:: *)
(*Overall prefactor (\[CapitalUpsilon]_2f)*)


(* ::Text:: *)
(*The expression for the 2-fluid prefactor \[CapitalUpsilon]2f is obtained from*)
(*\!\( *)
(*\*SubsuperscriptBox[\(\[Integral]\), \(1\), \(2\)]\(Simplify[*)
(*\*FractionBox[*)
(*SuperscriptBox[\(F[\[Tau] - 1]\), \(2\)], *)
(*SuperscriptBox[\(\[Tau]\), \(2\)]]] \[DifferentialD]\[Tau]\)\)*)


FArg[d\[Tau]H_,{\[Gamma]_,r_}]:=1+r E^(-2\[Gamma] d\[Tau]H)
F[d\[Tau]H_,{\[Gamma]_,r_}]:=FArg[d\[Tau]H,{\[Gamma],r}]/(1+r) (* separate the constant term *)


(*Assuming[{_\[Element]Reals,_>0},\!\(
\*SubsuperscriptBox[\(\[Integral]\), \(1\), \(1 + d\[Theta]\)]\(
\*FractionBox[
SuperscriptBox[\(F[\(-1\) + \[Theta], {\[Gamma], r}]\), \(2\)], 
SuperscriptBox[\(\[Theta]\), \(2\)]]\[DifferentialD]\[Theta]\)\)]*)
(*Assuming[{_\[Element]Reals,_>0},\!\(
\*SubsuperscriptBox[\(\[Integral]\), \(1\), \(1 + d\[Theta]\)]\(
\*FractionBox[
SuperscriptBox[\(\[ExponentialE]\), \(\(-\[Gamma]\)\((\[Theta] - 1)\)\)], 
SuperscriptBox[\(\[Theta]\), \(2\)]]\[DifferentialD]\[Theta]\)\)]*)
g[\[Delta]\[Tau]H_,\[Gamma]_]:=\[Gamma] E^\[Gamma] Gamma[-1,\[Gamma],\[Gamma](1+\[Delta]\[Tau]H)]
\[CapitalUpsilon]2f::usage="\[CapitalUpsilon]2f[\!\(\*SubscriptBox[\(\[Delta]\[Tau]\), \(fin\)]\)\!\(\*SubscriptBox[\(\[ScriptCapitalH]\), \(*\)]\),{\[Gamma],r}] gives the suppression factor for the 2-fluid system.";
\[CapitalUpsilon]2f[\[Delta]\[Tau]H_,{\[Gamma]_,r_}]:=1/(1+r)^2 (\[CapitalUpsilon][\[Delta]\[Tau]H]+2r g[\[Delta]\[Tau]H,2\[Gamma]]+r^2 g[\[Delta]\[Tau]H,4\[Gamma]]);
\[CapitalUpsilon]2fOver1f::usage="\[CapitalUpsilon]2fOver1f[\!\(\*SubscriptBox[\(\[Delta]\[Tau]\), \(fin\)]\)\!\(\*SubscriptBox[\(\[ScriptCapitalH]\), \(*\)]\),{\[Gamma],r}] gives the ratio of suppression factors \[CapitalUpsilon]2f/\[CapitalUpsilon]";
\[CapitalUpsilon]2fOver1f[\[Delta]\[Tau]H_,{\[Gamma]_,r_}]:=\[CapitalUpsilon]2f[\[Delta]\[Tau]H,{\[Gamma],r}]/\[CapitalUpsilon][\[Delta]\[Tau]H]//Simplify//Evaluate


(* ::Section:: *)
(*End*)


End[] (* Private` *)


EndPackage[]
