(* ::Package:: *)

(* ::Title:: *)
(*Sound Shell Model*)


(* ::Text:: *)
(*Implementation of the Sound shell model  (Hindmarsh 2016, Hindmarsh & Hijazi 2019)*)
(*The Sound shell model, in essence, aims to capture the linear regime of in a cosmological 1^st-order phase transition producing a stochastic gravitational wave background through compressional waves.*)
(*This notebook follows the formulation from Roper Pol, Procacci & Caprini 2024.*)


BeginPackage["SSModel`"]


ClearAll["SSModel`*"]


(* ::Subsection:: *)
(*Available public Symbols*)


OmegaGWPref;
OmegaGWIntegrand;
OmegaGWIntegrandP;
OmegaGW;
ComputeOmegaGW;
GetKineticSpectrum;
{\[CapitalDelta]Ci,\[CapitalDelta]Si,pmn,\[CapitalDelta]mn,\[CapitalDelta],\[CapitalDelta]0};
EddyTurnover;


\[Nu];
RBeta;
EkOverRstarFun;
\[ScriptCapitalC]Fun;
KineticSpectrum;
InterpolateKineticSpectrum;


\[CapitalUpsilon];
OmegaGWHighK;
OmegaGWLowK;
OmegaGWPrefHighK;
OmegaGWPrefLowK;


cs2;
cs;
\[CapitalGamma]meanAdiabIndex;


AutoComplete;
AppendAssociation;


(* file paths *)
$DataFolder;
$KineticFile;
$GWFile;


(* plot *)
PlotKineticSpectrum;
PlotGWData;
PlotTiming;


(* other *)
$PlotOptions;
$Verbose;
$KOmega;
FilterOptions;
LogSubdivide;


(* ::Section:: *)
(*Initialize*)


(* ::Subsection:: *)
(*Bubble profiler (v,w,\[Lambda])*)


$SSModelFolder=DirectoryName[$InputFileName];
$ProFile=FileNameJoin[$SSModelFolder,"BubbleFluidProfiles.wl"];


(* ::Text:: *)
(*Load the fluid profiler for the FOPT bubbles*)


(* ::Input::Initialization:: *)
Get[$ProFile];
Print["Loaded",$ProFile];


(* ::Subsection::Closed:: *)
(*Private*)


Begin["`Private`"];


(* ::Subsection::Closed:: *)
(*Messages*)


cs2::careful="Redefining \!\(\*
StyleBox[\"cs2\",\nFontWeight->\"Bold\"]\) doesn't modify the value of \!\(\*
StyleBox[\"cs\",\nFontWeight->\"Bold\"]\). Both are used in the implementation.";


msg["dataFolder"]=StringTemplate["$DataFolder folder set to ``.\nRedefine $DataFolder to set new data path."];


msg["notDataset"]=StringTemplate["Content of `` is not a Dataset. Association can't be appended."];
msg["samePars"][file_String,assoc_Association]:=Module[{num},
	num=Select[assoc,NumericQ];
	Sequence@@{StringTemplate["File `` already contains data with parameters:\n"][file],Dataset[num]}
	]
msg["saveKin"]=StringTemplate["Kinetic spectrum with \!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\)=``, \[Alpha]=`` and \!\(\*SubscriptBox[\(c\), \(s\)]\)=`` already present in ``."];


msg["interpolAmps"]=StringTemplate["Interpolating amplitudes f(z),l(z)\[Rule]\!\(\*SuperscriptBox[\(A\), \(2\)]\)(z) with parameters \[Xi]w=``, \[Alpha]=``, cs=``."];
msg["interpolKin"]=StringTemplate["Interpolating kinetic spectrum \!\(\*SubscriptBox[\(\[Zeta]\), \(kin\)]\)(K)."]
msg["noKin"]=StringTemplate["No kinetic spectrum with \!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\)=``, \[Alpha]=``, \!\(\*SubscriptBox[\(c\), \(s\)]\)=`` in file ``.\nInterpolating via \!\(\*SubscriptBox[\(v\), \(ip\)]\)(\[Xi]),\!\(\*SubscriptBox[\(\[Lambda]\), \(ip\)]\)(\[Xi]) \[Rule] f(z),l(z) \[Rule] \!\(\*SuperscriptBox[\(A\), \(2\)]\)(z) \[Rule] \!\(\*SubscriptBox[\(\[Zeta]\), \(kin\)]\)(K)."]
msg["computeKin"]=StringTemplate["Computing kinetic spectrum with parameters \[Xi]w=``, \[Alpha]=``, cs=``, and storing in ``."]
msg["computedKin"][assoc_Association]:=Sequence["Computed kinetic spectrum with parameters:\n",Dataset[Select[assoc,NumericQ]]]
msg["loadedKin"][assoc_Association]:=Sequence["Loaded kinetic spectrum with parameters:\n",Dataset[Select[assoc,NumericQ]]]
msg["multiKin"]=StringTemplate["`` kinetic spectra with parameters \[Xi]w=``, \[Alpha]=``, cs=`` in ``.\nSelecting first one."]


msg["nComputed"]=StringTemplate["`1`/`2` (`3`) points computed"]


lbl["plotGW"]=StringTemplate["\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\)=``   \[Alpha]=``   \!\(\*SubscriptBox[\(\[Delta]\[Tau]\), \(fin\)]\)\!\(\*SubscriptBox[\(\[ScriptCapitalH]\), \(*\)]\)=``   \!\(\*SubscriptBox[\(\[ScriptCapitalH]\), \(*\)]\)\!\(\*SubscriptBox[\(R\), \(*\)]\)=``"]


(* ::Subsection::Closed:: *)
(*Paths*)


$DataFolder=FileNameJoin[$SSModelFolder,"data"];
$KineticFile=FileNameJoin[$DataFolder,"kinetic.wl"];
$GWFile=FileNameJoin[$DataFolder,"gw.wl"];
Print[msg["dataFolder"][$DataFolder]];


(* ::Subsection::Closed:: *)
(*Autocomplete options*)


(* ::Text:: *)
(*A function to add "code completion" to any symbol. *)
(*This looks up the Options in the symbol's definition which key matches a String, and displays a drop-down menu with all the options, once " is typed (after its arguments).*)
(*Example usage:*)
(*	mySymbol[x_,"myOption"->True]:=..*)
(*	AutoComplete[mySymbol]*)
(*Then writing*)
(*	mySymbol[x,"*)
(*will show a drop-down menu with "myOption"*)


(* ::Input::Initialization:: *)
(* 
https://resources.wolframcloud.com/FunctionRepository/resources/\
AddCodeCompletion
https://mathematica.stackexchange.com/questions/56984
*)
addCodeCompletion[function_String][args___] := Module[
   	{processed},
   	processed = {args} /. {
      		None -> 0,
      		"AbsoluteFileName" -> 2,
      		"RelativeFileName" -> 3,
      		"Color" -> 4,
      		"PackageName" -> 7,
      		"DirectoryName" -> 8,
      		"InterpreterType" -> 9
      	};
   	(FE`Evaluate[FEPrivate`AddSpecialArgCompletion[#1]] &)[
    function -> processed]
   ];
SetAttributes[AutoComplete,HoldFirst];
AutoComplete[sym_]:=Module[{name,keys,keySeq,min,max},
Switch[sym,_Symbol,name=ToString@sym;
	keys=Keys@Options@sym;
	(*count min/max number of arguments*){min,max}=MinMax[Length/@List@@ReleaseHold[DownValues[sym][[All,1]]]]-1/. {\[Infinity],-\[Infinity]}->{0,0};
	keySeq=Sequence@@Join[Table[0,min],ConstantArray[keys,Length@keys+(max-min)]],_Association,name=ToString@HoldForm@sym;
	keySeq=Keys@sym,_Transition|_ActionFunction,name=ToString@HoldForm@sym;
	keySeq=Keys@sym[Association]];
	addCodeCompletion[name][keySeq]]
	(*autocomplete options from multiple functions*)
AutoComplete[sym_Symbol,functions_List]:=Module[{name,keys,keySeq},name=ToString@sym;
	keys=Keys[Union@@Options/@Union[{sym},functions]];
	keySeq=Sequence@@ConstantArray[keys,Length@keys];
	addCodeCompletion[name][keySeq]]


(* ::Subsection:: *)
(*Utilities*)


(* ::Text:: *)
(*Math functions*)


LogN[0]=-\[Infinity];
LogN[0.]=-\[Infinity];
LogN[x_]=Log[x];


(* ::Text:: *)
(*Filter options by function*)


FilterOptions[opts___,funs_List]:=FilterRules[{opts},Union@@Options/@funs]
FilterOptions[opts___,fun_Symbol]:=Module[{dep=$dependencies[fun]},
	If[MatchQ[dep,_List],
		FilterOptions[opts,dep], (* attempt to filter Options from $dependencies[fun] *)
		FilterRules[{opts},Options[fun]]
		]]


(* ::Text:: *)
(*Logarithmic uniform list*)


(* ::Input::Initialization:: *)
LogSubdivide[start_,stop_,n_]:=Exp[Subdivide[Log[start],Log[stop],n]]//N
LogSubdivide[start_,stop_]:=LogSubdivide[start,stop,30]


(* ::Text:: *)
(*Extend a spectrum with by extrapolating 2 points at "Left" and "Right", assuming a power-law behavior.*)


(* ::Input::Initialization:: *)
Options[ExtendSpectrum]={"Left"->10.^-10,"Right"->10.^10,"NPoints"->10};
ExtendSpectrum[data_List,OptionsPattern[]]:=Module[{
	x,vLeft,vRight,
	left=OptionValue["Left"],right=OptionValue["Right"],
	N=OptionValue["NPoints"]
	},
vLeft=E^LinearModelFit[Log@data[[;;10]],x,x][Log@left];
vRight=E^LinearModelFit[Log@data[[-10;;]],x,x][Log@right];
Join[{{left,vLeft}},data,{{right,vRight}}]
]


(* ::Input::Initialization:: *)
Options[SmoothenData]={"Plot"->False,"SplineDegree"->15,"NPoints"->100};
SmoothenData[data_List,OptionsPattern[]]:=Module[{dataSmooth,
	spDegree=OptionValue["SplineDegree"],
	npoints=OptionValue["NPoints"]},
	dataSmooth=Exp[BSplineFunction[Log@data,SplineDegree->spDegree]/@Subdivide[0,1,npoints]];
	If[OptionValue["Plot"],Print@ListLogLogPlot[dataSmooth,Joined->True]];
	dataSmooth
	]


$Verbose=True;
$SSPrint[expr__]:=If[$Verbose,Print[expr],None]
$SSPrint[expr__,verbose:(True|False)]:=If[$Verbose\[And]verbose,Print[expr],None]


(* ::Subsection:: *)
(*Storing data*)


(* ::Input::Initialization:: *)
Options[AppendAssociation]={"Overwrite"->False};
AppendAssociation[file_String,assoc_Association,OptionsPattern[]]:=Module[{
data,matching,crit,overwrite=OptionValue["Overwrite"]},
	If[!FileExistsQ[file],$SSPrint@StringTemplate["New file `` generated"][file];
	Put[Dataset[{assoc}],file];Return[]];
	data=Get[file];
	If[!MatchQ[data,_Dataset],$SSPrint[msg["notDataset"][file]];Return[]];
	(* find all the rows in data matching numerical values in assoc *)
	crit=Function[row,And@@KeyValueMap[(row[#1]==#2)&,Select[assoc,NumericQ]]];
	matching=data[Select[crit]];
	If[Length[matching]>0,$SSPrint[msg["samePars"][file,assoc]];
		If[overwrite===False,$SSPrint["Aborted"];Return[]];
		If[overwrite===True,$SSPrint["Overwriting"];
			(* remove previous data matching selection *)	
			data=data[Select[Not@*crit]],
			$SSPrint["Appending data with same parameters."]
		] 
		];
	AppendTo[data,assoc];
	Put[data,file]
	]
AutoComplete[AppendAssociation];


(* ::Subsection::Closed:: *)
(*Plotting*)


ClearAll[PlotKineticSpectrum]
Options[PlotKineticSpectrum]={"KSpectrumData"->True,"KRange"->{10^-3,10^4}};
PlotKineticSpectrum[kinSpectrum_Association,opts:OptionsPattern[]]:=Module[{spec,data,epilog,
	Km,KM,optPlot},
	optPlot=FilterRules[{opts},Options[LogLogPlot]];
	{Km,KM}=OptionValue["KRange"];
	spec=kinSpectrum["KSpectrum"];
	epilog=If[OptionValue["KSpectrumData"],
		data=InterpolationData[kinSpectrum["LogKSpectrum"]];
		{Black,Point[data]},
		{}];
	LogLogPlot[spec[K],{K,Km,KM},
		Frame->True,GridLines->Automatic,FrameLabel->{"\!\(\*SubscriptBox[\(R\), \(*\)]\)k", "\!\(\*SubscriptBox[\(\[Zeta]\), \(kin\)]\)"},
		Epilog->epilog,
		optPlot//Evaluate
		]
	]
PlotKineticSpectrum[kinSpectra_Dataset,opts:OptionsPattern[]]:=Module[{spec,data,epilog,
	Km,KM,optPlot},
	optPlot=FilterRules[{opts},Options[LogLogPlot]];
	{Km,KM}=OptionValue["KRange"];
	spec=kinSpectra[All,"KSpectrum"]//Normal;
	epilog=If[OptionValue["KSpectrumData"],
		data=Flatten[Log@kinSpectra[All,"KSpectrumData"]//Normal,1];
		{Black,Point[data]},
		{}];
	LogLogPlot[Through[spec[K]]//Evaluate,{K,Km,KM},
		PlotLegends->(StringTemplate["\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\)=``"]/@Normal[kinSpectra[All,"vw"]]),
		Frame->True,GridLines->Automatic,FrameLabel->{"\!\(\*SubscriptBox[\(R\), \(*\)]\)k", "\!\(\*SubscriptBox[\(\[Zeta]\), \(kin\)]\)"},
		Epilog->epilog,
		optPlot//Evaluate
		]
	]
AutoComplete[PlotKineticSpectrum];


(* ::Input::Initialization:: *)
ClearAll[PlotGWData]
Options[PlotGWData]={
	ScalingFunctions->{"Log","Log"},
	Joined->True,
	GridLines->Automatic,
	Frame->True,
	"LegendTemplate"->StringTemplate["vw=`vw`  \[Alpha]=`\[Alpha]`  \!\(\*SubscriptBox[\(\[Delta]\[Tau]\), \(fin\)]\)\!\(\*SubscriptBox[\(\[ScriptCapitalH]\), \(*\)]\)=`SourceDuration`   \!\(\*SubscriptBox[\(\[ScriptCapitalH]\), \(*\)]\)\!\(\*SubscriptBox[\(R\), \(*\)]\)=`HR`"],
	FrameLabel->{"\!\(\*SubscriptBox[\(R\), \(*\)]\)k","\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(GW\)]\) / \!\(\*SubscriptBox[\(\[ScriptCapitalT]\), \(GW\)]\)"},
	"gammas"->None
	};
PlotGWData[data_List,opts:OptionsPattern[$dependencies[PlotGWData]]]:=Module[{
dim=data//Dimensions,dataPlot,
	optLeg,gammas=OptionValue["gammas"],
	optPlot=FilterRules[{opts}~Join~Options[PlotGWData],Options[ListLogLogPlot]//Evaluate]},
	optLeg=PlotLegends->If[ListQ@gammas,
		Placed[StringTemplate["\[CapitalGamma]/\!\(\*SubscriptBox[\(\[ScriptCapitalH]\), \(*\)]\)=``"][Round[#]]&/@gammas,{Center,Bottom}],
	None
	];
	dataPlot=If[dim[[-1]]==3,Map[#[[;;2]]&,data,{-2}],data];
	ListPlot[dataPlot,optPlot,optLeg]
	]
PlotGWData[gw_Association,opts:OptionsPattern[$dependencies[PlotGWData]]]:=Module[{title,data=gw["KOmega"]//Normal},
title=lbl["plotGW"][Sequence@@Lookup[gw,{"vw","\[Alpha]","SourceDuration","HR"}]];
	PlotGWData[data,opts//Evaluate,PlotLabel->title]
	]
PlotGWData[dataset_Dataset,opts:OptionsPattern[$dependencies[PlotGWData]]]:=Module[{
legend=OptionValue["LegendTemplate"],data=dataset[All,"KOmega"]//Normal
	},
	legend=legend/@Normal[dataset];
	PlotGWData[data,opts//Evaluate,PlotLegends->legend]
	]
AutoComplete[PlotGWData];
$dependencies[PlotGWData]={PlotGWData,ListLogLogPlot};

(* default plotting options *)
$PlotOptions={Frame->True,FrameLabel->{"p"},GridLines->Automatic};

(* extract data from plot, in the format {{{x1,y1},{x2,y2},..},{..}} *)
DataFromPlot[plot_]:=Cases[Normal[plot],Line[p_]:>p,Infinity]


Options[PlotTiming]={
	ScalingFunctions->{"Log","Log"},
	Joined->True,
	GridLines->Automatic,
	Frame->True,
	"LegendTemplate"->StringTemplate["vw=`vw`  \[Alpha]=`\[Alpha]`  \!\(\*SubscriptBox[\(\[Delta]\[Tau]\), \(fin\)]\)\!\(\*SubscriptBox[\(\[ScriptCapitalH]\), \(*\)]\)=`SourceDuration`   \!\(\*SubscriptBox[\(\[ScriptCapitalH]\), \(*\)]\)\!\(\*SubscriptBox[\(R\), \(*\)]\)=`HR`"],
	FrameLabel->{"K", "t (s)"},
	"gammas"->None
	};
PlotTiming[KOmega_List,opts:OptionsPattern[]]:=Module[{KT,dim,optPlot},
	optPlot=FilterRules[{opts}~Join~Options[PlotTiming],Options[ListLogLogPlot]];
	dim=Dimensions[KOmega];
	KT=If[Length[dim]==2\[And]dim[[2]]==3,KOmega[[All,{1,3}]],KOmega]; (* option to pass (ki,\[CapitalOmega]i,Ti) data *)
	ListPlot[KT,optPlot//Evaluate]
	]
PlotTiming[GW_Association,opts:OptionsPattern[]]:=PlotTiming[GW["KTime"],opts//Evaluate]
PlotTiming[GW_Dataset,opts:OptionsPattern[]]:=PlotTiming[GW[All,"KTime"]//Normal,opts//Evaluate]
AutoComplete[PlotTiming];


(* ::Subsection::Closed:: *)
(*Physics*)


(* ::Input::Initialization:: *)
cs2=1/3;
cs=Sqrt[1.cs2];
num=<|HR->.001,\[Delta]\[Tau]H->1.,\[CapitalDelta]\[Tau]H->1.,K->1.,q->.1,\[ScriptQ]->.01,z->.3,r->10.,Qmn->1.,Pmn->1.|>;
\[CapitalGamma]meanAdiabIndex=1+cs^2; (* \[CapitalGamma]=Overscript[w, _]/Overscript[\[Rho], _] *)
Print["Auxiliary parameters set to:"];
Echo[cs2,"\!\(\*SubsuperscriptBox[\(c\), \(s\), \(2\)]\) \[Rule]"];
Echo[\[CapitalGamma]meanAdiabIndex,"\!\(\*SubscriptBox[\(\[CapitalGamma]\), \(ad\)]\) \[Rule]"];
Message[cs2::careful]


EddyTurnover::usage="EddyTurnOver[\!\(\*SubscriptBox[\(\[ScriptCapitalH]\), \(\[Star]\)]\)\!\(\*SubscriptBox[\(R\), \(\[Star]\)]\),\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(K\)]\)] gives the eddy turnover time \!\(\*SubscriptBox[\(\[Delta]\[Tau]\), \(nl\)]\)\!\(\*SubscriptBox[\(H\), \(*\)]\), an estimate of the onset of non-linearities, and thus of the source duration \!\(\*SubscriptBox[\(\[Delta]\[Tau]\), \(fin\)]\).";
EddyTurnover[HR_,\[CapitalOmega]K_]:=HR/Sqrt[\[CapitalOmega]K]


(* ::Text:: *)
(*2do: user-friendly way to modify these parameters.*)


(* ::Section:: *)
(*Sound waves amplitudes (f,l,A^2)*)


(* ::Text:: *)
(*Introduce z\[Congruent]Subscript[kT, n]*)


(* ::Input::Initialization:: *)
ClearAll[f,l,Asq]
f[vip_,z_,{\[Xi]m_,\[Xi]M_}]:=(4\[Pi])/z NIntegrate[vip[\[Xi]] Sin[z \[Xi]],{\[Xi],\[Xi]m,\[Xi]M},Method->"Oscillatory"]
l[\[Lambda]ip_,z_,{\[Xi]m_,\[Xi]M_}]:=(4\[Pi])/z NIntegrate[\[Lambda]ip[\[Xi]] Sin[z \[Xi]]\[Xi],{\[Xi],0.,\[Xi]M},Method->"Oscillatory"]
Asq[vip_,\[Lambda]ip_,z_]:=1/4 (f'[z]^2+cs2 l[z]^2)


(* ::Input::Initialization:: *)
ClearAll[InterpolateAmplitudes]
Options[InterpolateAmplitudes]={"zDomain"->{10.^-2,10.^3},"Nz"->3000};
InterpolateAmplitudes[vw_,\[Alpha]_,OptionsPattern[]]:=Module[{
	profile,vip,wip,\[Lambda]ip,shellDom,
	zs,fInt,lInt,
	domain=OptionValue["zDomain"]
	},
	(* compute bubble fluid profile *)
	profile=Profile[vw,\[Alpha],"Verbose"->False];
	$SSPrint["Bubble type \[Rule] ",profile["Type"]];
	{vip,wip,\[Lambda]ip,shellDom}=Lookup[profile,{"v","w","\[Lambda]","ShellDomain"}];
	f[z_]:=f[vip,z,shellDom];
	l[z_]:=l[\[Lambda]ip,z,shellDom];
	(* interpolate amplitudes *)
	$SSPrint[msg["interpolAmps"][vw,\[Alpha],cs]];zs=E^Subdivide[Sequence@@Log[domain],OptionValue["Nz"]];
	$SSPrint["f(z).."];
	fInt=Interpolation[Table[{z,f[z]},{z,zs}]]//EchoTiming;
	$SSPrint["l(z).."];
	lInt=Interpolation[Table[{z,l[z]},{z,zs}]]//EchoTiming;
	<|"f"->fInt,"l"->lInt|>
	]
AutoComplete[InterpolateAmplitudes];


(* ::Section:: *)
(*Kinetic spectrum (\[Zeta]_kin)*)


(* ::Text:: *)
(*Bubble lifetime distribution*)


(* ::Input::Initialization:: *)
(* bubble lifetime distribution *)
\[Nu][T_]=E^-T;  (* exponential nucleation *)
\[Nu][T_,"sim"]:=1/2 (T^2) (E^(-(T^3/6))) (* simultaneous nucleation *)


(* ::Text:: *)
(*Auxiliary functions*)


(* ::Input::Initialization:: *)
RBeta::usage="RBeta[\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\),\!\(\*SubscriptBox[\(c\), \(s\)]\)] gives \!\(\*SubscriptBox[\(R\), \(\[Star]\)]\)\!\(\*SubscriptBox[\(\[Beta]\), \(\[Star]\)]\)=(8\[Pi]\!\(\*SuperscriptBox[\()\), \(1/3\)]\)max(\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\),\!\(\*SubscriptBox[\(c\), \(s\)]\)).
RBeta[\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\)] gives \!\(\*SubscriptBox[\(R\), \(\[Star]\)]\)\!\(\*SubscriptBox[\(\[Beta]\), \(\[Star]\)]\)=(8\[Pi]\!\(\*SuperscriptBox[\()\), \(1/3\)]\)\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\).";
RBeta[vw_]:=(8\[Pi])^(1/3) vw;
RBeta[vw_,cs_]:=(8\[Pi])^(1/3) Max[vw,cs];
 (* approximation of R*\[Beta]* for fast wall expansions *)
EkOverRstarFun::usage="EkOverRstarFun[\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\),\!\(\*SubscriptBox[\(c\), \(s\)]\)] gives the factor Ek/(\!\(\*SubscriptBox[\(R\), \(\[Star]\)]\)\[Zeta]k)=1/(2\!\(\*SuperscriptBox[\(\[Pi]\), \(2\)]\)(R\[Beta]\!\(\*SuperscriptBox[\()\), \(6\)]\)).
EkOverRstarFun[\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\),\!\(\*SubscriptBox[\(c\), \(s\)]\),\"\[Zeta]kMax\"\[Rule]\!\(\*SubscriptBox[\(\[Zeta]\), \(max\)]\)] rescales by the \[Zeta]k maximum.";
Options[EkOverRstarFun]={"\[Zeta]kMax"->1.};
EkOverRstarFun[vw_,OptionsPattern[]]:=OptionValue["\[Zeta]kMax"]/(2\[Pi]^2 RBeta[vw]^6)(* such that Ekin/R* = Ek*/R* \[Zeta]kin *)
EkOverRstarFun[vw_,cs_,OptionsPattern[]]:=OptionValue["\[Zeta]kMax"]/(2\[Pi]^2 RBeta[vw,cs]^6)


(* ::Text:: *)
(*Spectrum \[Zeta]_kin(K)*)


(* ::Input::Initialization:: *)
KineticSpectrum[K_?NumericQ,R\[Beta]_?NumericQ,A2_,\[Nu]_:\[Nu]]:=K^2 NIntegrate[ \[Nu][T]T^6 A2[T K/R\[Beta]]//Evaluate,{T,0,\[Infinity](*R\[Beta]/Kzm,R\[Beta]/KzM*)},Method->"AdaptiveMonteCarlo"]


(* ::Text:: *)
(*\[ScriptCapitalC] square integral*)


\[ScriptCapitalC]Fun::usage="\[ScriptCapitalC]Fun[\[Zeta]k] gives the integral \[ScriptCapitalC] of the square kinetic energy spectrum: 16/15 \!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(\[Infinity]\)]\)\[DifferentialD]K \!\(\*SuperscriptBox[\(\[Zeta]\), \(2\)]\)(K)/\!\(\*SuperscriptBox[\(K\), \(2\)]\).";
\[ScriptCapitalC]Fun[\[Zeta]k_]:=16/15 NIntegrate[(\[Zeta]k[K]/K)^2,{K,0,\[Infinity]},Method->"GlobalAdaptive",MaxRecursion->\[Infinity],AccuracyGoal->\[Infinity],PrecisionGoal->4]


(* ::Subsection:: *)
(*Interpolation*)


ClearAll[InterpolateKineticSpectrum]
Options[InterpolateKineticSpectrum]={
	"KDomain"->{10.^-2,10^3},
	"NormalizePeak"->False,
	"ExtendKineticSpectrum"->True,
	"Verbose"->False,
	"Store"->False,
	"NK"->100
	};
InterpolateKineticSpectrum::usage="InterpolateKineticSpectrum[\[Zeta]k,\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\),\!\(\*SubscriptBox[\(c\), \(s\)]\)] interpolates the kinetic spectrum function \[Zeta]k, returning an Association.
InterpolateKineticSpectrum[{\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\),\[Alpha],\!\(\*SubscriptBox[\(c\), \(s\)]\)}] computes the square velocity amplitude \!\(\*SuperscriptBox[\(A\), \(2\)]\) for the given parameters.";
InterpolateKineticSpectrum[\[Zeta]KFun_,vw_?NumericQ,cs_?NumericQ,OptionsPattern[]]:=Module[{
	Km,KM,Ks,\[Zeta]kData,\[Zeta]kMax,Log\[Zeta]k,\[Zeta]k,\[ScriptCapitalK],\[CapitalOmega]K,EkOverRstar,\[ScriptCapitalC],KSpec,
	KDomain=OptionValue["KDomain"],N=OptionValue["NK"]
	},
	{Km,KM}=KDomain;
	Ks=LogSubdivide[Km,KM,N-1];
	\[Zeta]kData=Table[Check[{K,\[Zeta]KFun[K]},Nothing],{K,Ks}]//Quiet;
	(* peak and Subscript[E, k]/Subscript[R, *] *)
	\[Zeta]kMax=MaximalBy[\[Zeta]kData,Last]//First;
	If[OptionValue["Verbose"],
		Echo[\[Zeta]kMax,"{\!\(\*SuperscriptBox[\(k\), \(max\)]\),\!\(\*SubsuperscriptBox[\(\[Zeta]\), \(k\), \(max\)]\)} \[Rule]",MatrixForm]
		];
	(* normalize \[Zeta]k: \[Zeta]k,max\[Congruent]1 *)
	If[OptionValue["NormalizePeak"],
		\[Zeta]kData=MapAt[#/\[Zeta]kMax[[2]]&,\[Zeta]kData,{All,2}];
		EkOverRstar=EkOverRstarFun[vw,cs,"\[Zeta]kMax"->\[Zeta]kMax[[2]]],
		EkOverRstar=EkOverRstarFun[vw,cs]
		];
	(* extrapolate by power laws *)
	If[OptionValue["ExtendKineticSpectrum"],\[Zeta]kData=ExtendSpectrum[\[Zeta]kData]];
	Log\[Zeta]k=Interpolation[Log@\[Zeta]kData,InterpolationOrder->1];
	KDomain=E^Log\[Zeta]k["Domain"][[1]];
	\[Zeta]k=Function[K,E^Log\[Zeta]k[Log@K]//Evaluate];
	\[ScriptCapitalK]=NIntegrate[\[Zeta]k[K],{K,Km,KM},AccuracyGoal->\[Infinity],PrecisionGoal->3,MaxRecursion->\[Infinity]];
	\[CapitalOmega]K=EkOverRstar \[ScriptCapitalK];
	\[ScriptCapitalC]=\[ScriptCapitalC]Fun[\[Zeta]k];
	KSpec=<|
		"EkOverRstar"->EkOverRstar,
		"KSpectrumIntegral"->\[ScriptCapitalK],
		"KEnergyDensityFraction"->\[CapitalOmega]K,
		"\[ScriptCapitalC]"->\[ScriptCapitalC],
		"KDomain"->KDomain,
		"LogKSpectrum"->Log\[Zeta]k,
		"KSpectrum"->\[Zeta]k
		|>;
	Return[KSpec]
	]
InterpolateKineticSpectrum[{vw_?NumericQ,\[Alpha]_?NumericQ,cs_?NumericQ},opts:OptionsPattern[]]:=Module[{
	amps,fInt,lInt,A2Int,
	R\[Beta],\[Zeta]kFun,kinSpectrum,optsAmp,optsKin,store=OptionValue["Store"]
	},
	optsAmp=FilterRules[{opts},Options[InterpolateAmplitudes]]//Evaluate;
	optsKin=FilterRules[{opts},Options[InterpolateKineticSpectrum]]//Evaluate;
	(* interpolate amplitues and extract f, l, A^2 *)
	amps=InterpolateAmplitudes[vw,\[Alpha],optsAmp];
	{fInt,lInt}=Lookup[amps,{"f","l"}];
	A2Int[z_]=.25(fInt'[z]^2+cs^2 lInt[z]^2);
	R\[Beta]=RBeta[vw,cs];
	$SSPrint[msg["interpolKin"][]];
	\[Zeta]kFun[K_]:=KineticSpectrum[K,R\[Beta],A2Int];
	kinSpectrum=InterpolateKineticSpectrum[\[Zeta]kFun,vw,cs,optsKin];
	PrependTo[kinSpectrum,<|"BubbleType"->BubbleType[vw,\[Alpha]],"vw"->vw,"\[Alpha]"->\[Alpha],"cs"->cs,"R\[Beta]"->R\[Beta],"f"->fInt,"l"->lInt|>];
	If[MatchQ[store,_String],
		$SSPrint[StringTemplate["Appending kinetic spectrum to ``"][store]];
		AppendAssociation[store,kinSpectrum]];
	Return[kinSpectrum]
	]
AutoComplete[InterpolateKineticSpectrum];


(* ::Subsection:: *)
(*Get spectrum*)


(* ::Text:: *)
(*The following function *)
(*- checks if a file contains a kinetic spectrum with the given parameters (or the file does not exist)*)
(*- if it does, returns the corresponding spectrum (as association)*)
(*- if it doesn't, computes it through interpolation of f(z),l(z) -> A^2 -> \[Zeta]_kin(K),*)
(*	and stores it to the same file*)


(* ::Input::Initialization:: *)
GetKineticSpectrum::usage="GetKineticSpectrum[\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\),\[Alpha],\!\(\*SubscriptBox[\(c\), \(s\)]\)] gives an Association with the kinetic spectrum and relative parameters. If the corresponding spectrum does not alread exist in file, it is computed and stored.
GetKineticSpectrum[\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\),\[Alpha],\!\(\*SubscriptBox[\(c\), \(s\)]\),\"File\"\[Rule]file] sets the lookup file (defaults to \"$KineticFile\").";
Options[GetKineticSpectrum]={"File"->$KineticFile,"Verbose"->True};
GetKineticSpectrum[vw_?NumericQ,\[Alpha]_?NumericQ,cs_?NumericQ,opts:OptionsPattern[]]:=Module[{
	kinSpectrum,nSpectra,fileExist=True,optsKin,
file=OptionValue["File"],verbose=OptionValue["Verbose"]
	},
	optsKin=FilterRules[{opts},Options[InterpolateKineticSpectrum]]//Evaluate;
	kinSpectrum=Check[Get[file][Select[#vw==vw\[And]#\[Alpha]==\[Alpha]\[And]#cs==cs&]]//Normal,fileExist=False];
	nSpectra=Length@kinSpectrum;
	Which[
		nSpectra==0\[Or]Not[fileExist],
	$SSPrint[msg["noKin"][vw,\[Alpha],cs,file],fileExist\[And]verbose];
		kinSpectrum=InterpolateKineticSpectrum[{vw,\[Alpha],cs},"Store"->file,optsKin//Evaluate];
	$SSPrint[msg["computedKin"][kinSpectrum],verbose],
		nSpectra>0,If[nSpectra>1,$SSPrint[msg["multiKin"][nSpectra,vw,\[Alpha],cs,file]]];
		kinSpectrum=First[kinSpectrum];
	$SSPrint[msg["loadedKin"][kinSpectrum],verbose]
		];
	(* defining \[Zeta]k=\[ExponentialE]^Log\[Zeta]k[Log[K]] only when GetKineticSpectrum is called (for efficiency) *)
	kinSpectrum["KSpectrum"]=Function[K,E^kinSpectrum["LogKSpectrum"][Log[K]]//Evaluate];
	kinSpectrum
	]
$dependencies[GetKineticSpectrum]={InterpolateKineticSpectrum};
AutoComplete[GetKineticSpectrum,$dependencies[GetKineticSpectrum]];


(* ::Section::Closed:: *)
(*Time integrals (\[CapitalDelta] functions)*)


(* ::Input::Initialization:: *)
(* quantities are H-normalized: Subscript[H, *]Subscript[\[Tau], *]=1 \[Implies] \[Tau]H\[Congruent]Subscript[H, *]\[Tau], q\[Congruent]p/Subscript[H, *] *)
ClearAll[pmn,\[CapitalDelta]mn,\[CapitalDelta],\[CapitalDelta]DefH];
pmn::usage="pmn[{m,n},{k,p,\[ScriptP],\!\(\*SubscriptBox[\(c\), \(s\)]\)}] gives the pseudo-momentum \!\(\*SubscriptBox[\(p\), \(mn\)]\)\[Congruent](p+m\[ScriptP])\!\(\*SubscriptBox[\(c\), \(s\)]\)+nk (eq. 53-54 of RPC).";
pmn[{m_Integer,n_Integer},{k_,p_,\[ScriptP]_,cs_}]:=(p+m \[ScriptP])cs+n k
\[CapitalDelta]Ci[\[Delta]\[Theta]_,q_]:=CosIntegral[q(1+\[Delta]\[Theta])]-CosIntegral[q]
\[CapitalDelta]Si[\[Delta]\[Theta]_,q_]:=SinIntegral[q(1+\[Delta]\[Theta])]-SinIntegral[q]
\[CapitalDelta]mn::usage="\[CapitalDelta]mn[\[Delta]\[Tau]H,\!\(\*SubscriptBox[\(p\), \(mn\)]\)] gives the \!\(\*SubscriptBox[\(\[CapitalDelta]\), \(mn\)]\) time integral (eq. 56 of RPC).
\[CapitalDelta]mn[\[Delta]\[Tau]H,{m,n},{k,p,\[ScriptP],\!\(\*SubscriptBox[\(c\), \(s\)]\)}] computes \!\(\*SubscriptBox[\(p\), \(mn\)]\) from the given momenta and speed of sound.";
\[CapitalDelta]mn[\[Delta]\[Tau]H_,qmn_]:=1/4 Chop[\[CapitalDelta]Ci[\[Delta]\[Tau]H,qmn]^2+\[CapitalDelta]Si[\[Delta]\[Tau]H,qmn]^2]; (* qmn<0 gives 0.\[ImaginaryI] *)
\[CapitalDelta]mn[\[Delta]\[Tau]H_,{m_Integer,n_Integer},pmnArgs_List]:=\[CapitalDelta]mn[\[Delta]\[Tau]H,pmn[{m,n},pmnArgs]];
\[CapitalDelta]::usage="\[CapitalDelta][\[Delta]\[Tau]H,{k,p,\[ScriptP],\!\(\*SubscriptBox[\(c\), \(s\)]\)}] gives the \[CapitalDelta] time integral (eq. 55 of RPC).";
\[CapitalDelta][\[Delta]\[Tau]H_,{\[ScriptK]_,q_,\[ScriptQ]_,cs_}]:=Sum[\[CapitalDelta]mn[\[Delta]\[Tau]H,{m,n},{\[ScriptK],q,\[ScriptQ],cs}],{m,{-1,1}},{n,{-1,1}}]//Evaluate


(* flat-Universe expression for comparison *)
\[CapitalDelta]mnFlatOverH2[\[Delta]\[Tau]H_,qmn_]:=(1-Cos[qmn \[Delta]\[Tau]H])/(2qmn^2)//Evaluate;


(* ::Subsection:: *)
(*\[CapitalDelta]0*)


\[CapitalDelta]0[\[Delta]\[Tau]H_]:=Log[1+\[Delta]\[Tau]H]^2/2


(* ::Section:: *)
(*Gravitational wave spectrum (\[CapitalOmega]_GW)*)


(* ::Subsection:: *)
(*Compute spectrum (fixed K)*)


(* ::Subsubsection:: *)
(*\[CapitalOmega]_GW*)


(* ::Input::Initialization:: *)
ClearAll[OmegaGWPref]
\[ScriptCapitalT]GW[gstar_]:=1.6 10^-5 (100/gstar)^(1/3) (* EW approximation for the redshift factor *)
\[ScriptCapitalT]GW[None]=1; (* default to no redshift *)
Options[OmegaGWPref]={"DensityFraction"->1,"RelativisticDOF"->None};
OmegaGWPref[\[CapitalGamma]meanAdiabIndex_,EkOverRstar_,OptionsPattern[]]:= Module[{
	\[Rho]Ratio=OptionValue["DensityFraction"],gstar=OptionValue["RelativisticDOF"]
	},
	3\[CapitalGamma]meanAdiabIndex^2  EkOverRstar^2 \[Rho]Ratio^2 \[ScriptCapitalT]GW[gstar]
	]
(* 3(Overscript[w, _]/Overscript[\[Rho], _])^2(Subscript[E, k*]/Subscript[R, *])^2 Subscript[T, GW](Subscript[g, *]) *)


(* H*R*-normed GW spectrum *)
{Pm,PM}={0,\[Infinity]};
z[K_,P_,\[ScriptCapitalP]_]:=(K^2+P^2-\[ScriptCapitalP]^2)/(2 P K);
OmegaGWIntegrand::usage="OmegaGWIntegrand[\[Zeta]k,\[CapitalDelta]N,{K,P,z}] gives the GW spectrum integrand for a given kinetic spectrum \[Zeta]k[K], and a numerical time integral \[CapitalDelta]N[p,\[ScriptP]].";
OmegaGWIntegrand[\[Zeta]k_,\[CapitalDelta]N_,{K_?NumericQ,logP_?NumericQ,log\[ScriptCapitalP]_?NumericQ}]:=Module[{P=E^logP,\[ScriptCapitalP]=E^log\[ScriptCapitalP]},
	(1-z[K,P,\[ScriptCapitalP]]^2)^2/\[ScriptCapitalP]^2 \[Zeta]k[\[ScriptCapitalP]]\[CapitalDelta]N[P,\[ScriptCapitalP]]
	]
OmegaGWIntegrandP::usage="OmegaGWIntegrandP[\[Zeta]k,\[CapitalDelta]N,{K,z}] gives the GW spectrum integrated over P but not z, for a given kinetic spectrum \[Zeta]k[K], and a numerical time integral \[CapitalDelta]N[p,\[ScriptP]].";
Options[OmegaGWIntegrandP]={AccuracyGoal->\[Infinity],PrecisionGoal->2,MaxRecursion->\[Infinity]};
OmegaGWIntegrandP[\[Zeta]k_,\[CapitalDelta]N_,{K_?NumericQ,logP_?NumericQ},opts:OptionsPattern[]]:=NIntegrate[
	OmegaGWIntegrand[\[Zeta]k,\[CapitalDelta]N,{K,logP,log\[ScriptCapitalP]}],
	{log\[ScriptCapitalP],LogN[Abs[K-E^logP]],LogN[K+E^logP]},
	opts//Evaluate,
	AccuracyGoal->\[Infinity],
	PrecisionGoal->2,
	MaxRecursion->\[Infinity]
	]
OmegaGW::usage="OmegaGW[K,{\!\(\*SubscriptBox[\(c\), \(s\)]\),\[Delta]\[Tau]H,HR},kineticSpectrum] gives the GW spectrum at given K\[Congruent]R*k.
OmegaGW[K,{\!\(\*SubscriptBox[\(c\), \(s\)]\),\[Delta]\[Tau]H,HR},{\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\),\[Alpha]}] computes/loads the kinetic spectrum for given vw, \[Alpha] and cs parameters. 
";
Options[OmegaGW]={
	AccuracyGoal->\[Infinity],
	PrecisionGoal->2,
	MaxRecursion->\[Infinity],
	"RelativisticDOF"->None,
	"Verbose"->False
	};
OmegaGW[K_?NumericQ,{cs_?NumericQ,\[Delta]\[Tau]H_?NumericQ,HR_?NumericQ},kinSpectrum:(_Association|_List),opts:OptionsPattern[$dependencies[OmegaGW]]]:=Module[{
	pref,\[Zeta]k,EkOverRstar,\[CapitalDelta]N,optsPref,optsInt
	},
	optsPref=FilterRules[{opts}~Join~Options[OmegaGW],Options[OmegaGWPref]];
	optsInt=FilterRules[{opts}~Join~Options[OmegaGW],Options[NIntegrate]];
	{\[Zeta]k,EkOverRstar}=Lookup[kinSpectrum,{"KSpectrum","EkOverRstar"}];
	pref=OmegaGWPref[\[CapitalGamma]meanAdiabIndex,EkOverRstar,optsPref//Evaluate];
	\[CapitalDelta]N[p_,\[ScriptP]_]:=Re@\[CapitalDelta][\[Delta]\[Tau]H,{K/HR,p/HR,\[ScriptP]/HR,cs}]//Simplify//Evaluate;
	pref K^2
	NIntegrate[E^(2logP)\[Zeta]k[E^logP]
		OmegaGWIntegrandP[\[Zeta]k,\[CapitalDelta]N,{K,logP},optsInt//Evaluate],
		{logP,-\[Infinity],\[Infinity]},
		optsInt//Evaluate
		]
	]
(* alternative implementation: explicit vw, \[Alpha] dependence *)
OmegaGW[K_?NumericQ,{cs_?NumericQ,\[Delta]\[Tau]H_?NumericQ,HR_?NumericQ},{vw_?NumericQ,\[Alpha]_?NumericQ},opts:OptionsPattern[]]:=Module[{
	pref,kinSpectrum,\[Zeta]k,EkOverRstar,\[CapitalDelta]N,
	optsKin,optsPref,optsInt
	},
	optsKin=FilterRules[{opts},Options[GetKineticSpectrum]]//Evaluate;
	kinSpectrum=GetKineticSpectrum[vw,\[Alpha],cs,optsKin//Evaluate];
	OmegaGW[K,{cs,\[Delta]\[Tau]H,HR},kinSpectrum,opts//Evaluate]
	]
$dependencies[OmegaGW]={OmegaGW,OmegaGWPref,GetKineticSpectrum};
AutoComplete[OmegaGW,$dependencies[OmegaGW]];


(* ::Subsubsection::Closed:: *)
(*Old implementation (P,z integral)*)


(*(* H*R*-normed GW spectrum *)
{Pm,PM}={0,\[Infinity]};
OmegaGWIntegrand::usage="OmegaGWIntegrand[\[Zeta]k,\[CapitalDelta]N,{K,P,z}] gives the GW spectrum integrand for a given kinetic spectrum \[Zeta]k[K], and a numerical time integral \[CapitalDelta]N[p,\[ScriptP]].";
OmegaGWIntegrand[\[Zeta]k_,\[CapitalDelta]N_,{K_?NumericQ,P_?NumericQ,z_?NumericQ}]:=With[{\[ScriptCapitalP]=Sqrt[P^2+K^2-2P K z]},(1-z^2)^2 P^2/\[ScriptCapitalP]^4 \[Zeta]k[P]\[Zeta]k[\[ScriptCapitalP]]\[CapitalDelta]N[P,\[ScriptCapitalP]]]
(*OmegaGW[\[Delta]\[Tau]H_,K_,HR_,r_]:=OmegaGWPref[\[CapitalGamma]meanAdiabIndex,r,EkOverRstar] K^3NIntegrate[OmegaGWIntegrand[\[Delta]\[Tau]H,{K,P,z,cs},HR],{z,-1,1},{P,Pm,PM},Method->"AdaptiveMonteCarlo"]*)
OmegaGWIntegrandP::usage="OmegaGWIntegrandP[\[Zeta]k,\[CapitalDelta]N,{K,z}] gives the GW spectrum integrated over P but not z, for a given kinetic spectrum \[Zeta]k[K], and a numerical time integral \[CapitalDelta]N[p,\[ScriptP]].";
Options[OmegaGWIntegrandP]={AccuracyGoal->\[Infinity],PrecisionGoal->3,MaxRecursion->3};
OmegaGWIntegrandP[\[Zeta]k_,\[CapitalDelta]N_,{K_?NumericQ,z_?NumericQ},opts:OptionsPattern[]]:=NIntegrate[OmegaGWIntegrand[\[Zeta]k,\[CapitalDelta]N,{K,P,z}],{P,Pm,PM},
	Method->"GlobalAdaptive",
	opts//Evaluate
	]
OmegaGW::usage="OmegaGW[K,{\!\(\*SubscriptBox[\(c\), \(s\)]\),\[Delta]\[Tau]H,HR},kineticSpectrum] gives the GW spectrum at given K\[Congruent]R*k.
OmegaGW[K,{\!\(\*SubscriptBox[\(c\), \(s\)]\),\[Delta]\[Tau]H,HR},{\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\),\[Alpha]}] computes/loads the kinetic spectrum for given vw, \[Alpha] and cs parameters. 
";
Options[OmegaGW]={AccuracyGoal->\[Infinity],PrecisionGoal->3,MaxRecursion->3,"RelativisticDOF"->None};
OmegaGW[K_?NumericQ,{cs_?NumericQ,\[Delta]\[Tau]H_?NumericQ,HR_?NumericQ},kinSpectrum:(_Association|_List),opts:OptionsPattern[]]:=Module[{
	pref,\[Zeta]k,EkOverRstar,\[CapitalDelta]N,optsPref,optsInt
	},
	optsPref=FilterRules[{opts},Options[OmegaGWPref]];
	optsInt=FilterRules[{opts},Options[NIntegrate]];
	{\[Zeta]k,EkOverRstar}=Lookup[kinSpectrum,{"KSpectrum","EkOverRstar"}];
	pref=OmegaGWPref[\[CapitalGamma]meanAdiabIndex,EkOverRstar,optsPref//Evaluate];
	\[CapitalDelta]N[p_,\[ScriptP]_]:=Re@\[CapitalDelta][\[Delta]\[Tau]H,{K/HR,p/HR,\[ScriptP]/HR,cs}]//Simplify//Evaluate;
	pref K^3 
	NIntegrate[OmegaGWIntegrandP[\[Zeta]k,\[CapitalDelta]N,{K,z},optsInt],{z,-.999,.999},
		Method->"LocalAdaptive",
		optsInt//Evaluate
		]
	]
(* alternative implementation: explicit vw, \[Alpha] dependence *)
OmegaGW[K_?NumericQ,{cs_?NumericQ,\[Delta]\[Tau]H_?NumericQ,HR_?NumericQ},{vw_?NumericQ,\[Alpha]_?NumericQ},opts:OptionsPattern[]]:=Module[{
	pref,kinSpectrum,\[Zeta]k,EkOverRstar,\[CapitalDelta]N,
	optsKin,optsPref,optsInt
	},
	optsKin=FilterRules[{opts},Options[GetKineticSpectrum]]//Evaluate;
	kinSpectrum=GetKineticSpectrum[vw,\[Alpha],cs,optsKin//Evaluate];
	OmegaGW[K,{cs,\[Delta]\[Tau]H,HR},kinSpectrum,opts]
	]
$dependencies[OmegaGW]={OmegaGWPref,GetKineticSpectrum};
AutoComplete[OmegaGW,$dependencies[OmegaGW]];*)


(* ::Subsubsection::Closed:: *)
(*\[Delta]-approximations (HH19)*)


OmegaGWPrefHighK[\[CapitalGamma]meanAdiabIndex_,EkOverRstar_,cs_,\[Delta]\[Tau]H_,HR_,opts:OptionsPattern[]]:=\[Pi]/2 HR/cs \[CapitalUpsilon][\[Delta]\[Tau]H]OmegaGWPref[\[CapitalGamma]meanAdiabIndex,EkOverRstar,opts//Evaluate]


OmegaGWHighK::usage="OmegaGWHighK[k SubStar[R],{\!\(\*SubscriptBox[\(c\), \(s\)]\),\[Delta]\[Tau]H,HR},kineticSpectrum] gives the GW spectrum at given K\[Congruent]R*k.
OmegaGWHighK[k SubStar[R],{\!\(\*SubscriptBox[\(c\), \(s\)]\),\[Delta]\[Tau]H,HR},{\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\),\[Alpha]}] computes/loads the kinetic spectrum for given vw, \[Alpha] and cs parameters. 
";
Options[OmegaGWHighK]={
	AccuracyGoal->\[Infinity],
	PrecisionGoal->2,
	MaxRecursion->\[Infinity],
	"RelativisticDOF"->None,
	"Verbose"->False
	};
OmegaGWHighK[K_?NumericQ,{cs_?NumericQ,\[Delta]\[Tau]H_?NumericQ,HR_?NumericQ},kinSpectrum_Association,opts:OptionsPattern[]]:=Module[{
	pref,\[Zeta]k,EkOverRstar,optsPref,optsInt,
	Pm=K/2(1/cs-1),
	Pp=K/2(1/cs+1),
	z,\[ScriptCapitalP]
	},
	optsPref=FilterRules[{opts},Options[OmegaGWPref]];
	optsInt=FilterRules[{opts},Options[NIntegrate]];
	{\[Zeta]k,EkOverRstar}=Lookup[kinSpectrum,{"KSpectrum","EkOverRstar"}];
	z[P_]:=1/cs+K/(2P)(1-1/cs^2);
	\[ScriptCapitalP][P_]:=K/cs-P;
	pref=OmegaGWPrefHighK[\[CapitalGamma]meanAdiabIndex,EkOverRstar,cs,\[Delta]\[Tau]H,HR,optsPref//Evaluate];
	pref K^2 
	NIntegrate[(1-z[P]^2)^2 P/\[ScriptCapitalP][P]^3 \[Zeta]k[P]\[Zeta]k[\[ScriptCapitalP][P]],
		{P,Pm,Pp},
		optsInt//Evaluate
		]
	]
OmegaGWHighK[K_?NumericQ,{cs_?NumericQ,\[Delta]\[Tau]H_?NumericQ,HR_?NumericQ},{vw_?NumericQ,\[Alpha]_?NumericQ},opts:OptionsPattern[]]:=Module[{
	pref,kinSpectrum,\[Zeta]k,EkOverRstar,
	optsKin,optsPref,optsInt
	},
	optsKin=FilterRules[{opts},Options[GetKineticSpectrum]]//Evaluate;
	kinSpectrum=GetKineticSpectrum[vw,\[Alpha],cs,optsKin//Evaluate];
	OmegaGWHighK[K,{cs,\[Delta]\[Tau]H,HR},kinSpectrum,opts//Evaluate]
	]
$dependencies[OmegaGWHighK]={OmegaGWHighK,OmegaGWPref,GetKineticSpectrum};
AutoComplete[OmegaGWHighK,$dependencies[OmegaGWHighK]];


OmegaGWPrefLowK[\[CapitalGamma]meanAdiabIndex_,EkOverRstar_,\[ScriptCapitalC]_,\[Delta]\[Tau]H_,opts:OptionsPattern[]]:=\[ScriptCapitalC] \[CapitalDelta]0[\[Delta]\[Tau]H] OmegaGWPref[\[CapitalGamma]meanAdiabIndex,EkOverRstar,opts//Evaluate]


OmegaGWLowK::usage="OmegaGWLowK[k SubStar[R],\[Delta]\[Tau]H,,kineticSpectrum] gives the GW spectrum at given K\[Congruent]R*k.
OmegaGWLowK[k SubStar[R],{\!\(\*SubscriptBox[\(c\), \(s\)]\),\[Delta]\[Tau]H},{\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\),\[Alpha]}] computes/loads the kinetic spectrum for given vw, \[Alpha] and cs parameters. 
";
Options[OmegaGWLowK]={
	"RelativisticDOF"->None,
	"Verbose"->False
	};
OmegaGWLowK[K_?NumericQ,\[Delta]\[Tau]H_?NumericQ,kinSpectrum_Association,opts:OptionsPattern[]]:=Module[{
	pref,\[Zeta]k,EkOverRstar,\[ScriptCapitalC],optsPref
	},
	optsPref=FilterRules[Join[{opts},Options[OmegaGWLowK]],Options[OmegaGWPref]];
	{\[Zeta]k,EkOverRstar,\[ScriptCapitalC]}=Lookup[kinSpectrum,{"KSpectrum","EkOverRstar","\[ScriptCapitalC]"}];
	pref=OmegaGWPrefLowK[\[CapitalGamma]meanAdiabIndex,EkOverRstar,\[ScriptCapitalC],\[Delta]\[Tau]H,optsPref//Evaluate];
	pref K^3
	]
OmegaGWLowK[K_?NumericQ,{cs_?NumericQ,\[Delta]\[Tau]H_?NumericQ},{vw_?NumericQ,\[Alpha]_?NumericQ},opts:OptionsPattern[]]:=Module[{
	pref,kinSpectrum,\[Zeta]k,EkOverRstar,
	optsKin,optsPref
	},
	optsKin=FilterRules[Join[{opts},Options[OmegaGWLowK]],Options[GetKineticSpectrum]];
	kinSpectrum=GetKineticSpectrum[vw,\[Alpha],cs,optsKin//Evaluate];
	OmegaGWLowK[K,\[Delta]\[Tau]H,kinSpectrum,opts//Evaluate]
	]
$dependencies[OmegaGWLowK]={OmegaGWLowK,OmegaGWPref,GetKineticSpectrum};
AutoComplete[OmegaGWLowK,$dependencies[OmegaGWLowK]];


(* ::Subsection:: *)
(*Compute spectrum (K range)*)


ComputeOmegaGW::usage="ComputeOmegaGW[Function] computes the GW spectra for a given function of the momentum K\[Congruent]Rk.";
(*ComputeOmegaGW[{\!\(\*SubscriptBox[\(c\), \(s\)]\),d\[Tau]H,HR},{\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\),\[Alpha]}] computes the spectrum via OmegaGW[K,{\!\(\*SubscriptBox[\(c\), \(s\)]\),d\[Tau]H,HR},{\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\),\[Alpha]}], returning an Association.";*)
ClearAll[ComputeOmegaGW]
Options[ComputeOmegaGW]={
	"KRange"->{10.^-3,2 10.^3},
	"N"->30,"Plot"->False,
	"DynamicPlot"->True,
	"Verbose"->True,
	"Store"->False,
	"Metadata"->None
	};
ComputeOmegaGW[\[CapitalOmega]Function_,OptionsPattern[]]:=Module[{N=OptionValue["N"],Ncomp,Ks,K\[CapitalOmega]T,Km,KM,verbose=OptionValue["Verbose"]},
	{Km,KM}=OptionValue["KRange"];
	Ks=LogSubdivide[Km,KM,N-1];
	$KOmega={}; (* global variable for dynamic plot *)
	$KOmegaFail={}; (* global variable for data which failed to compute *)
	If[OptionValue["DynamicPlot"],Print@Dynamic[PlotGWData[Abs@$KOmega[[All,;;2]]]]];
	Do[K\[CapitalOmega]T=Flatten[{K,AbsoluteTiming[\[CapitalOmega]Function[K]]}][[{1,3,2}]];
		If[MatchQ[K\[CapitalOmega]T,{_?NumericQ,_?NumericQ,_?NumericQ}],
			AppendTo[$KOmega,K\[CapitalOmega]T],
			AppendTo[$KOmegaFail,K\[CapitalOmega]T]],
		{K,Ks}
		]//EchoTiming;
	Ncomp=Length[$KOmega];
	$SSPrint[msg["nComputed"][Ncomp,N,PercentForm[1.Ncomp/N]],verbose];
	If[OptionValue["Plot"],Print@PlotGWData[$KOmega[[All,;;2]]]];
	$KOmega
	]
$GWKeyOrder={"vw","\[Alpha]","cs","SourceDuration","HR","R\[Beta]","EkOverRstar","KSpectrumIntegral","KEnergyDensityFraction","\[ScriptCapitalC]","EddyTurnover","KOmega","KTime","Metadata"};
ComputeOmegaGW[{cs_,\[Delta]\[Tau]H_,HR_},{vw_,\[Alpha]_},opts:OptionsPattern[$dependencies[ComputeOmegaGW]]]:=Module[{
	K\[CapitalOmega],KTime,KSpectrum,GWSpectrum,\[Delta]\[Tau]Hnl,store,
	optsKin,optsGW,optsComp},
	optsKin=FilterOptions[opts,GetKineticSpectrum];
	optsGW=FilterOptions[opts,OmegaGW];
	optsComp=FilterOptions[opts,ComputeOmegaGW];
	store=OptionValue["Store"]/.True->$GWFile;
	(* extract kinetic spectrum (fixed vw and \[Alpha]) *)
	KSpectrum=GetKineticSpectrum[vw,\[Alpha],cs,optsKin//Evaluate];
	(* compute spectrum *)
	K\[CapitalOmega]=ComputeOmegaGW[Function[K,OmegaGW[K,{cs,\[Delta]\[Tau]H,HR},KSpectrum,optsGW//Evaluate]],
		optsComp//Evaluate];
	{K\[CapitalOmega],KTime}={K\[CapitalOmega][[All,;;2]],K\[CapitalOmega][[All,{1,3}]]};
	(* create spectrum object *)
	GWSpectrum=Select[KSpectrum,NumericQ];
	\[Delta]\[Tau]Hnl=EddyTurnover[HR,GWSpectrum["KEnergyDensityFraction"]];
	AppendTo[GWSpectrum,<|
		"SourceDuration"->\[Delta]\[Tau]H,
		"HR"->HR,
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
$dependencies[ComputeOmegaGW]={ComputeOmegaGW,OmegaGW,OmegaGWPref,GetKineticSpectrum};
AutoComplete[ComputeOmegaGW,$dependencies[ComputeOmegaGW]];


(* ::Section:: *)
(*Overall prefactor \[CapitalUpsilon](\[Delta]\[Tau]_fin \[ScriptCapitalH])*)


\[CapitalUpsilon]::usage="\[CapitalUpsilon][\!\(\*SubscriptBox[\(\[Delta]\[Tau]\), \(fin\)]\)\!\(\*SubscriptBox[\(\[ScriptCapitalH]\), \(*\)]\)] gives the suppression factor 1-1/\!\(\*SubscriptBox[\(\[Delta]\), \(fin\)]\)\!\(\*SubscriptBox[\(\[Tau]\[ScriptCapitalH]\), \(*\)]\)";
\[CapitalUpsilon][\[Delta]\[Tau]H_]:=\[Delta]\[Tau]H/(1+\[Delta]\[Tau]H);


(* ::Section:: *)
(*End*)


End[] (* Private` *)


(* Make SecondPackage symbols accessible in the user's session *)
If[! MemberQ[$ContextPath, "Bubble`"],
    AppendTo[$ContextPath, "Bubble`"]
	]


EndPackage[]
