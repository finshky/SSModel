(* ::Package:: *)

(* ::Title:: *)
(*Fluid profiles for*)
(*First-Order Phase Transitions Bubbles*)


(* ::Section:: *)
(*Initialize*)


BeginPackage["Bubble`"]


(* ::Subsection::Closed:: *)
(*Public functions*)


\[Gamma]sq;
\[Mu];
wRatio;
wShock;
wMinus;
\[Lambda]FromW;
\[Lambda]FromWSol;
shock;
ChapmanJouget;


BubbleType;
ProfileFromTargetAlpha;
Detonation;
Deflagration;
Hybrid;
Profile;


PlotProfile;
ManipulateVelocityProfile;


(* utilities *)
InterpolationData;
ExtendProfile;


(* ::Subsection::Closed:: *)
(*Private*)


Begin["`Private`"]


(* ::Section::Closed:: *)
(*Equations*)


(* ::Subsection::Closed:: *)
(*Utilities*)


(* 
https://resources.wolframcloud.com/FunctionRepository/resources/AddCodeCompletion
https://mathematica.stackexchange.com/questions/56984
*)
addCodeCompletion[function_String][args___]:=Module[
	{processed},
	processed={args}/.{
		None -> 0,
		"AbsoluteFileName" -> 2,
		"RelativeFileName" -> 3,
		"Color" -> 4,
		"PackageName" -> 7,
		"DirectoryName" -> 8,
		"InterpreterType" -> 9
	};
	(FE`Evaluate[FEPrivate`AddSpecialArgCompletion[#1]]&)[function -> processed]
];


SetAttributes[AutoComplete,HoldFirst];
AutoComplete[sym_]:=Module[{name,keys,keySeq,min,max},
	Switch[sym,
		_Symbol,
			name=ToString@sym;
			keys=Keys@Options@sym;
			(* count min/max number of arguments *)
			{min,max}=MinMax[Length/@List@@ReleaseHold[DownValues[sym][[All,1]]]] - 1/.{\[Infinity],-\[Infinity]}->{0,0};
			keySeq=Sequence@@Join[Table[0,min],ConstantArray[keys,Length@keys+(max-min)]],
		_Association,
			name=ToString@HoldForm@sym;
			keySeq=Keys@sym,
		_Transition|_ActionFunction,
			name=ToString@HoldForm@sym;
			keySeq=Keys@sym[Association]
		];
	addCodeCompletion[name][keySeq]
	]
(* autocomplete options from multiple functions *)
AutoComplete[sym_Symbol,functions_List]:=Module[
	{name,keys,keySeq},
	name=ToString@sym;
	keys=Keys[Union@@Options/@functions];
	keySeq=Sequence@@ConstantArray[keys,Length@keys];
	addCodeCompletion[name][keySeq]
	]


(* ::Subsection:: *)
(*Physics utilities*)


(* ::Input::Initialization:: *)
\[Gamma]sq[v_]:=1/(1-v^2) (* square Lorentz factor *)
\[Mu][v_,\[Xi]w_]:=(\[Xi]w-v)/(1-\[Xi]w v) (* Lorentz transformation of v from wall to universe frame *)
wRatio[v_,\[Xi]w_]:=\[Xi]w/v \[Gamma]sq[\[Xi]w]/\[Gamma]sq[v] (* enthalpy ratio w+/w- *)
wShock[\[Xi]_]:=(9\[Xi]^2-1)/(3(1-\[Xi]^2))
wMinus[wPlus_,wRatio_]:=wPlus/wRatio
(* energy perturbations from enthalpy. NB: \[Alpha] only contributes inside the bubble, where \[Xi]<\[Xi]w *)
\[Lambda]FromW[w_,\[Alpha]_,cs2_]:=1/(1+cs2) (w-1-\[Alpha])
\[Lambda]FromWSol[wSol_,\[Xi]w_,\[Alpha]_,cs2_]:=Module[{\[Xi]s=wSol["Coordinates"]//First},
	Interpolation[{\[Xi]s,\[Lambda]FromW[wSol["ValuesOnGrid"],\[Alpha] UnitStep[\[Xi]w-\[Xi]s],cs2]}\[Transpose],InterpolationOrder->1]
	]
shock[c_,\[Xi]_]:=v/.Solve[\[Xi] \[Mu][v,\[Xi]]==c^2,v]//First//Evaluate
ChapmanJouget::usage="ChapmanJouget[\[Alpha]] computes the Chapman-Jouget velocity froma given transition strength \[Alpha], in the Bag EoS.";
ChapmanJouget[\[Alpha]p_]:=(1+Sqrt[\[Alpha]p(2+3 \[Alpha]p)])/(Sqrt[3](1+\[Alpha]p)) (* Chapman-Jouget velocity, in the Bag EoS! *)


(* ::Subsection::Closed:: *)
(*Bubble characterization*)


(* Bubble expansion type *)
BubbleTypeRules={"Def"->Deflagration,"Det"->Detonation,"Hyb"->Hybrid};
Options[BubbleType]={"cs2"->cs2,"Input"->"\[Alpha]"};
BubbleType[\[Xi]w_,in_,OptionsPattern[{BubbleType}]]:=Module[{vCJ,cs=Sqrt[OptionValue["cs2"]]},
	If[OptionValue["Input"]==="\[Alpha]",
		vCJ=ChapmanJouget[in];
		Which[
		\[Xi]w<cs,"Def",
		\[Xi]w>vCJ,"Det",
		True,"Hyb" (* actually one should consider the possibility of no solution *)
		],
	Which[
		\[Xi]w<cs\[And]in<\[Xi]w,"Def",
		\[Xi]w>cs\[And]in<\[Mu][cs,\[Xi]w],"Det",
		\[Xi]w>cs\[And]in>shock[cs,\[Xi]w]\[And]in<\[Xi]w,"Hyb", (* actually one should consider the possibility of no solution *)
		True,None
		]
	]]


(* ::Subsection::Closed:: *)
(*Interpolation*)


(* extract data from InterpolatingFunction *)
InterpolationData[f_InterpolatingFunction]:={First[f["Coordinates"]],f["ValuesOnGrid"]}\[Transpose]
Options[InterpolationUnion]={InterpolationOrder->1,"LeftMultiplier"->1,"RightMultiplier"->1};
InterpolationUnion[left_InterpolatingFunction,right_InterpolatingFunction,x_,opt:OptionsPattern[]]:=Module[
	{lCoord,rCoord,coord,opts,lFac=OptionValue["LeftMultiplier"],rFac=OptionValue["RightMultiplier"]},
	opts=Evaluate@FilterRules[Join[{opt},Options[InterpolationUnion]],Options[Interpolation]];
	lCoord=MapAt[lFac #&,Select[InterpolationData[left],First[#]<x&],{All,2}];
	rCoord=MapAt[rFac #&,Select[InterpolationData[right],First[#]>=x&],{All,2}];
	Interpolation[Join[lCoord,rCoord],opts]
	]

(* extend solution for velocity, to include v=0 regions *)
Options[ExtendProfile]={"Type"->"v","wMinus"->None};
ExtendProfile[prof_InterpolatingFunction,OptionsPattern[]]:=Module[{\[Xi]s,fs,d\[Xi],wp,wm,type=OptionValue["Type"]},
	{\[Xi]s,fs}=InterpolationData[prof]\[Transpose];
	d\[Xi]=Min[Differences[\[Xi]s]];
	wm=OptionValue["wMinus"]/.None->fs[[1]];
	(* extend \[Xi] domain *)
	\[Xi]s=Join[{0.,\[Xi]s[[1]]-d\[Xi]},\[Xi]s,{\[Xi]s[[-1]]+d\[Xi],1.}];
	fs=Switch[type,
		"v",Join[{0.,0.},fs,{0.,0.}], (* extend velocity profile *)
		"w",Join[Table[wm,2],fs,{1.,1.}], (* extend enthalpy profile *)
		True,Abort[]
		];
	Interpolation[{\[Xi]s,fs}\[Transpose],InterpolationOrder->1]
	]


(* ::Text:: *)
(*Modify built-in functions (not safe!!)*)


Unprotect[NDSolveValue,Interpolation];
AppendTo[Options[NDSolveValue],"ExtrapolationHandler"->{Indeterminate&,"WarningMessage"->False}];
AppendTo[Options[Interpolation],"ExtrapolationHandler"->{Indeterminate&,"WarningMessage"->False}];
Protect[NDSolveValue,Interpolation];


(* ::Subsection::Closed:: *)
(*Fluid equations*)


(* ::Input::Initialization:: *)
fluidEqs=<|
	"d\[Xi]dv"->\[Xi]'[v]==(\[Gamma]sq[v]\[Xi][v](1-v \[Xi][v]))/(2v) (\[Mu][\[Xi][v],v]^2/("cs2")-1), (* v[\[Xi]] diverges, then we solve for \[Xi][v] and invert *)
	"dlnwdv"->(1+1/("cs2"))\[Gamma]sq[v]\[Mu][\[Xi],v],
	"dvd\[Xi]"->v'[\[Xi]]==(2v[\[Xi]])/(\[Gamma]sq[v[\[Xi]]]\[Xi](1-v[\[Xi]] \[Xi])) (\[Mu][\[Xi],v[\[Xi]]]^2/("cs2")-1)^-1,
	"dwd\[Xi]"->w'[\[Xi]]==w[\[Xi]](1+1/("cs2"))\[Gamma]sq[v[\[Xi]]]\[Mu][\[Xi],v[\[Xi]]]v'[\[Xi]]
|>;

(* integrate the enthalpy, Eq. 29 of Espinosa:2010hh *)
intW[\[CapitalXi]_,\[Xi]0_,v_,cs2_]:=NIntegrate[v'[\[Xi]](1+1/cs2)\[Gamma]sq[v[\[Xi]]]\[Mu][\[Xi],v[\[Xi]]],{\[Xi],\[Xi]0,\[CapitalXi]}]


(* ::Subsection::Closed:: *)
(*Integration of enthalpy*)


(* ::Input::Initialization:: *)
(* enthalpy integrand *)
integrandW[\[Xi]_, v_, cs2_] :=(1 + 1/cs2)\[Gamma]sq[v]\[Mu][v, \[Xi]]
integrateW[\[Xi]s_List,vs_List, cs2_] :=Module[
{ fVals, dvs, trapz, cumSum,\[Xi]sSort,vsSort},
{\[Xi]sSort,vsSort}=ReverseSortBy[{\[Xi]s,vs}\[Transpose],First,Order]\[Transpose];
	 (* Integrand evaluated at every (xi, v) pair *)
	 fVals =integrandW[\[Xi]sSort,vsSort, cs2] ;
	 dvs =Differences[vsSort]; (* Step sizes in v *)
	 trapz =dvs MovingAverage[fVals,2]; (* Trapezoid areas between consecutive pairs *)
	 cumSum=Prepend[Accumulate[trapz], 0.]; (* Cumulative sum; prepend 0 so w[[1]]=Exp[0]=1 *)
	Reverse@ Exp[cumSum]
	]
ComputeW::usage="ComputeW[v_,\[Xi]w_,vWFrame] computes the enthalpy w from a fluid velocity profile and the wall velocity.";
Options[ComputeW]={"cs2"->cs2};
ComputeW[v_InterpolatingFunction,\[Xi]w_Real,opts:OptionsPattern[]]:=Module[{\[Xi]s,vs,ws,norm,cs2},
	cs2=OptionValue["cs2"];
	{\[Xi]s,vs}=InterpolationData[v]\[Transpose];
	ws=integrateW[\[Xi]s,vs,cs2];
	{\[Xi]s,ws}\[Transpose]
	]


(* ::Subsection::Closed:: *)
(*Plotting*)


(* ::Input::Initialization:: *)
PlotProfile[profile_,opts:OptionsPattern[]]:=Module[{opt={opts,Frame->True,ImageSize->Medium}},
Column[{
	Plot[profile["v"][\[Xi]],{\[Xi],Sequence@@(profile["v"]["Domain"][[1]]+0.01{-1,1})},PlotStyle->Blue,ImagePadding->{{50,1},{0,1}},FrameLabel->{None,"v"},Evaluate[Sequence@@opt],PlotRange->All],
	Plot[profile["w"][\[Xi]],{\[Xi],Sequence@@(profile["v"]["Domain"][[1]]+0.01{-1,1})},PlotStyle->Orange,PlotRange->All,ImagePadding->{{50,1},{0,0}},FrameLabel->{None,"w/\!\(\*OverscriptBox[\(w\), \(_\)]\)"},Evaluate[Sequence@@opt]],
	Plot[profile["\[Lambda]"][\[Xi]],{\[Xi],Sequence@@(profile["v"]["Domain"][[1]]+0.01{-1,1})},PlotStyle->Purple,PlotRange->All,ImagePadding->{{50,1},{Automatic,0}},FrameLabel->{"\[Xi]","\[Lambda]"},Evaluate[Sequence@@opt]]
	},Spacings->0]]


ManipulateVelocityProfile[]:=Manipulate[Plot[{
	Profile[v[[1]],v[[2]],"cs2"->cs^2,"Input"->"v"]["v"][\[Xi]]//Evaluate,
	\[Mu][cs,\[Xi]],
	shock[cs,\[Xi]],
	\[Xi]
	},{\[Xi],0,1},
	PlotRange->{{0,1},{0,1}},
	PlotStyle->{
	BubbleType[v[[1]],v[[2]],"cs2"->cs^2,"Input"->"v"]/.{"Det"->Blue,"Hyb"->Orange,"Def"->Red},
	Directive[Dashed,Gray],Directive[Dotted,Gray],Directive[Gray,Mesh->None]},
	GridLines->Automatic,
	Frame->True,FrameLabel->{"\[Xi]","v"},
	Filling->{1->Axis,4->Top},
	PlotLabel->"Self-similar profile"
	],
	{{v,{\[Xi]wDet,vwDet}},Locator},
	{{cs,Sqrt[cs2]},0,1},
	ControlPlacement->Top
	]


(* ::Section::Closed:: *)
(*Boundary conditions*)


(* ::Subsubsection::Closed:: *)
(*Default values*)


(* ::Input::Initialization:: *)
(* from Hindmarsh:2019phv *)
cs2=1/3;

(* detonation *)
\[Xi]wDet=0.77;
\[Alpha]Det=0.091;
vwDet=0.33; (* Subscript[v, -] for a detonation *)

(* deflagration *)
\[Alpha]Def=0.263;
\[Xi]wDef=0.5;
vwDef=.45;

(* hybrid *)
\[Xi]wHyb=.7;
vwHyb=.45;
\[Alpha]Hyb=.052;


(* ::Subsubsection::Closed:: *)
(*Conversion of boundary conditions*)


(* ::Input::Initialization:: *)
VWFrame::"VWFrame[vw,\[Alpha]] returns the velocity in the wall frame from the wall velocity and \[Alpha].";
Options[VWFrame]={ "Plus" ->True, "Sign" ->"Plus"};
VWFrame[vw_, \[Alpha]_, opts : OptionsPattern[]] :=
	 Module[{plusBranch=OptionValue[ "Plus"], sg=OptionValue[ "Sign"],
	   a1, a2, aa},
	 If[sg === "Plus",
	  a1=1/(3 vw) + vw;
	  a2=Sqrt[(1/(3 vw) - vw)^2 + 4 \[Alpha]^2 + 8/3 \[Alpha]];
	  aa=0.5/(1 + \[Alpha]),
	 (* Sign === "minus" *)
	  a1=(1 + \[Alpha]) vw + (1 - 3 \[Alpha])/(3 vw);
	  a2=Sqrt[((1 + \[Alpha]) vw + (1 - 3 \[Alpha])/(3 vw))^2 - 4/3];
	  aa=0.5
	 ];
	 aa If[plusBranch, a1+a2, a1-a2]
	 ]

Options[VplusVminus]={ "Type" ->"Det", "cs2" ->cs2};
VplusVminus[\[Xi]w_,\[Alpha]_, opts : OptionsPattern[]] :=Module[{ 
ty=OptionValue["Type"],
cs2=OptionValue["cs2"], cs, vplus, vminus
},
cs=Sqrt[cs2];
Which[
	  ty === "Det",
	  vplus =\[Xi]w;
	  vminus=VWFrame[\[Xi]w, \[Alpha], "Plus" ->True, "Sign" ->"minus"],
	  ty === "Def",
	  vminus=\[Xi]w;
	  vplus =VWFrame[\[Xi]w, \[Alpha], "Plus" ->False, "Sign" ->"Plus"],
	 ty === "Hyb",
	 vminus=cs;
	 vplus =VWFrame[cs, \[Alpha], "Plus" ->False, "Sign" ->"Plus"]
	 ];
	 {vplus, vminus}
	 ]


(* ::Section::Closed:: *)
(*Profiles*)


(* ::Subsection::Closed:: *)
(*\[Alpha]/\[Alpha]+ shooting algorithm*)


(* ::Input::Initialization:: *)
Options[ProfileFromTargetAlpha]={"Type"->Automatic,"MaxIterations"->10,"Tolerance"->0.001,"Verbose"->False,"cs2"->cs2};
ProfileFromTargetAlpha[\[Xi]w_,\[Alpha]_,opts:OptionsPattern[{ProfileFromTargetAlpha,Deflagration,Hybrid}]]:=Module[{
	\[Alpha]Plus=\[Alpha],\[Alpha]Comp,profile,sol,type,niter=1,
	wp,res,optsProf,tol=OptionValue["Tolerance"],maxIter=OptionValue["MaxIterations"],conv=False},
	type=OptionValue["Type"]/.Automatic->BubbleType[\[Xi]w,\[Alpha]];
	profile=type/.BubbleTypeRules;
	optsProf=Evaluate@FilterRules[Join[{opts},Options[profile]],Options[profile]];
	If[OptionValue["Verbose"],Print[StringTemplate["Running shooting algorithm for target \[Alpha]=``"][\[Alpha]]]];
	res=NestWhileList[(
		sol=profile[\[Xi]w,#[[2]],"Input"->"\[Alpha]+",optsProf,"Verbose"->False];
		wp=sol["wp"];
		\[Alpha]Comp=\[Alpha]Plus wp;
		If[Abs[\[Alpha]Comp-\[Alpha]]/\[Alpha]<tol,conv=True];
		{\[Alpha]Comp,\[Alpha]Plus=\[Alpha]/wp}
	)&,{None,\[Alpha]},(niter<maxIter\[And]Not[conv])&
	];
	sol["\[Lambda]"]=\[Lambda]FromWSol[sol["w"],\[Xi]w,\[Alpha],OptionValue["cs2"]];
	sol["\[Alpha]Computed"]=\[Alpha]Comp;
	sol["\[Alpha]+"]=\[Alpha]Plus;
	sol["\[Alpha]ComputedList"]=res[[2;;,1]];
	sol["\[Alpha]+List"]=res[[2;;,2]];
	sol
	]


(* ::Subsection::Closed:: *)
(*Detonation*)


(* ::Input::Initialization:: *)
Options[Detonation]={"cs2"->cs2,"Verbose"->False,"Input"->"\[Alpha]","Normalize"->True};
$options["Det"]=Options[Detonation];
Detonation[\[Xi]w_,in_,opts:OptionsPattern[]]:=Module[{
	shellDom,cs2,\[Alpha],vm,vmWFrame,BCond,cs2rule,vSol,wSol,\[Lambda]Sol,wp=1.,wm,w0
	},
	cs2rule=FilterRules[{opts}~Join~Options[Detonation],"cs2"]//First;
	cs2=1.OptionValue["cs2"];
	If[OptionValue["Input"]==="\[Alpha]",
		\[Alpha]=in;
		vmWFrame=VplusVminus[\[Xi]w,\[Alpha],cs2rule,"Type"->"Det"][[2]];
		vm=\[Mu][vmWFrame,\[Xi]w],
		(* if input is Subscript[v, -] *)
		vm=in;
		vmWFrame=\[Mu][vm,\[Xi]w];
		\[Alpha]=1.
		];
	(* solve ODE *)
	BCond=(v[\[Xi]w]==vm);
	vSol=NDSolveValue[{fluidEqs["dvd\[Xi]"]/.cs2rule,BCond},v,{\[Xi],Sqrt[cs2],\[Xi]w}];
	wSol=ComputeW[vSol,\[Xi]w];
	wm=wRatio[vmWFrame,\[Xi]w];
	If[OptionValue["Normalize"],wSol[[All,2]]*=wm];
	wSol=Interpolation[wSol,InterpolationOrder->1];
	w0=wSol["ValuesOnGrid"][[1]];
	(* extend profiles *)
	shellDom=vSol["Domain"]//First;
	vSol=ExtendProfile[vSol];
	wSol=ExtendProfile[wSol,"Type"->"w"];
	(* fluid density perturbations from enthalpy. \[Alpha] is added only inside the bubble (0 otherwise) *)
	\[Lambda]Sol=\[Lambda]FromWSol[wSol,\[Xi]w,\[Alpha],cs2];
	If[OptionValue["Verbose"],Print[TableForm[ {{\[Alpha],\[Xi]w,vmWFrame,vm,wm,w0}}\[Transpose],TableHeadings->{{"\[Alpha]","\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\)","\!\(\*SubscriptBox[OverscriptBox[\(v\), \(~\)], \(-\)]\)","\!\(\*SubscriptBox[\(v\), \(-\)]\)","\!\(\*SubscriptBox[\(w\), \(-\)]\)","\!\(\*SubscriptBox[\(w\), \(0\)]\)"}}]]];
	AssociationThread[{"v","w","\[Lambda]","\[Xi]w","vm","vp","wp","wm","w0","\[Alpha]","ShellDomain"},{vSol,wSol,\[Lambda]Sol,\[Xi]w,vm,Sqrt[cs2],wp,wm,w0,\[Alpha],shellDom}]
	]


(* ::Subsection::Closed:: *)
(*Deflagration*)


(* ::Input::Initialization:: *)
Options[Deflagration]={"cs2"->cs2,"Verbose"->False,"Input"->"\[Alpha]","Normalize"->True,"\[Lambda]From\[Alpha]"->False};
$options["Def"]=Union[Options[Deflagration],Options[ProfileFromTargetAlpha]];
Deflagration[\[Xi]w_,in_,opts:OptionsPattern[{Deflagration,ProfileFromTargetAlpha}]]:=Module[{
	cs2,cs,\[Alpha],\[Alpha]Plus,vp,vpWFrame,wp,wm,BCond,cs2rule,vSol,wSol,\[Lambda]Sol=None,\[Xi]s,shellDom,optsProf,
	\[Xi]sh,wsh,input=OptionValue["Input"],\[Alpha]2=OptionValue["\[Lambda]From\[Alpha]"],verbose=OptionValue["Verbose"]
	},
	If[input=="\[Alpha]",
		\[Alpha]=in;
		(* shooting algorithm for \[Alpha]+ *)
		optsProf=Evaluate@FilterRules[Join[{opts},$options["Def"]],$options["Def"]];
		ProfileFromTargetAlpha[\[Xi]w,\[Alpha],optsProf,"Type"->"Def"],
		(* from given \[Alpha]+ or v+ *)
		cs2rule=FilterRules[{opts}~Join~Options[Deflagration],"cs2"]//First;
		cs2=1.OptionValue["cs2"];cs=Sqrt[cs2];
		If[input=="\[Alpha]+",
			\[Alpha]Plus=in;
			vpWFrame=VplusVminus[\[Xi]w,\[Alpha]Plus,cs2rule,"Type"->"Def"][[1]];
			vp=\[Mu][vpWFrame,\[Xi]w],
			(* if input is v+ *)
			vp=in;
			vpWFrame=\[Mu][vp,\[Xi]w];
			\[Alpha]Plus=1. (* to avoid breaking \[Lambda]Sol *)
			];
		(* solve ODE *)
		BCond=(v[\[Xi]w]==vp);
		vSol=NDSolveValue[{
			fluidEqs["dvd\[Xi]"]/.cs2rule,BCond,
			WhenEvent[v[\[Xi]]<shock[cs,\[Xi]],"StopIntegration"]
			},v,{\[Xi],\[Xi]w,1}];
		shellDom=vSol["Domain"]//First;
		(* integrate enthalpy *)
		wSol=ComputeW[vSol,\[Xi]w];
		\[Xi]sh=shellDom[[2]];
		wsh=wShock[\[Xi]sh];
		If[OptionValue["Normalize"],wSol[[All,2]]*=wsh/wSol[[-1,2]]]; (* normalize enthalpy on the shell *)
		wp=wSol[[1,2]];
		wm=wp/wRatio[vpWFrame,\[Xi]w];
		wSol=Interpolation[wSol,InterpolationOrder->1];
		(* extend profiles *)
		vSol=ExtendProfile[vSol];
		wSol=ExtendProfile[wSol,"Type"->"w","wMinus"->wm];
		(* fluid density perturbations from enthalpy. \[Alpha] is added only inside the bubble (0 otherwise) *)
		If[\[Alpha]2, \[Lambda]Sol=\[Lambda]FromWSol[wSol,\[Xi]w,\[Alpha]2,cs2]];
	If[verbose,Print[TableForm[ {{\[Alpha],\[Xi]w,vpWFrame,vp,wm,\[Xi]sh,wsh}}\[Transpose],
TableHeadings->{{"\[Alpha]+","\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\)","\!\(\*SubscriptBox[OverscriptBox[\(v\), \(~\)], \(+\)]\)","\!\(\*SubscriptBox[\(v\), \(+\)]\)","\!\(\*SubscriptBox[\(w\), \(-\)]\)","\!\(\*SubscriptBox[\(\[Xi]\), \(sh\)]\)","\!\(\*SubscriptBox[\(w\), \(sh\)]\)"}}]]];
		AssociationThread[{"v","w","\[Lambda]","\[Xi]w","vp","vm","wp","wm","wsh","\[Alpha]+","ShellDomain"},
		{vSol,wSol,\[Lambda]Sol,\[Xi]w,vp,0.,wp,wm,wsh,\[Alpha],shellDom}]
	]]


(* ::Subsection::Closed:: *)
(*Hybrid*)


(* ::Input::Initialization:: *)
ClearAll[Hybrid];
Options[Hybrid]={"cs2"->cs2,"Verbose"->False,"Input"->"\[Alpha]","\[Lambda]From\[Alpha]"->False};
$options["Hyb"]=Union[Options[Hybrid],Options[ProfileFromTargetAlpha]];
Hybrid[\[Xi]w_,in_,opts:OptionsPattern[{Hybrid,ProfileFromTargetAlpha}]]:=Module[{
cs2,cs,\[Alpha],vp,vpWFrame,vm,vmWFrame,
BCond,cs2rule,solDef,solDet,vSol,wSol,\[Lambda]Sol=None,shellDom,wp,wm,w0,
optsProf,input=OptionValue["Input"],verbose=OptionValue["Verbose"]},
	If[input=="\[Alpha]",
	\[Alpha]=in;
	(* shooting algorithm for \[Alpha]+ *)
	optsProf=Evaluate@FilterRules[Join[{opts},$options["Hyb"]],$options["Hyb"]];
	ProfileFromTargetAlpha[\[Xi]w,\[Alpha],optsProf,"Type"->"Hyb"],
	(* from given \[Alpha]+ or v+ *)
		cs2rule=FilterRules[{opts}~Join~Options[Detonation],"cs2"]//First;
		cs2=1.OptionValue["cs2"];cs=Sqrt[cs2];
		If[input=="\[Alpha]+",
			\[Alpha]=in;
			{vpWFrame,vmWFrame}=VplusVminus[\[Xi]w,\[Alpha],cs2rule,"Type"->"Hyb"];
			vp=\[Mu][vpWFrame,\[Xi]w],
			(* if input is Subscript[v, +] *)
			vp=in;
			vpWFrame=\[Mu][vp,\[Xi]w];
			vmWFrame=cs;
			\[Alpha]=1. (* to avoid breaking \[Lambda]Sol *)
			];
		vm=\[Mu][vmWFrame,\[Xi]w]; (* v- is the Lorentz boosted v-^~ = cs *)
		(* get def and det partial solutions *)
		solDef=Deflagration[\[Xi]w,vp,"Input"->"v","cs2"->cs2];
		solDet=Detonation[\[Xi]w,vm-10.^-3,"Input"->"v","cs2"->cs2,"Normalize"->False];
		shellDom={First@solDet["ShellDomain"],Last@solDef["ShellDomain"]};
		wp=solDef["wp"];
		wm=wp /wRatio[vpWFrame,cs];
		{vSol,wSol}=InterpolationUnion[solDet[#],solDef[#],\[Xi]w,"LeftMultiplier"->#/.{"v"->1,"w"->wm}]&/@{"v","w"};
	w0=wSol["ValuesOnGrid"][[1]];
		(* energy density perturbation *)
		If[OptionValue["\[Lambda]From\[Alpha]"],\[Lambda]Sol=\[Lambda]FromWSol[wSol,\[Xi]w,\[Alpha],cs2]];
		(* print and return results *)
		If[verbose,
			Print[TableForm[ {{\[Alpha],\[Xi]w,vpWFrame,vp,wp,wm}}\[Transpose],
			TableHeadings->{{"\[Alpha]","\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\)","\!\(\*SubscriptBox[OverscriptBox[\(v\), \(~\)], \(+\)]\)","\!\(\*SubscriptBox[\(v\), \(+\)]\)","\!\(\*SubscriptBox[\(w\), \(+\)]\)","\!\(\*SubscriptBox[\(w\), \(-\)]\)"}}]]];
		AssociationThread[{"v","w","\[Lambda]","\[Xi]w","vp","vm","wp","wm","w0","\[Alpha]","ShellDomain"},
			{vSol,wSol,\[Lambda]Sol,\[Xi]w,vp,vm,wp,wm,w0,\[Alpha],shellDom}]
	]
	]
$options["Hyb"]=Union[Options[Hybrid],Options[ProfileFromTargetAlpha]];


(* ::Subsection::Closed:: *)
(*Generic profile*)


(* ::Input::Initialization:: *)
Profile::"Profile[\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\),\[Alpha]] gives the fluid profile solution as a function the wall velocity and transition strength.
Profile[\!\(\*SubscriptBox[\(\[Xi]\), \(w\)]\),\!\(\*SubscriptBox[\(v\), \(w\)]\),\"Input\"\[Rule]\"v\" uses the fluid velocity at the wall (in the fluid frame) as second input.";
Options[Profile]={"cs2"->cs2,"Input"->"\[Alpha]","Verbose"->False};
Profile[\[Xi]w_,\[Alpha]_,opt:OptionsPattern[{Profile,Detonation,Deflagration,Hybrid,ProfileFromTargetAlpha}]]:=Module[{optsBub,optsProf,profile,profileFunc,type},
	optsBub=Evaluate@FilterRules[{opt}~Join~Options[Profile],Options[BubbleType]];
	type=BubbleType[\[Xi]w,\[Alpha],optsBub]//If[OptionValue["Verbose"],Echo[#,"Bubble type \[Rule]"]&,Identity];
	optsProf=FilterRules[Join[{opt},$options[type]],$options[type]];
	profileFunc=type/.BubbleTypeRules;
	profile=profileFunc[\[Xi]w,\[Alpha],optsProf];
	profile["Type"]=type;
	profile
	]
AutoComplete[Profile];


(* ::Section::Closed:: *)
(*Example usage*)


(* ::Text:: *)
(*Solve for a given profile*)


(*
(* get a detonation profile *)
detSol=Detonation[\[Xi]w,\[Alpha],"Verbose"->True]
(* get a deflagration profile *)
defSol=Deflagration[\[Xi]wDef,\[Alpha]Def,"Input"->"\[Alpha]","Verbose"->True]
(* get a hybrid solution *)
hybSol=Hybrid[\[Xi]wHyb,\[Alpha]Hyb]//EchoTiming
(* get a generic profile *)
prof=Profile[\[Alpha],\[Xi]w]
*)


(* ::Text:: *)
(*Manipulate the initial conditions (Subscript[\[Xi], w],v(Subscript[\[Xi], w])), solve and plot the profile*)


(*Manipulate[Plot[{
Profile[v[[1]],v[[2]],"cs2"->cs^2,"Input"->"v"]["v"][\[Xi]]//Evaluate,
\[Mu][cs,\[Xi]],
shock[cs,\[Xi]],
\[Xi]
},{\[Xi],0,1},
PlotRange->{{0,1},{0,1}},
PlotStyle->{
BubbleType[v[[1]],v[[2]],"cs2"->cs^2,"Input"->"v"]/.{"Det"->Blue,"Hyb"->Orange,"Def"->Red},
Directive[Dashed,Gray],Directive[Dotted,Gray],Directive[Gray,Mesh->None]},
GridLines->Automatic,
Frame->True,FrameLabel->{"\[Xi]","v"},
Filling->{1->Axis,4->Top},
PlotLabel->"Self-similar profile"
],
{{v,{\[Xi]w,vw}},Locator},
{{cs,Sqrt[cs2]},0,1},
ControlPlacement->Top
]*)



(* ::Section:: *)
(*End*)


End[] (* "Private`" *)


EndPackage[]
