(* ::Package:: *)

BeginPackage["utilities`"]

export::usage = "export[filename,expr] exports expr to filename in the Figures directory."
ChangePropensity::usage = "ChangePropensity[model,offset] changes the intercept of the model by offset."
modelParameterTable::usage = "modelParamterTable[model] displays information about the model parameters."
tickLabelsOnly::usage="Return tick labels without marks."

Begin["`Private`"]

export[filename_,expr_,opts:OptionsPattern[Export]] :=
    Module[ {fullname},
        fullname = FileNameJoin[{"Figures/",filename}];
        Export[fullname,expr,opts,Background->None,ImageResolution->1200];
        expr
    ]

ChangePropensity[model_FittedModel,offset_] :=
	ChangePropensity[model["BestFit"],offset]
    
ChangePropensity[fitFun_,offset_] :=
    fitFun /. Power[Plus[numerator_,Power[E,Plus[a_,r__]]],-1]:>Power[Plus[numerator,Power[E,Plus[a-offset,r]]],-1]

modelParameterTable[model_]:=Module[{paramNames,confidenceIntervalTable},
	(* assuming parameters are returned in the same order *)
	paramNames=Replace[model["BestFitParameters"],{Rule[Subscript[name_,index_],value_]:>Subscript[SymbolName[name],index],Rule[name_,value_]:>SymbolName[name]},1];
	confidenceIntervalTable=Join[Transpose[{paramNames}],model["ParameterConfidenceIntervalTableEntries",ConfidenceLevel->0.95],2];
    Column[{model["BestFit"],Grid[MapAll[If[MachineNumberQ[#],NumberForm[#,{6,2}],#]&,confidenceIntervalTable]]},Left,2]
]

tickLabelsOnly[max_Integer,n_Integer]:={#,#,0}&/@FindDivisions[{0,max},n]
tickLabelsOnly[max_Real,n_Integer]:={#,NumberForm[N@#,{2,1}],0}&/@FindDivisions[{0,max},n]

End[]

EndPackage[]






