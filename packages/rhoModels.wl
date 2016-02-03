(* ::Package:: *)

BeginPackage["rhoModels`", { "utilities`"}]

measuredMeans::usage="Averaged data from rho TALE series."
fitRho::usage="fitRho[outputGene] fits a model to outputGene."
getPredRules::usage="getPredRules[fit,actInc,repInc] extracts the prediction rules from fit, including for TALEA and TALER using the increments for activation and repression."
plotTALE::usage="Plot TALE data and prediction."
fitAndPlotTraining::usage="Fit and plot with training data."
plotActivity::usage="Plots activity of a formula in terms of the training data."


responsePlot::usage="Plots response of model."
positionPlot::usage="Plots position of response along DV axis."
combinedSigmoidPlot::usage="Plots region of response on sigmoidal curve."
mutantPrediction::usage="Plot mutant predictions."


Begin["`Private`"] (* Begin Private Context *) 


(* Import data *)
measuredImport=Import[FileNameJoin[{ParentDirectory[DirectoryName[$InputFileName]],"data/rho_NEE_TALE_series.xlsx"}]];
minX=0; maxX=0.65;
background=10;
adjData=Map[Max[0,#]&,measuredImport-background,{3}][[{5,6,1,2,3,4}]];
wtindex=3;
maxTALEAindex=6;
initialMeans=Map[Mean,adjData,{2}];
adjMeasured=adjData/Max[initialMeans];
measuredMeans=initialMeans/Max[initialMeans];

dvData=Import[FileNameJoin[{ParentDirectory[DirectoryName[$InputFileName]],"externalData/Papatensko_DV_DATA_8.xls"}]];
columnlabels={"dl","twi","sna","rho"};
importdata=Dataset[Association@Thread[columnlabels-># ]& /@dvData[[1,All,{15,33,51,81}]]]; (* This selects the smoothed data *)
importdatarefined=importdata[All,All,If[#=="",0,#]&];
data=importdatarefined[Range[1,648,2]];  (* This is a simpler, less general method of taking every second measurement *)
newdata=Join[Normal[data],Association/@Thread["rhoscaled"->Rescale[Normal[data[All,"rho"]]]*Max[measuredMeans[[wtindex]]]],Association/@Thread["measured"->measuredMeans[[wtindex]]],2]//Dataset;


(* Function definitions *)


colorsTALEA=Lighter/@Table[ColorData["BlueGreenYellow"][i],{i,{.1,.3,.7,.9}}];
colorsTALER=Table[ColorData["Aquamarine"][i],{i,{.5,.9}}];
colors=Join[colorsTALER,colorsTALEA]

listPlotVariation[embryoData_,index_,sdprop_:0.5]:=Module[{interval},
	interval=Transpose@Map[{Mean[#]-StandardDeviation[#]*sdprop,Mean[#]+StandardDeviation[#]*sdprop}&,embryoData];
	ListPlot[interval,PlotRange->{{minX,maxX}*100,{0,1.2}},Filling->{1->{2}},FillingStyle->Directive[{Opacity[0.5],colors[[index]]}],PlotStyle->{{Thin,Gray},{Thin,Gray}},Joined->True,DataRange->{minX,maxX}*100]
]

plotTALE[predTALE_,index_,sdprop_:1]:=Module[{},
	Show[listPlotVariation[adjMeasured[[#]],#,sdprop]&/@index,ListPlot[Transpose[predTALE[[All,index]]],PlotRange->{{minX,maxX}*100,{0,1.1}},DataRange->{minX,maxX}*100,Joined->True,PlotStyle->Thread[{Thickness[.004],Darker[#,0.2]&/@colors[[index]]}]],ImageSize->800,Ticks->None,Axes->False,AspectRatio->1/2]
]

Options[fitRho]={formula->"dltwicoop",scaling->False,subset->All};
fitRho::formula="`1` is not a defined model formula.";
fitRho[outputGene_,OptionsPattern[]]:=Module[{maxparamlist,fitformula,regs,geneNames,ssdata,vars,paramlist},
	maxparamlist={\[Alpha],Subscript[\[Beta], 0],Subscript[\[Beta], 1],Subscript[\[Beta], 2],Subscript[\[Beta], 3],Subscript[\[Gamma], 0],Subscript[\[Gamma], 1],Subscript[\[Gamma], 2],Subscript[\[Xi], 0],Subscript[\[Xi], 1],Subscript[\[Xi], 2]};
	fitformula=Switch[OptionValue[formula],
		"dltwicoop",\[Alpha]/(1 + Exp[-(Subscript[\[Beta], 0]+Subscript[\[Beta], 1]Symbol["sna"]+Subscript[\[Beta], 2](2/(1+E^(-Subscript[\[Gamma], 1](Symbol["dl"]+Symbol["twi"]))) -1))]),  
		"independent",\[Alpha]/(1 + Exp[-(Subscript[\[Beta], 0]+Subscript[\[Beta], 1]Symbol["sna"]+Subscript[\[Beta], 2] Symbol["twi"] +Subscript[\[Beta], 3] Symbol["dl"])]),
		_,Return[Message[fitRho::formula,OptionValue[formula]]]];
	fitformula=If[OptionValue[scaling],
		fitformula,
		fitformula /. {\[Alpha]->1}];
	paramlist=Pick[maxparamlist,!FreeQ[fitformula,#]&/@maxparamlist];
	regs={"dl","sna","twi"};
	vars=Symbol[#]& /@regs;
	geneNames=Append[regs,outputGene];
	ssdata=(newdata[OptionValue[subset],geneNames] //Normal //Values);
	NonlinearModelFit[ssdata,fitformula,paramlist,vars,
		PrecisionGoal->Automatic,AccuracyGoal->Automatic,MaxIterations->100,WorkingPrecision->Automatic,Method-> {"NMinimize",Method->"NelderMead"}]
]

getPredRules[fit_,actInc_ ,repInc_]:=Module[{wtpredrules,predTALEA,predTALER},
	wtpredrules=Normal[Normal[newdata[All,{"dl","sna","twi"}]]]/.{"dl"->Symbol["dl"],"sna"->Symbol["sna"],"twi"->Symbol["twi"]};
	predTALEA=Transpose[Table[ChangePropensity[fit,adjust] /. wtpredrules,{adjust,0,3 actInc,actInc}]];
	predTALER=Transpose[Table[ChangePropensity[fit,adjust] /. wtpredrules,{adjust,0,-2 repInc,-repInc}]];
	{wtpredrules,Join[predTALER[[All,{3,2}]],predTALEA,2]}
]


plotVsTraining[predTALE_,gene_,index_:{wtindex}]:=Module[{geneList},
	geneList=Append[{"dl","sna","twi"},gene];
	ListPlot[Join[Transpose[Values[Normal[newdata[All,geneList]]]],Transpose[predTALE[[All,index]]]],Joined->True,PlotStyle->Join[ConstantArray[Thickness[0.004],Length[geneList]],{Directive[Thickness[0.004],Black,Dashed],Directive[Thickness[0.004],Brown,Dashed]}], AspectRatio->1/2,PlotLegends->Placed[SwatchLegend[Join[geneList,{"prediction","TALEA X3"}]],{Right,Top}]]
]

Options[fitAndPlotTraining]=Options[fitRho];
fitAndPlotTraining[gene_,opts:OptionsPattern[]]:=Module[{nlfit,wtpredrules,predTALE},
	nlfit=fitRho[gene,FilterRules[{opts}, Options[fitRho]]];
	{wtpredrules,predTALE}=getPredRules[nlfit,0.85,1];
	plotVsTraining[predTALE,gene,{3,6}]
]


positionPlot[i_]:=ListPlot[{measuredMeans[[wtindex]],measuredMeans[[maxTALEAindex]]},Joined->True,Axes->False,PlotStyle->{Black,Gray},PlotRange->{Full,{0,1.1}},
	Epilog->{Orange,Arrowheads[Medium],Arrow[{{i,measuredMeans[[wtindex,i]]},{i,measuredMeans[[maxTALEAindex,i]]}}]}]

sigmoidPlot[nlfit_,wtpredrules_,i_,offset_]:=Module[{linearpred,posvalue,exprvalueStart,exprvalueEnd},
	linearpred=Replace[nlfit["BestFit"],1/(1+Exp[mu_]):>-mu];
	posvalue=linearpred/.wtpredrules[[i]];
	exprvalueStart=1/(1+Exp[-posvalue]);
	exprvalueEnd=1/(1+Exp[-posvalue-offset]);
	Show[{Plot[1/(1+Exp[-(posvalue+(x-1)offset)]),{x,-2,6},PlotRange->{{-2,6},{-0.1,1.2}},PlotStyle->{Gray,Thickness[.012]},Frame->True,Axes->False,FrameTicks->None]}]
]

responsePlot[nlfit_,wtpredrules_,i_,offset_]:=Show[{sigmoidPlot[nlfit,wtpredrules,i,offset],BoxWhiskerChart[adjMeasured[[{3,4,5,6},i]],{{"MedianMarker",1,Black}},
	PlotRange->{{-2,6},{0,1.1}},FrameTicks->None,AspectRatio->Full,ChartStyle->{EdgeForm[Black],colors}]}
]

combinedSigmoidPlot[nlfit_,wtpredrules_,i_,offset_]:=Module[{linearpred,posvalue,exprvalueStart,exprvalueEnd},
	linearpred=Replace[nlfit["BestFit"],1/(1+Exp[mu_]):>-mu];
	posvalue=linearpred/.wtpredrules[[i]];
	exprvalueStart=1/(1+Exp[-posvalue]);
	exprvalueEnd=1/(1+Exp[-posvalue-offset]);
	Show[Graphics[{FaceForm[None],Rectangle[{-15,0},{12,1.1}]},AspectRatio->1/5],Plot[1/(1+Exp[-x]),{x,-7,7},PlotRange->{{-7,7},{0,1}},PlotStyle->{Gray,Thickness[.015]}],
	Plot[1/(1+Exp[-x]),{x,posvalue,posvalue+offset},PlotRange->{{-7,7},{0,1}},PlotStyle->{Yellow,Thickness[.01]}]]
]


geneNames={"dl","sna","twi","measured"};

mutantPrediction[nlfit_,changeRules_,newName_String]:=Module[{predFunction,mutantData,mutantPred},
	predfunction=nlfit["BestFit"] /. reg_Symbol/;MemberQ[Normal[newdata[1,Keys]],ToString[reg]] :>Slot[ToString[reg]];
	mutantData=newdata[All,Function[Evaluate[<|#,changeRules|>]]];
	mutantPred=mutantData[All,Function[Evaluate[<|#,newName->predfunction|>]]];
	ListPlot[Transpose[mutantPred[All,Append[geneNames,newName]] //Normal //Values],Joined->True,AspectRatio->1/2,PlotRange->{Full,{0,1.1}},PlotLegends->Append[geneNames,newName],ImageSize->600]
]


plotActivity[formula_]:=Module[{activity},
	activity=Rescale@Normal@data[All,  formula /. # &];
	ListPlot[Join[Transpose[Values[Normal[data]]],{activity,1-activity}],Joined->True,Ticks->{tickLabelsOnly[Length[data],7],tickLabelsOnly[1.0,5]},AspectRatio->1/2,PlotLegends->Join[Normal[data[1,Keys]]/.{"dl"->"Dl","twi"->"Twi","sna"->"Sna"},{"Equivalent activation","Equivalent repression"}],ImageSize->300]
]


End[] (* End Private Context *)   
EndPackage[]
