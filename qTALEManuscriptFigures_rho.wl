(* ::Package:: *)

(* Mathematica script *)


PrependTo[$Path, FileNameJoin[{DirectoryName[$InputFileName], "packages"}]];

<< rhoModels`

nlfit=fitRho["measured",formula->"dltwicoop",scaling->False];
{wtpredrules,predTALE}=getPredRules[nlfit,0.85,1];
"Figure3/TALEDataAndPredictions.pdf"~export~ListPlot[Join[Transpose[predTALE[[All,1;;6]]],measuredMeans[[1;;6]]],Joined->True,ImageSize->800]
"Figure3/TALEPredictions.pdf"~export~ListPlot[Transpose[predTALE[[All,1;;6]]],Joined->True,ImageSize->800,AspectRatio->1/2]


"Figure3/measuredAndPrediction_WT1stdDev1.0.pdf"~export~plotTALE[predTALE,{3},1]
"Figure3/measuredAndPrediction_TALERx1stdDev1.0.pdf"~export~plotTALE[predTALE,{2,3},1]
"Figure3/measuredAndPrediction_TALERx2stdDev1.0.pdf"~export~plotTALE[predTALE,{1,3},1]
"Figure3/measuredAndPrediction_TALEAx1stdDev1.0.pdf"~export~plotTALE[predTALE,{4,3},1]
"Figure3/measuredAndPrediction_TALEAx2stdDev1.0.pdf"~export~plotTALE[predTALE,{5,3},1]
"Figure3/measuredAndPrediction_TALEAx3stdDev1.0.pdf"~export~plotTALE[predTALE,{6,3},1]

iterRange={10,60,95,115,135,165,220,275};
result=Reap[Partition[Table[responsePlot[nlfit,wtpredrules,i,.85],{i,iterRange}],10,10,1,{}]];
sigmoidal=Partition[Table[combinedSigmoidPlot[nlfit,wtpredrules,i,3*.85],{i,iterRange}],10,10,1,{}];
positions=Partition[Table[positionPlot[i],{i,iterRange}],10,10,1,{}];

"Figure3/sigmoidalResponse.pdf"~export~GraphicsGrid[Join[result[[1]],sigmoidal,positions],ImageSize->1800,ImageMargins->15]

"Figure3/WTPrediction.pdf"~export~mutantPrediction[nlfit,Unevaluated[Sequence[]],"WT prediction"]

"FigureS6/ZinzenIndependent.pdf"~export~fitAndPlotTraining["rhoscaled",formula->"independent",scaling->False]
"FigureS6/ZinzenSubmodule.pdf"~export~fitAndPlotTraining["rhoscaled",formula->"dltwicoop",scaling->False]

"FigureS7/rho_snaMutant.pdf"~export~mutantPrediction[nlfit,"sna"->0,"sna mutant"]
"FigureS7/rho_dlMutant.pdf"~export~mutantPrediction[nlfit,"dl"->0,"dl mutant"]
"FigureS7/rho_twiMutant.pdf"~export~mutantPrediction[nlfit,"twi"->0,"twi mutant"]



linearrep=Times[a_.,1 / (1+Exp[b_])] :> -b;
coefs=Replace[nlfit["BestFit"],linearrep] ;
dltwi=FirstCase[coefs,a_ (b_-1),"",Infinity,Heads->True];
"SupplementaryNote/rho_dltwiComponent.png"~export~Plot[Evaluate[dltwi /. dl+twi->x],{x,0,1},AxesLabel->{"\!\(\*
StyleBox[\"dl\",\nFontSlant->\"Italic\"]\)+\!\(\*
StyleBox[\"twi\",\nFontSlant->\"Italic\"]\)"},ImageSize->300]

dltwiStrings=dltwi /. {twi->"twi",dl->"dl"};
"SupplementaryNote/rho_dlTwiActivity.png"~export~plotActivity[dltwiStrings]
"SupplementaryNote/rho_parameterTable.pdf"~export~modelParameterTable[nlfit]

