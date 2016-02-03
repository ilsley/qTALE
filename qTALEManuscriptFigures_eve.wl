(* ::Package:: *)

(* Mathematica script *)

PrependTo[$Path, FileNameJoin[{DirectoryName[$InputFileName], "packages"}]];

<< eveModels`

(* Figure 1 *)

eve3NL = FlyModel[{3}, {"hbP", "kni", "tll"},
  modelMethod -> "nonlinear", trainingSelection -> "resetRest",
  trainingExcludeStripes -> {7}, SQregulators -> {"hbP"},
  mRNACohorts -> {3}, proteinCohorts -> {3},
  trainingFilter -> Abs[#z] < 40 \[And] (#eve > 0.2 \[Or] #eve < 0.1 ),
  scale -> True];

"Figure1/eve37_HbP.pdf"~export~PlotEmbryo[eve3NL["embryo"], "hbP"];
"Figure1/eve37_tll.pdf"~export~PlotEmbryo[eve3NL["embryo"],"tll"];
"Figure1/eve37_kni.pdf"~export~PlotEmbryo[eve3NL["embryo"],"kni"];

eve3fun = eve3NL["model"]["BestFit"];
"Figure1/eve37HbandKniTllZero.pdf"~export~
    ContourPlot[eve3fun /. tll -> 0, {hbP, 0, .45}, {kni, 0, .7},
      AxesLabel -> Automatic, PlotRange -> Full, PlotLegends -> Automatic,
      FrameLabel -> {"Hb", "kni"}, AspectRatio -> 1, ImageSize -> 300];

"Figure1/eve37NonlinearPrediction.pdf"~export~PlotModel[eve3NL, plotType -> "embryo" ];

"Figure1/eve37TALER_WT_1x_2x_4x.pdf"~export~PlotTALE[eve3NL, {0, -1, -2, -4}];


(* Figure 2 *)
eve4NL=FlyModel[{4},{"hb","kni","tll"},modelMethod->"nonlinear",trainingSelection->"resetRest",
	trainingExcludeStripes->{6},SQregulators->{"kni"},mRNACohorts->{1,2},proteinCohorts->{1,2},
	trainingFilter->Abs[#z]<40 \[And](#eve >0.2\[Or]#eve <0.1 ),scale->True];  

"Figure2/eve46NonlinearPrediction.pdf"~export~PlotModel[eve4NL];

plotFig2TALE[adj_]:=Module[{data},
	data=PlotModel[eve4NL["embryo"],ChangePropensity[eve4NL["model"],adj],plotType->"dataset"];
	Show[PlotEmbryo[data[Select[-10<#x\[And]Abs[#z]<40&]],"pred",pointSize2D->0.01],ImageSize->1000]]

"Figure2/eve46TALER_WT.pdf"~export~plotFig2TALE[0];
"Figure2/eve46TALER_1x.pdf"~export~plotFig2TALE[-2.2]
"Figure2/eve46TALER_2x.pdf"~export~plotFig2TALE[-4.4]

"Figure2/tllWithDVDifferences.pdf"~export~PlotEmbryoTll[eve4NL["embryo"][Select[#y>0&]],"tll",maxcolor->Red]

"Figure2/eve4_hbmRNA.pdf"~export~PlotEmbryo[eve4NL["embryo"],"hb"]


(* Supplementary Figures *)
"FigureS2/eve46_kni_0.75_Pred.pdf"~export~PlotModel[MutantEmbryo[eve4NL["embryo"],"kni",#kni*.75],eve4NL["model"]]
"FigureS2/eve46_kni_0_Pred.pdf"~export~PlotModel[MutantEmbryo[eve4NL["embryo"],"kni",0],eve4NL["model"]]

"FigureS3/eve46Ventral_WTPred.pdf"~export~PlotModel[eve4NL["embryo"],eve4NL["model"],dimensions->2,view->"ventral"]
"FigureS3/eve46Ventral_kniInput_0.4.pdf"~export~(PlotEmbryo[MutantEmbryo[eve4NL["embryo"],"kni",#kni+#sna*inc],"kni",dimensions->2,view->"ventral",pointSize2D->0.012,maxvalue->1.2] /. inc -> 0.4)
"FigureS3/eve46Ventral_snakniPred_0.2.pdf"~export~(PlotModel[MutantEmbryo[eve4NL["embryo"],"kni",#kni+#sna*inc],eve4NL["model"],dimensions->2,view->"ventral"] /. inc -> 0.2)
"FigureS3/eve46Ventral_snakniPred_0.4.pdf"~export~(PlotModel[MutantEmbryo[eve4NL["embryo"],"kni",#kni+#sna*inc],eve4NL["model"],dimensions->2,view->"ventral"] /. inc -> 0.4)

<< tllMutantData`

predfunction3=(eve3NL["model"]["BestFit"] /. hbP-> hb/. reg_Symbol/;MemberQ[probeNames,ToString[reg]] :>Slot[ToString[reg]]) /. Times[numerator_.,rest__]:>rest; 
predfunction4=(eve4NL["model"]["BestFit"] /.  reg_Symbol/;MemberQ[probeNames,ToString[reg]] :>Slot[ToString[reg]]) /. Times[numerator_.,rest__]:>rest;  

"FigureS5/JanssensDataWTPrediction.pdf"~export~Show[showPred[predfunction3,predfunction4,dataWT,7,0],ImageSize->600]
"FigureS5/JanssensDataTllMutantPrediction.pdf"~export~Show[showPred[predfunction3,predfunction4,dataMutant,7,0],ImageSize->600]


(* Supplementary Note *)
linearrep=Times[a_.,1 / (1+Exp[b_])] :> -b;
"SupplementaryNote/eve37_QuadraticComponent.pdf"~export~Plot[Evaluate[Plus@@Cases[eve3NL["model"]["BestFit"] /. linearrep ,Real_ hbP|Real_ hbP^2]],{hbP,0,.7},AxesLabel->{"Hb"}]
"SupplementaryNote/eve46_QuadraticComponent.pdf"~export~Plot[Evaluate[Plus@@Cases[eve4NL["model"]["BestFit"] /. linearrep ,Real_ kni|Real_ kni^2]],{kni,0,.6},AxesLabel->{"\!\(\*
StyleBox[\"kni\",\nFontSlant->\"Italic\"]\)"}]

eve4NLNoQuadNoTllFullDV=FlyModel[{4,6},{"hb","kni"},modelMethod->"nonlinear",trainingSelection->"resetRest",trainingExcludeStripes->{},SQregulators->{},mRNACohorts->{1,2},proteinCohorts->{1,2},trainingFilter->(#eve >0.2\[Or]#eve <0.1),scale->True];  
"SupplementaryNote/eve46_NoQuadPred.pdf"~export~PlotModel[eve4NLNoQuadNoTllFullDV]

"SupplementaryNote/eve37_parameterTable.pdf"~export~modelParameterTable[eve3NL["model"]]
"SupplementaryNote/eve46_parameterTable.pdf"~export~modelParameterTable[eve4NL["model"]]
"SupplementaryNote/eve46_NoQuad_parameterTable.pdf"~export~modelParameterTable[eve4NLNoQuadNoTllFullDV["model"]]


