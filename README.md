# MetaOdysseus

MetaOdysseus is an R-package designed for native MS, bottom-up and native top-down experiments aimed at:

- Determining metal-to-protein stoichiometry.
- Identification of metal binding sites and mapping of Cys residues coordinating metal ions.
- Analysis of protein structure and folding. 

MetaOdysseus address:
 -Charge-state deconvolution.
 -Statistical scoring.
 -Mass assignment of protein complexes.

The following sections will examplify some examples published in (doi:)

# Installation in R

`install.packages("devtools")`

`devtools::install_github("ManuelPerisDiaz/CysMpro")`

Download data set for next examples from:

https://github.com/ManuelPerisDiaz/Data_MetaOdysseus


# Examples

# Case 1
## Native ESI-MS (low resolution)
The first experimental approach listed above consists in the use of direct injection native ESI-MS to determine metal ion stoichiometry and the study of the charge distribution, which allows to get insights into protein conformation change upon metal ion binding.

1) Specify localization of experimental file. The experimental file used here was converted to xy file format.

The protein used here is apo-metallothionein 2. 

`setwd("Directory where is located the file")`

`LowRes_apoMT2<-read.table("0ZnMT2a.xy")`

We may plot the spectra:

`plot_Peak(LowRes_apoMT2,t1=900,t2=1800,norm = TRUE,dir="H:",ID=c("low"))`

2) Deconvolution

`lista<-ProcessESI_Decon(comb = "none",file=MT2,plotOriginal = FALSE, plotProcess =FALSE,t1=960,t2=1800,xy=TRUE,dir="C:",
                          SNR.Th=5,Wa=1, FWHM=0.7,top=20,interval=100,maxMass=9000,minMass=5000,chargemax = 10,chargemin = 1
  )`

3) Score the deconvolution results

`Scored<- Score_deconvolution(lista, intervaLo2=5,intervaLo=200,mz1=6000,mz2=7500)`

4) Mass assignment

`sequence<-("MDPNCSCAAGDSCTCAGSCKCKECKCTSCKKSCCSCCPVGCAKCAQGCICKGASDKCSCCA")`


`Assign_Mass(Scored,sequence, Metal=c("Zn"), NbindingMetal = 7, mzw=2,mono = FALSE)`

## Native ESI-MS (high resolution)

1) Specify localization of experimental file. The experimental file used here was converted to xy file format.

`HighRes_apoMT2<-read.table("20180703_ApoThionein000003.txt")`

`plot_Peak(HighRes_apoMT2,t1=1204,t2=1215,norm = TRUE,dir="C:",ID=c("high"))`

`MT2<-HighRes_apoMT2`

2a) Deconvolution with simulation algorithm

 ` lista<-ProcessESI_Decon(comb = "comb3",file=MT2,plotOriginal = FALSE, plotProcess =FALSE,t1=960,t2=1800,xy=TRUE,dir="H:/MS-DUAL",
                          SNR.Th=5,Wa=0.5, FWHM=0.8 ,top=50,interval=20,maxMass=9000,minMass=5000,chargemax = 10,chargemin = 1
  )`
  
2b) Deconvolution with zscore algorithm

`DeconvHigh<-Deconvolved_ESI_High(exp=MT2,maxMass = 8000,peak_picking = TRUE,SNR.Th = 6)`

3) Score the deconvolution results

`Scored<- Score_deconvolution(lista, intervaLo2=10,intervaLo=200,mz1=6000,mz2=7500)`

4) Mass assignment

`sequence<-("MDPNCSCAAGDSCTCAGSCKCKECKCTSCKKSCCSCCPVGCAKCAQGCICKGASDKCSCCA")`

`Assign_Mass(Scored,sequence, Metal=c("Zn"), NbindingMetal = 7, mzw=3,mono=TRUE)`



### Zn7_MT2 labelled with IAM "high resolution"


`MT2<-read.table("Zn7MT2_IAM_50mM.xy")`

`plot_Peak(MT2,t1=960,t2=1600,norm = TRUE,dir="H:",ID=c("iam"))`


1)   Deconvolution

`  lista<-ProcessESI_Decon(comb = "comb3",file=MT2,plotOriginal = FALSE, plotProcess =FALSE,t1=1000,t2=1600,xy=TRUE,dir="H:/MS-DUAL",
                          SNR.Th=5,Wa=0.5, FWHM=0.5,top=20,interval=200,maxMass=7500,minMass=6000,chargemax = 10,chargemin = 1
  )`

2)  Scoring deconvolution

`Scored<- Score_deconvolution(lista, intervaLo2=10,mz1=6000,mz2=7500,intervaLo=500,loess = FALSE)`


3)  Peak assignment

`sequence<-("MDPNCSCAAGDSCTCAGSCKCKECKCTSCKKSCCSCCPVGCAKCAQGCICKGASDKCSCCA")`

`Assign_MassB(Scored,external=NULL,sequence=sequence, Metal=c("Zn"), NbindingMetal = 7, mzw=2,mono=TRUE)`


# Case 2
## MALDI-MS Intact protein

In this example, we analyzed a recombinant metallothionein-2 protein labelled by a set of nucleophiles targeting towards Cys residues by MALDI-MS in order to determine the number of chemically-labelled Cys residues and therefore, the number of Cys residues participating in metal coordination.

1) Specify localization of experimental file and protein sequence. The experimental file used here was saved from MS instrument as .mzXML file format. Metallothionein 2 was labelled by ethyl iodoacetate and directly analyzed by MALDI-MS.


`setwd("Directory where is located the file")`

`file<-"apomt_50mM_ET_60min.mzXML"` 

`MT2a<-("MDPNCSCAAGDSCTCAGSCKCKECKCTSCKKSCCSCCPVGCAKCAQGCICKGASDKCSCCA")`


2) Theoretical MALDI-MS pattern spectra is constructed selecting the labelling reagent.

`TheoreticalMass<-Protein.to.List(MT2a, ET= TRUE)`

3) Pre-processing experimental MALDI-MS spectra

`ProcessMALDIMS<-ProcessMALDIMS(file,t1=6000,t2=7000,plotProcess=TRUE,plotDetection=TRUE,plotOriginal = TRUE)`

4) Annotate experimental spectra by cross-correlation with simulated MALDI-MS pattern. 

`maldiTOF<-Target(ProcessMALDIMS,TheoreticalMass, mzw = 2,orig=TRUE,write=FALSE)`

## MALDI-MS Peptide-mass fingerprinting

The peptide mixture analysed experimentally is identified by peptide matching with a comprehensive peptide database constructed in silico. Focused on metalloprotein studies, the software considers all combinations for labelling of Cys residues. We might determine protein regions where metal ion binds with analyzing a set of independent samples (e.g. Zn0MT2, Zn1MT2…Zn7MT2), enzymatically digested and analysed by MALDI-TOF-MS. Comparing how the number of modifications change upon metal ion addition, one might determine the number of Cys residues that participates in metal ion coordination and determine the region of the protein.

1) Specify localization of experimental file and protein sequence.


`setwd("Directory where is located the file")`

`MT2<-"0_ZnMT2_PMF.mzXML"`

`sequence<-("MDPNCSCAAGDSCTCAGSCKCKECKCTSCKKSCCSCCPVGCAKCAQGCICKGASDKCSCCA")`


2) Theoretical PMF MALDI-MS pattern spectrum:

`TheoreticalMass<-Protein.to.Peptide(sequence)`

3) Pre-processing experimental MALDI-MS spectra

`ProcessMALDIMS<-ProcessMALDIMS(MT2,t1=850,t2=3500,amp.Th=0.008,SNR.Th=2,plotProcess=FALSE,plotDetection=FALSE,plotOriginal = FALSE,comb = "comb1")`

4) Annotate experimental spectra by cross-correlation with simulated MALDI-MS pattern by a mass error of 1 Da and a threshold of 10 % relative intensity.

`maldiTOF<-annotationMALDI_Bot(save.name=c("0_MT2_PMF"),ProcessMALDIMS,TheoreticalMass,mzw = 2,plot=FALSE,plotIndi=FALSE,orig=TRUE,write=TRUE,In=1)`

5) Score the annotated PMF.

`maldiTOF.Amp(sequence,maldiTOF,missed=5)`

`ScorePMF<-ScoresPMF(sequence,maldiTOF,TheoreticalMass,Total=61,Mw=6.042,save.name=c("0_MT2_PMF"),write=TRUE,missed=5)`

`pmfscore<-as.numeric(ScorePMF[1,5])`
`threshold<-  pmfscore`
`Mw=6042`
`bootstrapping_PMF<-PMF_bootstrapping(data=maldiTOF,R=999,sequence,TheoreticalMass,Mw=6.042,ScorePMF,missed=5)`
`Prob_score<-Prob_score(maldiTOF,TheoreticalMass,mzw=2,t2=3500,t1=500)`

5) Validate the score

`Validate_score<-Validate_score(ProcessMALDIMS,TheoreticalMass,mzw=2,nAA=61,threshold, Mw=6.042,sequence=sequence,Npermu=100,
In=10,prob=FALSE,score_prob_threshold=Prob_score,missed=5)`

`Distribution_plot(pmfscore_totalb=Validate_score,threshold=threshold,binwidth=1,sd=7,name="Zn7MT3",do.plotB=FALSE)`


# Case 3
## Bottom-up MS/MS

1) Specify localization of experimental files and protein sequence.

`Download the files: 1547.mzML, 1561.mzML, 1567.mzML`
`setwd("Directory where is located the set of MS/MS files")`

`sequence<-("MDPNCSCAAGDSCTCAGSCKCKECKCTSCKKSCCSCCPVGCAKCAQGCICKGASDKCSCCA")`

2) Pre-processing and annotation of experimental spectra

`target<-msms.batch(sequence,mzML = TRUE,label="target",missed=2,tol=2, Score=1)`


3) Validation with a target-decoy approach 


`df <- data.frame(matrix(unlist(data.frame(target)), nrow=length(target), byrow=T), stringsAsFactors=FALSE)`
`colnames(df) <- c("Sequence","Mean Error (Da)","N. ions matched","Per. peaks matched","Per. ions matched","Int. ions matched","Score","Parent ion")`

`df_full<-cbind(df,label="target")`

`full_target_decoy<-decoy(sequence,Npermu=1)`


## Native top-down MS/MS


1) Specify localization of experimental files and protein sequence.

`Download the file: MT1E_7_B_CID000005.mzML`
`setwd("Directory where is located the set of MS/MS files")`

`sequence<-("MDPNCSCAAGDSCTCAGSCKCKECKCTSCKKSCCSCCPVGCAKCAQGCICKGASDKCSCCA")`

2) Pre-processing and annotation of a set of experimental spectra

`results<-Topdown_batch(mzML=TRUE,sequence=sequence,peak=FALSE,SNR.Th=0.1,t=3,b=0,mzw=3,Metal=c("Zn"),NbindingMetal=7,maxMass = 8000,Maxcharge = 6)`


3) Validation of a particular MS/MS spectra

`ev200<-Topdown(experimental = exp,tol=0.15, sequence=sequence,peak=TRUE,SNR.Th=10,t=2,b=0,parent=50,mzw=3,Metal=c("Zn"),NbindingMetal=7,maxMass = 7000,Maxcharge = 6)`

`ev200results<-ev200[[1]]`
`ev200_true<-ev200results[!duplicated(ev200results$expt_mz),]`

`Validation_topdown(ev200_true,sequence,Deconv_CysMpro = ev200[[2]],mzw=2, intervalmz=3000,Npermu=1000,  Nmetal=7,Metal=c("Zn"))`




Please, for debugging or any other question related contact us: manuel.perisdiaz@uwr.edu.pl





