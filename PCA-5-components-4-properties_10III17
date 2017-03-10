---
title: "Machine learning: 5 components 4 properties after PCA.rmd"
author: "Mikhail Orlov"
output: pdf_document
---

Clearing workspace, setting working directory and loading R packages
```{r}
rm(list=ls())

setwd('/home/mikhail/Documents/Script_2016_all/PCA/')

library(R.matlab)
library(seqinr)
library(factoextra)
library(data.table)
library(caret)
library(doMC)
library(Biostrings)

library(reldna)
```
Reading in E. coli K12 genome (GenBank Accesion U00096.2), creating reverse complement genome and tranformation genome into character form 
```{r}
e.coli_U00096.2<-unlist(read.fasta('e.coli_U00096.2.fasta', seqonly = T))
reverseComplement_e.coli_U00096.2<-as.character(reverseComplement(DNAString(e.coli_U00096.2)))

e.coli_U00096.2_char<-unlist(strsplit(e.coli_U00096.2, ''))
reverseComplement_e.coli_U00096.2_char<-unlist(strsplit(reverseComplement_e.coli_U00096.2, ''))

```

Dynamical properties calculation from within R using octave script translated according to Grinevich et al., 2013
Creating function for dynamical properties calculation

```{r}



dynchars<-function(seq, interval_size) {
  if (missing(seq))
    stop("Need to specify sequence (as a vector of chars)")
  
  if (missing(interval_size))
    stop("Need to specify interval size")
  
  if(!is.character(seq))
    stop("Sequence must be a character vector containing A, C, G, T letters only")
  
seq<-toupper(seq)
seq<-c(seq, seq[2:(interval_size)])
  
a<-3.4*10^(-10)
I<-c(7.6, 4.8, 8.2, 4.1)*10^(-44)
K<-c(227, 155, 220, 149)*10^(-20)
V<-c(2.09, 1.43, 3.12, 2.12)*10^(-20)
tau<-c(127, 99, 140, 84)
  
csA<-cumsum(seq=='A') 
csT<-cumsum(seq=='T')
csG<-cumsum(seq=='G')
csC<-cumsum(seq=='C')

countA = csA[interval_size:length(csA)]-c(0, csA[1:(length(csA)-interval_size)])
countT = csT[interval_size:length(csT)]-c(0, csT[1:(length(csT)-interval_size)])
countG = csG[interval_size:length(csG)]-c(0, csG[1:(length(csG)-interval_size)])
countC = csC[interval_size:length(csC)]-c(0, csC[1:(length(csC)-interval_size)])
  
M<-cbind(countA, countT, countG, countC)/interval_size
M_comp<-cbind(countT, countA, countC, countG)/interval_size
M_comp<-apply(t(M_comp),1,rev) 
Is<-as.numeric(M%*%I)#! numeric conversion
Ks<-as.numeric(M%*%K)
Vs<-as.numeric(M%*%V)
  
E01<-(8*(Ks*Vs)^0.5)* 6E23 / 4184
d1<-((Ks*a^2)/Vs)^(0.5)/a;
c1<-(Ks*a^2/Is)^0.5
m1<-E01/c1/6.011E-26
taus1<-as.numeric(M%*%tau) #!as.numeric conversion
gc1 = M[,3] + M[,4]
  
Is<-as.numeric(M%*%I)#! numeric conversion
Ks<-as.numeric(M%*%K)
Vs<-as.numeric(M%*%V)

  
E02<- 8*(Ks*Vs)^0.5  * 6E23 / 4184;
d2<-((Ks*a^2)/Vs)^(0.5)/a;
c2<-(Ks*a^2/Is)^0.5;
m2<-E02/c2/6.011E-26;
taus2<-as.numeric(M_comp%*%tau)
gc2 = M_comp[,3] + M_comp[,4]

dynchars_return<-list(E01=E01, d1=d1, c1=c1, m1=m1, taus1=taus1, gc1=gc1, E02=E02, d2=d2, c2=c2, m2=m2, taus2=taus2, gc2=gc2)
  
return(dynchars_return)
  
}
```
Calculation the properties for a given genome using sliding window 200 nt
```{r}

dynchars_output<-dynchars(e.coli_U00096.2_char, 200)
E01<-dynchars_output$E01
E02<-dynchars_output$E02

d1<-dynchars_output$d1
d2<-dynchars_output$d2

gc200matlab1<-dynchars_output$gc1 #name of the variable is from before
gc200matlab2<-dynchars_output$gc2 #name of the variable is from before
```
Loading data on sequences of different types (promoters, non-promoters, genes, islands, and lowscore) from .Rdata files (must be copied separetely)
```{r}

load('spline_dataset_pro.Rdata')
load('spline_dataset_notpro.Rdata')
load('spline_dataset_gen.Rdata')
load('spline_dataset_isl.Rdata')
load('dataset_lowscore.Rdata')

##dataset_notpro[[1488]]<-NULL # to close to the left flank
```


Check-up for reverse strand
```{r, echo=FALSE}
gregexpr(toupper(dataset_isl$isl2228$seq), reverseComplement_e.coli_U00096.2)
nchar(e.coli_U00096.2)-dataset_isl$isl2228$tss-270
```
Extracting data on all promoters and on experimentaly found ones - including previosely calculated electrostatic potential profiles
```{r}
pro_names<-names(dataset_pro)
tsss<-c()
seqs<-c()
exp_strands<-c()
exp_tsss<-c()
exp_mpots_check<-c() #for experimentally found one
exp_names<-c()

for (i in 1:length(dataset_pro)){
  tsss<-c(tsss, dataset_pro[[i]]$tss)
  seqs<-c(seqs, dataset_pro[[i]]$seq)
  if (dataset_pro[[i]]$evidence=='experimental'){ 
    #exp_proms<-rbind(exp_proms, c(pro_names[i], strands[i], tsss[i], seqs[i]))
    exp_names<-c(exp_names, pro_names[i])
    exp_mpots_check<-rbind(exp_mpots_check, dataset_pro[[i]]$mpot)  	
    #exp_promoters<-rbind(exp_promoters, c(pro_names[i], strands[i], tsss[i], seqs[i]))
    exp_tsss<-c(exp_tsss, dataset_pro[[i]]$tss) 
    exp_strands<-c(exp_strands, dataset_pro[[i]]$strand)}
}
```

# # # Changing genomic coordinate used for dataset_... - TSS is according to 5'-end for both strands 
```{r}

exp_tsss[which(exp_strands=='reverse')]<-nchar(e.coli_U00096.2)-exp_tsss[which(exp_strands=='reverse')]
#exp_tsss<--exp_tsss
```


Creating matrices for data on activation energy ('aeos1') and size ('aeos3'). The matrices is 699*201 - 699 promoters sequences, 201 value for a physical property profile
to 5'-end for both strands 


Creating matrices of dynamical properties and GC-content for intervals [-150; 50] nts
```{r} 
for (i in c('aeos1forward', 'aeos1reverse', 'aeos3forward', 'aeos3reverse', 'gc200forward', 'gc200reverse', 'mpotsforward', 'mpotsreverse')) {
  assign(i, c())
}
zout<- -480:239


for (i in 1:length(exp_tsss)) {
  if (exp_strands[i]=='forward') {
    aeos1forward<-rbind(aeos1forward, E01[as.numeric(exp_tsss[i]-150):(exp_tsss[i]+50)])
    #aeos2forward<-rbind(aeos2forward, matr1[(as.numericexp_tsss[i]-150):(as.numericexp_tsss[i]+50),2])
    aeos3forward<-rbind(aeos3forward, d1[(as.numeric(exp_tsss[i])-150):(exp_tsss[i]+50)])
    #aeos4forward<-rbind(aeos4forward, matr1[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50),4])            
    gc200forward<-rbind(gc200forward, gc200matlab1[(as.numeric(exp_tsss[i])-150):(as.numeric(exp_tsss[i])+50)])
    
    p<-lseqspline1D(substr(e.coli_U00096.2, exp_tsss[i]-250, exp_tsss[i]+150), bound=c(50, 350), ref=251 )
    mpotsforward<-rbind(mpotsforward, 
                        ##mpot_whole_genome_forward$mpot[mpot_whole_genome_forward$x %in% (exp_tsss[i]-150):(exp_tsss[i]+50)])
                        p$mpot[p$x %in% zout])
    
     } else {
    aeos1reverse<-rbind(aeos1reverse, E02[(as.numeric(exp_tsss[i])-50):(as.numeric(exp_tsss[i])+150)])
    #aeos2reverse<-rbind(aeos2reverse, matr1[(as.numeric(exp_tsss[i])-50):(as.numeric(exp_tsss[i])+150),2])
    aeos3reverse<-rbind(aeos3reverse, d2[(as.numeric(exp_tsss[i])-50):(as.numeric(exp_tsss[i])+150)])
    #aeos4reverse<-rbind(aeos4reverse, matr1[(as.numeric(exp_tsss[i])-50):(as.numeric(exp_tsss[i])+150),4])
    gc200reverse<-rbind(gc200reverse, (gc200matlab2[(as.numeric(exp_tsss[i])-50):(as.numeric(exp_tsss[i])+150)]))
    
    
    
    p<-lseqspline1D(substr(reverseComplement_e.coli_U00096.2, exp_tsss[i]-250, exp_tsss[i]+150), bound=c(50, 350), ref=251 )
    mpotsreverse<-rbind(mpotsreverse, 
                        ##mpot_whole_genome_forward$mpot[mpot_whole_genome_forward$x %in% (exp_tsss[i]-150):(exp_tsss[i]+50)])
                        p$mpot[p$x %in% zout])
    
    
    }
}


#merging matrices for forward and reverse strands  together
aeos1<-rbind(aeos1forward, aeos1reverse)
#aeos2<-rbind(aeos2forward, aeos2reverse)
aeos3<-rbind(aeos3forward, aeos3reverse)
#aeos4<-rbind(aeos4forward, aeos4reverse)
gc200<-rbind(gc200forward, gc200reverse)

mpots<-rbind(mpotsforward, mpotsreverse)
```


extracting data for non-promoters. 1488 th sequences remove since it is 188 nts from flank and EP cannot be calculated
```{r}


notpro_names<-names(dataset_notpro)
nottsss<-rep(0, length(dataset_notpro))
notseqs<-c()
notstrands<-c()
#notmpots<-c()

for (i in 1:length(dataset_notpro)){
  nottsss[i]<-dataset_notpro[[i]]$tss
  notseqs<-toupper(c(notseqs, dataset_notpro[[i]]$seq))
  notstrands<-c(notstrands, dataset_notpro[[i]]$strand)
 # notmpots<-rbind(notmpots,  dataset_notpro[[i]]$mpot)  
}

# # # no reverse strand!
#nottsss[which(notstrands=='reverse')]<-nchar(e.coli_U00096.2)-nottsss[which(notstrands=='reverse')]

#nottsss<-nchar(e.coli_U00096.2)-nottsss
```
Matrices for physical properties profiles - creating matrices for data on activation energy ('aeos1') and size ('aeos3')
```{r}

for (i in c('notaeos1forward', 'notaeos1reverse', 'notaeos3forward', 'notaeos3reverse', 'notgc200forward', 'notgc200reverse', 'notmpotsforward', 'notmpotsreverse')) {
  assign(i, c())
}

for (i in 1:length(nottsss)) {
  if (notstrands[i]=='forward') {
                if (nottsss[i]<250) {
                  notaeos1forward<-rbind(notaeos1forward, E01[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50)])
    #notaeos2forward<-rbind(notaeos2forward, matr1[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50),2])
    notaeos3forward<-rbind(notaeos3forward, d1[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50)])
    #notaeos4forward<-rbind(notaeos4forward, matr1[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50),4])			
    notgc200forward<-rbind(notgc200forward, gc200matlab1[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50)])
  
                  
                spec_for1488<-substr(
                paste0(substr(e.coli_U00096.2, nchar(e.coli_U00096.2)-199, nchar(e.coli_U00096.2)), e.coli_U00096.2), 
                nottsss[i]-250+200, nottsss[i]+150+200)
                
                
                p<-lseqspline1D(spec_for1488, bound=c(50, 350), ref=251 )
                notmpotsforward<-rbind(notmpotsforward, 
                p$mpot[p$x %in% zout])
  
                } else {
    notaeos1forward<-rbind(notaeos1forward, E01[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50)])
    #notaeos2forward<-rbind(notaeos2forward, matr1[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50),2])
    notaeos3forward<-rbind(notaeos3forward, d1[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50)])
    #notaeos4forward<-rbind(notaeos4forward, matr1[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50),4])			
    notgc200forward<-rbind(notgc200forward, gc200matlab1[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50)])
  
    p<-lseqspline1D(substr(e.coli_U00096.2, nottsss[i]-250, nottsss[i]+150), bound=c(50, 350), ref=251 )
    notmpotsforward<-rbind(notmpotsforward, 
                        ##mpot_whole_genome_forward$mpot[mpot_whole_genome_forward$x %in% (exp_tsss[i]-150):(exp_tsss[i]+50)])
                        p$mpot[p$x %in% zout])
  }
        #    } else {
  #  notaeos1reverse<-rbind(notaeos1reverse, E02[(as.numeric(nottsss[i])-50):(as.numeric(nottsss[i])+150)])
    #notaeos2reverse<-rbind(notaeos2reverse, matr1[(as.numeric(nottsss[i])-50):(as.numeric(nottsss[i])+150),2])
   # notaeos3reverse<-rbind(notaeos3reverse, d2[(as.numeric(nottsss[i])-50):(as.numeric(nottsss[i])+150)])
    #notaeos4reverse<-rbind(notaeos4reverse, matr1[(as.numeric(nottsss[i])-50):(as.numeric(nottsss[i])+150),4])
   # notgc200reverse<-rbind(notgc200reverse, (gc200matlab2[(as.numeric(nottsss[i])-50):(as.numeric(nottsss[i])+150)]))
    
    
  #  p<-lseqspline1D(substr(reverseComplement_e.coli_U00096.2, nottsss[i]-250, nottsss[i]+150), bound=c(50, 350), ref=251 )
   # notmpotsreverse-rbind(notmpotsreverse, 
                        ##mpot_whole_genome_forward$mpot[mpot_whole_genome_forward$x %in% (exp_tsss[i]-150):(exp_tsss[i]+50)])
                #        p$mpot[p$x %in% zout])
    
    
    
    }
}

#mergind data for strands is not neaded - reverse strand is empty
notaeos1<-notaeos1forward
#notaeos2<-notaeos2forward
notaeos3<-notaeos3forward
#notaeos4<-notaeos4forward
notgc200<-notgc200forward

notmpots<-notmpotsforward
```
Data extracion for genes

```{r}
gen_names<-names(dataset_gen)
gentsss<-c()
genseqs<-c
genstrands<-c()
#genmpots<-c()

for (i in 1:length(dataset_gen)){
  gentsss<-c(gentsss, dataset_gen[[i]]$tss)
  genseqs<-toupper(c(genseqs, dataset_gen[[i]]$seq))
  genstrands<-c(genstrands, dataset_gen[[i]]$strand)
 # genmpots<-rbind(genmpots, dataset_gen[[i]]$mpot)
}
```

```{r}
# # #
gentsss[which(genstrands=='reverse')]<-nchar(e.coli_U00096.2)-gentsss[which(genstrands=='reverse')]


for (i in c('genaeos1forward', 'genaeos1reverse', 'genaeos3forward', 'genaeos3reverse', 'gengc200forward', 'gengc200reverse', 'genmpotsforward', 'genmpotsreverse' )) {
  assign(i, c())
}

for (i in 1:length(dataset_gen)) {
  if (genstrands[i]=='forward') {
    genaeos1forward<-rbind(genaeos1forward, E01[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50)])
    #   genaeos2forward<-rbind(genaeos2forward, matr1[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50),2])
    genaeos3forward<-rbind(genaeos3forward, d1[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50)])
    #  genaeos4forward<-rbind(genaeos4forward, matr1[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50),4])            
    gengc200forward<-rbind(gengc200forward, gc200matlab1[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50)])
    
    
   p<-lseqspline1D(substr(e.coli_U00096.2, gentsss[i]-250, gentsss[i]+150), bound=c(50, 350), ref=251 )
    genmpotsforward<-rbind(genmpotsforward, 
                        ##mpot_whole_genome_forward$mpot[mpot_whole_genome_forward$x %in% (exp_tsss[i]-150):(exp_tsss[i]+50)])
                        p$mpot[p$x %in% zout])
    
        } else {
    genaeos1reverse<-rbind(genaeos1reverse, E02[(as.numeric(gentsss[i])-50):(as.numeric(gentsss[i])+150)])
    # genaeos2reverse<-rbind(genaeos2reverse, matr1[(as.numeric(gentsss[i])-50):(as.numeric(gentsss[i])+150),2])
    genaeos3reverse<-rbind(genaeos3reverse, d2[(as.numeric(gentsss[i])-50):(as.numeric(gentsss[i])+150)])
    #genaeos4reverse<-rbind(genaeos4reverse, matr1[(as.numeric(gentsss[i])-50):(as.numeric(gentsss[i])+150),4])
    gengc200reverse<-rbind(gengc200reverse, (gc200matlab2[(as.numeric(gentsss[i])-50):(as.numeric(gentsss[i])+150)]))
   
    
    p<-lseqspline1D(substr(reverseComplement_e.coli_U00096.2, gentsss[i]-250, gentsss[i]+150), bound=c(50, 350), ref=251 )
      genmpotsreverse<-rbind(genmpotsreverse, 
                        ##mpot_whole_genome_forward$mpot[mpot_whole_genome_forward$x %in% (exp_tsss[i]-150):(exp_tsss[i]+50)])
                        p$mpot[p$x %in% zout]) 
    }
}

genaeos1<-rbind(genaeos1forward, genaeos1reverse)
#genaeos2<-rbind(genaeos2forward, genaeos2reverse)
genaeos3<-rbind(genaeos3forward, genaeos3reverse)
#genaeos4<-rbind(genaeos4forward, genaeos4reverse)
gengc200<-rbind(gengc200forward, gengc200reverse)

genmpots<-rbind(genmpotsforward, genmpotsreverse)
```
Islands. data extracion
```{r}

isl_names<-names(dataset_isl)
isltsss<-c()
islseqs<-c
islstrands<-c()
#islmpots<-c()

for (i in 1:length(dataset_isl)){
  isltsss<-c(isltsss, dataset_isl[[i]]$tss)
  islseqs<-toupper(c(islseqs, dataset_isl[[i]]$seq))
  islstrands<-c(islstrands, dataset_isl[[i]]$strand)
 # islmpots<-rbind(islmpots, dataset_isl[[i]]$mpot)
}

# # #

isltsss[which(islstrands=='reverse')]<-nchar(e.coli_U00096.2)-isltsss[which(islstrands=='reverse')]
```

```{r}

for (i in c('islaeos1forward', 'islaeos1reverse', 'islaeos3forward', 'islaeos3reverse', 'islgc200forward', 'islgc200reverse', 'islmpotsforward', 'islmpotsreverse')) {
  assign(i, c())
}

for (i in 1:length(dataset_isl)) {
  if (islstrands[i]=='forward') {
    islaeos1forward<-rbind(islaeos1forward, E01[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50)])
    #   islaeos2forward<-rbind(islaeos2forward, matr1[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50),2])
    islaeos3forward<-rbind(islaeos3forward, d1[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50)])
    #  islaeos4forward<-rbind(islaeos4forward, matr1[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50),4])            
    islgc200forward<-rbind(islgc200forward, gc200matlab1[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50)])
    
    p<-lseqspline1D(substr(e.coli_U00096.2, isltsss[i]-250, isltsss[i]+150), bound=c(50, 350), ref=251 )
    islmpotsforward<-rbind(islmpotsforward, 
                        ##mpot_whole_genome_forward$mpot[mpot_whole_genome_forward$x %in% (exp_tsss[i]-150):(exp_tsss[i]+50)])
                        p$mpot[p$x %in% zout]) 
    
    } else {
    islaeos1reverse<-rbind(islaeos1reverse, E02[(as.numeric(isltsss[i])-50):(as.numeric(isltsss[i])+150)])
    # islaeos2reverse<-rbind(islaeos2reverse, matr1[(as.numeric(isltsss[i])-50):(as.numeric(isltsss[i])+150),2])
    islaeos3reverse<-rbind(islaeos3reverse, d2[(as.numeric(isltsss[i])-50):(as.numeric(isltsss[i])+150)])
    #islaeos4reverse<-rbind(islaeos4reverse, matr1[(as.numeric(isltsss[i])-50):(as.numeric(isltsss[i])+150),4])
    islgc200reverse<-rbind(islgc200reverse, (gc200matlab2[(as.numeric(isltsss[i])-50):(as.numeric(isltsss[i])+150)]))
  
    p<-lseqspline1D(substr(reverseComplement_e.coli_U00096.2, isltsss[i]-250, isltsss[i]+150), bound=c(50, 350), ref=251 )
    islmpotsreverse<-rbind(islmpotsreverse, 
                        ##mpot_whole_genome_forward$mpot[mpot_whole_genome_forward$x %in% (exp_tsss[i]-150):(exp_tsss[i]+50)])
                        p$mpot[p$x %in% zout]) 
    
    }
}

islaeos1<-rbind(islaeos1forward, islaeos1reverse)
#islaeos2<-rbind(islaeos2forward, islaeos2reverse)
islaeos3<-rbind(islaeos3forward, islaeos3reverse)
#islaeos4<-rbind(islaeos4forward, islaeos4reverse)
islgc200<-rbind(islgc200forward, islgc200reverse)

islmpots<-rbind(islmpotsforward, islmpotsreverse)
```
Lowscore. data extraction
only first 2000 among lowscore sequences are used!
```{r}
lowscore_names<-names(dataset_lowscore)
lowscoretsss<-c()
lowscoreseqs<-c
lowscorestrands<-c()
#lowscorempots<-c()

```

```{r}
for (i in 1:2000){
  lowscoretsss<-c(lowscoretsss, dataset_lowscore[[i]]$tss)
  lowscoreseqs<-toupper(c(lowscoreseqs, dataset_lowscore[[i]]$seq))
  lowscorestrands<-c(lowscorestrands, dataset_lowscore[[i]]$strand)
 # lowscorempots<-rbind(lowscorempots, dataset_lowscore[[i]]$mpot)
}

# # #

lowscoretsss[which(lowscorestrands=='reverse')]<-nchar(e.coli_U00096.2)-lowscoretsss[which(lowscorestrands=='reverse')]
```


```{r}

for (i in c('lowscoreaeos1forward', 'lowscoreaeos1reverse', 'lowscoreaeos3forward', 'lowscoreaeos3reverse', 'lowscoregc200forward', 'lowscoregc200reverse', 'lowscorempotsforward', 'lowscorempotsreverse' )) {
  assign(i, c())
}

```

only first 2000 among lowscore sequences are used!
```{r}

for (i in 1:2000) {
  if (lowscorestrands[i]=='forward') {
    lowscoreaeos1forward<-rbind(lowscoreaeos1forward, E01[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50)])
    #   lowscoreaeos2forward<-rbind(lowscoreaeos2forward, matr1[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50),2])
    lowscoreaeos3forward<-rbind(lowscoreaeos3forward, d1[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50)])
    #  lowscoreaeos4forward<-rbind(lowscoreaeos4forward, matr1[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50),4])            
    lowscoregc200forward<-rbind(lowscoregc200forward, gc200matlab1[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50)])

      p<-lseqspline1D(substr(e.coli_U00096.2, lowscoretsss[i]-250, lowscoretsss[i]+150), bound=c(50, 350), ref=251 )
    lowscorempotsforward<-rbind(lowscorempotsforward, 
                        ##mpot_whole_genome_forward$mpot[mpot_whole_genome_forward$x %in% (exp_tsss[i]-150):(exp_tsss[i]+50)])
                        p$mpot[p$x %in% zout]) 
    
      } else {
    lowscoreaeos1reverse<-rbind(lowscoreaeos1reverse, E02[(as.numeric(lowscoretsss[i])-50):(as.numeric(lowscoretsss[i])+150)])
    # lowscoreaeos2reverse<-rbind(lowscoreaeos2reverse, matr1[(as.numeric(lowscoretsss[i])-50):(as.numeric(lowscoretsss[i])+150),2])
    lowscoreaeos3reverse<-rbind(lowscoreaeos3reverse, d2[(as.numeric(lowscoretsss[i])-50):(as.numeric(lowscoretsss[i])+150)])
    #lowscoreaeos4reverse<-rbind(lowscoreaeos4reverse, matr1[(as.numeric(lowscoretsss[i])-50):(as.numeric(lowscoretsss[i])+150),4])
    lowscoregc200reverse<-rbind(lowscoregc200reverse, (gc200matlab2[(as.numeric(lowscoretsss[i])-50):(as.numeric(lowscoretsss[i])+150)]))
  
         p<-lseqspline1D(substr(reverseComplement_e.coli_U00096.2, lowscoretsss[i]-250, lowscoretsss[i]+150), bound=c(50, 350), ref=251 )
    lowscorempotsreverse<-rbind(lowscorempotsreverse, 
                        ##mpot_whole_genome_forward$mpot[mpot_whole_genome_forward$x %in% (exp_tsss[i]-150):(exp_tsss[i]+50)])
                        p$mpot[p$x %in% zout]) 
    
    }
}

lowscoreaeos1<-rbind(lowscoreaeos1forward, lowscoreaeos1reverse)
#lowscoreaeos2<-rbind(lowscoreaeos2forward, lowscoreaeos2reverse)
lowscoreaeos3<-rbind(lowscoreaeos3forward, lowscoreaeos3reverse)
#lowscoreaeos4<-rbind(lowscoreaeos4forward, lowscoreaeos4reverse)
lowscoregc200<-rbind(lowscoregc200forward, lowscoregc200reverse)

lowscorempots<-rbind(lowscorempotsforward, lowscorempotsreverse)

```
Setting names for created datasets 
```{r}

rownames(aeos1)<-exp_names
rownames(aeos3)<-exp_names
rownames(mpots)<-exp_names
rownames(gc200)<-exp_names

rownames(notaeos1)<-paste0('Non_promoter_', 1:nrow(notaeos1))
rownames(notaeos3)<-paste0('Non_promoter_', 1:nrow(notaeos1))
rownames(notmpots)<-paste0('Non_promoter_', 1:nrow(notaeos1))
rownames(notgc200)<-paste0('Non_promoter_', 1:nrow(notaeos1))

rownames(genaeos1)<-paste0('Gene_', 1:nrow(genaeos1))
rownames(genaeos3)<-paste0('Gene_', 1:nrow(genaeos1))
rownames(genmpots)<-paste0('Gene_', 1:nrow(genaeos1))
rownames(gengc200)<-paste0('Gene_', 1:nrow(genaeos1))

rownames(islaeos1)<-paste0('Islands_', 1:nrow(islaeos1))
rownames(islaeos3)<-paste0('Islands_', 1:nrow(islaeos1))
rownames(islmpots)<-paste0('Islands_', 1:nrow(islaeos1))
rownames(islgc200)<-paste0('Islands_', 1:nrow(islaeos1))

rownames(lowscoreaeos1)<-paste0('Lowscore_', 1:nrow(lowscoreaeos1))
rownames(lowscoreaeos3)<-paste0('Lowscore_', 1:nrow(lowscoreaeos1))
rownames(lowscorempots)<-paste0('Lowscore_', 1:nrow(lowscoreaeos1))
rownames(lowscoregc200)<-paste0('Lowscore_', 1:nrow(lowscoreaeos1))


```

#Merging datasets. Scaling data
```{r}


to_pca_5components_4props<-rbind(cbind((aeos1), (aeos3), (mpots), (gc200)),
                                 cbind((notaeos1), (notaeos3), (notmpots), (notgc200)),
                                 cbind((genaeos1), (genaeos3), (genmpots), (gengc200)),
                                 cbind((islaeos1), (islaeos3), (islmpots), (islgc200)),
                                 cbind((lowscoreaeos1), (lowscoreaeos3), (lowscorempots), (lowscoregc200))
)
```
# saving the initial matrix for 5 components, 4 properties
```{r}
save(to_pca_5components_4props, file='to_pca_5components_4props.Rdata')
```

#setting sequences groups and variables groups 
```{r}
habillage_5components_4props<-c(rep('Promoters', length(exp_names)), rep('Non-promoters', length(dataset_notpro)), rep('Genes', length(dataset_gen)), rep('Islands', length(dataset_isl)), rep('Lowscore', 2000))

habillage_5components_4props_cols<-c(rep('red', length(exp_names)), rep('blue', length(dataset_notpro)), rep('green', length(dataset_gen)), rep('yellow', length(dataset_isl)), rep('magenta', 2000))

colnames(to_pca_5components_4props)<-NULL

groups_variables <- c(rep('cyan', ncol(aeos1)), rep('navy', ncol(aeos3)), rep('darkgreen', ncol(mpots)), rep('darkred', ncol(gc200)))
```




#PRINCIPAL COMPONENTS ANALYSIS ITSELF
```{r}

princ.return.5comps.4props <- prcomp(to_pca_5components_4props, scale=T, center = T)


```

``` {r PCA variables plot using base graphics}
# Helper function : 
# Correlation between variables and principal components
var_cor_func <- function(var.loadings, comp.sdev){
  var.loadings*comp.sdev
  }
# Variable correlation/coordinates
loadings <- princ.return.5comps.4props$rotation
sdev <- princ.return.5comps.4props$sdev
var.coord <- t(apply(loadings, 1, var_cor_func, sdev))

#Graph of variables using R base graph



#INDIVIDUALS
# Contributions of individuals
ind.coord <- princ.return.5comps.4props$x
contrib <- function(ind.coord, comp.sdev, n.ind){
  100*(1/n.ind)*ind.coord^2/comp.sdev^2
}
ind.contrib <- t(apply(ind.coord,1, contrib, 
                       princ.return.5comps.4props$sdev, nrow(ind.coord)))
svg('/home/mikhail/Documents/Script_2016_all/PCA/custom_biplot_pca_5_comps.svg', height = 10, width = 10)
plot(ind.coord[,1], ind.coord[,2], cex=0.8, pch=1,
     xlab="PC1",ylab="PC2", col=habillage_5components_4props_cols)
abline(h=0, v=0, lty = 2)
#text(ind.coord[,1], ind.coord[,2], labels=rownames(ind.coord),
     #   cex=0.7, pos = 3)


# Plot the correlation circle
a <- seq(0, 2*pi, length = 100)
lines( cos(a)*50, sin(a)*50, type = 'l', col="gray", #*50 !
      xlab = "PC1",  ylab = "PC2")
abline(h = 0, v = 0, lty = 2)
# Add active variables
arrows(0, 0, var.coord[, 1]*50, var.coord[, 2]*50, #*50 so the size of arrows be proper
      length = 0.1, angle = 15, code = 2, col = groups_variables)
# Add labels
#text(var.coord, labels=rownames(var.coord), cex = 1, adj=1)
legend('bottomright', legend = as.character(c('E0', 'd', 'EP', 'GC')), fill = unique(groups_variables), cex = 1.5)
legend("topright", legend=c('Promoters', 'Non-promoters', 'Genes', 'Promoter islands', 'Lowscore'), fill=unique(habillage_5components_4props_cols))
dev.off()
```


```{r 3D variables and individuals vusualization}
library(rgl)
plot3d(princ.return.5comps.4props$x[,1:3], col=habillage_5components_4props_cols, size=2, xlab = paste0('PC1, approx. ', ceiling(variance[1]), '% of variance'), ylab = paste0('PC2, approx. ', ceiling(variance[2]), '% of variance'), zlab = paste0('PC3, approx. ', ceiling(variance[3]), '% of variance'))
coords <- NULL
for (i in 1:nrow(princ.return.5comps.4props$rotation)) {
  coords <- rbind(coords, rbind(c(0,0,0),princ.return.5comps.4props$rotation[i,1:3]))
}
lines3d(coords*1000, col=groups_variables, lwd=0.2)

M <- par3d("userMatrix") # #pca3d
movie3d( spin3d(), duration=15,  dir="/home/mikhail/Documents/Script_2016_all/PCA/pca_movie3d/", clean = T, convert = T ) 



```


```{r Variables visualization}
svg('/home/mikhail/Documents/Script_2016_all/PCA/fviz_pca_var_princ.return.5comps.4props.svg', height = 10, width = 10)
fviz_pca_var(princ.return.5comps.4props, geom = c('arrow', 'text'), col.var="contrib")+
scale_color_gradient2(low="white", mid="white", 
      high="black", midpoint=0.1) + theme_minimal()
dev.off()

svg('/home/mikhail/Documents/Script_2016_all/PCA/fviz_pca_var_princ.return.5comps.4props.svg', height = 10, width = 10)
fviz_pca_biplot(princ.return.5comps.4props, 
  habillage = as.factor(habillage_5components_4props), addEllipses = TRUE,
  col.var = rainbow(ncol(to_pca_5components_4props)), alpha.var = "contrib",
  label = "var"
  ) +
  scale_color_brewer(palette="Dark2")+
  theme_minimal()

dev.off()
library(ggbiplot)
#print(ggbiplot(princ.return.5comps.4props, obs.scale = 1, var.scale = 1, groups = (habillage_5components_4props), ellipse = TRUE, circle = TRUE))

#biplot(princ.return.5comps.4props, choices = 1:2, scale = 1, pc.biplot = TRUE, col=habillage_5components_4props_cols)

```

#saving initial prcomp output 
```{r}
save(princ.return.5comps.4props, file='princ.return.5comps.4props.Rdata')

```
prcomp output transformation
```{r}
eig.val <- get_eigenvalue(princ.return.5comps.4props)

# Eigenvalues
eig <- (princ.return.5comps.4props$sdev)^2

# Variances in percentage
variance <- eig*100/sum(eig)

# Cumulative variances
cumvar <- cumsum(variance)
#png('cumvar_5components_4props.png', height=1250, width=1250, res=130)
plot(cumvar, type='l', main='Cumulative variance for principal components \n on 4 properties for 5 components mixture', ylab='Cumulative variance (%)', xlab='Principal components')
v.5comps.4props=100 # 65 
h.5comps.4props=98
abline(h=h.5comps.4props, v=v.5comps.4props, col=12)
text(v.5comps.4props-15, h.5comps.4props ,paste('Number of \nPCs=', v.5comps.4props,'\n' ,h.5comps.4props, '% of \nvariance \nretained'),srt=0.2,pos=3)

#dev.off()

```
Futher transformation - chosen number of 'rotation' columns (i.e. eigenvectors are multiplied by inital large matrix of variables)
```{r}
#pcs.from.var<-c()
#for (i in 1:v.5comps.4props) {
#  load <- princ.return.5comps.4props$rotation[,i]
#  pr.cp <- to_pca_5components_4props %*% load
#  pr <- as.numeric(pr.cp)
#  pcs.from.var<-cbind(pcs.from.var, pr)
#}
```
#saving PCA data: variable are conversed to PCs
```{r}
#rownames(pcs.from.var)<-rownames(to_pca_5components_4props)
#colnames(pcs.from.var)<-1:ncol(pcs.from.var)
```
#renaming the variable, saving the transformated output
```{r}
#mixt_5comps_after_pca_on_4_props_ae_size_mpots_gc200<- pcs.from.var
#save(mixt_5comps_after_pca_on_4_props_ae_size_mpots_gc200, file='new_no_factor_8XII_mixt_5comps_after_pca_on_4_props_ae_size_mpots_gc200.Rdata')
```
#general PCA results visualization
habillage i.e.  partitions into groups
```{r}
svg(filename = 'PCA_5_components_4_properties_not_conversion_into_variables_axes_1_2.svg', width = 10, height = 7)
fviz_pca_ind(princ.return.5comps.4props, label="none", habillage=habillage_5components_4props, pointsize = 0.5, addEllipses = T, ellipse.level = 0.99)
dev.off()

svg(filename = 'PCA_5_components_4_properties_not_conversion_into_variables_axes_1_3.svg', width = 10, height = 7)
fviz_pca_ind(princ.return.5comps.4props, axes = c(1,3), label="none", habillage=habillage_5components_4props, pointsize = 0.5, addEllipses = T, ellipse.level = 0.99)
dev.off()


svg(filename = 'PCA_5_components_4_properties_not_conversion_into_variables_axes_2_3.svg', width = 10, height = 7)
fviz_pca_ind(princ.return.5comps.4props, axes = c(2,3), label="none", habillage=habillage_5components_4props, pointsize = 0.5, addEllipses = T, ellipse.level = 0.99)
dev.off()



library(rgl)

habillage.pro.not.isl.gen<-c(rep('red', length(exp_names)), rep ('blue', length(dataset_notpro)), rep('green', length(dataset_gen)), rep('orange', length(dataset_isl)), rep('magenta', 2000)) 
#4 components correspond to colors 1-4





open3d()
plot3d(princ.return.5comps.4props$x[,1:3], col=(habillage.pro.not.isl.gen), size=2.5)

M <- par3d("userMatrix") # #pca3d
movie3d( spin3d(), duration=15, dir="/home/mikhail/Documents/Script_2016_all/PCA/pca_movie3d/", clean = T, convert = T ) 

# #
library(pca3d)
pca3d(princ.return.5comps.4props, group = habillage.pro.not.isl.gen, col = habillage.pro.not.isl.gen, radius = 0.5)
#makeMoviePCA()

```


# # # SUPERVISED MACHINE LEARNING
Creating subsets for separate sequence types
```{r}
registerDoMC(cores = 3)
df <- as.data.frame(((princ.return.5comps.4props$x)[,1:100]))
#df<-cbind(habillage_5components_4props, df)
promoters <- df[1:699,]
non_promoters <- df[700:2579,]
genes <- df[2580:6006,]
islands <- df[6007:8234,]
lowscore <- df[8235:10234,]
```
Adding factor columns
```{r}
factor_to_promoters_vs_lowscore<-as.factor(c(rep('Promoter', nrow(promoters)), rep('Lowscore', nrow(lowscore))))
promoters_vs_lowscore <- cbind(factor_to_promoters_vs_lowscore, rbind(promoters, lowscore))

factor_to_promoters_vs_non_promoters<-as.factor(c(rep('Promoter', nrow(promoters)), rep('Non_promoter', nrow(non_promoters))))
promoters_vs_non_promoters <- cbind(factor_to_promoters_vs_non_promoters, rbind(promoters, non_promoters))

factor_to_promoters_vs_islands<-as.factor(c(rep('Promoter', nrow(promoters)), rep('Island', nrow(islands))))
promoters_vs_islands <- cbind(factor_to_promoters_vs_islands, rbind(promoters, islands))

factor_to_promoters_vs_genes<-as.factor(c(rep('Promoter', nrow(promoters)), rep('Gene', nrow(genes))))
promoters_vs_genes <- cbind(factor_to_promoters_vs_genes, rbind(promoters, genes))
```

#for promoters_vs_lowscore model training and evaluation

```{r}
set.seed(999)
inTraining <- createDataPartition(promoters_vs_lowscore$factor_to_promoters_vs_lowscore, p = 0.7, list = F)
training <- promoters_vs_lowscore[inTraining,]
testing  <- promoters_vs_lowscore[-inTraining,]
fitControl <- trainControl(method = "repeatedcv", 
                           number = 10, 
                           repeats = 15, 
                           allowParallel = T, 
                           classProbs = T, 
                           summaryFunction = twoClassSummary
)
```
#train_control <- trainControl(method="repeatedcv", number=10, repeats=3)
```{r}
fit_promoters_vs_lowscore <- train(factor_to_promoters_vs_lowscore ~ .,
                                   data = training,
                                   method = "rf",#nb, lda, nbDiscrete, nbSearch
   #                                #preProcess=c("center", "scale"),
                                   tuneLength = 15,
                                   trControl = fitControl,
                                   metric = "ROC"
)

predictionClasses_promoters_vs_lowscore <- predict(fit_promoters_vs_lowscore, newdata = testing)
predictionProb_promoters_vs_lowscore <- predict(fit_promoters_vs_lowscore, newdata = testing, type ="prob")
confusionMatrix_promoters_vs_lowscore <- confusionMatrix(data = predictionClasses_promoters_vs_lowscore, testing$factor_to_promoters_vs_lowscore)

```

For promoters_vs_non_promoters
```{r}
set.seed(999)
inTraining <- createDataPartition(promoters_vs_non_promoters$factor_to_promoters_vs_non_promoters, p = 0.7, list = F)
training <- promoters_vs_non_promoters[inTraining,]
testing  <- promoters_vs_non_promoters[-inTraining,]
fitControl <- trainControl(method = "repeatedcv", 
                           number = 10, 
                           repeats = 3, 
                           allowParallel = T, 
                           classProbs = T, 
                           summaryFunction = twoClassSummary
)

#train_control <- trainControl(method="repeatedcv", number=10, repeats=3)
fit_promoters_vs_non_promoters <- train(factor_to_promoters_vs_non_promoters ~ .,
                                        data = training,
                                        method = "nb",#nb, lda, nbDiscrete, nbSearch
                                        ##preProcess=c("center", "scale"),
                                        tuneLength = 7,
                                        trControl = fitControl,
                                        metric = "ROC"
)
predictionClasses_promoters_vs_non_promoters <- predict(fit_promoters_vs_non_promoters, newdata = testing)
predictionProb_promoters_vs_non_promoters <- predict(fit_promoters_vs_non_promoters, newdata = testing, type ="prob")
confusionMatrix_promoters_vs_non_promoters <- confusionMatrix(data = predictionClasses_promoters_vs_non_promoters, testing$factor_to_promoters_vs_non_promoters)
```
#for promoters_vs_islands

```{r}

set.seed(999)
inTraining <- createDataPartition(promoters_vs_islands$factor_to_promoters_vs_islands, p = 0.7, list = F)
training <- promoters_vs_islands[inTraining,]
testing  <- promoters_vs_islands[-inTraining,]
fitControl <- trainControl(method = "repeatedcv", 
                           number = 10, 
                           repeats = 3, 
                           allowParallel = T, 
                           classProbs = T, 
                           summaryFunction = twoClassSummary
)

#train_control <- trainControl(method="repeatedcv", number=10, repeats=3)
fit_promoters_vs_islands <- train(factor_to_promoters_vs_islands ~ .,
                                  data = training,
                                  method = "nb",#nb, lda, nbDiscrete, nbSearch
                                  ##preProcess=c("center", "scale"),
                                  tuneLength = 7,
                                  trControl = fitControl,
                                  metric = "ROC"
)
predictionClasses_promoters_vs_islands <- predict(fit_promoters_vs_islands, newdata = testing)
predictionProb_promoters_vs_islands <- predict(fit_promoters_vs_islands, newdata = testing, type ="prob")
confusionMatrix_promoters_vs_islands <- confusionMatrix(data = predictionClasses_promoters_vs_islands, testing$factor_to_promoters_vs_islands)
```

#for promoters_vs_genes
```{r}
set.seed(999)
inTraining <- createDataPartition(promoters_vs_genes$factor_to_promoters_vs_genes, p = 0.7, list = F)
training <- promoters_vs_genes[inTraining,]
testing  <- promoters_vs_genes[-inTraining,]
fitControl <- trainControl(method = "repeatedcv", 
                           number = 10, 
                           repeats = 3, 
                           allowParallel = T, 
                           classProbs = T, 
                           summaryFunction = twoClassSummary
)

#train_control <- trainControl(method="repeatedcv", number=10, repeats=3)
fit_promoters_vs_genes <- train(factor_to_promoters_vs_genes ~ .,
                                data = training,
                                method = "nb",#nb, lda, nbDiscrete, nbSearch
                                #preProcess=c("center", "scale"),
                                tuneLength = 7,
                                trControl = fitControl,
                                metric = "ROC"
)
predictionClasses_promoters_vs_genes <- predict(fit_promoters_vs_genes, newdata = testing)
predictionProb_promoters_vs_genes <- predict(fit_promoters_vs_genes, newdata = testing, type ="prob")
confusionMatrix_promoters_vs_genes <- confusionMatrix(data = predictionClasses_promoters_vs_genes, testing$factor_to_promoters_vs_genes)

```

Mean values for used properties for sequences sets

```{r}
aes<-grep('*aeos1$', ls(), value = T)[c(1,5,2,3,4)]

meanaes<-c()

for (i in aes) {
  meanaes<-c(meanaes, mean(get(i)))

}

names(meanaes)<-c('Promoters', 'Non-promoters', 'Genes', 'Islands', 'Lowscore')
names(medianaes)<-c('Promoters', 'Non-promoters', 'Genes', 'Islands', 'Lowscore')


ds<-grep('*aeos3$', ls(), value = T)[c(1,5,2,3,4)]

meands<-c()
for (i in ds) {
  meands<-c(meands, mean(get(i)))
}

names(meands)<-c('Promoters', 'Non-promoters', 'Genes', 'Islands', 'Lowscore')


eps<-grep('*mpots$', ls(), value = T)[c(1,5,2,3,4)]

meaneps<-c()
for (i in eps) {
  meaneps<-c(meaneps, mean(get(i)))
}

names(meaneps)<-c('Promoters', 'Non-promoters', 'Genes', 'Islands', 'Lowscore')

gcs<-grep('*gc200$', ls(), value = T)[c(1,5,2,3,4)]

meangcs<-c()
for (i in gcs) {
  meangcs<-c(meangcs, mean(get(i)))
}

names(meangcs)<-c('Promoters', 'Non-promoters', 'Genes', 'Islands', 'Lowscore')



```
Islands cluster analisys
```{r}
library(fastcluster)
islclust<-hclust.vector(princ.return.5comps.4props$x[6007:8234,], method = 'ward')

plot(islclust, labels=F)
cutisl<-cutree(islclust, k=2)

habillage2isl.clusts<-c(rep(1, length(exp_names)), rep (1, length(dataset_notpro)), rep(1, length(dataset_gen)), rep('red', length(which(cutisl==1))), rep('blue', length(which(cutisl==2))), rep(1, 2000)) 


fviz_pca_ind(princ.return.5comps.4props, label="none", habillage=habillage2isl.clusts, pointsize = 0.5, addEllipses = T, ellipse.level = 0.99)


fviz_pca_ind(princ.return.5comps.4props, axes = c(1,3), label="none", habillage=habillage2isl.clusts, pointsize = 0.5, addEllipses = T, ellipse.level = 0.99)

fviz_pca_ind(princ.return.5comps.4props, axes = c(2,3), label="none", habillage=habillage2isl.clusts, pointsize = 0.5, addEllipses = T, ellipse.level = 0.99)

plot3d(princ.return.5comps.4props$x[,1:3], col=(habillage2isl.clusts), size=2, cex.lab=1.8)

summary(as.vector(islaeos1[which(cutisl==1),]))
summary(as.vector(islaeos1[which(cutisl==2),]))

isl_props<-grep('isl', ls(), value=T)

lmts <- range(as.vector(islaeos1[which(cutisl==1),]),
              as.vector(islaeos1[which(cutisl==2),]))

par(mfrow = c(1, 2))
boxplot(as.vector(islaeos1[which(cutisl==1),]),ylim=lmts)
boxplot(as.vector(islaeos1[which(cutisl==2),]),ylim=lmts)

lmts <- range(as.vector(islaeos3[which(cutisl==1),]),
              as.vector(islaeos3[which(cutisl==2),]))

par(mfrow = c(1, 2))
boxplot(as.vector(islaeos3[which(cutisl==1),]),ylim=lmts)
boxplot(as.vector(islaeos3[which(cutisl==2),]),ylim=lmts)

summary(as.vector(islaeos3[which(cutisl==1),]))
summary(as.vector(islaeos3[which(cutisl==2),]))

summary(as.vector(islmpots[which(cutisl==1),]))
summary(as.vector(islmpots[which(cutisl==2),]))

summary(as.vector(islgc200[which(cutisl==1),]))
summary(as.vector(islgc200[which(cutisl==2),]))

summary(isltsss[which(cutisl==1)])
summary(isltsss[which(cutisl==2)])

table(islstrands[which(cutisl==1)])
table(islstrands[which(cutisl==2)])

lmts <- range(as.vector(islaeos1[which(cutisl==1),]),
              as.vector(islaeos1[which(cutisl==2),]))
par(mfrow = c(1, 2))
matplot(islaeos1[which(cutisl==1),], type='l', lwd=0.03, col='red', ylim = lmts)
matplot(islaeos1[which(cutisl==2),], type='l', lwd=0.03, col='blue', ylim = lmts)


lmts <- range(as.vector(islaeos3[which(cutisl==1),]),
              as.vector(islaeos3[which(cutisl==2),]))
par(mfrow = c(1, 2))
matplot(islaeos3[which(cutisl==1),], type='l', lwd=0.03, col='red', ylim = lmts)
matplot(islaeos3[which(cutisl==2),], type='l', lwd=0.03, col='blue', ylim = lmts)


lmts <- range(as.vector(islmpots[which(cutisl==1),]),
              as.vector(islmpots[which(cutisl==2),]))
par(mfrow = c(1, 2))
matplot(islmpots[which(cutisl==1),], type='l', lwd=0.03, col='red', ylim = lmts)
matplot(islmpots[which(cutisl==2),], type='l', lwd=0.03, col='blue', ylim = lmts)


lmts <- range(as.vector(islgc200[which(cutisl==1),]),
              as.vector(islgc200[which(cutisl==2),]))
par(mfrow = c(1, 2))

matplot(islgc200[which(cutisl==1),], type='l', lwd=0.03, col='red', ylim = lmts)
matplot(islgc200[which(cutisl==2),], type='l', lwd=0.03, col='blue', ylim = lmts)



```


# # # UNSUPERVISED MACHINE LEARNING - CLUSTER ANALYSIS

```{r}
library(fastcluster)

ward.princ.return.5comps.4props<-hclust.vector(princ.return.5comps.4props$x[,1:100], method = 'ward')

library(NbClust)

indices_that_works<-c("kl", "ch", "hartigan", "cindex", "db", "silhouette", "ratkowsky", "ball", "ptbiserial", "frey", "mcclain",  "dunn", "sdbw")

all_for_5comp<-c()

for (i in indices_that_works) {
    print(i)
    siddres<-NbClust(data = princ.return.5comps.4props$x[,1:100], diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 22, method = 'ward.D2', index = i)
    assign(paste0('sidd_nbres_',i), siddres)
    all_for_5comp<-c(all_for_5comp, siddres) 
    print(siddres$Best.nc[1])
}



library(dendextend)


d1=color_branches(as.dendrogram(ward.princ.return.5comps.4props), k=nrow(princ.return.5comps.4props$x[,1:100]), col = habillage.pro.not.isl.gen[order.dendrogram(as.dendrogram(ward.princ.return.5comps.4props))])

d1 %>% rect.dendrogram(k = 4, horiz = TRUE, border = 8, lty = 5, lwd = 2, col='red')
png('ward_on_pca_4_prop_5_components_circlized.png', height=3000, width=3000, res=300)
raised_d1<-raise.dendrogram(d1, 7000)

#raised_d1 %>% rect.dendrogram(k = 4, horiz = TRUE, border = 8, lty = 5, lwd = 2, col='red')

circlize_dendrogram(raised_d1, labels=F, dend_track_height = 0.95)
dev.off()


png('ward_on_pca_5_prop_4_components_circlized_raised.png', height=3000, width=3000, res=300)
circlize_dendrogram(raised_d1, labels=F, dend_track_height = 0.95)
##plot(raised_d1, main="Ward's method clustering for PCA results \n(promoters, non-promoters, genes, promoter islands, and lowscore sequences", labels=F)
  legend("topright", legend=c('Promoters', 'Non-promoters', 'Genes', 'Promoter islands', 'Lowscore'), fill=unique(habillage.pro.not.isl.gen))

dev.off()

```


separately for properties
```{r}
library(fastcluster)

ward.princ.return.5comps.4props<-hclust.vector(aeos$x[,1:100], method = 'ward')

library(NbClust)

indices_that_works<-c("kl", "ch", "hartigan", "cindex", "db", "silhouette", "ratkowsky", "ball", "ptbiserial", "frey", "mcclain",  "dunn", "sdbw")

all_for_5comp<-c()

for (i in indices_that_works) {
    print(i)
    siddres<-NbClust(data = princ.return.5comps.4props$x[,1:100], diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 22, method = 'ward.D2', index = i)
    assign(paste0('sidd_nbres_',i), siddres)
    all_for_5comp<-c(all_for_5comp, siddres) 
    print(siddres$Best.nc[1])
}



library(dendextend)


d1=color_branches(as.dendrogram(ward.princ.return.5comps.4props), k=nrow(princ.return.5comps.4props$x[,1:100]), col = habillage.pro.not.isl.gen[order.dendrogram(as.dendrogram(ward.princ.return.5comps.4props))])

d1 %>% rect.dendrogram(k = 4, horiz = TRUE, border = 8, lty = 5, lwd = 2, col='red')
png('ward_on_pca_4_prop_5_components_circlized.png', height=3000, width=3000, res=300)
raised_d1<-raise.dendrogram(d1, 14000)

##plot(hang.dendrogram(raised_d1), horiz = TRUE)

#raised_d1 %>% rect.dendrogram(k = 4, horiz = TRUE, border = 8, lty = 5, lwd = 2, col='red')

circlize_dendrogram(raised_d1, labels=F, dend_track_height = 0.95)
dev.off()


png('ward_on_pca_5_prop_4_components_circlized_raised.png', height=3000, width=3000, res=300)
circlize_dendrogram(raised_d1, labels=F, dend_track_height = 0.95)
##plot(raised_d1, main="Ward's method clustering for PCA results \n(promoters, non-promoters, genes, promoter islands, and lowscore sequences", labels=F)
  legend("topright", legend=c('Promoters', 'Non-promoters', 'Genes', 'Promoter islands', 'Lowscore'), fill=unique(habillage.pro.not.isl.gen))

dev.off()

```


```{r}



#poster
XX<-(-540):(179)


plot(XX, mpots[1,], ty='l', lwd=1.5, ylab='EP value', xlab='Sequence (angstrom)', cex.lab=1.4 )
abline(v=0, lty=2)

X<--150:50
plot(X, aeos1[1,], ty='l', lwd=1.5, ylab='Activation energy (kcal/mol)', xlab='Sequence (nts)', cex.lab=1.4 )
abline(v=0, lty=2)

plot(X, aeos3[1,], ty='l', lwd=1.5, ylab='Size of open state (nts)', xlab='Sequence (nts)', cex.lab=1.4 )
abline(v=0, lty=2)

plot(X, gc200[1,], ty='l', lwd=1.5, ylab=' GC-content for 200 b.p.', xlab='Sequence (nts)', cex.lab=1.4 )
abline(v=0, lty=2)

plot(1:100, princ.return.5comps.4props$x[1,1:100], ty='l', lwd=2, xlab='PCs most contibuting to variance', ylab='Relative units', cex.lab=1.4, col='red' )

#plot 3d for 3 main PCs colored according to cluster analysis

plot3d(princ.return.5comps.4props$x[,1:3], col=habillage.pro.not.isl.gen[order.dendrogram(as.dendrogram(ward.princ.return.5comps.4props))], size=2.5)


#d1 %>% rect.dendrogram(k = 4, horiz = TRUE, border = 8, lty = 5, lwd = 2, col='red')

#confusionMatrix visualization

par(mfrow=c(2,2))
fourfoldplot(confusionMatrix_promoters_vs_lowscore$table, color = c("#CC6666", "#99CC99"),
              conf.level = 0, margin = 1, main = "Promoters vs lowscore")

fourfoldplot(confusionMatrix_promoters_vs_non_promoters$table, color = c("#CC6666", "#99CC99"),
              conf.level = 0, margin = 1, main = "Promoters vs non-promoters")

fourfoldplot(confusionMatrix_promoters_vs_genes$table, color = c("#CC6666", "#99CC99"),
              conf.level = 0, margin = 1, main = "Promoters vs genes")

fourfoldplot(confusionMatrix_promoters_vs_islands$table, color = c("#CC6666", "#99CC99"),
              conf.level = 0, margin = 1, main = "Promoters vs islands")



ward_all_aeos1<-hclust.vector(to_pca_5components_4props[, 1:ncol(aeos1)], method = 'ward')
ward_all_aeos3<-hclust.vector(to_pca_5components_4props[, 202:(202+201)], method = 'ward')
ward_all_mpots<-hclust.vector(to_pca_5components_4props[, 403:(403+720)], method = 'ward')

ward_all_gc200<-hclust.vector(to_pca_5components_4props[,(402+721):ncol(to_pca_5components_4props)], method = 'ward')



library(dendextend)

d1aeos1=color_branches(as.dendrogram(ward_all_aeos1), k=nrow(to_pca_5components_4props), col = habillage.pro.not.isl.gen[order.dendrogram(as.dendrogram(ward_all_aeos1))])
raised_d1aeos1<-raise.dendrogram(d1aeos1, 8000)

d1aeos3=color_branches(as.dendrogram(ward_all_aeos3), k=nrow(to_pca_5components_4props), col = habillage.pro.not.isl.gen[order.dendrogram(as.dendrogram(ward_all_aeos3))])
raised_d1aeos3<-raise.dendrogram(d1aeos3, 700)

d1mpots=color_branches(as.dendrogram(ward_all_mpots), k=nrow(to_pca_5components_4props), col = habillage.pro.not.isl.gen[order.dendrogram(as.dendrogram(ward_all_mpots))])
raised_d1mpots<-raise.dendrogram(d1mpots, 700)

d1gc200=color_branches(as.dendrogram(ward_all_gc200), k=nrow(to_pca_5components_4props), col = habillage.pro.not.isl.gen[order.dendrogram(as.dendrogram(ward_all_gc200))])
raised_d1gc200<-raise.dendrogram(d1gc200, 700)

svg('ward_on_pca_4_prop_5_components_circlized_raised_all_4.png', height=15, width=15)
par(mfrow=c(2,2))
circlize_dendrogram(raised_d1aeos1, labels=F, dend_track_height = 0.95)
circlize_dendrogram(raised_d1aeos3, labels=F, dend_track_height = 0.95)
circlize_dendrogram(raised_d1mpots, labels=F, dend_track_height = 0.95)
circlize_dendrogram(raised_d1gc200, labels=F, dend_track_height = 0.95)
dev.off()


png('ward_on_pca_4_prop_4_components_circlized.png', height=3000, width=3000, res=300)
```




# # # GC3 calculation

#function to calculate GC-3
library(seqinr)


GC1<-function (s) {
  GC_for_1s<-GC(s[seq(1, length(s), 3)])
  return(GC_for_1s)
}

GC2<-function (s) {
  GC_for_2s<-GC(s[seq(2, length(s), 3)])
  return(GC_for_2s)
}

GC3<-function (s) {
  GC_for_3s<-GC(s[seq(3, length(s),3)])
  return(GC_for_3s)
}

GC1_2_3<-function (s) {
  GC_for_1s<-GC(s[seq(1, length(s), 3)])
  GC_for_2s<-GC(s[seq(2, length(s), 3)])
  GC_for_3s<-GC(s[seq(3, length(s),3)])
  return(list(GC_for_1s, GC_for_2s, GC_for_3s))
}

#GC by codons spread calculation



GC_spread <- function (seq) {
  library(zoo)
  library(seqinr)
  return(range(GC1(seq), GC2(seq), GC3(seq))[1]-range(GC1(seq), GC2(seq), GC3(seq))[2])
}


rollGCspread<-function (seq, interval) {
  library(zoo)
  library(seqinr)
  
  #GC_spread <- function (seq) {
   # return(range(GC1(seq), GC2(seq), GC3(seq))[1]-range(GC1(seq), GC2(seq), GC3(seq))[2])
  #}
  
  return(rollapply(seq, interval, GC_spread))
}



#sapply(dataset_gen, function(x) {return(GC1_2_3(strsplit(x$seq, '')))})
gen_allseqs<-sapply(dataset_gen, function(x) {return((strsplit(x$seq, '')))})
#gc1_gen_allseqs<-sapply(gen_allseqs, GC1)
#gc2_gen_allseqs<-sapply(gen_allseqs, GC2)
#gc3_gen_allseqs<-sapply(gen_allseqs, GC3)

spread_gc1_2_3_gen<-c()
for (i in seq_along(gen_allseqs)){
 spread_gc1_2_3_gen<-c(spread_gc1_2_3_gen, (range(GC1_2_3(gen_allseqs[[i]]))[2] - range(GC1_2_3(gen_allseqs[[i]]))[1]))
}
  

#for seqs - ALL promoters


#sapply(dataset_pro, function(x) {return(GC1_2_3(strsplit(x$seq, '')))})
pro_allseqs<-sapply(dataset_pro, function(x) {return((strsplit(x$seq, '')))})

spread_gc1_2_3_pro<-c()
for (i in seq_along(pro_allseqs)){
  spread_gc1_2_3_pro<-c(spread_gc1_2_3_pro, (range(GC1_2_3(pro_allseqs[[i]]))[2] - range(GC1_2_3(pro_allseqs[[i]]))[1]))
}

#for isl_seqs - islands


#sapply(dataset_pro, function(x) {return(GC1_2_3(strsplit(x$seq, '')))})
isl_allseqs<-sapply(dataset_isl, function(x) {return((strsplit(x$seq, '')))})

spread_gc1_2_3_isl<-c()
for (i in seq_along(isl_allseqs)){
  spread_gc1_2_3_isl<-c(spread_gc1_2_3_isl, (range(GC1_2_3(isl_allseqs[[i]]))[2] - range(GC1_2_3(isl_allseqs[[i]]))[1]))
}

#for lowscore_seqs - lowscore


#sapply(dataset_pro, function(x) {return(GC1_2_3(strsplit(x$seq, '')))})
lowscore_allseqs<-sapply(dataset_lowscore, function(x) {return((strsplit(x$seq, '')))})

spread_gc1_2_3_lowscore<-c()
for (i in seq_along(lowscore_allseqs)){
  spread_gc1_2_3_lowscore<-c(spread_gc1_2_3_lowscore, (range(GC1_2_3(lowscore_allseqs[[i]]))[2] - range(GC1_2_3(lowscore_allseqs[[i]]))[1]))
}

#for isl_seqs - islands


#sapply(dataset_pro, function(x) {return(GC1_2_3(strsplit(x$seq, '')))})
isl_allseqs<-sapply(dataset_isl, function(x) {return((strsplit(x$seq, '')))})

spread_gc1_2_3_isl<-c()
for (i in seq_along(isl_allseqs)){
  spread_gc1_2_3_isl<-c(spread_gc1_2_3_isl, (range(GC1_2_3(isl_allseqs[[i]]))[2] - range(GC1_2_3(isl_allseqs[[i]]))[1]))
}


#for notpro


#sapply(dataset_pro, function(x) {return(GC1_2_3(strsplit(x$seq, '')))})
not_allseqs<-sapply(dataset_notpro , function(x) {return((strsplit(x$seq, '')))})

spread_gc1_2_3_not<-c()
for (i in seq_along(not_allseqs)){
  spread_gc1_2_3_not<-c(spread_gc1_2_3_not, (range(GC1_2_3(not_allseqs[[i]]))[2] - range(GC1_2_3(not_allseqs[[i]]))[1]))
}


#for complete genome
GC3_total_genome<-range(GC1_2_3(e.coli_U00096.2_char))[2] - range(GC1_2_3(e.coli_U00096.2_char))[1]

lmts<-range(spread_gc1_2_3_pro, spread_gc1_2_3_isl, spread_gc1_2_3_gen, spread_gc1_2_3_lowscore, spread_gc1_2_3_not )
#lmts<-c(0,0.45)

par(mfrow=c(1,5))
boxplot(spread_gc1_2_3_pro, ylim = lmts)
boxplot(spread_gc1_2_3_isl, ylim = lmts)
boxplot(spread_gc1_2_3_gen, ylim = lmts)
boxplot(spread_gc1_2_3_lowscore, ylim = lmts)
boxplot(spread_gc1_2_3_not, ylim = lmts)

mtext('GC spread by triplets for promoters, islands, genes, and lowscore sequences', at=-1.5, line=1, cex=1.5 )

