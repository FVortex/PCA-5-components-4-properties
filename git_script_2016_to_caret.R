#setting working directory where the script, 2 octave scripts and .Rdata files are stored
setwd('/home/mikhail/Documents/Script_2016_all/PCA/')

# removing saved working space
rm(list=ls())

#loading R libraries
library(R.matlab)
library(seqinr)

## calculation of dynamical properties (activation energy and size of open states) from within R using 2 octave script
#writing in R E.coli K12 MG1655 genome from fasta file (must be copied to the working directory)

e.coli_U00096.2<-unlist(read.fasta('e.coli_U00096.2.fasta', seqonly = T))
#saving the genome version to .txt for further usage by octave scripts
writeLines(e.coli_U00096.2, 'e.coli_U00096.2.txt')
# check-up (whether the genome has right length and only 4 nucleotides composition)
check<-readLines('e.coli_U00096.2.txt')

##running octave scripts; the second one is for reading into octave genome in .txt format
system('octave dynamic_characteristics_vec_cumsum_read_txt.m')

system('octave dynamic_characteristics_vec_cumsum_read_txt_additional.m')

#reading in files created by octave scripts for activation energy and size of open states, and Gc-content (200 bp-long sliding window)

E01<-read.table('E01.mat')$V1
E02<-read.table('E02.mat')$V1
d1<-read.table('d1.mat')$V1
d2<-read.table('d2.mat')$V1
gc200matlab<-read.table('gc.mat')$V1

#loading data on sequences of different types (promoters, non-promoters, genes, islands, and lowscore) from .Rdata files (must be copied separetely)
load('/home/mikhail/Documents/Script_2016_all/spline_dataset_pro.Rdata')
load('/home/mikhail/Documents/Script_2016_all/spline_dataset_notpro.Rdata')
load('/home/mikhail/Documents/Script_2016_all/spline_dataset_gen.Rdata')
load('spline_dataset_isl.Rdata')
load('dataset_lowscore.Rdata')

# extracting data on all promoters and on experimentaly found ones - including previosely calculated electrostatic potential profiles

pro_names<-names(dataset_pro)
tsss<-c()
seqs<-c()
exp_strands<-c()
exp_tsss<-c()
mpots<-c() #for experimentally found one
exp_names<-c()

for (i in 1:length(dataset_pro)){
  tsss<-c(tsss, dataset_pro[[i]]$tss)
  seqs<-c(seqs, dataset_pro[[i]]$seq)
  if (dataset_pro[[i]]$evidence=='experimental'){ 
    #exp_proms<-rbind(exp_proms, c(pro_names[i], strands[i], tsss[i], seqs[i]))
    exp_names<-c(exp_names, pro_names[i])
    mpots<-rbind(mpots, dataset_pro[[i]]$mpot)  	
    #exp_promoters<-rbind(exp_promoters, c(pro_names[i], strands[i], tsss[i], seqs[i]))
    exp_tsss<-c(exp_tsss, dataset_pro[[i]]$tss) 
    exp_strands<-c(exp_strands, dataset_pro[[i]]$strand)}
}


#creating matrices for data on activation energy ('aeos1') and size ('aeos3'). The matrices is 699*201 - 699 promoters sequences, 201 value for a physical property profile
for (i in c('aeos1forward', 'aeos1reverse', 'aeos3forward', 'aeos3reverse', 'gc200forward', 'gc200reverse')) {
  assign(i, c())
}

for (i in 1:length(exp_tsss)) {
  if (exp_strands[i]=='forward') {
    aeos1forward<-rbind(aeos1forward, E01[as.numeric(exp_tsss[i]-150):(exp_tsss[i]+50)])
    #aeos2forward<-rbind(aeos2forward, matr1[(as.numericexp_tsss[i]-150):(as.numericexp_tsss[i]+50),2])
    aeos3forward<-rbind(aeos3forward, d1[(as.numeric(exp_tsss[i])-150):(exp_tsss[i]+50)])
    #aeos4forward<-rbind(aeos4forward, matr1[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50),4])            
    gc200forward<-rbind(gc200forward, gc200matlab[(as.numeric(exp_tsss[i])-150):(as.numeric(exp_tsss[i])+50)])
     } else {
    aeos1reverse<-rbind(aeos1reverse, E02[(as.numeric(exp_tsss[i])-150):(as.numeric(exp_tsss[i])+50)])
    #aeos2reverse<-rbind(aeos2reverse, matr1[(as.numeric(exp_tsss[i])-150):(as.numeric(exp_tsss[i])+50),2])
    aeos3reverse<-rbind(aeos3reverse, d2[(as.numeric(exp_tsss[i])-150):(as.numeric(exp_tsss[i])+50)])
    #aeos4reverse<-rbind(aeos4reverse, matr1[(as.numeric(exp_tsss[i])-150):(as.numeric(exp_tsss[i])+50),4])
    gc200reverse<-rbind(gc200reverse, rev(gc200matlab[(as.numeric(exp_tsss[i])-150):(as.numeric(exp_tsss[i])+50)]))
    }
}
#merging matrices for forward and reverse strands  together
aeos1<-rbind(aeos1forward, aeos1reverse)
#aeos2<-rbind(aeos2forward, aeos2reverse)
aeos3<-rbind(aeos3forward, aeos3reverse)
#aeos4<-rbind(aeos4forward, aeos4reverse)
gc200<-rbind(gc200forward, gc200reverse)



#extracting data for non-promoters

notpro_names<-names(dataset_notpro)
nottsss<-c()
notseqs<-c
notstrands<-c()
notmpots<-c()

for (i in 1:length(dataset_notpro)){
  nottsss<-c(nottsss, dataset_notpro[[i]]$tss)
  notseqs<-toupper(c(notseqs, dataset_notpro[[i]]$seq))
  notstrands<-c(notstrands, dataset_notpro[[i]]$strand)
  notmpots<-rbind(notmpots,  dataset_notpro[[i]]$mpot)  
}


#matrices for physical properties profiles
#creating matrices for data on activation energy ('aeos1') and size ('aeos3'). 
for (i in c('notaeos1forward', 'notaeos1reverse', 'notaeos3forward', 'notaeos3reverse', 'notgc200forward', 'notgc200reverse')) {
  assign(i, c())
}
for (i in 1:length(nottsss)) {
  if (notstrands[i]=='forward') {
    notaeos1forward<-rbind(notaeos1forward, E01[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50)])
    #notaeos2forward<-rbind(notaeos2forward, matr1[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50),2])
    notaeos3forward<-rbind(notaeos3forward, d1[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50)])
    #notaeos4forward<-rbind(notaeos4forward, matr1[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50),4])			
    notgc200forward<-rbind(notgc200forward, gc200matlab[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50)])
          } else {
    notaeos1reverse<-rbind(notaeos1reverse, E02[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50)])
    #notaeos2reverse<-rbind(notaeos2reverse, matr1[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50),2])
    notaeos3reverse<-rbind(notaeos3reverse, d2[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50)])
    #notaeos4reverse<-rbind(notaeos4reverse, matr1[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50),4])
    notgc200reverse<-rbind(notgc200reverse, rev(gc200matlab[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50)]))
    }
}

#mergind data for strands is not neaded - reverse strand is empty
notaeos1<-notaeos1forward
#notaeos2<-notaeos2forward
notaeos3<-notaeos3forward
#notaeos4<-notaeos4forward
notgc200<-notgc200forward

# # #genes. data axtracion
gen_names<-names(dataset_gen)
gentsss<-c()
genseqs<-c
genstrands<-c()
genmpots<-c()

for (i in 1:length(dataset_gen)){
  gentsss<-c(gentsss, dataset_gen[[i]]$tss)
  genseqs<-toupper(c(genseqs, dataset_gen[[i]]$seq))
  genstrands<-c(genstrands, dataset_gen[[i]]$strand)
  genmpots<-rbind(genmpots, dataset_gen[[i]]$mpot)
}

for (i in c('genaeos1forward', 'genaeos1reverse', 'genaeos3forward', 'genaeos3reverse', 'gengc200forward', 'gengc200reverse' )) {
  assign(i, c())
}

for (i in 1:length(dataset_gen)) {
  if (genstrands[i]=='forward') {
    genaeos1forward<-rbind(genaeos1forward, E01[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50)])
    #   genaeos2forward<-rbind(genaeos2forward, matr1[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50),2])
    genaeos3forward<-rbind(genaeos3forward, d1[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50)])
    #  genaeos4forward<-rbind(genaeos4forward, matr1[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50),4])            
    gengc200forward<-rbind(gengc200forward, gc200matlab[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50)])
        } else {
    genaeos1reverse<-rbind(genaeos1reverse, E02[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50)])
    # genaeos2reverse<-rbind(genaeos2reverse, matr1[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50),2])
    genaeos3reverse<-rbind(genaeos3reverse, d2[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50)])
    #genaeos4reverse<-rbind(genaeos4reverse, matr1[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50),4])
    gengc200reverse<-rbind(gengc200reverse, rev(gc200matlab[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50)]))
    }
}

genaeos1<-rbind(genaeos1forward, genaeos1reverse)
#genaeos2<-rbind(genaeos2forward, genaeos2reverse)
genaeos3<-rbind(genaeos3forward, genaeos3reverse)
#genaeos4<-rbind(genaeos4forward, genaeos4reverse)
gengc200<-rbind(gengc200forward, gengc200reverse)

## islands. data axtracion
isl_names<-names(dataset_isl)
isltsss<-c()
islseqs<-c
islstrands<-c()
islmpots<-c()

for (i in 1:length(dataset_isl)){
  isltsss<-c(isltsss, dataset_isl[[i]]$tss)
  islseqs<-toupper(c(islseqs, dataset_isl[[i]]$seq))
  islstrands<-c(islstrands, dataset_isl[[i]]$strand)
  islmpots<-rbind(islmpots, dataset_isl[[i]]$mpot)
}

for (i in c('islaeos1forward', 'islaeos1reverse', 'islaeos3forward', 'islaeos3reverse', 'islgc200forward', 'islgc200reverse' )) {
  assign(i, c())
}

for (i in 1:length(dataset_isl)) {
  if (islstrands[i]=='forward') {
    islaeos1forward<-rbind(islaeos1forward, E01[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50)])
    #   islaeos2forward<-rbind(islaeos2forward, matr1[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50),2])
    islaeos3forward<-rbind(islaeos3forward, d1[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50)])
    #  islaeos4forward<-rbind(islaeos4forward, matr1[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50),4])            
    islgc200forward<-rbind(islgc200forward, gc200matlab[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50)])
  } else {
    islaeos1reverse<-rbind(islaeos1reverse, E02[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50)])
    # islaeos2reverse<-rbind(islaeos2reverse, matr1[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50),2])
    islaeos3reverse<-rbind(islaeos3reverse, d2[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50)])
    #islaeos4reverse<-rbind(islaeos4reverse, matr1[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50),4])
    islgc200reverse<-rbind(islgc200reverse, rev(gc200matlab[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50)]))
  }
}

islaeos1<-rbind(islaeos1forward, islaeos1reverse)
#islaeos2<-rbind(islaeos2forward, islaeos2reverse)
islaeos3<-rbind(islaeos3forward, islaeos3reverse)
#islaeos4<-rbind(islaeos4forward, islaeos4reverse)
islgc200<-rbind(islgc200forward, islgc200reverse)

## lowscore. data axtracion
lowscore_names<-names(dataset_lowscore)
lowscoretsss<-c()
lowscoreseqs<-c
lowscorestrands<-c()
lowscorempots<-c()

#only first 2000 among lowscore sequences are used
for (i in 1:2000){
  lowscoretsss<-c(lowscoretsss, dataset_lowscore[[i]]$tss)
  lowscoreseqs<-toupper(c(lowscoreseqs, dataset_lowscore[[i]]$seq))
  lowscorestrands<-c(lowscorestrands, dataset_lowscore[[i]]$strand)
  lowscorempots<-rbind(lowscorempots, dataset_lowscore[[i]]$mpot)
}

for (i in c('lowscoreaeos1forward', 'lowscoreaeos1reverse', 'lowscoreaeos3forward', 'lowscoreaeos3reverse', 'lowscoregc200forward', 'lowscoregc200reverse' )) {
  assign(i, c())
}

#only first 2000 among lowscore sequences are used
for (i in 1:2000) {
  if (lowscorestrands[i]=='forward') {
    lowscoreaeos1forward<-rbind(lowscoreaeos1forward, E01[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50)])
    #   lowscoreaeos2forward<-rbind(lowscoreaeos2forward, matr1[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50),2])
    lowscoreaeos3forward<-rbind(lowscoreaeos3forward, d1[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50)])
    #  lowscoreaeos4forward<-rbind(lowscoreaeos4forward, matr1[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50),4])            
    lowscoregc200forward<-rbind(lowscoregc200forward, gc200matlab[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50)])
  } else {
    lowscoreaeos1reverse<-rbind(lowscoreaeos1reverse, E02[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50)])
    # lowscoreaeos2reverse<-rbind(lowscoreaeos2reverse, matr1[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50),2])
    lowscoreaeos3reverse<-rbind(lowscoreaeos3reverse, d2[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50)])
    #lowscoreaeos4reverse<-rbind(lowscoreaeos4reverse, matr1[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50),4])
    lowscoregc200reverse<-rbind(lowscoregc200reverse, rev(gc200matlab[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50)]))
  }
}

lowscoreaeos1<-rbind(lowscoreaeos1forward, lowscoreaeos1reverse)
#lowscoreaeos2<-rbind(lowscoreaeos2forward, lowscoreaeos2reverse)
lowscoreaeos3<-rbind(lowscoreaeos3forward, lowscoreaeos3reverse)
#lowscoreaeos4<-rbind(lowscoreaeos4forward, lowscoreaeos4reverse)
lowscoregc200<-rbind(lowscoregc200forward, lowscoregc200reverse)

#setting names

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


#Merging datasets. Normalizing data

to_pca_5components_4props<-rbind(cbind(scale(aeos1), scale(aeos3), scale(mpots), scale(gc200)),
                                 cbind(scale(notaeos1), scale(notaeos3), scale(notmpots), scale(notgc200)),
                                 cbind(scale(genaeos1), scale(genaeos3), scale(genmpots), scale(gengc200)),
                                 cbind(scale(islaeos1), scale(islaeos3), scale(islmpots), scale(islgc200)),
                                 cbind(scale(lowscoreaeos1), scale(lowscoreaeos3), scale(lowscorempots), scale(lowscoregc200))
)
# saving the initial matrix for 5 components, 4 properties
save(to_pca_5components_4props, file='to_pca_5components_4props.Rdata')


#setting sequences groups
habillage_5components_4props<-c(rep('P', length(exp_names)), rep('N', length(dataset_notpro)), rep('G', length(dataset_gen)), rep('I', length(dataset_isl)), rep('L', 2000))

colnames(to_pca_5components_4props)<-NULL


#############3 #scale =F?
princ.return.5comps.4props <- prcomp(to_pca_5components_4props, scale=T)

#saving prcomp output 
save(princ.return.5comps.4props, file='princ.return.5comps.4props.Rdata')

library("factoextra")
eig.val <- get_eigenvalue(princ.return.5comps.4props)

# Eigenvalues
eig <- (princ.return.5comps.4props$sdev)^2

# Variances in percentage
variance <- eig*100/sum(eig)

# Cumulative variances
cumvar <- cumsum(variance)
png('cumvar_5components_4props.png', height=1250, width=1250, res=130)
plot(cumvar, type='l', main='Cumulative variance for principal components \n on 4 properties for 5 components mixture', ylab='Cumulative variance (%)', xlab='Principal components')
v.5comps.4props=100 # 65 
h.5comps.4props=98
abline(h=h.5comps.4props, v=v.5comps.4props, col=12)
text(v.5comps.4props-15, h.5comps.4props ,paste('Number of \nPCs=', v.5comps.4props,'\n' ,h.5comps.4props, '% of \nvariance \nretained'),srt=0.2,pos=3)

dev.off()



pcs.from.var<-c()
for (i in 1:v.5comps.4props) {
  load <- princ.return.5comps.4props$rotation[,i]
  pr.cp <- to_pca_5components_4props %*% load
  pr <- as.numeric(pr.cp)
  pcs.from.var<-cbind(pcs.from.var, pr)
}

#saving PCA data: variable are conversed to PCs

rownames(pcs.from.var)<-rownames(to_pca_5components_4props)
colnames(pcs.from.var)<-1:ncol(pcs.from.var)
#renaming the variable
mixt_5comps_after_pca_on_4_props_ae_size_mpots_gc200<- pcs.from.var
save(mixt_5comps_after_pca_on_4_props_ae_size_mpots_gc200, file='new_no_factor_8XII_mixt_5comps_after_pca_on_4_props_ae_size_mpots_gc200.Rdata')
#geneal PCA results visualization

library(factoextra)

fviz_pca_ind(princ.return.5comps.4props, label="none", habillage=habillage_5components_4props, pointsize = 0.5, addEllipses = T, ellipse.level = 0.99)


library(rgl)

habillage.pro.not.isl.gen<-c(rep('red', length(exp_names)), rep ('blue', length(dataset_notpro)), rep('green', length(dataset_gen)), rep('orange', length(dataset_isl)), rep('magenta', 2000)) 
#4 components correspond to colors 1-4
plot3d(princ.return.5comps.4props$x[,1:3], col=(habillage.pro.not.isl.gen), size=1)
plot3d(princ.return.pro.not.isl.gen.aeos3$x[,1:3], col=(habillage.pro.not.isl.gen), size=1)
# #pca3d
# #
library(pca3d)
pca3d(princ.return.5comps.4props, group = habillage.pro.not.isl.gen, col = habillage.pro.not.isl.gen, radius = 0.5)



# # # SUPERVISED MACHINE LEARNING
library(data.table)
library(caret)
library(doMC)
#library(DMwR)

registerDoMC(cores = 3)
#load("~/work/2016/iteb/promoters_nonpromoters_after_pca_on_4_props_ae_size_mpots_gc200.Rdata")
#load("~/work/2016/iteb/promoters_nonpromoters_initial_ae_data.Rdata")
#load("new_8XII_mixt_5comps_after_pca_on_4_props_ae_size_mpots_gc200.Rdata")
#load("new_no_factor_8XII_mixt_5comps_after_pca_on_4_props_ae_size_mpots_gc200.Rdata")
df <- as.data.frame(mixt_5comps_after_pca_on_4_props_ae_size_mpots_gc200)
#df<-cbind(habillage_5components_4props, df)
promoters <- df[1:699,]
non_promoters <- df[700:2579,]
genes <- df[2580:6006,]
islands <- df[6007:8234,]
lowscore <- df[8235:10234,]


factor_to_promoters_vs_lowscore<-as.factor(c(rep('Promoter', nrow(promoters)), rep('Lowscore', nrow(lowscore))))
promoters_vs_lowscore <- cbind(factor_to_promoters_vs_lowscore, rbind(promoters, lowscore))

factor_to_promoters_vs_non_promoters<-as.factor(c(rep('Promoter', nrow(promoters)), rep('Non_promoter', nrow(non_promoters))))
promoters_vs_non_promoters <- cbind(factor_to_promoters_vs_non_promoters, rbind(promoters, non_promoters))

factor_to_promoters_vs_islands<-as.factor(c(rep('Promoter', nrow(promoters)), rep('Island', nrow(islands))))
promoters_vs_islands <- cbind(factor_to_promoters_vs_islands, rbind(promoters, islands))

factor_to_promoters_vs_genes<-as.factor(c(rep('Promoter', nrow(promoters)), rep('Gene', nrow(genes))))
promoters_vs_genes <- cbind(factor_to_promoters_vs_genes, rbind(promoters, genes))

#for (data_set in c(promoters_vs_lowscore, promoters_vs_non_promoters, promoters_vs_islands, promoters_vs_genes)) {


#for promoters_vs_lowscore
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

#train_control <- trainControl(method="repeatedcv", number=10, repeats=3)
fit_promoters_vs_lowscore <- train(factor_to_promoters_vs_lowscore ~ .,
                                   data = training,
                                   method = "nb",#nb, lda, nbDiscrete, nbSearch
                                   preProcess=c("center", "scale"),
                                   tuneLength = 15,
                                   trControl = fitControl,
                                   metric = "Kappa"
)

predictionClasses_promoters_vs_lowscore <- predict(fit_promoters_vs_lowscore, newdata = testing)
predictionProb_promoters_vs_lowscore <- predict(fit_promoters_vs_lowscore, newdata = testing, type ="prob")
confusionMatrix_promoters_vs_lowscore <- confusionMatrix(data = predictionClasses_promoters_vs_lowscore, testing$factor_to_promoters_vs_lowscore)



#for promoters_vs_non_promoters
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
                                        method = "rf",#nb, lda, nbDiscrete, nbSearch
                                        preProcess=c("center", "scale"),
                                        tuneLength = 7,
                                        trControl = fitControl,
                                        metric = "ROC"
)
predictionClasses_promoters_vs_non_promoters <- predict(fit_promoters_vs_non_promoters, newdata = testing)
predictionProb_promoters_vs_non_promoters <- predict(fit_promoters_vs_non_promoters, newdata = testing, type ="prob")
confusionMatrix_promoters_vs_non_promoters <- confusionMatrix(data = predictionClasses_promoters_vs_non_promoters, testing$factor_to_promoters_vs_non_promoters)

#for promoters_vs_islands
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
                                  method = "rf",#nb, lda, nbDiscrete, nbSearch
                                  preProcess=c("center", "scale"),
                                  tuneLength = 7,
                                  trControl = fitControl,
                                  metric = "ROC"
)
predictionClasses_promoters_vs_islands <- predict(fit_promoters_vs_islands, newdata = testing)
predictionProb_promoters_vs_islands <- predict(fit_promoters_vs_islands, newdata = testing, type ="prob")
confusionMatrix_promoters_vs_islands <- confusionMatrix(data = predictionClasses_promoters_vs_islands, testing$factor_to_promoters_vs_islands)

#for promoters_vs_genes
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
                                method = "rf",#nb, lda, nbDiscrete, nbSearch
                                preProcess=c("center", "scale"),
                                tuneLength = 7,
                                trControl = fitControl,
                                metric = "ROC"
)
predictionClasses_promoters_vs_genes <- predict(fit_promoters_vs_genes, newdata = testing)
predictionProb_promoters_vs_genes <- predict(fit_promoters_vs_genes, newdata = testing, type ="prob")
confusionMatrix_promoters_vs_genes <- confusionMatrix(data = predictionClasses_promoters_vs_genes, testing$factor_to_promoters_vs_genes)


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


sapply(dataset_gen, function(x) {return(GC1_2_3(strsplit(x$seq, '')))})
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

#bendability profiles
#library(reldna)


#for (i in c('bendforward', 'bendreverse')) {
#  assign(i, c())
#}

#for (i in 1:length(exp_names)) {
#  if (exp_strands[i]=='forward') {
#    aeos1forward<-rbind(aeos1forward, E01[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50)])
#    #aeos2forward<-rbind(aeos2forward, matr1[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50),2])
#    aeos3forward<-rbind(aeos3forward, d1[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50)])
    #aeos4forward<-rbind(aeos4forward, matr1[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50),4])            
#  } else {
#    aeos1reverse<-rbind(aeos1reverse, E02[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50)])
    #aeos2reverse<-rbind(aeos2reverse, matr1[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50),2])
#    aeos3reverse<-rbind(aeos3reverse, d2[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50)])
    #aeos4reverse<-rbind(aeos4reverse, matr1[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50),4])
    
#  }
  
#}

#aeos1forward<-aeos1forward[-1,] #удаление первой строки с нулями
#aeos1reverse<-aeos1reverse[-1,]
#dimnames(aeos1reverse)<-NULL #удаление аттрибутов
#dimnames(aeos1forward)<-NULL

#aeos2forward<-aeos2forward[-1,] #СѓРґР°Р»РµРЅРёРµ РёСЃС…РѕРґРЅС‹С… РЅСѓР»РµР№
#aeos2reverse<-aeos2reverse[-1,]
#dimnames(aeos2reverse)<-NULL
#dimnames(aeos2forward)<-NULL

#aeos3forward<-aeos3forward[-1,] #СѓРґР°Р»РµРЅРёРµ РёСЃС…РѕРґРЅС‹С… РЅСѓР»РµР№
#aeos3reverse<-aeos3reverse[-1,]
#dimnames(aeos3reverse)<-NULL
#dimnames(aeos3forward)<-NULL

#aeos4forward<-aeos4forward[-1,] #СѓРґР°Р»РµРЅРёРµ РёСЃС…РѕРґРЅС‹С… РЅСѓР»РµР№
#aeos4reverse<-aeos4reverse[-1,]
#dimnames(aeos4reverse)<-NULL
#dimnames(aeos4forward)<-NULL

#РѕР±СЉРµРґРёРЅРµРЅРёРµ РґР°РЅРЅС‹С… РґР»СЏ РѕР±РѕРёС… СЃС‚СЂСЌРЅРґРѕРІ

#aeos1<-rbind(aeos1forward, aeos1reverse)
#aeos2<-rbind(aeos2forward, aeos2reverse)
#aeos3<-rbind(aeos3forward, aeos3reverse)
#aeos4<-rbind(aeos4forward, aeos4reverse)



