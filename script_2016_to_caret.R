setwd('/home/mikhail/Documents/Script_2016_all/PCA/')

library(R.matlab)
library(seqinr)
# #mat_all<-readMat('ecoli_full_data_set_dynamic_chars.mat')

#mat1<-readMat('M1_ecoli.mat')

#mat2<-readMat('M2_ecoli.mat')

#E01<-as.matrix(mat1$M1)

#E02<-as.matrix(mat2$M2)

## alternatively - dynamical properties calculation from within R
e.coli_U00096.2<-unlist(read.fasta('e.coli_U00096.2.fasta', seqonly = T))
writeLines(e.coli_U00096.2, 'e.coli_U00096.2.txt')
# check-up 
check<-readLines('e.coli_U00096.2.txt')

system('octave dynamic_characteristics_vec_cumsum_read_txt.m')

system('octave dynamic_characteristics_vec_cumsum_read_txt_additional.m')

#reading in files created by octave

E01<-read.table('E01.mat')$V1
E02<-read.table('E02.mat')$V1
d1<-read.table('d1.mat')$V1
d2<-read.table('d2.mat')$V1
gc200matlab<-read.table('gc.mat')$V1
#system('octave [E01, E02, d1, d2, m1, m2, c1, c2, taus1, taus2, gc, gc_comp] = dynamic_characteristics_vec_cumsum_read_txt("e.coli_U00096.2.fasta", 200, 0)')
# #gc200matlab<-as.numeric(readMat('gc_ecoli_200.mat')$gc)

# repeated 200 b.p.-long sliding window GC-calculation in R
#library(zoo)
#library(seqinr)
# #totgenome<-read.fasta("/home/mikhail/Documents/Script_2016_all/Genome_to_dataset_pro_U00096.2.txt",seqtype = "DNA", as.string = TRUE, forceDNAtolower = TRUE, set.attributes = FALSE, seqonly=TRUE) 
#e.coli_U00096.2_char<-unlist(strsplit(e.coli_U00096.2, ''))

#totgenome_rollapply_GC<-rollapply(e.coli_U00096.2_char, 200, GC)
#susbtitution of calculated in octave GC-content
#gc200matlab<-c(rep(mean(totgenome_rollapply_GC), 100), totgenome_rollapply_GC, rep(mean(totgenome_rollapply_GC), 100))

# #alternatively from the new file
# #E01<-as.numeric(mat_all$E01)
# #E02<-as.numeric(mat_all$E02)

# #d1<-as.numeric(mat_all$d1)
# #d2<-as.numeric(mat_all$d2)

# #gc200matlab<-as.numeric(mat_all$gc)

#loading data on sequences
load('/home/mikhail/Documents/Script_2016_all/spline_dataset_pro.Rdata')
load('/home/mikhail/Documents/Script_2016_all/spline_dataset_notpro.Rdata')
load('/home/mikhail/Documents/Script_2016_all/spline_dataset_gen.Rdata')
load('/home/mikhail/Documents/Script_2016_all/spline_dataset_isl.Rdata')
#load('spline_dataset_rand.Rdata')
load('dataset_lowscore.Rdata')
#функция для реверсии строк (strings)

strReverse <- function(x)
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")

#извлечение данных о промоторах

pro_names<-names(dataset_pro)
tsss<-c()
seqs<-c()
exp_strands<-c()
exp_tsss<-c()
exp_promoters<-c()
exp_proms<-c()
mpots<-c()

for (i in 1:length(dataset_pro)){
  tsss<-c(tsss, dataset_pro[[i]]$tss)
  seqs<-c(seqs, dataset_pro[[i]]$seq)
#  strands<-c(strands, dataset_pro[[i]]$strand)
  if (dataset_pro[[i]]$evidence=='experimental'){ 
    exp_proms<-rbind(exp_proms, c(pro_names[i], strands[i], tsss[i], seqs[i]))
    mpots<-rbind(mpots, dataset_pro[[i]]$mpot)  	
    exp_promoters<-rbind(exp_promoters, c(pro_names[i], strands[i], tsss[i], seqs[i]))
    exp_tsss<-c(exp_tsss, dataset_pro[[i]]$tss) 
    exp_strands<-c(exp_strands, dataset_pro[[i]]$strand)}
}



for (i in c('aeos1forward', 'aeos1reverse', 'aeos3forward', 'aeos3reverse')) {
  assign(i, c())
}
#aeos1forward<-c()
#aeos1reverse<-c()
#aeos2forward<-rep(0,201)
#aeos2reverse<-rep(0,201)
#aeos3forward<-c()
#aeos3reverse<-c()
#aeos4forward<-rep(0,201)
#aeos4reverse<-rep(0,201)
for (i in 1:length(exp_proms[,1])) {
  if (exp_strands[i]=='forward') {
    aeos1forward<-rbind(aeos1forward, E01[(as.numeric(exp_tsss[i])-150):(as.numeric(exp_proms[i,3])+50)])
    #aeos2forward<-rbind(aeos2forward, matr1[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50),2])
    aeos3forward<-rbind(aeos3forward, d1[(as.numeric(exp_tsss[i])-150):(as.numeric(exp_proms[i,3])+50)])
    #aeos4forward<-rbind(aeos4forward, matr1[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50),4])            
  } else {
    aeos1reverse<-rbind(aeos1reverse, E02[(as.numeric(exp_tsss[i])-150):(as.numeric(exp_proms[i,3])+50)])
    #aeos2reverse<-rbind(aeos2reverse, matr1[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50),2])
    aeos3reverse<-rbind(aeos3reverse, d2[(as.numeric(exp_tsss[i])-150):(as.numeric(exp_proms[i,3])+50)])
    #aeos4reverse<-rbind(aeos4reverse, matr1[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50),4])
    
  }
  
}

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

aeos1<-rbind(aeos1forward, aeos1reverse)
#aeos2<-rbind(aeos2forward, aeos2reverse)
aeos3<-rbind(aeos3forward, aeos3reverse)
#aeos4<-rbind(aeos4forward, aeos4reverse)



gc200forward<-c()
gc200reverse<-c()

#forward strand only


for (i in 1:length(exp_tsss)) {
  if (exp_strands[i]=='forward') {
    gc200forward<-rbind(gc200forward, gc200matlab[(as.numeric(exp_tsss[i])-150):(as.numeric(exp_tsss[i])+50)])
    #y<-gc200matlab[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50)]
    #x<-1:length(y)
    #lo<- loess(y~x)
    #lines(predict(lo), col='red', lwd=0.2)    
  } else {
    gc200reverse<-rbind(gc200reverse, rev(gc200matlab[(as.numeric(exp_tsss[i])-150):(as.numeric(exp_tsss[i])+50)]))
    #y<-rev(gc200matlab[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50)])
    #x<-1:length(y)
    #lo<-loess(y~x)
    #lines(predict(lo), col='red', lwd=0.2)    
  }
}


#gc200forward<-gc200forward[-1,] #удаление первой строки с нулями
#gc200reverse<-gc200reverse[-1,]
#dimnames(gc200reverse)<-NULL #удаление аттрибутов
#dimnames(gc200forward)<-NULL

gc200<-rbind(gc200forward, gc200reverse)



#извлечение данных о НЕпромоторах

load("spline_dataset_notpro.Rdata")

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


#извлечение интервалов физических характеристик [-150, 50] п.о. от ТСТ для всех 4 параметров динамической модели (случай непромоторных областей)

for (i in c('notaeos1forward', 'notaeos1reverse', 'notaeos3forward', 'notaeos3reverse')) {
  assign(i, c())
}
#notaeos1forward<-rep(0,201)
#notaeos1reverse<-rep(0,201)
#notaeos2forward<-rep(0,201)
#notaeos2reverse<-rep(0,201)
#notaeos3forward<-rep(0,201)
#notaeos3reverse<-rep(0,201)
#notaeos4forward<-rep(0,201)
#notaeos4reverse<-rep(0,201)
for (i in 1:length(nottsss)) {
  if (notstrands[i]=='forward') {
    notaeos1forward<-rbind(notaeos1forward, E01[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50)])
    #notaeos2forward<-rbind(notaeos2forward, matr1[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50),2])
    notaeos3forward<-rbind(notaeos3forward, d1[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50)])
    #notaeos4forward<-rbind(notaeos4forward, matr1[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50),4])			
  } else {
    notaeos1reverse<-rbind(notaeos1reverse, E02[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50)])
    #notaeos2reverse<-rbind(notaeos2reverse, matr1[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50),2])
    notaeos3reverse<-rbind(notaeos3reverse, d2[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50)])
    #notaeos4reverse<-rbind(notaeos4reverse, matr1[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50),4])
  }
}

#notaeos1forward<-notaeos1forward[-1,] #удаление исходных нулей
#notaeos1reverse<-notaeos1reverse[-1,]
#dimnames(notaeos1reverse)<-NULL
#dimnames(notaeos1forward)<-NULL

#notaeos2forward<-notaeos2forward[-1,] #удаление исходных нулей
#notaeos2reverse<-notaeos2reverse[-1,]
#dimnames(notaeos2reverse)<-NULL
#dimnames(notaeos2forward)<-NULL

#notaeos3forward<-notaeos3forward[-1,] #удаление исходных нулей
#notaeos3reverse<-notaeos3reverse[-1,]
#dimnames(notaeos3reverse)<-NULL
#dimnames(notaeos3forward)<-NULL

#notaeos4forward<-notaeos4forward[-1,] #удаление исходных нулей
#notaeos4reverse<-notaeos4reverse[-1,]
#dimnames(notaeos4reverse)<-NULL
#dimnames(notaeos4forward)<-NULL

#объединение данных для обоих стрэндов не требуется- непромоторные области на 2м стрэнде не взяты


notaeos1<-notaeos1forward
#notaeos2<-notaeos2forward
notaeos3<-notaeos3forward
#notaeos4<-notaeos4forward


#forward strand only


notgc200forward<-c()
notgc200reverse<-c()

for (i in 1:length(dataset_notpro)) {
  if (notstrands[i]=='forward') {
    notgc200forward<-rbind(notgc200forward, gc200matlab[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50)])
    #y<-gc200matlab[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50)]
    #x<-1:length(y)
    #lo<- loess(y~x)
    #lines(predict(lo), col='red', lwd=0.2)    
  } else {
    notgc200reverse<-rbind(notgc200reverse, rev(gc200matlab[(as.numeric(nottsss[i])-150):(as.numeric(nottsss[i])+50)]))
    #y<-rev(gc200matlab[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50)])
    #x<-1:length(y)
    #lo<-loess(y~x)
    #lines(predict(lo), col='red', lwd=0.2)    
  }
}


#notgc200forward<-notgc200forward[-1,] #удаление первой строки с нулями
#notgc200reverse<-notgc200reverse[-1,] no reverse strand sequences
#dimnames(notgc200reverse)<-NULL #удаление аттрибутов
#dimnames(notgc200forward)<-NULL

notgc200<-notgc200forward

#dimnames(notgc200)<-NULL

# # #genes
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

for (i in c('genaeos1forward', 'genaeos1reverse', 'genaeos3forward', 'genaeos3reverse')) {
  assign(i, c())
}
#genaeos1forward<-rep(0,201)
#genaeos1reverse<-rep(0,201)
#genaeos2forward<-rep(0,201)
#genaeos2reverse<-rep(0,201)
#genaeos3forward<-rep(0,201)
#genaeos3reverse<-rep(0,201)
#genaeos4forward<-rep(0,201)
#genaeos4reverse<-rep(0,201)
for (i in 1:length(dataset_gen)) {
  if (genstrands[i]=='forward') {
    genaeos1forward<-rbind(genaeos1forward, E01[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50)])
    #   genaeos2forward<-rbind(genaeos2forward, matr1[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50),2])
    genaeos3forward<-rbind(genaeos3forward, d1[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50)])
    #  genaeos4forward<-rbind(genaeos4forward, matr1[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50),4])            
  } else {
    genaeos1reverse<-rbind(genaeos1reverse, E02[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50)])
    # genaeos2reverse<-rbind(genaeos2reverse, matr1[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50),2])
    genaeos3reverse<-rbind(genaeos3reverse, d2[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50)])
    #genaeos4reverse<-rbind(genaeos4reverse, matr1[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50),4])
  }
}

#genaeos1forward<-genaeos1forward[-1,] #удаление первой строки с нулями
#genaeos1reverse<-genaeos1reverse[-1,]
#dimnames(genaeos1reverse)<-NULL #удаление аттрибутов
#dimnames(genaeos1forward)<-NULL

#genaeos2forward<-genaeos2forward[-1,] #СѓРґР°Р»РµРЅРёРµ РёСЃС…РѕРґРЅС‹С… РЅСѓР»РµР№
#genaeos2reverse<-genaeos2reverse[-1,]
#dimnames(genaeos2reverse)<-NULL
#dimnames(genaeos2forward)<-NULL

#genaeos3forward<-genaeos3forward[-1,] #СѓРґР°Р»РµРЅРёРµ РёСЃС…РѕРґРЅС‹С… РЅСѓР»РµР№
#genaeos3reverse<-genaeos3reverse[-1,]
#dimnames(genaeos3reverse)<-NULL
#dimnames(genaeos3forward)<-NULL

#genaeos4forward<-genaeos4forward[-1,] #СѓРґР°Р»РµРЅРёРµ РёСЃС…РѕРґРЅС‹С… РЅСѓР»РµР№
#genaeos4reverse<-genaeos4reverse[-1,]
#dimnames(genaeos4reverse)<-NULL
#dimnames(genaeos4forward)<-NULL

#РѕР±СЉРµРґРёРЅРµРЅРёРµ РґР°РЅРЅС‹С… РґР»СЏ РѕР±РѕРёС… СЃС‚СЂСЌРЅРґРѕРІ

genaeos1<-rbind(genaeos1forward, genaeos1reverse)
#genaeos2<-rbind(genaeos2forward, genaeos2reverse)
genaeos3<-rbind(genaeos3forward, genaeos3reverse)
#genaeos4<-rbind(genaeos4forward, genaeos4reverse)



gengc200forward<-c()
gengc200reverse<-c()

#forward strand only


for (i in 1:length(gentsss)) {
  if (genstrands[i]=='forward') {
    gengc200forward<-rbind(gengc200forward, gc200matlab[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50)])
    #y<-gc200matlab[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50)]
    #x<-1:length(y)
    #lo<- loess(y~x)
    #lines(predict(lo), col='red', lwd=0.2)    
  } else {
    gengc200reverse<-rbind(gengc200reverse, rev(gc200matlab[(as.numeric(gentsss[i])-150):(as.numeric(gentsss[i])+50)]))
    #y<-rev(gc200matlab[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50)])
    #x<-1:length(y)
    #lo<-loess(y~x)
    #lines(predict(lo), col='red', lwd=0.2)    
  }
  
}


#gengc200forward<-gengc200forward[-1,] #удаление первой строки с нулями
#gengc200reverse<-gengc200reverse[-1,]
#dimnames(gengc200reverse)<-NULL #удаление аттрибутов
#dimnames(gengc200forward)<-NULL

gengc200<-rbind(gengc200forward, gengc200reverse)



## ISL

isl_names<-names(dataset_isl)
isltsss<-c()
islseqs<-c()
islstrands<-c()
islmpots<-c()

for (i in 1:length(dataset_isl)){
  isltsss<-c(isltsss, dataset_isl[[i]]$tss)
  islseqs<-toupper(c(islseqs, dataset_isl[[i]]$seq))
  islstrands<-c(islstrands, dataset_isl[[i]]$strand)
  islmpots<-rbind(islmpots, dataset_isl[[i]]$mpot)
}


for (i in c('islaeos1forward', 'islaeos1reverse', 'islaeos3forward', 'islaeos3reverse')) {
  assign(i, c())
}
#islaeos1forward<-rep(0,201)
#islaeos1reverse<-rep(0,201)
#islaeos2forward<-rep(0,201)
#islaeos2reverse<-rep(0,201)
#islaeos3forward<-rep(0,201)
#islaeos3reverse<-rep(0,201)
#islaeos4forward<-rep(0,201)
#islaeos4reverse<-rep(0,201)

for (i in 1:length(dataset_isl)) {
  if (islstrands[i]=='forward') {
    islaeos1forward<-rbind(islaeos1forward, E01[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50)])
    #   islaeos2forward<-rbind(islaeos2forward, matr1[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50),2])
    islaeos3forward<-rbind(islaeos3forward, d1[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50)])
    #  islaeos4forward<-rbind(islaeos4forward, matr1[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50),4])            
  } else {
    islaeos1reverse<-rbind(islaeos1reverse, E02[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50)])
    # islaeos2reverse<-rbind(islaeos2reverse, matr1[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50),2])
    islaeos3reverse<-rbind(islaeos3reverse, d2[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50)])
    #islaeos4reverse<-rbind(islaeos4reverse, matr1[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50),4])
  }
}

#islaeos1forward<-islaeos1forward[-1,] #удаление первой строки с нулями
#islaeos1reverse<-islaeos1reverse[-1,]
#dimnames(islaeos1reverse)<-NULL #удаление аттрибутов
#dimnames(islaeos1forward)<-NULL

#islaeos2forward<-islaeos2forward[-1,] #СѓРґР°Р»РµРЅРёРµ РёСЃС…РѕРґРЅС‹С… РЅСѓР»РµР№
#islaeos2reverse<-islaeos2reverse[-1,]
#dimnames(islaeos2reverse)<-NULL
#dimnames(islaeos2forward)<-NULL

#islaeos3forward<-islaeos3forward[-1,] #СѓРґР°Р»РµРЅРёРµ РёСЃС…РѕРґРЅС‹С… РЅСѓР»РµР№
#islaeos3reverse<-islaeos3reverse[-1,]
#dimnames(islaeos3reverse)<-NULL
#dimnames(islaeos3forward)<-NULL

#islaeos4forward<-islaeos4forward[-1,] #СѓРґР°Р»РµРЅРёРµ РёСЃС…РѕРґРЅС‹С… РЅСѓР»РµР№
#islaeos4reverse<-islaeos4reverse[-1,]
#dimnames(islaeos4reverse)<-NULL
#dimnames(islaeos4forward)<-NULL

#РѕР±СЉРµРґРёРЅРµРЅРёРµ РґР°РЅРЅС‹С… РґР»СЏ РѕР±РѕРёС… СЃС‚СЂСЌРЅРґРѕРІ

islaeos1<-rbind(islaeos1forward, islaeos1reverse)
#islaeos2<-rbind(islaeos2forward, islaeos2reverse)
islaeos3<-rbind(islaeos3forward, islaeos3reverse)
#islaeos4<-rbind(islaeos4forward, islaeos4reverse)
#dimnames(islaeos1)<-NULL
#dimnames(islaeos3)<-NULL

islgc200forward<-c()
islgc200reverse<-c()

#forward strand only


for (i in 1:length(isltsss)) {
  if (islstrands[i]=='forward') {
    islgc200forward<-rbind(islgc200forward, gc200matlab[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50)])
    #y<-gc200matlab[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50)]
    #x<-1:length(y)
    #lo<- loess(y~x)
    #lines(predict(lo), col='red', lwd=0.2)    
  } else {
    islgc200reverse<-rbind(islgc200reverse, rev(gc200matlab[(as.numeric(isltsss[i])-150):(as.numeric(isltsss[i])+50)]))
    #y<-rev(gc200matlab[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50)])
    #x<-1:length(y)
    #lo<-loess(y~x)
    #lines(predict(lo), col='red', lwd=0.2)    
  }
  
}


#islgc200forward<-islgc200forward[-1,] #удаление первой строки с нулями
#islgc200reverse<-islgc200reverse[-1,]
#dimnames(islgc200reverse)<-NULL #удаление аттрибутов
#dimnames(islgc200forward)<-NULL

islgc200<-rbind(islgc200forward, islgc200reverse)
#dimnames(islgc200)<-NULL


## lowscore

lowscore_names<-names(dataset_lowscore)
lowscoretsss<-c()
lowscoreseqs<-c
lowscorestrands<-c()
lowscorempots<-c()

#only a part of lowscore sequences is used
for (i in 1:2000){
  lowscoretsss<-c(lowscoretsss, dataset_lowscore[[i]]$tss)
  lowscoreseqs<-toupper(c(lowscoreseqs, dataset_lowscore[[i]]$seq))
  lowscorestrands<-c(lowscorestrands, dataset_lowscore[[i]]$strand)
  lowscorempots<-rbind(lowscorempots, dataset_lowscore[[i]]$mpot)
}

for (i in c('lowscoreaeos1forward', 'lowscoreaeos1reverse', 'lowscoreaeos3forward', 'lowscoreaeos3reverse')) {
  assign(i, c())
}
#lowscoreaeos1forward<-rep(0,201)
#lowscoreaeos1reverse<-rep(0,201)
#lowscoreaeos2forward<-rep(0,201)
#lowscoreaeos2reverse<-rep(0,201)
#lowscoreaeos3forward<-rep(0,201)
#lowscoreaeos3reverse<-rep(0,201)
#lowscoreaeos4forward<-rep(0,201)
#lowscoreaeos4reverse<-rep(0,201)
for (i in 1:2000) {
  if (lowscorestrands[i]=='forward') {
    lowscoreaeos1forward<-rbind(lowscoreaeos1forward, E01[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50)])
    #   lowscoreaeos2forward<-rbind(lowscoreaeos2forward, matr1[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50),2])
    lowscoreaeos3forward<-rbind(lowscoreaeos3forward, d1[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50)])
    #  lowscoreaeos4forward<-rbind(lowscoreaeos4forward, matr1[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50),4])            
  } else {
    lowscoreaeos1reverse<-rbind(lowscoreaeos1reverse, E02[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50)])
    # lowscoreaeos2reverse<-rbind(lowscoreaeos2reverse, matr1[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50),2])
    lowscoreaeos3reverse<-rbind(lowscoreaeos3reverse, d2[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50)])
    #lowscoreaeos4reverse<-rbind(lowscoreaeos4reverse, matr1[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50),4])
  }
}

#lowscoreaeos1forward<-lowscoreaeos1forward[-1,] #удаление первой строки с нулями
#lowscoreaeos1reverse<-lowscoreaeos1reverse[-1,]
#dimnames(lowscoreaeos1reverse)<-NULL #удаление аттрибутов
#dimnames(lowscoreaeos1forward)<-NULL

#lowscoreaeos2forward<-lowscoreaeos2forward[-1,] #СѓРґР°Р»РµРЅРёРµ РёСЃС…РѕРґРЅС‹С… РЅСѓР»РµР№
#lowscoreaeos2reverse<-lowscoreaeos2reverse[-1,]
#dimnames(lowscoreaeos2reverse)<-NULL
#dimnames(lowscoreaeos2forward)<-NULL

#lowscoreaeos3forward<-lowscoreaeos3forward[-1,] #СѓРґР°Р»РµРЅРёРµ РёСЃС…РѕРґРЅС‹С… РЅСѓР»РµР№
#lowscoreaeos3reverse<-lowscoreaeos3reverse[-1,]
#dimnames(lowscoreaeos3reverse)<-NULL
#dimnames(lowscoreaeos3forward)<-NULL

#lowscoreaeos4forward<-lowscoreaeos4forward[-1,] #СѓРґР°Р»РµРЅРёРµ РёСЃС…РѕРґРЅС‹С… РЅСѓР»РµР№
#lowscoreaeos4reverse<-lowscoreaeos4reverse[-1,]
#dimnames(lowscoreaeos4reverse)<-NULL
#dimnames(lowscoreaeos4forward)<-NULL

#РѕР±СЉРµРґРёРЅРµРЅРёРµ РґР°РЅРЅС‹С… РґР»СЏ РѕР±РѕРёС… СЃС‚СЂСЌРЅРґРѕРІ

lowscoreaeos1<-lowscoreaeos1forward
#lowscoreaeos2<-rbind(lowscoreaeos2forward, lowscoreaeos2reverse)
lowscoreaeos3<-lowscoreaeos3forward
#lowscoreaeos4<-rbind(lowscoreaeos4forward, lowscoreaeos4reverse)


lowscoregc200forward<-c()
lowscoregc200reverse<-c()

#forward stlowscore only


for (i in 1:length(lowscoretsss)) {
  if (lowscorestrands[i]=='forward') {
    lowscoregc200forward<-rbind(lowscoregc200forward, gc200matlab[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50)])
    #y<-gc200matlab[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50)]
    #x<-1:length(y)
    #lo<- loess(y~x)
    #lines(predict(lo), col='red', lwd=0.2)    
  } else {
    lowscoregc200reverse<-rbind(lowscoregc200reverse, rev(gc200matlab[(as.numeric(lowscoretsss[i])-150):(as.numeric(lowscoretsss[i])+50)]))
    #y<-rev(gc200matlab[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50)])
    #x<-1:length(y)
    #lo<-loess(y~x)
    #lines(predict(lo), col='red', lwd=0.2)    
  }
  
}


#lowscoregc200forward<-lowscoregc200forward[-1,] #удаление первой строки с нулями
#lowscoregc200reverse<-lowscoregc200reverse[-1,]
#dimnames(lowscoregc200reverse)<-NULL #удаление аттрибутов
#dimnames(lowscoregc200forward)<-NULL

lowscoregc200<-lowscoregc200forward




#setting names


rownames(aeos1)<-exp_proms[,1]
rownames(aeos3)<-exp_proms[,1]
rownames(mpots)<-exp_proms[,1]
rownames(gc200)<-exp_proms[,1]

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

habillage_5components_4props<-c(rep('P', nrow(exp_proms)), rep('N', length(dataset_notpro)), rep('G', length(dataset_gen)), rep('I', length(dataset_isl)), rep('L', 2000))

colnames(to_pca_5components_4props)<-NULL


#############3 #scale =F?
princ.return.5comps.4props <- prcomp(to_pca_5components_4props, scale=F)

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

#plot(pcs.from.var[1,], type='l', lwd=0.2)
#for (i in 1:nrow(pcs.from.var)){
#  lines(pcs.from.var[i,], type='l', lwd=0.2)
#}


rownames(pcs.from.var)<-rownames(to_pca_5components_4props)
colnames(pcs.from.var)<-1:ncol(pcs.from.var)

#factor_to_promoters_vs_lowscore<-as.factor(c(rep('Promoter', nrow(promoters)), rep('Lowscore', nrow(lowscore))))
#promoters_vs_lowscore <- cbind(factor_to_promoters_vs_lowscore, rbind(promoters, lowscore))



#habillage_5components_4props<-as.factor(c(rep('Promoter', nrow(exp_proms)), rep('Non-promoter', length(dataset_notpro)), rep('Gene', length(dataset_gen)), rep('Islands', length(dataset_isl)), rep('Lowscore', 2000)))
#mixt_5comps_after_pca_on_4_props_ae_size_mpots_gc200<-cbind(habillage_5components_4props, pcs.from.var)
mixt_5comps_after_pca_on_4_props_ae_size_mpots_gc200<- pcs.from.var
save(mixt_5comps_after_pca_on_4_props_ae_size_mpots_gc200, file='new_no_factor_8XII_mixt_5comps_after_pca_on_4_props_ae_size_mpots_gc200.Rdata')


eig.5comps.4props <- data.frame(eig = eig, variance = variance,
                                cumvariance = cumvar)

#saving initial profiles sets as pairs 'PROMOTERS - SMT ELSE"
pro.not.all.props<-rbind(cbind(aeos1, aeos3, mpots, gc200),
                         cbind(notaeos1, notaeos3, notmpots, notgc200))
rownames(pro.not.all.props)<-c(exp_proms[,1], rownames(notaeos1))
colnames(pro.not.all.props)<-NULL
save(pro.not.all.props, file='pro.not.all.props.Rdata')

pro.gen.all.props<-rbind(cbind(aeos1, aeos3, mpots, gc200),
                         cbind(genaeos1, genaeos3, genmpots, gengc200))
rownames(pro.gen.all.props)<-c(exp_proms[,1], rownames(genaeos1))
colnames(pro.gen.all.props)<-NULL
save(pro.gen.all.props, file='pro.gen.all.props.Rdata')

pro.isl.all.props<-rbind(cbind(aeos1, aeos3, mpots, gc200),
                         cbind(islaeos1, islaeos3, islmpots, islgc200))
rownames(pro.isl.all.props)<-c(exp_proms[,1], rownames(islaeos1))
colnames(pro.isl.all.props)<-NULL
save(pro.isl.all.props, file='pro.isl.all.props.Rdata')

pro.lowscore.all.props<-rbind(cbind(aeos1, aeos3, mpots, gc200),
                              cbind(lowscoreaeos1, lowscoreaeos3, lowscorempots, lowscoregc200))
rownames(pro.lowscore.all.props)<-c(exp_proms[,1], rownames(lowscoreaeos1))
colnames(pro.lowscore.all.props)<-NULL
save(pro.lowscore.all.props, file='pro.lowscore.all.props.Rdata')

#barplot(eig.5comps.4props [, 2], names.arg=1:nrow(eig.5comps.4props), 
 #       main = "Variances",
  #      xlab = "Principal Components",
   #     ylab = "Percentage of variances",
    #    col ="steelblue")
# Add connected line segments to the plot
#lines(x = 1:nrow(eig.5comps.4props ), 
 #     eig.5comps.4props [, 2], 
  #    type="b", pch=19, col = "red")

library(factoextra)

fviz_pca_ind(princ.return.5comps.4props, label="none", habillage=habillage_5components_4props, pointsize = 0.5, addEllipses = T, ellipse.level = 0.99)


library(rgl)

habillage.pro.not.isl.gen<-c(rep('red', nrow(exp_proms)), rep ('blue', length(dataset_notpro)), rep('green', length(dataset_gen)), rep('orange', length(dataset_isl)), rep('magenta', 2000)) 
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
load("new_no_factor_8XII_mixt_5comps_after_pca_on_4_props_ae_size_mpots_gc200.Rdata")
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
library(reldna)


for (i in c('bendforward', 'bendreverse')) {
  assign(i, c())
}

for (i in 1:length(exp_proms[,1])) {
  if (exp_strands[i]=='forward') {
    aeos1forward<-rbind(aeos1forward, E01[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50)])
    #aeos2forward<-rbind(aeos2forward, matr1[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50),2])
    aeos3forward<-rbind(aeos3forward, d1[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50)])
    #aeos4forward<-rbind(aeos4forward, matr1[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50),4])            
  } else {
    aeos1reverse<-rbind(aeos1reverse, E02[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50)])
    #aeos2reverse<-rbind(aeos2reverse, matr1[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50),2])
    aeos3reverse<-rbind(aeos3reverse, d2[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50)])
    #aeos4reverse<-rbind(aeos4reverse, matr1[(as.numeric(exp_proms[i,3])-150):(as.numeric(exp_proms[i,3])+50),4])
    
  }
  
}

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

aeos1<-rbind(aeos1forward, aeos1reverse)
#aeos2<-rbind(aeos2forward, aeos2reverse)
aeos3<-rbind(aeos3forward, aeos3reverse)
#aeos4<-rbind(aeos4forward, aeos4reverse)



