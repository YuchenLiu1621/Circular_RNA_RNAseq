setwd('/mnt/00D65E61D65E56CE/Deep_seq_22.3.15/1')
library.package <- function(){
  #library(readxl)
  library(stringr)
  library(Biostrings)
  library(GenomicAlignments)
}

#解压，fq转fa，不需要返回####################################
fq.2.fa <- function(){
  
  #fq->fa########################################################################################################
  
  # print("please enter *Path*")
  # path<-readline()
  # setwd(path)
  
  files.gz<-list.files(pattern = ".fq.gz")
  for (i in 1:length(files.gz)) {
    shell_cmd<-paste0("gzip -d ",files.gz[i])
    system(shell_cmd,intern = T)
    print(paste0("gzip",i))
  }
  files.fq<-list.files(pattern = ".fq")
  for (i in 1:length(files.fq)) {
    shell_cmd<-paste0("grep -A 1 \"^@\" ",files.fq[i],"  | grep -v \"^-\" | sed \'s/@/>/\' > ",
                      paste0(str_sub(files.fq[i],1,-2),"a"))
    system(shell_cmd,intern = T)
    print(paste0("fq->fa",i))
  }
}
#############################################################

#bwa序列比对，得到bam文件，不需要返回
bwa <- function(){
  #ref.info######################################################################################################
  
  ref.info <- read.csv(list.files(pattern = "ref.info.csv"),header = F,stringsAsFactors = F)
  ref.info <- na.omit(ref.info)
  colnames(ref.info) <- c("names","refseq","guideseq")
  names <- ref.info[,1]
  refseq <- ref.info[,2]
  guideseq <- ref.info[,3]
  ref.seq<-file("ref.seq.fa","w")
  writeLines(paste0(">",ref.info[,1],"\n",ref.info[,2]),ref.seq)
  close(ref.seq)
  s_e<-matrix(as.numeric(str_locate(ref.info[,2],ref.info[,3])),ncol = 2)
  colnames(s_e)<-c("Guide_start","Guide_end")
  ref.info<-data.frame(ref.info,s_e)
  
  #bwa###################################################################################################
  
  # print("please enter *Max length of index*")
  # index.max.length<-readline()
  # index.max.length<-as.numeric(index.max.length)
  file.index<-list.files(pattern = "index.info.csv")
  
  index.info<-read.csv(file.index,blank.lines.skip = T,stringsAsFactors = F,header = F)
  index.info<-na.omit(index.info)
  colnames(index.info)<-c("sample","id","r1_index","r2_index")
  index.max.length=max(nchar(c(index.info$r1_index,index.info$r2_index)))
  #r1
  file.r1<-list.files(pattern = "1.fq")
  r1_fq<-file("r1_del_index.fq","w")
  deep_clean_r1<-file(file.r1,"r")
  line_clean_r1<-readLines(deep_clean_r1,n=1)
  i=1
  while((length(line_clean_r1) != 0)){
    a<-(i+2) %% 4
    b<-i %% 4
    if(a == 0 || b == 0){
      
      del.index<-substr(line_clean_r1,(index.max.length+1),nchar(line_clean_r1))
      #cat(c(del.index,"\n"),file = "r1_del_index.fq",append = T)
      writeLines(del.index,r1_fq)
      
    }else{
      
      #cat(c(line_clean_r1,"\n"),file = "r1_del_index.fq",append = T)
      writeLines(line_clean_r1,r1_fq)
      
    }
    line_clean_r1<-readLines(deep_clean_r1,n=1)
    if(i %% 100000 == 0){
      print(i)
    }
    i<-i+1
  }
  close(deep_clean_r1)
  close(r1_fq)
  #r2
  file.r2<-list.files(pattern = "2.fq")
  r2_fq<-file("r2_del_index.fq","w")
  deep_clean_r2<-file(file.r2,"r")
  line_clean_r2<-readLines(deep_clean_r2,n=1)
  i=1
  while((length(line_clean_r2) != 0)){
    a<-(i+2) %% 4
    b<-i %% 4
    if(a == 0 || b == 0){
      
      del.index<-substr(line_clean_r2,(index.max.length+1),nchar(line_clean_r2))
      #cat(c(del.index,"\n"),file = "r2_del_index.fq",append = T)
      writeLines(del.index,r2_fq)
      
    }else{
      
      #cat(c(line_clean_r2,"\n"),file = "r2_del_index.fq",append = T)
      writeLines(line_clean_r2,r2_fq)
      
    }
    line_clean_r2<-readLines(deep_clean_r2,n=1)
    if(i %% 100000 == 0){
      print(i)
    }
    i<-i+1
  }
  close(deep_clean_r2)
  close(r2_fq)
  
  #file.ref<-list.files(pattern = "ref.seq.fa")
  shell_cmd<-"bwa index ref.seq.fa"
  system(shell_cmd,intern = T)
  print("bwa_index created")
  
  shell_cmd<-"bwa mem ref.seq.fa r1_del_index.fq r2_del_index.fq > bwa-mem-pairend.sam"
  system(shell_cmd,intern = T)
  print("sequences are aligned")
  shell_cmd<-"samtools view -S -b bwa-mem-pairend.sam > bwa-mem-pairend.bam"
  system(shell_cmd,intern = T)
  shell_cmd<-"samtools sort bwa-mem-pairend.bam -o bwa-mem-pairend-sorted.bam"
  system(shell_cmd,intern = T)
  print("BWA complete")
}

indel.and.base.editing <- function(){
  #ref.info######################################################################################################
  
  ref.info <- read.csv(list.files(pattern = "ref.info.csv"),header = F,stringsAsFactors = F)
  ref.info<-na.omit(ref.info)
  # ref.info$V2[which(nchar(ref.info$V2)+12 > 300)] <- paste0(substr(ref.info$V2[which(nchar(ref.info$V2)+12 > 300)],1,144),
  #                                                           substr(ref.info$V2[which(nchar(ref.info$V2)+12 > 300)],nchar(ref.info$V2[which(nchar(ref.info$V2)+12 > 300)])-144+1,nchar(ref.info$V2[which(nchar(ref.info$V2)+12 > 300)])))
  
  colnames(ref.info) <- c("names","refseq","guideseq","seq-condition","compare-to-300","guide-in")
  names <- ref.info[,1]
  refseq <- ref.info[,2]
  guideseq <- ref.info[,3]
  # ref.seq<-file("ref.seq.fa","w")
  # writeLines(paste0(">",ref.info[,1],"\n",ref.info[,2]),ref.seq)
  # close(ref.seq)
  s_e<-matrix(as.numeric(str_locate(ref.info[,2],ref.info[,3])),ncol = 2)
  colnames(s_e)<-c("Guide_start","Guide_end")
  ref.info<-data.frame(ref.info,s_e)
  ref.info$refseq_length <- nchar(ref.info$refseq)
  
  #index.info##########################################################################################################
  
  print("Creating index.info")
  file.index<-list.files(pattern = "index.info.csv")
  #index.info<-read_excel(file.index,col_names = F)
  #write.csv(index.info,"index.info.csv",row.names = F,col.names = F)
  index.info<-read.csv(file.index,blank.lines.skip = T,stringsAsFactors = F,header = F)
  index.info<-na.omit(index.info)
  colnames(index.info)<-c("sample","id","r1_index","r2_index")
  #index.info<-index.info[-which(is.na(index.info$sample)==T),]
  #index.info$sample[13:24]<-"IL1RN"
  #unique(index.info$...1)
  rev_pair<-function(x){
    r<-reverse(x)
    a<-unlist(strsplit(r,""))
    b<-cbind(c("A","T","G","C"),c("T","A","C","G"))
    c<-b[,2][match(a,b[,1])]
    d<-c()
    for (i in 1:length(c)) {
      d<-paste0(d,c[i])
    }
    return(d)
  }
  
  r1.exp<-c()
  for (i in 1:length(ref.info$names)) {
    sam.index.i<-which(index.info$sample == ref.info$names[i])
    r1.exp.i<-paste0(index.info$r1_index[sam.index.i],rep(substr(refseq[i],1,4),length(sam.index.i)))
    r1.exp<-c(r1.exp,r1.exp.i)
  }
  r2.exp<-c()
  for (i in 1:length(ref.info$names)) {
    sam.index.i<-which(index.info$sample == ref.info$names[i])
    r2.exp.i<-paste0(index.info$r2_index[sam.index.i],rep(rev_pair(substr(refseq[i],nchar(refseq[i])-3,nchar(refseq[i]))),length(sam.index.i)))
    r2.exp<-c(r2.exp,r2.exp.i)
  }
  index.info$r1_exp<-r1.exp
  index.info$r2_exp<-r2.exp
  index.info$r2_index_rev_pair<-sapply(index.info$r2_index,function(x){rev_pair(x)})
  r2.rev.pair.exp<-c()
  for (i in 1:length(ref.info$names)) {
    sam.index.i<-which(index.info$sample == ref.info$names[i])
    r2.rev.pair.exp.i<-paste0(index.info$r2_index_rev_pair[sam.index.i],rep(rev_pair(substr(refseq[i],nchar(refseq[i])-3,nchar(refseq[i]))),length(sam.index.i)))
    r2.rev.pair.exp<-c(r2.rev.pair.exp,r2.rev.pair.exp.i)
  }
  index.info$r2_index_rev_pair_exp<-r2.rev.pair.exp
  index.info$order<-c(1:nrow(index.info))
  
  #bam of r1 r2##########################################################################################################
  
  print("Reading bam of r1&r2")
  r1.seq<-readDNAStringSet(list.files(pattern = "1.fa$"))
  r2.seq<-readDNAStringSet(list.files(pattern = "2.fa$"))
  flag1<-scanBamFlag(isFirstMateRead = T,
                     isSecondMateRead = F,
                     isDuplicate = F,
                     isNotPassingQualityControls = F)
  param1<-ScanBamParam(flag = flag1,what = "seq")
  
  r1_bam<-readGAlignments("bwa-mem-pairend-sorted.bam",use.names = T,param = param1)
  
  flag2<-scanBamFlag(isFirstMateRead = F,
                     isSecondMateRead = T,
                     isDuplicate = F,
                     isNotPassingQualityControls = F)
  param2<-ScanBamParam(flag = flag2,what = "seq")
  
  r2_bam<-readGAlignments("bwa-mem-pairend-sorted.bam",use.names = T,param = param2)
  
  #mis== 1 | 0 #######################################################################################################
  
  print("please wait for *mis0and1*")
  index.in.seq.exp<-paste0(subseq(r1.seq,1,min(nchar(as.character(index.info$r1_exp)))),
                           subseq(r2.seq,1,min(nchar(as.character(index.info$r2_exp)))))
  index.in.seq.rev.exp<-paste0(subseq(r1.seq,1,min(nchar(as.character(index.info$r1_exp)))),
                               subseq(r2.seq,1,min(nchar(as.character(index.info$r2_index_rev_pair_exp)))))
  index.ori.exp<-paste0(substr(index.info$r1_exp,1,min(nchar(as.character(index.info$r1_exp)))),
                        substr(index.info$r2_exp,1,min(nchar(as.character(index.info$r2_exp)))))
  index.ori.rev.exp<-paste0(substr(index.info$r1_exp,1,min(nchar(as.character(index.info$r1_exp)))),
                            substr(index.info$r2_index_rev_pair_exp,1,min(nchar(as.character(index.info$r2_index_rev_pair_exp)))))
  #min(nchar(as.character(index.info$r1_exp)))
  
  sam.index.mis0.exp<-match(index.in.seq.exp,index.ori.exp)
  match.percent.exp=length(na.omit(sam.index.mis0.exp))/length(sam.index.mis0.exp)
  
  sam.index.mis0.rev.exp<-match(index.in.seq.rev.exp,index.ori.rev.exp)
  match.percent.rev.exp=length(na.omit(sam.index.mis0.rev.exp))/length(sam.index.mis0.rev.exp)
  
  if(match.percent.exp > match.percent.rev.exp){
    index.in.seq = index.in.seq.exp
    index.ori = index.ori.exp
    sam.index.mis0 = sam.index.mis0.exp
    match.percent = match.percent.exp
  }else{
    index.in.seq = index.in.seq.rev.exp
    index.ori = index.ori.rev.exp
    sam.index.mis0 = sam.index.mis0.rev.exp
    match.percent = match.percent.rev.exp
  }
  
  print(paste0('match percent of sam.index.mis0 = ',match.percent))
  sam.index.NA.mis1<-sapply(index.in.seq[is.na(sam.index.mis0)],function(x){res=as.matrix(stringDist(c(as.character(x),index.ori),method = "hamming"));
  searchres=res[,1];
  searchres=searchres[-1];
  minv=min(searchres);
  if(minv<2){return(which.min(searchres))}else{return(NA)}})
  
  sam.index.mis1<-sam.index.mis0
  sam.index.mis1[is.na(sam.index.mis1)]<-sam.index.NA.mis1
  print("mis0and1 complete")
  
  #cigar######################################################################################################
  
  #x:names of seqs
  extrace_cigar_r1<-function(x){
    readid<-sapply(x,function(x){strsplit(x," ")[[1]][1]})
    readid_bam.index<-match(names(r1_bam),readid)
    cigar<-as.character(cigar(r1_bam)[!is.na(readid_bam.index)])
    #need.index<-grep(names(sort(table(cigar),decreasing = T))[1],cigar,invert = T)#返回没有匹配的下标，是我们需要的indel
    #cigar.filter<-cigar[need.index]
    return(cigar)
  }
  extrace_cigar_r2<-function(x){
    readid<-sapply(x,function(x){strsplit(x," ")[[1]][1]})
    readid_bam.index<-match(names(r2_bam),readid)
    cigar<-as.character(cigar(r2_bam)[!is.na(readid_bam.index)])
    #need.index<-grep(names(sort(table(cigar),decreasing = T))[1],cigar,invert = T)#返回没有匹配的下标，是我们需要的indel
    #cigar.filter<-cigar[need.index]
    return(cigar)
  }
  #x:cigar
  d.i.num<-function(x){
    d.num.i<-sapply(as.character(x),function(x){d=unlist(str_extract_all(x,"\\d*D"));d.numi<-sum(na.omit(as.numeric(unlist(str_extract_all(d,"\\d*")))));return(d.numi)})
    d.num<-sum(d.num.i)
    i.num.i<-sapply(as.character(x),function(x){i=unlist(str_extract_all(x,"\\d*I"));i.numi<-sum(na.omit(as.numeric(unlist(str_extract_all(i,"\\d*")))));return(i.numi)})
    i.num<-sum(i.num.i)
    return(c(d.num,i.num))
  }
  # all.num<-function(x){
  #   all.num.i<-sapply(as.character(x),function(x){all=unlist(str_extract_all(x,"\\d*"));all.numi<-sum(na.omit(as.numeric(all)));return(all.numi)})
  #   all.num<-sum(all.num.i)
  #   return(all.num)
  # }
  
  ####################################################################################################
  
  #str_locate(refseq[1],guide_seq_1)
  #x<-("19M125S")
  #输入paste后的r1 cigar和r2 cigar（以空格分隔），起始、终止位置，对应的参考序列长度
  cigar.need<-function(x,n1,n2,length.of.ref){
    if(length.of.ref <= 300){
      r1 <- unlist(str_split(x," "))[1]
      a<-unlist(str_extract_all(r1,"\\d*"))
      a[which(a=="")]<-NA
      b<-as.numeric(na.omit(a))
      
      c<-unlist(str_extract_all(r1,"[A-Z]*"))
      c[which(c=="")]<-NA
      d<-as.character(na.omit(c))
      
      r2 <- unlist(str_split(x," "))[2]
      e<-unlist(str_extract_all(r2,"\\d*"))
      e[which(e=="")]<-NA
      f<-rev(as.numeric(na.omit(e)))
      
      g<-unlist(str_extract_all(r2,"[A-Z]*"))
      g[which(g=="")]<-NA
      h<-rev(as.character(na.omit(g)))
      
      cigar.want.r1 <- character()
      for (i in 1:length(b)) {
        cigar.seq.i <- paste0(rep(d[i],b[i]),collapse = "")
        cigar.want.r1 <- paste0(cigar.want.r1,cigar.seq.i)
      }
      
      cigar.want.r2 <- character()
      for (i in 1:length(f)) {
        cigar.seq.i <- paste0(rep(h[i],f[i]),collapse = "")
        cigar.want.r2 <- paste0(cigar.want.r2,cigar.seq.i)
      }
      
      cigar.r2.del.length <- nchar(cigar.want.r1)+nchar(cigar.want.r2)-length.of.ref
      cigar.union <- paste0(cigar.want.r1,substr(cigar.want.r2,cigar.r2.del.length+1,nchar(cigar.want.r2)))
      
      cigar.union <- substr(cigar.union,n1,n2)
      
      cigar.all.chara <- unlist(str_split(cigar.union,""))
      cigar.all.uni.chara <- cigar.all.chara[1]
      cigar.all.uni.chara.num <- c()
      for (i in 2:length(cigar.all.chara)) {
        if(cigar.all.chara[i] == cigar.all.chara[i-1]){
          cigar.all.uni.chara <- cigar.all.uni.chara
        }else{
          cigar.all.uni.chara <- c(cigar.all.uni.chara,cigar.all.chara[i])
          cigar.all.uni.chara.num <- c(cigar.all.uni.chara.num,i)
        }
      }
      cigar.all.uni.chara.num <- c(1,cigar.all.uni.chara.num,length(cigar.all.chara)+1)
      cigar.all.uni.chara.num.1 <- c()
      for (i in 2:length(cigar.all.uni.chara.num)) {
        cigar.all.uni.chara.num.1 <- c(cigar.all.uni.chara.num.1,cigar.all.uni.chara.num[i]-cigar.all.uni.chara.num[i-1])
      }
      
      cigar.union.return <- paste0(paste0(cigar.all.uni.chara.num.1,cigar.all.uni.chara),collapse = "")
      
      # cigar.need<-c()
      # for (i in 1:length(b)) {
      #   sum<-sum(b[1:i])
      #   if(sum>n){##
      #     num.i<-(b[i]-(n-sum(b[1:i-1])))##
      #     cigar.need<-paste0(na.omit(c(num.i,b[i+1:length(b)])),d[i:length(d)])
      #     break()
      #   }
      # }
      return(cigar.union.return)
    }else{
      r1 <- unlist(str_split(x," "))[1]
      a<-unlist(str_extract_all(r1,"\\d*"))
      a[which(a=="")]<-NA
      b<-as.numeric(na.omit(a))
      
      c<-unlist(str_extract_all(r1,"[A-Z]*"))
      c[which(c=="")]<-NA
      d<-as.character(na.omit(c))
      
      r2 <- unlist(str_split(x," "))[2]
      e<-unlist(str_extract_all(r2,"\\d*"))
      e[which(e=="")]<-NA
      f<-rev(as.numeric(na.omit(e)))
      
      g<-unlist(str_extract_all(r2,"[A-Z]*"))
      g[which(g=="")]<-NA
      h<-rev(as.character(na.omit(g)))
      
      cigar.want.r1 <- character()
      for (i in 1:length(b)) {
        cigar.seq.i <- paste0(rep(d[i],b[i]),collapse = "")
        cigar.want.r1 <- paste0(cigar.want.r1,cigar.seq.i)
      }
      
      cigar.want.r2 <- character()
      for (i in 1:length(f)) {
        cigar.seq.i <- paste0(rep(h[i],f[i]),collapse = "")
        cigar.want.r2 <- paste0(cigar.want.r2,cigar.seq.i)
      }
      cigar.length <- sum(b,f)
      if(length.of.ref-cigar.length >= 0){
        cigar.in <- paste0(paste0(rep("G",length.of.ref-cigar.length),collapse = ""),collapse = "")#在中间加G不影响统计I、D
        cigar.union <- paste0(cigar.want.r1,cigar.in,cigar.want.r2,collapse = "")
        
        cigar.union <- substr(cigar.union,n1,n2)
        
        cigar.all.chara <- unlist(str_split(cigar.union,""))
        cigar.all.uni.chara <- cigar.all.chara[1]
        cigar.all.uni.chara.num <- c()
        for (i in 2:length(cigar.all.chara)) {
          if(cigar.all.chara[i] == cigar.all.chara[i-1]){
            cigar.all.uni.chara <- cigar.all.uni.chara
          }else{
            cigar.all.uni.chara <- c(cigar.all.uni.chara,cigar.all.chara[i])
            cigar.all.uni.chara.num <- c(cigar.all.uni.chara.num,i)
          }
        }
        cigar.all.uni.chara.num <- c(1,cigar.all.uni.chara.num,length(cigar.all.chara)+1)
        cigar.all.uni.chara.num.1 <- c()
        for (i in 2:length(cigar.all.uni.chara.num)) {
          cigar.all.uni.chara.num.1 <- c(cigar.all.uni.chara.num.1,cigar.all.uni.chara.num[i]-cigar.all.uni.chara.num[i-1])
        }
        
        cigar.union.return <- paste0(paste0(cigar.all.uni.chara.num.1,cigar.all.uni.chara),collapse = "")
        return(cigar.union.return)
      }
      
    }
    
  }
  
  print("please wait for *Count d.i.all*")
  d.i.all<-c(0,0,0)
  sam.name<-c()
  for (j in 1:nrow(ref.info)) {
    sample.exist <- sort(unique(sam.index.mis1))
    min.i <- min(which(index.info$sample==ref.info$names[j]))
    max.i <- max(which(index.info$sample==ref.info$names[j]))
    for (i in min.i:max.i) {
      ref.length.of.i <- ref.info$refseq_length[j]
      guide.start.of.i <- ref.info$Guide_start[j]
      guide.end.of.i <- ref.info$Guide_end[j]
      index.length.of.i <- nchar(index.info$r1_index[i])
      seq.name.r1 <- names(r1.seq)[which(sam.index.mis1==index.info$order[i])]
      seq.name.r2 <- names(r2.seq)[which(sam.index.mis1==index.info$order[i])]
      cigar.i.r1 <- extrace_cigar_r1(seq.name.r1)
      cigar.i.r2 <- extrace_cigar_r2(seq.name.r2)
      cigar.i.all <- paste(cigar.i.r1,cigar.i.r2," ")
      n1 <- guide.start.of.i - 25
      n2 <- guide.end.of.i + 25
      if(length(cigar.i.r1) != 0){
        cigar.del<-sapply(cigar.i.all,function(x){cigar.need(x,n1,n2,ref.length.of.i)})
        d.i.all<-rbind(d.i.all,c(d.i.num(cigar.del),(ref.length.of.i + nchar(index.info$r1_index[i]) + nchar(index.info$r2_index[i])) * length(cigar.i.all)))
        sam.name<-c(sam.name,paste0(index.info$sample[i],"-",index.info$id[i]))
      }else{
        d.i.all<-rbind(d.i.all,c(0,0,0))
        sam.name<-c(sam.name,paste0(index.info$sample[i],"-",index.info$id[i]))
      }
      print(i)
    }
  }
  d.i.all.final<-d.i.all
  d.i.all.final<-d.i.all.final[-1,]
  rownames(d.i.all.final)<-sam.name
  colnames(d.i.all.final)<-c("Deletion","Insertion","All")
  write.csv(d.i.all.final,"Result.csv")
  
  ###############################################################################base editing
  
  #####1.数据按样本拆分
  #R1
  pass.ori <- getwd()
  
  r1.fq <- file(list.files(pattern = "1.fq"),"r")
  line.r1.fq <- readLines(r1.fq,n=1)
  i=1
  
  shell_cmd <- paste0("mkdir"," ",getwd(),"/Separate.R1")
  system(shell_cmd)
  setwd(paste0(pass.ori,"/Separate.R1"))
  all.sample.r1.file.to.write <- list()
  for (j in 1:nrow(index.info)) {
    all.sample.r1.file.to.write[[j]] <- file(paste0(index.info$sample[j],"-",index.info$id[j],".R1.fq"),"w")
  }
  while(length(line.r1.fq) != 0){
    seq.number <- ceiling(i/4)
    sample.number <- sam.index.mis1[seq.number]
    if(is.na(sample.number) == F){
      setwd(paste0(pass.ori,"/Separate.R1"))
      writeLines(line.r1.fq,all.sample.r1.file.to.write[[sample.number]])
    }
    setwd(pass.ori)
    line.r1.fq <- readLines(r1.fq,n=1)
    if(i%%100000 == 0){
      print(i)
    }
    i=i+1
  }
  print("##")
  for (j in 1:nrow(index.info)) {
    close(all.sample.r1.file.to.write[[j]])
  }
  close(r1.fq)
  
  #R2
  r2.fq <- file(list.files(pattern = "2.fq"),"r")
  line.r2.fq <- readLines(r2.fq,n=1)
  i=1
  
  shell_cmd <- paste0("mkdir"," ",getwd(),"/Separate.R2")
  system(shell_cmd)
  setwd(paste0(pass.ori,"/Separate.R2"))
  all.sample.r2.file.to.write <- list()
  for (j in 1:nrow(index.info)) {
    all.sample.r2.file.to.write[[j]] <- file(paste0(index.info$sample[j],"-",index.info$id[j],".R2.fq"),"w")
  }
  while(length(line.r2.fq) != 0){
    seq.number <- ceiling(i/4)
    sample.number <- sam.index.mis1[seq.number]
    if(is.na(sample.number) == F){
      setwd(paste0(pass.ori,"/Separate.R2"))
      writeLines(line.r2.fq,all.sample.r2.file.to.write[[sample.number]])
    }
    setwd(pass.ori)
    line.r2.fq <- readLines(r2.fq,n=1)
    if(i%%100000 == 0){
      print(i)
    }
    i=i+1
  }
  for (j in 1:nrow(index.info)) {
    close(all.sample.r2.file.to.write[[j]])
  }
  close(r2.fq)
  
  for (i in 1:nrow(index.info)) {
    if(ref.info$compare.to.300[match(index.info$sample[i],ref.info$names)] == ">" & ref.info$guide.in[match(index.info$sample[i],ref.info$names)] == "R2"){
      setwd(paste0(pass.ori,"/Separate.R2"))
      shell_cmd <- paste0("mkdir ",pass.ori,"/Separate.R2","/",index.info$sample[i],"-",index.info$id[i],"-R2-2-R1")
      system(shell_cmd)
      
      file.r2.i <- file(paste0(index.info$sample[i],"-",index.info$id[i],".R2.fq"),"r")
      file.r1.2.r1.i <- file(paste0(pass.ori,"/Separate.R2","/",index.info$sample[i],"-",index.info$id[i],"-R2-2-R1","/",index.info$sample[i],"-",index.info$id[i],"-R2-2-R1.fq"),"w")
      lines.r2 <- readLines(file.r2.i,n=1)
      h=1
      
      while(length(lines.r2 != 0)){
        if(h%%4 == 1){
          writeLines(lines.r2,file.r1.2.r1.i)
        }else if(h%%4 ==2){
          writeLines(rev_pair(lines.r2),file.r1.2.r1.i)
        }else{
          writeLines(reverse(lines.r2),file.r1.2.r1.i)
        }
        h=h+1
        lines.r2 <- readLines(file.r2.i,n=1)
      }
      
      close(file.r2.i)
      close(file.r1.2.r1.i)
    }
    print(i)
  }
  
  print("###")
  
  shell_cmd <- paste0("mkdir"," ",pass.ori,"/CRISPResso-res")
  system(shell_cmd)
  
  setwd(paste0(pass.ori,"/CRISPResso-res"))
  
  for (i in 1:nrow(index.info)) {
    if(ref.info$compare.to.300[match(index.info$sample[i],ref.info$names)] == "<"){
      r1.fq.i <- paste0(pass.ori,"/Separate.R1","/",index.info$sample[i],"-",index.info$id[i],".R1.fq")
      r2.fq.i <- paste0(pass.ori,"/Separate.R2","/",index.info$sample[i],"-",index.info$id[i],".R2.fq")
      shell_cmd <- paste0("CRISPResso"," --fastq_r1 ",r1.fq.i," --fastq_r2 ",r2.fq.i," --amplicon_seq ",
                          ref.info$refseq[which(ref.info$names==index.info$sample[i])],
                          " --guide_seq ",ref.info$guideseq[which(ref.info$names==index.info$sample[i])],
                          " --quantification_window_size 10 --quantification_window_center -10 --base_editor_output --min_paired_end_reads_overlap 1 --max_paired_end_reads_overlap 150")
      system(shell_cmd,intern = T)
    }else if(ref.info$compare.to.300[match(index.info$sample[i],ref.info$names)] == ">" & ref.info$guide.in[match(index.info$sample[i],ref.info$names)] == "R1"){
      r1.fq.i <- paste0(pass.ori,"/Separate.R1","/",index.info$sample[i],"-",index.info$id[i],".R1.fq")
      ref.seq.i <- substr(ref.info$refseq[which(ref.info$names==index.info$sample[i])],1,(150-nchar(index.info$r1_index[i])))
      shell_cmd <- paste0("CRISPResso"," --fastq_r1 ",r1.fq.i," --amplicon_seq ",
                          ref.seq.i,
                          " --guide_seq ",ref.info$guideseq[which(ref.info$names==index.info$sample[i])],
                          " --quantification_window_size 10 --quantification_window_center -10 --base_editor_output --min_paired_end_reads_overlap 1 --max_paired_end_reads_overlap 150")
      system(shell_cmd,intern = T)
    }else if(ref.info$compare.to.300[match(index.info$sample[i],ref.info$names)] == ">" & ref.info$guide.in[match(index.info$sample[i],ref.info$names)] == "R2"){
      r2.2.r1.fq <- paste0(pass.ori,"/Separate.R2","/",index.info$sample[i],"-",index.info$id[i],"-R2-2-R1","/",index.info$sample[i],"-",index.info$id[i],"-R2-2-R1.fq")
      ref.seq.i <- substr(ref.info$refseq[which(ref.info$names==index.info$sample[i])],
                          nchar(ref.info$refseq[which(ref.info$names==index.info$sample[i])])-(150-nchar(index.info$r2_index[i]))+1,
                          nchar(ref.info$refseq[which(ref.info$names==index.info$sample[i])]))
      shell_cmd <- paste0("CRISPResso"," --fastq_r1 ",r2.2.r1.fq," --amplicon_seq ",
                          ref.seq.i,
                          " --guide_seq ",ref.info$guideseq[which(ref.info$names==index.info$sample[i])],
                          " --quantification_window_size 10 --quantification_window_center -10 --base_editor_output --min_paired_end_reads_overlap 1 --max_paired_end_reads_overlap 150")
      system(shell_cmd,intern = T)
    }
    print(i)
  }
  
  #合并结果
  
  for (j in 1:nrow(ref.info)) {
    min.i <- min(which(index.info$sample==ref.info$names[j]))
    max.i <- max(which(index.info$sample==ref.info$names[j]))
    
    if(ref.info$compare.to.300[j] == "<"){
      condition.j <- c()
      have.i <- c()
      for (i in min.i : max.i) {
        files.i <- list.files(path = paste0(pass.ori,"/CRISPResso-res","/","CRISPResso_on_",
                                            index.info$sample[i],"-",index.info$id[i],".R1_",index.info$sample[i],"-",index.info$id[i],".R2"),
                              pattern = "Nucleotide_frequency_table.txt")
        condition.j <- c(condition.j,files.i)
        if(length(files.i) != 0){
          have.i <- c(have.i,i)
        }
        
      }
      
      if(length(condition.j) != 0){
        count.table.j <- matrix(rep(0,ref.info$refseq_length[j]*2),nrow = 2)
        percent.table.j <- matrix(rep(0,ref.info$refseq_length[j]*2),nrow = 2)
        
        colnames(count.table.j) <- colnames(read.table(paste0(pass.ori,"/CRISPResso-res","/","CRISPResso_on_",
                                                              index.info$sample[have.i[1]],"-",index.info$id[have.i[1]],".R1_",index.info$sample[have.i[1]],"-",index.info$id[have.i[1]],".R2",
                                                              "/","Nucleotide_frequency_table.txt")))
        colnames(percent.table.j) <- colnames(read.table(paste0(pass.ori,"/CRISPResso-res","/","CRISPResso_on_",
                                                                index.info$sample[have.i[1]],"-",index.info$id[have.i[1]],".R1_",index.info$sample[have.i[1]],"-",index.info$id[have.i[1]],".R2",
                                                                "/","Nucleotide_frequency_table.txt")))
        
        if(ref.info$seq.condition[j] == "positive"){
          namepart1 <- paste0("Site_",-(ref.info$Guide_start[j]-1):-1,"_",unlist(str_split(ref.info$refseq[j],""))[1:(ref.info$Guide_start[j]-1)])
          namepart2 <- paste0("PAM_",1:4,"_",unlist(str_split(ref.info$refseq[j],""))[ref.info$Guide_start[j]:(ref.info$Guide_start[j]+3)])
          namepart3 <- paste0("Guide_",1:(nchar(ref.info$guideseq[j])-4),"_",unlist(str_split(ref.info$refseq[j],""))[(ref.info$Guide_start[j]+4):ref.info$Guide_end[j]])
          namepart4 <- paste0("site_",1:(ref.info$refseq_length[j]-ref.info$Guide_end[j]),"_",unlist(str_split(ref.info$refseq[j],""))[(ref.info$Guide_end[j]+1):ref.info$refseq_length[j]])
          
        }else{
          namepart1 <- paste0("Site_",-(ref.info$refseq_length[j]-ref.info$Guide_end[j]):-1,"_",rev(unlist(str_split(ref.info$refseq[j],""))[(ref.info$Guide_end[j]+1):ref.info$refseq_length[j]]))
          namepart2 <- paste0("PAM_",1:4,"_",rev(unlist(str_split(ref.info$refseq[j],""))[(ref.info$Guide_end[j]-3):ref.info$Guide_end[j]]))
          namepart3 <- paste0("Guide_",1:(nchar(ref.info$guideseq[j])-4),"_",rev(unlist(str_split(ref.info$refseq[j],""))[ref.info$Guide_start[j]:(ref.info$Guide_end[j]-4)]))
          namepart4 <- paste0("site_",1:(ref.info$Guide_start[j]-1),"_",rev(unlist(str_split(ref.info$refseq[j],""))[1:(ref.info$Guide_start[j]-1)]))
          
        }
        rownames.of.table <- c()
        for (i in min.i : max.i) {
          condition.i <- list.files(path = paste0(pass.ori,"/CRISPResso-res","/","CRISPResso_on_",
                                                  index.info$sample[i],"-",index.info$id[i],".R1_",index.info$sample[i],"-",index.info$id[i],".R2"),
                                    pattern = "Nucleotide_frequency_table.txt")
          if(length(condition.i) != 0){
            count.table.j <- rbind(count.table.j,read.table(paste0(pass.ori,"/CRISPResso-res","/","CRISPResso_on_",
                                                                   index.info$sample[i],"-",index.info$id[i],".R1_",index.info$sample[i],"-",index.info$id[i],".R2",
                                                                   "/","Nucleotide_frequency_table.txt"),
                                                            header = T,sep = "\t",row.names = 1))
            percent.table.j <- rbind(percent.table.j,read.table(paste0(pass.ori,"/CRISPResso-res","/","CRISPResso_on_",
                                                                       index.info$sample[i],"-",index.info$id[i],".R1_",index.info$sample[i],"-",index.info$id[i],".R2",
                                                                       "/","Nucleotide_percentage_table.txt"),
                                                                header = T,sep = "\t",row.names = 1))
            rownames.of.table <- c(rownames.of.table,paste0(paste0(index.info$sample[i],"_",index.info$id[i]),"_",c("A","C","G","T","N","-")))
          }
          
          
        }
        if(ref.info$seq.condition[j] == "negative"){
          count.table.j <- count.table.j[c(-1,-2),ncol(count.table.j):1]
          percent.table.j <- percent.table.j[c(-1,-2),ncol(percent.table.j):1]
        }else{
          count.table.j <- count.table.j[c(-1,-2),]
          percent.table.j <- percent.table.j[c(-1,-2),]
        }
        colnames(count.table.j) <- c(namepart1,namepart2,namepart3,namepart4)
        colnames(percent.table.j) <- c(namepart1,namepart2,namepart3,namepart4)
        
        rownames(count.table.j) <- rownames.of.table
        rownames(percent.table.j) <- rownames.of.table
        
        shell_cmd <- paste0("mkdir ",pass.ori,"/",ref.info$names[j])
        system(shell_cmd)
        setwd(paste0(pass.ori,"/",ref.info$names[j]))
        write.csv(count.table.j,paste0(ref.info$names[j],"_frequency.csv"))
        write.csv(percent.table.j,paste0(ref.info$names[j],"_percentage.csv"))
        print(j)
      }
    }else if(ref.info$compare.to.300[j] == ">" & ref.info$guide.in[j] == "R2"){
      condition.j <- c()
      have.i <- c()
      for (i in min.i : max.i) {
        files.i <- list.files(path = paste0(pass.ori,"/CRISPResso-res","/","CRISPResso_on_",
                                            index.info$sample[i],"-",index.info$id[i],".R2"),
                              pattern = "Nucleotide_frequency_table.txt")
        condition.j <- c(condition.j,files.i)
        if(length(files.i) != 0){
          have.i <- c(have.i,i)
        }
        
      }
      if(length(condition.j) != 0){
        count.table.j <- matrix(rep(0,(150-nchar(index.info$r2_index[i]))*2),nrow = 2)
        percent.table.j <- matrix(rep(0,(150-nchar(index.info$r2_index[i]))*2),nrow = 2)
        
        colnames(count.table.j) <- colnames(read.table(paste0(pass.ori,"/CRISPResso-res","/","CRISPResso_on_",
                                                              index.info$sample[have.i[1]],"-",index.info$id[have.i[1]],"-R2-2-R1",
                                                              "/","Nucleotide_frequency_table.txt")))
        colnames(percent.table.j) <- colnames(read.table(paste0(pass.ori,"/CRISPResso-res","/","CRISPResso_on_",
                                                                index.info$sample[have.i[1]],"-",index.info$id[have.i[1]],"-R2-2-R1",
                                                                "/","Nucleotide_frequency_table.txt")))
        ref.seq.j <- paste0(str_extract(colnames(count.table.j),"[A-Z]"),collapse = "")
        guide.seq.j <- ref.info$guideseq[j]
        ref.length <- nchar(ref.seq.j)
        guide.start <- str_locate(ref.seq.j,guide.seq.j)[1]
        guide.end <- str_locate(ref.seq.j,guide.seq.j)[2]
        if(ref.info$seq.condition[j] == "positive"){
          namepart1 <- paste0("Site_",-(guide.start-1):-1,"_",unlist(str_split(ref.seq.j,""))[1:(guide.start-1)])
          namepart2 <- paste0("PAM_",1:4,"_",unlist(str_split(ref.seq.j,""))[guide.start:(guide.start+3)])
          namepart3 <- paste0("Guide_",1:(nchar(guide.seq.j)-4),"_",unlist(str_split(ref.seq.j,""))[(guide.start+4):guide.end])
          namepart4 <- paste0("site_",1:(ref.length-guide.end),"_",unlist(str_split(ref.seq.j,""))[(guide.end+1):ref.length])
          
        }else{
          namepart1 <- paste0("Site_",-(ref.length-guide.end):-1,"_",rev(unlist(str_split(ref.seq.j,""))[(guide.end+1):ref.length]))
          namepart2 <- paste0("PAM_",1:4,"_",rev(unlist(str_split(ref.seq.j,""))[(guide.end-3):guide.end]))
          namepart3 <- paste0("Guide_",1:(nchar(guide.seq.j)-4),"_",rev(unlist(str_split(ref.seq.j,""))[guide.start:(guide.end-4)]))
          namepart4 <- paste0("site_",1:(guide.start-1),"_",rev(unlist(str_split(ref.seq.j,""))[1:(guide.start-1)]))
          
        }
        rownames.of.table <- c()
        for (i in min.i : max.i) {
          condition.i <- list.files(path = paste0(pass.ori,"/CRISPResso-res","/","CRISPResso_on_",
                                                  index.info$sample[i],"-",index.info$id[i],"-R2-2-R1"),
                                    pattern = "Nucleotide_frequency_table.txt")
          if(length(condition.i) != 0){
            count.table.j <- rbind(count.table.j,read.table(paste0(pass.ori,"/CRISPResso-res","/","CRISPResso_on_",
                                                                   index.info$sample[i],"-",index.info$id[i],"-R2-2-R1",
                                                                   "/","Nucleotide_frequency_table.txt"),
                                                            header = T,sep = "\t",row.names = 1))
            percent.table.j <- rbind(percent.table.j,read.table(paste0(pass.ori,"/CRISPResso-res","/","CRISPResso_on_",
                                                                       index.info$sample[i],"-",index.info$id[i],"-R2-2-R1",
                                                                       "/","Nucleotide_percentage_table.txt"),
                                                                header = T,sep = "\t",row.names = 1))
            rownames.of.table <- c(rownames.of.table,paste0(paste0(index.info$sample[i],"_",index.info$id[i]),"_",c("A","C","G","T","N","-")))
          }
          
          
        }
        if(ref.info$seq.condition[j] == "negative"){
          count.table.j <- count.table.j[c(-1,-2),ncol(count.table.j):1]
          percent.table.j <- percent.table.j[c(-1,-2),ncol(percent.table.j):1]
        }else{
          count.table.j <- count.table.j[c(-1,-2),]
          percent.table.j <- percent.table.j[c(-1,-2),]
        }
        colnames(count.table.j) <- c(namepart1,namepart2,namepart3,namepart4)
        colnames(percent.table.j) <- c(namepart1,namepart2,namepart3,namepart4)
        
        rownames(count.table.j) <- rownames.of.table
        rownames(percent.table.j) <- rownames.of.table
        
        shell_cmd <- paste0("mkdir ",pass.ori,"/",ref.info$names[j])
        system(shell_cmd)
        setwd(paste0(pass.ori,"/",ref.info$names[j]))
        write.csv(count.table.j,paste0(ref.info$names[j],"_frequency.csv"))
        write.csv(percent.table.j,paste0(ref.info$names[j],"_percentage.csv"))
        print(j)
      }
    }else{
      condition.j <- c()
      have.i <- c()
      for (i in min.i : max.i) {
        files.i <- list.files(path = paste0(pass.ori,"/CRISPResso-res","/","CRISPResso_on_",
                                            index.info$sample[i],"-",index.info$id[i],".R1"),
                              pattern = "Nucleotide_frequency_table.txt")
        condition.j <- c(condition.j,files.i)
        if(length(files.i) != 0){
          have.i <- c(have.i,i)
        }
        
      }
      if(length(condition.j) != 0){
        count.table.j <- matrix(rep(0,(150-nchar(index.info$r1_index[i]))*2),nrow = 2)
        percent.table.j <- matrix(rep(0,(150-nchar(index.info$r1_index[i]))*2),nrow = 2)
        
        colnames(count.table.j) <- colnames(read.table(paste0(pass.ori,"/CRISPResso-res","/","CRISPResso_on_",
                                                              index.info$sample[have.i[1]],"-",index.info$id[have.i[1]],".R1",
                                                              "/","Nucleotide_frequency_table.txt")))
        colnames(percent.table.j) <- colnames(read.table(paste0(pass.ori,"/CRISPResso-res","/","CRISPResso_on_",
                                                                index.info$sample[have.i[1]],"-",index.info$id[have.i[1]],".R1",
                                                                "/","Nucleotide_frequency_table.txt")))
        ref.seq.j <- paste0(str_extract(colnames(count.table.j),"[A-Z]"),collapse = "")
        guide.seq.j <- ref.info$guideseq[j]
        ref.length <- nchar(ref.seq.j)
        guide.start <- str_locate(ref.seq.j,guide.seq.j)[1]
        guide.end <- str_locate(ref.seq.j,guide.seq.j)[2]
        if(ref.info$seq.condition[j] == "positive"){
          namepart1 <- paste0("Site_",-(guide.start-1):-1,"_",unlist(str_split(ref.seq.j,""))[1:(guide.start-1)])
          namepart2 <- paste0("PAM_",1:4,"_",unlist(str_split(ref.seq.j,""))[guide.start:(guide.start+3)])
          namepart3 <- paste0("Guide_",1:(nchar(guide.seq.j)-4),"_",unlist(str_split(ref.seq.j,""))[(guide.start+4):guide.end])
          namepart4 <- paste0("site_",1:(ref.length-guide.end),"_",unlist(str_split(ref.seq.j,""))[(guide.end+1):ref.length])
          
        }else{
          namepart1 <- paste0("Site_",-(ref.length-guide.end):-1,"_",rev(unlist(str_split(ref.seq.j,""))[(guide.end+1):ref.length]))
          namepart2 <- paste0("PAM_",1:4,"_",rev(unlist(str_split(ref.seq.j,""))[(guide.end-3):guide.end]))
          namepart3 <- paste0("Guide_",1:(nchar(guide.seq.j)-4),"_",rev(unlist(str_split(ref.seq.j,""))[guide.start:(guide.end-4)]))
          namepart4 <- paste0("site_",1:(guide.start-1),"_",rev(unlist(str_split(ref.seq.j,""))[1:(guide.start-1)]))
          
        }
        rownames.of.table <- c()
        for (i in min.i : max.i) {
          condition.i <- list.files(path = paste0(pass.ori,"/CRISPResso-res","/","CRISPResso_on_",
                                                  index.info$sample[i],"-",index.info$id[i],".R1"),
                                    pattern = "Nucleotide_frequency_table.txt")
          if(length(condition.i) != 0){
            count.table.j <- rbind(count.table.j,read.table(paste0(pass.ori,"/CRISPResso-res","/","CRISPResso_on_",
                                                                   index.info$sample[i],"-",index.info$id[i],".R1",
                                                                   "/","Nucleotide_frequency_table.txt"),
                                                            header = T,sep = "\t",row.names = 1))
            percent.table.j <- rbind(percent.table.j,read.table(paste0(pass.ori,"/CRISPResso-res","/","CRISPResso_on_",
                                                                       index.info$sample[i],"-",index.info$id[i],".R1",
                                                                       "/","Nucleotide_percentage_table.txt"),
                                                                header = T,sep = "\t",row.names = 1))
            rownames.of.table <- c(rownames.of.table,paste0(paste0(index.info$sample[i],"_",index.info$id[i]),"_",c("A","C","G","T","N","-")))
          }
          
          
        }
        if(ref.info$seq.condition[j] == "negative"){
          count.table.j <- count.table.j[c(-1,-2),ncol(count.table.j):1]
          percent.table.j <- percent.table.j[c(-1,-2),ncol(percent.table.j):1]
        }else{
          count.table.j <- count.table.j[c(-1,-2),]
          percent.table.j <- percent.table.j[c(-1,-2),]
        }
        colnames(count.table.j) <- c(namepart1,namepart2,namepart3,namepart4)
        colnames(percent.table.j) <- c(namepart1,namepart2,namepart3,namepart4)
        
        rownames(count.table.j) <- rownames.of.table
        rownames(percent.table.j) <- rownames.of.table
        
        shell_cmd <- paste0("mkdir ",pass.ori,"/",ref.info$names[j])
        system(shell_cmd)
        setwd(paste0(pass.ori,"/",ref.info$names[j]))
        write.csv(count.table.j,paste0(ref.info$names[j],"_frequency.csv"))
        write.csv(percent.table.j,paste0(ref.info$names[j],"_percentage.csv"))
        print(j)
      }  
    }
    
    
    
  }
  
  # #####indel
  # insertion.table <- c()
  # deletion.table <- c()
  # indel.table.rowname <- c()
  # for (j in 1:nrow(ref.info)){
  #   min.i <- min(which(index.info$sample==ref.info$names[j]))
  #   max.i <- max(which(index.info$sample==ref.info$names[j]))
  #   
  #   if(ref.info$compare.to.300[j] == "<"){
  #     condition.j <- c()
  #     for (i in min.i : max.i) {
  #       condition.j <- c(condition.j,list.files(path = paste0(pass.ori,"/CRISPResso-res","/","CRISPResso_on_",
  #                                                             index.info$sample[i],"-",index.info$id[i],".R1_",index.info$sample[i],"-",index.info$id[i],".R2"),
  #                                               pattern = "Nucleotide_frequency_table.txt"))
  #     }
  #     
  #     if(length(condition.j) != 0){
  #       indel.table.i <- read.table(paste0(pass.ori,"/CRISPResso-res","/","CRISPResso_on_",
  #                                          index.info$sample[i],"-",index.info$id[i],".R1_",index.info$sample[i],"-",index.info$id[i],".R2",
  #                                          "/","Indel_histogram.txt"),header = T,sep = "\t")
  #       deletion.table.i <- indel.table.i[indel.table.i$indel_size<0 & indel.table.i$fq != 0,]
  #       deletion.num.i <- -sum(deletion.table.i[,1]*deletion.table.i[,2])
  #       deletion.table <- c(deletion.table,deletion.num.i)
  #       
  #       insertion.table.i <- indel.table.i[indel.table.i$indel_size>0 & indel.table.i$fq != 0,]
  #       deletion.num.i <- sum(insertion.table.i[,1]*insertion.table.i[,2])
  #       deletion.table <- c(deletion.table,deletion.num.i)
  #     }
  #   }
  #   
  # }
  
}

only.indel=function(){
  ref.info <- read.csv(list.files(pattern = "ref.info.csv"),header = F,stringsAsFactors = F)
  ref.info<-na.omit(ref.info)
  # ref.info$V2[which(nchar(ref.info$V2)+12 > 300)] <- paste0(substr(ref.info$V2[which(nchar(ref.info$V2)+12 > 300)],1,144),
  #                                                           substr(ref.info$V2[which(nchar(ref.info$V2)+12 > 300)],nchar(ref.info$V2[which(nchar(ref.info$V2)+12 > 300)])-144+1,nchar(ref.info$V2[which(nchar(ref.info$V2)+12 > 300)])))
  
  colnames(ref.info) <- c("names","refseq","guideseq","seq-condition","compare-to-300","guide-in")
  names <- ref.info[,1]
  refseq <- ref.info[,2]
  guideseq <- ref.info[,3]
  # ref.seq<-file("ref.seq.fa","w")
  # writeLines(paste0(">",ref.info[,1],"\n",ref.info[,2]),ref.seq)
  # close(ref.seq)
  s_e<-matrix(as.numeric(str_locate(ref.info[,2],ref.info[,3])),ncol = 2)
  colnames(s_e)<-c("Guide_start","Guide_end")
  ref.info<-data.frame(ref.info,s_e)
  ref.info$refseq_length <- nchar(ref.info$refseq)
  
  #index.info##########################################################################################################
  
  print("Creating index.info")
  file.index<-list.files(pattern = "index.info.csv")
  #index.info<-read_excel(file.index,col_names = F)
  #write.csv(index.info,"index.info.csv",row.names = F,col.names = F)
  index.info<-read.csv(file.index,blank.lines.skip = T,stringsAsFactors = F,header = F)
  index.info<-na.omit(index.info)
  colnames(index.info)<-c("sample","id","r1_index","r2_index")
  #index.info<-index.info[-which(is.na(index.info$sample)==T),]
  #index.info$sample[13:24]<-"IL1RN"
  #unique(index.info$...1)
  rev_pair<-function(x){
    r<-reverse(x)
    a<-unlist(strsplit(r,""))
    b<-cbind(c("A","T","G","C"),c("T","A","C","G"))
    c<-b[,2][match(a,b[,1])]
    d<-c()
    for (i in 1:length(c)) {
      d<-paste0(d,c[i])
    }
    return(d)
  }
  
  r1.exp<-c()
  for (i in 1:length(ref.info$names)) {
    sam.index.i<-which(index.info$sample == ref.info$names[i])
    r1.exp.i<-paste0(index.info$r1_index[sam.index.i],rep(substr(refseq[i],1,4),length(sam.index.i)))
    r1.exp<-c(r1.exp,r1.exp.i)
  }
  r2.exp<-c()
  for (i in 1:length(ref.info$names)) {
    sam.index.i<-which(index.info$sample == ref.info$names[i])
    r2.exp.i<-paste0(index.info$r2_index[sam.index.i],rep(rev_pair(substr(refseq[i],nchar(refseq[i])-3,nchar(refseq[i]))),length(sam.index.i)))
    r2.exp<-c(r2.exp,r2.exp.i)
  }
  index.info$r1_exp<-r1.exp
  index.info$r2_exp<-r2.exp
  index.info$r2_index_rev_pair<-sapply(index.info$r2_index,function(x){rev_pair(x)})
  r2.rev.pair.exp<-c()
  for (i in 1:length(ref.info$names)) {
    sam.index.i<-which(index.info$sample == ref.info$names[i])
    r2.rev.pair.exp.i<-paste0(index.info$r2_index_rev_pair[sam.index.i],rep(rev_pair(substr(refseq[i],nchar(refseq[i])-3,nchar(refseq[i]))),length(sam.index.i)))
    r2.rev.pair.exp<-c(r2.rev.pair.exp,r2.rev.pair.exp.i)
  }
  index.info$r2_index_rev_pair_exp<-r2.rev.pair.exp
  index.info$order<-c(1:nrow(index.info))
  
  #bam of r1 r2##########################################################################################################
  
  print("Reading bam of r1&r2")
  r1.seq<-readDNAStringSet(list.files(pattern = "1.fa$"))
  r2.seq<-readDNAStringSet(list.files(pattern = "2.fa$"))
  flag1<-scanBamFlag(isFirstMateRead = T,
                     isSecondMateRead = F,
                     isDuplicate = F,
                     isNotPassingQualityControls = F)
  param1<-ScanBamParam(flag = flag1,what = "seq")
  
  r1_bam<-readGAlignments("bwa-mem-pairend-sorted.bam",use.names = T,param = param1)
  
  flag2<-scanBamFlag(isFirstMateRead = F,
                     isSecondMateRead = T,
                     isDuplicate = F,
                     isNotPassingQualityControls = F)
  param2<-ScanBamParam(flag = flag2,what = "seq")
  
  r2_bam<-readGAlignments("bwa-mem-pairend-sorted.bam",use.names = T,param = param2)
  
  #mis== 1 | 0 #######################################################################################################
  
  print("please wait for *mis0and1*")
  index.in.seq.exp<-paste0(subseq(r1.seq,1,min(nchar(as.character(index.info$r1_exp)))),
                           subseq(r2.seq,1,min(nchar(as.character(index.info$r2_exp)))))
  index.in.seq.rev.exp<-paste0(subseq(r1.seq,1,min(nchar(as.character(index.info$r1_exp)))),
                               subseq(r2.seq,1,min(nchar(as.character(index.info$r2_index_rev_pair_exp)))))
  index.ori.exp<-paste0(substr(index.info$r1_exp,1,min(nchar(as.character(index.info$r1_exp)))),
                        substr(index.info$r2_exp,1,min(nchar(as.character(index.info$r2_exp)))))
  index.ori.rev.exp<-paste0(substr(index.info$r1_exp,1,min(nchar(as.character(index.info$r1_exp)))),
                            substr(index.info$r2_index_rev_pair_exp,1,min(nchar(as.character(index.info$r2_index_rev_pair_exp)))))
  #min(nchar(as.character(index.info$r1_exp)))
  
  sam.index.mis0.exp<-match(index.in.seq.exp,index.ori.exp)
  match.percent.exp=length(na.omit(sam.index.mis0.exp))/length(sam.index.mis0.exp)
  
  sam.index.mis0.rev.exp<-match(index.in.seq.rev.exp,index.ori.rev.exp)
  match.percent.rev.exp=length(na.omit(sam.index.mis0.rev.exp))/length(sam.index.mis0.rev.exp)
  
  if(match.percent.exp > match.percent.rev.exp){
    index.in.seq = index.in.seq.exp
    index.ori = index.ori.exp
    sam.index.mis0 = sam.index.mis0.exp
    match.percent = match.percent.exp
  }else{
    index.in.seq = index.in.seq.rev.exp
    index.ori = index.ori.rev.exp
    sam.index.mis0 = sam.index.mis0.rev.exp
    match.percent = match.percent.rev.exp
  }
  
  print(paste0('match percent of sam.index.mis0 = ',match.percent))
  sam.index.NA.mis1<-sapply(index.in.seq[is.na(sam.index.mis0)],function(x){res=as.matrix(stringDist(c(as.character(x),index.ori),method = "hamming"));
  searchres=res[,1];
  searchres=searchres[-1];
  minv=min(searchres);
  if(minv<2){return(which.min(searchres))}else{return(NA)}})
  
  sam.index.mis1<-sam.index.mis0
  sam.index.mis1[is.na(sam.index.mis1)]<-sam.index.NA.mis1
  print("mis0and1 complete")
  
  #cigar######################################################################################################
  
  #x:names of seqs
  extrace_cigar_r1<-function(x){
    readid<-sapply(x,function(x){strsplit(x," ")[[1]][1]})
    readid_bam.index<-match(names(r1_bam),readid)
    cigar<-as.character(cigar(r1_bam)[!is.na(readid_bam.index)])
    #need.index<-grep(names(sort(table(cigar),decreasing = T))[1],cigar,invert = T)#返回没有匹配的下标，是我们需要的indel
    #cigar.filter<-cigar[need.index]
    return(cigar)
  }
  extrace_cigar_r2<-function(x){
    readid<-sapply(x,function(x){strsplit(x," ")[[1]][1]})
    readid_bam.index<-match(names(r2_bam),readid)
    cigar<-as.character(cigar(r2_bam)[!is.na(readid_bam.index)])
    #need.index<-grep(names(sort(table(cigar),decreasing = T))[1],cigar,invert = T)#返回没有匹配的下标，是我们需要的indel
    #cigar.filter<-cigar[need.index]
    return(cigar)
  }
  #x:cigar
  d.i.num<-function(x){
    d.num.i<-sapply(as.character(x),function(x){d=unlist(str_extract_all(x,"\\d*D"));d.numi<-sum(na.omit(as.numeric(unlist(str_extract_all(d,"\\d*")))));return(d.numi)})
    d.num<-sum(d.num.i)
    i.num.i<-sapply(as.character(x),function(x){i=unlist(str_extract_all(x,"\\d*I"));i.numi<-sum(na.omit(as.numeric(unlist(str_extract_all(i,"\\d*")))));return(i.numi)})
    i.num<-sum(i.num.i)
    return(c(d.num,i.num))
  }
  # all.num<-function(x){
  #   all.num.i<-sapply(as.character(x),function(x){all=unlist(str_extract_all(x,"\\d*"));all.numi<-sum(na.omit(as.numeric(all)));return(all.numi)})
  #   all.num<-sum(all.num.i)
  #   return(all.num)
  # }
  
  ####################################################################################################
  
  #str_locate(refseq[1],guide_seq_1)
  #x<-("19M125S")
  #输入paste后的r1 cigar和r2 cigar（以空格分隔），起始、终止位置，对应的参考序列长度
  cigar.need<-function(x,n1,n2,length.of.ref){
    if(length.of.ref <= 300){
      r1 <- unlist(str_split(x," "))[1]
      a<-unlist(str_extract_all(r1,"\\d*"))
      a[which(a=="")]<-NA
      b<-as.numeric(na.omit(a))
      
      c<-unlist(str_extract_all(r1,"[A-Z]*"))
      c[which(c=="")]<-NA
      d<-as.character(na.omit(c))
      
      r2 <- unlist(str_split(x," "))[2]
      e<-unlist(str_extract_all(r2,"\\d*"))
      e[which(e=="")]<-NA
      f<-rev(as.numeric(na.omit(e)))
      
      g<-unlist(str_extract_all(r2,"[A-Z]*"))
      g[which(g=="")]<-NA
      h<-rev(as.character(na.omit(g)))
      
      cigar.want.r1 <- character()
      for (i in 1:length(b)) {
        cigar.seq.i <- paste0(rep(d[i],b[i]),collapse = "")
        cigar.want.r1 <- paste0(cigar.want.r1,cigar.seq.i)
      }
      
      cigar.want.r2 <- character()
      for (i in 1:length(f)) {
        cigar.seq.i <- paste0(rep(h[i],f[i]),collapse = "")
        cigar.want.r2 <- paste0(cigar.want.r2,cigar.seq.i)
      }
      
      cigar.r2.del.length <- nchar(cigar.want.r1)+nchar(cigar.want.r2)-length.of.ref
      cigar.union <- paste0(cigar.want.r1,substr(cigar.want.r2,cigar.r2.del.length+1,nchar(cigar.want.r2)))
      
      cigar.union <- substr(cigar.union,n1,n2)
      
      cigar.all.chara <- unlist(str_split(cigar.union,""))
      cigar.all.uni.chara <- cigar.all.chara[1]
      cigar.all.uni.chara.num <- c()
      for (i in 2:length(cigar.all.chara)) {
        if(cigar.all.chara[i] == cigar.all.chara[i-1]){
          cigar.all.uni.chara <- cigar.all.uni.chara
        }else{
          cigar.all.uni.chara <- c(cigar.all.uni.chara,cigar.all.chara[i])
          cigar.all.uni.chara.num <- c(cigar.all.uni.chara.num,i)
        }
      }
      cigar.all.uni.chara.num <- c(1,cigar.all.uni.chara.num,length(cigar.all.chara)+1)
      cigar.all.uni.chara.num.1 <- c()
      for (i in 2:length(cigar.all.uni.chara.num)) {
        cigar.all.uni.chara.num.1 <- c(cigar.all.uni.chara.num.1,cigar.all.uni.chara.num[i]-cigar.all.uni.chara.num[i-1])
      }
      
      cigar.union.return <- paste0(paste0(cigar.all.uni.chara.num.1,cigar.all.uni.chara),collapse = "")
      
      # cigar.need<-c()
      # for (i in 1:length(b)) {
      #   sum<-sum(b[1:i])
      #   if(sum>n){##
      #     num.i<-(b[i]-(n-sum(b[1:i-1])))##
      #     cigar.need<-paste0(na.omit(c(num.i,b[i+1:length(b)])),d[i:length(d)])
      #     break()
      #   }
      # }
      return(cigar.union.return)
    }else{
      r1 <- unlist(str_split(x," "))[1]
      a<-unlist(str_extract_all(r1,"\\d*"))
      a[which(a=="")]<-NA
      b<-as.numeric(na.omit(a))
      
      c<-unlist(str_extract_all(r1,"[A-Z]*"))
      c[which(c=="")]<-NA
      d<-as.character(na.omit(c))
      
      r2 <- unlist(str_split(x," "))[2]
      e<-unlist(str_extract_all(r2,"\\d*"))
      e[which(e=="")]<-NA
      f<-rev(as.numeric(na.omit(e)))
      
      g<-unlist(str_extract_all(r2,"[A-Z]*"))
      g[which(g=="")]<-NA
      h<-rev(as.character(na.omit(g)))
      
      cigar.want.r1 <- character()
      for (i in 1:length(b)) {
        cigar.seq.i <- paste0(rep(d[i],b[i]),collapse = "")
        cigar.want.r1 <- paste0(cigar.want.r1,cigar.seq.i)
      }
      
      cigar.want.r2 <- character()
      for (i in 1:length(f)) {
        cigar.seq.i <- paste0(rep(h[i],f[i]),collapse = "")
        cigar.want.r2 <- paste0(cigar.want.r2,cigar.seq.i)
      }
      cigar.length <- sum(b,f)
      if(length.of.ref-cigar.length >= 0){
        cigar.in <- paste0(paste0(rep("G",length.of.ref-cigar.length),collapse = ""),collapse = "")#在中间加G不影响统计I、D
        cigar.union <- paste0(cigar.want.r1,cigar.in,cigar.want.r2,collapse = "")
        
        cigar.union <- substr(cigar.union,n1,n2)
        
        cigar.all.chara <- unlist(str_split(cigar.union,""))
        cigar.all.uni.chara <- cigar.all.chara[1]
        cigar.all.uni.chara.num <- c()
        for (i in 2:length(cigar.all.chara)) {
          if(cigar.all.chara[i] == cigar.all.chara[i-1]){
            cigar.all.uni.chara <- cigar.all.uni.chara
          }else{
            cigar.all.uni.chara <- c(cigar.all.uni.chara,cigar.all.chara[i])
            cigar.all.uni.chara.num <- c(cigar.all.uni.chara.num,i)
          }
        }
        cigar.all.uni.chara.num <- c(1,cigar.all.uni.chara.num,length(cigar.all.chara)+1)
        cigar.all.uni.chara.num.1 <- c()
        for (i in 2:length(cigar.all.uni.chara.num)) {
          cigar.all.uni.chara.num.1 <- c(cigar.all.uni.chara.num.1,cigar.all.uni.chara.num[i]-cigar.all.uni.chara.num[i-1])
        }
        
        cigar.union.return <- paste0(paste0(cigar.all.uni.chara.num.1,cigar.all.uni.chara),collapse = "")
        return(cigar.union.return)
      }
      
    }
    
  }
  
  print("please wait for *Count d.i.all*")
  d.i.all<-c(0,0,0)
  sam.name<-c()
  for (j in 1:nrow(ref.info)) {
    sample.exist <- sort(unique(sam.index.mis1))
    min.i <- min(which(index.info$sample==ref.info$names[j]))
    max.i <- max(which(index.info$sample==ref.info$names[j]))
    for (i in min.i:max.i) {
      ref.length.of.i <- ref.info$refseq_length[j]
      guide.start.of.i <- ref.info$Guide_start[j]
      guide.end.of.i <- ref.info$Guide_end[j]
      index.length.of.i <- nchar(index.info$r1_index[i])
      seq.name.r1 <- names(r1.seq)[which(sam.index.mis1==index.info$order[i])]
      seq.name.r2 <- names(r2.seq)[which(sam.index.mis1==index.info$order[i])]
      cigar.i.r1 <- extrace_cigar_r1(seq.name.r1)
      cigar.i.r2 <- extrace_cigar_r2(seq.name.r2)
      cigar.i.all <- paste(cigar.i.r1,cigar.i.r2," ")
      n1 <- guide.start.of.i - 25
      n2 <- guide.end.of.i + 25
      if(length(cigar.i.r1) != 0){
        cigar.del<-sapply(cigar.i.all,function(x){cigar.need(x,n1,n2,ref.length.of.i)})
        d.i.all<-rbind(d.i.all,c(d.i.num(cigar.del),(ref.length.of.i + nchar(index.info$r1_index[i]) + nchar(index.info$r2_index[i])) * length(cigar.i.all)))
        sam.name<-c(sam.name,paste0(index.info$sample[i],"-",index.info$id[i]))
      }else{
        d.i.all<-rbind(d.i.all,c(0,0,0))
        sam.name<-c(sam.name,paste0(index.info$sample[i],"-",index.info$id[i]))
      }
      print(i)
    }
  }
  d.i.all.final<-d.i.all
  d.i.all.final<-d.i.all.final[-1,]
  rownames(d.i.all.final)<-sam.name
  colnames(d.i.all.final)<-c("Deletion","Insertion","All")
  write.csv(d.i.all.final,"Result.csv")
}

library.package()
fq.2.fa()

bwa()
indel.and.base.editing()


for (i in 1:4) {
  if(i==1){
    setwd(paste0('/mnt/00D65E61D65E56CE/Deep_seq_22.3.15/',i))
    
    
    
    only.indel()
  }else{
    setwd(paste0('/mnt/00D65E61D65E56CE/Deep_seq_22.3.15/',i))
    fq.2.fa()
    
    bwa()
    
    only.indel()
  }
  
  print(i)
}



















