######################################
# the reference code of program2 
######################################

######################################
# initial
######################################

library("Biostrings",verbose=F,quietly=T)

# read parameters
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("USAGE: Rscript hw2_105753035.R --input test.fasta --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output result.fasta", call.=FALSE)
}

# parse parameters
i<-1 
while(i < length(args))
{
  if(args[i] == "--input"){
    i_f<-args[i+1]
    i<-i+1
  }else if(args[i] == "--score"){
    s_f<-args[i+1]
    i<-i+1
  }else if(args[i] == "--aln"){
    aln_mode <- args[i+1]
    i<-i+1
  }else if(args[i] == "--gap_open"){
    g_o<- as.integer(args[i+1]) 
    i<-i+1
  }else if(args[i] == "--gap_extend"){
    g_e<-as.integer(args[i+1])
    i<-i+1    
  }else if(args[i] == "--output"){
    o_f<-args[i+1]
    i<-i+1
  }else{
    stop(paste("Unknown flag", args[i]), call.=FALSE)
  }
  i<-i+1
}

print("PARAMETERS")
print(paste("input file         :", i_f))
print(paste("output file        :", o_f))
print(paste("score file         :", s_f))
print(paste("gap open penalty   :", g_o))
print(paste("gap extend penalty :", g_e))

######################################
# main
######################################
# read fasta file
ff <- readAAStringSet(i_f)
seq_name = names(ff)
sequence = paste(ff)


# aln length
aln_length <- nchar(sequence[1])

# read score file
s_m<-read.table(s_f)
s_m<-as.matrix(s_m)

#bulid matrix

nchx = nchar(ff[1])
repx = rep(0,nchx)
for(i in 1:nchx) repx[i] = substr(ff[1],start=i,stop=i)

nchy = nchar(ff[2])
repy = rep(0,nchy)
for(i in 1:nchy) repy[i] = substr(ff[2],start=i,stop=i)

F = matrix(0,nchx+1,nchy+1)
dimnames(F) = list(c("",repx),c("",repy))

for (i in 1:nchx) F[i+1,1] = i*g_e
for (j in 1:nchy) F[1,j+1] = j*g_e

for (i in 1:nchx)
{
  for (j in 1:nchy)
  {
    if (repx[i] == repy[j])
    {
      
      t1 = F[i,j] + s_m[repx[i],repy[j]]
    } else {
      t1 = F[i,j] + g_o
    }
    
    t2 = F[i,j+1] + g_e
    t3 = F[i+1,j] + g_e
    
    F[i+1,j+1] = max(t1,t2,t3)
    
  }
}

# trace back 

AlignmentA = ""
AlignmentB = ""
Final_a <- 2
Final_b <- 2

while ( Final_a < nchx + 1 || Final_b < nchy + 1 )
{
  Score <- F[nchx + 1,nchy + 1]
  ScoreDiag <- F[nchx , nchy ]
  ScoreUp <- F[nchx + 1, nchy]
  ScoreLeft <- F[nchx , nchy + 1]
  if (Score == ScoreDiag + s_m[repx[nchx],repy[nchy]])
  {
    AlignmentA <- paste(AlignmentA,repx[nchx], sep="") 
    AlignmentB <- paste(AlignmentB,repy[nchy], sep="")
    nchx <- nchx - 1
    nchy <- nchy - 1
    Final_a <- Final_a + 1
    Final_b <- Final_b + 1
  }
  else if (Score == ScoreLeft + as.numeric(g_e) )
  {
    AlignmentA <- paste(AlignmentA,repx[nchx], sep="") 
    AlignmentB <- paste(AlignmentB,"-", sep="")
    nchx <- nchx - 1
    Final_a <- Final_a + 1
  }
  else if (Score == ScoreUp + as.numeric(g_e) )
  {
    AlignmentA <- paste(AlignmentA,"-", sep="") 
    AlignmentB <- paste(AlignmentB,repy[nchy], sep="")
    nchy <- nchy - 1
    Final_b <- Final_b + 1
  }
  
}




aln_score<-0
for(i in 1:aln_length)
{
  a<-substring(sequence[1], i, i)
  b<-substring(sequence[2], i, i)
  
  if((a != "-")&&(b != "-"))
  {
    print(paste(a, "-", b, "=", s_m[a,b]))
    aln_score = aln_score + s_m[a,b]
  }
  else{
    aln_score = aln_score + g_o
  }
}

print(aln_score)

# output
writeXStringSet(ff, o_f)
write(AlignmentA,file=o_f,append=FALSE)
write(AlignmentB,file=o_f,append=TRUE)
print(AlignmentA)
print(AlignmentB)
