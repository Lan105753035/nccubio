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

nx = nchar(ff[1])
xx = rep(0,nx)
for(i in 1:nx) xx[i] = substr(ff[1],start=i,stop=i)

ny = nchar(ff[2])
yy = rep(0,ny)
for(i in 1:ny) yy[i] = substr(ff[2],start=i,stop=i)

F = matrix(0,nx+1,ny+1)
dimnames(F) = list(c("",xx),c("",yy))

for (i in 1:nx) F[i+1,1] = i*g_e
for (j in 1:ny) F[1,j+1] = j*g_e

for (i in 1:nx)
{
  for (j in 1:ny)
  {
    if (xx[i] == yy[j])
    {
      
      t1 = F[i,j] + s_m[xx[i],yy[j]]
    } else {
      t1 = F[i,j] + g_o
    }
    
    t2 = F[i,j+1] + g_e
    t3 = F[i+1,j] + g_e
    
    F[i+1,j+1] = max(t1,t2,t3)
    
  }
}
xx[nx]
yy[ny]
s_m[xx[nx],yy[ny]]
F[nx + 1,ny + 1]





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
