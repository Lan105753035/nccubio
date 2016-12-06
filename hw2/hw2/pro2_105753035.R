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
  stop("USAGE: Rscript pro2_<your student ID>.R --input test.fasta --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output result.fasta", call.=FALSE)
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
print(paste("aln_mode           :", aln_mode))
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

# mapping matrix
x=nchar(sequence[1])
repx = rep(0,x)
for(i in 1:x) repx[i] = substr(sequence[1],start=i,stop=i)
y=nchar(sequence[2])
repy = rep(0,y)
for(j in 1:y) repy[j] = substr(sequence[2],start=j,stop=j)
M = matrix(0,x+1,y+1)

for (i in 1:x)
{
  for (j in 1:y){
    if(aln_mode=="global"){
      for (a in 1:x) M[a+1,1] = g_o+(a-1)*g_e
      for (b in 1:y) M[1,b+1] = g_o+(b-1)*g_e
    if(repx[i]==repy[j])
      {aboveleft=M[i,j]+s_m[repx[i],repy[j]]}
    else{aboveleft=M[i,j]+g_o}
    left=M[i,j+1]+g_e
    above=M[i+1,j]+g_e
    M[i+1,j+1]=max(aboveleft,left,above)
    }
    if(aln_mode=="local"){
      if(repx[i]==repy[j])
      {aboveleft=M[i,j]+s_m[repx[i],repy[j]]}
      else{aboveleft=M[i,j]+g_o}
      left=M[i,j+1]+g_e
      above=M[i+1,j]+g_e
      M[i+1,j+1]=max(aboveleft,left,above,0)
    }
  }
}

write.table(M,file="matrix.txt")
# trace back

align1Buf = ""
align2Buf = ""
xaxis<-1
yaxis<-1

while ( xaxis < x + 1 || yaxis < y + 1 )
{
  currentCell <- M[x + 1,y + 1]
  saboveleft <- M[x , y ]
  sabove <- M[x + 1, y]
  sleft <- M[x , y + 1]
  if (currentCell == saboveleft + s_m[repx[x],repy[y]])
  {
    align1Buf <- paste(align1Buf,repx[x]) 
    align2Buf <- paste(align2Buf,repy[y])
    x <- x - 1
    y <- y - 1
  }
  else if (currentCell == sleft + as.numeric(g_e) )
  {
    align1Buf <- paste(align1Buf,repx[x]) 
    align2Buf <- paste(align2Buf,"-")
    x <- x - 1
    yaxis<-yaxis+1
  }
  else if (currentCell == sabove + as.numeric(g_e) )
  {
    align1Buf <- paste(align1Buf,"-") 
    align2Buf <- paste(align2Buf,repy[y])
    y <- y - 1
    xaxis<-xaxis+1
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
    aln_score = aln_score + g_O
  }
}

print(aln_score)

# output
writeXStringSet(ff, o_f)
write(align1Buf,file=o_f,append=FALSE)
write(align2Buf,file=o_f,append=TRUE)
print(align1Buf)
print(align2Buf)
