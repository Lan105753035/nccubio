# read PAM1 from data
pam1<-read.table(file="pam1.txt",header = TRUE)

# check PAM1 data
dim(pam1)
str(pam1)

# construct PAM250 from PAM1
pam250<-as.matrix(pam1)/10000
pam<-as.matrix(pam1)/10000

i<-0
while (i<=250){
  pam250<-as.matrix(pam250)%*%as.matrix(pam)
  i<-i+1
}

pam250<-pam250*100
pam250<-round(pam250,0)

# output PAM250 as a file
write.table(pam250,file = "pam250.txt")
