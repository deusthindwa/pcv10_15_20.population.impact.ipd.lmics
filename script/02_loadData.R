#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================


```{r setupdata, echo=FALSE}
d1<-read.csv("C:/Users/dmw63/Desktop/My documents h/SSI/data jan 2018/data_jan2018.csv")

range(d1$year) #range of years

d1$agec<-NA
d1$agec[0<=d1$agey_round & d1$agey_round<5] <-1
d1$agec[5<=d1$agey_round & d1$agey_round<18] <-2
d1$agec[18<=d1$agey_round & d1$agey_round<40] <-3
d1$agec[40<=d1$agey_round & d1$agey_round<65] <-4
d1$agec[65<=d1$agey_round & d1$agey_round<80] <-5
d1$agec[80<=d1$agey_round & d1$agey_round<110] <-6

agelabs<-c('<5y','5-17y','18-39y','40-64y','65-79y','80+y')

pcv7st<-c('4','6A','6B','9V','18C','19F','23F') #group 6A w pcv7
pcv13st<-c('1','3','5','6C','7F','19A' ) #group 6C with PCV13
pcv15st<-c( '22F','33F')
pcv20st<-c( '8','10A','11A','12F','15B','22F','33F') #Patent by Pfizer (Watson)

d1$pcv7<-0
d1$pcv7[d1$serotype_1 %in% pcv7st]<-1
d1$pcv13<-0
d1$pcv13[d1$serotype_1 %in% pcv13st]<-1
d1$pcv15<-0
d1$pcv15[d1$serotype_1 %in% pcv15st]<-1
d1$pcv20<-0
d1$pcv20[d1$serotype_1 %in% pcv20st]<-1

d1$pcv<-0
d1$pcv[d1$pcv7==1]<-7
d1$pcv[d1$pcv13==1]<-13
d1$pcv[d1$pcv15==1]<-15
d1$pcv[d1$pcv20==1]<-20
```

```{r tabulate}

#table without PCV15
simp.table<-as.array(table(d1$year,d1$pcv, d1$agec))
plot.cols<-c('black','#66c2a5','#fc8d62','#8da0cb')
for(i in 1:dim(simp.table)[3]){
  #pct of IPD covered by pcv20
  pct.pcv7<-paste0(round(100*t(simp.table[,2,i]) /colSums(t(simp.table[,,i]))),"%")
  pct.pcv13<-paste0(round(100*t(simp.table[,3,i]) /colSums(t(simp.table[,,i]))),"%")
  pct.pcv20<-paste0(round(100*t(simp.table[,4,i]) /colSums(t(simp.table[,,i]))),"%")
  
  xx<-barplot( t(simp.table[,,i]), main=paste0("Possible coverage ", agelabs[i]),  xlab="Year", col=plot.cols, ylim=c(0, max(rowSums(simp.table[,,i]))*1.5))
  legend("topright", legend=rev(c("NVT", "PCV7", '+PCV13', '+PCV20')),
         fill=rev(plot.cols),  cex=0.8,  box.lty=0)
  plotmax<-max(rowSums(simp.table[,,i]))
  text(x = xx, y =rowSums(simp.table[,,i]), label =pct.pcv7, pos = 3, cex = 0.5, col = '#66c2a5')
  text(x = xx, y =rowSums(simp.table[,,i]+plotmax*0.015), label =pct.pcv13, pos = 3, cex = 0.5, col = '#fc8d62')
  text(x = xx, y =rowSums(simp.table[,,i]+plotmax*0.03), label =pct.pcv20, pos = 3, cex = 0.5, col = '#8da0cb')
}