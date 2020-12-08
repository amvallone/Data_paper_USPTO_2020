setwd("/Volumes/GoogleDrive/Mi unidad/Patents 2020/UPSTO/Markov")

#Require pacakgaes

pkg <- c("magrittr", "ggplot2","directlabels","ggalt","ggrepel","dplyr","plyr","tidyverse","PAFit","igraph","rSHAPE","rshape","reshape2")
sapply(pkg,library,character.only = TRUE)


#loading the data.
lapply( paste("./base/",dir("./base")[grep("rdata",dir("./base"))],sep=""),load,.GlobalEnv) #the rdata file must be stored into a folder called "base" in the working directory.
rm(list=ls()[2:length(ls())]) #Reduce the size of the Global Enviroment to preserve the RAM


#cleansing ip_section from NAs
for(i in 1:length(data)){
  data[[i]] <- data[[i]][!is.na(data[[i]][,"ipc_section"]),]
}

#fixing typos
for(x in 1:length(data)) {
  data[[x]][,"ipc_section"]<- toupper(data[[x]][,"ipc_section"])
}

#let only the 8 ipc section

for (k in 1:length(data)){
  data[[k]] <- subset(data[[k]], subset=data[[k]][,"ipc_section"] %in% LETTERS[1:8] )
}
 
dd<-select(Reduce(rbind.fill,list(as.data.frame(t2007),as.data.frame(t2008),as.data.frame(t2009),as.data.frame(t2010),as.data.frame(t2011),as.data.frame(t2012),as.data.frame(t2013),as.data.frame(t2014),as.data.frame(t2015),as.data.frame(t2016),as.data.frame(t2017),as.data.frame(t2018),as.data.frame(t2019))),patent_id,inventor_id,patent_year,ipc_section)

dd$patent_id<-as.character(dd$patent_id)
dd$inventor_id<-as.character(dd$inventor_id)
dd$ipc_section<-as.character(dd$ipc_section)
dd$patent_year<-as.character(dd$patent_year)
dd$patent_year<-as.numeric(dd$patent_year)


pyear<-select(dd,patent_id,patent_year)
pyear<-distinct(pyear)

psection<-select(dd,patent_id,ipc_section)
psection<-distinct(psection)

pinventor<-select(dd,patent_id,inventor_id)
pinventor<-distinct(pinventor)


pinventor$inventor <- sequence(rle(pinventor$patent_id)$lengths)
pinv<-dcast(pinventor,patent_id ~ inventor,value.var = "inventor_id")

pinv<-filter(pinv,!is.na(pinv$`2`)) #delete the isolated patent


base<-list()
for(k in 2:ncol(pinv)){
  expo<-list()
  for (j in 2:ncol(pinv)){
    data <- pinv[c(1,k,j)]
    expo[[j]]<-data%>%drop_na()
  }
  base[[k]]<-expo
}

base<-base[!sapply(base,is.null)]

base2<-list()
for(n in  1:length(base)){
li<-base[[n]]
li<-compact(li)
li<-lapply(li,setNames,c("patent_id","inv_i","inv_j"))
ll<-do.call("rbind.fill",Filter(is.data.frame,li))
base2[[n]]<-ll
}

base2<-Reduce(rbind.fill,base2)
base2<-distinct(base2)
base2<-filter(base2,inv_i!=inv_j)

#modificar section
psection$ipc_section[psection$ipc_section=="e"]<-"E"
psection$ipc_section[psection$ipc_section=="g"]<-"G"
psection$ipc_section[psection$ipc_section=="h"]<-"H"
#solo se considera section A-H
pA<-psection %>% filter(ipc_section=="A")
pB<-psection %>% filter(ipc_section=="B")
pC<-psection %>% filter(ipc_section=="C")
pD<-psection %>% filter(ipc_section=="D")
pE<-psection %>% filter(ipc_section=="E")
pF<-psection %>% filter(ipc_section=="F")
pG<-psection %>% filter(ipc_section=="G")
pH<-psection %>% filter(ipc_section=="H")

pAH<-Reduce(rbind.fill,list(pA,pB,pC,pD,pE,pF,pG,pH))
pAH$num<-1
base2<-merge(base2,pAH,by.x="patent_id",by.y="patent_id")
base2<-select(base2,-num)

aux1<-select(base2,inv_j)
aux2<-select(base2,inv_i)
names(aux1)[1]<-"inv"
names(aux2)[1]<-"inv"
aux<-rbind.fill(aux1,aux2)
aux<-unique(aux)
aux<-arrange(aux,inv)
aux$num<-ave(aux$inv,FUN=seq_along)
aux$num<-as.numeric(aux$num)


names(aux)[1]<-"inv_i"
base2<-merge(base2,aux,by.x="inv_i",by.y="inv_i")
names(base2)[4]<-"ff"
names(aux)[1]<-"inv_j"
base2<-merge(base2,aux,by.x="inv_j",by.y="inv_j")
names(base2)[5]<-"tt"

base<-select(base2,patent_id,ff,tt)

base<-merge(base,pyear,by.x="patent_id",by.y="patent_id")

#########################

#network construction

b1<-base %>% filter(patent_year<=2007)
z <- data.frame(b1$ff, b1$tt)
n1 <- graph.data.frame(z, directed = F)
n1<-simplify(n1)
n1b<-decompose(n1,min.vertices = 65000)
n1b<-n1b[[1]]


b2<-base %>% filter(patent_year<=2008)
z <- data.frame(b2$ff, b2$tt)
n2 <- graph.data.frame(z, directed = F)
n2<-simplify(n2)
n2b<-decompose(n2,min.vertices = 65062)
n2b<-n2b[[1]]


b3<-base %>% filter(patent_year<=2009)
z <- data.frame(b3$ff, b3$tt)
n3 <- graph.data.frame(z, directed = F)
n3<-simplify(n3)
n3b<-decompose(n3,min.vertices = 65062)
n3b<-n3b[[1]]


b4<-base %>% filter(patent_year<=2010)
z <- data.frame(b4$ff, b4$tt)
n4 <- graph.data.frame(z, directed = F)
n4<-simplify(n4)
n4b<-decompose(n4,min.vertices = 65062)
n4b<-n4b[[1]]


b5<-base %>% filter(patent_year<=2011)
z <- data.frame(b5$ff, b5$tt)
n5 <- graph.data.frame(z, directed = F)
n5<-simplify(n5)
n5b<-decompose(n5,min.vertices = 65062)
n5b<-n5b[[1]]

b6<-base %>% filter(patent_year<=2012)
z <- data.frame(b6$ff, b6$tt)
n6 <- graph.data.frame(z, directed = F)
n6<-simplify(n6)
n6b<-decompose(n6,min.vertices = 65062)
n6b<-n6b[[1]]


b7<-base %>% filter(patent_year<=2013)
z <- data.frame(b7$ff, b7$tt)
n7 <- graph.data.frame(z, directed = F)
n7<-simplify(n7)
n7b<-decompose(n7,min.vertices = 65062)
n7b<-n7b[[1]]


b8<-base %>% filter(patent_year<=2014)
z <- data.frame(b8$ff, b8$tt)
n8 <- graph.data.frame(z, directed = F)
n8<-simplify(n8)
n8b<-decompose(n8,min.vertices = 65062)
n8b<-n8b[[1]]


b9<-base %>% filter(patent_year<=2015)
z <- data.frame(b9$ff, b9$tt)
n9 <- graph.data.frame(z, directed = F)
n9<-simplify(n9)
n9b<-decompose(n9,min.vertices = 65062)
n9b<-n9b[[1]]


b10<-base %>% filter(patent_year<=2016)
z <- data.frame(b10$ff, b10$tt)
n10 <- graph.data.frame(z, directed = F)
n10<-simplify(n10)
n10b<-decompose(n10,min.vertices = 65062)
n10b<-n10b[[1]]


b11<-base %>% filter(patent_year<=2017)
z <- data.frame(b11$ff, b11$tt)
n11 <- graph.data.frame(z, directed = F)
n11<-simplify(n11)
n11b<-decompose(n11,min.vertices = 65062)
n11b<-n11b[[1]]


b12<-base %>% filter(patent_year<=2018)
z <- data.frame(b12$ff, b12$tt)
n12 <- graph.data.frame(z, directed = F)
n12<-simplify(n12)
n12b<-decompose(n12,min.vertices = 65062)
n12b<-n12b[[1]]


b13<-base %>% filter(patent_year<=2019)
z <- data.frame(b13$ff, b13$tt)
n13 <- graph.data.frame(z, directed = F)
n13<-simplify(n13)
n13b<-decompose(n13,min.vertices = 65062)
n13b<-n13b[[1]]

##PL
data<-degree(n13b)
data<-data[data>0]
m_pl<-dispois$new(data)
est_pl<-estimate_xmin(m_pl)
est_pl$xmin
est_pl$pars
est_pl$gof
m_pl$setXmin(est_pl)
bs_pl<-bootstrap_p(m_pl,no_of_sims = 1000,threads = 8,seed=123)



# diameter 
diam1<-diameter(n1b)
diam2<-diameter(n2b)
diam3<-diameter(n3b)
diam4<-diameter(n4b)
diam5<-diameter(n5b)
diam6<-diameter(n6b)
diam7<-diameter(n7b)
diam8<-diameter(n8b)
diam9<-diameter(n9b)
diam10<-diameter(n10b)
diam11<-diameter(n11b)
diam12<-diameter(n12b)
diam13<-diameter(n13b)


sp1<-shortest.paths(n1b)
sp2<-shortest.paths(n2b)
sp3<-shortest.paths(n3b)
sp4<-shortest.paths(n4b)
sp5<-shortest.paths(n5b)
sp6<-shortest.paths(n6b)
sp7<-shortest.paths(n7b)
sp8<-shortest.paths(n8b)
sp9<-shortest.paths(n9b)
sp10<-shortest.paths(n10b)
sp11<-shortest.paths(n11b)
sp12<-shortest.paths(n12b)
sp13<-shortest.paths(n13b)


### pfit most importat component
big <- decompose(net13, min.vertices=1485782)
big<-big[[1]]
bb<-as_data_frame(big)
bb<-PAFit::to_igraph(Pnet)

bb <- decompose(bb, min.vertices=1485782)
bb<-bb[[1]]
Pbig<-as.PAFit_net(bb,type="undirected")
View(bb)
com<-components(Pnet)
ba<-as_data_frame(bb)
View(ba)
z<-data.frame(ba$from,ba$to,ba$time)
names(z)[1]<-"ff_inv"
names(z)[2]<-"to_inv"
names(z)[3]<-"t_pat"
Pbig<-as.PAFit_net(z,type="undirected")
sPbig<-get_statistics(Pbig)
z$ff_inv<-as.numeric(z$ff_inv)
Z$to_inv<-as.numeric(z$to_inv)
Pbig<-as.PAFit_net(z,type="undirected")
sPbig<-get_statistics(Pbig)
result_Pbig<-only_A_estimate(Pbig,sPbig)
plot(result_Pbig,sPbig,line="TRUE",cex=2,cex.axis=2,cex.lab=2)
summary(sPbig)


## agregar diam, gp y ne.

##diam
#diam13<-diameter(n13b) para obtener tabla año-diametro, realizar para acumulado de cada año y section, GCC

##gp
#gp<-components(n13b) realizar para el acumulado de cada año y section 3 principales components y dividir nodes sobre el total,  GCC

#ne node vs edge en el acumulado de cada año.


########################


n.ipc <- function(x){
  p <- unique(data[[x]][,1:3])
  out <- table(p[,2],p[,3])
  out
}

ipc <- lapply(1:length(data),function(x) n.ipc(x))

sec <- lapply(1:length(ipc),function(x)tcrossprod(t(ipc[[x]])))
sec1 <- lapply(1:length(sec), function(x) sec[[x]]/rowSums(sec[[x]]))
secT <- Reduce('+',sec)
secT1 <- secT/rowSums(secT)


p.sec <- sapply(1:length(sec1),function(x) 1-diag(sec1[[x]]))
colnames(p.sec) <- 2007:2019
prais <- reshape2::melt(p.sec)
p.sec1<-p.sec/p.sec[,1]
prais1 <- reshape2::melt(p.sec1)
colnames(prais1) <- c("Section","Year","Index")
colnames(prais) <- c("Section","Year","Index")

diam1 <- diam
diam1[,2:ncol(diam1)] <- diam[,2:ncol(diam1)] /diam[,2]
tca.diam <- (diam[2:ncol(diam)]-diam[,2])/matrix(0:12,8,13,byrow=T)
tca.diam[,1]<- LETTERS[1:8]
colnames(tca.diam)[1] <- "X1"
tcr.diam <- (log(diam[2:ncol(diam)])-log(diam[,2]))/matrix(0:12,8,13,byrow=T)
tcr.diam[,1]<- LETTERS[1:8]
colnames(tcr.diam)[1] <- "X1"


diam <- reshape2::melt(diam)
diam1 <- reshape2::melt(diam1)
diam[,2] %<>% as.character %>% as.numeric
diam1[,2] %<>% as.character %>% as.numeric
colnames(diam1) <- c("Section","Year","Index")
colnames(diam) <- c("Section","Year","Index")

tcr.div <- (log(p.sec)-log(p.sec[,1]))/matrix(0:12,8,13,byrow=T)
tcr.div <- tcr.div[,2:13]

tcr.div1 <- reshape2::melt(tcr.div)
colnames(tcr.div1) <- c("Section","Year","Index")

tcr.diam1 <- reshape2::melt(tcr.diam)
tcr.diam1[,2] %<>% as.character %>% as.numeric
colnames(tcr.diam1) <- c("Section","Year","Index")


#Figure 4a
ggplot(diam,aes(x=Year,y=Index,color=Section))+geom_xspline(aes(linetype=Section),size=.9)+geom_point(aes(shape=Section),size=2)+theme_bw()+theme(axis.title=element_text(size=16,face="bold"),axis.text = element_text(size=14),legend.text =element_text(size=14))+scale_x_continuous("Year", 2007:2019, 2007:2019)+labs(y="Diameter")+scale_colour_grey(name="IPC Section")+ scale_shape_manual("IPC Section",values=c(15, 16, 17,18,19,3,4,8))+geom_dl(aes(label = Section), method = list("last.bumpup",cex=1.1,hjust=-1,vjust=-.1), color="black")+scale_linetype_manual("IPC Section",values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "F1", "12345678"))+ scale_y_continuous(breaks=seq(0,80,5))

# Figure 4b
ggplot(tcr.diam1,aes(x=Year,y=Index,color=Section))+geom_xspline(aes(linetype=Section),size=.9)+geom_point(aes(shape=Section),size=2)+theme_bw()+theme(axis.title=element_text(size=16,face="bold"),axis.text = element_text(size=14),legend.text =element_text(size=14))+scale_x_continuous("Year", 2008:2019, 2008:2019)+labs(y="Diameter rate of relative growth")+scale_colour_grey(name="IPC Section")+ scale_shape_manual("IPC Section",values=c(15, 16, 17,18,19,3,4,8))+geom_dl(aes(label = Section), method = list("last.bumpup",cex=1.1,hjust=-1,vjust=-.1), color="black")+scale_linetype_manual("IPC Section",values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "F1", "12345678"))#+ scale_y_continuous(breaks=seq(0,80,5))

# Figure 7a
ggplot(prais,aes(x=Year,y=Index,color=Section))+geom_xspline(aes(linetype=Section),size=.9)+geom_point(aes(shape=Section),size=2)+theme_bw()+theme(axis.title=element_text(size=16,face="bold"),axis.text = element_text(size=14),legend.text =element_text(size=14))+scale_x_continuous("Year", 2007:2019, 2007:2019)+labs(y="Diversificaction (Pr)")+ scale_y_continuous(breaks=seq(0,1,.05))+scale_colour_grey(name="IPC Section")+ scale_shape_manual("IPC Section",values=c(15, 16, 17,18,19,3,4,8))+geom_dl(aes(label = Section), method = list("last.bumpup",cex=1.1,hjust=-1,vjust=-.1), color="black")+scale_linetype_manual("IPC Section",values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "F1", "12345678"))

# Figure 7b
ggplot(tcr.div1,aes(x=Year,y=Index,color=Section))+geom_xspline(aes(linetype=Section),size=.9)+geom_point(aes(shape=Section),size=2)+theme_bw()+theme(axis.title=element_text(size=16,face="bold"),axis.text = element_text(size=14),legend.text =element_text(size=14))+scale_x_continuous("Year", 2008:2019, 2008:2019)+labs(y="Diversificaction rate of relative growth")+scale_colour_grey(name="IPC Section")+ scale_shape_manual("IPC Section",values=c(15, 16, 17,18,19,3,4,8))+geom_dl(aes(label = Section), method = list("last.bumpup",cex=1.1,hjust=-1,vjust=-.1), color="black")+scale_linetype_manual("IPC Section",values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "F1", "12345678"))#+ scale_y_continuous(breaks=seq(0,80,5))

#Figure 8

data2 <- diam1[,1:3]
colnames(data2)[3] <- "Diameter"
data2["Diversification"] <- prais1[,3]
ggplot(data2, aes(x=Diversification, y=Diameter,label=Year))+ geom_point(size=2)+geom_text_repel(aes(label=Year),hjust=0, vjust=0,force=1,color="black")+geom_hline(yintercept=1,linetype="dashed")+geom_vline(xintercept=1,linetype="dashed")+theme_bw()+facet_wrap(.~Sector)+xlab("Diversification Index (base 2007=1.0)")+ylab("Diameter Index (base 2007=1.0)")+scale_y_continuous(breaks=seq(0,7,1))+scale_x_continuous(breaks=seq(0,2.5,.25))



# Figure

ggplot(neL, aes(x=node, y=edge)) + 
  geom_line(size=1.5)+
  scale_y_continuous("Number of edges", trans = "log10")+
  scale_x_continuous("Number of nodes", trans = "log10")+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major=element_blank(),panel.grid.minor=element_blank())

# Figure

ggplot(gR, aes(x=year, y=size, col=gp,group=gp)) + 
  geom_line(aes(linetype=gp,colour=gp))+
  scale_y_continuous("CC Size", trans = "log10")+
  scale_x_continuous("Year",breaks = seq(from = 2007, to = 2019, by = 1))+
  scale_linetype_manual(values=c("solid", "solid", "solid", "dotdash"))+
  theme_bw()+scale_color_manual(values = c("gray45", "gray75", "gray25","gray25"))+
  theme(legend.position = "none",panel.grid.major=element_blank(),panel.grid.minor=element_blank())



# Figure 

p<-ggplot(gR, aes(x=year, y=size, col=gp,group=gp)) + 
  geom_line(aes(linetype=gp,colour=gp))+
  scale_y_continuous("CC Size", trans = "log10")+
  scale_x_continuous("Year",breaks = seq(from = 2007, to = 2019, by = 1))+
  scale_linetype_manual(values=c("solid", "dotdash", "solid", "solid"))+
  theme_bw()+scale_color_manual(values = c("gray25", "gray25", "gray45","gray75"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

p + facet_wrap( ~ class, nrow = 3) + 
  theme(legend.title= element_blank(),legend.position = "bottom")










