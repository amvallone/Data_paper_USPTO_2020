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

# preparing the data for sna anañysis
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


#section A-H
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

##Total
b1<-base %>% filter(patent_year<=2007)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n1b<-decompose(n1,min.vertices =gcc)
n1b<-n1b[[1]]
diam1<-diameter(n1b)
v1<-as.numeric(length(V(n1b)))
e1<-gsize(n1b)
gcc1<-c1
names(gcc1)[1]<-"size"
gcc1$year<-2007
gcc1$gp<-c("GCC","CC2","CC3")

b1<-base %>% filter(patent_year<=2008)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n2b<-decompose(n1,min.vertices =gcc)
n2b<-n2b[[1]]
diam2<-diameter(n2b)
v2<-as.numeric(length(V(n2b)))
e2<-gsize(n2b)
gcc2<-c1
names(gcc2)[1]<-"size"
gcc2$year<-2008
gcc2$gp<-c("GCC","CC2","CC3")

b1<-base %>% filter(patent_year<=2009)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n3b<-decompose(n1,min.vertices =gcc)
n3b<-n3b[[1]]
diam3<-diameter(n3b)
v3<-as.numeric(length(V(n3b)))
e3<-gsize(n3b)
gcc3<-c1
names(gcc3)[1]<-"size"
gcc3$year<-2009
gcc3$gp<-c("GCC","CC2","CC3")


b1<-base %>% filter(patent_year<=2010)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n4b<-decompose(n1,min.vertices =gcc)
n4b<-n4b[[1]]
diam4<-diameter(n4b)
v4<-as.numeric(length(V(n4b)))
e4<-gsize(n4b)
gcc4<-c1
names(gcc4)[1]<-"size"
gcc4$year<-2010
gcc4$gp<-c("GCC","CC2","CC3")


b1<-base %>% filter(patent_year<=2011)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n5b<-decompose(n1,min.vertices =gcc)
n5b<-n5b[[1]]
diam5<-diameter(n5b)
v5<-as.numeric(length(V(n5b)))
e5<-gsize(n5b)
gcc5<-c1
names(gcc5)[1]<-"size"
gcc5$year<-2011
gcc5$gp<-c("GCC","CC2","CC3")


b1<-base %>% filter(patent_year<=2012)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n6b<-decompose(n1,min.vertices =gcc)
n6b<-n6b[[1]]
diam6<-diameter(n6b)
v6<-as.numeric(length(V(n6b)))
e6<-gsize(n6b)
gcc6<-c1
names(gcc6)[1]<-"size"
gcc6$year<-2012
gcc6$gp<-c("GCC","CC2","CC3")


b1<-base %>% filter(patent_year<=2013)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n7b<-decompose(n1,min.vertices =gcc)
n7b<-n7b[[1]]
diam7<-diameter(n7b)
v7<-as.numeric(length(V(n7b)))
e7<-gsize(n7b)
gcc7<-c1
names(gcc7)[1]<-"size"
gcc7$year<-2013
gcc7$gp<-c("GCC","CC2","CC3")


b1<-base %>% filter(patent_year<=2014)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n8b<-decompose(n1,min.vertices =gcc)
n8b<-n8b[[1]]
diam8<-diameter(n8b)
v8<-as.numeric(length(V(n8b)))
e8<-gsize(n8b)
gcc8<-c1
names(gcc8)[1]<-"size"
gcc8$year<-2014
gcc8$gp<-c("GCC","CC2","CC3")

b1<-base %>% filter(patent_year<=2015)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n9b<-decompose(n1,min.vertices =gcc)
n9b<-n9b[[1]]
diam9<-diameter(n9b)
v9<-as.numeric(length(V(n9b)))
e9<-gsize(n9b)
gcc9<-c1
names(gcc9)[1]<-"size"
gcc9$year<-2015
gcc9$gp<-c("GCC","CC2","CC3")


b1<-base %>% filter(patent_year<=2016)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n10b<-decompose(n1,min.vertices =gcc)
n10b<-n10b[[1]]
diam10<-diameter(n10b)
v10<-as.numeric(length(V(n10b)))
e10<-gsize(n10b)
gcc10<-c1
names(gcc10)[1]<-"size"
gcc10$year<-2016
gcc10$gp<-c("GCC","CC2","CC3")

b1<-base %>% filter(patent_year<=2017)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n11b<-decompose(n1,min.vertices =gcc)
n11b<-n11b[[1]]
diam11<-diameter(n11b)
v11<-as.numeric(length(V(n11b)))
e11<-gsize(n11b)
gcc11<-c1
names(gcc11)[1]<-"size"
gcc11$year<-2017
gcc11$gp<-c("GCC","CC2","CC3")


b1<-base %>% filter(patent_year<=2018)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n12b<-decompose(n1,min.vertices =gcc)
n12b<-n12b[[1]]
diam12<-diameter(n12b)
v12<-as.numeric(length(V(n12b)))
e12<-gsize(n12b)
gcc12<-c1
names(gcc12)[1]<-"size"
gcc12$year<-2018
gcc12$gp<-c("GCC","CC2","CC3")

b1<-base %>% filter(patent_year<=2019)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n13b<-decompose(n1,min.vertices =gcc)
n13b<-n13b[[1]]
diam13<-diameter(n13b)
v13<-as.numeric(length(V(n13b)))
e13<-gsize(n13b)
gcc13<-c1
names(gcc13)[1]<-"size"
gcc13$year<-2019
gcc13$gp<-c("GCC","CC2","CC3")

year<-c(2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)
L<-c(diam1,diam2,diam3,diam4,diam5,diam6,diam7,diam8,diam9,diam10,diam11,diam12,diam13)
node<-c(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13)
edge<-c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13)
diamL<-data.frame(year,L)
neL<-data.frame(year,node,edge)

gpL<-Reduce(rbind.fill,list(gcc1,gcc2,gcc3,gcc4,gcc5,gcc6,gcc7,gcc8,gcc9,gcc10,gcc11,gcc12,gcc13))

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

### pfit most important component Total
big <- decompose(n13b, min.vertices=gcc)
big<-big[[1]]
bb<-as_data_frame(big)
bb<-PAFit::to_igraph(Pnet)

bb <- decompose(bb, min.vertices=gcc)
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

deg1<-table(degree(n1b))
deg1<-as.data.frame(deg1)
deg1$cat<-2007
deg2<-table(degree(n7b))
deg2<-as.data.frame(deg2)
deg2$cat<-2013
deg3<-table(degree(n13b))
deg3<-as.data.frame(deg3)
deg3$cat<-2019
degL<-Reduce(rbind.fill,list(deg1,deg2,deg3))
degL$cat<-as.factor(degL$cat)
names(degL)[1]<-"degree"

b1<-base %>% filter(patent_year>=2010 & patent_year<=2010)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc1<-c1
gcc1$year<-2010
gcc1$gp<-"GCC Post 2010 subgraph, no past"


b1<-base %>% filter(patent_year>=2010 & patent_year<=2011)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc2<-c1
gcc2$year<-2011
gcc2$gp<-"GCC Post 2010 subgraph, no past"


b1<-base %>% filter(patent_year>=2010 & patent_year<=2012)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc3<-c1
gcc3$year<-2012
gcc3$gp<-"GCC Post 2010 subgraph, no past"

b1<-base %>% filter(patent_year>=2010 & patent_year<=2013)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc4<-c1
gcc4$year<-2013
gcc4$gp<-"GCC Post 2010 subgraph, no past"

b1<-base %>% filter(patent_year>=2010 & patent_year<=2014)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc5<-c1
gcc5$year<-2014
gcc5$gp<-"GCC Post 2010 subgraph, no past"

b1<-base %>% filter(patent_year>=2010 & patent_year<=2015)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc6<-c1
gcc6$year<-2015
gcc6$gp<-"GCC Post 2010 subgraph, no past"

b1<-base %>% filter(patent_year>=2010 & patent_year<=2016)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc7<-c1
gcc7$year<-2016
gcc7$gp<-"GCC Post 2010 subgraph, no past"


b1<-base %>% filter(patent_year>=2010 & patent_year<=2017)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc8<-c1
gcc8$year<-2017
gcc8$gp<-"GCC Post 2010 subgraph, no past"


b1<-base %>% filter(patent_year>=2010 & patent_year<=2018)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc9<-c1
gcc9$year<-2018
gcc9$gp<-"GCC Post 2010 subgraph, no past"

b1<-base %>% filter(patent_year>=2010 & patent_year<=2019)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc10<-c1
gcc10$year<-2019
gcc10$gp<-"GCC Post 2010 subgraph, no past"

gpL<-Reduce(rbind.fill,list(gpL,gcc1,gcc2,gcc3,gcc4,gcc5,gcc6,gcc7,gcc8,gcc9,gcc10))


###filter for section
baseA<-base %>% filter(ipc_section=="A")
baseB<-base %>% filter(ipc_section=="B")
baseC<-base %>% filter(ipc_section=="C")
baseD<-base %>% filter(ipc_section=="D")
baseE<-base %>% filter(ipc_section=="E")
baseF<-base %>% filter(ipc_section=="F")
baseG<-base %>% filter(ipc_section=="G")
baseH<-base %>% filter(ipc_section=="H")

##Section A
b1<-baseA %>% filter(patent_year<=2007)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n1b<-decompose(n1,min.vertices =gcc)
n1b<-n1b[[1]]
diam1<-diameter(n1b)
v1<-as.numeric(length(V(n1b)))
e1<-gsize(n1b)
gcc1<-c1
names(gcc1)[1]<-"size"
gcc1$year<-2007
gcc1$gp<-c("GCC","CC2","CC3")

b1<-baseA %>% filter(patent_year<=2008)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n2b<-decompose(n1,min.vertices =gcc)
n2b<-n2b[[1]]
diam2<-diameter(n2b)
v2<-as.numeric(length(V(n2b)))
e2<-gsize(n2b)
gcc2<-c1
names(gcc2)[1]<-"size"
gcc2$year<-2008
gcc2$gp<-c("GCC","CC2","CC3")

b1<-baseA %>% filter(patent_year<=2009)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n3b<-decompose(n1,min.vertices =gcc)
n3b<-n3b[[1]]
diam3<-diameter(n3b)
v3<-as.numeric(length(V(n3b)))
e3<-gsize(n3b)
gcc3<-c1
names(gcc3)[1]<-"size"
gcc3$year<-2009
gcc3$gp<-c("GCC","CC2","CC3")


b1<-baseA %>% filter(patent_year<=2010)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n4b<-decompose(n1,min.vertices =gcc)
n4b<-n4b[[1]]
diam4<-diameter(n4b)
v4<-as.numeric(length(V(n4b)))
e4<-gsize(n4b)
gcc4<-c1
names(gcc4)[1]<-"size"
gcc4$year<-2010
gcc4$gp<-c("GCC","CC2","CC3")


b1<-baseA %>% filter(patent_year<=2011)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n5b<-decompose(n1,min.vertices =gcc)
n5b<-n5b[[1]]
diam5<-diameter(n5b)
v5<-as.numeric(length(V(n5b)))
e5<-gsize(n5b)
gcc5<-c1
names(gcc5)[1]<-"size"
gcc5$year<-2011
gcc5$gp<-c("GCC","CC2","CC3")


b1<-baseA %>% filter(patent_year<=2012)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n6b<-decompose(n1,min.vertices =gcc)
n6b<-n6b[[1]]
diam6<-diameter(n6b)
v6<-as.numeric(length(V(n6b)))
e6<-gsize(n6b)
gcc6<-c1
names(gcc6)[1]<-"size"
gcc6$year<-2012
gcc6$gp<-c("GCC","CC2","CC3")


b1<-baseA %>% filter(patent_year<=2013)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n7b<-decompose(n1,min.vertices =gcc)
n7b<-n7b[[1]]
diam7<-diameter(n7b)
v7<-as.numeric(length(V(n7b)))
e7<-gsize(n7b)
gcc7<-c1
names(gcc7)[1]<-"size"
gcc7$year<-2013
gcc7$gp<-c("GCC","CC2","CC3")


b1<-baseA %>% filter(patent_year<=2014)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n8b<-decompose(n1,min.vertices =gcc)
n8b<-n8b[[1]]
diam8<-diameter(n8b)
v8<-as.numeric(length(V(n8b)))
e8<-gsize(n8b)
gcc8<-c1
names(gcc8)[1]<-"size"
gcc8$year<-2014
gcc8$gp<-c("GCC","CC2","CC3")

b1<-baseA %>% filter(patent_year<=2015)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n9b<-decompose(n1,min.vertices =gcc)
n9b<-n9b[[1]]
diam9<-diameter(n9b)
v9<-as.numeric(length(V(n9b)))
e9<-gsize(n9b)
gcc9<-c1
names(gcc9)[1]<-"size"
gcc9$year<-2015
gcc9$gp<-c("GCC","CC2","CC3")


b1<-baseA %>% filter(patent_year<=2016)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n10b<-decompose(n1,min.vertices =gcc)
n10b<-n10b[[1]]
diam10<-diameter(n10b)
v10<-as.numeric(length(V(n10b)))
e10<-gsize(n10b)
gcc10<-c1
names(gcc10)[1]<-"size"
gcc10$year<-2016
gcc10$gp<-c("GCC","CC2","CC3")

b1<-baseA %>% filter(patent_year<=2017)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n11b<-decompose(n1,min.vertices =gcc)
n11b<-n11b[[1]]
diam11<-diameter(n11b)
v11<-as.numeric(length(V(n11b)))
e11<-gsize(n11b)
gcc11<-c1
names(gcc11)[1]<-"size"
gcc11$year<-2017
gcc11$gp<-c("GCC","CC2","CC3")


b1<-baseA %>% filter(patent_year<=2018)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n12b<-decompose(n1,min.vertices =gcc)
n12b<-n12b[[1]]
diam12<-diameter(n12b)
v12<-as.numeric(length(V(n12b)))
e12<-gsize(n12b)
gcc12<-c1
names(gcc12)[1]<-"size"
gcc12$year<-2018
gcc12$gp<-c("GCC","CC2","CC3")

b1<-baseA %>% filter(patent_year<=2019)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n13b<-decompose(n1,min.vertices =gcc)
n13b<-n13b[[1]]
diam13<-diameter(n13b)
v13<-as.numeric(length(V(n13b)))
e13<-gsize(n13b)
gcc13<-c1
names(gcc13)[1]<-"size"
gcc13$year<-2019
gcc13$gp<-c("GCC","CC2","CC3")


year<-c(2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)
A<-c(diam1,diam2,diam3,diam4,diam5,diam6,diam7,diam8,diam9,diam10,diam11,diam12,diam13)
node<-c(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13)
edge<-c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13)
diamA<-data.frame(year,A)
neA<-data.frame(year,node,edge)

gpA<-Reduce(rbind.fill,list(gcc1,gcc2,gcc3,gcc4,gcc5,gcc6,gcc7,gcc8,gcc9,gcc10,gcc11,gcc12,gcc13))

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

### pfit most important component A
big <- decompose(n13b, min.vertices=gcc)
big<-big[[1]]
bb<-as_data_frame(big)
bb<-PAFit::to_igraph(Pnet)

bb <- decompose(bb, min.vertices=gcc)
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

deg1<-table(degree(n1b))
deg1<-as.data.frame(deg1)
deg1$cat<-2007
deg2<-table(degree(n7b))
deg2<-as.data.frame(deg2)
deg2$cat<-2013
deg3<-table(degree(n13b))
deg3<-as.data.frame(deg3)
deg3$cat<-2019
degA<-Reduce(rbind.fill,list(deg1,deg2,deg3))
degA$cat<-as.factor(degA$cat)
names(degA)[1]<-"degree"

b1<-baseA %>% filter(patent_year>=2010 & patent_year<=2010)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc1<-c1
gcc1$year<-2010
gcc1$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseA %>% filter(patent_year>=2010 & patent_year<=2011)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc2<-c1
gcc2$year<-2011
gcc2$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseA %>% filter(patent_year>=2010 & patent_year<=2012)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc3<-c1
gcc3$year<-2012
gcc3$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseA %>% filter(patent_year>=2010 & patent_year<=2013)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc4<-c1
gcc4$year<-2013
gcc4$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseA %>% filter(patent_year>=2010 & patent_year<=2014)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc5<-c1
gcc5$year<-2014
gcc5$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseA %>% filter(patent_year>=2010 & patent_year<=2015)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc6<-c1
gcc6$year<-2015
gcc6$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseA %>% filter(patent_year>=2010 & patent_year<=2016)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc7<-c1
gcc7$year<-2016
gcc7$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseA %>% filter(patent_year>=2010 & patent_year<=2017)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc8<-c1
gcc8$year<-2017
gcc8$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseA %>% filter(patent_year>=2010 & patent_year<=2018)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc9<-c1
gcc9$year<-2018
gcc9$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseA %>% filter(patent_year>=2010 & patent_year<=2019)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc10<-c1
gcc10$year<-2019
gcc10$gp<-"GCC Post 2010 subgraph, no past"


gpA<-Reduce(rbind.fill,list(gpA,gcc1,gcc2,gcc3,gcc4,gcc5,gcc6,gcc7,gcc8,gcc9,gcc10))

##Section B
b1<-baseB %>% filter(patent_year<=2007)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n1b<-decompose(n1,min.vertices =gcc)
n1b<-n1b[[1]]
diam1<-diameter(n1b)
v1<-as.numeric(length(V(n1b)))
e1<-gsize(n1b)
gcc1<-c1
names(gcc1)[1]<-"size"
gcc1$year<-2007
gcc1$gp<-c("GCC","CC2","CC3")

b1<-baseB %>% filter(patent_year<=2008)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n2b<-decompose(n1,min.vertices =gcc)
n2b<-n2b[[1]]
diam2<-diameter(n2b)
v2<-as.numeric(length(V(n2b)))
e2<-gsize(n2b)
gcc2<-c1
names(gcc2)[1]<-"size"
gcc2$year<-2008
gcc2$gp<-c("GCC","CC2","CC3")

b1<-baseB %>% filter(patent_year<=2009)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n3b<-decompose(n1,min.vertices =gcc)
n3b<-n3b[[1]]
diam3<-diameter(n3b)
v3<-as.numeric(length(V(n3b)))
e3<-gsize(n3b)
gcc3<-c1
names(gcc3)[1]<-"size"
gcc3$year<-2009
gcc3$gp<-c("GCC","CC2","CC3")


b1<-baseB %>% filter(patent_year<=2010)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n4b<-decompose(n1,min.vertices =gcc)
n4b<-n4b[[1]]
diam4<-diameter(n4b)
v4<-as.numeric(length(V(n4b)))
e4<-gsize(n4b)
gcc4<-c1
names(gcc4)[1]<-"size"
gcc4$year<-2010
gcc4$gp<-c("GCC","CC2","CC3")


b1<-baseB %>% filter(patent_year<=2011)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n5b<-decompose(n1,min.vertices =gcc)
n5b<-n5b[[1]]
diam5<-diameter(n5b)
v5<-as.numeric(length(V(n5b)))
e5<-gsize(n5b)
gcc5<-c1
names(gcc5)[1]<-"size"
gcc5$year<-2011
gcc5$gp<-c("GCC","CC2","CC3")


b1<-baseB %>% filter(patent_year<=2012)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n6b<-decompose(n1,min.vertices =gcc)
n6b<-n6b[[1]]
diam6<-diameter(n6b)
v6<-as.numeric(length(V(n6b)))
e6<-gsize(n6b)
gcc6<-c1
names(gcc6)[1]<-"size"
gcc6$year<-2012
gcc6$gp<-c("GCC","CC2","CC3")


b1<-baseB %>% filter(patent_year<=2013)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n7b<-decompose(n1,min.vertices =gcc)
n7b<-n7b[[1]]
diam7<-diameter(n7b)
v7<-as.numeric(length(V(n7b)))
e7<-gsize(n7b)
gcc7<-c1
names(gcc7)[1]<-"size"
gcc7$year<-2013
gcc7$gp<-c("GCC","CC2","CC3")


b1<-baseB %>% filter(patent_year<=2014)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n8b<-decompose(n1,min.vertices =gcc)
n8b<-n8b[[1]]
diam8<-diameter(n8b)
v8<-as.numeric(length(V(n8b)))
e8<-gsize(n8b)
gcc8<-c1
names(gcc8)[1]<-"size"
gcc8$year<-2014
gcc8$gp<-c("GCC","CC2","CC3")

b1<-baseB %>% filter(patent_year<=2015)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n9b<-decompose(n1,min.vertices =gcc)
n9b<-n9b[[1]]
diam9<-diameter(n9b)
v9<-as.numeric(length(V(n9b)))
e9<-gsize(n9b)
gcc9<-c1
names(gcc9)[1]<-"size"
gcc9$year<-2015
gcc9$gp<-c("GCC","CC2","CC3")


b1<-baseB %>% filter(patent_year<=2016)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n10b<-decompose(n1,min.vertices =gcc)
n10b<-n10b[[1]]
diam10<-diameter(n10b)
v10<-as.numeric(length(V(n10b)))
e10<-gsize(n10b)
gcc10<-c1
names(gcc10)[1]<-"size"
gcc10$year<-2016
gcc10$gp<-c("GCC","CC2","CC3")

b1<-baseB %>% filter(patent_year<=2017)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n11b<-decompose(n1,min.vertices =gcc)
n11b<-n11b[[1]]
diam11<-diameter(n11b)
v11<-as.numeric(length(V(n11b)))
e11<-gsize(n11b)
gcc11<-c1
names(gcc11)[1]<-"size"
gcc11$year<-2017
gcc11$gp<-c("GCC","CC2","CC3")


b1<-baseB %>% filter(patent_year<=2018)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n12b<-decompose(n1,min.vertices =gcc)
n12b<-n12b[[1]]
diam12<-diameter(n12b)
v12<-as.numeric(length(V(n12b)))
e12<-gsize(n12b)
gcc12<-c1
names(gcc12)[1]<-"size"
gcc12$year<-2018
gcc12$gp<-c("GCC","CC2","CC3")

b1<-baseB %>% filter(patent_year<=2019)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n13b<-decompose(n1,min.vertices =gcc)
n13b<-n13b[[1]]
diam13<-diameter(n13b)
v13<-as.numeric(length(V(n13b)))
e13<-gsize(n13b)
gcc13<-c1
names(gcc13)[1]<-"size"
gcc13$year<-2019
gcc13$gp<-c("GCC","CC2","CC3")


year<-c(2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)
B<-c(diam1,diam2,diam3,diam4,diam5,diam6,diam7,diam8,diam9,diam10,diam11,diam12,diam13)
node<-c(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13)
edge<-c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13)
diamB<-data.frame(year,B)
neB<-data.frame(year,node,edge)

gpB<-Reduce(rbind.fill,list(gcc1,gcc2,gcc3,gcc4,gcc5,gcc6,gcc7,gcc8,gcc9,gcc10,gcc11,gcc12,gcc13))

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

### pfit most important component B
big <- decompose(n13b, min.vertices=gcc)
big<-big[[1]]
bb<-as_data_frame(big)
bb<-PAFit::to_igraph(Pnet)

bb <- decompose(bb, min.vertices=gcc)
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

deg1<-table(degree(n1b))
deg1<-as.data.frame(deg1)
deg1$cat<-2007
deg2<-table(degree(n7b))
deg2<-as.data.frame(deg2)
deg2$cat<-2013
deg3<-table(degree(n13b))
deg3<-as.data.frame(deg3)
deg3$cat<-2019
degB<-Reduce(rbind.fill,list(deg1,deg2,deg3))
degB$cat<-as.factor(degB$cat)
names(degB)[1]<-"degree"

b1<-baseB %>% filter(patent_year>=2010 & patent_year<=2010)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc1<-c1
gcc1$year<-2010
gcc1$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseB %>% filter(patent_year>=2010 & patent_year<=2011)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc2<-c1
gcc2$year<-2011
gcc2$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseB %>% filter(patent_year>=2010 & patent_year<=2012)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc3<-c1
gcc3$year<-2012
gcc3$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseB %>% filter(patent_year>=2010 & patent_year<=2013)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc4<-c1
gcc4$year<-2013
gcc4$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseB %>% filter(patent_year>=2010 & patent_year<=2014)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc5<-c1
gcc5$year<-2014
gcc5$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseB %>% filter(patent_year>=2010 & patent_year<=2015)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc6<-c1
gcc6$year<-2015
gcc6$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseB %>% filter(patent_year>=2010 & patent_year<=2016)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc7<-c1
gcc7$year<-2016
gcc7$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseB %>% filter(patent_year>=2010 & patent_year<=2017)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc8<-c1
gcc8$year<-2017
gcc8$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseB %>% filter(patent_year>=2010 & patent_year<=2018)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc9<-c1
gcc9$year<-2018
gcc9$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseB %>% filter(patent_year>=2010 & patent_year<=2019)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc10<-c1
gcc10$year<-2019
gcc10$gp<-"GCC Post 2010 subgraph, no past"


gpB<-Reduce(rbind.fill,list(gpB,gcc1,gcc2,gcc3,gcc4,gcc5,gcc6,gcc7,gcc8,gcc9,gcc10))


##Section C
b1<-baseC %>% filter(patent_year<=2007)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n1b<-decompose(n1,min.vertices =gcc)
n1b<-n1b[[1]]
diam1<-diameter(n1b)
v1<-as.numeric(length(V(n1b)))
e1<-gsize(n1b)
gcc1<-c1
names(gcc1)[1]<-"size"
gcc1$year<-2007
gcc1$gp<-c("GCC","CC2","CC3")

b1<-baseC %>% filter(patent_year<=2008)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n2b<-decompose(n1,min.vertices =gcc)
n2b<-n2b[[1]]
diam2<-diameter(n2b)
v2<-as.numeric(length(V(n2b)))
e2<-gsize(n2b)
gcc2<-c1
names(gcc2)[1]<-"size"
gcc2$year<-2008
gcc2$gp<-c("GCC","CC2","CC3")

b1<-baseC %>% filter(patent_year<=2009)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n3b<-decompose(n1,min.vertices =gcc)
n3b<-n3b[[1]]
diam3<-diameter(n3b)
v3<-as.numeric(length(V(n3b)))
e3<-gsize(n3b)
gcc3<-c1
names(gcc3)[1]<-"size"
gcc3$year<-2009
gcc3$gp<-c("GCC","CC2","CC3")


b1<-baseC %>% filter(patent_year<=2010)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n4b<-decompose(n1,min.vertices =gcc)
n4b<-n4b[[1]]
diam4<-diameter(n4b)
v4<-as.numeric(length(V(n4b)))
e4<-gsize(n4b)
gcc4<-c1
names(gcc4)[1]<-"size"
gcc4$year<-2010
gcc4$gp<-c("GCC","CC2","CC3")


b1<-baseC %>% filter(patent_year<=2011)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n5b<-decompose(n1,min.vertices =gcc)
n5b<-n5b[[1]]
diam5<-diameter(n5b)
v5<-as.numeric(length(V(n5b)))
e5<-gsize(n5b)
gcc5<-c1
names(gcc5)[1]<-"size"
gcc5$year<-2011
gcc5$gp<-c("GCC","CC2","CC3")


b1<-baseC %>% filter(patent_year<=2012)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n6b<-decompose(n1,min.vertices =gcc)
n6b<-n6b[[1]]
diam6<-diameter(n6b)
v6<-as.numeric(length(V(n6b)))
e6<-gsize(n6b)
gcc6<-c1
names(gcc6)[1]<-"size"
gcc6$year<-2012
gcc6$gp<-c("GCC","CC2","CC3")


b1<-baseC %>% filter(patent_year<=2013)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n7b<-decompose(n1,min.vertices =gcc)
n7b<-n7b[[1]]
diam7<-diameter(n7b)
v7<-as.numeric(length(V(n7b)))
e7<-gsize(n7b)
gcc7<-c1
names(gcc7)[1]<-"size"
gcc7$year<-2013
gcc7$gp<-c("GCC","CC2","CC3")


b1<-baseC %>% filter(patent_year<=2014)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n8b<-decompose(n1,min.vertices =gcc)
n8b<-n8b[[1]]
diam8<-diameter(n8b)
v8<-as.numeric(length(V(n8b)))
e8<-gsize(n8b)
gcc8<-c1
names(gcc8)[1]<-"size"
gcc8$year<-2014
gcc8$gp<-c("GCC","CC2","CC3")

b1<-baseC %>% filter(patent_year<=2015)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n9b<-decompose(n1,min.vertices =gcc)
n9b<-n9b[[1]]
diam9<-diameter(n9b)
v9<-as.numeric(length(V(n9b)))
e9<-gsize(n9b)
gcc9<-c1
names(gcc9)[1]<-"size"
gcc9$year<-2015
gcc9$gp<-c("GCC","CC2","CC3")


b1<-baseC %>% filter(patent_year<=2016)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n10b<-decompose(n1,min.vertices =gcc)
n10b<-n10b[[1]]
diam10<-diameter(n10b)
v10<-as.numeric(length(V(n10b)))
e10<-gsize(n10b)
gcc10<-c1
names(gcc10)[1]<-"size"
gcc10$year<-2016
gcc10$gp<-c("GCC","CC2","CC3")

b1<-baseC %>% filter(patent_year<=2017)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n11b<-decompose(n1,min.vertices =gcc)
n11b<-n11b[[1]]
diam11<-diameter(n11b)
v11<-as.numeric(length(V(n11b)))
e11<-gsize(n11b)
gcc11<-c1
names(gcc11)[1]<-"size"
gcc11$year<-2017
gcc11$gp<-c("GCC","CC2","CC3")


b1<-baseC %>% filter(patent_year<=2018)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n12b<-decompose(n1,min.vertices =gcc)
n12b<-n12b[[1]]
diam12<-diameter(n12b)
v12<-as.numeric(length(V(n12b)))
e12<-gsize(n12b)
gcc12<-c1
names(gcc12)[1]<-"size"
gcc12$year<-2018
gcc12$gp<-c("GCC","CC2","CC3")

b1<-baseC %>% filter(patent_year<=2019)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n13b<-decompose(n1,min.vertices =gcc)
n13b<-n13b[[1]]
diam13<-diameter(n13b)
v13<-as.numeric(length(V(n13b)))
e13<-gsize(n13b)
gcc13<-c1
names(gcc13)[1]<-"size"
gcc13$year<-2019
gcc13$gp<-c("GCC","CC2","CC3")


year<-c(2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)
C<-c(diam1,diam2,diam3,diam4,diam5,diam6,diam7,diam8,diam9,diam10,diam11,diam12,diam13)
node<-c(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13)
edge<-c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13)
diamC<-data.frame(year,C)
neC<-data.frame(year,node,edge)

gpC<-Reduce(rbind.fill,list(gcc1,gcc2,gcc3,gcc4,gcc5,gcc6,gcc7,gcc8,gcc9,gcc10,gcc11,gcc12,gcc13))

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

### pfit most important component C
big <- decompose(n13b, min.vertices=gcc)
big<-big[[1]]
bb<-as_data_frame(big)
bb<-PAFit::to_igraph(Pnet)

bb <- decompose(bb, min.vertices=gcc)
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

deg1<-table(degree(n1b))
deg1<-as.data.frame(deg1)
deg1$cat<-2007
deg2<-table(degree(n7b))
deg2<-as.data.frame(deg2)
deg2$cat<-2013
deg3<-table(degree(n13b))
deg3<-as.data.frame(deg3)
deg3$cat<-2019
degC<-Reduce(rbind.fill,list(deg1,deg2,deg3))
degC$cat<-as.factor(degC$cat)
names(degC)[1]<-"degree"

b1<-baseC %>% filter(patent_year>=2010 & patent_year<=2010)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc1<-c1
gcc1$year<-2010
gcc1$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseC %>% filter(patent_year>=2010 & patent_year<=2011)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc2<-c1
gcc2$year<-2011
gcc2$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseC %>% filter(patent_year>=2010 & patent_year<=2012)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc3<-c1
gcc3$year<-2012
gcc3$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseC %>% filter(patent_year>=2010 & patent_year<=2013)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc4<-c1
gcc4$year<-2013
gcc4$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseC %>% filter(patent_year>=2010 & patent_year<=2014)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc5<-c1
gcc5$year<-2014
gcc5$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseC %>% filter(patent_year>=2010 & patent_year<=2015)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc6<-c1
gcc6$year<-2015
gcc6$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseC %>% filter(patent_year>=2010 & patent_year<=2016)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc7<-c1
gcc7$year<-2016
gcc7$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseC %>% filter(patent_year>=2010 & patent_year<=2017)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc8<-c1
gcc8$year<-2017
gcc8$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseC %>% filter(patent_year>=2010 & patent_year<=2018)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc9<-c1
gcc9$year<-2018
gcc9$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseC %>% filter(patent_year>=2010 & patent_year<=2019)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc10<-c1
gcc10$year<-2019
gcc10$gp<-"GCC Post 2010 subgraph, no past"


gpC<-Reduce(rbind.fill,list(gpC,gcc1,gcc2,gcc3,gcc4,gcc5,gcc6,gcc7,gcc8,gcc9,gcc10))

##Section D
b1<-baseD %>% filter(patent_year<=2007)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n1b<-decompose(n1,min.vertices =gcc)
n1b<-n1b[[1]]
diam1<-diameter(n1b)
v1<-as.numeric(length(V(n1b)))
e1<-gsize(n1b)
gcc1<-c1
names(gcc1)[1]<-"size"
gcc1$year<-2007
gcc1$gp<-c("GCC","CC2","CC3")

b1<-baseD %>% filter(patent_year<=2008)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n2b<-decompose(n1,min.vertices =gcc)
n2b<-n2b[[1]]
diam2<-diameter(n2b)
v2<-as.numeric(length(V(n2b)))
e2<-gsize(n2b)
gcc2<-c1
names(gcc2)[1]<-"size"
gcc2$year<-2008
gcc2$gp<-c("GCC","CC2","CC3")

b1<-baseD %>% filter(patent_year<=2009)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n3b<-decompose(n1,min.vertices =gcc)
n3b<-n3b[[1]]
diam3<-diameter(n3b)
v3<-as.numeric(length(V(n3b)))
e3<-gsize(n3b)
gcc3<-c1
names(gcc3)[1]<-"size"
gcc3$year<-2009
gcc3$gp<-c("GCC","CC2","CC3")


b1<-baseD %>% filter(patent_year<=2010)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n4b<-decompose(n1,min.vertices =gcc)
n4b<-n4b[[1]]
diam4<-diameter(n4b)
v4<-as.numeric(length(V(n4b)))
e4<-gsize(n4b)
gcc4<-c1
names(gcc4)[1]<-"size"
gcc4$year<-2010
gcc4$gp<-c("GCC","CC2","CC3")


b1<-baseD %>% filter(patent_year<=2011)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n5b<-decompose(n1,min.vertices =gcc)
n5b<-n5b[[1]]
diam5<-diameter(n5b)
v5<-as.numeric(length(V(n5b)))
e5<-gsize(n5b)
gcc5<-c1
names(gcc5)[1]<-"size"
gcc5$year<-2011
gcc5$gp<-c("GCC","CC2","CC3")


b1<-baseD %>% filter(patent_year<=2012)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n6b<-decompose(n1,min.vertices =gcc)
n6b<-n6b[[1]]
diam6<-diameter(n6b)
v6<-as.numeric(length(V(n6b)))
e6<-gsize(n6b)
gcc6<-c1
names(gcc6)[1]<-"size"
gcc6$year<-2012
gcc6$gp<-c("GCC","CC2","CC3")


b1<-baseD %>% filter(patent_year<=2013)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n7b<-decompose(n1,min.vertices =gcc)
n7b<-n7b[[1]]
diam7<-diameter(n7b)
v7<-as.numeric(length(V(n7b)))
e7<-gsize(n7b)
gcc7<-c1
names(gcc7)[1]<-"size"
gcc7$year<-2013
gcc7$gp<-c("GCC","CC2","CC3")


b1<-baseD %>% filter(patent_year<=2014)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n8b<-decompose(n1,min.vertices =gcc)
n8b<-n8b[[1]]
diam8<-diameter(n8b)
v8<-as.numeric(length(V(n8b)))
e8<-gsize(n8b)
gcc8<-c1
names(gcc8)[1]<-"size"
gcc8$year<-2014
gcc8$gp<-c("GCC","CC2","CC3")

b1<-baseD %>% filter(patent_year<=2015)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n9b<-decompose(n1,min.vertices =gcc)
n9b<-n9b[[1]]
diam9<-diameter(n9b)
v9<-as.numeric(length(V(n9b)))
e9<-gsize(n9b)
gcc9<-c1
names(gcc9)[1]<-"size"
gcc9$year<-2015
gcc9$gp<-c("GCC","CC2","CC3")


b1<-baseD %>% filter(patent_year<=2016)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n10b<-decompose(n1,min.vertices =gcc)
n10b<-n10b[[1]]
diam10<-diameter(n10b)
v10<-as.numeric(length(V(n10b)))
e10<-gsize(n10b)
gcc10<-c1
names(gcc10)[1]<-"size"
gcc10$year<-2016
gcc10$gp<-c("GCC","CC2","CC3")

b1<-baseD %>% filter(patent_year<=2017)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n11b<-decompose(n1,min.vertices =gcc)
n11b<-n11b[[1]]
diam11<-diameter(n11b)
v11<-as.numeric(length(V(n11b)))
e11<-gsize(n11b)
gcc11<-c1
names(gcc11)[1]<-"size"
gcc11$year<-2017
gcc11$gp<-c("GCC","CC2","CC3")


b1<-baseD %>% filter(patent_year<=2018)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n12b<-decompose(n1,min.vertices =gcc)
n12b<-n12b[[1]]
diam12<-diameter(n12b)
v12<-as.numeric(length(V(n12b)))
e12<-gsize(n12b)
gcc12<-c1
names(gcc12)[1]<-"size"
gcc12$year<-2018
gcc12$gp<-c("GCC","CC2","CC3")

b1<-baseD %>% filter(patent_year<=2019)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n13b<-decompose(n1,min.vertices =gcc)
n13b<-n13b[[1]]
diam13<-diameter(n13b)
v13<-as.numeric(length(V(n13b)))
e13<-gsize(n13b)
gcc13<-c1
names(gcc13)[1]<-"size"
gcc13$year<-2019
gcc13$gp<-c("GCC","CC2","CC3")


year<-c(2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)
D<-c(diam1,diam2,diam3,diam4,diam5,diam6,diam7,diam8,diam9,diam10,diam11,diam12,diam13)
node<-c(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13)
edge<-c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13)
diamD<-data.frame(year,D)
neD<-data.frame(year,node,edge)

gpD<-Reduce(rbind.fill,list(gcc1,gcc2,gcc3,gcc4,gcc5,gcc6,gcc7,gcc8,gcc9,gcc10,gcc11,gcc12,gcc13))

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

### pfit most important component D
big <- decompose(n13b, min.vertices=gcc)
big<-big[[1]]
bb<-as_data_frame(big)
bb<-PAFit::to_igraph(Pnet)

bb <- decompose(bb, min.vertices=gcc)
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

deg1<-table(degree(n1b))
deg1<-as.data.frame(deg1)
deg1$cat<-2007
deg2<-table(degree(n7b))
deg2<-as.data.frame(deg2)
deg2$cat<-2013
deg3<-table(degree(n13b))
deg3<-as.data.frame(deg3)
deg3$cat<-2019
degD<-Reduce(rbind.fill,list(deg1,deg2,deg3))
degD$cat<-as.factor(degD$cat)
names(degD)[1]<-"degree"

b1<-baseD %>% filter(patent_year>=2010 & patent_year<=2010)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc1<-c1
gcc1$year<-2010
gcc1$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseD %>% filter(patent_year>=2010 & patent_year<=2011)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc2<-c1
gcc2$year<-2011
gcc2$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseD %>% filter(patent_year>=2010 & patent_year<=2012)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc3<-c1
gcc3$year<-2012
gcc3$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseD %>% filter(patent_year>=2010 & patent_year<=2013)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc4<-c1
gcc4$year<-2013
gcc4$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseD %>% filter(patent_year>=2010 & patent_year<=2014)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc5<-c1
gcc5$year<-2014
gcc5$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseD %>% filter(patent_year>=2010 & patent_year<=2015)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc6<-c1
gcc6$year<-2015
gcc6$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseD %>% filter(patent_year>=2010 & patent_year<=2016)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc7<-c1
gcc7$year<-2016
gcc7$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseD %>% filter(patent_year>=2010 & patent_year<=2017)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc8<-c1
gcc8$year<-2017
gcc8$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseD %>% filter(patent_year>=2010 & patent_year<=2018)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc9<-c1
gcc9$year<-2018
gcc9$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseD %>% filter(patent_year>=2010 & patent_year<=2019)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc10<-c1
gcc10$year<-2019
gcc10$gp<-"GCC Post 2010 subgraph, no past"


gpD<-Reduce(rbind.fill,list(gpD,gcc1,gcc2,gcc3,gcc4,gcc5,gcc6,gcc7,gcc8,gcc9,gcc10))

##Section E
b1<-baseE %>% filter(patent_year<=2007)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n1b<-decompose(n1,min.vertices =gcc)
n1b<-n1b[[1]]
diam1<-diameter(n1b)
v1<-as.numeric(length(V(n1b)))
e1<-gsize(n1b)
gcc1<-c1
names(gcc1)[1]<-"size"
gcc1$year<-2007
gcc1$gp<-c("GCC","CC2","CC3")

b1<-baseE %>% filter(patent_year<=2008)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n2b<-decompose(n1,min.vertices =gcc)
n2b<-n2b[[1]]
diam2<-diameter(n2b)
v2<-as.numeric(length(V(n2b)))
e2<-gsize(n2b)
gcc2<-c1
names(gcc2)[1]<-"size"
gcc2$year<-2008
gcc2$gp<-c("GCC","CC2","CC3")

b1<-baseE %>% filter(patent_year<=2009)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n3b<-decompose(n1,min.vertices =gcc)
n3b<-n3b[[1]]
diam3<-diameter(n3b)
v3<-as.numeric(length(V(n3b)))
e3<-gsize(n3b)
gcc3<-c1
names(gcc3)[1]<-"size"
gcc3$year<-2009
gcc3$gp<-c("GCC","CC2","CC3")


b1<-baseE %>% filter(patent_year<=2010)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n4b<-decompose(n1,min.vertices =gcc)
n4b<-n4b[[1]]
diam4<-diameter(n4b)
v4<-as.numeric(length(V(n4b)))
e4<-gsize(n4b)
gcc4<-c1
names(gcc4)[1]<-"size"
gcc4$year<-2010
gcc4$gp<-c("GCC","CC2","CC3")


b1<-baseE %>% filter(patent_year<=2011)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n5b<-decompose(n1,min.vertices =gcc)
n5b<-n5b[[1]]
diam5<-diameter(n5b)
v5<-as.numeric(length(V(n5b)))
e5<-gsize(n5b)
gcc5<-c1
names(gcc5)[1]<-"size"
gcc5$year<-2011
gcc5$gp<-c("GCC","CC2","CC3")


b1<-baseE %>% filter(patent_year<=2012)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n6b<-decompose(n1,min.vertices =gcc)
n6b<-n6b[[1]]
diam6<-diameter(n6b)
v6<-as.numeric(length(V(n6b)))
e6<-gsize(n6b)
gcc6<-c1
names(gcc6)[1]<-"size"
gcc6$year<-2012
gcc6$gp<-c("GCC","CC2","CC3")


b1<-baseE %>% filter(patent_year<=2013)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n7b<-decompose(n1,min.vertices =gcc)
n7b<-n7b[[1]]
diam7<-diameter(n7b)
v7<-as.numeric(length(V(n7b)))
e7<-gsize(n7b)
gcc7<-c1
names(gcc7)[1]<-"size"
gcc7$year<-2013
gcc7$gp<-c("GCC","CC2","CC3")


b1<-baseE %>% filter(patent_year<=2014)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n8b<-decompose(n1,min.vertices =gcc)
n8b<-n8b[[1]]
diam8<-diameter(n8b)
v8<-as.numeric(length(V(n8b)))
e8<-gsize(n8b)
gcc8<-c1
names(gcc8)[1]<-"size"
gcc8$year<-2014
gcc8$gp<-c("GCC","CC2","CC3")

b1<-baseE %>% filter(patent_year<=2015)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n9b<-decompose(n1,min.vertices =gcc)
n9b<-n9b[[1]]
diam9<-diameter(n9b)
v9<-as.numeric(length(V(n9b)))
e9<-gsize(n9b)
gcc9<-c1
names(gcc9)[1]<-"size"
gcc9$year<-2015
gcc9$gp<-c("GCC","CC2","CC3")


b1<-baseE %>% filter(patent_year<=2016)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n10b<-decompose(n1,min.vertices =gcc)
n10b<-n10b[[1]]
diam10<-diameter(n10b)
v10<-as.numeric(length(V(n10b)))
e10<-gsize(n10b)
gcc10<-c1
names(gcc10)[1]<-"size"
gcc10$year<-2016
gcc10$gp<-c("GCC","CC2","CC3")

b1<-baseE %>% filter(patent_year<=2017)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n11b<-decompose(n1,min.vertices =gcc)
n11b<-n11b[[1]]
diam11<-diameter(n11b)
v11<-as.numeric(length(V(n11b)))
e11<-gsize(n11b)
gcc11<-c1
names(gcc11)[1]<-"size"
gcc11$year<-2017
gcc11$gp<-c("GCC","CC2","CC3")


b1<-baseE %>% filter(patent_year<=2018)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n12b<-decompose(n1,min.vertices =gcc)
n12b<-n12b[[1]]
diam12<-diameter(n12b)
v12<-as.numeric(length(V(n12b)))
e12<-gsize(n12b)
gcc12<-c1
names(gcc12)[1]<-"size"
gcc12$year<-2018
gcc12$gp<-c("GCC","CC2","CC3")

b1<-baseE %>% filter(patent_year<=2019)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n13b<-decompose(n1,min.vertices =gcc)
n13b<-n13b[[1]]
diam13<-diameter(n13b)
v13<-as.numeric(length(V(n13b)))
e13<-gsize(n13b)
gcc13<-c1
names(gcc13)[1]<-"size"
gcc13$year<-2019
gcc13$gp<-c("GCC","CC2","CC3")

year<-c(2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)
E<-c(diam1,diam2,diam3,diam4,diam5,diam6,diam7,diam8,diam9,diam10,diam11,diam12,diam13)
node<-c(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13)
edge<-c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13)
diamE<-data.frame(year,E)
neE<-data.frame(year,node,edge)

gpE<-Reduce(rbind.fill,list(gcc1,gcc2,gcc3,gcc4,gcc5,gcc6,gcc7,gcc8,gcc9,gcc10,gcc11,gcc12,gcc13))

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

### pfit most important component E
big <- decompose(n13b, min.vertices=gcc)
big<-big[[1]]
bb<-as_data_frame(big)
bb<-PAFit::to_igraph(Pnet)

bb <- decompose(bb, min.vertices=gcc)
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

deg1<-table(degree(n1b))
deg1<-as.data.frame(deg1)
deg1$cat<-2007
deg2<-table(degree(n7b))
deg2<-as.data.frame(deg2)
deg2$cat<-2013
deg3<-table(degree(n13b))
deg3<-as.data.frame(deg3)
deg3$cat<-2019
degE<-Reduce(rbind.fill,list(deg1,deg2,deg3))
degE$cat<-as.factor(degE$cat)
names(degE)[1]<-"degree"

b1<-baseE %>% filter(patent_year>=2010 & patent_year<=2010)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc1<-c1
gcc1$year<-2010
gcc1$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseE %>% filter(patent_year>=2010 & patent_year<=2011)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc2<-c1
gcc2$year<-2011
gcc2$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseE %>% filter(patent_year>=2010 & patent_year<=2012)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc3<-c1
gcc3$year<-2012
gcc3$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseE %>% filter(patent_year>=2010 & patent_year<=2013)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc4<-c1
gcc4$year<-2013
gcc4$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseE %>% filter(patent_year>=2010 & patent_year<=2014)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc5<-c1
gcc5$year<-2014
gcc5$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseE %>% filter(patent_year>=2010 & patent_year<=2015)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc6<-c1
gcc6$year<-2015
gcc6$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseE %>% filter(patent_year>=2010 & patent_year<=2016)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc7<-c1
gcc7$year<-2016
gcc7$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseE %>% filter(patent_year>=2010 & patent_year<=2017)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc8<-c1
gcc8$year<-2017
gcc8$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseE %>% filter(patent_year>=2010 & patent_year<=2018)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc9<-c1
gcc9$year<-2018
gcc9$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseE %>% filter(patent_year>=2010 & patent_year<=2019)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc10<-c1
gcc10$year<-2019
gcc10$gp<-"GCC Post 2010 subgraph, no past"


gpE<-Reduce(rbind.fill,list(gpE,gcc1,gcc2,gcc3,gcc4,gcc5,gcc6,gcc7,gcc8,gcc9,gcc10))

##Section F
b1<-baseF %>% filter(patent_year<=2007)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n1b<-decompose(n1,min.vertices =gcc)
n1b<-n1b[[1]]
diam1<-diameter(n1b)
v1<-as.numeric(length(V(n1b)))
e1<-gsize(n1b)
gcc1<-c1
names(gcc1)[1]<-"size"
gcc1$year<-2007
gcc1$gp<-c("GCC","CC2","CC3")

b1<-baseF %>% filter(patent_year<=2008)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n2b<-decompose(n1,min.vertices =gcc)
n2b<-n2b[[1]]
diam2<-diameter(n2b)
v2<-as.numeric(length(V(n2b)))
e2<-gsize(n2b)
gcc2<-c1
names(gcc2)[1]<-"size"
gcc2$year<-2008
gcc2$gp<-c("GCC","CC2","CC3")

b1<-baseF %>% filter(patent_year<=2009)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n3b<-decompose(n1,min.vertices =gcc)
n3b<-n3b[[1]]
diam3<-diameter(n3b)
v3<-as.numeric(length(V(n3b)))
e3<-gsize(n3b)
gcc3<-c1
names(gcc3)[1]<-"size"
gcc3$year<-2009
gcc3$gp<-c("GCC","CC2","CC3")


b1<-baseF %>% filter(patent_year<=2010)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n4b<-decompose(n1,min.vertices =gcc)
n4b<-n4b[[1]]
diam4<-diameter(n4b)
v4<-as.numeric(length(V(n4b)))
e4<-gsize(n4b)
gcc4<-c1
names(gcc4)[1]<-"size"
gcc4$year<-2010
gcc4$gp<-c("GCC","CC2","CC3")


b1<-baseF %>% filter(patent_year<=2011)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n5b<-decompose(n1,min.vertices =gcc)
n5b<-n5b[[1]]
diam5<-diameter(n5b)
v5<-as.numeric(length(V(n5b)))
e5<-gsize(n5b)
gcc5<-c1
names(gcc5)[1]<-"size"
gcc5$year<-2011
gcc5$gp<-c("GCC","CC2","CC3")


b1<-baseF %>% filter(patent_year<=2012)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n6b<-decompose(n1,min.vertices =gcc)
n6b<-n6b[[1]]
diam6<-diameter(n6b)
v6<-as.numeric(length(V(n6b)))
e6<-gsize(n6b)
gcc6<-c1
names(gcc6)[1]<-"size"
gcc6$year<-2012
gcc6$gp<-c("GCC","CC2","CC3")


b1<-baseF %>% filter(patent_year<=2013)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n7b<-decompose(n1,min.vertices =gcc)
n7b<-n7b[[1]]
diam7<-diameter(n7b)
v7<-as.numeric(length(V(n7b)))
e7<-gsize(n7b)
gcc7<-c1
names(gcc7)[1]<-"size"
gcc7$year<-2013
gcc7$gp<-c("GCC","CC2","CC3")


b1<-baseF %>% filter(patent_year<=2014)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n8b<-decompose(n1,min.vertices =gcc)
n8b<-n8b[[1]]
diam8<-diameter(n8b)
v8<-as.numeric(length(V(n8b)))
e8<-gsize(n8b)
gcc8<-c1
names(gcc8)[1]<-"size"
gcc8$year<-2014
gcc8$gp<-c("GCC","CC2","CC3")

b1<-baseF %>% filter(patent_year<=2015)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n9b<-decompose(n1,min.vertices =gcc)
n9b<-n9b[[1]]
diam9<-diameter(n9b)
v9<-as.numeric(length(V(n9b)))
e9<-gsize(n9b)
gcc9<-c1
names(gcc9)[1]<-"size"
gcc9$year<-2015
gcc9$gp<-c("GCC","CC2","CC3")


b1<-baseF %>% filter(patent_year<=2016)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n10b<-decompose(n1,min.vertices =gcc)
n10b<-n10b[[1]]
diam10<-diameter(n10b)
v10<-as.numeric(length(V(n10b)))
e10<-gsize(n10b)
gcc10<-c1
names(gcc10)[1]<-"size"
gcc10$year<-2016
gcc10$gp<-c("GCC","CC2","CC3")

b1<-baseF %>% filter(patent_year<=2017)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n11b<-decompose(n1,min.vertices =gcc)
n11b<-n11b[[1]]
diam11<-diameter(n11b)
v11<-as.numeric(length(V(n11b)))
e11<-gsize(n11b)
gcc11<-c1
names(gcc11)[1]<-"size"
gcc11$year<-2017
gcc11$gp<-c("GCC","CC2","CC3")

b1<-baseF %>% filter(patent_year<=2018)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n12b<-decompose(n1,min.vertices =gcc)
n12b<-n12b[[1]]
diam12<-diameter(n12b)
v12<-as.numeric(length(V(n12b)))
e12<-gsize(n12b)
gcc12<-c1
names(gcc12)[1]<-"size"
gcc12$year<-2018
gcc12$gp<-c("GCC","CC2","CC3")

b1<-baseF %>% filter(patent_year<=2019)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n13b<-decompose(n1,min.vertices =gcc)
n13b<-n13b[[1]]
diam13<-diameter(n13b)
v13<-as.numeric(length(V(n13b)))
e13<-gsize(n13b)
gcc13<-c1
names(gcc13)[1]<-"size"
gcc13$year<-2019
gcc13$gp<-c("GCC","CC2","CC3")

year<-c(2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)
FF<-c(diam1,diam2,diam3,diam4,diam5,diam6,diam7,diam8,diam9,diam10,diam11,diam12,diam13)
node<-c(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13)
edge<-c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13)
diamF<-data.frame(year,FF)
neF<-data.frame(year,node,edge)

gpF<-Reduce(rbind.fill,list(gcc1,gcc2,gcc3,gcc4,gcc5,gcc6,gcc7,gcc8,gcc9,gcc10,gcc11,gcc12,gcc13))

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

### pfit most important component F
big <- decompose(n13b, min.vertices=gcc)
big<-big[[1]]
bb<-as_data_frame(big)
bb<-PAFit::to_igraph(Pnet)

bb <- decompose(bb, min.vertices=gcc)
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

deg1<-table(degree(n1b))
deg1<-as.data.frame(deg1)
deg1$cat<-2007
deg2<-table(degree(n7b))
deg2<-as.data.frame(deg2)
deg2$cat<-2013
deg3<-table(degree(n13b))
deg3<-as.data.frame(deg3)
deg3$cat<-2019
degF<-Reduce(rbind.fill,list(deg1,deg2,deg3))
degF$cat<-as.factor(degF$cat)
names(degF)[1]<-"degree"

b1<-baseF %>% filter(patent_year>=2010 & patent_year<=2010)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc1<-c1
gcc1$year<-2010
gcc1$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseF %>% filter(patent_year>=2010 & patent_year<=2011)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc2<-c1
gcc2$year<-2011
gcc2$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseF %>% filter(patent_year>=2010 & patent_year<=2012)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc3<-c1
gcc3$year<-2012
gcc3$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseF %>% filter(patent_year>=2010 & patent_year<=2013)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc4<-c1
gcc4$year<-2013
gcc4$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseF %>% filter(patent_year>=2010 & patent_year<=2014)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc5<-c1
gcc5$year<-2014
gcc5$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseF %>% filter(patent_year>=2010 & patent_year<=2015)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc6<-c1
gcc6$year<-2015
gcc6$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseF %>% filter(patent_year>=2010 & patent_year<=2016)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc7<-c1
gcc7$year<-2016
gcc7$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseF %>% filter(patent_year>=2010 & patent_year<=2017)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc8<-c1
gcc8$year<-2017
gcc8$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseF %>% filter(patent_year>=2010 & patent_year<=2018)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc9<-c1
gcc9$year<-2018
gcc9$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseF %>% filter(patent_year>=2010 & patent_year<=2019)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc10<-c1
gcc10$year<-2019
gcc10$gp<-"GCC Post 2010 subgraph, no past"


gpF<-Reduce(rbind.fill,list(gpF,gcc1,gcc2,gcc3,gcc4,gcc5,gcc6,gcc7,gcc8,gcc9,gcc10))

##Section G
b1<-baseG %>% filter(patent_year<=2007)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n1b<-decompose(n1,min.vertices =gcc)
n1b<-n1b[[1]]
diam1<-diameter(n1b)
v1<-as.numeric(length(V(n1b)))
e1<-gsize(n1b)
gcc1<-c1
names(gcc1)[1]<-"size"
gcc1$year<-2007
gcc1$gp<-c("GCC","CC2","CC3")

b1<-baseG %>% filter(patent_year<=2008)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n2b<-decompose(n1,min.vertices =gcc)
n2b<-n2b[[1]]
diam2<-diameter(n2b)
v2<-as.numeric(length(V(n2b)))
e2<-gsize(n2b)
gcc2<-c1
names(gcc2)[1]<-"size"
gcc2$year<-2008
gcc2$gp<-c("GCC","CC2","CC3")

b1<-baseG %>% filter(patent_year<=2009)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n3b<-decompose(n1,min.vertices =gcc)
n3b<-n3b[[1]]
diam3<-diameter(n3b)
v3<-as.numeric(length(V(n3b)))
e3<-gsize(n3b)
gcc3<-c1
names(gcc3)[1]<-"size"
gcc3$year<-2009
gcc3$gp<-c("GCC","CC2","CC3")


b1<-baseG %>% filter(patent_year<=2010)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n4b<-decompose(n1,min.vertices =gcc)
n4b<-n4b[[1]]
diam4<-diameter(n4b)
v4<-as.numeric(length(V(n4b)))
e4<-gsize(n4b)
gcc4<-c1
names(gcc4)[1]<-"size"
gcc4$year<-2010
gcc4$gp<-c("GCC","CC2","CC3")


b1<-baseG %>% filter(patent_year<=2011)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n5b<-decompose(n1,min.vertices =gcc)
n5b<-n5b[[1]]
diam5<-diameter(n5b)
v5<-as.numeric(length(V(n5b)))
e5<-gsize(n5b)
gcc5<-c1
names(gcc5)[1]<-"size"
gcc5$year<-2011
gcc5$gp<-c("GCC","CC2","CC3")


b1<-baseG %>% filter(patent_year<=2012)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n6b<-decompose(n1,min.vertices =gcc)
n6b<-n6b[[1]]
diam6<-diameter(n6b)
v6<-as.numeric(length(V(n6b)))
e6<-gsize(n6b)
gcc6<-c1
names(gcc6)[1]<-"size"
gcc6$year<-2012
gcc6$gp<-c("GCC","CC2","CC3")


b1<-baseG %>% filter(patent_year<=2013)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n7b<-decompose(n1,min.vertices =gcc)
n7b<-n7b[[1]]
diam7<-diameter(n7b)
v7<-as.numeric(length(V(n7b)))
e7<-gsize(n7b)
gcc7<-c1
names(gcc7)[1]<-"size"
gcc7$year<-2013
gcc7$gp<-c("GCC","CC2","CC3")


b1<-baseG %>% filter(patent_year<=2014)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n8b<-decompose(n1,min.vertices =gcc)
n8b<-n8b[[1]]
diam8<-diameter(n8b)
v8<-as.numeric(length(V(n8b)))
e8<-gsize(n8b)
gcc8<-c1
names(gcc8)[1]<-"size"
gcc8$year<-2014
gcc8$gp<-c("GCC","CC2","CC3")

b1<-baseG %>% filter(patent_year<=2015)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n9b<-decompose(n1,min.vertices =gcc)
n9b<-n9b[[1]]
diam9<-diameter(n9b)
v9<-as.numeric(length(V(n9b)))
e9<-gsize(n9b)
gcc9<-c1
names(gcc9)[1]<-"size"
gcc9$year<-2015
gcc9$gp<-c("GCC","CC2","CC3")


b1<-baseG %>% filter(patent_year<=2016)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n10b<-decompose(n1,min.vertices =gcc)
n10b<-n10b[[1]]
diam10<-diameter(n10b)
v10<-as.numeric(length(V(n10b)))
e10<-gsize(n10b)
gcc10<-c1
names(gcc10)[1]<-"size"
gcc10$year<-2016
gcc10$gp<-c("GCC","CC2","CC3")

b1<-baseG %>% filter(patent_year<=2017)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n11b<-decompose(n1,min.vertices =gcc)
n11b<-n11b[[1]]
diam11<-diameter(n11b)
v11<-as.numeric(length(V(n11b)))
e11<-gsize(n11b)
gcc11<-c1
names(gcc11)[1]<-"size"
gcc11$year<-2017
gcc11$gp<-c("GCC","CC2","CC3")


b1<-baseG %>% filter(patent_year<=2018)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n12b<-decompose(n1,min.vertices =gcc)
n12b<-n12b[[1]]
diam12<-diameter(n12b)
v12<-as.numeric(length(V(n12b)))
e12<-gsize(n12b)
gcc12<-c1
names(gcc12)[1]<-"size"
gcc12$year<-2018
gcc12$gp<-c("GCC","CC2","CC3")

b1<-baseG %>% filter(patent_year<=2019)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n13b<-decompose(n1,min.vertices =gcc)
n13b<-n13b[[1]]
diam13<-diameter(n13b)
v13<-as.numeric(length(V(n13b)))
e13<-gsize(n13b)
gcc13<-c1
names(gcc13)[1]<-"size"
gcc13$year<-2019
gcc13$gp<-c("GCC","CC2","CC3")

year<-c(2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)
G<-c(diam1,diam2,diam3,diam4,diam5,diam6,diam7,diam8,diam9,diam10,diam11,diam12,diam13)
node<-c(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13)
edge<-c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13)
diamG<-data.frame(year,G)
neG<-data.frame(year,node,edge)

gpG<-Reduce(rbind.fill,list(gcc1,gcc2,gcc3,gcc4,gcc5,gcc6,gcc7,gcc8,gcc9,gcc10,gcc11,gcc12,gcc13))

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

### pfit most important component G
big <- decompose(n13b, min.vertices=gcc)
big<-big[[1]]
bb<-as_data_frame(big)
bb<-PAFit::to_igraph(Pnet)

bb <- decompose(bb, min.vertices=gcc)
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

deg1<-table(degree(n1b))
deg1<-as.data.frame(deg1)
deg1$cat<-2007
deg2<-table(degree(n7b))
deg2<-as.data.frame(deg2)
deg2$cat<-2013
deg3<-table(degree(n13b))
deg3<-as.data.frame(deg3)
deg3$cat<-2019
degG<-Reduce(rbind.fill,list(deg1,deg2,deg3))
degG$cat<-as.factor(degG$cat)
names(degG)[1]<-"degree"

b1<-baseG %>% filter(patent_year>=2010 & patent_year<=2010)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc1<-c1
gcc1$year<-2010
gcc1$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseG %>% filter(patent_year>=2010 & patent_year<=2011)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc2<-c1
gcc2$year<-2011
gcc2$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseG %>% filter(patent_year>=2010 & patent_year<=2012)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc3<-c1
gcc3$year<-2012
gcc3$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseG %>% filter(patent_year>=2010 & patent_year<=2013)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc4<-c1
gcc4$year<-2013
gcc4$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseG %>% filter(patent_year>=2010 & patent_year<=2014)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc5<-c1
gcc5$year<-2014
gcc5$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseG %>% filter(patent_year>=2010 & patent_year<=2015)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc6<-c1
gcc6$year<-2015
gcc6$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseG %>% filter(patent_year>=2010 & patent_year<=2016)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc7<-c1
gcc7$year<-2016
gcc7$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseG %>% filter(patent_year>=2010 & patent_year<=2017)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc8<-c1
gcc8$year<-2017
gcc8$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseG %>% filter(patent_year>=2010 & patent_year<=2018)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc9<-c1
gcc9$year<-2018
gcc9$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseG %>% filter(patent_year>=2010 & patent_year<=2019)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc10<-c1
gcc10$year<-2019
gcc10$gp<-"GCC Post 2010 subgraph, no past"


gpG<-Reduce(rbind.fill,list(gpG,gcc1,gcc2,gcc3,gcc4,gcc5,gcc6,gcc7,gcc8,gcc9,gcc10))

##Section H
b1<-baseH %>% filter(patent_year<=2007)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n1b<-decompose(n1,min.vertices =gcc)
n1b<-n1b[[1]]
diam1<-diameter(n1b)
v1<-as.numeric(length(V(n1b)))
e1<-gsize(n1b)
gcc1<-c1
names(gcc1)[1]<-"size"
gcc1$year<-2007
gcc1$gp<-c("GCC","CC2","CC3")

b1<-baseH %>% filter(patent_year<=2008)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n2b<-decompose(n1,min.vertices =gcc)
n2b<-n2b[[1]]
diam2<-diameter(n2b)
v2<-as.numeric(length(V(n2b)))
e2<-gsize(n2b)
gcc2<-c1
names(gcc2)[1]<-"size"
gcc2$year<-2008
gcc2$gp<-c("GCC","CC2","CC3")

b1<-baseH %>% filter(patent_year<=2009)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n3b<-decompose(n1,min.vertices =gcc)
n3b<-n3b[[1]]
diam3<-diameter(n3b)
v3<-as.numeric(length(V(n3b)))
e3<-gsize(n3b)
gcc3<-c1
names(gcc3)[1]<-"size"
gcc3$year<-2009
gcc3$gp<-c("GCC","CC2","CC3")


b1<-baseH %>% filter(patent_year<=2010)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n4b<-decompose(n1,min.vertices =gcc)
n4b<-n4b[[1]]
diam4<-diameter(n4b)
v4<-as.numeric(length(V(n4b)))
e4<-gsize(n4b)
gcc4<-c1
names(gcc4)[1]<-"size"
gcc4$year<-2010
gcc4$gp<-c("GCC","CC2","CC3")


b1<-baseH %>% filter(patent_year<=2011)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n5b<-decompose(n1,min.vertices =gcc)
n5b<-n5b[[1]]
diam5<-diameter(n5b)
v5<-as.numeric(length(V(n5b)))
e5<-gsize(n5b)
gcc5<-c1
names(gcc5)[1]<-"size"
gcc5$year<-2011
gcc5$gp<-c("GCC","CC2","CC3")


b1<-baseH %>% filter(patent_year<=2012)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n6b<-decompose(n1,min.vertices =gcc)
n6b<-n6b[[1]]
diam6<-diameter(n6b)
v6<-as.numeric(length(V(n6b)))
e6<-gsize(n6b)
gcc6<-c1
names(gcc6)[1]<-"size"
gcc6$year<-2012
gcc6$gp<-c("GCC","CC2","CC3")


b1<-baseH %>% filter(patent_year<=2013)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n7b<-decompose(n1,min.vertices =gcc)
n7b<-n7b[[1]]
diam7<-diameter(n7b)
v7<-as.numeric(length(V(n7b)))
e7<-gsize(n7b)
gcc7<-c1
names(gcc7)[1]<-"size"
gcc7$year<-2013
gcc7$gp<-c("GCC","CC2","CC3")


b1<-baseH %>% filter(patent_year<=2014)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n8b<-decompose(n1,min.vertices =gcc)
n8b<-n8b[[1]]
diam8<-diameter(n8b)
v8<-as.numeric(length(V(n8b)))
e8<-gsize(n8b)
gcc8<-c1
names(gcc8)[1]<-"size"
gcc8$year<-2014
gcc8$gp<-c("GCC","CC2","CC3")

b1<-baseH %>% filter(patent_year<=2015)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n9b<-decompose(n1,min.vertices =gcc)
n9b<-n9b[[1]]
diam9<-diameter(n9b)
v9<-as.numeric(length(V(n9b)))
e9<-gsize(n9b)
gcc9<-c1
names(gcc9)[1]<-"size"
gcc9$year<-2015
gcc9$gp<-c("GCC","CC2","CC3")


b1<-baseH %>% filter(patent_year<=2016)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n10b<-decompose(n1,min.vertices =gcc)
n10b<-n10b[[1]]
diam10<-diameter(n10b)
v10<-as.numeric(length(V(n10b)))
e10<-gsize(n10b)
gcc10<-c1
names(gcc10)[1]<-"size"
gcc10$year<-2016
gcc10$gp<-c("GCC","CC2","CC3")

b1<-baseH %>% filter(patent_year<=2017)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n11b<-decompose(n1,min.vertices =gcc)
n11b<-n11b[[1]]
diam11<-diameter(n11b)
v11<-as.numeric(length(V(n11b)))
e11<-gsize(n11b)
gcc11<-c1
names(gcc11)[1]<-"size"
gcc11$year<-2017
gcc11$gp<-c("GCC","CC2","CC3")


b1<-baseH %>% filter(patent_year<=2018)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n12b<-decompose(n1,min.vertices =gcc)
n12b<-n12b[[1]]
diam12<-diameter(n12b)
v12<-as.numeric(length(V(n12b)))
e12<-gsize(n12b)
gcc12<-c1
names(gcc12)[1]<-"size"
gcc12$year<-2018
gcc12$gp<-c("GCC","CC2","CC3")

b1<-baseH %>% filter(patent_year<=2019)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1:3)
c1<-as.data.frame(c1)
gcc<-c1[1,1]
n13b<-decompose(n1,min.vertices =gcc)
n13b<-n13b[[1]]
diam13<-diameter(n13b)
v13<-as.numeric(length(V(n13b)))
e13<-gsize(n13b)
gcc13<-c1
names(gcc13)[1]<-"size"
gcc13$year<-2019
gcc13$gp<-c("GCC","CC2","CC3")

year<-c(2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019)
H<-c(diam1,diam2,diam3,diam4,diam5,diam6,diam7,diam8,diam9,diam10,diam11,diam12,diam13)
node<-c(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13)
edge<-c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13)
diamH<-data.frame(year,H)
neH<-data.frame(year,node,edge)

gpH<-Reduce(rbind.fill,list(gcc1,gcc2,gcc3,gcc4,gcc5,gcc6,gcc7,gcc8,gcc9,gcc10,gcc11,gcc12,gcc13))

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

### pfit most important component H
big <- decompose(n13b, min.vertices=gcc)
big<-big[[1]]
bb<-as_data_frame(big)
bb<-PAFit::to_igraph(Pnet)

bb <- decompose(bb, min.vertices=gcc)
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

deg1<-table(degree(n1b))
deg1<-as.data.frame(deg1)
deg1$cat<-2007
deg2<-table(degree(n7b))
deg2<-as.data.frame(deg2)
deg2$cat<-2013
deg3<-table(degree(n13b))
deg3<-as.data.frame(deg3)
deg3$cat<-2019
degH<-Reduce(rbind.fill,list(deg1,deg2,deg3))
degH$cat<-as.factor(degH$cat)
names(degH)[1]<-"degree"

b1<-baseH %>% filter(patent_year>=2010 & patent_year<=2010)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc1<-c1
gcc1$year<-2010
gcc1$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseH %>% filter(patent_year>=2010 & patent_year<=2011)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc2<-c1
gcc2$year<-2011
gcc2$gp<-"GCC Post 2010 subgraph, no past"


b1<-baseH %>% filter(patent_year>=2010 & patent_year<=2012)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc3<-c1
gcc3$year<-2012
gcc3$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseH %>% filter(patent_year>=2010 & patent_year<=2013)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc4<-c1
gcc4$year<-2013
gcc4$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseH %>% filter(patent_year>=2010 & patent_year<=2014)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc5<-c1
gcc5$year<-2014
gcc5$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseH %>% filter(patent_year>=2010 & patent_year<=2015)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc6<-c1
gcc6$year<-2015
gcc6$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseH %>% filter(patent_year>=2010 & patent_year<=2016)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc7<-c1
gcc7$year<-2016
gcc7$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseH %>% filter(patent_year>=2010 & patent_year<=2017)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc8<-c1
gcc8$year<-2017
gcc8$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseH %>% filter(patent_year>=2010 & patent_year<=2018)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc9<-c1
gcc9$year<-2018
gcc9$gp<-"GCC Post 2010 subgraph, no past"

b1<-baseH %>% filter(patent_year>=2010 & patent_year<=2019)
z<-data.frame(b1$ff,b1$tt)
n1<-graph.data.frame(z,directed = F)
n1<-simplify(n1)
c1<-components(n1)
c1<-table(c1$csize)
c1<-as.data.frame(c1)
c1<-select(c1,Var1)
c1$Var1<-as.character(c1$Var1)
c1$Var1<-as.numeric(c1$Var1)
c1<-rev(c1$Var1)
c1<-as.data.frame(c1)
c1<-c1 %>% slice(1)
c1<-as.data.frame(c1)
names(c1)[1]<-"size"
gcc10<-c1
gcc10$year<-2019
gcc10$gp<-"GCC Post 2010 subgraph, no past"

gpH<-Reduce(rbind.fill,list(gpH,gcc1,gcc2,gcc3,gcc4,gcc5,gcc6,gcc7,gcc8,gcc9,gcc10))

##diameter,section
AA<-t(diamA)
AA<-as.data.frame(AA)
AA<-AA %>% slice(2)
dA<-as.data.frame("A")
names(dA)[1]<-"section"
dA<-cbind(dA,AA)
dA$section<-as.character(dA$section)

BB<-t(diamB)
BB<-as.data.frame(BB)
BB<-BB %>% slice(2)
dB<-as.data.frame("B")
names(dB)[1]<-"section"
dB<-cbind(dB,BB)
dB$section<-as.character(dB$section)

CC<-t(diamC)
CC<-as.data.frame(CC)
CC<-CC %>% slice(2)
dC<-as.data.frame("C")
names(dC)[1]<-"section"
dC<-cbind(dC,CC)
dC$section<-as.character(dC$section)

DD<-t(diamD)
DD<-as.data.frame(DD)
DD<-DD %>% slice(2)
dD<-as.data.frame("D")
names(dD)[1]<-"section"
dD<-cbind(dD,DD)
dD$section<-as.character(dD$section)

EE<-t(diamE)
EE<-as.data.frame(EE)
EE<-EE %>% slice(2)
dE<-as.data.frame("E")
names(dE)[1]<-"section"
dE<-cbind(dE,EE)
dE$section<-as.character(dE$section)

FF<-t(diamF)
FF<-as.data.frame(FF)
FF<-FF %>% slice(2)
dF<-as.data.frame("F")
names(dF)[1]<-"section"
dF<-cbind(dF,FF)
dF$section<-as.character(dF$section)

GG<-t(diamG)
GG<-as.data.frame(GG)
GG<-GG %>% slice(2)
dG<-as.data.frame("G")
names(dG)[1]<-"section"
dG<-cbind(dG,GG)
dG$section<-as.character(dG$section)

HH<-t(diamH)
HH<-as.data.frame(HH)
HH<-HH %>% slice(2)
dH<-as.data.frame("H")
names(dH)[1]<-"section"
dH<-cbind(dH,HH)
dH$section<-as.character(dH$section)

diam<-Reduce(rbind.fill,list(dA,dB,dC,dD,dE,dF,dG,dH))

names(diam)[2]<-"2007"
names(diam)[3]<-"2008"
names(diam)[4]<-"2009"
names(diam)[5]<-"2010"
names(diam)[6]<-"2011"
names(diam)[7]<-"2012"
names(diam)[8]<-"2013"
names(diam)[9]<-"2014"
names(diam)[10]<-"2015"
names(diam)[11]<-"2016"
names(diam)[12]<-"2017"
names(diam)[13]<-"2018"
names(diam)[14]<-"2019"

gpL$class<-"Total"
gpA$class<-"A"
gpB$class<-"B"
gpC$class<-"C"
gpD$class<-"D"
gpE$class<-"E"
gpF$class<-"F"
gpG$class<-"G"
gpH$class<-"H"

gR<-Reduce(rbind.fill,list(gpL,gpA,gpB,gpC,gpD,gpE,gpF,gpG,gpH))

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

#Figure 2

ggplot(degL, aes(x=degree, y=Freq, col=cat)) + 
  geom_point(aes(colour=cat))+
  scale_y_continuous("Frequency", trans = "log10",labels = trans_format("log10", math_format(10^.x)))+
  scale_x_continuous("Degree", trans = "log10",labels = trans_format("log10", math_format(10^.x)))+
  scale_color_manual(values = c("gray75", "gray45", "gray10"))+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major=element_blank(),panel.grid.minor=element_blank())

ggplot(degA, aes(x=degree, y=Freq, col=cat)) + 
  geom_point(aes(colour=cat))+
  scale_y_continuous("Frequency", trans = "log10",labels = trans_format("log10", math_format(10^.x)))+
  scale_x_continuous("Degree", trans = "log10",labels = trans_format("log10", math_format(10^.x)))+
  scale_color_manual(values = c("gray75", "gray45", "gray10"))+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major=element_blank(),panel.grid.minor=element_blank())

ggplot(degB, aes(x=degree, y=Freq, col=cat)) + 
  geom_point(aes(colour=cat))+
  scale_y_continuous("Frequency", trans = "log10",labels = trans_format("log10", math_format(10^.x)))+
  scale_x_continuous("Degree", trans = "log10",labels = trans_format("log10", math_format(10^.x)))+
  scale_color_manual(values = c("gray75", "gray45", "gray10"))+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major=element_blank(),panel.grid.minor=element_blank())

ggplot(degC, aes(x=degree, y=Freq, col=cat)) + 
  geom_point(aes(colour=cat))+
  scale_y_continuous("Frequency", trans = "log10",labels = trans_format("log10", math_format(10^.x)))+
  scale_x_continuous("Degree", trans = "log10",labels = trans_format("log10", math_format(10^.x)))+
  scale_color_manual(values = c("gray75", "gray45", "gray10"))+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major=element_blank(),panel.grid.minor=element_blank())

ggplot(degD, aes(x=degree, y=Freq, col=cat)) + 
  geom_point(aes(colour=cat))+
  scale_y_continuous("Frequency", trans = "log10",labels = trans_format("log10", math_format(10^.x)))+
  scale_x_continuous("Degree", trans = "log10",labels = trans_format("log10", math_format(10^.x)))+
  scale_color_manual(values = c("gray75", "gray45", "gray10"))+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major=element_blank(),panel.grid.minor=element_blank())

ggplot(degE, aes(x=degree, y=Freq, col=cat)) + 
  geom_point(aes(colour=cat))+
  scale_y_continuous("Frequency", trans = "log10",labels = trans_format("log10", math_format(10^.x)))+
  scale_x_continuous("Degree", trans = "log10",labels = trans_format("log10", math_format(10^.x)))+
  scale_color_manual(values = c("gray75", "gray45", "gray10"))+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major=element_blank(),panel.grid.minor=element_blank())

ggplot(degF, aes(x=degree, y=Freq, col=cat)) + 
  geom_point(aes(colour=cat))+
  scale_y_continuous("Frequency", trans = "log10",labels = trans_format("log10", math_format(10^.x)))+
  scale_x_continuous("Degree", trans = "log10",labels = trans_format("log10", math_format(10^.x)))+
  scale_color_manual(values = c("gray75", "gray45", "gray10"))+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major=element_blank(),panel.grid.minor=element_blank())

ggplot(degG, aes(x=degree, y=Freq, col=cat)) + 
  geom_point(aes(colour=cat))+
  scale_y_continuous("Frequency", trans = "log10",labels = trans_format("log10", math_format(10^.x)))+
  scale_x_continuous("Degree", trans = "log10",labels = trans_format("log10", math_format(10^.x)))+
  scale_color_manual(values = c("gray75", "gray45", "gray10"))+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major=element_blank(),panel.grid.minor=element_blank())

ggplot(degH, aes(x=degree, y=Freq, col=cat)) + 
  geom_point(aes(colour=cat))+
  scale_y_continuous("Frequency", trans = "log10",labels = trans_format("log10", math_format(10^.x)))+
  scale_x_continuous("Degree", trans = "log10",labels = trans_format("log10", math_format(10^.x)))+
  scale_color_manual(values = c("gray75", "gray45", "gray10"))+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major=element_blank(),panel.grid.minor=element_blank())


#Figure 4a
ggplot(diam,aes(x=Year,y=Index,color=Section))+geom_xspline(aes(linetype=Section),size=.9)+geom_point(aes(shape=Section),size=2)+theme_bw()+theme(axis.title=element_text(size=16,face="bold"),axis.text = element_text(size=14),legend.text =element_text(size=14))+scale_x_continuous("Year", 2007:2019, 2007:2019)+labs(y="Diameter")+scale_colour_grey(name="IPC Section")+ scale_shape_manual("IPC Section",values=c(15, 16, 17,18,19,3,4,8))+geom_dl(aes(label = Section), method = list("last.bumpup",cex=1.1,hjust=-1,vjust=-.1), color="black")+scale_linetype_manual("IPC Section",values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "F1", "12345678"))+ scale_y_continuous(breaks=seq(0,80,5))

# Figure 4b
ggplot(tcr.diam1,aes(x=Year,y=Index,color=Section))+geom_xspline(aes(linetype=Section),size=.9)+geom_point(aes(shape=Section),size=2)+theme_bw()+theme(axis.title=element_text(size=16,face="bold"),axis.text = element_text(size=14),legend.text =element_text(size=14))+scale_x_continuous("Year", 2008:2019, 2008:2019)+labs(y="Diameter rate of relative growth")+scale_colour_grey(name="IPC Section")+ scale_shape_manual("IPC Section",values=c(15, 16, 17,18,19,3,4,8))+geom_dl(aes(label = Section), method = list("last.bumpup",cex=1.1,hjust=-1,vjust=-.1), color="black")+scale_linetype_manual("IPC Section",values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "F1", "12345678"))#+ scale_y_continuous(breaks=seq(0,80,5))

# Figure 5

ggplot(neL, aes(x=node, y=edge)) + 
  geom_line(size=1.5)+
  scale_y_continuous("Number of edges", trans = "log10")+
  scale_x_continuous("Number of nodes", trans = "log10")+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major=element_blank(),panel.grid.minor=element_blank())

ggplot(neA, aes(x=node, y=edge)) + 
  geom_line(size=1.5)+
  scale_y_continuous("Number of edges", trans = "log10")+
  scale_x_continuous("Number of nodes", trans = "log10")+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major=element_blank(),panel.grid.minor=element_blank())

ggplot(neB, aes(x=node, y=edge)) + 
  geom_line(size=1.5)+
  scale_y_continuous("Number of edges", trans = "log10")+
  scale_x_continuous("Number of nodes", trans = "log10")+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major=element_blank(),panel.grid.minor=element_blank())

ggplot(neC, aes(x=node, y=edge)) + 
  geom_line(size=1.5)+
  scale_y_continuous("Number of edges", trans = "log10")+
  scale_x_continuous("Number of nodes", trans = "log10")+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major=element_blank(),panel.grid.minor=element_blank())

ggplot(neD, aes(x=node, y=edge)) + 
  geom_line(size=1.5)+
  scale_y_continuous("Number of edges", trans = "log10")+
  scale_x_continuous("Number of nodes", trans = "log10")+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major=element_blank(),panel.grid.minor=element_blank())

ggplot(neE, aes(x=node, y=edge)) + 
  geom_line(size=1.5)+
  scale_y_continuous("Number of edges", trans = "log10")+
  scale_x_continuous("Number of nodes", trans = "log10")+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major=element_blank(),panel.grid.minor=element_blank())

ggplot(neF, aes(x=node, y=edge)) + 
  geom_line(size=1.5)+
  scale_y_continuous("Number of edges", trans = "log10")+
  scale_x_continuous("Number of nodes", trans = "log10")+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major=element_blank(),panel.grid.minor=element_blank())

ggplot(neG, aes(x=node, y=edge)) + 
  geom_line(size=1.5)+
  scale_y_continuous("Number of edges", trans = "log10")+
  scale_x_continuous("Number of nodes", trans = "log10")+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major=element_blank(),panel.grid.minor=element_blank())

ggplot(neH, aes(x=node, y=edge)) + 
  geom_line(size=1.5)+
  scale_y_continuous("Number of edges", trans = "log10")+
  scale_x_continuous("Number of nodes", trans = "log10")+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major=element_blank(),panel.grid.minor=element_blank())


# Figure 6
p<-ggplot(gR, aes(x=year, y=size, col=gp,group=gp)) + 
  geom_line(aes(linetype=gp,colour=gp))+
  scale_y_continuous("CC Size", trans = "log10")+
  scale_x_continuous("Year",breaks = seq(from = 2007, to = 2019, by = 1))+
  scale_linetype_manual(values=c("solid", "dotdash", "solid", "solid"))+
  theme_bw()+scale_color_manual(values = c("gray25", "gray25", "gray45","gray75"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

p + facet_wrap( ~ class, nrow = 3) + 
  theme(legend.title= element_blank(),legend.position = "bottom")

# Figure 7a
ggplot(prais,aes(x=Year,y=Index,color=Section))+geom_xspline(aes(linetype=Section),size=.9)+geom_point(aes(shape=Section),size=2)+theme_bw()+theme(axis.title=element_text(size=16,face="bold"),axis.text = element_text(size=14),legend.text =element_text(size=14))+scale_x_continuous("Year", 2007:2019, 2007:2019)+labs(y="Diversificaction (Pr)")+ scale_y_continuous(breaks=seq(0,1,.05))+scale_colour_grey(name="IPC Section")+ scale_shape_manual("IPC Section",values=c(15, 16, 17,18,19,3,4,8))+geom_dl(aes(label = Section), method = list("last.bumpup",cex=1.1,hjust=-1,vjust=-.1), color="black")+scale_linetype_manual("IPC Section",values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "F1", "12345678"))

# Figure 7b
ggplot(tcr.div1,aes(x=Year,y=Index,color=Section))+geom_xspline(aes(linetype=Section),size=.9)+geom_point(aes(shape=Section),size=2)+theme_bw()+theme(axis.title=element_text(size=16,face="bold"),axis.text = element_text(size=14),legend.text =element_text(size=14))+scale_x_continuous("Year", 2008:2019, 2008:2019)+labs(y="Diversificaction rate of relative growth")+scale_colour_grey(name="IPC Section")+ scale_shape_manual("IPC Section",values=c(15, 16, 17,18,19,3,4,8))+geom_dl(aes(label = Section), method = list("last.bumpup",cex=1.1,hjust=-1,vjust=-.1), color="black")+scale_linetype_manual("IPC Section",values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "F1", "12345678"))#+ scale_y_continuous(breaks=seq(0,80,5))

#Figure 8

data2 <- diam1[,1:3]
colnames(data2)[3] <- "Diameter"
data2["Diversification"] <- prais1[,3]
ggplot(data2, aes(x=Diversification, y=Diameter,label=Year))+ geom_point(size=2)+geom_text_repel(aes(label=Year),hjust=0, vjust=0,force=1,color="black")+geom_hline(yintercept=1,linetype="dashed")+geom_vline(xintercept=1,linetype="dashed")+theme_bw()+facet_wrap(.~Sector)+xlab("Diversification Index (base 2007=1.0)")+ylab("Diameter Index (base 2007=1.0)")+scale_y_continuous(breaks=seq(0,7,1))+scale_x_continuous(breaks=seq(0,2.5,.25))









