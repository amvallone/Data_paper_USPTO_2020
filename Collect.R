
setwd("")
library(patentsview)

pull.data <-function(year,class){
  m.querry <- function(year){
    y <- year + 1
    m1 <- qry_funs$and(qry_funs$gte(patent_date = paste(year, "-01-01",sep="")),qry_funs$lt(patent_date = paste(year, "-02-01",sep="")))
    m2 <- qry_funs$and(qry_funs$gte(patent_date = paste(year, "-02-01",sep="")),qry_funs$lt(patent_date = paste(year, "-03-01",sep="")))
    m3 <- qry_funs$and(qry_funs$gte(patent_date = paste(year, "-03-01",sep="")),qry_funs$lt(patent_date = paste(year, "-04-01",sep="")))
    m4 <- qry_funs$and(qry_funs$gte(patent_date = paste(year, "-04-01",sep="")),qry_funs$lt(patent_date = paste(year, "-05-01",sep="")))
    m5 <- qry_funs$and(qry_funs$gte(patent_date = paste(year, "-05-01",sep="")),qry_funs$lt(patent_date = paste(year, "-06-01",sep="")))
    m6 <- qry_funs$and(qry_funs$gte(patent_date = paste(year, "-06-01",sep="")),qry_funs$lt(patent_date = paste(year, "-07-01",sep="")))
    m7 <- qry_funs$and(qry_funs$gte(patent_date = paste(year, "-07-01",sep="")),qry_funs$lt(patent_date = paste(year, "-08-01",sep="")))
    m8 <- qry_funs$and(qry_funs$gte(patent_date = paste(year, "-08-01",sep="")),qry_funs$lt(patent_date = paste(year, "-09-01",sep="")))
    m9 <- qry_funs$and(qry_funs$gte(patent_date = paste(year, "-09-01",sep="")),qry_funs$lt(patent_date = paste(year, "-10-01",sep="")))
    m10 <- qry_funs$and(qry_funs$gte(patent_date = paste(year, "-10-01",sep="")),qry_funs$lt(patent_date = paste(year, "-11-01",sep="")))
    m11 <- qry_funs$and(qry_funs$gte(patent_date = paste(year, "-11-01",sep="")),qry_funs$lt(patent_date = paste(year, "-12-01",sep="")))
    m12 <- qry_funs$and(qry_funs$gte(patent_date = paste(year, "-12-01",sep="")),qry_funs$lt(patent_date = paste(y, "-01-01",sep="")))
    out <- list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12)
    out
  }	
  
  get <- m.querry(year)
  get.info <- function(class){
    data <- list()
    for (i in 1:length(get)){
      fields <- c("patent_id",get_fields("patents",c(class)))
      cat("it",i,"de",length(get),"\n")
      data[[i]] <- search_pv(get[[i]],  fields = fields, all_pages=TRUE)
      save(data, file = paste(year,"_",class,".rdata",sep=""))
    }
    invisible(data)
  }
  out <- get.info(class)
  invisible(out)
}


compile <- function(year){
  path <- paste("./",year,sep="")
  files <- dir(path)[grep(pattern="rdata",dir(path))]
  full.path <- paste(path,files,sep="/")
  full.data <- list()
  for (i in 1:length(full.path)){
    store <-list()
    load(full.path[i])
    for (j in 1:length(data)){
      z <-unnest_pv_data(data[[j]]$data)
      store[[j]] <-z[[1]]
    }
    full.data[[i]] <-do.call(rbind.data.frame,store)
  }
  out <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "patent_id", all.x = TRUE),full.data) 
  invisible(out)
}


inventor <- lapply(2007:2029,function(x) pull.data(x,"inventors"))
assigness <- lapply(2011:2019,function(x) pull.data(x,"assignees"))
patents <- lapply(2007:2019,function(x) pull.data(x,"patents"))
ipc  <- lapply(2007:2019,function(x) pull.data(x,"ipcs"))


t2007 <- compile("2007")
t2007<-cbind(t2007$patent_id,t2007$inventor_id,t2007$ipc_section,t2007$ipc_class,t2007$ipc_subclass,t2007$patent_year)

t2008 <- compile("2008")
t2008<-cbind(t2008$patent_id,t2008$inventor_id,t2008$ipc_section,t2008$ipc_class,t2008$ipc_subclass,t2008$patent_year)

t2009 <- compile("2009")
t2009<-cbind(t2008$patent_id,t2009$inventor_id,t2009$ipc_section,t2009$ipc_class,t2009$ipc_subclass,t2009$patent_year)

t2010 <- compile("2010")
t2010<-cbind(t2010$patent_id,t2010$inventor_id,t2010$ipc_section,t2010$ipc_class,t2010$ipc_subclass,t2010$patent_year)

t2011 <- compile("2011")
t2011<-cbind(t2011$patent_id,t2011$inventor_id,t2011$ipc_section,t2011$ipc_class,t2011$ipc_subclass,t2011$patent_year)

t2012 <- compile("2012")
t2012<-cbind(t2012$patent_id,t2012$inventor_id,t2012$ipc_section,t2012$ipc_class,t2012$ipc_subclass,t2012$patent_year)

t2013 <- compile("2013")
t2013<-cbind(t2013$patent_id,t2013$inventor_id,t2013$ipc_section,t2013$ipc_class,t2013$ipc_subclass,t2013$patent_year)

t2014 <- compile("2014")
t2014<-cbind(t2014$patent_id,t2014$inventor_id,t2014$ipc_section,t2014$ipc_class,t2014$ipc_subclass,t2014$patent_year)

t2015 <- compile("2015")
t2015<-cbind(t2015$patent_id,t2015$inventor_id,t2015$ipc_section,t2015$ipc_class,t2015$ipc_subclass,t2015$patent_year)

t2016 <- compile("2016")
t2016<-cbind(t2016$patent_id,t2016$inventor_id,t2016$ipc_section,t2016$ipc_class,t2016$ipc_subclass,t2016$patent_year)

t2017 <- compile("2017")
t2017<-cbind(t2017$patent_id,t2017$inventor_id,t2017$ipc_section,t2017$ipc_class,t2017$ipc_subclass,t2017$patent_year)

t2018 <- compile("2018")
t2018<-cbind(t2018$patent_id,t2018$inventor_id,t2018$ipc_section,t2018$ipc_class,t2018$ipc_subclass,t2018$patent_year)

t2019 <- compile("2019")
t2019<-cbind(t2019$patent_id,t2019$inventor_id,t2019$ipc_section,t2019$ipc_class,t2019$ipc_subclass,t2019$patent_year)



