## ----init1, echo=T, warning=F, message=FALSE, cache=F-------------------------
#load the R packages we need. 
library("pander")
library("coin")
library("knitr")
library("rmarkdown")
library("summarytools")
library("parallel")
library("DescTools")
library("ggplot2")
library("rcompanion")
library("ggpubr")
library("ggforce")

source("statfuns.r")

options(stringsAsFactors=F)
`%nin%` = Negate(`%in%`)

#myfile = "Mix_and_Match_master_20220428.csv"
#myfile = "Mix_and_Match_master_20220503.csv"
myfile = "Mix_and_Match_master_20220510.csv"


#load the booster dataset
x = read.csv(myfile)

#clean up the header names
h1 = read.csv(myfile, header=F, nrow=1)
#h1 = make.names(gsub(h1[1,], pat="\\(.*\\)", rep=""))
h1 = make.names(gsub(h1, pat="\\.\\.", rep="."))
h1 = make.names(gsub(h1, pat="\\.$", rep=""))
colnames(x) = make.unique(h1)

#trim whitespace and set to lowercase from every character
for (v in colnames(x)){
    if (is.character(x[,v])){
        x[,v] = tolower(trimws(x[,v]))
        #get rid of wierd white spaces
        x[,v] = trimws(x[, v], whitespace = "[\\h\\v]")
    }
}

#exclude patients from first column
table(x[,1] == 1)
x = x[is.na(x[,1]),]

#exclude any row without an MRN
x = x[!is.na(x$MRN),]



#outcome variables

#responders vs non-responders (column CM)
table(x$Outcome)
#compare change in spike ab before and after 4th dose (column BI and CL)
x$spike.baseline = as.numeric(gsub(x$Baseline.Study.Spike.Ab.Result..AU.mL, pat="<", rep=""))
x$spike.4 = as.numeric(gsub(x$X4.week.Spike.Ab.Result..AU.mL, pat="<", rep=""))
x$spike.diff = x$spike.4 - x$spike.baseline


#compare change in spike ab titer before and after 4th dose in responder vs non-responder (column BI and CM)
#1:pos, 0:neg
x$spike.res.baseline = as.numeric(as.factor(x[,"Baseline.Study.Spike.Ab.Result..Positive.or.Negative"])) - 1
x$spike.res.4 = as.numeric(as.factor(x[,"X4.week.Spike.Antibody.Result..Positive.or.Negative"])) - 1
x$spike.res.diff = x$spike.res.4 - x$spike.res.baseline
table(x[,"Baseline.Study.Spike.Ab.Result..Positive.or.Negative"])
table(x$spike.res.baseline)
table(x[,"X4.week.Spike.Antibody.Result..Positive.or.Negative"])
table(x$spike.res.4)
table(x$spike.res.diff)


#compare change in T-cell response (value) before and after 4 th dose (column BK and DJ)
x$tcell.baseline = as.numeric(gsub(x[,"Baseline.T.cell.assay.result.mIU.mL" ], pat=">", rep=""))
x$tcell.4 = as.numeric(gsub(x[,"X4.week.T.cell.assay.value.mIU.mL" ], pat=">", rep=""))
x$tcell.diff = x$tcell.4-x$tcell.baseline


#compare change in T-cell response (pos/neg) before and after 4th dose (column BJ and DI)
#1:pos, 0:neg
x$tcell.res.baseline = as.numeric(factor(x[,"Baseline.T.cell.Assay.Result"], levels=c("neg", "pos")))-1
x$tcell.res.4 = as.numeric(factor(x[,"X4.week.T.cell.assay.result"], levels=c("neg", "pos")))-1
x$tcell.res.diff = x$tcell.res.4 - x$tcell.res.baseline
table(x[,"Baseline.T.cell.Assay.Result"])
table(x$tcell.res.baseline)
table(x[,"X4.week.T.cell.assay.result"])
table(x$tcell.res.4)
table(x$tcell.res.diff) 


#dichotomize variables at 180
x[,"Time.since.last.cd20.and.4th.dose.over.180"] = x[,"Time.since.last.cd20.and.4th.dose",] > 180
x[,"Time.since.last.btki.and.4th.dose.over.180"] = x[,"Time.since.last.btki.and.4th.dose"] > 180
x[,"Time.since.CAR.T.over.180"] = x[,"Time.since.CAR.T"] > 180
x[,"Time.since.transplant.over.180"] = x[,"Time.since.transplant"] > 180




#numerical variables
num.vars = c(
    "Age",
    "Time.since.last.cd20.and.4th.dose",
    "Time.since.last.btki.and.4th.dose",
    "Time.since.CAR.T",
    "Time.since.transplant",
    "Baseline.total.WBC",
    "Baseline.ANC",
    "Baseline.ALC",
    "Baseline.AMC",
    "CD3.Abs",
    "CD4.abs",
    "CD8.abs",
    "CD.19.abs",
    "CD.16.56.abs",
    "Baseline.IgG",
    "Baseline.IGA",
    "Baseline.IgM",
    "X4.week.total.WBC",
    "X4.week.ANC",
    "X4.week.ALC",
    "X4.week.AMC",
    "CD3.Abs.1",
    "CD4.abs.1",
    "CD8.abs.1",
    "CD.19.abs.1",
    "CD.16.56.abs.1",
    "X4.wk.IgG",
    "X4.wk.IGA",
    "X4.wk.IgM"
)

#categorical
cat.vars = c(
    "Previous.mRNA.Vaccine.Given",
    "chemo",
    "cd20",
    "BTKi",
    "imid",
    "ADC",
    "surgery",
    "SCT",
    "XRT",
    "CAR.T",
    "hormone",
    "PI3Ki",
    "supportive",
    "immunotherapy",
    "proteasome",
    "Trial",
    "Steroids",
    "none",
    "BCL2",
    "Surgery",
    "Time.since.last.cd20.and.4th.dose.over.180",
    "Time.since.last.btki.and.4th.dose.over.180",
    "Time.since.CAR.T.over.180",
    "Time.since.transplant.over.180"
)

#summarize a variable include # of NAs
mysummary <- function(x){
    c(summary(x), n = sum(!is.na(x)), na=sum(is.na(x)))
}




## ----results="asis"-----------------------------------------------------------


#run a quick summary for each variable
for (v in num.vars){
    cat("##", v, "\n")
    cat("\n")
    cat(pander(mysummary(x[,v])))
    cat("\n")
}

for (v in cat.vars){
    cat("##", v, "\n")
    cat("\n")
    cat(pander(table(x[,v], useNA="ifany")))
    cat("\n")
}





## ---- warning=F---------------------------------------------------------------

#run kruskal wallis on a numeric vs categorical variable
#var1 the name of the numeric var
#var2 the name of the categorical var
runDiffAssociationTest <- function(x, var1, var2){
    #skip everything if its the same variable
    if (var1 == var2){
        return(F)
    }

    cat(paste0("\n###", var1, " by ", var2, "\n"))
    #make sure vars are factors so levels aren't dropped
    x2 = x[,c(var1, var2)]
    x2[,var2] = factor(x2[,var2])

    vals = x2[,var1]
    group = x2[,var2]
    try({
        res = kruskal.test(vals~group)
        if (res$p.value < .1){

            cat(paste0("#### * \n"))
        }
        cat(pander(kruskal.test(vals~group)))
        cat(pander(by(vals, group, FUN=mysummary)))
    }, silent=T)
}

#simple fisher test of proportions
runPropTest <- function(x, var1, var2){
    #skip everything if its the same variable
    if (var1 == var2){
        return (F)
    }


    cat(paste0("###", var1, " by ", var2, "\n"))
    #make sure vars are factors so levels aren't dropped
    x2 = x[,c(var1, var2)]
    x2[,var1] = factor(x2[,var1])
    x2[,var2] = factor(x2[,var2])

    #split into 2x2 table
    ctab1 = ctable(x2[,var1], x2[,var2], dnn=c(var1, var2), useNA="no")
    cat(pander(ctab1[[1]]))

    tab1 = table(x2[,var1], x2[,var2], dnn=c(var1, var2))

    try({
        res = fisher.test(tab1)
        if (res$p.value < .1){

            cat(paste0("#### * \n"))
        }
        #run fisher's exact test 
        cat(pander(fisher.test(tab1)))
        cat(paste0("\n"))
        knitr::normal_print(fisher.test(tab1))
        cat(paste0("\n"))
    }, silent=T)
}


#look for correclation between two variables
runCorTest <- function(x, var1, var2){

    cat(paste0("###", var1, " by ", var2, "\n"))
    
    vals1 = x[,var1]
    vals2 = x[,var2]

    #simple scatter plot
    df1 = data.frame(vals1, vals2)
    #p1 = ggplot(df1, aes(x = vals1, y = vals2)) + geom_point() +
    #    xlab(var1) + ylab(var2)
    #print(p1)
    #cat(paste0("\n"))

    cat(pander(cor.test(vals1, vals2, method="kendall")))

    #if significant show a star
    try({
        res = cor.test(vals1, vals2, method="kendall")
            if (res$p.value < .10){
                cat(paste0("#### *\n"))
    }}, silent=T)
    cat(paste0("\n"))
}



plotTitersByVar <- function(t1, t2, group, ylim=NA, timepoint = "4w", ylab="titers", alpha=.50){
    df1 = data.frame(rbind(
        data.frame(time=0, samp = 1:length(t1), titer=t1, group=group),
        data.frame(time=timepoint, samp = 1:length(t1), titer=t2, group=group)
    ))
    df1$time = as.factor(df1$time) 

    p1 = ggplot(df1, aes(x=time, y=titer, group=samp, color=group)) + geom_jitter(width=.05, size=.75) + geom_line(alpha=alpha) + ylab(ylab) + theme_classic()
    if (!is.na(ylim)){
        p1 = p1 + ylim(ylim)
    }

    #calculate CIs
    #https://www.r-graph-gallery.com/4-barplot-with-error-bar.html
    ##add mean trends for each level
    df2 = lapply(unique(group), function(g){
        alpha=0.05
        t=qt((1-alpha)/2 + .5, sum(group==g)-1)   # tend to 1.96 if sample size is big enough
        rbind(
            data.frame(time=0, samp = g, titer=mean(t1[group==g], na.rm=T), titer.ci = sd(t1[group==g], na.rm=T)/sqrt(sum(group==g))*t, group=g),
            data.frame(time=timepoint, samp = g, titer=mean(t2[group==g], na.rm=T), titer.ci = sd(t2[group==g], na.rm=T)/sqrt(sum(group==g))*t, group=g)
        )
    })
    df2 = do.call(rbind, df2)
    df2$time = as.factor(df2$time) 
    p2 = p1 + geom_line(alpha=1, data=df2, size=1.5) 
    p3 = p2 + geom_errorbar(aes(ymin = titer - titer.ci, ymax = titer + titer.ci), data=df2, size=.75, width=.15)
    print(p2)
    print(p3)


    df3 = data.frame(rbind(
        data.frame(time=0, titer=t1, group=group),
        data.frame(time=timepoint, titer=t2, group=group)
    ))

    #b1 = ggplot(df3, aes(x=time, y=titer, fill=group)) +  geom_boxplot(outlier.shape=NA) +  geom_point(position=position_jitterdodge(jitter.width=.1), color="black") + ylab(ylab) + theme_classic()

    b2 = ggplot(df3, aes(x=time, y=titer, fill=group)) +
        stat_summary(position=position_dodge(), fun.data ="mean_cl_boot", geom="crossbar")  + 
        geom_point(position=position_jitterdodge(jitter.width=.1), color="black") + 
        ylab(ylab) + 
        theme_classic()

    print(b2)


}






## ---- results="asis"----------------------------------------------------------


for (v in cat.vars){
    runPropTest(x, "Outcome", v)
}
for (v in num.vars){
    runDiffAssociationTest(x, v, "Outcome")
}


## -----------------------------------------------------------------------------



nullgroup = x$spike.baseline*0
plotTitersByVar(x$spike.baseline, x$spike.4, nullgroup)
wilcox.test(x$spike.diff)
pander(mysummary(x$spike.baseline))
pander(mysummary(x$spike.4))



## ---- results="asis", warning=F-----------------------------------------------


for (v in cat.vars){

    runDiffAssociationTest(x, "spike.diff", v)
    cat("\n")
    plotTitersByVar(x$spike.baseline, x$spike.4, as.factor(x[,v]))
    cat("\n")
}
for (v in num.vars){
    runCorTest(x, v, "spike.diff")
}


## -----------------------------------------------------------------------------

tab.spikeres = ctable(x$spike.res.baseline, x$spike.res.4, dnn=c("spike.base", "spike.4w"), useNA="no")
print(tab.spikeres)
pander(tab.spikeres)


## ---- results="asis"----------------------------------------------------------


for (v in cat.vars){
    runPropTest(x, "spike.res.diff", v)
}
for (v in num.vars){
    runDiffAssociationTest(x, v, "spike.res.diff")
}


## -----------------------------------------------------------------------------
wilcox.test(x$tcell.diff)
pander(mysummary(x$tcell.baseline))
pander(mysummary(x$tcell.4))



## ---- results="asis", warning=F-----------------------------------------------


for (v in cat.vars){
    runDiffAssociationTest(x, "tcell.diff", v)
    cat("\n")
    plotTitersByVar(x$tcell.baseline, x$tcell.4, as.factor(x[,v]))
    cat("\n")
}
for (v in num.vars){
    runCorTest(x, v, "tcell.diff")
}


## -----------------------------------------------------------------------------

tab.tcellres = ctable(x$tcell.res.baseline, x$tcell.res.4, dnn=c("tcell.base", "tcell.4w"), useNA="no")
print(tab.tcellres)
pander(tab.tcellres)


## ---- results="asis"----------------------------------------------------------


for (v in cat.vars){
    runPropTest(x, "tcell.res.diff", v)
}
for (v in num.vars){
    runDiffAssociationTest(x, v, "tcell.res.diff")
}


## -----------------------------------------------------------------------------

asis_output("##Outcome vs Baseline Igm\n")

ggplot(x, aes(x=Outcome, y=Baseline.IgM, fill=Outcome)) + 
    geom_boxplot(outlier.shape = NA) + 
    coord_cartesian(ylim = quantile(x$Baseline.IgM, c(0.05, 0.95), na.rm=T)) +
    geom_jitter(width = 0.2)

#ggplot(x, aes(x=Outcome, y=Baseline.IgM, color=Outcome)) + geom_point()
#ggplot(x, aes(x=Outcome, y=Baseline.IgM, color=Outcome)) + geom_point() + ylim(c(0,250))

asis_output("##Outcome vs 4 week Igm\n")
ggplot(x, aes(x=Outcome, y=X4.wk.IgM, fill=Outcome)) + 
    geom_boxplot(outlier.shape = NA) + 
    coord_cartesian(ylim = quantile(x$X4.wk.IgM, c(0.05, 0.95), na.rm=T)) +
    geom_jitter(width = 0.2)

#ggplot(x, aes(x=Outcome, y=X4.wk.IgM, color=Outcome)) + geom_point() 
#ggplot(x, aes(x=Outcome, y=X4.wk.IgM, color=Outcome)) + geom_point() + ylim(c(0, 250))

asis_output("##Outcome vs Baseline CD19\n")
ggplot(x, aes(x=Outcome, y=CD.19.abs, fill=Outcome)) + 
    geom_boxplot(outlier.shape = NA) + 
    coord_cartesian(ylim = quantile(x$CD.19.abs, c(0.05, 0.95), na.rm=T)) +
    geom_jitter(width = 0.2)

#ggplot(x, aes(x=Outcome, y=CD.19.abs, color=Outcome)) + geom_point() 
#ggplot(x, aes(x=Outcome, y=CD.19.abs, color=Outcome)) + geom_point() + ylim(c(0,3000))

asis_output("##Outcome vs 4 week CD19\n")
ggplot(x, aes(x=Outcome, y=CD.19.abs.1, fill=Outcome)) + 
    geom_boxplot(outlier.shape = NA) + 
    coord_cartesian(ylim = quantile(x$CD.19.abs.1, c(0.05, 0.95), na.rm=T)) +
    geom_jitter(width = 0.2)

#ggplot(x, aes(x=Outcome, y=CD.19.abs.1, color=Outcome)) + geom_point() 
#ggplot(x, aes(x=Outcome, y=CD.19.abs.1, color=Outcome)) + geom_point() + ylim(c(0,3000))




## -----------------------------------------------------------------------------
x[,"BASELINE.WT.ID50"] = as.numeric(x[,"BASELINE.WT.ID50"]) 
x[,"OMICRON.B.1.1.529.ID50"] = as.numeric(x[,"OMICRON.B.1.1.529.ID50"])
x[,"X4.WEEK.WT.ID50"] = as.numeric(x[,"X4.WEEK.WT.ID50"])
x[,"X4.WEEK.B.1.1.529.ID50"] = as.numeric(x[,"X4.WEEK.B.1.1.529.ID50"])


asis_output("##spike.baseline vs BASELINE.WT.ID50\n")
cor.test(x$spike.baseline, x$BASELINE.WT.ID50, method="spearman", use="pairwise.complete.obs")
#plot(x$spike.baseline, x$BASELINE.WT.ID50)
#ggplot(x, aes(x=spike.baseline, y=BASELINE.WT.ID50)) + geom_point() + geom_smooth(method=lm, se=T)
ggscatter(x, x = "spike.baseline", y = "BASELINE.WT.ID50",
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
    ) + stat_cor(method = "spearman")


asis_output("##spike.baseline vs OMICRON.B.1.1.529.ID50\n")
cor.test(x$spike.baseline, x$OMICRON.B.1.1.529.ID50, method="spearman", use="pairwise.complete.obs")
#plot(x$spike.baseline, x$OMICRON.B.1.1.529.ID50)
#ggplot(x, aes(x=spike.baseline, y=OMICRON.B.1.1.529.ID50)) + geom_point() + geom_smooth(method=lm, se=T)
ggscatter(x, x = "spike.baseline", y = "OMICRON.B.1.1.529.ID50",
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
    ) + stat_cor(method = "spearman")



asis_output("##spike.4week vs X4.WEEK.WT.ID50\n")
cor.test(x$spike.4, x$X4.WEEK.WT.ID50, method="spearman", use="pairwise.complete.obs")
#plot(x$spike.4, x$X4.WEEK.WT.ID50)
#ggplot(x, aes(x=spike.4, y=X4.WEEK.WT.ID50)) + geom_point() + geom_smooth(method=lm, se=T)
ggscatter(x, x = "spike.4", y = "X4.WEEK.WT.ID50",
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
    ) + stat_cor(method = "spearman")


asis_output("##spike.4week vs X4.WEEK.B.1.1.529.ID50\n")
cor.test(x$spike.4, x$X4.WEEK.B.1.1.529.ID50, method="spearman", use="pairwise.complete.obs")
#plot(x$spike.4, x$X4.WEEK.B.1.1.529.ID50)
#ggplot(x, aes(x=spike.4, y=X4.WEEK.B.1.1.529.ID50)) + geom_point() + geom_smooth(method=lm, se=T)
ggscatter(x, x = "spike.4", y = "X4.WEEK.B.1.1.529.ID50",
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
    ) + stat_cor(method = "spearman")




asis_output("##combined  spike vs WT.ID50\n")
df1 = rbind(
    data.frame(ab.spike = x$spike.baseline, wt.id50 = x[,"BASELINE.WT.ID50"], time="baseline"),
    data.frame(ab.spike = x$spike.4, wt.id50 = x[,"X4.WEEK.WT.ID50"], time="four week")
)

ggplot(df1, aes(x=ab.spike, y=wt.id50, color=time)) +
    geom_point() + 
    geom_smooth(method = "lm") +
    stat_cor(method="spearman", label.x.npc = 0.1, show.legend = FALSE) +
    #coord_cartesian(xlim =c(0, 1000), ylim= c(0,200)) +
    facet_zoom(
        xlim =c(50, 1000),
        ylim= c(0,200), 
        horizontal=F,
        shrink=F,
        zoom.size=.5
    )+
    #theme_classic()
    theme_bw() +
    #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme()


asis_output("##combined  spike vs 529.ID50\n")

df2 = rbind(
    data.frame(ab.spike = x$spike.baseline, o529.id50 = x[,"OMICRON.B.1.1.529.ID50"], time="baseline"),
    data.frame(ab.spike = x$spike.4, o529.id50 = x[,"X4.WEEK.B.1.1.529.ID50"], time="four week")
)

ggplot(df2, aes(x=ab.spike, y=o529.id50, color=time)) +
    geom_point() + 
    geom_smooth(method = "lm") +
    stat_cor(method="spearman", label.x.npc = 0.1, show.legend = FALSE) +
    facet_zoom(
        xlim =c(50, 1000),
        ylim= c(0,10), 
        horizontal=F,
        shrink=F,
        zoom.size=.5
    )+
    #theme_classic()
    theme_bw() +
    #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme()





## -----------------------------------------------------------------------------


asis_output("##spike.baseline vs tcell baseline\n")
cor.test(x$spike.baseline, x$tcell.baseline, method="spearman", use="pairwise.complete.obs")
#plot(x$spike.baseline, x$BASELINE.WT.ID50)
#ggplot(x, aes(x=spike.baseline, y=BASELINE.WT.ID50)) + geom_point() + geom_smooth(method=lm, se=T)
ggscatter(x, x = "spike.baseline", y = "tcell.baseline",
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
    ) + stat_cor(method = "spearman")


asis_output("##spike.baseline.4w vs tcell 4w\n")
cor.test(x$spike.4, x$tcell.4, method="spearman", use="pairwise.complete.obs")
#plot(x$spike.baseline, x$BASELINE.WT.ID50)
#ggplot(x, aes(x=spike.baseline, y=BASELINE.WT.ID50)) + geom_point() + geom_smooth(method=lm, se=T)
ggscatter(x, x = "spike.4", y = "tcell.4",
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
    ) + stat_cor(method = "spearman")


asis_output("##combined spike.baseline vs tcell\n")


df3 = rbind(
    data.frame(ab.spike = x$spike.baseline, tcell = x[,"tcell.baseline"], time="baseline"),
    data.frame(ab.spike = x$spike.4, tcell = x[,"tcell.4"], time="four week")
)

ggplot(df3, aes(x=ab.spike, y=tcell, color=time)) +
    geom_point() + 
    geom_smooth(method = "lm") +
    stat_cor(method="spearman", label.x.npc = 0.1, show.legend = FALSE) +
    #coord_cartesian(xlim =c(0, 1000), ylim= c(0,200)) +
    facet_zoom(
        xlim =c(50, 1000),
        ylim= c(0,2000), 
        horizontal=F,
        shrink=T,
        show.area=T,
        zoom.size=.5
    )+
    #theme_classic()
    theme_bw() +
    #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme()

ggplot(df3, aes(x=ab.spike, y=tcell, color=time)) +
    geom_point() + 
    geom_smooth(method = "lm") +
    stat_cor(method="spearman", label.x.npc = 0.1, show.legend = FALSE) +
    #coord_cartesian(xlim =c(0, 1000), ylim= c(0,200)) +
    facet_zoom(
        xlim =c(50, 1000),
        ylim= c(0,2000), 
        horizontal=F,
        shrink=T,
        show.area=F,
        zoom.size=.5
    )+
    #theme_classic()
    theme_bw() +
    #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme()




