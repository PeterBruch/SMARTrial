removeCorrelatedFeatures <- function(x, cutoff = 0.8, distance = "pearson", cluster_method = "ward.D2") {
    # calculate distance matrix
    if (distance == "binary") {
        #maybe also useful is the input is a sparse matrix
        distMat <- stats::dist(t(x), method = "binary")
    } else if (distance == "pearson") {
        #otherwise, using pearson correlation
        distMat <- stats::as.dist(1-cor(x))
    } else if (distance == "euclidean") {
        distMat <- stats::dist(t(x), method = "euclidean")
    } else if (distance == "cosine") {
        # cosine similarity maybe preferred for sparse matrix
        cosineSimi <- function(x){
            x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))
        }
        distMat <- stats::as.dist(1-cosineSimi(t(x)))
    } else if (distance == "canberra") {
        distMat <- stats::as.dist(as.matrix(dist(t(x), method = "canberra"))/nrow(x))
    }
    
    #hierarchical clustering
    hc <- stats::hclust(distMat, method = cluster_method)
    clusters <- stats::cutree(hc, h = 1-cutoff)
    x.re <- x[,!duplicated(clusters)]
    
    #record the removed features
    mapList <- lapply(colnames(x.re), function(i) {
        members <- names(clusters[clusters == clusters[i]])
        members[members != i]
    })
    names(mapList) <- colnames(x.re)
    
    return(list(reduced = x.re,
                mapReduce = mapList))
}


testPartition <- function(y, ratio) {
    #balanced sampling of test set
    ord <- seq_along(y)
    testIdx <- lapply(unique(y),function(n) {
        subY <- ord[y == n]
        sample(subY, size = as.integer(length(subY)  * ratio)) 
    }) %>% do.call(c,.) %>% sort()
    return(testIdx)
}

#Function for multi-variant binomial regression
runGlm.bin <- function(X, y, alpha=1, repeats=20, folds = 3, testRatio = NULL, lambda = "lambda.1se") {
    modelList <- list()
    lambdaList <- c()
    aucCV <- c()
    aucTest <- c()
    rocTest <- list()
    coefMat  <- matrix(NA, ncol(X), repeats)
    rownames(coefMat) <- colnames(X)
    
    for (i in seq(repeats)) {
        if (!is.null(testRatio)) {
            testIdx <- testPartition(y, testRatio)
            X.test <- X[testIdx,,drop=FALSE]
            X.train <- X[-testIdx,, drop=FALSE]
            y.test <- y[testIdx]
            y.train <- y[-testIdx]
        } else {
            X.train <- X
            y.train <- y
        }
        
        vecFold <- mltools::folds(y.train,nfolds = folds, stratified = TRUE)
        
        #train model
        res <- cv.glmnet(X.train,y.train, type.measure = "auc",
                         foldid = vecFold, alpha = alpha, standardize = FALSE,
                         intercept = TRUE, family = "binomial")
        lambdaList <- c(lambdaList, res[[lambda]])
        aucCV <- c(aucCV, res$cvm[res$lambda == res[[lambda]]])
        modelList[[i]] <- res
        coefMat[,i] <- coef(res, s = lambda)[-1]
        
        #test model if testRatio is speficied
        if(!is.null(testRatio)) {
            rocRes <- plotROC(res, X.test, y.test, lambda)
            aucTest <- c(aucTest, rocRes$auc)
            rocTest[[i]] <- rocRes$plot
        }
    }
    list(modelList = modelList, lambdaList = lambdaList, aucCV = aucCV, coefMat = coefMat,
         aucTest = aucTest, rocTest = rocTest)
}





#function to plot ROC curves
plotROC <- function(glmModel, X, y, lambda = "lambda.1se") {
    lambdaChose <- glmModel[[lambda]]
    glmPred <- prediction(predict(glmModel, type = "response", newx = X, s=lambdaChose), y)
    glmPerform <- performance(glmPred,"tpr","fpr")
    aucVal <- performance(glmPred, measure = "auc")@y.values[[1]]
    xname <- glmPerform@x.name
    yname <- glmPerform@y.name
    plotTab <- tibble(x= glmPerform@x.values[[1]],
                      y = glmPerform@y.values[[1]])
    p <- ggplot(plotTab, aes(x=x, y =y )) + geom_line(color = "red") +
        xlab(xname) + ylab(yname) + theme(panel.grid = element_blank())
    
    if (!is.null(aucVal)) {
        p <- p + annotate("text", x= 0.75, y = 0.25, label = sprintf("AUC: %1.2f", aucVal))
    }
    list(plot=p, auc = aucVal)
}


# Function to run random forest
runForestRun <- function(x, Y, n=1, testRatio = 0.2, ntreeTry=1000, ntree = 10000, plotROC = FALSE) {
    # x is the matrix for explanaotry vector and Y is for response vector, currently, y should be a catagorical variable
    # n is the number of time to run random forest
    # bin is (the number of test set samples) / (the number of training set samples)
    
    #recording the results
    mtryList <- c()
    aucList <- c()
    matMDA <- c()
    matMDG <- c()
    
    #combine dataset
    stopifnot(nrow(x) == length(Y))
    allTab <- data.frame(cbind(x,Y))
    allTab$Y <- as.factor(allTab$Y)
    
    for (i in seq(n)) {
        #seperate training set and test set
        
        train_id <- testPartition(y, 1-testRatio)
        #train_id <- sample(seq_len(nrow(allTab)),size = as.integer(nrow(allTab)*bin),replace = FALSE, prob = NULL)
        dataTrain <- allTab[train_id,]
        dataTest <- allTab[-train_id,]
        
        #choosing a best value for mtry 
        bestmtry <- tuneRF(dataTrain[,-ncol(dataTrain)],dataTrain$Y, ntreeTry=ntreeTry, stepFactor=1.5,improve=0.01, trace=FALSE, plot=FALSE, dobest=FALSE)
        n_mtry <- bestmtry[order(bestmtry[,2]),][1,1] #get the mtry value with least OOBError
        mtryList <- c(mtryList,n_mtry)
        
        #run random forest
        resRF <-randomForest(Y~.,data=dataTrain, mtry=n_mtry, ntree=ntree, keep.forest=TRUE, importance=TRUE, test=dataTest)
        
        #evaluate the model by AUC of ROC
        resRF.pr <- predict(resRF, type = "prob", newdata = dataTest)[,2] 
        resRF.pred <- prediction(resRF.pr, dataTest$Y)
        resRF.auc <- slot(performance(resRF.pred,"auc"),"y.values")[[1]]
        aucList <- c(aucList,resRF.auc)
        
        #plot ROC curve
        if (plotROC) {
            resRF.perf <- performance(resRF.pred,"tpr","fpr")
            plot(resRF.perf, main="ROC Curve for Random Forest", col=2, lwd=2)
            abline(a=0,b=1,col="gray",lwd=2,lty=2)
            legend("bottomright",sprintf("AUC = %1.2f",resRF.auc),bty = "n")
        }
        
        #get the importance value of each variable
        resImp <- importance(resRF)
        matMDA <- cbind(matMDA,resImp[,3])
        matMDG <- cbind(matMDG,resImp[,4])
        #colnames(matMDA) <- colnames(matMDG) <- seq(n)
    }
    res <- list(mtry = mtryList, auc = aucList, MDA = matMDA, MDG = matMDG)
}

#function for scaling predictors
dataScale <- function(x, censor = NULL, robust = FALSE) {
    #function to scale different variables
    if (length(unique(na.omit(x))) <=3){
        #a binary variable, change to -0.5 and 0.5 for 1 and 2
        x - 0.5
    } else {
        if (robust) {
            #continuous variable, centered by median and divied by 2*mad
            mScore <- (x-median(x,na.rm=TRUE))/mad(x,na.rm=TRUE)
            if (!is.null(censor)) {
                mScore[mScore > censor] <- censor
                mScore[mScore < -censor] <- -censor
            }
            mScore/2
        } else {
            mScore <- (x-mean(x,na.rm=TRUE))/(sd(x,na.rm=TRUE))
            if (!is.null(censor)) {
                mScore[mScore > censor] <- censor
                mScore[mScore < -censor] <- -censor
            }
            mScore/2
        }
    }
}
#functions to cmpare R2
compareR2 <- function(matList, y, testRatio = 0.3, repeats = 1000, model = "linear", boots=TRUE) {
    #matList is a list of feature matrix to be evaluated. Sample names in rows and feature names in columns
    overSample <- names(y)
    for (eachMat in matList) {
        overSample <- intersect(overSample, rownames(eachMat))
    }
    
    resp <- y[overSample]
    extList <- matList
    extList[["Combine"]] <- lapply(matList, function(x) x[overSample,,drop=FALSE]) %>% do.call(cbind,.)
    for (eachSet in names(extList)) {
        x <- data.frame(extList[[eachSet]])
        x <- x[overSample,,drop=FALSE]
        #x[["response"]] <- resp
        extList[[eachSet]] <- x
    }
    
    # for each dataset, build regression models with repeats
    allResTab <- lapply(names(extList), function(nn) {
        
        eachMat <- extList[[nn]]
        eachRes <- lapply(seq(repeats), function(i) {
            smpTrain <- sample(overSample, as.integer((1-testRatio)*length(overSample)), replace = TRUE)
            smpTest <- overSample[!overSample %in% smpTrain]
            
            yTrain <- y[smpTrain]
            yTest <- y[smpTest]
            xTrain <- eachMat[smpTrain,,drop=FALSE]
            xTest <- eachMat[smpTest,,drop=FALSE]
            dataSet <- xTrain
            dataSet[["y"]] <- yTrain
            if (model == "linear") {
                mm <- glm(y ~ ., data = dataSet, family = "gaussian")
                yPredTrain <- predict(mm)
                yPredTest <- predict(mm, xTest)
                pTrain <- cor(yPredTrain, yTrain)^2
                pTest <- cor(yPredTest, yTest)^2
            } else if (model == "logistic") {
                mm <- glm(y ~ ., data = dataSet, family = "binomial")
                yPredTrain <- predict(mm, type = "response")
                yPredTest <- predict(mm, xTest, type = "response")
                pTrain <- mean((yPredTrain>0.5) == yTrain)
                pTest <- mean((yPredTest>0.5) == yTest)
            }
            
            tibble(name = nn, term = c("train","test"), performance = c(pTrain, pTest), rep=i)
            
        }) %>% bind_rows()
        
    }) %>% bind_rows()
    
    sumRes <- group_by(allResTab, name, term) %>%
        nest() %>%
        mutate(CI = map(data, ~calcCI(.$performance, boots = boots)),
               meanPerform = map(data, ~mean(.$performance))) %>%
        unnest(c(CI,meanPerform)) %>% ungroup() %>%
        mutate(type = rep(c("lower","upper"), nrow(.)/2)) %>%
        select(-data) %>%
        pivot_wider(names_from = type, values_from = CI)
    
    xLab <- ifelse(model == "linear", "Variance explained", "Accuracy")
    p <- ggplot(sumRes, aes(x=name, y=meanPerform)) +
        geom_bar(aes(fill = name), stat="identity") +
        scale_fill_manual(values = c("#7fc97f","#beaed4", "#fdc086")) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width =0.3) +
        facet_wrap(~term) +
        ylab(xLab) + xlab("") +
        coord_flip(ylim = c(0.5,0.8)) +
        theme_bw() +
        theme(legend.position = "none")
    
    return(list(table = sumRes, plot = p))
}

#Function to clean and combine 
prepareInput <- function(y, mutMatrix, drugMatrix, onlyCombine = FALSE, censor = NULL, geneIsEln = FALSE) {
    
    #generate explainatory variable table
    expTab <- list()
    
    #change gene matrix name if it's eln
    geneName <- ifelse(geneIsEln, "ELN22", "gene")
    
    #cnv/snv matrix
    overSample1 <- intersect(names(y), rownames(mutMatrix))
    geneTab <- mutMatrix[overSample1,,drop=FALSE]
    #at least 3 mutated sample after subsetting
    if (ncol(geneTab) > 1) {
        geneTab <- geneTab[, colSums(geneTab) >= 3]
    }
    expTab[[geneName]] <- list(y = y[overSample1],
                             X = apply(geneTab,2,dataScale))
    
    #drug responses
    overSample2 <- intersect(rownames(drugMatrix), names(y))
    drugMatrix <- drugMatrix[overSample2,]
    expTab[["drug"]] <- list(y=y[overSample2],
                             X = apply(drugMatrix, 2, dataScale, censor))
    
    
    #combine three tables
    overSample <- intersect(overSample1, overSample2)
    
    expTab[["combine"]] <- list(y = y[overSample],
                                X = cbind(expTab[["drug"]][["X"]][overSample,,drop=FALSE], 
                                          expTab[[geneName]][["X"]][overSample,,drop=FALSE]))
    
    if (onlyCombine) {
        #only return combined results, for feature selection
        expTab <-expTab[length(expTab)]
    }
    
    return(expTab)
    
}

#Function for plotting variance explained for each measurement
plotVar <- function(glmResult, set = "train") {
    
    plotTab <- lapply(names(glmResult), function(x) {
        if (set == "train") {
            tibble(variable = x, value = glmResult[[x]]$aucCV)
        } else {
            tibble(variable = x, value = glmResult[[x]]$aucTest)
        }
        
    }) %>% bind_rows() %>% group_by(variable) %>% 
        summarise(mean=mean(value, na.rm = TRUE),sd=sd(value, na.rm=TRUE))
    
    
    g <- ggplot(plotTab,aes(x=variable, y=mean, fill= variable)) + 
        geom_bar(position=position_dodge(), stat="identity", width = 0.8, col="black") + 
        geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, width = 0.3), position=position_dodge(.9)) + 
        theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =12), 
                                plot.title = element_text(hjust =0.5),
                                axis.text.y = element_text(size=12),
                                axis.title = element_text(size=15),
                                legend.position = "none") + 
        scale_fill_brewer("Set1",type = "qual") + coord_cartesian(ylim = c(0,1)) +
        ylab("AUC") + xlab("") + ggtitle(paste0("cross-validated area under ROC"))
    
    return(list(plot = g, varTab = plotTab))
}


# For bootstraping
#function for using bootstraping for evaluating model performance
calcCI <- function(stats, alphaCI=0.95, boots = TRUE) {
    if (boots) {
        p <- ((1.0-alphaCI)/2.0) 
        lower <- as.numeric(quantile(stats, p))
        p <- (alphaCI+((1.0-alphaCI)/2.0))
        upper <- as.numeric(quantile(stats, p))
    } else {
        lower <- mean(stats) - 1.96*sd(stats)/sqrt(length(stats))
        upper <- mean(stats) + 1.96*sd(stats)/sqrt(length(stats))
    }
   
    return(c(lower=lower, upper=upper))
}

#function for run bootstrapped glmnet
bootstrapGlm.bin <- function(X, y, alpha=1, repeats=1000, folds = 3, testRatio = 0.3, lambda = "lambda.1se") {
    
    cPerform <- rep(NA, repeats)
    coefMat  <- matrix(NA, ncol(X), repeats)
    rownames(coefMat) <- colnames(X)
    
    #build evaluate model for each bootstrapped sample
    for (i in seq(repeats)) {
        
        trainIdx <- sample(seq_along(y), size = as.integer(length(y)*(1-testRatio)), replace = TRUE)
        X.train <- X[trainIdx,]
        X.test <- X[-trainIdx,]
        y.train <- y[trainIdx]
        y.test <- y[-trainIdx]
        
        vecFold <- mltools::folds(y.train,nfolds = folds, stratified = TRUE)
        
        #train model
        res <- cv.glmnet(X.train,y.train, type.measure = "class",
                         foldid = vecFold, alpha = alpha, standardize = FALSE,
                         intercept = TRUE, family = "binomial")
        
        coefMat[,i] <- coef(res, s = lambda)[-1]
        
        y.pred <- as.numeric(predict(res, newx = X.test, s = lambda, type = "class"))
        cPerform[i] <- mean(y.pred == y.test)
    }
    
    #calculate bootstraped statistics for model performance and feature coefficient. 
    list(coefMat = coefMat, performance = cPerform)
}

#function to plot bootstrapped CI for model performance
plotModelBoots <- function(glmResult, alphaCI = 0.95) {
    
    plotTab <- lapply(names(glmResult), function(x) {
        ci <- calcCI(glmResult[[x]]$performance, alphaCI = 0.95)
        tibble(variable = x, mean = mean(glmResult[[x]]$performance),
               lower = ci[1], upper = ci[2])
    }) %>% bind_rows() 
    
    g <- ggplot(plotTab,aes(x=variable, y=mean, fill= variable)) + 
        geom_bar(position=position_dodge(), stat="identity", width = 0.8, col="black") + 
        geom_errorbar(aes(ymin=lower, ymax=upper, width = 0.3), position=position_dodge(.9)) + 
        theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =12), 
                                plot.title = element_text(hjust =0.5),
                                axis.text.y = element_text(size=12),
                                axis.title = element_text(size=15),
                                legend.position = "none") + 
        scale_fill_brewer("Set1",type = "qual") + coord_cartesian(ylim = c(0,1)) +
        ylab("Accuracy") + xlab("")
    
    return(list(plot = g, varTab = plotTab))
}

#function to plot bootstrapped CI for feature coefficient
plotFeatureBoots <- function(glmResult, alphaCI = 0.95, freqCut = 0.1) {
    
    lapply(glmResult, function(eachRes) {
        coefMat <- eachRes$coefMat
        ciTab <- apply(coefMat, 1, calcCI, alphaCI) %>%
            t() %>% as_tibble(rownames = "feature") %>%
            mutate(meanCoef = rowMeans(coefMat, na.rm=TRUE),
                   freq = rowMeans(coefMat!=0, na.rm = TRUE)) %>%
            filter(freq >= freqCut) %>%
            arrange(desc(abs(meanCoef))) %>% mutate(feature = factor(feature, levels= feature))
        
        ggplot(ciTab, aes(x=feature, y=meanCoef, fill = freq)) +
            geom_bar(stat = "identity") +
            geom_errorbar(aes(ymax = upper, ymin = lower), width =0.2) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) 
    })
}


# For survival analysis

# Defien a color scheme, based on ggsci_NEJM panel, for the paper
colList <- c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF","#FFDC91FF","#EE4C97FF")

library(survival)
library(survminer)
library(maxstat)
#function for cox regression
com <- function(response, time, endpoint, scale =FALSE) {
    
    if (scale) {
        #calculate z-score
        response <- (response - mean(response, na.rm = TRUE))/sd(response, na.rm=TRUE)
    }
    surv <- coxph(Surv(time, endpoint) ~ response)
    
    
    tibble(p = summary(surv)[[7]][,5],
           HR = summary(surv)[[7]][,2],
           lower = summary(surv)[[8]][,3],
           higher = summary(surv)[[8]][,4])
}

# Function for Kaplan-Meier plot
km <- function(response, time, endpoint, titlePlot = "KM plot", pval = NULL,
               stat = "median", maxTime =NULL, showP = TRUE, showTable = FALSE,
               ylab = "Fraction", xlab = "Time (months)",
               table_ratio = c(0.7,0.3), yLabelAdjust = 0) {
    #function for km plot
    survS <- tibble(time = time,
                    endpoint = endpoint)
    
    if (!is.null(maxTime))
        survS <- mutate(survS, endpoint = ifelse(time > maxTime, FALSE, endpoint),
                        time = ifelse(time > maxTime, maxTime, time))
    
    if (stat == "maxstat") {
        ms <- maxstat.test(Surv(time, endpoint)  ~ response,
                           data = survS,
                           smethod = "LogRank",
                           minprop = 0.2,
                           maxprop = 0.8,
                           alpha = NULL)
        
        survS$group <- factor(ifelse(response >= ms$estimate, "less sensitive", "sensitive"))
        p <- com(survS$group, survS$time, survS$endpoint)$p
        
    } else if (stat == "median") {
        med <- median(response, na.rm = TRUE)
        survS$group <- factor(ifelse(response >= med, "less sensitive", "sensitive"))
        p <- com(survS$group, survS$time, survS$endpoint)$p
        
    } else if (stat == "binary") {
        survS$group <- factor(response)
        if (nlevels(survS$group) > 2) {
            sdf <- survdiff(Surv(survS$time,survS$endpoint) ~ survS$group)
            p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
        } else {
            p <- com(survS$group, survS$time, survS$endpoint)$p
        }
    }
    
    if (is.null(pval)) {
        if(p< 1e-16) {
            pAnno <- bquote(italic("P")~"< 1e-16")
        } else {
            pval <- formatNum(p, digits = 1)
            pAnno <- bquote(italic("P")~"="~.(pval))
        }
        
    } else {
        pval <- formatNum(pval, digits = 1)
        pAnno <- bquote(italic("P")~"="~.(pval))
    }
    
    if (!showP) pAnno <- ""
    
    colorPal <- rev(colList[1:length(unique(survS$group))])
    p <- ggsurvplot(survfit(Surv(time, endpoint) ~ group, data = survS),
                    data = survS, pval = FALSE,  conf.int = FALSE, palette = colorPal,
                    legend = ifelse(showTable, "none","bottom"),
                    legend.title=bquote(italic("ex vivo")~"sensitivity"),
                    surv.median.line = "hv",
                    ylab = "Event-free survival probability", xlab = "Time (in years)", title = titlePlot,
                    risk.table = showTable, legend.labs = sort(unique(survS$group)),
                    ggtheme = theme_half + theme(panel.border = element_blank(),
                                                 axis.title.y = element_text(vjust =yLabelAdjust)))
    if (!showTable) {
        p <- p$plot + annotate("text",label=pAnno, x = 3, y=0.8, hjust =0, size =5)
        return(p)
    } else {
        #construct a gtable
        pp <- p$plot + annotate("text",label=pAnno, x = 3, y=0.8, hjust =0, size=5)
        pt <- p$table + ylab("") + xlab("") 
        p <- plot_grid(pp,pt, rel_heights = table_ratio, nrow =2, align = "v")
        return(p)
    }
}

#function for plot hazard ratio
plotHazard <- function(survRes, title = "") {
    sumTab <- summary(survRes)$coefficients
    confTab <- summary(survRes)$conf.int
    #correct feature name
    nameOri <- rownames(sumTab)
    nameMod <- substr(nameOri, 1, nchar(nameOri) -1)
    plotTab <- tibble(feature = rownames(sumTab),
                      nameMod = substr(nameOri, 1, nchar(nameOri) -1),
                      HR = sumTab[,2],
                      p = sumTab[,5],
                      Upper = confTab[,4],
                      Lower = confTab[,3]) %>%
        mutate(feature = ifelse(nameMod %in% names(survRes$xlevels), nameMod, feature)) %>%
        mutate(feature = str_replace(feature, "[.]","/")) %>%
        mutate(feature = str_replace(feature, "[_]","-")) %>%
        arrange(desc(abs(p))) %>% mutate(feature = factor(feature, levels = feature)) %>%
        mutate(type = ifelse(HR >1 ,"up","down")) %>%
        mutate(Upper = ifelse(Upper > 10, 10, Upper))
    
    ggplot(plotTab, aes(x=feature, y = HR, color = type)) +
        geom_hline(yintercept = 1, linetype = "dotted", color = "grey50") +
        geom_point(position = position_dodge(width=0.8), size=3, color = "black") +
        geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.3, size=1,color = "grey20") +
        geom_text(position = position_nudge(x = 0.3),
                  aes(y = HR, label =  sprintf("italic(P)~'='~'%s'",
                                               formatNum(p, digits = 1))),
                  color = "black", size =3, parse = TRUE) +
        expand_limits(y=c(-0.5,0))+
        #scale_color_manual(values = c(up = colList[1], down = colList[2])) +
        ggtitle(title) + scale_y_log10() +
        ylab("Hazard ratio") +
        coord_flip() +
        theme_full +
        theme(legend.position = "none", axis.title.y = element_blank())
}

#Function to run multivariate Cox model on test table
runCox <- function(testTab) {
    surv1 <- coxph(
        Surv(time, endpoint) ~
            .,
        data = testTab)
    return(surv1)
}

#Function to calculate C-index
calcConcord <- function(testTab) {
    surv1 <- coxph(
        Surv(time, endpoint) ~
            .,
        data = testTab)
    cc <- summary(surv1)$concordance
    return(cc)
}

#compare concord
compareCindex <- function(eachTab, drugName) {
   ccFull <- calcConcord(eachTab)
   ccDrug <- calcConcord(select(eachTab,-ELN22))
   ccELN <- calcConcord(select(eachTab, -auc))
   plotTab <- tibble(model = c(paste0("ELN22 + ",drugName), drugName, "ELN22"),
                     cIndex = c(ccFull[1], ccDrug[1], ccELN[1]),
                     se = c(ccFull[2], ccDrug[2], ccELN[2])) %>%
       mutate(model = factor(model, levels = model))
   g <- ggplot(plotTab, aes(x=model, y=cIndex, fill = model)) +
       geom_bar(stat = "identity") +
       geom_errorbar(aes(max = cIndex + se, min=cIndex-se), width =0.5) +
       theme_bw() +
       theme(legend.position = "none") +
       coord_flip() +
       xlab("") + ylab("C-index")
   return(g)
}

#Function to format floats
formatNum <- function(i, limit = 0.01, digits =1, format="e") {
    r <- sapply(i, function(n) {
        if (n < limit) {
            formatC(n, digits = digits, format = format)
        } else {
            format(n, digits = digits)
        }
    })
    return(r)
}

#customized pieDonuat function
PieDonutCustom <- function (data, mapping, start = getOption("PieDonut.start", 
                                                             0), addPieLabel = TRUE, addDonutLabel = TRUE, showRatioDonut = TRUE, 
                            showRatioPie = TRUE, ratioByGroup = TRUE, showRatioThreshold = getOption("PieDonut.showRatioThreshold", 
                                                                                                     0.02), labelposition = getOption("PieDonut.labelposition", 
                                                                                                                                      2), labelpositionThreshold = 0.1, r0 = getOption("PieDonut.r0", 
                                                                                                                                                                                       0.3), r1 = getOption("PieDonut.r1", 1), r2 = getOption("PieDonut.r2", 
                                                                                                                                                                                                                                              1.2), explode = NULL, selected = NULL, explodePos = 0.1, 
                            color = "white", pieAlpha = 0.8, donutAlpha = 1, maxx = NULL, 
                            showPieName = TRUE, showDonutName = FALSE, title = NULL, 
                            pieLabelSize = 4, donutLabelSize = 3, titlesize = 5, explodePie = TRUE, 
                            explodeDonut = FALSE, use.label = TRUE, use.labels = TRUE, 
                            family = getOption("PieDonut.family", ""), piLabel = "", donutLabel = "")
{
    colList <- c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF","#FFDC91FF","#EE4C97FF")
    (cols = colnames(data))
    if (use.labels) 
        data = moonBook::addLabelDf(data, mapping)
    count <- NULL
    if ("count" %in% names(mapping)) 
        count <- moonBook::getMapping(mapping, "count")
    count
    pies <- donuts <- NULL
    (pies = moonBook::getMapping(mapping, "pies"))
    if (is.null(pies)) 
        (pies = moonBook::getMapping(mapping, "pie"))
    if (is.null(pies)) 
        (pies = moonBook::getMapping(mapping, "x"))
    (donuts = moonBook::getMapping(mapping, "donuts"))
    if (is.null(donuts)) 
        (donuts = moonBook::getMapping(mapping, "donut"))
    if (is.null(donuts)) 
        (donuts = moonBook::getMapping(mapping, "y"))
    if (!is.null(count)) {
        df <- data %>% group_by(.data[[pies]]) %>% dplyr::summarize(Freq = sum(.data[[count]]))
        df
    }
    else {
        df = data.frame(table(data[[pies]]))
    }
    colnames(df)[1] = pies
    df$end = cumsum(df$Freq)
    df$start = dplyr::lag(df$end)
    df$start[1] = 0
    total = sum(df$Freq)
    df$start1 = df$start * 2 * pi/total
    df$end1 = df$end * 2 * pi/total
    df$start1 = df$start1 + start
    df$end1 = df$end1 + start
    df$focus = 0
    if (explodePie) 
        df$focus[explode] = explodePos
    df$mid = (df$start1 + df$end1)/2
    df$x = ifelse(df$focus == 0, 0, df$focus * sin(df$mid))
    df$y = ifelse(df$focus == 0, 0, df$focus * cos(df$mid))
    df$label = df[[pies]]
    df$ratio = df$Freq/sum(df$Freq)
    if (showRatioPie) {
        df$label = ifelse(df$ratio >= showRatioThreshold, paste0(df$label, 
                                                                 "\n(", scales::percent(df$ratio), ")"), 
                          as.character(df$label))
    }
    df$labelx = (r0 + r1)/2 * sin(df$mid) + df$x
    df$labely = (r0 + r1)/2 * cos(df$mid) + df$y
    if (!is.factor(df[[pies]])) 
        df[[pies]] <- factor(df[[pies]])
    df
    mainCol = colList[1:nrow(df)]
    df$radius = r1
    df$radius[df$focus != 0] = df$radius[df$focus != 0] + df$focus[df$focus != 
                                                                       0]
    df$hjust = ifelse((df$mid%%(2 * pi)) > pi, 1, 0)
    df$vjust = ifelse(((df$mid%%(2 * pi)) < (pi/2)) | (df$mid%%(2 * 
                                                                    pi) > (pi * 3/2)), 0, 1)
    df$segx = df$radius * sin(df$mid)
    df$segy = df$radius * cos(df$mid)
    df$segxend = (df$radius + 0.05) * sin(df$mid)
    df$segyend = (df$radius + 0.05) * cos(df$mid)
    df
    if (!is.null(donuts)) {
        subColor = makeSubColor(mainCol, no = length(unique(data[[donuts]])))
        subColor
        data
        if (!is.null(count)) {
            df3 <- as.data.frame(data[c(donuts, pies, count)])
            colnames(df3) = c("donut", "pie", "Freq")
            df3
            df3 <- eval(parse(text = "complete(df3,donut,pie)"))
            df3$Freq[is.na(df3$Freq)] = 0
            if (!is.factor(df3[[1]])) 
                df3[[1]] = factor(df3[[1]])
            if (!is.factor(df3[[2]])) 
                df3[[2]] = factor(df3[[2]])
            df3 <- df3 %>% arrange(.data$pie, .data$donut)
            a <- df3 %>% spread(.data$pie, value = .data$Freq)
            a = as.data.frame(a)
            a
            rownames(a) = a[[1]]
            a = a[-1]
            a
            colnames(df3)[1:2] = c(donuts, pies)
        }
        else {
            df3 = data.frame(table(data[[donuts]], data[[pies]]), 
                             stringsAsFactors = FALSE)
            colnames(df3)[1:2] = c(donuts, pies)
            a = table(data[[donuts]], data[[pies]])
            a
        }
        a
        df3
        df3$group = rep(colSums(a), each = nrow(a))
        df3$pie = rep(1:ncol(a), each = nrow(a))
        total = sum(df3$Freq)
        total
        df3$ratio1 = df3$Freq/total
        df3
        if (ratioByGroup) {
            df3$ratio = scales::percent(df3$Freq/df3$group)
        }
        else {
            df3$ratio <- scales::percent(df3$ratio1)
        }
        df3$end = cumsum(df3$Freq)
        df3
        df3$start = dplyr::lag(df3$end)
        df3$start[1] = 0
        df3$start1 = df3$start * 2 * pi/total
        df3$end1 = df3$end * 2 * pi/total
        df3$start1 = df3$start1 + start
        df3$end1 = df3$end1 + start
        df3$mid = (df3$start1 + df3$end1)/2
        df3$focus = 0
        if (!is.null(selected)) {
            df3$focus[selected] = explodePos
        }
        else if (!is.null(explode)) {
            selected = c()
            for (i in 1:length(explode)) {
                start = 1 + nrow(a) * (explode[i] - 1)
                selected = c(selected, start:(start + nrow(a) - 
                                                  1))
            }
            selected
            df3$focus[selected] = explodePos
        }
        df3
        df3$x = 0
        df3$y = 0
        df
        if (!is.null(explode)) {
            explode
            for (i in 1:length(explode)) {
                xpos = df$focus[explode[i]] * sin(df$mid[explode[i]])
                ypos = df$focus[explode[i]] * cos(df$mid[explode[i]])
                df3$x[df3$pie == explode[i]] = xpos
                df3$y[df3$pie == explode[i]] = ypos
            }
        }
        df3$no = 1:nrow(df3)
        df3$label = df3[[donuts]]
        if (showRatioDonut) {
            if (max(nchar(levels(df3$label))) <= 2) 
                df3$label = paste0(df3$label, "(", df3$ratio, 
                                   ")")
            else df3$label = paste0(df3$label, "\n(", df3$ratio, 
                                    ")")
        }
        df3$label[df3$ratio1 == 0] = ""
        df3$label[df3$ratio1 < showRatioThreshold] = ""
        df3$hjust = ifelse((df3$mid%%(2 * pi)) > pi, 1, 0)
        df3$vjust = ifelse(((df3$mid%%(2 * pi)) < (pi/2)) | (df3$mid%%(2 * 
                                                                           pi) > (pi * 3/2)), 0, 1)
        df3$no = factor(df3$no)
        df3
        labelposition
        if (labelposition > 0) {
            df3$radius = r2
            if (explodeDonut) 
                df3$radius[df3$focus != 0] = df3$radius[df3$focus != 
                                                            0] + df3$focus[df3$focus != 0]
            df3$segx = df3$radius * sin(df3$mid) + df3$x
            df3$segy = df3$radius * cos(df3$mid) + df3$y
            df3$segxend = (df3$radius + 0.05) * sin(df3$mid) + 
                df3$x
            df3$segyend = (df3$radius + 0.05) * cos(df3$mid) + 
                df3$y
            if (labelposition == 2) 
                df3$radius = (r1 + r2)/2
            df3$labelx = (df3$radius) * sin(df3$mid) + df3$x
            df3$labely = (df3$radius) * cos(df3$mid) + df3$y
        }
        else {
            df3$radius = (r1 + r2)/2
            if (explodeDonut) 
                df3$radius[df3$focus != 0] = df3$radius[df3$focus != 
                                                            0] + df3$focus[df3$focus != 0]
            df3$labelx = df3$radius * sin(df3$mid) + df3$x
            df3$labely = df3$radius * cos(df3$mid) + df3$y
        }
        df3$segx[df3$ratio1 == 0] = 0
        df3$segxend[df3$ratio1 == 0] = 0
        df3$segy[df3$ratio1 == 0] = 0
        df3$segyend[df3$ratio1 == 0] = 0
        if (labelposition == 0) {
            df3$segx[df3$ratio1 < showRatioThreshold] = 0
            df3$segxend[df3$ratio1 < showRatioThreshold] = 0
            df3$segy[df3$ratio1 < showRatioThreshold] = 0
            df3$segyend[df3$ratio1 < showRatioThreshold] = 0
        }
        df3
        del = which(df3$Freq == 0)
        del
        if (length(del) > 0) 
            subColor <- subColor[-del]
        subColor
    }
    p <- ggplot() + ggforce::theme_no_axes() + coord_fixed()
    if (is.null(maxx)) {
        r3 = r2 + 0.3
    }
    else {
        r3 = maxx
    }
    p1 <- p + ggforce::geom_arc_bar(aes_string(x0 = "x", y0 = "y", 
                                               r0 = as.character(r0), r = as.character(r1), start = "start1", 
                                               end = "end1", fill = pies), alpha = pieAlpha, color = color, 
                                    data = df) + transparent() + scale_fill_manual(values = mainCol) + 
        xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
    if ((labelposition == 1) & (is.null(donuts))) {
        p1 <- p1 + geom_segment(aes_string(x = "segx", 
                                           y = "segy", xend = "segxend", yend = "segyend"), 
                                data = df) + geom_text(aes_string(x = "segxend", 
                                                                  y = "segyend", label = "label", hjust = "hjust", 
                                                                  vjust = "vjust"), size = pieLabelSize, data = df, 
                                                       family = family)
    }
    else if ((labelposition == 2) & (is.null(donuts))) {
        p1 <- p1 + geom_segment(aes_string(x = "segx", 
                                           y = "segy", xend = "segxend", yend = "segyend"), 
                                data = df[df$ratio < labelpositionThreshold, ]) + 
            geom_text(aes_string(x = "segxend", y = "segyend", 
                                 label = "label", hjust = "hjust", 
                                 vjust = "vjust"), size = pieLabelSize, 
                      data = df[df$ratio < labelpositionThreshold, 
                      ], family = family) + geom_text(aes_string(x = "labelx", 
                                                                 y = "labely", label = "label"), size = pieLabelSize, 
                                                      data = df[df$ratio >= labelpositionThreshold, ], 
                                                      family = family)
    }
    else {
        p1 <- p1 + geom_text(aes_string(x = "labelx", y = "labely", 
                                        label = "label"), size = pieLabelSize, data = df, 
                             family = family)
    }
    if (showPieName) 
        p1 <- p1 + annotate("text", x = 0, y = 0, label = pies, 
                            size = titlesize, family = family)
    p1 <- p1 + theme(text = element_text(family = family))
    if (!is.null(donuts)) {
        if (explodeDonut) {
            p3 <- p + ggforce::geom_arc_bar(aes_string(x0 = "x", 
                                                       y0 = "y", r0 = as.character(r1), r = as.character(r2), 
                                                       start = "start1", end = "end1", fill = "no", 
                                                       explode = "focus"), alpha = donutAlpha, 
                                            color = color, data = df3)
        }
        else {
            p3 <- p + ggforce::geom_arc_bar(aes_string(x0 = "x", 
                                                       y0 = "y", r0 = as.character(r1), r = as.character(r2), 
                                                       start = "start1", end = "end1", fill = "no"), 
                                            alpha = donutAlpha, color = color, data = df3)
        }
        p3 <- p3 + transparent() + scale_fill_manual(values = subColor) + 
            xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
        p3
        if (labelposition == 1) {
            p3 <- p3 + geom_segment(aes_string(x = "segx", 
                                               y = "segy", xend = "segxend", yend = "segyend"), 
                                    data = df3) + geom_text(aes_string(x = "segxend", 
                                                                       y = "segyend", label = "label", hjust = "hjust", 
                                                                       vjust = "vjust"), size = donutLabelSize, 
                                                            data = df3, family = family)
        }
        else if (labelposition == 0) {
            p3 <- p3 + geom_text(aes_string(x = "labelx", 
                                            y = "labely", label = "label"), size = donutLabelSize, 
                                 data = df3, family = family)
        }
        else {
            p3 <- p3 + geom_segment(aes_string(x = "segx", 
                                               y = "segy", xend = "segxend", yend = "segyend"), 
                                    data = df3[df3$ratio1 < labelpositionThreshold, 
                                    ]) + geom_text(aes_string(x = "segxend", 
                                                              y = "segyend", label = "label", hjust = "hjust", 
                                                              vjust = "vjust"), size = donutLabelSize, 
                                                   data = df3[df3$ratio1 < labelpositionThreshold, 
                                                   ], family = family) + geom_text(aes_string(x = "labelx", 
                                                                                              y = "labely", label = "label"), size = donutLabelSize, 
                                                                                   data = df3[df3$ratio1 >= labelpositionThreshold, 
                                                                                   ], family = family)
        }
        if (!is.null(title)) 
            p3 <- p3 + annotate("text", x = 0, y = r3, 
                                label = title, size = titlesize, family = family)
        else if (showDonutName) 
            p3 <- p3 + annotate("text", x = (-1) * r3, 
                                y = r3, label = donuts, hjust = 0, size = titlesize, 
                                family = family)
        p3 <- p3 + theme(text = element_text(family = family))
        grid::grid.newpage()
        print(p1 + annotate(geom = "text", x=0, y=0, label = piLabel, fontface= "bold" ), vp = grid::viewport(height = 1, width = 1))
        print(p3 + annotate(geom = "text", x=-0.8, y=0.8, fontface = "bold", label = donutLabel), vp = grid::viewport(height = 1, width = 1))
    }
    else {
        p1
    }
}


mscale <- function(x, center = TRUE, scale = TRUE, censor = NULL, useMad = FALSE){
    if (scale & center) {
        if (useMad) {
            x.scaled <- apply(x, 1, function(y) (y-median(y,na.rm = T))/meanAD(y))
        } else {
            x.scaled <- apply(x, 1, function(y) (y-mean(y,na.rm=T))/sd(y,na.rm = T))
        }
    } else if (center & !scale) {
        if (useMad) {
            x.scaled <- apply(x, 1, function(y) (y-median(y,na.rm=T)))
        } else {
            x.scaled <- apply(x, 1, function(y) (y-mean(y,na.rm=T)))
        }
    } else if (!center & scale) {
        if (useMad) {
            x.scaled <- apply(x, 1, function(y) y/meanAD(y))
        } else {
            x.scaled <- apply(x, 1, function(y) y/sd(y,na.rm = T))
        }
    } else {
        x.scaled <- t(x)
    }
    
    if (!is.null(censor)) {
        if (length(censor) == 1) {
            x.scaled[x.scaled > censor] <- censor
            x.scaled[x.scaled < -censor] <- -censor
        } else {
            x.scaled[x.scaled > censor[2]] <- censor[2] #higher limit
            x.scaled[x.scaled < censor[1]] <- censor[1] #lower limit
        }
    }
    return(t(as.matrix(x.scaled)))
}

# Function for IC50 curves

fitIC50 <- function(formula, data = NULL, weights = NULL, logDose = NULL, ...) {
    if (!is.null(data)) {
        modelFrame <- model.frame(formula, data)
    } else {
        modelFrame <- model.frame(formula)
    }
    
    if (!is.null(logDose)) {
        modelFrame[, 2] <- logDose^modelFrame[, 2]
    }
    
    parm_fit <- dr4pl::dr4pl(modelFrame[, 2], modelFrame[, 1], ...)
    newModel <- list(model = modelFrame, formula = formula, parm_fit = parm_fit, logDose = logDose)
    
    class(newModel) <- "fitIC50"
    return(newModel)
}


predict.fitIC50 <- function(object, newdata = NULL, se.fit = FALSE, level = 0.95,
                            interval = c("none", "confidence", "prediction"), ...) {
    
    if (is.null(newdata))
        newdata <- object$model else newdata <- newdata
        
        params <- object$parm_fit$parameters
        logDose <- object$logDose
        
        a <- min(params[1], params[4])
        d <- max(params[4], params[1])
        c <- params[2]
        b <- params[3]
        
        predY <- function(x) {
            y = d + (a - d)/(1 + (x/c)^b)
            names(y) <- NULL
            return(y)
        }
        
        if (is.vector(newdata)) {
            conc <- newdata
        } else if (is.data.frame(newdata)) {
            if (ncol(newdata) == 1) {
                conc <- newdata[, 1]
            } else {
                concName <- as.character(object$formula)[3]
                conc <- newdata[, concName]
            }
        }
        
        if (!is.null(logDose))
            conc <- logDose^conc
        res <- vapply(conc, predY, numeric(1))
        
        return(res)
}


#set the global ggplot theme
theme_full <- theme_bw() + theme(axis.line = element_blank(),
                                 panel.border = element_rect(size = 1.5),
                                 axis.ticks = element_line(),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 plot.title = element_text(size=15, face = "bold"),
                                 axis.title = element_text(size=15),
                                 axis.text = element_text(size=15),
                                 legend.text = element_text(size=14),
                                 legend.title = element_text(size=14))

theme_half <- theme_bw() + theme(axis.line =element_line(),
                                 panel.border = element_rect(size=1.5),
                                 axis.ticks = element_line(),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 plot.title = element_text(size=15, face = "bold"),
                                 axis.title = element_text(size=15),
                                 axis.text = element_text(size=15),
                                 legend.text = element_text(size=14),
                                 legend.title = element_text(size=14))
