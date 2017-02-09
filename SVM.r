##main function
#default parameter values: fold is for cross-val, para=FALSE, it will not output SVM parameter values, n is the times cross-val will be done
model.test=function(training.data,pred.data,fold=5,model="svm",method="radial", gamma, cost, output.file,para="FALSE",n=100, m=20){
    library(e1071)
    library(ROCR)
    library(verification)
    # progress bar, indicates the progress of the program
    #library(tcltk)
    par(mfrow = c(3,3),mar=c(5, 6, 4, 2))
    ##check arguments##
    if (is.numeric(model)||is.numeric(training.data)||is.numeric(pred.data))
    stop("invalid argument!")
    
    input=training.data
    ##save input data## change 1 to 2, 0 to 1
    input.num=input
    input[,1][input[,1]==1]="2"
    input[,1][input[,1]==0]="1"
    # change numeric to factor
    input[,1]=factor(input[,1])
    
    
    ##############################################
    #             parameter testing              #
    ##############################################
    
    ###########search for best parameter combination###########
    best.gamma=0
    best.cost=0
    if(model=="svm"){
        cat("\n####################################\n")
        cat("#Detecting parameters for SVM model#\n")
        cat("####################################\n")
        cat("Searching, please wait...\n")
        # txtProgressBar is used for counting the number of parameter combinations
        tot=length(gamma)*length(cost)
        #pb <- tkProgressBar(title = "Searching for parameters, please wait..", min = 0, max =tot , width = 500)
        k=0
        best.perf=999
        sum.para=c()
        for(i in gamma){
            for(j in cost){
                k=k+1
                #setTkProgressBar(pb, k, label=paste( round(k/tot*100, 0),"% done"))
                obj <- tune.svm(CellType ~ ., data = input, gamma = i, cost = j,kernel =method,cross=fold)
                # the smaller best.performance, the smaller error
                obj.out=c(i,j,obj$best.performance)
                sum.para=rbind(sum.para,obj.out)
                if(best.perf>obj$best.performance){
                    best.gamma=obj$best.parameters["gamma"]
                    best.cost=obj$best.parameters["cost"]
                    best.perf=obj$best.performance
                }
            }
        }
        #close(pb)
        cat("\nBest parameter values:\n")
        cat("Gamma: ",best.gamma[1,1],"\n")
        cat("Cost: ",best.cost[1,1],"\n\n")
        cat("best.performance=error: ",best.perf,"\n\n")
    }
    
    
    ##############################################
    #          prediction and testing            #
    ##############################################
    cat("\n####################################\n")
    cat(  "#     prediction and testing       #\n")
    cat(  "####################################\n")
    cat("randomizing , please wait...\n")
    ###########prediction###########
    if(model=="svm"){
        svm.model <- svm(CellType ~ .,kernel = method, data = input, gamma = best.gamma, cost = best.cost, probability = TRUE)
        ##test using training data itself##
        pred.v.test <- predict(svm.model, input[, -1], decision.values = TRUE, probability = TRUE)
        pred.prob.real.test=round(attributes(pred.v.test)$probabilities[,2],3)
        ##pred using genome-wide data##
        pred.v <- predict(svm.model, pred.data, decision.values = TRUE, probability = TRUE)
        pred.prob.real=round(attributes(pred.v)$probabilities[,2],3)
    }
    if(model == "logistic"){
        glm.model <- glm(CellType ~ ., family='binomial', input)
        ##test##
        glm.pred.test<-predict.glm(glm.model, input[, -1])
        pred.prob.real.test=round(exp(glm.pred.test)/(1+exp(glm.pred.test)),3)
        ##pred##
        glm.pred<-predict.glm(glm.model, pred.data)
        pred.prob.real=round(exp(glm.pred)/(1+exp(glm.pred)),3) #same as: glm.pred.prob <- predict(glm.model, newdata=input[, -1], type="response");
    }
    ##hist1
    hist(pred.prob.real, col="lightblue",xlim=c(0,1),main="Genome Wide Prediction, real", xlab="Probality", ylab="#Number Of Gene")
    ##hist2
    hist(pred.prob.real.test, col="lightblue",xlim=c(0,1),main="training data Prediction, real", xlab="Probality", ylab="#Number Of Gene")
    ##hist3
    tr.data=cbind(pred.prob.real.test,input[,1])
    tr.pos=tr.data[tr.data[,2]==2,]
    tr.neg=tr.data[tr.data[,2]==1,]
    hist(tr.pos[,1], col="red", xlim=c(0,1), ylim=c(0,120),main="Training Data, Prediction vs real ", xlab="Probality", ylab="#Number Of Gene")
    hist(tr.neg[,1], col=rgb(0, 1, 0, 0.5), xlim=c(0,1), add=T)
    ##randomization##
    sum.pred.prob.rand=c()
    #pb <- tkProgressBar(title = "randomization, please wait..", min = 0, max =m , width = 500)
    for (i in 1:m){
        #setTkProgressBar(pb, i, label=paste( round(i/m*100, 0),"% done"))
        input3=training.data
        CellType=sample(input3[,1],length(input3[,1]))
        input.rand=cbind(CellType,input3[,-1])
        input.rand[,1][input.rand[,1]==1]="1"
        input.rand[,1][input.rand[,1]==0]="0"
        # change numeric to factor
        input.rand[,1]=factor(input.rand[,1])
        
        if(model=="svm"){
            svm.model <- svm(CellType ~ .,kernel = method, data = input.rand, gamma = best.gamma, cost = best.cost, probability = TRUE)
            ##test##
            pred.v.svm.test <- predict(svm.model, input[, -1], decision.values = TRUE, probability = TRUE)
            pred.prob.rand.test=round(attributes(pred.v.svm.test)$probabilities[,2],3)
            # need to define sum.pred.prob.rand.test as a new variable
            if(i==1){
                sum.pred.prob.rand.test=pred.prob.rand.test
            }
            if(i>1){
                sum.pred.prob.rand.test=sum.pred.prob.rand.test+pred.prob.rand.test
            }
            ##pred##
            pred.v.svm <- predict(svm.model, pred.data, decision.values = TRUE, probability = TRUE)
            pred.prob.rand=round(attributes(pred.v.svm)$probabilities[,2],3)
            if(i==1){
                sum.pred.prob.rand=pred.prob.rand
            }
            if(i>1){
                sum.pred.prob.rand=sum.pred.prob.rand+pred.prob.rand
            }
        }
        if(model == "logistic"){
            glm.model <- glm(CellType ~ ., family='binomial', input.rand)
            ##test##
            glm.pred <- predict.glm(glm.model, input[, -1])
            pred.prob.rand.test = round(exp(glm.pred)/(1+exp(glm.pred)),3) #same as: glm.pred.prob <- predict(glm.model, newdata=input[, -1], type="response");
            if(i==1){
                sum.pred.prob.rand.test=pred.prob.rand.test
            }
            if(i>1){
                sum.pred.prob.rand.test=sum.pred.prob.rand.test+pred.prob.rand.test
            }
            ##pred##
            glm.pred <- predict.glm(glm.model, pred.data)
            pred.prob.rand = round(exp(glm.pred)/(1+exp(glm.pred)),3)
            if(i==1){
                sum.pred.prob.rand=pred.prob.rand
            }
            if(i>1){
                sum.pred.prob.rand=sum.pred.prob.rand+pred.prob.rand
            }
        }
    }
    # close(pb)
    ##test##
    avg.pred.prob.rand.test=sum.pred.prob.rand.test/m
    ##pred##
    avg.pred.prob.rand=sum.pred.prob.rand/m
    
    ##classification cutoff##
    #if the predicted prob > 0.5 AND the predicted prob > averaged randomized predicted prob, we consider it a predicted 1
    ##test using training data##
    pred.v.rand.final.test <- c(0,0,1)[((pred.prob.real.test>0.5) + (pred.prob.real.test > avg.pred.prob.rand.test)) + 1]
    ##pred using entire data##
    pred.v.rand.final <- c(0,0,1)[((pred.prob.real>0.5) + (pred.prob.real > avg.pred.prob.rand)) + 1]
    
    ##hist1
    hist(avg.pred.prob.rand, col="lightblue",main="Genome Wide Prediction, random", xlab="Probality", ylab="#Number Of Gene")
    ##hist2
    hist(avg.pred.prob.rand.test, col="lightblue",main="training data Prediction, random", xlab="Probality", ylab="#Number Of Gene")
    ##hist3
    tr.data=cbind(avg.pred.prob.rand.test,input[,1])
    tr.pos=tr.data[tr.data[,2]==2,]
    tr.neg=tr.data[tr.data[,2]==1,]
    hist(tr.pos[,1], col="red", ylim=c(0,200),main="training data Prediction, random", xlab="Probality", ylab="#Number Of Gene")
    hist(tr.neg[,1], col=rgb(0, 1, 0, 0.5), add=T)
    
    ##################### predicted ####################
    #                     NOT USE!
    #
    #        tr.data=cbind(pred.prob.real.test,pred.v.rand.final.test)
    #        tr.pos=tr.data[tr.data[,2]==1,]
    #        tr.neg=tr.data[tr.data[,2]==0,]
    #        hist(tr.pos[,1], col="red", xlim=c(0,1), ylim=c(0,150),main="Training Data, Prediction vs random cutoff", xlab="Probality", ylab="#Number Of Gene")
    #        hist(tr.neg[,1], col=rgb(0, 1, 0, 0.5), xlim=c(0,1), add=T)
    #        hist(avg.pred.prob.rand.test, col=rgb(0,0,1,1/4),xlim=c(0,1), add=T)        #Genome Wide Prediction(Rand)
    
    ##output##
    output=cbind(pred.prob.real, pred.v.rand.final, avg.pred.prob.rand, pred.data);
    colnames(output[,1:4])=c("logist.prob","logist.pred","svm.prob","svm.pred");
    write.table(output, file = output.file,  sep = "\t", row.names = TRUE,quote = FALSE)
    ####################################################
    
    ###########model testing with entire data###########
    
    ###### prob>0.5 AND obs>randomization as cut off  #######
    cat("\n#randomization as cut off:\n")
    table.s=table(pred = pred.v.rand.final.test, true = input[, 1])
    FP=table.s[2,1]/sum(table.s[2,])
    FN=table.s[1,2]/sum(table.s[1,])
    cat("Prediction (Based on whole data set):\n")
    cat("False Positive Rate:\t",FP,"\n")
    cat("False Negative Rate:\t",FN,"\n\n")
    #	cat(table.s[1,1],"\n");
    #	cat(table.s[1,2],"\n");
    #	cat(table.s[2,1],"\n");
    #	cat(table.s[2,2],"\n\n");
    
    ###### prob=0.5 as cut off #######
    cat("#prob=0.5 as cut off:\n")
    ##test using training data##
    pred.v.rand.final.test.5 <- c(0,1)[(pred.prob.real.test>0.5) + 1]
    
    table.s=table(pred = pred.v.rand.final.test.5, true = input[, 1])
    FP=table.s[2,1]/sum(table.s[2,])
    FN=table.s[1,2]/sum(table.s[1,])
    cat("Prediction (Based on whole data set):\n")
    cat("False Positive Rate:\t",FP,"\n")
    cat("False Negative Rate:\t",FN,"\n\n")
    #	cat(table.s[1,1],"\n");
    #	cat(table.s[1,2],"\n");
    #	cat(table.s[2,1],"\n");
    #	cat(table.s[2,2],"\n\n");
    
    
    
    
    
    ##############################################
    #        cross-classification testing        #
    ##############################################
    cat("\n####################################\n")
    cat("#   cross-classification testing   #\n")
    cat("####################################\n")
    cat(n," round, cross-validation, please wait..\n")
    #pb <- tkProgressBar(title = "cross-validation, please wait..", min = 0, max = n, width = 500)
    # get index for input data
    index <- 1:nrow(input)
    num_fold=trunc(length(index)/fold)
    sum.data=c()
    sum.table=matrix(0,2,2)
    
    for(ii in 1:n){
        #setTkProgressBar(pb, ii, label=paste( round(ii/n*100, 0),"% done"))
        data.sum=c()
        table.sum=matrix(0,2,2)
        # randomize the order of training examples
        index.rand = sample(index, length(index))
        i=0
        while (i < fold){
            i=i+1
            ##select training and testing data
            if(i < fold){
                testindex=index.rand[(num_fold*(i-1)+1):(num_fold*i)]
            }else if(i == fold) {
                testindex=index.rand[(num_fold*(i-1)+1):length(index)]
            }
            testset <- input[testindex, ]
            trainset <- input[-testindex, ]
            
            ##predict
            if(model == "svm"){
                svm.model <- svm(CellType ~ .,kernel = method, data = trainset, gamma = best.gamma, cost = best.cost, probability = TRUE)
                pred.v <- predict(svm.model, testset[, -1], decision.values = TRUE, probability = TRUE)
                # if the 2X2 prediction matrix become 2X1 matrix (for example when using linear kernel), stop calculation, redo the randomization of the entire training data
                if((sum(pred.v=="2")==length(pred.v)) || (sum(pred.v=="1")==length(pred.v))){
                    i=0
                    data.sum=c()
                    table.sum=matrix(0,2,2)
                    index.rand = sample(index, length(index))
                    cat("restart!\n");
                    # go back to while loop, starting i=0
                    next;
                }
                # decision.values: continuous values, need to choose a cutoff value to determine two classes.
                # probabilities: continuous between 0 and 1, also need to choose a cutoff value to determine two classes.
                # ROC results are same for decision.values and probilities
                predict.data=cbind(input.num[testindex, 1],testset[, 1],pred.v,attributes(pred.v)$decision.values,attributes(pred.v)$probabilities[,2])
                # real.num=numeric response value of input, real=factor response value of input, predict=binary value
                colnames(predict.data)=c("real.num","real","predict","decisition.values","prob")
            }
            if(model == "logistic"){
                glm.model <- glm(CellType ~ ., family='binomial', data = trainset)
                glm.pred<-predict.glm(glm.model,testset[, -1])
                glm.pred.prob=exp(glm.pred)/(1+exp(glm.pred))
                pred.v <- c(0,1)[(glm.pred.prob>0.5)+1]
                if((sum(pred.v=="2")==length(pred.v)) || (sum(pred.v=="1")==length(pred.v))){
                    i=0
                    data.sum=c()
                    table.sum=matrix(0,2,2)
                    index.rand = sample(index, length(index))
                    cat("restart!\n");
                    next;
                }
                predict.data=cbind(input.num[testindex, 1],testset[, 1], glm.pred.prob)
                colnames(predict.data)=c("real.num","real","prob")
            }
            
            data.sum=rbind(data.sum,predict.data)
            table.s=table(pred = pred.v, true = testset[, 1])
            table.sum=table.sum+table.s
            #    cat(table.s,"\n")
            #    cat(table.sum,"\n")
        }
        sum.data=rbind(sum.data,data.sum)
        sum.table=sum.table+table.sum
    }
    # close(pb)
    
    ##output
    FP=sum.table[2,1]/sum(sum.table[2,])
    FN=sum.table[1,2]/sum(sum.table[1,])
    cat("Prediction (Based on ",fold," fold cross-validation):\n")
    cat("False Positive Rate:\t",FP,"\n")
    cat("False Negative Rate:\t",FN,"\n\n")
    #	cat(sum.table[1,1],"\n");
    #	cat(sum.table[1,2],"\n");
    #	cat(sum.table[2,1],"\n");
    #	cat(sum.table[2,2],"\n");
    
    ###########ROC plot###########
    cat("\n####################################\n")
    cat("#    plotting ROC and PC Curve     #\n")
    cat("####################################\n")
    cat("drawing graph, please wait...\n")
    # use ROCR package
    # real is 0 or 0, prob is continuous
    svm.roc <- prediction(sum.data[,"prob"], sum.data[,"real"])
    svm.pr <- performance(svm.roc, "prec", "rec")
    plot(svm.pr,colorize = TRUE,xaxis.cex.axis=1.6,yaxis.cex.axis=1.6,ylab="",xlab="")
    mtext("Precision-Recall", 3, line = 1, font = 2,cex = 1.7)
    mtext("Recall", 1, line = 2.8, font = 2,cex = 1.3)
    mtext("Precision", 2, line = 2.8, font = 2,cex = 1.3)
    
    svm.auc <- performance(svm.roc, 'tpr', 'fpr')
    plot(svm.auc,colorize = TRUE,xaxis.cex.axis=1.6,yaxis.cex.axis=1.6,ylab="",xlab="")
    mtext("ROC curve", 3, line = 1, font = 2,cex = 1.7)
    mtext("False Positive Rate", 1, line = 2.8, font = 2,cex = 1.3)
    mtext("True Positive Rate", 2, line = 2.8, font = 2,cex = 1.3)
    
    # verification package
    area=roc.plot(sum.data[,"real.num"], sum.data[,"prob"],threshold = seq(0.1,0.9,0.1),xlab = "False positive rate",ylab = "True positive rate", main = NULL, CI = T, n.boot = 100, plot = "both", binormal = TRUE)
    mtext("ROC curve", 3, line = 1, font = 2,cex = 1.7)
    cat("Plotting ROC...\n")
    cat("AUC: ",area$roc.vol[,2],"\n")
    
    cat("\n--END--\n")
}






