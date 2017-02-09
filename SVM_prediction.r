################################################################################################
#Goal: prediction by using SVM and Logistic model on cell type classifiation
#Author: Debjit Ray
#Date: 04/26/2011
################################################################################################
##### Female Human #####
##Input Data
trainning.data.hf=read.table("trainingData_Hf_final.txt",row.names=1,header = TRUE, sep="\t")
pred.data.hf=read.table("embryonic_human-h133_plus_2_Hf_only_GeneName.na30.txt",row.names=1, sep="\t")
colnames(pred.data.hf)=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13");
pred.data.hf=round(pred.data.hf,3)

##Example Command
model.test(training.data=trainning.data.hf, pred.data=pred.data.hf,fold=5, gamma = 2^(-8:2), cost = 2^(-3:15),output.file="Hf.griswold.txt")



##### Male Mouse #####
##Input Data
trainning.data.hf=read.table("trainingData_Mm_final.txt",row.names=1,header = TRUE, sep="\t")
pred.data.hf=read.table("embryonic_mouse-m430_2_new_mm_only_GeneName.na30.txt",row.names=1, sep="\t")
colnames(pred.data.hf)=c("TP1","TP2","TP3","TP4");
pred.data.hf=round(pred.data.hf,3)

##Example Command
model.test(training.data=trainning.data.hf, pred.data=pred.data.hf,fold=5, gamma = 2^(-7:0), cost = 2^(10:20),output.file="Mm.griswold.txt")



##### Female Mouse #####
##Input Data
trainning.data.hf=read.table("trainingData_Mf_griswold(MAS5)_final.txt",row.names=1,header = TRUE, sep="\t")
pred.data.hf=read.table("embryonic_mouse-m430_2_new_fm_only_GeneName.na30.txt",row.names=1, sep="\t")
colnames(pred.data.hf)=c("TP1","TP2","TP3","TP4");
pred.data.hf=round(pred.data.hf,3)

##Example Command
model.test(training.data=trainning.data.hf, pred.data=pred.data.hf,fold=5, gamma = 2^(-7:0), cost = 2^(10:20),,output.file="Mf.griswold.txt")
################################################################################################
##main function

#default parameter values: fold is for cross-val, para=FALSE, it will not output SVM parameter values, n is the times cross-val will be done
model.test=function(training.data,pred.data,fold=5, gamma = 2^(-3:0), cost = 2^(0:4),para="FALSE",output.file){
	library(e1071)
	library(ROCR)
	library(verification)
# progress bar, indicates the progress of the program
	library(tcltk)
	input=training.data
#        file.name<-paste(input,".pdf", sep = "")
#        pdf(file=file.name,width=10,height=5,pointsize=6);
##save input data## change 1 to 2, 0 to 1
	input.num=input
	input[,1][input[,1]==1]="1"
	input[,1][input[,1]==0]="0"
# change numeric to factor
	input[,1]=factor(input[,1])
#	str(input)
#	str(input.num)
##check arguments##
	if (is.numeric(model)||is.numeric(input)) 
          stop("invalid argument!")
	
##search for best parameter combination##
	best.gamma=0
	best.cost=0
		  cat("\nSearching for parameters, please wait..\n")	
                                        # txtProgressBar is used for counting the number of parameter combinations
          tot=length(gamma)*length(cost)
          pb <- tkProgressBar(title = "Searching for parameters, please wait..", min = 0, max =tot , width = 500)
          k=0
          best.perf=999
          sum.para=c()
          for(i in gamma){
            for(j in cost){
              k=k+1
              setTkProgressBar(pb, k, label=paste( round(k/tot*100, 0),"% done"))
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
          close(pb)
          cat("\nBest parameter values:\n")
          cat("Gamma: ",best.gamma[1,1],"\n")
          cat("Cost: ",best.cost[1,1],"\n\n")
		  cat("best.performance=error: ",best.perf,"\n\n")
	
	##whole data model##
		svm.model <- svm(CellType ~ .,kernel = method, data = input, gamma = best.gamma, cost = best.cost, probability = TRUE)
		pred.v.svm <- predict(svm.model, pred.data, decision.values = TRUE, probability = TRUE)
		pred.prob.svm=round(attributes(pred.v.svm)$probabilities[,2],3)

		glm.model <- glm(CellType ~ ., family='binomial', input)
		glm.pred<-predict.glm(glm.model, pred.data)
		pred.prob.glm=round(exp(glm.pred)/(1+exp(glm.pred)),3) #same as: glm.pred.prob <- predict(glm.model, newdata=input[, -1], type="response");
		pred.v.glm <- c(0,1)[(pred.prob.glm>0.5)+1]			
	
	output=cbind(pred.prob.glm,pred.v.glm,pred.prob.svm,pred.v.svm,pred.data);
	colnames(output[,1:4])=c("logist.prob","logist.pred","svm.prob","svm.pred");	
#dev.off();
	write.table(output, file = output.file,  sep = "\t", row.names = TRUE,quote = FALSE)
# output SVM parameter values
	if(para=="TRUE"){
		colnames(sum.para)=c("gamma","cost","model.error");
		cat("parameter test:\n")
		return(sum.para);
	}
}






