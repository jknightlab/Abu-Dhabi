# Number of GEx PCs to include in eQTL model

n.eqtl <- data.frame("PCs"=c(0, 5, 10, 15, 20, 25, 30, 35, 40), 
                     "eGenes"=c(4111, 5998, 6746, 7149, 7380, 7492, 7605, 7657, 7706))
n.eqtl$Diff <- NA
for(i in 2:nrow(n.eqtl)){
  n.eqtl$Diff[i] <- n.eqtl$eGenes[i] - n.eqtl$eGenes[i-1]
}

p1 <- ggplot(n.eqtl, aes(PCs, eGenes)) + geom_point() + geom_line() + theme_bw()
p2 <- ggplot(n.eqtl, aes(PCs, Diff)) + geom_point() + geom_line() + theme_bw()

pdf("U:/Abu-Dhabi/RNASeq/eQTL/Increasing_GEx_PCs_results.pdf")
p1
p2
dev.off()


# Regress out PCs and save adjusted GEx data
load("eQTL.RData")

PCs <- cbind(gex.pcs[, 1:25], geno.pcs[, 1:4])
expression.set <- exprs
expression.set <- t(expression.set)
expression.pc <- PCs
num.PC <- ncol(expression.pc)

# regress out expression
regressed.data <- matrix(nrow=nrow(expression.set), ncol=ncol(expression.set))
for(i in 1:(dim(expression.set)[2])){
  model <- lm(as.matrix(expression.set[, i]) ~ as.matrix(expression.pc[, 1:(num.PC)]))
  regressed.data[, i] <- expression.set[, i] - 
    rowSums(sapply(1:(num.PC), function(i)model$coefficients[i+1]*expression.pc[, 1:(num.PC)][,i]))
}
rownames(regressed.data) <- rownames(expression.set)
colnames(regressed.data) <- colnames(expression.set)
regressed.data <- t(regressed.data)

save(list=c("regressed.data", "geno", "info"), file="../../eQTL.25PCs.RData")