### R code from vignette source 'diffHicUG.Rnw'

###################################################
### code chunk number 1: diffHicUG.Rnw:77-80
###################################################
picdir <- "plots-ug"
if (file.exists(picdir)) { unlink(picdir, recursive=TRUE) }
dir.create(picdir)


###################################################
### code chunk number 2: diffHicUG.Rnw:123-126
###################################################
input <- c("merged_flox_1.h5", "merged_flox_2.h5", "merged_ko_1.h5", "merged_ko_2.h5")
fragments <- readRDS("mm10-hindIII.rds")
design <- model.matrix(~factor(c("flox", "flox", "ko", "ko")))


###################################################
### code chunk number 3: diffHicUG.Rnw:134-137
###################################################
require(diffHic)
param <- pairParam(fragments=fragments)
data <- squareCounts(input, width=1e6, param=param)


###################################################
### code chunk number 4: diffHicUG.Rnw:140-143
###################################################
require(edgeR)
keep <- aveLogCPM(asDGEList(data)) > 0
data <- data[keep,]


###################################################
### code chunk number 5: diffHicUG.Rnw:146-148
###################################################
y <- asDGEList(data)
y$offset <- normalize(data, type="loess")


###################################################
### code chunk number 6: diffHicUG.Rnw:151-153
###################################################
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)


###################################################
### code chunk number 7: diffHicUG.Rnw:156-157
###################################################
result <- glmQLFTest(fit)


###################################################
### code chunk number 8: diffHicUG.Rnw:190-191
###################################################
system.file("python", "presplit_map.py", package="diffHic", mustWork=TRUE)


###################################################
### code chunk number 9: diffHicUG.Rnw:204-207
###################################################
require(BSgenome.Hsapiens.UCSC.hg19)
hs.frag <- cutGenome(BSgenome.Hsapiens.UCSC.hg19, "AAGCTT", 4)
hs.frag


###################################################
### code chunk number 10: diffHicUG.Rnw:215-217
###################################################
hs.param <- pairParam(hs.frag)
hs.param


###################################################
### code chunk number 11: diffHicUG.Rnw:227-230
###################################################
diagnostics <- preparePairs("SRR027957.bam", hs.param, 
    file="SRR027957.h5", dedup=TRUE, minq=10)
diagnostics


###################################################
### code chunk number 12: diffHicUG.Rnw:257-258
###################################################
diagnostics$chimeras[["invalid"]]/diagnostics$chimeras[["multi"]]


###################################################
### code chunk number 13: diffHicUG.Rnw:273-277
###################################################
min.inward <- 1000
min.outward <- 25000
prunePairs("SRR027957.h5", hs.param, file.out="SRR027957_trimmed.h5", 
    max.frag=600, min.inward=min.inward, min.outward=min.outward)


###################################################
### code chunk number 14: fraglenhist (eval = FALSE)
###################################################
## diags <- getPairData("SRR027957.h5", hs.param)
## hist(diags$length[diags$length < 1000], ylab="Frequency", xlab="Spacing (bp)")


###################################################
### code chunk number 15: diffHicUG.Rnw:293-294
###################################################
diags <- getPairData("SRR027957.h5", hs.param)
hist(diags$length[diags$length < 1000], ylab="Frequency", xlab="Spacing (bp)")


###################################################
### code chunk number 16: diffHicUG.Rnw:313-315
###################################################
intra <- !is.na(diags$insert)
table(diags$orientation[!intra])


###################################################
### code chunk number 17: strorient (eval = FALSE)
###################################################
## llinsert <- log2(diags$insert + 1L)
## intra <- !is.na(llinsert)
## breaks <- seq(min(llinsert[intra]), max(llinsert[intra]), length.out=30)
## inward <- hist(llinsert[diags$orientation==1L], plot=FALSE, breaks=breaks)
## outward <- hist(llinsert[diags$orientation==2L] ,plot=FALSE, breaks=breaks)
## samestr <- hist(llinsert[diags$orientation==0L | diags$orientation==3L], 
##    plot=FALSE, breaks=breaks)
## samestr$counts <- samestr$counts/2
## #
## ymax <- max(inward$counts, outward$counts, samestr$counts)/1e6
## xmax <- max(inward$mids, outward$mids, samestr$mids)
## xmin <- min(inward$mids, outward$mids, samestr$mids)
## #
## plot(0,0,type="n", xlim=c(xmin, xmax), ylim=c(0, ymax),
##     xlab=expression(log[2]~"[insert size (bp)]"), ylab="Frequency (millions)")
## lines(inward$mids, inward$counts/1e6, col="darkgreen", lwd=2)
## abline(v=log2(min.inward), col="darkgrey")
## lines(outward$mids, outward$counts/1e6, col="red", lwd=2)
## abline(v=log2(min.outward), col="darkgrey", lty=2)
## lines(samestr$mids, samestr$counts/1e6, col="blue", lwd=2)
## legend("topright", c("inward", "outward", "same"), 
##     col=c("darkgreen", "red", "blue"), lwd=2)


###################################################
### code chunk number 18: diffHicUG.Rnw:352-353
###################################################
llinsert <- log2(diags$insert + 1L)
intra <- !is.na(llinsert)
breaks <- seq(min(llinsert[intra]), max(llinsert[intra]), length.out=30)
inward <- hist(llinsert[diags$orientation==1L], plot=FALSE, breaks=breaks)
outward <- hist(llinsert[diags$orientation==2L] ,plot=FALSE, breaks=breaks)
samestr <- hist(llinsert[diags$orientation==0L | diags$orientation==3L], 
   plot=FALSE, breaks=breaks)
samestr$counts <- samestr$counts/2
#
ymax <- max(inward$counts, outward$counts, samestr$counts)/1e6
xmax <- max(inward$mids, outward$mids, samestr$mids)
xmin <- min(inward$mids, outward$mids, samestr$mids)
#
plot(0,0,type="n", xlim=c(xmin, xmax), ylim=c(0, ymax),
    xlab=expression(log[2]~"[insert size (bp)]"), ylab="Frequency (millions)")
lines(inward$mids, inward$counts/1e6, col="darkgreen", lwd=2)
abline(v=log2(min.inward), col="darkgrey")
lines(outward$mids, outward$counts/1e6, col="red", lwd=2)
abline(v=log2(min.outward), col="darkgrey", lty=2)
lines(samestr$mids, samestr$counts/1e6, col="blue", lwd=2)
legend("topright", c("inward", "outward", "same"), 
    col=c("darkgreen", "red", "blue"), lwd=2)


###################################################
### code chunk number 19: diffHicUG.Rnw:377-382
###################################################
prepped <- preparePairs("SRR027958.bam", hs.param, file="SRR027958.h5", 
    dedup=TRUE, minq=10)
counted <- prunePairs("SRR027958.h5", hs.param, file.out="SRR027958_trimmed.h5", 
    max.frag=600, min.inward=min.inward, min.outward=min.outward)
mergePairs(files=c("SRR027957_trimmed.h5", "SRR027958_trimmed.h5"), "merged.h5")


###################################################
### code chunk number 20: diffHicUG.Rnw:390-394
###################################################
anchor.id <- as.integer(runif(100, 1, length(hs.param$fragments)))
target.id <- as.integer(runif(100, 1, length(hs.param$fragments)))
dummy <- data.frame(anchor.id, target.id, other.data=as.integer(runif(100, 1, 100)))
savePairs(dummy, "example.h5", hs.param)


###################################################
### code chunk number 21: diffHicUG.Rnw:429-432
###################################################
require(BSgenome.Mmusculus.UCSC.mm10)
mm.frag <- cutGenome(BSgenome.Mmusculus.UCSC.mm10, "AAGCTT", 4)
input <- c("merged_flox_1.h5", "merged_flox_2.h5", "merged_ko_1.h5", "merged_ko_2.h5")


###################################################
### code chunk number 22: diffHicUG.Rnw:443-450
###################################################
bin.size <- 1e6
mm.param <- pairParam(mm.frag)
data <- squareCounts(input, mm.param, width=bin.size, filter=1)
data
head(counts(data))
head(anchors(data))
head(targets(data))


###################################################
### code chunk number 23: diffHicUG.Rnw:470-471
###################################################
head(regions(data))


###################################################
### code chunk number 24: diffHicUG.Rnw:490-492
###################################################
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
gene.body <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)


###################################################
### code chunk number 25: diffHicUG.Rnw:495-496
###################################################
strand(gene.body) <- "*"


###################################################
### code chunk number 26: diffHicUG.Rnw:504-506
###################################################
redata <- connectCounts(input, mm.param, regions=gene.body)
redata


###################################################
### code chunk number 27: diffHicUG.Rnw:514-515
###################################################
head(regions(redata))


###################################################
### code chunk number 28: diffHicUG.Rnw:532-534
###################################################
margin.data <- marginCounts(input, mm.param, width=bin.size)
margin.data


###################################################
### code chunk number 29: mamargin (eval = FALSE)
###################################################
## adjc <- cpm(asDGEList(margin.data), log=TRUE, prior.count=5)
## smoothScatter(0.5*(adjc[,1]+adjc[,3]), adjc[,1]-adjc[,3],
##     xlab="A", ylab="M", main="Flox (1) vs. Ko (1)")


###################################################
### code chunk number 30: diffHicUG.Rnw:555-556
###################################################
adjc <- cpm(asDGEList(margin.data), log=TRUE, prior.count=5)
smoothScatter(0.5*(adjc[,1]+adjc[,3]), adjc[,1]-adjc[,3],
    xlab="A", ylab="M", main="Flox (1) vs. Ko (1)")


###################################################
### code chunk number 31: diffHicUG.Rnw:580-582
###################################################
new.param <- reform(mm.param, restrict=c("chr1", "chr2"))
new.param


###################################################
### code chunk number 32: diffHicUG.Rnw:589-591
###################################################
new.param <- reform(mm.param, restrict=cbind("chr2", "chr19"))
new.param


###################################################
### code chunk number 33: diffHicUG.Rnw:602-605
###################################################
dummy.repeat <- GRanges("chr1", IRanges(10000, 1000000))
new.param <- reform(mm.param, discard=dummy.repeat)
new.param


###################################################
### code chunk number 34: diffHicUG.Rnw:623-625
###################################################
new.param <- reform(mm.param, cap=5)
new.param


###################################################
### code chunk number 35: avehistplot (eval = FALSE)
###################################################
## require(edgeR)
## ave.ab <- aveLogCPM(asDGEList(data))
## hist(ave.ab, xlab="Average abundance")


###################################################
### code chunk number 36: diffHicUG.Rnw:671-672
###################################################
require(edgeR)
ave.ab <- aveLogCPM(asDGEList(data))
hist(ave.ab, xlab="Average abundance")


###################################################
### code chunk number 37: diffHicUG.Rnw:683-685
###################################################
count.keep <- ave.ab >= aveLogCPM(5, lib.size=mean(data$totals))
summary(count.keep)


###################################################
### code chunk number 38: diffHicUG.Rnw:697-699
###################################################
dummy <- data[count.keep,]
dummy


###################################################
### code chunk number 39: diffHicUG.Rnw:714-719
###################################################
direct <- filterDirect(data)
direct$threshold
direct.keep <- direct$abundances > log2(5) + direct$threshold
dummy <- data[direct.keep,]
summary(direct.keep)


###################################################
### code chunk number 40: diffHicUG.Rnw:726-728
###################################################
direct.keep2 <- direct.keep & count.keep
summary(direct.keep2)


###################################################
### code chunk number 41: diffHicUG.Rnw:746-747
###################################################
trended <- filterTrended(data)


###################################################
### code chunk number 42: avedistplot (eval = FALSE)
###################################################
## smoothScatter(trended$log.distance, trended$abundances, 
##     xlab="Log-Distance", ylab="Normalized abundance")
## o <- order(trended$log.distance)
## lines(trended$log.distance[o], trended$threshold[o], col="red", lwd=2)


###################################################
### code chunk number 43: diffHicUG.Rnw:763-764
###################################################
smoothScatter(trended$log.distance, trended$abundances, 
    xlab="Log-Distance", ylab="Normalized abundance")
o <- order(trended$log.distance)
lines(trended$log.distance[o], trended$threshold[o], col="red", lwd=2)


###################################################
### code chunk number 44: diffHicUG.Rnw:773-776
###################################################
trend.keep <- trended$abundances > trended$threshold
dummy <- data[trend.keep,]
summary(trend.keep)


###################################################
### code chunk number 45: diffHicUG.Rnw:788-790
###################################################
trend.keep2 <- trend.keep & count.keep
summary(trend.keep2)


###################################################
### code chunk number 46: diffHicUG.Rnw:806-809
###################################################
chr1.data <- squareCounts(input, reform(mm.param, restrict="chr1"), width=1e5, filter=1)
enrich.chr1 <- filterPeaks(chr1.data)
summary(enrich.chr1)   


###################################################
### code chunk number 47: diffHicUG.Rnw:829-832
###################################################
chr1.ab <- aveLogCPM(asDGEList(chr1.data))
keep <- enrich.chr1 > 1 & chr1.ab > aveLogCPM(5, lib.size=mean(chr1.data$totals)) 
sum(keep)


###################################################
### code chunk number 48: diffHicUG.Rnw:841-863
###################################################
#which(keep)[order(enrich.chr1[keep], decreasing=TRUE)]
#plotPlaid(input[3], param=mm.param, width=2.5e4, cap=20, col="black",
#	anchor=resize(anchors(chr1.data[chosen]), fix="center", width=1e6),
#	target=resize(targets(chr1.data[chosen]), fix="center", width=1e6))
a.list <- GRanges("chr1", IRanges(
			c(150402687, 179710378),
			c(150501589, 179796841)))
t.list <- GRanges("chr1", IRanges(
			c(139400813, 150402687),
			c(139499958, 150501589)))
capped <- c(20, 10)
par(mfrow=c(1,2))
for (chosen in 1:length(a.list)) { 
	if (!any(anchors(chr1.data[keep])==a.list[chosen] & targets(chr1.data[keep])==t.list[chosen])) { stop("can't find example for peak-calling") }
	plotPlaid(input[3], param=mm.param, width=2.5e4, max.count=capped[chosen], col="black",
			anchor=resize(a.list[chosen], fix="center", width=1e6),
			target=resize(t.list[chosen], fix="center", width=1e6),
			main=paste("Example", chosen))
	box()
	rect(start(a.list[chosen]), start(t.list[chosen]),
		end(a.list[chosen]), end(t.list[chosen]), border="red")
}


###################################################
### code chunk number 49: diffHicUG.Rnw:898-902
###################################################
new.bin.size <- 1e5
smaller.data <- squareCounts(input, mm.param, width=new.bin.size, filter=20)
smaller.ab <- aveLogCPM(asDGEList(smaller.data))
summary(smaller.ab)


###################################################
### code chunk number 50: diffHicUG.Rnw:914-919
###################################################
direct <- filterDirect(data, scale=(bin.size/new.bin.size)^2)
direct$threshold
small.keep <- smaller.ab > direct$threshold + log2(5)
smaller.filtered <- smaller.data[small.keep,]
summary(small.keep)


###################################################
### code chunk number 51: diffHicUG.Rnw:927-932
###################################################
trended <- filterTrended(data, scale=(bin.size/new.bin.size)^2)
s.dist <- log10(getDistance(smaller.data) + new.bin.size)
s.threshold <- approx(x=trended$log.distance, y=trended$threshold, xout=s.dist, rule=2)$y
s.threshold[is.na(s.dist)] <- trended$threshold[is.na(trended$log.distance)][1]
summary(smaller.ab > s.threshold)


###################################################
### code chunk number 52: diffHicUG.Rnw:954-955
###################################################
summary(aveLogCPM(asDGEList(redata)))


###################################################
### code chunk number 53: diffHicUG.Rnw:996-1000
###################################################
dist <- getDistance(data)
diag.keep <- is.na(dist) | dist > 0L  
dummy <- data[diag.keep & direct.keep2,]
summary(diag.keep)


###################################################
### code chunk number 54: diffHicUG.Rnw:1039-1041
###################################################
original <- data
data <- data[direct.keep2,]


###################################################
### code chunk number 55: diffHicUG.Rnw:1049-1050
###################################################
smaller.data <- smaller.filtered


###################################################
### code chunk number 56: diffHicUG.Rnw:1066-1068
###################################################
normfacs <- normalize(data)
normfacs


###################################################
### code chunk number 57: makema (eval = FALSE)
###################################################
## mval <- adj.counts[,3]-adj.counts[,2]
## smoothScatter(ab, mval, xlab="A", ylab="M", main="KO vs. Flox")
## fit <- loessFit(x=ab, y=mval)
## lines(ab[o], fit$fitted[o], col="red")


###################################################
### code chunk number 58: maunnorm (eval = FALSE)
###################################################
## ab <- aveLogCPM(asDGEList(data))
## o <- order(ab)
## adj.counts <- cpm(asDGEList(data), log=TRUE)
## mval <- adj.counts[,3]-adj.counts[,2]
## smoothScatter(ab, mval, xlab="A", ylab="M", main="KO vs. Flox")
## fit <- loessFit(x=ab, y=mval)
## lines(ab[o], fit$fitted[o], col="red")


###################################################
### code chunk number 59: diffHicUG.Rnw:1097-1098
###################################################
ab <- aveLogCPM(asDGEList(data))
o <- order(ab)
adj.counts <- cpm(asDGEList(data), log=TRUE)
mval <- adj.counts[,3]-adj.counts[,2]
smoothScatter(ab, mval, xlab="A", ylab="M", main="KO vs. Flox")
fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")


###################################################
### code chunk number 60: diffHicUG.Rnw:1109-1111
###################################################
nb.off <- normalize(data, type="loess")
head(nb.off)


###################################################
### code chunk number 61: manorm (eval = FALSE)
###################################################
## adj.counts <- log2(counts(data) + 0.5) - nb.off/log(2)
## mval <- adj.counts[,3]-adj.counts[,2]
## smoothScatter(ab, mval, xlab="A", ylab="M", main="KO vs. Flox")
## fit <- loessFit(x=ab, y=mval)
## lines(ab[o], fit$fitted[o], col="red")


###################################################
### code chunk number 62: diffHicUG.Rnw:1127-1128
###################################################
adj.counts <- log2(counts(data) + 0.5) - nb.off/log(2)
mval <- adj.counts[,3]-adj.counts[,2]
smoothScatter(ab, mval, xlab="A", ylab="M", main="KO vs. Flox")
fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")


###################################################
### code chunk number 63: diffHicUG.Rnw:1145-1147
###################################################
corrected <- correctedContact(original, winsor.high=0.02, ignore.low=0.02)
head(corrected$truth)


###################################################
### code chunk number 64: diffHicUG.Rnw:1157-1158
###################################################
corrected$max


###################################################
### code chunk number 65: diffHicUG.Rnw:1229-1231
###################################################
cnv.offs <- normalizeCNV(data, margin.data)
head(cnv.offs)


###################################################
### code chunk number 66: diffHicUG.Rnw:1273-1277
###################################################
count.files <- c("merged_erg.h5", "merged_gfp.h5")
rick.data <- squareCounts(count.files, hs.param, width=1e6)
rick.marg <- marginCounts(count.files, hs.param, width=1e6)
rick.data <- rick.data[aveLogCPM(asDGEList(rick.data)) > 0,]


###################################################
### code chunk number 67: diffHicUG.Rnw:1284-1288
###################################################
matched <- matchMargins(rick.data, rick.marg)
m.adjc <- cpm(asDGEList(rick.marg), log=TRUE)
sum.madjc <- m.adjc[matched$amatch,] + m.adjc[matched$tmatch,]
margin.lr <- sum.madjc[,1] - sum.madjc[,2]


###################################################
### code chunk number 68: cnvplotter (eval = FALSE)
###################################################
## before <- cpm(asDGEList(rick.data), log=TRUE)
## after <- log2(counts(rick.data)+0.5) - normalizeCNV(rick.data, rick.marg, maxk=1000)/log(2)
## par(mfrow=c(1,2), cex.axis=1.2, cex.lab=1.4)
## smoothScatter(margin.lr, before[,1]-before[,2], ylim=c(-4, 4), main="Before",
##     xlab="Sum of marginal log-ratios", ylab="Interaction log-ratio")
## smoothScatter(margin.lr, after[,1]-after[,2], ylim=c(-4, 4), main="After",
##     xlab="Sum of marginal log-ratios", ylab="Interaction log-ratio")


###################################################
### code chunk number 69: diffHicUG.Rnw:1308-1309
###################################################
before <- cpm(asDGEList(rick.data), log=TRUE)
after <- log2(counts(rick.data)+0.5) - normalizeCNV(rick.data, rick.marg, maxk=1000)/log(2)
par(mfrow=c(1,2), cex.axis=1.2, cex.lab=1.4)
smoothScatter(margin.lr, before[,1]-before[,2], ylim=c(-4, 4), main="Before",
    xlab="Sum of marginal log-ratios", ylab="Interaction log-ratio")
smoothScatter(margin.lr, after[,1]-after[,2], ylim=c(-4, 4), main="After",
    xlab="Sum of marginal log-ratios", ylab="Interaction log-ratio")


###################################################
### code chunk number 70: diffHicUG.Rnw:1393-1396
###################################################
design <- model.matrix(~factor(c("flox", "flox", "ko", "ko")))
colnames(design) <- c("Intercept", "KO")
design


###################################################
### code chunk number 71: diffHicUG.Rnw:1405-1407
###################################################
y <- asDGEList(data)
y$offset <- nb.off


###################################################
### code chunk number 72: diffHicUG.Rnw:1415-1417
###################################################
y <- estimateDisp(y, design)
y$common.dispersion


###################################################
### code chunk number 73: bcvplot (eval = FALSE)
###################################################
## plotBCV(y)


###################################################
### code chunk number 74: diffHicUG.Rnw:1430-1431
###################################################
plotBCV(y)


###################################################
### code chunk number 75: diffHicUG.Rnw:1447-1448
###################################################
fit <- glmQLFit(y, design, robust=TRUE)


###################################################
### code chunk number 76: qlplot (eval = FALSE)
###################################################
## plotQLDisp(fit)


###################################################
### code chunk number 77: diffHicUG.Rnw:1462-1463
###################################################
plotQLDisp(fit)


###################################################
### code chunk number 78: diffHicUG.Rnw:1473-1474
###################################################
summary(fit$df.prior)


###################################################
### code chunk number 79: diffHicUG.Rnw:1505-1507
###################################################
result <- glmQLFTest(fit, coef=2)
topTags(result)


###################################################
### code chunk number 80: diffHicUG.Rnw:1529-1531
###################################################
adj.p <- p.adjust(result$table$PValue, method="BH")
sum(adj.p <= 0.05)


###################################################
### code chunk number 81: diffHicUG.Rnw:1537-1544
###################################################
ax <- anchors(data)
tx <- targets(data)
final <- data.frame(anchor.chr=seqnames(ax), anchor.start=start(ax), anchor.end=end(ax),
    target.chr=seqnames(tx), target.start=start(tx), target.end=end(tx), 
    result$table, FDR=adj.p)
o <- order(final$PValue)
write.table(final[o,], file="results.tsv", sep="\t", quote=FALSE, row.names=FALSE)


###################################################
### code chunk number 82: diffHicUG.Rnw:1595-1600
###################################################
y.small <- asDGEList(smaller.data)
y.small$offset <- normalize(smaller.data, type="loess")
y.small <- estimateDisp(y.small, design)
fit.small <- glmQLFit(y.small, design, robust=TRUE)
result.small <- glmQLFTest(fit.small)


###################################################
### code chunk number 83: diffHicUG.Rnw:1608-1612
###################################################
cons <- consolidatePairs(list(larger=data, smaller=smaller.data), 
    list(result$table, result.small$table))
head(cons$id$larger)
head(cons$id$smaller)


###################################################
### code chunk number 84: diffHicUG.Rnw:1622-1624
###################################################
head(cons$table)
sum(cons$table$FDR <= 0.05)


###################################################
### code chunk number 85: diffHicUG.Rnw:1635-1642
###################################################
ax.2 <- anchors(cons$pairs)
tx.2 <- targets(cons$pairs)
final.2 <- data.frame(anchor.chr=seqnames(ax.2), anchor.start=start(ax.2), 
    anchor.end=end(ax.2), target.chr=seqnames(tx.2), target.start=start(tx.2), 
    target.end=end(tx.2), cons$table)
o2 <- order(final.2$PValue)
write.table(final.2[o2,], file="results.2.tsv", sep="\t", quote=FALSE, row.names=FALSE)


###################################################
### code chunk number 86: diffHicUG.Rnw:1658-1670
###################################################
inside <- csaw::getBestTest(cons$id$smaller, result.small$table)
best.vec <- rep(NA, nrow(cons$pairs))
best.vec[as.integer(rownames(inside))] <- inside$best
head(best.vec)
ax.3 <- as.data.frame(anchors(smaller.data))[best.vec,]
tx.3 <- as.data.frame(targets(smaller.data))[best.vec,]
nested <- data.frame(anchor.start=ax.3$start, anchor.end=ax.3$end, 
    target.start=tx.3$start, target.end=tx.3$end, 
    result.small$table[best.vec,c("logFC", "F")])
head(nested)
final.3 <- data.frame(final.2, nest=nested)
write.table(final.3[o2,], file="results.3.tsv", sep="\t", quote=FALSE, row.names=FALSE)


###################################################
### code chunk number 87: plaid1 (eval = FALSE)
###################################################
## plotPlaid(input[1], anchor=expanded.a, target=expanded.t, max.count=bound1, 
##     width=5e4, col="red", param=mm.param, main="Flox")
## rect(start(ax.2[chosen]), start(tx.2[chosen]), end(ax.2[chosen]), end(tx.2[chosen]))


###################################################
### code chunk number 88: plaid3 (eval = FALSE)
###################################################
## plotPlaid(input[3], anchor=expanded.a, target=expanded.t, max.count=bound3, 
##     width=5e4, col="blue", param=mm.param, main="KO")
## rect(start(ax.2[chosen]), start(tx.2[chosen]), end(ax.2[chosen]), end(tx.2[chosen]))


###################################################
### code chunk number 89: setup (eval = FALSE)
###################################################
## chosen <- o2[1]
## expanded.a <- resize(ax.2[chosen], fix="center", width=bin.size*5)
## expanded.t <- resize(tx.2[chosen], fix="center", width=bin.size*5)
## bound1 <- 100
## bound3 <- bound1*data$totals[3]/data$totals[1]


###################################################
### code chunk number 90: diffHicUG.Rnw:1707-1710 (eval = FALSE)
###################################################
## chosen <- o2[1]
## expanded.a <- resize(ax.2[chosen], fix="center", width=bin.size*5)
## expanded.t <- resize(tx.2[chosen], fix="center", width=bin.size*5)
## bound1 <- 100
## bound3 <- bound1*data$totals[3]/data$totals[1]
## plotPlaid(input[1], anchor=expanded.a, target=expanded.t, max.count=bound1, 
##     width=5e4, col="red", param=mm.param, main="Flox")
## rect(start(ax.2[chosen]), start(tx.2[chosen]), end(ax.2[chosen]), end(tx.2[chosen]))
## plotPlaid(input[3], anchor=expanded.a, target=expanded.t, max.count=bound3, 
##     width=5e4, col="blue", param=mm.param, main="KO")
## rect(start(ax.2[chosen]), start(tx.2[chosen]), end(ax.2[chosen]), end(tx.2[chosen]))


###################################################
### code chunk number 91: diffHicUG.Rnw:1713-1718
###################################################
chosen <- o2[1]
expanded.a <- resize(ax.2[chosen], fix="center", width=bin.size*5)
expanded.t <- resize(tx.2[chosen], fix="center", width=bin.size*5)
bound1 <- 100
bound3 <- bound1*data$totals[3]/data$totals[1]
if (as.character(seqnames(ax.2[chosen]))!="chr11" || start(ax.2[chosen])!=30001298 || end(ax.2[chosen])!=31000010 ||
    as.character(seqnames(tx.2[chosen]))!="chr11" || start(tx.2[chosen])!=29000942 || end(tx.2[chosen])!=30001301) {
	warning("first plaid plot hits the wrong spot")
}


###################################################
### code chunk number 92: diffHicUG.Rnw:1728-1729
###################################################
plotPlaid(input[1], anchor=expanded.a, target=expanded.t, max.count=bound1, 
    width=5e4, col="red", param=mm.param, main="Flox")
rect(start(ax.2[chosen]), start(tx.2[chosen]), end(ax.2[chosen]), end(tx.2[chosen]))


###################################################
### code chunk number 93: diffHicUG.Rnw:1731-1732
###################################################
plotPlaid(input[3], anchor=expanded.a, target=expanded.t, max.count=bound3, 
    width=5e4, col="blue", param=mm.param, main="KO")
rect(start(ax.2[chosen]), start(tx.2[chosen]), end(ax.2[chosen]), end(tx.2[chosen]))


###################################################
### code chunk number 94: plaid1b (eval = FALSE)
###################################################
## plotPlaid(input[1], anchor=expanded.a, target=expanded.t, max.count=bound1, 
##     width=5e4, param=mm.param, col="red", main="Flox")
## rect(start(ax[chosen]), start(tx[chosen]), end(ax[chosen]), end(tx[chosen]))


###################################################
### code chunk number 95: plaid3b (eval = FALSE)
###################################################
## plotPlaid(input[3], anchor=expanded.a, target=expanded.t, max.count=bound3, 
##     width=5e4, col="blue", param=mm.param, main="KO")
## rect(start(ax[chosen]), start(tx[chosen]), end(ax[chosen]), end(tx[chosen]))


###################################################
### code chunk number 96: setup2 (eval = FALSE)
###################################################
## chosen <- o[1]
## expanded.a <- resize(ax[chosen], fix="center", width=bin.size*5)
## expanded.t <- resize(tx[chosen], fix="center", width=bin.size*5)
## bound1 <- 30
## bound3 <- bound1*data$totals[3]/data$totals[1]


###################################################
### code chunk number 97: diffHicUG.Rnw:1773-1776 (eval = FALSE)
###################################################
## chosen <- o[1]
## expanded.a <- resize(ax[chosen], fix="center", width=bin.size*5)
## expanded.t <- resize(tx[chosen], fix="center", width=bin.size*5)
## bound1 <- 30
## bound3 <- bound1*data$totals[3]/data$totals[1]
## plotPlaid(input[1], anchor=expanded.a, target=expanded.t, max.count=bound1, 
##     width=5e4, param=mm.param, col="red", main="Flox")
## rect(start(ax[chosen]), start(tx[chosen]), end(ax[chosen]), end(tx[chosen]))
## plotPlaid(input[3], anchor=expanded.a, target=expanded.t, max.count=bound3, 
##     width=5e4, col="blue", param=mm.param, main="KO")
## rect(start(ax[chosen]), start(tx[chosen]), end(ax[chosen]), end(tx[chosen]))


###################################################
### code chunk number 98: diffHicUG.Rnw:1779-1784
###################################################
chosen <- o[1]
expanded.a <- resize(ax[chosen], fix="center", width=bin.size*5)
expanded.t <- resize(tx[chosen], fix="center", width=bin.size*5)
bound1 <- 30
bound3 <- bound1*data$totals[3]/data$totals[1]
if (as.character(seqnames(ax[chosen]))!="chr3" || start(ax[chosen])!=151001261 || end(ax[chosen])!=152002835 ||
    as.character(seqnames(tx[chosen]))!="chr3" || start(tx[chosen])!=145998892 || end(tx[chosen])!=147005847) {
	warning("second plaid plot hits the wrong spot")
}


###################################################
### code chunk number 99: diffHicUG.Rnw:1789-1790
###################################################
plotPlaid(input[1], anchor=expanded.a, target=expanded.t, max.count=bound1, 
    width=5e4, param=mm.param, col="red", main="Flox")
rect(start(ax[chosen]), start(tx[chosen]), end(ax[chosen]), end(tx[chosen]))


###################################################
### code chunk number 100: diffHicUG.Rnw:1792-1793
###################################################
plotPlaid(input[3], anchor=expanded.a, target=expanded.t, max.count=bound3, 
    width=5e4, col="blue", param=mm.param, main="KO")
rect(start(ax[chosen]), start(tx[chosen]), end(ax[chosen]), end(tx[chosen]))


###################################################
### code chunk number 101: rot1 (eval = FALSE)
###################################################
## rotPlaid(input[1], mm.param, region=example, width=2e4, 
##     main="Flox", col="Red", max.count=bound1)


###################################################
### code chunk number 102: rot2 (eval = FALSE)
###################################################
## rotPlaid(input[3], mm.param, region=example, width=2e4,
##     main="KO", col="blue", max.count=bound3)


###################################################
### code chunk number 103: rotset
###################################################
chosen <- o2[3]
example <- tx.2[chosen]
end(example) <- end(ax.2[chosen])
nest.mid.a <- (ax.3$start[chosen]+ax.3$end[chosen])/2
nest.mid.t <- (tx.3$start[chosen]+tx.3$end[chosen])/2
nest.mid <- (nest.mid.a + nest.mid.t)/2
nest.gap <- nest.mid.a - nest.mid.t


###################################################
### code chunk number 104: diffHicUG.Rnw:1826-1830
###################################################
if (as.character(seqnames(ax.2[chosen]))!="chr2" || start(ax.2[chosen])!=136998216 || end(ax.2[chosen])!=138001392 ||
    as.character(seqnames(tx.2[chosen]))!="chr2" || start(tx.2[chosen])!=136001779 || end(tx.2[chosen])!=136998219) {
	warning("third plaid plot hits the wrong spot")
}


###################################################
### code chunk number 105: enbox (eval = FALSE)
###################################################
## points(nest.mid, nest.gap, cex=7) 


###################################################
### code chunk number 106: diffHicUG.Rnw:1837-1842 (eval = FALSE)
###################################################
## chosen <- o2[3]
## example <- tx.2[chosen]
## end(example) <- end(ax.2[chosen])
## nest.mid.a <- (ax.3$start[chosen]+ax.3$end[chosen])/2
## nest.mid.t <- (tx.3$start[chosen]+tx.3$end[chosen])/2
## nest.mid <- (nest.mid.a + nest.mid.t)/2
## nest.gap <- nest.mid.a - nest.mid.t
## rotPlaid(input[1], mm.param, region=example, width=2e4, 
##     main="Flox", col="Red", max.count=bound1)
## points(nest.mid, nest.gap, cex=7) 
## rotPlaid(input[3], mm.param, region=example, width=2e4,
##     main="KO", col="blue", max.count=bound3)
## points(nest.mid, nest.gap, cex=7) 


###################################################
### code chunk number 107: diffHicUG.Rnw:1847-1849
###################################################
rotPlaid(input[1], mm.param, region=example, width=2e4, 
    main="Flox", col="Red", max.count=bound1)
points(nest.mid, nest.gap, cex=7) 


###################################################
### code chunk number 108: diffHicUG.Rnw:1852-1854
###################################################
rotPlaid(input[3], mm.param, region=example, width=2e4,
    main="KO", col="blue", max.count=bound3)
points(nest.mid, nest.gap, cex=7) 


###################################################
### code chunk number 109: diplot (eval = FALSE)
###################################################
## chosen <- o[2]
## expanded.a <- resize(ax[chosen], fix="center", width=5e7)
## expanded.t <- resize(tx[chosen], fix="center", width=5e7)
## colfun <- plotDI(data, result$table$logFC, expanded.a, expanded.t, diag=FALSE)


###################################################
### code chunk number 110: diffHicUG.Rnw:1879-1880
###################################################
chosen <- o[2]
expanded.a <- resize(ax[chosen], fix="center", width=5e7)
expanded.t <- resize(tx[chosen], fix="center", width=5e7)
colfun <- plotDI(data, result$table$logFC, expanded.a, expanded.t, diag=FALSE)


###################################################
### code chunk number 111: diffHicUG.Rnw:1885-1889
###################################################
if (as.character(seqnames(ax[chosen]))!="chr16" || start(ax[chosen])!=70000341 || end(ax[chosen])!=70999295 ||
    as.character(seqnames(tx[chosen]))!="chr16" || start(tx[chosen])!=65999261 || end(tx[chosen])!=67001066) {
	warning("first DI plot hits the wrong spot")
}


###################################################
### code chunk number 112: colorbar (eval = FALSE)
###################################################
## logfc <- -20:20/10
## plot(0,0,type="n", axes=FALSE, xlab="", ylab="", xlim=range(logfc), ylim=c(0,1))
## rect(logfc - 0.05, 0, logfc + 0.05, 1, col=colfun(logfc), border=NA)
## axis(1, cex.axis=1.2)
## mtext("logFC", side=1, line=3, cex=1.4)


###################################################
### code chunk number 113: diffHicUG.Rnw:1912-1913
###################################################
logfc <- -20:20/10
plot(0,0,type="n", axes=FALSE, xlab="", ylab="", xlim=range(logfc), ylim=c(0,1))
rect(logfc - 0.05, 0, logfc + 0.05, 1, col=colfun(logfc), border=NA)
axis(1, cex.axis=1.2)
mtext("logFC", side=1, line=3, cex=1.4)


###################################################
### code chunk number 114: diffHicUG.Rnw:1942-1943
###################################################
sessionInfo()


