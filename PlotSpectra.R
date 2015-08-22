polyCurve <- function(data, from, to, n = 5000, miny=0, col = "red", border = col, lty=1)
{
     x <- data[ ,1]
     y <- data[ ,2]
     drawPoly <- function(fun, from, to, n = 5000, miny=0, col, border, lty) {
          Sq <- seq(from = from, to = to, length = n)
          polygon(x = c(Sq[1], Sq, Sq[n]),
                  y = c(miny, fun(Sq), miny),
                  col = col, border = border, lty = lty)
     }
     lf <- length(from)
     stopifnot(identical(lf, length(to)))
     if(length(col) != lf)
          col <- rep(col, length.out = lf)
     if(length(border) != lf)
          border <- rep(border, length.out = lf)
     if(missing(miny))
          miny <- min(y)
     interp <- approxfun(x = x, y = y)
     mapply(drawPoly, from = from, to = to, col = col, border = border, lty = lty,
            MoreArgs = list(fun = interp, n = n, miny = miny))
     invisible()
}


ex <- list()
# Excitation
ex$t390 <- read.table(file='~/Public/DropBox/GitHub/R-Microscopy/Spectra/390.txt', header=TRUE)
#ex$t438 <- read.table(file='~/Public/DropBox/GitHub/R-Microscopy/Spectra/438.txt', header=TRUE)
ex$t485 <- read.table(file='~/Public/DropBox/GitHub/R-Microscopy/Spectra/485.txt', header=TRUE)
ex$t560 <- read.table(file='~/Public/DropBox/GitHub/R-Microscopy/Spectra/560.txt', header=TRUE)
#ex$t585 <- read.table(file='~/Public/DropBox/GitHub/R-Microscopy/Spectra/585.txt', header=TRUE)
ex$t648 <- read.table(file='~/Public/DropBox/GitHub/R-Microscopy/Spectra/648.txt', header=TRUE)
ex$t740 <- read.table(file='~/Public/DropBox/GitHub/R-Microscopy/Spectra/740.txt', header=TRUE)

# Emission
em <- list()
em$t440 <- read.table(file='~/Public/DropBox/GitHub/R-Microscopy/Spectra/440.txt', header=TRUE)
em$t440$X.T <- em$t440$X.T * 100
em$t525 <- read.table(file='~/Public/DropBox/GitHub/R-Microscopy/Spectra/525.txt', header=TRUE)
em$t525$X.T <- em$t525$X.T * 100
em$t607 <- read.table(file='~/Public/DropBox/GitHub/R-Microscopy/Spectra/607.txt', header=TRUE)
em$t607$X.T <- em$t607$X.T * 100
em$t684 <- read.table(file='~/Public/DropBox/GitHub/R-Microscopy/Spectra/684.txt', header=TRUE)
em$t684$X.T <- em$t684$X.T * 100
em$t809 <- read.table(file='~/Public/DropBox/GitHub/R-Microscopy/Spectra/809.txt', header=TRUE)
em$t809$X.T <- em$t809$X.T * 100

flNames <- list(Alexa488 = 'Alexa 488',
Alexa555 = 'Alexa 555',
Alexa647 = 'Alexa 647',
Alexa750 = 'Alexa 750',
Calcein = 'Calcein',
CellTrackerRed = 'Cell Tracker Red CMTPX',
Cy3 = 'Cy3',
Cy5.5 = 'Cy5.5',
Cy5 = 'Cy5',
Cy7 = 'Cy7',
DAPI = 'DAPI',
DiI = 'DiI',
DyLight680 = 'DyLight 680',
DyLight800 = 'DyLight 800',
eGFP = 'eGFP',
Eosin = 'Eosin',
EtBr = 'Ethidium Bromide',
EtHD = 'Ethidium Homodimer',
FITC = 'FITC',
Hoechst = 'Hoechst 33342',
IRDye680 = 'LI-CORE IRDye 680',
IRDye800 = 'LI-CORE IRDye 800',
mCherry = 'mCherry',
MitoFarRed = 'MitoTracker Far Red',
MitoGreen = 'MitoTracker Green',
MitoOrange = 'MitoTracker Orange',
MitoRed = 'MitoTracker Red',
PE = 'R-Phycoerythrin (PE)',
PI = 'Propidium Iodide',
RFP = 'RFP',
SYBR = 'SYBR Green',
TRITC = 'TRITC',
TxRed = 'Texas Red')

fl <- list()
for(name in names(flNames))
{
     fl[[paste0(name, 'EX')]] <- read.table(file=paste0('~/Public/DropBox/GitHub/R-Microscopy/Spectra/', name, 'EX.txt'), header=TRUE)
     if(name != 'Hoechst' && name != 'CellTrackerRed')
     {
          fl[[paste0(name, 'EX')]]$X.T <- fl[[paste0(name, 'EX')]]$X.T * 100
     }
     lastWL <- fl[[paste0(name, 'EX')]]$wl[nrow(fl[[paste0(name, 'EX')]])]
     fl[[paste0(name, 'EX')]] <- rbind(data.frame(wl=fl[[paste0(name, 'EX')]]$wl[1]-1.5, X.T=0), fl[[paste0(name, 'EX')]])
     fl[[paste0(name, 'EX')]] <- rbind(fl[[paste0(name, 'EX')]], data.frame(wl=lastWL+1.5, X.T=0))
     fl[[paste0(name, 'EX')]] <- rbind(data.frame(wl=fl[[paste0(name, 'EX')]]$wl[1]-2.5, X.T=0), fl[[paste0(name, 'EX')]])
     fl[[paste0(name, 'EX')]] <- rbind(fl[[paste0(name, 'EX')]], data.frame(wl=lastWL+2.5, X.T=0))
     fl[[paste0(name, 'EM')]] <- read.table(file=paste0('~/Public/DropBox/GitHub/R-Microscopy/Spectra/', name, 'EM.txt'), header=TRUE)
     if(name != 'Hoechst' && name != 'CellTrackerRed')
     {
          fl[[paste0(name, 'EM')]]$X.T <- fl[[paste0(name, 'EM')]]$X.T * 100
     }
     lastWL <- fl[[paste0(name, 'EM')]]$wl[nrow(fl[[paste0(name, 'EM')]])]
     fl[[paste0(name, 'EM')]] <- rbind(data.frame(wl=fl[[paste0(name, 'EM')]]$wl[1]-1.5, X.T=0), fl[[paste0(name, 'EM')]])
     fl[[paste0(name, 'EM')]] <- rbind(fl[[paste0(name, 'EM')]], data.frame(wl=lastWL+1.5, X.T=0))
     fl[[paste0(name, 'EM')]] <- rbind(data.frame(wl=fl[[paste0(name, 'EM')]]$wl[1]-2.5, X.T=0), fl[[paste0(name, 'EM')]])
     fl[[paste0(name, 'EM')]] <- rbind(fl[[paste0(name, 'EM')]], data.frame(wl=lastWL+2.5, X.T=0))
}

for(name in names(flNames))
{
     pdf(file=paste0('~/Public/DropBox/GitHub/R-Microscopy/Spectra/Plots/',name,'.pdf'), width=9, height=6.5)
     plot(c(), c(), xlim=c(350,850), ylim=c(0,100), xlab='Wavelength [nm]', ylab='% Transmission', main=flNames[[name]], cex=1.5)
     polyCurve(fl[[paste0(name, 'EX')]], from=350, to=1200, n=1000, miny=0, col=rgb(0,0,0,0.15), bor=gray(0.5), lty=1)
     polyCurve(fl[[paste0(name, 'EM')]], from=350, to=1200, n=1000, miny=0, col=rgb(0,0,0,0.3), bor='black', lty=1)
     for(x in names(ex))
     {
          lines(ex[[x]]$wl, ex[[x]]$X.T, col=gray(0.5), lty=1)
     }
     for(m in names(em))
     {
          lines(em[[m]]$wl, em[[m]]$X.T, col='black', lty=1)
     }
     dev.off()
}

# Dichroic
tx$cube <- read.table(file='~/Public/DropBox/GitHub/R-Microscopy/Spectra/cube.txt', header=TRUE)
tx$cube$X.T <- tx$cube$X.T * 100

# Fluorophores
tx$DAPI <-

plotSpectra <- function(ex, em, name='t390')
par(mar=c(4.1,4.1,2.5,0.5))
plot(t390$wl,t390$X.T, type='l', ylab='% Transmission', xlab='Wavelength [nm]', main='390 X 440')
polyCurve(em$t390,from=200,to=1200, n=nrow(em$t390), col='red')

duh <- read.csv(file='~/Desktop/Scope Data/Hoechst.csv', header=TRUE)
duhEX <- duh[ ,c(1,2)]
names(duhEX) <- c('wl','X.T')
write.table(duhEX, file='~/Desktop/Scope Data/HoechstEX.txt', row.names=F)
duhEM <- duh[ ,c(1,3)]
names(duhEM) <- c('wl','X.T')
write.table(duhEM, file='~/Desktop/Scope Data/HoechstEM.txt', row.names=F)

duh <- read.csv(file='~/Desktop/Scope Data/CellTrackerRedCMTPX.csv', header=TRUE)
duhEX <- duh[ ,c(1,2)]
names(duhEX) <- c('wl','X.T')
write.table(duhEX, file='~/Desktop/Scope Data/CellTrackerRedEX.txt', row.names=F)
duhEM <- duh[ ,c(1,3)]
names(duhEM) <- c('wl','X.T')
write.table(duhEM, file='~/Desktop/Scope Data/CellTrackerRedEM.txt', row.names=F)
