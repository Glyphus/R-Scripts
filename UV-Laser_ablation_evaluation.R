#segmented package for the piecewise linear regression
library("segmented")
library("zoo")
library("strucchange")
library("outliers")
library("MASS")
library("robustbase")

killbreaks <- function(breaks) {
	if ((breaks[2]-breaks[1]) < 30) {
		breaks <- breaks[-2];
		breaks <- killbreaks(breaks);
	}
	return(breaks)
}

#user selects the parent directory for the files
dir <- choose.dir(default = "", caption = "Select folder")
#looking for all csv files from knime
name.files <- list.files(path = dir, pattern = "csv", all.files = FALSE,
           full.names = FALSE, recursive = TRUE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
k <- length(name.files)
#starting a loop that loads each csv file after another and starts the regression
for (i in 1:k) {
	try({
		#loading the data
		df.cut <- read.table(paste(dir,"\\",name.files[i[1]],sep=""), sep=";", header=TRUE);
		offset <- df.cut[1,1]-1
		df.cut <- data.frame((df.cut$x-df.cut[1,1]+1),(df.cut$Y));
		df.cut
		names(df.cut) <- c("x", "Y");
		df.cut$Y <- c(head(df.cut$Y, 2),rollmean(df.cut$Y,5),tail(df.cut$Y, 2))
		#generate ln of x
		df.Logx <- data.frame(log(df.cut[1]),df.cut[2]);
		bp <- breakpoints(Y ~ x, data = df.Logx, breaks = 4);
		reg <- segmented.lm(obj = lm(Y ~ x, data = df.Logx),
			seg.Z = ~ x, psi = NA,
			seg.control (toll = 1e-5, it.max = 100, n.boot = 100, h = 1, K = length(bp$breakpoints), display = FALSE));
		#saving breakpoints to Variable
		breaks <- exp(reg$psi[(length(bp$breakpoints)+1):(length(bp$breakpoints)*2)])
		#breaks <- breaks[which(breaks < (145-offset) | breaks > (161-offset))]
		if (breaks[1] < 25) {
			breaks <- killbreaks(breaks);
		}
		span <- max(df.cut$Y)-min(df.cut$Y[1:100])
		#regression of the whole graph
		span.checks <- c()
		##for	(a in 1:length(breaks)) {
		##		if (((df.cut$Y[which(abs(df.cut$x-breaks[a])==min(abs(df.cut$x-breaks[a])))]-min(df.cut$Y[1:100])) < 0.85*span) & ((df.cut$Y[which(abs(df.cut$x-breaks[a])==min(abs(df.cut$x-breaks[a])))]-min(df.cut$Y[1:100])) > 0.15*span) ) {
		##			span.checks <- c(span.checks, a);
		##		};
		##	};
		##	if (length(span.checks) > 0) {
		##		breaks <- breaks[-span.checks];
		##}
		for (n in 1:length(breaks)) {
			nam <-  paste ("X",n-1, sep = "");
			assign(nam, breaks[n]);
		};
		if (breaks[1] < 40) {
			if (length(breaks) == 4) {
				kinetics <- nls (Y ~ (x < X0)* 	(N0) +
					(x >= X0 & x < X1)* 		(N0 + N1*(1 - exp(-kon * (x-X0)))) +
					(x >= X1 & x < X2)* 		(N0 + N1) +
					(x >= X2 & x < X3)*		(N2 + (N0 + N1 - N2)*(exp(-koff * (x-X2)))) +
					(x >= X3)* 				(N2),
					start = list (N0 = head(df.cut$Y,1), N1 = max(df.cut$Y), N2 = tail(df.cut$Y,1), kon = 0.025, koff = 0.004),
					data = df.cut,
					control=nls.control(minFactor=1e-5, maxiter=10000, warnOnly = TRUE));
				Structure <- c("TRUE","TRUE","TRUE");
				Structure <- setNames(Structure, c("Lag","Plateau","Tail"));
			};
			if (length(breaks) == 2) {
				kinetics <- nls (Y ~ (x < X0)* 	(N0) +
					(x >= X0 & x < X1)* 		(N0 + N1*(1 - exp(-kon * (x-X0)))) +
					(x >= X1 )*				(N2 + (N0 + N1 - N2)*(exp(-koff * (x-X1)))),
					start = list (N0 = head(df.cut$Y,1), N1 = max(df.cut$Y), N2 = tail(df.cut$Y,1), kon = 0.025, koff = 0.004),
					data = df.cut,
					control=nls.control(minFactor=1e-5, maxiter=10000, warnOnly = TRUE));
				Structure <- c("TRUE","FALSE","FALSE");
				Structure <- setNames(Structure, c("Lag","Plateau","Tail"));
			};
			if (length(breaks) == 3) {
				if ((mean(df.cut$Y[(which(abs(df.cut$x-breaks[2])==min(abs(df.cut$x-breaks[2])))):(which(abs(df.cut$x-breaks[3])==min(abs(df.cut$x-breaks[3]))))])-min(df.cut$Y[1:100])) < 0.85*span) {
					kinetics <- nls (Y ~ (x < X0)* 	(N0) +
						(x >= X0 & x < X1)* 		(N0 + N1*(1 - exp(-kon * (x-X0)))) +
						(x >= X1 & x < X2)*		(N2 + (N0 + N1 - N2)*(exp(-koff * (x-X1)))) +
						(x >= X2)* 				(N2),
						start = list (N0 = head(df.cut$Y,1), N1 = max(df.cut$Y), N2 = tail(df.cut$Y,1), kon = 0.025, koff = 0.004),
						data = df.cut,
						control=nls.control(minFactor=1e-5, maxiter=10000, warnOnly = TRUE));
					Structure <- c("TRUE","FALSE","TRUE");
					Structure <- setNames(Structure, c("Lag","Plateau","Tail"));
				} else {
					kinetics <- nls (Y ~ (x < X0)* 	(N0) +
						(x >= X0 & x < X1)* 		(N0 + N1*(1 - exp(-kon * (x-X0)))) +
						(x >= X1 & x < X2)* 		(N0 + N1) +
						(x >= X2)*				(N2 + (N0 + N1 - N2)*(exp(-koff * (x-X2)))),
						start = list (N0 = head(df.cut$Y,1), N1 = max(df.cut$Y), N2 = tail(df.cut$Y,1), kon = 0.025, koff = 0.004),
						data = df.cut,
						control=nls.control(minFactor=1e-5, maxiter=10000, warnOnly = TRUE));
					Structure <- c("TRUE","TRUE","FALSE");
					Structure <- setNames(Structure, c("Lag","Plateau","Tail"));
				};
			};
		} else {
			if (length(breaks) == 3) {
				kinetics <- nls (Y ~ (x < X0)* 	(N0 + N1*(1 - exp(-kon * (x)))) +
					(x >= X0 & x < X1)* 		(N0 + N1) +
					(x >= X1 & x < X2)* 		(N2 + (N0 + N1 - N2)*(exp(-koff * (x-X1)))) +
					(x >= X2)*				(N2),
					start = list (N0 = head(df.cut$Y,1), N1 = max(df.cut$Y), N2 = tail(df.cut$Y,1), kon = 0.025, koff = 0.004),
					data = df.cut,
					control=nls.control(minFactor=1e-5, maxiter=10000, warnOnly = TRUE));
				Structure <- c("FALSE","TRUE","TRUE");
				Structure <- setNames(Structure, c("Lag","Plateau","Tail"));
			};
			if (length(breaks) == 1) {
				kinetics <- nls (Y ~ (x < X0)* 	(N0 + N1*(1 - exp(-kon * (x)))) +
					(x >= X0)*		 		(N2 + (N0 + N1 - N2)*(exp(-koff * (x-X0)))),
					start = list (N0 = head(df.cut$Y,1), N1 = max(df.cut$Y), N2 = tail(df.cut$Y,1), kon = 0.025, koff = 0.004),
					data = df.cut,
					control=nls.control(minFactor=1e-5, maxiter=10000, warnOnly = TRUE));
				Structure <- c("FALSE","FALSE","FALSE");
				Structure <- setNames(Structure, c("Lag","Plateau","Tail"));
			};
			if (length(breaks) == 2) {
				if ((mean(df.cut$Y[(which(abs(df.cut$x-breaks[1])==min(abs(df.cut$x-breaks[1])))):(which(abs(df.cut$x-breaks[2])==min(abs(df.cut$x-breaks[2]))))])-min(df.cut$Y)) < 0.85*span) {
						kinetics <- nls (Y ~ (x < X0)* 	(N0 + N1*(1 - exp(-kon * (x)))) +
						(x >= X0 & x < X1)* 		(N2 + (N0 + N1 - N2)*(exp(-koff * (x-X0)))) +
						(x >= X1)*				(N2),
						start = list (N0 = head(df.cut$Y,1), N1 = max(df.cut$Y), N2 = tail(df.cut$Y,1), kon = 0.025, koff = 0.004),
						data = df.cut,
						control=nls.control(minFactor=1e-5, maxiter=10000, warnOnly = TRUE));
				Structure <- c("FALSE","FALSE","TRUE");
				Structure <- setNames(Structure, c("Lag","Plateau","Tail"));
				} else {
					kinetics <- nls (Y ~ (x < X0)* 	(N0 + N1*(1 - exp(-kon * (x)))) +
						(x >= X0 & x < X1)* 		(N0 + N1) +
						(x >= X1)*		 		(N2 + (N0 + N1 - N2)*(exp(-koff * (x-X1)))),
						start = list (N0 = head(df.cut$Y,1), N1 = max(df.cut$Y), N2 = tail(df.cut$Y,1), kon = 0.025, koff = 0.004),
						data = df.cut,
						control=nls.control(minFactor=1e-5, maxiter=10000, warnOnly = TRUE));
					Structure <- c("FALSE","TRUE","FALSE");
					Structure <- setNames(Structure, c("Lag","Plateau","Tail"));
				};
			};
		};
		#saving the relevant parameters
		names <- c()
		for (n in 1:length(breaks)) {
			names <-  c(names,paste ("X",n-1, sep = ""));
		};
		names <- c(names, "R2");
		breaks <- c(breaks, (1-(kinetics$m$deviance()/(sum((df.cut$Y-mean(df.cut$Y))^2)))))
		out.breaks <- setNames(breaks,names);
		output <- c(coefficients(kinetics), out.breaks,Structure)
		write.table(output, paste(dir,"\\",name.files[i[1]],".txt",sep=""), sep="\t")
		#saving a plot of the fitted curves and raw data
		jpeg(paste(dir,"\\",name.files[i[1]],".jpg",sep=""));
		plot(df.cut$x,kinetics$m$fitted(),xlab = "time in sec", ylab = "mean bg corrected intensity");
		points(df.cut$x,df.cut$Y,col="red");
		dev.off();
	});
}
