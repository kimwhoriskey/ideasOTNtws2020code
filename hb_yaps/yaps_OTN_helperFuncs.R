plotNobs <- function(toa){
	nobs <- apply(toa, 1, function(k) sum(!is.na(k)))
	plot(nobs)
	lines(caTools::runmean(nobs, k=50), col="red")
	lines(caTools::runmean(nobs, k=100), col="blue")
}

getNobs <- function(toa){
	nobs <- apply(toa, 1, function(k) sum(!is.na(k)))
	return(nobs)

}

plotToa <- function(toa){
	best_h <- which.max(apply(toa, 2, function(k) sum(!is.na(k))))
	suppressWarnings(matplot(toa[,best_h], (toa[,best_h] - toa) * 1450))
}



downSampleToa <- function(toa, max_hydro){
	toa_down <- toa
	nobs <- apply(toa, 1, function(k) sum(!is.na(k)))
	# plot(nobs)
	# lines(runmean(nobs, k=100), col="red")

	while(sum(nobs > max_hydro) > 0){
		toa_down[nobs > max_hydro, ] <- toa_down[nobs >  max_hydro] * rbinom(length(toa_down[nobs >  max_hydro]), 1, .95)
		toa_down[toa_down == 0] <- NA
		nobs <- apply(toa_down, 1, function(k) sum(!is.na(k)))
		# plot(nobs)
		# lines(runmean(nobs, k=100), col="red")
	}
	nobs <- apply(toa_down, 1, function(k) sum(!is.na(k)))
	plot(nobs)
	lines(runmean(nobs, k=100), col="red")

	return(toa_down)
}


getSS <- function(temp, sal=0, depth=5){
	ss <- 1410 + 4.21 * temp - 0.037 * temp^2 + 1.10 * sal + 0.018 * depth
	return(ss)
}

buildToaKnownSeq <- function(seq, aligned_dat, hydros){
	toa <- matrix(nrow=length(seq), ncol=nrow(hydros))
	inp <- as.matrix(aligned_dat[,c('seq_ping_idx', 'hydro_idx', 'eposync')])
	toa[inp[, 1:2]] <- inp[,3]
	
	nobs <- apply(toa, 1, function(k) sum(!is.na(k)))
	
	first_det <- head(which(nobs >= 1), n=1)
	last_det <- tail(which(nobs >= 1), n=1)
	toa <- toa[first_det:last_det, ]
	seq <- 	seq[first_det:last_det]

	return(list(toa=toa, seq=seq))
}
