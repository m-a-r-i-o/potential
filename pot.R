workingdir <- "./"

center_var <- function(u)
{	#Centering a variable by median subtraction (you can explore other possibilites...)
	#Used to center velocities and positions to the center of mass (which would be the mean, but still)
	return(u - median(u))
}

cum_mass <- function(snap, epsi)
{	#Calculating the sum of the masses of stars (Mbin) within increasing radii (rbin)
	#epsi determines the size of bins indirectly by selecting the fraction of stars
	#in each bin. The bins are quantiles in radius. Notice that in principle the masses
	#can be different, so we are actually summing them instead of just assuming that all
	#bins have the same mass because they have different numbers.
	#This function is later used to calculate the gravitational field and potential.
	x <- center_var(snap$V1)
	y <- center_var(snap$V2)
	z <- center_var(snap$V3)

	r <- sqrt(x*x + y*y + z*z)
	rq <- seq(from = 0.0, to = 1.0, by = epsi)
	rbin <- quantile(r, rq)
	Mbin <- sapply(rbin, function(rbin) length(r[r < rbin])/length(r))
	return(data.frame(rbin, Mbin))
}

integrate_g <- function(x, y)
{	#trapezoidal method integration of a function y of x (both x and y are arrays of numbers)
	#the array that is returned is shorter (by 1) than x and y
	if((length(x) != length(y)) | (length(x) < 2)) {stop("x and y should have the same length and more than 1 point!")}
	dx <- diff(x)
	yright <- y[-1]
	yleft <- y[-length(y)]
	inte <- cumsum((yright + yleft)*dx/2.0)
	return(inte)
}

potential <- function(snap, epsi)
{	#In the following I refer to potential and force but mean potential per unit mass and acceleration.
	#The function calculates potential per unit mass by integrating minus the gravitational acceleration.
	#Calculate the mean potential by integrating the gravitational acceleration. Remember that the potential
	#is MINUS the work I have to expend to unbind a star from the cluster (bringing it at rest to infinity)
	#so it is the integral of MINUS g (I am pushing up with force g a star that is moving up) and the
	#additive constant has to be forced to 0 at infinity. Lacking infinity, I set it to 0 in the last bin. 
	rMbin <- cum_mass(snap, epsi)
	Mbin <- rMbin$Mbin
	rbin <- rMbin$rbin
	gbin <- -Mbin/(rbin*rbin)

	Ubin <- integrate_g(rbin, -gbin)

	Ubin <- Ubin - Ubin[length(Ubin)]
	Ubin <- c(Ubin, 0.0)
	return(data.frame(rbin, Mbin, gbin, Ubin))
}

energy <- function(snap, epsi)
{	#Calculating the kinetic, potential, and total energy (per unit mass) of stars.
	x <- center_var(snap$V1)
	y <- center_var(snap$V2)
	z <- center_var(snap$V3)
	vx <- center_var(snap$V4)
	vy <- center_var(snap$V5)
	vz <- center_var(snap$V6)
	r <- sqrt(x*x + y*y + z*z)
	potsnap <- potential(snap, epsi)
	U <- approx(x = potsnap$rbin, y = potsnap$Ubin, xout = r, rule = 2)$y
	K <- 0.5*(vx*vx + vy*vy + vz*vz)
	E <- K + U
	return(data.frame(r, U, K, E))	
}

angular_momentum <- function(snap)
{	#Calculating angular momentum of stars (per unit mass)
	x <- center_var(snap$V1)
	y <- center_var(snap$V2)
	z <- center_var(snap$V3)
	vx <- center_var(snap$V4)
	vy <- center_var(snap$V5)
	vz <- center_var(snap$V6)
	r <- sqrt(x*x + y*y + z*z)

	lx <- y*vz - z*vy
	ly <- z*vx - x*vz
	lz <- x*vy - y*vx

	l2 <- lx*lx + ly*ly + lz*lz

	return(data.frame(r, lx, ly, lz, l2))
}

effective_potential <- function(snap, epsi)
{
	potsnap <- potential(snap, epsi)
	angmomsnap <- angular_momentum(snap)
	l2 <- angmomsnap$l2
	r <- angmomsnap$r
	K <- energy(snap, epsi)$K
	radialK <- K - l2/(2.0*r*r)
	Ueffmin <- sapply(l2, function(l2i) min(potsnap$Ubin + 0.5*l2i/(potsnap$rbin*potsnap$rbin)))
	rcirc <- sapply(l2, function(l2i) potsnap$rbin[which.min(potsnap$Ubin + 0.5*l2i/(potsnap$rbin*potsnap$rbin))])
	U <- approx(x = potsnap$rbin, y = potsnap$Ubin, xout = r, rule = 2)$y
	Ueff <- U + 0.5*l2/(r*r)
	rmax <- sapply(1:length(l2), function(i) {
		rright <- potsnap$rbin[(potsnap$rbin > rcirc[i])]
		Uright <- potsnap$Ubin[(potsnap$rbin > rcirc[i])]
		radialE <- radialK[i] + Ueff[i] #K[i] + U[i]
		rmaxbin <- rright[which.min(abs(Uright + 0.5*l2[i]/(rright*rright) - radialE))]
		return(rmaxbin)
		})
	rmin <- sapply(1:length(l2), function(i) {
		rleft <- potsnap$rbin[!(potsnap$rbin > rcirc[i])]
		Uleft <- potsnap$Ubin[!(potsnap$rbin > rcirc[i])]
		radialE <- radialK[i] + Ueff[i] #K[i] + U[i]
		rminbin <- rleft[which.min(abs(Uleft + 0.5*l2[i]/(rleft*rleft) - radialE))]
		return(rminbin)
		})

	return(data.frame(r, U, l2, Ueff, Ueffmin, rcirc, radialK, rmax, rmin))
	#r is the radius of an actual star, U its actual potential, l2 its specific angular momentum squared
	#Ueff the current value of the associated effective potential (must always be greater than Ueffmin)
	#Ueffmin the minimum of the associated effective potential (the value for which motion would be circular)
	#rcirc is the radius associated to Ueffmin for a given particle (its circular motion radius that would
	#be reached if as much energy as possible were removed while keeping angular momentum constant)
}

plotpot <- function(snap, epsi, col)
{
	potsnap <- potential(snap, epsi)
	lines(potsnap$rbin, potsnap$Ubin, col = col)
}

plotthissnap <- function(snapnum, epsi)
{	#This function makes a bunch of plots for each snapshot
	snapname <- paste(workingdir, "c_", snapnum, ".dat", sep = "")
	snap <- read.table(snapname, header = F)
	epposnap <- effective_potential(snap, epsi)
	enesnap <- energy(snap, epsi)
	pdfname <- paste(snapnum, ".pdf", sep = "")
	pdf(pdfname)

	#Actual radius of stars, versus circular radius (the radius at the minimum of the effective potential)
	#The actual radius can be larger or smaller than the circular radius, a 1:1 line is plotted for guidance
	limz <- c(0.0, 3.0)
	plot(epposnap$r, epposnap$rcirc, xlim = limz, ylim = limz, pch = ".", xlab = "r", ylab = expression(r[c*i*r*c]))
	lines(limz, limz)

	#Energy (total: kinetic + potential) of stars as a function of radius. Energy should be negative, except
	#for very few stars which are going to leave the cluster unless they are scattered back by other stars
	#It can be seen that for a given radius there is a minimum energy to be there. The stars that have zero
	#kinetic energy occupy the boundary region (i.e. have the minimum energy required to be at a given radius)
	#They should be on very elongated (high eccentricity) orbits, exactly linear orbits if exactly K=0.
	#Here I plot the low 5% of stars in kinetic energy (stars slower than the remaining 95%) in orange, so they
	#trace the boundary. Also all stars on the boundary should be orange.
	plot(enesnap$r, enesnap$E, pch = ".", xlab = "r", ylab = "E")
	lowK <- (enesnap$K < quantile(enesnap$K, 0.05))
	points(enesnap$r[lowK], enesnap$E[lowK], pch = ".", col = "orange")

	#Potential energy as a function of kinetic energy. The K + U = 0 line (U = - K line) represents stars
	#that are unbound from the cluster. There should be next to no stars above this line. The "virial"
	#2K + U = 0 line is also plotted, and it should roughly go through the middle of the stars.
	plot(enesnap$U, enesnap$K, pch = ".", xlab = "U", ylab = "K")
	rUra <- range(enesnap$U)
	lines(rUra, -0.5*rUra, lty = 2)
	lines(rUra, -rUra)

	#Ratio between kinetic and potential energy (sometimes called virial ratio or double/half the virial ratio)
	#as a function of radius. Plotting window trimmed due to extreme outliers.
	plot(enesnap$r, enesnap$K/enesnap$U, pch = ".", xlim = c(0.0, quantile(enesnap$r, 0.95)), ylim = c(-1.0, 0.0), xlab = "r", ylab = "K/U")

	#Distribution of energy (total energy per unit mass) obtained by KDE			
	plot(density(enesnap$E), main = "", xlab = "E")

	#Effective potential energy per unit mass (potential energy in the equivalent one-dimensional problem)
	#versus kinetic energy of the radial motion. The potential energy is shiftd so its minimum is 0.
	#Note that the effective potential is different for every star due to its different angular momentum
	plot(epposnap$Ueff - epposnap$Ueffmin, enesnap$K - epposnap$l2/(2.0*epposnap$r*epposnap$r), pch = ".", xlab = expression(U[e*f*f]), ylab = expression(K[r]))

	#Maximum radius attainable by a star on an orbit in the mean potential versus actual radius.
	#This is calculated by computing the effective potential (mean potential corrected with the centrifugal term,
	#turning the two-dimensional problem into a one-dimensional one, read the Landau 1) and finding the
	#points where the kinetic energy would be 0 in the one-dimensional problem. The outermost of the two is rmax.
	#For all stars r < rmax, because the actual radius is smaller than the maximum radius (always by construction)
	#in the circular potential problem, but in practice we can see some (few!) r > rmax because of approximations
	#in calculating the mean potential (which is obtained by binning and numerical integration) or because the
	#mean potential approximation breaks down (a star is interacting with another star or with an overdensity) 
	plot(epposnap$r, epposnap$rmax, pch = ".", xlab = "r", ylab = expression(r[m*a*x]))
	lines(c(0, 100), c(0, 100))

	#Now rmax is plotted against the radius at the minimum effective potential, the radius that the orbit would
	#have if its energy were reduced, at constant angular momentum, until the orbit circularizes. The maximum
	#radius is always larger than the circular radius, and when the two are similar, the orbit is almost circular
	plot(epposnap$rcirc, epposnap$rmax, pch = ".", xlab = expression(r[c*i*r*c]), ylab = expression(r[m*a*x]))
	lines(c(0, 100), c(0, 100))

	#Same as above but for the minimum radius. r > rmin always except for breakdowns... like before
	plot(epposnap$r, epposnap$rmin, pch = ".", xlab = "r", ylab = expression(r[m*i*n]))
	lines(c(0, 100), c(0, 100))

	#Same as above but for the minimum radius. rcirc > rmin always except for breakdowns... like before	
	plot(epposnap$rcirc, epposnap$rmin, pch = ".", xlab = expression(r[c*i*r*c]), ylab = expression(r[m*i*n]))
	lines(c(0, 100), c(0, 100))

	#Effective eccentricity plotted as a function of the circular radius. This would be the actual eccentricity
	#in a Keplerian motion, as computed using the perihelion and the aphelion. For this non-keplerian potential
	#the effective eccentricity still measures the shape of the orbit, in the sense that it is 0 if the orbit is
	#circular, and tends to 1 if the ratio between rmax and rmin tends to infinity. In the case of radial
	#anisotropy, this should increase with radius (it does not seem to do this)
	ecce <- (epposnap$rmax - epposnap$rmin)/(epposnap$rmax + epposnap$rmin)
	plot(epposnap$rcirc, ecce,  pch = ".", xlab = expression(r[c*i*r*c]), ylab = expression(e[e*f*f]))

	#Energy as a function of angular momentum squared (both per unit mass). If you want to have a given
	#angular momentum, there is a minimum energy that you must have (that energy corresponds to the circular
	#orbit, with more eccentric orbits having more energy for the same angular momentum). So we plot in orange
	#the orbits that are in the lowest 5 percentiles of eccentricity (those that are almost circular), and they
	#are actually along the boundary.
	plot(epposnap$l2, enesnap$E, pch = ".", xlab = expression(L^2), ylab = "E")
	almostcirc <- ecce < quantile(ecce, 0.05)
	points(epposnap$l2[almostcirc], enesnap$E[almostcirc], pch = ".", col = "orange")
	dev.off()

	#Return something so that we can make plots as a function of time with it
	#return(epposnap$rcirc[1:2000]) #have to limit it, otherwise the length is different
	qrcirc <- seq(from = 0.05, to = 0.95, by = 0.05)
	#return(quantile(epposnap$rcirc, qrcirc))
	return(quantile(epposnap$rcirc, qrcirc))
}

eccetime <- function(snapnum, epsi)
{
	snapname <- paste(workingdir, "c_", snapnum, ".dat", sep = "")
	snap <- read.table(snapname, header = F)
	epposnap <- effective_potential(snap, epsi)
	ecce <- (epposnap$rmax - epposnap$rmin)/(epposnap$rmax + epposnap$rmin)
	return(quantile(ecce, seq(from = 0.0, to = 1.0, by = 0.02)))
}

rcirctime <- function(snapnum, epsi)
{
	snapname <- paste(workingdir, "c_", snapnum, ".dat", sep = "")
	snap <- read.table(snapname, header = F)
	epposnap <- effective_potential(snap, epsi)
	return(quantile(epposnap$rcirc, seq(from = 0.0, to = 1.0, by = 0.02)))
}

radeccemass <- function(snapnum, epsi)
{
	snapname <- paste(workingdir, "c_", snapnum, ".dat", sep = "")
	snap <- read.table(snapname, header = F)
	m <- snap$V7
	epposnap <- effective_potential(snap, epsi)
	ecce <- (epposnap$rmax - epposnap$rmin)/(epposnap$rmax + epposnap$rmin)
	rcirc <- epposnap$rcirc
	rmin <- epposnap$rmin
	rmax <- epposnap$rmax
	return(data.frame(m, ecce, rcirc, rmin, rmax))
}

#snappi <- c("0000", "0100", "0200", "0300", "0400", "0500", "0600", "0700", "0800", "0900", "1000", "1100", "1200", "1300", "1400", "1500", "1600", "1700", "1800", "1900")
snappi <- c("0000", "0300", "0600", "0900", "1200", "1500", "1800")

plomat <- sapply(snappi, function(snapnum) plotthissnap(snapnum, 0.0001))

plomat

q()

time <- (1:ncol(plomat))-1
pdf("rcirctime.pdf")
plot(range(time),c(0, 1), type = "n", xlab = "t", ylab = "eccentricity")
sapply(1:nrow(plomat), function(i) lines(time, plomat[i,])) 
dev.off()

#radeccemass("0000", 0.01)

warnings()

