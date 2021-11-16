#!/applications/R/R-4.0.0/bin/Rscript

xs <- seq(from = -2*pi, to = 2*pi, by = pi/100)
wave.1 <- sin(3*xs)
wave.2 <- sin(10*xs)

wave.3 <- 0.5 * wave.1 + 0.25 * wave.2

# Plot trajectory given a Fourier series
plot.fourier <- function(fourier.series, f.0, ts) {
  w <- 2*pi*f.0
  trajectory <- sapply(ts, function(t) fourier.series(t,w))
  plot(ts, trajectory, type = "l", xlab = "time", ylab = "f(t)")
  abline(h = 0, lty = 3)
}

pdf("plots/fs_trajectory.pdf")
plot.fourier(fourier.series = function(t,w) { sin(w*t) }, f.0 = 1, ts = seq(0, 1, 1/100))
dev.off()


# Plot the signal for equation f(t) = 0.5 * sin(3wt) + 0.25 * sin(10wt)
# data aquisition frequency (Hz)
acq.freq <- 100
# measuring time interval (seconds)
time <- 6
# vector of sampling time points (s)
ts <- seq(0, time, 1/acq.freq)
# fundamental frequency
f.0 <- 1/time

dc.component <- 0
# frequency of signal components (Hz)
component.freqs <- c(3, 10)
# delay of signal components (radians)
component.delay <- c(0, 0)
# strength of signal components
component.strength <- c(.5, .25)

f <- function(t,w) {
  dc.component +
  sum( component.strength * sin(component.freqs*w*t + component.delay))
}

pdf("plots/signal.pdf")
plot.fourier(f, f.0, ts)
dev.off()

# Phase shifts - translations over the x-axis; "phase angle" or "starting point"
component.delay <- c(pi/2,0)
pdf("plots/signal_phase_shift_wave.1_by_halfpi.pdf")
plot.fourier(f, f.0, ts)
dev.off()

# DC components - translations over the y-axis; constant amplitude for the 0Hz cycle
dc.component <- -2
pdf("plots/signal_dc_component_constant_amp_minus2.pdf")
plot.fourier(f, f.0, ts)
dev.off()

fft(c(1, 1, 1, 1) / 4)
fft(1:4) / 4

# Function to convert the fft() output to real numbers
# cs is the vector of complex points to convert
convert.fft <- function(cs, sample.rate=1) {
  cs <- cs / length(cs) # normalize

  distance.center <- function(c) signif( Mod(c), 4)
  angle           <- function(c) signif( 180*Arg(c)/pi, 3)

  df <- data.frame(cycle    = 0:(length(cs)-1),
                   freq     = 0:(length(cs)-1) * sample.rate / length(cs),
                   strength = sapply(cs, distance.center),
                   delay    = sapply(cs, angle)) 

  df
}

convert.fft(fft(1:4))

# Function to plot a frequency spectrum of a given fft()-transformed trajectory or signal 
plot.frequency.spectrum <- function(X.k, xlimits = c(0, length(X.k))) {
  plot.data <- cbind(0:(length(X.k)-1), Mod(X.k))

  # TODO: why is this scaling necessary?
  plot.data[2:length(X.k), 2] <- 2 * plot.data[2:length(X.k), 2]

  plot(plot.data, t = "h", lwd = 2, main = "",
       xlab = "Frequency (Hz)", ylab = "Strength",
       xlim = xlimits, ylim = c(0, max(Mod(plot.data[,2]))))
}


acq.freq <- 100                    # data acquisition (sample) frequency (Hz)
time     <- 6                      # measuring time interval (seconds)
ts       <- seq(0,time-1/acq.freq,1/acq.freq) # vector of sampling time-points (s) 
f.0 <- 1/time

dc.component <- 1
component.freqs <- c(3,7,10)        # frequency of signal components (Hz)
component.delay <- c(0,0,0)         # delay of signal components (radians)
component.strength <- c(1.5,.5,.75) # strength of signal components

f   <- function(t,w) { 
  dc.component + 
  sum( component.strength * sin(component.freqs*w*t + component.delay)) 
}

w <- 2*pi*f.0
trajectory <- sapply(ts, function(t) f(t,w))
head(trajectory,n=30)

X.k <- fft(trajectory)

pdf("plots/freq_spectrum2.pdf")
plot.frequency.spectrum(X.k = X.k, xlimits = c(0, 20))
dev.off()

y <- trajectory
x <- seq(0, length(trajectory)-1, 1)
FT <- spec.fft(y, x, center = F)


s.rate<-44100                        # sampling frequency
t <- 2                               # seconds, for my situation, I've got 1000s of 1 - 5 minute files to go through
ind <- seq(s.rate*t)/s.rate          # time indices for each step
                                     # let's add two sin waves together to make the sound wave
f1 <- 600                            # Hz: freq of sound wave 1
y <- 100*sin(2*pi*f1*ind)            # sine wave 1
f2 <- 1500                           # Hz: freq of sound wave 2
z <- 500*sin(2*pi*f2*ind+1)          # sine wave 2
s <- y+z                             # the sound wave: my data isn't this nice, but I think this is an OK example

pdf("plots/test_sine.pdf", height = 20, width = 200)
plot(x = 1:length(s), y = s, type = "l")
dev.off()



x <- seq(0, 1, length.out = 1e3)
y <- sin(4 * 2 * pi * x) + 0.5 * sin(20 * 2 * pi * x)
X.k <- Mod(fft(trajectory))

FT <- spec.fft(y, x, center = T)
pdf("plots/spectral_freq_spectrum.pdf")
par(mfrow = c(2, 1))
plot(x, y, type = "l", main = "Signal")
plot(
FT,
ylab = "Amplitude",
xlab = "Frequency",
type = "l",
#xlim = c(-30, 30),
main = "Spectrum"
)
dev.off()

summary(FT)

pdf("plots/spectral_freq_spectrum.pdf")
plot(FT)
dev.off()





