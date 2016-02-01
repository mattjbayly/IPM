##################################################
#### PERIODIC MATRIX MODEL TUTORIAL
##################################################
# SOURCE
#http://www.r-bloggers.com/periodic-matrix-model-for-annual-plant-demography/

#### TRANSITION VALUES
seed  <- 0.9^4  # Seed surviving rate; annual
germ  <- 0.3    # Germination rate; spring
plant <- 0.05   # Plant surviving rate; from germination to mature
yield <- 100    # Seed production number; per matured plant
mature <- 0.5  # pollination rate, per mature plant

##################################################
### MATRIX FUNCTION VALUES

seedq             <- function() seed^(1/4)

p11               <- function() plant / mature
q11               <- function() mature
r11               <- function() yield
s11               <- function() germ * seedq()
s21               <- function() (1 - germ) * seedq()
p22 <- q22 <- r12 <- function() seedq()

##################################################
