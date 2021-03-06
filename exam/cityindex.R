# ===========================================================================
#
#      Test for the Zipf index of city size distributions
#
# ===========================================================================
rm (list = ls(all=TRUE))
graphics.off()
library(MASS)

############
# Part (a) #
############
#Since the sample is ordered, the observations are no longer independent and the optimality properties of OLS
#vanish. In particular, the ordinary t-test does not work correctly anymore. 


############
# Part (b) #
############
R <- 1000
n <- 20
y <- 1:n
Z <- rep(NA,R)

for (r in 1:R) {
  x <- sort(exp(rexp(n)),decreasing=TRUE)
  obj <- lm(log(y) ~ log(x))
  Z[r] <- coefficients(obj)[2]
}
truehist(-Z)

############
# Part (c) #
############
R <- 1000
n <- 20
y <- 1:n
Tstat <- rep(NA,R)
# Hypothetical value, Nullhypothesis
alpha0 <- 1

for (r in 1:R) {
  x <- sort(exp(rexp(n)),decreasing=TRUE)
  obj <- lm(log(y) ~ log(x))
  Tstat[r] <- (-coefficients(obj)[2]-alpha0)/ sqrt(vcov(obj)[2,2])
}
truehist(Tstat,xlim=c(-40,40))
curve(dt(x,df=n-2),add=T)

# Economic interpretation of the Hypothesis:
# The second largest city is one half the size of the largest city, the third largest city is one-third the size of the largest city
# The rank i of a city is proportional to the number of cities larger than i
#Die �konomische Intuition hinter der Hypothese ist das sogenannte Zipfsche Gesetz (vergleiche Gabaix und Ioannides (2004, S.2344)): 
# Dieses besagt, dass die zweitgr��te Stadt halb so gro� ist wie die gr��te, die drittgr��te Stadt ein Drittel der Gr��e der gr��ten besitzt usw..
# Ordnet man n St�dte also hinsichtlich ihrer Gr��e X_1,...,X_n vom gr��ten (Rang 1) zum kleinsten (Rang n),
# so ist der Rang i einer Stadt der Gr��e X_i proportional zum Anteil der St�dte die gr��er als i sind. 
# Es gilt also: X_i proportional zu k/i f�r eine Konstante k. Das Bilden von Logarithmen f�hrt zu einem linearen Modell der Form: 
# ln(X_i) = ln(k) - ln(i) bzw. ln(i) = ln(k) - ln(X_i). Ein Wert von alpha=1 im angegebenen Regressionsmodell w�rde also bedeuten,
# dass das Zipfsche Gesetz G�ltigkeit hat.
# nte-gr��te Stadt in einem Land hat 1/n mal so viele Einwohner wie die gr��te Stadt. Ist alpha=1 spricht dies f�r ein ausgewogenes
# St�dtesystem, f�r alpha>1 sind St�dte tendenziell gr��er im Vergleich zur gr��ten Stadt und f�r alpha <1 tendenziell kleiner.

############
# Part (d) #
############
# Note that in part (b) we already did a parametric LM-Type bootstrap
Tsharp = Tstat
# Sort the Tsharp values
Tsharp <- sort(Tsharp)
truehist(Tsharp)
critlow <- Tsharp[0.025*R]
crithigh <- Tsharp[0.975*R]
print(c(critlow,crithigh))