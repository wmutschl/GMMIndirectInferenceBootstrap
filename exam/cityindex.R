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
#Die Ökonomische Intuition hinter der Hypothese ist das sogenannte Zipfsche Gesetz (vergleiche Gabaix und Ioannides (2004, S.2344)): 
# Dieses besagt, dass die zweitgrößte Stadt halb so groß ist wie die größte, die drittgrößte Stadt ein Drittel der Größe der größten besitzt usw..
# Ordnet man n Städte also hinsichtlich ihrer Größe X_1,...,X_n vom größten (Rang 1) zum kleinsten (Rang n),
# so ist der Rang i einer Stadt der Größe X_i proportional zum Anteil der Städte die größer als i sind. 
# Es gilt also: X_i proportional zu k/i für eine Konstante k. Das Bilden von Logarithmen führt zu einem linearen Modell der Form: 
# ln(X_i) = ln(k) - ln(i) bzw. ln(i) = ln(k) - ln(X_i). Ein Wert von alpha=1 im angegebenen Regressionsmodell würde also bedeuten,
# dass das Zipfsche Gesetz Gültigkeit hat.
# nte-größte Stadt in einem Land hat 1/n mal so viele Einwohner wie die größte Stadt. Ist alpha=1 spricht dies für ein ausgewogenes
# Städtesystem, für alpha>1 sind Städte tendenziell größer im Vergleich zur größten Stadt und für alpha <1 tendenziell kleiner.

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