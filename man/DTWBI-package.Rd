\name{DTWBI-package}
\alias{DTWBI-package}
\alias{DTWBI}
\docType{package}
\title{
\packageTitle{DTWBI}
}
\description{
\packageDescription{DTWBI}
}
\details{
%\packageDESCRIPTION{DTWBI}
\packageIndices{DTWBI}
}
\author{
\packageAuthor{DTWBI}

Maintainer: \packageMaintainer{DTWBI}
}
\references{
Thi-Thu-Hong Phan, Emilie Poisson-Caillault, Alain Lefebvre, Andre Bigand. Dynamic time warping-
based imputation for univariate time series data.  Pattern Recognition Letters, Elsevier, 2017,
<DOI:10.1016/j.patrec.2017.08.019>. <hal-01609256>
}
\keyword{package}
\keyword{DTW}
\keyword{times series}
\keyword{imputation}
\examples{
# Load package dataset
data(dataDTWBI)

# Create a query and a reference signal
query <- dataDTWBI$query
ref <- dataDTWBI$query

# Create a gap within query (10\% of signal size)
query <- gapCreation(query, rate = 0.1)
t <- query$begin_gap
T <- query$gap_size

# Fill gap using DTWBI algorithm
results_DTWBI <- DTWBI_univariate(query$output_vector, t, T)

# Plot
plot(ref, type = "l")
lines(results_DTWBI$output_vector, col = "red", lty = "dashed")

# Compute the similarity of imputed vector and reference
compute.sim(ref, results_DTWBI$output_vector)
}
