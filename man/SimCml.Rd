% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimCml.R
\docType{data}
\name{SimCml}
\alias{SimCml}
\title{Simulated data for CML patients}
\format{A data frame of 1777 rows (patients) on the following 14 variables:
\describe{
\item{cp1}{Time in months to chronic phase}
\item{cp1.s}{Indicator of the occurrence of chronic phase; 1=occur, 0=censored}
\item{ap}{Time in months to accelerated phase of CML}
\item{ap.s}{Indicator of the occurrence of accelerated phase; 1=occur, 0=censored}
\item{cp2}{Time in months to chronic phase after progression to advanced stage}
\item{cp2.s}{Indicator whether subject come back to chronic phase after progression to an advanced stage; 1=occur, 0=censored}
\item{bc}{Time in months to blast crisis phase of CML}
\item{bc.s}{Indicator of the occurrence of blast crisis phase; 1=occur, 0=censored}
\item{sct}{Time in months to receive transplant}
\item{sct.s}{Indicator whether patient receive transplant; 1=transplant, 0=no transplant}
\item{death}{Time in months to death}
\item{death.s}{Death indicator; 1=death, 0=censored}
\item{age}{0=age less than 50; 1=age larger than 50}
\item{sex}{1=female; 2=male}
}}
\description{
A dataset containing information of CML patients who have received transplant.
Both states before and after transplant are included in this data set.
This data set is used for illustration, so the magnitude of event times are out of clinical consideration.
}
\examples{
data(SimCml)
}
