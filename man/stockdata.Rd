% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/huge-package.R
\docType{data}
\name{stockdata}
\alias{stockdata}
\title{Stock price of S&P 500 companies from 2003 to 2008}
\format{
The format is a list containing contains two matrices.
1. data - 1258x452, represents the 452 stocks' close prices for 1258 trading days.
  2. info - 452x3:
  The 1st column: the query symbol for each company.
  The 2nd column: the category for each company.
  The 3rd column: the full name of each company.
}
\source{
It was publicly available at finance.yahoo, which is now out of date
}
\usage{
data(stockdata)
}
\description{
This data set consists of stock price and company information.
}
\details{
This data set can be used to perform high-dimensional graph estimation to analyze the relationships between S&P 500 companies.
}
\examples{
data(stockdata)
image(stockdata$data)
stockdata$info
}
\keyword{datasets}
