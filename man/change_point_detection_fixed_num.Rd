\name{change_point_detection_fixed_num}
\alias{change_point_detection_fixed_num}
\title{Fixed-Segment Change-Point Detection with Gurobi}
\usage{
change_point_detection_fixed_num(Length, X, Delta, output_flag = 1L)
}
\arguments{
  \item{Length}{Integer. Number of segments. Must be at least 2.}
  \item{X}{Numeric vector containing the sequential observations.}
  \item{Delta}{Integer. Minimum number of observations assigned to each segment.}
  \item{output_flag}{Integer. Use 1 to display Gurobi output and 0 to suppress it.}
}
\value{
A numeric vector of cumulative segment ending positions.
}
\description{
Solves a mixed-integer optimization formulation for segmented empirical
cumulative distribution functions. The number of segments is fixed by the
user. The returned vector gives the cumulative ending positions of the
detected segments.
}
\examples{
\dontrun{
x <- c(1, 1, 2, 2, 10, 10)
change_point_detection_fixed_num(Length = 2, X = x, Delta = 2)
}
}
