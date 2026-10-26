\name{s3cr_census}
\alias{s3cr_census}
\alias{census}
\alias{species}
\title{s3cR census and species trait data}
\description{
Four-site, three-date community counts and species-level trait values from
the example data distributed with the \pkg{s3cR} package. The data are a
lightweight subset suitable for demonstrating functional change methods.
}
\usage{
data(census)
data(species)
}
\format{
\code{census} is a data frame with 48 rows and columns \code{site},
\code{date}, \code{species}, and \code{n}. \code{species} is a data frame
with four rows and columns \code{species}, \code{size}, and \code{size_v}.
}
\source{
\url{https://github.com/pgauzere/s3cR}; original data objects
\code{census} and \code{species}.
}
\seealso{\code{\link{make_fundecomp_input}}}
