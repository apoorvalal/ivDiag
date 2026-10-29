\docType{data}
\name{mcclendon}
\alias{mcclendon}
\title{Data from McClendon (2014)}
\format{
A data frame with 3647 rows and 4 columns:
\describe{
  \item{intended}{1 if the person said they would attend the rally (RSVP), 0 otherwise.}
  \item{attended}{1 if the person attended the rally, 0 otherwise.}
  \item{openedesteem}{1 if the person opened an email promising public recognition for attending (the treatment), 0 otherwise.}
  \item{esteem}{1 if the person was randomly assigned to receive such an email (the instrument), 0 if assigned to the control email.}
}
}
\description{
  Data from a field experiment in McClendon (2014) AJPS.
}
\details{
McClendon (2014) worked with an LGBT advocacy organization to invite its members to a Pride rally. Members were randomly assigned to a control email or to one of two emails that promised to recognize attendees publicly, in the organization's newsletter or on its Facebook page. Random assignment to a recognition email (\code{esteem}) is an instrument for opening such an email (\code{openedesteem}); only people assigned to a recognition email could open one. The data keep the 3,647 members with complete records. In Lal et al. (2024) the outcome is \code{intended}, as in Table 2, column 1, of the original article.

The replication data are in the public domain (CC0 1.0), Harvard Dataverse, \doi{10.7910/DVN/LVOSAF}.
}
\source{
  Replication data for McClendon (2014), Harvard Dataverse, \doi{10.7910/DVN/LVOSAF} (CC0 1.0).
}
\references{
  McClendon, Gwyneth H. 2014. "Social Esteem and Participation in Contentious Politics: A Field Experiment at an LGBT Pride Rally." American Journal of Political Science 58 (2): 279--290.
}
\keyword{datasets}
