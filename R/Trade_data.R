#' Determinants of International Trade in the European Union
#'
#' @description
#' The dataset contains bilateral trade and its economic determinants
#' for 26 European Union countries over the period 1995–2015.
#' Each observation represents a country pair, resulting in 325
#' unique trading pairs.
#'
#' The countries included are: Austria, Belgium, Bulgaria, Cyprus,
#' Czechia, Denmark, Estonia, Finland, France, Germany, Greece,
#' Hungary, Ireland, Italy, Latvia, Lithuania, Luxembourg,
#' the Netherlands, Poland, Portugal, Romania, Slovakia, Slovenia,
#' Spain, Sweden, and the United Kingdom.
#'
#' The dataset was used in:
#' Beck (2020), *What drives international trade? Robust analysis
#' for the European Union*, Journal of International Studies,
#' 13(3), 68–84.
#'
#' @format A data frame with 325 rows and 18 variables:
#' \describe{
#'   \item{LNTRADE}{Natural logarithm of bilateral trade between two countries.}
#'   \item{B}{Border dummy (1 if countries share a common border; 0 otherwise).}
#'   \item{LNDGEO}{Natural logarithm of the shortest distance between capital cities.}
#'   \item{L}{Common language dummy (1 if countries share at least one official language; 0 otherwise).}
#'   \item{LNRGDPPROD}{Natural logarithm of the product of real GDPs of the two countries.}
#'   \item{RGDPpcDIFF}{Absolute difference in real GDP per capita.}
#'   \item{GOV}{Absolute difference in government expenditure shares of GDP.}
#'   \item{HUMAN}{Absolute difference in human capital indicators.}
#'   \item{CPW}{Absolute difference in capital per worker.}
#'   \item{INFVAR}{Absolute difference in the standard deviation of inflation rates.}
#'   \item{ARABLE}{Absolute difference in arable land.}
#'   \item{ARABLEpw}{Absolute difference in arable land per worker.}
#'   \item{LAND}{Absolute difference in total land area.}
#'   \item{LANDpc}{Absolute difference in land area per capita.}
#'   \item{EPCpc}{Absolute difference in electricity consumption per capita.}
#'   \item{FDI}{Absolute difference in foreign direct investment flows.}
#'   \item{KSI}{Krugman specialization index (value added, 35 sectors).}
#'   \item{BCIDIFF}{Absolute difference in the Bayesian Corruption Index.}
#' }
#'
#' @source
#' Harvard Dataverse.
#' \doi{10.7910/DVN/JMMOEA}
#'
#' @usage data(Trade_data)
#'
"Trade_data"
