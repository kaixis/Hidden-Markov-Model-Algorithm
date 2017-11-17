#' kaixisHW4.
#'
#' @name kaixisHW4
#' @docType package
#' @author Kaixi Song
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib kaixisHW4
#' 
#' 
#' 
#' This function implement a Viterbi algorithm, given obeserved outcomes, transition and emission
#' matrices, and initial probabilities. 
#'
#'@param obs an m*n matrix representing probability from a given state to an observed outcome. 
#'@param theta a parameter which determines transition probability between states. 
#'@param ts a vector contains time point
#'
#'@return a vector represent the path that can lead to maximum likelihood
#'


