#'
#' @title Predict with given DataSHIELD GLM
#' @description This function calculates the prediction scores for a given GLM
#' @param params_char Character containing all the parameter
#' @param data Character containing the name of the dataset
#' @param formula The formula used to fit the model
#' @return Prediction scores
#' @author Stefan B., Daniel S.
predictStringGLM = function (params_char, data, formula)
{
  if (missing(data))
    stop("Need name of dataset used to train 'ds_mod'")

  if ( (length(params_char) != 1) && (! is.character(param_char)) )
    stop("'param_char' has to be a character of length 1")

  fm = as.formula(formula)

  params = unlist(lapply(strsplit(params_char, "<n>")[[1]], FUN = function (p) {
    sp = strsplit(p, "<=>")
    params = vapply(sp, FUN.VALUE = numeric(1L), function (s) as.numeric(s[2]))
    names(params) = vapply(sp, FUN.VALUE = character(1L), function (s) s[1])

    return (params)
  }))

  X = eval(parse(text = paste0("model.matrix(fm, data = ", data, ")")))
  return (X %*% params)
}

#'
#' @title Compute Placement Values on Server
#' @description This function calculates the placement values required to calculate the ROC-GLM.
#' @param score Prediction scores
#' @param truth True zero-one values as integer
#' @param tset Set of thresholds
#' @return Placement values
#' @author Stefan B., Daniel S.
computePlacementValues = function (score, truth, tset = NULL)
{
  if (is.null(tset)) {
    tset = score[truth == 1]
  }
  tset = sort(tset)

  out = vapply(X = tset, FUN.VALUE = numeric(1L), FUN = function (t) {
    # FPR for each threshold
    pred = ifelse(score >= t, 1, 0)
    return (mean((pred != truth)[truth == 0]))
  })
  return (out)
}

#'
#' @title Calculate U Matrix for ROC-GLM
#' @description This function calculates U matrix which is used as target variable for the ROC-GLM.
#' @param tset Set of thresholds
#' @param pv Placement values
#' @return Matrix of zeros and ones that are used as target variable for Probit regression
#' @author Stefan B., Daniel S.
calcU = function (tset, pv)
{
  tset_sorted = sort(tset)
  out = vapply(X = tset, FUN.VALUE = integer(length(pv)), FUN = function (th) {
    ifelse(pv <= th, 1L, 0L)
  })
  return (out)
}

#'
#' @title Get Data for ROC-GLM
#' @description This function calculates the data used for the ROC-GLM
#' @param U Response matrix for ROC-GLM
#' @param tset Set of thresholds
#' @return Data used for the ROC-GLM
#' @author Stefan B., Daniel S.
rocGLMData = function (U, tset)
{
  roc_glm_data = data.frame(
    y = rep(c(0, 1), times = length(tset)),
    x = rep(qnorm(tset), each = 2L),
    w = as.vector(apply(U, 2, function (x) c(sum(x == 0), sum(x == 1))))
  )
  return (roc_glm_data)
  #mod = glm(y ~ x, data = roc_glm_data, family = binomial(link = "probit"), weights = roc_glm_data$w)
  #return (mod)
}

#'
#' @title Calculate data for ROC-GLM
#' @description This function stores the data on the server used for the ROC-GLM#' @param params_char Character containing all the parameter
#' @param params_char Character containing all the parameter
#' @param data Character containing the name of the dataset
#' @param formula The formula used to fit the model
#' @return Data.frame used to calculate ROC-GLM
#' @author Stefan B., Daniel S.
#' @export
rocGLMFrame = function (params_char, data, formula)
{
  scores = predictStringGLM(params_char, data, formula)
  target = strsplit(ds_mod$formula, " ~ ")[[1]][1]
  truth = eval(parse(text = paste0(data, "[[", target, "]]")))

  tset = scores[truth == 1]
  pv = computePlacementValues(scores, truth, tset)
  U  = calcU(tset, pv)
  roc_clm_data = rocGLMData(U, tset)

  return (roc_glm_data)
}

