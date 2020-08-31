dsProbitRegr = function (connections, formula, data, w = NULL, stop_tol = 1e-8, iter_max = 25L, trace = FALSE)
{
  iter = 0L
  dev_old = Inf

  if (trace) cat("\n")

  while (iter <= iter_max) {

    if (iter == 0L) {
      params_char = "xxx"
    } else {
      params_char = paste(paste0(seq_along(params), "xex", params), collapse = "xnx")
    }
    if (is.null(w)) {
      call = paste0("calcDistrParts(formula = ", formula, ", data = '", data, "', params_char = '", params_char, "')")
    } else {
      call = paste0("calcDistrParts(formula = ", formula, ", data = '", data, "', w = '", w, "', params_char = '", params_char, "')")
    }
    eval(parse(text = paste0("cq = quote(", call, ")")))

    update = DSI::datashield.aggregate(conns = connections, cq)

    lh = Reduce("+", lapply(update, function (x) x$likelihood))
    dev = -2 * log(lh)
    XtX = Reduce("+", lapply(update, function (x) x$XtX))
    Xty = Reduce("+", lapply(update, function (x) x$Xy))

    if (iter == 0L) params = rep(0, length(update[[1]]$Xy))
    params = params + solve(XtX) %*% Xty

    iter = iter + 1L

    if (trace) cat("Deviance of iter", iter, "=", round(dev, digits = 4L), "\n")

    stop_crit = abs(dev - dev_old) / (abs(dev) + 0.1)
    if (stop_crit < stop_tol) { if (trace) { cat("\n")}; break; }
    dev_old = dev
  }
  out = list(iter = iter, parameter = params, deviance = dev)
  return (out)
}

dsROCGLM = function (connections, data, formula, ds_model)
{
  formula = format(formula)
  params = ds_model$coefficients[,1]
  params_char = paste(paste0(seq_along(params), "xex", params), collapse = "xnx")

  call = paste0("rocGLMFrame('", params_char, "', '", data, "', ", formula, ")")
  eval(parse(text = paste0("cq = quote(", call, ")")))
  DSI::datashield.assign(connections, "roc_data", cq)

  roc_glm = dsProbitRegr(connections, "y ~ x", "roc_data", w = "w", trace = TRUE)
}

calculateAUC = function (roc_glm)
{
  params = roc_glm$parameter
  temp = function (x) pnorm(params[1] + params[2] * qnorm(x))
  int = integrate(f = temp, lower = 0, upper = 1)
  return (1 - int$value)
}

plotROCGLM = function (roc_glm)
{
  params = roc_glm$parameter

  x = seq(0, 1, 0.002)
  y = pnorm(params[1] + params[2]*qnorm(x))

  plot(1 - x, 1 - y, type = "l", xlab = "FPR", ylab = "TPR")
  polygon(x = c(1 - x, rev(1 - x)), y = c(1 - y, rep(0, length(y))),
    col = rgb(54,100, 139, 100, maxColorValue = 255), border = NA)
  lines(1 - x, 1 - y, lwd = 2)
}
