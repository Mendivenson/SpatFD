compile_model <- function(df){

  # Función para fijar los argumentos de una función dada
  fix_args <- function(f, args_to_fix) {
    new_f <- f
    body(new_f) <- do.call(substitute, list(body(f), args_to_fix))
    formals(new_f) <- formals(f)[!names(formals(f)) %in% names(args_to_fix)]
    environment(new_f) <- globalenv()
    new_f
  }

  # Función para forzar la suma de funciones como una nueva función
  combine_models <- function(componentes) {
    llamadas <- lapply(componentes, function(fn) as.call(list(fn, quote(h))))
    expr_suma <- Reduce(function(a, b) bquote(.(a) + .(b)), llamadas)

    f <- function(h) NULL
    body(f) <- expr_suma
    environment(f) <- globalenv()
    f
  }

  all_models = list()
  for (v in unique(df$var)){
    model <- df[df$var == v,]

    componentes <- lapply(seq_len(nrow(model)), function(i) {
      x <- as.list(model[i, ])
      m <- x$model
      params <- x[setdiff(names(x), c('model', 'var'))]
      params <- params[!vapply(params, is.na, logical(1))]  # quitar NA (params que ese modelo no usa)
      fix_args(models[[m]], params)
    })

    all_models[[v]] = combine_models(componentes)
  }
  all_models
}

# par(mfrow = c(2,2))
# a <- compile_model()
# invisible(
#   lapply(names(a),
#        function(f){
#          curve(a[[f]](x), from = 0, to = 2000, main = f,
#                          xlab = 'h', ylab = 'variograma')})
# )
