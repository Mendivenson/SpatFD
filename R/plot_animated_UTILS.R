# ===================== FUNCIONES INTERNAS DE plot_animated() ==================
# Solo las usa plot_animated(). La grilla de celdas está en utils.R, porque
# también la usa el mapa estático de plot.SpatFDpred().

# TIEMPOS DEL SLIDER, EN DÍAS: Un número es la cantidad de tiempos
# equidistantes entre el primero y el último; un vector, los tiempos en el
# formato de los datos
slider_days <- function(times, range, template){
  if (is.numeric(times) && length(times) == 1 && !inherits(times, c('Date', 'POSIXt'))) {
    stopifnot("'times' must be a whole number of at least 2, or a vector of times" =
                times >= 2 && times == round(times))
    return(seq(range[1], range[2], length.out = times))
  }

  days <- as_days(times)
  stopifnot(
    "'times' must be in the same format as the times of the data" = is.numeric(days) && !anyNA(days),
    "'times' must be within the times of the curves" = all(days >= range[1] & days <= range[2])
  )
  sort(days)
}

# ETIQUETAS DEL SLIDER: Los tiempos en el formato de los datos
slider_labels <- function(days, template){
  if (is.numeric(template)) return(format(round(days, 2), trim = TRUE, drop0trailing = TRUE))
  labels <- format(from_days(days, template))
  make.unique(labels)                                    # Dos tiempos en el mismo día, por ejemplo
}

# CONTORNO DEL MAPA: Coordenadas de los bordes de los polígonos, con NA entre
# un borde y el siguiente para que plotly no los una
map_border <- function(map){
  border <- sf::st_coordinates(sf::st_cast(sf::st_boundary(sf::st_geometry(map)), 'MULTILINESTRING'))
  pieces <- split(as.data.frame(border[, 1:2]), interaction(border[, 'L1'], border[, 'L2'], drop = TRUE))
  do.call(rbind, lapply(pieces, function(piece) rbind(piece, NA)))
}

# MAPAS ANIMADOS: Dos paneles lado a lado. El primer trazo es el mapa de la
# predicción, y es el único que cambia: cada cuadro de la animación solo
# reemplaza su matriz de valores. El contorno, los sitios observados y el mapa
# de la varianza quedan fijos.
animated_maps <- function(cells, prediction, variance, labels, border, stations, limits, main){
  hover <- 'x: %{x}<br>y: %{y}<br>value: %{z:.4g}<extra></extra>'

  fig <- plotly::plot_ly()

  # PANEL IZQUIERDO: Predicción (trazo 0, el que se anima)
  fig <- plotly::add_heatmap(fig, x = cells$x, y = cells$y, z = prediction[[1]], zmin = limits[1], zmax = limits[2],
                             colorscale = 'Viridis', reversescale = TRUE, hovertemplate = hover,
                             colorbar = list(title = 'Prediction', x = 0.44, len = 0.75), xaxis = 'x', yaxis = 'y')

  # PANEL DERECHO: Varianza
  fig <- plotly::add_heatmap(fig, x = cells$x, y = cells$y, z = variance,
                             colorscale = 'Viridis', reversescale = TRUE, hovertemplate = hover,
                             colorbar = list(title = 'Variance', x = 1, len = 0.75), xaxis = 'x2', yaxis = 'y2')

  # Contorno y sitios observados, en los dos paneles
  for (axes in list(c('x', 'y'), c('x2', 'y2'))) {
    fig <- plotly::add_paths(fig, x = border$X, y = border$Y, line = list(color = 'black', width = 0.6),
                             hoverinfo = 'skip', showlegend = FALSE, xaxis = axes[1], yaxis = axes[2])
    fig <- plotly::add_markers(fig, x = stations[, 1], y = stations[, 2], text = rownames(stations),
                               hoverinfo = 'text', showlegend = FALSE, xaxis = axes[1], yaxis = axes[2],
                               marker = list(color = 'white', size = 7, line = list(color = 'black', width = 1)))
  }

  # SLIDER Y BOTONES: Cada paso muestra el cuadro de su tiempo
  step <- function(label) {
    list(label = label, method = 'animate',
         args = list(list(label), list(mode = 'immediate', frame = list(duration = 0, redraw = TRUE),
                                       transition = list(duration = 0))))
  }
  play  <- list(label = 'Play', method = 'animate',
                args = list(NULL, list(fromcurrent = TRUE, frame = list(duration = 200, redraw = TRUE),
                                       transition = list(duration = 0))))
  pause <- list(label = 'Pause', method = 'animate',
                args = list(list(NULL), list(mode = 'immediate', frame = list(duration = 0, redraw = FALSE))))

  axis <- list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, ticks = '', title = '')
  fig <- plotly::layout(fig,
    title       = list(text = main),
    xaxis       = c(axis, list(domain = c(0, 0.42))),
    yaxis       = c(axis, list(scaleanchor = 'x')),
    xaxis2      = c(axis, list(domain = c(0.55, 0.97))),
    yaxis2      = c(axis, list(scaleanchor = 'x2', anchor = 'x2')),
    annotations = list(list(text = 'Prediction', x = 0.21, y = -0.02, xref = 'paper', yref = 'paper',
                            showarrow = FALSE, font = list(size = 15)),
                       list(text = 'Variance', x = 0.76, y = -0.02, xref = 'paper', yref = 'paper',
                            showarrow = FALSE, font = list(size = 15))),
    sliders     = list(list(active = 0, currentvalue = list(prefix = 'Time: '), pad = list(t = 30),
                            steps = lapply(labels, step))),
    updatemenus = list(list(type = 'buttons', showactive = FALSE, x = 0, y = -0.08, xanchor = 'right',
                            yanchor = 'top', direction = 'left', buttons = list(play, pause))))

  # CUADROS: Uno por tiempo, solo con la matriz de la predicción
  fig <- plotly::plotly_build(fig)
  fig$x$frames <- Map(function(z, label) list(name = label, data = list(list(z = z)), traces = list(0)),
                      prediction, labels)
  fig
}
