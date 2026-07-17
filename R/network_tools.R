# =============================================================================
# network_tools.R
#
# Interactive relatedness-network visualization, built on visNetwork.
# =============================================================================

#' Plot an interactive relatedness/kinship network
#'
#' Builds an interactive network (via \pkg{visNetwork}) of individuals
#' arranged in two semicircles by sex (females on one side, males on the
#' other), connected by edges for pairs whose kinship falls at or below a
#' chosen threshold. Intended for identifying and exploring the
#' \emph{least}-related pairs/candidates (e.g. for selecting low-relatedness
#' breeding groups) — edges are only drawn for pairs below the threshold,
#' not above it.
#'
#' Clicking a node highlights it, its directly connected nodes, and the
#' edges between them, and fades everything else out. Clicking empty space
#' (deselecting) restores the original edge colours.
#'
#' @param krdf A data fram of pairwise kinship/relatedness values, with
#'   columns \code{from}, \code{to}, \code{kinship|relatedness}, for example
#'   the output of \code{swingeR:::QGrel()})or \code{SNPRelate::snpgdsIBDKING()}
#'   reshaped to data frame format.
#' @param meta A data frame with at least columns \code{ID} and \code{Sex}
#'   (\code{"F"}/\code{"M"}), used to lay out and colour the nodes.
#' @param td Threshold used to decide which edges to draw. Interpreted
#'   differently depending on \code{method} (see below). Default \code{0.1}.
#' @param method Character string specifying how \code{"td"} is interpreted.
#'   \code{"fixed"} uses \code{"td"} directly as threshold, drawing
#'   edges where the relationship value is less than or equal to \code{"td"}.
#'   \code{"relative"} (default); interprets \code{"td"} as lowest quantile of the observed
#'   relationship values; for example, \code{"td = 0.1"} draws the least-related
#'   10\% of pairs.
#' @param kr Character string specifying the name of the column containing
#'   the pairwise relationship values used to construct the network
#'   (e.g. \code{"kinship"} or \code{"relatedness"}). Default is \code{"kinship"}.
#`
#' @return A \code{visNetwork} htmlwidget object. Print it (or let it
#'   auto-print at the top level) to display the interactive network; save
#'   it to a standalone HTML file with
#'   \code{visNetwork::visSave(net, "relatedness_network.html")}.
#'
#' @details
#' Requires the \pkg{dplyr} and \pkg{visNetwork} packages (declared as
#' package dependencies). Node positions are fixed (not physics-based) in a
#' semicircular layout by sex, so the graph stays readable for larger
#' individual counts rather than re-arranging on every interaction.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' net <- networkY(kinship_df, meta, td = 0.1, method = "fixed")
#' net
#' visNetwork::visSave(net, "relatedness_network.html")
#' }
networkY <- function(krdf, meta, td = 0.1, method = "relative", kr = "kinship") {

  females <- meta$ID[meta$Sex == "F"]
  males   <- meta$ID[meta$Sex == "M"]
  n_f <- length(females)
  n_m <- length(males)
  cat("Females:", n_f, "Males:", n_m, "\n")

  radius <- 400
  female_angles <- seq(pi / 2, 3 * pi / 2, length.out = n_f + 2)[c(-1, -(n_f + 2))]
  male_angles   <- seq(3 * pi / 2, 5 * pi / 2, length.out = n_m + 2)[c(-1, -(n_m + 2))]

  nodes_f <- data.frame(
    id    = females,
    label = females,
    sex   = "F",
    color = "#E87D7D",
    origColor = "#E87D7D",
    shape = "dot",
    x     = as.numeric(-radius * cos(female_angles - pi / 2)),
    y     = as.numeric(radius  * sin(female_angles - pi / 2))
  )
  nodes_m <- data.frame(
    id    = males,
    label = males,
    sex   = "M",
    color = "#7DB3E8",
    origColor = "#7DB3E8",
    shape = "square",
    x     = as.numeric(-radius * cos(male_angles - pi / 2)),
    y     = as.numeric(radius  * sin(male_angles - pi / 2))
  )
  nodes <- rbind(nodes_f, nodes_m)

  edges <- krdf %>%
    dplyr::mutate(
      value = .data[[kr]],
      title = paste0("\u03c6 = ", round(.data[[kr]], 4)),
      color = dplyr::case_when(
        .data[[kr]] >= 0.25   ~ "#E84040",
        .data[[kr]] >= 0.125  ~ "#E8A040",
        .data[[kr]] >= 0.0625 ~ "#E8D840",
        .data[[kr]] >= 0.02   ~ "#33CC00",
        .data[[kr]] >= 0.01   ~ "#99CC00",
        .data[[kr]] < -0.02   ~ "#009900",
        TRUE                  ~ "#AAAAAA"
      )
    )
  cat("Edges shown:", nrow(edges), "\n")

  edges <- edges %>%
    dplyr::filter(as.character(from) < as.character(to))

  if (method == "fixed") {
    threshold <- td
  } else {
    threshold <- stats::quantile(edges[[kr]], td)
  }

  edges_filtered <- edges %>%
    dplyr::filter(.data[[kr]] <= threshold) %>%
    dplyr::mutate(
      value  = 1,
      title  = paste0("\u03c6 = ", round(.data[[kr]], 4)),
      color  = "#848484",
      hidden = FALSE
    )
  cat("Edges shown:", nrow(edges_filtered), "\n")

  centerX <- mean(nodes$x)
  centerY <- mean(nodes$y)

  edges_with_curve <- edges_filtered %>%
    dplyr::left_join(nodes %>% dplyr::select(id, x, y), by = c("from" = "id")) %>%
    dplyr::rename(fx = x, fy = y) %>%
    dplyr::left_join(nodes %>% dplyr::select(id, x, y), by = c("to" = "id")) %>%
    dplyr::rename(tx = x, ty = y) %>%
    dplyr::mutate(
      midX      = (fx + tx) / 2,
      midY      = (fy + ty) / 2,
      toCenterX = centerX - midX,
      toCenterY = centerY - midY,
      # Edge vector based on actual from/to direction
      ex        = tx - fx,
      ey        = ty - fy,
      cross     = ex * toCenterY - ey * toCenterX,
      # so each edge gets the correct curve for ITS direction
      smooth.type      = ifelse(cross > 0, "curvedCCW", "curvedCW"),
      smooth.enabled   = TRUE,
      smooth.roundness = 0.5
    ) %>%
    dplyr::select(-fx, -fy, -tx, -ty, -midX, -midY,
                  -ex, -ey, -toCenterX, -toCenterY, -cross)

  net <- visNetwork::visNetwork(
    nodes = nodes,
    edges = edges_with_curve,
    width  = "100%",
    height = "800px",
    main   = list(
      text  = paste0(kr, " network (", kr, " < ", round(threshold,4), ")"),
      style = "font-size:16px; font-weight:bold; text-align:center;"
    )
  ) %>%
    visNetwork::visPhysics(enabled = FALSE, stabilization = FALSE) %>%
    visNetwork::visNodes(
      fixed       = list(x = TRUE, y = TRUE),
      size        = 10,
      font        = list(size = 11, bold = TRUE),
      borderWidth = 2
    ) %>%
    visNetwork::visEdges(
      color   = list(color = "#848484", opacity = 0.3),
      scaling = list(min = 1, max = 6)
    ) %>%
    visNetwork::visOptions(
      highlightNearest = FALSE,
      nodesIdSelection = list(
        enabled = TRUE,
        main    = "Select individual",
        style   = "width:200px; font-size:14px;"
      )
    ) %>%
    visNetwork::visInteraction(hover = TRUE, tooltipDelay = 100) %>%
    visNetwork::visLegend(
      addNodes = list(
        list(label = "Female", shape = "dot",    color = "#E87D7D", size = 15),
        list(label = "Male",   shape = "square", color = "#7DB3E8", size = 15)
      ),
      useGroups = FALSE,
      position  = "left",
      width     = 0.15
    ) %>%
    visNetwork::visEvents(
      selectNode = "function(params) {
        var selectedId     = params.nodes[0];
        var allNodes       = this.body.data.nodes.get();
        var connectedNodes = this.getConnectedNodes(selectedId);
        connectedNodes.push(selectedId);
        var nodeUpdates = allNodes.map(function(n) {
          return {
            id:    n.id,
            color: connectedNodes.includes(n.id) ? n.origColor : '#DDDDDD'
          };
        });
        this.body.data.nodes.update(nodeUpdates);
        var connectedEdges = this.getConnectedEdges(selectedId);
        var allEdges       = this.body.data.edges.get();
        var updates        = allEdges.map(function(e) {
          return {
            id:    e.id,
            color: connectedEdges.includes(e.id)
              ? { color: '#E84040', opacity: 1.0 }
              : { color: '#848484', opacity: 0.03 }
          };
        });
        this.body.data.edges.update(updates);
      }",
      deselectNode = "function() {
        var allEdges = this.body.data.edges.get();
        var updates  = allEdges.map(function(e) {
          return { id: e.id, color: '#848484' };
        });
        this.body.data.edges.update(updates);
      }"
    )

  return(net)
}
