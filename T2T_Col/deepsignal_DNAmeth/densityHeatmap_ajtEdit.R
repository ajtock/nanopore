densityHeatmap_ajtEdit <-
function (data, density_param = list(na.rm = TRUE), col = rev(brewer.pal(11, 
    "Spectral")), color_space = "LAB", ylab = deparse(substitute(data)), 
    column_title = paste0("Density heatmap of ", deparse(substitute(data))), 
    title = column_title, ylim = NULL, range = ylim, title_gp = gpar(fontsize = 14), 
    ylab_gp = gpar(fontsize = 12), tick_label_gp = gpar(fontsize = 10), 
    quantile_gp = gpar(fontsize = 10), show_quantiles = TRUE, 
    column_order = NULL, column_names_side = "bottom", show_column_names = TRUE, 
    column_names_max_height = unit(6, "cm"), column_names_gp = gpar(fontsize = 12), 
    column_names_rot = 90, cluster_columns = FALSE, clustering_distance_columns = "ks", 
    clustering_method_columns = "complete", mc.cores = 1, ...) 
{
    arg_list = list(...)
    if (length(arg_list)) {
        if (any(c("row_km", "row_split", "split", "km") %in% 
            names(arg_list))) {
            stop_wrap("density heatmaps do not allow row splitting.")
        }
        if (any(grepl("row", names(arg_list)))) {
            stop_wrap("density heatmaps do not allow to set rows.")
        }
        if ("anno" %in% names(arg_list)) {
            stop_wrap("`anno` is removed from the argument. Please directly construct a `HeatmapAnnotation` object and set to `top_annotation` or `bottom_annotation`.")
        }
    }
    ylab = ylab
    column_title = column_title
    density_param$na.rm = TRUE
    if (!is.matrix(data) && !is.data.frame(data) && !is.list(data)) {
        stop_wrap("only matrix and list are allowed.")
    }
    if (is.matrix(data)) {
        data2 = as.list(as.data.frame(data))
        names(data2) = colnames(data)
        data = data2
    }
    density_list = lapply(data, function(x) do.call(density, 
        c(list(x = x), density_param)))
    quantile_list = sapply(data, quantile, na.rm = TRUE)
    mean_value = sapply(data, mean, na.rm = TRUE)
    n = length(density_list)
    nm = names(density_list)
    max_x = quantile(unlist(lapply(density_list, function(x) x$x)), 
        0.98999999999999999)
    min_x = quantile(unlist(lapply(density_list, function(x) x$x)), 
        0.01)
    if (!is.null(range)) {
        max_x = range[2]
        min_x = range[1]
    }
    x = seq(min_x, max_x, length = 500)
    mat = lapply(density_list, function(r) {
        f = approxfun(r$x, r$y)
        res = f(x)
        res[is.na(res)] = 0
        rev(res)
    })
    mat = as.matrix(as.data.frame(mat))
    colnames(mat) = nm
    if (cluster_columns) {
        if (clustering_distance_columns == "ks") {
            d = ks_dist(mat, mc.cores = mc.cores)
            hc = hclust(d, clustering_method_columns)
            cluster_columns = hc
        }
    }
    if (inherits(col, "function")) {
        col = col(mat)
    }
    else {
        col = colorRamp2(seq(0, quantile(mat, 0.98999999999999999, 
            na.rm = TRUE), length = length(col)), col, space = color_space)
    }
    bb = grid.pretty(c(min_x, max_x))
    ht = Heatmap(mat, col = col, name = NULL, column_title = title, 
        column_title_gp = title_gp, cluster_rows = FALSE, cluster_columns = cluster_columns, 
        clustering_distance_columns = clustering_distance_columns, 
        clustering_method_columns = clustering_method_columns, 
        column_dend_reorder = mean_value, column_names_side = column_names_side, 
        show_column_names = show_column_names, column_names_max_height = column_names_max_height, 
        column_names_gp = column_names_gp, column_names_rot = column_names_rot, 
        column_order = column_order, left_annotation = rowAnnotation(axis = anno_empty(border = FALSE, 
            width = grobHeight(textGrob(ylab, gp = ylab_gp)) * 
                2 + max_text_width(bb, gp = tick_label_gp) + 
                unit(4, "mm")), show_annotation_name = FALSE), 
        right_annotation = {
            if (show_quantiles) {
                rowAnnotation(quantile = anno_empty(border = FALSE, 
                  width = grobWidth(textGrob("100%", gp = quantile_gp)) + 
                    unit(6, "mm")), show_annotation_name = FALSE)
            }
            else NULL
        }, ...)
    random_str = paste(sample(c(letters, LETTERS, 0:9), 8), collapse = "")
    ht@name = paste0(ht@name, "_", random_str)
    names(ht@left_annotation) = paste0(names(ht@left_annotation), 
        "_", random_str)
    if (show_quantiles) {
        names(ht@right_annotation) = paste0(names(ht@right_annotation), 
            "_", random_str)
    }
    post_fun = function(ht) {
        column_order = column_order(ht)
        if (!is.list(column_order)) {
            column_order = list(column_order)
        }
        n_slice = length(column_order)
        decorate_annotation(paste0("axis_", random_str), {
            grid.text(ylab, x = grobHeight(textGrob(ylab, gp = ylab_gp)), 
                rot = 90)
        }, slice = 1)
        if (!is.null(ht@right_annotation)) {
            for (i_slice in 1:n_slice) {
                decorate_heatmap_body(paste0("density_", random_str), 
                  {
                    n = length(column_order[[i_slice]])
                    pushViewport(viewport(xscale = c(0.5, n + 
                      0.5), yscale = c(min_x, max_x), clip = TRUE))
                    for (i in seq_len(5)) {
                      grid.lines(1:n, quantile_list[i, column_order[[i_slice]]], 
                        default.units = "native", gp = gpar(lty = 2))
                    }
                    grid.lines(1:n, mean_value[column_order[[i_slice]]], 
                      default.units = "native", gp = gpar(lty = 2, 
                        col = "darkred"))
                    upViewport()
                  }, column_slice = i_slice)
            }
        }
        decorate_heatmap_body(paste0("density_", random_str), 
            {
                pushViewport(viewport(yscale = c(min_x, max_x), 
                  clip = FALSE))
                grid.rect(gp = gpar(fill = NA))
                grid.yaxis(gp = tick_label_gp)
                upViewport()
            }, column_slice = 1)
        if (!is.null(ht@right_annotation)) {
            decorate_heatmap_body(paste0("density_", random_str), 
                {
                  n = length(column_order[[n_slice]])
                  lq = !apply(quantile_list, 1, function(x) all(x > 
                    max_x) || all(x < min_x))
                  lq = c(lq, !(all(mean_value > max_x) || all(mean_value < 
                    min_x)))
                  if (sum(lq) == 0) {
                    return(NULL)
                  }
                  labels = c(rownames(quantile_list), "mean")
                  y = c(quantile_list[, column_order[[n_slice]][n]], 
                    mean_value[column_order[[n_slice]][n]])
                  labels = labels[lq]
                  y = y[lq]
                  od = order(y)
                  y = y[od]
                  labels = labels[od]
                  pushViewport(viewport(xscale = c(0.5, n + 0.5), 
                    yscale = c(min_x, max_x), clip = FALSE))
                  text_height = convertHeight(grobHeight(textGrob(labels[1])) * 
                    2, "native", valueOnly = TRUE)
                  h1 = y - text_height * 0.5
                  h2 = y + text_height * 0.5
                  pos = rev(smartAlign(h1, h2, c(min_x, max_x)))
                  h = (pos[, 1] + pos[, 2])/2
                  link_width = unit(6, "mm")
                  n2 = length(labels)
                  grid.text(labels, unit(1, "npc") + rep(link_width, 
                    n2), h, default.units = "native", just = "left", 
                    gp = quantile_gp)
                  link_width = link_width - unit(1, "mm")
                  ly = y <= max_x & y >= min_x
                  if (sum(ly)) {
                    grid.segments(unit(rep(1, n2), "npc")[ly], 
                      y[ly], unit(1, "npc") + rep(link_width * 
                        (1/3), n2)[ly], y[ly], default.units = "native")
                    grid.segments(unit(1, "npc") + rep(link_width * 
                      (1/3), n2)[ly], y[ly], unit(1, "npc") + 
                      rep(link_width * (2/3), n2)[ly], h[ly], 
                      default.units = "native")
                    grid.segments(unit(1, "npc") + rep(link_width * 
                      (2/3), n2)[ly], h[ly], unit(1, "npc") + 
                      rep(link_width, n2)[ly], h[ly], default.units = "native")
                  }
                  upViewport()
                }, column_slice = n_slice)
        }
    }
    ht@heatmap_param$post_fun = post_fun
    ht_list = ht %v% NULL
    return(ht_list)
}
