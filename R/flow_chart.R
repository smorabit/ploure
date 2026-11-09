#---------------------------------------------------------#
# QC flow chart
#
# Generalized from PLEURA's v2/overview/qc_flow_chart.Rmd: a left-to-right chart in plain ggplot (boxes are tiles,
# edges are segments). A shared trunk of steps fans out into one row per branch (e.g. lineage); each step box shows
# the cells kept and its criterion, with the excluded cells (and why) in a box below it.
#---------------------------------------------------------#

#' QC flow chart
#'
#' Draws a cell-filtering pipeline from a table of steps: a trunk of shared steps, which can fan out into one row per
#' branch (e.g. one per lineage). Each step box shows the cells kept and the step's criterion; steps that remove cells
#' get an exclusion box below them, with the number removed and, from `removals`, the count per category. Branch steps
#' are aligned in columns by step name, and a branch that skips a step gets a dashed placeholder box.
#'
#' Text uses the Unicode minus sign, so write the figure with `cairo_pdf()`.
#'
#' @param steps A data frame with one row per step, in pipeline order within each branch:
#'   * `branch` (or the column named by `branch_col`): the trunk value (`trunk`) or a branch name;
#'   * `step`: the step name (box title), unique within a branch;
#'   * `n_out`: cells kept; `n_in` (cells entering, needed for the branch entry boxes) and/or `n_removed` (default
#'     `n_in - n_out`);
#'   * optional `label`: text under the cell count (e.g. the criterion; use `\n` for line breaks);
#'   * optional `excluded`: text in the exclusion box when `removals` has no rows for the step;
#'   * optional `body`: the whole box text, replacing the cell count and `label` (e.g. for a step that flags cells
#'     without removing them).
#' @param removals Optional data frame of removed cells by category: `branch`, `step`, `category`, `n_removed`
#'   (several rows per category are summed).
#' @param branch_col Name of the branch column in `steps` and `removals`.
#' @param trunk Value of the branch column that marks the trunk steps.
#' @param branch_labels Optional named vector of branch display names.
#' @param branch_colors Optional named colors (by branch) for each branch's entry and final boxes.
#' @param entry_boxes Start each branch with a box holding its name and the cells entering it.
#' @param step_order Optional order of the branch step columns; default: the order the steps appear in.
#' @param skip_label Text of the placeholder box of a skipped step.
#' @param category_labels Optional named vector to rename removal categories (e.g.
#'   `c('Damaged / low quality' = 'Damaged / LQ')`).
#' @param titles Optional section titles above the trunk and the branches (length 2; NA to leave one out).
#' @param box_width Width of a step box (plot units; a column is `box_width + col_gap` wide).
#' @param col_gap Gap between columns, where the arrows run.
#' @param title_size,text_size Box title and body text sizes.
#' @param colors Box colors: `box`, `border`, `excluded`, `excluded_border` and `edge`.
#'
#' @return A ggplot. A good PDF size is about 0.65 inches per column wide and 2 inches per branch row high.
#' @export
#' @examples
#' QCFlowChart(sim_qc_steps, sim_qc_removals,
#'             branch_colors = c(Lymphoid = '#1f78b4', Myeloid = '#e31a1c', Stromal = '#33a02c'),
#'             category_labels = c('Damaged / low quality' = 'Damaged / LQ'),
#'             titles = c('All cells', 'Per lineage'))
QCFlowChart <- function(
    steps,
    removals = NULL,
    branch_col = 'branch',
    trunk = 'all',
    branch_labels = NULL,
    branch_colors = NULL,
    entry_boxes = TRUE,
    step_order = NULL,
    skip_label = 'skipped',
    category_labels = NULL,
    titles = NULL,
    box_width = 2.55,
    col_gap = 0.3,
    title_size = 3.1,
    text_size = 2.7,
    colors = list(box = 'grey97', border = 'grey30', excluded = '#fbe3e3', excluded_border = '#b2182b',
                  edge = 'grey30')
){
    st <- .FlowSteps(steps, branch_col)
    rm_df <- if(is.null(removals)) NULL else .FlowRemovals(removals, branch_col)
    pitch <- box_width + col_gap
    excl_gap <- 0.35

    trunk_df <- st[st$branch == trunk, ]
    branches <- unique(st$branch[st$branch != trunk])
    branch_steps <- lapply(branches, function(b) st$step[st$branch == b])
    if(is.null(step_order)) step_order <- .MergeOrder(branch_steps)
    missing <- setdiff(unlist(branch_steps), step_order)
    if(length(missing) > 0) stop('steps missing from step_order: ', paste(missing, collapse = ', '))

    #---------------------------------------------------------#
    # nodes: x, row (0 = trunk, 1.. = branches), title, body, type
    #---------------------------------------------------------#

    step_node <- function(cur, x, row, type){
        body <- if(!is.na(cur$body)) cur$body else
            paste0(.Fmt(cur$n_out), ' cells', ifelse(is.na(cur$label), '', paste0('\n', cur$label)))
        data.frame(x = x, row = row, title = cur$step, body = body, type = type, branch = cur$branch, w = box_width)
    }
    excl_node <- function(cur, x, row){
        if(is.na(cur$n_removed) || cur$n_removed <= 0) return(NULL)
        data.frame(x = x, row = row, title = NA, type = 'excl', branch = cur$branch, w = box_width,
                   body = .ExcludedText(cur, rm_df, category_labels))
    }

    nodes <- list()
    excl <- list()
    trunk_x <- (seq_len(nrow(trunk_df)) - 1) * pitch
    for(i in seq_len(nrow(trunk_df))){
        nodes[[length(nodes) + 1]] <- step_node(trunk_df[i, ], trunk_x[i], 0, 'main')
        excl[[length(excl) + 1]] <- excl_node(trunk_df[i, ], trunk_x[i], 0)
    }

    # branch columns start after the fan-out
    x0 <- if(nrow(trunk_df) > 0) max(trunk_x) + pitch + col_gap else 0
    col_x <- x0 + (seq_len(length(step_order) + entry_boxes) - 1) * pitch
    names(col_x) <- c(if(entry_boxes) '.entry', step_order)
    for(j in seq_along(branches)){
        cur_branch <- branches[j]
        cur_df <- st[st$branch == cur_branch, ]
        if(entry_boxes){
            n_entry <- if(is.na(cur_df$n_in[1])) cur_df$n_out[1] + max(0, cur_df$n_removed[1], na.rm = TRUE) else
                cur_df$n_in[1]
            nodes[[length(nodes) + 1]] <- data.frame(
                x = col_x[['.entry']], row = j, title = .MapLabels(cur_branch, branch_labels),
                body = paste0(.Fmt(n_entry), ' cells'), type = 'entry', branch = cur_branch, w = box_width
            )
        }
        # steps between the branch's first and last step are placeholders when skipped
        present <- match(cur_df$step, step_order)
        for(k in seq(min(present), max(present))){
            cur_step <- step_order[k]
            if(cur_step %in% cur_df$step){
                cur <- cur_df[cur_df$step == cur_step, ]
                nodes[[length(nodes) + 1]] <- step_node(cur, col_x[[cur_step]], j,
                                                        if(k == max(present)) 'final' else 'main')
                excl[[length(excl) + 1]] <- excl_node(cur, col_x[[cur_step]], j)
            } else {
                nodes[[length(nodes) + 1]] <- data.frame(x = col_x[[cur_step]], row = j, title = NA, body = skip_label,
                                                         type = 'skip', branch = cur_branch, w = box_width * 0.8)
            }
        }
    }
    nodes <- do.call(rbind, nodes)
    excl <- do.call(rbind, excl)
    all_nodes <- rbind(nodes, excl)
    all_nodes$h <- 0.2 + 0.2 * (lengths(regmatches(all_nodes$body, gregexpr('\n', all_nodes$body))) + 1 +
                                    ifelse(is.na(all_nodes$title), 0, 1.3))

    #---------------------------------------------------------#
    # rows: spacing from the tallest boxes; trunk at y = 0, centered on the branch rows
    #---------------------------------------------------------#

    main <- all_nodes[all_nodes$type != 'excl', ]
    exc <- all_nodes[all_nodes$type == 'excl', ]
    row_main_h <- tapply(main$h, main$row, max)
    row_excl_h <- tapply(exc$h, exc$row, max)
    max_main <- max(row_main_h)
    max_excl <- if(length(row_excl_h) > 0) max(row_excl_h) else 0
    row_gap <- max_main + excl_gap + max_excl + 0.5
    n_branch <- length(branches)
    row_y <- c(0, ((n_branch - 1) / 2 - seq_len(n_branch) + 1) * row_gap)
    all_nodes$y <- row_y[all_nodes$row + 1]
    # exclusion boxes hang below their row, tops aligned
    is_excl <- all_nodes$type == 'excl'
    all_nodes$y[is_excl] <- all_nodes$y[is_excl] - row_main_h[as.character(all_nodes$row[is_excl])] / 2 - excl_gap -
        all_nodes$h[is_excl] / 2

    #---------------------------------------------------------#
    # edges
    #---------------------------------------------------------#

    main <- all_nodes[!is_excl, ]
    arrows <- list()
    lines <- list()
    # consecutive boxes in each row
    for(r in unique(main$row)){
        cur <- main[main$row == r, ]
        cur <- cur[order(cur$x), ]
        if(nrow(cur) > 1){
            i <- seq_len(nrow(cur) - 1)
            arrows[[length(arrows) + 1]] <- data.frame(x = cur$x[i] + cur$w[i] / 2, xend = cur$x[i + 1] - cur$w[i + 1] / 2,
                                                       y = cur$y[i], yend = cur$y[i])
        }
    }
    # fan-out from the last trunk step to the first box of each branch
    if(nrow(trunk_df) > 0 && n_branch > 0){
        fan_x <- max(trunk_x) + box_width / 2 + (x0 - box_width / 2 - max(trunk_x) - box_width / 2) / 2
        first <- do.call(rbind, lapply(seq_len(n_branch), function(j){
            cur <- main[main$row == j, ]
            cur[which.min(cur$x), ]
        }))
        lines[[length(lines) + 1]] <- data.frame(x = max(trunk_x) + box_width / 2, xend = fan_x, y = 0, yend = 0)
        if(n_branch > 1){
            lines[[length(lines) + 1]] <- data.frame(x = fan_x, xend = fan_x, y = min(first$y), yend = max(first$y))
        }
        arrows[[length(arrows) + 1]] <- data.frame(x = fan_x, xend = first$x - first$w / 2, y = first$y, yend = first$y)
    }
    arrows <- do.call(rbind, arrows)
    lines <- do.call(rbind, lines)
    excl_edges <- NULL
    if(any(is_excl)){
        ex <- all_nodes[is_excl, ]
        parent <- main[match(paste(ex$row, ex$x), paste(main$row, main$x)), ]
        excl_edges <- data.frame(x = ex$x, xend = ex$x, y = parent$y - parent$h / 2, yend = ex$y + ex$h / 2)
    }

    #---------------------------------------------------------#
    # draw
    #---------------------------------------------------------#

    all_nodes$fill <- colors$box
    all_nodes$border <- colors$border
    all_nodes$linetype <- 'solid'
    all_nodes$fill[is_excl] <- colors$excluded
    all_nodes$border[is_excl] <- colors$excluded_border
    skip <- all_nodes$type == 'skip'
    all_nodes$fill[skip] <- 'white'
    all_nodes$border[skip] <- 'grey60'
    all_nodes$linetype[skip] <- 'dashed'
    if(!is.null(branch_colors)){
        colored <- all_nodes$type %in% c('entry', 'final') & all_nodes$branch %in% names(branch_colors)
        all_nodes$fill[colored] <- scales::alpha(branch_colors[all_nodes$branch[colored]], 0.35)
    }
    all_nodes$title_y <- all_nodes$y + all_nodes$h / 2 - 0.22
    all_nodes$body_y <- ifelse(is.na(all_nodes$title), all_nodes$y, all_nodes$y - 0.13)
    titled <- all_nodes[!is.na(all_nodes$title), ]
    head <- grid::arrow(length = grid::unit(0.08, 'inches'), type = 'closed')

    p <- ggplot()
    if(!is.null(arrows)){
        p <- p + geom_segment(data = arrows, aes(x = x, xend = xend, y = y, yend = yend), arrow = head,
                              linewidth = 0.4, color = colors$edge)
    }
    if(!is.null(lines)){
        p <- p + geom_segment(data = lines, aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.4,
                              color = colors$edge)
    }
    if(!is.null(excl_edges)){
        p <- p + geom_segment(data = excl_edges, aes(x = x, xend = xend, y = y, yend = yend), arrow = head,
                              linewidth = 0.3, color = colors$excluded_border, linetype = 'dashed')
    }
    p <- p +
        geom_tile(data = all_nodes, aes(x = x, y = y, width = w, height = h, fill = fill, color = border,
                                        linetype = linetype), linewidth = 0.4) +
        geom_text(data = titled, aes(x = x, y = title_y, label = title), fontface = 'bold', size = title_size) +
        geom_text(data = all_nodes, aes(x = x, y = body_y, label = body), size = text_size, lineheight = 0.95) +
        scale_fill_identity() +
        scale_color_identity() +
        scale_linetype_identity() +
        coord_cartesian(clip = 'off') +
        theme_void()

    if(!is.null(titles)){
        top <- max(all_nodes$y + all_nodes$h / 2) + 0.6
        title_df <- data.frame(
            x = c(if(nrow(trunk_df) > 0) mean(trunk_x) else NA, if(n_branch > 0) mean(range(main$x[main$row > 0])) else NA),
            label = titles
        )
        title_df <- title_df[!is.na(title_df$x) & !is.na(title_df$label), ]
        p <- p + geom_text(data = title_df, aes(x = x, y = top, label = label), fontface = 'italic',
                           size = title_size + 0.3)
    }
    p
}

# steps table with branch, step, n_in, n_out, n_removed, label, excluded
.FlowSteps <- function(steps, branch_col){
    st <- as.data.frame(steps)
    if(!(branch_col %in% colnames(st))) stop("steps needs a '", branch_col, "' column.")
    missing <- setdiff(c('step', 'n_out'), colnames(st))
    if(length(missing) > 0) stop('steps needs the columns: ', paste(missing, collapse = ', '))
    st$branch <- as.character(st[[branch_col]])
    for(cc in c('n_in', 'n_removed')) if(!(cc %in% colnames(st))) st[[cc]] <- NA_real_
    for(cc in c('label', 'excluded', 'body')) if(!(cc %in% colnames(st))) st[[cc]] <- NA_character_
    st$n_removed <- ifelse(is.na(st$n_removed), st$n_in - st$n_out, st$n_removed)
    if(any(duplicated(st[c('branch', 'step')]))) stop('step names must be unique within a branch.')
    st
}

.FlowRemovals <- function(removals, branch_col){
    rm_df <- as.data.frame(removals)
    missing <- setdiff(c(branch_col, 'step', 'category', 'n_removed'), colnames(rm_df))
    if(length(missing) > 0) stop('removals needs the columns: ', paste(missing, collapse = ', '))
    rm_df$branch <- as.character(rm_df[[branch_col]])
    rm_df
}

# "-N cells", then one line per category (largest first), or the step's `excluded` text
.ExcludedText <- function(cur, rm_df, category_labels){
    header <- paste0('\u2212', .Fmt(cur$n_removed), ' cells')
    rows <- if(is.null(rm_df)) NULL else rm_df[rm_df$branch == cur$branch & rm_df$step == cur$step, ]
    if(!is.null(rows) && nrow(rows) > 0){
        n <- tapply(rows$n_removed, rows$category, sum)
        n <- sort(n, decreasing = TRUE)
        cats <- .MapLabels(names(n), category_labels)
        return(paste0(header, '\n', paste0(cats, ': ', .Fmt(n), collapse = '\n')))
    }
    if(!is.na(cur$excluded)) return(paste0(header, '\n', cur$excluded))
    header
}

.Fmt <- function(x){
    scales::comma(x, accuracy = 1)
}

# merge several step sequences into one order that keeps each sequence's order (ties: first appearance)
.MergeOrder <- function(seqs){
    nodes <- unique(unlist(seqs))
    if(length(nodes) == 0) return(character(0))
    edges <- unique(do.call(rbind, lapply(seqs, function(s){
        if(length(s) < 2) return(NULL)
        cbind(s[-length(s)], s[-1])
    })))
    out <- character(0)
    while(length(nodes) > 0){
        blocked <- if(is.null(edges) || nrow(edges) == 0) character(0) else edges[edges[, 1] %in% nodes, 2]
        ready <- nodes[!(nodes %in% blocked)]
        if(length(ready) == 0) stop('the step orders of the branches conflict; give step_order.')
        out <- c(out, ready[1])
        nodes <- setdiff(nodes, ready[1])
    }
    out
}
