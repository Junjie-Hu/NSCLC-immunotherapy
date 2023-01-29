My_plot_gene_pseudotime <- function (cds_subset, min_expr = NULL, cell_size = 0.75, nrow = NULL, 
    ncol = 1, panel_order = NULL, color_by = "State", trend_formula = "~ sm.ns(Pseudotime, df=3)", 
    label_by_short_name = TRUE, relative_expr = TRUE, vertical_jitter = NULL, 
    horizontal_jitter = NULL) 
{
    f_id <- NA
    Cell <- NA
    if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", 
        "negbinomial.size")) {
        integer_expression <- TRUE
    }
    else {
        integer_expression <- FALSE
        relative_expr <- TRUE
    }
    if (integer_expression) {
        cds_exprs <- exprs(cds_subset)
        if (relative_expr) {
            if (is.null(sizeFactors(cds_subset))) {
                stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
            }
            cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
        }
        cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
    }
    else {
        cds_exprs <- reshape2::melt(as.matrix(exprs(cds_subset)))
    }
    if (is.null(min_expr)) {
        min_expr <- cds_subset@lowerDetectionLimit
    }
    colnames(cds_exprs) <- c("f_id", "Cell", "expression")
    cds_pData <- pData(cds_subset)
    cds_fData <- fData(cds_subset)
    cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
    cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
    if (integer_expression) {
        cds_exprs$adjusted_expression <- cds_exprs$expression
    }
    else {
        cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
    }
    if (label_by_short_name == TRUE) {
        if (is.null(cds_exprs$gene_short_name) == FALSE) {
            cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
            cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
        }
        else {
            cds_exprs$feature_label <- cds_exprs$f_id
        }
    }
    else {
        cds_exprs$feature_label <- cds_exprs$f_id
    }
    cds_exprs$f_id <- as.character(cds_exprs$f_id)
    cds_exprs$feature_label <- factor(cds_exprs$feature_label)
    new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime)
    model_expectation <- genSmoothCurves(cds_subset, cores = 1, 
        trend_formula = trend_formula, relative_expr = T, new_data = new_data)
    colnames(model_expectation) <- colnames(cds_subset)
    library(plyr)
    expectation <- ddply(cds_exprs, .(f_id, Cell), function(x) data.frame(expectation = model_expectation[x$f_id, 
        x$Cell]))
    cds_exprs <- merge(cds_exprs, expectation)
    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
    if (is.null(panel_order) == FALSE) {
        cds_exprs$feature_label <- factor(cds_exprs$feature_label, 
            levels = panel_order)
    }
    q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)
    # if (is.null(color_by) == FALSE) {
    #     q <- q + geom_point(aes_string(color = color_by), size = I(cell_size), 
    #         position = position_jitter(horizontal_jitter, vertical_jitter))
    # }
    # else {
    #     q <- q + geom_point(size = I(cell_size), position = position_jitter(horizontal_jitter, 
    #         vertical_jitter))
    # }
    q <- q + geom_line(aes(x = Pseudotime, y = expectation), 
        data = cds_exprs)
    q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow = nrow, 
        ncol = ncol, scales = "free_y")
    if (min_expr < 1) {
        q <- q + expand_limits(y = c(min_expr, 1))
    }
    if (relative_expr) {
        q <- q + ylab("Relative Expression")
    }
    else {
        q <- q + ylab("Absolute Expression")
    }
    q <- q + xlab("Pseudo-time")
    q <- q + monocle_theme_opts()
    q
}