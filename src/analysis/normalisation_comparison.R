
stop("Finish the script dickhead!")

references

calc_norm_factors <- function(chain) {
    params <- chain@parameters
    if (is.null(params$phi)) {
        params$nu
    } else {
        params$nu * params$phi
    }
}

lapply(seq_len(nrow(df)),
    function(i) {
        cat(i, "/", nrow(df), "\n")
        file <- df[i, "file"]
        chain <- readRDS(file)
        match <- which(references$data == data)
        ref_norm <- calc_norm_factors(references[[match[[1]], "chain"]])
        my_norm <- calc_norm_factors(chain)
    }
)
