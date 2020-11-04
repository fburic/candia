library("data.table")
library("tidyverse")
library("feather")
library("filesstrings")
library("optparse")
library("yaml")


main <- function() {
    config <- get_config()
    peak_counts <- as.data.table(read_feather(file.path(config$root_dir, config$time_modes_values)))
    model_index <- read_feather(file.path(config$root_dir, config$model_index))

    print("INFO: Selecting best models according to criteria: unimodality")

    models_with_unimodality_fraction <-
        peak_counts[
            peak_counts[
                npeaks == 1,
                .(n_unimodal_components = .N),
                by = .(model_id)
                ],
            on = .(model_id)
            ][
                ,
                .(unimodal_fraction = n_unimodal_components / .N),
                by = .(model_id)
                ] %>% unique
    print("INFO: Calculated unimodality fraction for all models")

    models_with_mostly_unimodal_elution_profiles <-
        models_with_unimodality_fraction %>%
        inner_join(model_index, by = "model_id") %>%
        group_by(swath_start, rt_window) %>%
        filter(unimodal_fraction == max(unimodal_fraction)) %>%
        ungroup()
    print("INFO: Selected models with predominantly unimodal time modes")

    models_with_mostly_unimodal_elution_profiles %>%
        mutate(swath_start = formatC(swath_start / 100, format = "f", digits = 2)) %>%
        write_csv(file.path(config$root_dir, config$best_models))
    print(paste0("INFO: Wrote ", file.path(config$root_dir, config$best_models)))
}


get_config <- function() {
    option_list = list(
        make_option(c("-c", "--config"),
                    type = "character", default = NULL,
                    help = "YAML experiment config file")
    )
    opt_parser = OptionParser(option_list = option_list)
    options = parse_args(opt_parser)
    if (is.null(options$config)) {
        print_help(opt_parser)
        stop("Missing arguments", call. = F)
    } else {
        config <- read_yaml(options$config)
    }
    return(config)
}


main()
