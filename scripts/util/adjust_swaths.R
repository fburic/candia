library(tidyverse)
library(optparse)
library(filesstrings)

option_list = list(
    make_option(c("-i", "--input_map"),
                type="character", default=NULL,
                help="input file name",
                metavar="character"),
    make_option(c("-o", "--output_map"),
                type="character", default=NULL,
                help="Output file name",
                metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
options = parse_args(opt_parser)

if (is.null(options$input_map) | is.null(options$output_map)){
    print_help(opt_parser)
    stop("Missing arguments", call.=F)
}

input_map = read_csv(options$input_map)

swaths =
    input_map %>%
    distinct(prec_isolation_window_start, prec_isolation_window_end) %>%
    arrange(prec_isolation_window_start)

non_overlapping_swaths =
    swaths %>%
    mutate(lower = prec_isolation_window_start,
           upper = prec_isolation_window_end,
           lead_lower = lead(lower)) %>%
    group_by(lower) %>%
    mutate(swath_upper_adjusted = round(mean(c(lead_lower, upper), na.rm = T), digits=2)) %>%
    ungroup() %>%
    mutate(swath_lower_adjusted = lag(swath_upper_adjusted))

non_overlapping_swaths$swath_lower_adjusted[1] = non_overlapping_swaths$lower[1]

map_with_adjusted_swaths <- input_map %>% left_join(non_overlapping_swaths)

map_with_adjusted_swaths <-
    map_with_adjusted_swaths %>%
    select(-c(lead_lower,
              lower,
              upper,
              prec_isolation_window_start,
              prec_isolation_window_end))

readr::write_csv(map_with_adjusted_swaths, options$output_map)

# Save intervals for quick reference
readr::write_csv(map_with_adjusted_swaths %>%
                     select(swath_lower_adjusted, swath_upper_adjusted) %>%
                     distinct(swath_lower_adjusted, swath_upper_adjusted) %>%
                     arrange(swath_lower_adjusted),
                 paste0(options$output_map, '.intervals'))
