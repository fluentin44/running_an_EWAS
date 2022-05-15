gen_reportee <- function(voi, raw_results, report_template_filename, output_file_name, covariate_names, results_summary, date_started, p_data){
report_path <- paste0("./reports/", report_template_filename, ".Rmd")
output_path <- paste0("./reports/dma_reports/", Sys.Date(), "/")
rmarkdown::render(
  report_path,
  output_dir = output_path,
  output_file = paste0(output_file_name, ".html"),
  quiet=T,
  params = list(name = str_replace(output_file_name, ".html", ""), # orange in the rmd
                raw_results = raw_results,
                covariates = covariate_names,
                results_summary = results_summary,
                date_started = date_started,
                voi = voi,
                p_data = p_data
  )
)
}
