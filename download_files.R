# download data, run from project base directory

dir.create(file.path("data2", "aggregated_data", "outs", "count"),
           recursive = TRUE,
           showWarnings = FALSE)

download.file("https://www.dropbox.com/s/b9k6bm1idwk8xpe/cloupe.cloupe?dl=1",
              file.path("data2", "aggregated_data", "outs", "count", "cloupe.cloupe"))

download.file("https://www.dropbox.com/s/7mkg6a2jdhwa0ei/so.rds?dl=1",
              file.path("data2", "so.rds"))
