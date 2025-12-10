install.packages("googledrive")
library(googledrive)

folder_url <- "https://drive.google.com/drive/folders/1RMBPLA1wzeKKzQLEUs8drSFNixa4W2fk"
folder_drib <- drive_get(as_id(folder_url))
files_list <- drive_ls(folder_drib)

# Download each file into data/
for (i in seq_len(nrow(files))) {
  drive_download(
    file = files$id[i],
    path = file.path("data", files$name[i]),
    overwrite = TRUE
  )
}