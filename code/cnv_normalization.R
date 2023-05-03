library(ggplot2)

main_folder = '../../HRD_score/data/gene_cnv/'

new_folder = '../../HRD_score/data/gene_cnv2/'
dir.create(new_folder)

txt_files = list.files(path = main_folder, pattern = "*.tsv", recursive = TRUE, full.names = TRUE)

file.copy(txt_files, new_folder)



files = list.files(path = "../../HRD_score/data/allele_specific_cnv/formatted_asc_files/", pattern = NULL, all.files = FALSE,
                   full.names = FALSE, recursive = FALSE,
                   ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

path = '../../HRD_score/data/allele_specific_cnv/formatted_asc_files/'
test_file = files[1]
test = read.csv(paste(path,test_file,sep = ''), sep = '\t')


path_gene = '../../HRD_score/data/gene_cnv2/gene_cnv_tsv/'

files = list.files(path = path_gene, pattern = NULL, all.files = FALSE,
                   full.names = FALSE, recursive = FALSE,
                   ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)




#test$Copy_Number = test$Copy_Number - 2
test$Copy_Number <- ((test$Copy_Number - 2) / (max(test$Copy_Number) - 2))

splitter = function(string){
  split_string <- strsplit(string, split = '\\.')[[1]]
  split_string <- split_string[-length(split_string)]
  remaining_string <- paste(split_string, collapse = '.')
  return (remaining_string)
}

savepath = '../../HRD_score/data/allele_specific_cnv/normalizedCNV/'

savepath_gene = '../../HRD_score/data/gene_cnv2/gene_cnv_norm/'


i = 0
for (file in files){
  df = read.csv(paste(path,file,sep = ''), sep = '\t')
  df$Copy_Number = df$Copy_Number - 2
  name = splitter(file)
  write.csv(df,paste(savepath,name,'.csv',sep = ''))
  if (i%%500 == 0){
    print(paste('Processed',as.character(i),sep=' '))
  }
  i = i + 1
}


i = 0
for (file in files){
  df = read.csv(paste(path_gene,file,sep = ''), sep = '\t')
  df$copy_number = df$copy_number - 2
  name = splitter(file)
  write.csv(df,paste(savepath_gene,name,'.csv',sep = ''))
  if (i%%500 == 0){
    print(paste('Processed',as.character(i),sep=' '))
  }
  i = i + 1
}
