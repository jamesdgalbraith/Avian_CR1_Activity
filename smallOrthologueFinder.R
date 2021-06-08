#! /usr/bin/env Rscript

# import libraries
library(tidyverse)
library(GenomicRanges)
library(BSgenome)
library(plyranges)

species_list <- read_tsv("~/Genomes/Birds/bird_genomes.tsv", col_names = c("species_name", "genome_name"))
genome_dir <- "~/Genomes/Birds"

for(j in 1:nrow(species_list)){
  
  # set species and genome names
  species_name <- species_list$species_name[j]
  genome_name <- species_list$genome_name[j]
  print(species_name)
  # set path to query
  if(!file.exists(paste0("divergence/best_hits/", species_name, "_small_best_hit_coords.tsv"))){next()}
  # import genome
  print(paste0("Importing ", species_name))
  genome_seq <- readDNAStringSet(paste0("~/Genomes/Birds/", species_name, "/", genome_name))
  names(genome_seq) <- sub(" .*", "", names(genome_seq))
  genome_idx <- tibble(seqnames = names(genome_seq), scaffold_length = width(genome_seq))
  
  best_hits <- read_tsv(paste0("divergence/best_hits/", species_name, "_small_best_hit_coords.tsv"),
                        col_names = c("seqnames", "start", "end", "strand", "jc_dist", "class")) %>%
    mutate(names = paste0(seqnames, ":", start, "-", end, "(", strand, ")#", jc_dist, "#", class)) %>%
    filter(end - start  + 1 >= 100) %>%
    dplyr::rename(ostart = start, oend = end) %>%
    inner_join(genome_idx)
  
  best_hit_ranges <- best_hits %>%
    mutate(start = ostart - 600, end = oend + 600) %>%
    filter(start > 0, end <= scaffold_length) %>%
    plyranges::as_granges()
  
   # get 3' repeat seq and name
  best_hits_seq <- getSeq(genome_seq, best_hit_ranges)
  names(best_hits_seq) <- paste0(best_hit_ranges$names)
  
  # determine suitibility
  best_hits_acceptable <- alphabetFrequency(best_hits_seq) %>%
    as_tibble() %>%
    mutate(seqnames = names(best_hits_seq)) %>%
    filter(N < 300) %>%
    dplyr::rename(names = seqnames)
  
  best_hits_final_ranges <- best_hits %>%
    filter(names %in% best_hits_acceptable$names) %>%
    mutate(start = ostart - 600, end = oend + 600) %>%
    as_granges()
  
  best_hits_final_seq <- getSeq(genome_seq, best_hits_final_ranges)
  names(best_hits_final_seq) <- best_hits_final_ranges$names
  
  # write hit seq to file
  writeXStringSet(best_hits_final_seq, paste0("extended_fastas/", species_name, "_extended_small.fasta"))
  
}


for(j in c(9)){
  # set query species name
  species_name <- species_list$species_name[j]
  
  print(paste0("echo ", species_name))
  
  
  for(i in c(10:13)){

    # set subject species and genome names
  
    s_species <- species_list$species_name[i]
    s_genome <- species_list$genome_name[i]
  
    if(species_name != s_species){
      print(paste0("echo ", species_name, " in ", s_species))
      # print(paste0("Blasting ", species_name, " in ", s_species))
      print(paste0("blastn -query extended_fastas/", species_name, "_extended_small.fasta -db ~/Genomes/Birds/", s_species, "/", s_genome,
                   " -outfmt \"6 qseqid sseqid qstart qend sstart send bitscore length qlen slen pident\" -task dc-megablast -num_threads 3 -out extended_fastas/blast_out/", species_name, "_in_", s_species, "_small.tsv -max_target_seqs 1 -evalue 1e-50"))
    }
  }
}

 
# # Species lists
# Amazona_collaria c(7, 6, 96, 5, 72)
# Anas_platyrhynchos c(9:13)
# Apteryx_australis c(16, 19, 18, 17, 35, 49)
# Bubo_bubo c(24, 25, 109, 84, 22, 116)
# Calidris_pugnax c(29, 30, 67, 104)
# Charadrius vociferus c(38, 39, 101, 64)
# Gallus_gallus c(58, 87, 115, 69, 71, 46)
# Taeniopygia_guttata c(112, 105, 86, 107, 118, 60, 112, 115, 107, 119)
 
best_hits_transformed <- best_hits %>%
  dplyr::select(class, length, pident) %>%
  mutate(length = ceiling(length/50)*50, pident = ceiling(pident*2)/2)



best_hits_transformed <- data.frame(base::table(best_hits_transformed$pident, best_hits_transformed$length)) %>%
  as_tibble() %>%
  dplyr::mutate(Var1 = as.double(as.character(Var1)), Var2 = as.double(as.character(Var2))) %>%
  dplyr::rename(pident = Var1, length= Var2) %>%
  dplyr::arrange(length, pident)

ggplot(best_hits_transformed, aes(pident, length, fill = Freq)) + geom_tile()

ggplot(best_hits, aes(pident)) + geom_density()

ggplot() + geom_point(data = best_hits, mapping = aes(x = pident, y = length)) + scale_x_continuous(limits = c(90, 100))
ggplot() + geom_density(data = best_hits, mapping = aes(x = 100 - pident)) + scale_x_continuous(limits = c(0, 15))

