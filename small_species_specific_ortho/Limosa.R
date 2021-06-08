library(tidyverse)
library(GenomicRanges)
library(BSgenome)
library(plyranges)
library(reshape2)
library(RColorBrewer)

species_list <- read_tsv("~/Genomes/Birds/bird_genomes.tsv", col_names = c("species_name", "genome_name"))
work_dir <- "~/temp"
query_dir <- "~/Birds/OrthologueRegions/cluster"
genome_dir <- "~/Genomes/Birds"
repbase_rep <- "RepBase_for_class.fasta"

subfamily_colours <- tibble(subfamily = c("CR1-C", "CR1-Croc", "CR1-E", "CR1-J", "CR1-W", "CR1-X", "CR1-Y", "CR1-Z"),
                            class_colour = c("#FEFA70", "#9D9D9D", "#FED16B", "#98FF6F", "#D363FE", "#6F87FF", "#B56FFF", "#FE6363"))

subfamilies <- read_tsv("extended_fastas/recip_blast/classes.tsv") %>%
  inner_join(subfamily_colours)

# set query species
species_name <- species_list$species_name[67]

query_info <- tibble(names = names(readDNAStringSet(paste0("extended_fastas/", species_name, "_extended_small.fasta")))) %>%
  separate(names, into = c("qseqid", "pident", "class"), sep = "#") %>%
  mutate(qseqid_2 = qseqid,
         pident = as.double(pident)) %>%
  tidyr::separate(qseqid_2, into = c("qseqnames", "ranges"), ":") %>%
  tidyr::separate(ranges, into = c("ranges", "ostrand"), "\\(") %>%
  tidyr::separate(ranges, into = c("ostart", "oend"), "-") %>%
  mutate(ostrand = sub(")", "", ostrand),
         ostart = as.integer(ostart),
         oend = as.integer(oend),
         olength = oend - ostart + 1)

collated_info <- query_info %>% dplyr::select(qseqid)

queries <- query_info %>% dplyr::select(qseqid, olength)

for(i in c(29, 30, 67, 102)){
  
  # set subject species
  s_species <- species_list$species_name[i]
  
  # skip if subject species is query species
  if(species_name == s_species){
    collated_info <- collated_info %>% mutate(!! s_species := 0)
    next()
  }
  
  # read in blast, remove small hits, label flank, and determine start, end and strand
  other_blast <- read_tsv(paste0("extended_fastas/blast_out/", species_name, "_in_", s_species, "_small.tsv") ,col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "bitscore", "length", "qlen", "slen", "pident")) %>%
    mutate(qseqid = sub("#.*", "", qseqid))
  
  old <- other_blast %>%
    filter(qstart <= 450, qend >= qlen - 450) %>%
    dplyr::group_by(qseqid) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  other_filtered <- other_blast %>%
    filter(!qseqid %in% old$qseqid) %>%
    mutate(olen = qlen - 1200) %>%
    filter(qstart <= 450 | qend >= qlen - 450) %>%
    mutate(sstrand = ifelse(sstart < send, "+", "-"),
           start = ifelse(sstart < send, sstart, send),
           end = ifelse(sstart < send, send, sstart),
           strand = case_when(qstart <= 450 ~ "5",
                              qend >= qlen - 450 ~ "3"),
           strand_filter = paste0(qseqid, "#", sstrand)) %>%
    arrange(qseqid, start) %>%
    mutate(sstart = start, send = end) %>%
    dplyr::select(-start, -end)
  
  # select single hits
  hit_no <- tibble(strand_filter = names(table(other_filtered$strand_filter)), n = as.integer(table(other_filtered$strand_filter))) %>%
    filter(n == 2)
  hit_no_s <- tibble(qseqid = names(table(sub("#.", "", hit_no$strand_filter))), n = as.integer(table(sub("#.", "", hit_no$strand_filter)))) %>%
    filter(n == 1)
  hit_no <- hit_no %>%
    filter(sub("#.", "", strand_filter) %in% hit_no_s$qseqid)
  
  # work out distances between query stop start and subject stop start
  other_further_filtered_5_3 <- other_filtered %>%
    filter(strand_filter %in% hit_no$strand_filter, sstrand == "+") %>%
    dplyr::select(-strand_filter) %>%
    filter((qseqid == lead(qseqid) & strand == "5" & lead(strand == "3")) |
             (qseqid == lag(qseqid) & strand == "3" & lag(strand == "5"))) %>%
    mutate(sdist = case_when(qseqid == lead(qseqid) ~ lead(sstart) - send + 1,
                             qseqid == lag(qseqid) ~ sstart - lag(send) + 1),
           qdist = case_when(qseqid == lead(qseqid) ~ lead(qstart) - qend + 1,
                             qseqid == lag(qseqid) ~ qstart - lag(qend) + 1),
           diff = qdist - sdist)
  other_further_filtered_3_5 <- other_filtered %>%
    filter(strand_filter %in% hit_no$strand_filter, sstrand == "-") %>%
    dplyr::select(-strand_filter) %>%
    filter((qseqid == lead(qseqid) & strand == "3" & lead(strand == "5")) |
             (qseqid == lag(qseqid) & strand == "5" & lag(strand == "3"))) %>%
    mutate(sdist = case_when(qseqid == lead(qseqid) ~ lead(sstart) - send + 1,
                             qseqid == lag(qseqid) ~ sstart - lag(send) + 1),
           qdist = case_when(qseqid == lead(qseqid) ~ qstart - lead(qend) + 1,
                             qseqid == lag(qseqid) ~ lag(qstart) - qend + 1),
           diff = qdist - sdist)
  
  # find old by weeding out overlaps (potential deletions)
  other_further_filtered_old_5_3 <- other_further_filtered_5_3 %>%
    filter((strand == 5 & qend > 616 & lead(qseqid) == qseqid) |
             (lag(qseqid) == qseqid & lag(strand) == 5 & lag(qend) > 616) |
             (strand == 3 & qstart < qlen - 616 & lag(qseqid) == qseqid) |
             (lead(qseqid) == qseqid & lead(strand) == 3 & lead(qstart) < qlen - 616)) %>%
    mutate(start = ifelse(qseqid == lead(qseqid), sstart, lag(sstart)),
           end = ifelse(qseqid == lead(qseqid), lead(send), send)) %>%
    group_by(qseqid) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    mutate(state = 0) %>%
    dplyr::select(sseqid, start, end, qseqid, sstrand, state)
  
  # find old by weeding out overlaps (potential deletions)
  other_further_filtered_old_3_5 <- other_further_filtered_3_5 %>%
    filter((strand == 5 & qend > 616 & lag(qseqid) == qseqid) |
             (strand == 3 & lead(qend) > 616 & lead(qseqid) == qseqid) |
             (strand == 3 & qstart < qlen - 616 & lead(qseqid) == qseqid) |
             (strand == 5 & lag(qstart) < qlen - 616 & lag(qseqid) == qseqid)) %>%
    mutate(start = ifelse(qseqid == lead(qseqid), sstart, lag(sstart)),
           end = ifelse(qseqid == lead(qseqid), lead(send), send)) %>%
    group_by(qseqid) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    mutate(state = 0) %>%
    dplyr::select(sseqid, start, end, qseqid, sstrand, state)
  
  # determine old vs new
  other_further_filtered <- rbind(other_further_filtered_5_3, other_further_filtered_3_5) %>%
    filter(sdist > -1000, sdist < 1000, !qseqid %in% other_further_filtered_old_3_5$qseqid, !qseqid %in% other_further_filtered_old_5_3$qseqid) %>%
    mutate(state = case_when((sdist <= 16 & qdist >= olen - 24) ~ 1,
                             (qdist >= -16 & qdist <= 16 & sdist >= 90) ~ 0,
                             TRUE ~ -1)) %>%
    mutate(start = ifelse(qseqid == lead(qseqid), sstart, lag(sstart)),
           end = ifelse(qseqid == lead(qseqid), lead(send), send)) %>%
    group_by(qseqid) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    dplyr::select(sseqid, start, end, qseqid, sstrand, state)
  
  old <- old %>%
    mutate(sstrand = ifelse(sstart < send, "+", "-"),
           start = ifelse(sstart < send, sstart, send),
           end = ifelse(sstart < send, send, sstart),
           state = 0) %>%
    dplyr::select(sseqid, start, end, qseqid, sstrand, state)
  
  present_bound_resolved <- rbind(old, other_further_filtered_old_3_5, other_further_filtered_old_5_3, other_further_filtered) %>%
    arrange(qseqid)
  
  # write resolved locations to file
  write_tsv(x = present_bound_resolved, file = paste0("extended_fastas/resolved_locations/small_", species_name, "_in_", s_species, ".tsv"), col_names = F)
  
  # select relevant info for further resolution
  pre_info <- present_bound_resolved %>%
    select(qseqid, state) %>%
    full_join(queries) %>%
    mutate(state = ifelse(is.na(state), -1, state)) %>%
    dplyr::rename(!! s_species := state) %>%
    dplyr::select(-olength)
  
  # collate with previous info
  collated_info <- full_join(pre_info, collated_info)
  
}

# missing = -1, new = 1, ancestral = 0
# join all data
collated_info <- inner_join(collated_info, query_info)

# write to file
write_tsv(x = collated_info, file = paste0("extended_fastas/collated_results/", species_name, "_small_collated_results.tsv"), col_names = T)

# plotting
species_name <- species_list$species_name[67]
collated_info <- read_tsv(paste0("extended_fastas/collated_results/", species_name, "_small_collated_results.tsv")) %>%
  mutate(round_jc = floor(pident * 250)/250) %>%
  filter(olength >= 100)

# determine whether ancestral or recent (needs to be adjusted for each species)
collated_info_2 <- collated_info %>%
  filter(Calidris_pugnax %in% c(0, 1) & Scolopax_mira %in% c(0, 1) & Limosa_lapponica_baueri %in% c(0, 1)) %>%
  mutate(state = case_when(
    Calidris_pugnax == 0 & Scolopax_mira == 0 & Calidris_pygmaea == 0 ~ 0, # 
    Calidris_pugnax == 0 & Scolopax_mira == 0 & Calidris_pygmaea == 1 ~ 1, #
    Calidris_pugnax == 0 & Scolopax_mira == 1 & Calidris_pygmaea == 0 ~ 2, # 
    Calidris_pugnax == 1 & Scolopax_mira == 0 & Calidris_pygmaea == 0 ~ 3, #
    Calidris_pugnax == 0 & Scolopax_mira == 1 & Calidris_pygmaea == 1 ~ 4, # 
    Calidris_pugnax == 1 & Scolopax_mira == 0 & Calidris_pygmaea == 1 ~ 5, # 
    Calidris_pugnax == 1 & Scolopax_mira == 1 & Calidris_pygmaea == 0 ~ 6, # 
    Calidris_pugnax == 1 & Scolopax_mira == 1 & Calidris_pygmaea == 1 ~ 7, # 
    TRUE ~ -1 # unclear
  )
  ) %>%
  arrange(round_jc)

data.frame(table(collated_info_2$state))

# count number of each state
counts <- data.frame(table(collated_info_2$state, collated_info_2$round_jc))

# fix up counts
counts_2 <- counts %>%
  as_tibble() %>%
  dplyr::mutate(Var2 = as.double(as.character(Var2)), Var1 = as.character(Var1)) %>%
  dplyr::arrange(Var1, Var2) %>%
  dplyr::rename(state = Var1, round_jc = Var2)

# # select only resolved
counts <- counts_2 %>% filter(state != -1)

counts_pie <- counts_2 %>%
  group_by(state) %>%
  mutate(sum_Freq = sum(Freq)) %>%
  ungroup() %>%
  dplyr::select(state, sum_Freq) %>%
  base::unique()

ggplot(counts_pie, aes(x="", y=sum_Freq, fill=state)) + geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0)

counts_max <- counts %>%
  group_by(round_jc) %>%
  mutate(sumFreq = sum(Freq)) %>%
  ungroup()

divergence_plot <- ggplot(data = counts, aes(round_jc, Freq, fill = state)) + geom_bar(stat = "identity") +
  scale_x_continuous(name = "Jukes-Cantor distance", expand = c(0,0), limits = c(-0.002, 0.402)) + theme_bw() +
  scale_fill_manual(labels = c("Ancestral", "Since Limosa",  "Since Scolopax", "Since Scolopax and Limosa", "Since Calidris pugnax"),
                    values = brewer.pal(n = 6, name = 'BrBG')[c(1:3, 5:6)]) +
  scale_y_continuous(limits = c(0, 350), expand = c(0,0), name = "Insertions") +
  ggtitle(label = gsub("_", " ", species_name)) +
  theme(plot.title = element_text(family = "Arial", face = "italic", hjust = 0.5, size = 14),
        axis.title.x = element_text(family = "Arial", size = 12),
        axis.title.y = element_text(family = "Arial", size = 12),
        axis.text.x = element_text(family = "Arial", size = 11),
        axis.text.y = element_text(family = "Arial", size = 11),
        legend.text = element_text(family = "Arial", size = 11),
        legend.title = element_blank())

divergence_plot

ggsave(divergence_plot, filename = paste0("extended_fastas/plots/", species_name, "_small_insertion_timeline.svg"), device = "svg", height = 10, width = 18, units = "cm")
