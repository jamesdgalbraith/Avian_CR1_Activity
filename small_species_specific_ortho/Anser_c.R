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
species_name <- species_list$species_name[12]

query_info <- tibble(qseqid = names(readDNAStringSet(paste0("extended_fastas/", species_name, "_extended_small.fasta")))) %>%
  mutate(qseqid_2 = qseqid) %>%
  tidyr::separate(qseqid_2, into = c("qseqnames", "ranges"), ":") %>%
  tidyr::separate(ranges, into = c("ranges", "ostrand"), "\\(") %>%
  tidyr::separate(ranges, into = c("ostart", "oend"), "-") %>%
  mutate(ostrand = sub(")", "", ostrand),
         ostart = as.integer(ostart),
         oend = as.integer(oend),
         olength = oend - ostart + 1)

collated_info <- query_info %>% dplyr::select(qseqid)

queries <- query_info %>% dplyr::select(qseqid, olength)

for(i in c(9:13, 58)){
  
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
    filter(abs(diff) < 10000) %>%
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
    filter(abs(diff) < 10000) %>%
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
species_name <- species_list$species_name[12]
collated_info <- read_tsv(paste0("extended_fastas/collated_results/", species_name, "_small_collated_results.tsv")) %>%
  filter(olength >= 100)

subfamilies <- read_tsv("extended_fastas/recip_blast/classes.tsv")

reciprocal_blast_out <- read_tsv(paste0("extended_fastas/recip_blast/", species_name, ".out"),
                                 col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) %>%
  mutate(d = mismatch/length) %>%
  mutate(jc_dist = (-3 / 4) * log(1 - (4 * d / 3))) %>%
  dplyr::group_by(qseqid) %>% # group using query sequence
  dplyr::slice(1) %>% # select top hits (highest bitscore)
  dplyr::ungroup() %>%
  mutate(sseqid = sub(";.*", "", sub(".*#", "", as.character(sseqid)))) %>%
  inner_join(subfamilies) %>%
  dplyr::select(qseqid, jc_dist, subfamily, sseqid)

# determine whether ancestral or recent (needs to be adjusted for each species)
collated_info_2 <- collated_info %>%
  mutate(state = case_when(
    (Anas_zonorhyncha == 0 & Anser_indicus == 0 & Anser_brachyrhynchus == 0 & Gallus_gallus == 0) ~ 0, # very ancestral
    (Anas_zonorhyncha == 0 & Anser_indicus == 0 & Anser_brachyrhynchus == 0 & Gallus_gallus == 1) ~ 1, # ancestral
    (Anas_zonorhyncha == 1 & Anser_indicus == 0 & Anser_brachyrhynchus == 0 & Gallus_gallus == 1) ~ 2, # since Anas
    (Anas_zonorhyncha == 1 & Anser_indicus == 0 & Anser_brachyrhynchus == 1) ~ 3, # indicus x brachyrhynchus hybrid
    (Anas_zonorhyncha == 1 & Anser_indicus == 1 & Anser_brachyrhynchus == 0) ~ 4, # cygnoides x brachyrhynchus hybrid
    (Anas_zonorhyncha == 1 & Anser_indicus == 1 & Anser_brachyrhynchus == 1) ~ 5, # since cygnoides+indicus
    TRUE ~ -1 # unclear
  )
  ) %>%
  inner_join(reciprocal_blast_out) %>%
  # filter(state != -1) %>%
  base::unique()

table(collated_info_2$state, collated_info_2$subfamily)


divergence_plot <- ggplot(data = collated_info_2, aes(jc_dist, fill = as.factor(state))) + geom_histogram(binwidth = 0.005) +
  scale_x_continuous(name = "Jukes-Cantor distance", expand = c(0,0), limits = c(-0.005, 0.405)) + theme_bw() +
  scale_fill_manual(labels = c("Ancestral",  "Since Anas", "Anser i. + Anser c.", "Since Anser indicus", "Since Anser brachyrhynchus"),
                    values = brewer.pal(n = 6, name = 'BrBG')[c(1:2, 4:6)]) +
  scale_y_continuous(limits = c(0, 1200), expand = c(0,0), name = "Insertions") +
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

collated_info_3 <- collated_info_2 %>%
  filter(state %in% c(1, 2, 3))

table(collated_info_3$state)

collated_info_3_tbl <- as_tibble(as.data.frame(table(collated_info_3$subfamily))) %>%
  mutate(Var1 = as.character(Var1)) %>%
  filter(Freq/sum(Freq) > 0.01, Freq > 1)

collated_info_3 <- collated_info_3 %>%
  filter(subfamily %in% collated_info_3_tbl$Var1)

subfamily_plot <- ggplot(data = collated_info_3, aes(jc_dist, subfamily, fill = as.factor(subfamily))) + geom_violin(scale = "count") + theme_bw() +
  scale_fill_manual(values = subfamily_colours$class_colour[subfamily_colours$subfamily %in% collated_info_3$subfamily]) +
  scale_y_discrete(expand = c(0,0), limits = rev(levels(as.factor(collated_info_3$subfamily)))) +
  scale_x_continuous(expand = c(0,0), name = "Jukes-Cantor distance", limits = c(0, 0.26)) +
  ggtitle(label = sub("_", " ", species_name), subtitle = paste0(nrow(collated_info_3), " CR1s")) +
  theme(legend.position = "none",
        plot.title = element_text(family = "Arial", face = "italic", hjust = 0.5, size = 14),
        plot.subtitle = element_text(family = "Arial", hjust = 0.5, size = 12),
        axis.title.x = element_text(family = "Arial", size = 12),
        axis.title.y = element_text(family = "Arial", size = 12),
        axis.text.x = element_text(family = "Arial", size = 11),
        axis.text.y = element_text(family = "Arial", size = 11))

subfamily_plot

ggsave(subfamily_plot, filename = paste0("extended_fastas/plots/", species_name, "_subfamily_timeline.svg"), device = "svg", height = 10, width = 18, units = "cm")

for(i in 1:nrow(collated_info_3_tbl)){
  
  to_align_ranges <- collated_info_3_ranges %>%
    filter(sseqid == collated_info_3_tbl$Var1[i])
  
  to_align_seq <- getSeq(genome_seq, to_align_ranges)
  
  names(to_align_seq) <- paste0(seqnames(to_align_ranges), ":", ranges(to_align_ranges), "(", strand(to_align_ranges), ")")
  
  writeXStringSet(to_align_seq, "temp/temp.fa")
  
  system(paste0("mafft --localpair --thread 12 temp/temp.fa > extended_fastas/new_class/alignments/", species_name, "#", collated_info_3_tbl$subfamily[i],
                "#", collated_info_3_tbl$Var1[i],".fasta"))
  
}
