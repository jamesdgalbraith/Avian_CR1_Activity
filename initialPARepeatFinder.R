#! /usr/bin/env Rscript
Sys.setenv(PATH = "/apps/software/blast+/2.7.1/bin:/bin:/usr/bin:/opt/apps/lmod/lmod/libexec:/home/james/bin:/home/james/bin/mummer:/home/james/bin/RepeatMasker:/home/james/bin/minimap2:/home/james/bin/RepeatModeler:/usr/local/go/bin:/home/james/go/bin:/home/james/bin/vmatch:/usr/local/bin/")

# import libraries
library(tidyverse)
library(GenomicRanges)
library(BSgenome)
library(plyranges)

# import species list
species_list <- read_tsv("~/Genomes/Birds/bird_genomes.tsv", col_names = c("species_name", "genome_name"))

for(j in 1:nrow(species_list)){
  
  # set species and genome names
  species_name <- species_list$species_name[j]
  genome_name <- species_list$genome_name[j]
  
  # import genome
  print(paste0("Importing ", species_name))
  genome_seq <- readDNAStringSet(paste0("~/Genomes/Birds/", species_name, "/", genome_name))
  names(genome_seq) <- sub(" .*", "", names(genome_seq))
  genome_idx <- tibble(seqnames = names(genome_seq), scaffold_length = width(genome_seq))
  
  # find CR1s
  # search for already found repeats repeats
  print(paste0("Initial blast ", species_name))
  timestamp()
  blast_out <- read_tsv(system(paste0("blastn -query ~/Birds/OrthologueRegions/cr1Set.fa -db ~/Genomes/Birds/", species_name, "/", genome_name, " -outfmt \"6 qseqid sseqid qstart qend sstart send length qlen slen\" -num_threads 12"), intern = T), col_names = c("qseqid", "seqnames", "qstart", "qend", "sstart", "send", "length", "qlen", "slen"))
  timestamp()
  
  # filter based on size, reduce, further filter on size
  initial_ranges <- blast_out %>%
    filter(length > 1000, qend >= qlen - 18) %>%
    mutate(strand = ifelse(sstart < send, "+", "-"),
           start = ifelse(sstart < send, sstart, send),
           end = ifelse(sstart > send, sstart, send))  %>%
    dplyr::select(seqnames, start, end, strand, slen) %>%
    as_granges() %>%
    reduce_ranges_directed(slen = mean(slen)) %>% # as all collapsed are on the same scaffold mean scaffold length is equal to the true scaffold length
    as_tibble() %>%
    filter(width > 2300, width < 5100) %>%
    arrange(-width) %>%
    mutate(names = paste0(seqnames, ":", start, "-", end, "(", strand, ")")) %>%
    as_granges()
  
  # getting seq to identify domains in
  initial_seq <- getSeq(genome_seq, initial_ranges)
  names(initial_seq) <- paste0(seqnames(initial_ranges), ":", ranges(initial_ranges), "(", strand(initial_ranges), ")")
  writeXStringSet(initial_seq, paste0("orthologueFinder/", species_name, "_initial.fa"))
  
  # run rpstblast to identify domains
  print("rpstblastning")
  timestamp()
  rpst_blast_out <- read_tsv(system(paste0("rpstblastn -query ~/Birds/OrthologueRegions/orthologueFinder/", species_name, "_initial.fa -db ~/Databases/localrpsb/transposons/transposons -outfmt \"6 qseqid sseqid qstart qend sstart send length qlen slen pident stitle\" -num_threads 12"), intern = T), col_names = c("qseqid", "seqnames", "qstart", "qend", "sstart", "send", "length", "qlen", "slen", "pident", "stitle")) %>%
    separate(stitle, into = c("code", "name", "description"), sep = ", ")
  timestamp()
  
  # pull out domains
  rpst_blast_out_RT <- rpst_blast_out %>%
    filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.1 * slen)
  rpst_blast_out_EN <- rpst_blast_out %>%
    filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.1 * slen)
  
  # identify seq containing RT and EN
  intact_cr1_ranges <- rpst_blast_out_EN %>%
    filter(qseqid %in% rpst_blast_out_RT$qseqid) %>%
    dplyr::select(qseqid) %>%
    base::unique() %>%
    tidyr::separate(qseqid, into = c("seqnames", "ranges"), sep = ":") %>%
    tidyr::separate(ranges, into = c("ranges", "strand"), sep = "\\(") %>%
    tidyr::separate(ranges, into = c("start", "end"), sep = "-") %>%
    dplyr::mutate(strand = sub(")", "", strand), start = as.integer(start), end = as.integer(end)) %>%
    as_granges()
  
  if(length(intact_cr1_ranges) <= 500){
    
    # get seq containing RT and EN
    intact_seq <- getSeq(genome_seq, intact_cr1_ranges)
    names(intact_seq) <- paste0(seqnames(intact_cr1_ranges), ":", ranges(intact_cr1_ranges), "(", strand(intact_cr1_ranges), ")")
    writeXStringSet(intact_seq, paste0("orthologueFinder/intact_queries/", species_name, "_intact.fa"))
    
  } else {
    
    # get intact sequences and write to file
    intact_seq <- getSeq(genome_seq, intact_cr1_ranges)
    names(intact_seq) <- paste0(seqnames(intact_cr1_ranges), ":", ranges(intact_cr1_ranges), "(", strand(intact_cr1_ranges), ")")
    writeXStringSet(intact_seq, paste0("orthologueFinder/", species_name, "_intact_initial.fa"))
    
    # cluster intact sequences with vsearch at 99% id
    print(paste0("Clustering intact ", species_name, " sequences"))
    timestamp()
    system(paste0("vsearch --cluster_fast orthologueFinder/", species_name , "_intact_initial.fa --uc orthologueFinder/", species_name,"_intact_full_length.uc --consout orthologueFinder/", species_name,"_intact_full_length.fasta --id 0.94 --threads 0"))
    timestamp()
    
    # read in cluster data
    uc_out <- read_tsv(paste0("orthologueFinder/", species_name, "_intact_full_length.uc"), col_names = c("type", "cluster_n", "length", "pident", "h_strand", "X1", "X2", "Compressed_alignment", "qseqid", "sseqid"))
    
    # select clusters with 3 or more members
    uc_out <- uc_out %>%
      filter(type %in% c("C"), length >= 3)
    
    if(nrow(uc_out) < 2){
      
      intact_cr1_ranges <- rpst_blast_out_EN %>%
        filter(qseqid %in% rpst_blast_out_RT$qseqid) %>%
        dplyr::select(qseqid) %>%
        base::unique() %>%
        tidyr::separate(qseqid, into = c("seqnames", "ranges"), sep = ":") %>%
        tidyr::separate(ranges, into = c("ranges", "strand"), sep = "\\(") %>%
        tidyr::separate(ranges, into = c("start", "end"), sep = "-") %>%
        dplyr::mutate(strand = sub(")", "", strand), start = as.integer(start), end = as.integer(end)) %>%
        dplyr::slice(1:500) %>%
        as_granges()
      
      intact_seq <- getSeq(genome_seq, intact_cr1_ranges)
      names(intact_seq) <- paste0(seqnames(intact_cr1_ranges), ":", ranges(intact_cr1_ranges), "(", strand(intact_cr1_ranges), ")")
      writeXStringSet(intact_seq, paste0("orthologueFinder/intact_queries/", species_name, "_intact.fa"))
      
    } else {
      
      # read in cluster_consensuses
      uc_seq <- readDNAStringSet(paste0("orthologueFinder/", species_name, "_intact_full_length.fasta"))
      
      # create ranges object of clusters of 3 or more seqs
      uc_ranges <- tibble(seqnames = names(uc_seq), start = 1, end = width(uc_seq)) %>%
        dplyr::mutate(qseqid = sub("centroid=", "", seqnames), qseqid = sub(";.*", "", qseqid)) %>%
        inner_join(uc_out) %>%
        as_granges()
      
      # get seq of clusters over 3 and write to file
      uc_cons_seq <- getSeq(uc_seq, uc_ranges)
      names(uc_cons_seq) <- paste0(species_name, "_cluster_", uc_ranges$cluster_n)
      writeXStringSet(uc_cons_seq, paste0("orthologueFinder/intact_queries/", species_name, "_intact.fa"))
      
    }
    
  }
  
  # search for seq containing RT and EN
  print(paste0("Second blast ", species_name))
  timestamp()
  second_blast_out <- read_tsv(system(paste0("blastn -query ~/Birds/OrthologueRegions/orthologueFinder/intact_queries/", species_name, "_intact.fa -db ~/Genomes/Birds/", species_name, "/", genome_name, " -outfmt \"6 qseqid sseqid qstart qend sstart send length pident qlen slen\" -num_threads 12"), intern = T), col_names = c("qseqid", "seqnames", "qstart", "qend", "sstart", "send", "length", "pident", "qlen", "slen")) %>%
    mutate(strand = ifelse(sstart < send, "+", "-"),
           start = ifelse(sstart < send, sstart, send),
           end = ifelse(sstart > send, sstart, send)) %>%
    dplyr::select(-sstart, -send) %>%
    filter(qseqid != paste0(seqnames, ":", start, "-", end, "(", strand, ")"))
  timestamp()
  
  # pull out seqs to search other species
  # must have high id to query, be distant from contig ends, be very close to 3' and be short (50 - 500bp)
  repeats_to_search <- second_blast_out %>%
    filter(pident >= 99, start > 1000, end <= slen - 1000, qend >= qlen - 25) %>% 
    as_granges() %>%
    GenomicRanges::reduce(min.gapwidth = 50) %>%
    as_tibble() %>%
    inner_join(genome_idx) %>%
    arrange(-width) %>%
    filter(width >= 50, width <= 500)
  
  if(nrow(repeats_to_search) == 0){
    print(paste0("None to search in ", species_name))
    next()
  }
  
  # extend flanks for repeat to search and make ranges object
  complete_tbl <- repeats_to_search %>%
    dplyr::mutate(start = start - 1000, end = end + 1000, repeat_name = paste0(species_name, "_repeat_", 1:nrow(repeats_to_search)))
  complete_ranges <- complete_tbl %>%
    select(-width) %>%
    as_granges()
  
  # get sequence of extended queries
  complete_seq <- getSeq(genome_seq, complete_ranges)
  
  # remove queries with high 'N' content, if over 500 remain select top 500
  acceptable_N_freq_ranges <- alphabetFrequency(complete_seq) %>%
    as_tibble() %>%
    mutate(seqnames = complete_tbl$seqnames, start = complete_tbl$start, end = complete_tbl$end, strand = complete_tbl$strand,
           repeat_name = paste0(seqnames, ":", start, "-", end, "(", strand, ")")) %>%
    filter(N < 250) %>%
    as_granges()
  
  if(length(acceptable_N_freq_ranges) == 0){
    print(paste0("Too much N in ", species_name))
    next()
  }
  
  # get sequence of acceptable extended queries
  complete_seq <- getSeq(genome_seq, acceptable_N_freq_ranges)
  names(complete_seq) <- acceptable_N_freq_ranges$repeat_name
  writeXStringSet(complete_seq, paste0("orthologueFinder/", species_name, "_complete.fa"))
  
  # check for TinTs
  print(paste0("TinT blast ", species_name))
  timestamp()
  system(paste0("blastn -task dc-megablast -query orthologueFinder/", species_name, "_complete.fa -db ~/Databases/RepBase/RepBase24.07.fasta -outfmt \"6 qseqid sseqid qstart qend sstart send length pident qlen slen\" -num_threads 12 -out orthologueFinder/", species_name, "_tint.tsv"))
  tint_blast_out <- read_tsv(paste0("orthologueFinder/", species_name, "_tint.tsv"), col_names = c("qseqid", "seqnames", "qstart", "qend", "sstart", "send", "length", "pident", "qlen", "slen")) %>%
    filter(length > qlen - 1800) # allows for 200bp of flanks to be repeat
  timestamp()
  
  # filter out TinTs  
  no_tint_tbl <- acceptable_N_freq_ranges %>%
    filter(!repeat_name %in% tint_blast_out$qseqid) %>%
    as_tibble()
  
  # if none without TinTs skip to next
  if(nrow(no_tint_tbl) == 0){
    print(paste0("Only TinTs in ", species_name))
    next()
  }
  
  # get ranges of TinTless, slice to restrict to top 500
  no_tint_ranges <- no_tint_tbl %>%
    mutate(repeat_name = paste0(species_name, "_", 1:nrow(no_tint_tbl))) %>%
    dplyr::slice(1:500) %>%
    as_granges()
  
  # get sequences of TinT free
  complete_seq <- getSeq(genome_seq, no_tint_ranges)
  names(complete_seq) <- no_tint_ranges$repeat_name
  writeXStringSet(complete_seq, paste0("orthologueFinder/", species_name, "_complete.fa"))
  
}
