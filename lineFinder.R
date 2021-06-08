library(tidyverse)
library(GenomicRanges)
library(BSgenome)
library(plyranges)

species_list <- read_tsv("~/Genomes/Birds/bird_genomes.tsv", col_names = c("species_name", "genome_name"))

for(i in 1:nrow(species_list)){
  
  print(species_list$species_name[i])
  
  # read in consensus of CARP output
  consensus_seq <- Biostrings::readDNAStringSet(filepath = paste0("CARP/", species_list$species_name[i], "/", species_list$species_name[i], "_consensus.fasta"))
  names(consensus_seq) <- sub(" .*", "", names(consensus_seq))
  
  # extract sequences large enough to contain both EN and RT
  big_ranges <- tibble(seqnames = names(consensus_seq), start = 1, end = width(consensus_seq)) %>%
    filter(end > 1800, end < 6000)
  
  if(nrow(big_ranges) == 0){
    next()
  }
  
  big_ranges <- as_granges(big_ranges)
  
  big_seq <-getSeq(consensus_seq, big_ranges)
  names(big_seq) <- seqnames(big_ranges)
  
  Biostrings::writeXStringSet(big_seq, paste0("CARP/", species_list$species_name[i], "/", species_list$species_name[i], "_big_consensus.fasta"))
  if(!file.exists(paste0("CARP/", species_list$species_name[i], "/", species_list$species_name[i], "_rps.out"))){
    system(command = paste0("rpstblastn -query CARP/", species_list$species_name[i], "/", species_list$species_name[i], "_big_consensus.fasta -db ~/Databases/localrpsb/transposons/transposons -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs stitle\" -evalue 0.01 -out CARP/", species_list$species_name[i], "/", species_list$species_name[i], "_rps.out -num_threads 12"))
  }
  
  rps_out <- read_tsv(paste0("CARP/", species_list$species_name[i], "/", species_list$species_name[i], "_rps.out"),
                      col_names = c("seqnames", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                                    "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs", "stitle"))
  
  if(!file.exists(paste0("CARP/", species_list$species_name[i], "/", species_list$species_name[i], "_repbase.out"))){
    system(command = paste0("blastn -task dc-megablast -query CARP/", species_list$species_name[i], "/", species_list$species_name[i], "_big_consensus.fasta -db ~/Databases/RepBase/Separate15_2_20/LINEs.fa -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs stitle\" -evalue 0.01 -out CARP/", species_list$species_name[i], "/", species_list$species_name[i], "_repbase.out -num_threads 12"))
  }
  
  
  repbase_out <- read_tsv(paste0("CARP/", species_list$species_name[i], "/", species_list$species_name[i], "_repbase.out"),
                          col_names = c("seqnames", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                                        "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs", "stitle"))
  
  if(nrow(rps_out) == 0){
    next
  }
  
  rps_out <- rps_out%>%
    separate(stitle, into = c("code", "name", "description"), sep = ", ") %>%
    mutate(strand = if_else(qstart < qend, "+", "-"))
  
  RT <- rps_out %>%
    filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_3", "TERT"), (send - sstart + 1) >= 0.3 * slen)
  EN <- rps_out %>%
    filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"), (send - sstart + 1) >= 0.3 * slen)
  
  potential_LINEs <- rps_out %>%
    dplyr::filter(seqnames %in% RT$seqnames, seqnames %in% EN$seqnames) %>%
    dplyr::select(seqnames, qlen, strand) %>%
    dplyr::mutate(start = 1) %>%
    dplyr::rename(end = qlen) %>%
    base::unique()
  
  if(nrow(potential_LINEs) > 1){
    
    potential_LINEs <- potential_LINEs %>%
      as_granges()
    
    potential_LINEs_seq <- getSeq(consensus_seq, potential_LINEs)
    names(potential_LINEs_seq) <- seqnames(potential_LINEs)
    writeXStringSet(potential_LINEs_seq, paste0("CARP/", species_list$species_name[i], "/potential_LINEs.fa"))
    system(paste0("makeblastdb -dbtype nucl -in CARP/", species_list$species_name[i], "/potential_LINEs.fa -out CARP/temp"))
    
    potential_LINEs_blast <- read_tsv(system(paste0("blastn -dust yes -num_threads 12 -db CARP/temp -query CARP/", species_list$species_name[i], "/potential_LINEs.fa -outfmt \"6 qseqid sseqid qstart qend sstart send pident qcovs bitscore length mismatch evalue qlen slen\" -task blastn"), intern = TRUE), col_names = c("qseqid", "seqnames", "qstart", "qend", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue", "qlen", "slen"))
    
    redundant_potential_LINEs <- potential_LINEs_blast %>%
      dplyr::filter(seqnames != qseqid, pident >= 94, qlen < slen, qcovs > 50) %>%
      dplyr::select(seqnames) %>%
      base::unique()
    
    potential_LINEs <- potential_LINEs %>%
      filter(!seqnames %in% redundant_potential_LINEs$seqnames)
    potential_LINEs_seq <- getSeq(consensus_seq, potential_LINEs)
    names(potential_LINEs_seq) <- seqnames(potential_LINEs)
    writeXStringSet(potential_LINEs_seq, paste0("CARP/", species_list$species_name[i], "/", species_list$species_name[i], "_LINEs.fa"))
    
  } else if(nrow(potential_LINEs) == 1){
    
    potential_LINEs <- potential_LINEs %>%
      as_granges()
    
    potential_LINEs_seq <- getSeq(consensus_seq, potential_LINEs)
    names(potential_LINEs_seq) <- seqnames(potential_LINEs)
    writeXStringSet(potential_LINEs_seq, paste0("CARP/", species_list$species_name[i], "/", species_list$species_name[i], "_LINEs.fa"))
    
  }
  
}

absent <- tibble(species = character())

for(i in 1:nrow(species_list)){
  
  print(species_list$species_name[i])
  
  if(!file.exists(paste0("CARP/", species_list$species_name[i], "/", species_list$species_name[i], "_LINEs.fa"))){
    a <- tibble(species = species_list$species_name[i])
    absent <- rbind(absent, a)
    next()
  }
  
  system(paste0("blastn -query CARP/", species_list$species_name[i], "/", species_list$species_name[i], "_LINEs.fa -db ~/Genomes/Birds/", species_list$species_name[i], "/", species_list$genome_name[i], " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen\" -out CARP/", species_list$species_name[i], "/", species_list$species_name[i], "_LINEs.out -num_threads 12"))
  
  blast_search <- read_tsv(file = paste0("CARP/", species_list$species_name[i], "/", species_list$species_name[i], "_LINEs.out"), col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs", "qlen", "slen"))
  
  filtered_blast_search <- blast_search %>%
    dplyr::filter(pident >= 94) %>%
    dplyr::filter(length >= 1000) %>%
    dplyr::mutate(start = ifelse(sstart < send, sstart, send),
                  end = ifelse(sstart > send, sstart, send),
                  strand = ifelse(sstart < send, "+", "-")) %>%
    dplyr::rename(seqnames = sseqid)
  
  genome_seq <- readDNAStringSet(paste0("~/Genomes/Birds/", species_list$species_name[i], "/", species_list$genome_name[i]))
  names(genome_seq) <- sub( " .*", "", names(genome_seq))
  
  hits <- filtered_blast_search %>%
    dplyr::select(qseqid, qlen) %>%
    dplyr::mutate(start = 1) %>%
    dplyr::rename(end = qlen, seqnames = qseqid) %>%
    base::unique()
  
  hits_ranges <- plyranges::as_granges(hits)
  
  consensus_seq <- readDNAStringSet(paste0("CARP/", species_list$species_name[i], "/", species_list$species_name[i], "_LINEs.fa"))
  
  consensus_seq <- getSeq(consensus_seq, hits_ranges)
  
  names(consensus_seq) <- seqnames(hits_ranges)
  
  for(j in 1:nrow(hits)){
    
    select_filtered_blast_search <- filtered_blast_search %>%
      dplyr::filter(qseqid == hits$seqnames[j]) %>%
      dplyr::arrange(-bitscore) %>%
      dplyr::select(seqnames, start, end, strand, slen) %>%
      dplyr::mutate(start = ifelse(strand == "+", start - 2500, start - 1200),
                    end = ifelse(strand == "+", end + 1200, end + 2500),
                    start = ifelse(start < 1, 1, start),
                    end = ifelse(end > slen, slen, end))
    
    if(nrow(select_filtered_blast_search) > 30){
      select_filtered_blast_search <- select_filtered_blast_search %>%
        dplyr::slice(1:30)
    } else if(nrow(select_filtered_blast_search) < 10){
      
      select_filtered_blast_search <- blast_search %>%
        dplyr::filter(qseqid == hits$seqnames[j]) %>%
        dplyr::filter(pident >= 94, length > 100) %>%
        dplyr::mutate(start = ifelse(sstart < send, sstart, send),
                      end = ifelse(sstart > send, sstart, send),
                      strand = ifelse(sstart < send, "+", "-")) %>%
        dplyr::rename(seqnames = sseqid) %>%
        dplyr::arrange(-bitscore) %>%
        dplyr::select(seqnames, start, end, strand, slen) %>%
        dplyr::mutate(start = ifelse(strand == "+", start - 2500, start - 1200),
                      end = ifelse(strand == "+", end + 1200, end + 2500),
                      start = ifelse(start < 1, 1, start),
                      end = ifelse(end > slen, slen, end)) %>%
        dplyr::slice(1:10)
      
    }
    
    select_filtered_blast_search_ranges <- as_granges(select_filtered_blast_search) %>%
      plyranges::reduce_ranges_directed()
    
    select_filtered_blast_search_seq <- Biostrings::getSeq(genome_seq, select_filtered_blast_search_ranges)
    names(select_filtered_blast_search_seq) <- paste0(seqnames(select_filtered_blast_search_ranges), ":", ranges(select_filtered_blast_search_ranges),
                                                      "(", strand(select_filtered_blast_search_ranges), ")")
    
    
    select_filtered_blast_search_seq <- c(consensus_seq[j], select_filtered_blast_search_seq)
    writeXStringSet(select_filtered_blast_search_seq, "CARP/temp.fa")
    
    system(paste0("mafft --localpair --adjustdirection --thread 12 CARP/temp.fa > CARP/curation/", species_list$species_name[i], "_", sub("_consensus","",names(consensus_seq)[j]), ".fasta"))
    
  }
  
}


