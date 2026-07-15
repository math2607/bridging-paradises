#=###=#=###=#=###=#=###=#=###=
# Functions
#=###=#=###=#=###=#=###=#=###=

rename.locations <- function(align) {
  align <- gsub("\\|+", "_", align)
  align <- gsub(":+", "_", align)
  align <- gsub("\\(+", "_", align)
  align <- gsub("\\)+", "_", align)
  align <- gsub("\\.+", "_", align)
  align <- gsub("\\;+", "_", align)
  align <- gsub("\\t+", "_", align)
  align <- gsub(" +", "_", align)
  align <- gsub("=+", "_", align)
  align <- gsub(",+", "_", align)
  return(align)
}

# genbank sequences by gene
get_genbank_safe <- function(accessions, chunk_size = 50, genename = "Not-Specified") {
  n <- length(accessions)
  # Split accessions into a list of smaller vectors
  indices <- split(1:n, ceiling(seq_along(1:n) / chunk_size))
  
  # Download each chunk and store in a list
  results_list <- lapply(indices, function(idx) {
    base::message(paste("Downloading ", genename, " chunk...", min(idx), "to", max(idx)))
    Sys.sleep(0.5) # A tiny pause to keep NCBI happy
    return(read.GenBank(accessions[idx], species.names = TRUE))
  })
  
  # Combine the DNAbin list objects into one
  return(do.call(base::c, results_list))
}

rename.tips.align <- function(align) {
  align <- gsub("\\|", "_", align)
  align <- gsub(":", "__", align)
  align <- gsub("\\(", "___", align)
  align <- gsub("\\)", "____", align)
  align <- gsub("\\.", "", align)
  return(align)
}


retrieve_names_wtlocation <- function(x) {
  # Detecta a posição do último "|"
  last_pipe <- str_locate_all(x, fixed("|")) |> 
    lapply(function(mat) if (nrow(mat) > 0) mat[nrow(mat), "start"] else NA_integer_)
  
  # Para cada string, processa a remoção do ":" após o último "|"
  mapply(function(label, pipe_pos) {
    if (is.na(pipe_pos)) return(label)  # Nenhum "|"
    
    after_pipe <- str_sub(label, 1, pipe_pos-1)
    return(after_pipe)
  }, x, last_pipe, USE.NAMES = FALSE)
}


all_locations_clean <- function(x) {
  # Detecta a posição do último "|"
  last_pipe <- str_locate_all(x, fixed("|")) |> 
    lapply(function(mat) if (nrow(mat) > 0) mat[nrow(mat), "start"] else NA_integer_)
  
  # Para cada string, processa a remoção do ":" após o último "|"
  mapply(function(label, pipe_pos) {
    if (is.na(pipe_pos)) return(label)  # Nenhum "|"
    
    after_pipe <- str_sub(label, pipe_pos + 1)
    return(after_pipe)
  }, x, last_pipe, USE.NAMES = FALSE)
}

remove_colon_after_last_pipe <- function(x) {
  # Detecta a posição do último "|"
  last_pipe <- str_locate_all(x, fixed("|")) |> 
    lapply(function(mat) if (nrow(mat) > 0) mat[nrow(mat), "start"] else NA_integer_)
  
  # Para cada string, processa a remoção do ":" após o último "|"
  mapply(function(label, pipe_pos) {
    if (is.na(pipe_pos)) return(label)  # Nenhum "|"
    
    after_pipe <- str_sub(label, pipe_pos + 1)
    colon_index <- str_locate(after_pipe, ":")[, "start"]
    
    if (is.na(colon_index)) {
      return(label)  # Nenhum ":" após o último "|"
    } else {
      cut_pos <- pipe_pos + colon_index  # Posição do ":" (vamos cortar ANTES dele)
      return(str_sub(label, 1, cut_pos - 1))
    }
  }, x, last_pipe, USE.NAMES = FALSE)
}


# Search function for one gene
search_by_gene <- function(genes, species, retmax = 10000, db = "nuccore", delay = 0.5) {
  # Construct a single species term (all species combined with OR)
  species_term <- paste0("(", paste(species, collapse = "[Organism] OR "), "[Organism])")
  
  cat(
    "Retriving GenBank seach (db =", db, ") for:\n",
    "- Genes:", paste(genes, collapse = ", "), "\n",
    "- Species:", paste(species, collapse = ", "), "\n",
    "- Max records:", retmax, "\n\n",
    "Constructed search terms:\n", paste(species_term, collapse = "\n"), "\n\n"
  )
  
  
  # Map over genes (not species)
  results <- purrr::map(genes, ~ {
    term <- paste0(species_term, " AND ", .x, "[gene]")
    cat("Searching for gene:", .x, "\n")
    
    Sys.sleep(delay)  # Avoid NCBI rate limits
    tryCatch(
      entrez_search(db = db, term = term, use_history = TRUE, retmax = retmax),
      error = function(e) {
        base::message("Query failed for gene: ", .x, "\nError: ", e$message)
        NULL
      }
    )
  })
  
  # Name the results with the genes
  names(results) <- genes
  return(results)
}

# Process each gene's results
process_gene_results <- function(gene_result, db = "nuccore", rettype = "gb", retmode = "text", delay = 0.5) {
  tryCatch({
    
    # Fetch records (moved outside of map_dfr)
    cat(
      "Searching GenBank genes (db =", db, ") for:\n",
      "- Genes:", paste(names(gene_result), collapse = ", "), "\n",
      "- Species:", paste(paste(unique(str_extract_all(gene_result[[1]]$QueryTranslation, '(?<=\\\")[A-Z][a-z]+(?=\\\"\\[Organism\\])')[[1]], collapse = ", ")), collapse = ", "), "\n",
      "- Rettype:", rettype, "\n\n"
    )
    
    gene_df <- data.frame()
    
    for (i in 1:length(gene_result)){
      
      cat("Searching for gene:", names(gene_result)[i] , "\n")
      
      gb_text <- entrez_fetch(
        db = db,
        web_history = gene_result[[i]]$web_history,
        rettype = rettype,
        retmode = retmode
      )
      
      # Split records
      record <- unlist(strsplit(gb_text, "//\n"))
      
      # Extract fields function
      accession <- str_match(record, "ACCESSION\\s+(\\S+)")[,2]
      organism <- str_match(record, "ORGANISM\\s+([^\n]+)")[,2]
      title <- str_match(record, "TITLE\\s+([^\n]+)")[,2]
      geo_loc <- str_match(record, "/geo_loc_name=\"([^\"]+)\"")[,2]
      gene <- str_match(record, "/gene=\"([^\"]+)\"")[,2]
      voucher <- str_match(record, "/specimen_voucher=\"([^\"]+)\"")[,2]  
      
      ph <- data.frame(
        ORGANISM = organism,
        VOUCHER = voucher,
        GENE = gene,
        GEO_LOC = geo_loc,
        ACCESSION = accession,
        TITLE = title,
        stringsAsFactors = FALSE
      )
      
      if (i == 1){
        gene_df <- ph
      } else {
        gene_df <- rbind(gene_df, ph)
      }
    }
    
    # Process all records
    gene_df[] <- lapply(gene_df, function(x) if(is.character(x)) gsub("\n", "", x) else x)
    
    return(gene_df)
  }, error = function(e) {
    base::message("Error processing gene: ", e)
    return(NULL)
  })
}

# 1. Function to process each element in the list bstringset = a
process_sequences <- function(bstringset) {
  
  cat(
    "Aligning sequences for:\n",
    "- Genes:", paste(str_extract_all(names(bstringset), "(?<=_)[^.]+(?=\\.)"), collapse = ", ")[1], "\n",
    "- Vouchers:", paste(unique(unlist(lapply(bstringset, function(x) {as.data.frame(x@ranges)[,"names"]})))[unique(unlist(lapply(bstringset, function(x) {as.data.frame(x@ranges)[,"names"]}))) %in% "1" == FALSE], collapse = ", "), "\n\n"
  )
  
  ph <- list()
  
  for (i in seq_along(bstringset)) {
    
    if (i == length(seq_along(bstringset))){
      is_next <- ", this was the last gene! \n\n Finishing run!\n\n"
    } else {
      is_next <- ", continuing with next... \n\n"
    }
    
    cat("Aligning gene:", paste0(str_extract_all(names(bstringset), "(?<=_)[^.]+(?=\\.)")[i]), "\n")
    current_gene <- names(bstringset)[i]
    current_seqs <- bstringset[[i]]
    
    # Skip if only contains NA (width=2)
    if(all(BiocGenerics::width(current_seqs) == 2 & grepl("NA", as.character(current_seqs)))) {
      cat("NA only found for gene ", paste0(str_extract_all(names(bstringset), "(?<=_)[^.]+(?=\\.)")[i]), is_next)
      ph[[current_gene]] <- NULL
      next
    }
    
    # Remove invalid sequences
    valid_seqs <- current_seqs[BiocGenerics::width(current_seqs) > 2 & !is.na(current_seqs)]
    
    # Skip if no valid sequences left
    if(length(valid_seqs) == 0) {
      cat("No valid sequences found for gene ", paste0(str_extract_all(names(bstringset), "(?<=_)[^.]+(?=\\.)")[i]), is_next)
      ph[[current_gene]] <- NULL
      next
    }
    
    # Clean sequences
    clean_seqs <- gsub(" ", "", valid_seqs)  # Remove spaces
    clean_seqs <- toupper(clean_seqs)       # Convert to uppercase
    clean_seqs <- gsub("[^ACGT]", "N", clean_seqs)  # Replace non-standard bases
    
    # Convert to DNAStringSet
    dna_seqs <- tryCatch(
      DNAStringSet(clean_seqs),
      error = function(e) {
        warning("Failed to convert sequences for ", current_gene, ": ", e$message)
        return(NULL)
      }
    )
    
    # Handle single sequence case
    if(length(dna_seqs) == 1) {
      warning("Only one valid sequence found for gene ", current_gene)
      cat("Finished gene:", paste0(str_extract_all(names(bstringset), "(?<=_)[^.]+(?=\\.)")[i]), is_next)
      ph[[current_gene]] <- dna_seqs
      next
    }
    
    # Multiple sequences - perform alignment
    aligned <- tryCatch(
      msa(dna_seqs, type = "dna"),
      error = function(e) {
        warning("Alignment failed for ", current_gene, ": ", e$message)
        return(NULL)
      }
    )
    
    ph[[current_gene]] <- if(!is.null(aligned)) as(aligned, "DNAStringSet") else dna_seqs
    cat("Finished gene:", paste0(str_extract_all(names(bstringset), "(?<=_)[^.]+(?=\\.)")[i]) , is_next)
  }
  
  return(ph)
}

# 4. Merge sequences by name - corrected version
merge_alignments <- function(processed_list) {
  # Get all unique sequence names across all genes
  all_names <- unique(unlist(lapply(processed_list, names)))
  
  # Get the width of each gene alignment
  gene_widths <- sapply(processed_list, function(x) BiocGenerics::width(x)[1])
  
  # Initialize a list to store merged sequences
  merged_sequences <- list()
  
  for (name in all_names) {
    seq_parts <- character(length(processed_list))
    
    for (i in seq_along(processed_list)) {
      gene_alignment <- processed_list[[i]]
      if (name %in% names(gene_alignment)) {
        seq_parts[i] <- as.character(gene_alignment[name])
      } else {
        # Fill with gaps of the correct gene length
        seq_parts[i] <- paste(rep("-", gene_widths[i]), collapse = "")
      }
    }
    
    # Combine all gene sequences for this individual
    merged_sequences[[name]] <- paste(seq_parts, collapse = "")
  }
  
  # Convert to DNAStringSet and ensure equal lengths
  merged_alignment <- DNAStringSet(unlist(merged_sequences))
  return(merged_alignment)
}

# boundaries between delimitation groups, for the bars beside the tree
get_boundaries <- function(df, col_name) {
  df %>%
    dplyr::group_by(!!sym(col_name)) %>%
    dplyr::summarise(
      ymin = min(y) - 0.5, # -0.5 to cover the full tip height
      ymax = max(y) + 0.5,
      .groups = "drop"
    )
}