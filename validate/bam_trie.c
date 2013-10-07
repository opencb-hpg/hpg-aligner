#include "bam_trie.h"

//--------------------------------------------------------------------------------------

void generate_logs(const char *file, cp_trie *trie);
int add_region_to_trie(cp_trie *trie, bam1_t *bam_line, bam_header_t *header);

int dna_map_region_equal_margin_soft(dna_map_region_t *region1, dna_map_region_t* region2,
				     int margin);

int rna_map_region_equal_margin_soft(dna_map_region_t *region1, dna_map_region_t* region2,
				     int margin);

//--------------------------------------------------------------------------------------

char *create_sequence_string(bam1_t *bam1) {
  char *bam_seq = bam1_seq(bam1);
  int seq_len = bam1->core.l_qseq;
  
  char *seq = (char *) malloc(seq_len * sizeof(char));

  // nucleotide content                                                                              
  for (int i = 0; i < seq_len; i++) {
    switch (bam1_seqi(bam_seq, i)) {
    case 1:
      seq[i] = 'A';
      break;
    case 2:
      seq[i] = 'C';
      break;
    case 4:
      seq[i] = 'G';
      break;
    case 8:
      seq[i] = 'T';
      break;
    case 15:
      seq[i] = 'N';
      printf("N");
      break;
    default:
      seq[i] = 'N';
      break;
    }
  }
  return seq;
}

//--------------------------------------------------------------------------------------

char *create_quality_string(bam1_t *bam1, int base_quality) {
  char *bam_qual = bam1_qual(bam1);
  int qual_len = bam1->core.l_qseq;

  char *qual = (char *) malloc(qual_len * sizeof(char));
  for (int i = 0; i < qual_len; i++) {
    qual[i] = base_quality + bam_qual[i];
  }

  return qual;
}



cp_trie *dna_dataset_to_trie(char * file, trie_result_t* result) {
  char buffer[MAX_DATASET_LINE_LENGTH];
  const char delimiters[] = "_";
  char* token = NULL;
  char* id="";
  char s;
  int pos;
  short int tam;
  FILE* fd = fopen(file,"r");
  
  // Create trie
  cp_trie *trie = cp_trie_create_trie(COLLECTION_MODE_DEEP, 0, (cp_destructor_fn) trie_node_free);

  dna_map_region_t *region;
  
  if(fd == NULL) { // Mejorar esta gestión de errores
    printf("Fallo al abrir el fichero");
    return NULL;
  }
  
  while(fgets(buffer,MAX_DATASET_LINE_LENGTH,fd))  {
    if (buffer[0] != '@') continue;
    
    id = strdup(&buffer[1]);
    id[strlen(id) - 1] = 0;

    // insert tot the list
    array_list_insert(id, id_list);

    if (strstr(id, "rand")) {
      region = NULL;
    } else {
      pos = 0;
      token = strtok(buffer, delimiters);
      region = (dna_map_region_t*) calloc(1, sizeof(dna_map_region_t)); 
      
      while (token != NULL) {
	switch (pos) {
	case 0:
	  region->chromosome = strdup(&token[1]);
	  break;
	case 1:
	  region->start_position = atoi(token);
	  break;
	case 2:
	  region->end_position = atoi(token);
	  break;
	case 3:
	  region->strand = atoi(token);
	  break;
	}
	token = strtok(NULL, delimiters);  
	pos++;  
	
	if (pos > 3) break;
      }
    }
    
    //    printf("id: %s\t\t", id); dna_print_region(region); printf("\n");
   
    cp_trie_add(trie, id, trie_node_new((void *)region));      
  }
  fclose(fd);
  return trie;
}

//--------------------------------------------------------------------------------------

cp_trie *rna_dataset_to_trie(char * file, trie_result_t* result) {
  char buffer[MAX_DATASET_LINE_LENGTH];
  const char delimiters[] = "@";
  char* token = NULL;
  char* id="";
  char s;
  int pos;
  short int tam;
  FILE* fd = fopen(file,"r");
  
  // Create trie
  cp_trie *trie = cp_trie_create_trie(COLLECTION_MODE_DEEP, 0, (cp_destructor_fn) trie_node_free);

  rna_map_region_t *region;
  
  if(fd == NULL) { // Mejorar esta gestión de errores
    printf("Fallo al abrir el fichero");
    return NULL;
  }
  
  while(fgets(buffer,MAX_DATASET_LINE_LENGTH,fd))  {
    if (buffer[0] != '@') continue;
    
    id = strdup(&buffer[1]);
    //printf("1.INSERT TO TRIE (%lu): %s\n",  strlen(id), id);
    id[strlen(id) - 3] = 0;    
    // insert tot the list
    //printf("2.INSERT TO TRIE: %s\n", id);
    array_list_insert(id, id_list);
    
    if (strstr(id, "rand")) {
      region = NULL;
    } else {
      pos = 0;
      token = strtok(buffer, delimiters);
      region = (rna_map_region_t*) calloc(1, sizeof(dna_map_region_t)); 
      
      while (token != NULL) {
	switch (pos) {
	case 3:
	  region->chromosome = strdup(&token[0]);
	  break;
	case 4:
	  region->start_position = atoi(token);
	  break;
	case 5:
	  region->end_position = atoi(token);
	  break;
	}
	token = strtok(NULL, delimiters);  
	pos++;  
	
	if (pos > 5) break;
      }
    }    
    //    printf("id: %s\t\t", id); dna_print_region(region); printf("\n");   
    cp_trie_add(trie, id, trie_node_new((void *)region));      
  }

  fclose(fd);

  return trie;
}

//--------------------------------------------------------------------------------------
exon_coords_t *new_exon_coords(int start, int end) {
  exon_coords_t *p = (exon_coords_t *) malloc(sizeof(exon_coords_t));

  p->start = start;
  p->end = end;

  return p;
}

void exon_coords_free(exon_coords_t *p) {
  if (p) { free(p); }
}

exon_data_t *exon_data_new(int strand, 
			   int start, int end,
			   char *transcript_id) {
  exon_data_t *p = (exon_data_t *) malloc(sizeof(exon_data_t));

  p->strand = strand;
  p->start = start;
  p->end = end;
  p->transcript_id = transcript_id;

  return p;
}

void exon_data_free(exon_data_t *p) {
  if (p) {
    if (p->transcript_id) free(p->transcript_id);
    free(p);
  }
}

//--------------------------------------------------------------------------------------

static inline char *parse_attribute(char *name, char *attrs) {
  int name_len = strlen(name);
  char *p1 = strstr(attrs, name);
  char * p2 = strstr(p1 + name_len, "\";");
  int len = p2 - p1 - name_len;
  char *res = (char *) malloc(len + 2);
  strncpy(res, p1 + name_len, len);
  res[len] = 0;
  return res;
}

//--------------------------------------------------------------------------------------

static inline exon_data_t *parse_exon_line(FILE *f) {
  const int MAX_LENGTH = 8192;

  int field, found;
  char line[MAX_LENGTH], *token, *str, *p1, *p2;
  
  int chr, strand, start, end, exon_number;
  char *chr_name, *gene_id, *transcript_id, *exon_id, *exon_number_str;

  while (fgets(line, MAX_LENGTH, f) != NULL) {
    //printf("PARSE: %s", line);
    str = strdup(line);
    token = strtok(str, "\t");
    field = 0;
    found = 1;
    while (token != NULL) {
      if (field == 2) { // feature type name: exon
	if (strcmp(token, "exon") != 0) {
	  found = 0;
	  break;
	}
      } else if (field == 3) { // start position of the feature (starting at 1)
	start = atoi(token);
      } else if (field == 4) { // end position of the feature (starting at 1)  
        end = atoi(token);
      } else if (field == 6) { // strand + (forward) or - (reverse)
	strand = ((strcmp(token, "-") == 0) ? 1 : 0);
      } else if (field == 8) { // attributes, a semicolon-separated list of tag-value pairs
	// gene_id, transcript_id, exon_id
	gene_id = parse_attribute("gene_id \"", token);
	transcript_id = parse_attribute("transcript_id \"", token);
	//exon_id = parse_attribute("exon_id \"", token);	
	free(gene_id);
      }      
      token = strtok(NULL, "\t");
      field++;
    }

    free(str);

    if (found) {
      return exon_data_new(strand, start, end, transcript_id);
    }

  }

  return NULL;

}

//--------------------------------------------------------------------------------------


cp_hashtable *load_transcriptome_validate(char *file) {
  FILE* f = fopen(file, "r");

  int pos, direction, count = 0;
  exon_data_t *exon = NULL, *exon1 = NULL, *exon2 = NULL;

  size_t g_start, g_end;
  int type, splice_strand, strand;
  char nt_start[2], nt_end[2];
  char *transcript_id;

  array_list_t *list = array_list_new(100, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *list_aux = array_list_new(100, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *list_new = array_list_new(100,
					  1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  
  cp_hashtable *t = cp_hashtable_create_by_option(COLLECTION_MODE_NOSYNC |
						  COLLECTION_MODE_COPY |
						  COLLECTION_MODE_DEEP,
						  10000,
						  cp_hash_istring,
						  (cp_compare_fn) strcasecmp,
						  (cp_copy_fn) strdup,
						  (cp_destructor_fn) free,
						  (cp_copy_fn) array_list_dup,
						  (cp_destructor_fn) array_list_free);
  while (1) {
    // read the first exon
    if (!exon) {
      exon = parse_exon_line(f);  
      if (!exon) {
	break;
      }
    }
    
    strand = exon->strand;
    array_list_insert(exon, list);
    transcript_id = exon->transcript_id;

    while ((exon = parse_exon_line(f)) && strcmp(exon->transcript_id, transcript_id) == 0) {
      array_list_insert(exon, list);
    }
    
    // process the whole list
    if (strand) {
      // - strand
      for (int i = array_list_size(list) - 1; i >= 0 ; i--) {
	exon1 = array_list_get(i, list);
	array_list_insert(new_exon_coords(exon1->start, exon1->end), list_new);
      }      
    } else {
      // + strand
      for (int i = 0; i < array_list_size(list); i++) {
	exon1 = array_list_get(i, list);
	//printf("\t EXON%i: %i-%i (%s)\n", i, exon1->start, exon1->end, exon1->transcript_id);
	array_list_insert(new_exon_coords(exon1->start, exon1->end), list_new);
      }
    }
    
    cp_hashtable_put(t, transcript_id, list_new);

    array_list_clear(list_new, NULL);    
    array_list_clear(list, exon_data_free);

  } // end while

  fclose(f);
  array_list_free(list, NULL);

  //list = cp_hashtable_get(t, "ENST00000606142");

  /*
    if (array_list_size(list) > 500 || 
      array_list_size(list) < 0) { 
    printf("1.ERROR OCURRED %i | %lu\n", array_list_size(list), array_list_size(list)); 
    exit(-1);
  } else {
    printf("1.OK %lu\n", array_list_size(list));
    exit(-1);
    }
  */

  return t;
}

//--------------------------------------------------------------------------------------

int limit_wrong_not_mapped = 0;
int limit_wrong_mapped = 0;

void dna_intersection(cp_trie *trie, int margin, char *filename, trie_result_t *result) {
  int read_bytes; 
  
  char *id, *prev_id = NULL;
  trie_node_t *node, *prev_node = NULL;

  dna_map_region_t region;
  dna_map_region_t* region_trie = NULL;

  result->margin = margin / 2;

  bam1_t* bam_line = bam_init1();

  // Open bam file for read
  bam_file_t* bam_file_p =  bam_fopen(filename);
  
  //header for BAM file has been done in the opening
  bam_header_t* bam_header_p = bam_file_p->bam_header_p;

  while ((read_bytes = bam_read1(bam_file_p->bam_fd, bam_line)) > 0) {
    id = strdup(bam1_qname(bam_line));

    if ((node = cp_trie_exact_match(trie, bam1_qname(bam_line))) == NULL)  {
      printf("id %s not found !!\n", id);
    }

    if (bam_line->core.flag == 4) {
      // unmapped read
      node->not_mapped++;      
      if (strstr(id, "rand")) {
	node->right_not_mapped++;
      } else {
	node->wrong_not_mapped++; 
      }
    } else {
      // mapped read
      node->mapped++;      
	
      // set region
      region_trie = (dna_map_region_t *) node->info;
      region.chromosome = strdup(bam_header_p->target_name[bam_line->core.tid]);;
      region.start_position = bam_line->core.pos + 1;
      region.end_position = region.start_position + bam_line->core.l_qseq + 1; // TODO mirar si usar el CIGAR
      region.strand = (((uint32_t) bam_line->core.flag) & BAM_FREVERSE) ? 1 : 0;
      
      // check the position
      if (dna_map_region_equal_margin_soft(region_trie, &region, margin)) {
	node->right_mapped++;
      } else {
	char *seq = create_sequence_string(bam_line);
	node->wrong_mapped++;
	if (strcmp(region.chromosome, region_trie->chromosome) != 0) {
	  node->wrong_mapped_chromosome++;
	} else	if (region.strand != region_trie->strand) {
	  node->wrong_mapped_strand++;
	}
	free(seq);
      }
      free(region.chromosome);
    }
    free(id);
  }
  
  //  result->filename= strdup(filename);
  
  // free memory and close file
  bam_destroy1(bam_line);
  bam_fclose(bam_file_p);
}


void rna_intersection(cp_trie *trie, int margin, char *filename, 
		      trie_result_t *result, cp_hashtable *trans) {
  int read_bytes; 
  
  char id[5000], *prev_id = NULL;
  trie_node_t *node, *prev_node = NULL;

  rna_map_region_t region;
  rna_map_region_t* region_trie = NULL;

  result->margin = margin / 2;

  bam1_t* bam_line = bam_init1();

  // Open bam file for read
  bam_file_t* bam_file_p =  bam_fopen(filename);
  
  //header for BAM file has been done in the opening
  bam_header_t* bam_header_p = bam_file_p->bam_header_p;
  char trans_id[1024];
  char *id_aux;
  int num_lines;

  while ((read_bytes = bam_read1(bam_file_p->bam_fd, bam_line)) > 0) {
    //printf("Lines %i\n", num_lines);
    char *id_aux2 = bam1_qname(bam_line);
    id_aux = strdup(id_aux2);
    num_lines++;
    //printf("--> %s\n", id_aux);
    int c = 0;
    while (c < strlen(id_aux) && id_aux[c] != '/') {
      id[c] = id_aux[c++];
    }
    id[c] = '\0';

    if ((node = cp_trie_exact_match(trie, id/*bam1_qname(bam_line)*/)) == NULL)  {
      printf("id %s not found !!\n", id);
      exit(-1);
    }

    if (bam_line->core.flag == 4) {
      // unmapped read
      node->not_mapped++;      
      if (strstr(id, "rand")) {
	node->right_not_mapped++;
      } else {
	node->wrong_not_mapped++; 
      }
    } else {
      // mapped read
      node->mapped++;      
	
      // set region
      region_trie = (rna_map_region_t *) node->info;
      char *chr_aux = bam_header_p->target_name[bam_line->core.tid];
      if (chr_aux[0] == 'c') {
	//PATCH for mapSplice v2
	region.chromosome = strdup(&chr_aux[3]);
      } else {
	region.chromosome = strdup(bam_header_p->target_name[bam_line->core.tid]);;
      }
      //printf("CHROMOSOME :%s\n", region.chromosome);
      //exit(-1);      

      region.start_position = bam_line->core.pos + 1;
      region.end_position = region.start_position + bam_line->core.l_qseq + 1;
      region.strand = (((uint32_t) bam_line->core.flag) & BAM_FREVERSE) ? 1 : 0;
      
      // check the position
      if (rna_map_region_equal_margin_soft(region_trie, &region, margin)) {
	node->right_mapped++;
	uint32_t k;
	int c = 0;
	uint32_t *cigar_bam = bam1_cigar(bam_line);
	char *cigar = convert_to_cigar_string(cigar_bam, bam_line->core.n_cigar);
	
	while (id[c] != '@') { trans_id[c] = id[c++]; }
	trans_id[c] = '\0';
	
	//printf("Search trans_id: %s\n", trans_id);
	array_list_t *list = cp_hashtable_get(trans, trans_id);
	
	//if (array_list_size(list) > 100) { printf("1.ERROR OCURRED %i | %lu\n", array_list_size(list), array_list_size(list)); exit(-1); }
	//printf("READ %s :\n", bam_line->data);
	//printf("\t %s : CIGAR(%i): %s\n", trans_id, bam_line->core.n_cigar, cigar);

	//Search exon location
	int found = 0;
	cigar_code_t *cigar_code = cigar_code_new_by_string(cigar);
	int genome_start = region.start_position;
	int genome_end   = region.start_position;

	//printf("First OPs\n");
	for (int i = 0; i < cigar_code->ops->size; i++) {
	  cigar_op_t *op = array_list_get(i, cigar_code->ops);
	  if (op->name == 'M' || op->name == 'D' || 
	      op->name == 'N' || op->name == 'P' ||
	      op->name == '=' ) {
	    genome_end += op->number;
	  }
	}
	//printf("First OPs End %lu\n", array_list_size(list));
	int found_start = 0, found_end = 0;
	int i;

	//Search start position in exon
	for (i = 0; i < array_list_size(list); i++) {
	  exon_coords_t *coords = array_list_get(i, list); 
	  //printf("\tSTART: %i :[%lu-%lu]\n", genome_start, coords->start, coords->end);
	  if (genome_start >= coords->start && genome_start <= coords->end) {
	    //Found exon!
	    found_start = coords->start;
	    break;
	  }
	}

	//printf("Genome Start %i: %s\n", genome_start, found_start == 0 ? "NO" : "YES");

	for (; i < array_list_size(list); i++) {
	  exon_coords_t *coords = array_list_get(i, list); 
	  //printf("\tEND: %i :[%lu-%lu]\n", genome_end, coords->start, coords->end);
	  if (genome_end >= coords->start && genome_end <= coords->end) {
	    found_end = coords->end;
	    break;
	  }
	}

	//printf("Genome End %i: %s\n", genome_end, found_end == 0 ? "NO" : "YES");
	if (found_end && found_start) {
	  cigar_op_t *first_op = array_list_get(0, cigar_code->ops);
	  if (first_op->name == 'S' || 
	      first_op->name == 'H') {
	    genome_start -= first_op->number;
	    if (genome_start < found_start) {
	      found_start = 0;
	    }
	  }
	  
	  //printf("Cigar Precision START 'S' & 'H' : %s\n", found_start == 0 ? "NO" : "YES");

	  if (found_start) {
	    //Search exact start and end position
	    cigar_op_t *last_op = array_list_get(cigar_code->ops->size - 1, cigar_code->ops);	    
	    if (last_op->name == 'S' || 
		last_op->name == 'H') {
	      genome_end += last_op->number;
	      if (genome_end > found_end) {
		found_end = 0;
	      }
	    }
	  }

	  //printf("Cigar Precision END 'S' & 'H' : %s\n", found_end == 0 ? "NO" : "YES");

	}
	    
	if (!found_end || !found_start) {
	  node->wrong_sj++;
	  //printf("@@@@@ ERROR SPLICE : %s ****\n", bam_line->data);
	} else {
	  node->right_sj++; 
	  //printf("@@@@@ CORRECT SPLICE : %s ****\n", bam_line->data);
	}	

	free(cigar);
      } else {
	char *seq = create_sequence_string(bam_line);
	node->wrong_mapped++;
	if (strcmp(region.chromosome, region_trie->chromosome) != 0) {
	  node->wrong_mapped_chromosome++;
	}
	free(seq);
      }
      free(region.chromosome);
    }
    free(id_aux);
  }
  
  // free memory and close file
  //bam_destroy1(bam_line);
  //bam_fclose(bam_file_p);
}


void print_result(trie_result_t *result, int log) {

  int is_rand, num_reads = array_list_size(id_list);

  char *id;
  trie_node_t *node;
  for (int i = 0; i < num_reads; i++) {
    is_rand = 0;
    id = array_list_get(i, id_list);
    node = cp_trie_exact_match(trie, id);
    if (strstr(id, "rand")) {
      is_rand = 1;
      result->rand_reads++;
    }

    if (node->mapped) {
      result->mapped++;
      
      if (node->right_mapped) {
	result->right_mapped++;
	if (node->right_sj) {
	  result->right_sj++;
	} else {
	  result->wrong_sj++;
	}

      } else {
	result->wrong_mapped++;
	if (node->wrong_mapped_chromosome) {
	  if (++limit_wrong_mapped < 10) {
	    //printf("\twrong mapped (chromosome): %s\n", id);
	  }
	  result->wrong_mapped_chromosome++;
	} else if (node->wrong_mapped_strand) {
	  if (++limit_wrong_mapped < 10) {
	    //printf("\twrong mapped (strand): %s\n", id);
	  } 
	  result->wrong_mapped_strand++;
	}
      }
      result->multi_right_mapped += node->right_mapped;
      result->multi_wrong_mapped += node->wrong_mapped;

      result->multi_mapped += node->mapped;
    } else {
      result->not_mapped++;
      if (is_rand) {
	result->right_not_mapped++;
      } else {
	result->wrong_not_mapped++;
	if (++limit_wrong_not_mapped < 5) {
	  //printf("\twrong not mapped: %s\n", id);
	}
      }
    }
  }


  int total = result->mapped + result->not_mapped;

  printf("\nU N I Q U E   A L I G N M E N T\n");
  printf("-------------------------------------\n");
  printf("\tMargin length: %d\n", (result->margin * 2));
  printf("\tNum. reads   : %d\n", num_reads);
  printf("\tNum. mappings: %d\n\n", total);
  printf("\tMapped: %d (%0.2f %%)\n",result->mapped, 100.0f * result->mapped / total);
  printf("\t\tRight mapped: %d (%0.2f %%) -> removing random reads: %0.2f %%\n", 
	 result->right_mapped, 100.0f * result->right_mapped / total, 
	 100.0f * result->right_mapped / (total - result->rand_reads));
  printf("\t\tWrong mapped: %d (%0.2f %%): chromosome mismatch: %d (%0.2f %%), strand mismatch: %d (%0.2f %%)\n", 
	 result->wrong_mapped, 100.0f * result->wrong_mapped / total,
	 result->wrong_mapped_chromosome, 100.0f * result->wrong_mapped_chromosome / total,
	 result->wrong_mapped_strand, 100.0f * result->wrong_mapped_strand / total);
  printf("\n");
  printf("\tNot mapped: %d (%0.2f %%)\n", result->not_mapped, 100.0f * result->not_mapped / total);
  printf("\t\tRight not mapped: %d (%0.2f %%) @rand reads: %i (%0.2f %%)\n", 
	 result->right_not_mapped, 100.0f * result->right_not_mapped / total, 
	 result->rand_reads, 100.0f * result->right_not_mapped / result->rand_reads);
  printf("\t\tWrong not mapped: %d (%0.2f %%)\n", result->wrong_not_mapped, 100.0f * result->wrong_not_mapped / total);

  total = result->multi_mapped + result->not_mapped;
  printf("\nM U L T I   -   A L I G N M E N T S\n");
  printf("-------------------------------------\n");
  printf("\tMargin length      : %d\n", (result->margin * 2));
  printf("\tNum. reads         : %d\n", num_reads);
  printf("\tNum. multi mappings: %d\n\n", total);
  printf("\tMapped: %d (%0.2f %%)\n", result->multi_mapped, 100.0f * result->multi_mapped / total);
  printf("\t\tRight mapped: %d (%0.2f %%)\n", result->multi_right_mapped, 100.0f * result->multi_right_mapped / total);
  printf("\t\tWrong mapped: %d (%0.2f %%)\n", result->multi_wrong_mapped, 100.0f * result->multi_wrong_mapped / total);
  printf("\n");

  printf("\nS P L I C E    -    S E N S I T I V I T Y\n");
  printf("-------------------------------------------\n");
  printf("\Total Reads Right Mapped : %d (%0.2f)\n", result->right_mapped, (100.0f * result->right_mapped) / num_reads);
  printf("\tRight SJ alignments: %d (%0.2f %%)\n", result->right_sj, (100.0f * result->right_sj) / result->right_mapped);
  printf("\tWrong SJ alignments: %d (%0.2f %%)\n", result->wrong_sj, (100.0f * result->wrong_sj) / result->right_mapped);

  int total_right = result->right_mapped - result->wrong_sj;

  printf("\tTOTAL: %d/%d (%0.2f %%)\n", total_right, num_reads, (100.0f * total_right) / num_reads);
}

//--------------------------------------------------------------------------------------

int dna_map_region_equal_margin_soft(dna_map_region_t *region1, dna_map_region_t* region2,
				     int margin) {
  int result = 0;

  if (region1->strand != region2->strand) {
    if (strcmp(region1->chromosome, region2->chromosome) == 0) { 
      result = ((region1->start_position - margin <= region2->end_position && 
		 region1->start_position + margin >= region2->end_position)
		||
		(region1->end_position - margin <= region2->start_position && 
		 region1->end_position + margin >= region2->start_position));

      //      result = ((region1->start_position == region2->end_position) || 
      //		(region2->start_position == region1->end_position));
    }
  } else {
    if (strcmp(region1->chromosome, region2->chromosome) == 0) { 
      result = (region1->start_position - margin <= region2->start_position && 
		region1->start_position + margin >= region2->start_position);
      //      result = (region1->start_position == region2->start_position);
      
    }
  }
  
  return result;
}

int rna_map_region_equal_margin_soft(dna_map_region_t *region1, dna_map_region_t* region2,
				     int margin) {
  int result = 0;
  
  if (strcmp(region1->chromosome, region2->chromosome) == 0) { 
    result = (region1->start_position - margin <= region2->start_position && 
	      region1->end_position + margin >= region2->end_position);    
    
  }
  
  return result;
}

//--------------------------------------------------------------------------------------

cp_trie* dna_bam_to_trie(char * file) {
  // Create trie
  cp_trie *trie = cp_trie_create(0);
  
  int read_bytes;  
  //		char* bam_string;
  int cont=0;
  bam1_t* bam_p = bam_init1();
  
  // Open bam file for read
  bam_file_t* bam_file_p =  bam_fopen(file);
  
  //header for BAM file has been done in the opening
  bam_header_t* bam_header_p = bam_file_p->bam_header_p;
  
  while ((read_bytes = bam_read1(bam_file_p->bam_fd, bam_p)) > 0) {
    add_region_to_trie(trie,bam_p,bam_header_p);
    cont++;
  }
  
  // Free Memory
  bam_fclose(bam_file_p);
  bam_destroy1(bam_p);
  return trie; 
}

//--------------------------------------------------------------------------------------

int add_region_to_trie(cp_trie *trie, bam1_t *bam_line, bam_header_t *header) {  
  // Create new region
  dna_map_region_t * region = (dna_map_region_t *) calloc(1,sizeof(dna_map_region_t));
  
  // Add chromosome info
  int tid = bam_line->core.tid;
  
  if (tid == -1) {
    // region->chromosome = (char*)calloc(2,sizeof(char));
    // strcpy(region->chromosome,"-");
    region->chromosome = "-";
  } else {
    region->chromosome = (char*) calloc(strlen(header->target_name[tid])+1,sizeof(char));
    strcpy(region->chromosome,header->target_name[tid] );
  }
  
  // Add start position info
  region->start_position = bam_line->core.pos + 1; // TODO mirar el FLAG si es 1-base o 0-base
  
  // Add end position info
  region->end_position = region->start_position + bam_line->core.l_qseq+1; // TODO mirar si usamos el CIGAR
  
  // Add strand info
  uint32_t flag = (uint32_t) bam_line->core.flag;
  region->strand = (flag & BAM_FREVERSE) ? 1 : 0;
  
  // Get the Seq name Id
  char * id= (char*) calloc(bam_line->core.l_qname,sizeof(char));
  strcpy(id,bam1_qname(bam_line));
  
  cp_trie_add(trie, id, region);
  //		free(id);
  return 1;
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
