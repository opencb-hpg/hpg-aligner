#include "bam_trie.h"

//--------------------------------------------------------------------------------------

void generate_logs(const char *file, cp_trie *trie);
int add_region_to_trie(cp_trie *trie, bam1_t *bam_line, bam_header_t *header);

int dna_map_region_equal_margin_soft(dna_map_region_t *region1, dna_map_region_t* region2,
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
      } else {
	result->wrong_mapped++;
	if (node->wrong_mapped_chromosome) {
	  if (++limit_wrong_mapped < 10) {
	    printf("\twrong mapped (chromosome): %s\n", id);
	  }
	  result->wrong_mapped_chromosome++;
	} else if (node->wrong_mapped_strand) {
	  if (++limit_wrong_mapped < 10) {
	    printf("\twrong mapped (strand): %s\n", id);
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
	  printf("\twrong not mapped: %s\n", id);
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
  printf("\tNum. multi mappings: %d\n\n",total);
  printf("\tMapped: %d (%0.2f %%)\n", result->multi_mapped, 100.0f * result->multi_mapped / total);
  printf("\t\tRight mapped: %d (%0.2f %%)\n", result->multi_right_mapped, 100.0f * result->multi_right_mapped / total);
  printf("\t\tWrong mapped: %d (%0.2f %%)\n", result->multi_wrong_mapped, 100.0f * result->multi_wrong_mapped / total);
  printf("\n");

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
