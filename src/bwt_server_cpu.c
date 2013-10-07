#include "bwt_server.h"

#define SIMPLE_SW 1
#define SP_SW 2

#define SW_NORMAL 0 
#define SW_FINAL 1 

//====================================================================================
// bwt_server_input functions: init
//====================================================================================

void bwt_server_input_init(list_t* read_list_p, unsigned int batch_size, bwt_optarg_t *bwt_optarg_p, 
			   bwt_index_t *bwt_index_p, list_t* write_list_p, unsigned int write_size, 
			   list_t* unmapped_read_list_p, metaexons_t *metaexons, sw_optarg_t *sw_optarg,
			   genome_t *genome, bwt_server_input_t* input_p) {
  input_p->read_list_p = read_list_p;
  input_p->batch_size = batch_size;
  input_p->bwt_optarg_p = bwt_optarg_p;
  input_p->write_list_p = write_list_p;
  input_p->write_size = write_size;
  input_p->bwt_index_p = bwt_index_p;
  input_p->unmapped_read_list_p = unmapped_read_list_p;
  input_p->metaexons = metaexons;
  input_p->sw_optarg = sw_optarg;
  input_p->genome = genome;
}

//====================================================================================
// apply_bwt
//====================================================================================

cal_t *convert_bwt_anchor_to_CAL(bwt_anchor_t *bwt_anchor, size_t read_start, size_t read_end) {
  linked_list_t *linked_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
  seed_region_t *seed_region = seed_region_new(read_start, read_end,
					       bwt_anchor->start, bwt_anchor->end, 0);

  linked_list_insert_first(seed_region, linked_list);

  cal_t *cal = cal_new(bwt_anchor->chromosome + 1, bwt_anchor->strand,
		       bwt_anchor->start, bwt_anchor->end,
		       1, linked_list,
		       linked_list_new(COLLECTION_MODE_ASYNCHRONIZED));

  return cal;

}

size_t bwt_search_pair_anchors(array_list_t *list, unsigned int read_length) {
  bwt_anchor_t *bwt_anchor;
  int max_double_anchor = 0, max_anchor_length = 0;
  
  bwt_anchor_t *max_anchor = NULL;
  bwt_anchor_t *bwt_anchor_back, *bwt_anchor_forw;
  int anchor_length_tmp, anchor_back, anchor_forw;
  int strand, type;
  int found_anchor = 0, found_double_anchor = 0;

  const int MIN_ANCHOR = 25;
  const int MIN_SINGLE_ANCHOR = 40;

  //const int MIN_DOUBLE_ANCHOR = MIN_ANCHOR*2;
  const int MAX_BWT_REGIONS = 50;
  const int MAX_BWT_ANCHOR_DISTANCE = 500000;

  array_list_t *anchor_list_tmp, *forward_anchor_list, *backward_anchor_list;
  cal_t *cal;
  int seed_size, gap_read, gap_genome;

  array_list_t *backward_anchor_list_0 = array_list_new(MAX_BWT_REGIONS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *forward_anchor_list_0 = array_list_new(MAX_BWT_REGIONS, 1.25f , COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *backward_anchor_list_1 = array_list_new(MAX_BWT_REGIONS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *forward_anchor_list_1 = array_list_new(MAX_BWT_REGIONS, 1.25f , COLLECTION_MODE_ASYNCHRONIZED);

  array_list_t *big_anchor_list = array_list_new(MAX_BWT_REGIONS, 1.25f , COLLECTION_MODE_ASYNCHRONIZED);

  //printf("Tot Anchors %i\n", array_list_size(list));
  for (int i = 0; i < array_list_size(list); i++) {
    bwt_anchor = array_list_get(i, list);
    if (bwt_anchor->strand == 1) {
      //printf("(-)bwt anchor %i:%lu-%lu (%i): \n", bwt_anchor->chromosome + 1, bwt_anchor->start, bwt_anchor->end, bwt_anchor->end - bwt_anchor->start + 1);
      if (bwt_anchor->type == FORWARD_ANCHOR) {
	array_list_insert(bwt_anchor, forward_anchor_list_1);
	//printf("FORW\n");
      } else {
	array_list_insert(bwt_anchor, backward_anchor_list_1);
	//printf("BACK\n");
      }
    } else {
      //printf("(+)bwt anchor %i:%lu-%lu (%i): \n", bwt_anchor->chromosome + 1, bwt_anchor->start, bwt_anchor->end, bwt_anchor->end - bwt_anchor->start + 1);
      if (bwt_anchor->type == FORWARD_ANCHOR) {
	array_list_insert(bwt_anchor, forward_anchor_list_0);
	//printf("FORW\n");
      } else {
	array_list_insert(bwt_anchor, backward_anchor_list_0);
	//printf("BACK\n");
      }
    }

    anchor_length_tmp = bwt_anchor->end - bwt_anchor->start + 1;
    if (anchor_length_tmp > MIN_SINGLE_ANCHOR && anchor_length_tmp > max_anchor_length) {
      max_anchor_length = anchor_length_tmp;
      found_anchor = 1;
      strand = bwt_anchor->strand;
      type = bwt_anchor->type;
    }
    
    if (read_length - anchor_length_tmp < 16) {
      array_list_insert(bwt_anchor, big_anchor_list);
    } 
    
  }
  
  array_list_clear(list, NULL);

  if (array_list_size(big_anchor_list) > 0) {
    for (int i = array_list_size(big_anchor_list) - 1; i >= 0; i--) {
      //printf("Insert cal %i\n", i);
      bwt_anchor = array_list_remove_at(i, big_anchor_list);
      size_t seed_size = bwt_anchor->end - bwt_anchor->start;

      if (bwt_anchor->type == FORWARD_ANCHOR) {
	cal = convert_bwt_anchor_to_CAL(bwt_anchor, 0, seed_size);
      } else {
	cal = convert_bwt_anchor_to_CAL(bwt_anchor, read_length - seed_size - 1, read_length - 1);
      }
      
      array_list_insert(cal, list);
    }
    array_list_set_flag(SINGLE_ANCHORS, list);
    
    goto exit;
  }

  for (int type = 1; type >= 0; type--) {
    if (!type) {
      forward_anchor_list = forward_anchor_list_1;
      backward_anchor_list = backward_anchor_list_1;
      //printf("Strand (+): %i-%i\n", array_list_size(forward_anchor_list), array_list_size(backward_anchor_list));
    } else { 
      forward_anchor_list = forward_anchor_list_0;
      backward_anchor_list = backward_anchor_list_0;
      //printf("Strand (-): %i-%i\n", array_list_size(forward_anchor_list), array_list_size(backward_anchor_list));
    }

    //Associate Anchors (+)/(-)
    for (int i = 0; i < array_list_size(forward_anchor_list); i++) { 
      bwt_anchor_forw = array_list_get(i, forward_anchor_list);
      for (int j = 0; j < array_list_size(backward_anchor_list); j++) { 
	bwt_anchor_back = array_list_get(j, backward_anchor_list);
	anchor_forw = (bwt_anchor_forw->end - bwt_anchor_forw->start + 1);
	anchor_back = (bwt_anchor_back->end - bwt_anchor_back->start + 1); 

	anchor_length_tmp = anchor_forw + anchor_back;

	//printf("\tCommpare %i:%lu-%lu with %i:%lu-%lu\n", bwt_anchor_forw->chromosome + 1, 
	//     bwt_anchor_forw->start, bwt_anchor_forw->end, bwt_anchor_back->chromosome + 1, 
	//     bwt_anchor_back->start, bwt_anchor_back->end);
	if (bwt_anchor_forw->chromosome == bwt_anchor_back->chromosome &&
	    abs(bwt_anchor_back->start - bwt_anchor_forw->end) <= MAX_BWT_ANCHOR_DISTANCE && 
	    anchor_forw >= MIN_ANCHOR && anchor_back >= MIN_ANCHOR) {
	  
	  if (bwt_anchor_back->start < bwt_anchor_forw->end) { continue; }
	  
	  gap_read = read_length - (anchor_forw + anchor_back);
	  gap_genome = bwt_anchor_back->start - bwt_anchor_forw->end;

	  //printf("anchor_forw = %i, anchor_back = %i, gap_read = %i, gap_genome = %i\n",
	  //	 anchor_forw, anchor_back, gap_read, gap_genome);
	  	  
	  int apply_flank = 0;
	  if (gap_read < 2 || gap_genome < 2) {
	    int gap;
	    if (gap_read < 0 && gap_genome < 0) {
	      gap = abs(gap_read) > abs(gap_genome) ? abs(gap_read) : abs(gap_genome);
	    } else if (gap_read < 0) {
	      gap = abs(gap_read);
	    } else if (gap_genome < 0) {
	      gap = abs(gap_genome);
	    } else {
	      gap = 2;
	    }
	    
	    int flank  = 5;
	    apply_flank = 1;
	    
	    if (abs(gap) >= flank*2) {
	      //Solve read overlap
	      flank = abs(gap)/2 + flank/2;
	    }
	    //printf("\tgap = %i, flank = %i\n", gap, flank);
	    if (flank >= anchor_forw) {
	      bwt_anchor_forw->end -= anchor_forw/2;	      
	    } else {
	      bwt_anchor_forw->end -= flank;
	    }

	    if (flank >= anchor_back) {
	      bwt_anchor_back->start += anchor_back/2;	    
	    } else {
	      bwt_anchor_back->start += flank;
	    }
	  } 
	  	  
	  cal = convert_bwt_anchor_to_CAL(bwt_anchor_forw, 0, bwt_anchor_forw->end - bwt_anchor_forw->start);
	  //printf("INSERT-1 (%i)[%i:%lu-%lu]\n", cal->strand, cal->chromosome_id, cal->start, cal->end);
	  array_list_insert(cal, list);
	  seed_size = bwt_anchor_back->end - bwt_anchor_back->start + 1;
	  //if (bwt_anchor_forw->end + read_length >= bwt_anchor_back->start) {	    
	  //seed_region_t *seed_region = seed_region_new(read_length - seed_size, read_length - 1,
	  //bwt_anchor_back->start, bwt_anchor_back->end, 1);
	  //cal->end = bwt_anchor_back->end;
	  //linked_list_insert_last(seed_region, cal->sr_list);	
	  //} else {
	  cal = convert_bwt_anchor_to_CAL(bwt_anchor_back, read_length - seed_size, read_length - 1);
	  //printf("INSERT-2 (%i)[%i:%lu-%lu]\n", cal->strand, cal->chromosome_id, cal->start, cal->end);
	  array_list_insert(cal, list);
	  if (array_list_size(list) > 5) { goto exit; }
	  //}

	  array_list_set_flag(DOUBLE_ANCHORS, list);
	  found_double_anchor = 1;
	}                                                                                                                      
      }         
    }
  }

  if (!found_double_anchor && found_anchor) { 
    //Not Double anchor found but one Yes!!
    if (strand == 1) {
      if (type == FORWARD_ANCHOR) {
	anchor_list_tmp = forward_anchor_list_1;
      } else {
	anchor_list_tmp =  backward_anchor_list_1;
      }
    } else {
      if (type == FORWARD_ANCHOR) {
	anchor_list_tmp =  forward_anchor_list_0;
      } else {
	anchor_list_tmp =  backward_anchor_list_0;
      }
    }

    //printf("LIST SIZE %i\n", array_list_size(anchor_list_tmp));
    for (int i = 0; i < array_list_size(anchor_list_tmp); i++) {
      bwt_anchor = array_list_get(i, anchor_list_tmp);
      size_t seed_size = bwt_anchor->end - bwt_anchor->start;
      //array_list_insert(bwt_anchor_new(bwt_anchor->strand, bwt_anchor->chromosome, 
      //			       bwt_anchor->start, bwt_anchor->end, bwt_anchor->type), anchor_list);
      if (bwt_anchor->type == FORWARD_ANCHOR) {
	//printf("------------------------> start %i\n", 0);
	cal = convert_bwt_anchor_to_CAL(bwt_anchor, 0, seed_size);
      } else {
	//printf("------------------------> start %i\n", read_length - seed_size);
	cal = convert_bwt_anchor_to_CAL(bwt_anchor, read_length - seed_size - 1, read_length - 1);
      }
      array_list_insert(cal, list);
    }
    array_list_set_flag(SINGLE_ANCHORS, list);
  } 

 exit:
  array_list_free(forward_anchor_list_1, (void *)bwt_anchor_free);
  array_list_free(backward_anchor_list_1,  (void *)bwt_anchor_free);
  array_list_free(forward_anchor_list_0,  (void *)bwt_anchor_free);
  array_list_free(backward_anchor_list_0,  (void *)bwt_anchor_free);
  array_list_free(big_anchor_list,  (void *)bwt_anchor_free);

  return array_list_size(list);
  
}

//------------------------------------------------------------------------------------
// DNA
//------------------------------------------------------------------------------------

int apply_bwt(bwt_server_input_t* input, batch_t *batch) {
  //printf("APPLY BWT SERVER...\n");
  struct timeval start, end;
  double time;
  metaexons_t *metaexons = input->metaexons;
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  size_t num_reads = array_list_size(mapping_batch->fq_batch);
  size_t num_mappings;
  array_list_t *list;
  size_t *unmapped_indices = mapping_batch->targets;
  size_t num_unmapped = 0;
  size_t num_anchors;

  for (int i = 0; i < num_reads; i++) {
    fastq_read_t *read = array_list_get(i, mapping_batch->fq_batch);
    //printf("BWT: %s\n", read->id);
    list = mapping_batch->mapping_lists[i];    
    array_list_set_flag(0, list);
    num_mappings = bwt_map_inexact_read(read,
					input->bwt_optarg_p,
					input->bwt_index_p,
					list);
    if (array_list_get_flag(list) != 2) { //If flag 2, the read exceded the max number of mappings
      if (array_list_get_flag(list) == 1) {
	if (num_mappings > 0) {
	  num_anchors = bwt_search_pair_anchors(list, read->length);	
	  if (num_anchors == 0) {
	    array_list_set_flag(NOT_ANCHORS, list);
	  } 
	} else {
	  array_list_set_flag(NOT_ANCHORS, list);
	}
	//printf("tot anchors found %i %s\n", num_anchors, read->id);
	unmapped_indices[num_unmapped++] = i;
      } else if (num_mappings <= 0) {
	array_list_set_flag(0, list);
	//printf("Read NO Mapped %i %s\n", num_anchors, read->id);
	unmapped_indices[num_unmapped++] = i;
      }
    } else {
      array_list_set_flag(ALIGNMENTS_EXCEEDED, list);
    }
  }
  // array_list flag: 0 -> Not  BWT Anchors found
  //                  1 -> One  BWT Anchors found
  //                  2 -> Pair BWT Anchors found
  //                  3 -> Alignments found
  //                  4 -> Alignments exceded
  
  mapping_batch->num_targets = num_unmapped;
    
  if (batch->mapping_batch->num_targets > 0) {
    return SEEDING_STAGE;
  }

  if (batch->mapping_mode == RNA_MODE) {
    return RNA_POST_PAIR_STAGE;
  } else {
    return DNA_POST_PAIR_STAGE;
  }    
}

//------------------------------------------------------------------------------------
// RNA
//------------------------------------------------------------------------------------

int apply_bwt_rna(bwt_server_input_t* input, batch_t *batch) {

  LOG_DEBUG("========= APPLY BWT RNA =========\n");

  struct timeval start, end;
  double time;
  metaexons_t *metaexons = input->metaexons;
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  size_t num_reads = array_list_size(mapping_batch->fq_batch);
  size_t num_mappings;
  array_list_t *list;
  size_t *unmapped_indices = mapping_batch->targets;
  size_t num_unmapped = 0;
  size_t num_anchors;

  for (int i = 0; i < num_reads; i++) {
    fastq_read_t *read = array_list_get(i, mapping_batch->fq_batch);
    //printf("BWT: %s\n", read->id);
    list = mapping_batch->mapping_lists[i];    
    array_list_set_flag(1, list);
    num_mappings = bwt_map_inexact_read(read,
					input->bwt_optarg_p,
					input->bwt_index_p,
					list);
    if (array_list_get_flag(list) != 2) { //If flag 2, the read exceded the max number of mappings
      if (array_list_get_flag(list) == 1) {
	if (num_mappings > 0) {
	  num_anchors = bwt_search_pair_anchors(list, read->length);	
	  if (num_anchors == 0) {
	    array_list_set_flag(NOT_ANCHORS, list);
	    unmapped_indices[num_unmapped++] = i;
	  }
	} else {
	  array_list_set_flag(NOT_ANCHORS, list);
	  unmapped_indices[num_unmapped++] = i;
	} 
	//printf("tot anchors found %i %s\n", num_anchors, read->id);
      } else if (num_mappings <= 0) {
	array_list_set_flag(0, list);
	//printf("Read NO Mapped %i %s\n", num_anchors, read->id);
	unmapped_indices[num_unmapped++] = i;
      }else {
	//Read Map, Metaexon Actualization
	array_list_set_flag(ALIGNMENTS_FOUND, list);
	for (int i = 0; i < num_mappings; i++) {
	  alignment_t *alignment = array_list_get(i, list);
	  metaexon_insert(0/*alignment->seq_strand*/, alignment->chromosome,
			  alignment->position, alignment->position + read->length, 40,
			  METAEXON_NORMAL, NULL,
			  metaexons);
	}
      }
    } else {
      array_list_set_flag(ALIGNMENTS_EXCEEDED, list);
    }
  
    if (array_list_get_flag(list) == DOUBLE_ANCHORS) {
      //printf("DOUBLE ANCHORS\n");
      for (int j = 0; j < array_list_size(list); j++) {
	//bwt_anchor_t *bwt_anchor_prev = array_list_get(j, list);
	cal_t *cal = array_list_get(j, list);
	metaexon_insert(0/*cal->strand*/, cal->chromosome_id - 1,
			cal->start, cal->end, 40,
			METAEXON_NORMAL, NULL,
			metaexons);
      }    
    } else if (array_list_get_flag(list) == SINGLE_ANCHORS) {
      for (int j = 0; j < array_list_size(list); j++) {
	//bwt_anchor_t *bwt_anchor = array_list_get(j, list);
	cal_t *cal = array_list_get(j, list);
	metaexon_t *metaexon;
	if (metaexon_search(0/*cal->strand*/, cal->chromosome_id - 1,
			    cal->start, cal->end, &metaexon,
			    metaexons)) {
	  metaexon_insert(0/*cal->strand*/, cal->chromosome_id - 1,
			  cal->start, cal->end, 40,
			  METAEXON_NORMAL, NULL,
			  metaexons);
	}
      }
    }
  }
  // array_list flag: 0 -> Not  BWT Anchors found
  //                  1 -> One  BWT Anchors found
  //                  2 -> Pair BWT Anchors found
  //                  3 -> Alignments found
  //                  4 -> Alignments exceded
  
  mapping_batch->num_targets = num_unmapped;

  LOG_DEBUG("========= APPLY BWT RNA END =========\n");
    
  if (batch->mapping_batch->num_targets > 0) {
    return RNA_CAL_STAGE;
  } else {
    return RNA_STAGE;
  }
    
}

//------------------------------------------------------------------------------------
/*
int fill_anchor_gap(bwt_anchor_t *bwt_anchor_prev, bwt_anchor_t *bwt_anchor_next, 
		    genome_t *genome, fastq_read_t *read, char *sequence, 
		    const int max_err, int *l_extend, int *r_extend) {

  int anchor_len_prev = (bwt_anchor_prev->end - bwt_anchor_prev->start) + 1;
  int anchor_len_next = (bwt_anchor_next->end - bwt_anchor_next->start) + 1;
  size_t genome_start = bwt_anchor_prev->end + 1;
  size_t genome_end = bwt_anchor_next->start - 1;

  int read_gap = read->length - (anchor_len_prev + anchor_len_next);
  int genome_gap = genome_end - genome_start + 1;

  //printf("len_prev=%i / len_next=%i \n", anchor_len_prev, anchor_len_next);

  if (genome_gap <= 1 || read_gap <= 1) { 
    *l_extend = 0;
    *r_extend = 0;
    return -1; 
  }
  
  int num_err = 0;
  char *genome_ref = (char *)malloc(sizeof(char)*(genome_gap + 1));
  char *read_ref = (char *)calloc(read_gap + 1, sizeof(char));

  //printf("%lu-%lu = %lu, chr = %i\n", genome_start, genome_end, genome_gap, bwt_anchor_prev->chromosome);

  genome_read_sequence_by_chr_index(genome_ref, 0,
				    bwt_anchor_prev->chromosome, &genome_start, &genome_end, genome);
  memcpy(read_ref, sequence + anchor_len_prev, read_gap);

  LOG_DEBUG_F("REFERENCE : %s (%i)\n", genome_ref, read_gap);
  LOG_DEBUG_F("SEQUENCE  : %s (%i)\n", read_ref, genome_gap);

  int read_l = 0, read_r = read_gap - 1;
  int gen_l = 0, gen_r = genome_gap - 1;
  int not_close = 0;

  while (1) {
    printf("Distance gen=%i/%i, read=%i/%i, %i (%c == %c) (%c == %c)\n", gen_l, gen_r, read_l, read_r, num_err, 
	   genome_ref[gen_l], read_ref[read_l], genome_ref[gen_r], read_ref[read_r]);
    if (genome_ref[gen_l++] != read_ref[read_l++]) { num_err++; }
    if (num_err > max_err ) { not_close = 1; break; }

    if (genome_ref[gen_r--] != read_ref[read_r--]) { num_err++; }
    if (num_err > max_err) { not_close = 1; break; }
    
    if (gen_l >= gen_r || read_l >= read_r) { break; }
  }

  free(genome_ref);
  free(read_ref);

  *l_extend = read_l;
  *r_extend = read_r;

  if (not_close) {
    return -1;
  } else {
    printf("Close gap %i == %i\n", genome_gap, read_gap);
    if (genome_gap == read_gap) { return num_err; }
    else { return -1; }
  }  
}

int fill_extrem_gap (char *genome_ref, char *read_ref) {
  int len = strlen(genome_ref);
  int num_err = 0;
  for (int i = 0; i < len; i++) {
    if (genome_ref[i] != read_ref[i]) { num_err++; }
  }
  return num_err;  
}

typedef struct sp_info_sw {
  void *info;
} sp_info_sw;

typedef struct sw_item {
  int type_sw;
  int read_id;
  bwt_anchor_t *anchor_prev;
  bwt_anchor_t *anchor_next;
  sp_info_sw *info;
} sw_item_t;

sw_item_t *sw_item_new(int type_sw, int read_id, 
		       bwt_anchor_t *anchor_prev, 
		       bwt_anchor_t *anchor_next, 
		       sp_info_sw *info) {
  sw_item_t *sw_item = (sw_item_t *)malloc(sizeof(sw_item_t));
  
  sw_item->type_sw = type_sw;
  sw_item->read_id = read_id;
  sw_item->anchor_prev = anchor_prev;
  sw_item->anchor_next = anchor_next;
  sw_item->info = info;
  
  return sw_item;
}

void sw_item_free(sw_item_t *sw_item) {
  free(sw_item);
}

typedef struct sw_depth {
  char *q[MAX_DEPTH];
  char *r[MAX_DEPTH];
  sw_item_t *items[MAX_DEPTH];
  int depth;
} sw_depth_t;


void sw_depth_process(sw_optarg_t *sw_optarg, sw_multi_output_t *output, 
		      array_list_t **alignments_list, array_list_t *fq_batch, 
		      sw_depth_t *sw_depth, int step) {
  int distance, len;
  float norm_score;
  float match = sw_optarg->subst_matrix['A']['A'];
  if (sw_depth->depth == MAX_DEPTH || 
      (step == SW_FINAL && sw_depth->depth > 0)) {
    smith_waterman_mqmr(sw_depth->q, sw_depth->r, sw_depth->depth, sw_optarg, 1, output);
    for (int i = 0; i < sw_depth->depth; i++) {
      sw_item_t *sw_item = sw_depth->items[i];     
      fastq_read_t *read = array_list_get(sw_item->read_id, fq_batch);
      LOG_DEBUG_F("Read : %s | %s\n", read->id, read->sequence);
      printf("%i-QUE: %s\n", i, output->query_map_p[i]);
      printf("%i-REF: %s\n", i, output->ref_map_p[i]);

      if (sw_item->type_sw == SIMPLE_SW) {
	output->score_p[i] += (sw_item->anchor_prev->end - sw_item->anchor_prev->start)*match;
	output->score_p[i] += (sw_item->anchor_next->end - sw_item->anchor_next->start)*match;

	norm_score = NORM_SCORE(output->score_p[i], read->length, match);
	LOG_DEBUG_F("SW REF: %s (%f)\n", output->ref_map_p[i], norm_score);
	LOG_DEBUG_F("SW SEQ: %s\n", output->query_map_p[i]);
	cigar_code_t *cigar_code = generate_cigar_code(output->query_map_p[i], output->ref_map_p[i], 
						       strlen(output->query_map_p[i]),
						       output->query_start_p[i], output->ref_start_p[i],
						       strlen(sw_depth->q[i]), strlen(sw_depth->r[i]),
						       &distance, MIDDLE_SW);
	len = sw_item->anchor_prev->end - sw_item->anchor_prev->start;
	cigar_code_insert_first_op(cigar_op_new(len, 'M'), cigar_code);
	
	len = sw_item->anchor_next->end - sw_item->anchor_next->start;
	cigar_code_append_op(cigar_op_new(len, 'M'), cigar_code);

	printf("read len = %i, cigar len = %i\n", read->length, cigar_code_validate(read->length, cigar_code));
	if (!cigar_code_validate(read->length, cigar_code)) {
	  LOG_FATAL_F("   %s : %s\n", read->id, new_cigar_code_string(cigar_code));
	}

	alignment_t *alignment = alignment_new();
	alignment_init_single_end(strdup(&read->id[1]), strdup(read->sequence), strdup(read->quality),
				  sw_item->anchor_prev->strand,
				  sw_item->anchor_prev->chromosome - 1, sw_item->anchor_prev->start - 1,
				  strdup(new_cigar_code_string(cigar_code)), 
				  cigar_code_get_num_ops(cigar_code), norm_score * 254, 1, 1,
				  0, NULL, 0, alignment);
	array_list_insert(alignment, alignments_list[sw_item->read_id]);
	array_list_clear(cigar_code->ops, cigar_op_free);
	cigar_code_free(cigar_code);
      } else {
	continue;
      }
    }
    sw_depth->depth = 0;
  }
}

void sw_depth_insert(char *query, char *reference, sw_item_t *sw_item, 
		     sw_optarg_t *sw_optarg, sw_multi_output_t *output, 
		     array_list_t **alignments_list, array_list_t *fq_batch, 
		     sw_depth_t *sw_depth) {
  sw_depth->q[sw_depth->depth]     = strdup(query);
  sw_depth->r[sw_depth->depth]     = strdup(reference);
  sw_depth->items[sw_depth->depth++] = sw_item;
  printf("Insert SW depth %i\n", sw_depth->depth);
  printf("QUE: %s\n", query);
  printf("REF: %s\n", reference);
  sw_depth_process(sw_optarg, output, 
		   alignments_list, fq_batch, 
		   sw_depth, SW_NORMAL);
}



int first_phase(bwt_server_input_t* input, batch_t *batch) {
  struct timeval start, end;
  double time;
  metaexons_t *metaexons = input->metaexons;
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  size_t num_reads = array_list_size(mapping_batch->fq_batch);
  size_t num_mappings;
  array_list_t *list;
  size_t *unmapped_indices = mapping_batch->targets;
  size_t num_unmapped = 0;
  size_t num_anchors;
  sw_optarg_t *sw_optarg = input->sw_optarg;
  float mismatch_penalty = sw_optarg->subst_matrix['A']['C'];
  genome_t *genome = input->genome;
  sw_depth_t sw_depth;
  
  sw_depth.depth = 0;

  sw_multi_output_t *output = sw_multi_output_new(MAX_DEPTH);
  array_list_t **alignments_list = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
  alignment_t *alignment;
  const int flank = 20;
  for (int i = 0; i < num_reads; i++) {
    alignments_list[i] = array_list_new(20, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  }

  for (int i = 0; i < num_reads; i++) {
    fastq_read_t *read = array_list_get(i, mapping_batch->fq_batch);
    printf("%s\n", read->id);
    list = mapping_batch->mapping_lists[i];    
    array_list_set_flag(1, list);
    num_mappings = bwt_map_inexact_read(read,
					input->bwt_optarg_p,
					input->bwt_index_p,
					list);
    if (array_list_get_flag(list) != 2) { //If flag 2, the read exceded the max number of mappings
      if (array_list_get_flag(list) == 1) {
	if (num_mappings > 0) {
	  num_anchors = bwt_search_pair_anchors(list, read->length);	
	  if (num_anchors == 0) {
	    array_list_set_flag(0, list);
	  }
	} else {
	  array_list_set_flag(0, list);
	}
	//printf("tot anchors found %i %s\n", num_anchors, read->id);
	unmapped_indices[num_unmapped++] = i;
      } else if (num_mappings <= 0) {
	array_list_set_flag(0, list);
	//printf("Read NO Mapped %i %s\n", num_anchors, read->id);
	unmapped_indices[num_unmapped++] = i;
      } else {
	//Read Map, Store position 
	array_list_set_flag(3, list);
	for (int i = 0; i < num_mappings; i++) {
	  alignment = array_list_get(i, list);
	  metaexon_insert(alignment->seq_strand, alignment->chromosome,
			  alignment->position, alignment->position + read->length, 40,
			  METAEXON_NORMAL, NULL,
			  metaexons);
	}
      }
    }
    

    if (array_list_get_flag(list) == DOUBLE_ANCHORS) {
      //Double anchor found
      int found_prev;
      int found_next;
      metaexon_t *metaexon_prev, *metaexon_next;
      size_t search_prev, search_next;
      size_t anchor_len_prev, anchor_len_next;
      int read_gap, genome_gap;
      char q[2048];
      char r[2048];

      size_t start_exon_l, end_exon_l;
      size_t start_exon_r, end_exon_r;

      size_t genome_start, genome_end;
      sw_item_t *sw_item;
      alignment_t *alignment;
      
      int flank_ex = 10;
      char *sequence;
      char *rev_comp = NULL;

      for (int j = 0; j < array_list_size(list); j += 2) {
	bwt_anchor_t *bwt_anchor_prev = array_list_get(j, list);
	bwt_anchor_t *bwt_anchor_next = array_list_get(j + 1, list);
	if (bwt_anchor_prev->strand == 1) {
	  rev_comp = (char *) calloc(read->length + 1, sizeof(char));
          strcpy(rev_comp, read->sequence);
          seq_reverse_complementary(rev_comp, read->length);
	  sequence = rev_comp;
	} else {
	  sequence = read->sequence;
	}
	//1st: Search if exist one Metaexon with finals completed
	anchor_len_prev = (bwt_anchor_prev->end - bwt_anchor_prev->start) + 1;
	anchor_len_next = (bwt_anchor_next->end - bwt_anchor_next->start) + 1;
	genome_start = bwt_anchor_prev->end + 1;
	genome_end   = bwt_anchor_next->start - 1;
	
	read_gap = read->length - (anchor_len_prev + anchor_len_next);
	genome_gap = genome_end - genome_start + 1;

	if (read_gap <= 0 || genome_gap <= 0) { 
	  bwt_anchor_prev->end -= flank_ex;
	  bwt_anchor_next->start += flank_ex;
	  anchor_len_prev -= flank_ex;
	  anchor_len_next -= flank_ex;
	  
	  read_gap = read->length - (anchor_len_prev + anchor_len_next);
	  genome_gap = bwt_anchor_next->start - bwt_anchor_prev->end;
	}

	LOG_DEBUG_F(" ANCHOR FORWARD (%i)[%i:%lu-%lu] = %i\n", bwt_anchor_prev->strand, bwt_anchor_prev->chromosome,
		    bwt_anchor_prev->start, bwt_anchor_prev->end, anchor_len_prev);

	LOG_DEBUG_F(" ANCHOR BACKWARD (%i)[%i:%lu-%lu] = %i\n", bwt_anchor_next->strand, bwt_anchor_next->chromosome,
		    bwt_anchor_next->start, bwt_anchor_next->end, anchor_len_next);
	
	LOG_DEBUG_F(" READ_GAP = %i, GENOME_GAP = %i\n", read_gap, genome_gap);
	
	if (abs(genome_gap - read_gap) <= 40) { //min_intron_size
	  //Exonic Read, Complete metaexons structure and process
	  metaexon_insert(bwt_anchor_prev->strand, bwt_anchor_prev->chromosome,
			  bwt_anchor_prev->start, bwt_anchor_next->end, 40,
			  METAEXON_NORMAL, NULL,
			  metaexons);	 
	  int l_extend, r_extend;
	  int num_err;	    

	  num_err = fill_anchor_gap(bwt_anchor_prev, bwt_anchor_next, 
				    genome, read, sequence, read_gap / 3, 
				    &l_extend, &r_extend);
	  if (num_err == -1) {
	    //We need smith-waterman to close gap
	    sw_item = sw_item_new(SIMPLE_SW, i, 
				  bwt_anchor_prev, 
				  bwt_anchor_next, 
				  NULL);
	    memcpy(q, sequence + anchor_len_prev, read_gap);
	    q[read_gap] = '\0';

	    genome_read_sequence_by_chr_index(r, 0,
					      bwt_anchor_prev->chromosome, &genome_start, &genome_end, genome);
	    // printf("sw-REF : %s\n", r);
	    // printf("SW-SEQ : %s\n", q);
	    
	    sw_depth_insert(q, r, sw_item, 
			    sw_optarg, output, 
			    alignments_list,
			    mapping_batch->fq_batch,
			    &sw_depth);
	  } else {
	    //Not need smith-waterman gap closed, report alignment
	    //printf("Insert ok %i, %f, %f\n", num_err, mismatch_penalty, 254 - (num_err*mismatch_penalty*254)/(read->length*mismatch_penalty));
	    alignment = alignment_new();
	    char *cigar_str[128];
	    sprintf(cigar_str, "%iM", strlen(read->sequence));
	    alignment_init_single_end(strdup(&read->id[1]), strdup(read->sequence), strdup(read->quality),
				      bwt_anchor_prev->strand,
				      bwt_anchor_prev->chromosome - 1, bwt_anchor_prev->start - 1,
				      strdup(cigar_str), 1, 
				      254 - (num_err*mismatch_penalty*254)/(read->length*mismatch_penalty), 
				      1, 1,
				      0, NULL, 0, alignment);
	    array_list_insert(alignment, alignments_list[i]);
	  }
	} else {
	  //SP Read
	  continue;
	  search_prev = bwt_anchor_prev->start + anchor_len_prev - (anchor_len_prev / 2);
	  search_next = bwt_anchor_next->start + (anchor_len_next / 2);

	  printf("Position Search Prev(%i)[%i:%lu]: ", bwt_anchor_prev->strand,
		 bwt_anchor_prev->chromosome, search_prev);
	  
	  found_prev = metaexon_search(bwt_anchor_prev->strand, bwt_anchor_prev->chromosome, 
				       search_prev, &metaexon_prev, metaexons);
	  if (found_prev) {
	    found_prev = metaexon_prev->right_closed;
	  }
	  
	  printf("Position Search Next(%i)[%i:%lu]: ", 
		 bwt_anchor_next->strand, bwt_anchor_next->chromosome, 
		 search_next);

	  found_next = metaexon_search(bwt_anchor_next->strand, bwt_anchor_next->chromosome, 
				       search_next, &metaexon_next, metaexons);
	  
	  if (found_next) {
	    found_next = metaexon_next->left_closed;
	  }
	  	  
	  if (found_prev && !found_next) {
	    //Found Forward, but not Backward... We need smith-waterman to found SP
	    start_exon_l = metaexon_prev->end - flank;
	    end_exon_l = metaexon_prev->end + flank;
	    
	    
	  } else if (!found_prev && found_next) {
	    //Found Backward, but not Forward... We need smith-waterman to found SP
	    start_exon_r = metaexon_next->end - flank;
	    end_exon_r = metaexon_next->end + flank;
	    
	    
	  } else if (!found_prev && !found_next) {
	    //Normal Case ¿¿Seeding?? ¿¿Caling??...
	    
	  } else { 
	    //Know splice junction, 
	    int report_alig = 0;
	    if (anchor_len_prev + anchor_len_next > read->length) {
	      //Some nucleotides extend when intron starts, Recalculate anchors coords, and Report Alignment!
	      bwt_anchor_prev->end = metaexon_prev->end;
	      bwt_anchor_next->start = metaexon_next->start;
	      retport_alig = 1;
	    } else if (anchor_len_prev + anchor_len_next == read->length) {
	      //Exact anchors found with not errors! Report Alignment!
	      retport_alig = 1;
	    } else {
	      //Extend... and close gap, if not close... We need smith-waterman
	      int l_extend = metaexon_prev->end - bwt_anchor_prev->end;
	      int found_l = 0, found_r = 0;
	      if (l_extend < 0) {
		l_extend = 0;
		found_l = 1;
	      } else {
		//For to close gap to end Exon [|-->]
		int num_err;
		char genome_ref[l_extend + 1];
		genome_start = bwt_anchor_prev->end;
		genome_end = metaexon_prev->end;
		genome_read_sequence_by_chr_index(genome_ref, 0,
						  bwt_anchor_prev->chromosome - 1, &genome_start, &genome_end, genome);
		memcpy(read_ref, read->sequence + anchor_len_prev, l_extend);		
		num_err = fill_extrem_gap(genome_ref, read_ref);		
		if (num_err > (l_extend / 3) + 1) {
		  //We need smith-waterman to resolve final left Exon
		  
		} else {		  
		  bwt_anchor_prev->end  = metaexon_prev->end;
		  found_l = 1;
		}
	      }

	      int r_extend = bwt_anchor_next->start - metaexon_next->start;
	      if (r_extend < 0) {
		r_extend = 0;
		found_r = 1;
	      } else {
		//For to close gap to start Exon [<--|]
		int num_err = 0;
		char genome_ref[r_extend + 1];
		genome_start = metaexon_next->start;
		genome_end = bwt_anchor_next->start;
		genome_read_sequence_by_chr_index(genome_ref, 0,
						  bwt_anchor_prev->chromosome - 1, &genome_start, &genome_end, genome);
		memcpy(read_ref, (read->sequence + read->length) - r_extend, r_extend);		
		num_err = fill_extrem_gap(genome_ref, read_ref);		
		if (num_err > (l_extend / 3) + 1) {
		  //We need smith-waterman to resolve start right Exon
		  
		} else {
		  bwt_anchor_next->start = metaexon_next->start;
		  found_r = 1;
		}
	      }

	      if (found_l && found_r) {
		//Report Alignment
		report_alig = 1;
	      }	      
	    }
	  
	    if (report_alig) {
	      alignment = alignment_new();
	      char *cigar_str[128];
	      sprintf(cigar_str, "%iM%iN%iM", bwt_anchor_prev->end - bwt_anchor_prev->start, 
		      bwt_anchor_next->start - bwt_anchor_prev->end,  bwt_anchor_next->end - bwt_anchor_next->start);
	      alignment_init_single_end(strdup(&read->id[1]), strdup(read->sequence), strdup(read->quality),
					cal->strand,
					cal->chromosome_id - 1, cal->start - 1,
					strdup(cigar_str), 1, 254, 1, 1,
					0, NULL, 0, alignment);
	      array_list_insert(alignment, alignments_list[i]);	      	      
	    }

	  } //Else know splice junction	
	  
	} // Else Exonic Read or SP Read
      } //End loop anchors

      if (rev_comp != NULL) { free(rev_comp); }
 
    } //End if array_list_flag    
  } //End loop reads

  sw_depth_process(sw_optarg, output, 
		   alignments_list, mapping_batch->fq_batch, 
		   &sw_depth, SW_FINAL);
  // array_list flag: 0 -> Not  BWT Anchors found
  //                  1 -> One  BWT Anchors found
  //                  2 -> Pair BWT Anchors found


  for (int i = 0; i < num_reads; i++) {
    list = mapping_batch->mapping_lists[i];    
    if (array_list_get_flag(list) == 1 || 
	array_list_get_flag(list) == 2) {
      array_list_clear(list, bwt_anchor_free);
    }
    if (array_list_get_flag(list) == 2) {
      for (int j = 0; j < array_list_size(alignments_list[i]); j++) {
	alignment = array_list_get(j, alignments_list[i]);
	array_list_insert(alignment, list);
      }
    }
    array_list_free(alignments_list[i], NULL);
  }

  
  free(alignments_list);

  //show_metaexons(metaexons);
  
  mapping_batch->num_targets = num_unmapped;

  return CONSUMER_STAGE;

}
*/
//------------------------------------------------------------------------------------

