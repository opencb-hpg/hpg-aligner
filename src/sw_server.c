#include "sw_server.h"

//====================================================================================
//  Input structure for Smith-Waterman server
//====================================================================================

void sw_server_input_init(list_t* sw_list, list_t* alignment_list, unsigned int write_size, 
			  float match, float mismatch, float gap_open, float gap_extend, 
			  float min_score, unsigned int flank_length, genome_t* genome, 
			  size_t max_intron_size, int min_intron_size, 
			  size_t seed_max_distance, bwt_optarg_t* bwt_optarg_p, 
			  allocate_splice_elements_t *chromosome_avls_p, 
			  sw_server_input_t* input) {
  
  input->sw_list_p = sw_list;
  input->alignment_list_p = alignment_list;
  input->write_size = write_size;
  input->genome_p = genome;
  input->max_intron_size = max_intron_size;
  input->min_intron_size = min_intron_size;
  input->seed_max_distance = seed_max_distance;
  input->bwt_optarg_p =  bwt_optarg_p; 

  // Smith-Waterman parameters
  input->match = match;
  input->mismatch = mismatch;
  input->gap_open = gap_open;
  input->gap_extend = gap_extend;
  input->min_score = min_score;

  input->sw_optarg.gap_open = gap_open;
  input->sw_optarg.gap_extend = gap_extend;

  input->sw_optarg.subst_matrix['A']['A'] = input->match;
  input->sw_optarg.subst_matrix['C']['A'] = input->mismatch;
  input->sw_optarg.subst_matrix['T']['A'] = input->mismatch;
  input->sw_optarg.subst_matrix['G']['A'] = input->mismatch;
  input->sw_optarg.subst_matrix['N']['A'] = input->mismatch;

  input->sw_optarg.subst_matrix['A']['C'] = input->mismatch;
  input->sw_optarg.subst_matrix['C']['C'] = input->match;
  input->sw_optarg.subst_matrix['T']['C'] = input->mismatch;
  input->sw_optarg.subst_matrix['G']['C'] = input->mismatch;
  input->sw_optarg.subst_matrix['N']['C'] = input->mismatch;

  input->sw_optarg.subst_matrix['A']['T'] = input->mismatch;
  input->sw_optarg.subst_matrix['C']['T'] = input->mismatch;
  input->sw_optarg.subst_matrix['T']['T'] = input->match;
  input->sw_optarg.subst_matrix['G']['T'] = input->mismatch;
  input->sw_optarg.subst_matrix['N']['T'] = input->mismatch;

  input->sw_optarg.subst_matrix['A']['G'] = input->mismatch;
  input->sw_optarg.subst_matrix['C']['G'] = input->mismatch;
  input->sw_optarg.subst_matrix['T']['G'] = input->mismatch;
  input->sw_optarg.subst_matrix['G']['G'] = input->match;
  input->sw_optarg.subst_matrix['N']['G'] = input->mismatch;

  input->sw_optarg.subst_matrix['A']['N'] = input->mismatch;
  input->sw_optarg.subst_matrix['C']['N'] = input->mismatch;
  input->sw_optarg.subst_matrix['T']['N'] = input->mismatch;
  input->sw_optarg.subst_matrix['G']['N'] = input->mismatch;
  input->sw_optarg.subst_matrix['N']['N'] = input->match;

  // CAL
  input->flank_length = flank_length;
  input->chromosome_avls_p = chromosome_avls_p;
}

//====================================================================================
//  Smith-Waterman channel for SIMD implementation
//====================================================================================

inline void sw_channel_allocate_ref(unsigned int length, sw_channel_t* channel_p) {
     if (channel_p == NULL) return;
     
     if (channel_p->allocated_ref_size < length) {
	  if (channel_p->ref_p == NULL) {
	       channel_p->ref_p = (char*) calloc(length, sizeof(char));
	  } else {
	       free(channel_p->ref_p);
	       channel_p->ref_p = (char*) calloc(length, sizeof(char));
	  }
	  channel_p->allocated_ref_size = length;
     }
}

//------------------------------------------------------------------------------------

inline void sw_channel_update(size_t read_index, unsigned int cal_index, unsigned int read_len,
			      unsigned int header_len, unsigned int ref_len, sw_channel_t *channel_p) {
     channel_p->read_index = read_index;
     channel_p->cal_index = cal_index;
     channel_p->read_len = read_len;
     channel_p->header_len = header_len;
     channel_p->ref_len = ref_len;
}

//====================================================================================
// main sw function
//====================================================================================

void set_sw_sequences(char **q, char **r, size_t sw_count, char *sequence, 
		      genome_t *genome, int chromosome, seed_region_t *sr) {
  int gap_len;

  // get query sequence, revcomp if necessary
  gap_len = sr->read_end - sr->read_start + 1;
  q[sw_count] = (char *) malloc((gap_len + 1) * sizeof(char));
  memcpy(q[sw_count], &sequence[sr->read_start], gap_len);
  q[sw_count][gap_len] = '\0';

  // get ref. sequence
  gap_len = sr->genome_end - sr->genome_start + 1;
  r[sw_count] = (char *) malloc((gap_len + 1) * sizeof(char));
  genome_read_sequence_by_chr_index(r[sw_count], 0, chromosome, 
  				    &sr->genome_start, &sr->genome_end, genome);
  r[sw_count][gap_len] = '\0';
}

//------------------------------------------------------------------------------------
// apply_sw
//------------------------------------------------------------------------------------

int apply_sw(sw_server_input_t* input, batch_t *batch) {

  int tid = omp_get_thread_num();
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  cal_t *cal = NULL;
  array_list_t *cal_list = NULL, *mapping_list = NULL;//, *old_list = NULL, *new_list = NULL;

  array_list_t *fq_batch = mapping_batch->fq_batch;
  fastq_read_t *fq_read;

  size_t start, end;
  genome_t *genome = input->genome_p;
     
  size_t flank_length = input->flank_length;

  // SIMD support for Smith-Waterman
  float score, norm_score, min_score = input->min_score;
  
  size_t read_index, num_cals;

  size_t num_targets = mapping_batch->num_targets;
  size_t new_num_targets = 0;

  size_t sw_total = mapping_batch->num_to_do;

  // set to zero
  mapping_batch->num_to_do = 0;

  sw_optarg_t *sw_optarg = &input->sw_optarg;
  sw_multi_output_t *output = sw_multi_output_new(sw_total);
  char *q[sw_total], *r[sw_total];
  int sw_count = 0, run_sw;
  int read_len, ref_len, max_ref_len;
  int gap_read_len, gap_genome_len;
  size_t gap_read_start, gap_read_end, gap_genome_start, gap_genome_end;

  seed_region_t *s, *prev_s, *new_s;
  linked_list_iterator_t* itr;

  char *revcomp_seq;

  //  LOG_DEBUG("\n\n P R E   -   P R O C E S S\n");

  // initialize query and reference sequences to Smith-Waterman
  for (size_t i = 0; i < num_targets; i++) {

    read_index = mapping_batch->targets[i];
    fq_read = (fastq_read_t *) array_list_get(read_index, fq_batch);

    cal_list = mapping_batch->mapping_lists[read_index];
    num_cals = array_list_size(cal_list);

    read_len = fq_read->length;

    revcomp_seq = strdup(fq_read->sequence);
    seq_reverse_complementary(revcomp_seq, read_len);

    //    LOG_DEBUG_F("%s: seq : %s\n", fq_read->id, fq_read->sequence);

    // processing each CAL from this read
    for(size_t j = 0; j < num_cals; j++) {

      // get cal and read index
      cal = array_list_get(j, cal_list);

      //      LOG_DEBUG_F("cal %i of %i (strand %i)\n", j, num_cals, cal->strand);

      prev_s = NULL;

      itr = linked_list_iterator_new(cal->sr_list);
      s = (seed_region_t *) linked_list_iterator_curr(itr);
      while (s != NULL) {
	//	LOG_DEBUG_F("\tseed: [%i|%i - %i|%i]\n", 
	//		    s->genome_start, s->read_start, s->read_end, s->genome_end);

	if ((prev_s == NULL && s->read_start != 0) || (prev_s != NULL)) {
	  mapping_batch->num_gaps++;
	  run_sw = 0;
	  if (prev_s == NULL) {
	    // gap at the first position
	    gap_read_start = 0;
	    gap_read_end = s->read_start - 1;

	    gap_genome_start = s->genome_start - s->read_start;
	    gap_genome_end = s->genome_start - 1;

	    gap_read_len = gap_read_end - gap_read_start + 1;
	    gap_genome_len = gap_genome_end - gap_genome_start + 1;

	    // always run a SW in this case
	    run_sw = 1;
	    mapping_batch->num_ext_sws++;
	  } else {
	    assert(prev_s->read_end < s->read_start);

	    // gap in a middle position
	    gap_read_start = prev_s->read_end + 1;
	    gap_read_end = s->read_start - 1;

	    gap_genome_start = prev_s->genome_end + 1;
	    gap_genome_end = s->genome_start - 1;

	    gap_read_len = gap_read_end - gap_read_start + 1;
	    gap_genome_len = gap_genome_end - gap_genome_start + 1;

	    // only run a SW if gaps have different lengths
	    if (gap_read_len != gap_genome_len) {
	      run_sw = 1;
	    }
	  }
	  
	  //	  LOG_DEBUG_F("\t\tgap_read_len = %i, gap_genome_len = %i\n", gap_read_len, gap_genome_len);
	  //	  LOG_DEBUG_F("\t\t%i of %i: [%lu|%lu - %lu|%lu]\n", 
	  //		      sw_count, sw_total, gap_genome_start, gap_read_start, gap_read_end, gap_genome_end);

	  // insert gap in the list
	  new_s = seed_region_new(gap_read_start, gap_read_end, gap_genome_start, gap_genome_end, 0);
	  new_s->gap = 1;
	  linked_list_iterator_insert(new_s, itr);

	  if (run_sw == 1) {
	    mapping_batch->num_sws++;
	    new_s->run_sw = 1;
	    
	    // generate sequences for sw
	    if (cal->strand == 1) {
	      set_sw_sequences(q, r, sw_count, revcomp_seq, genome, cal->chromosome_id - 1, new_s);
	    } else {
	      set_sw_sequences(q, r, sw_count, fq_read->sequence, genome, cal->chromosome_id - 1, new_s);
	    }
	    
	    //	  LOG_DEBUG_F("\t\tquery-seq : %s\n", q[sw_count]);
	    //	  LOG_DEBUG_F("\t\tref-seq   : %s\n", r[sw_count]);
	    
	    // increase counter
	    sw_count++;	  
	  }
	}

	//continue loop...
	prev_s = s;
	linked_list_iterator_next(itr);
	s = linked_list_iterator_curr(itr);
      }

      if (prev_s != NULL && prev_s->read_end < read_len - 1) { 
	//	LOG_DEBUG("\tgap at the last position...");
	mapping_batch->num_gaps++;
	mapping_batch->num_sws++;
	mapping_batch->num_ext_sws++;

	// gap at the last position
	gap_read_start = prev_s->read_end + 1;
	gap_read_end = read_len - 1;
	gap_read_len = gap_read_end - gap_read_start + 1;
	
	gap_genome_len = gap_read_len;
	gap_genome_start = prev_s->genome_end + 1;
	gap_genome_end = gap_genome_start + gap_genome_len - 1;

	//	LOG_DEBUG_F("\t\tgap_read_len = %i, gap_genome_len = %i\n", gap_read_len, gap_genome_len);
	//	LOG_DEBUG_F("\t\t%i of %i: [%lu|%lu - %lu|%lu]\n", 
	//		    sw_count, sw_total, gap_genome_start, gap_read_start, gap_read_end, gap_genome_end);

	// insert gap in the list
	new_s = seed_region_new(gap_read_start, gap_read_end, gap_genome_start, gap_genome_end, 0);
	new_s->gap = 1;
	new_s->run_sw = 1;
	linked_list_insert_last(new_s, cal->sr_list);
	
	// generate sequences for sw
	if (cal->strand == 1) {
	  set_sw_sequences(q, r, sw_count, revcomp_seq, genome, cal->chromosome_id - 1, new_s);
	} else {
	  set_sw_sequences(q, r, sw_count, fq_read->sequence, genome, cal->chromosome_id - 1, new_s);
	}
	
	//	LOG_DEBUG_F("\t\tquery-seq : %s\n", q[sw_count]);
	//	LOG_DEBUG_F("\t\tref-seq   : %s\n", r[sw_count]);
	
	// increase counter
	sw_count++;	  
      }
      linked_list_iterator_free(itr);

      /*      
      //      LOG_DEBUG("\tUpdated list:");
      itr = linked_list_iterator_new(cal->sr_list);
      s = (seed_region_t *) linked_list_iterator_curr(itr);
      while (s != NULL) {
	//	LOG_DEBUG_F("\t\t%s [%i|%i - %i|%i]\n", (s->gap ? "gap   " : "region"),
	//		    s->genome_start, s->read_start, s->read_end, s->genome_end);
	linked_list_iterator_next(itr);
	s = linked_list_iterator_curr(itr);
      }
      linked_list_iterator_free(itr);
      */
    }
    // free memory
    free(revcomp_seq);
  }

  //  LOG_DEBUG_F("R U N   S W (sw_total = %i, sw_count = %i)\n", sw_total, sw_count);
  assert(mapping_batch->num_sws == sw_count);
  //  mapping_batch->num_gaps += sw_total;

  // run Smith-Waterman
  //  smith_waterman_mqmr(q, r, sw_count, sw_optarg, 1, output);
  smith_waterman_mqmr(q, r, sw_count, sw_optarg, 1, output);

  //  LOG_DEBUG("\n\n P O S T   -   P R O C E S S\n");

  // re-construct the CAL score from region gaps
  sw_count = 0;
  int header_len, mquery_start, mref_start;
  char *header_match, *read_match, *quality_match;
  char *p, *optional_fields;
  int pos, optional_fields_length, distance, AS, first, cigar_number;

  cigar_op_t* cigar_op;
  cigar_code_t *cigar_code;

  alignment_t *alignment;
  array_list_t *alignment_list;

  for (size_t i = 0; i < num_targets; i++) {
    
    //    LOG_DEBUG_F("targets %i of %i\n", i, num_targets);
    
    read_index = mapping_batch->targets[i];
    fq_read = (fastq_read_t *) array_list_get(read_index, fq_batch);
    
    cal_list = mapping_batch->mapping_lists[read_index];
    num_cals = array_list_size(cal_list);
    
    //    LOG_DEBUG_F("read_index =  %i, num_cals = %i\n", read_index, num_cals);

    read_len = fq_read->length;
    header_match = (char *) calloc(strlen(fq_read->id) + 1, sizeof(char));
    header_len = get_to_first_blank(fq_read->id, strlen(fq_read->id), header_match);
    header_match[header_len] = '\0';    

    alignment_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    array_list_set_flag(1, alignment_list);

    // processing each CAL from this read
    score = 0.0f;
    for(size_t j = 0; j < num_cals; j++) {
      // get cal and read index
      cal = array_list_get(j, cal_list);

      //      LOG_DEBUG_F("cal %i of %i\n", j, num_cals);

      cigar_code = cigar_code_new();
      distance = 0;

      read_match = (char *) calloc((read_len * 2), sizeof(char));
      quality_match = (char *) calloc((read_len * 2), sizeof(char));
      /*
      first = 1;
      for (linked_list_item_t *list_item = cal->sr_list->first; list_item != NULL; list_item = list_item->next) {
	s = list_item->item;
	//	LOG_DEBUG_F("is gap ? %s -> [%i|%i - %i|%i]\n", 
	//		    (s->gap ? "yes" : "no"), s->genome_start, s->read_start, s->read_end, s->genome_end);
	if (s->gap) {	  
	  // gap found !!
	  if (s->run_sw) {
	    int distance;
	    int read_gap_len = s->read_end - s->read_start + 1;
	    int genome_gap_len = s->genome_end - s->genome_start + 1;
	
	    LOG_DEBUG_F("\tgap (read %lu-%lu, genome %lu-%lu) = (%i, %i)\n", 
			s->read_start, s->read_end, s->genome_start, s->genome_end,
			read_gap_len, genome_gap_len);
	    LOG_DEBUG_F("\tmquery: %s (start %i)\n", output->query_map_p[sw_count], output->query_start_p[sw_count]);
	    LOG_DEBUG_F("\tmref  : %s (start %i)\n", output->ref_map_p[sw_count], output->ref_start_p[sw_count]);
	
	    cigar_code_t *cigar_c = generate_cigar_code(output->query_map_p[sw_count], output->ref_map_p[sw_count],
							strlen(output->query_map_p[sw_count]), output->ref_start_p[sw_count],
							read_gap_len, &distance);

	    //	    LOG_DEBUG_F("\tscore : %0.2f, cigar: %s (distance = %i)\n", 
	    //		output->score_p[sw_count], cigar_code_get_string(cigar_c), distance);
	    cigar_code_free(cigar_c);
	
	    score += output->score_p[sw_count];
	    if (first) {
	      mquery_start = output->query_start_p[sw_count];
	      mref_start = output->ref_start_p[sw_count];
	    }
	    // free query and reference
	    free(q[sw_count]);
	    free(r[sw_count]);
	    sw_count++;
	  } else {
	    
	  }
	} else {
	  // no gap: exact region
	  cigar_number = s->read_end - s->read_start + 1;
	  cigar_op = cigar_code_get_last_op(cigar_code);
	  if (cigar_op && cigar_op->name == 'M') {
	    cigar_op->number += cigar_number;
	  } else {
	    cigar_code_append_op(cigar_op_new(cigar_number, 'M'), cigar_code);
	  }
	  score += (cigar_number * input->match);
	  if (first) {
	    mquery_start = s->read_start;
	    mref_start = mquery_start;
	  }
	}
	first = 0;
      }

      memcpy(read_match, fq_read->quality, read_len);
      read_match[read_len] = '\0';

      memcpy(quality_match, fq_read->quality, read_len);
      quality_match[read_len] = '\0';

      // filter by SW score
      norm_score = NORM_SCORE(score, read_len, input->match);

      //      LOG_DEBUG_F("score: %0.2f -> (%0.2f, min. %0.2f), mquery_start = %i, mref_start = %i\n",
      //		  score, norm_score, min_score, mquery_start, mref_start);
      */
      if (0) {
	//      if (norm_score >= min_score) {
	
        pos = cal->start + mref_start - 1;

	// set optional fields
	optional_fields_length = 100;
	optional_fields = (char *) calloc(optional_fields_length, sizeof(char));
	
	p = optional_fields;
	AS = (int) score;
	
	sprintf(p, "ASi");
	p += 3;
	memcpy(p, &AS, sizeof(int));
	p += sizeof(int);
	
	sprintf(p, "NHi");
	p += 3;
	memcpy(p, &num_cals, sizeof(int));
	p += sizeof(int);
	
	sprintf(p, "NMi");
	p += 3;
	memcpy(p, &distance, sizeof(int));
	p += sizeof(int);

	// create an alignment and insert it into the list
	alignment = alignment_new();
	alignment_init_single_end(strdup(header_match), read_match, quality_match, 
				  cal->strand, cal->chromosome_id - 1, pos,
				  cigar_code_get_string(cigar_code), cigar_code_get_num_ops(cigar_code), 
				  norm_score * 254, 1, (num_cals > 1),
				  optional_fields_length, optional_fields, 0, alignment);
	array_list_insert(alignment, alignment_list);
      } else {
	// free memory
	if (read_match) free(read_match);
	if (quality_match) free(quality_match);
      }
      
      // free memory
      cigar_code_free(cigar_code);
    }
    
    // free the cal list, and update the mapping list with the alignment list
    free(header_match);
    array_list_free(cal_list, (void *) cal_free);
    mapping_batch->mapping_lists[read_index] = alignment_list;
  }

  // free
  sw_multi_output_free(output);

  // go to the next stage
  return POST_PAIR_STAGE;

  //  printf("END: apply_sw, (%d Smith-Waterman, %d valids)\n", total, valids);
}

//--------------------------------------------------------------------------------------

sw_output_t *sw_output_new(int strand, size_t chrom, size_t ref_start, size_t ref_len,
			   size_t mref_len, size_t mquery_start, size_t mref_start,
			   float score, float norm_score, char* mquery, char* mref) {

  sw_output_t *p = (sw_output_t *) calloc(1, sizeof(sw_output_t));

  p->strand = strand;
  p->chromosome = chrom;
  p->ref_start = ref_start;
  p->ref_len = ref_len;
  p->mref_len = mref_len;
  p->mquery_start = mquery_start;
  p->mref_start = mref_start;
  p->score = score;
  p->norm_score = norm_score;
  p->mquery = strdup(mquery);
  p->mref = strdup(mref);
  //p->mquery = NULL;
  //p->mref = NULL;

  return p;
}

//--------------------------------------------------------------------------------------

void sw_output_free(sw_output_t *p) {
  if (p == NULL) return;

  if (p->mquery != NULL) free(p->mquery);
  if (p->mref != NULL) free(p->mref);

  free(p);
}

//--------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

