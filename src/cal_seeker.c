#include "cal_seeker.h"

//------------------------------------------------------------------------------------
// functions to handle gaps between mapped seeds
//    - fill_gaps
//    - merge_seed_regions
//------------------------------------------------------------------------------------

#define NONE_POS   0
#define BEGIN_POS  1
#define END_POS    2

typedef struct sw_prepare {
  int left_flank;
  int right_flank;  
  char *query;
  char *ref;
  seed_region_t *seed_region;
  cal_t *cal;
} sw_prepare_t;

sw_prepare_t *sw_prepare_new(char *query, char *ref, int left_flank, int right_flank) {
  sw_prepare_t *p = (sw_prepare_t *) malloc(sizeof(sw_prepare_t));
  p->query = query;
  p->ref = ref;
  p->left_flank = left_flank;
  p->right_flank = right_flank;
  p->seed_region = NULL;
  p->cal = NULL;
  return p;
}

void sw_prepare_free(sw_prepare_t *p) {
  if (p) free(p);
}

//------------------------------------------------------------------------------------

void fill_gaps(mapping_batch_t *mapping_batch, sw_optarg_t *sw_optarg, 
	       genome_t *genome, int min_gap, int min_distance) {

  int sw_count = 0;

  fastq_read_t *fq_read;
  array_list_t *fq_batch = mapping_batch->fq_batch;

  size_t read_index, read_len;

  cal_t *cal;
  array_list_t *cal_list = NULL;
  size_t num_cals, num_targets = mapping_batch->num_targets;

  char *revcomp_seq = NULL;

  seed_region_t *s, *prev_s, *new_s;
  linked_list_iterator_t* itr;

  cigar_code_t *cigar_code;

  size_t gap_read_start, gap_read_end, gap_read_len;
  size_t gap_genome_start, gap_genome_end, gap_genome_len;

  int left_flank, right_flank;
  sw_prepare_t *sw_prepare;
  array_list_t *sw_prepare_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  //  LOG_DEBUG("\n\n P R E   -   P R O C E S S\n");

  // initialize query and reference sequences to Smith-Waterman
  for (size_t i = 0; i < num_targets; i++) {

    read_index = mapping_batch->targets[i];
    fq_read = (fastq_read_t *) array_list_get(read_index, fq_batch);

    cal_list = mapping_batch->mapping_lists[read_index];
    num_cals = array_list_size(cal_list);
    
    if (num_cals <= 0) continue;

    read_len = fq_read->length;

    // processing each CAL from this read
    for(size_t j = 0; j < num_cals; j++) {

      // get cal and read index
      cal = array_list_get(j, cal_list);
      LOG_DEBUG_F("CAL #%i of %i (strand %i), sr_list size = %i\n", j, num_cals, cal->strand, cal->sr_list->size);

      prev_s = NULL;
      itr = linked_list_iterator_new(cal->sr_list);
      s = (seed_region_t *) linked_list_iterator_curr(itr);
      while (s != NULL) {
	LOG_DEBUG_F("\tseed: [%i|%i - %i|%i]\n", 
		    s->genome_start, s->read_start, s->read_end, s->genome_end);

	// set the cigar for the current region
	gap_read_len = s->read_end - s->read_start + 1;
	cigar_code = cigar_code_new();
	cigar_code_append_op(cigar_op_new(gap_read_len, 'M'), cigar_code);
	s->info = (void *) cigar_code;

	cigar_code = NULL;
	sw_prepare = NULL;

	if ((prev_s == NULL && s->read_start != 0) || (prev_s != NULL)) {
	  mapping_batch->num_gaps++;
	  if (prev_s == NULL) {
	    // gap at the first position
	    gap_read_start = 0;
	    gap_read_end = s->read_start - 1;

	    gap_genome_start = s->genome_start - s->read_start;
	    gap_genome_end = s->genome_start - 1;

	    gap_read_len = gap_read_end - gap_read_start + 1;
	    gap_genome_len = gap_genome_end - gap_genome_start + 1;

	    if (gap_read_len > min_gap) {
	      // the gap is too big, may be there's another CAL to cover it
	      cigar_code = cigar_code_new();
	      cigar_code_append_op(cigar_op_new(gap_read_len, 'H'), cigar_code);	      
	    } else {
	      left_flank = 0;
	      right_flank = 4;
	    }
	  } else {
	    assert(prev_s->read_end < s->read_start);

	    // gap in a middle position
	    gap_read_start = prev_s->read_end + 1;
	    gap_read_end = s->read_start - 1;

	    gap_genome_start = prev_s->genome_end + 1;
	    gap_genome_end = s->genome_start - 1;

	    gap_read_len = gap_read_end - gap_read_start + 1;
	    gap_genome_len = gap_genome_end - gap_genome_start + 1;

	    left_flank = 2;
	    right_flank = 2;
	  }

	  if (!cigar_code) {
	    // we have to try to fill this gap and get a cigar
	    if (gap_read_len == gap_genome_len) {

	      //    1) first, for from  begin -> end, and begin <- end
	      int distance, first = -1, last = -1;
	      char *query;
	      char *ref = (char *) malloc((gap_genome_len + 1) * sizeof(char));;
	      genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1, 
						&gap_genome_start, &gap_genome_end, genome);
	      // handle strand -
	      if (cal->strand) {
		if (revcomp_seq == NULL) {
		  revcomp_seq = strdup(fq_read->sequence);
		  seq_reverse_complementary(revcomp_seq, read_len);
		}
		query = &revcomp_seq[gap_read_start];
	      } else {
		query = &fq_read->sequence[gap_read_start];
	      }
	      
	      //LOG_DEBUG_F("query: %s\n", query);
	      //LOG_DEBUG_F("ref  : %s\n", ref);
	      distance = 0;
	      for (int k = 0; k < gap_read_len; k++) {
		if (query[k] != ref[k]) {
		  distance++;
		  if (first == -1) first = k;
		  last = k;
		}
	      }
	      LOG_DEBUG_F("dist.: %i of %i (first = %i, last = %i)\n", distance, gap_read_len, first, last);
	      if (distance <= min_distance) {
		cigar_code = cigar_code_new();
		cigar_code_append_op(cigar_op_new(gap_read_len, 'M'), cigar_code);
		cigar_code_inc_distance(distance, cigar_code);
	      }
	      
	      // free memory
	      free(ref);
	    }
	    if (!cigar_code) {
	      //    2) second, prepare SW to run

	      // get query sequence, revcomp if necessary
	      size_t read_start = gap_read_start - left_flank;
	      size_t read_end = gap_read_end + right_flank;
	      int gap_read_len_ex = read_end - read_start + 1;
	      char *query = (char *) malloc((gap_read_len_ex + 1) * sizeof(char));
	      // handle strand -
	      if (cal->strand) {
		if (revcomp_seq == NULL) {
		  revcomp_seq = strdup(fq_read->sequence);
		  seq_reverse_complementary(revcomp_seq, read_len);
		}
		memcpy(query, &revcomp_seq[read_start], gap_read_len_ex);
	      } else {
		memcpy(query, &fq_read->sequence[read_start], gap_read_len_ex);
	      }
	      query[gap_read_len_ex] = '\0';
	      
	      // get ref. sequence
	      size_t genome_start = gap_genome_start - left_flank;
	      size_t genome_end = gap_genome_end + right_flank;
	      int gap_genome_len_ex = genome_end - genome_start + 1;
	      char *ref = (char *) malloc((gap_genome_len_ex + 1) * sizeof(char));;
	      genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1, 
						&genome_start, &genome_end, genome);	      
	      ref[gap_genome_len_ex] = '\0';

	      sw_prepare = sw_prepare_new(query, ref, left_flank, right_flank);
	      array_list_insert(sw_prepare, sw_prepare_list);
	      
	      // increase counter
	      sw_count++;	  
	    }
	  }
	  
	  // insert gap in the list
	  new_s = seed_region_new(gap_read_start, gap_read_end, gap_genome_start, gap_genome_end, 0);
	  new_s->info = (void *) cigar_code;
	  linked_list_iterator_insert(new_s, itr);

	  if (sw_prepare) {
	    sw_prepare->seed_region = new_s;
	  }
	}

	// continue loop...
	prev_s = s;
	linked_list_iterator_next(itr);
	s = linked_list_iterator_curr(itr);
      }

      // check for a gap at the last position
      sw_prepare = NULL;
      if (prev_s != NULL && prev_s->read_end < read_len - 1) { 
	cigar_code = NULL;
	mapping_batch->num_gaps++;
	//	mapping_batch->num_sws++;
	//	mapping_batch->num_ext_sws++;

	// gap at the last position
	gap_read_start = prev_s->read_end + 1;
	gap_read_end = read_len - 1;
	gap_read_len = gap_read_end - gap_read_start + 1;
	
	gap_genome_len = gap_read_len;
	gap_genome_start = prev_s->genome_end + 1;
	gap_genome_end = gap_genome_start + gap_genome_len - 1;

	LOG_DEBUG_F("\t\tgap_read_len = %i, gap_genome_len = %i\n", gap_read_len, gap_genome_len);
	LOG_DEBUG_F("\t\t%i : [%lu|%lu - %lu|%lu]\n", 
		    sw_count, gap_genome_start, gap_read_start, gap_read_end, gap_genome_end);

	if (gap_read_len > min_gap) {
	  // the gap is too big, may be there's another CAL to cover it
	  cigar_code = cigar_code_new();
	  cigar_code_append_op(cigar_op_new(gap_read_len, 'H'), cigar_code);	      
	} else {
	  // we have to try to fill this gap and get a cigar
	  
	  //    1) first, for from  begin -> end, and begin <- end
	  int distance, first = -1, last = -1;
	  char *query;
	  char *ref = (char *) malloc((gap_genome_len + 1) * sizeof(char));;
	  genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1, 
					    &gap_genome_start, &gap_genome_end, genome);
	  // handle strand -
	  if (cal->strand) {
	    if (revcomp_seq == NULL) {
	      revcomp_seq = strdup(fq_read->sequence);
	      seq_reverse_complementary(revcomp_seq, read_len);
	    }
	    query = &revcomp_seq[gap_read_start];
	  } else {
	    query = &fq_read->sequence[gap_read_start];
	  }
	  
	  //	  LOG_DEBUG_F("query: %s\n", query);
	  //	  LOG_DEBUG_F("ref  : %s\n", ref);
	  distance = 0;
	  for (int k = 0; k < gap_read_len; k++) {
	    if (query[k] != ref[k]) {
	      distance++;
	      if (first == -1) first = k;
	      last = k;
	    }
	  }
	  LOG_DEBUG_F("dist.: %i of %i (first = %i, last = %i)\n", distance, gap_read_len, first, last);
	  if (distance <= min_distance) {
	    cigar_code = cigar_code_new();
	    cigar_code_append_op(cigar_op_new(gap_read_len, 'M'), cigar_code);
	    cigar_code_inc_distance(distance, cigar_code);
	  } else {
	    //    2) second, prepare SW to run

	    left_flank = 4;
	    right_flank = 0;
	    
	    // get query sequence, revcomp if necessary
	    size_t read_start = gap_read_start - left_flank;
	    size_t read_end = gap_read_end + right_flank;
	    int gap_read_len_ex = read_end - read_start + 1;
	    char *query = (char *) malloc((gap_read_len_ex + 1) * sizeof(char));
	    // handle strand -
	    if (cal->strand) {
	      if (revcomp_seq == NULL) {
		revcomp_seq = strdup(fq_read->sequence);
		seq_reverse_complementary(revcomp_seq, read_len);
	      }
	      memcpy(query, &revcomp_seq[read_start], gap_read_len_ex);
	    } else {
	      memcpy(query, &fq_read->sequence[read_start], gap_read_len_ex);
	    }
	    query[gap_read_len_ex] = '\0';
	    
	    // get ref. sequence
	    size_t genome_start = gap_genome_start - left_flank;
	    size_t genome_end = gap_genome_end + right_flank;
	    int gap_genome_len_ex = genome_end - genome_start + 1;
	    char *ref = (char *) malloc((gap_genome_len_ex + 1) * sizeof(char));;
	    genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1, 
					      &genome_start, &genome_end, genome);
	    query[gap_genome_len_ex] = '\0';

	    sw_prepare = sw_prepare_new(query, ref, left_flank, right_flank);
	    array_list_insert(sw_prepare, sw_prepare_list);
	    
	    // increase counter
	    sw_count++;	  
	  }
	}
	
	// insert gap in the list
	new_s = seed_region_new(gap_read_start, gap_read_end, gap_genome_start, gap_genome_end, 0);
	new_s->info = (void *) cigar_code;
	linked_list_insert_last(new_s, cal->sr_list);

	if (sw_prepare) {
	  sw_prepare->seed_region = new_s;
	}
      }
      linked_list_iterator_free(itr);      
    }

    // free memory
    if (revcomp_seq) {
      free(revcomp_seq);
      revcomp_seq = NULL;
    }
  }

  //  LOG_DEBUG_F("R U N   S W (sw_count = %i, sw_prepare_list size = %i)\n", sw_count, array_list_size(sw_prepare_list));
  assert(sw_count == array_list_size(sw_prepare_list));

  char *q[sw_count], *r[sw_count];
  for (int i = 0; i < sw_count; i++) {
    sw_prepare = array_list_get(i, sw_prepare_list);
    q[i] = sw_prepare->query;
    r[i] = sw_prepare->ref;
  }
  sw_multi_output_t *output = sw_multi_output_new(sw_count);

  // run Smith-Waterman
  smith_waterman_mqmr(q, r, sw_count, sw_optarg, 1, output);
  
  //  LOG_DEBUG("\n\n P O S T   -   P R O C E S S\n");
  cigar_op_t* cigar_op;
  cigar_code_t *cigar_c;
  for (int i = 0; i < sw_count; i++) {
    sw_prepare = array_list_get(i, sw_prepare_list);
    s = sw_prepare->seed_region;
    int distance;
    int read_gap_len = s->read_end - s->read_start + 1;
    int genome_gap_len = s->genome_end - s->genome_start + 1;

    int read_gap_len_ex = read_gap_len_ex + sw_prepare->left_flank + sw_prepare->right_flank;
    int genome_gap_len_ex = genome_gap_len_ex + sw_prepare->left_flank + sw_prepare->right_flank;

    LOG_DEBUG_F("\tgap (read %lu-%lu, genome %lu-%lu) = (%i, %i)\n", 
		s->read_start, s->read_end, s->genome_start, s->genome_end,
		read_gap_len, genome_gap_len);
    LOG_DEBUG_F("\tflanks (left, right) = (%i, %i)\n", sw_prepare->left_flank, sw_prepare->right_flank);
    LOG_DEBUG_F("\tquery : %s\n", sw_prepare->query);
    LOG_DEBUG_F("\tref   : %s\n", sw_prepare->ref);
    LOG_DEBUG_F("\tmquery: %s (start %i)\n", output->query_map_p[i], output->query_start_p[i]);
    LOG_DEBUG_F("\tmref  : %s (start %i)\n", output->ref_map_p[i], output->ref_start_p[i]);

    cigar_code_t *cigar_c = generate_cigar_code(output->query_map_p[i], output->ref_map_p[i],
						strlen(output->query_map_p[i]), output->ref_start_p[i],
						read_gap_len, &distance);
    LOG_DEBUG_F("\tscore : %0.2f, cigar: %s (distance = %i)\n", 
		output->score_p[i], new_cigar_code_string(cigar_c), distance);

    assert(output->query_start_p[i] == 0);
    //    assert(output->ref_start_p[i] == 0);

    if (cigar_code_get_num_ops(cigar_c) > 2) {
      if (sw_prepare->left_flank > 0) {
	cigar_op = cigar_code_get_op(0, cigar_c);
	assert(cigar_op->number >= sw_prepare->left_flank && cigar_op->name == 'M');
	cigar_op->number -= sw_prepare->left_flank;
      }
      if (sw_prepare->right_flank > 0) {
	cigar_op = cigar_code_get_last_op(cigar_c);
	assert(cigar_op->number >= sw_prepare->right_flank && cigar_op->name == 'M');
	cigar_op->number -= sw_prepare->right_flank;
      }
      init_cigar_string(cigar_c);
      LOG_DEBUG_F("\tnew cigar: %s\n", new_cigar_code_string(cigar_c));
    } else {
      assert(cigar_code_get_num_ops(cigar_c) == 1);
      if (sw_prepare->right_flank > 0) {
	cigar_op = cigar_code_get_last_op(cigar_c);
	assert(cigar_op->number >= sw_prepare->right_flank && cigar_op->name == 'M');
	cigar_op->number -= (sw_prepare->left_flank + sw_prepare->right_flank);
	if (cigar_op->number > read_gap_len) {
	  cigar_code_append_op(cigar_op_new(cigar_op->number - read_gap_len, 'D'), cigar_c);
	} else if (cigar_op->number < read_gap_len) {
	  cigar_code_append_op(cigar_op_new(read_gap_len - cigar_op->number, 'I'), cigar_c);
	} else{
	  init_cigar_string(cigar_c);
	}
	//	LOG_DEBUG_F("\tnew cigar: %s\n", new_cigar_code_string(cigar_c));
      }
    }

    // and now set the cigar for this gap
    s->info = (void *) cigar_c;

    // free
    sw_prepare_free(sw_prepare);
  }
  
  // debugging....
  for (size_t i = 0; i < num_targets; i++) {
    read_index = mapping_batch->targets[i];
    fq_read = (fastq_read_t *) array_list_get(read_index, fq_batch);

    LOG_DEBUG_F("Read %s\n", fq_read->id);

    cal_list = mapping_batch->mapping_lists[read_index];
    num_cals = array_list_size(cal_list);
    
    if (num_cals <= 0) continue;

    // processing each CAL from this read
    for(size_t j = 0; j < num_cals; j++) {

      // get cal and read index
      cal = array_list_get(j, cal_list);
      LOG_DEBUG_F("\tCAL #%i of %i (strand %i), sr_list size = %i\n", j, num_cals, cal->strand, cal->sr_list->size);
      itr = linked_list_iterator_new(cal->sr_list);
      s = (seed_region_t *) linked_list_iterator_curr(itr);
      while (s != NULL) {
	LOG_DEBUG_F("\t\t%s (dist. %i)\t[%i|%i - %i|%i]\n", 
		    (s->info ? new_cigar_code_string((cigar_code_t *) s->info) : ">>>>>> gap"),
		    (s->info ? ((cigar_code_t *) s->info)->distance : -1),
		    s->genome_start, s->read_start, s->read_end, s->genome_end);
	linked_list_iterator_next(itr);
	s = linked_list_iterator_curr(itr);
      }
      linked_list_iterator_free(itr);
    }
  }


  // free memory
  sw_multi_output_free(output);
  array_list_free(sw_prepare_list, (void *) NULL);
}

//------------------------------------------------------------------------------------

void fill_end_gaps(mapping_batch_t *mapping_batch, sw_optarg_t *sw_optarg, 
		   genome_t *genome, int min_H, int min_distance) {

  int sw_count = 0;

  fastq_read_t *fq_read;
  array_list_t *fq_batch = mapping_batch->fq_batch;

  size_t read_index, read_len;

  cal_t *cal;
  array_list_t *cal_list = NULL;
  size_t num_cals, num_targets = mapping_batch->num_targets;

  char *seq, *revcomp_seq = NULL;

  seed_region_t *s;

  cigar_op_t *cigar_op;
  cigar_code_t *cigar_code;

  size_t start, end;
  size_t gap_read_start, gap_read_end, gap_read_len;
  size_t gap_genome_start, gap_genome_end, gap_genome_len;

  int first, last, mode, distance, flank = 5;
  sw_prepare_t *sw_prepare;
  array_list_t *sw_prepare_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  // initialize query and reference sequences to Smith-Waterman
  for (size_t i = 0; i < num_targets; i++) {

    read_index = mapping_batch->targets[i];
    fq_read = (fastq_read_t *) array_list_get(read_index, fq_batch);

    cal_list = mapping_batch->mapping_lists[read_index];
    num_cals = array_list_size(cal_list);
    
    if (num_cals <= 0) continue;

    read_len = fq_read->length;
    revcomp_seq = NULL;

    // processing each CAL from this read
    for(size_t j = 0; j < num_cals; j++) {

      // get cal and read index
      cal = array_list_get(j, cal_list);
      if (cal->sr_list->size == 0) continue;

      sw_prepare = NULL;
      s = (seed_region_t *) linked_list_get_first(cal->sr_list);
      cigar_code = (cigar_code_t *) s->info;
      LOG_DEBUG_F("CAL #%i of %i (strand %i), sr_list size = %i, cigar = %s\n", 
		  j, num_cals, cal->strand, cal->sr_list->size, new_cigar_code_string(cigar_code));
      
      for (int k = 0; k < 2; k++) {
	mode = NONE_POS;
	if (k == 0) {
	  if ((cigar_op = cigar_code_get_op(0, cigar_code)) &&
	      cigar_op->name == 'H' && cigar_op->number > min_H) {
	    LOG_DEBUG_F("%i%c\n", cigar_op->number, cigar_op->name);
	    mode = BEGIN_POS;
	    gap_read_start = 0;
	    gap_read_end = cigar_op->number - 1;
	    gap_genome_start = s->genome_start;
	    gap_genome_end = gap_genome_start + cigar_op->number - 1;
	  }
	} else {
	  if ((cigar_op = cigar_code_get_last_op(cigar_code)) &&
	      cigar_op->name == 'H' && cigar_op->number > min_H) {
	    LOG_DEBUG_F("%i%c\n", cigar_op->number, cigar_op->name);
	    sw_count++;
	    mode = END_POS;
	  }
	}
	    
	if (mode == NONE_POS || mode == END_POS) continue;

	// get query sequence, revcomp if necessary
	if (cal->strand) {
	  if (revcomp_seq == NULL) {
	    revcomp_seq = strdup(fq_read->sequence);
	    seq_reverse_complementary(revcomp_seq, read_len);
	  }
	  seq = revcomp_seq;
	} else {
	  seq = fq_read->sequence;
	}

	gap_read_len = gap_read_end - gap_read_start;
	/*	
	char *query = (char *) malloc((gap_read_len + 1) * sizeof(char));
	memcpy(query, seq, gap_len);
	query[gap_read_len] = '\0';
	*/

	// get ref. sequence
	start = gap_genome_start;
	end = gap_genome_end;
	gap_genome_len = end - start + 1;
	char *ref = (char *) malloc((gap_genome_len + 1) * sizeof(char));
	genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1, 
					  &start, &end, genome);
	ref[gap_genome_len] = '\0';
	
	first = -1; 
	last = -1;
	distance = 0;
	for (int k = 0; k < gap_read_len; k++) {
	  if (seq[k] != ref[k]) {
	    distance++;
	    if (first == -1) first = k;
	    last = k;
	  }
	}
	//	LOG_DEBUG_F("query: %s\n", seq);
	//	LOG_DEBUG_F("ref. : %s\n", ref);
	LOG_DEBUG_F("distance = %i: first = %i, last = %i\n", distance, first, last);
	if (distance < min_distance) {
	  cigar_op->name = 'M';
	  continue;
	}

	// we must run the SW algorithm
	

	//	sw_prepare = sw_prepare_new(0, 0, 0, 0);
	//	sw_prepare_sequences( cal, genome, sw_prepare);
	//	array_list_insert(sw_prepare, sw_prepare_list);
	//	sw_count++;
      }
      /*

      //      sw_prepare = sw_prepare_new(query, ref, left_flank, right_flank);
      */
    }
  }
  LOG_DEBUG_F("sw_count = %i\n", sw_count);


  // debugging....
  for (size_t i = 0; i < num_targets; i++) {
    read_index = mapping_batch->targets[i];
    fq_read = (fastq_read_t *) array_list_get(read_index, fq_batch);

    LOG_DEBUG_F("Read %s\n", fq_read->id);
    
    cal_list = mapping_batch->mapping_lists[read_index];
    num_cals = array_list_size(cal_list);
    
    if (num_cals <= 0) continue;

    for(size_t j = 0; j < num_cals; j++) {

      // get cal and read index
      cal = array_list_get(j, cal_list);
      if (cal->sr_list->size == 0) continue;

      sw_prepare = NULL;
      s = (seed_region_t *) linked_list_get_first(cal->sr_list);
      cigar_code = (cigar_code_t *) s->info;
      LOG_DEBUG_F("CAL #%i of %i (strand %i), sr_list size = %i, cigar = %s\n", 
		  j, num_cals, cal->strand, cal->sr_list->size, new_cigar_code_string(cigar_code));
    }
  }
}

//------------------------------------------------------------------------------------

void merge_seed_regions(mapping_batch_t *mapping_batch) {
  linked_list_item_t *list_item;
  cal_t *cal;
  seed_region_t *s, *s_first;
  cigar_code_t *cigar_code, *cigar_code_prev;
  cigar_op_t *cigar_op, *cigar_op_prev;
  int num_ops;
  int op;
  array_list_t *cals_list;
  fastq_read_t *fq_read;
  linked_list_iterator_t itr;

  register size_t num_targets = mapping_batch->num_targets;
  register size_t num_cals;
  register int i;
  register size_t t;

  for (t = 0; t < num_targets; t++) {
    cals_list = mapping_batch->mapping_lists[mapping_batch->targets[t]];
    num_cals = array_list_size(cals_list);
    fq_read = array_list_get(mapping_batch->targets[t], mapping_batch->fq_batch);

    LOG_DEBUG_F("Read %s\n", fq_read->id);

    for (i = 0; i < num_cals; i++) {
      cal = array_list_get(i, cals_list);
      LOG_DEBUG_F("\tCAL %i:\n", i);
      linked_list_iterator_init(cal->sr_list, &itr);

      s_first = linked_list_iterator_curr(&itr);      
      if (!s_first) { LOG_DEBUG("\t\tLINKED LIST EMPTY"); continue; }

      cigar_code_prev = (cigar_code_t *)s_first->info;
      LOG_DEBUG_F("\t\tpart. cigar: %s\n", new_cigar_code_string(cigar_code_prev));

      s = linked_list_iterator_next(&itr);
      while (s) {
	//LOG_DEBUG_F("\t\tItem [%lu|%i - %i|%lu]: \n", s->genome_start, s->read_start, s->read_end, s->genome_end);
	cigar_code = (cigar_code_t *)s->info;
	if (cigar_code) { //TODO: delete
	  LOG_DEBUG_F("\t\tpart. cigar: %s\n", new_cigar_code_string(cigar_code));
	  num_ops = array_list_size(cigar_code->ops);
	  for (op = 0, cigar_op = array_list_get(op, cigar_code->ops); 
	       op < num_ops;
	       op++, cigar_op = array_list_get(op, cigar_code->ops)) {
	    cigar_code_append_op(cigar_op, cigar_code_prev);	    
	  }
	  cigar_code_prev->distance += cigar_code->distance;
	} else {
	  LOG_DEBUG("\t\tpart. cigar: NULL <<<<<<<");
	} 

	s_first->read_end = s->read_end;
	s_first->genome_end = s->genome_end;

	linked_list_iterator_remove(&itr);
	s = linked_list_iterator_curr(&itr);
      }

      LOG_DEBUG_F("\tFusion Result %i\n", linked_list_size(cal->sr_list));

      linked_list_iterator_init(cal->sr_list, &itr);
      s = linked_list_iterator_curr(&itr);
      while (s) {
	cigar_code = (cigar_code_t *)s->info;	
	LOG_DEBUG_F("\t\tItem [%lu|%i - %i|%lu]: Distance(%i) %s\n", s->genome_start, s->read_start, s->read_end, s->genome_end, cigar_code->distance, new_cigar_code_string(cigar_code));
	s = linked_list_iterator_next(&itr);
      }
    }
  }
}



int apply_caling_rna(cal_seeker_input_t* input, batch_t *batch) {
  //printf("APPLY CALING ...\n");
  struct timeval start, end;
  double time;
  if (time_on) { start_timer(start); }

  mapping_batch_t *mapping_batch = batch->mapping_batch;
  array_list_t *allocate_cals;
  size_t num_cals, select_cals, total_cals = 0;
  size_t num_batches = 0, num_reads_unmapped = 0, num_without_cals = 0;
  size_t total_reads = 0;
  size_t num_targets, target_pos, total_targets, extra_target_pos;
  fastq_read_t *read;
  genome_t *genome = input->genome;
  size_t *targets_aux;
  size_t min_seeds, max_seeds;
  int seed_size = 16;

  num_targets = mapping_batch->num_targets;
  total_targets = 0;
  extra_target_pos = 0;
  total_reads += num_targets;
  target_pos = 0;

  mapping_batch->extra_stage_do = 1;

  for (size_t i = 0; i < num_targets; i++) {
    allocate_cals = array_list_new(1000, 
				   1.25f, 
				   COLLECTION_MODE_ASYNCHRONIZED);

    read = array_list_get(mapping_batch->targets[i], mapping_batch->fq_batch); 

    //printf("%i mini-mappings\n", array_list_size(mapping_batch->mapping_lists[mapping_batch->targets[i]]));
    max_seeds = (read->length / 15)*2 + 10;
    num_cals = bwt_generate_cal_list_linked_list(mapping_batch->mapping_lists[mapping_batch->targets[i]], 
						 input->cal_optarg,
						 &min_seeds, &max_seeds,
						 genome->num_chromosomes + 1,
						 allocate_cals, read->length);

    //printf("\t Target %i: %i\n", mapping_batch->targets[i], num_cals);
    array_list_free(mapping_batch->mapping_lists[mapping_batch->targets[i]], region_bwt_free);
    mapping_batch->mapping_lists[mapping_batch->targets[i]] = allocate_cals;

    if (num_cals > MAX_RNA_CALS) {
      if (!mapping_batch->extra_stage_do) {
	mapping_batch->extra_targets[extra_target_pos] = mapping_batch->targets[i];
	mapping_batch->extra_stage_id[extra_target_pos++] = MAX_CALS;
	array_list_clear(mapping_batch->mapping_lists[mapping_batch->targets[i]], cal_free);	
      } else {
	select_cals = num_cals - MAX_RNA_CALS;
	for(size_t j = num_cals - 1; j >= MAX_RNA_CALS; j--) {
	  cal_free(array_list_remove_at(j, mapping_batch->mapping_lists[mapping_batch->targets[i]]));
	}
	mapping_batch->targets[target_pos++] = mapping_batch->targets[i];
      }
    }else if (num_cals > 0) {
      mapping_batch->targets[target_pos++] = mapping_batch->targets[i];
    }else {
      if (!mapping_batch->extra_stage_do) {
	mapping_batch->extra_targets[extra_target_pos] = mapping_batch->targets[i];
	mapping_batch->extra_stage_id[extra_target_pos++] = NO_CALS;
      }
    }
  }
  
  mapping_batch->num_targets = target_pos;


  if (time_on) { stop_timer(start, end, time); timing_add(time, CAL_SEEKER, timing); }

  //printf("APPLY CAL SEEKER DONE!\n");

  return RNA_STAGE;




  // this code does not offer great improvements !!!
  // should we comment it ?
  if (mapping_batch->extra_stage_do) {
    //printf("Go to original targets & Fusion...\n");
    targets_aux = mapping_batch->targets;
    mapping_batch->targets = mapping_batch->extra_targets;
    mapping_batch->extra_targets = targets_aux;
    mapping_batch->num_targets = mapping_batch->num_extra_targets;
   
    for (size_t i = 0; i < target_pos; i++) {
      mapping_batch->targets[mapping_batch->num_targets++] = mapping_batch->extra_targets[i];
    }
    //mapping_batch->num_targets += target_pos;
    mapping_batch->num_extra_targets = 0;
    mapping_batch->extra_stage_do = 0;
    /*printf("Total Targets %i\n", mapping_batch->num_targets);
    printf("\t--->TARGETS: ", mapping_batch->extra_stage);
    for (int i = 0; i < mapping_batch->num_targets; i++){
      printf("%i,", mapping_batch->targets[i]);
    }
    printf("\n");
    */
  }else if (extra_target_pos) {
    targets_aux = mapping_batch->targets;
    mapping_batch->targets = mapping_batch->extra_targets;
    mapping_batch->extra_targets = targets_aux;
    
    mapping_batch->num_targets = extra_target_pos;
    mapping_batch->num_extra_targets = target_pos;
    //printf("Go back to stage Seeding. Extra Targets = %i, Targets = %i\n", mapping_batch->num_extra_targets, mapping_batch->num_targets );
    mapping_batch->extra_stage_do = 1;

    return SEEDING_STAGE;
  } 

  //return RNA_PREPROCESS_STAGE;

}

//====================================================================================
// apply_caling
//====================================================================================
int apply_caling(cal_seeker_input_t* input, batch_t *batch) {
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  array_list_t *list = NULL;
  size_t read_index, num_cals, min_seeds, max_seeds;
  int min_limit;

  cal_t *cal;
  array_list_t *cal_list;

  fastq_read_t *read;

  size_t num_chromosomes = input->genome->num_chromosomes + 1;
  size_t num_targets = mapping_batch->num_targets;
  size_t *targets = mapping_batch->targets;
  size_t new_num_targets = 0;
  max_seeds = input->cal_optarg->num_seeds;
  
  //  size_t *new_targets = (size_t *) calloc(num_targets, sizeof(size_t));
  
  // set to zero
  mapping_batch->num_to_do = 0;

  for (size_t i = 0; i < num_targets; i++) {

    read_index = targets[i];
    read = array_list_get(read_index, mapping_batch->fq_batch); 

    // for debugging
    //    LOG_DEBUG_F("%s\n", ((fastq_read_t *) array_list_get(read_index, mapping_batch->fq_batch))->id);
    
    if (!list) {
      list = array_list_new(1000, 
			    1.25f, 
			    COLLECTION_MODE_ASYNCHRONIZED);
    }

    // optimized version
    max_seeds = (read->length / 15)*2 + 10;
    num_cals = bwt_generate_cal_list_linked_list(mapping_batch->mapping_lists[read_index], 
						 input->cal_optarg,
						 &min_seeds, &max_seeds,
						 num_chromosomes,
						 list,
						 read->length);
    /*
    // for debugging
    LOG_DEBUG_F("num. cals = %i, min. seeds = %i, max. seeds = %i\n", num_cals, min_seeds, max_seeds);

    for (size_t j = 0; j < num_cals; j++) {
      cal = array_list_get(j, list);
      LOG_DEBUG_F("\tchr: %i, strand: %i, start: %lu, end: %lu, num_seeds = %i, num. regions = %lu\n", 
		  cal->chromosome_id, cal->strand, cal->start, cal->end, cal->num_seeds, cal->sr_list->size);
    }
    */
    //    printf("min_seeds = %i, max_seeds = %i, min_limit = %i, num_cals = %i\n", 
    //	   min_seeds, max_seeds, min_limit, array_list_size(list));

    // filter incoherent CALs
    int founds[num_cals], found = 0;
    for (size_t j = 0; j < num_cals; j++) {
      cal = array_list_get(j, list);
      if (cal->sr_list->size > 0) {
	int start = 0;
	for (linked_list_item_t *list_item = cal->sr_list->first; list_item != NULL; list_item = list_item->next) {
	  seed_region_t *s = list_item->item;
	  if (start > s->read_start) {
	    found++;
	    founds[j] = 1;
	  }
	  start = s->read_end + 1;
	}
      }
    }
    if (found) {
      cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
      for (size_t j = 0; j < num_cals; j++) {
	if (!founds[j]) {
	  cal = array_list_get(j, list);
	  array_list_insert(cal, cal_list);
	  array_list_set(j, NULL, list);
	}
      }
      array_list_free(list, (void *) cal_free);
      num_cals = array_list_size(cal_list);
      list = cal_list;
    }

    // filter CALs by the number of seeds
    int min_limit = input->cal_optarg->min_num_seeds_in_cal;
    if (min_limit < 0) min_limit = max_seeds;
    
    if (min_seeds == max_seeds || min_limit <= min_seeds) {
      cal_list = list;
      list = NULL;
    } else {
      cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
      for (size_t j = 0; j < num_cals; j++) {
	cal = array_list_get(j, list);
	if (cal->num_seeds >= min_limit) {
	  array_list_insert(cal, cal_list);
	  array_list_set(j, NULL, list);
	}
      }
      array_list_clear(list, (void *) cal_free);
      num_cals = array_list_size(cal_list);
    }

    if (num_cals > MAX_CALS) {
      for (size_t j = num_cals - 1; j >= MAX_CALS; j--) {
	cal = (cal_t *) array_list_remove_at(j, cal_list);
	cal_free(cal);
      }
      num_cals = array_list_size(cal_list);
    }
    
    if (num_cals > 0 && num_cals <= MAX_CALS) {
      array_list_set_flag(2, cal_list);
      targets[new_num_targets++] = read_index;

      int count1 = 0, count2 = 0;
      // count number of sw to do

      // method #1
      //      printf("method #1\n");
      seed_region_t *s, *prev_s;
      linked_list_iterator_t* itr;
      for (size_t j = 0; j < num_cals; j++) {
	prev_s = NULL;
	cal = array_list_get(j, cal_list);
	itr = linked_list_iterator_new(cal->sr_list);
	s = (seed_region_t *) linked_list_iterator_curr(itr);
	while (s != NULL) {
	  if ((prev_s == NULL && s->read_start != 0) || (prev_s != NULL)) {
	    //	    printf("\t\t\tcase 1\n");
	    count1++;
	  }
	  prev_s = s;
	  linked_list_iterator_next(itr);
	  s = linked_list_iterator_curr(itr);
	}
	if (prev_s != NULL && prev_s->read_end < read->length - 1) { 
	  count1++;
	  //	  printf("\t\t\tcase 2 (%i < %i)\n", prev_s->read_end, read->length - 1);
	}
	linked_list_iterator_free(itr);
      }
      /*
      // method #2
      printf("method #2\n");
      for (size_t j = 0; j < num_cals; j++) {
	cal = array_list_get(j, cal_list);
	printf("\t: %i\n", j);
	if (cal->sr_list->size > 0) {
	  int start = 0;
	  for (linked_list_item_t *list_item = cal->sr_list->first; list_item != NULL; list_item = list_item->next) {
	    seed_region_t *s = list_item->item;
	    printf("\t\t[%i|%i - %i|%i]\n", s->genome_start, s->read_start, s->read_end, s->genome_end);
	    if (s->read_start != start) {
	      count2++;
	    }
	    start = s->read_end + 1;
	  }
	  if (start < read->length) { 
	    count2++;
	  }
	}
      }
      printf("count #1 = %i, count #2 = %i\n", count1, count2);
      assert(count1 == count2);
*/
      mapping_batch->num_to_do += count1;

      // we have to free the region list
      array_list_free(mapping_batch->mapping_lists[read_index], (void *) region_bwt_free);
      mapping_batch->mapping_lists[read_index] = cal_list;
    } else {
      array_list_set_flag(0, mapping_batch->mapping_lists[read_index]);
      // we have to free the region list
      array_list_clear(mapping_batch->mapping_lists[read_index], (void *) region_bwt_free);
      if (cal_list) array_list_free(cal_list, (void *) cal_free);
      if (list) array_list_clear(list, (void *) cal_free);
    }

    /*    
    cal_list = list;
    list = NULL;
    array_list_set_flag(2, cal_list);
    //    mapping_batch->num_to_do += num_cals;
    targets[new_num_targets++] = read_index;
    
    // we have to free the region list
    array_list_free(mapping_batch->mapping_lists[read_index], (void *) region_bwt_free);
    mapping_batch->mapping_lists[read_index] = cal_list;
    */
    /*
    // filter CALs by the number of seeds
    int min_limit = input->cal_optarg->min_num_seeds_in_cal;
    if (min_limit < 0) min_limit = max_seeds;

    printf("min_seeds = %i, max_seeds = %i, min_limit = %i, num_cals = %i\n", 
	   min_seeds, max_seeds, min_limit, array_list_size(list));
    
    if (min_seeds == max_seeds || min_limit <= min_seeds) {
      cal_list = list;
      list = NULL;
    } else {
      cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
      for (size_t j = 0; j < num_cals; j++) {
	cal = array_list_get(j, list);
	if (cal->num_seeds >= min_limit) {
	  array_list_insert(cal, cal_list);
	  array_list_set(j, NULL, list);
	}
      }
      array_list_clear(list, (void *) cal_free);
      num_cals = array_list_size(cal_list);
      printf("************, num_cals = %i\n", num_cals);
    }

    if (num_cals > MAX_CALS) {
      for (size_t j = num_cals - 1; j >= MAX_CALS; j--) {
	cal = (cal_t *) array_list_remove_at(j, cal_list);
	cal_free(cal);
      }
      num_cals = array_list_size(cal_list);
    }

    if (num_cals > 0 && num_cals <= MAX_CALS) {
      array_list_set_flag(2, cal_list);
      mapping_batch->num_to_do += num_cals;
      targets[new_num_targets++] = read_index;
      
      // we have to free the region list
      array_list_free(mapping_batch->mapping_lists[read_index], (void *) region_bwt_free);
      mapping_batch->mapping_lists[read_index] = cal_list;
    } else {
      array_list_set_flag(0, mapping_batch->mapping_lists[read_index]);
      // we have to free the region list
      array_list_clear(mapping_batch->mapping_lists[read_index], (void *) region_bwt_free);
      if (cal_list) array_list_free(cal_list, (void *) cal_free);
      if (list) array_list_clear(list, (void *) cal_free);
    }
    */
  } // end for 0 ... num_targets

  // update batch
  mapping_batch->num_targets = new_num_targets;

  //  LOG_DEBUG_F("num. SW to do: %i\n", 	mapping_batch->num_to_do);

  //  exit(-1);

  // free memory
  if (list) array_list_free(list, NULL);

  if (batch->pair_input->pair_mng->pair_mode != SINGLE_END_MODE) {
    return PRE_PAIR_STAGE;
  } else if (batch->mapping_batch->num_targets > 0) {
    return SW_STAGE;
  }
  
  return DNA_POST_PAIR_STAGE;
}

//------------------------------------------------------------------------------------
// cal_seeker_input functions: init
//------------------------------------------------------------------------------------

void cal_seeker_input_init(list_t *regions_list, cal_optarg_t *cal_optarg, 
			   list_t* write_list, unsigned int write_size, 
			   list_t *sw_list, list_t *pair_list, 
			   genome_t *genome, cal_seeker_input_t *input) {
  input->regions_list = regions_list;
  input->cal_optarg = cal_optarg;
  input->batch_size = write_size;
  input->sw_list = sw_list;
  input->pair_list = pair_list;
  input->write_list = write_list;
  input->genome = genome;
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
