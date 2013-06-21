#include "cal_seeker.h"

//------------------------------------------------------------------------------------
// cal_seeker_input functions: init
//------------------------------------------------------------------------------------

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

  num_targets = mapping_batch->num_targets;
  total_targets = 0;
  extra_target_pos = 0;
  total_reads += num_targets;
  target_pos = 0;

  /*
  printf("(STEP %i)TARGETS: ", mapping_batch->extra_stage);
  for (int i = 0; i < num_targets; i++){
    printf("%i,", mapping_batch->targets[i]);
  }
  printf("\n");
  */

  mapping_batch->extra_stage_do = 1;

  for (size_t i = 0; i < num_targets; i++) {
    allocate_cals = array_list_new(1000, 
				   1.25f, 
				   COLLECTION_MODE_ASYNCHRONIZED);

    read = array_list_get(mapping_batch->targets[i], mapping_batch->fq_batch); 

    //printf("%i mini-mappings\n", array_list_size(mapping_batch->mapping_lists[mapping_batch->targets[i]]));
    max_seeds = (read->length / 15)*2 + 10;
    num_cals = bwt_generate_cal_rna_list_linked_list(mapping_batch->mapping_lists[mapping_batch->targets[i]], 
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

  return RNA_PREPROCESS_STAGE;

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

  return RNA_PREPROCESS_STAGE;

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

    printf("min_seeds = %i, max_seeds = %i, min_limit = %i, num_cals = %i\n", 
	   min_seeds, max_seeds, min_limit, array_list_size(list));
    */
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
  
  return POST_PAIR_STAGE;
}

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
