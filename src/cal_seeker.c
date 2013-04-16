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

    num_cals = bwt_generate_cal_list_rna_linked_list(mapping_batch->mapping_lists[mapping_batch->targets[i]], 
						     input->cal_optarg, 
						     allocate_cals, 
						     read->length, genome->num_chromosomes + 1);  
						     /*
    num_cals = bwt_generate_cal_list_rna_linkedlist(mapping_batch->mapping_lists[mapping_batch->targets[i]], 
						     input->cal_optarg, 
						     allocate_cals, 
						     read->length, genome->num_chromosomes + 1);  
						     */

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

  size_t num_chromosomes = input->genome->num_chromosomes + 1;
  size_t num_targets = mapping_batch->num_targets;
  size_t *targets = mapping_batch->targets;
  size_t new_num_targets = 0;
  //  size_t *new_targets = (size_t *) calloc(num_targets, sizeof(size_t));
  
  // set to zero
  mapping_batch->num_to_do = 0;

  for (size_t i = 0; i < num_targets; i++) {

    read_index = targets[i];

    // for debugging
    LOG_DEBUG_F("%s\n", ((fastq_read_t *) array_list_get(read_index, mapping_batch->fq_batch))->id);
    
    if (!list) {
      list = array_list_new(1000, 
			    1.25f, 
			    COLLECTION_MODE_ASYNCHRONIZED);
    }
    // optimized version
    num_cals = bwt_generate_cal_list_linkedlist(mapping_batch->mapping_lists[read_index], 
						input->cal_optarg,
						&min_seeds, &max_seeds,
						num_chromosomes,
						list);

    // for debugging
    LOG_DEBUG_F("num. cals = %i, min. seeds = %i, max. seeds = %i\n", num_cals, min_seeds, max_seeds);
    for (size_t j = 0; j < num_cals; j++) {
      cal = array_list_get(j, list);
      LOG_DEBUG_F("\tchr: %i, strand: %i, start: %lu, end: %lu, num. seeds = %lu, flank: (start, end) = (%lu, %lu)\n", 
		  cal->chromosome_id, cal->strand, cal->start, cal->end, cal->num_seeds, cal->flank_start, cal->flank_end);
    }

    // filter CALs by the number of seeds
    if (min_seeds == max_seeds) {
      cal_list = list;
      list = NULL;
    } else {
      int seed_filter = ceil(0.3f * max_seeds);
      cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
      for (size_t j = 0; j < num_cals; j++) {
	cal = array_list_get(j, list);
	if (cal->num_seeds >= seed_filter) {
	  LOG_DEBUG_F("\t\tchr: %i, strand: %i, start: %lu, end: %lu, num. seeds = %lu (filter %i), flank: (start, end) = (%lu, %lu)\n", 
		      cal->chromosome_id, cal->strand, cal->start, cal->end, cal->num_seeds, seed_filter, cal->flank_start, cal->flank_end);
	  
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
  } // end for 0 ... num_seqs

  // update batch
  mapping_batch->num_targets = new_num_targets;

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
			   genome_t *genome, cal_seeker_input_t *input){
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
