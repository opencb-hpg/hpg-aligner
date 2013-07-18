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

    //printf("%i mini-mappings\n", array_list_size(mapping_batch->mapping_lists[mapping_batch->targets[i]]));

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
    //    LOG_DEBUG_F("%s\n", ((fastq_read_t *) array_list_get(read_index, mapping_batch->fq_batch))->id);
    
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

    /*
    // for debugging
    LOG_DEBUG_F("num. cals = %i, min. seeds = %i, max. seeds = %i\n", num_cals, min_seeds, max_seeds);
    for (size_t j = 0; j < num_cals; j++) {
      cal = array_list_get(j, list);
      LOG_DEBUG_F("\tchr: %i, strand: %i, start: %lu, end: %lu, num. seeds = %lu, flank: (start, end) = (%lu, %lu)\n", 
		  cal->chromosome_id, cal->strand, cal->start, cal->end, cal->num_seeds, cal->flank_start, cal->flank_end);
    }
    */

    // filter CALs by the number of seeds
    //int min_limit = 0; //input->cal_optarg->min_num_seeds_in_cal;
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
	  //	  LOG_DEBUG_F("\t\tchr: %i, strand: %i, start: %lu, end: %lu, num. seeds = %lu (min. limit %i seeds), flank: (start, end) = (%lu, %lu)\n", 
	  //		      cal->chromosome_id, cal->strand, cal->start, cal->end, cal->num_seeds, min_limit, cal->flank_start, cal->flank_end);
	  
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

//====================================================================================
// apply_caling bs
//====================================================================================
int apply_caling_bs(cal_seeker_input_t* input, batch_t *batch) {

  //printf("APPLY CALLING BS...\n");

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
  
  // new variables for bisulfite
  size_t num_targets2 = mapping_batch->num_targets2;
  size_t *targets2 = mapping_batch->targets2;
  size_t new_num_targets2 = 0;
  array_list_t *list2 = NULL;
  size_t num_cals2;
  array_list_t *cal_list2;
  size_t read_index2;

  // set to zero
  mapping_batch->num_to_do = 0;
  mapping_batch->num_to_do2 = 0;


  ////////////////////////////////
  /*
  size_t reads_mapp  = 0;
  size_t reads_mapp2 = 0;
  size_t reads_no_mapp  = 0;
  size_t reads_no_mapp2 = 0;
  size_t reads_discard  = 0;
  size_t reads_discard2 = 0;

  size_t reads_cals = 0;
  size_t reads_cals_acum = 0;
  */
  ////////////////////////////////


  //printf("targets 1 %lu\ntargets 2 %lu\n", num_targets, num_targets2);

  //printf("----->primera lista\n");
  for (size_t i = 0; i < num_targets; i++) {

    read_index = targets[i];

    // for debugging
    //    LOG_DEBUG_F("%s\n", ((fastq_read_t *) array_list_get(read_index, mapping_batch->fq_batch))->id);
    
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
    //printf("read %lu\tcals1 = %lu\n", read_index, num_cals);

    /*
    // for debugging
    LOG_DEBUG_F("num. cals = %i, min. seeds = %i, max. seeds = %i\n", num_cals, min_seeds, max_seeds);
    for (size_t j = 0; j < num_cals; j++) {
      cal = array_list_get(j, list);
      LOG_DEBUG_F("\tchr: %i, strand: %i, start: %lu, end: %lu, num. seeds = %lu, flank: (start, end) = (%lu, %lu)\n", 
		  cal->chromosome_id, cal->strand, cal->start, cal->end, cal->num_seeds, cal->flank_start, cal->flank_end);
    }
    */

    // filter CALs by the number of seeds
    int min_limit = input->cal_optarg->min_num_seeds_in_cal;

    if (min_limit < 0) min_limit = max_seeds;

    if (min_seeds == max_seeds || min_limit <= min_seeds) {
      cal_list = list;
      list = NULL;

    } else {
      cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

      /////////////////
      /*
      if (num_cals > 0)
	reads_cals++;
      */
      /////////////////

      for (size_t j = 0; j < num_cals; j++) {
	cal = array_list_get(j, list);
	//filter cals with few seeds
	if (cal->num_seeds >= min_limit) {

	  /////////////////
	  //reads_cals_acum++;
	  /////////////////


	  //	  LOG_DEBUG_F("\t\tchr: %i, strand: %i, start: %lu, end: %lu, num. seeds = %lu (min. limit %i seeds), flank: (start, end) = (%lu, %lu)\n", 
	  //		      cal->chromosome_id, cal->strand, cal->start, cal->end, cal->num_seeds, min_limit, cal->flank_start, cal->flank_end);
	  
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
      /////////////
      //reads_mapp++;
      /////////////

      array_list_set_flag(2, cal_list);
      mapping_batch->num_to_do += num_cals;
      targets[new_num_targets++] = read_index;
      
      // we have to free the region list
      array_list_free(mapping_batch->mapping_lists[read_index], (void *) region_bwt_free);
      mapping_batch->mapping_lists[read_index] = cal_list;
    } else {
      /////////////
      //reads_discard++;
      /////////////
      array_list_set_flag(0, mapping_batch->mapping_lists[read_index]);
      // we have to free the region list
      array_list_clear(mapping_batch->mapping_lists[read_index], (void *) region_bwt_free);
      if (cal_list) array_list_free(cal_list, (void *) cal_free);
      if (list) array_list_clear(list, (void *) cal_free);
    }

    //printf("read %lu\tcals  = %lu\n", read_index, num_cals);
  } // end for 0 ... num_seqs


  //printf("----->segunda lista\n");
  // add for bisulfite
  for (size_t i = 0; i < num_targets2; i++) {

    read_index2 = targets2[i];
    //printf("\nread %lu\n", read_index2);

    // for debugging
    //    LOG_DEBUG_F("%s\n", ((fastq_read_t *) array_list_get(read_index, mapping_batch->fq_batch))->id);
    
    if (!list2) {
      list2 = array_list_new(1000, 
			     1.25f, 
			     COLLECTION_MODE_ASYNCHRONIZED);
    }

    // optimized version
    num_cals2 = bwt_generate_cal_list_linkedlist(mapping_batch->mapping_lists2[read_index2], 
						 input->cal_optarg,
						 &min_seeds, &max_seeds,
						 num_chromosomes,
						 list2);
    //printf("read %lu\tcals2 = %lu\n", read_index2, num_cals2);

    /*
    // for debugging
    LOG_DEBUG_F("num. cals = %i, min. seeds = %i, max. seeds = %i\n", num_cals2, min_seeds, max_seeds);
    for (size_t j = 0; j < num_cals2; j++) {
      cal = array_list_get(j, list);
      LOG_DEBUG_F("\tchr: %i, strand: %i, start: %lu, end: %lu, num. seeds = %lu, flank: (start, end) = (%lu, %lu)\n", 
		  cal->chromosome_id, cal->strand, cal->start, cal->end, cal->num_seeds, cal->flank_start, cal->flank_end);
    }
    */

    // filter CALs by the number of seeds
    int min_limit = input->cal_optarg->min_num_seeds_in_cal;

    if (min_limit < 0) min_limit = max_seeds;

    if (min_seeds == max_seeds || min_limit <= min_seeds) {
      //printf("read %lu\tborrar listas\n", read_index2);
      cal_list2 = list2;
      list2 = NULL;
      //printf("read %lu\tborrar listas\n", read_index2);

    } else {
      //printf("read %lu\t1recortar lista\n", read_index2);
      cal_list2 = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

      /////////////////
      /*
      if (num_cals2 > 0)
	reads_cals++;
      */
      /////////////////

      for (size_t j = 0; j < num_cals2; j++) {
	cal = array_list_get(j, list2);
	//printf("read %lu\tnum_seeds %lu\tmin_limit %i\n", read_index2, cal->num_seeds, min_limit);
	//filter cals with few seeds
	if (cal->num_seeds >= min_limit) {

	  /////////////////
	  //reads_cals_acum++;
	  /////////////////

	  //printf("keep cal %lu, with %lu seed (of %lu needed)\n", j, cal->num_seeds, min_limit);
	  //	  LOG_DEBUG_F("\t\tchr: %i, strand: %i, start: %lu, end: %lu, num. seeds = %lu (min. limit %i seeds), flank: (start, end) = (%lu, %lu)\n", 
	  //		      cal->chromosome_id, cal->strand, cal->start, cal->end, cal->num_seeds, min_limit, cal->flank_start, cal->flank_end);
	  
	  //printf("pre  callist2 size %lu\n", array_list_size(cal_list2));
	  array_list_insert(cal, cal_list2);
	  //printf("post callist2 size %lu\n", array_list_size(cal_list2));
	  array_list_set(j, NULL, list2);
	  //printf("list2 size %lu\n", array_list_size(list2));
	}
      }
      array_list_clear(list2, (void *) cal_free);
      num_cals2 = array_list_size(cal_list2);
      //printf("read %lu\t2recortar lista\n", read_index2);
    }

    if (num_cals2 > MAX_CALS) {
      //printf("read %lu\tacortar lista\n", read_index2);
      for (size_t j = num_cals2 - 1; j >= MAX_CALS; j--) {
	cal = (cal_t *) array_list_remove_at(j, cal_list2);
	cal_free(cal);
      }
      num_cals2 = array_list_size(cal_list2);
      //printf("read %lu\tacortar lista\n", read_index2);
    }

    //printf("read %lu\tcals2 = %lu\n", read_index2, num_cals2);

    if (num_cals2 > 0 && num_cals2 <= MAX_CALS) {
      /////////////
      //reads_mapp2++;
      /////////////

      //printf("read %lu\tactualizar lista\n", read_index2);
      array_list_set_flag(2, cal_list2);
      mapping_batch->num_to_do2 += num_cals2;
      targets2[new_num_targets2++] = read_index2;
      
      // we have to free the region list
      array_list_free(mapping_batch->mapping_lists2[read_index2], (void *) region_bwt_free);
      mapping_batch->mapping_lists2[read_index2] = cal_list2;
      //printf("read %lu\tactualizar lista\n", read_index2);
    } else {
      /////////////
      //reads_discard2++;
      /////////////

      //printf("read %lu\tdescartar lista\n", read_index2);
      array_list_set_flag(0, mapping_batch->mapping_lists2[read_index2]);
      // we have to free the region list
      array_list_clear(mapping_batch->mapping_lists2[read_index2], (void *) region_bwt_free);
      if (cal_list2) array_list_free(cal_list2, (void *) cal_free);
      if (list2) array_list_clear(list2, (void *) cal_free);
      //printf("read %lu\tdescartar lista\n", read_index2);
    }

    //printf("read %lu\tcals2 = %lu\n", read_index2, num_cals2);

    //printf("read %lu\tcals1 = %lu\tsize mapps  %i\n", read_index, num_cals, array_list_size(mapping_batch->mapping_lists[read_index]));
    //printf("read %lu\tcals2 = %lu\tsize mapps2 %lu\n", read_index2, num_cals2, array_list_size(mapping_batch->mapping_lists2[read_index2]));

    
  } // end for 0 ... num_seqs
  // end add for bisulfite

  /*
  printf("CAL_seek1   \t%3lu\thave CAL (to SW)\t%3lu\thave no CALs     \t%3lu\n", 
	 num_targets, reads_mapp, reads_discard);
  printf("CAL_seek2   \t%3lu\thave CAL (to SW)\t%3lu\thave no CALs     \t%3lu\n", 
	 num_targets2, reads_mapp2, reads_discard2);
  */
  /*
  printf("2 reads CAL   \t%3lu\ttotal CALs      \t%3lu\tCALs promedio    \t%6.2f\n", 
	 reads_cals, reads_cals_acum, 1.0 * reads_cals_acum / reads_cals);
  */

  //printf("new targets1 = %lu, new targets2 = %lu\n", new_num_targets, new_num_targets2);

  // update batch
  mapping_batch->num_targets = new_num_targets;

  // free memory
  if (list) array_list_free(list, NULL);

  // added for bisulfite
  // update batch
  mapping_batch->num_targets2 = new_num_targets2;

  // free memory
  if (list2) array_list_free(list2, NULL);
  // end added for bisulfite

  //printf("End cal stage\nSW_STAGE = %lu\n", SW_STAGE);

  //return SW_STAGE;

  if (batch->pair_input->pair_mng->pair_mode != SINGLE_END_MODE) {
    //printf("return PRE_PAIR_STAGE\n");
    return PRE_PAIR_STAGE;
  } else if (batch->mapping_batch->num_targets > 0 || batch->mapping_batch->num_targets2 > 0) {
    //printf("return SW_STAGE\n");
    return SW_STAGE;
  }

  //printf("return POST_PAIR_STAGE\n");
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
