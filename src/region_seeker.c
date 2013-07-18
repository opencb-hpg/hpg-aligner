#include "region_seeker.h"


void region_seeker_input_init(list_t *unmapped_read_list, cal_optarg_t *cal_optarg, 
			      bwt_optarg_t *bwt_optarg, bwt_index_t *bwt_index, 
			      list_t* region_list, unsigned int region_threads, 
			      unsigned int gpu_enable, int padding_left, int padding_right,
			      region_seeker_input_t *input_p) {

  input_p->unmapped_read_list_p = unmapped_read_list;
  input_p->cal_optarg_p = cal_optarg;
  input_p->bwt_optarg_p = bwt_optarg;
  input_p->bwt_index_p = bwt_index;
  input_p->region_list_p = region_list;
  input_p->region_threads = region_threads;
  input_p->gpu_enable = gpu_enable;
  input_p->padding_left = padding_left;
  input_p->padding_right = padding_right;
}

//====================================================================================
// apply_seeding
//====================================================================================

int apply_seeding(region_seeker_input_t* input, batch_t *batch) {
  //printf("APPLY SEEDING...\n");
  struct timeval start, end;
  double time;
  if (time_on) { start_timer(start); }

  mapping_batch_t *mapping_batch = batch->mapping_batch;
  array_list_t *list = NULL;
  size_t read_index, num_mappings;

  size_t num_seeds = input->cal_optarg_p->num_seeds;
  size_t seed_size = input->cal_optarg_p->seed_size;
  size_t min_seed_size = input->cal_optarg_p->min_seed_size;
  int padding_left = input->padding_left;
  int padding_right = input->padding_right;

  fastq_read_t *fq_read;
  array_list_t *fq_batch = mapping_batch->fq_batch;

  size_t num_targets = mapping_batch->num_targets;
  size_t *targets = mapping_batch->targets;
  size_t new_num_targets = 0;
  fastq_read_t *read;
  // set to zero
  mapping_batch->num_to_do = 0;
  
  //TODO: omp parallel for !!
  if (batch->mapping_mode == RNA_MODE) {
    for (size_t i = 0; i < num_targets; i++) {
      //printf("Seq (i=%i)(target=%i): %s\n", i, targets[i], read->sequence);
      read = array_list_get(targets[i], mapping_batch->fq_batch);
      num_mappings = bwt_map_exact_seeds_seq(padding_left,
					     padding_right,
					     read->sequence,
					     seed_size,
					     min_seed_size,
					     input->bwt_optarg_p, 
					     input->bwt_index_p, 
					     mapping_batch->mapping_lists[targets[i]],
					     mapping_batch->extra_stage_id[targets[i]]);
      
      //printf("Num mappings %i\n", num_mappings);
      if (num_mappings > 0) {
	array_list_set_flag(2, mapping_batch->mapping_lists[targets[i]]);
	targets[new_num_targets++] = targets[i];
	mapping_batch->num_to_do += num_mappings;
      }
    }
  } else {
    for (size_t i = 0; i < num_targets; i++) {
      read = array_list_get(targets[i], mapping_batch->fq_batch);
      //      printf("region_seeker: seeds for %s\n", read->id);
      num_mappings = bwt_map_exact_seeds_seq_by_num(read->sequence, num_seeds,
						    seed_size, min_seed_size,
						    input->bwt_optarg_p, input->bwt_index_p, 
						    mapping_batch->mapping_lists[targets[i]]);
      if (num_mappings > 0) {
	array_list_set_flag(2, mapping_batch->mapping_lists[targets[i]]);
	targets[new_num_targets++] = targets[i];
	mapping_batch->num_to_do += num_mappings;
      }
    }
  }

  // update batch targets
  mapping_batch->num_targets = new_num_targets;

  if (time_on) { stop_timer(start, end, time); timing_add(time, REGION_SEEKER, timing); }

  //printf("APPLY SEEDING DONE!\n");

  return CAL_STAGE;

}

//------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

int apply_seeding_bs(region_seeker_input_t* input, batch_t *batch) {

  //printf("APPLY SEEDING BS...\n");
  struct timeval start, end;
  double time;
  if (time_on) { start_timer(start); }

  mapping_batch_t *mapping_batch = batch->mapping_batch;
  array_list_t *list = NULL;
  size_t read_index;
  size_t num_mapps1 = 0, num_mapps2 = 0,
    num_mapps3 = 0, num_mapps4 = 0;

  size_t num_seeds = input->cal_optarg_p->num_seeds;
  size_t seed_size = input->cal_optarg_p->seed_size;
  size_t min_seed_size = input->cal_optarg_p->min_seed_size;
  int padding_left = input->padding_left;
  int padding_right = input->padding_right;

  fastq_read_t *fq_read;
  array_list_t *fq_batch = mapping_batch->fq_batch;

  size_t num_targets = mapping_batch->num_targets;
  size_t *targets = mapping_batch->targets;
  size_t *targets2 = mapping_batch->targets2;
  size_t new_num_targets = 0;
  size_t new_num_targets2 = 0;
  fastq_read_t *read;

  // set to zero
  mapping_batch->num_to_do = 0;
  mapping_batch->num_to_do2 = 0;
  

  // create 
  size_t num_reads = array_list_size(mapping_batch->fq_batch);


  /*
  // mostrar las reads
  {
    fastq_read_t* fq_read_src;

    for (size_t i = 0; i < num_targets; i++) {
      fq_read_src  = (fastq_read_t *) array_list_get(targets[i], mapping_batch->fq_batch);
      printf("\nId = %lu\tOrig:   %s\n", i, fq_read_src->sequence);
      fq_read_src  = (fastq_read_t *) array_list_get(targets[i], mapping_batch->CT_fq_batch);
      printf("Id = %lu\tCT:     %s\n", i, fq_read_src->sequence);
      fq_read_src  = (fastq_read_t *) array_list_get(targets[i], mapping_batch->CT_rev_fq_batch);
      printf("Id = %lu\tCT_rev: %s\n", i, fq_read_src->sequence);
      fq_read_src  = (fastq_read_t *) array_list_get(targets[i], mapping_batch->GA_fq_batch);
      printf("Id = %lu\tGA:     %s\n", i, fq_read_src->sequence);
      fq_read_src  = (fastq_read_t *) array_list_get(targets[i], mapping_batch->GA_rev_fq_batch);
      printf("Id = %lu\tGA_rev: %s\n", i, fq_read_src->sequence);
    }
  }
  */


  ////////////////////////////////
  /*
  size_t reads_mapp = 0;
  size_t reads_mapp2 = 0;
  size_t reads_no_mapp = 0;
  size_t reads_no_mapp2 = 0;
  size_t reads_discard = 0;
  */
  ////////////////////////////////


  //TODO: omp parallel for !!
  //if (batch->mapping_mode == BS_MODE) {
  for (size_t i = 0; i < num_targets; i++) {
    // original call

    /*
    read = array_list_get(targets[i], mapping_batch->fq_batch);
    //printf("Seq (i=%i)(target=%i): %s\n", i, targets[i], read->sequence);
    //printf("region_seeker: seeds for %s\n", read->id);
    num_mappings = bwt_map_exact_seeds_seq_by_num_bs(read->sequence, num_seeds,
						     seed_size, min_seed_size,
						     input->bwt_optarg_p, input->bwt_index_p, 
						     mapping_batch->mapping_lists[targets[i]]);
    */

    num_mapps1 = 0;
    num_mapps2 = 0;
    num_mapps3 = 0;
    num_mapps4 = 0;

    read = array_list_get(targets[i], mapping_batch->GA_rev_fq_batch);
    //printf("Seq (i=%i)(target=%i): %s\n", i, targets[i], read->sequence);
    //printf("region_seeker: seeds for %s\n", read->id);
    num_mapps2 = bwt_map_exact_seeds_seq_by_num_bs(read->sequence, num_seeds,
						   seed_size, min_seed_size,
						   input->bwt_optarg_p, input->bwt_index2_p, 
						   mapping_batch->mapping_lists[targets[i]]);
    // transform the reads from the search 2 to the reverse strand
    if (num_mapps2 > 0) {
      //printf("transform maps1\n");
      transform_regions(mapping_batch->mapping_lists[targets[i]]);
    }

    //    printf("---->is null ? %i, size = %i\n", (mapping_batch->mapping_lists[targets[i]] == NULL),
    //	   mapping_batch->mapping_lists[targets[i]]->size);

    read = array_list_get(targets[i], mapping_batch->GA_fq_batch);
    //printf("Seq (i=%i)(target=%i): %s\n", i, targets[i], read->sequence);
    //printf("region_seeker: seeds for %s\n", read->id);
    num_mapps1 = bwt_map_exact_seeds_seq_by_num_bs(read->sequence, num_seeds,
						   seed_size, min_seed_size,
						   input->bwt_optarg_p, input->bwt_index_p, 
						   mapping_batch->mapping_lists[targets[i]]);

    //    printf("---->is null ? %i, size = %i\n", (mapping_batch->mapping_lists[targets[i]] == NULL),
    //	   mapping_batch->mapping_lists[targets[i]]->size);


    //printf("----<is null ? %i, size = %i\n", (mapping_batch->mapping_lists2[targets2[i]] == NULL),
    //	   mapping_batch->mapping_lists2[targets2[i]]->size);

    read = array_list_get(targets[i], mapping_batch->CT_rev_fq_batch);
    //printf("Seq (i=%i)(target=%i): %s\n", i, targets[i], read->sequence);
    //printf("region_seeker: seeds for %s\n", read->id);
    num_mapps4 = bwt_map_exact_seeds_seq_by_num_bs(read->sequence, num_seeds,
						   seed_size, min_seed_size,
						   input->bwt_optarg_p, input->bwt_index_p, 
						   mapping_batch->mapping_lists2[targets[i]]);

    // transform the reads from the search 4 to the reverse strand
    if (num_mapps4 > 0) {
      //printf("transform maps2\n");
      transform_regions(mapping_batch->mapping_lists2[targets[i]]);
    }

    //    printf("----<is null ? %i, size = %i\n", (mapping_batch->mapping_lists2[targets2[i]] == NULL),
    //	   mapping_batch->mapping_lists2[targets2[i]]->size);

    read = array_list_get(targets[i], mapping_batch->CT_fq_batch);
    //printf("Seq (i=%i)(target=%i): %s\n", i, targets[i], read->sequence);
    //printf("region_seeker: seeds for %s\n", read->id);
    num_mapps3 = bwt_map_exact_seeds_seq_by_num_bs(read->sequence, num_seeds,
						   seed_size, min_seed_size,
						   input->bwt_optarg_p, input->bwt_index2_p, 
						   mapping_batch->mapping_lists2[targets[i]]);

    //    printf("----<is null ? %i, size = %i\n", (mapping_batch->mapping_lists2[targets2[i]] == NULL),
    //	   mapping_batch->mapping_lists2[targets2[i]]->size);


    //printf("Num mappings for read %lu\ns1: %lu\ns2: %lu\ns3: %lu\ns4: %lu\n\n", targets[i], num_mapps1, num_mapps2, num_mapps3, num_mapps4);

    if (num_mapps1 > 0) {
      //printf("set flags\n");
      // si hay semillas en la read, se marca como elemento a procesar y se guarda la lista
      array_list_set_flag(2, mapping_batch->mapping_lists[targets[i]]);
      targets[new_num_targets++] = targets[i];
      mapping_batch->num_to_do += num_mapps1;

      //////////////
      //reads_mapp++;
      //////////////
    } else {
      // si no hay semillas en la read, se borra la lista
      array_list_clear(mapping_batch->mapping_lists[targets[i]], NULL);

      //////////////
      //reads_no_mapp++;
      //////////////
    }

    if (num_mapps3 > 0) {
      //printf("set flags\n");
      // si hay semillas en la read, se marca como elemento a procesar y se guarda la lista
      array_list_set_flag(2, mapping_batch->mapping_lists2[targets[i]]);
      targets2[new_num_targets2++] = targets[i];
      mapping_batch->num_to_do2 += num_mapps3;

      //////////////
      //reads_mapp2++;
      //////////////
    } else {
      // si no hay semillas en la read, se borra la lista
      array_list_clear(mapping_batch->mapping_lists2[targets[i]], NULL);

      //////////////
      //reads_no_mapp2++;
      //////////////
    }

  }
  //} // end if MODE_BS

  // update batch targets
  mapping_batch->num_targets = new_num_targets;
  mapping_batch->num_targets2 = new_num_targets2;

  /*
  printf("BWT_exact1  \t%3lu\thave seeds (to CAL)\t%3lu\thave no seed     \t%3lu\n", 
	 num_targets, reads_mapp, reads_no_mapp);
  printf("BWT_exact2  \t%3lu\thave seeds (to CAL)\t%3lu\thave no seed     \t%3lu\n", 
	 num_targets, reads_mapp2, reads_no_mapp2);
  */

  if (time_on) { stop_timer(start, end, time); timing_add(time, REGION_SEEKER, timing); }

  //printf("new targets1 = %lu, new targets2 = %lu\n", new_num_targets, new_num_targets2);
  //printf("APPLY SEEDING DONE!\n");

  return CAL_STAGE;
  //return CONSUMER_STAGE;

}

//------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
