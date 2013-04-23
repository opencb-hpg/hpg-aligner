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
