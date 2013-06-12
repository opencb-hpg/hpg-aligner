#include "bwt_server.h"

//====================================================================================
// bwt_server_input functions: init
//====================================================================================

void bwt_server_input_init(list_t* read_list_p, unsigned int batch_size, bwt_optarg_t *bwt_optarg_p, 
			   bwt_index_t *bwt_index_p, list_t* write_list_p, unsigned int write_size, 
			   list_t* unmapped_read_list_p, bwt_server_input_t* input_p) {
  input_p->read_list_p = read_list_p;
  input_p->batch_size = batch_size;
  input_p->bwt_optarg_p = bwt_optarg_p;
  input_p->write_list_p = write_list_p;
  input_p->write_size = write_size;
  input_p->bwt_index_p = bwt_index_p;
  input_p->unmapped_read_list_p = unmapped_read_list_p;

}

//====================================================================================
// apply_bwt
//====================================================================================

int apply_bwt(bwt_server_input_t* input, batch_t *batch) {
  // run bwt _by_filter
  //printf("APPLY BWT SERVER...\n");
  struct timeval start, end;
  double time;

  if (time_on) { start_timer(start); }

  mapping_batch_t *mapping_batch = batch->mapping_batch;
  
  bwt_map_inexact_array_list_by_filter(mapping_batch->fq_batch, input->bwt_optarg_p,
				       input->bwt_index_p, 
				       mapping_batch->mapping_lists,
				       &mapping_batch->num_targets, mapping_batch->targets);
  
  size_t num_mapped_reads = array_list_size(mapping_batch->fq_batch) - mapping_batch->num_targets;
  mapping_batch->num_to_do = num_mapped_reads;

  if (time_on) { stop_timer(start, end, time); timing_add(time, BWT_SERVER, timing); }
  
  //printf("APPLY BWT SERVER DONE!\n");

  if (batch->mapping_batch->num_targets > 0) {
    //TODO: DELETE
    //printf("Web have targets\n");
    for (int i = 0; i < batch->mapping_batch->num_targets; i++) {
      batch->mapping_batch->bwt_mappings[batch->mapping_batch->targets[i]] = 1;
    }
    return SEEDING_STAGE;
  }
  //printf("Reads are mapped\n");
  return POST_PAIR_STAGE;
}

//====================================================================================
// apply_bwt
//====================================================================================

int apply_bwt_bs(bwt_server_input_t* input, batch_t *batch) {
  // run bwt _by_filter
  //printf("APPLY BWT_BS SERVER...\n");
  struct timeval start, end;
  double time;

  if (time_on) { start_timer(start); }

  mapping_batch_t *mapping_batch = batch->mapping_batch;
  array_list_t *array_tmp;
  
  //copy the batch reads
  size_t num_reads = array_list_size(mapping_batch->fq_batch);

  mapping_batch->CT_fq_batch     = array_list_new(num_reads + 2, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  mapping_batch->CT_rev_fq_batch = array_list_new(num_reads + 2, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  mapping_batch->GA_fq_batch     = array_list_new(num_reads + 2, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  mapping_batch->GA_rev_fq_batch = array_list_new(num_reads + 2, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  cpy_array_bs(mapping_batch->fq_batch, 
	       mapping_batch->CT_fq_batch, mapping_batch->CT_rev_fq_batch, 
	       mapping_batch->GA_fq_batch, mapping_batch->GA_rev_fq_batch);

  //printf("******end copy array\n");

  //transform the batch reads
  replace_array(mapping_batch->GA_fq_batch, ACT);
  //printf("******end replace G->A\n");
  rev_comp_array(mapping_batch->GA_rev_fq_batch, mapping_batch->GA_fq_batch);
  //printf("******end rev_comp G->A\n");
  replace_array(mapping_batch->CT_fq_batch, AGT);
  //printf("******end replace C->T\n");
  rev_comp_array(mapping_batch->CT_rev_fq_batch, mapping_batch->CT_fq_batch);
  //printf("******end rev_comp C->T\n");

  /*
  // input->bwt_index_p  = index ACT
  // input->bwt_index2_p = index AGT
  // check index nucleotides and reads conversion
  printf("index_p  %s\n", input->bwt_index_p->nucleotides);
  printf("index2_p %s\n", input->bwt_index2_p->nucleotides);
  fastq_read_t* fq_read_src;
  for (size_t i = 0; i < num_reads; i++) {
    printf("\nread = %lu\n", i);
    fq_read_src  = (fastq_read_t *) array_list_get(i, mapping_batch->fq_batch);
    printf("orig   = %s\n", fq_read_src->sequence);
    fq_read_src  = (fastq_read_t *) array_list_get(i, mapping_batch->CT_fq_batch);
    printf("CT     = %s\n", fq_read_src->sequence);
    fq_read_src  = (fastq_read_t *) array_list_get(i, mapping_batch->CT_rev_fq_batch);
    printf("CT_rev = %s\n", fq_read_src->sequence);
    fq_read_src  = (fastq_read_t *) array_list_get(i, mapping_batch->GA_fq_batch);
    printf("GA     = %s\n", fq_read_src->sequence);
    fq_read_src  = (fastq_read_t *) array_list_get(i, mapping_batch->GA_rev_fq_batch);
    printf("GA_rev = %s\n", fq_read_src->sequence);
  }
  */

  /*
  // original call
  bwt_map_inexact_array_list_by_filter_bs(mapping_batch->fq_batch, input->bwt_optarg_p,
					  input->bwt_index_p, 
					  mapping_batch->mapping_lists,
					  &mapping_batch->num_targets, mapping_batch->targets);
  */

  array_tmp = mapping_batch->fq_batch;
  // mapping against input->bwt_index_p
  mapping_batch->fq_batch = mapping_batch->GA_fq_batch;
  //mapping_batch->fq_batch = mapping_batch->CT_rev_fq_batch;
  // mapping against input->bwt_index2_p
  //mapping_batch->fq_batch = mapping_batch->CT_fq_batch;
  //mapping_batch->fq_batch = mapping_batch->GA_rev_fq_batch;

  /*
  fastq_read_t* fq_read_src;
  for (size_t i = 0; i < num_reads; i++) {
    printf("\nread = %lu\n", i);
    fq_read_src  = (fastq_read_t *) array_list_get(i, mapping_batch->fq_batch);
    printf("orig   = %s\n", fq_read_src->sequence);
  }
  */
  bwt_map_inexact_array_list_by_filter_bs(mapping_batch->fq_batch, input->bwt_optarg_p,
					  input->bwt_index_p, 
					  mapping_batch->mapping_lists,
					  &mapping_batch->num_targets, mapping_batch->targets);

  mapping_batch->fq_batch = array_tmp;

  /*
  mapping_batch->fq_batch = mapping_batch->GA_fq_batch;
  bwt_map_inexact_array_list_by_filter_bs(mapping_batch->fq_batch, input->bwt_optarg_p,
					  input->bwt_index2_p, 
					  mapping_batch->mapping_lists,
					  &mapping_batch->num_targets, mapping_batch->targets);

  */
  //make the four searches
  /*
  printf("******search 1 on index \t%s\n", input->bwt_index_p->nucleotides);

  bwt_map_inexact_array_list_by_filter_bs(mapping_batch->CT_fq_batch, input->bwt_optarg_p,
					  input->bwt_index_p, 
					  mapping_batch->mapping_lists,
					  &mapping_batch->num_targets, mapping_batch->targets);

  printf("******end search 1 on index \t%s\n\n", input->bwt_index_p->nucleotides);
  */
  /*
  printf("******search 2 on index \t%s\n", input->bwt_index2_p->nucleotides);

  bwt_map_inexact_array_list_by_filter_bs(mapping_batch->CT_rev_fq_batch, input->bwt_optarg_p,
					  input->bwt_index2_p, 
					  mapping_batch->mapping_lists,
					  &mapping_batch->num_targets, mapping_batch->targets);

  printf("******end search 2 on index \t%s\n\n", input->bwt_index2_p->nucleotides);
  */
  /*
  printf("******search 3 on index \t%s\n", input->bwt_index2_p->nucleotides);

  bwt_map_inexact_array_list_by_filter_bs(mapping_batch->GA_fq_batch, input->bwt_optarg_p,
					  input->bwt_index2_p, 
					  mapping_batch->mapping_lists,
					  &mapping_batch->num_targets, mapping_batch->targets);

  printf("******end search 3 on index \t%s\n\n", input->bwt_index2_p->nucleotides);
  */
  /*
  printf("******search 4 on index \t%s\n", input->bwt_index_p->nucleotides);

  bwt_map_inexact_array_list_by_filter_bs(mapping_batch->GA_rev_fq_batch, input->bwt_optarg_p,
					  input->bwt_index_p, 
					  mapping_batch->mapping_lists,
					  &mapping_batch->num_targets, mapping_batch->targets);

  printf("******end search 4 on index \t%s\n\n", input->bwt_index_p->nucleotides);
  */

  size_t num_mapped_reads = array_list_size(mapping_batch->fq_batch) - mapping_batch->num_targets;
  mapping_batch->num_to_do = num_mapped_reads;

  if (time_on) { stop_timer(start, end, time); timing_add(time, BWT_SERVER, timing); }
  
  //printf("APPLY BWT SERVER DONE!\n");
  return CONSUMER_STAGE;

  /*
  if (batch->mapping_batch->num_targets > 0) {
    //TODO: DELETE
    //printf("Web have targets\n");
    for (int i = 0; i < batch->mapping_batch->num_targets; i++) {
      batch->mapping_batch->bwt_mappings[batch->mapping_batch->targets[i]] = 1;
    }
    return SEEDING_STAGE;
  }
  //printf("Reads are mapped\n");
  return POST_PAIR_STAGE;
  */
}

//------------------------------------------------------------------------------------

/*
void bwt_server_cpu(bwt_server_input_t* input, pair_mng_t *pair_mng) {    
    LOG_DEBUG_F("bwt_server_cpu(%d): START\n", omp_get_thread_num());
    list_item_t *item = NULL;
    array_list_t *fq_batch;
    mapping_batch_t *mapping_batch;
    //list_t *write_list = input->write_list_p;
    list_t *unmapped_read_list = input->unmapped_read_list_p;
    size_t num_mappings, num_mappings_tot = 0;
    size_t total_reads = 0;
    size_t reads_no_mapped = 0; 
    //write_batch_t* write_batch_p = write_batch_new(write_size, MATCH_FLAG);
    size_t num_batches = 0;

    while ( (item = list_remove_item(input->read_list_p)) != NULL ) {
      //printf("BWT Server extract one item...\n");
	num_batches++;
	
	if (time_on) { timing_start(BWT_SERVER, 0, timing_p); }	
	fq_batch = (array_list_t *) item->data_p;
	mapping_batch = mapping_batch_new(fq_batch, pair_mng); //TODO:¿? set pair_mng from ¿?¿?¿?

	//printf("\tCall function process\n");
	num_mappings = bwt_map_inexact_array_list(fq_batch,
						  input->bwt_optarg_p,
						  input->bwt_index_p, 
						  mapping_batch->mapping_lists,
						  &mapping_batch->num_targets, 
						  mapping_batch->targets);
	//printf("\tEnd call\n");
	num_mappings_tot += num_mappings;
	total_reads += array_list_size(mapping_batch->fq_batch);
	reads_no_mapped += mapping_batch->num_targets;
	//Results
	//printf("Process Batch (bwt_server_cpu): (%d)Mappings - (%d)Unmappings\n", num_mappings, unmapped_batch_p->num_reads);
	//printf("BWT Processing batch finish!reads no mapped\n");
	list_item_free(item);
	item = list_item_new(0, WRITE_ITEM, mapping_batch);

	if (time_on) { timing_stop(BWT_SERVER, 0, timing_p); }
	list_insert_item(item, unmapped_read_list);
    }
    
    //list_decr_writers(write_list);
    list_decr_writers(unmapped_read_list);
    LOG_DEBUG_F("bwt_server_cpu (Total reads process %lu, Reads unmapped %lu): END\n", 
	   total_reads, reads_no_mapped); 
}
*/
/*
void bwt_map_inexact_batch_by_filter(fastq_batch_t *batch,
				     bwt_optarg_t *bwt_optarg, 
				     bwt_index_t *index,
				     size_t num_selected, size_t *selected,
				     size_t num_lists, array_list_t **lists,
				     size_t *num_mapped, size_t* mapped,
				     size_t *num_unmapped, size_t *unmapped);
*/
//------------------------------------------------------------------------------------
//int mapped_by_bwt[100];
