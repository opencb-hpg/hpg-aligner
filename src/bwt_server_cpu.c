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
  
  // variables for histogram
  float histogram[256];


  if (time_on) { start_timer(start); }

  mapping_batch_t *mapping_batch = batch->mapping_batch;
  // bs variables
  //array_list_t *array_tmp;
  array_list_t 
    *mapping_list1,
    *mapping_list2;

  //copy the batch reads
  size_t num_reads = array_list_size(mapping_batch->fq_batch);
  //printf("reads %lu\n", num_reads);

  //intialize the new mappings and indices
  //printf("Init new variables - mappings\n");
  mapping_batch->mapping_lists2 = (array_list_t **) calloc(num_reads, sizeof(array_list_t*));
  mapping_list1 = array_list_new(500,
				 1.25f,
				 COLLECTION_MODE_ASYNCHRONIZED);
  mapping_list2 = array_list_new(500,
				 1.25f,
				 COLLECTION_MODE_ASYNCHRONIZED);
  for (size_t i = 0; i < num_reads; i++) {
    mapping_batch->mapping_lists2[i] = array_list_new(500,
						      1.25f,
						      COLLECTION_MODE_ASYNCHRONIZED);
  }

  //printf("Init new variables - reads transformations\n");
  mapping_batch->CT_fq_batch     = array_list_new(num_reads + 2, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  mapping_batch->CT_rev_fq_batch = array_list_new(num_reads + 2, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  mapping_batch->GA_fq_batch     = array_list_new(num_reads + 2, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  mapping_batch->GA_rev_fq_batch = array_list_new(num_reads + 2, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  // copy and transform the reads simultaneously
  cpy_transform_array_bs(mapping_batch->fq_batch, 
			 mapping_batch->CT_fq_batch, mapping_batch->CT_rev_fq_batch, 
			 mapping_batch->GA_fq_batch, mapping_batch->GA_rev_fq_batch);

  /*
  // mostrar las reads
  fastq_read_t* fq_read_src;

  for (size_t i = 0; i < num_reads; i++) {
    fq_read_src  = (fastq_read_t *) array_list_get(i, mapping_batch->fq_batch);
    printf("\nId = %lu\tOrig:   %s\n", i, fq_read_src->sequence);
    fq_read_src  = (fastq_read_t *) array_list_get(i, mapping_batch->CT_fq_batch);
    printf("Id = %lu\tCT:     %s\n", i, fq_read_src->sequence);
    fq_read_src  = (fastq_read_t *) array_list_get(i, mapping_batch->CT_rev_fq_batch);
    printf("Id = %lu\tCT_rev: %s\n", i, fq_read_src->sequence);
    fq_read_src  = (fastq_read_t *) array_list_get(i, mapping_batch->GA_fq_batch);
    printf("Id = %lu\tGA:     %s\n", i, fq_read_src->sequence);
    fq_read_src  = (fastq_read_t *) array_list_get(i, mapping_batch->GA_rev_fq_batch);
    printf("Id = %lu\tGA_rev: %s\n", i, fq_read_src->sequence);
  }
  */

// make the four searches
  alignment_t *alignment;
  size_t header_len, num_mapps1 = 0, num_mapps2 = 0, num_mapps3 = 0, num_mapps4 = 0;
  size_t num_threads = input->bwt_optarg_p->num_threads;
  num_reads = array_list_size(mapping_batch->fq_batch);
  size_t chunk = MAX(1, num_reads/(num_threads*10));
  fastq_read_t* fq_read;
  mapping_batch->num_targets = 0;

  for (size_t i = 0; i < num_reads; i++) {
    //printf("\n********** read %lu **********\n", i);

    array_list_set_flag(1, mapping_batch->mapping_lists[i]);
    array_list_set_flag(1, mapping_batch->mapping_lists2[i]);

    /*
    // make histogram of the original read
    fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->GA_fq_batch);
    // initialize the counter
    histogram['A'] = 0.0;
    histogram['C'] = 0.0;
    histogram['G'] = 0.0;
    histogram['T'] = 0.0;

    for (size_t j = 0; j < fq_read->length; j++) {
      histogram[fq_read->sequence[j]]++;
    }

    //printf("A=%.0f\tC=%.0f\tG=%.0f\tT=%.0f\n",
    //	   histogram['A'], histogram['C'],histogram['G'], histogram['T']);
    //printf("A=%.2f\tC=%.2f\tG=%.2f\tT=%.2f\n",
    //	   histogram['A'] / fq_read->length, histogram['C'] / fq_read->length,
    //	   histogram['G'] / fq_read->length, histogram['T'] / fq_read->length);

    if (histogram['C'] / fq_read->length < LIMIT_INF &&
	histogram['T'] / fq_read->length > LIMIT_SUP) {
      //printf("1\n");

      fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->CT_fq_batch);
      //printf("index2 %s\n", input->bwt_index2_p->nucleotides);
      //printf("read3   %s\n", fq_read->sequence);
      num_mapps3 = bwt_map_forward_inexact_seq(fq_read->sequence,
					       input->bwt_optarg_p, input->bwt_index2_p,
					       mapping_batch->mapping_lists2[i]);
      
      fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->CT_rev_fq_batch);
      //printf("index1 %s\n", input->bwt_index_p->nucleotides);
      //printf("read4   %s\n", fq_read->sequence);
      num_mapps4 = bwt_map_forward_inexact_seq(fq_read->sequence,
					       input->bwt_optarg_p, input->bwt_index_p,
					       mapping_list2);

      if (num_mapps3 + num_mapps4 == 0) {
	fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->GA_fq_batch);
	//printf("index1 %s\n", input->bwt_index_p->nucleotides);
	//printf("read1   %s\n", fq_read->sequence);
	num_mapps1 = bwt_map_forward_inexact_seq(fq_read->sequence,
						 input->bwt_optarg_p, input->bwt_index_p,
						 mapping_batch->mapping_lists[i]);

	fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->GA_rev_fq_batch);
	//printf("index2 %s\n", input->bwt_index2_p->nucleotides);
	//printf("read2   %s\n", fq_read->sequence);
	num_mapps2 = bwt_map_forward_inexact_seq(fq_read->sequence,
						 input->bwt_optarg_p, input->bwt_index2_p,
						 mapping_list1);
      }
    }
    else {
      if (histogram['G'] / fq_read->length < LIMIT_INF &&
	  histogram['A'] / fq_read->length > LIMIT_SUP) {
	//printf("2\n");

	fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->GA_fq_batch);
	//printf("index1 %s\n", input->bwt_index_p->nucleotides);
	//printf("read1   %s\n", fq_read->sequence);
	num_mapps1 = bwt_map_forward_inexact_seq(fq_read->sequence,
						 input->bwt_optarg_p, input->bwt_index_p,
						 mapping_batch->mapping_lists[i]);

	fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->GA_rev_fq_batch);
	//printf("index2 %s\n", input->bwt_index2_p->nucleotides);
	//printf("read2   %s\n", fq_read->sequence);
	num_mapps2 = bwt_map_forward_inexact_seq(fq_read->sequence,
						 input->bwt_optarg_p, input->bwt_index2_p,
						 mapping_list1);
	if (num_mapps1 + num_mapps2 == 0) {
	  fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->CT_fq_batch);
	  //printf("index2 %s\n", input->bwt_index2_p->nucleotides);
	  //printf("read3   %s\n", fq_read->sequence);
	  num_mapps3 = bwt_map_forward_inexact_seq(fq_read->sequence,
						   input->bwt_optarg_p, input->bwt_index2_p,
						   mapping_batch->mapping_lists2[i]);

	  fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->CT_rev_fq_batch);
	  //printf("index1 %s\n", input->bwt_index_p->nucleotides);
	  //printf("read4   %s\n", fq_read->sequence);
	  num_mapps4 = bwt_map_forward_inexact_seq(fq_read->sequence,
						   input->bwt_optarg_p, input->bwt_index_p,
						   mapping_list2);
	}
      }
      else {
	//printf("3\n");

	fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->GA_fq_batch);
	//printf("index1 %s\n", input->bwt_index_p->nucleotides);
	//printf("read1   %s\n", fq_read->sequence);
	num_mapps1 = bwt_map_forward_inexact_seq(fq_read->sequence,
						 input->bwt_optarg_p, input->bwt_index_p,
						 mapping_batch->mapping_lists[i]);

	fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->GA_rev_fq_batch);
	//printf("index2 %s\n", input->bwt_index2_p->nucleotides);
	//printf("read2   %s\n", fq_read->sequence);
	num_mapps2 = bwt_map_forward_inexact_seq(fq_read->sequence,
						 input->bwt_optarg_p, input->bwt_index2_p,
						 mapping_list1);

	fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->CT_fq_batch);
	//printf("index2 %s\n", input->bwt_index2_p->nucleotides);
	//printf("read3   %s\n", fq_read->sequence);
	num_mapps3 = bwt_map_forward_inexact_seq(fq_read->sequence,
						 input->bwt_optarg_p, input->bwt_index2_p,
						 mapping_batch->mapping_lists2[i]);

	fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->CT_rev_fq_batch);
	//printf("index1 %s\n", input->bwt_index_p->nucleotides);
	//printf("read4   %s\n", fq_read->sequence);
	num_mapps4 = bwt_map_forward_inexact_seq(fq_read->sequence,
						 input->bwt_optarg_p, input->bwt_index_p,
						 mapping_list2);
      }
    }
    */

    fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->GA_fq_batch);
    //printf("index1 %s\n", input->bwt_index_p->nucleotides);
    //printf("read1   %s\n", fq_read->sequence);
    num_mapps1 = bwt_map_forward_inexact_seq(fq_read->sequence,
					     input->bwt_optarg_p, input->bwt_index_p,
					     mapping_batch->mapping_lists[i]);
    //printf("mapps  %lu\n", num_mapps1);
    //printf("end    %s\n", fq_read->sequence);

    fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->GA_rev_fq_batch);
    //printf("index2 %s\n", input->bwt_index2_p->nucleotides);
    //printf("read2   %s\n", fq_read->sequence);
    num_mapps2 = bwt_map_forward_inexact_seq(fq_read->sequence,
					     input->bwt_optarg_p, input->bwt_index2_p,
					     mapping_list1);
    //printf("mapps  %lu\n", num_mapps2);
    //printf("end    %s\n", fq_read->sequence);

    fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->CT_fq_batch);
    //printf("index2 %s\n", input->bwt_index2_p->nucleotides);
    //printf("read3   %s\n", fq_read->sequence);
    num_mapps3 = bwt_map_forward_inexact_seq(fq_read->sequence,
					     input->bwt_optarg_p, input->bwt_index2_p,
					     mapping_batch->mapping_lists2[i]);
    //printf("mapps  %lu\n", num_mapps3);
    //printf("end    %s\n", fq_read->sequence);

    fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->CT_rev_fq_batch);
    //printf("index1 %s\n", input->bwt_index_p->nucleotides);
    //printf("read4   %s\n", fq_read->sequence);
    num_mapps4 = bwt_map_forward_inexact_seq(fq_read->sequence,
					     input->bwt_optarg_p, input->bwt_index_p,
					     mapping_list2);
    //printf("mapps  %lu\n", num_mapps4);
    //printf("end    %s\n", fq_read->sequence);


    if (num_mapps1 + num_mapps2 + num_mapps3 + num_mapps4 > 0) {
      array_list_set_flag(1, mapping_batch->mapping_lists[i]);
      array_list_set_flag(1, mapping_batch->mapping_lists2[i]);

      // transform and unify the mappings of search 2 in the mappings of search 1
      if (num_mapps2 > 0) {
	transform_mappings(mapping_list1);
	insert_mappings(mapping_batch->mapping_lists[i], mapping_list1);
      }

      // transform and unify the mappings of search 4 in the mappings of search 3
      if (num_mapps4 > 0) {
	transform_mappings(mapping_list2);
	insert_mappings(mapping_batch->mapping_lists2[i], mapping_list2);
      }

      //printf("num_reads = %lu\tnum_mapps1 = %lu\tnum_mapps2 = %lu\tnum_mapps3 = %lu\tnum_mapps4 = %lu\n",
      //     num_reads, num_mapps1, num_mapps2, num_mapps3, num_mapps4);

      num_mapps1 = array_list_size(mapping_batch->mapping_lists[i]);
      num_mapps3 = array_list_size(mapping_batch->mapping_lists2[i]);

      //printf("num_reads = %lu\tnum_mapps1 = %lu\tnum_mapps2 = %lu\tnum_mapps3 = %lu\tnum_mapps4 = %lu\n\n",
      //     num_reads, num_mapps1, num_mapps2, num_mapps3, num_mapps4);

      //printf("Update alignments1\n");
      for (size_t j = 0; j < num_mapps1; j++) {
	alignment = (alignment_t *) array_list_get(j, mapping_batch->mapping_lists[i]);
	header_len = strlen(fq_read->id);
	alignment->query_name = (char *) malloc(sizeof(char) * (header_len + 1));
	
	get_to_first_blank(fq_read->id, header_len, alignment->query_name);
	bwt_cigar_cpy(alignment, fq_read->quality);
	
	// ************************* OPTIONAL FIELDS ***************************
	alignment = add_optional_fields(alignment, num_mapps1);
	// *********************** OPTIONAL FIELDS END *************************
      }

      //printf("Update alignments2\n");
      for (size_t j = 0; j < num_mapps3; j++) {
	alignment = (alignment_t *) array_list_get(j, mapping_batch->mapping_lists2[i]);
	header_len = strlen(fq_read->id);
	alignment->query_name = (char *) malloc(sizeof(char) * (header_len + 1));
	
	get_to_first_blank(fq_read->id, header_len, alignment->query_name);
	bwt_cigar_cpy(alignment, fq_read->quality);
	
	// ************************* OPTIONAL FIELDS ***************************
	alignment = add_optional_fields(alignment, num_mapps3);
	// *********************** OPTIONAL FIELDS END *************************
      }

    } else {
      //imprimir (o guardar) las reads no mapeadas exactas (para depuración)
      //fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->fq_batch);
      //printf("read no mapp %s\n", fq_read->sequence);

      if (array_list_get_flag(mapping_batch->mapping_lists[i]) != 2 && array_list_get_flag(mapping_batch->mapping_lists2[i]) != 2) {
	mapping_batch->targets[(mapping_batch->num_targets)++] = i;
	array_list_set_flag(0, mapping_batch->mapping_lists[i]);
	array_list_set_flag(0, mapping_batch->mapping_lists2[i]);
      } else {
	array_list_set_flag(2, mapping_batch->mapping_lists[i]);
	array_list_set_flag(2, mapping_batch->mapping_lists2[i]);
      }
    }
    array_list_clear(mapping_list1, NULL);
    array_list_clear(mapping_list2, NULL);
  }

  array_list_free(mapping_list1, NULL);
  array_list_free(mapping_list2, NULL);

  //printf("End postprocess\n");

  size_t num_mapped_reads = array_list_size(mapping_batch->fq_batch) - mapping_batch->num_targets;
  mapping_batch->num_to_do = num_mapped_reads;

  //printf("mapped_reads = %lu\treads = %lu\ttargets = %lu\n",
  //	 num_mapped_reads, array_list_size(mapping_batch->fq_batch), mapping_batch->num_targets);
  //printf("reads to process = %lu\n", mapping_batch->num_to_do);


  //revert_mappings_seqs(mapping_batch->mapping_lists,  mapping_batch->fq_batch);
  //revert_mappings_seqs(mapping_batch->mapping_lists2, mapping_batch->fq_batch);

  if (time_on) { stop_timer(start, end, time); timing_add(time, BWT_SERVER, timing); }
  
  //printf("APPLY BWT SERVER DONE!\n");
  return CONSUMER_STAGE;

  if (batch->mapping_batch->num_targets > 0) {
    //TODO: DELETE
    printf("We have targets\n");
    for (int i = 0; i < batch->mapping_batch->num_targets; i++) {
      batch->mapping_batch->bwt_mappings[batch->mapping_batch->targets[i]] = 1;
    }
    return SEEDING_STAGE;
  }

  //printf("Reads are mapped\n");
  //return POST_PAIR_STAGE;
  return CONSUMER_STAGE;
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
