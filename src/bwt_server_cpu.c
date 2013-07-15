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

size_t bwt_search_pair_anchors(array_list_t *anchor_list) {
  bwt_anchor_t *bwt_anchor;
  int max_double_anchor = 0, max_anchor = 0;
  bwt_anchor_t *max_backward = NULL, *max_forward = NULL;                                                                                                         
  bwt_anchor_t *bwt_anchor_back, *bwt_anchor_forw;
  int anchor_tmp, anchor_back, anchor_forw;
  int strand, type;
  int found_anchor = 0, found_double_anchor = 0;
  const int MIN_ANCHOR = 16;
  const int MIN_DOUBLE_ANCHOR = MIN_ANCHOR*2;
  const int MAX_BWT_REGIONS = 50;
  const int MAX_BWT_ANCHOR_DISTANCE = 1000000;

  array_list_t *anchor_list_tmp, *forward_anchor_list, *backward_anchor_list;
  array_list_t *backward_anchor_list_0 = array_list_new(MAX_BWT_REGIONS + 1, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *forward_anchor_list_0 = array_list_new(MAX_BWT_REGIONS + 1, 1.25f , COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *backward_anchor_list_1 = array_list_new(MAX_BWT_REGIONS + 1, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *forward_anchor_list_1 = array_list_new(MAX_BWT_REGIONS + 1, 1.25f , COLLECTION_MODE_ASYNCHRONIZED);

  //printf("Tot Anchors %i\n", array_list_size(anchor_list));

  for (int i = 0; i < array_list_size(anchor_list); i++) {
    bwt_anchor = array_list_get(i, anchor_list);
    if (bwt_anchor->strand == 1) {
      //printf("(-)bwt anchor %i:%lu-%lu (%i): ", bwt_anchor->chromosome + 1, bwt_anchor->start, bwt_anchor->end, bwt_anchor->end - bwt_anchor->start);
      if (bwt_anchor->type == FORWARD_ANCHOR) {
	array_list_insert(bwt_anchor, forward_anchor_list_1);
	//printf("FORW\n");
      } else {
	array_list_insert(bwt_anchor, backward_anchor_list_1);
	//printf("BACK\n");
      }
    } else {
      //printf("(+)bwt anchor %i:%lu-%lu (%i): ", bwt_anchor->chromosome + 1, bwt_anchor->start, bwt_anchor->end, bwt_anchor->end - bwt_anchor->start);
      if (bwt_anchor->type == FORWARD_ANCHOR) {
	array_list_insert(bwt_anchor, forward_anchor_list_0);
	//printf("FORW\n");
      } else {
	array_list_insert(bwt_anchor, backward_anchor_list_0);
	//printf("BACK\n");
      }      
    }
    anchor_tmp = bwt_anchor->end - bwt_anchor->start;
    if (anchor_tmp > MIN_ANCHOR && anchor_tmp > max_anchor) {
      max_anchor = anchor_tmp;
      found_anchor = 1;
      strand = bwt_anchor->strand;
      type = bwt_anchor->type;
    }
  }
  
  array_list_clear(anchor_list, NULL);

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

    //max_double_anchor = 0;//MIN_DOUBLE_ANCHOR;
    //Associate Anchors (+)/(-)
    for (int i = 0; i < array_list_size(forward_anchor_list); i++) { 
      bwt_anchor_forw = array_list_get(i, forward_anchor_list);
      for (int j = 0; j < array_list_size(backward_anchor_list); j++) { 
	bwt_anchor_back = array_list_get(j, backward_anchor_list);
	anchor_forw = (bwt_anchor_forw->end - bwt_anchor_forw->start);
	anchor_back = (bwt_anchor_back->end - bwt_anchor_back->start); 

	anchor_tmp = anchor_forw + anchor_back;

	//printf("\tCommpare %i:%lu-%lu with %i:%lu-%lu\n", bwt_anchor_forw->chromosome + 1, 
	//     bwt_anchor_forw->start, bwt_anchor_forw->end, bwt_anchor_back->chromosome + 1, 
	//     bwt_anchor_back->start, bwt_anchor_back->end);
	if (bwt_anchor_forw->chromosome == bwt_anchor_back->chromosome &&
	    abs(bwt_anchor_back->start - bwt_anchor_forw->end) <= MAX_BWT_ANCHOR_DISTANCE && 
	    anchor_forw >= MIN_ANCHOR && anchor_back >= MIN_ANCHOR) {

	  if (bwt_anchor_back->start < bwt_anchor_forw->end) { continue; }

	  array_list_insert(bwt_anchor_new(bwt_anchor_forw->strand, bwt_anchor_forw->chromosome, 
					   bwt_anchor_forw->start,  bwt_anchor_forw->end, bwt_anchor_forw->type), anchor_list);
	  
	  array_list_insert(bwt_anchor_new(bwt_anchor_back->strand, bwt_anchor_back->chromosome, 
					   bwt_anchor_back->start,  bwt_anchor_back->end, bwt_anchor_back->type), anchor_list);
	  
	  array_list_set_flag(2, anchor_list);
	  //max_backward = bwt_anchor_back;
	  //max_forward = bwt_anchor_forw;
	  //max_double_anchor = anchor_tmp;
	  found_double_anchor = 1;
	  //printf("\t Found!\n");
	}                                                                                                                      
      }         
    }
  }

  if (!found_double_anchor && 
      found_anchor) { 
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
    for (int i = 0; i < array_list_size(anchor_list_tmp); i++) {
      bwt_anchor = array_list_get(i, anchor_list_tmp);
      array_list_insert(bwt_anchor_new(bwt_anchor->strand, bwt_anchor->chromosome, 
				       bwt_anchor->start, bwt_anchor->end, bwt_anchor->type), anchor_list);
      array_list_set_flag(1, anchor_list);
    }
  }

  array_list_free(forward_anchor_list_1, bwt_anchor_free);
  array_list_free(backward_anchor_list_1, bwt_anchor_free);
  array_list_free(forward_anchor_list_0, bwt_anchor_free);
  array_list_free(backward_anchor_list_0, bwt_anchor_free);
  
  return array_list_size(anchor_list);
  
}


int apply_bwt(bwt_server_input_t* input, batch_t *batch) {
  //printf("APPLY BWT SERVER...\n");
  struct timeval start, end;
  double time;

  mapping_batch_t *mapping_batch = batch->mapping_batch;
  size_t num_reads = array_list_size(mapping_batch->fq_batch);
  size_t num_mappings;
  array_list_t *list;
  size_t *unmapped_indices = mapping_batch->targets;
  size_t num_unmapped = 0;
  size_t num_anchors;
  for (int i = 0; i < num_reads; i++) {
    fastq_read_t *read = array_list_get(i, mapping_batch->fq_batch);
    //printf("%s\n", read->id);
    list = mapping_batch->mapping_lists[i];    
    array_list_set_flag(1, list);
    num_mappings = bwt_map_inexact_read(read,
					input->bwt_optarg_p,
					input->bwt_index_p,
					list);
    if (array_list_get_flag(list) != 2) { //If flag 2, the read exceded the max number of mappings
      if (array_list_get_flag(list) == 1) {
	if (num_mappings > 0) {
	  num_anchors = bwt_search_pair_anchors(list);	
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
      } //else {
	//printf("#Read Mapped %s\n", read->id);
      //}
    }
  }
  // array_list flag: 0 -> Not  BWT Anchors found
  //                  1 -> One  BWT Anchors found
  //                  2 -> Pair BWT Anchors found
  
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
