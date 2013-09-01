#include "workflow_functions.h"

int limit_of_reads = 0;
int n_insert = 0;

//====================================================================================

wf_input_t *wf_input_new(fastq_batch_reader_input_t *fq_reader_input,
                         batch_t *batch) {

  wf_input_t *wfi = (wf_input_t *) calloc(1, sizeof(wf_input_t));
  wfi->fq_reader_input = fq_reader_input;
  wfi->batch = batch;

  return wfi;
}

wf_input_buffer_t *wf_input_buffer_new(linked_list_t *buffer,
				       batch_t *batch) {
  wf_input_buffer_t *wfi = (wf_input_buffer_t *) calloc(1, sizeof(wf_input_buffer_t));
  wfi->buffer = buffer;
  wfi->batch = batch;

  return wfi;
}

void wf_input_free(wf_input_t *wfi) {
  if (wfi) free(wfi);
}

//--------------------------------------------------------------------
// workflow producer
//--------------------------------------------------------------------

void *fastq_reader(void *input) {
     struct timeval start, end;
     double time;

     if (time_on) { start_timer(start); }

     wf_input_t *wf_input = (wf_input_t *) input;
     batch_t *new_batch = NULL;
     batch_t *batch = wf_input->batch;
     fastq_batch_reader_input_t *fq_reader_input = wf_input->fq_reader_input;
     array_list_t *reads = array_list_new(10000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

     if (fq_reader_input->flags == SINGLE_END_MODE) {
       fastq_fread_bytes_se(reads, fq_reader_input->batch_size, fq_reader_input->fq_file1);
     } else {
       fastq_fread_bytes_aligner_pe(reads, fq_reader_input->batch_size, 
				    fq_reader_input->fq_file1, fq_reader_input->fq_file2);
     }

     size_t num_reads = array_list_size(reads);

     if (num_reads == 0) {
	  array_list_free(reads, (void *)fastq_read_free);
     } else {
	  mapping_batch_t *mapping_batch = mapping_batch_new(reads, 
							     batch->pair_input->pair_mng);

	  new_batch = batch_new(batch->bwt_input, batch->region_input, batch->cal_input, 
				batch->pair_input, batch->preprocess_rna, batch->sw_input, batch->writer_input, 
				batch->mapping_mode, mapping_batch);
     }

     if (time_on) { stop_timer(start, end, time); timing_add(time, FASTQ_READER, timing); }

     return new_batch;
}

void *buffer_reader(void *input) {
  wf_input_buffer_t *wf_input = (wf_input_t *) input;

  linked_list_t *buffer = wf_input->buffer;
  batch_t *batch = wf_input->batch;
  buffer_item_t *buffer_item;
  const int MAX_READS = 100;
  int num_reads = 0;
  batch_t *new_batch = NULL;

  if (linked_list_size(buffer) > 0) {
    array_list_t *reads = array_list_new(MAX_READS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    mapping_batch_t *mapping_batch = mapping_batch_new_2(MAX_READS, 
							 reads,
							 batch->pair_input->pair_mng);
    while (num_reads < MAX_READS) {
      buffer_item = linked_list_remove_last(buffer);
      if (buffer_item == NULL) { break; }
      fastq_read_t *read = buffer_item->read;
      array_list_insert(buffer_item->read, reads);
      mapping_batch->mapping_lists[num_reads] = array_list_new(50,
							       1.25f, 
							       COLLECTION_MODE_ASYNCHRONIZED);

      for (int i = 0; i < array_list_size(buffer_item->items_list); i++) {
	void *item = array_list_get(i, buffer_item->items_list);
	array_list_insert(item, mapping_batch->mapping_lists[num_reads]);
      }

      array_list_set_flag(array_list_get_flag(buffer_item->items_list),
			  mapping_batch->mapping_lists[num_reads]);
      num_reads++;
      //printf("TOTAL READS %i\n", num_reads);
      buffer_item_free(buffer_item);
    }
    
    mapping_batch->num_allocated_targets = num_reads;
    new_batch = batch_new(batch->bwt_input, batch->region_input, batch->cal_input, 
			  batch->pair_input, batch->preprocess_rna, batch->sw_input,
			  batch->writer_input, batch->mapping_mode, mapping_batch); 
  }


  return new_batch;
}

//--------------------------------------------------------------------
// workflow consumer
//--------------------------------------------------------------------

int search_hard_clipping(array_list_t *array_list);
void write_mapped_read(array_list_t *array_list, bam_file_t *bam_file);
void write_unmapped_read(fastq_read_t *fq_read, bam_file_t *bam_file);

//--------------------------------------------------------------------

int num_sws = 0;
int num_ext_sws = 0;
int num_gaps = 0;

int bam_writer(void *data) {
     struct timeval start, end;
     double time;

     if (time_on) { start_timer(start); }

     batch_t *batch = (batch_t *) data;
     fastq_read_t *fq_read;
     array_list_t *array_list;
     size_t num_items;

     //bam1_t *bam1;
     //alignment_t *alig;

     mapping_batch_t *mapping_batch = (mapping_batch_t *) batch->mapping_batch;

     num_sws += mapping_batch->num_sws;
     num_ext_sws += mapping_batch->num_ext_sws;
     num_gaps += mapping_batch->num_gaps;

     batch_writer_input_t *writer_input = batch->writer_input;
     bam_file_t *bam_file = writer_input->bam_file;     
     linked_list_t *linked_list = writer_input->list_p;
     size_t num_reads = array_list_size(mapping_batch->fq_batch);
     size_t num_mapped_reads = 0;
     size_t total_mappings = 0;
     unsigned char found_p1 = 0;
     unsigned char found_p2 = 0;
     int i = 0;

     extern size_t bwt_correct;
     extern size_t bwt_error;
     extern pthread_mutex_t bwt_mutex;

     writer_input->total_batches++;

     extern size_t *histogram_sw;

     if (batch->mapping_mode == DNA_MODE || batch->mapping_mode == RNA_MODE) {
       //for (int i = 0; i < 1024; i++) {
       //histogram_sw[i] += mapping_batch->histogram_sw[i];
       //}
       free(mapping_batch->histogram_sw);
       //
       // DNA mode
       //
       for (size_t i = 0; i < num_reads; i++) {
	 num_items = array_list_size(mapping_batch->mapping_lists[i]);
	 total_mappings += num_items;
	 fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->fq_batch);

	 // mapped or not mapped ?	 
	 if (num_items == 0) {
	   total_mappings++;
	   write_unmapped_read(fq_read, bam_file);
	   if (mapping_batch->mapping_lists[i]) {
	     array_list_free(mapping_batch->mapping_lists[i], NULL);
	   }	 
	 } else {
	   num_mapped_reads++;
	   write_mapped_read(mapping_batch->mapping_lists[i], bam_file);
	 }
       }
     } else {
       //TODO: This section needs supracals implementatation, that it is not implemented yet.
       //
       // RNA mode
       //
       i = 0;
       while (i < num_reads) {
	 num_items = array_list_size(mapping_batch->mapping_lists[i]);
	 total_mappings += num_items;
	 fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->fq_batch);
	 // mapped or not mapped ?
	 if (!num_items) {
	   total_mappings++;
	   write_unmapped_read(fq_read, bam_file);
	 } else {
	   //Search if some read contains a left/right hard clipping
	   if (batch->pair_input->pair_mng->pair_mode == SINGLE_END_MODE) {
	     if (search_hard_clipping(mapping_batch->mapping_lists[i])) {
	       //Store to Buffer fq_read and alignments
	       buffer_item_t *buffer_item = buffer_item_new(fastq_read_dup(buffer_item_new(fq_read)), 
							    mapping_batch->mapping_lists[i]);

	       linked_list_insert((void *)buffer_item, linked_list);
	     } else {
	       //Not Found Hard Clipping! Write to Diks all alignments
	       num_mapped_reads++;
	       write_mapped_read(mapping_batch->mapping_lists[i], bam_file);
	     }
	   } else {
	      if (i%2 == 0) {
		found_p1 = search_hard_clipping(mapping_batch->mapping_lists[i]);
		if (!mapping_batch->mapping_lists[i + 1]) { printf("Pair-End without pair.\n"); break;  }
		found_p2 = search_hard_clipping(mapping_batch->mapping_lists[i + 1]);
		if (found_p1 || found_p2) {
		  //Store to Buffer fq_read and fq_read + 1 and alignments
		  buffer_pair_item_t *buffer_item = buffer_item_pair_new(fastq_read_dup(buffer_item_new(fq_read)), 
									 mapping_batch->mapping_lists[i], 
									 fastq_read_dup(array_list_get(i + 1, 
												       mapping_batch->fq_batch)),
									 mapping_batch->mapping_lists[i + 1]);
		  linked_list_insert((void *)buffer_item, linked_list);
		} else {
		  //Not Found Hard Clipping in any Pair! Write to Diks all alignments
		  if (num_items) { num_mapped_reads++; }
		  num_items = array_list_size(mapping_batch->mapping_lists[i + 1]);
		  if (num_items) { num_mapped_reads++;total_mappings += num_items; }
		  
		  write_mapped_read(mapping_batch->mapping_lists[i], bam_file);
		  write_mapped_read(mapping_batch->mapping_lists[i + 1], bam_file);
		}
		i++;
	      }
	   }
	 } //else num_items
	 i++;
       } //end of while

       if (global_status == WORKFLOW_STATUS_FINISHED) {
	 pthread_cond_signal(&cond_sp);
       }
     }
     
     //fprintf(stderr, "TOTAL READS PROCESS: %lu\n", basic_st->total_reads);
     if (basic_st->total_reads >= writer_input->limit_print) {
       //LOG_DEBUG_F("TOTAL READS PROCESS: %lu\n", basic_st->total_reads);
       //LOG_DEBUG_F("\tTotal Reads Mapped: %lu(%.2f%)\n", 
       //	   basic_st->num_mapped_reads, 
       //	   (float) (basic_st->num_mapped_reads*100)/(float)(basic_st->total_reads));
       //writer_input->limit_print += 1000000;
       fprintf(stderr, "TOTAL READS PROCESS: %lu\n", basic_st->total_reads);
       printf("\tTotal Reads Mapped: %lu(%.2f%)\n", 
		   basic_st->num_mapped_reads, 
		   (float) (basic_st->num_mapped_reads*100)/(float)(basic_st->total_reads));
       
       writer_input->limit_print += 100000;
     }
     
     //printf("Batch Write OK!\n");     
     
     if (mapping_batch) {
       mapping_batch_free(mapping_batch);
     }
     
     if (batch) batch_free(batch);
     
     basic_statistics_add(num_reads, num_mapped_reads, total_mappings, basic_st);
     
     if (time_on) { stop_timer(start, end, time); timing_add(time, BAM_WRITER, timing); }
}

//--------------------------------------------------------------------

int search_hard_clipping(array_list_t *array_list) {
  size_t num_items = array_list_size(array_list);
  alignment_t *alig;
  int found = 0;

  for (size_t j = 0; j < num_items; j++) {
    alig = (alignment_t *) array_list_get(j, array_list);
    if (alig->large_hard_clipping) {
      return 0;
    }
  }
  return 0;
}

//--------------------------------------------------------------------

void write_mapped_read(array_list_t *array_list, bam_file_t *bam_file) {
  size_t num_items = array_list_size(array_list);
  alignment_t *alig;
  bam1_t *bam1;
  for (size_t j = 0; j < num_items; j++) {
    alig = (alignment_t *) array_list_get(j, array_list);

    if (alig != NULL) {
      bam1 = convert_to_bam(alig, 33);
      bam_fwrite(bam1, bam_file);
      bam_destroy1(bam1);	 
      alignment_free(alig);
    } else {
      LOG_FATAL_F("alig is NULL, num_items = %lu\n", num_items)
    }
  }
  if (array_list) { array_list_free(array_list, NULL); }
}

//--------------------------------------------------------------------

void write_unmapped_read(fastq_read_t *fq_read, bam_file_t *bam_file) {
  static char aux[1024];
  alignment_t *alig;
  size_t header_len;
  char *id;
  bam1_t *bam1;

  // calculating cigar
  sprintf(aux, "%luX", fq_read->length);	       
  alig = alignment_new();	       
  header_len = strlen(fq_read->id);
  id = (char *) malloc(sizeof(char) * (header_len + 1));
  get_to_first_blank(fq_read->id, header_len, id);
  //free(fq_read->id);
  alignment_init_single_end(id, fq_read->sequence, fq_read->quality, 
			    0, -1, -1, aux, 1, 0, 0, 0, 0, NULL, 0, alig);

  bam1 = convert_to_bam(alig, 33);
  bam_fwrite(bam1, bam_file);
  bam_destroy1(bam1);
	       
  alig->sequence = NULL;
  alig->quality = NULL;
  alig->cigar = NULL;
  alignment_free(alig);	       
  //printf("\tWRITE : read %i (%d items): unmapped...done !!\n", i, num_items);
  
}

//--------------------------------------------------------------------
// stage functions
//--------------------------------------------------------------------

int bwt_stage(void *data) {
     batch_t *batch = (batch_t *) data;

     if (batch->mapping_mode == DNA_MODE) {
       return apply_bwt(batch->bwt_input, batch);     
     } else {
       return apply_bwt_rna(batch->bwt_input, batch);     
     }
}

//--------------------------------------------------------------------

int seeding_stage(void *data) {
     batch_t *batch = (batch_t *) data;

     return apply_seeding(batch->region_input, batch);
}

//--------------------------------------------------------------------

int cal_stage(void *data) {
     batch_t *batch = (batch_t *) data;

     if (batch->mapping_mode == DNA_MODE) {
       return apply_caling(batch->cal_input, batch);
     } else {
       return apply_caling_rna(batch->cal_input, batch);
     }
}

//--------------------------------------------------------------------

int rna_preprocess_stage(void *data) {
  batch_t *batch = (batch_t *) data;

  return apply_preprocess_rna(batch->preprocess_rna, batch);  
}

//---------------------------------------------------------------------

int pre_pair_stage(void *data) {
     batch_t *batch = (batch_t *) data;
     return apply_pair(batch->pair_input, batch);
}

//--------------------------------------------------------------------

int sw_stage(void *data) {
     batch_t *batch = (batch_t *) data;
     
     if (batch->mapping_mode == RNA_MODE) {
       return apply_sw_rna(batch->sw_input, batch);
     } else {
       return apply_sw(batch->sw_input, batch);
     }
}

//--------------------------------------------------------------------

int rna_last_stage(void *data) {
   batch_t *batch = (batch_t *) data;
   return apply_rna_last(batch->sw_input, batch);
}

//--------------------------------------------------------------------

int rna_last_hc_stage(void *data) {
   batch_t *batch = (batch_t *) data;
   return apply_rna_last_hc(batch->sw_input, batch);
}

//--------------------------------------------------------------------

int post_pair_stage(void *data) {
     batch_t *batch = (batch_t *) data;
     return prepare_alignments(batch->pair_input, batch);
}

//--------------------------------------------------------------------
