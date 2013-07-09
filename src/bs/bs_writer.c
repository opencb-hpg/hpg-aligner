#include "bs_writer.h"
#include "methylation.h"

int bs_writer(void *data) {

  printf("\t-----> bs_writer\n");

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

  printf("\t-----> bs_writer 0000\n");
     // set the sequences of the mapping to the original
  {
    size_t num_reads = array_list_size(mapping_batch->fq_batch);
    printf("\t-----> num_reads %i\n", num_reads);
  }
  revert_mappings_seqs(mapping_batch->mapping_lists, mapping_batch->mapping_lists2, mapping_batch->fq_batch);

  printf("\t-----> bs_writer 1\n");

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

     array_list_t **mapping_lists;
     int *found = (int *) calloc(num_reads, sizeof(int));
     
  printf("\t-----> bs_writer 2\n");

     // process mapping_lists and mapping_lists2
     for (int k = 0; k < 2; k++) {
       
       mapping_lists = (k == 0) ? mapping_batch->mapping_lists : mapping_batch->mapping_lists2;
       
       for (size_t i = 0; i < num_reads; i++) {
	 num_items = array_list_size(mapping_lists[i]);
	 total_mappings += num_items;
	 fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->fq_batch);
	 
	 printf("\t-----> num_items (%i) = %i\n", k, num_items);

	 // mapped or not mapped ?	 
	 if (num_items == 0) {
	   //total_mappings++;
	   //	   write_unmapped_read(fq_read, bam_file); //
	   if (mapping_lists[i]) {
	     array_list_free(mapping_lists[i], NULL);
	   }
	 } else {
	   found[i] = 1;
	   //	   num_mapped_reads++; //
	   write_mapped_read(mapping_lists[i], bam_file);
	 }
       }
     }

     for (size_t i = 0; i < num_reads; i++) {
       if (found[i]) {
	 num_mapped_reads++; //
       } else {
	 total_mappings++;
	 fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->fq_batch);
	 write_unmapped_read(fq_read, bam_file);
       }
     }

     if (basic_st->total_reads >= writer_input->limit_print) {
       LOG_DEBUG_F("TOTAL READS PROCESS: %lu\n", basic_st->total_reads);
       LOG_DEBUG_F("\tTotal Reads Mapped: %lu(%.2f%)\n", 
		   basic_st->num_mapped_reads, 
		   (float) (basic_st->num_mapped_reads*100)/(float)(basic_st->total_reads));
       writer_input->limit_print += 1000000;
     }
     
     //printf("Batch Write OK!\n");     
     
     if (mapping_batch) {
       mapping_batch_free(mapping_batch);
     }
     
     if (batch) batch_free(batch);

     if (found) free(found);
     
     basic_statistics_add(num_reads, num_mapped_reads, total_mappings, basic_st);
     
     if (time_on) { stop_timer(start, end, time); timing_add(time, BAM_WRITER, timing); }
}
