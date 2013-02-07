#include "dna_aligner.h"

//--------------------------------------------------------------------
bam_header_t *create_bam_header_by_genome(genome_t *genome);
//--------------------------------------------------------------------
  
#define BWT_STAGE         0
#define SEEDING_STAGE     1
#define CAL_STAGE         2
#define PRE_PAIR_STAGE    3
#define SW_STAGE          4
#define POST_PAIR_STAGE   5
#define CONSUMER_STAGE   -1

//--------------------------------------------------------------------
// structure between the different workflow stages
//--------------------------------------------------------------------

typedef struct batch {
     bwt_server_input_t *bwt_input;
     region_seeker_input_t *region_input;
     cal_seeker_input_t *cal_input;
     pair_server_input_t *pair_input;
     sw_server_input_t *sw_input;
     batch_writer_input_t *writer_input;
     mapping_batch_t *mapping_batch;
} batch_t;

batch_t *batch_new(bwt_server_input_t *bwt_input,
		   region_seeker_input_t *region_input,
		   cal_seeker_input_t *cal_input,
		   pair_server_input_t *pair_input,
		   sw_server_input_t *sw_input,
		   batch_writer_input_t *writer_input,
		   mapping_batch_t *mapping_batch) {

     batch_t *b = (batch_t *) calloc(1, sizeof(batch_t));
     b->bwt_input = bwt_input;
     b->region_input = region_input;
     b->cal_input = cal_input;
     b->pair_input = pair_input;
     b->sw_input = sw_input;
     b->writer_input = writer_input;
     b->mapping_batch = mapping_batch;
     return b;
}

void batch_free(batch_t *b) {
     if (b) free(b);
}

//--------------------------------------------------------------------
// workflow input
//--------------------------------------------------------------------

typedef struct wf_input {
     fastq_batch_reader_input_t *fq_reader_input;
     batch_t *batch;
} wf_input_t;

wf_input_t *wf_input_new(fastq_batch_reader_input_t *fq_reader_input,
			 batch_t *batch) {
     
     wf_input_t *wfi = (wf_input_t *) calloc(1, sizeof(wf_input_t));
     wfi->fq_reader_input = fq_reader_input;
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
				batch->pair_input, batch->sw_input, batch->writer_input, 
				mapping_batch);
     }
     

     return new_batch;
}

//--------------------------------------------------------------------
// workflow consumer
//--------------------------------------------------------------------

int bam_writer(void *data) {
     batch_t *batch = (batch_t *) data;

     static char aux[1024];

     char *id;
     fastq_read_t *fq_read;
     array_list_t *array_list;
     size_t num_items, total_mappings, header_len;

     bam1_t *bam1;
     alignment_t *alig;

     mapping_batch_t *mapping_batch = (mapping_batch_t *) batch->mapping_batch;

     batch_writer_input_t *writer_input = batch->writer_input;
     bam_file_t *bam_file = writer_input->bam_file;
     
     
     size_t num_reads = array_list_size(mapping_batch->fq_batch);

     writer_input->total_batches++;
     writer_input->total_reads += num_reads;

     //printf("Read %i reads\n", num_reads);
     for (size_t i = 0; i < num_reads; i++) {
	  array_list = mapping_batch->mapping_lists[i];
	  num_items = array_list_size(array_list);
	  writer_input->total_mappings += num_items;
	  // mapped or not mapped ?
	  if (num_items == 0) {
	       writer_input->total_mappings++;

	       fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->fq_batch);

	       // calculating cigar
	       sprintf(aux, "%luX", fq_read->length);
	       
	       alig = alignment_new();
	       
	       header_len = strlen(fq_read->id);
	       id = (char *) malloc(sizeof(char) * (header_len + 1));
	       get_to_first_blank(fq_read->id, header_len, id);
	       //free(fq_read->id);
	       alignment_init_single_end(id, fq_read->sequence, fq_read->quality, 
					 0, -1, -1, aux, 1, 0, 0, 0, 0, NULL, alig);
	       
	       bam1 = convert_to_bam(alig, 33);
	       bam_fwrite(bam1, bam_file);
	       bam_destroy1(bam1);
	       
	       // some cosmetic stuff before freeing the alignment,
	       // (in order to not free twice some fields)
	       //alig->query_name = NULL;
	       alig->sequence = NULL;
	       alig->quality = NULL;
	       alig->cigar = NULL;
	       alignment_free(alig);
	       
	       //	printf("\tWRITE : read %i (%d items): unmapped...done !!\n", i, num_items);
	       
	  } else {
	       writer_input->num_mapped_reads++;
	       //printf("\tWRITE : read %d (%d items): mapped...\n", i, num_items);
	       if (num_items > 1000) {
		    alig = (alignment_t *) array_list_get(0, array_list); 
		    printf("%s\n", alig->sequence); 
	       }
	       for (size_t j = 0; j < num_items; j++) {
		    alig = (alignment_t *) array_list_get(j, array_list);
		    //printf("\t%s\n", alig->cigar);
		    if (alig != NULL) {			 
			 bam1 = convert_to_bam(alig, 33);
			 bam_fwrite(bam1, bam_file);
			 bam_destroy1(bam1);			 
			 alignment_free(alig);
		    }
	       }
	       //	printf("\tWRITE : read %d (%d items): mapped...done !!\n", i, num_items);
	  }
	  if (array_list) array_list_free(array_list, NULL);
     }
     
     if (writer_input->total_reads >= writer_input->limit_print) {
	  LOG_DEBUG_F("TOTAL READS PROCESS: %lu\n", writer_input->total_reads);
	  LOG_DEBUG_F("\tTotal Reads Mapped: %lu(%.2f%)\n", 
		      writer_input->num_mapped_reads, 
		      (float) (writer_input->num_mapped_reads*100)/(float)(writer_input->total_reads));
	  writer_input->limit_print += 1000000;
     }
     
     //printf("Batch Write OK!\n");
     
     if (mapping_batch) {
	  mapping_batch_free(mapping_batch);
     }
     if (batch) batch_free(batch);

}

//--------------------------------------------------------------------
// stage functions
//--------------------------------------------------------------------

int bwt_stage(void *data) {
     batch_t *batch = (batch_t *) data;
//     apply_bwt(batch->bwt_input, batch->mapping_batch);

     // go to the next stage
     if (batch->mapping_batch->num_targets > 0) {
	  return SEEDING_STAGE;
     } else if (batch->pair_input->pair_mng->pair_mode != SINGLE_END_MODE) {      
	  return PRE_PAIR_STAGE;
     }
     return CONSUMER_STAGE;
}

int seeding_stage(void *data) {
     batch_t *batch = (batch_t *) data;
     apply_seeding(batch->region_input, batch->mapping_batch);

     // go to the next stage
     return CAL_STAGE;
}

int cal_stage(void *data) {
     batch_t *batch = (batch_t *) data;
     apply_caling(batch->cal_input, batch->mapping_batch);

     // go to the next stage
     if (batch->pair_input->pair_mng->pair_mode != SINGLE_END_MODE) {
	  return PRE_PAIR_STAGE;
     } else if (batch->mapping_batch->num_targets > 0) {
	  return SW_STAGE;
     }
     return CONSUMER_STAGE;
}

int pre_pair_stage(void *data) {
     batch_t *batch = (batch_t *) data;
     apply_pair(batch->pair_input, batch->mapping_batch);

     // go to the next stage
     if (batch->mapping_batch->num_targets > 0) {
	  return SW_STAGE;
     }
     return CONSUMER_STAGE;
}

int sw_stage(void *data) {
     batch_t *batch = (batch_t *) data;
     apply_sw(batch->sw_input, batch->mapping_batch);

     // go to the next stage
     if (batch->mapping_batch->num_targets > 0) {
	  return POST_PAIR_STAGE;
     }
     return CONSUMER_STAGE;
}

int post_pair_stage(void *data) {
     batch_t *batch = (batch_t *) data;
     prepare_alignments(batch->pair_input, batch->mapping_batch);
     
     return CONSUMER_STAGE;
}

//--------------------------------------------------------------------

void run_dna_aligner(genome_t *genome, bwt_index_t *bwt_index, 
		     bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg, 
		     pair_mng_t *pair_mng, options_t *options) {

     int path_length = strlen(options->output_name);
     int extend_length = 0;
     if (options->extend_name) {
	  extend_length = strlen(options->extend_name);
     }
     
     char *reads_results = (char *) calloc((60 + extend_length), sizeof(char));
     char *output_filename = (char *) calloc((path_length + extend_length + 60), sizeof(char));
     
     if (options->extend_name) {
	  strcat(reads_results, "/");
	  strcat(reads_results, options->extend_name);
	  strcat(reads_results, "_alignments.bam");  
     } else {
	  strcat(reads_results, "/alignments.bam");
     }
     
     strcat(output_filename, options->output_name);
     strcat(output_filename, reads_results);
     free(reads_results);
     
     // display selected options
     LOG_DEBUG("Displaying options...\n");
     options_display(options);
     
     omp_set_nested(1);  
     
     // preparing input FastQ file
     fastq_batch_reader_input_t reader_input;
     fastq_batch_reader_input_init(options->in_filename, options->in_filename2, 
				   options->pair_mode, options->batch_size, 
				   NULL, &reader_input);
     if (options->pair_mode == SINGLE_END_MODE) {
	  reader_input.fq_file1 = fastq_fopen(options->in_filename);
     } else {
	  reader_input.fq_file1 = fastq_fopen(options->in_filename);
	  reader_input.fq_file2 = fastq_fopen(options->in_filename2);
     }
     
     // preparing output BAM file
     batch_writer_input_t writer_input;
     batch_writer_input_init(output_filename, NULL, NULL, NULL, genome, &writer_input);
     bam_header_t *bam_header = create_bam_header_by_genome(genome);
     writer_input.bam_file = bam_fopen_mode(output_filename, bam_header, "w");
     bam_fwrite_header(bam_header, writer_input.bam_file);

     // preparing workflow stage input
     bwt_server_input_t bwt_input;
     bwt_server_input_init(NULL, 0, bwt_optarg, bwt_index, 
			   NULL, 0, NULL, &bwt_input);
     
     region_seeker_input_t region_input;
     region_seeker_input_init(NULL, cal_optarg, bwt_optarg, 
			      bwt_index, NULL, 0, options->gpu_process, 
			      &region_input);
     
     cal_seeker_input_t cal_input;
     cal_seeker_input_init(NULL, cal_optarg, NULL, 0, 
			   NULL, NULL, &cal_input);
     
     pair_server_input_t pair_input;
     pair_server_input_init(pair_mng, bwt_optarg->report_best, bwt_optarg->report_n_hits, 
			    bwt_optarg->report_all, NULL, NULL, NULL, &pair_input);
     
     sw_server_input_t sw_input;
     sw_server_input_init(NULL, NULL, 0, options->match, options->mismatch, 
			  options->gap_open, options->gap_extend, options->min_score, 
			  options->flank_length, genome, 0, 0, 0,  bwt_optarg, &sw_input);

     if (time_on) { 
	  timing_start(MAIN_INDEX, 0, timing_p);
     }
     
     //--------------------------------------------------------------------------------------
     // workflow management
     //
     //
     batch_t *batch = batch_new(&bwt_input, &region_input, &cal_input, 
				&pair_input, &sw_input, &writer_input, NULL);
     wf_input_t *wf_input = wf_input_new(&reader_input, batch);

     // create and initialize workflow
     workflow_t *wf = workflow_new();
     
     workflow_stage_function_t stage_functions[] = {bwt_stage, seeding_stage, cal_stage, 
						    pre_pair_stage, sw_stage, post_pair_stage};
     char *stage_labels[] = {"BWT", "SEEDING", "CAL", "PRE PAIR", "SW", "POST PAIR"};
     workflow_set_stages(6, &stage_functions, stage_labels, wf);
     
     // optional producer and consumer functions
     workflow_set_producer(fastq_reader, "FastQ reader", wf);
     workflow_set_consumer(bam_writer, "BAM writer", wf);
     
     workflow_run_with(options->num_cpu_threads, wf_input, wf);
     
     // free memory
     workflow_free(wf);
     wf_input_free(wf_input);
     batch_free(batch);
     //
     //
     // end of workflow management
     //--------------------------------------------------------------------------------------


     //closing files
     if (options->pair_mode == SINGLE_END_MODE) {
	  fastq_fclose(reader_input.fq_file1);
     } else {
	  fastq_fclose(reader_input.fq_file1);
	  fastq_fclose(reader_input.fq_file2);
     }
     bam_fclose(writer_input.bam_file);


     if (time_on) { 
	  timing_stop(MAIN_INDEX, 0, timing_p);
     }
     
     if (statistics_on) {
	  size_t total_item = 0;
	  double max_time = 0, total_throughput = 0;
	  printf("\nBWT time:\n");
	  for (int i = 0; i < options->num_cpu_threads; i++) {
	       printf("\tThread %d: %0.4f s (%d batches, %d items) -> %0.2f BWT/s (reads)\n", 
		      i, bwt_time[i] / 1e6, thr_batches[i], thr_bwt_items[i], 1e6 * thr_bwt_items[i] / bwt_time[i]);
	       total_item += thr_bwt_items[i];
	       total_throughput += (1e6 * thr_bwt_items[i] / bwt_time[i]);
	       if (max_time < bwt_time[i]) max_time = bwt_time[i];
	  }
	  printf("\n\tTotal BWTs: %lu, Max time = %0.4f, Throughput = %0.2f BWT/s\n", 
		 total_item, max_time / 1e6, total_throughput);
	  
	  total_item = 0; max_time = 0; total_throughput = 0;
	  printf("\nSeeding time:\n");
	  for (int i = 0; i < options->num_cpu_threads; i++) {
	       printf("\tThread %d: %0.4f s (%d batches, %d items) -> %0.2f BWT/s (seeds)\n", 
		      i, seeding_time[i] / 1e6, thr_batches[i], thr_seeding_items[i], 
		      1e6 * thr_seeding_items[i] / seeding_time[i]);
	       total_item += thr_seeding_items[i];
	       total_throughput += (1e6 * thr_seeding_items[i] / seeding_time[i]);
	       if (max_time < seeding_time[i]) max_time = seeding_time[i];
	  }
	  printf("\n\tTotal BWTs: %lu, Max time = %0.4f, Throughput = %0.2f BWT/s\n", 
		 total_item, max_time / 1e6, total_throughput);
	  
	  total_item = 0; max_time = 0; total_throughput = 0;
	  printf("\nCAL time:\n");
	  for (int i = 0; i < options->num_cpu_threads; i++) {
	       printf("\tThread %d: %0.4f s (%d batches, %d items) -> %0.2f CAL/s)\n", 
		      i, cal_time[i] / 1e6, thr_batches[i], thr_cal_items[i], 1e6 * thr_cal_items[i] / cal_time[i]);
	       total_item += thr_cal_items[i];
	       total_throughput += (1e6 * thr_cal_items[i] / cal_time[i]);
	       if (max_time < cal_time[i]) max_time = cal_time[i];
	  }
	  printf("\n\tTotal CALs: %lu, Max time = %0.4f, Throughput = %0.2f CAL/s\n", 
		 total_item, max_time / 1e6, total_throughput);
	  
	  total_item = 0; max_time = 0; total_throughput = 0;
	  printf("\nSW time:\n");
	  for (int i = 0; i < options->num_cpu_threads; i++) {
	       printf("\tThread %d: %0.4f s (%d batches, %d items) -> %0.2f SW/s)\n", 
		      i, sw_time[i] / 1e6, thr_batches[i], thr_sw_items[i], 1e6 * thr_sw_items[i] / sw_time[i]);
	       total_item += thr_sw_items[i];
	       total_throughput += (1e6 * thr_sw_items[i] / sw_time[i]);
	       if (max_time < sw_time[i]) max_time = sw_time[i];
	  }
	  printf("\n\tTotal SWs: %lu, Max time = %0.4f, Throughput = %0.2f SW/s\n", 
		 total_item, max_time / 1e6, total_throughput);
  }  
}
