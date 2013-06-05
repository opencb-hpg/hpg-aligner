#include "bs_aligner.h"

//--------------------------------------------------------------------
// run dna aligner
//--------------------------------------------------------------------

void run_bs_aligner(genome_t *genome1, genome_t *genome2, 
		    bwt_index_t *bwt_index1, bwt_index_t *bwt_index2, 
		    bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg, 
		    pair_mng_t *pair_mng, report_optarg_t *report_optarg, 
		    options_t *options) {

     int path_length = strlen(options->output_name);
     int prefix_length = 0;
     if (options->prefix_name) {
	  prefix_length = strlen(options->prefix_name);
     }
     
     char *reads_results = (char *) calloc((60 + prefix_length), sizeof(char));
     char *output_filename = (char *) calloc((path_length + prefix_length + 60), sizeof(char));
     
     if (options->prefix_name) {
	  strcat(reads_results, "/");
	  strcat(reads_results, options->prefix_name);
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
     batch_writer_input_init(output_filename, NULL, NULL, NULL, genome1, &writer_input);

     bam_header_t *bam_header = create_bam_header_by_genome(genome1);
     writer_input.bam_file = bam_fopen_mode(output_filename, bam_header, "w");
     bam_fwrite_header(bam_header, writer_input.bam_file);

     // preparing workflow stage input
     bwt_server_input_t bwt_input;
     bwt_server_input_init(NULL, 0, bwt_optarg, bwt_index1, 
			   NULL, 0, NULL, &bwt_input);
     
     region_seeker_input_t region_input;
     region_seeker_input_init(NULL, cal_optarg, bwt_optarg, 
			      bwt_index1, NULL, 0, options->gpu_process, 0, 0, 
			      &region_input);
     
     cal_seeker_input_t cal_input;
     cal_seeker_input_init(NULL, cal_optarg, NULL, 0, 
			   NULL, NULL, genome1, &cal_input);
     
     pair_server_input_t pair_input;
     pair_server_input_init(pair_mng, report_optarg, NULL, NULL, NULL, &pair_input);
     
     sw_server_input_t sw_input;
     sw_server_input_init(NULL, NULL, 0, options->match, options->mismatch, 
			  options->gap_open, options->gap_extend, options->min_score, 
			  options->flank_length, genome1, 0, 0, 0,  bwt_optarg, NULL, &sw_input);


     //--------------------------------------------------------------------------------------
     // workflow management
     //
     //
     // timing
     struct timeval start, end;
     extern double main_time;

     batch_t *batch = batch_new(&bwt_input, &region_input, &cal_input, 
				&pair_input, NULL, &sw_input, &writer_input, DNA_MODE, NULL);

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
     
     //if (time_on) {
     start_timer(start);
       //}

     workflow_run_with(options->num_cpu_threads, wf_input, wf);

     //if (time_on) {
     stop_timer(start, end, main_time);
       //printf("Total Time: %4.04f sec\n", time / 1000000);
       //}
     
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
     
     free(output_filename);

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
