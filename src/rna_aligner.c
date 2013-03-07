#include "rna_aligner.h"

#define NUM_SECTIONS_TIME 		8

//--------------------------------------------------------------------
// workflow input                                                                                                                                                     
//--------------------------------------------------------------------  


/*
buffer_item_t *buffer_item_new(fastq_read_t *fq_read, array_list_t *array_list) {
  buffer_item_t *item = (buffer_item_t *)malloc(sizeof(buffer_item_t));
  item->read = fq_read;
  item->alignments_list = array_list;

  return item;
}

buffer_pair_item_t *buffer_pair_item_new(fastq_read_t *fq_read_1, array_list_t *array_list_1, 
					 fastq_read_t *fq_read_2, array_list_t *array_list_2) {
  buffer_pair_item_t *item = (buffer_pair_item_t *)malloc(sizeof(buffer_pair_item_t));
  item->read_1 = fq_read_1;
  item->alignments_list_1 = array_list_1;

  item->read_2 = fq_read_2;
  item->alignments_list_2 = array_list_2;

  return item;
}

int convert_alignments_to_CAL(array_list_t *alignments_list, array_list_t *cal_list) {
  alignment_t *alginment;
  size_t num_items = array_list_size(alignments_list);
  int found = 0;

  for (int i = num_items - 1; i >= 0; i--) {
    alignment = array_list_get_at(i, alignments_list);
    if (alignment->large_hard_clipping) {
      found = 1;
      array_list_remove_at(i, alignments_list);
      //Extract strand, chr, start and search in AVL. Nex, we will form CALs and insert to cal_list
    } else {
      alignment->primary_alignment = 0;
    }
  }
  return found;
}

void thread_function(extra_stage_t *extra_stage_input) {
  linked_list_t *list = extra_stage_input->align_list;
  linked_list_iterator_init(align_list, itr);
  workflow_t *wf = extra_stage_input->workflow;
  pair_mng_t *pair_mng = extra_stage_input->pair_mng;
  batch_t *batch;
  size_t num_alignments;
  alignment_t *align;
  void *buffer_item;
  buffer_pair_item_t *item_pair;
  buffer_item_t *item;
  mapping_batch_t *mapping_batch = mapping_batch_new_by_num(MAX_READS_RNA + 1, pair_mng);
  array_list_t *fq_batch = array_list_new(MAX_READS_RNA, 
					  1.25f, 
					  COLLECTION_MODE_ASYNCHRONIZED);
  size_t num_reads = 0;
  size_t num_targets = 0;
  int found;
  global_status = WORKFLOW_STATUS_RUNNING;
  wf->complete_extra_stage = 0;

  printf("Extrem search STARTs...\n");
  
  while(workflow_get_simple_status(wf) == WORKFLOW_STATUS_RUNNING) {
    pthread_cond_wait(&cond_sp, &mutex_sp);
    while(buffer_item = linked_list_remove_last(list) {
	if (linked_list_get_flag(list) != SINGLE_END_MODE) {
	  buffer_item_t *item = (buffer_item_t *)buffer_item;
	  array_list_insert(item->read, fq_batch);
	  convert_alignments_to_CAL(item->alignment_list, mapping_batch->mapping_lists[num_reads]);
	  mapping_batch->old_mapping_list[num_reads] = item->alignment_list;
	  mapping_batch->targets[num_targets++] = 1;
	} else {
	  buffer_pair_item_t *item = (buffer_pair_item_t *)buffer_item;
	  array_list_insert(item->read_1, fq_batch); 
	  if (convert_alignments_to_CAL(item->alignment_list_1, mapping_batch->mapping_lists[num_reads])) {
	    mapping_batch->targets[num_targets++] = 1;
	  }
	  mapping_batch->old_mapping_list[num_reads++] = item->alignment_list_1;

	  array_list_insert(item->read_2, fq_batch);
	  if (convert_alignments_to_CAL(item->alignment_list_2, mapping_batch->mapping_lists[num_reads])) {
	    mapping_batch->targets[num_targets++] = 1;
	  }
	  mapping_batch->old_mapping_list[num_reads++] = item->alignment_list_2;
	}
	if (array_list_size(fq_batch) >= MAX_READS_RNA ) {
	  mapping_batch->num_targets = num_targets;
	  //Insert in the corrected site
	  mapping_batch = mapping_batch_new_by_num(MAX_READS_RNA + 1, pair_mng);
	}
      }
  }

  wf->complete_extra_stage = 1;
  
  printf("Finish search!\n");
}
*/
void run_rna_aligner(genome_t *genome, bwt_index_t *bwt_index, pair_mng_t *pair_mng,
		     bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg, 
		     options_t *options) {
  int path_length = strlen(options->output_name);
  int extend_length = 0;

  if (options->extend_name) {
    extend_length = strlen(options->extend_name);
  }

  char *reads_results = (char *)calloc((60 + extend_length), sizeof(char));
  char *extend_junctions = (char *)calloc((60 + extend_length), sizeof(char));
  char *exact_junctions = (char *)calloc((60 + extend_length), sizeof(char));

  char *output_filename = (char *)calloc((path_length + extend_length + 60), sizeof(char));
  char *extend_filename = (char *)calloc((path_length + extend_length + 60), sizeof(char));
  char *exact_filename = (char *)calloc((path_length + extend_length + 60), sizeof(char));

  if (options->extend_name) {
    strcat(reads_results, "/");
    strcat(reads_results, options->extend_name);
    strcat(reads_results, "_alignments.bam");  

    strcat(extend_junctions, "/");
    strcat(extend_junctions, options->extend_name);
    strcat(extend_junctions, "_extend_junctions.bed");

    strcat(exact_junctions, "/");
    strcat(exact_junctions, options->extend_name);
    strcat(exact_junctions, "_exact_junctions.bed");
 
  } else {
    strcat(reads_results, "/alignments.bam");
    strcat(extend_junctions, "/extend_junctions.bed");
    strcat(exact_junctions, "/exact_junctions.bed");
  } 

  strcat(output_filename, options->output_name);
  strcat(output_filename, reads_results);
  free(reads_results);

  strcat(extend_filename, options->output_name);
  strcat(extend_filename, extend_junctions);
  free(extend_junctions);

  strcat(exact_filename, options->output_name);
  strcat(exact_filename, exact_junctions);
  free(exact_junctions);


  LOG_DEBUG("Auto Thread Configuration Done !");

  // timing
  struct timeval start, end;
  double time;

  if (time_on) { 
    char* labels_time[NUM_SECTIONS_TIME] = {"FASTQ Reader               ", 
                                            "BWT Server                 ", 
					    "REGION Seeker              ", 
					    "CAL Seeker                 ", 
					    "RNA Preprocess             ", 
					    "RNA Server                 ",
					    "BAM Writer                 ", 
					    "TOTAL Time                 "};    
    timing = timing_new((char**) labels_time, NUM_SECTIONS_TIME);
  }

  // display selected options
  LOG_DEBUG("Displaying options...\n");
  options_display(options);

  //============================= INPUT INITIALIZATIONS =========================
  allocate_splice_elements_t chromosome_avls[genome->num_chromosomes];
  init_allocate_splice_elements(chromosome_avls, genome->num_chromosomes);
  load_intron_file(genome, options->intron_filename, chromosome_avls);


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

  linked_list_t *alignments_list = linked_list_new(COLLECTION_MODE_SYNCHRONIZED);
  linked_list_set_flag(options->pair_mode, alignments_list);

  bwt_server_input_t bwt_input;
  bwt_server_input_init(NULL,  0,  bwt_optarg, 
			bwt_index, NULL,  0, 
			NULL, &bwt_input);
      
  region_seeker_input_t region_input;
  region_seeker_input_init(NULL, cal_optarg, 
			   bwt_optarg, bwt_index, NULL, 
			   0, options->gpu_process, options->min_seed_padding_left, 
			   options->min_seed_padding_right, &region_input);
  
  cal_seeker_input_t cal_input;
  cal_seeker_input_init(NULL, cal_optarg, NULL, 0, NULL, NULL, genome, &cal_input);
  
  preprocess_rna_input_t preprocess_rna;
  
  preprocess_rna_input_init(options->max_intron_length, options->flank_length, 
			    options->seeds_max_distance, options->seed_size, genome, &preprocess_rna);

  
  sw_server_input_t sw_input;
  sw_server_input_init(NULL, NULL, 0,  options->match,  
		       options->mismatch,  options->gap_open, options->gap_extend,  
		       options->min_score,  options->flank_length, genome,  
		       options->max_intron_length, options->min_intron_length,  
		       options->seeds_max_distance,  bwt_optarg, chromosome_avls, &sw_input);
  

  pair_server_input_t pair_input;
  pair_server_input_init(pair_mng, bwt_optarg->report_best, bwt_optarg->report_n_hits, 
			 bwt_optarg->report_all, NULL, NULL, NULL, &pair_input);

  batch_writer_input_t writer_input;
  batch_writer_input_init( output_filename,
			   exact_filename, 
			   extend_filename, 
			   alignments_list, genome, &writer_input);

  bam_header_t *bam_header = create_bam_header_by_genome(genome);
  writer_input.bam_file = bam_fopen_mode(output_filename, bam_header, "w");
  bam_fwrite_header(bam_header, writer_input.bam_file);
  
  extra_stage_t extra_stage_input;
  //===================================================================================
  //-----------------------------------------------------------------------------------
  // workflow management
  //
  //
  batch_t *batch = batch_new(&bwt_input, &region_input, &cal_input, 
			     &pair_input, &preprocess_rna, &sw_input, &writer_input, RNA_MODE, NULL);

  wf_input_t *wf_input = wf_input_new(&reader_input, batch);
  
  //create and initialize workflow
  workflow_t *wf = workflow_new();
     
  workflow_stage_function_t stage_functions[] = {bwt_stage, seeding_stage, cal_stage, 
						 rna_preprocess_stage, sw_stage, post_pair_stage};

  char *stage_labels[] = {"BWT", "SEEDING", "CAL", "PRE PAIR", "SW", "POST PAIR"};
  workflow_set_stages(6, &stage_functions, stage_labels, wf);
     
  // optional producer and consumer functions
  workflow_set_producer(fastq_reader, "FastQ reader", wf);
  workflow_set_consumer(bam_writer, "BAM writer", wf);

  // Create new thread POSIX for search extra Splice Junctions
  //============================================================
  pthread_attr_t attr;
  pthread_t thread;
  void *status;
  int ret;

  extra_stage_input.align_list = alignments_list;
  extra_stage_input.workflow = wf;
  //extra_stage_input.pair_mng = pair_mng_new(pair_mng->pair_mode, pair_mng->min_distance, pair_mng->max_distance);
  
  if (time_on) {
    start_timer(start);
  }

  /*pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
     
  //Create thread extra Splice Junctions
  if (ret = pthread_create(&thread, &attr, thread_function, &extra_stage_input)) {
    printf("ERROR; return code from pthread_create() is %d\n", ret);
    exit(-1);
  }
  */
  //Run workflow
  workflow_run_with(options->num_cpu_threads, wf_input, wf);
  //printf("Finish workflow\n");
  pthread_cond_signal(&cond_sp);

  //Write chromosome avls
  write_chromosome_avls(chromosome_avls, NULL,  
			extend_filename, 
			exact_filename,
			options->write_size, genome->num_chromosomes);

  /*  pthread_attr_destroy(&attr);
  if (ret = pthread_join(thread, &status)) {
    printf("ERROR; return code from pthread_join() is %d\n", ret);
    exit(-1);
  } 
  */
  if (time_on) { 
    stop_timer(start, end, time);
    timing_add(time, TOTAL_TIME, timing);
  }
    
  //closing files
  if (options->pair_mode == SINGLE_END_MODE) {
    fastq_fclose(reader_input.fq_file1);
  } else {
    fastq_fclose(reader_input.fq_file1);
    fastq_fclose(reader_input.fq_file2);
  }

  bam_fclose(writer_input.bam_file);
    
  // free memory
  workflow_free(wf);
  wf_input_free(wf_input);
  batch_free(batch);

  //
  //
  // end of workflow management
  //--------------------------------------------------------------------------------------


  free(output_filename);
  free(exact_filename);
  free(extend_filename);
  linked_list_free(alignments_list, NULL);
}

//--------------------------------------------------------------------
