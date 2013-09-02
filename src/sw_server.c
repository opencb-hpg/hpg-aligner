#include "sw_server.h"

//====================================================================================
//  Input structure for Smith-Waterman server
//====================================================================================

void sw_server_input_init(list_t* sw_list, list_t* alignment_list, unsigned int write_size, 
			  float match, float mismatch, float gap_open, float gap_extend, 
			  float min_score, unsigned int flank_length, genome_t* genome, 
			  size_t max_intron_size, int min_intron_size, 
			  size_t seed_max_distance, bwt_optarg_t* bwt_optarg_p, 
			  allocate_splice_elements_t *chromosome_avls_p, 
			  sw_server_input_t* input) {
  
  input->sw_list_p = sw_list;
  input->alignment_list_p = alignment_list;
  input->write_size = write_size;
  input->genome_p = genome;
  input->max_intron_size = max_intron_size;
  input->min_intron_size = min_intron_size;
  input->seed_max_distance = seed_max_distance;
  input->bwt_optarg_p =  bwt_optarg_p; 

  // Smith-Waterman parameters
  input->match = match;
  input->mismatch = mismatch;
  input->gap_open = gap_open;
  input->gap_extend = gap_extend;
  input->min_score = min_score;

  input->sw_optarg.gap_open = gap_open;
  input->sw_optarg.gap_extend = gap_extend;

  input->sw_optarg.subst_matrix['A']['A'] = input->match;
  input->sw_optarg.subst_matrix['C']['A'] = input->mismatch;
  input->sw_optarg.subst_matrix['T']['A'] = input->mismatch;
  input->sw_optarg.subst_matrix['G']['A'] = input->mismatch;
  input->sw_optarg.subst_matrix['N']['A'] = input->mismatch;

  input->sw_optarg.subst_matrix['A']['C'] = input->mismatch;
  input->sw_optarg.subst_matrix['C']['C'] = input->match;
  input->sw_optarg.subst_matrix['T']['C'] = input->mismatch;
  input->sw_optarg.subst_matrix['G']['C'] = input->mismatch;
  input->sw_optarg.subst_matrix['N']['C'] = input->mismatch;

  input->sw_optarg.subst_matrix['A']['T'] = input->mismatch;
  input->sw_optarg.subst_matrix['C']['T'] = input->mismatch;
  input->sw_optarg.subst_matrix['T']['T'] = input->match;
  input->sw_optarg.subst_matrix['G']['T'] = input->mismatch;
  input->sw_optarg.subst_matrix['N']['T'] = input->mismatch;


  input->sw_optarg.subst_matrix['A']['G'] = input->mismatch;
  input->sw_optarg.subst_matrix['C']['G'] = input->mismatch;
  input->sw_optarg.subst_matrix['T']['G'] = input->mismatch;
  input->sw_optarg.subst_matrix['G']['G'] = input->match;
  input->sw_optarg.subst_matrix['N']['G'] = input->mismatch;

  input->sw_optarg.subst_matrix['A']['N'] = input->mismatch;
  input->sw_optarg.subst_matrix['C']['N'] = input->mismatch;
  input->sw_optarg.subst_matrix['T']['N'] = input->mismatch;
  input->sw_optarg.subst_matrix['G']['N'] = input->mismatch;
  input->sw_optarg.subst_matrix['N']['N'] = input->match;

  // CAL
  input->flank_length = flank_length;
  input->chromosome_avls_p = chromosome_avls_p;
}

//====================================================================================
//  Smith-Waterman server main function
//====================================================================================
/*
void sw_server(sw_server_input_t *input_p) {
     int sw_id = omp_get_thread_num();
     printf("sw_server (%i): START\n", sw_id);

     unsigned int i, num_reads, num_cals, bytes;
     char header_id[1024];
     // lists and items
     list_t *sw_list_p = input_p->sw_list_p;
     list_item_t *sw_item_p = NULL;
     sw_batch_t *sw_batch_p = NULL;
     
     list_item_t *item_p = NULL;
     
     list_t* write_list_p = input_p->write_list_p;
     unsigned int write_size = input_p->write_size;
     write_batch_t *write_batch_p = write_batch_new(write_size, MATCH_FLAG);
          
     // genome
     char *ref_p;
     unsigned int ref_len;
     unsigned long int start, end;
     genome_t *genome_p = input_p->genome_p;
     
     unsigned int flank_length = input_p->flank_length;
     
     // SIMD support for Smith-Waterman
     float min_score = input_p->min_score;
     unsigned int curr_depth = 0;
     sw_simd_input_t *sw_input_p = sw_simd_input_new(SIMD_DEPTH);
     sw_simd_output_t *sw_output_p = sw_simd_output_new(SIMD_DEPTH);
     
     sw_simd_context_t *context_p = sw_simd_context_new(input_p->match, input_p->mismatch, 
							input_p->gap_open, input_p->gap_extend); 
     
     // for tracking what reads are being mapped successfully
     unsigned int allocated_mapping_reads = 10000;
     unsigned char *mapping_reads_p = (unsigned char*) calloc(allocated_mapping_reads, sizeof(unsigned char));
     
     // for tracking the current read, cal being processed using sw_channel_t
     sw_channel_t *channel_p, *sw_channels_p = (sw_channel_t*) calloc(SIMD_DEPTH, sizeof(sw_channel_t));
     
     unsigned int header_len, read_len;
     unsigned int total = 0, total_valids = 0, total_reads = 0;
     char *header_id_match_p, *search_match_p, *quality_match_p, *cigar_p;
     unsigned int len_id_header;
     alignment_t *alignment_p;
     cal_t *cal_p;
     // main loop
     while ( (sw_item_p = list_remove_item(sw_list_p)) != NULL ) {

	  curr_depth = 0;
	  
	  if (timing_p != NULL) { timing_start(SW_SERVER, sw_id, timing_p); }
	  
	  sw_batch_p = (sw_batch_t*) sw_item_p->data_p;
	  num_reads = sw_batch_p->num_reads;
	  total_reads += num_reads;
	  
	  if (num_reads > allocated_mapping_reads) {
	       allocated_mapping_reads = num_reads;
	       mapping_reads_p = (unsigned char *) realloc(mapping_reads_p, 
							   allocated_mapping_reads * sizeof(unsigned char));
	  }
	  memset(mapping_reads_p, 0, allocated_mapping_reads * sizeof(unsigned char));
	  
	  // for each read
	  for (int i = 0; i < num_reads; i++) {
	       read_len = strlen(sw_batch_p->allocate_reads_p[i]->sequence);
	       header_len = strlen(sw_batch_p->allocate_reads_p[i]->id);

	       // for each cal
	       num_cals = array_list_size(sw_batch_p->allocate_cals_p[i]);
	       for(int j = 0; j < num_cals; j++) {
		    cal_p = array_list_get(j, sw_batch_p->allocate_cals_p[i]);
		    start = cal_p->start - flank_length;
		    end = cal_p->end + flank_length;
		    printf("Id:%s\n", sw_batch_p->allocate_reads_p[i]->id);
		    printf("\tProcess smith-waterman in coordenates: [%d - %d] - [chromosome %d] - [strand %d]\n", start, end, cal_p->chromosome_id, cal_p->strand);
		    channel_p = &sw_channels_p[curr_depth];
		    sw_channel_allocate_ref((unsigned int) end - start + 2, channel_p);
		    
		    genome_read_sequence_by_chr_index(channel_p->ref_p, cal_p->strand,
						      cal_p->chromosome_id - 1, &start, &end, genome_p);
		    
		    printf("Read to mapped: %s\n", sw_batch_p->allocate_reads_p[i]->sequence);
		    printf("Mapped to reference : %s\n", channel_p->ref_p);
		    
		    sw_channel_update(i, j, read_len, header_len, end - start + 1, channel_p);
		    
		    sw_simd_input_add(sw_batch_p->allocate_reads_p[i]->sequence, read_len,
				      channel_p->ref_p, channel_p->ref_len, 
				      curr_depth, sw_input_p);

		    // if depth is full, run SMID Smith-Waterman
		    if ((++curr_depth) == SIMD_DEPTH) {
			 smith_waterman_simd(sw_input_p, sw_output_p, context_p);
			 write_batch_p = process_sw_output(sw_output_p, sw_input_p, min_score, curr_depth, sw_channels_p, sw_batch_p, write_list_p, write_batch_p, write_size, sw_id, &total_valids, mapping_reads_p, genome_p);
			 
			 curr_depth = 0;
		    }
	       } // end of for 0..num_cals
	  } // end of for 0..num_reads
	  
	  if (curr_depth > 0) {
	       //printf("remaining smith-watermans = %i\n", curr_depth);
	       i = channel_p->read_index;
	       ref_p = channel_p->ref_p;
	       ref_len = channel_p->ref_len;
	       read_len = channel_p->read_len;
	       
	       for (int k = curr_depth; k < SIMD_DEPTH; k++) {
		    sw_simd_input_add(sw_batch_p->allocate_reads_p[i]->sequence, read_len,
				      ref_p, ref_len, 
				      k, sw_input_p);
	       }
	       smith_waterman_simd(sw_input_p, sw_output_p, context_p);
	       
	       write_batch_p = process_sw_output(sw_output_p, sw_input_p, min_score, curr_depth, sw_channels_p, sw_batch_p, write_list_p, write_batch_p, write_size, sw_id, &total_valids, mapping_reads_p, genome_p);
	       curr_depth = 0;
	  }
	  	  
	  for (i = 0; i < num_reads; i++) {
	       if (mapping_reads_p[i] == 0) {
		    printf("****************** read %i NO MAPPED !!!\n", i);
		    read_len = strlen(sw_batch_p->allocate_reads_p[i]->sequence);
		    header_len = strlen(sw_batch_p->allocate_reads_p[i]->id);

		    if ( write_batch_p->size > write_batch_p->allocated_size ) {
		      item_p = list_item_new(0, WRITE_ITEM, write_batch_p);
		      //if (time_on) { timing_stop(EXACT_SEEKER_INDEX, 0, timing_p); }
		      list_insert_item(item_p, write_list_p);
		      //if (time_on) { timing_start(EXACT_SEEKER_INDEX, 0, timing_p); }
		      
		      write_batch_p = write_batch_new(write_size, MATCH_FLAG);
		    }
		    
		    search_match_p = (char *)malloc(sizeof(char)*(read_len + 1));
		    memcpy(search_match_p, sw_batch_p->allocate_reads_p[i]->sequence, read_len);
		    search_match_p[read_len] = '\0';
		    
		    quality_match_p = (char *)malloc(sizeof(char)*(read_len + 1));
		    memcpy(quality_match_p, sw_batch_p->allocate_reads_p[i]->quality, read_len);
		    quality_match_p[read_len] = '\0';
		
		    len_id_header = get_to_first_blank(sw_batch_p->allocate_reads_p[i]->id, header_len, &header_id);
		    
		    header_id_match_p = (char *)malloc(sizeof(char)*len_id_header);
		    memcpy(header_id_match_p, &header_id, len_id_header);
		    
		    alignment_p = alignment_new();
		    cigar_p = (char *)malloc(sizeof(char)*10);
		    sprintf(cigar_p, "%d%c\0", read_len, 'X');
		    //TODO:chromosome 0??
		    
		    alignment_init_single_end(header_id_match_p, search_match_p, quality_match_p, 0, 0, 0, cigar_p, 1, 255, 0, 0, alignment_p);

		    //printf("seq: %s\n", alignment_p->sequence);
		    ((alignment_t **)write_batch_p->buffer_p)[write_batch_p->size] = alignment_p;
		    write_batch_p->size++;
		      
		    mapping_reads_p[i] = 2;
	       } // end of if mapping_reads_p[i] == 0
	  } // end of for 0..num_reads
	  
	  sw_batch_free(sw_batch_p);
	  list_item_free(sw_item_p);
	  
	  if (timing_p != NULL) { timing_stop(SW_SERVER, sw_id, timing_p); }
     } // end of while 
     
     // insert or free memory
     if (write_batch_p != NULL) {
	  if (write_batch_p->size > 0) {
	       item_p = list_item_new(0, WRITE_ITEM, write_batch_p);
	       list_insert_item(item_p, write_list_p);
	  } else {
	       write_batch_free(write_batch_p);
	  }
     }
     
     for(int k = 0; k < SIMD_DEPTH; k++) {
	  free(sw_channels_p[k].ref_p);
     }
     free(sw_channels_p);
     
     sw_simd_input_free(sw_input_p);
     sw_simd_output_free(sw_output_p);
     sw_simd_context_free(context_p); 
     
     list_decr_writers(write_list_p);
     
     printf("sw_server: END (%i reads, %i smith-waterman -> %i valids)\n", total_reads, total, total_valids);
}

//------------------------------------------------------------------------------------

write_batch_t* process_sw_output(sw_simd_output_t* sw_output_p, sw_simd_input_t* sw_input_p,
				 float min_score, unsigned int depth, sw_channel_t* sw_channels_p,
				 sw_batch_t* sw_batch_p, list_t* write_list_p, write_batch_t* write_batch_p, 
				 unsigned int write_size, unsigned int sw_id, unsigned int* total_valids_p,
				 unsigned char* mapping_reads_p, genome_t* genome_p) {

     unsigned int i, j;
     unsigned int header_len, read_len, mapped_len, bytes;
     short int primary_alignment;
     char header_id[1024];
     char *header_match_p, *read_match_p, *quality_match_p;
     char *cigar_p;
     short int number_cigar_op;
     alignment_t *alignment_p;
     unsigned int chromosome;
     unsigned int deletion_n;
     list_item_t* item_p = NULL;
     
     printf(" ======================== Process Output SW =========================\n");
     sw_simd_input_display(depth, sw_input_p);
     sw_simd_output_display(depth, sw_output_p);
     printf("======================================================================\n");
     for (int k = 0; k < depth; k++) {
	  
	  read_len = sw_channels_p[k].read_len;
	  mapped_len = sw_output_p->mapped_len_p[k];
	  
	  // is it a valid mapping ? && TODO:mapped_len >= read_len ??
	  if (sw_output_p->norm_score_p[k] > min_score ) {
	       
	       i = sw_channels_p[k].read_index;
	       j = sw_channels_p[k].cal_index;
	       
	       
	       (*total_valids_p)++;
	       
	       // process valid alignment: sam format and save into batch to write
	        header_len = sw_channels_p[k].header_len;
	        if ( write_batch_p->size > write_batch_p->allocated_size ) {
		  item_p = list_item_new(0, WRITE_ITEM, write_batch_p);
		  //if (time_on) { timing_stop(EXACT_SEEKER_INDEX, 0, timing_p); }
		  list_insert_item(item_p, write_list_p);
		  //if (time_on) { timing_start(EXACT_SEEKER_INDEX, 0, timing_p); }
		  
		  write_batch_p = write_batch_new(write_size, MATCH_FLAG);
		  
		}
		
		deletion_n = 0;
		for(int i = 0; i < mapped_len; i++){
		  if( sw_output_p->mapped_seq_p[k][i] == '-'){
		    deletion_n++;
		  }
		}
		
		if(mapping_reads_p[i] == 1){
		    primary_alignment = 1;
		}else{
		    primary_alignment = 0;
		}
		alignment_p = alignment_new();
		
		
		header_len = get_to_first_blank(sw_batch_p->allocate_reads_p[i]->id, header_len, header_id);
		header_match_p = (char *)malloc(sizeof(char)*(header_len + 1));
		memcpy(header_match_p, &header_id, header_len);
		header_match_p[header_len] = '\0';
		printf("%d deletion\n", deletion_n);
		
		
		read_match_p = (char *)malloc(sizeof(char)*(mapped_len + 1));
		memcpy(read_match_p, sw_batch_p->allocate_reads_p[i]->sequence + sw_output_p->start_seq_p[k], mapped_len - deletion_n);
		read_match_p[mapped_len - deletion_n ] = '\0';
		printf("%s\n", read_match_p);
		
		quality_match_p = (char *)malloc(sizeof(char)*(mapped_len + 1));
		memcpy(quality_match_p, sw_batch_p->allocate_reads_p[i]->quality + sw_output_p->start_seq_p[k], mapped_len - deletion_n);
		quality_match_p[mapped_len - deletion_n] = '\0';
		//printf("%s\n", quality_match_p);
		
		cigar_p =  generate_cigar_str(sw_output_p->mapped_seq_p[k], sw_output_p->mapped_ref_p[k], sw_output_p->start_seq_p[k], sw_input_p->seq_len_p[k], sw_output_p->mapped_len_p[k], &number_cigar_op);
		printf("Generate cigar end\n");
		printf("chromosome %d\n", ((cal_t *)array_list_get(j, sw_batch_p->allocate_cals_p[i]))->chromosome_id );
		printf("cigar_p(%d)=%s\n", number_cigar_op, cigar_p);
		
		
		alignment_init_single_end(header_match_p, read_match_p, quality_match_p, ((cal_t *)array_list_get(j, sw_batch_p->allocate_cals_p[i]))->strand, ((cal_t *)array_list_get(j, sw_batch_p->allocate_cals_p[i]))->chromosome_id - 1, ((cal_t *)array_list_get(j, sw_batch_p->allocate_cals_p[i]))->start, cigar_p, number_cigar_op, 255, primary_alignment, 1, alignment_p);
		//printf("seq: %s\n", alignment_p->sequence);
		
		((alignment_t **)write_batch_p->buffer_p)[write_batch_p->size] = alignment_p;
		write_batch_p->size++;
		mapping_reads_p[i] = 1;
	           
	  } // end of if norm_score
	  
	  // free mapped sequence and reference
	  free(sw_output_p->mapped_seq_p[k]);
	  free(sw_output_p->mapped_ref_p[k]);
     } // end of for 0..depth
     
     return write_batch_p;
}
*/
//====================================================================================
//  Smith-Waterman channel for SIMD implementation
//====================================================================================

inline void sw_channel_allocate_ref(unsigned int length, sw_channel_t* channel_p) {
     if (channel_p == NULL) return;
     
     if (channel_p->allocated_ref_size < length) {
	  if (channel_p->ref_p == NULL) {
	       channel_p->ref_p = (char*) calloc(length, sizeof(char));
	  } else {
	       free(channel_p->ref_p);
	       channel_p->ref_p = (char*) calloc(length, sizeof(char));
	  }
	  channel_p->allocated_ref_size = length;
     }
}

//------------------------------------------------------------------------------------

inline void sw_channel_update(size_t read_index, unsigned int cal_index, unsigned int read_len,
			      unsigned int header_len, unsigned int ref_len, sw_channel_t *channel_p) {
     channel_p->read_index = read_index;
     channel_p->cal_index = cal_index;
     channel_p->read_len = read_len;
     channel_p->header_len = header_len;
     channel_p->ref_len = ref_len;
}

//====================================================================================
// apply_sw
//====================================================================================
//int unmapped_by_score_counter[100];

//FILE *fd_ref = NULL, *fd_query = NULL;

int apply_sw(sw_server_input_t* input, batch_t *batch) {

  //  if (fd_ref == NULL) { fd_ref = fopen("sw_ref2.txt", "w"); }
  //  if (fd_query == NULL) { fd_query = fopen("sw_query2.txt", "w"); }

  //  printf("START: apply_sw\n"); 
  int tid = omp_get_thread_num();
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  cal_t *cal = NULL;
  array_list_t *cal_list = NULL, *mapping_list = NULL;//, *old_list = NULL, *new_list = NULL;

  array_list_t *fq_batch = mapping_batch->fq_batch;
  fastq_read_t *fq_read;

  size_t start, end;
  genome_t *genome = input->genome_p;
     
  size_t flank_length = input->flank_length;

  // SIMD support for Smith-Waterman
  float score, min_score = input->min_score;
  //  size_t curr_depth = 0;
  sw_output_t *sw_output;
  //  sw_simd_input_t *sw_sinput = sw_simd_input_new(SIMD_DEPTH);
  //  sw_simd_output_t *sw_soutput = sw_simd_output_new(SIMD_DEPTH);
  //sw_simd_context_t *context = sw_simd_context_new(input->match, input->mismatch, 
  //						    input->gap_open, input->gap_extend); 

  // for tracking the current read, cal being processed using sw_channel_t
  //sw_channel_t *channel;
  //sw_channel_t sw_channels[SIMD_DEPTH];
  //memset(sw_channels, 0, sizeof(sw_channels));
  
  //size_t header_len, read_len;
  //size_t strands[SIMD_DEPTH], chromosomes[SIMD_DEPTH], starts[SIMD_DEPTH];
  
  size_t read_index, num_cals;
  //  size_t total = 0, valids = 0;

  size_t num_targets = mapping_batch->num_targets;
  size_t new_num_targets = 0;

  size_t sw_total = mapping_batch->num_to_do;

  // set to zero
  mapping_batch->num_to_do = 0;

  /*
  // for all seqs pending to process !!
  size_t sw_total = 0;
  for (size_t i = 0; i < num_seqs; i++) {
    sw_total += array_list_size(batch->mapping_lists[batch->targets[i]]);
  }
  printf("number of sw to run: %d (vs num_done = %d)\n", sw_total, batch->num_done);
  */

  sw_optarg_t *sw_optarg = &input->sw_optarg;
    /*
  sw_optarg_t sw_optarg; //= sw_optarg_new(gap_open, gap_extend, matrix_filename);
  sw_optarg.gap_open = input->gap_open;
  sw_optarg.gap_extend = input->gap_extend;
  sw_optarg.subst_matrix['A']['A'] = input->match;    sw_optarg.subst_matrix['C']['A'] = input->mismatch; sw_optarg.subst_matrix['T']['A'] = input->mismatch; sw_optarg.subst_matrix['G']['A'] = input->mismatch;
  sw_optarg.subst_matrix['A']['C'] = input->mismatch; sw_optarg.subst_matrix['C']['C'] = input->match;    sw_optarg.subst_matrix['T']['C'] = input->mismatch; sw_optarg.subst_matrix['G']['C'] = input->mismatch;
  sw_optarg.subst_matrix['A']['G'] = input->mismatch; sw_optarg.subst_matrix['C']['T'] = input->mismatch; sw_optarg.subst_matrix['T']['T'] = input->match;    sw_optarg.subst_matrix['G']['T'] = input->mismatch;
  sw_optarg.subst_matrix['A']['T'] = input->mismatch; sw_optarg.subst_matrix['C']['G'] = input->mismatch; sw_optarg.subst_matrix['T']['G'] = input->mismatch; sw_optarg.subst_matrix['G']['G'] = input->match;
    */
  sw_multi_output_t *output = sw_multi_output_new(sw_total);
  char *q[sw_total], *r[sw_total];
  uint8_t strands[sw_total], chromosomes[sw_total];
  size_t starts[sw_total];
  size_t sw_count = 0, read_indices[sw_total];
  int read_len, ref_len, max_ref_len;

  // debugging: to kown how many reads are not mapped by SW score
  //  int unmapped_by_score[fq_batch->num_reads];
  //  memset(unmapped_by_score, 0, fq_batch->num_reads * sizeof(int));

  //  printf("num of sw to do: %i\n", sw_total);

  // initialize query and reference sequences to Smith-Waterman
  for (size_t i = 0; i < num_targets; i++) {
    //    printf("sw_server: target #%i of %i\n", i, num_seqs);

    read_index = mapping_batch->targets[i];
    fq_read = (fastq_read_t *) array_list_get(read_index, fq_batch);

    //    printf("sw_server: read #%i\n", read_index);

    cal_list = mapping_batch->mapping_lists[read_index];
    num_cals = array_list_size(cal_list);

    read_len = fq_read->length;
    //    max_ref_len = read_len + (read_len / 2);

    //    printf("sw_server: num_cals = %i cals\n", num_cals);

    // processing each CAL from this read
    for(size_t j = 0; j < num_cals; j++) {

      // get cal and read index
      cal = array_list_get(j, cal_list);
      read_indices[sw_count] = read_index;


      if (flank_length >= cal->start) {
        start = 0;
      } else {
        start = cal->start - flank_length;
      }

      end = cal->end + flank_length;
      if (end >= genome->chr_size[cal->chromosome_id - 1]) {
        end = genome->chr_size[cal->chromosome_id - 1] - 1;
      }

      ref_len = end - start + 2;
      //      if (ref_len < max_ref_len) {

	// query sequence, revcomp if necessary
	q[sw_count] = (char *) calloc((read_len + 1), sizeof(char));
	memcpy(q[sw_count], fq_read->sequence, read_len);
	if (cal->strand == 1) {
	  seq_reverse_complementary(q[sw_count], read_len);
	}
	//q[sw_count] = &(fq_batch->seq[fq_batch->data_indices[index]]);
	
	// reference sequence
	//printf("\tSW: %d.[chromosome:%d]-[strand:%d]-[start:%d, end:%d]\n", j, cal->chromosome_id, cal->strand, cal->start, cal->end);
	
	r[sw_count] = calloc(1, end - start + 2);
	genome_read_sequence_by_chr_index(r[sw_count], cal->strand,
					  cal->chromosome_id - 1, &start, &end, genome);
	
	// save some stuff, we'll use them after...
	strands[sw_count] = cal->strand;
	chromosomes[sw_count] = cal->chromosome_id;
	starts[sw_count] = start;
	
	//      if ( strchr(r[sw_count], 'N') == NULL && strlen(r[sw_count]) < (strlen(q[sw_count]) + (strlen(q[sw_count])/2)) ) {
	//	fprintf(fd_ref, "%s\n", r[sw_count]);
	//	fprintf(fd_query, "%s\n", q[sw_count]);
	//      }
	
	/*
	  printf("\tread #%i (sw #%i of %i):\n", index, sw_count, sw_total);
	  printf("\t\tquery: %s (%i)\n\t\tref  : %s (%i)\n\n", 
	  q[sw_count], strlen(q[sw_count]), r[sw_count], strlen(r[sw_count]));
	*/
	// increase counter
	sw_count++;
	//      } else {
	//	printf("ref_len = %i (max. %i)\n", ref_len, max_ref_len);
	//      }
    }

    // free cal_list
    array_list_clear(cal_list, (void *) cal_free);
    //    batch->mapping_lists[index] = NULL;
  }

  // run Smith-Waterman
  //  printf("before smith_waterman: sw_total = %i, sw_count = %i\n", sw_total, sw_count);
  smith_waterman_mqmr(q, r, sw_count, sw_optarg, 1, output);
  //  printf("after smith_waterman\n");

  /*
  // debugging
  {
    FILE *fd = fopen("sw.out", "w");
    sw_multi_output_save(sw_total, output, fd);
    fclose(fd);
  }
  */
  double norm_score;
  // filter alignments by min_score
  for (size_t i = 0; i < sw_count; i++) {

    //    score = output->score_p[i] / (strlen(output->query_map_p[i]) * input->match);
    //    if (score >= min_score) {
    /*
    printf("--------------------------------------------------------------\n");
    printf("Smith-Waterman results: read_indices = %i\n", read_indices[i]);
    //    printf("id\t%s\n", &(batch->fq_batch->header[batch->fq_batch->header_indices[read_indices[i]]]));
    printf("ref\n%s\n", r[i]);
    printf("query\n%s\n", q[i]);
    printf("map\n%s\n", output->ref_map_p[i]);
    printf("ref: chr = %d, strand = %d, start = %d, len = %d\n", chromosomes[i], strands[i], starts[i], strlen(r[i]));
    printf("query-map-start = %d, ref-map-start = %d\n", 
	   output->query_start_p[i], output->ref_start_p[i]);
    printf("score = %0.2f (min. score = %0.2f)\n", output->score_p[i], min_score);
    printf("--------------------------------------------------------------\n");
    */

    read_index = read_indices[i];
    fq_read = (fastq_read_t *) array_list_get(read_index, fq_batch);

    read_len = fq_read->length;
    norm_score = NORM_SCORE(output->score_p[i], read_len, input->match);

    /*
    LOG_DEBUG_F("ref. map : %s (start: %i)\n", 
		output->ref_map_p[i], output->ref_start_p[i]);
    LOG_DEBUG_F("query map: %s (start: %i)\n", 
		output->query_map_p[i], output->query_start_p[i]);
    LOG_DEBUG("\n");
    */

    if (norm_score >= min_score) {
      // valid mappings, 
      //insert in the list for further processing
      mapping_list = mapping_batch->mapping_lists[read_index];
      array_list_set_flag(0, mapping_list);

      if (array_list_size(mapping_list) == 0) {
      //	mapping_list = array_list_new(1000, 
      //				      1.25f, 
      //				      COLLECTION_MODE_ASYNCHRONIZED);
	
      //	batch->mapping_lists[index] = mapping_list;
	mapping_batch->targets[new_num_targets++] = read_index;
      }

      sw_output = sw_output_new(strands[i],
				chromosomes[i],
				starts[i],
				strlen(r[i]),
				strlen(output->query_map_p[i]),
				output->query_start_p[i],
				output->ref_start_p[i],
				output->score_p[i],
				norm_score,
				output->query_map_p[i],
				output->ref_map_p[i]);
      array_list_insert(sw_output, mapping_list);

      mapping_batch->num_to_do++;

      // debugging
      //unmapped_by_score[index] = 1;
    }

    // free query and reference
    free(q[i]);
    free(r[i]);
  }
  mapping_batch->num_targets = new_num_targets;
  /*
  // debugging
  for (size_t i = 0; i < fq_batch->num_reads; i++) {
    if (unmapped_by_score[i] == 0) {
	unmapped_by_score_counter[tid]++;
	//printf("by score: %s\n", &(batch->fq_batch->header[batch->fq_batch->header_indices[index]]));
      }
  }
  */

  // update counter
  //  thr_sw_items[tid] += sw_count;

  // free
  sw_multi_output_free(output);

  // go to the next stage

  return POST_PAIR_STAGE;

  //  printf("END: apply_sw, (%d Smith-Waterman, %d valids)\n", total, valids);
}

//--------------------------------------------------------------------------------------

int apply_sw_bs(sw_server_input_t* input, batch_t *batch) {

  int sw_3_nucleotides = 0;

  /*
  sw_optarg_t *sw_optarg2 = &input->sw_optarg;

  printf("Matrix Table\n\tA\tC\tG\tT\tN\nA\t%+.0f\t%+.0f\t%+.0f\t%+.0f\t%+.0f\nC\t%+.0f\t%+.0f\t%+.0f\t%+.0f\t%+.0f\nG\t%+.0f\t%+.0f\t%+.0f\t%+.0f\t%+.0f\nT\t%+.0f\t%+.0f\t%+.0f\t%+.0f\t%+.0f\nN\t%+.0f\t%+.0f\t%+.0f\t%+.0f\t%+.0f\n\n",
	 sw_optarg2->subst_matrix['A']['A'],
	 sw_optarg2->subst_matrix['C']['A'],
	 sw_optarg2->subst_matrix['G']['A'],
	 sw_optarg2->subst_matrix['T']['A'],
	 sw_optarg2->subst_matrix['N']['A'],

	 sw_optarg2->subst_matrix['A']['C'],
	 sw_optarg2->subst_matrix['C']['C'],
	 sw_optarg2->subst_matrix['G']['C'],
	 sw_optarg2->subst_matrix['T']['C'],
	 sw_optarg2->subst_matrix['N']['C'],

	 sw_optarg2->subst_matrix['A']['G'],
	 sw_optarg2->subst_matrix['C']['G'],
	 sw_optarg2->subst_matrix['G']['G'],
	 sw_optarg2->subst_matrix['T']['G'],
	 sw_optarg2->subst_matrix['N']['G'],

	 sw_optarg2->subst_matrix['A']['T'],
	 sw_optarg2->subst_matrix['C']['T'],
	 sw_optarg2->subst_matrix['G']['T'],
	 sw_optarg2->subst_matrix['T']['T'],
	 sw_optarg2->subst_matrix['N']['T'],

	 sw_optarg2->subst_matrix['A']['N'],
	 sw_optarg2->subst_matrix['C']['N'],
	 sw_optarg2->subst_matrix['G']['N'],
	 sw_optarg2->subst_matrix['T']['N'],
	 sw_optarg2->subst_matrix['N']['N']
	 );
  */

  if (sw_3_nucleotides == 0) {
    apply_sw_bs_4nt(input, batch);
  } else {
    
    //printf("START: apply_sw\n"); 
    int tid = omp_get_thread_num();
    mapping_batch_t *mapping_batch = batch->mapping_batch;
    cal_t *cal = NULL;
    array_list_t *cal_list = NULL, *mapping_list = NULL;
    
    array_list_t *fq_batch = mapping_batch->fq_batch;
    fastq_read_t *fq_read;
    
    // added by PP for bisulfite
    array_list_t *CT_fq_batch = mapping_batch->CT_fq_batch;
    array_list_t *GA_fq_batch = mapping_batch->GA_fq_batch;
    array_list_t *CT_rev_fq_batch = mapping_batch->CT_rev_fq_batch;
    array_list_t *GA_rev_fq_batch = mapping_batch->GA_rev_fq_batch;
    fastq_read_t *fq_read2;
    // end added by PP for bisulfite
    
    size_t start, end;
    size_t start2, end2;
    /*
    genome_t *genome = input->genome_p;
    */
    // added by PP for bisulfite
    genome_t *genome1 = input->genome1_p;
    genome_t *genome2 = input->genome2_p;
    // end added by PP for bisulfite
    
    size_t flank_length = input->flank_length;
    
    // SIMD support for Smith-Waterman
    float score, min_score = input->min_score;
    
    sw_output_t *sw_output;
    
    size_t read_index, num_cals;
    
    size_t num_targets = mapping_batch->num_targets;
    size_t new_num_targets = 0;
    // added by PP for bisulfite
    size_t num_targets2 = mapping_batch->num_targets2;
    size_t new_num_targets2 = 0;
    // added by PP for bisulfite
    
    // added by PP for bisulfite
    size_t sw_total1 = mapping_batch->num_to_do;
    size_t sw_total2 = mapping_batch->num_to_do2;
    size_t sw_total = sw_total1 + sw_total2;
    // end added by PP for bisulfite
    
    // set to zero
    mapping_batch->num_to_do = 0;
    // added by PP for bisulfite
    mapping_batch->num_to_do2 = 0;
    int g[sw_total];
    // end added by PP for bisulfite
    
    sw_optarg_t *sw_optarg = &input->sw_optarg;
    
    sw_multi_output_t *output = sw_multi_output_new(sw_total);
    char *q[sw_total], *r[sw_total];
    uint8_t strands[sw_total], chromosomes[sw_total];
    size_t starts[sw_total];
    size_t sw_count = 0, read_indices[sw_total], sw_count2 = 0;
    int read_len, ref_len, max_ref_len;
    
    //printf("num of sw to do: %i\n", sw_total);
    
    // initialize query and reference sequences to Smith-Waterman
    for (size_t i = 0; i < num_targets; i++) {
      //    printf("sw_server: target #%i of %i\n", i, num_seqs);
      read_index = mapping_batch->targets[i];
      
      // to use with the three nucleotides searches
      fq_read  = (fastq_read_t *) array_list_get(read_index, GA_fq_batch);
      fq_read2 = (fastq_read_t *) array_list_get(read_index, GA_rev_fq_batch);
      
      //printf("read %lu = %s\n", read_index, fq_read->sequence);
      //printf("read %lu = %s\n", read_index, fq_read2->sequence);
      
      //    printf("sw_server: read #%i\n", read_index);
      
      cal_list = mapping_batch->mapping_lists[read_index];
      num_cals = array_list_size(cal_list);
      
      read_len = fq_read->length;
      //    max_ref_len = read_len + (read_len / 2);
      
      //printf("sw_server: num_cals = %i cals\n", num_cals);
      
      // processing each CAL from this read
      for(size_t j = 0; j < num_cals; j++) {
	
	// get cal and read index
	cal = array_list_get(j, cal_list);
	read_indices[sw_count] = read_index;
	
	if (flank_length >= cal->start) {
	  start = 0;
	} else {
	  start = cal->start - flank_length;
	}
	
	end = cal->end + flank_length;
	if (end >= genome1->chr_size[cal->chromosome_id - 1]) {
	  end = genome1->chr_size[cal->chromosome_id - 1] - 1;
	}
	
	ref_len = end - start + 2;
	//      if (ref_len < max_ref_len) {
	
	// query sequence, revcomp if necessary
	q[sw_count] = (char *) calloc((read_len + 1), sizeof(char));
	
	// to use with the three nucleotides searches
	if (cal->strand == 0) {
	  memcpy(q[sw_count], fq_read->sequence, read_len);
	  //seq_reverse_complementary(q[sw_count], read_len);
	} else {
	  memcpy(q[sw_count], fq_read2->sequence, read_len);
	}
	
	//q[sw_count] = &(fq_batch->seq[fq_batch->data_indices[index]]);
	
	// reference sequence
	//printf("\tSW: %d.[chromosome:%d]-[strand:%d]-[start:%d, end:%d]\n", j, cal->chromosome_id, cal->strand, cal->start, cal->end);
	
	r[sw_count] = calloc(1, end - start + 2);
	
	// to use with the three nucleotides searches

	if (cal->strand == 0) {
	  genome_read_sequence_by_chr_index(r[sw_count], 0,
					    cal->chromosome_id - 1, &start, &end, genome1);
	} else {

	  genome_read_sequence_by_chr_index(r[sw_count], 0,
					    cal->chromosome_id - 1, &start, &end, genome2);

	  /*
	  start2 = genome1->chr_size[cal->chromosome_id - 1] - 1 - end;
	  end2   = genome1->chr_size[cal->chromosome_id - 1] - 1 - start;
	  genome_read_sequence_by_chr_index(r[sw_count], 0,
					    cal->chromosome_id - 1, &start2, &end2, genome2);
	  */
	}

	/*
	genome_read_sequence_by_chr_index(r[sw_count], cal->strand,
					  cal->chromosome_id - 1, &start, &end, genome1);
	*/
	
	// save some stuff, we'll use them after...
	strands[sw_count] = cal->strand;
	chromosomes[sw_count] = cal->chromosome_id;
	starts[sw_count] = start;

	/*
	printf("st = %lu\tend = %lu\n", cal->start, cal->end);
	printf("1\nseq %s\ngen %s\nstrand %2lu chromo %lu start %lu end %lu\n",
	       q[sw_count], r[sw_count], cal->strand, cal->chromosome_id, start, end);
	*/
	// increase counter
	sw_count++;
      }
      
      // free cal_list
      array_list_clear(cal_list, (void *) cal_free);
      //    batch->mapping_lists[index] = NULL;
    }
    ////////////////
    sw_count2 = sw_count;
    
    for (size_t i = 0; i < num_targets2; i++) {
      //    printf("sw_server: target #%i of %i\n", i, num_seqs);
      read_index = mapping_batch->targets2[i];
      
      // to use with the three nucleotides searches
      fq_read  = (fastq_read_t *) array_list_get(read_index, CT_fq_batch);
      fq_read2 = (fastq_read_t *) array_list_get(read_index, CT_rev_fq_batch);
      
      //printf("read %lu = %s\n", read_index, fq_read->sequence);
      //printf("read %lu = %s\n", read_index, fq_read2->sequence);
      
      //    printf("sw_server: read #%i\n", read_index);
      
      cal_list = mapping_batch->mapping_lists2[read_index];
      num_cals = array_list_size(cal_list);
      
      read_len = fq_read->length;
      //    max_ref_len = read_len + (read_len / 2);
      
      //printf("sw_server: num_cals = %i cals\n", num_cals);
      
      // processing each CAL from this read
      for(size_t j = 0; j < num_cals; j++) {
	
	// get cal and read index
	cal = array_list_get(j, cal_list);
	read_indices[sw_count] = read_index;
	
	if (flank_length >= cal->start) {
	  start = 0;
	} else {
	  start = cal->start - flank_length;
	}
	
	end = cal->end + flank_length;
	if (end >= genome1->chr_size[cal->chromosome_id - 1]) {
	  end = genome1->chr_size[cal->chromosome_id - 1] - 1;
	}
	
	ref_len = end - start + 2;
	//      if (ref_len < max_ref_len) {
	
	// query sequence, revcomp if necessary
	q[sw_count] = (char *) calloc((read_len + 1), sizeof(char));
	
	// to use with the three nucleotides searches
	if (cal->strand == 0) {
	  memcpy(q[sw_count], fq_read->sequence, read_len);
	  //seq_reverse_complementary(q[sw_count], read_len);
	} else {
	  memcpy(q[sw_count], fq_read2->sequence, read_len);
	}
	
	//q[sw_count] = &(fq_batch->seq[fq_batch->data_indices[index]]);
	
	// reference sequence
	//printf("\tSW: %d.[chromosome:%d]-[strand:%d]-[start:%d, end:%d]\n", j, cal->chromosome_id, cal->strand, cal->start, cal->end);
	
	r[sw_count] = calloc(1, end - start + 2);
	
	// to use with the three nucleotides searches

	if (cal->strand == 0) {
	  genome_read_sequence_by_chr_index(r[sw_count], 0,
					    cal->chromosome_id - 1, &start, &end, genome2);
	} else {

	  genome_read_sequence_by_chr_index(r[sw_count], 0,
					    cal->chromosome_id - 1, &start, &end, genome1);

	  /*
	  start2 = genome1->chr_size[cal->chromosome_id - 1] - 1 - end;
	  end2   = genome1->chr_size[cal->chromosome_id - 1] - 1 - start;
	  genome_read_sequence_by_chr_index(r[sw_count], 0,
					    cal->chromosome_id - 1, &start2, &end2, genome1);
	  */
	}

	/*
	genome_read_sequence_by_chr_index(r[sw_count], cal->strand,
					  cal->chromosome_id - 1, &start, &end, genome2);
	*/
	
	// save some stuff, we'll use them after...
	strands[sw_count] = cal->strand;
	chromosomes[sw_count] = cal->chromosome_id;
	starts[sw_count] = start;
	
	//printf("2\nseq %s\ngen %s\nstrand %2lu chromo %lu start %lu end %lu\n",
	//     q[sw_count], r[sw_count], cal->strand, cal->chromosome_id, start, end);

	// increase counter
	sw_count++;
      }
      
      // free cal_list
      array_list_clear(cal_list, (void *) cal_free);
      //    batch->mapping_lists[index] = NULL;
    }
    
    //printf("before smith_waterman: sw_total = %i, sw_count = %i, sw_count2 = %i\n", sw_total, sw_count, sw_count2);
    
    // run Smith-Waterman
    //  printf("before smith_waterman: sw_total = %i, sw_count = %i\n", sw_total, sw_count);
    smith_waterman_mqmr(q, r, sw_count, sw_optarg, 1, output);
    //  printf("after smith_waterman\n");
    
    
    for (size_t i = 0; i < sw_count; i++) {
      LOG_DEBUG_F("cal: start = %lu, strand = %i\n", starts[i], strands[i]);
      LOG_DEBUG_F("\tquery : %s\n", q[i]); 
      LOG_DEBUG_F("\tref.  : %s\n", r[i]); 
      LOG_DEBUG_F("\tquery map: %s (start: %i)\n", 
		  output->query_map_p[i], output->query_start_p[i]);
      LOG_DEBUG_F("\tref. map : %s (start: %i)\n", 
		  output->ref_map_p[i], output->ref_start_p[i]);
      LOG_DEBUG("\n");
    }
    
    
    //size_t mapp = 0, mapp2 = 0;
    
    
    
    double norm_score;
    // filter alignments by min_score
    for (size_t i = 0; i < sw_count2; i++) {
      
      read_index = read_indices[i];
      fq_read = (fastq_read_t *) array_list_get(read_index, GA_fq_batch);
      fq_read2 = (fastq_read_t *) array_list_get(read_index, GA_rev_fq_batch);
      
      read_len = fq_read->length;
      norm_score = NORM_SCORE(output->score_p[i], read_len, input->match);
      
      if (norm_score >= min_score) {
	// valid mappings, 
	//insert in the list for further processing
	mapping_list = mapping_batch->mapping_lists[read_index];
	array_list_set_flag(0, mapping_list);
	
	if (array_list_size(mapping_list) == 0) {
	  mapping_batch->targets[new_num_targets++] = read_index;
	  
	  //mapp++;
	}
	
	sw_output = sw_output_new(strands[i],
				  chromosomes[i],
				  starts[i],
				  strlen(r[i]),
				  strlen(output->query_map_p[i]),
				  output->query_start_p[i],
				  output->ref_start_p[i],
				  output->score_p[i],
				  norm_score,
				  output->query_map_p[i],
				  output->ref_map_p[i]);
	array_list_insert(sw_output, mapping_list);
	
	mapping_batch->num_to_do++;
	
      }
      
      // free query and reference
      free(q[i]);
      free(r[i]);
    }
    mapping_batch->num_targets = new_num_targets;
    
    for (size_t i = sw_count2; i < sw_count; i++) {
      
      read_index = read_indices[i];
      fq_read = (fastq_read_t *) array_list_get(read_index, CT_fq_batch);
      fq_read2 = (fastq_read_t *) array_list_get(read_index, CT_rev_fq_batch);
      
      read_len = fq_read->length;
      norm_score = NORM_SCORE(output->score_p[i], read_len, input->match);
      
      if (norm_score >= min_score) {
	// valid mappings, 
	//insert in the list for further processing
	mapping_list = mapping_batch->mapping_lists2[read_index];
	array_list_set_flag(0, mapping_list);
	
	if (array_list_size(mapping_list) == 0) {
	  mapping_batch->targets2[new_num_targets2++] = read_index;
	  
	  //mapp2++;
	}
	
	sw_output = sw_output_new(strands[i],
				  chromosomes[i],
				  starts[i],
				  strlen(r[i]),
				  strlen(output->query_map_p[i]),
				  output->query_start_p[i],
				  output->ref_start_p[i],
				  output->score_p[i],
				  norm_score,
				  output->query_map_p[i],
				  output->ref_map_p[i]);
	array_list_insert(sw_output, mapping_list);
	
	mapping_batch->num_to_do2++;
	
      }
      
      // free query and reference
      free(q[i]);
      free(r[i]);
    }
    mapping_batch->num_targets2 = new_num_targets2;
    
    // update counter
    //  thr_sw_items[tid] += sw_count;
    
    // free
    sw_multi_output_free(output);
    
    // go to the next stage
    
    /*
      printf("3 SW1         \t%3lu\tmapp               \t%3lu\tno map (discard) \t%3lu\n", 
      num_targets, mapp, num_targets - mapp);
      printf("3 SW2         \t%3lu\tmapp               \t%3lu\tno map (discard) \t%3lu\n", 
      num_targets2, mapp2, num_targets2 - mapp2);
    */
    
    //printf("END: apply_sw, (%d Smith-Waterman)\n", sw_total);
    
  }
  
  //return CONSUMER_STAGE;
  return POST_PAIR_STAGE;
  
  //  printf("END: apply_sw, (%d Smith-Waterman, %d valids)\n", total, valids);
}

//--------------------------------------------------------------------------------------

void fill_matrix(subst_matrix_t subst_matrix, float match, float mismatch, int type, float factor_match, float factor_mismatch) {

  subst_matrix['A']['A'] = match;
  subst_matrix['C']['A'] = mismatch;
  subst_matrix['T']['A'] = mismatch;
  subst_matrix['N']['A'] = mismatch;

  subst_matrix['A']['C'] = mismatch;
  subst_matrix['G']['C'] = mismatch;
  subst_matrix['N']['C'] = mismatch;

  subst_matrix['A']['T'] = mismatch;
  subst_matrix['T']['T'] = match;
  subst_matrix['G']['T'] = mismatch;
  subst_matrix['N']['T'] = mismatch;

  subst_matrix['C']['G'] = mismatch;
  subst_matrix['T']['G'] = mismatch;
  subst_matrix['N']['G'] = mismatch;

  subst_matrix['A']['N'] = mismatch;
  subst_matrix['C']['N'] = mismatch;
  subst_matrix['T']['N'] = mismatch;
  subst_matrix['G']['N'] = mismatch;
  subst_matrix['N']['N'] = match;

  /*
  float a = match * factor_match;
  float b = mismatch / factor_mismatch;
  float c = mismatch / factor_mismatch;
  */
  float x = 5;
  float y = 5;
  float z = 5;

  if (type == 1) {
    subst_matrix['C']['C'] = x;
    subst_matrix['T']['C'] = y;
    subst_matrix['C']['T'] = z;
    subst_matrix['G']['G'] = match;
    subst_matrix['A']['G'] = mismatch;
    subst_matrix['G']['A'] = mismatch;
  } else {
    subst_matrix['C']['C'] = match;
    subst_matrix['T']['C'] = mismatch;
    subst_matrix['C']['T'] = mismatch;
    subst_matrix['G']['G'] = x;
    subst_matrix['A']['G'] = y;
    subst_matrix['G']['A'] = z;
  }

}

//--------------------------------------------------------------------------------------

void apply_sw_bs_4nt(sw_server_input_t* input, batch_t *batch) {

  //printf("\nSTART: apply_sw\n"); 
  int tid = omp_get_thread_num();
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  cal_t *cal = NULL;
  array_list_t *cal_list = NULL, *mapping_list = NULL;//, *old_list = NULL, *new_list = NULL;

  array_list_t *fq_batch = mapping_batch->fq_batch;
  fastq_read_t *fq_read;

  size_t start, end;
  genome_t *genome = input->genome_p;
     
  size_t flank_length = input->flank_length;

  // SIMD support for Smith-Waterman
  float score, min_score = input->min_score;

  sw_output_t *sw_output;

  size_t read_index, num_cals;

  size_t num_targets = mapping_batch->num_targets;
  size_t new_num_targets = 0;
  size_t num_targets2 = mapping_batch->num_targets2;
  size_t new_num_targets2 = 0;

  size_t sw_total1 = mapping_batch->num_to_do;
  size_t sw_total2 = mapping_batch->num_to_do2;
  size_t sw_total = sw_total1 + sw_total2;

  // set to zero
  mapping_batch->num_to_do = 0;
  mapping_batch->num_to_do2 = 0;
  int g[sw_total];

  sw_optarg_t *sw_optarg = &input->sw_optarg;

  //sw_multi_output_t *output  = sw_multi_output_new(sw_total);
  sw_multi_output_t *output1 = sw_multi_output_new(sw_total);
  sw_multi_output_t *output2 = sw_multi_output_new(sw_total);

  //char *q[sw_total],  *r[sw_total];
  char *q1[sw_total], *r1[sw_total];
  char *q2[sw_total], *r2[sw_total];

  //uint8_t strands[sw_total], chromosomes[sw_total];
  uint8_t strands1[sw_total], chromosomes1[sw_total];
  uint8_t strands2[sw_total], chromosomes2[sw_total];

  //size_t starts[sw_total];
  size_t starts1[sw_total];
  size_t starts2[sw_total];

  size_t sw_count = 0, read_indices[sw_total], sw_count2 = 0;

  int read_len, ref_len, max_ref_len;

  //printf("matrix 1\n");
  //sw_optarg_t *sw_optarg1 = sw_optarg_new(sw_optarg->gap_open, sw_optarg->gap_extend, sw_optarg->subst_matrix_name);
  sw_optarg_t sw_optarg1;
  //printf("matrix 2\n");
  //sw_optarg_t *sw_optarg2 = sw_optarg_new(sw_optarg->gap_open, sw_optarg->gap_extend, sw_optarg->subst_matrix_name);
  sw_optarg_t sw_optarg2;
 
  //printf("data matrix 1\n");
  float match = sw_optarg->subst_matrix['A']['A'];
  //printf("data matrix 2\n");
  float missm = sw_optarg->subst_matrix['C']['A'];

  sw_optarg1.gap_open   = sw_optarg->gap_open;
  sw_optarg1.gap_extend = sw_optarg->gap_extend;
  sw_optarg2.gap_open   = sw_optarg->gap_open;
  sw_optarg2.gap_extend = sw_optarg->gap_extend;

  //printf("open   %f, %f, %f\n", sw_optarg->gap_open, sw_optarg1.gap_open, sw_optarg2.gap_open);
  //printf("extend %f, %f, %f\n", sw_optarg->gap_extend, sw_optarg1.gap_extend, sw_optarg2.gap_extend);

  //printf("fill matrix 1\n");
  fill_matrix(sw_optarg1.subst_matrix, match, missm, 0, 8, 2);
  //printf("fill matrix 2\n");
  fill_matrix(sw_optarg2.subst_matrix, match, missm, 1, 8, 2);

  /*
  printf("Matrix Table1\n\tA\tC\tG\tT\tN\nA\t%+.0f\t%+.0f\t%+.0f\t%+.0f\t%+.0f\nC\t%+.0f\t%+.0f\t%+.0f\t%+.0f\t%+.0f\nG\t%+.0f\t%+.0f\t%+.0f\t%+.0f\t%+.0f\nT\t%+.0f\t%+.0f\t%+.0f\t%+.0f\t%+.0f\nN\t%+.0f\t%+.0f\t%+.0f\t%+.0f\t%+.0f\n",
	 sw_optarg1.subst_matrix['A']['A'],
	 sw_optarg1.subst_matrix['C']['A'],
	 sw_optarg1.subst_matrix['G']['A'],
	 sw_optarg1.subst_matrix['T']['A'],
	 sw_optarg1.subst_matrix['N']['A'],

	 sw_optarg1.subst_matrix['A']['C'],
	 sw_optarg1.subst_matrix['C']['C'],
	 sw_optarg1.subst_matrix['G']['C'],
	 sw_optarg1.subst_matrix['T']['C'],
	 sw_optarg1.subst_matrix['N']['C'],

	 sw_optarg1.subst_matrix['A']['G'],
	 sw_optarg1.subst_matrix['C']['G'],
	 sw_optarg1.subst_matrix['G']['G'],
	 sw_optarg1.subst_matrix['T']['G'],
	 sw_optarg1.subst_matrix['N']['G'],

	 sw_optarg1.subst_matrix['A']['T'],
	 sw_optarg1.subst_matrix['C']['T'],
	 sw_optarg1.subst_matrix['G']['T'],
	 sw_optarg1.subst_matrix['T']['T'],
	 sw_optarg1.subst_matrix['N']['T'],

	 sw_optarg1.subst_matrix['A']['N'],
	 sw_optarg1.subst_matrix['C']['N'],
	 sw_optarg1.subst_matrix['G']['N'],
	 sw_optarg1.subst_matrix['T']['N'],
	 sw_optarg1.subst_matrix['N']['N']
	 );

  printf("Matrix Table2\n\tA\tC\tG\tT\tN\nA\t%+.0f\t%+.0f\t%+.0f\t%+.0f\t%+.0f\nC\t%+.0f\t%+.0f\t%+.0f\t%+.0f\t%+.0f\nG\t%+.0f\t%+.0f\t%+.0f\t%+.0f\t%+.0f\nT\t%+.0f\t%+.0f\t%+.0f\t%+.0f\t%+.0f\nN\t%+.0f\t%+.0f\t%+.0f\t%+.0f\t%+.0f\n\n",
	 sw_optarg2.subst_matrix['A']['A'],
	 sw_optarg2.subst_matrix['C']['A'],
	 sw_optarg2.subst_matrix['G']['A'],
	 sw_optarg2.subst_matrix['T']['A'],
	 sw_optarg2.subst_matrix['N']['A'],

	 sw_optarg2.subst_matrix['A']['C'],
	 sw_optarg2.subst_matrix['C']['C'],
	 sw_optarg2.subst_matrix['G']['C'],
	 sw_optarg2.subst_matrix['T']['C'],
	 sw_optarg2.subst_matrix['N']['C'],

	 sw_optarg2.subst_matrix['A']['G'],
	 sw_optarg2.subst_matrix['C']['G'],
	 sw_optarg2.subst_matrix['G']['G'],
	 sw_optarg2.subst_matrix['T']['G'],
	 sw_optarg2.subst_matrix['N']['G'],

	 sw_optarg2.subst_matrix['A']['T'],
	 sw_optarg2.subst_matrix['C']['T'],
	 sw_optarg2.subst_matrix['G']['T'],
	 sw_optarg2.subst_matrix['T']['T'],
	 sw_optarg2.subst_matrix['N']['T'],

	 sw_optarg2.subst_matrix['A']['N'],
	 sw_optarg2.subst_matrix['C']['N'],
	 sw_optarg2.subst_matrix['G']['N'],
	 sw_optarg2.subst_matrix['T']['N'],
	 sw_optarg2.subst_matrix['N']['N']
	 );
  */

  size_t elem_1 = 0, elem_2 = 0;
  
  //printf("num of sw to do: %lu\n", sw_total);

  // initialize query and reference sequences to Smith-Waterman
  for (size_t i = 0; i < num_targets; i++) {
    //    printf("sw_server: target #%i of %i\n", i, num_seqs);
    read_index = mapping_batch->targets[i];
    
    // to use with the four nucleotides searches
    fq_read  = (fastq_read_t *) array_list_get(read_index, fq_batch);
    
    //    printf("sw_server: read #%i\n", read_index);
    
    cal_list = mapping_batch->mapping_lists[read_index];
    num_cals = array_list_size(cal_list);
    
    read_len = fq_read->length;
    
    //printf("sw_server: num_cals = %i cals\n", num_cals);
    
    // processing each CAL from this read
    for(size_t j = 0; j < num_cals; j++) {
      
      // get cal and read index
      cal = array_list_get(j, cal_list);
      read_indices[sw_count] = read_index;

      if (flank_length >= cal->start) {
        start = 0;
      } else {
        start = cal->start - flank_length;
      }

      end = cal->end + flank_length;
      if (end >= genome->chr_size[cal->chromosome_id - 1]) {
        end = genome->chr_size[cal->chromosome_id - 1] - 1;
      }

      ref_len = end - start + 2;
      //      if (ref_len < max_ref_len) {
      
      if (cal->strand == 0) {
	g[sw_count] = 0;
	
	q1[elem_1] = (char *) calloc((read_len + 1), sizeof(char));
	
	// to use with the four nucleotides searches
	memcpy(q1[elem_1], fq_read->sequence, read_len);
	/*
	if (cal->strand == 1) {
	  seq_reverse_complementary(q1[elem_1], read_len);
	}
	*/
	
	r1[elem_1] = calloc(1, end - start + 2);
	
	// to use with the four nucleotides searches
	genome_read_sequence_by_chr_index(r1[elem_1], cal->strand,
					  cal->chromosome_id - 1, &start, &end, genome);
	// to use with the four nucleotides searches
	
	// save some stuff, we'll use them after...
	strands1[elem_1] = cal->strand;
	chromosomes1[elem_1] = cal->chromosome_id;
	starts1[elem_1] = start;

	//printf("11\nseq %s\ngen %s\nstrand %2lu chromo %lu start %lu end %lu\n",
	//     q1[elem_1], r1[elem_1], cal->strand, cal->chromosome_id, start, end);

	elem_1++;
      } else {
	g[sw_count] = 1;
	
	q2[elem_2] = (char *) calloc((read_len + 1), sizeof(char));
	
	// to use with the four nucleotides searches
	memcpy(q2[elem_2], fq_read->sequence, read_len);
	//	if (cal->strand == 1) {
	seq_reverse_complementary(q2[elem_2], read_len);
	//	}
	
	r2[elem_2] = calloc(1, end - start + 2);
	
	// to use with the four nucleotides searches
	genome_read_sequence_by_chr_index(r2[elem_2], cal->strand,
					  cal->chromosome_id - 1, &start, &end, genome);
	// to use with the four nucleotides searches
	
	// save some stuff, we'll use them after...
	strands2[elem_2] = cal->strand;
	chromosomes2[elem_2] = cal->chromosome_id;
	starts2[elem_2] = start;
	/*
	printf("st = %lu\tend = %lu\n", cal->start, cal->end);
	printf("12\nseq %s\ngen %s\nstrand %2lu chromo %lu start %lu end %lu\n",
	       q2[elem_2], r2[elem_2], cal->strand, cal->chromosome_id, start, end);
	*/
	elem_2++;
      }
      
      // query sequence, revcomp if necessary
      
      // increase counter
      sw_count++;
    }
    
    // free cal_list
    array_list_clear(cal_list, (void *) cal_free);
    //    batch->mapping_lists[index] = NULL;
  }
  ////////////////
  sw_count2 = sw_count;
  /*
  printf("num of sw to do from first list: %lu\n", sw_count2);
  printf("elem_1 %lu, elem_2 %lu\n", elem_1, elem_2);
  */

  for (size_t i = 0; i < num_targets2; i++) {
    //    printf("sw_server: target #%i of %i\n", i, num_seqs);
    read_index = mapping_batch->targets2[i];
    
    // to use with the four nucleotides searches
    fq_read  = (fastq_read_t *) array_list_get(read_index, fq_batch);
    
    //    printf("sw_server: read #%i\n", read_index);
    
    cal_list = mapping_batch->mapping_lists2[read_index];
    num_cals = array_list_size(cal_list);
    
    read_len = fq_read->length;
    //    max_ref_len = read_len + (read_len / 2);
    
    //printf("sw_server: num_cals = %i cals\n", num_cals);
    
    // processing each CAL from this read
    for(size_t j = 0; j < num_cals; j++) {
      
      // get cal and read index
      cal = array_list_get(j, cal_list);
      read_indices[sw_count] = read_index;
      
      if (flank_length >= cal->start) {
        start = 0;
      } else {
        start = cal->start - flank_length;
      }
      
      end = cal->end + flank_length;
      if (end >= genome->chr_size[cal->chromosome_id - 1]) {
        end = genome->chr_size[cal->chromosome_id - 1] - 1;
      }
      
      ref_len = end - start + 2;
      //      if (ref_len < max_ref_len) {
      
      if (cal->strand == 0) {
	g[sw_count] = 1;
	// query sequence, revcomp if necessary
	q2[elem_2] = (char *) calloc((read_len + 1), sizeof(char));
	
	// to use with the four nucleotides searches
	memcpy(q2[elem_2], fq_read->sequence, read_len);
	/*
	if (cal->strand == 1) {
	  seq_reverse_complementary(q2[elem_2], read_len);
	}
	*/
	
	r2[elem_2] = calloc(1, end - start + 2);
	
	// to use with the four nucleotides searches
	genome_read_sequence_by_chr_index(r2[elem_2], cal->strand,
					  cal->chromosome_id - 1, &start, &end, genome);
	// to use with the four nucleotides searches
	
	// save some stuff, we'll use them after...
	strands2[elem_2] = cal->strand;
	chromosomes2[elem_2] = cal->chromosome_id;
	starts2[elem_2] = start;

	//printf("22\nseq %s\ngen %s\nstrand %2lu chromo %lu start %lu end %lu\n",
	//      q2[elem_2], r2[elem_2], cal->strand, cal->chromosome_id, start, end);

	elem_2++;
      } else {
	g[sw_count] = 0;
	// query sequence, revcomp if necessary
	q1[elem_1] = (char *) calloc((read_len + 1), sizeof(char));
	
	// to use with the four nucleotides searches
	memcpy(q1[elem_1], fq_read->sequence, read_len);
	//	if (cal->strand == 1) {
	seq_reverse_complementary(q1[elem_1], read_len);
	//	}
	
	r1[elem_1] = calloc(1, end - start + 2);
	
	// to use with the four nucleotides searches
	genome_read_sequence_by_chr_index(r1[elem_1], cal->strand,
					  cal->chromosome_id - 1, &start, &end, genome);
	// to use with the four nucleotides searches
	
	// save some stuff, we'll use them after...
	strands1[elem_1] = cal->strand;
	chromosomes1[elem_1] = cal->chromosome_id;
	starts1[elem_1] = start;

	//printf("21\nseq %s\ngen %s\nstrand %2lu chromo %lu start %lu end %lu\n",
	//     q1[elem_1], r1[elem_1], cal->strand, cal->chromosome_id, start, end);

	elem_1++;
      }
      
      
      // increase counter
      sw_count++;
    }
    
    // free cal_list
    array_list_clear(cal_list, (void *) cal_free);
    //    batch->mapping_lists[index] = NULL;
  }
  /*
  printf("num of sw to do from second list %lu\n", sw_count - sw_count2);
  printf("num of sw to do %lu\n", sw_count);
  printf("elem_1 %lu, elem_2 %lu\n", elem_1, elem_2);
  */

  //printf("before smith_waterman: sw_total = %i, sw_count = %i, sw_count2 = %i\n", sw_total, sw_count, sw_count2);
  
  // run Smith-Waterman
  //  printf("before smith_waterman: sw_total = %i, sw_count = %i\n", sw_total, sw_count);
  //smith_waterman_mqmr(q, r, sw_count, sw_optarg, 1, output);
  if (elem_1 > 0)
    smith_waterman_mqmr(q1, r1, elem_1, &sw_optarg1, 1, output1);
  if (elem_2 > 0)
    smith_waterman_mqmr(q2, r2, elem_2, &sw_optarg2, 1, output2);
  //printf("after smith_waterman\n");
  
  /*
  for (size_t i = 0; i < sw_count; i++) {
    LOG_DEBUG_F("cal: start = %lu, strand = %i\n", starts[i], strands[i]);
    LOG_DEBUG_F("\tquery : %s\n", q[i]); 
    LOG_DEBUG_F("\tref.  : %s\n", r[i]); 
    LOG_DEBUG_F("\tquery map: %s (start: %i)\n", 
		output->query_map_p[i], output->query_start_p[i]);
    LOG_DEBUG_F("\tref. map : %s (start: %i)\n", 
		output->ref_map_p[i], output->ref_start_p[i]);
    LOG_DEBUG("\n");
  }
  */
  
  //size_t mapp = 0, mapp2 = 0;
  
  double norm_score;
  // filter alignments by min_score
  
  size_t el_1 = 0, el_2 = 0;
  
  for (size_t i = 0; i < sw_count2; i++) {
    
    read_index = read_indices[i];
    fq_read = (fastq_read_t *) array_list_get(read_index, fq_batch);
    
    read_len = fq_read->length;
    
    if (g[i] == 0) {
      norm_score = NORM_SCORE(output1->score_p[el_1], read_len, input->match);
      
      if (norm_score >= min_score) {
	// valid mappings, 
	//insert in the list for further processing
	mapping_list = mapping_batch->mapping_lists[read_index];
	array_list_set_flag(0, mapping_list);
	
	if (array_list_size(mapping_list) == 0) {
	  mapping_batch->targets[new_num_targets++] = read_index;
	  
	  //mapp++;
	}
	
	sw_output = sw_output_new(strands1[el_1],
				  chromosomes1[el_1],
				  starts1[el_1],
				  strlen(r1[el_1]),
				  strlen(output1->query_map_p[el_1]),
				  output1->query_start_p[el_1],
				  output1->ref_start_p[el_1],
				  output1->score_p[el_1],
				  norm_score,
				  output1->query_map_p[el_1],
				  output1->ref_map_p[el_1]);
	array_list_insert(sw_output, mapping_list);
	
	mapping_batch->num_to_do++;
      }
      // free query and reference
      free(q1[el_1]);
      free(r1[el_1]);
      
      el_1++;
    } else {
      norm_score = NORM_SCORE(output2->score_p[el_2], read_len, input->match);
      
      if (norm_score >= min_score) {
	// valid mappings, 
	//insert in the list for further processing
	mapping_list = mapping_batch->mapping_lists[read_index];
	array_list_set_flag(0, mapping_list);
	
	if (array_list_size(mapping_list) == 0) {
	  mapping_batch->targets[new_num_targets++] = read_index;
	  
	  //mapp++;
	}
	
	sw_output = sw_output_new(strands2[el_2],
				  chromosomes2[el_2],
				  starts2[el_2],
				  strlen(r2[el_2]),
				  strlen(output2->query_map_p[el_2]),
				  output2->query_start_p[el_2],
				  output2->ref_start_p[el_2],
				  output2->score_p[el_2],
				  norm_score,
				  output2->query_map_p[el_2],
				  output2->ref_map_p[el_2]);
	array_list_insert(sw_output, mapping_list);
	
	mapping_batch->num_to_do++;
      }
      // free query and reference
      free(q2[el_2]);
      free(r2[el_2]);
      
      el_2++;
    }
    
  }
  mapping_batch->num_targets = new_num_targets;

  /*
  printf("elem from first list:  %lu\n", el_1);
  printf("elem from second list: %lu\n", el_2);
  printf("num of sw do from first list: %lu\n", sw_count2);

  printf("from sw_count2: %lu\tto sw_count: %lu\n", sw_count2, sw_count);
  */
  
  for (size_t i = sw_count2; i < sw_count; i++) {
    
    read_index = read_indices[i];
    fq_read = (fastq_read_t *) array_list_get(read_index, fq_batch);
    
    //printf("g[%lu]: %i\n", i, g[i]);
    
    if (g[i] == 0) {
      
      read_len = fq_read->length;
      norm_score = NORM_SCORE(output1->score_p[el_1], read_len, input->match);
      
      if (norm_score >= min_score) {
	// valid mappings, 
	//insert in the list for further processing
	mapping_list = mapping_batch->mapping_lists2[read_index];
	array_list_set_flag(0, mapping_list);
	
	if (array_list_size(mapping_list) == 0) {
	  mapping_batch->targets2[new_num_targets2++] = read_index;
	  
	  //mapp2++;
	}
	
	sw_output = sw_output_new(strands1[el_1],
				  chromosomes1[el_1],
				  starts1[el_1],
				  strlen(r1[el_1]),
				  strlen(output1->query_map_p[el_1]),
				  output1->query_start_p[el_1],
				  output1->ref_start_p[el_1],
				  output1->score_p[el_1],
				  norm_score,
				  output1->query_map_p[el_1],
				  output1->ref_map_p[el_1]);
	array_list_insert(sw_output, mapping_list);
	
	mapping_batch->num_to_do2++;
	
      }
      
      // free query and reference
      free(q1[el_1]);
      free(r1[el_1]);
      el_1++;
    } else {
      //printf("g[i] = 1\n");
      //printf("elem_2 %lu, el_2 %lu\n", elem_2, el_2);
      
      read_len = fq_read->length;
      norm_score = NORM_SCORE(output2->score_p[el_2], read_len, input->match);
      
      //printf("g[i] = 1\n");

      if (norm_score >= min_score) {
	// valid mappings, 
	//insert in the list for further processing
	mapping_list = mapping_batch->mapping_lists2[read_index];
	array_list_set_flag(0, mapping_list);
	
	//printf("g[i] = 1\n");

	if (array_list_size(mapping_list) == 0) {
	  mapping_batch->targets2[new_num_targets2++] = read_index;
	  
	  //mapp2++;
	}
	
	sw_output = sw_output_new(strands2[el_2],
				  chromosomes2[el_2],
				  starts2[el_2],
				  strlen(r2[el_2]),
				  strlen(output2->query_map_p[el_2]),
				  output2->query_start_p[el_2],
				  output2->ref_start_p[el_2],
				  output2->score_p[el_2],
				  norm_score,
				  output2->query_map_p[el_2],
				  output2->ref_map_p[el_2]);
	//printf("g[i] = 1\n");
	array_list_insert(sw_output, mapping_list);
	//printf("g[i] = 1\n");
	
	mapping_batch->num_to_do2++;
	
      }
      
      // free query and reference
      free(q2[el_2]);
      free(r2[el_2]);
      el_2++;
    }
  }
  mapping_batch->num_targets2 = new_num_targets2;

  /*  
  printf("elem from first list:  %lu\n", el_1);
  printf("elem from second list: %lu\n", el_2);
  printf("num of sw do from second list: %lu\n", sw_count - sw_count2);
  */

  // update counter
  //  thr_sw_items[tid] += sw_count;
  
  // free
  //sw_multi_output_free(output);
  sw_multi_output_free(output1);
  sw_multi_output_free(output2);
  
  // go to the next stage
  //printf("END: apply_sw, (%d Smith-Waterman)\n", sw_total);
  //  printf("END: apply_sw, (%d Smith-Waterman, %d valids)\n", total, valids);
}

//--------------------------------------------------------------------------------------

sw_output_t *sw_output_new(int strand, size_t chrom, size_t ref_start, size_t ref_len,
			   size_t mref_len, size_t mquery_start, size_t mref_start,
			   float score, float norm_score, char* mquery, char* mref) {

  sw_output_t *p = (sw_output_t *) calloc(1, sizeof(sw_output_t));

  p->strand = strand;
  p->chromosome = chrom;
  p->ref_start = ref_start;
  p->ref_len = ref_len;
  p->mref_len = mref_len;
  p->mquery_start = mquery_start;
  p->mref_start = mref_start;
  p->score = score;
  p->norm_score = norm_score;
  p->mquery = strdup(mquery);
  p->mref = strdup(mref);
  //p->mquery = NULL;
  //p->mref = NULL;

  return p;
}

//--------------------------------------------------------------------------------------

void sw_output_free(sw_output_t *p) {
  if (p == NULL) return;

  if (p->mquery != NULL) free(p->mquery);
  if (p->mref != NULL) free(p->mref);

  free(p);
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

