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
  float score, norm_score, min_score = input->min_score;
  
  size_t read_index, num_cals;

  size_t num_targets = mapping_batch->num_targets;
  size_t new_num_targets = 0;

  size_t sw_total = mapping_batch->num_to_do;

  // set to zero
  mapping_batch->num_to_do = 0;

  sw_optarg_t *sw_optarg = &input->sw_optarg;
  sw_multi_output_t *output = sw_multi_output_new(sw_total);
  char *q[sw_total], *r[sw_total];
  uint8_t strands[sw_total], chromosomes[sw_total];
  size_t starts[sw_total];
  size_t sw_count = 0, read_indices[sw_total];
  int read_len, ref_len, max_ref_len, read_start, pos, gap_len;

  seed_region_t *s, *new_s;
  linked_list_iterator_t* itr;

  // initialize query and reference sequences to Smith-Waterman
  for (size_t i = 0; i < num_targets; i++) {

    read_index = mapping_batch->targets[i];
    fq_read = (fastq_read_t *) array_list_get(read_index, fq_batch);

    cal_list = mapping_batch->mapping_lists[read_index];
    num_cals = array_list_size(cal_list);

    read_len = fq_read->length;

    // processing each CAL from this read
    for(size_t j = 0; j < num_cals; j++) {

      // get cal and read index
      cal = array_list_get(j, cal_list);
      read_indices[sw_count] = read_index;

      read_start = 1;
      itr = linked_list_iterator_new(cal->sr_list);
      s = (seed_region_t *) linked_list_iterator_curr(itr);
      while (s != NULL) {
	LOG_DEBUG_F("[%i|%i - %i|%i] (read_start %i)\n", 
		    s->genome_start, s->read_start, s->read_end, s->genome_end, read_start);
	if (s->read_start != read_start) {
	  // prepare sw
	  // query sequence, revcomp if necessary
	  int gap_len = s->read_start - read_start + 1;
	  LOG_DEBUG_F("\tgap-length = %i\n", gap_len);
	  
	  q[sw_count] = (char *) calloc(gap_len, sizeof(char));
	  memcpy(q[sw_count], fq_read->sequence + read_start - 1, gap_len);
	  LOG_DEBUG_F("\tgap-seq = %s\n", q[sw_count]);
	  if (cal->strand == 1) {
	    seq_reverse_complementary(q[sw_count], gap_len);
	    LOG_DEBUG_F("\tgap-seq-rev = %s\n", q[sw_count]);
	  }
	
	  // reference sequence
	  if (read_start == 1) {
	    end = s->genome_end;
	    start = end - gap_len + 1;
	  } else {
	    start = pos;
	    end = start + gap_len - 1;
	  }
	  r[sw_count] = calloc(1, gap_len);
	  genome_read_sequence_by_chr_index(r[sw_count], cal->strand,
					    cal->chromosome_id - 1, &start, &end, genome);
	  LOG_DEBUG_F("\tref-seq = (%lu - %lu) : %s\n", start, end, r[sw_count]);
	  
	  // save some stuff, we'll use them after...
	  strands[sw_count] = cal->strand;
	  chromosomes[sw_count] = cal->chromosome_id;
	  starts[sw_count] = start;
	  
	  // increase counter
	  sw_count++;	  

	  // insert gap in the list
	  new_s = seed_region_new(read_start, s->read_end - 1, 0, 0, 0);
	  new_s->gap = 1;
	  linked_list_iterator_insert(new_s, itr);
	  linked_list_iterator_prev(itr);
	}
	read_start = s->read_end + 1;
	pos = s->genome_end + 1;

	//continue loop...
	linked_list_iterator_next(itr);
	s = linked_list_iterator_curr(itr);
      }
      if (read_start < read_len) { 
	// prepare sw
	// query sequence, revcomp if necessary
	int gap_len = read_len - read_start + 1;
	LOG_DEBUG_F("seq : %s\n", fq_read->sequence);
	LOG_DEBUG_F("final gap-length = %i + %i\n", gap_len, (2 * flank_length));
	
	q[sw_count] = (char *) calloc(gap_len + flank_length, sizeof(char));
	memcpy(q[sw_count], fq_read->sequence + read_start - 1 - flank_length, gap_len + flank_length);
	LOG_DEBUG_F("final gap-seq     = %s\n", q[sw_count]);
	if (cal->strand == 1) {
	  seq_reverse_complementary(q[sw_count], gap_len + flank_length);
	  LOG_DEBUG_F("final gap-seq-rev = %s\n", q[sw_count]);
	}
	
	// reference sequence
	start = pos - flank_length;
	end = start + gap_len - 1 + flank_length + flank_length;
	r[sw_count] = calloc(1, gap_len + flank_length + flank_length);
	genome_read_sequence_by_chr_index(r[sw_count], cal->strand,
					  cal->chromosome_id - 1, &start, &end, genome);
	LOG_DEBUG_F("final ref-seq = (%lu - %lu), gap_len = %i + %i: %s\n", 
		    start, end, gap_len, (2 * flank_length), r[sw_count]);
	
	// save some stuff, we'll use them after...
	strands[sw_count] = cal->strand;
	chromosomes[sw_count] = cal->chromosome_id;
	starts[sw_count] = start;
	
	// increase counter
	sw_count++;

	// insert gap at the end of the list
	new_s = seed_region_new(read_start - flank_length - 1, read_start + gap_len - 1, 0, 0, 0);
	new_s->gap = 1;
	linked_list_insert_last(new_s, cal->sr_list);
      }
    }
  }

  assert(sw_total == sw_count);
  
  // run Smith-Waterman
  smith_waterman_mqmr(q, r, sw_count, sw_optarg, 1, output);
  LOG_DEBUG_F("\tmquery: %s (start %i)\n", output->query_map_p[0], output->query_start_p[0]);
  LOG_DEBUG_F("\tmref  : %s (start %i)\n", output->ref_map_p[0], output->ref_start_p[0]);
  LOG_DEBUG_F("\tscore : %0.2f\n", output->score_p[0]);

  // re-construct the CAL score from region gaps
  sw_count = 0;
  int header_len, mapped_len, mquery_start, mref_start;
  array_list_t *alignment_list;
  size_t deletion_n, num_cigar_ops, shift;
  char *cigar, *header_match, *read_match, *quality_match, *quality, *header_clean;
  char *p, *optional_fields;
  alignment_t *alignment;
  int c, optional_fields_length, distance, AS, first;

  for (size_t i = 0; i < num_targets; i++) {

    LOG_DEBUG_F("targets %i of %i\n", i, num_targets);

    read_index = mapping_batch->targets[i];
    fq_read = (fastq_read_t *) array_list_get(read_index, fq_batch);
    
    cal_list = mapping_batch->mapping_lists[read_index];
    num_cals = array_list_size(cal_list);
    
    LOG_DEBUG_F("read_index =  %i, num_cals = %i\n", read_index, num_cals);

    //    read_len = fq_read->length;
    //    header_len = get_to_first_blank(fq_read->id, strlen(fq_read->id), header_clean);

    alignment_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    array_list_set_flag(1, alignment_list);

    // processing each CAL from this read
    score = 0.0f;
    for(size_t j = 0; j < num_cals; j++) {
      // get cal and read index
      cal = array_list_get(j, cal_list);

      LOG_DEBUG_F("cal %i of %i\n", j, num_cals);

      first = 1;
      for (linked_list_item_t *list_item = cal->sr_list->first; list_item != NULL; list_item = list_item->next) {
	s = list_item->item;
	LOG_DEBUG_F("is gap ? %s -> [%i|%i - %i|%i]\n", 
		    (s->gap ? "yes" : "no"), s->genome_start, s->read_start, s->read_end, s->genome_end);
	if (s->gap) {	  
	  LOG_DEBUG_F("\tmquery: %s (start %i)\n", output->query_map_p[sw_count], output->query_start_p[sw_count]);
	  LOG_DEBUG_F("\tmref  : %s (start %i)\n", output->ref_map_p[sw_count], output->ref_start_p[sw_count]);
	  LOG_DEBUG_F("\tscore : %0.2f\n", output->score_p[sw_count]);
	  score += output->score_p[sw_count];
	  if (first) {
	    mquery_start = output->query_start_p[sw_count];
	    mref_start = output->ref_start_p[sw_count];
	  }
	  // free query and reference
	  free(q[sw_count]);
	  free(r[sw_count]);
	  sw_count++;
	} else {
	  score += ((s->read_end - s->read_start + 1) * input->match);
	  if (first) {
	    mquery_start = s->read_start - 1;
	    mref_start = mquery_start;
	  }
	}
	first = 0;
      }
      LOG_DEBUG_F("score: %0.2f, mquery_start = %i, mref_start = %i\n",
		  score, mquery_start, mref_start);

      exit(-1);

      // filter by SW score
      norm_score = NORM_SCORE(score, read_len, input->match);
      
      if (norm_score >= min_score) {
	/*
	header_match = (char *) malloc(sizeof(char) * (header_len + 1));
	memcpy(header_match, header_clean, header_len);
	header_match[header_len] = '\0';
      
	mquery_start = sw_output->mquery_start;
	mref_start = sw_output->mref_start;

	//printf("%s\n", header_match);
	//printf("\tstrand = %i, chromosome = %i\n", sw_output->strand, sw_output->chromosome - 1);
	//printf("\tmquery_start = %i, mref_start = %i\n", mquery_start, mref_start);
	
	pos = cal->start + mref_start - 1;

	read_match = (char *) malloc(sizeof(char) * (mapped_len + 1));
	quality_match = (char *) malloc(sizeof(char) * (mapped_len + 1));
	c = 0;
	quality = fq_read->quality + mquery_start;
	for (int ii = 0; ii < mapped_len; ii++){
	  if (sw_output->mquery[ii] != '-') {
	    read_match[c] = sw_output->mquery[ii]; 
	    quality_match[c] = quality[c];
	    c++;
	  }
	}
	read_match[c] = '\0'; 
	quality_match[c] = '\0'; 
	
	// set optional fields
	optional_fields_length = 100;
	optional_fields = (char *) calloc(optional_fields_length, sizeof(char));
	
	p = optional_fields;
	AS = (int) score;
	
	sprintf(p, "ASi");
	p += 3;
	memcpy(p, &AS, sizeof(int));
	p += sizeof(int);
	
	sprintf(p, "NHi");
	p += 3;
	memcpy(p, &num_cals, sizeof(int));
	p += sizeof(int);
	
	sprintf(p, "NMi");
	p += 3;
	memcpy(p, &distance, sizeof(int));
	p += sizeof(int);

	// create an alignment and insert it into the list
	alignment = alignment_new();
	alignment_init_single_end(header_match, read_match, quality_match, 
				  cal->strand, 
				  cal->chromosome_id - 1, 
				  pos,
				  cigar, num_cigar_ops, norm_score * 254, 1, (num_cals > 1),
				  optional_fields_length, optional_fields, 0, alignment);
	
	array_list_insert(alignment, alignment_list);
	*/
      }
    }

    // free the cal list, and update the mapping list with the alignment list
    array_list_free(cal_list, (void *) cal_free);
    mapping_batch->mapping_lists[read_index] = alignment_list;
  }

/*
  double norm_score;
  // filter alignments by min_score
  for (size_t i = 0; i < sw_count; i++) {

    read_index = read_indices[i];
    fq_read = (fastq_read_t *) array_list_get(read_index, fq_batch);

    read_len = fq_read->length;
    //    norm_score = 0.0f;
    norm_score = NORM_SCORE(output->score_p[i], read_len, input->match);

    if (norm_score >= min_score) {
      // valid mappings, 
      //insert in the list for further processing
      mapping_list = mapping_batch->mapping_lists[read_index];
      array_list_set_flag(0, mapping_list);

      if (array_list_size(mapping_list) == 0) {
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
    }

    // free query and reference
    free(q[i]);
    free(r[i]);
  }
  mapping_batch->num_targets = new_num_targets;
*/

  // free
  sw_multi_output_free(output);

  // go to the next stage
  return POST_PAIR_STAGE;

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
//------------------------------------------------------------------------------------

