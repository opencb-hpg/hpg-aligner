
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "buffers.h"


//#define MAXLINE 2048

//===================================================================================
batch_t *batch_new(bwt_server_input_t *bwt_input,
                   region_seeker_input_t *region_input,
                   cal_seeker_input_t *cal_input,
                   pair_server_input_t *pair_input,
		   preprocess_rna_input_t *preprocess_rna,
                   sw_server_input_t *sw_input,
                   batch_writer_input_t *writer_input,
		   int mapping_mode,
                   mapping_batch_t *mapping_batch) {

  batch_t *b = (batch_t *) calloc(1, sizeof(batch_t));
  b->bwt_input = bwt_input;
  b->region_input = region_input;
  b->cal_input = cal_input;
  b->pair_input = pair_input;
  b->sw_input = sw_input;
  b->writer_input = writer_input;
  b->mapping_batch = mapping_batch;
  b->mapping_mode = mapping_mode;
  b->preprocess_rna = preprocess_rna;

  return b;
}

void batch_free(batch_t *b) {
  if (b) free(b);
}


//====================================================================================

void region_batch_init(array_list_t **allocate_mapping_p, fastq_batch_t *unmapped_batch_p, region_batch_t *region_batch_p){
  region_batch_p->allocate_mapping_p = allocate_mapping_p;
  region_batch_p->unmapped_batch_p = unmapped_batch_p;
}

void region_batch_free(region_batch_t *region_batch_p){
  for(int i = 0; i < region_batch_p->unmapped_batch_p->num_reads; i++){
    array_list_free(region_batch_p->allocate_mapping_p[i], (void *)region_bwt_free);
  }
  free(region_batch_p->allocate_mapping_p);
  fastq_batch_free(region_batch_p->unmapped_batch_p);
  free(region_batch_p);
  
}

//====================================================================================

sw_batch_t* sw_batch_new(unsigned int num_reads, array_list_t **allocate_cals_p, 
			 fastq_read_t **allocate_reads_p) {
  sw_batch_t* sw_batch_p = (sw_batch_t *)malloc(sizeof(sw_batch_t));

  sw_batch_p->num_reads = num_reads;
  sw_batch_p->allocate_reads_p = allocate_reads_p;
  sw_batch_p->allocate_cals_p = allocate_cals_p;
  
  return sw_batch_p;
}

void sw_batch_free(sw_batch_t *sw_batch_p) {
  
  for(int i = 0; i < sw_batch_p->num_reads; i++){
    array_list_free(sw_batch_p->allocate_cals_p[i], (void *)cal_free);
    fastq_read_free(sw_batch_p->allocate_reads_p[i]);
  }

  free(sw_batch_p->allocate_cals_p);
  free(sw_batch_p->allocate_reads_p);
  free(sw_batch_p);
}
void sw_batch_init(unsigned int num_reads, array_list_t **allocate_cals_p, 
		   fastq_read_t **allocate_reads_p, sw_batch_t *sw_batch_p) {
  sw_batch_p->num_reads = num_reads;
  sw_batch_p->allocate_reads_p = allocate_reads_p;
  sw_batch_p->allocate_cals_p = allocate_cals_p;
}

//====================================================================================
//  write_batch functions
//====================================================================================

write_batch_t* write_batch_new(unsigned int allocate_size, unsigned char flag) {
  write_batch_t* write_batch_p = (write_batch_t*) calloc(1, sizeof(write_batch_t));
  
  write_batch_p->flag = flag;
  write_batch_p->size = 0;
  
  if(flag != MATCH_FLAG){
    write_batch_p->allocated_size = allocate_size;
    write_batch_p->buffer_p = (void *) calloc(allocate_size, sizeof(char));
  }else{
    write_batch_p->allocated_size = allocate_size/sizeof(alignment_t *);
    write_batch_p->buffer_p = (void *) calloc(write_batch_p->allocated_size, sizeof(alignment_t *));
  }
  return write_batch_p;
}

//------------------------------------------------------------------------------------

void write_batch_free(write_batch_t* write_batch_p) {
 if (write_batch_p == NULL) return;
 
 if (write_batch_p->buffer_p != NULL) free(write_batch_p->buffer_p);
 
 free(write_batch_p);
}

//====================================================================================

report_optarg_t *report_optarg_new(int all, int n_best, int n_hits, int only_paired, int best) {
  report_optarg_t *p = (report_optarg_t*) calloc(1, sizeof(report_optarg_t));

  p->all = all;
  p->n_best = n_best;
  p->n_hits = n_hits;
  p->only_paired = only_paired;
  p->best = best;

  return p;
}

//------------------------------------------------------------------------------------

void report_optarg_free(report_optarg_t *p) {
  if (p != NULL)
    free(p);
}

//====================================================================================

pair_mng_t *pair_mng_new(int pair_mode, size_t min_distance, 
			 size_t max_distance, int report_only_paired) {
  pair_mng_t *p = (pair_mng_t*) calloc(1, sizeof(pair_mng_t));

  p->pair_mode = pair_mode;
  p->min_distance = min_distance;
  p->max_distance = max_distance;
  p->report_only_paired = report_only_paired;

  return p;
}

//------------------------------------------------------------------------------------

void pair_mng_free(pair_mng_t *p) {
  if (p != NULL)
    free(p);
}

//====================================================================================

cal_batch_t* cal_batch_new(array_list_t **allocate_mapping, fastq_batch_t *unmapped_batch){
  cal_batch_t* cal_batch = (cal_batch_t *)malloc(sizeof(cal_batch_t));
  
  cal_batch->allocate_mapping = allocate_mapping;
  cal_batch->unmapped_batch = unmapped_batch;
  
  return cal_batch;
}

void cal_batch_free(cal_batch_t *cal_batch){
  for(int i = 0; i < cal_batch->unmapped_batch->num_reads; i++){
    array_list_free(cal_batch->allocate_mapping[i], (void *)region_bwt_free);
  }
  free(cal_batch->allocate_mapping);
  fastq_batch_free(cal_batch->unmapped_batch);
  free(cal_batch);
  
}

//====================================================================================

unsigned int pack_junction(unsigned int chromosome, unsigned int strand, size_t start, size_t end, size_t junction_id, size_t num_reads, char* buffer_p){
  int len;
  char str[1024];
  char *chr_p, *p = buffer_p;
  char strand_char[2] = {'+', '-'};

  if (chromosome == 23) { sprintf(str, "%c\0", 'X'); }
  else if (chromosome == 24) { sprintf(str, "%c\0", 'Y'); }
  else if (chromosome == 25) { sprintf(str, "%s\0", "MT"); }
  else { sprintf(str, "%i\0", chromosome); }
 
  len = strlen(str);
  memcpy(p, str, len);
  p += len;
  *p = '\t';
  p++;
  
  sprintf(str, "%lu", start);
  len = strlen(str);
  memcpy(p, str, len); 
  p += len;
  *p = '\t'; 
  p++;
  
  sprintf(str, "%lu", end);
  len = strlen(str);
  memcpy(p, str, len); 
  p += len;
  *p = '\t'; 
  p++;
  
  sprintf(str, "JUNCTION_%lu", junction_id);
  len = strlen(str);
  memcpy(p, str, len); 
  p += len;
  *p = '\t'; 
  p++;
  
  sprintf(str, "%lu", num_reads);
  len = strlen(str);
  memcpy(p, str, len); 
  p += len;
  *p = '\t'; 
  p++;
  
  *p = strand_char[strand]; 
  p++;
  *p = '\n'; 
  p++;
  

  return (p - buffer_p);
}

//=====================================================================================
//=====================================================================================

mapping_batch_t *mapping_batch_new(array_list_t *fq_batch, pair_mng_t *pair_mng) {

  mapping_batch_t *p = (mapping_batch_t *) calloc(1, sizeof(mapping_batch_t));
  size_t num_reads = array_list_size(fq_batch);

  p->action = BWT_ACTION;
  p->num_targets = 0;
  p->num_extra_targets = 0;
  p->num_allocated_targets = num_reads;
  p->extra_stage_do = 0;

  if (!pair_mng) { 
    p->pair_mng = pair_mng_new(SINGLE_END_MODE, 0, 0, 0); 
  } else {
    p->pair_mng = pair_mng_new(pair_mng->pair_mode, pair_mng->min_distance, 
			       pair_mng->max_distance, pair_mng->report_only_paired); 
  }

  p->num_to_do = 0;
  p->fq_batch = fq_batch;
  p->targets = (size_t *) calloc(num_reads, sizeof(size_t));
  p->extra_targets = (size_t *) calloc(num_reads, sizeof(size_t));
  p->extra_stage_id = (unsigned char *) calloc(num_reads, sizeof(unsigned char));
  p->mapping_lists = (array_list_t **) calloc(num_reads, sizeof(array_list_t*));

  //for debug. TODO:delete
  p->bwt_mappings = (unsigned char *)calloc(num_reads, sizeof(unsigned char));

  for (size_t i = 0; i < num_reads; i++) {
    p->mapping_lists[i] = array_list_new(500, 
					 1.25f, 
					 COLLECTION_MODE_ASYNCHRONIZED); 
  }
    
  // added by PP for bisulfite
  p->num_targets2 = 0;
  p->num_to_do2 = 0;
  p->targets2 = (size_t *) calloc(num_reads, sizeof(size_t));

  return p;
}

//------------------------------------------------------------------------------------

mapping_batch_t *mapping_batch_new_by_num(size_t num_reads, pair_mng_t *pair_mng) {

  mapping_batch_t *p = (mapping_batch_t *) calloc(1, sizeof(mapping_batch_t));

  p->action = BWT_ACTION;
  p->num_targets = 0;
  p->num_extra_targets = 0;
  p->num_allocated_targets = num_reads;
  p->extra_stage_do = 0;

  if (!pair_mng) { 
    p->pair_mng = pair_mng_new(SINGLE_END_MODE, 0, 0, 0); 
  } else {
    p->pair_mng = pair_mng_new(pair_mng->pair_mode, pair_mng->min_distance, 
			       pair_mng->max_distance, pair_mng->report_only_paired); 
  }

  p->num_to_do = 0;
  p->targets = (size_t *) calloc(num_reads, sizeof(size_t));
  p->extra_targets = (size_t *) calloc(num_reads, sizeof(size_t));
  p->extra_stage_id = (unsigned char *) calloc(num_reads, sizeof(unsigned char));
  p->mapping_lists = (array_list_t **) calloc(num_reads, sizeof(array_list_t*));
  p->old_mapping_lists = (array_list_t **) calloc(num_reads, sizeof(array_list_t*));

  for (size_t i = 0; i < num_reads; i++) {
    p->mapping_lists[i] = array_list_new(10, 
					 1.25f, 
					 COLLECTION_MODE_ASYNCHRONIZED);
  }
    
  return p;
}


//------------------------------------------------------------------------------------

void mapping_batch_free(mapping_batch_t *p) {
  if (p == NULL) return;
  
  if (p->fq_batch) { array_list_free(p->fq_batch, (void *) fastq_read_free); }
  if (p->targets) { free(p->targets); }
  if (p->mapping_lists) { free(p->mapping_lists); }
  if (p->pair_mng) { free(p->pair_mng); }
  if (p->extra_stage_id) { free(p->extra_stage_id); }
  if (p->extra_targets) { free(p->extra_targets); }
  if (p->old_mapping_lists) {
    free(p->old_mapping_lists);
  }
  if (p->bwt_mappings) { free(p->bwt_mappings); }
  
  // added by PP
  if (p->CT_fq_batch) { array_list_free(p->CT_fq_batch, (void *) fastq_read_free); }
  if (p->CT_rev_fq_batch) { array_list_free(p->CT_rev_fq_batch, (void *) fastq_read_free); }
  if (p->GA_fq_batch) { array_list_free(p->GA_fq_batch, (void *) fastq_read_free); }
  if (p->GA_rev_fq_batch) { array_list_free(p->GA_rev_fq_batch, (void *) fastq_read_free); }
  if (p->mapping_lists2) { free(p->mapping_lists2); }
  if (p->targets2) { free(p->targets2); }
  if (p->bs_status) {free(p->bs_status); }
  
  free(p);
}

//------------------------------------------------------------------------------------
/*
rna_batch_t *rna_batch_new(array_list_t *fq_batch) {
  rna_batch_t *p = (rna_batch_t *) calloc(1, sizeof(rna_batch_t));

  size_t num_reads = array_list_size(fq_batch);

  p->action = BWT_ACTION;
  p->all_targets = 1;
  p->num_targets = 0;
  p->num_allocated_targets = num_reads;
  p->num_mapping_lists = num_reads;

  p->num_done = 0;
  p->num_to_do = 0;

  p->fq_batch = fq_batch;
  p->targets = (size_t *) calloc(num_reads, sizeof(size_t));
  p->mapping_lists = (array_list_t **) calloc(num_reads, sizeof(array_list_t*));
  for (size_t i = 0; i < num_reads; i++) {
    p->mapping_lists[i] = array_list_new(500, 
					 1.25f, 
					 COLLECTION_MODE_ASYNCHRONIZED); 
  }
    
  return p;
}

//------------------------------------------------------------------------------------

void rna_batch_free(rna_batch_t *p) {
  if (p == NULL) return;
  
  if (p->fq_batch != NULL) array_list_free(p->fq_batch, (void *)fastq_read_free);
  if (p->targets != NULL) free(p->targets);
  if (p->mapping_lists != NULL) { free(p->mapping_lists); }
  
  free(p);
}
*/

//------------------------------------------------------------------------------------

void bs_context_init(bs_context_t *bs_context, size_t num_reads) {
  bs_context->context_CpG = array_list_new(num_reads * 4, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  bs_context->context_CHG = array_list_new(num_reads * 4, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  bs_context->context_CHH = array_list_new(num_reads * 4, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  bs_context->context_MUT = array_list_new(num_reads * 4, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  
  bs_context->CpG_methyl   = 0;
  bs_context->CpG_unmethyl = 0;
  bs_context->CHG_methyl   = 0;
  bs_context->CHG_unmethyl = 0;
  bs_context->CHH_methyl   = 0;
  bs_context->CHH_unmethyl = 0;
  bs_context->MUT_methyl   = 0;
  bs_context->num_bases    = 0;
}

//------------------------------------------------------------------------------------
