#ifndef BATCH_WRITER_H
#define BATCH_WRITER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "cprops/hashtable.h"

#include "commons/commons.h"
#include "commons/system_utils.h"
#include "containers/list.h"

#include "bioformats/fastq/fastq_file.h"
#include "bioformats/fastq/fastq_batch.h"
#include "bioformats/bam-sam/bam_file.h"

#include "buffers.h"
#include "timing.h"


//====================================================================================
//  structures and prototypes
//====================================================================================

struct batch_writer_input {
  char* match_filename;
  char* mismatch_filename;

  char* splice_exact_filename;
  char* splice_extend_filename;
  
  //  char* header_filename;
  genome_t* genome;

  linked_list_t* list_p;

  // internal
  bam_file_t *bam_file;
  size_t total_batches;
  size_t total_reads;
  size_t total_mappings;
  size_t num_mapped_reads;
  size_t limit_print;
};

//------------------------------------------------------------------------------------

void batch_writer_input_init(char* match_filename, char* splice_exact_filename, 
			     char* splice_extend_filename, linked_list_t* list_p, 
			     genome_t* genome, batch_writer_input_t* input);

//====================================================================================

void batch_writer(batch_writer_input_t* input_p);
void batch_writer2(batch_writer_input_t* input_p);

bam_header_t *create_bam_header_by_genome(genome_t *genome);

//void batch_writer_splice(batch_writer_splice_input_t* input_p);
//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

#endif // BATCH_WRITER_H
