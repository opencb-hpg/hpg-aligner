#ifndef RNA_ALIGNER_H
#define RNA_ALIGNER_H


#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>

#include "commons/log.h"
#include "commons/file_utils.h"

#include "containers/workflow_scheduler.h"

#include "bioformats/fastq/fastq_batch_reader.h"

#include "error.h"
#include "timing.h"
#include "buffers.h"
#include "bwt_server.h"
#include "batch_writer.h"
#include "region_seeker.h"
#include "cal_seeker.h"
#include "sw_server.h"
#include "pair_server.h"
#include "preprocess_rna.h"
#include "options.h"
#include "statistics.h"
#include "workflow_functions.h"


typedef struct extra_stage {
  workflow_t *workflow;
  linked_list_t *align_list;
  pair_mng_t *pair_mng;
  char *intron_filename;
} extra_stage_t;

typedef struct buffer_item {
  fastq_read_t *read;
  array_list_t *alignments_list;
}buffer_item_t;

typedef struct buffer_pair_item {
  fastq_read_t *read_1;
  array_list_t *alignments_list_1;
  fastq_read_t *read_2;
  array_list_t *alignments_list_2;
}buffer_pair_item_t;

void run_rna_aligner(genome_t *genome, bwt_index_t *bwt_index, 
		     pair_mng_t *pair_mng,
		     bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg, 
		     options_t *options);

#endif
