#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

//#include "bam_commons.h"
//#include "commons.h"
//#include "bam.h"

#include "breakpoint.h"
#include "containers/array_list.h"
#include "containers/cprops/trie.h"

#include "bioformats/bam/bam_file.h"
#include "bioformats/bam/alignment.h"

#include "bioformats/features/region/dna_map_region.h"
#include "bioformats/features/region/rna_map_region.h"

//--------------------------------------------------------------------

#define MAX_DATASET_LINE_LENGTH		4096

//--------------------------------------------------------------------

array_list_t *id_list;
cp_trie *trie;

//--------------------------------------------------------------------


typedef struct trie_node {
  uint32_t  mapped;
  uint32_t  right_mapped;
  uint32_t  wrong_mapped;
  uint32_t  wrong_mapped_chromosome;
  uint32_t  wrong_mapped_strand;
  uint32_t  not_mapped;
  uint32_t  right_not_mapped;
  uint32_t  wrong_not_mapped;

  uint32_t  right_sj;
  uint32_t  wrong_sj;

  size_t num_mappings;
  void *info;
} trie_node_t;

static inline trie_node_t *trie_node_new(void *info) {
  trie_node_t *p = (trie_node_t *) malloc(sizeof(trie_node_t));
  p->info = info;
  p->num_mappings = 0;

  p->mapped = 0;
  p->right_mapped = 0;
  p->wrong_mapped = 0;
  p->wrong_mapped_chromosome = 0;
  p->wrong_mapped_strand = 0;
  p->not_mapped = 0;
  p->right_not_mapped = 0;
  p->wrong_not_mapped = 0;

  p->right_sj = 0;
  p->wrong_sj = 0;

  return p;
}

static inline void trie_node_free(trie_node_t *p) {
  if (p) {
    if (p->info) {
      dna_map_region_t *region = p->info;
      if (region->chromosome) free(region->chromosome);
      free(region);
    }
    free(p);
  }
}

//--------------------------------------------------------------------
 
typedef struct trie_result {
  char* filename;
  uint32_t  mapped;
  uint32_t  right_mapped;
  uint32_t  wrong_mapped;
  uint32_t  wrong_mapped_chromosome;
  uint32_t  wrong_mapped_strand;
  uint32_t  not_mapped;
  uint32_t  right_not_mapped;
  uint32_t  wrong_not_mapped;

  uint32_t  rand_reads;
  uint32_t  multi_mapped;
  uint32_t  multi_right_mapped;
  uint32_t  multi_wrong_mapped;

  int margin;
  int right_sj;
  int wrong_sj;
} trie_result_t;


static inline void trie_result_init(trie_result_t *p) {
  if (p) {
    p->mapped = 0;
    p->right_mapped = 0;
    p->wrong_mapped = 0;
    p->wrong_mapped_chromosome = 0;
    p->wrong_mapped_strand = 0;
    p->not_mapped = 0;
    p->right_not_mapped = 0;
    p->wrong_not_mapped = 0;
    
    p->rand_reads = 0;
    p->multi_mapped = 0;
    p->multi_right_mapped = 0;
    p->multi_wrong_mapped = 0;
    p->right_sj = 0;
    p->wrong_sj = 0;
    //    memset(p, sizeof(trie_result_t), 0);
  }
}

static inline void trie_result_add(trie_result_t *src, trie_result_t *dst) {
  if (src && dst) {
    dst->mapped += src->mapped;
    dst->right_mapped += src->right_mapped;
    dst->wrong_mapped += src->wrong_mapped;
    dst->not_mapped += src->not_mapped;
    dst->right_not_mapped += src->right_not_mapped;
    dst->wrong_not_mapped += src->wrong_not_mapped;
    
    dst->rand_reads += src->rand_reads;
    dst->mapped += src->multi_mapped;
    dst->multi_right_mapped += src->multi_right_mapped;
    dst->multi_wrong_mapped += src->multi_wrong_mapped;
  }
}

static inline trie_result_t *trie_result_new() {
  trie_result_t *p = (trie_result_t *) calloc(1, sizeof(trie_node_t));
  trie_result_init(p);
  return p;
}

static inline void trie_result_free(trie_result_t *p) {
  if (p) free(p);
}

//--------------------------------------------------------------------

/**
 *	@brief Convert bam file into trie
 *	@param filename. Filename containing the bam
 *	@return cp_trie
 *
 **/
cp_trie* dna_bam_to_trie(char* filename);

//--------------------------------------------------------------------

/**
 *	@brief Convert simulated file (dwgsim) into trie
 *	@param filename. Filename containing the simulated dataset
 *	@return cp_trie
 *
**/
cp_trie* dna_dataset_to_trie(char* filename, trie_result_t* result);

//--------------------------------------------------------------------

/**
 * 	@brief Calculate the intersection between the trie and the bam_file
 * 	@param trie. The Trie calculated with dataset_to_trie or bam_to_trie
 * 	@param filename. The Bam filename to calculate the intersection
 * 	@result trie_result_t
 *
 *
 **/
void dna_intersection(cp_trie* trie, int margin, char* filename, trie_result_t* result);

//--------------------------------------------------------------------

void rna_intersection(cp_trie* trie, int margin, char* filename, 
		      trie_result_t* result, cp_hashtable *t);

//--------------------------------------------------------------------

void print_result(trie_result_t* result, int log);

//--------------------------------------------------------------------

cp_hashtable *load_transcriptome_validate(char *f);

cp_trie *rna_dataset_to_trie(char * file, trie_result_t* result);

typedef struct exon_coords {
  int start;
  int end;
} exon_coords_t;

typedef struct exon_data {
  int strand;
  int start;
  int end;
  char *transcript_id;
} exon_data_t;

exon_coords_t *new_exon_coords(int start, int end);
void exon_coords_free(exon_coords_t *p);
exon_data_t *exon_data_new(int strand, 
			   int start, int end,
			   char *transcript_id);
void exon_data_free(exon_data_t *p);

//--------------------------------------------------------------------

