#ifndef METHYLATION_H
#define METHYLATION_H

#include <stdio.h>
#include <stdlib.h>
#include "aligners/bwt/bwt.h"
#include "aligners/sw/smith_waterman.h"
#include "bioformats/fastq/fastq_file.h"
#include "aligners/bwt/genome.h"
#include "buffers.h"

#include "containers/array_list.h"
#include "containers/list.h"


#define ACGT 0
#define AGT  1
#define ACT  2
#define AT   3


/**
 * @brief  Transfor the sequence to the metilation case
 * @param  refs sequence of reference
 * @param  len length of the sequence
 * @param  type type of transformation
 * 
 * Replace the bases of the sequence depending the type.
 */
void replace(char *refs, int len, int type);

/**
 * @brief  Obtain the complementary base of the base.
 * @param  c base of a sequence.
 * @return complementary base.
 * 
 * Obtain the complementary base of the input.
 */
char complement (char c);

/**
 * @brief  
 * @param  refs 
 * @param  len 
 * @param  type 
 * 
 * 
 */
void rev_comp(char *query, char *seq, int len);

/**
 * @brief  
 * @param  seq 
 * @param  len 
 * 
 * 
 */
void comp(char *seq, int len);

/**
 * @brief  
 * @param  src 
 * @param  dest 
 * 
 * 
 */
void cpy_array(array_list_t *src, array_list_t *dest);

/**
 * @brief  
 * @param  src 
 * @param  dest1 
 * @param  dest2 
 * @param  dest3 
 * @param  dest4 
 * 
 * 
 */
void cpy_array_bs(array_list_t *src, array_list_t *dest1, array_list_t *dest2, array_list_t *dest3, array_list_t *dest4);

/**
 * @brief  
 * @param  seq 
 * @param  len 
 * 
 * 
 */
void replace_array(array_list_t *reads, int type);

/**
 * @brief  
 * @param  dest 
 * @param  src 
 * 
 * 
 */
void rev_comp_array(array_list_t *dest, array_list_t *src);

/**
 * @brief  
 * @param  dest 
 * @param  src 
 * 
 * 
 */
void copy_array(array_list_t *dest, array_list_t *src);

/**
 * @brief  
 * @param  dest 
 * @param  src 
 * 
 * 
 */
void insert_mappings_array(array_list_t **dest, array_list_t **src);

/**
 * @brief  Include the new mappings existing in the "src" into "dest"
 * @param  dest Mappings that 
 * @param  src 
 * 
 * 
 */
void insert_mappings(array_list_t *dest, array_list_t *src);

/**
 * @brief  Set the sequence strand of the mappings on the array to 1
 * @param  src 
 * 
 * 
 */
void transform_mappings_array(array_list_t **src);

/**
 * @brief  Set the sequence strand of the mappings to 1
 * @param  src Array of mappings
 * 
 * 
 */
void transform_mappings(array_list_t *src);

/**
 * @brief  
 * @param  num_unmapped 
 * @param  unmapped_indices 
 * @param  unmapped1,unmapped2,unmapped3,unmapped4 
 * @param  indices1,indices2,indices3,indices4 
 * @param  lists,lists2
 * 
 * 
 */
void select_targets_bs(size_t *num_unmapped, size_t *unmapped_indices,
		       size_t *indices_1, size_t num_reads,
		       array_list_t **lists, array_list_t **lists2);

/**
 * @brief  
 * @param  num_unmapped 
 * @param  unmapped_indices 
 * @param  unmapped1,unmapped2,unmapped3,unmapped4 
 * @param  indices1,indices2,indices3,indices4 
 * @param  lists,lists2
 * 
 * 
 */
void update_targets_bs(size_t num_reads, size_t *unmapped_indices,
		       size_t unmapped, size_t *indices);

/**
 * @brief  Set the original read in the sequences of each mapping
 * @param  src1  Array with the mappings of the reads transformed
 * @param  src1  Array with the mappings of the reads transformed
 * @param  orig Fastq batch with the original reads
 * 
 * Go on the array of mappings, and set the sequence to the original
 */
void revert_mappings_seqs(array_list_t **src1, array_list_t **src2, array_list_t *orig);

#endif // METHYLATION_H
