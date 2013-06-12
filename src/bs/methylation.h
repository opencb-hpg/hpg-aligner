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
void insert_mappings(array_list_t **dest, array_list_t **src);

/**
 * @brief  
 * @param  src 
 * 
 * 
 */
void transform_mappings(array_list_t **src);

#endif // METHYLATION_H
