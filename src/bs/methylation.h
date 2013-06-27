#ifndef METHYLATION_H
#define METHYLATION_H

//====================================================================================

#include <stdio.h>
#include <stdlib.h>
#include "aligners/bwt/bwt.h"
#include "aligners/sw/smith_waterman.h"
#include "bioformats/fastq/fastq_file.h"
#include "aligners/bwt/genome.h"
#include "buffers.h"

#include "containers/array_list.h"
#include "containers/list.h"

//====================================================================================

#define ACGT 0
#define AGT  1
#define ACT  2
#define AT   3

#define LIMIT_INF   0.20
#define LIMIT_SUP   0.30

//====================================================================================

/**
 * @brief  Transfor the sequence to the metilation case
 * @param  refs sequence of reference
 * @param  len length of the sequence
 * @param  type type of transformation
 * 
 * Replace the bases of the sequence depending the type.
 */
void replace(char *refs, int len, int type);

//====================================================================================

/**
 * @brief  Obtain the complementary base of the base.
 * @param  c base of a sequence.
 * @return complementary base.
 * 
 * Obtain the complementary base of the input.
 */
char complement (char c);

//====================================================================================

/**
 * @brief  Obtain the reverse complementary of a dna sequence
 * @param  orig sequence to transform
 * @param  dest sequence transformed
 * @param  len  length of the sequence
 * 
 * Transform the @a orig sequence of DNA, to the reverse complementary
 */
void rev_comp(char *query, char *seq, int len);

//====================================================================================

/**
 * @brief  Transform the sequence into the complementary
 * @param  seq sequence of nucleotides
 * @param  len length of the sequence
 * 
 * Complement the sequence element by element
 */
void comp(char *seq, int len);

//====================================================================================

/**
 * @brief  Make a copy of an array_list
 * @param  dest pointer to the new structure that is gona store the reads
 * @param  src  pointer to the structure that is gone be copied
 * 
 * Copy an array of reads into a new one
 */
void copy_array(array_list_t *dest, array_list_t *src);

//====================================================================================

/**
 * @brief  Make four copies of the array of reads for the bisulfite case
 * @param  src original data array
 * @param  dest1 
 * @param  dest2 
 * @param  dest3 
 * @param  dest4 
 * 
 * Copy the reads and insert in the four arrays (old version, no use)
 */
void cpy_array_bs(array_list_t *src, array_list_t *dest1, array_list_t *dest2, array_list_t *dest3, array_list_t *dest4);

//====================================================================================

/**
 * @brief  Make four copies of the array of reads, and transform for the bisulfite case
 * @param  src original data array
 * @param  dest_ct     
 * @param  dest_ct_rev 
 * @param  dest_ga     
 * @param  dest_ga_rev  
 * 
 * Copy the reads and transform it before insert in the new array
 */
void cpy_transform_array_bs(array_list_t *src, array_list_t *dest_ct, array_list_t *dest_ct_rev, array_list_t *dest_ga, array_list_t *dest_ga_rev);

//====================================================================================

/**
 * @brief  Combine the alignments from both arrays
 * @param  dest array is gona store the combination
 * @param  src array is going to unify the other
 * 
 * Insert the alignments from @a src into the mappings on @a dest
 */
void insert_mappings_array(array_list_t **dest, array_list_t **src);

//====================================================================================

/**
 * @brief  Include the new mappings existing in the "src" into "dest"
 * @param  dest Mappings that are gone store the combination
 * @param  src mapping are goin to be joined
 * 
 * Unify two mapping_list into one
 */
void insert_mappings(array_list_t *dest, array_list_t *src);

//====================================================================================

/**
 * @brief  Include the new regions existing in the "src" into "dest"
 * @param  dest Regions that are gone store the combination
 * @param  src regions are goin to be joined
 * 
 * Unify two region_list into one
 */
void insert_regions(array_list_t *dest, array_list_t *src);

//====================================================================================

/**
 * @brief  Set the sequence strand of the mappings on the array to 1
 * @param  src mappings to transform
 * 
 * Set all the mappins from all the reads in @a src on the reverse strand
 */
void transform_mappings_array(array_list_t **src);

//====================================================================================

/**
 * @brief  Set the sequence strand of the mappings to 1
 * @param  src Array of mappings
 * 
 * Set all the mappins in @a src on the reverse strand
 */
void transform_mappings(array_list_t *src);

//====================================================================================

/**
 * @brief  Set the region strand to 1
 * @param  src Array of regions
 * 
 * Set all the regions in @a src on the reverse strand
 */
void transform_regions(array_list_t *src);

//====================================================================================

/**
 * @brief  Select the index values where no read is mapped
 * @param  num_unmapped  store the number of unmapped reads
 * @param  unmapped_indices store the index of unmapped reads
 * @param  indices_1 vector with the count of unmapped reads
 * @param  num_reads number of reads analyzed
 * @param  lists,lists2 arrays with the information relative to the mappings
 * 
 * Updates the unmapped index for the next stage
 */
void select_targets_bs(size_t *num_unmapped, size_t *unmapped_indices,
		       size_t *indices_1, size_t num_reads,
		       array_list_t **lists, array_list_t **lists2);

//====================================================================================

/**
 * @brief  
 * @param  num_reads 
 * @param  unmapped_indices 
 * @param  unmapped 
 * @param  indices1 
 * 
 * 
 */
void update_targets_bs(size_t num_reads, size_t *unmapped_indices,
		       size_t unmapped, size_t *indices);

//====================================================================================

/**
 * @brief  Set the original read in the sequences of each mapping
 * @param  src1  Array with the mappings of the reads transformed
 * @param  src1  Array with the mappings of the reads transformed
 * @param  orig Fastq batch with the original reads
 * 
 * Go on the array of mappings, and set the sequence to the original
 */
void revert_mappings_seqs(array_list_t **src1, array_list_t **src2, array_list_t *orig);

//====================================================================================

/**
 * @brief  Obtain the reverse complementary reads stored in the array
 * @param  dest 
 * @param  src 
 * 
 * 
 */
void rev_comp_array(array_list_t *dest, array_list_t *src);

//====================================================================================

/**
 * @brief  Make an histogram of the read
 * @param  seq  
 * @param  len  
 * @param  type 
 * 
 * Return true if the number of certain base is greater than the filter
 */
int histogram_seq(char *seq, size_t len, int type);

//====================================================================================

#endif // METHYLATION_H
