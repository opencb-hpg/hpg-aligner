#include "bam_trie.h"


#define DEFAULT_MAX_DISTANCE_SIZE		500


#define BAM_TRIE_USAGE_HELP "Usage: bam_trie_test --ref-align <fastq_file> --bam <bam_file> --mode <dna|rna> --transcriptome <file> --pair-mode <0|1>\n"

/* **********************************************
 *
 *    		Global variables  		*
 * **********************************************/

int main(int argc, char **argv) {

  int c;
  int pos = 0;
  int option_index = 0;
  
  static struct option long_options[]={
    {"ref-align",		required_argument, 0, 'a'},
    {"transcriptome",		required_argument, 0, 't'},
    {"ref-bam", 		required_argument, 0, 'b'},
    {"bam",			required_argument, 0, 'c'},
    {"margin-length",		required_argument, 0, 'm'},
    {"log-file",		required_argument, 0, 'd'},
    {"mode",		        required_argument, 0, 'e'},
    {"pair-mode",		required_argument, 0, 'p'},
    {0,0,0,0}
  };
  
  int argc_with_file_options = 0;
  char ** argv_with_file_options = NULL;
  //		char ** argv_from_file_options = NULL;
  char *token[100];
  char *transcriptome_file;
  char *ref_file;
  char *bam_file;
  int log = 1;
  int margin = 0;
  trie_result_t *result;
  int align_bam = 0; // align_bam = 0 => Align, align_bam=1 => BAM
  int mode;
  int pair_mode = 0;

  if(argc < 5) {
    printf(BAM_TRIE_USAGE_HELP);	
    exit(0);
  }

  argc_with_file_options = argc;
  argv_with_file_options = argv;
  
  while((c = getopt_long(argc_with_file_options, argv_with_file_options, "", long_options, &option_index)) != -1) {
    switch ( c ) {
    case 'a': // ref-align	
      ref_file = (char*) calloc(strlen(optarg) + 1, sizeof(char));
      strcpy(ref_file,optarg);
      align_bam = 0;
      break;
      
    case 'b': // ref-bam	
      ref_file = (char*) calloc(strlen(optarg) + 1, sizeof(char));
      strcpy(ref_file, optarg);
      align_bam = 1;
      break;
      
    case 'c': // bam	
      bam_file = (char*) calloc(strlen(optarg) + 1, sizeof(char));
      strcpy(bam_file, optarg);
      
      token[pos] = strtok(bam_file, ",");
      while (token[pos]!= NULL) {
	pos++;
	token[pos] = strtok(NULL,",");
      }
      break;
      
    case 'm':  
      margin = atoi(optarg);
      break;

    case 'e':
      if (strcmp(optarg, "dna") == 0) {
	mode = 0;
      } else if (strcmp(optarg, "rna") == 0) {
	mode = 1;
      } else {
	printf("Unknown mode %s\n", optarg);
	printf(BAM_TRIE_USAGE_HELP);
	exit(-1);
      }
      break;
    case 't':
      transcriptome_file = (char*) calloc(strlen(optarg) + 1, sizeof(char));
      strcpy(transcriptome_file, optarg);
      break;
    case 'p':
      pair_mode = atoi(optarg);
      break;
    default:	
      break;
    }				/* -----  end switch  ----- */

  }
  
  result = (trie_result_t*) calloc(1, sizeof(trie_result_t));
  id_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  if (mode == 1) {
    printf("Loading transcriptome file '%s' ...\n", transcriptome_file);
    cp_hashtable *t = load_transcriptome_validate(transcriptome_file);
    printf("Load transcriptome file done!\n", transcriptome_file);
    
    printf("Loading FASTQ file '%s' ...\n", ref_file);
    if (!pair_mode) {
      trie = rna_dataset_to_trie(ref_file, result);
    } else {
      trie = rna_dataset_to_trie_pair(ref_file, result);
    }
    printf("Load FASTQ file done!\n", transcriptome_file);

    //for(int i = 0; i < pos; i++) {
    printf("Validating BAM file '%s' ...\n\n", token[0]);
    if (!pair_mode) {
      rna_intersection(trie, margin, token[0], result, t);
      print_result(result, log);
    } else {
      rna_intersection_pair(trie, margin, token[0], result, t);
      print_result_pair(result, log);
    }

    printf("\nValidate BAM file done!\nDONE!\n");
  } else {
    printf("Loading FASTQ file '%s' ...\n", ref_file);
    trie = dna_dataset_to_trie(ref_file, result);
    printf("Load FASTQ file done!\n", transcriptome_file);

    for(int i = 0; i < pos; i++) {
      dna_intersection(trie, margin, token[i], result);
      print_result(result,log);
    }

  }

  /*
  array_list_t *list = cp_hashtable_get(t, "ENST00000539153");
  printf("SEARCH RESULTS 'ENST00000539153': \n");
  for (int i = 0; i < list->size; i++) {
    exon_coords_t *coords = array_list_get(i, list);
    printf("\t[%lu-%lu]\n", coords->start, coords->end);
  }
  */  


  // free memory
  free(result);
  free(ref_file);
  free(bam_file);
  cp_trie_destroy(trie);
  array_list_free(id_list, (void *) free);

  return 0;
}

