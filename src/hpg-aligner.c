#include "dna/dna_aligner.h"
#include "rna/rna_aligner.h"
#include "bs/bs_aligner.h"
#include "build-index/index_builder.h"


double emboss_matrix_t = 0.0f, emboss_tracking_t = 0.0f;
double sse_matrix_t = 0.0f, sse_tracking_t = 0.0f;
double sse1_matrix_t = 0.0f, sse1_tracking_t = 0.0f;
double avx_matrix_t = 0.0f, avx_tracking_t = 0.0f;
double avx1_matrix_t = 0.0f, avx1_tracking_t = 0.0f;

//--------------------------------------------------------------------
// constants
//--------------------------------------------------------------------

#define OPTIONS 			30
#define MIN_ARGC  			5
#define NUM_SECTIONS_STATISTICS 	5
#define NUM_SECTIONS_STATISTICS_SB	21
//#define REQUIRED 1
//#define NO_REQUIRED 0

//--------------------------------------------------------------------
// global variables for log functions
//--------------------------------------------------------------------

//int log_level = DEFAULT_LOG_LEVEL;
//bool verbose = true;

//--------------------------------------------------------------------
// global variables for timing and capacity meausures
//--------------------------------------------------------------------

char time_on = 0;
char statistics_on = 0;
timing_t* timing = NULL;
basic_statistics_t *basic_st;
cal_st_t cal_st;
double kl_time;

pthread_cond_t cond_sp;
pthread_mutex_t mutex_sp;

size_t bwt_correct = 0;
size_t bwt_error = 0;
pthread_mutex_t bwt_mutex;

// timing
double main_time;


void run_dna_aligner(genome_t *genome, bwt_index_t *bwt_index, 
		     bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg, 
		     pair_mng_t *pair_mng, report_optarg_t *report_optarg,
		     options_t *options);

void run_rna_aligner(genome_t *genome, bwt_index_t *bwt_index, pair_mng_t *pair_mng,
		     bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg,
		     report_optarg_t *report_optarg, options_t *options);


//void run_bs_aligner(genome_t *genome, genome_t *genome1, genome_t *genome2, 
void run_bs_aligner(genome_t *genome2, genome_t *genome1,
		    bwt_index_t *bwt_index2, bwt_index_t *bwt_index1,
		    bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg,
		    pair_mng_t *pair_mng, report_optarg_t *report_optarg,
		    options_t *options);


//--------------------------------------------------------------------
// main parameters support
//--------------------------------------------------------------------
int main(int argc, char* argv[]) {
  pthread_mutex_init(&cal_st.mutex, NULL);
  pthread_mutex_init(&mutex_sp, NULL);
  pthread_mutex_init(&bwt_mutex, NULL);

  const char HEADER_FILE[1024] = "Human_NCBI37.hbam\0";
  basic_st = basic_statistics_new();

  // init logs, after parsing the command-line
  // logs will be re-set according to the command-line
  log_level = LOG_FATAL_LEVEL;
  log_verbose = 1;
  log_file = NULL;

  if (argc <= 1) {
    LOG_FATAL("Missing command.\nValid commands are:\n\tdna: to map DNA sequences\n\trna: to map RNA sequences\n\tbs: to map BS sequences\n\tbuild-index: to create the genome index.\nUse -h or --help to display hpg-aligner options.");
  }

  if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
    usage_cli();
  }

  char *command = argv[1];  
  // We need to consume command: {dna | rna | bs | build-index}
  argc -= 1;
  argv += 1;

  if(strcmp(command, "dna") != 0 && 
     strcmp(command, "rna") != 0 &&
     strcmp(command, "bs" ) != 0 && 
     strcmp(command, "build-index") != 0) {
    LOG_FATAL("Command unknown.\nValid commands are:\n\tdna: to map DNA sequences\n\trna: to map RNA sequences\n\tbs: to map BS sequences\n\tbuild-index: to create the genome index.\nUse -h or --help to display hpg-aligner options.");
  }

  // parsing options
  options_t *options = parse_options(argc, argv);

  // now, we can set logs according to the command-line
  init_log_custom(options->log_level, 1, "hpg-aligner.log", "w");

  validate_options(options, command);
  LOG_DEBUG_F("Command Mode: %s\n", command);

  if (!strcmp(command, "build-index")) {
    if (options->bs_index == 0) { 
      //printf("Regular index generation\n");
      run_index_builder(options->genome_filename, options->bwt_dirname, options->index_ratio);
      LOG_DEBUG("Done !!\n");
      exit(0);
    }
    else { // bisulphite index generation
      printf("\nBisulphite index generation\n");
      
      /** **************************************************************************	*
       * 										*
       * Generates the genome transform from the input and builds the index		*
       * 										*
       * The genome transformed are stored in the directory give by the user,		*
       * and the index are stored in subfolders				       		*
       * 										*
       * ***************************************************************************	*/

      // generate binary code for original genome
      char binary_filename[strlen(options->bwt_dirname) + 128];
      sprintf(binary_filename, "%s/dna_compression.bin", options->bwt_dirname);
      generate_codes(binary_filename, options->genome_filename);

      char bs_dir1[256];
      sprintf(bs_dir1, "%s/AGT_index", options->bwt_dirname);
      //if (is_directory(bs_dir1) == 0) {
      create_directory(bs_dir1);
      //}

      LOG_DEBUG("Generation of AGT index\n");
      char genome_1[256];
      sprintf(genome_1, "%s/AGT_genome.fa", options->bwt_dirname);
      char gen1[256];
      sprintf(gen1, "sed 's/C/T/g' %s > %s",options->genome_filename, genome_1);
      system(gen1);

      run_index_builder_bs(genome_1, bs_dir1, options->index_ratio, "AGT");
      LOG_DEBUG("AGT index Done !!\n");

      LOG_DEBUG("Generation of ACT index\n");
      char bs_dir2[256];
      sprintf(bs_dir2, "%s/ACT_index", options->bwt_dirname);
      //if (is_directory(bs_dir2) == 0) {
      create_directory(bs_dir2);
      //}

      char genome_2[256];
      sprintf(genome_2, "%s/ACT_genome.fa", options->bwt_dirname);
      char gen2[256];
      sprintf(gen2, "sed 's/G/A/g' %s > %s",options->genome_filename, genome_2);
      system(gen2);

      run_index_builder_bs(genome_2, bs_dir2, options->index_ratio, "ACT");
      LOG_DEBUG("ACT index Done !!\n");
      exit(0);
    }
  }

  time_on =  (unsigned int) options->timming;
  statistics_on =  (unsigned int) options->statistics;

  genome_t *genome, *genome1, *genome2;
  bwt_index_t *bwt_index, *bwt_index1, *bwt_index2;

  // load index for dna/rna or for bisulfite case

  if (strcmp(command, "bs" ) != 0)
  {
    // genome parameters
    LOG_DEBUG("Reading genome...");
    //genome_t* genome = genome_new("dna_compression.bin", options->bwt_dirname);
    genome = genome_new("dna_compression.bin", options->bwt_dirname);
    LOG_DEBUG("Done !!");
    
    // BWT index
    LOG_DEBUG("Reading bwt index...");
    //if (time_on) { timing_start(INIT_BWT_INDEX, 0, timing_p); }
    //bwt_index_t *bwt_index = bwt_index_new(options->bwt_dirname);
    
    bwt_index = bwt_index_new(options->bwt_dirname);
    LOG_DEBUG("Reading bwt index done !!");
  }
  else {
    // change the filter_*_mappings values for bisulfite alignments
    options->filter_read_mappings = DEFAULT_FILTER_READ_MAPPINGS_BS;
    options->filter_seed_mappings = DEFAULT_FILTER_SEED_MAPPINGS_BS;

    char bs_dir1[256];
    sprintf(bs_dir1, "%s/AGT_index", options->bwt_dirname);
    char bs_dir2[256];
    sprintf(bs_dir2, "%s/ACT_index", options->bwt_dirname);

    // genome parameters
    LOG_DEBUG("Reading genomes...");
    //genome  = genome_new("dna_compression.bin", options->bwt_dirname);
    genome1 = genome_new("dna_compression.bin", bs_dir1);
    genome2 = genome_new("dna_compression.bin", bs_dir2);
    LOG_DEBUG("Done !!");
    
    // BWT index
    //if (time_on) { timing_start(INIT_BWT_INDEX, 0, timing_p); }

    LOG_DEBUG("Loading AGT index...");
    bwt_index1 = bwt_index_new(bs_dir1);
    /*
    printf("+++bwt_index1 -> %s\n" ,bwt_index1->nucleotides);
    bwt_index1->nucleotides = strdup(readNucleotide(bs_dir1, "Nucleotide"));
    bwt_init_replace_table(bwt_index1->nucleotides, bwt_index1->table, bwt_index1->rev_table);
    printf("***bwt_index1 -> %s\n" ,bwt_index1->nucleotides);
    */
    LOG_DEBUG("Loading AGT index done !!");

    LOG_DEBUG("Loading ACT index...");
    bwt_index2 = bwt_index_new(bs_dir2);
    /*
    printf("+++bwt_index2 -> %s\n", bwt_index2->nucleotides);
    bwt_index2->nucleotides = strdup(readNucleotide(bs_dir2, "Nucleotide"));
    bwt_init_replace_table(bwt_index2->nucleotides, bwt_index2->table, bwt_index2->rev_table);
    printf("***bwt_index2 -> %s\n", bwt_index2->nucleotides);
    */
    LOG_DEBUG("Loading ACT index done !!");

    //exit(-1);
  }

  
  //BWT parameters
  bwt_optarg_t *bwt_optarg = bwt_optarg_new(1, 0,
					    options->filter_read_mappings, 
					    options->filter_seed_mappings);
  
  // CAL parameters
  cal_optarg_t *cal_optarg = cal_optarg_new(options->min_cal_size, options->seeds_max_distance, 
					    options->num_seeds, options->min_num_seeds_in_cal,
					    options->seed_size, options->min_seed_size, 
					    options->cal_seeker_errors);
  
  // paired mode parameters
  pair_mng_t *pair_mng = pair_mng_new(options->pair_mode, options->pair_min_distance, 
				      options->pair_max_distance, options->report_only_paired);
  
  // report parameters
  report_optarg_t *report_optarg = report_optarg_new(options->report_all,
						     options->report_n_best,
						     options->report_n_hits, 
						     options->report_only_paired,
						     options->report_best);
  /*
  //*********************************
  {
    char *seq = strdup("CCTAACCAACATAATAAAACCCCATCTCTACTAAAAATACAAAAAAATTAACAAAACATAATAACAAATACCTATAATCCCAACTACTCAAAAAACTAAA");
    alignment_t *alig;
    array_list_t *mapping_list = array_list_new(100000, 1.25f, 
						COLLECTION_MODE_SYNCHRONIZED);
    
    size_t num_mappings;
    
    num_mappings = bwt_map_forward_inexact_seq(seq, bwt_optarg, 
					       bwt_index2, mapping_list);
    printf("aux seq: %s\n", seq);
    printf("num_mappings = %lu\n", num_mappings);
    for (size_t i = 0; i < num_mappings; i++) {
      alig = array_list_get(i, mapping_list);
      printf("%lu\t---------------------\n", i);
      printf("\tstrand = %i, chromosome = %i, position = %i\n", 
	     alig->seq_strand, alig->chromosome, alig->position);
    }
  }
  //*********************************
  */

  if (!strcmp(command, "rna")) {
    LOG_DEBUG("init table...");
    initTable();
    LOG_DEBUG("init table done !!");
    // RNA version
    run_rna_aligner(genome, bwt_index, pair_mng, bwt_optarg, cal_optarg, report_optarg, options);
  } else if (!strcmp(command, "dna")) {
    LOG_DEBUG("init table...");
    initTable();
    LOG_DEBUG("init table done !!");
    // DNA version
    run_dna_aligner(genome, bwt_index, bwt_optarg, cal_optarg, pair_mng, report_optarg, options);
  } else { // if (!strcmp(command, "bs")) {
    // BS version
    run_bs_aligner(genome1, genome2, bwt_index1, bwt_index2,
    //run_bs_aligner(genome, genome1, genome2, bwt_index1, bwt_index2,
		   bwt_optarg, cal_optarg, pair_mng, report_optarg, options);
  }

  LOG_DEBUG("main done !!");

  // Free memory
  
  if (strcmp(command, "bs" ) != 0)
  {
    bwt_index_free(bwt_index);
    genome_free(genome);
  } else {
    //genome_free(genome);
    bwt_index_free(bwt_index1);
    genome_free(genome1);
    bwt_index_free(bwt_index2);
    genome_free(genome2);
  }

  bwt_optarg_free(bwt_optarg);
  cal_optarg_free(cal_optarg);
  pair_mng_free(pair_mng);
  report_optarg_free(report_optarg);

  if (time_on) { 
    //timing_display(timing);
  }

  basic_statistics_display(basic_st, !strcmp(command, "rna"), main_time / 1000000);

  if (time_on){ timing_free(timing); }

  options_free(options);

  //printf("Reads mapped in BWT Section %lu (%lu correct(%f) and %lu misses(%f))\n", bwt_correct + bwt_error, bwt_correct, (float)(bwt_correct * 100)/(float)(bwt_correct + bwt_error), bwt_error, (float)(bwt_error * 100)/(float)(bwt_correct + bwt_error));
  
  return 0;
}

//--------------------------------------------------------------------
