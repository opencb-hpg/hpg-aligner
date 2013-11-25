gcc -g -std=c99 -D_XOPEN_SOURCE=600 -D_GNU_SOURCE -fopenmp -o bam_trie_test bam_trie.c bam_trie_test.c ../src/breakpoint.c -I../src/ -I ../lib/common-libs/ -I ../lib/bioinfo-libs -L ../lib/bioinfo-libs/ -L ../lib/common-libs/ -lbioinfo -lcommon -lm -lz 

gcc -g -std=c99 -D_XOPEN_SOURCE=600 -D_GNU_SOURCE -fopenmp -o bam_trie_real bam_trie.c  bam_trie_real.c ../src/breakpoint.c -I../src/ -I ../lib/common-libs/ -I ../lib/bioinfo-libs -L ../lib/bioinfo-libs/ -L ../lib/common-libs/ -lbioinfo -lcommon -lm -lz 
