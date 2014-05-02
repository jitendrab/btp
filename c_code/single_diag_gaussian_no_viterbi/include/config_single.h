#ifndef CONFIG_SINGLE
#define CONFIG_SINGLE

#define Max_Segs 50000 // max number of voiced speech segments in .scp file 
#define DIM 19
#define MAX_NUM_FEATURES 500000
#define MAX_NUM_STATES 100
#define MAX_NUM_MIX 10
#define INITIAL_NUM_MIX 1 //single gaussian

#define LAMBDA 1.1
#define INITIAL_STATES 20

#endif
