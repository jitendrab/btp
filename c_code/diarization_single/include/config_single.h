#ifndef CONFIG_SINGLE
#define CONFIG_SINGLE

#define Max_Segs 5000 // max number of voiced speech segments in .scp file 
#define DIM 19
#define MAX_NUM_FEATURES 200000
#define MAX_NUM_STATES 16
#define MAX_NUM_MIX 100
#define INITIAL_NUM_MIX 1 //single gaussian
#define LAMBDA 8
#define INITIAL_STATES 10
#define MIN_DUR 250
#endif
