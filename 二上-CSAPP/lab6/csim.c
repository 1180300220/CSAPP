#include <getopt.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <errno.h>
#include "cachelab.h"

//#define DEBUG_ON 
#define ADDRESS_LENGTH 64

/* Type: Memory address */
typedef unsigned long long int mem_addr_t;

/* Type: Cache line
   LRU is a counter used to implement LRU replacement policy  */
typedef struct cache_line {
    char valid;
    mem_addr_t tag;
    unsigned long long int lru;
} cache_line_t;

typedef cache_line_t* cache_set_t;
typedef cache_set_t* cache_t;

/* Globals set by command line args */
int verbosity = 0; /* print trace if set */
int s = 0; /* set index bits */
int b = 0; /* block offset bits */
int E = 0; /* associativity */
char* trace_file = NULL;

/* Derived from command line args */
int S; /* number of sets */
int B; /* block size (bytes) */

/* Counters used to record cache statistics */
int miss_count = 0;
int hit_count = 0;
int eviction_count = 0;
unsigned long long int lru_counter = 1;

/* The cache we are simulating */
cache_t cache;  
mem_addr_t set_index_mask;


/* 
 * initCache - Allocate memory, write 0's for valid and tag and LRU
 * also computes the set_index_mask
 */
void initCache()
{
 //todo...
    //为组申请空间
    cache = (cache_set_t *)malloc(sizeof(cache_set_t) * S);
    for(int i = 0;i < S;i++){
        //为组中每行申请空间
        cache[i] = (cache_line_t*)malloc(sizeof(cache_line_t)*E);
        for(int j = 0;j < E;j++){
            //初始化每行中的有效位，标记值，lru值。
            cache[i][j].valid = 0;
            cache[i][j].tag = 0;
            cache[i][j].lru = 0;
        }
    }
    //计算组索引
    set_index_mask = (mem_addr_t)(S-1);
}


/* 
 * freeCache - free allocated memory
 */
void freeCache()
{
 //todo...
    for(int i = 0;i < S;i++)    free(cache[i]);
    free(cache);
}


/* 
 * accessData - Access data at memory address addr.
 *   If it is already in cache, increast hit_count
 *   If it is not in cache, bring it in cache, increase miss count.
 *   Also increase eviction_count if a line is evicted.
 */
void accessData(mem_addr_t addr)
{
 //todo...
    //计算组偏移量，标记值
    mem_addr_t set_index = (addr>>b) & set_index_mask;
    mem_addr_t tag = addr>>(s+b);
    cache_set_t cache_set = cache[set_index];
    
    //命中情况
    for(int i = 0;i < E;i++){
        //条件：有效值为1，标记值相同
        if(cache_set[i].valid == 1 && cache_set[i].tag == tag ){
            if(verbosity == 1)   printf("hit \n");
            hit_count++;
            //将其他的行的lru值加1
            for(int j = 0;j < E; j++){
                if(cache_set[j].valid && cache_set[j].lru < cache_set[i].lru){
                    cache_set[j].lru++;
                }
            }
            //将本行的lru置0
            cache_set[i].lru = 0;
            //如果命中，函数直接结束
            return;
        }
    }
    
    //未命中情况
    if(verbosity == 1)  printf("miss\n");
    miss_count++;
    
    int j;
    //初始化最大距离行标号，最大距离。
    int maxIndex = 0;
    unsigned long long maxLru = 0;
    //逐行检测，找到最大距离行的标号和距离
    for(j = 0;j < E && cache_set[j].valid;j++){
        if(cache_set[j].lru >= maxLru){
            maxIndex = j;
            maxLru = cache_set[j].lru;
        }
    }
    
    //如果cache中有空行
    if(j != E){
        //所有行的lru值+1，将本行lru置0，将空行的标记值设置为目标地址的标记值，有效值改为1，即将数据写入。
        for(int k = 0;k < E;k++){
            if(cache_set[k].valid)  cache_set[k].lru++;
        }
        cache_set[j].lru = 0;
        cache_set[j].valid = 1;
        cache_set[j].tag = tag;
    }
    //如果没有空行
    else{
        //需要替换，将替换次数加一
        if(verbosity)   printf("eviction \n");
        eviction_count++;
        //所有行的lru值+1，之后将本行lru置0，将标记值写为目标地址的标记值，有效值改为1，即将数据写入
        for(int k = 0;k < E;k++){
            cache_set[k].lru++;
        }
        cache_set[maxIndex].lru = 0;
        cache_set[maxIndex].tag = tag;
    }
}


/*
 * replayTrace - replays the given trace file against the cache 
 */
void replayTrace(char* trace_fn)
{
    char buf[1000];
    mem_addr_t addr=0;
    unsigned int len=0;
    FILE* trace_fp = fopen(trace_fn, "r");

    if(!trace_fp){
        fprintf(stderr, "%s: %s\n", trace_fn, strerror(errno));
        exit(1);
    }

    while( fgets(buf, 1000, trace_fp) != NULL) {
        if(buf[1]=='S' || buf[1]=='L' || buf[1]=='M') {
            sscanf(buf+3, "%llx,%u", &addr, &len);
      
            if(verbosity)
                printf("%c %llx,%u ", buf[1], addr, len);

            accessData(addr);

            /* If the instruction is R/W then access again */
            if(buf[1]=='M')
                accessData(addr);
            
            if (verbosity)
                printf("\n");
        }
    }
    
    fclose(trace_fp);
}

/*
 * printUsage - Print usage info
 */
void printUsage(char* argv[])
{
    printf("Usage: %s [-hv] -s <num> -E <num> -b <num> -t <file>\n", argv[0]);
    printf("Options:\n");
    printf("  -h         Print this help message.\n");
    printf("  -v         Optional verbose flag.\n");
    printf("  -s <num>   Number of set index bits.\n");
    printf("  -E <num>   Number of lines per set.\n");
    printf("  -b <num>   Number of block offset bits.\n");
    printf("  -t <file>  Trace file.\n");
    printf("\nExamples:\n");
    printf("  linux>  %s -s 4 -E 1 -b 4 -t traces/yi.trace\n", argv[0]);
    printf("  linux>  %s -v -s 8 -E 2 -b 4 -t traces/yi.trace\n", argv[0]);
    exit(0);
}

/*
 * main - Main routine 
 */
int main(int argc, char* argv[])
{
    char c;

    while( (c=getopt(argc,argv,"s:E:b:t:vh")) != -1){
        switch(c){
        case 's':
            s = atoi(optarg);
            break;
        case 'E':
            E = atoi(optarg);
            break;
        case 'b':
            b = atoi(optarg);
            break;
        case 't':
            trace_file = optarg;
            break;
        case 'v':
            verbosity = 1;
            break;
        case 'h':
            printUsage(argv);
            exit(0);
        default:
            printUsage(argv);
            exit(1);
        }
    }

    /* Make sure that all required command line args were specified */
    if (s == 0 || E == 0 || b == 0 || trace_file == NULL) {
        printf("%s: Missing required command line argument\n", argv[0]);
        printUsage(argv);
        exit(1);
    }

    /* Compute S, E and B from command line args */
    S = 1<<s;
    B = 1<<b;
 
    /* Initialize cache */
    initCache();

#ifdef DEBUG_ON
    printf("DEBUG: S:%u E:%u B:%u trace:%s\n", S, E, B, trace_file);
    printf("DEBUG: set_index_mask: %llu\n", set_index_mask);
#endif
 
    replayTrace(trace_file);

    /* Free allocated memory */
    freeCache();

    /* Output the hit and miss statistics for the autograder */
    printSummary(hit_count, miss_count, eviction_count);
    return 0;
}
