#ifndef UTIL_H_INCLUDED
#define UTIL_H_INCLUDED
#include <cstddef>
#include <cstdio>
#define MB                          1024*1024


void die(void);
void * checkalloc(size_t count, size_t size);
void * checkrealloc(void * ptr, size_t count, size_t size);
void trim_string(char * s);
void fill_string(char * s,size_t size);
char yesno(int x);
void strlower(char * s);
char * read_multiline(FILE * input);
unsigned int digital_crc32(const unsigned char *buf, size_t len);
#if defined(PARALLEL) || defined(EXCHANGE)
void parallel_init(int * mynod, int * numnod, char * outputfmt);
void parallel_finish();
#endif
#if defined(__unix__)
double convtime(struct timeval time);
#endif

#ifdef TIMERS

#define TIMER_NONE           -1
#define TIMER_MC_MOVE        0
#define TIMER_NT_BONDS       1
#define TIMER_NT_ANGLES      2
#define TIMER_NT_DIHEDRALS   3
#define TIMER_NT_IMPROPERS   4
#define TIMER_NT_VDW_ELEC    5
#define TIMER_COV_TABLES     6
#define TIMER_CHECK_CUTOFF   7
#define TIMER_INT_EXACT_PREP 8
#define TIMER_INT_EXACT      9
#define TIMER_INT_PREP       10
#define TIMER_INT_ORIENT     11
#define TIMER_INT_TRANS      12
#define TIMER_INT_INDEX      13
#define TIMER_INT_LOOKUP     14
#define TIMER_INT_OTHER      15
#define TIMER_NB_LIST        16
#define TIMER_OTHER          17
#define NTIMERS              18

struct timer {
    unsigned long long int stop_count;
    unsigned long long int last_started;
    unsigned long long int total_ticks;
};

static struct timer timers[NTIMERS];
static volatile int current_timer;
static const char * timer_names[NTIMERS] = {"MC moves","Bonds","Angles","Dihedrals","Impropers","Non tab. VDW/Elec","Backbone tables","Check cutoff",
    "Exact prep","Exact interaction","Table prep","Table orientational","Table translational","Table index calc","Table lookup","Other interaction","Nonbond list update","Other"};
static double overhead;

void switch_timer(int timer);
void init_timers(void);
void stop_current_timer(void);
void print_timers(void);

#endif //TIMERS
#endif // UTIL_H_INCLUDED
