#ifndef UTIL_H_INCLUDED
#define UTIL_H_INCLUDED
#include <cstddef>
#include <cstdio>
#include <stdint.h>
#include <string.h>
//#include <cctype>
#define MB                          1024*1024


void die(void);
void * checkalloc(size_t count, size_t size);
void * checkrealloc(void * ptr, size_t count, size_t size);
void trim_string(char * s);
void fill_string(char * s,size_t size);
char yesno(int x);
void strlower(char * s);
char * read_multiline(FILE * input);
int count_words(char * str);

unsigned int digital_crc32(const unsigned char *buf, size_t len);
#if defined(PARALLEL) || defined(EXCHANGE)
void parallel_init(int * mynod, int * numnod, char * outputfmt);
void parallel_finish();
#endif
#if defined(__unix__)
double convtime(struct timeval time);
#endif

//type and size of bitmap entries in "subset"
#define BIT_ENTRY      uint64_t
//the size of a BIT_ENTRY must be 2^SHIFT bits
#define SHIFT          6
//#define ENTRY_SIZE     1<<SHIFT
#define BIT_ENTRY_SIZE sizeof(BIT_ENTRY)
#define ENTRY_SIZE     BIT_ENTRY_SIZE*8
//this allows us to divide by the number of bits using &
#define MASK           ENTRY_SIZE-1

class subset {
    unsigned long int n, nentries;
    BIT_ENTRY * bits;
public:
    subset() {
        n=0;
        nentries=0;
        bits=NULL;
    }
    subset(unsigned long int _n) {
        init(_n);
    }
    subset(const subset& other) {
        n=0;
        nentries=0;
        bits=NULL;
        (*this)=other;  //use the overloaded operator = below to do a deep copy
    }
    ~subset() {
        if (bits!=NULL) free(bits);
    }
    inline void init(unsigned long int _n) {
        n=_n;
        nentries=(n>>SHIFT)+1;
        bits = (BIT_ENTRY *) checkalloc(nentries,BIT_ENTRY_SIZE);
        clear();
    }
    inline void clear(void) {
        memset(bits,0,nentries*BIT_ENTRY_SIZE);
    }
    //copy constructor
    inline subset& operator=(const subset& other) {
        //check for self assignment to avoid unnecessarily deallocating memory
        if (this!=&other) {
            if (n!=other.n) {
                if (bits!=NULL) free(bits);
                init(other.n);
            }
            memcpy(bits,other.bits,nentries*BIT_ENTRY_SIZE);
        }
        return *this;
    }
    //set membership
    inline bool operator [](const unsigned long int idx) {
        if (idx>=n) return false;
        long int entry=idx>>SHIFT;
        //This must be typed BIT_ENTRY so that the expression 1<<bit comes out with the right type.
        BIT_ENTRY bit=idx&MASK;
        bit=((BIT_ENTRY) 1)<<bit;
        if ((bits[entry])&bit) return true;
        return false;
    }
    //add member
    inline subset& operator +=(const unsigned long int idx) {
        if (idx>=n) return *this;
        long int entry=idx>>SHIFT;
        BIT_ENTRY bit=idx&MASK;
        bit=((BIT_ENTRY) 1)<<bit;
                //if ((bits[b])&(1<<bit)) return; //ensures we don't double add edge to connections
        bits[entry]|=bit;
        return *this;
    }
    //remove member
    inline subset& operator -=(const unsigned long int idx) {
        if (idx>=n) return *this;
        long int entry=idx>>SHIFT;
        BIT_ENTRY bit=idx&MASK;
        bit=((BIT_ENTRY) 1)<<bit;
        //if ((bits[b])&(1<<bit)) return; //ensures we don't double add edge to connections
        bits[entry]&=(~bit);
        return *this;
    }
    //set intersection (&=), set removal (/=)
    inline subset& operator &=(const subset& other) {
        //long int minentries = min(nentries,other.nentries);
        long int minentries=nentries;
        if (other.nentries<minentries) minentries=other.nentries;
        for (long int i=0; i<minentries; i++) {
            bits[i]&=other.bits[i];
        }
        return *this;
    }
    //set union (|=)
    inline subset& operator |=(const subset& other) {
        //long int minentries = min(nentries,other.nentries);
        long int minentries=nentries;
        if (other.nentries<minentries) minentries=other.nentries;
        for (long int i=0; i<minentries; i++) {
            bits[i]|=other.bits[i];
        }
        return *this;
    }
    //set removal (/=)
    inline subset& operator /=(const subset& other) {
        //long int minentries = min(nentries,other.nentries);
        long int minentries=nentries;
        if (other.nentries<minentries) minentries=other.nentries;
        for (long int i=0; i<minentries; i++) {
            bits[i]&=~other.bits[i];
        }
        return *this;
    }
    inline bool is_empty(void){
        long int i;
        for (i=0; i<nentries; i++) if (bits[i]!=0) return false;
        return true;
    }
    inline void print(char * name) {
        long int i;
        printf("%s = ",name);
        for (i=0; i<n; i++) if ((*this)[i]) printf("%ld,",i+1);
        printf("\n");
    }
    /*inline void to_clique_info(clique_info * info) {
        int i;
        info->count=0;
        info->vertices=NULL;
        for (i=0; i<n; i++) if ((*this)[i]) {
            info->vertices=(long int *) checkrealloc(info->vertices,info->count+1,sizeof(long int));
            info->vertices[info->count]=i;
            info->count++;
        }
    }*/
    void parse_int_list(char * str);
};





#ifdef TIMERS

#define TIMER_NONE           -1
#define TIMER_MC_MOVE        0
#define TIMER_NT_BONDS       1
#define TIMER_NT_ANGLES      2
#define TIMER_NT_DIHEDRALS   3
#define TIMER_NT_IMPROPERS   4
#define TIMER_NT_VDW_ELEC    5
//#define TIMER_COV_TABLES     6
#define TIMER_GO             6
/*#define TIMER_CHECK_CUTOFF   8
#define TIMER_INT_EXACT_PREP
#define TIMER_INT_EXACT      9
#define TIMER_INT_PREP       10
#define TIMER_INT_ORIENT     11
#define TIMER_INT_TRANS      12
#define TIMER_INT_INDEX      13
#define TIMER_INT_LOOKUP     14*/
#define TIMER_INT_OTHER      7
#define TIMER_NB_LIST        8
#define TIMER_OTHER          9
#define NTIMERS              10

struct timer {
    unsigned long long int stop_count;
    unsigned long long int last_started;
    unsigned long long int total_ticks;
};

static struct timer timers[NTIMERS];
static volatile int current_timer;
static const char * timer_names[NTIMERS] = {"MC moves","Bonds","Angles","Dihedrals","Impropers","Non tab. VDW/Elec","Go potential","Other interaction","Nonbond list update","Other"};

/*"Backbone tables","Check cutoff",
    "Exact prep","Exact interaction","Table prep","Table orientational","Table translational","Table index calc","Table lookup","Other interaction",*/

static double overhead;

void switch_timer(int timer);
void init_timers(void);
void stop_current_timer(void);
void print_timers(void);

#endif //TIMERS
#endif // UTIL_H_INCLUDED
