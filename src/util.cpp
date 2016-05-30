#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include "util.h"
#if defined(__unix__) && defined(TIMERS)
//needed for getpid and setaffinity
#include <sched.h>
#include <unistd.h>
#endif
#if defined(PARALLEL) || defined(EXCHANGE)
#include <mpi.h>
#endif

//Unified error-handling exit routine.
#if defined(PARALLEL) || defined(EXCHANGE)
void parallel_init(int * mynod, int * numnod,char * outputfmt)
{
    int is_parallel;
    char outfname[255];
    MPI_Initialized(&is_parallel);
    if (!is_parallel) {
        MPI_Init(NULL,NULL);
        MPI_Comm_rank(MPI_COMM_WORLD,mynod);
        MPI_Comm_size(MPI_COMM_WORLD,numnod);
        snprintf(outfname,sizeof(outfname),outputfmt,(*mynod)+1);
        trim_string(outfname);
        freopen(outfname,"w",stdout);
        printf("Node %d of %d\n",(*mynod)+1,*numnod);
        fflush(stdout);
    }
}

void parallel_finish(void)
{
    int is_parallel;
    MPI_Initialized(&is_parallel);
    if (is_parallel) MPI_Finalize();
}

void die(void)
{
    int is_parallel;
    MPI_Initialized(&is_parallel);
    if (is_parallel) {
        MPI_Abort(MPI_COMM_WORLD,1);
    } else {
        exit(1);
    }
}
#else
void die(void)
{
    exit(1);
}
#endif

#ifdef __unix__
double convtime(struct timeval time)
{
    return (double) time.tv_sec + ((double) time.tv_usec) * 1e-6;
}
#endif

void * checkalloc(size_t count, size_t size) {
    void * p;
    p=calloc(count,size);
    if (p==NULL) {
        printf("Memory allocation error.\n");
        die();
    }
    return p;
}

void * checkrealloc(void * ptr, size_t count, size_t size) {
    void * ptr2;
    ptr2=realloc(ptr,count*size);
    if (ptr2==NULL) {
        printf("Memory allocation error.\n");
        die();
    }
    return ptr2;
}

void trim_string(char * s)
{
    int i,l;
    for (i=0; i<strlen(s); i++)
        if (isspace(s[i])) s[i]='\0';
}

void fill_string(char * s,size_t size)
{
    int i,j,l;
    l=strlen(s);
    if (l>size) l=size;
    for (i=0; i<l; i++)
        if (isspace(s[i])) break;
    for (j=i; j<size; j++)
        s[j]=' ';
}

char yesno(int x)
{
    if (x) return 'Y'; else return 'N';
}

void strlower(char * s)
{
    int i;
    for (i=0; i<strlen(s); i++) s[i]=tolower(s[i]);
}
/*void trim_leading_whitespace(char * s1, char * s2, char * n2)
{
    char * p;
    p=s1;
    while (isspace(*p)) p++;
    //p now points to the first non-whitespace character
    strncpy(*/

//Read multiple lines.  Used for reading sequence from input file
char * read_multiline(FILE * input)
{
    char * seq;
    char * pch;
    bool endseq;
    int lseq;
    char buffer[255];
    memset(buffer,0,sizeof(buffer));
    seq=(char *) checkalloc(1,sizeof(char));
    *seq='\0';
    lseq=0;
    endseq=false;
    do {
        fgets(buffer,sizeof(buffer),input);
        lseq+=strlen(buffer)+1;
        pch=strstr(buffer,"\n");
        strncpy(pch," \0",2); //this could theoretically overflow if "pch" is too closse to the end of the buffer
        pch=strstr(buffer,"END");
        if (pch!=NULL) {
            endseq=true;
            *pch='\0';
        }
        seq=(char *) checkrealloc(seq,lseq,sizeof(char));
        strncat(seq,buffer,strlen(buffer)+1);
    } while (!endseq);
    return seq;
}



unsigned long read_random_seed(void)
{
    FILE * random;     
    unsigned long seed;           
    //use urandom to avoid blocking
    random=fopen("/dev/urandom","rb");
    fread(&seed,sizeof(seed),1,random);     
    fclose(random);
    return seed;
}

//used to pre-count the number of residues in order to instantiate the all-atom region early
//I got this from stackoverflow.com
int count_words(char * str) {
    if (str==NULL) return 0;
    bool inSpaces = true;
   int numWords = 0;
   //this is potentially unsafe, could run off the end of the string
   while (*str != '\0')
   {
      if (std::isspace(*str))
      {
         inSpaces = true;
      }
      else if (inSpaces)
      {
         numWords++;
         inSpaces = false;
      }

      ++str;
   }

   return numWords;
}


//Parse a string containing something of the form 1,3,5,7-9 and set the corresponding flags.
//Integers mentioned go from 1 to count. flags is zero-based, so 1 corresponds to flags[0], etc.
void subset::parse_int_list(char * str)
{
    char * token;
    char * tokenend;
    char * hyphen;
    //here "n" is the number of possible items in the subset
    int start,end,i,m;
    //for (i=0; i<count; i++) flags[i]=false;
    init(n); //zero out any previous set
    if ((strlen(str)==0) || (strstr(str,"NONE")!=NULL) || (strstr(str,"none")!=NULL)) return; //empty set 
    if ((strstr(str,"ALL")!=NULL) || (strstr(str,"all")!=NULL)) { //universal set
        for (i=0; i<n; i++) (*this)+=i;
        return;
    }
    token=str;
    while (true) {
        tokenend=strchr(token,',');
        hyphen=strchr(token,'-');
        if ((hyphen!=NULL) && ((tokenend==NULL) || (hyphen<tokenend))) { //this token is of the form m-n
            m=sscanf(token,"%d-%d",&start,&end);
            if ((m<2) || (start<1) || (start>n) || (end<1) || (end>n)) {
                printf("Invalid selection.\n");
                exit(1);
            }
            start--; //convert to zero based
            end--;
            for (i=start; i<=end; i++) (*this)+=i;
        } else { //it's a single integer
            m=sscanf(token,"%d",&i);
            if ((m<1) || (i<1) || (i>n)) {
                printf("Invalid selection.\n");
                exit(1);
            }
            i--; //zero-based
            (*this)+=i;
        }
        if (tokenend==NULL) break; //no more tokens
        token=tokenend+1; //advance past the comma
    }
}

//This code from http://gnuradio.org/redmine/projects/gnuradio/repository/revisions/1cb52da49230c64c3719b4ab944ba1cf5a9abb92/entry/gr-digital/lib/digital_crc32.cc


// Automatically generated CRC function
// polynomial: 0x104C11DB7
unsigned int digital_update_crc32(unsigned int crc, const unsigned char *data, size_t len)
{
    static const unsigned int table[256] = {
    0x00000000U,0x04C11DB7U,0x09823B6EU,0x0D4326D9U,
    0x130476DCU,0x17C56B6BU,0x1A864DB2U,0x1E475005U,
    0x2608EDB8U,0x22C9F00FU,0x2F8AD6D6U,0x2B4BCB61U,
    0x350C9B64U,0x31CD86D3U,0x3C8EA00AU,0x384FBDBDU,
    0x4C11DB70U,0x48D0C6C7U,0x4593E01EU,0x4152FDA9U,
    0x5F15ADACU,0x5BD4B01BU,0x569796C2U,0x52568B75U,
    0x6A1936C8U,0x6ED82B7FU,0x639B0DA6U,0x675A1011U,
    0x791D4014U,0x7DDC5DA3U,0x709F7B7AU,0x745E66CDU,
    0x9823B6E0U,0x9CE2AB57U,0x91A18D8EU,0x95609039U,
    0x8B27C03CU,0x8FE6DD8BU,0x82A5FB52U,0x8664E6E5U,
    0xBE2B5B58U,0xBAEA46EFU,0xB7A96036U,0xB3687D81U,
    0xAD2F2D84U,0xA9EE3033U,0xA4AD16EAU,0xA06C0B5DU,
    0xD4326D90U,0xD0F37027U,0xDDB056FEU,0xD9714B49U,
    0xC7361B4CU,0xC3F706FBU,0xCEB42022U,0xCA753D95U,
    0xF23A8028U,0xF6FB9D9FU,0xFBB8BB46U,0xFF79A6F1U,
    0xE13EF6F4U,0xE5FFEB43U,0xE8BCCD9AU,0xEC7DD02DU,
    0x34867077U,0x30476DC0U,0x3D044B19U,0x39C556AEU,
    0x278206ABU,0x23431B1CU,0x2E003DC5U,0x2AC12072U,
    0x128E9DCFU,0x164F8078U,0x1B0CA6A1U,0x1FCDBB16U,
    0x018AEB13U,0x054BF6A4U,0x0808D07DU,0x0CC9CDCAU,
    0x7897AB07U,0x7C56B6B0U,0x71159069U,0x75D48DDEU,
    0x6B93DDDBU,0x6F52C06CU,0x6211E6B5U,0x66D0FB02U,
    0x5E9F46BFU,0x5A5E5B08U,0x571D7DD1U,0x53DC6066U,
    0x4D9B3063U,0x495A2DD4U,0x44190B0DU,0x40D816BAU,
    0xACA5C697U,0xA864DB20U,0xA527FDF9U,0xA1E6E04EU,
    0xBFA1B04BU,0xBB60ADFCU,0xB6238B25U,0xB2E29692U,
    0x8AAD2B2FU,0x8E6C3698U,0x832F1041U,0x87EE0DF6U,
    0x99A95DF3U,0x9D684044U,0x902B669DU,0x94EA7B2AU,
    0xE0B41DE7U,0xE4750050U,0xE9362689U,0xEDF73B3EU,
    0xF3B06B3BU,0xF771768CU,0xFA325055U,0xFEF34DE2U,
    0xC6BCF05FU,0xC27DEDE8U,0xCF3ECB31U,0xCBFFD686U,
    0xD5B88683U,0xD1799B34U,0xDC3ABDEDU,0xD8FBA05AU,
    0x690CE0EEU,0x6DCDFD59U,0x608EDB80U,0x644FC637U,
    0x7A089632U,0x7EC98B85U,0x738AAD5CU,0x774BB0EBU,
    0x4F040D56U,0x4BC510E1U,0x46863638U,0x42472B8FU,
    0x5C007B8AU,0x58C1663DU,0x558240E4U,0x51435D53U,
    0x251D3B9EU,0x21DC2629U,0x2C9F00F0U,0x285E1D47U,
    0x36194D42U,0x32D850F5U,0x3F9B762CU,0x3B5A6B9BU,
    0x0315D626U,0x07D4CB91U,0x0A97ED48U,0x0E56F0FFU,
    0x1011A0FAU,0x14D0BD4DU,0x19939B94U,0x1D528623U,
    0xF12F560EU,0xF5EE4BB9U,0xF8AD6D60U,0xFC6C70D7U,
    0xE22B20D2U,0xE6EA3D65U,0xEBA91BBCU,0xEF68060BU,
    0xD727BBB6U,0xD3E6A601U,0xDEA580D8U,0xDA649D6FU,
    0xC423CD6AU,0xC0E2D0DDU,0xCDA1F604U,0xC960EBB3U,
    0xBD3E8D7EU,0xB9FF90C9U,0xB4BCB610U,0xB07DABA7U,
    0xAE3AFBA2U,0xAAFBE615U,0xA7B8C0CCU,0xA379DD7BU,
    0x9B3660C6U,0x9FF77D71U,0x92B45BA8U,0x9675461FU,
    0x8832161AU,0x8CF30BADU,0x81B02D74U,0x857130C3U,
    0x5D8A9099U,0x594B8D2EU,0x5408ABF7U,0x50C9B640U,
    0x4E8EE645U,0x4A4FFBF2U,0x470CDD2BU,0x43CDC09CU,
    0x7B827D21U,0x7F436096U,0x7200464FU,0x76C15BF8U,
    0x68860BFDU,0x6C47164AU,0x61043093U,0x65C52D24U,
    0x119B4BE9U,0x155A565EU,0x18197087U,0x1CD86D30U,
    0x029F3D35U,0x065E2082U,0x0B1D065BU,0x0FDC1BECU,
    0x3793A651U,0x3352BBE6U,0x3E119D3FU,0x3AD08088U,
    0x2497D08DU,0x2056CD3AU,0x2D15EBE3U,0x29D4F654U,
    0xC5A92679U,0xC1683BCEU,0xCC2B1D17U,0xC8EA00A0U,
    0xD6AD50A5U,0xD26C4D12U,0xDF2F6BCBU,0xDBEE767CU,
    0xE3A1CBC1U,0xE760D676U,0xEA23F0AFU,0xEEE2ED18U,
    0xF0A5BD1DU,0xF464A0AAU,0xF9278673U,0xFDE69BC4U,
    0x89B8FD09U,0x8D79E0BEU,0x803AC667U,0x84FBDBD0U,
    0x9ABC8BD5U,0x9E7D9662U,0x933EB0BBU,0x97FFAD0CU,
    0xAFB010B1U,0xAB710D06U,0xA6322BDFU,0xA2F33668U,
    0xBCB4666DU,0xB8757BDAU,0xB5365D03U,0xB1F740B4U,
    };

    while (len > 0)
    {
      crc = table[*data ^ ((crc >> 24) & 0xff)] ^ (crc << 8);
      data++;
      len--;
    }
    return crc;
}

unsigned int digital_crc32(const unsigned char *buf, size_t len)
{
  return digital_update_crc32(0xffffffff, buf, len) ^ 0xffffffff;
}

#ifdef TIMERS


#if defined(__i386__)

static __inline__ unsigned long long rdtsc(void)
{
  unsigned long long int x;
     __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
     return x;
}
#elif defined(__x86_64__)


static __inline__ unsigned long long rdtsc(void)
{
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

#elif defined(__powerpc__)


static __inline__ unsigned long long rdtsc(void)
{
  unsigned long long int result=0;
  unsigned long int upper, lower,tmp;
  __asm__ volatile(
                "0:                  \n"
                "\tmftbu   %0           \n"
                "\tmftb    %1           \n"
                "\tmftbu   %2           \n"
                "\tcmpw    %2,%0        \n"
                "\tbne     0b         \n"
                : "=r"(upper),"=r"(lower),"=r"(tmp)
                );
  result = upper;
  result = result<<32;
  result = result|lower;

  return(result);
}

#endif //processor type


void switch_timer(int timer)
{
    unsigned long long int now = rdtsc();
    //printf("Switching to timer %s at %lld\n",timer_names[timer],now);
    if (current_timer>=0) { //stop current timer
        timers[current_timer].total_ticks+=(now-timers[current_timer].last_started);
        //printf("Stopping timer %s, time %lld, total time %lld, count %lld\n",timer_names[current_timer],(now-timers[current_timer].last_started),
        //     timers[current_timer].total_ticks,timers[current_timer].stop_count);
        timers[current_timer].stop_count++;
    }
    current_timer=timer;
    if (timer>=0) timers[timer].last_started=now;
    //printf("Started timer %s at %lld\n",timer_names[current_timer],timers[current_timer].last_started);
}

#define OVERHEAD_TEST 1000000
void init_timers(void)
{
    int itimer,itest;
#ifdef __unix__
//This prevents inaccuracies from switching cores by fixing us to a single cpu.
//we are a little nicer in that we try to fix to the cpu we were running on.
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(sched_getcpu(), &cpuset);
    sched_setaffinity (getpid(), sizeof(cpuset), &cpuset);
    printf("Timer initialization: Fixed process to cpu %d\n",sched_getcpu());
#endif
    for (itimer=0; itimer<NTIMERS; itimer++) {
        timers[itimer].stop_count=0;
        timers[itimer].total_ticks=0;
    }
    current_timer=TIMER_NONE;
    //overhead test
    for (itest=1; itest<=OVERHEAD_TEST; itest++) switch_timer(TIMER_OTHER);
    switch_timer(TIMER_NONE);
    overhead=((double)timers[TIMER_OTHER].total_ticks)/((double)OVERHEAD_TEST);
    printf("Timers initialized. Overhead %.2f cycles\n",overhead);
    timers[TIMER_OTHER].stop_count=0;
    timers[TIMER_OTHER].total_ticks=0;
}




void print_timers(void)
{
    unsigned long long int grand_total;
    double avg_ticks, frac;
    int itimer;
    grand_total=0;
    for (itimer=0; itimer<NTIMERS; itimer++) {
        //overhead correction
        timers[itimer].total_ticks-=(unsigned long long int) (overhead*timers[itimer].stop_count);
        grand_total+=timers[itimer].total_ticks;
    }
    printf("Total cycles (overhead corrected): %lld\n",grand_total);
    //printf("Name                         Count         Total cycles (10^9)   Avg. cycles %%\n");
    printf("%-22s %15s %15s %15s %7s\n","Name","Count","Total cycles","Avg. cycles","Percent");
    for (itimer=0; itimer<NTIMERS; itimer++) if (timers[itimer].total_ticks>0) {
        avg_ticks=((double) timers[itimer].total_ticks)/((double) timers[itimer].stop_count);
        frac=((double) timers[itimer].total_ticks)/((double) grand_total);
        printf("%-22s %15lld %15lld %15.2f %7.2f\n",timer_names[itimer],timers[itimer].stop_count,timers[itimer].total_ticks,
            avg_ticks,frac*100);
    }
}

#endif //TIMERS

