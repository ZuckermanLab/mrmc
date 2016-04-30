#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "string.h"
#include "mc.h"
//#include "tables.h"
//#include "covalent_tables.h"
#include "rotations.h"
//#include "fragments.h"
#include "util.h"
#include "ffield.h"

#ifdef __unix__
#define _GNU_SOURCE   // gives us feenableexcept on older gcc's
#define __USE_GNU     // gives us feenableexcept on newer gcc's
#include <fenv.h>
#include <execinfo.h>
#include <signal.h>
#include <unistd.h>
#else
#include <float.h>
//#ifndef _EM_OVERFLOW
//#define _EM_OVERFLOW EM_OVERFLOW
//#endif
#endif

#ifdef __unix__
#define MAX_BACKTRACE 20
//Signal handler to generate backtrace if we get a segmentation fault or floating point exception
void signal_handler(int sig) {
  void *array[MAX_BACKTRACE];
  size_t size;
  char msg[255];
  // get void*'s for all entries on the stack
  size = backtrace(array, MAX_BACKTRACE);

  // print out all the frames to stderr
  fprintf(stderr, "%s:\n", sys_siglist[sig]);
  //psignal(sig,msg);
  //fprintf(stderr,"%s\n",msg);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  fflush(stderr);
  exit(1);
}
#endif


using namespace std;

int main(int argc, char * argv[])
{
    int imoved,i,mynod,numnod;
    //double en,en2,de;
    char fname[255],buffer[255];
    time_t now;

    /*table * newtable;
    table * newdectable;
    covalent_table * newcovtable;*/
    simulation * sim;
#ifdef __unix__
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    signal(SIGSEGV,signal_handler);
    signal(SIGFPE,signal_handler);
#else
//    _controlfp(0, _EM_ZERODIVIDE | _EM_INVALID | _EM_OVERFLOW);
#endif
        if (argc<2) {
        printf("Syntax: tablemc input_file\n");
        die();
    }
#if defined(PARALLEL) || defined(EXCHANGE)
    parallel_init(&mynod,&numnod,argv[2]);
#endif
    time(&now);
    strncpy(buffer,ctime(&now),sizeof(buffer));
    printf("Justin's MRMC starting at %s\n",buffer);
    printf("Executable name: %s\n",argv[0]);
    fflush(stdout);
    sim=new simulation();
    sim->process_commands(argv[1]);
    delete sim;
    time(&now);
    strncpy(buffer,ctime(&now),sizeof(buffer));
    printf("Finished at %s\n",buffer);
    fflush(stdout);
#if defined(PARALLEL) || defined(EXCHANGE)
    parallel_finish();
#endif
    return 0;
}
