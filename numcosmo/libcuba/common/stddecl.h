/*
	stddecl.h
		declarations common to all Cuba routines
		last modified 23 Apr 15 th
*/


#ifndef _stddecl_h_
#define _stddecl_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/*
#define _BSD_SOURCE
#define _SVID_SOURCE
*/
#define _DEFAULT_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <unistd.h>
#include <assert.h>
#include <fcntl.h>
#include <setjmp.h>
#include <sys/stat.h>
#include <sys/types.h>
#ifdef HAVE_FORK
#include <sys/wait.h>
#include <sys/socket.h>
#include <signal.h>
#ifdef HAVE_SHMGET
#include <sys/ipc.h>
#include <sys/shm.h>
#endif
#endif

#ifdef HAVE_ALLOCA_H
#include <alloca.h>
#elif defined __GNUC__
#define alloca __builtin_alloca
#elif defined _AIX
#define alloca __alloca
#elif defined _MSC_VER
#include <malloc.h>
#define alloca _alloca
#else
#include <stddef.h>
#ifdef  __cplusplus
extern "C"
#endif
void *alloca (size_t);
#endif

#ifndef NDIM
#define NDIM t->ndim
#define MAXDIM 1024
#else
#define MAXDIM NDIM
#endif

#ifndef NCOMP
#define NCOMP t->ncomp
#define MAXCOMP 1024
#else
#define MAXCOMP NCOMP
#endif

#if defined(VEGAS) || defined(SUAVE)
#define VES_ONLY(...) __VA_ARGS__
#define NW 1
#else
#define VES_ONLY(...)
#define NW 0
#endif

#ifdef DIVONNE
#define DIV_ONLY(...) __VA_ARGS__
#else
#define DIV_ONLY(...)
#endif

#define SAMPLESIZE (NW + t->ndim + t->ncomp)*sizeof(real)


enum { uninitialized = 0x61627563 };

#define EnvInit(var, name, default) \
  if( var == uninitialized ) { \
    cchar *env = getenv(name); \
    if( env == NULL ) var = default; \
    else { \
      var = atoi(env); \
      if( cubaverb_ ) { \
        char out[64]; \
        sprintf(out, "env " name " = %d", (int)var); \
        Print(out); \
      } \
    } \
  }

#define VerboseInit() EnvInit(cubaverb_, "CUBAVERBOSE", 0)
#define MaxVerbose(flags) (flags + IDim(IMin(cubaverb_, 3) - ((flags) & 3)))

#define VERBOSE (t->flags & 3)
#define LAST (t->flags & 4)
#define SHARPEDGES (t->flags & 8)
#define KEEPFILE (t->flags & 16)
#define ZAPSTATE (t->flags & 32)
#define REGIONS (t->flags & 128)
#define RNG (t->flags >> 8)

#define INFTY DBL_MAX

#if __STDC_VERSION__ >= 199901L
#define POW2(n) 0x1p-##n
#else
#define POW2(n) ldexp(1., -n)
#endif

#define NOTZERO POW2(300)

#define ABORT -999

#define Elements(x) (sizeof(x)/sizeof(*x))

#define Copy(d, s, n) memcpy(d, s, (n)*sizeof(*(d)))

#define Move(d, s, n) memmove(d, s, (n)*sizeof(*(d)))

#define XCopy(d, s) Copy(d, s, t->ndim)

#define FCopy(d, s) Copy(d, s, t->ncomp)

#define Clear(d, n) memset(d, 0, (n)*sizeof(*(d)))

#define XClear(d) Clear(d, t->ndim)

#define FClear(d) Clear(d, t->ncomp)

#define Zap(d) memset(d, 0, sizeof(d))

#define MaxErr(avg) Max(t->epsrel*fabsx(avg), t->epsabs)

#ifdef __cplusplus
#define mallocset(p, n) (*(void **)&p = malloc(n))
#define reallocset(p, n) (*(void **)&p = realloc(p, n))
#else
#define mallocset(p, n) (p = malloc(n))
#define reallocset(p, n) (p = realloc(p, n))
#endif

#define Abort(s) abort1(s, __LINE__)
#define abort1(s, line) abort2(s, line)
#define abort2(s, line) { perror(s " " __FILE__ "(" #line ")"); exit(1); }

#define Die(p) if( (p) == NULL ) Abort("malloc")

#define MemAlloc(p, n) Die(mallocset(p, n))
#define ReAlloc(p, n) Die(reallocset(p, n))
#define Alloc(p, n) MemAlloc(p, (n)*sizeof(*p))

#if __STDC_VERSION__ >= 199901L
#define Sized(type, var, size) char var##_[size]; type *var = (type *)var##_
#define Vector(type, var, n1) type var[n1]
#define Array(type, var, n1, n2) type var[n1][n2]
#else
#define Sized(type, var, size) type *var = alloca(size)
#define Vector(type, var, n1) type *var = alloca((n1)*sizeof(type))
#define Array(type, var, n1, n2) type (*var)[n2] = alloca((n1)*(n2)*sizeof(type))
#endif

#define FORK_ONLY(...)
#define SHM_ONLY(...)
#define ShmAlloc(...)
#define ShmFree(...)

#ifdef MLVERSION
#define ML_ONLY(...) __VA_ARGS__
#define ML_NOT(...)
#else
#define ML_ONLY(...)
#define ML_NOT(...) __VA_ARGS__

#define CORE_MASTER (int []){32768}
#define MasterInit() do if( !cubafun_.init ) { \
  cubafun_.init = true; \
  if( cubafun_.initfun ) cubafun_.initfun(cubafun_.initarg, CORE_MASTER); \
} while( 0 )
#define MasterExit() do if( cubafun_.init ) { \
  cubafun_.init = false; \
  if( cubafun_.exitfun ) cubafun_.exitfun(cubafun_.exitarg, CORE_MASTER); \
} while( 0 )
#define Invalid(s) ((s) == NULL || *(int *)(s) == -1)

#ifdef HAVE_FORK
#undef FORK_ONLY
#define FORK_ONLY(...) __VA_ARGS__

#ifdef HAVE_SHMGET
#undef SHM_ONLY
#define SHM_ONLY(...) __VA_ARGS__

#define MasterAlloc(t) \
  t->shmid = shmget(IPC_PRIVATE, t->nframe*SAMPLESIZE, IPC_CREAT | 0600)
#define MasterFree(t) shmctl(t->shmid, IPC_RMID, NULL)
#define WorkerAlloc(t)
#define WorkerFree(r)

#undef ShmAlloc
#define ShmAlloc(t, who) \
  who##Alloc(t); \
  if( t->shmid != -1 ) { \
    t->frame = shmat(t->shmid, NULL, 0); \
    if( t->frame == (void *)-1 ) Abort("shmat"); \
  }

#undef ShmFree
#define ShmFree(t, who) \
  if( t->shmid != -1 ) { \
    shmdt(t->frame); \
    who##Free(t); \
  }

#endif
#endif
#endif
  
#define FrameAlloc(t, who) \
  SHM_ONLY(ShmAlloc(t, who) else) \
  MemAlloc(t->frame, t->nframe*SAMPLESIZE);

#define FrameFree(t, who) \
  DIV_ONLY(if( t->nframe )) { \
    SHM_ONLY(ShmFree(t, who) else) \
    free(t->frame); \
  }


#define StateDecl \
char *statefile_tmp = NULL, *statefile_XXXXXX = NULL; \
int statemsg = VERBOSE; \
ssize_t ini = 1; \
struct stat st

#define StateSetup(t) if( (t)->statefile ) { \
  if( *(t)->statefile == 0 ) (t)->statefile = NULL; \
  else { \
    ccount len = strlen((t)->statefile); \
    statefile_tmp = alloca(len + 8); \
    strcpy(statefile_tmp, (t)->statefile); \
    statefile_XXXXXX = statefile_tmp + len; \
  } \
}

typedef long long int signature_t;

enum { signature = 0x41425543 };

#define StateSignature(t, i) (signature + \
  ((signature_t)(i) << 60) + \
  ((signature_t)(t)->ncomp << 48) + \
  ((signature_t)(t)->ndim << 32))

#define StateReadTest(t) (t)->statefile && \
  stat((t)->statefile, &st) == 0 && (st.st_mode & 0400)

#define StateReadOpen(t, fd) do { \
  int fd; \
  if( (fd = open((t)->statefile, O_RDONLY)) != -1 ) { \
    do

#define StateRead(fd, buf, size) \
  ini += size - read(fd, buf, size)

#define StateReadClose(t, fd) \
    while( (--ini, 0) ); \
    close(fd); \
  } \
  if( ini | statemsg ) { \
    char s[512]; \
    sprintf(s, ini ? \
      "\nError restoring state from %s, starting from scratch." : \
      "\nRestored state from %s.", (t)->statefile); \
    Print(s); \
  } \
} while( 0 )


#define StateWriteTest(t) ((t)->statefile)

#define StateWriteOpen(t, fd) do { \
  ssize_t fail = 1; \
  int fd; \
  strcpy(statefile_XXXXXX, "-XXXXXX"); \
  if( (fd = mkstemp(statefile_tmp)) != -1 ) { \
    do

#define StateWrite(fd, buf, size) \
  fail += size - write(fd, buf, size)

#define StateWriteClose(t, fd) \
    while( (--fail, 0) ); \
    close(fd); \
    if( fail == 0 ) fail |= rename(statefile_tmp, (t)->statefile); \
  } \
  if( fail | statemsg ) { \
    char s[512]; \
    sprintf(s, fail ? \
      "\nError saving state to %s." : \
      "\nSaved state to %s.", (t)->statefile); \
    Print(s); \
    statemsg &= fail & -2; \
  } \
} while( 0 )


#define StateRemove(t) \
if( fail == 0 && (t)->statefile && KEEPFILE == 0 ) unlink((t)->statefile)


#ifdef __cplusplus
#define Extern extern "C"
#else
#define Extern extern
typedef enum { false, true } bool;
#endif

typedef const char cchar;

typedef const bool cbool;

typedef const int cint;

typedef const long clong;

typedef const size_t csize_t;

#define COUNT "%d"
typedef /*unsigned*/ int count;
typedef const count ccount;

#ifdef LONGLONGINT
#define PREFIX(s) ll##s
#define NUMBER "%lld"
#define NUMBER7 "%7lld"
#define NUMBER_MAX LLONG_MAX
typedef long long int number;
#else
#define PREFIX(s) s
#define NUMBER "%d"
#define NUMBER7 "%7d"
#define NUMBER_MAX INT_MAX
typedef int number;
#endif
typedef const number cnumber;

#define REAL "%g"
#define REALF "%f"
#define SHOW(r) (double)(r)
	/* floating-point numbers are printed with SHOW */

#if REALSIZE == 16
#include <quadmath.h>
typedef __float128 real;
#define RC(x) x##Q
#define sqrtx sqrtq
#define expx expq
#define powx powq
#define erfx erfq
#define fabsx fabsq
#define ldexpx ldexpq
#define REAL_MAX_EXP FLT128_MAX_EXP
#define REAL_MAX FLT128_MAX
#elif REALSIZE == 10
typedef long double real;
#define RC(x) x##L
#define sqrtx sqrtl
#define expx expl
#define powx powl
#define erfx erfl
#define fabsx fabsl
#define ldexpx ldexpl
#define REAL_MAX_EXP LDBL_MAX_EXP
#define REAL_MAX LDBL_MAX
#define MLPutRealxList MLPutReal128List
#define MLGetRealxList MLGetReal128List
#define MLReleaseRealxList MLReleaseReal128List
#else
typedef double real;
#define RC(x) x
#define sqrtx sqrt
#define expx exp
#define powx pow
#define erfx erf
#define fabsx fabs
#define ldexpx ldexp
#define REAL_MAX_EXP DBL_MAX_EXP
#define REAL_MAX DBL_MAX
#define MLPutRealxList MLPutReal64List
#define MLGetRealxList MLGetReal64List
#define MLReleaseRealxList MLReleaseReal64List
#endif

typedef const real creal;

typedef void (*subroutine)(void *, cint *);

typedef struct {
  subroutine initfun;
  void *initarg;
  subroutine exitfun;
  void *exitarg;
  bool init;
} coreinit;

typedef struct {
  int ncores, naccel;
  int pcores, paccel;
} corespec;

typedef struct {
  int fd, pid;
} fdpid;

typedef struct {
  corespec spec;
  fdpid fp[];
} Spin;


struct _this;

typedef struct {
  void (*worker)(struct _this *, csize_t, cint, cint);
  struct _this *thisptr;
  size_t thissize;
} dispatch;


typedef unsigned int state_t;

#define SOBOL_MINDIM 1
#define SOBOL_MAXDIM 40

/* length of state vector */
#define MERSENNE_N 624

/* period parameter */
#define MERSENNE_M 397

typedef struct {
  void (*getrandom)(struct _this *t, real *x);
  void (*skiprandom)(struct _this *t, cnumber n);
  union {
    struct {
      real norm;
      number v[SOBOL_MAXDIM][30], prev[SOBOL_MAXDIM];
      number seq;
    } sobol;
    struct {
      state_t state[MERSENNE_N];
      count next;
    } mersenne;
    struct {
      count n24, i24, j24, nskip;
      int carry, state[24];
    } ranlux;
  };
} RNGState;


#if NOUNDERSCORE
#define SUFFIX(s) s
#else
#define SUFFIX(s) s##_
#endif

#define EXPORT(s) EXPORT_(PREFIX(s))
#define EXPORT_(s) SUFFIX(s)


#define CString(cs, fs, len) { \
  char *_s = NULL; \
  if( fs ) { \
    int _l = len; \
    while( _l > 0 && fs[_l - 1] == ' ' ) --_l; \
    if( _l > 0 && (_s = alloca(_l + 1)) ) { \
      memcpy(_s, fs, _l); \
      _s[_l] = 0; \
    } \
  } \
  cs = _s; \
}

static inline real Sq(creal x) {
  return x*x;
}

static inline real Min(creal a, creal b) {
  return (a < b) ? a : b;
}

static inline real Max(creal a, creal b) {
  return (a > b) ? a : b;
}

static inline real Weight(creal sum, creal sqsum, cnumber n) {
  creal w = sqrtx(sqsum*n);
  return (n - 1)/Max((w + sum)*(w - sum), NOTZERO);
}


/* (a < 0) ? -1 : 0 */
#define NegQ(a) ((a) >> (sizeof(a)*8 - 1))

/* (a < 0) ? -1 : 1 */
#define Sign(a) (1 + 2*NegQ(a))

/* (a < 0) ? 0 : a */
#define IDim(a) ((a) & NegQ(-(a)))

/* (a < b) ? a : b */
#define IMin(a, b) ((a) - IDim((a) - (b)))

/* (a > b) ? a : b */
#define IMax(a, b) ((b) + IDim((a) - (b)))

/* (a == 0) ? 0 : -1 */
#define TrueQ(a) NegQ((a) | (-a))

/* a + (a == 0) */
#define Min1(a) ((a) + 1 + TrueQ(a))

/* abs(a) + (a == 0) */
#define Abs1(a) (((a) ^ NegQ(a)) - NegQ((a) - 1))


#ifdef MLVERSION

static inline void Print(MLCONST char *s)
{
  MLPutFunction(stdlink, "EvaluatePacket", 1);
  MLPutFunction(stdlink, "Print", 1);
  MLPutString(stdlink, s);
  MLEndPacket(stdlink);

  MLNextPacket(stdlink);
  MLNewPacket(stdlink);
}

#else

#define Print(s) do { puts(s); fflush(stdout); } while (0)

#endif

#endif

