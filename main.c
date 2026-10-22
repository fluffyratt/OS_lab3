#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdatomic.h>
#include <pthread.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include <inttypes.h>

static volatile uint64_t sink = 0;
static int quick = 0;
static FILE *out;
static double now_s(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + t.tv_nsec * 1e-9; }
static void *xalloc(size_t n) { void *p=NULL; if(posix_memalign(&p,64,n)!=0) {perror("allocation");exit(1);} return p; }
static uint64_t rng=0x9e3779b97f4a7c15ULL;
static uint64_t rand64(void){rng^=rng>>12;rng^=rng<<25;rng^=rng>>27;return rng*2685821657736338717ULL;}
static void record(const char *experiment,const char *variant,size_t size,int threads,int repeat,double seconds,double metric,const char *unit){
 fprintf(out,"%s,%s,%zu,%d,%d,%.9f,%.6f,%s\n",experiment,variant,size,threads,repeat,seconds,metric,unit);
 printf("%-16s %-21s size=%-9zu threads=%d run=%d time=%.5fs  %s=%.3f\n",experiment,variant,size,threads,repeat,seconds,unit,metric);
}
static int reps(void){return quick?2:5;}
static void benchmark_access(void){
 const size_t sizes[]={4096,32768,262144,2097152,16777216,67108864};
 for(size_t k=0;k<sizeof(sizes)/sizeof(sizes[0]);k++){
  size_t n=sizes[k]/sizeof(uint32_t);uint32_t *a=xalloc(n*sizeof(*a));size_t *next=xalloc(n*sizeof(*next));
  for(size_t i=0;i<n;i++){a[i]=(uint32_t)i;next[i]=i;}
  for(size_t i=n-1;i>0;i--){size_t j=rand64()%(i+1);size_t z=next[i];next[i]=next[j];next[j]=z;}
  size_t *links=xalloc(n*sizeof(*links));for(size_t i=0;i<n;i++)links[next[i]]=next[(i+1)%n];
  size_t steps=quick?200000:2000000;
  for(int r=0;r<reps();r++){
   uint64_t acc=0;double t=now_s();for(size_t j=0;j<steps;j++) acc+=a[j%n];double d=now_s()-t;sink^=acc;
   record("working_set","sequential",sizes[k],1,r,d,1e9*d/steps,"ns_per_access");
   size_t p=next[0];t=now_s();for(size_t j=0;j<steps;j++){p=links[p];}d=now_s()-t;sink^=p;
   record("working_set","random_pointer",sizes[k],1,r,d,1e9*d/steps,"ns_per_access");
  }
  free(a);free(next);free(links);
 }
}
static void benchmark_stride(void){
 const size_t bytes=quick?8*1024*1024:32*1024*1024;
 size_t n=bytes/sizeof(uint64_t);uint64_t *a=xalloc(bytes);for(size_t i=0;i<n;i++)a[i]=i+1;
 const size_t strides[]={1,2,4,8,16,32,64,128,256};
 size_t accesses=quick?300000:3000000;
 for(size_t x=0;x<sizeof(strides)/sizeof(strides[0]);x++)for(int r=0;r<reps();r++){
  size_t stride=strides[x],p=0;uint64_t acc=0;double t=now_s();
  for(size_t j=0;j<accesses;j++){acc+=a[p];p+=stride;if(p>=n)p-=n;}
  double d=now_s()-t;sink^=acc;
  char name[32];snprintf(name,sizeof(name),"stride_%zu_B",stride*sizeof(uint64_t));
  record("stride",name,bytes,1,r,d,1e9*d/accesses,"ns_per_access");
 }
 free(a);
}
static void benchmark_matrix(void){
 const int dims[]={64,128,256,384};int cnt=quick?3:4;
 for(int z=0;z<cnt;z++){
  int n=dims[z];size_t cells=(size_t)n*n;double *a=xalloc(cells*sizeof(double)),*b=xalloc(cells*sizeof(double)),*bt=xalloc(cells*sizeof(double)),*c=xalloc(cells*sizeof(double));
  for(size_t i=0;i<cells;i++){a[i]=(i%17)*0.1;b[i]=(i%13)*0.1;}
  for(int i=0;i<n;i++)for(int j=0;j<n;j++)bt[(size_t)j*n+i]=b[(size_t)i*n+j];
  for(int variant=0;variant<2;variant++)for(int r=0;r<reps();r++){
   memset(c,0,cells*sizeof(double));double t=now_s();
   if(!variant){for(int i=0;i<n;i++)for(int j=0;j<n;j++){double s=0;for(int k=0;k<n;k++)s+=a[(size_t)i*n+k]*b[(size_t)k*n+j];c[(size_t)i*n+j]=s;}}
   else {for(int i=0;i<n;i++)for(int j=0;j<n;j++){double s=0;for(int k=0;k<n;k++)s+=a[(size_t)i*n+k]*bt[(size_t)j*n+k];c[(size_t)i*n+j]=s;}}
   double d=now_s()-t;sink^=(uint64_t)(c[cells/2]*1000);
   record("matrix",variant?"transposed_B":"naive_ijk",n,1,r,d,2.0*n*n*n/(d*1e9),"GFLOP_per_s");
  }
  free(a);free(b);free(bt);free(c);
 }
}

typedef struct { _Alignas(64) atomic_uint_fast64_t counter; } padded_counter;
typedef struct {atomic_uint_fast64_t counter;} packed_counter;
static packed_counter packed[16];static padded_counter padded[16];
typedef struct {int id;size_t loops;int mode;int threads;} false_arg;
static void *false_worker(void *ptr){false_arg *a=ptr;atomic_uint_fast64_t *c=a->mode?&padded[a->id].counter:&packed[a->id].counter;
 for(size_t i=0;i<a->loops;i++){atomic_fetch_add_explicit(c,1,memory_order_relaxed);}
 return NULL;}
static void benchmark_sharing(void){
 const int counts[]={1,2,4,8};size_t loops=quick?150000:1500000;
 for(int ix=0;ix<4;ix++){int nt=counts[ix];for(int mode=0;mode<2;mode++)for(int r=0;r<reps();r++){
  pthread_t ts[16];false_arg args[16];for(int i=0;i<nt;i++){atomic_store(&packed[i].counter,0);atomic_store(&padded[i].counter,0);}
  double t=now_s();for(int i=0;i<nt;i++){args[i]=(false_arg){i,loops,mode,nt};if(pthread_create(&ts[i],NULL,false_worker,&args[i]))exit(2);}
  for(int i=0;i<nt;i++){pthread_join(ts[i],NULL);}
  double d=now_s()-t;
  uint64_t total=0;for(int i=0;i<nt;i++)total+=atomic_load(mode?&padded[i].counter:&packed[i].counter);
  if(total!=loops*(size_t)nt){fprintf(stderr,"counter mismatch\n");exit(3);}
  record("false_sharing",mode?"padded_64B":"adjacent",loops,nt,r,d,total/d/1e6,"Mops_per_s");
 }}
}
static atomic_uint_fast64_t shared_atomic;static pthread_mutex_t mu=PTHREAD_MUTEX_INITIALIZER;static uint64_t shared_mutex;
typedef struct {int mode;size_t loops;uint64_t local;} sync_arg;
static void *sync_worker(void *ptr){sync_arg *a=ptr;uint64_t local=0;
 for(size_t i=0;i<a->loops;i++){
  switch(a->mode){
   case 0:local++;break;
   case 1:atomic_fetch_add_explicit(&shared_atomic,1,memory_order_relaxed);break;
   case 2:pthread_mutex_lock(&mu);shared_mutex++;pthread_mutex_unlock(&mu);break;
   case 3:{uint64_t v=atomic_load_explicit(&shared_atomic,memory_order_relaxed);atomic_store_explicit(&shared_atomic,v+1,memory_order_relaxed);break;}
  }
 }
 a->local=local;return NULL;}
static void benchmark_sync(void){
 const int nts[]={1,2,4,8};const char *names[]={"thread_local","atomic_fetch_add","mutex","load_then_store"};size_t loops=quick?25000:250000;
 for(int x=0;x<4;x++){int nt=nts[x];for(int mode=0;mode<4;mode++)for(int r=0;r<reps();r++){
  pthread_t ts[8];sync_arg args[8];atomic_store(&shared_atomic,0);shared_mutex=0;
  double t=now_s();for(int i=0;i<nt;i++){args[i]=(sync_arg){mode,loops,0};if(pthread_create(&ts[i],NULL,sync_worker,&args[i]))exit(2);}
  uint64_t total=0;for(int i=0;i<nt;i++){pthread_join(ts[i],NULL);total+=args[i].local;}
  double d=now_s()-t;
  if(mode==1||mode==3){total=atomic_load(&shared_atomic);}
  if(mode==2){total=shared_mutex;}
  record("sync",names[mode],loops,nt,r,d,100.0*total/(loops*(size_t)nt),"correct_percent");
 }}
}
static void benchmark_volatile(void){
 size_t loops=quick?100000:1000000;
 for(int r=0;r<reps();r++){
  volatile uint64_t x=0;double t=now_s();for(size_t i=0;i<loops;i++)x++;double d=now_s()-t;sink^=x;
  record("volatile","single_thread_only",loops,1,r,d,1e9*d/loops,"ns_per_increment");
 }
}
int main(int argc,char **argv){
 for(int i=1;i<argc;i++){if(strcmp(argv[i],"--quick")==0)quick=1;else{fprintf(stderr,"Usage: %s [--quick]\n",argv[0]);return 1;}}
 out=fopen("results/raw.csv","w");if(!out){perror("results/raw.csv (run from project directory)");return 1;}
 fprintf(out,"experiment,variant,size,threads,repeat,seconds,metric,unit\n");
 benchmark_access();benchmark_stride();benchmark_matrix();benchmark_sharing();benchmark_sync();benchmark_volatile();
 fclose(out);printf("Done. CSV: results/raw.csv; guard=%" PRIu64 "\n",(uint64_t)sink);return 0;
}
