
#ifdef _CRAY
  #include <mpp/shmem.h>
#else
  #include <shmem.h> 
#endif


int main(int argc, char **argv)
{
  start_pes(0);
  return 0;
}
