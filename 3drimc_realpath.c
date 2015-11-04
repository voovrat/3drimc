#include <limits.h>
#include <stdlib.h>

main(int argc,char *argv[])
{
 
char fullpath[2048];

if(argc<=1)
{
  puts("usage:  3drimc_realpath path");
  return 1;
}

realpath(argv[1],fullpath);
puts(fullpath);

 
return 0;
}
