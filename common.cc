#include <string>
#include <stdexcept>
#include <fitsio.h>

#include "common.hh"

void check_fitsio_status(int status)
{
  if(status != 0)
    {
      std::string comb;
      char msg[81];
      msg[0] = '\0';

      for(;;)
        {
          int ret = fits_read_errmsg(msg);
          if(ret==0) break;
          comb += std::string(msg) + '\n';
        }
      
      throw std::runtime_error(comb);
    }
}
