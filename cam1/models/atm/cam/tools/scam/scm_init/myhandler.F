      integer function myhandler( sig, code, context )
#include <params.h>
#if ( defined sun )
      integer sig, code, context(5)
      myhandler = abort()
#elif ( defined OSF1 )
      stop
#else
      call abort()
#endif
      return
      end
