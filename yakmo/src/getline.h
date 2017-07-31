#ifndef _getline_h_
#define _getline_h_ 1

#define __STDC__ 1

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

#if defined (__GNUC__) || (defined (__STDC__) && __STDC__)
#define __PROTO(args) args
#else
#define __PROTO(args) ()
#endif  /* GCC.  */

#define GETLINE_NO_LIMIT -1

int
  getline __PROTO ((char **_lineptr, size_t *_n, FILE *_stream));
int
  getline_safe __PROTO ((char **_lineptr, size_t *_n, FILE *_stream,
                         int limit));
int
  getstr __PROTO ((char **_lineptr, size_t *_n, FILE *_stream,
		   int _terminator, int _offset, int limit));

#ifdef __cplusplus
}
#endif

#endif /* _getline_h_ */
