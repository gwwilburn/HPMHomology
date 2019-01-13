/* hpm_trace.h */

#include "hmmer.h"
#include "p7_config.h"

#include "easel.h"
#include "esl_vectorops.h"


#include "hpm_trace.h"

int
hpm_trace_Copy(const P7_TRACE *src, P7_TRACE *dst)
{

	return eslOK;
}


P7_TRACE *
hpm_trace_Clone(const P7_TRACE *tr)
{
	P7_TRACE *tr2      = p7_trace_Create();
	int       status;

	if ((status = hpm_trace_Copy(tr, tr2)) != eslOK) goto ERROR;

	return tr2;

	ERROR:
		p7_trace_Destroy(tr2);
		return NULL;
}
