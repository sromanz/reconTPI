#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <omp.h> 
#include "siemens_meas.h"

/* [traces dirs] = mexLoadSiemensTraces(filename) */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char * errmsg = 
    "usage: [traces dirs] = mexLoadSiemensTraces(filename)\n"
    "       [traces]      = mexLoadSiemensTraces(filename)\n";

  if (nrhs != 1 || 
      (nlhs != 1 && nlhs != 2) ||
      !mxIsChar(prhs[0]))
    mexErrMsgTxt(errmsg);

  measinfo_t measinfo;
  char measfile[512];
  if (mxGetString(prhs[0], measfile, sizeof(measfile)))
    mexErrMsgTxt(errmsg);

  if (siemens_get_measinfo(measfile, &measinfo)) {
    mexErrMsgTxt("cant get measinfo\n");
  }

  if (measinfo.np * measinfo.ntraces == 0) {
    mexErrMsgTxt("zero-sized data\n");
  }

  /* Construct outdata */
  mwSize outdims[3];
  outdims[0] = measinfo.np;
  outdims[1] = measinfo.ntraces;
  outdims[2] = measinfo.nchan;

  plhs[0] = mxCreateNumericArray(3, outdims, mxSINGLE_CLASS, mxCOMPLEX);
  float *outdataR = mxGetData(plhs[0]);
  float *outdataI = mxGetImagData(plhs[0]);

  float *outdata = calloc(measinfo.np * measinfo.ntraces * measinfo.nchan, 2 * sizeof(float));
  if (!outdata)
    mexErrMsgTxt("mexLoadSiemensTraces() unable to calloc\n");

  float *dirs = NULL;
  if (nlhs == 2) {
    outdims[0] = 2;
    outdims[1] = measinfo.ntraces;
    plhs[1] = mxCreateNumericArray(2, outdims, mxSINGLE_CLASS, mxREAL);
    dirs = mxGetData(plhs[1]);
  }

  /* load scandata */
  if (!siemens_get_mdh_traces(measfile, &measinfo, outdata, dirs))
    mexErrMsgTxt("error loading mdh traces\n");

  // uninterleave
  {
    size_t j;
    for (j = 0; j < (measinfo.np * measinfo.ntraces * measinfo.nchan); j++) {
      outdataR[j] = outdata[j*2];
      outdataI[j] = outdata[j*2+1];
    }
  }
  free(outdata);
}
