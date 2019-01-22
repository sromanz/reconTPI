#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <omp.h> 
#include "siemens_meas.h"

/* [traces] = mexLoadSiemensScandata(filename) */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char * errmsg = "usage: [traces] = mexLoadSiemensScandata(filename)";

  if (nrhs != 1 || nlhs != 1 || !mxIsChar(prhs[0]))
    mexErrMsgTxt(errmsg);

  measinfo_t measinfo;
  char measfile[512];
  if (mxGetString(prhs[0], measfile, sizeof(measfile)))
    mexErrMsgTxt(errmsg);

  if (siemens_get_measinfo(measfile, &measinfo)) {
    mexErrMsgTxt("cant get measinfo\n");
  }

  /* Construct outdata */
  mwSize outdims[2];
  outdims[0] = measinfo.max_aux_len / sizeof(float);
  outdims[1] = measinfo.num_aux;

  if (outdims[0] * outdims[1] == 0) {
    mexErrMsgTxt("zero-sized data\n");
  }

  //printf("size %d %d\n", outdims[0], outdims[1]);

  plhs[0] = mxCreateNumericArray(2, outdims, mxSINGLE_CLASS, mxREAL);
  float *outdata = mxGetData(plhs[0]);

  /* load scandata */
  if (!outdata || !siemens_get_mdh_scandata(measfile, outdata))
    mexErrMsgTxt(errmsg);

}
