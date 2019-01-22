#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <sys/types.h>

#include "siemens_meas.h"

static int _verbose = 0;

typedef struct _mdh {
  u_int32_t DMALength;
  int32_t MeasUID;
  u_int32_t ScanCounter;
  u_int32_t TimeStamp;
  u_int32_t PMUTimeStamp;
  u_int32_t EvalInfoMask[2];
  u_int16_t SamplesInScan;
  u_int16_t UsedChannels;
  u_int16_t LoopCounter[14];
  u_int16_t Pre;
  u_int16_t Post;
  u_int16_t KSpaceCenterColumn;
  u_int16_t Dummy;
  float ReadOutOffCenter;
  u_int32_t TimeSinceLastRF;
  u_int16_t KSpaceCenterLineNum;
  u_int16_t KSpaceCenterPartitionNum;
  u_int16_t IceProgramParams[4];
  u_int16_t FreeParams[4];
  float SliceData[7];
  u_int16_t ChannelId;
  u_int16_t PTABPosNeg;
} mdh_t;


#define MDH_ACQEND      (0x00000001)
#define MDH_SCANDATA    (0x00000020)

typedef struct _mdh_scandata {
  u_int32_t datalen;
  char name[8];
  u_int32_t tmp[12]; // dunno??  must be something
} mdh_scandata_t;

#define PERROR(s) do {                                                  \
    char message[256];                                                  \
    snprintf(message, sizeof(message), "%s:%d %s", __FILE__, __LINE__, s); \
    perror(message);                                                    \
  } while (0)


static FILE * siemens_open_to_mdh(const char *filename)
{
  u_int32_t hlen;

  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "file'%s'\n", filename);
    PERROR("siemens_open_to_mdh");
    return NULL;
  }
  if (1 != fread(&hlen, sizeof(hlen), 1, fp)) {
    PERROR("fread");
    fclose(fp);
    return NULL;
  }
  if (fseek(fp, hlen, SEEK_SET)) {
    PERROR("fseek");
    fclose(fp);
    return NULL;
  }

  return fp;
}

int siemens_get_measinfo(const char *filename, measinfo_t * measinfo)
{
  int ii = 0;
  mdh_t mdh;

  FILE *fp = siemens_open_to_mdh(filename);
  if (!fp)
    return -1;

  measinfo->np = 0;
  measinfo->ntraces = 0;
  measinfo->nchan = 1;
  measinfo->max_aux_len = 0;
  measinfo->num_aux = 0;

  while (1) {
    if (1 != fread(&mdh, sizeof(mdh), 1, fp)) {
      break;
    }
    ii++;
    mdh.DMALength &= 0x1ffff;
    int seeklen = mdh.DMALength - sizeof(mdh);
    if ((mdh.EvalInfoMask[0] & MDH_SCANDATA)) {
      mdh_scandata_t scandata_header;
      if (1 != fread(&scandata_header, sizeof(scandata_header), 1, fp)) {
        PERROR("fread");
        fclose(fp);
        return -1;
      }
      measinfo->num_aux++;
      if (!measinfo->max_aux_len) {
        measinfo->max_aux_len = scandata_header.datalen;
      } else {
        if (measinfo->max_aux_len != scandata_header.datalen)
          printf("odd scandata len: %d != %d, %d @ %d EvalInfoMask:[0x%x 0x%x]\n",
                 measinfo->max_aux_len, scandata_header.datalen,
                 mdh.DMALength, ii, mdh.EvalInfoMask[1], mdh.EvalInfoMask[0]);
      }
      if (_verbose) {
        /*fprintf(stdout, "scand: 0x%p\n", (void *)scandata_header.datalen);*/
      }
      if (fseek(fp, -sizeof(scandata_header), SEEK_CUR)) {
        PERROR("fseek");
        fclose(fp);
        return -1;
      }
    }
    // valid readout
    else if (!(mdh.EvalInfoMask[1] == 0 && mdh.EvalInfoMask[0] == MDH_ACQEND)) {
      if (mdh.ScanCounter >= measinfo->ntraces)
        measinfo->ntraces = mdh.ScanCounter;
      if (mdh.ChannelId >= measinfo->nchan)
        measinfo->nchan = mdh.ChannelId + 1;
      if (!measinfo->np) {
        measinfo->np = mdh.SamplesInScan;
      } else {
        if (measinfo->np != mdh.SamplesInScan)
          printf("odd mdhlen: %d @ %d EvalInfoMask:[0x%x 0x%x]\n", mdh.DMALength,
                 ii, mdh.EvalInfoMask[1], mdh.EvalInfoMask[0]);
      }
      // in this case, DMALength lies!
      seeklen = mdh.SamplesInScan * 2 * sizeof(float);
      if (_verbose) {
      /*  fprintf(stdout, "reado: ch%d %p samples: %d cnt: %d  f:%ld\n",
                mdh.ChannelId, (void *)mdh.DMALength, mdh.SamplesInScan, 
                mdh.ScanCounter, (long)ftell(fp));*/
      }
    } else {
      fprintf(stdout, "unknown mdh %d, 0x%x : 0x%x, dmalen: %d  chan%d\n", ii,
              mdh.EvalInfoMask[0], mdh.EvalInfoMask[1], mdh.DMALength, mdh.ChannelId);
      //fprintf(stdout, "ftell: %ld\n", ftell(fp));
    }
    // move forward
    if (fseek(fp, seeklen, SEEK_CUR)) {
      PERROR("fseek");
      fclose(fp);
      return -1;
    }
  }
  //printf("read %d mdhs\n", ii);
  fclose(fp);
  return 0;
}

size_t siemens_get_mdh_scandata(const char *filename, void * data)
{
  size_t len = 0;
  FILE *fp = siemens_open_to_mdh(filename);
  mdh_t mdh;
  size_t nread = 0;

  if (!fp)
    return 0;

  while (1) {
    if (1 != fread(&mdh, sizeof(mdh), 1, fp)) {
      break;
    }
    mdh.DMALength &= 0x1ffff;
    if ((mdh.EvalInfoMask[0] & MDH_SCANDATA)) {
      mdh_scandata_t scandata_header;
      if (1 != fread(&scandata_header, sizeof(scandata_header), 1, fp)) {
        PERROR("fread");
        fclose(fp);
        return 0;
      }
      if (1 != fread(data + len, scandata_header.datalen, 1, fp)) {
        PERROR("fread");
        fclose(fp);
        return 0;
      }
      len += scandata_header.datalen;
      nread++;
      //printf("scandata len: %d\n", scandata_header.datalen);
      if (fseek(fp, mdh.DMALength - sizeof(mdh) - sizeof(scandata_header) - scandata_header.datalen, SEEK_CUR)) {
        PERROR("fseek");
        fclose(fp);
        return 0;
      }
    } else {
      if (fseek(fp, mdh.DMALength - sizeof(mdh), SEEK_CUR)) {
        PERROR("fseek");
        fclose(fp);
        return 0;
      }
    }
  }
  fclose(fp);
  return nread;
}

static float * BASE(measinfo_t * measinfo, void * data, int chanid)
{
  float * fdata = data;
  if (chanid >= measinfo->nchan) {
    fprintf(stdout, "badchan: %d\n", chanid);
    assert(0);
  }
  return &fdata[measinfo->np * measinfo->ntraces * 2 * chanid];
}

size_t siemens_get_mdh_traces(const char *filename, measinfo_t * measinfo, void * data, float * mdh_fields)
{
  size_t len[measinfo->nchan];
  FILE *fp = siemens_open_to_mdh(filename);
  mdh_t mdh;
  size_t nread = 0;

  if (!fp)
    return 0;

  for (int i = 0; i < measinfo->nchan; i++)
    len[i] = 0;

  while (1) {
    if (1 != fread(&mdh, sizeof(mdh), 1, fp)) {
      break;
    }
    mdh.DMALength &= 0x1ffff;
    // is a readout?
    if (!(mdh.EvalInfoMask[0] & MDH_SCANDATA) &&
        !(mdh.EvalInfoMask[1] == 0 && mdh.EvalInfoMask[0] == MDH_ACQEND)) {
      size_t tlen = mdh.SamplesInScan * 2 * sizeof(float);
      size_t chanid = mdh.ChannelId;
      float *dest = BASE(measinfo, data, chanid) + (len[chanid] / sizeof(float));
      if (1 != fread((void *)dest, tlen, 1, fp))
        return 0;
      // verify fp values
      if (0) {
        int kk;
        for (kk = 0; kk < mdh.SamplesInScan * 2; kk++)
          if (!isfinite(dest[kk]))
            fprintf(stderr, "bad fp val, ch%ld @ %ld\n", chanid, len[chanid]);
      }

      // advance
      len[chanid] += tlen;

      if (chanid == measinfo->nchan - 1) {
        if (mdh_fields) {
          mdh_fields[nread * NUM_MDH_FIELDS + 0] = mdh.LoopCounter[0];
          mdh_fields[nread * NUM_MDH_FIELDS + 1] = (float)mdh.ReadOutOffCenter;
          
          
        }
        nread++;
      }
    } else {
      if (_verbose)
        fprintf(stdout, "skipping fwd %ld  @ %ld\n",
                (long)mdh.DMALength - sizeof(mdh), (long)nread);
      if (fseek(fp, mdh.DMALength - sizeof(mdh), SEEK_CUR))
        return 0;
    }
  }
  for (int i = 0; i < measinfo->nchan; i++) {
    if (len[i] != len[measinfo->nchan - 1]) {
      fprintf(stderr, "chan data incomplete: len[nchan-1]:%ld len[%d]:%ld\n",
              len[measinfo->nchan-1], i, len[i]);
    }
  }
  fclose(fp);
  return nread;
}

#ifdef TEST

#include <math.h>
#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif
#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif


/*! \brief      Calculate several statistics (min,max,avg,stddev) on an array.
 *
 * \param[in]           n       Array length.
 * \param[in]           a       Pointer to real_t array to inspect.
 * \param[out]          min     Gets min value.
 * \param[out]          max     Gets max value.
 * \param[out]          avg     Gets avg value.
 * \param[out]          stddev  Gets stddev.
 */
typedef float real_t;
void array1r_stats(size_t n, const real_t *a, real_t *min, real_t *max, real_t *avg, real_t *stddev)
{
  size_t i;

  double _min, _max, _avg, _stddev;
  _avg = 0;
  _min = _max = a[0];
  for (i = 0; i < n; i++) {
    _min = MIN(_min, a[i]);
    _max = MAX(_max, a[i]);
    _avg += a[i];
  }
  _avg /= n;
  *min = _min;
  *max = _max;
  *avg = _avg;

  // standard deviation
  _stddev = 0.;
  for (i = 0; i < n; i++) {
    _stddev += (a[i] - _avg) * (a[i] - _avg);
  }
  *stddev = sqrt(_stddev / n);
}

static int usage(char *argv[])
{
  printf("usage: %s <.meas file>\n", argv[0]);
  return -1;
}

int main(int argc, char *argv[])
{
  if (argc < 2) {
    usage(argv);
    exit(-1);
  }

  int ff;
  for (ff = 1; ff < argc; ff++) {
    char * measfile = argv[ff];
    measinfo_t measinfo;

    if (siemens_get_measinfo(measfile, &measinfo)) {
      //printf("unable to read: %s\n", measfile);
      continue;
    }
    printf("%s\n", argv[ff]);
    printf("measinfo: np          %d\n"
           "measinfo: ntraces     %d\n"
           "measinfo: nchan       %d\n"
           "measinfo: num_aux     %d\n"
           "measinfo: max_aux_len %d\n",
           measinfo.np, measinfo.ntraces, measinfo.nchan,
           measinfo.num_aux, measinfo.max_aux_len );

    void * scandata = calloc(measinfo.num_aux * measinfo.max_aux_len, 1);
    assert(scandata);
    siemens_get_mdh_scandata(measfile, scandata);
    free(scandata);

    void * tracedata = calloc(measinfo.np * measinfo.ntraces * measinfo.nchan, 2 * sizeof(float));
    assert(tracedata);
    siemens_get_mdh_traces(measfile, &measinfo, tracedata, NULL);

    float * traces = tracedata;

    for (int cc = 0; cc < measinfo.nchan; cc++) {
      float min, max, avg, stddev;
      float * traces2 = BASE(&measinfo, traces, cc);
      array1r_stats(measinfo.np * measinfo.ntraces * 2,
                    traces2, &min, &max, &avg, &stddev);
      float re = 0, im = 0;
      for (size_t i = 0; i < measinfo.np * measinfo.ntraces; i++) {
        re += traces2[i * 2];
        im += traces2[i * 2 + 1];
      }
      printf("ch%d: min %f  max %f  avg %f  stddev %f  |arg| %f\n", 
             cc, min, max, avg, stddev, atan2(re, im));
    }
    free(tracedata);
  }
  return 0;
}
#endif
