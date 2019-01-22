#ifdef __cplusplus
extern "C" {
#endif

typedef struct _measinfo {
  u_int32_t np;                 //! 
  u_int32_t ntraces;            //! 
  u_int32_t nchan;              //! number of channels
  u_int32_t max_aux_len;        //! max len of scandata
  u_int32_t num_aux;            //! number of scandatas
} measinfo_t;

#define NUM_MDH_FIELDS  (2)     //! currently x,y,z

int siemens_get_measinfo(const char *filename, measinfo_t * measinfo);
size_t siemens_get_mdh_scandata(const char *filename, void * data);
size_t siemens_get_mdh_traces(const char *filename, measinfo_t * measinfo, void * data, float * mdh_fields);

#ifdef __cplusplus
}
#endif
