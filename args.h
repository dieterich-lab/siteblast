bool get_argval(int, int *);
bool get_fargval(int, double *);
bool get_long_fargval(char *, double *);
bool get_cargval(int, char **);
void ckargs(const char *, int , char **, int );
void fprintf_argv(FILE* );
void ck_argc(const char *);

void get_argval_min(int c, int *v, int d, int min, const char *msg);
void get_argval_max(int c, int *v, int d, int max, const char *msg);
void get_argval_nonneg(int ch, int *val, int dflt);
void get_fargval_nonneg(int ch, double *val, double dflt);
void get_fargval_limit(int ch, double lower, double higher, double *val, double dflt);
void get_long_fargval_limit(char *id, double lower, double higher, double *val, double dflt);
void get_argval_pos(int ch, int *val, int dflt);
void get_argfile(int ch, FILE **fp, FILE *dflt);
char* get_argfilename(int ch);

extern char *argv0;
