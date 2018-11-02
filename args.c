#include "util.h"
#include "args.h"

#ifndef __lint
static const char rcsid[] =
"$Id: args.c,v 4.0 2004/08/26 09:49:25 mmichael Exp $";
#endif

static int argc;
static char **argv;
char *argv0;

/* ckargs  --  check that only certain parameters are set on the command line */
void ckargs(const char *options, int argcx, char **argvx, int non_options)
{
	int i;

	argc = argcx;
	argv = argvx;
	argv0 = argv0 ? argv0 : argv[0];
	for (i = non_options+1; i < argc; ++i)
	  /*
		if (argv[i][1] != '=')
			fatalf("Improper command option: '%s'.", argv[i]);
		else if (!strchr(options, argv[i][0]))
			fatalf("Available options: %s\n", options);
	  */
	        if (argv[i][1] == '=' && !strchr(options, argv[i][0]))
			fatalf("Available options: %s\n", options);
}

/* get_argval  --------------------- get the value of a command-line argument */
bool get_argval(int c, int *val_ptr)
{
	int i;

	ck_argc("get_argval");
	for (i = 0; i < argc; ++i)
		if (argv[i][0] == c && argv[i][1] == '=') {
			*val_ptr = atoi(argv[i]+2);
			return 1;
		}
	return 0;
}

/* get_fargval  -------------- get the float value of a command-line argument */
bool get_fargval(int c, double *val_ptr)
{
        int i;

        ck_argc("get_fargval");
        for (i = 0; i < argc; ++i)
                if (argv[i][0] == c && argv[i][1] == '=') {
                        *val_ptr = atof(argv[i]+2);
                        return 1;
                }
        return 0;
}
bool get_long_fargval(char *id, double *val_ptr){
  int i;

  ck_argc("get_long_fargval");
  for (i = 0; i < argc; ++i)
    if (!strncasecmp(argv[i], id, strlen(id))&& argv[i][strlen(id)] == '='){
      *val_ptr = atof(argv[i]+strlen(id)+1);
      return 1;
    }
  return 0;
}


bool get_cargval(int c, char **valp)
{
        int i;

        ck_argc("get_cargval");
        for (i = 0; i < argc; ++i)
                if (argv[i][0] == c && argv[i][1] == '=') {
                        *valp = argv[i]+2;
                        return 1;
                }
        return 0;
}

/* ck_argc - die if argc is unknown */
void ck_argc(const char *proc_name)
{
	if (argc == 0)
		fatalf("Call ckargs() before %s.\n", proc_name);
}

void fprintf_argv(FILE* fp)
{
	int i;
	fprintf(fp, "%s", argv0);
	for (i = 1; i < argc; ++i)
		(void)fprintf(fp, " %s", argv[i]);
}


void get_argval_min(int c, int *v, int d, int min, const char *msg)
{
	if (get_argval(c, v)) {
		if (*v < min)
			fatalf(msg, c);
	} else {
		*v = d;
	}
}	

void get_argval_max(int c, int *v, int d, int max, const char *msg)
{
	if (get_argval(c, v)) {
		if (*v > max)
			fatalf(msg, c);
	} else {
		*v = d;
	}
}

void get_argval_nonneg(int ch, int *val, int dflt)
{
	get_argval_min(ch, val, dflt, 0, "%c must be non-negative.");
}

void get_fargval_nonneg(int ch, double *val, double dflt)
{
  if (get_fargval(ch, val)) {
    if (*val < 0.)
      fatalf("%c must be non-negative.", ch);
  } else {
    *val = dflt;
  }
}

void get_fargval_limit(int ch, double lower, double higher, double *val, double dflt)
{
  if (get_fargval(ch, val)) {
    if (*val <= lower||*val >= higher)
	fatalf("%c must be in ]%f,%f[.",ch,lower,higher);
  } else {
    *val = dflt;
  }
}

void get_long_fargval_limit(char *id, double lower, double higher, double *val, double dflt){
   if (get_long_fargval(id, val)) {
    if (*val <= lower||*val >= higher)
      fatalf("%s must be in ]%f,%f[.",id,lower,higher);
  } else {
    *val = dflt;
  }
} 

void get_argval_pos(int ch, int *val, int dflt)
{
	get_argval_min(ch, val, dflt, 1, "%c must be positive.");
}

void get_argfile(int ch, FILE **fp, FILE *dflt){
  int i;

  ck_argc("get_argfile");
  for(i=0;i<argc;i++)
    if(argv[i][0]==ch && argv[i][1] == '=') {
      *fp=ckopen(argv[i]+2,"r");
      return;
    }
  *fp=dflt;
}

char* get_argfilename(int ch){
  int i;
  ck_argc("get_argfilename");  
  for(i=0;i<argc;i++)
    if(argv[i][0]==ch && argv[i][1] == '=')
      return argv[i]+2;
  return NULL;
}

