/*****************************************************************************/
/* pexec.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Simple program for shell command execution in parallel		     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Version: PEXEC_VERSION (see pexec.h), Last modification: see also pexec.h */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Copyright (C) 2007, 2008-2009; Pal, A. (apal@szofi.elte.hu)		     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Development of this program has been supported by the HATNet project	     */
/* (see also http://hatnet.hu).						     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*  This program is free software: you can redistribute it and/or modify     */
/*  it under the terms of the GNU General Public License as published by     */
/*  the Free Software Foundation, either version 3 of the License, or	     */
/*  (at your option) any later version.					     */
/*									     */
/*  This program is distributed in the hope that it will be useful,	     */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of	     */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	     */
/*  GNU General Public License for more details.			     */
/*									     */
/*  You should have received a copy of the GNU General Public License	     */
/*  along with this program.  If not, see <http://www.gnu.org/licenses/>.    */
/*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdarg.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/wait.h>
#include <sys/select.h>
#include <sys/resource.h>
#include <sys/socket.h>
#include <fcntl.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <netdb.h>
#include <signal.h>
#include <ctype.h>
#include <errno.h>
#ifdef HAVE_SYS_LOADAVG_H
#include <sys/loadavg.h>
#endif
#ifdef HAVE_LIBINTL_H
#include <libintl.h>
#define	_(string)	gettext(string)
#else
#define	_(string)	string
#endif

#include "tokenize.h"
#include "iof.h"
#include "numhash.h"
#include "fifo.h"
#include "format.h"
#include "str.h"
#include "linebuffer.h"
#include "longhelp.h"

#define	__LIST_MALLOC	xmalloc

#include "list.h"

#include "pexec.h"

#define	PEXEC_ABORT_ON_MALLOC
#ifdef	PEXEC_ABORT_ON_MALLOC
#define	realloc_check(ptr,size) \
	do \
	 {	if ( ptr==NULL && size>0 ) \
		 {	fprintf(stderr,"pexec.c: %s.\n",_("memory exhausted"));\
			abort(); \
		 } \
	 } while(0)
#else
#define	realloc_check(ptr,size) 
#endif
#define	malloc_check(ptr)	realloc_check(ptr,1)

static void *xmalloc(size_t size)
{
 void	*ret;
 ret=malloc(size);
 malloc_check(ret); 
 return(ret);
}

static void *xrealloc(void *ptr,size_t size)
{
 void	*ret;
 ret=realloc(ptr,size);
 realloc_check(ret,size); 
 return(ret);
}
static char *xstrdup(const char *s)
{
 char	*ret;
 ret=strdup(s);
 malloc_check(ret);
 return(ret);
}

/*****************************************************************************/

/* global variables used throughout this program: */

char	*progbasename;	/* basename of the program, normally it is 'pexec'   */
int	sig_pipe[2];	/* file descriptors of pipes used in signal handlers */

/*****************************************************************************/

#define		PEXEC_ISTR_AUTO		"auto"
#define		PEXEC_ISTR_MANAGED	"managed"
#define		PEXEC_ISTR_NCPU		"ncpu"

/*****************************************************************************/

char *logmsg_submit_task[] = 
 {	NULL,
	"unable to open " PEXEC_DEFAULT_NULLFILE " special file for reading",
	"unable to redirect from input file: file cannot be opened",
	"unable to open " PEXEC_DEFAULT_NULLFILE " special file for writing",
	"unable to create internal pipes",
	"unable to redirect to output file: file cannot be created",
	"unable to redirect to error file: file cannot be created",
	"unknown error",
	NULL
 };

int log_message(logdata *log,int loglevel,parameter *p,char *msg,...)
{
 va_list	ap;

 if ( log->fwlog==NULL || log->loglevel<loglevel || msg==NULL )
	return(0);

 if ( p != NULL )
	fprintf(log->fwlog,"[%s] ",p->name);

 va_start(ap,msg);
 vfprintf(log->fwlog,msg,ap);
 va_end(ap);

 fprintf(log->fwlog,"\n");
 
 fflush(log->fwlog);

 return(0);
}

/*****************************************************************************/

/* an implementation of the daemon() function for non-linux/GNU/BSD systems  */
int background(int nochdir,int noclose)
{
 int	child_pid;	/* PID of child process */
 int	nulldev_fd;	/* file descriptor for null device */

 child_pid = fork();

 if ( child_pid < 0 )	return(-1);

 else if ( child_pid == 0 )
  {
	#ifdef TIOCNOTTY
	if ( ( i=open("/dev/tty",O_RDWR) ) >= 0 )
	 {	ioctl(i, TIOCNOTTY, 0) ;
		close(i);
	 }
	#endif

	if ( ! nochdir )
		chdir("/");

	if ( ! noclose )
	 {	nulldev_fd=open("/dev/null",O_RDONLY);
		if ( nulldev_fd>=0 )
		 {	if ( isatty(0) )
			 {	close(0);
				dup2(nulldev_fd,0);
			 }
			close(nulldev_fd);
		 }
		nulldev_fd=open("/dev/null",O_WRONLY);
		if ( nulldev_fd>=0 )
		 {	if ( isatty(1) )
			 {	close(1);
				dup2(nulldev_fd,1);
			 }
			if ( isatty(2) )
			 {	close(2);
				dup2(nulldev_fd,2);
			 }
			close(nulldev_fd);
		 }
	 }
  }
 else
 	exit(0);

 return(0);
}

/*****************************************************************************/

void sig_act_child(int signum)
{
 int		pid,ret,status;
 signalinfo	sci;

 if ( sig_pipe[1]<0 )
	return;

 while ( (pid=waitpid(-1,&status,WNOHANG))>0 )
  {	sci.signal=SIGCHLD;
	sci.pid=pid;
	sci.exitstatus=WEXITSTATUS(status);
	if ( WIFSIGNALED(status) )
		sci.exitsignal=WTERMSIG(status);
	else
		sci.exitsignal=-1;

	ret=write(sig_pipe[1],&sci,sizeof(signalinfo));
  }

}

void sig_act_interrupt(int signum)
{
 signalinfo	sci;
 int		ret;

 if ( sig_pipe[1]<0 )
	return;

 sci.signal=signum;
 sci.pid=0;
 sci.exitstatus=0;
 sci.exitsignal=!0;

 ret=write(sig_pipe[1],&sci,sizeof(signalinfo));
}

/*****************************************************************************/

int fdwait(int fd,int timeout)
{
 fd_set		set;
 struct	timeval	tv;
 int		ret;

 FD_ZERO(&set);
 FD_SET(fd,&set);

 if ( timeout>0 )
  {	tv.tv_sec=timeout;
	tv.tv_usec=0;
	ret=select(fd+1,&set,NULL,NULL,&tv);
  }
 else
	ret=select(fd+1,&set,NULL,NULL,NULL);

 return(ret);
}

/*****************************************************************************/

#ifdef	HAVE_SETENV

int env_export(char *name,char *value)
{
 int	ret;
 ret=setenv(name,value,!0);
 return(ret); 
}

#else

int env_export(char *name,char *value)
{
 int	ret;
 char	*envstr;

 envstr=NULL;
 strappendf(&envstr,"%s=%s",name,value);
 ret=putenv(envstr); 
/* free(envstr); */

 return(ret); 
}

#endif

/*****************************************************************************/

int is_nasty_char(int c)
{
 if ( isalnum(c) )
	return(0);
 else if ( c=='-' || c=='+' || c=='@' || c=='^' )
	return(0);
 else if ( c=='_' || c=='=' || c==':' || c=='/' )
	return(0);
 else if ( c==',' || c=='.' )
	return(0);
 else
	return(1);
}

int get_nasty_char_number(char *arg)
{
 int	n;
 for ( n=0 ; *arg ; arg++ )
  {	if ( is_nasty_char(*arg) )
		n++;
  }
 return(n);
}

/*
char *make_escaped_shell_command(char **argv,int argc)
{
 char	*command;
 int	len,i,j,x,l;

 len=0;
 command=NULL;

 for ( i=0 ; i<argc ; i++ )
  {	x=get_nasty_char_number(argv[i]);
	l=strlen(argv[i]);
	if ( command==NULL )
		command=(char *)xmalloc(l+x+1);
	else
	 {	command=(char *)xrealloc(command,len+1+l+x+1);
		strcpy(command+len," ");
		len++;
	 }
	for ( j=0 ; argv[i][j] ; j++ )
	 {	command[len]=argv[i][j];
		len++;
	 }
  }
 command[len]=0;
	
 return(command);
}
*/

char *concatenate_arguments(int argc,char **argv)
{
 char	*command;
 int	len,i,l;

 len=0;
 command=NULL;

 for ( i=0 ; i<argc ; i++ )
  {	l=strlen(argv[i]);
	if ( command==NULL )
		command=(char *)xmalloc(l+1);
	else
	 {	command=(char *)xrealloc(command,len+1+l+1);
		strcpy(command+len," ");
		len++;
	 }
	strcpy(command+len,argv[i]);
	len+=l;
  }
 if ( command==NULL )
	command=(char *)xmalloc(16);

 command[len]=0;
	
 return(command);
}

/*****************************************************************************/

static int daemon_commandtoken_is_nasty(int c)
{
 if ( c<=32 || c>=127 || c=='\"' || c=='\'' || c=='\\' )
	return(1);
 else
	return(0);
}

static int hex_digit(int c)
{
 if ( '0' <= c && c <= '9' )
	return(c-'0');
 else if ( 'a' <= c && c <= 'f' )
	return(c-'a'+10);
 else if ( 'A' <= c && c <= 'F' )
	return(c-'A'+10);
 else
	return(-1);
}

char * daemon_commandtoken_escape(char *buff,int size)
{
 int	l,s,c;
 char	*p,*q,*out;

 if ( buff==NULL )
	return(NULL);
 
 l=0;
 for ( p=buff,s=size ; s>0 ; p++,s-- )
  {	if ( daemon_commandtoken_is_nasty((int)(*p)) )
		l+=3;
	else
		l++;
  }
 out=(char *)xmalloc(l+1);
 for ( p=buff,q=out,s=size ; s>0 ; p++,s-- )
  {	if ( daemon_commandtoken_is_nasty((int)(*p)) )
	 {	c=(*(unsigned char *)p);
		sprintf(q,"\\%.2X",c);
		q+=3;
	 }
	else
	 {	*q=*p;
		q++;
	 }
  }
 *q=0;

 return(out);
}

char * daemon_commandtoken_escape_string(char *buff)
{
 int	len;
 if ( buff==NULL || (len=strlen(buff))<=0 )
	return(NULL);
 else
	return(daemon_commandtoken_escape(buff,len));
}

int daemon_commandtoken_unescape(char *buff)
{
 char	*out;
 int	h1,h2,ret;
 
 ret=0;
 for ( out=buff ; *buff ; )
  {	if ( *buff=='\\' && 
	(h1=hex_digit((int)(*(buff+1))))>=0 && 
	(h2=hex_digit((int)(*(buff+2))))>=0 )
	 {	*out=h1*16+h2;
		buff+=3;
		out++;
	 }
	else
	 {	*out=*buff;
		buff++;
		out++;
	 }
	ret++;
  }

 *out=0;

 return(ret);
}

/*****************************************************************************/

#define	HPRINT_BUFFER_SIZE		256

int hprintf(int handle,char *msg,...)
{
 char		buff[HPRINT_BUFFER_SIZE],*tbuff;
 va_list	ap;
 int		n;

 va_start(ap,msg);
 n=vsnprintf(buff,HPRINT_BUFFER_SIZE,msg,ap);
 va_end(ap);
 if ( n<HPRINT_BUFFER_SIZE )
  {	write(handle,buff,n);
	return(n);
  }
 else
  {	tbuff=NULL;
	va_start(ap,msg);
	vstrappendf(&tbuff,msg,ap);
	va_end(ap);
	if ( tbuff != NULL )
	 {	n=strlen(tbuff);
		write(handle,tbuff,n);
		free(tbuff);
	 }
	else
		n=0;

	return(n);
  }

}

/*****************************************************************************/

#define		BUFF_LEN	1024

int get_number_of_cpus_proccpuinfo(void)
{
 FILE	*fr;
 char	buff[BUFF_LEN];
 int	ncpu,i,j;

 fr=fopen("/proc/cpuinfo","rb");
 if ( fr != NULL )
  {	ncpu=0;
	while ( ! feof(fr) ) 
	 {	if ( fgets(buff,BUFF_LEN-1,fr)==NULL )  break;
		for ( i=0 ; i<strlen(buff) ; )
		 {      if ( buff[i]==10 || buff[i]==9 || 
			buff[i]==13 || buff[i]==32 )
				memmove(buff+i,buff+i+1,strlen(buff+i));
			else
				i++;
		 }
		for ( j=0 ; j<strlen(buff) ; j++ )
		 {	if ( buff[j]==':' )	break;		}
		if ( j==strlen(buff) )  continue;
		buff[j]=0;
		if ( strcmp(buff,"processor")==0 )
			ncpu++;
	 };
	fclose(fr);
  }
 else
	ncpu=0;

 return(ncpu);
}

int get_number_of_cpus_sysconf(void)
{
 int	ncpu;
#ifdef _SC_NPROCESSORS_ONLN
 ncpu=sysconf(_SC_NPROCESSORS_ONLN);
#else
 ncpu=0;
#endif
 return(ncpu);
}

int get_number_of_cpus(void)
{
 int	ncpu,w;

 /* try /proc/cpuinfo: (assumed to work on linux kernels) */
 ncpu=get_number_of_cpus_proccpuinfo();
 /* try sysconf(): */
 w=get_number_of_cpus_sysconf();
 if ( w>0 && ncpu>0 )	ncpu=(w<ncpu?w:ncpu);
 
 /* otherwise, we have at least one: */
 if ( ncpu<=0 )	ncpu=1;

 return(ncpu);
}

#undef		BUFF_LEN

/*****************************************************************************/

int get_bit_size(int n)
{
 int	r;
 for ( r=0 ; n>0 ; n/=2 )
  {	r++;			}
 return(r);
}

/*****************************************************************************/

child *	get_child_by_pid(child *cc,int pid)
{
 for ( ; cc != NULL ; cc=cc->next )
  {	if ( cc->pid==pid )
		return(cc);
  }
 return(NULL);
}

child *	get_child_by_id(child *cc,int id)
{
 for ( ; cc != NULL ; cc=cc->next )
  {	if ( cc->id==id )
		return(cc);
  }
 return(NULL);
}

imutex * get_imutex_by_name(imutex *mx,char *name)
{
 for ( ; mx != NULL ; mx=mx->next )
  {	if ( strcmp(mx->name,name)==0 )
		return(mx);
  }
 return(NULL);
}

int fd_avail(int fd)
{
 fd_set		set;
 struct	timeval	tv;

 FD_ZERO(&set);
 FD_SET(fd,&set);
 tv.tv_sec=0;
 tv.tv_usec=0;
 select(fd+1,&set,NULL,NULL,&tv);

 if ( FD_ISSET(fd,&set) )
	return(1); 
 else
	return(0);
}

char *fileformat_replace(char *format,parameter *par)
{
 char *name;

 name=format_replace(format,0,
	's',FORMAT_STRING,par->name,
	'k',FORMAT_STRING,par->id,
	'd',FORMAT_STRING,par->id+1,
	0);

 return(name);
}

int submit_task(paralleldata *p,parameter *par,child *c,int no_format_replace,
	parallelstatus *ps)
{
 int	stdfd[3];
 int	pipeout[2],pipeerr[2];
 int	pid;

 if ( par->no_touch_std )
  {	stdfd[0]=-1;
	stdfd[1]=-1;
	stdfd[2]=-1;
	pipeout[0]=-1;
	pipeout[1]=-1;
	pipeerr[0]=-1;
	pipeerr[1]=-1;
  }

 else
  {	if ( p->in == NULL )
	 {	stdfd[0]=open(PEXEC_DEFAULT_NULLFILE,O_RDONLY);
		if ( stdfd[0]<0 )
			return(-1);
	 }
	else
	 {	char	*inname;
		if ( no_format_replace )
	 		stdfd[0]=open(p->in,O_RDONLY);
		else
		 {	inname=fileformat_replace(p->in,par);
			stdfd[0]=open(inname,O_RDONLY);
			free(inname);
		 }
		if ( stdfd[0]<0 )
			return(-2);
	 }
 
	if ( p->out == NULL )
	 {	stdfd[1]=open(PEXEC_DEFAULT_NULLFILE,O_WRONLY);
		if ( stdfd[1]<0 )
			return(-3);
		pipeout[0]=-1;
		pipeout[1]=-1;
	 }
	else if ( p->fwout != NULL )
	 {	if ( pipe(pipeout) )
			return(-4);
		stdfd[1]=pipeout[1];
	 }
	else
	 {	char	*outname;
		if ( no_format_replace )
			stdfd[1]=open(p->out,O_WRONLY|O_CREAT|O_TRUNC,0644);
		else
		 {	outname=fileformat_replace(p->out,par);
			stdfd[1]=open(outname,O_WRONLY|O_CREAT|O_TRUNC,0644);
			free(outname);
		 }
		if ( stdfd[1]<0 )
			return(-5);
		pipeout[0]=-1;
		pipeout[1]=-1;
	 }

	if ( p->err == NULL )
	 {	stdfd[2]=open(PEXEC_DEFAULT_NULLFILE,O_WRONLY);
		if ( stdfd[2]<0 )
			return(-3);
		pipeerr[0]=-1;
		pipeerr[1]=-1;
	 }
	else if ( p->fwerr != NULL )
	 {	if ( pipe(pipeerr) )
			return(-4);
		stdfd[2]=pipeerr[1];
	 }
	else
	 {	char	*errname;
		if ( no_format_replace )
			stdfd[2]=open(p->err,O_WRONLY|O_CREAT|O_TRUNC,0644);
		else
		 {	errname=fileformat_replace(p->err,par);
			stdfd[2]=open(errname,O_WRONLY|O_CREAT|O_TRUNC,0644);
			free(errname);
		 }
		if ( stdfd[2]<0 )
			return(-6);
		pipeerr[0]=-1;
		pipeerr[1]=-1;
	 }

  }

 pid=fork();

 if ( pid<0 )
	return(-4);

 else if ( pid>0 )
  {	if ( pipeout[0]>=0 )
	 {	c->fdstdout=pipeout[0];
		close(pipeout[1]);
	 }
	else
	 {	c->fdstdout=-1;
		if ( stdfd[1]>=0 )	close(stdfd[1]);
	 }
	if ( pipeerr[0]>=0 )
	 {	c->fdstderr=pipeerr[0];
		close(pipeerr[1]);
	 }
	else
	 {	c->fdstderr=-1;
		if ( stdfd[2]>=0 )	close(stdfd[2]);
	 }
	if ( stdfd[0]>=0 )	close(stdfd[0]);
	return(pid);
  }
 else
  {	if ( pipeout[0]>=0 )	close(pipeout[0]);
	if ( pipeerr[0]>=0 )	close(pipeerr[0]);

	/* is ps (parallelstatus) is defined: close opened filedesc's */
	if ( ps != NULL )
	 {	client	*cl;
		child	*cc;
		for ( cl=ps->clientlist ; cl != NULL ; cl=cl->next )
		 {	if ( cl->peer>=0 )	close(cl->peer);
		 }
		if ( ps->sock >= 0 )
			close(ps->sock);
		if ( ps->hsck >= 0 )
			close(ps->hsck);
		for ( cc=ps->childlist ; cc != NULL ; cc=cc->next )
	 	 {	if ( cc->fdstdout>=0 )	close(cc->fdstdout);
			if ( cc->fdstderr>=0 )	close(cc->fdstderr);
		 }
	 }

	if ( stdfd[0]>=0 )
	 {	close(0);
		dup2(stdfd[0],0);
		close(stdfd[0]);
	 }

	if ( stdfd[1]>=0 )
	 {	close(1);
		dup2(stdfd[1],1);
		close(stdfd[1]);
	 }

	if ( stdfd[2]>=0 )
	 {	close(2);
		dup2(stdfd[2],2);
		close(stdfd[2]);
	 }

	close(sig_pipe[0]);
	close(sig_pipe[1]);
	
	if ( p->envvarname != NULL && par->name != NULL )
		env_export(p->envvarname,par->name);

	/* this is the normal shell execution: */
	if ( par->c.is_shell > 0 )
	 {	char	*argv[4];

		argv[0]=p->shell;	/* /bin/sh	*/
		argv[1]="-c";
		argv[2]=concatenate_arguments(par->c.argc,par->c.argv);
		argv[3]=NULL;
		
		execv(p->shell,argv);

		fprintf(stderr,_("%s: unable to execute the shell '%s'.\n"),
		progbasename,p->shell);

		exit(2);
	 }
	/* this is the remote shell execution: */
	else if ( par->c.is_shell < 0 )
	 {	char	**argv;
		int	i,argc,n;

		n=0;
		for ( n=0 ; p->rshargs[n] != NULL ; )	n++;
		argc=n+1+par->c.argc+1; /* rsh <arg> <host> <command> NULL */
		argv=(char **)xmalloc(sizeof(char *)*(argc+1));
		argc=0;
		for ( i=0 ; i<n ; i++ )
		 {	argv[argc]=p->rshargs[i];
			argc++;
	 	 }
		if ( par->name != NULL )
		 {	argv[argc]=par->name;
			argc++;
		 }
		else
		 {	argv[argc]="localhost";
			argc++;
		 }
		for ( i=0 ; i<par->c.argc ; i++ )
		 {	argv[argc]=par->c.argv[i];
			argc++;
		 }
		argv[argc]=NULL;

		execvp(argv[0],argv);

		fprintf(stderr,_("%s: unable to execute the remote shell command '%s'.\n"),
		progbasename,argv[0]);

		exit(2);
			
	 }
	/* this is the normal command execution: */
	else
	 {	char	**argv;
		int	i,argc;

		argc=par->c.argc;
		argv=(char **)xmalloc(sizeof(char *)*(argc+1));
		for ( i=0 ; i<argc ; i++ )
		 {	argv[i]=par->c.argv[i];		}
		argv[argc]=NULL;

		execvp(argv[0],argv);

		fprintf(stderr,_("%s: unable to execute the command '%s'.\n"),
		progbasename,argv[0]);

		exit(2);
	 }	
  }

 return(-7);	/* unreachable point, however */
}

int send_task(paralleldata *p,parameter *par,child *cc)
{
 char	**args,*token;
 int	a,i;

 args=(char **)xmalloc(sizeof(char *)*16);
 for ( a=0 ; a<16 ; a++ )
  {	args[a]=NULL;		}

 a=0;
 strappendf(&args[a],"execute");a++;
 strappendf(&args[a],"identifier=%d",par->id);a++;
 if ( p->envvarname != NULL && par->name != NULL )
  {	strappendf(&args[a],"envname=%s",p->envvarname);a++;
	strappendf(&args[a],"envvalue=%s",par->name);a++;
  }
 if ( p->in != NULL )
  {	char	*inname;
	inname=fileformat_replace(p->in,par);
	strappendf(&args[a],"in=%s",inname);a++;
	free(inname);
  }
 if ( p->out != NULL && p->fwout==NULL )
  {	char	*outname;
	outname=fileformat_replace(p->out,par);
	strappendf(&args[a],"out=%s",outname);a++;
	free(outname);
  }
 else if ( p->fwout != NULL )
  {	strappendf(&args[a],"out=-");a++;
  }
 if ( p->err != NULL && p->fwerr==NULL )
  {	char	*errname;
	errname=fileformat_replace(p->err,par);
	strappendf(&args[a],"err=%s",errname);a++;
	free(errname);
  }
 else if ( p->fwerr != NULL )
  {	strappendf(&args[a],"err=-");a++;
  }

 if ( par->c.is_shell )
  {	strappendf(&args[a],"shell=%s",p->shell);a++;
	strappendf(&args[a],"-");a++;
	args[a]=concatenate_arguments(par->c.argc,par->c.argv);a++;
  }
 else
  {	strappendf(&args[a],"-");a++;
	args=(char **)xrealloc(args,sizeof(char *)*(16+par->c.argc));
	for ( i=0 ; i<par->c.argc ; i++ )
	 {	args[a]=xstrdup(par->c.argv[i]);
		a++;
	 }
  }

 args[a]=NULL;

 for ( i=0 ; i<a ; i++ )
  {	token=daemon_commandtoken_escape_string(args[i]);
	if ( i<a-1 )
		hprintf(cc->rs->fhsend,"%s ",token);
	else
		hprintf(cc->rs->fhsend,"%s\n",token);

	/* fprintf(stderr,"send_task(): token: %s\n",token); */

	free(token);
	free(args[i]);
  }

 free(args);

 return(0);
}

size_t read_fd_block(int fd,void *data,size_t size)
{
 char		*cdata;
 ssize_t	rd;
 size_t		tr;

 cdata=data;
 tr=0;

 while ( size>0 )
  {	rd=read(fd,cdata,size);
	/*
	fprintf(stderr,"read_fd_block(): read:%d size=%d\n",(int)rd,(int)size);
	*/
	if ( rd<0 && errno==EINTR )	/* can be interrupted by SIGCHLD     */
		continue;
	else if ( rd<=0 )		/* otherwise, something is wrong, or */
	 {	/*
		fprintf(stderr,"error: errno=%d '%s'\n",errno,strerror(errno));
		*/
		return(0);
	 }
	else				/* or... everything is fine...	     */
	 {	size-=rd;
		cdata+=rd;
		tr+=rd;
	 }
  };

 return(tr);
}

int read_signalinfo(int fd,signalinfo *sci)
{
 return((int)read_fd_block(fd,sci,sizeof(signalinfo)));
}

int write_out_line(FILE *fw,char *format,char *lbuffer,
	parameter *p,int omit_newlines)
{
 char	*outline;
 size_t	len;

 outline=format_replace(format,1,
	'k',FORMAT_INTEGER,p->id,
	'd',FORMAT_INTEGER,p->id+1,
	's',FORMAT_STRING,p->name,
	'l',FORMAT_STRING,lbuffer,
	NULL);

 if ( outline != NULL )
  {	len=strlen(outline);
	if ( omit_newlines )
		fwrite(outline,1,len,fw);
	else
		fprintf(fw,"%s\n",outline);
	fflush(fw);
	free(outline);
  }

 return(0);
}

int write_output(char *buff,size_t n,FILE *fw,
char *format,int omit_newlines,linebuffer *lb,parameter *p)
{
 char	*line; 

 if ( fw==NULL )	
	return(0);

 if ( buff==NULL )	n=0;

 if ( format==NULL || lb==NULL )
  {	if ( n>0 && buff != NULL )
	 {	fwrite(buff,1,n,fw);
		fflush(fw);
	 }
	return(0);
  }
 else if ( format != NULL && lb != NULL )
  {	linebuffer_concatenate(lb,buff,n);
	while ( (line=linebuffer_fetch(lb)) != NULL || (buff==NULL && (line=linebuffer_flush(lb)) != NULL) )
	 {	write_out_line(fw,format,line,p,omit_newlines);
		free(line);
	 }
		
	return(0);
  }
 else
	return(-1);
}

int imutex_cleanup(imutex *mx,int qid,void *data)
{
 int	i,j,r;

 for ( ; mx != NULL ; mx=mx->next )
  { 	if ( mx->pclients==NULL || mx->npclient<=0 )
		continue;
	r=0;
	for ( i=0 ; i<mx->npclient ; i++ )
	 {	if ( mx->pclients[i].qid==qid && mx->pclients[i].data==data )
		 {	mx->pclients[i].qid=-1;
			mx->pclients[i].data=NULL;
			r++;
		 }
	 }
	if ( r==mx->npclient )
	 {	free(mx->pclients);
		mx->pclients=NULL;
		mx->npclient=0;
	 }
	else if ( r )
	 {	for ( i=0,j=0 ; i<mx->npclient ; i++ )
		 {	if ( ! ( mx->pclients[i].qid<0 && mx->pclients[i].data==NULL ) )
			 {	mx->pclients[j]=mx->pclients[i];
				j++;
			 }
		 }
		mx->npclient=j;
	 }
  }

 return(0);
}

int imutex_lock_get_params(char **cmd,char **rname,int *rmaxnum)
{
 char	*name;
 int	i,maxnum;

 name=NULL;
 maxnum=1;
 for ( i=1 ; cmd[i] != NULL ; i++ )
  {	if ( ( strcmp(cmd[i],"-m")==0 || strcmp(cmd[i],"--maximum")==0 ) && cmd[i+1] != NULL )
	 {	i++;
		sscanf(cmd[i],"%d",&maxnum);
	 }
	else if ( cmd[i][0]=='-' )
		return(1);
	else if ( name==NULL )
		name=cmd[i];
	else
		return(1);
  }

 *rname=name;
 *rmaxnum=maxnum;

 return(0);
}

int imutex_unlock_get_params(char **cmd,char **rname)
{
 char	*name;
 int	i;

 name=NULL;
 for ( i=1 ; cmd[i] != NULL ; i++ )
  {	if ( cmd[i][0]=='-' )
		return(1);
	else if ( name==NULL )
		name=cmd[i];
	else
		return(1);
  }

 *rname=name;

 return(0);
}

int imutex_mutex_get_params(char **cmd,char **rname)
{
 char	*name;
 int	i;

 name=NULL;
 for ( i=1 ; cmd[i] != NULL ; i++ )
  {	if ( cmd[i][0]=='-' )
		return(1);
	else if ( name==NULL )
		name=cmd[i];
	else
		return(1);
  }

 *rname=name;

 return(0);
}

int remote_client_master_parse_command_tokens(char **cmd,parallelstatus *ps,int qid,void *dd)
{ 
 time_t		tc;
 int		i,ret; 
 remoteshell	*rs;
 client		*cl;

 if ( qid>=0 ) 
  {	rs=(remoteshell *)dd;
	cl=NULL;
  }
 else
  {	rs=NULL;
	cl=(client *)dd;
  }

 if ( cmd==NULL || cmd[0]==NULL )
	return(0);

 else if ( strcmp(cmd[0],"status")==0 )
  {	time(&tc);
	if ( qid >=0 && rs != NULL )
	 {	hprintf(rs->fhsend,"remote %d status %d/%d/%d %d %d\n",
			qid,
			ps->nfinished,ps->npending,ps->nparam,
			(int)ps->t0,(int)tc);
		/* fprintf(stderr,"qid=%d: status sent [%d,%d,%d]\n",qid,ps->nfinished,ps->npending,ps->nparam); */
	 }
	else if ( cl != NULL )
	 {	hprintf(cl->peer,"status %d/%d/%d %d %d\n",
			ps->nfinished,ps->npending,ps->nparam,
			(int)ps->t0,(int)tc);
	 }
	ret=0;
  }
 else if ( strcmp(cmd[0],"lock")==0 )
  {	char	*name;
	int	maxnum;
	imutex	*mx;

	name=NULL;
	maxnum=0;
	i=imutex_lock_get_params(cmd,&name,&maxnum);

	/* fprintf(stderr,"pexec: %d: lock request: '%s'.\n",(int)getpid(),name); */

	if ( i || name==NULL )
		hprintf(cl->peer,"error unexpected arguments for '%s'\n",cmd[0]);
	else if ( (mx=get_imutex_by_name(ps->imutexlist,name)) != NULL )
	 {	mx->pclients=(pendingclient *)xrealloc(mx->pclients,sizeof(pendingclient)*(mx->npclient+1));
		
		mx->pclients[mx->npclient].qid=qid;
		mx->pclients[mx->npclient].data=dd;
		mx->npclient++;
		mx->state++;
	 }	
	else	
	 {	mx=list_new(imutex);
		list_insert_first(ps->imutexlist,mx);
		mx->name=xstrdup(name);
		mx->state=1;
		mx->pclients=NULL;
		mx->npclient=0;
		if ( qid>=0 && rs != NULL )
			hprintf(rs->fhsend,"remote %d locked \"%s\"\n",qid,mx->name);
		else if ( cl != NULL )
			hprintf(cl->peer,"locked \"%s\"\n",mx->name);
	 }
	ret=0;
  }
 else if ( strcmp(cmd[0],"unlock")==0 )
  {	char		*name;
	imutex		*mx;
	pendingclient	pcw;
	int		is_pcw_set;

	name=NULL;
	i=imutex_unlock_get_params(cmd,&name);

	if ( (!i) && name != NULL && (mx=get_imutex_by_name(ps->imutexlist,name)) != NULL )
	 {	
		mx->state--;

		if ( mx->pclients != NULL && mx->npclient>0 )
		 {	memcpy(&pcw,&mx->pclients[0],sizeof(pendingclient));
			is_pcw_set=1;
			if ( mx->npclient>1 )
				memmove(mx->pclients,mx->pclients+1,sizeof(pendingclient)*(mx->npclient-1));
			else
			 {	free(mx->pclients);
				mx->pclients=NULL;
			 }
			mx->npclient--;
		 }
		else
			is_pcw_set=0;

		if ( is_pcw_set )
		 {	if ( pcw.qid>=0 && pcw.data != NULL )
				hprintf(((remoteshell *)pcw.data)->fhsend,"remote %d locked \"%s\"\n",pcw.qid,mx->name);
			else if ( pcw.data != NULL )
				hprintf(((client *)pcw.data)->peer,"locked \"%s\"\n",mx->name);
		 }
	
		if ( mx->state <= 0 )
		 {	free(mx->name);
			list_remove(ps->imutexlist,mx);
			free(mx);
		 }
	 }

	ret=0;
  }
 else if ( strcmp(cmd[0],"mutex")==0 )
  {	imutex	*mx;
	char	*name;
	int	lcnt;

	name=NULL;
	i=imutex_mutex_get_params(cmd,&name);
	if ( i || name==NULL )
		hprintf(cl->peer,"error unexpected arguments for '%s'\n",cmd[0]);
	else 
	 {	if ( (mx=get_imutex_by_name(ps->imutexlist,name)) != NULL )
			lcnt=mx->npclient;
		else
			lcnt=0;
		if ( qid>=0 && rs != NULL )
			hprintf(rs->fhsend,"remote %d mutex \"%s\" %d\n",qid,name,lcnt);
		else if ( cl != NULL )
			hprintf(cl->peer,"mutex \"%s\" %d\n",name,lcnt);
	 }

	ret=0;
  }

 else if ( strcmp(cmd[0],"exit")==0 )
  {	hprintf(cl->peer,"bye\n");
	ret=1;
  }
 else
  {	hprintf(cl->peer,"error unexpected command '%s'\n",cmd[0]);
	ret=0;
  }

 return(ret);
}

int remote_client_master_parse_command(char *line,parallelstatus *ps,int qid,void *dd)
{ 
 int	ret; 
 char	**cmd;

 cmd=tokenize_spaces_dyn(line);

 if ( cmd==NULL )
	return(0);
 else if ( cmd[0]==NULL )
  {	free(cmd);
	return(0);
  }

 ret=remote_client_master_parse_command_tokens(cmd,ps,qid,dd);

 free(cmd);

 return(ret);
}

int pexec_submit(paralleldata *p,parallelstatus *ps,
	numhashtable *ntp,numhashtable *ntf,
	remoteshell *rs,parameter *cp,int n)
{
 child	*cc;

 cc=list_new(child);
 cc->fdstdout=-1;
 cc->fdstderr=-1;
 list_insert_first(ps->childlist,cc);
 ps->achild++;
 rs->achild++;
 cc->rs=rs;

 if ( ! rs->num_processes )
	rs->estatus=0;

 linebuffer_reset(&cc->lout);
 linebuffer_reset(&cc->lerr);

 /* remote submission: */
 if ( rs->pid >= 0 )
  {	cc->id=n;
	cc->pid=-1;

	send_task(p,cp,cc);

	ps->npending++;
	numhash_add(ntp,n,NULL);
  }

 /* local submission: */
 else
  {	cc->id=n;
	cc->pid=submit_task(p,cp,cc,0,ps);
	ps->npending++;
	numhash_add(ntp,n,NULL);
	if ( cc->pid < 0 )
	 {	ps->nfinished++;
		numhash_add(ntf,n,NULL);
		log_message(p->log,1,cp,"error while invoking task:"
			" %s",logmsg_submit_task[-cc->pid]);
		list_remove(ps->childlist,cc);
		free(cc);
		ps->achild--;
		rs->achild--;
 	 }
  }

 return(0);
}

remoteshell *pexec_get_free_remote_shell(remoteshell *rshells,int nrshell)
{
 remoteshell	*rs;
 int		r;

 /* check if max number of processes has been reached... */
 for ( r=0,rs=rshells ; r<nrshell ; r++,rs++ )
  {	if ( ( 0 < rs->num_processes && rs->achild < rs->num_processes ) || rs->estatus>=2 )
		return(rs);
  }

 return(NULL);
}

int pexec_do_parallelized_execution(paralleldata *p,
	parameter *params,int nparam,
	remoteshell *rshells,int nrshell,
	int sock,int hsck)
{
 child			*cc;
 struct	sigaction	chldact;
 signalinfo		sci;
 numhashtable		ntp,ntf;	/* pending & finished jobs */
 int			i,n,r;
 parameter		*cp;
 fd_set			set;
 int			spipe,max,status;
 char			*buff;
 int			buffsize;
 remoteshell		*rs;
 client			*cl;
 parallelstatus		ps;
 linebuffer		lhyp;
 
 if ( pipe(sig_pipe)<0 )
	return(-1);

 spipe=sig_pipe[0];

 chldact.sa_handler=sig_act_child;
 sigemptyset(&chldact.sa_mask);
 chldact.sa_flags=(SA_NOCLDSTOP|SA_RESTART);
 sigaction(SIGCHLD,&chldact,NULL);

 ps.childlist=NULL;
 ps.achild=0;
 /* just to be sure, [*].achild must be init'ed by remote_shell_init(): */
 for ( r=0,rs=rshells ; r<nrshell ; r++,rs++ )
  {	rs->achild=0;
	rs->estatus=0;
  }

 ps.nfinished=0;
 numhash_init(&ntf,get_bit_size(nparam),4);
 ps.npending =0;
 numhash_init(&ntp,get_bit_size(nparam),4);
 ps.nparam=nparam;

 buffsize=getpagesize();	/* some natural buffer size 	*/
 buff=(char *)xmalloc(buffsize);

 ps.sock=sock;			/* inherit the same socket fd	*/
 ps.hsck=hsck;			/* hypervisor client socket fd  */
 ps.clientlist=NULL;		/* no clients			*/
 ps.imutexlist=NULL;		/* no mutexes			*/

 linebuffer_reset(&lhyp);

 time(&ps.t0);

 while ( ps.nfinished < ps.nparam )
  {	while ( 1 )
	 {	
		n=numhash_get_smallest_free(&ntp);
		/* fprintf(stderr,"numhash_get_smallest_free()=%d\n",n); */

		/* nothing to do... */
		if ( n<0 || n>=nparam )
			break;

		rs=pexec_get_free_remote_shell(rshells,nrshell);
		if ( rs==NULL )
			break;

		/* here we know that the remote shell 'rs' is	*/ 
		/* free for submission a new submission:	*/
		cp=&params[n];

		pexec_submit(p,&ps,&ntp,&ntf,rs,cp,n);

	 }

	for ( rs=rshells,r=0 ; r<nrshell ; r++,rs++ )
	 {	if ( ! rs->num_processes && ! rs->estatus )
		 {	rs->estatus=1;
			if ( rs->pid >=0 )
				hprintf(rs->fhsend,"request\n");
			else if ( hsck>=0 )
				hprintf(hsck,"request\n");
		 }
	 }

	/* This would happen if many of the submit_task()'s fail yielding    */
	/* an unexpected increase of "finished" tasks. If the increment was  */
	/* too high and all children have been exited yet, the select()	     */
	/* would block forever (i.e. we have to stop the main loop here).    */
	if ( ps.nfinished >= nparam )
		break;

	FD_ZERO(&set);
	FD_SET(spipe,&set);
	max=spipe;

	for ( cc=ps.childlist ; cc != NULL ; cc=cc->next )
	 {	if ( cc->id<0 || cc->pid<0 )
			continue;  /* ignore invalid or non-local childeren */
		if ( cc->fdstdout >= 0 )
		 {	FD_SET(cc->fdstdout,&set);
			if ( max<cc->fdstdout )	max=cc->fdstdout;
		 }
		if ( cc->fdstderr >= 0 )
		 {	FD_SET(cc->fdstderr,&set);
			if ( max<cc->fdstderr )	max=cc->fdstderr;
		 }
	 }

	for ( rs=rshells,r=0 ; r<nrshell ; r++,rs++ )
	 {	if ( rs->pid < 0 || rs->fhsend<0 || rs->fhrecv<0 )
			continue; /* ignore local "remote shells" */
		FD_SET(rs->fhrecv,&set);
		if ( max<rs->fhrecv )	max=rs->fhrecv;
	 }

	if ( sock>=0 )
	 {	FD_SET(sock,&set);
		if ( max<sock )		max=sock;
	 }
	if ( hsck>=0 )
	 {	FD_SET(hsck,&set);
		if ( max<hsck )		max=hsck;
	 }

	for ( cl=ps.clientlist ; cl != NULL ; cl=cl->next )
	 {	FD_SET(cl->peer,&set);
		if ( max<cl->peer )	max=cl->peer;
	 }

	i=select(max+1,&set,NULL,NULL,NULL);

	if ( i<0 && errno==EINTR )
		continue;

	/* somebody is dying: */
	if ( FD_ISSET(spipe,&set) )
	 {	if ( ! (read_signalinfo(spipe,&sci)>0) )
		 {	fprintf(stderr,_("%s: fatal error: read_signalinfo() failed.\n"),progbasename);
			exit(1);
		 }

		/*
		fprintf(stderr,"read_signalinfo(): pid=%d signal=%d exit=%d\n",sci.pid,sci.exitsignal,sci.exitstatus);
		*/

		/* normal termination, successful: */
		if ( sci.exitsignal<0 )
		 {	cc=get_child_by_pid(ps.childlist,sci.pid);
			n=cc->id;
			cc->id =-1;
			cc->pid=-1;
			cp=&params[n];
			numhash_add(&ntf,n,NULL);
			ps.nfinished++;
			/* LOG */
			if ( cc->fdstdout>=0 && p->fwout != NULL )
			 {	while ( fd_avail(cc->fdstdout) )
				 {	n=read(cc->fdstdout,buff,buffsize);	
					if ( n>0 )
						write_output(buff,n,p->fwout,p->formatout,p->omit_newlines,&cc->lout,cp);
					else if ( ! ( n<0 && errno==EINTR ) )
						break;
				 }
				write_output(NULL,0,p->fwout,p->formatout,p->omit_newlines,&cc->lout,cp);
				close(cc->fdstdout);
				cc->fdstdout=-1;
			 }	
			if ( cc->fdstderr>=0 && p->fwerr != NULL )
			 {	while ( fd_avail(cc->fdstderr) )
				 {	n=read(cc->fdstderr,buff,buffsize);
					if ( n>0 )
						write_output(buff,n,p->fwerr,p->formaterr,p->omit_newlines,&cc->lerr,cp);
					else if ( ! ( n<0 && errno==EINTR ) )
						break;
				 }
				write_output(NULL,0,p->fwerr,p->formaterr,p->omit_newlines,&cc->lerr,cp);
				close(cc->fdstderr);
				cc->fdstderr=-1;
			 }
			list_remove(ps.childlist,cc);
			cc->rs->achild--;
			ps.achild--;
			if ( ! cc->rs->num_processes && hsck>=0 )
				hprintf(hsck,"ready\n");
			free(cc);
		 }
		/* abnormal termination, caused by some discrepancy: */
		else if ( p->fallback_to_die )
		 {	close(sig_pipe[0]);
			sig_pipe[0]=-1;
			close(sig_pipe[1]);
			sig_pipe[1]=-1;
			for ( cc=ps.childlist ; cc != NULL ; cc=cc->next )
			 {	if ( cc->pid < 0 )
					continue;	/* ignore nonlocal */
				kill(cc->pid,SIGKILL);
				waitpid(cc->pid,&status,0);
				free(cc);
			 }
			ps.childlist=NULL;
			ps.achild=0;
			return(1);
		 }
		/* resubmit after the abnormal termination: */	
		else
		 {	cc=get_child_by_pid(ps.childlist,sci.pid);
			numhash_remove(&ntp,n);
			ps.npending--;
			list_remove(ps.childlist,cc);
			cc->rs->achild--;
			ps.achild--;
			free(cc);
		 }
	 } /* if ( FD_ISSET(spipe,&set) ) */

	for ( cc=ps.childlist ; cc != NULL ; cc=cc->next )
 	 {	if ( cc->id<0 || cc->pid<0 )
			continue;	/* skip invalid or non-local children */
		cp=&params[cc->id];
		if ( cc->fdstdout >= 0 && p->fwout != NULL && FD_ISSET(cc->fdstdout,&set) )
		 {	n=read(cc->fdstdout,buff,buffsize);
			if ( n>0 )
				write_output(buff,n,p->fwout,p->formatout,p->omit_newlines,&cc->lout,cp);
			else if ( ! ( n<0 && errno==EINTR ) )
			 {	close(cc->fdstdout);
				cc->fdstdout=-1;
			 }
		 }
		if ( cc->fdstderr >= 0 && p->fwerr != NULL && FD_ISSET(cc->fdstderr,&set) )
		 {	n=read(cc->fdstderr,buff,buffsize);
			if ( n>0 )
				write_output(buff,n,p->fwerr,p->formaterr,p->omit_newlines,&cc->lerr,cp);
			else if ( ! ( n<0 && errno==EINTR ) )
			 {	close(cc->fdstderr);
				cc->fdstderr=-1;
			 }
		 }
	 } /* for ( cc=ps.childlist ; cc != NULL ; cc=cc->next ) */

	/* some data has been received from the tunnel of the remote pexec: */
	for ( rs=rshells,r=0 ; r<nrshell ; r++,rs++ )
	 {	int	len;
		char	*line,**cmd;

		if ( rs->pid < 0 || rs->fhsend<0 || rs->fhrecv<0 )
			continue; /* ignore local "remote shells" */
		if ( ! FD_ISSET(rs->fhrecv,&set) )
			continue;

		n=read(rs->fhrecv,buff,buffsize);
		/* the remote shell process terminated unexpectedly: */
		if ( n==0 )
		 {	/* TBD */ 
			break;
		 }
		/* nothing, just a signal interrupted the reading, go on...*/
		else if ( n<0 && errno==EINTR )
			continue;

		linebuffer_concatenate(&rs->lrsh,buff,n);

		while ( (line=linebuffer_fetch(&rs->lrsh)) != NULL )
		 {	/* fprintf(stderr,"pexec: %d: remote shell: received: '%s'\n",(int)getpid(),line);  */
			cmd=tokenize_spaces_dyn(line);
			if ( cmd != NULL && cmd[0] != NULL )
			 {	if ( strcmp(cmd[0],"output")==0 && cmd[1] != NULL && cmd[2] != NULL )
				 {	sscanf(cmd[1],"%d",&n);
					cc=get_child_by_id(ps.childlist,n);
					cp=&params[n];
					len=daemon_commandtoken_unescape(cmd[2]);
					write_output(cmd[2],len,p->fwout,p->formatout,p->omit_newlines,&cc->lout,cp);
				 }
				else if ( strcmp(cmd[0],"error")==0 && cmd[1] != NULL && cmd[2] != NULL )
				 {	sscanf(cmd[1],"%d",&n);
					cc=get_child_by_id(ps.childlist,n);
					cp=&params[n];
					len=daemon_commandtoken_unescape(cmd[2]);
					write_output(cmd[2],len,p->fwerr,p->formaterr,p->omit_newlines,&cc->lerr,cp);
				 }
				else if ( strcmp(cmd[0],"finish")==0 && cmd[1] != NULL )
				 {	sscanf(cmd[1],"%d",&n);
					cc=get_child_by_id(ps.childlist,n);
					numhash_add(&ntf,n,NULL);
					ps.nfinished++;
					list_remove(ps.childlist,cc);
					cc->rs->achild--;
					ps.achild--;
					if ( ! cc->rs->num_processes )
						hprintf(cc->rs->fhsend,"ready\n");
					free(cc);
				 }
				else if ( strcmp(cmd[0],"remote")==0 && cmd[1] != NULL && cmd[2] != NULL )
				 {	int	qid;
					void	*dd;
					if ( sscanf(cmd[1],"%d",&qid)==1 )
					 {	if ( qid>=0 )	dd=(void *)rs;
						else		dd=NULL;
						remote_client_master_parse_command_tokens(cmd+2,&ps,qid,dd);
					 }
				 }
				else if ( strcmp(cmd[0],"acknowledged")==0 )
					rs->estatus=2;
				/* TBD: handle unexpected data: */
				/*
				else
				 {
				 }
				*/
			 }
			if ( cmd != NULL )				
				free(cmd);
			free(line);
		 }

	 } /* for ( rs=rshells,r=0 ; r<nrshell ; r++,rs++ ) */

	for ( cl=ps.clientlist ; cl != NULL ; cl=cl->next )
	 {	char	*line;
		int	closepeer;

		if ( ! FD_ISSET(cl->peer,&set) )
			continue;

		n=read(cl->peer,buff,buffsize);
		if ( n==0 )
		 {	close(cl->peer);
			linebuffer_free(&cl->lcli);
			imutex_cleanup(ps.imutexlist,-1,(void *)cl);
			list_remove(ps.clientlist,cl);
			free(cl);
			break;
		 }
		else if ( n<0 && errno==EINTR )
			continue;

		linebuffer_concatenate(&cl->lcli,buff,n);
		closepeer=0;
		
		while ( (line=linebuffer_fetch(&cl->lcli)) != NULL )
		 {	closepeer=remote_client_master_parse_command(line,&ps,-1,(void *)cl);
			free(line);
			if ( closepeer )	break;
		 };

		if ( closepeer )
		 {	close(cl->peer);
			linebuffer_free(&cl->lcli);
			imutex_cleanup(ps.imutexlist,-1,(void *)cl);
			list_remove(ps.clientlist,cl);
			free(cl);
			break;
		 }
	 } /* for ( cl=ps.clientlist ; cl != NULL ; cl=cl->next ) */

	if ( hsck>=0 && FD_ISSET(hsck,&set) )
	 {	char	*line,**cmd;

		n=read(hsck,buff,buffsize);
		if ( n>0 )
			linebuffer_concatenate(&lhyp,buff,n);

		while ( (line=linebuffer_fetch(&lhyp)) != NULL )
		 {	cmd=tokenize_spaces_dyn(line);
			if ( cmd != NULL )
			 {	if ( strcmp(cmd[0],"acknowledged")==0 )
				 {	for ( rs=rshells,r=0 ; r<nrshell ; r++,rs++ )
					 {	if ( ! rs->num_processes )
							rs->estatus=2;
					 }
				 }
				free(cmd);
			 }
			free(line);
		 }

 	 } /* if ( hsck>=0 && FD_ISSET(hsck,&set) ) */

	if ( sock>=0 && FD_ISSET(sock,&set) ) 
	 {	cl=list_new(client);
		cl->peer=accept(sock,NULL,NULL);
		linebuffer_reset(&cl->lcli);
		list_insert_first(ps.clientlist,cl);
		cl->status=0;
	 } /* if ( sock>=0 && FD_ISSET(sock,&set) )  */
	
  }

 for ( rs=rshells,r=0 ; r<nrshell ; r++,rs++ )
  {	if ( ! rs->num_processes )
	 {	if ( rs->pid >=0 )
			hprintf(rs->fhsend,"completed\n");
		else if ( hsck>=0 )
			hprintf(hsck,"completed\n");
	 }
  }

 free(buff);
 
 numhash_free(&ntp);
 numhash_free(&ntf);

 signal(SIGCHLD,SIG_DFL);

 close(sig_pipe[0]);
 close(sig_pipe[1]);

 return(0); 
}

int remote_client_daemon_parse_command(char *line,client *cl,parallelstatus *ps,int fhsend)
{ 
 int	i,ret; 
 char	**cmd;

 cmd=tokenize_spaces_dyn(line);

 if ( cmd==NULL )
	return(0);
 else if ( cmd[0]==NULL )
  {	free(cmd);
	return(0);
  }
 else if ( strcmp(cmd[0],"status")==0 )
  {	dqueue	*dq;

	dq=list_new(dqueue);
	dq->id=ps->iqueue;
	dq->qclient=cl;
	list_insert_first(ps->dqueuelist,dq);

	hprintf(fhsend,"remote %d status\n",dq->id);

	ps->iqueue++;

	ret=0;
  }
 else if ( strcmp(cmd[0],"lock")==0 )
  {	char	*name;
	int	maxnum;
	dqueue	*dq;

	name=NULL;
	maxnum=0;
	i=imutex_lock_get_params(cmd,&name,&maxnum);
	if ( i || name==NULL )
		hprintf(cl->peer,"error unexpected arguments for '%s'\n",cmd[0]);
	else
	 {	dq=list_new(dqueue);
		dq->id=ps->iqueue;
		dq->qclient=cl;
		list_insert_first(ps->dqueuelist,dq);
		hprintf(fhsend,"remote %d lock \"%s\"\n",dq->id,name);
		ps->iqueue++;
	 }

	ret=0;
  }
 else if ( strcmp(cmd[0],"unlock")==0 )
  {	char	*name;

	name=NULL;
	i=imutex_unlock_get_params(cmd,&name);
	if ( ! ( i || name==NULL ) )
		hprintf(fhsend,"remote -1 unlock \"%s\"\n",name);

	ret=0;
  }
 else if ( strcmp(cmd[0],"mutex")==0 )
  {	char	*name;
	dqueue	*dq;

	name=NULL;
	i=imutex_mutex_get_params(cmd,&name);
	if ( i || name==NULL )
		hprintf(cl->peer,"error unexpected arguments for '%s'\n",cmd[0]);
	else
	 {	dq=list_new(dqueue);
		dq->id=ps->iqueue;
		dq->qclient=cl;
		list_insert_first(ps->dqueuelist,dq);
		hprintf(fhsend,"remote %d mutex \"%s\"\n",dq->id,name);
		ps->iqueue++;
	 }

	ret=0;
  }

 else if ( strcmp(cmd[0],"exit")==0 )
  {	hprintf(cl->peer,"bye\n");
	ret=1;
  }
 else
  {	hprintf(cl->peer,"error unexpected command '%s'\n",cmd[0]);
	ret=0;
  }

 free(cmd);

 return(ret);
}

dqueue *dqueue_get_queue_by_id(parallelstatus *ps,int id)
{
 dqueue *dq;
 for ( dq=ps->dqueuelist ; dq != NULL ; dq=dq->next )
  {	if ( dq->id == id )
		return(dq);
  }
 return(NULL);
}

int daemon_process_command(parallelstatus *ps,int fhsend,char **cmd)
{
 int		i,n,id,zeroarg;
 char		*shell,*in,*out,*err,*envname,*envvalue,*eoc;
 paralleldata	p;
 parameter	par;
 child		*cc;
 dqueue		*dq;

 if ( cmd==NULL || cmd[0]==NULL )
	return(0);

 for ( n=0 ; cmd[n] != NULL ; )	n++;

 if ( strcmp(cmd[0],"exit")==0 )
	return(1);

 else if ( strcmp(cmd[0],"remote")==0 && n>=3 )
  {	int	qid,j;
	char	*buff;

	if ( sscanf(cmd[1],"%d",&qid)<1 )
		qid=-1;

	dq=dqueue_get_queue_by_id(ps,qid);
	if ( dq != NULL )
	 {	buff=NULL;

		strappendf(&buff,"%s",cmd[2]);

		if ( strcmp(cmd[2],"locked")==0 || strcmp(cmd[2],"mutex")==0 )
		 {	strappendf(&buff," \"%s\"",(cmd[3]!=NULL?cmd[3]:"-"));
			if ( cmd[3] != NULL )	j=4;
			else			j=3;
		 }
		else
			j=3;
		for ( ; cmd[j] != NULL ; j++ )
			strappendf(&buff," %s",cmd[j]);

		hprintf(dq->qclient->peer,"%s\n",buff);
		if ( buff != NULL )	free(buff);
	
		list_remove(ps->dqueuelist,dq);
		free(dq);
	 }
	/*
	else
		fprintf(stderr,"daemon pexec: %d: unexpected queue id: %d\n",(int)getpid(),qid);
	*/

	return(0);
  }

 else if ( strcmp(cmd[0],"execute")==0 )
  {	shell=NULL;
	in=NULL;
	out=NULL;
	err=NULL;
	envname=NULL;
	envvalue=NULL;
	id=-1;

	zeroarg=-1;
	for ( i=1 ; i<n ; i++ )
	 {	if ( cmd[i][0]=='-' )
		 {	zeroarg=i+1;
			break;
		 }
		eoc=strchr(cmd[i],'=');
		if ( eoc==NULL )
		 {	hprintf(fhsend,"message \"invalid argument '%s' for '%s'\"\n",cmd[i],cmd[0]);
			return(0);
		 }
		*eoc=0;
		eoc++;

		if ( strcmp(cmd[i],"id")==0 || strcmp(cmd[i],"identifier")==0 )
			sscanf(eoc,"%d",&id);	
		else if ( strcmp(cmd[i],"shell")==0 )
			shell=eoc;
		else if ( strcmp(cmd[i],"in")==0 )
			in=eoc;
		else if ( strcmp(cmd[i],"out")==0 )
			out=eoc;
		else if ( strcmp(cmd[i],"err")==0 )
			err=eoc;
		else if ( strcmp(cmd[i],"envname")==0 )
			envname=eoc;
		else if ( strcmp(cmd[i],"envvalue")==0 )
			envvalue=eoc;
		else 
		 {	hprintf(fhsend,"message \"invalid argument '%s' for '%s'\"\n",cmd[i],cmd[0]);
			return(0);
		 }
	 }
	if ( zeroarg<0 || id<0 )
	 {	hprintf(fhsend,"message \"invalid syntax near '%s'\"\n",cmd[0]);
		return(0);
	 }
	else if ( zeroarg>=n )
	 {	hprintf(fhsend,"message \"command specification is missing near '%s'\"\n",cmd[0]);
		return(0);
	 }

	if ( in==NULL )
		p.in=NULL;
	else
		p.in=in;

	if ( out==NULL )
	 {	p.out=NULL;
		p.fwout=NULL;
	 }
	else if ( strcmp(out,"-")==0 )
	 {	p.out=out;
		p.fwout=stdout;	/* just something which is definitely not NULL */
	 }
	else
	 {	p.out=out;
		p.fwout=NULL;
	 }

	if ( err==NULL )
	 {	p.err=NULL;
		p.fwerr=NULL;
	 }
	else if ( strcmp(err,"-")==0 )
	 {	p.err=err;
		p.fwerr=stderr;	/* just something which is definitely not NULL */
	 }
	else
	 {	p.err=err;
		p.fwerr=NULL;
	 }

	if ( envname != NULL && envvalue != NULL )
	 {	p.envvarname=envname;
		par.name=envvalue;
	 }
	else
	 {	p.envvarname=NULL;
		par.name=NULL;
	 }

	if ( shell != NULL )
	 {	p.shell=shell;
		par.c.is_shell=1;
		par.c.argc=n-zeroarg;
		par.c.argv=cmd+zeroarg;
	 }
	else
	 {	p.shell=NULL;
		par.c.is_shell=0;
		par.c.argc=n-zeroarg;
		par.c.argv=cmd+zeroarg;
	 }
	par.no_touch_std=0;

	cc=list_new(child);
	cc->fdstdout=-1;
	cc->fdstderr=-1;
	cc->id=id;
	/* it was: submit_task(&p,&par,cc,!0,NULL); */
	cc->pid=submit_task(&p,&par,cc,!0,ps);
	if ( cc->pid<0 )
	 {	hprintf(fhsend,"message \"unable to execute, reason code: %d\"\n",cc->pid);
		free(cc);
		return(0);
	 }
	else
	 {	list_insert_first(ps->childlist,cc);
		return(0);
	 }
	
  }
 
 else if ( strcmp(cmd[0],"request")==0 )
  {	if ( ps->hsck < 0 )
	 {	hprintf(fhsend,"message \"unexpected 'request': hypervisor has not been connected to daemon\"\n");	}
	else
	 {	hprintf(ps->hsck,"request\n");		}
	return(0);
  }

 else if ( strcmp(cmd[0],"ready")==0 )
  {	if ( ps->hsck < 0 )
	 {	hprintf(fhsend,"message \"unexpected 'ready': hypervisor has not been connected to daemon\"\n");	}
	else
	 {	hprintf(ps->hsck,"ready\n");		}
	return(0);
  }

 else if ( strcmp(cmd[0],"completed")==0 )
  {	if ( ps->hsck < 0 )
	 {	hprintf(fhsend,"message \"unexpected 'completed': hypervisor has not been connected to daemon\"\n");	}
	else
	 {	hprintf(ps->hsck,"completed\n");		}
	return(0);
  }

 else
  {	hprintf(fhsend,"message \"invalid command: %s\"\n",cmd[0]);
	return(0);
  }

}

int daemon_send_data(int fhsend,char *streamname,int id,char *buff,int size)
{
 char	*out;

 if ( buff==NULL )
	return(0);

 out=daemon_commandtoken_escape(buff,size);

 hprintf(fhsend,"%s %d ",streamname,id);
 write(fhsend,out,strlen(out));
 free(out);
 hprintf(fhsend,"\n");

 return(0);
}

int daemon_send_output(int fhsend,int id,char *buff,int size)
{
 return(daemon_send_data(fhsend,"output",id,buff,size));
}
int daemon_send_error(int fhsend,int id,char *buff,int size)
{
 return(daemon_send_data(fhsend,"error",id,buff,size));
}

int pexec_daemon_main_loop(int fhrecv,int fhsend,
	int num_processes,
	int sock,int hsck)
{
 fd_set			set;
 parallelstatus		ps;
 child			*cc;
 struct	sigaction	chldact;
 signalinfo		sci;
 int			spipe,max,status;
 int			buffsize,i,n;
 char			*buff,**cmd,*line;
 int			is_in_loop,ret;
 linebuffer		lcmd,lhyp;
 client			*cl;

 hprintf(fhsend,"initialization num_processes=%d\n",num_processes);

 buffsize=getpagesize();	/* some natural buffer size */
 buff=(char *)xmalloc(buffsize);

 if ( pipe(sig_pipe)<0 )
	return(-1);

 spipe=sig_pipe[0];

 chldact.sa_handler=sig_act_child;
 sigemptyset(&chldact.sa_mask);
 chldact.sa_flags=(SA_NOCLDSTOP|SA_RESTART);
 sigaction(SIGCHLD,&chldact,NULL);

 ps.childlist=NULL;
 ps.achild=0;

 ps.clientlist=NULL;
 ps.imutexlist=NULL;
 ps.dqueuelist=NULL;
 ps.iqueue=0;

 ps.sock=sock;		/* remote control server socket */
 ps.hsck=hsck;		/* hypervisor client socket	*/
 
 linebuffer_reset(&lcmd);
 linebuffer_reset(&lhyp);

 is_in_loop=1;

 while ( is_in_loop ) 
  {	
	FD_ZERO(&set);

	FD_SET(spipe,&set);
	max=spipe;
	FD_SET(fhrecv,&set);
	if ( max<fhrecv )	max=fhrecv;

	for ( cc=ps.childlist ; cc != NULL ; cc=cc->next )
	 {	if ( cc->id<0 || cc->pid<0 )	
			continue;	/* just to be sure... */
		if ( cc->fdstdout >= 0 )
		 {	FD_SET(cc->fdstdout,&set);
			if ( max<cc->fdstdout )	max=cc->fdstdout;
		 }
		if ( cc->fdstderr >= 0 )
		 {	FD_SET(cc->fdstderr,&set);
			if ( max<cc->fdstderr )	max=cc->fdstderr;
		 }
	 }

	if ( sock>=0 )
	 {	FD_SET(sock,&set);
		if ( max<sock )		max=sock;
	 }
	if ( hsck>=0 )
	 {	FD_SET(hsck,&set);
		if ( max<hsck )		max=hsck;
	 }

	for ( cl=ps.clientlist ; cl != NULL ; cl=cl->next )
	 {	FD_SET(cl->peer,&set);
		if ( max<cl->peer )	max=cl->peer;
	 }

	/* hprintf(1,"select(): waiting...\n"); */
	i=select(max+1,&set,NULL,NULL,NULL);
	/* hprintf(1,"select(): done.\n"); */

	/* somebody is dying: */
	if ( FD_ISSET(spipe,&set) )
	 {	/* hprintf(1,"read_signalinfo(): waiting...\n"); */
		if ( ! (read_signalinfo(spipe,&sci)>0) )
		 {	fprintf(stderr,_("%s: fatal error: read_signalinfo() failed.\n"),progbasename);
			exit(1);
		 }
		/* hprintf(1,"read_signalinfo(): done.\n"); */

		/*
		fprintf(stderr,"read_signalinfo(): pid=%d signal=%d exit=%d\n",sci.pid,sci.exitsignal,sci.exitstatus);
		*/

		/* normal termination, successful: */
		if ( sci.exitsignal<0 )
		 {	cc=get_child_by_pid(ps.childlist,sci.pid);
			/* LOG */
			if ( cc->fdstdout>=0 )
			 {	while ( fd_avail(cc->fdstdout) )
				 {	n=read(cc->fdstdout,buff,buffsize);	
					if ( n>0 )
						daemon_send_output(fhsend,cc->id,buff,n);
					else if ( ! ( n<0 && errno==EINTR ) )
						break;
				 }
				close(cc->fdstdout);
				cc->fdstdout=-1;
			 }	
			if ( cc->fdstderr>=0 )
			 {	while ( fd_avail(cc->fdstderr) )
				 {	n=read(cc->fdstderr,buff,buffsize);
					if ( n>0 )
						daemon_send_error(fhsend,cc->id,buff,n);
					else if ( ! ( n<0 && errno==EINTR ) )
						break;
				 }
				close(cc->fdstderr);
				cc->fdstderr=-1;
			 }
		 }

		hprintf(fhsend,"finish %d signal=%d status=%d\n",cc->id,sci.exitsignal,sci.exitstatus);

		list_remove(ps.childlist,cc);
		free(cc);
		ps.achild--;

		continue;

	 } /* if ( FD_ISSET(spipe,&set) ) */

	/* some data has been received from the tunnel of the parent pexec: */
	if ( FD_ISSET(fhrecv,&set) )
	 {	/* hprintf(1,"read(fhrecv,...): waiting...\n"); */
		n=read(fhrecv,buff,buffsize);
		if ( n==0 )
		 {	is_in_loop=0;
			break;
		 }
		else if ( n<0 && errno==EINTR )
			continue;

		linebuffer_concatenate(&lcmd,buff,n);

		while ( (line=linebuffer_fetch(&lcmd)) != NULL )
		 {	cmd=tokenize_spaces_dyn(line);
			if ( cmd != NULL )
			 {	for ( i=0 ; cmd[i] != NULL ; i++ )
					daemon_commandtoken_unescape(cmd[i]);
				ret=daemon_process_command(&ps,fhsend,cmd);
				free(cmd);
				free(line);
				if ( ret>0 )
				 {	is_in_loop=0;
					break;
				 }
			 }
			else
				free(line);
		 }

		continue;

	 } /* if ( FD_ISSET(fhrecv,&set) ) */

	for ( cc=ps.childlist ; cc != NULL ; cc=cc->next )
 	 {	if ( cc->id<0 || cc->pid<0 )
			continue;
		if ( cc->fdstdout >= 0 && FD_ISSET(cc->fdstdout,&set) )
		 {	n=read(cc->fdstdout,buff,buffsize);
			if ( n>0 )
				daemon_send_output(fhsend,cc->id,buff,n);
			else if ( ! ( n<0 && errno==EINTR ) )
			 {	close(cc->fdstdout);
				cc->fdstdout=-1;
			 }
		 }
		if ( cc->fdstderr >= 0 && FD_ISSET(cc->fdstderr,&set) )
		 {	n=read(cc->fdstderr,buff,buffsize);
			if ( n>0 )
				daemon_send_error(fhsend,cc->id,buff,n);
			else if ( ! ( n<0 && errno==EINTR ) )
			 {	close(cc->fdstderr);
				cc->fdstderr=-1;
			 }
		 }
	 } /* for ( cc in ps.childlist ) */

	for ( cl=ps.clientlist ; cl != NULL ; cl=cl->next )
	 {	char	*line;
		int	closepeer;

		if ( ! FD_ISSET(cl->peer,&set) )
			continue;

		n=read(cl->peer,buff,buffsize);
		if ( n==0 )
		 {	close(cl->peer);
			linebuffer_free(&cl->lcli);
			list_remove(ps.clientlist,cl);
			free(cl);
			break;
		 }
		else if ( n<0 && errno==EINTR )
			continue;

		linebuffer_concatenate(&cl->lcli,buff,n);
		closepeer=0;
		
		while ( (line=linebuffer_fetch(&cl->lcli)) != NULL )
		 {	/* fprintf(stderr,"daemon pexec: %d: line[%d]='%s'\n",(int)getpid(),cl->peer,line); */
			closepeer=remote_client_daemon_parse_command(line,cl,&ps,fhsend);
			free(line);
			if ( closepeer )	break;
		 };

		if ( closepeer )
		 {	close(cl->peer);
			linebuffer_free(&cl->lcli);
			list_remove(ps.clientlist,cl);
			free(cl);
			break;
		 }
	 } /* for ( cl in ps.clientlist ) */

	if ( hsck>=0 && FD_ISSET(hsck,&set) )
	 {
		n=read(hsck,buff,buffsize);
		if ( n>0 )
			linebuffer_concatenate(&lhyp,buff,n);

		while ( (line=linebuffer_fetch(&lhyp)) != NULL )
		 {	cmd=tokenize_spaces_dyn(line);
			if ( cmd != NULL )
			 {	if ( strcmp(cmd[0],"acknowledged")==0 )
				 {	hprintf(fhsend,"acknowledged\n");	}
				free(cmd);
			 }
			free(line);
		 }

 	 } /* if ( hsck>=0 && FD_ISSET(hsck,&set) ) */

	if ( sock>=0 && FD_ISSET(sock,&set) ) 
	 {	cl=list_new(client);
		cl->peer=accept(sock,NULL,NULL);
		linebuffer_reset(&cl->lcli);
		list_insert_first(ps.clientlist,cl);
		cl->status=0;
	 } /* is sock>=0 && FD_ISSET(sock) */

  }

 linebuffer_free(&lcmd);

 free(buff);

 signal(SIGCHLD,SIG_DFL);

 close(sig_pipe[0]); sig_pipe[0]=-1;
 close(sig_pipe[1]); sig_pipe[1]=-1;

 for ( cc=ps.childlist ; cc != NULL ; cc=cc->next )
  {	if ( cc->pid < 0 )
		continue;
	kill(cc->pid,SIGKILL);
	waitpid(cc->pid,&status,0);
	free(cc);
  }
 ps.childlist=NULL;
 ps.achild=0;

 return(0);
}

/*****************************************************************************/

int parse_host_data(char *arg,remotehost **rrhosts,int *rnrhost)
{
 char 		*hostlist,**hosts,*eoc;
 int		j,n,r;
 remotehost	*rhosts;
 int		nrhost;

 if ( arg==NULL || *arg==0 )
	return(-1);
		
 hostlist=xstrdup(arg);
 hosts=tokenize_char_dyn(hostlist,',');
 rhosts=NULL;
 nrhost=0;

 for ( j=0 ; hosts != NULL && hosts[j] != NULL ; j++ )
  {	eoc=strchr(hosts[j],':');
	rhosts=(remotehost *)xrealloc(rhosts,sizeof(remotehost)*(nrhost+1));
	if ( eoc==NULL )
	 {	if ( strcmp(hosts[j],PEXEC_ISTR_AUTO)==0 )
		 {	rhosts[nrhost].hostspec=NULL;
			rhosts[nrhost].num_processes=PEXEC_MNP_AUTO;
		 }
		else if ( strcmp(hosts[j],PEXEC_ISTR_MANAGED)==0 )
		 {	rhosts[nrhost].hostspec=NULL;
			rhosts[nrhost].num_processes=PEXEC_MNP_MANAGED;
		 }
		else if ( strcmp(hosts[j],PEXEC_ISTR_NCPU)==0 )
		 {	rhosts[nrhost].hostspec=NULL;
			rhosts[nrhost].num_processes=PEXEC_MNP_NCPU;
		 }
		else if ( sscanf(hosts[j],"%d",&n)==1 && n>0 )
		 {	rhosts[nrhost].hostspec=NULL;
			rhosts[nrhost].num_processes=n;
		 }
		else
		 {	rhosts[nrhost].hostspec=xstrdup(hosts[j]);
			rhosts[nrhost].num_processes=PEXEC_MNP_AUTO;
		 }
		/*
		else 
		 {	free(rhosts);
			return(1);
		 }
		*/

	 }
	else if ( eoc==hosts[j] )
	 {	free(rhosts);
		return(1);
	 }
	else
	 {	*eoc=0;
		if ( strcmp(eoc+1,PEXEC_ISTR_AUTO)==0 )
			n=PEXEC_MNP_AUTO;
		else if ( strcmp(eoc+1,PEXEC_ISTR_MANAGED)==0 )
			n=PEXEC_MNP_MANAGED;
		else if ( strcmp(eoc+1,PEXEC_ISTR_NCPU)==0 )
			n=PEXEC_MNP_NCPU;
		else if ( (r=sscanf(eoc+1,"%d",&n))==0 || (r==1 && n<=0) )
		 {	free(rhosts);
			free(hostlist);
			return(1);
		 }

		rhosts[nrhost].hostspec=xstrdup(hosts[j]);
		rhosts[nrhost].num_processes=n;
	 }
	nrhost++;

  }

 free(hosts); 
 free(hostlist);

 *rrhosts=rhosts;
 *rnrhost=nrhost;

 return(0);
}

int remote_shell_init(char *rsh,char **rshargs,char *pexec_self,char *ctrlport,
	int timeout,remotehost *rh,remoteshell *rs,int prio)
{
 int	pipesend[2],piperecv[2];
 int	pid;

 if ( rh->hostspec==NULL )
  {	rs->pid=-1;
	rs->fhsend=rs->fhrecv=-1;
	rs->num_processes=rh->num_processes;
	rs->achild=0;
	linebuffer_reset(&rs->lrsh);
	return(0);
  }

 if ( pipe(pipesend)<0 )
	return(-1);	/* creation of pipe failed */
 if ( pipe(piperecv)<0 )
	return(-1);	/* creation of pipe failed */

 pid=fork();

 if ( pid<0 )
	return(-2);	/* fork failed */

 else if ( pid>0 )
  {	char	*line,**tokens;

	rs->pid=pid;

	rs->fhsend=pipesend[1];
	close(pipesend[0]);
	rs->fhrecv=piperecv[0];
	close(piperecv[1]);
	
	rs->achild=0;

	linebuffer_reset(&rs->lrsh);

	line=linebuffer_read_line(rs->fhrecv,&rs->lrsh,timeout);

	if ( line==NULL )
		return(-3);	/* data retrieval failed */

	/* fprintf(stderr,"received: %s\n",line); */

	tokens=tokenize_spaces_dyn(line);

	if ( tokens != NULL && 
	tokens[0] != NULL && strcmp(tokens[0],"initialization")==0 && 
	tokens[1] != NULL && sscanf(tokens[1],"num_processes=%d",&rs->num_processes)==1 &&
	rs->num_processes >= 0 )
	 {	free(tokens);
		free(line);
		return(0);	/* successful */
	 }

	else
	 {	if ( tokens != NULL )	free(tokens);
		free(line);
		return(-4);	/* unexpected data received */
	 }
  }

 else
  {	int	i,n;
	char	**argv,num_buff[32],pri_buff[32];

	close(pipesend[1]);
	close(piperecv[0]);

	close(0);
	dup2(pipesend[0],0);
	close(pipesend[0]);
	close(1);
	dup2(piperecv[1],1);
	close(piperecv[1]);
	close(2);
	dup2(1,2);

	for ( n=0 ; rshargs[n] != NULL ; )	n++;

	argv=(char **)xmalloc(sizeof(char *)*(n+32));
	for ( i=0 ; i<n ; )
	 {	argv[i]=rshargs[i];
		i++;
	 }
	argv[i++]=rh->hostspec;
	argv[i++]=pexec_self;
	argv[i++]="--tunnel";
	if ( rh->num_processes>0 )
	 {	argv[i++]="--number";
		sprintf(num_buff,"%d",rh->num_processes);
		argv[i++]=num_buff;
	 }
	else if ( rh->num_processes == PEXEC_MNP_AUTO )
	 {	argv[i++]="--number";
		argv[i++]=PEXEC_ISTR_AUTO;
	 }
	else if ( rh->num_processes == PEXEC_MNP_MANAGED )
	 {	argv[i++]="--number";
		argv[i++]=PEXEC_ISTR_MANAGED;
	 }
	else if ( rh->num_processes == PEXEC_MNP_NCPU )
	 {	argv[i++]="--number";
		argv[i++]=PEXEC_ISTR_NCPU;
	 }

	if ( prio > 0 )
	 {	argv[i++]="--nice";
		sprintf(pri_buff,"%d",prio-128);
		argv[i++]=pri_buff;
	 }
	if ( ctrlport != NULL )
	 {	argv[i++]="--bind";
		argv[i++]=ctrlport;
	 }
	argv[i++]=NULL;	
	
	execvp(rsh,argv);

	fprintf(stdout,"initialization execution=failed\n");

	exit(1);

	return(0);
  }

}

/*****************************************************************************/

#ifndef	UNIX_PATH_MAX
#define UNIX_PATH_MAX    108
#endif

int is_unix_socket_name(char *name)
{
 if ( strchr(name,'/') != NULL )
	return(1);
 else
	return(0);
}

/* <port> or <hostname>:<port>, works only for ipv4: */
int is_inet_socket_name(char *name)
{
 int	port;
 char	*colon;

 if ( sscanf(name,"%d",&port)==1 )
	return(1);
 else if ( sscanf(name,"*:%d",&port)==1 )
	return(1);
 else if ( (colon=strchr(name,':')) != NULL && sscanf(colon+1,"%d",&port)==1 )
	return(2);
 else
	return(0);
}

int remote_control_port_bind(char *ctrlport,char **rctrlport,
	int allow_auto,int fail_on_existing)
{
 int			sock,port,ret,pid;
 char			buff[256],*colon;
 struct	stat		st;

 if ( ctrlport==NULL )
	return(-1);

 else if ( allow_auto && strcmp(ctrlport,"inet")==0 )
  {	struct	sockaddr_in	inaddr;

	sock=socket(PF_INET,SOCK_STREAM,0);

	if ( sock<0 )
		return(-1);
	
	for ( port=11228,ret=-1 ; port<16384 ; port++ )
	 {	inaddr.sin_family=AF_INET;
		inaddr.sin_port=htons(port);
		inaddr.sin_addr.s_addr=INADDR_ANY;
		ret=bind(sock,(struct sockaddr *)(&inaddr),(socklen_t)(sizeof(inaddr)));
		if ( ! ret )
			break;
	 }
	if ( ret<0 )
	 {	close(sock);
		return(-1);
	 }

	sprintf(buff,"%d",port);
	if ( rctrlport != NULL )	*rctrlport=xstrdup(buff);

	return(sock);
  }

 else if ( sscanf(ctrlport,"%d",&port)==1 || sscanf(ctrlport,"*:%d",&port)==1 )
  {	struct	sockaddr_in	inaddr;

 	sock=socket(PF_INET,SOCK_STREAM,0);

	if ( sock<0 )
		return(-1);

	inaddr.sin_family=AF_INET;
	inaddr.sin_port=htons(port);
	inaddr.sin_addr.s_addr=INADDR_ANY;

	ret=bind(sock,(struct sockaddr *)(&inaddr),(socklen_t)(sizeof(inaddr)));
	if ( ret<0 )
	 {	close(sock);
		return(-1);
	 }

	sprintf(buff,"%d",port);
	if ( rctrlport != NULL )	*rctrlport=xstrdup(buff);

	return(sock);
  }

 else if ( (colon=strchr(ctrlport,':')) != NULL && sscanf(colon+1,"%d",&port)==1 )
  {	struct	sockaddr_in	inaddr;

	char		*hostname;
	struct  hostent *peer;

	hostname=xmalloc((size_t)(colon-ctrlport)+1);
	memcpy(hostname,ctrlport,(colon-ctrlport));
	hostname[(colon-ctrlport)]=0;
	peer=gethostbyname(hostname);
	free(hostname);

	if ( peer==NULL || peer->h_addrtype != AF_INET )
		return(-1);

	sock=socket(PF_INET,SOCK_STREAM,0);

	if ( sock<0 )
		return(-1);

	inaddr.sin_family=AF_INET;
	inaddr.sin_port=htons(port);
	memcpy(&inaddr.sin_addr.s_addr,peer->h_addr,peer->h_length);

	ret=bind(sock,(struct sockaddr *)(&inaddr),(socklen_t)(sizeof(inaddr)));
	if ( ret<0 )
	 {	close(sock);
		return(-1);
	 }

	if ( rctrlport != NULL )	*rctrlport=xstrdup(ctrlport);

	return(sock);
  }

 else if ( allow_auto && strcmp(ctrlport,"unix")==0 )
  {	struct	sockaddr_un	unaddr;

	sock=socket(PF_UNIX,SOCK_STREAM,0);

	if ( sock<0 )
		return(-1);

	unaddr.sun_family=AF_UNIX;

	pid=(int)getpid();
	sprintf(unaddr.sun_path,"/tmp/pexec.%d.sock",pid);

	if ( ! stat(unaddr.sun_path,&st) && S_ISSOCK(st.st_mode) )
		unlink(unaddr.sun_path);

	ret=bind(sock,(struct sockaddr *)(&unaddr),(socklen_t)(sizeof(unaddr)));
	if ( ret<0 )
	 {	close(sock);
		return(-1);
	 }

	if ( rctrlport != NULL )	*rctrlport=xstrdup(unaddr.sun_path);

	return(sock);
  }

 else if ( is_unix_socket_name(ctrlport) )
  {	struct	sockaddr_un	unaddr;

	sock=socket(PF_UNIX,SOCK_STREAM,0);

	if ( sock<0 )
		return(-1);

	unaddr.sun_family=AF_UNIX;
	strncpy(unaddr.sun_path,ctrlport,UNIX_PATH_MAX-1);
	unaddr.sun_path[UNIX_PATH_MAX-1]=0;

	if ( fail_on_existing && (! stat(unaddr.sun_path,&st)) )
	 {	close(sock);
		return(-1);
	 }
	else if ( (! stat(unaddr.sun_path,&st)) && S_ISSOCK(st.st_mode) )
		unlink(unaddr.sun_path);

	ret=bind(sock,(struct sockaddr *)(&unaddr),(socklen_t)(sizeof(unaddr)));
	if ( ret<0 )
	 {	close(sock);
		return(-1);
	 }

	if ( rctrlport != NULL )	*rctrlport=xstrdup(unaddr.sun_path);

	return(sock);
  }
 else
 	return(-1);

}

int bind_variable_export(char *envvar,char *ctrlport)
{
 if ( ctrlport==NULL )
	return(-1);

 env_export(envvar,ctrlport);

 return(0);
}

int remote_control_port_connect(char *ctrlport)
{
 int			sock,port,ret;
 struct	sockaddr_in	inaddr;
 struct	sockaddr_un	unaddr;
 char			*sc;

 if ( ctrlport==NULL || (! ctrlport[0]) )
	return(-1);

 else if ( is_unix_socket_name(ctrlport) )
  {	if ( (sock=socket(PF_UNIX,SOCK_STREAM,0))<0 )
		return(-1);

	unaddr.sun_family=AF_UNIX;
	strncpy(unaddr.sun_path,ctrlport,UNIX_PATH_MAX-1);
	unaddr.sun_path[UNIX_PATH_MAX-1]=0;

	ret=connect(sock,(struct sockaddr *)(&unaddr),(socklen_t)(sizeof(unaddr)));
	if ( ret<0 )
	 {	close(sock);
		return(-1);
	 }

	return(sock);
  }

 else
  {	inaddr.sin_family=AF_INET;

	sc=strchr(ctrlport,':');
	if ( sc==NULL )
	 {	if ( sscanf(ctrlport,"%d",&port)<1 )
		 	return(-1);
		inaddr.sin_addr.s_addr=htonl(INADDR_LOOPBACK);
		inaddr.sin_port=htons(port);
	 }
	else
	 {	char		*hostname;
		int		len;
		struct	hostent	*peer;

		if ( sscanf(sc+1,"%d",&port)<1 )
			return(-1);

		len=(int)(sc-ctrlport);
		hostname=(char *)xmalloc(len+1);
		memcpy(hostname,ctrlport,len);
		hostname[len]=0;
		peer=gethostbyname(hostname);
		free(hostname);

		if ( peer==NULL)
			return(-1);

		memcpy(&inaddr.sin_addr.s_addr,peer->h_addr,peer->h_length);
		inaddr.sin_port=htons(port);
	 }

	if ( (sock=socket(PF_INET,SOCK_STREAM,0))<0 )
		return(-1);

	ret=connect(sock,(struct sockaddr *)(&inaddr),(socklen_t)(sizeof(inaddr)));
	if ( ret<0 )
	 {	close(sock);
		return(-1);
	 }

	return(sock);
  }	
	
}

/*****************************************************************************/

int remote_status(int sock,FILE *fw)
{
 linebuffer	lrcv;
 char		*line;

 hprintf(sock,"status --all\n");

 linebuffer_reset(&lrcv);
 line=linebuffer_read_line(sock,&lrcv,0);
 if ( line != NULL )
  {	remove_newlines_and_comments(line);
	fprintf(fw,"%s\n",line);
	free(line);
  }
 
 return(0);
}

int remote_lock(int sock,FILE *fw,char *name)
{
 linebuffer	lrcv;
 char		*line;

 hprintf(sock,"lock \"%s\"\n",name);

 linebuffer_reset(&lrcv);
 line=linebuffer_read_line(sock,&lrcv,0);
 if ( line != NULL )
  {	remove_newlines_and_comments(line);
	if ( fw != NULL )	fprintf(fw,"%s\n",line);
	free(line);
  }

 return(0);
}

int remote_unlock(int sock,FILE *fw,char *name)
{
 hprintf(sock,"unlock \"%s\"\n",name);
 return(0);
}

int remote_copy(int sock,char *name,FILE *frin,FILE *fwout)
{
 char	*buff;
 int	blksize,r,fd;
 fd_set	set;

 FD_ZERO(&set);
 fd=fileno(frin);
 FD_SET(fd,&set);
 select(fd+1,&set,NULL,NULL,NULL);

 if ( sock>=0 && name != NULL )	remote_lock(sock,NULL,name);

 blksize=getpagesize();
 buff=(char *)xmalloc(blksize);

 while ( ! feof(frin) )
  {	r=fread(buff,1,blksize,frin);
	if ( r<=0 )
		break;
	else
		fwrite(buff,1,r,fwout);
  }

 free(buff);
 
 if ( sock>=0 && name != NULL )	remote_unlock(sock,NULL,name);

 return(0);
}

int remote_atomic_execute(int is_shell_commands,char *shell,int argc,char **argv,int sock)
{
 int	pid,status;

 pid=fork();

 if ( pid>0 )
  {	waitpid(pid,&status,0);
	return(WEXITSTATUS(status));
  }
 else
  {	
	if ( sock>=0 )	close(sock);
	
	if ( is_shell_commands )
	 {	char	*largv[4];

		largv[0]=shell;	/* /bin/sh	*/
		largv[1]="-c";
		largv[2]=argv[0];
		largv[3]=NULL;
		
		execv(shell,largv);

		fprintf(stderr,_("%s: unable to execute the script '%s'.\n"),
			shell,progbasename);

		exit(2);
	 }
	else
	 {	char	**largv;
		int	i,largc;

		largc=argc;
		largv=(char **)xmalloc(sizeof(char *)*(argc+1));
		for ( i=0 ; i<argc ; i++ )
		 {	largv[i]=argv[i];		}
		largv[argc]=NULL;

		execvp(largv[0],largv);

		fprintf(stderr,_("%s: unable to execute the command '%s'.\n"),
			progbasename,argv[0]);

		exit(2);
	 }
  }

 return(-1);	/* unreachable */
}

int remote_disconnect(int sock)
{
 char	buff[16];
 int	n;

 hprintf(sock,"exit\n");
 n=read(sock,buff,16);

 return(n);
}

/*****************************************************************************/

int pexec_hypervisor_check_load(int loadtype)
{
 double	loadavg[3];

 if ( ! ( 0<=loadtype && loadtype<3 ) )
	return(0);

 getloadavg(loadavg,3);
 return((int)loadavg[loadtype]);
}

int pexec_hypervisor_request_cleanup(hypervisorstatus *hs,client *cl)
{
 int		i;

 for ( i=0 ; i<hs->nrequest ; )
  {	if ( hs->requests[i].cl == cl )
	 {	if ( i<hs->nrequest-1 )
		 {	memmove(hs->requests+i,hs->requests+i+1,sizeof(request)*(hs->nrequest-i-1));	}
		hs->nrequest--;
	 }
	else
		i++;
  }

 if ( hs->nrequest <= 0 )
  {	free(hs->requests);
	hs->requests=NULL;
	hs->nrequest=0; 
  }

 return(0);
}

int pexec_hypervisor_acknowledge_pending(hypervisorstatus *hs)
{
 request	*rw;

 while ( hs->nrequest>0 && hs->nrunning < hs->num_processes && pexec_hypervisor_check_load(hs->use_load) < hs->num_processes )
  {	
	if ( hs->use_fifo )
		rw=&hs->requests[0];
	else
		rw=&hs->requests[hs->nrequest-1];

	hprintf(rw->cl->peer,"acknowledged\n");
	hs->nrunning++;
	rw->cl->status++;

	if ( hs->nrequest>1 && hs->use_fifo )
	 {	memmove(hs->requests,hs->requests+1,sizeof(request)*(hs->nrequest-1));	}

	hs->nrequest--;
	if ( hs->nrequest <= 0 )
	 {	free(hs->requests);
		hs->requests=NULL;
		hs->nrequest=0; 
	 }
  }
 return(0);
}

int pexec_hypervisor_client_parse_command(char *line,hypervisorstatus *hs,client *cl)
{
 char		**cmd;
 request	*rw;
 int		ret;

 cmd=tokenize_spaces_dyn(line);
 if ( cmd==NULL )
	return(0);
 else if ( cmd[0]==NULL )
  {	free(cmd);
	return(0);
  }

 if ( strcmp(cmd[0],"request")==0 )
  {	/*
	if ( hs->nrunning < hs->num_processes && pexec_hypervisor_check_load(hs->use_load) < hs->num_processes )
	 {	hprintf(cl->peer,"acknowledged\n");
		hs->nrunning++;
		cl->status++;
	 }
	else
	 {
	*/
	hs->requests=(request *)xrealloc(hs->requests,sizeof(request)*(hs->nrequest+1));
	rw=&hs->requests[hs->nrequest];
	rw->cl=cl;
	hs->nrequest++;
	pexec_hypervisor_acknowledge_pending(hs);
	/*
	 }
	*/
	ret=0;
  }
 else if ( strcmp(cmd[0],"ready")==0 )
  {	hs->nrunning--;
	cl->status--;
	pexec_hypervisor_acknowledge_pending(hs);
	ret=0;
  }
 else if ( strcmp(cmd[0],"completed")==0 )
  { 	
	ret=0;
  }
 else if ( strcmp(cmd[0],"status")==0 )
  {	double	loadavg[3];

	getloadavg(loadavg,3);

	hprintf(cl->peer,"status num_processes=%d use_load=%d use_fifo=%d nrunning=%d nrequest=%d load=%.2f,%.2f,%.2f\n",
		hs->num_processes,hs->use_load,hs->use_fifo,hs->nrunning,hs->nrequest,
		loadavg[0],loadavg[1],loadavg[2]);

	ret=0;
  }
 else if ( strcmp(cmd[0],"set")==0 )
  {	int	i,w;
	char	*invvar;
	invvar=NULL;
	for ( i=1 ; cmd[i] != NULL ; i++ )
	 {	if ( sscanf(cmd[i],"num_processes=%d",&w)==1 )
		 {	if ( w<0 )	w=0;
			hs->num_processes=w;
		 }
		else if ( sscanf(cmd[i],"use_load=%d",&w)==1 )
		 {	if ( w<0 )	w=-1;
			else if ( w>2 )	w=2;
			hs->use_load=w;
		 }	
		else if ( sscanf(cmd[i],"use_fifo=%d",&w)==1 )
		 {	if ( ! w )	w=0;
			else		w=1;
			hs->use_fifo=w;
		 }
		else
			invvar=cmd[i];
	 }
	if ( invvar != NULL )
	 {	hprintf(cl->peer,"message \"invalid variable alternation '%s'\"\n",invvar);	}
	ret=0;
  }
 else if ( strcmp(cmd[0],"close")==0 )
	ret=1;
 else if ( strcmp(cmd[0],"terminate")==0 )
	ret=-1;
 else
  {	hprintf(cl->peer,"message \"invalid command '%s'\"\n",cmd[0]);
	ret=0;
  }

 free(cmd);

 return(ret);
}

int pexec_hypervisor_main_loop(int sock,int num_processes,int use_load,int use_fifo)
{ 
 client			*cl,*cnext;
 hypervisorstatus	hs;
 fd_set			set;
 int			i,max,buffsize,n;
 char			*buff;
 int			is_in_loop;
 struct sigaction	intact;
 signalinfo		sci;
 int			spipe;
 struct timeval		tv;

 if ( sock<0 )		/* strange */
	return(-1);

 hs.clientlist=NULL;
 hs.num_processes=num_processes;
 hs.use_load=use_load;
 hs.use_fifo=use_fifo;
 hs.requests=NULL;
 hs.nrequest=0;
 hs.nrunning=0;

 buffsize=getpagesize();
 buff=(char *)xmalloc((size_t)buffsize);

 is_in_loop=1;

 if ( pipe(sig_pipe)<0 )
	return(-1);

 intact.sa_handler=sig_act_interrupt;
 sigemptyset(&intact.sa_mask);
 intact.sa_flags=(SA_NOCLDSTOP|SA_RESTART);
 sigaction(SIGINT,&intact,NULL);

 spipe=sig_pipe[0];

 while ( is_in_loop )
  {	
	FD_ZERO(&set);
	FD_SET(sock,&set);
	max=sock;

	for ( cl=hs.clientlist ; cl != NULL ; cl=cl->next )
	 {	FD_SET(cl->peer,&set);
		if ( max<cl->peer )	max=cl->peer;
	 }
	if ( spipe>=0 )
	 {	FD_SET(spipe,&set);
		if ( max<spipe )	max=spipe;
	 }

	if ( use_load>=0 )
	 {	tv.tv_sec=PEXEC_LOAD_CHECK_PERIOD;
		tv.tv_usec=0;
		i=select(max+1,&set,NULL,NULL,&tv);
	 }	
	else
		i=select(max+1,&set,NULL,NULL,NULL);


	if ( i<0 && errno==EINTR )
		continue;

	pexec_hypervisor_acknowledge_pending(&hs);

	/* check clients: */
	for ( cl=hs.clientlist ; cl != NULL ; cl=cnext )
	 {	char	*line;
		int	closepeer;

		cnext=cl->next;

		if ( ! FD_ISSET(cl->peer,&set) )
			continue;

		n=read(cl->peer,buff,buffsize);

		if ( n<=0 )
			closepeer=1;
		else
			closepeer=0;

		if ( n>0 && (! closepeer) )
			linebuffer_concatenate(&cl->lcli,buff,n);
	
		while ( (line=linebuffer_fetch(&cl->lcli)) != NULL && ! closepeer )
		 {	closepeer=pexec_hypervisor_client_parse_command(line,&hs,cl);
			free(line);
			if ( closepeer )	break;
		 };

		if ( closepeer )
		 {	close(cl->peer);
			if ( cl->status>0 )
				hs.nrunning -= cl->status;
			linebuffer_free(&cl->lcli);
			pexec_hypervisor_request_cleanup(&hs,cl);
			list_remove(hs.clientlist,cl);
			free(cl);
			pexec_hypervisor_acknowledge_pending(&hs);
		 }
		if ( closepeer<0 )	/* ask for hypervisor termination too */
			is_in_loop=0;

	 } /* for ( cl=hs.clientlist ; cl != NULL ; cl=cnext ) */

 	/* multiplexed signal handling: */
	if ( FD_ISSET(spipe,&set) )
	 {	if ( ! (read_signalinfo(spipe,&sci)>0) )
		 {	fprintf(stderr,_("%s: fatal error: read_signalinfo() failed.\n"),progbasename);
			exit(1);
		 }
		/* interrupt from keyboard: */
		if ( sci.signal==SIGINT && sci.exitsignal )
			is_in_loop=0;
	 } /* if ( FD_ISSET(spipe,&set) ) */

	if ( FD_ISSET(sock,&set) )
	 {	cl=list_new(client);
		cl->peer=accept(sock,NULL,NULL);
		linebuffer_reset(&cl->lcli);
		list_insert_first(hs.clientlist,cl);
		cl->status=0;
	 } /* if ( FD_ISSET(sock,&set) ) */

  } /* pexec_hypervisor_main_loop(): while ( is_in_loop ) */

 signal(SIGINT,SIG_DFL);

 close(sig_pipe[0]); sig_pipe[0]=-1;
 close(sig_pipe[1]); sig_pipe[1]=-1;

 free(buff);

 return(0);
}

int pexec_hypervisor_stop(int sock)
{
 hprintf(sock,"terminate\n");
 return(0);
}

/*****************************************************************************/

int fprint_parameters(FILE *fw,parameter *params,int nparam)
{
 int	i,j;

 for ( i=0 ; i<nparam ; i++ )
  {	fprintf(fw,"name='%s' is_shell=%d ",params[i].name,params[i].c.is_shell);

	fprintf(fw,"args:");
	for ( j=0 ; j<params[i].c.argc ; j++ )
	 {	fprintf(fw," '%s'",params[i].c.argv[j]);	}

	fprintf(fw,"\n");
  }
 return(0);
}

/*****************************************************************************/

longhelp_entry	pexec_long_help[]=
{
 LONGHELP_OPTIONS,
 { "General options:", NULL },
 { "-h, --help", 
	"Gives general summary about the command line options." },
 { "--long-help", 
	"Gives a detailed list of command line options." } ,
 { "--version", 
	"Gives some version information about the program." } ,
 { "-s, --shell <shell>", 
	"Full path (e.g. /bin/sh) of the shell or interpreter to be used "
	"for script execution." },
 { "-c, --shell-command", 
	"Use the specified shell to interpret the command(s) instead of "
	"direct execution." },
 { "-m, --multiple-command",
	"Allow multiple individual shell command scripts to be executed "
	"in parallel with the variation of the parameters." },
 { "-e, --environment <variable>",
	"Name of an environmental variable which is set to the respective "
	"parameter before each execution." },
 { "-n, --number <number>",
	"The maximal number of processes running simultaneously. The <number> "
	"itself can even be a complex specification of remote hosts (see "
	"documentation for more details)." }, 
 { "-C, --control <port>",
	"The control port of a hypervisor daemon (full path of a UNIX "
	"socket or an INET host specification)." },
 { "-p, --list <list>",
	"The single-argument form of main parameter list." },
 { "-r, --parameters <list>",
	"The multiple-argument form of the main parameter list." },
 { "-f, --listfile <file>",
	"The main parameter list file." },
 { "-w, --column <index>",
	"The column index from where the parameters should be taken if "
	"they are read from a parameter file." },
 { "-t, --complete",
	"Threat the whole line as a single parameter if the "
	"parameters are read from a file." },
 { "-z, --nice",
	"Sets the scheduling priority of pexec and all children "
	"(executed processes) to the priority defined by this nice value." },
 { "--",
	"A marker after which the command to execute begins." } ,

 { "Redirecting standard input, output and error:", NULL },
 { "-i, --input <input>",
	"The (optionally formatted) name of the input file which is used "
	"for redirecting the standard input." },
 { "-o, --output <output>",
	"The (optionally formatted) name of the output "
	"file which is used for redirecting the standard output." },
 { "-u, --error <output>",
	"The (optionally formatted) name of the output "
	"error file, which is used for redirecting the standard error." },
 { "-R, --normal-redirection",
	"Equivalent to specifying --output -, --error - and --input /dev/null." },
 { "-a, --output-format <format>",
	"The format of the final standard output redirection if the output "
	"of all of the processes are gathered into the same file." },
 { "-b, --error-format <format>",
	"The same final redirection format for the standard error." },
 { "-x, --omit-newlines",
	"Disable automatic newlines after the output and error formats." }, 

 { "Execution using remote hosts:", NULL },
 { "-g, --remote-shell <remote_shell>",
	"The name or full path of the remote shell to be used for building "
	"the tunnel between the local and the peer host(s). "
	"Default: ``/usr/bin/ssh''." },
 { "-P, --pexec <pexec>",
	"The full path of the pexec program on the remote hosts. "
	"If this option is omitted, pexec tries to figure out from the "
	"invoking syntax and/or the current path." },
 { "-T, --tunnel",
	"Internal use only (pexec will start in tunnel daemon mode)." },

 { "Remote control, mutual exclusions and atomic command execution:", NULL },
 { "-y, --bind <port>",
	"This option lets pexec to be remote controlled via INET or "
	"UNIX domain sockets." },
 { "-E, --pexec-connection-variable <env>",
	"This option overrides the default environment name PEXEC_REMOTE_PORT "
	"to the specified value, which is used by the ``-p|--connect auto'' "
	"combination to determine the control socket with which the running "
	"pexec instance can be controlled." },
 { "-j, --remote",
	"Used to remote control and/or poll the status of other running "
	"instances of pexec." },
 { "-p, --connect <port>",
	"Remote control port to connect to." },
 { "-t, --status",
	"Prints the actual status of the running jobs in a "
	"human-readable form." },
 { "-l, --lock <mutex>",
	"Locks the specified mutex (if the mutex is not locked by "
	"someone else, otherwise it will block until the mutex is released)." },
 { "-u, --unlock <mutex>",
	"Unlocks the specified mutex." },
 { "-m, --mutex <mutex>",
	"Name of the mutex." },
 { "-d, --dump <filename>",
	"Dump the content of the given file to standard output, "
	"if ``-m|--mutex'' is given, this will be atomic." },
 { "-s, --save <filename>",
	"Save the content of standard input to the given file, "
	"if ``-m|--mutex'' is given, this will be atomic." },
 { "-a, --atomic <command>",
	"Execute the given command. If ``-m|--mutex'' is given, the "
	"exectution is going to be atomic with respect to that mutex." },

 { "Hypervisor mode:",NULL },
 { "-H, --hypervisor",
	"Starts pexec in hypervisor mode." },
 { "-C, --control <port>",
	"The control port used by the hypervisor." },
 { "-l, --load <window>",
	"Use load also to limit the number of simultaneous processes "
	"with the specified load average interval (0, 1 or 2, or 1min, "
	"5min or 15min, respectively)." },
 { "-f, --fifo",
	"First in first out queue processing." },
 { "-s, --lifo",
	"Last in first out (stack) queue processing (default)." },

 { "Logging:",NULL },
 { "-L, --log <file>",
	"The name of the log file." },
 { "-W, --log-level <level>",
	"The logging level." },
 { "-V, --verbose",
	"Increase the log level by one." }, 

 { NULL, NULL }
};
 

int fprint_pexec_long_help(FILE *fw)
{
 fprintf(fw,
"Usage:\tpexec [options] [-c|-m] [--] command [arguments] | 'compound command'\n"
"Execute commands or shell scripts in parallel on a single host or\n"
"on remote hosts using a remote shell.\n\n");

 longhelp_fprint(fw,pexec_long_help,0,-1);

 fprintf(fw,"\n");
 fprintf(fw,"Report bugs to <apal@szofi.elte.hu>\n");

 return(0);
}

/* all flags: -abcdefghijklmnopqrstuvwxyz-ABCDEFGHIJKLMNOPQRSTUVWXYZ- */
/* available: ----d-----------------v-----AB-DEFG-IJK--NO-QRS-U--XYZ- */

int fprint_pexec_usage(FILE *fw)
{
 fprintf(fw,
"Usage:\tpexec"
"\t[options] [--] command [arguments]\n"
"\t\t[options] -c [--] 'compound command'\n"
"\t\t[options] -m [--] 'compound command 1' 'compound command 2' ...\n");
 fprintf(fw,
"General Options:\n"
"\t[-h|--help|--long-help] [--version]\n"
"\t[-c|--shell-command] [-m|--multiple-command] [-s|--shell <shell>]\n"
"\t[-e|--environment <environmental_variable_name>]\n");
 fprintf(fw,
"\t[-p|--list <list_of_parametes> [-p ...] [-p ...]]\n"
"\t[-r|--parameters <list_of_parametes> {--|-<option} [-r ...] [-r ...]]\n"
"\t[-f|--listfile <parameter_file> [-w|--column <column>|-t|--complete]]\n"
"\t[-n|--number auto|<num>|managed|ncpu [-C|--control {<port>|</path>}]]\n"
"\t[-l|--load|--use-load <load>] [-z|--nice <nice>]\n");
 fprintf(fw,
"Redirecting input and output:\n"
"\t[-i|--input <format_for_standard_input_file>]\n"
"\t[-o|--output <format_for_standard_output_file>]\n"
"\t[-u|--error <format_for_standard_error_file>]\n"
"\t[-a|--output-format <format_for_stdout_redirection> [-x]]\n"
"\t[-b|--error-format <format_for_stderr_redirection> [-x|--omit-newlines]]\n"
"\t[-R|--normal-redirection]\n");
 fprintf(fw,
"Parallelization using remote hosts:\n"
"\t[-g|--remote-shell \"<remote_shell [options]>\"] [-q|--timeout <sec>]\n"
"\t-n|--number [<host>:]{auto|<num>|managed|ncpu}[,...],[auto|<num>|...]\n"
"\t[-P|--pexec <full_pexec_path_on_the_remote_host(s)>] [-k|--local-files]\n");
 fprintf(fw,
"Running as a tunnel daemon (only for internal use, see also the manual):\n"
"\t-T|--tunnel [-z|--nice <nice>]\n"
"\t[-n|--number auto|<num>|managed|ncpu [-C|--control {<port>|</path>}]]\n");
 fprintf(fw,
"Remote control, mutual exclusions and atomic command execution:\n"
"\t[-y|--bind inet|unix|<port>|/<path>]\n"
"\t-j|--remote [-p|--connect auto|/<path>|[host:]<port>] [-t|--status]\n"
"\t[-E|--pexec-connection-variable <env_variable_name>]\n"
"\t[{-l|--lock|--mutex-lock | -u|--unlock|--mutex-unlock} <name>]\n"
"\t[-m|--mutex <name> {-d|--dump | -s|--save} <filename> ]\n"
"\t[-m|--mutex <name> -a|--atomic [-c|--shell-command] [--] command [...]]\n");
 fprintf(fw,
"Hypervisor mode and operations:\n"
"\t-H|--hypervisor [-C|--control {<port>|</path>}] [start|stop]\n"
"\t[-n|--number auto|<num>] [-l|--load|--use-load <load>]\n"
"\t[-f|--fifo | -s|--lifo|--stack]\n"
"\t(default hypervisor socket: %s)\n",PEXEC_DEFAULT_HYPERVISOR_SOCKET);
 fprintf(fw,
"Logging (for normal and hypervisor modes):\n"
"\t[-L|--log <log_file>] [-W|--log-level <log_level> | -V|--verbose [...]]\n"
"Command specifications:\n"	
"\t[--] { command [args] | 'compound' | 'compound 1' ['compound 2'...] }\n");
 fprintf(fw,
"For more information, see --long-help or the full texinfo documentation.\n"
"Examples can be found in the ``Examples'' section of the documentation.\n");
 return(0);
}

int fprint_version(FILE *fw)
{
 fprintf(fw,"%s %s (%s)\n",
	progbasename,PEXEC_VERSION,PEXEC_LAST_MODIFICATION);

 fprintf(fw,"Copyright (C) 2007, 2008-2009; Pal, Andras <apal@szofi.elte.hu>\n\n");

 fprintf(fw,
"This is free software. You may redistribute copies of it under the terms of\n"
"the GNU General Public License <http://www.gnu.org/licenses/gpl.html>.\n"
"There is NO WARRANTY, to the extent permitted by law. \n\n");

 fprintf(fw,
"This software was written by Andras Pal. The core part was written \n"
"while working for the Hungarian-made Automated Telescope (HAT) project \n"
"to make the data processing more easier and therefore find many-many \n"
"extrasolar planets. See more information about this project: \n"
"http://hatnet.hu. Another internal libraries (e.g. numhash.[ch]) were \n"
"primarily written for other projects.\n");

 return(0);
}

int fprint_err_invarg0(char *arg)
{
 fprintf(stderr,_("%s: error: invalid command line argument '%s'.\n"),progbasename,arg);
 return(0);
}
int fprint_err_invarg1(char *arg)
{
 fprintf(stderr,_("%s: error: invalid or missing argument near '%s'.\n"),progbasename,arg);
 return(0);
}
int fprint_err_invarg2(char *arg)
{
 fprintf(stderr,_("%s: error: special command line argument '%s' must be the first in the list.\n"),progbasename,arg);
 return(0);
}

int main(int argc,char *argv[])
{
 int		is_shell_commands,		/* -c			*/
		is_multi_commands,		/* -m			*/
		zeroarg;			/* --			*/
 remotehost	*rhosts;			/* -n			*/
 int		nrhost;
 char		*list,				/* -l			*/
		*listfile,			/* -f			*/
		*logfile;			/* -l			*/
 int		listcolumn;			/* -w			*/
 paralleldata	p;				/* -[ab] -[iou] -n	*/
 logdata	log;

 int		ncmd;

 int		llen;
 int		i,ret,prio;
 parameter	*params;
 int		nparam;
 char		**pnames;
 
 char		*pexec_self,*ctrlport,*ctrlenv,*hypcport;
 int		timeout,use_load;

 int		mode;				/* pexec execution mode	*/

 is_shell_commands=0;
 is_multi_commands=0;
 zeroarg=-1;
 list=NULL;
 pnames=NULL;
 listfile=logfile=NULL;
 llen=0;

 rhosts=NULL;
 nrhost=PEXEC_MNP_AUTO;

 p.envvarname=NULL;
 p.out=p.err=p.in=NULL;
 p.shell=PEXEC_DEFAULT_SHELL;
 p.rshcmd=NULL;

 prio=0;

 p.formatout=NULL;
 p.formaterr=NULL;
 p.omit_newlines=0;

 log.loglevel=-1;
 log.fwlog=NULL;
 
 listcolumn=0;

 progbasename=strrchr(argv[0],'/');
 if ( progbasename != NULL )	progbasename++;
 else				progbasename=argv[0];

 mode=PEXEC_MODE_DEFAULT;

 if ( strchr(argv[0],'/') != NULL )
	pexec_self=argv[0];
 else
  {	/* pexec_self="/home/apal/bin/pexec"; */
	pexec_self=progbasename;
  }

 timeout=60;

 ctrlport=NULL;
 ctrlenv=PEXEC_DEFAULT_ENVVARIABLE;

 hypcport=NULL;

 use_load=-1;

 for ( i=1 ; i<argc ; i++ )
  {	/* -h, --help: give some general usage and exit successfully */
	if ( strcmp(argv[i],"-h")==0 || strcmp(argv[i],"--short-help")==0 || strcmp(argv[i],"--help")==0 || strcmp(argv[i],"--help-short")==0 )
	 {	fprint_pexec_usage(stdout);
		return(0);
	 }
	/* --long-help: give some more detailed usage and exit successfully */
	if ( strcmp(argv[i],"--long-help")==0 || strcmp(argv[i],"--help-long")==0 )
	 {	fprint_pexec_long_help(stdout);
		return(0);
	 }
	/* --version: give some version information about pexec: */
	else if ( strcmp(argv[i],"--version")==0 )
	 {	fprint_version(stdout);
		return(0);
	 }
	/* -c, --shell-command: execute commands using a shell */
	else if ( strcmp(argv[i],"-c")==0 || strcmp(argv[i],"--shell-command")==0 )
		is_shell_commands=1;
	/* -m, --multiple-command: execute multiple commands, implies -c */
	else if ( strcmp(argv[i],"-m")==0 || strcmp(argv[i],"--multiple-command")==0 )
		is_multi_commands=1;
	/* -e, --environment: name of an environmental variable to pass the current parameter */
	else if ( strcmp(argv[i],"-e")==0 || strcmp(argv[i],"--environment")==0 || strcmp(argv[i],"--setenv")==0 )
	 {	if ( (++i)>=argc )
		 {	fprint_err_invarg1(argv[i-1]);return(1);	}
		p.envvarname=argv[i];
	 }
	/* -s, --shell: full path of the shell which should understand '-c' */
	else if ( strcmp(argv[i],"-s")==0 || strcmp(argv[i],"--shell")==0 )
	 {	if ( (++i)>=argc )
		 {	fprint_err_invarg1(argv[i-1]);return(1);	}
		p.shell=argv[i];
	 }
	/* -g, --remote-shell: name (and optional arguments) for the remote shell */
	else if ( strcmp(argv[i],"-g")==0 || strcmp(argv[i],"--remote-shell")==0 )
 	 {	if ( (i++)>=argc )
		 {	fprint_err_invarg1(argv[i-1]);return(1);	}
		p.rshcmd=xstrdup(argv[i]);
	 }
	/* -p, --list: list of parameters: space separated single arg */
	else if ( strcmp(argv[i],"-p")==0 || strcmp(argv[i],"--list")==0 )
	 {	if ( (++i)>=argc )
		 {	fprint_err_invarg1(argv[i-1]);return(1);	}
		else if ( list==NULL )
		 {	list=xstrdup(argv[i]);
			llen=strlen(list);
		 }
		else
		 {	int	ilen;
			ilen=strlen(argv[i]);
			list=(char *)xrealloc(list,llen+ilen+2);
			strcpy(list+llen," ");llen++;
			strcpy(list+llen,argv[i]);llen+=ilen;
		 }
	 }
	/* -r, --parameters: list of parameters: multiple arguments */
	else if ( strcmp(argv[i],"-r")==0 || strcmp(argv[i],"--parameters")==0 )
	 {	i++;
		for ( ; i<argc && argv[i][0] != '-' ; i++ )
		 {	if ( list==NULL )
			 {	list=xstrdup(argv[i]);
				llen=strlen(list);
			 }
			else
			 {	int	ilen;
				ilen=strlen(argv[i]);
				list=(char *)xrealloc(list,llen+ilen+2);
				strcpy(list+llen," ");llen++;
				strcpy(list+llen,argv[i]);llen+=ilen;
			 }
		 }
		i--;
	 }
	/* -f, --listfile: parameter list file */
	else if ( strcmp(argv[i],"-f")==0 || strcmp(argv[i],"--listfile")==0 )
	 {	if ( (++i)>=argc )
		 {	fprint_err_invarg1(argv[i-1]);return(1);	}
		listfile=argv[i];
	 }
	/* -w, --column: parameter column in the list file */
	else if ( strcmp(argv[i],"-w")==0 || strcmp(argv[i],"--column")==0 )
	 {	if ( (++i)>=argc )
		 {	fprint_err_invarg1(argv[i-1]);return(1);	}
		else if ( sscanf(argv[i],"%d",&listcolumn)<1 ||
		listcolumn<=0 )
		 {	fprint_err_invarg1(argv[i-1]);
			return(1);
		 }
		listcolumn--;
	 }
	/* -t, --complete: use the whole line of the list file as a param */
	else if ( strcmp(argv[i],"-t")==0 || strcmp(argv[i],"--complete")==0 || strcmp(argv[i],"--complete-line")==0 )
		listcolumn=-1;
	/* -n, --number: number of processes to execute simultaneously */
	else if ( strcmp(argv[i],"-n")==0 || strcmp(argv[i],"--number")==0 )
	 {	if ( (++i)>=argc )
		 {	fprint_err_invarg1(argv[i-1]);return(1);	}
		else if ( strcmp(argv[i],PEXEC_ISTR_AUTO)==0 )
		 {	rhosts=NULL;
			nrhost=PEXEC_MNP_AUTO;
		 }
		else if ( strcmp(argv[i],PEXEC_ISTR_MANAGED)==0 )
		 {	rhosts=NULL;
			nrhost=PEXEC_MNP_MANAGED;
		 }
		else if ( strcmp(argv[i],PEXEC_ISTR_NCPU)==0 )
		 {	rhosts=NULL;
			nrhost=PEXEC_MNP_NCPU;
		 }
		else if ( parse_host_data(argv[i],&rhosts,&nrhost) )
		 {	fprint_err_invarg1(argv[i-1]);
			return(1);
		 }
	 }
	/* -M, --managed: equivalent to `-n|--number managed`: */
	else if ( strcmp(argv[i],"-M")==0 || strcmp(argv[i],"--managed")==0 )
	 {	rhosts=NULL;
		nrhost=PEXEC_MNP_MANAGED;
	 }
	/* -l, --load: use load */
	else if ( strcmp(argv[i],"-l")==0 || strcmp(argv[i],"--load")==0 || strcmp(argv[i],"--use-load")==0 )
	 {	if ( (++i)>=argc )
		 {	fprint_err_invarg1(argv[i-1]);return(1);	}
		else if ( strcmp(argv[i],"0")==0 || strcmp(argv[i], "1m")==0 || strcmp(argv[i], "1min")==0 )
			use_load=0;
		else if ( strcmp(argv[i],"1")==0 || strcmp(argv[i], "5m")==0 || strcmp(argv[i], "5min")==0 )
			use_load=1;
		else if ( strcmp(argv[i],"2")==0 || strcmp(argv[i],"15m")==0 || strcmp(argv[i],"15min")==0 )
			use_load=2;
		else
		 {	fprint_err_invarg1(argv[i-1]);
			return(1);
		 }
	 }
	/* -z, --nice: priority of pexec and the (initial) priority of the processes */
	else if ( strcmp(argv[i],"-z")==0 || strcmp(argv[i],"--nice")==0 )
	 {	if ( (++i)>=argc )
		 {	fprint_err_invarg1(argv[i-1]);return(1);	}
		else if ( sscanf(argv[i],"%d",&prio)<1 )
		 {	fprint_err_invarg1(argv[i-1]);
			return(1);
		 }
		prio+=128;
	 }
	/* -q, --timeout: timeout in seconds for remote shell execution */
	else if ( strcmp(argv[i],"-q")==0 || strcmp(argv[i],"--timeout")==0 )
	 {	if ( (++i)>=argc )
		 {	fprint_err_invarg1(argv[i-1]);return(1);	}
		else if ( sscanf(argv[i],"%d",&timeout)<1 )
		 {	fprint_err_invarg1(argv[i-1]);return(1);	}
		if ( timeout<0 )
			timeout=0;
	 }
	/* -P, --pexec: full path to pexec on the remote hosts: */
	else if ( strcmp(argv[i],"-P")==0 || strcmp(argv[i],"--pexec")==0 )
	 {	if ( (++i)>=argc )
		 {	fprint_err_invarg1(argv[i-1]);return(1);	}
		pexec_self=argv[i];
	 }
	/* -E, --pexec-connection-variable <environment variable name>: */
	else if ( strcmp(argv[i],"-E")==0 || strcmp(argv[i],"--pexec-connection-variable")==0 )
	 {	if ( (++i)>=argc )
		 {	fprint_err_invarg1(argv[i-1]);return(1);	}
		ctrlenv=argv[i];
	 }
	/* -i, --input: input file name(s) */
	else if ( strcmp(argv[i],"-i")==0 || strcmp(argv[i],"--input")==0 )
	 {	if ( (++i)>=argc )
		 {	fprint_err_invarg1(argv[i-1]);return(1);	}
		p.in=argv[i];
	 }
	/* -o, --output: output file name(s) or '-' or '-[12]' */
	else if ( strcmp(argv[i],"-o")==0 || strcmp(argv[i],"--output")==0 )
	 {	if ( (++i)>=argc )
		 {	fprint_err_invarg1(argv[i-1]);return(1);	}
		p.out=argv[i];
	 }
	/* -u, --error: error file name(s) or '-' or '-[12]' */
	else if ( strcmp(argv[i],"-u")==0 || strcmp(argv[i],"--error")==0 || strcmp(argv[i],"--output-error")==0 )
	 {	if ( (++i)>=argc )
		 {	fprint_err_invarg1(argv[i-1]);return(1);	}
		p.err=argv[i];
	 }
	/* -R, --normal-redirection: */
	else if ( strcmp(argv[i],"-R")==0 || strcmp(argv[i],"--normal-redirection")==0 )
	 {	p.in=NULL;
		p.out="-";
		p.err="-";
	 }
	/* -a, --output-format: post-format of std output lines */
	else if ( strcmp(argv[i],"-a")==0 || strcmp(argv[i],"--output-format")==0 )
	 {	if ( (++i)>=argc )
		 {	fprint_err_invarg1(argv[i-1]);return(1);	}
		p.formatout=argv[i];
	 }
	/* -b, --error-format: post-format of error lines */
	else if ( strcmp(argv[i],"-b")==0 || strcmp(argv[i],"--error-format")==0 )
	 {	if ( (++i)>=argc )
		 {	fprint_err_invarg1(argv[i-1]);return(1);	}
		p.formaterr=argv[i];
	 }
	/* -x, --omit-newlines: omit newlines from the end of post-format */
	else if ( strcmp(argv[i],"-x")==0 || strcmp(argv[i],"--omit-newlines")==0 )
		p.omit_newlines=1;
	/* -y, --bind: control socket binding port or path */
	else if ( strcmp(argv[i],"-y")==0 || strcmp(argv[i],"--bind")==0 )
	 {	if ( (++i)>=argc )
		 {	fprint_err_invarg1(argv[i-1]);return(1);	}
		ctrlport=argv[i];
	 }
	/* -C, --control: hypervisor socket binding port or path */
	else if ( strcmp(argv[i],"-C")==0 || strcmp(argv[i],"--control")==0 )
	 {	if ( (++i)>=argc )
		 {	fprint_err_invarg1(argv[i-1]);return(1);	}
		hypcport=argv[i];
	 }
	/* -L, --log: name of output log file or '-' or '-[12]' */
	else if ( strcmp(argv[i],"-L")==0 || strcmp(argv[i],"--log")==0 )
	 {	if ( (++i)>=argc )
		 {	fprint_err_invarg1(argv[i-1]);return(1);	}
		logfile=argv[i];
	 }
	/* -V, --verbose: increase log level by one: */
	else if ( strcmp(argv[i],"-V")==0 || strcmp(argv[i],"--verbose")==0 )
		log.loglevel++;
	/* -W, --log-level: log level, some non-negative integer */
	else if ( strcmp(argv[i],"-W")==0 || strcmp(argv[i],"--log-level")==0 )
	 {	if ( (++i)>=argc )
		 {	fprint_err_invarg1(argv[i-1]);return(1);	}
		else if ( sscanf(argv[i],"%d",&log.loglevel)<1 || log.loglevel<0 )
		 {	fprint_err_invarg1(argv[i-1]);
			return(1);
		 }
	 } 
	/* -T, --tunnel: run pexec as a tunnel daemon: */
	else if ( strcmp(argv[i],"-T")==0 || strcmp(argv[i],"--tunnel")==0 )
	 {	if ( i>1 )
		 {	fprint_err_invarg2(argv[i]);
			return(1);
		 }
		zeroarg=i+1;
		mode=PEXEC_MODE_DAEMON;
		break;
	 }
	/* -j, --remote: remote control */
	else if ( strcmp(argv[i],"-j")==0 || strcmp(argv[i],"--remote")==0 )
	 {	if ( i>1 )
		 {	fprint_err_invarg2(argv[i]);
			return(1);
		 }
		zeroarg=i+1;
		mode=PEXEC_MODE_REMOTECONTROL;
		break;
	 }
	/* -H, --hypervisor: hypervisor mode */
	else if ( strcmp(argv[i],"-H")==0 || strcmp(argv[i],"--hypervisor")==0 )
	 {	if ( i>1 )
		 {	fprint_err_invarg2(argv[i]);
			return(1);
		 }
		zeroarg=i+1;
		mode=PEXEC_MODE_HYPERVISOR;
		break;
	 }
	/* --: option delimiter */
	else if ( strcmp(argv[i],"--")==0 )
	 {	zeroarg=i+1;
		if ( zeroarg>=argc )
			zeroarg=-1;
		break;
	 }
	/* invalid argument */
	else if ( argv[i][0]=='-' )
	 {	fprint_err_invarg0(argv[i]);	
		return(1);
	 }
	/* the first real command */
	else
	 {	zeroarg=i;
		break;
	 }
  }

 /* command has not been supplied, exiting. */
 if ( zeroarg<0 )
	return(0);

 /* The default mode for pexec: execution of the command(s): */
 if ( mode == PEXEC_MODE_DEFAULT )
  {	remoteshell	*rshells,*rs;
	int		nrshell,r;
	int		status,sock,hsck;
	char		*pctrlport;
	int		is_rhosts_defined,no_touch_std;

	if ( rhosts==NULL )
	 {	rhosts=(remotehost *)xmalloc(sizeof(remotehost)*1);
		rhosts[0].hostspec=NULL;
		rhosts[0].num_processes=nrhost;
		nrhost=1;
		is_rhosts_defined=0;
	 }
	else
		is_rhosts_defined=!0;

	if ( p.rshcmd==NULL )
		p.rshcmd=xstrdup(PEXEC_DEFAULT_RSH);

	p.rshargs=tokenize_spaces_dyn(p.rshcmd);
	p.rsh=p.rshargs[0];

	/* only one parameter source can be defined: */
	if ( list != NULL && listfile != NULL )
	 {	fprintf(stderr,_("%s: error: both parameter list and list file are defined.\n"),progbasename);
		return(1);
	 }
	 /* parameters are taken from the command line argument(s): */
	else if ( list != NULL )
	 {	pnames=tokenize_spaces_dyn(list);
		 for ( nparam=0 ; pnames != NULL && pnames[nparam] != NULL ; )
		nparam++;
		params=(parameter *)xmalloc(sizeof(parameter)*nparam);
		for ( i=0 ; i<nparam ; i++ )
		 {	params[i].name=pnames[i];		}
	 }
	/* parameters are taken from a file: */
	else if ( listfile != NULL )
	 {	FILE	*fr;
		char	*line,**tokens;
		int	ntoken;

		if ( strcmp(listfile,"-")==0 )
			fr=stdin;
		else
			fr=fopen(listfile,"rb");
		if ( fr==NULL )
		 {	fprintf(stderr,_("%s: error: unable to open list file.\n"),progbasename);
			return(1);
		 }
		params=NULL;
		nparam=0;
		while ( ! feof(fr) )
		 {	line=freadline(fr);
			if ( line==NULL )
				break;
			remove_newlines_and_comments(line);
			if ( ! ( strlen(line)>0 ) )
			 {	free(line);
				continue;
			 }

			else if ( listcolumn<0 )
			 {	params=(parameter *)xrealloc(params,sizeof(parameter)*(nparam+1));
				params[nparam].name=xstrdup(line);
				free(line);
				nparam++;
			 }
			else
			 {	tokens=tokenize_spaces_dyn(line);
				if ( tokens==NULL )
				 {	free(line);
					continue;
				 }
				for ( ntoken=0 ; tokens[ntoken] != NULL ; )
					ntoken++;
				if ( listcolumn>=ntoken )
				 {	free(tokens);
					free(line);
					continue;
				 }
				params=(parameter *)xrealloc(params,sizeof(parameter)*(nparam+1));
				params[nparam].name=xstrdup(tokens[listcolumn]);
				free(tokens);
				free(line);
				nparam++;
			 }
		 }
		/*
		fprintf(stderr,"nparam=%d\n",nparam);
		for ( i=0 ; i<nparam ; i++ )
		 {	fprintf(stderr,"%s\n",params[i].name);		}
		*/
		if ( fileno(fr) != fileno(stdin) )
			fclose(fr);
	 }
	/* no parameters: */
	else
	 {	params=NULL;
		nparam=0;
	 }	

	/* some initialization: */
	for ( i=0 ; i<nparam ; i++ )
	 {	params[i].c.is_shell=0;
		params[i].c.argv=NULL;
		params[i].c.argc=0;
	 }

	if ( is_multi_commands )
		ncmd=argc-zeroarg;
	else
		ncmd=1;

	no_touch_std=0;

	if ( nparam>0 && is_multi_commands && nparam != ncmd )
	 {	fprintf(stderr,_("%s: error: number of parameters and commands mismatch.\n"),progbasename);
		return(1);
	 }
	else if ( is_multi_commands && nparam>0 )
	 {	for ( i=0 ; i<nparam ; i++ )
		 {	params[i].c.is_shell=1;
			params[i].c.argc=1;
			params[i].c.argv=argv+zeroarg+i;
		 }
	 }
	else if ( is_multi_commands )
	 {	nparam=ncmd;
		params=(parameter *)xmalloc(sizeof(parameter)*nparam);
		for ( i=0 ; i<ncmd ; i++ )
		 {	params[i].name=NULL;
		 	params[i].c.is_shell=1;
			params[i].c.argc=1;
			params[i].c.argv=argv+zeroarg+i;
		 }
	 }
	else if ( nparam>0 && is_shell_commands )
	 {	for ( i=0 ; i<nparam ; i++ )
		 {	params[i].c.is_shell=1;
			params[i].c.argc=argc-zeroarg;
			params[i].c.argv=argv+zeroarg;
		 }
	 }
	else if ( nparam>0 )
	 {	for ( i=0 ; i<nparam ; i++ )
		 {	params[i].c.is_shell=0;
			params[i].c.argc=argc-zeroarg;
			params[i].c.argv=argv+zeroarg;
		 }
	 }
	else if ( ! is_rhosts_defined )
	 {	nparam=1;
		params=(parameter *)xmalloc(sizeof(parameter)*nparam);
		params[0].name=NULL;
		params[0].c.is_shell=is_shell_commands;
		params[0].c.argv=argv+zeroarg;
		params[0].c.argc=argc-zeroarg;
		no_touch_std=!0;
	 }
	else
	 {	nparam=nrhost;
		params=(parameter *)xmalloc(sizeof(parameter)*nparam);
		for ( i=0 ; i<nparam ; i++ )
	 	 {	params[i].name=rhosts[i].hostspec;
			params[i].c.is_shell=-1;
			params[i].c.argv=argv+zeroarg;
			params[i].c.argc=argc-zeroarg;
		 }
		rhosts[0].hostspec=NULL;
		rhosts[0].num_processes=nparam;
		nrhost=1;
		no_touch_std=!0;
	 }

	for ( i=0 ; i<nparam ; i++ )
	 { 	params[i].status=0;
		params[i].no_touch_std=0;
		params[i].id=i;
	 }

	if ( ( no_touch_std || nparam<=1 ) && p.in==NULL && p.out==NULL && p.err==NULL )
	 {	for ( i=0 ; i<nparam ; i++ )
		 {	params[i].no_touch_std=1;		}
	 }
 
	/*
	 fprint_parameters(stderr,params,nparam);
	*/ /* debug: print execution parameter list */

	if ( p.out != NULL && (! format_check_if_formatted(p.out,"skd")) )
	 {	if ( strcmp(p.out,"-")==0 || strcmp(p.out,"-1")==0 )
			p.fwout=stdout;
		else if ( strcmp(p.out,"-2")==0 )	
			p.fwout=stderr;
		else
			p.fwout=fopen(p.out,"wb");
		if ( p.fwout==NULL )
		 {	fprintf(stderr,_("%s: error: unable to create collective output file '%s'.\n"),progbasename,p.out);
			return(1);
		 }
	  }
	else
		p.fwout=NULL;

	if ( p.err != NULL && (! format_check_if_formatted(p.err,"skd")) )
	 {	if ( strcmp(p.err,"-")==0 || strcmp(p.err,"-2")==0 )
			p.fwerr=stderr;
		else if ( strcmp(p.err,"-1")==0 )
			p.fwerr=stdout;
		else
			p.fwerr=fopen(p.err,"wb");
		if ( p.fwerr==NULL )
		 {	fprintf(stderr,_("%s: error: unable to create collective error file '%s'.\n"),progbasename,p.err);
			return(1);
		 }
	 }
	else
		p.fwerr=NULL;

	if ( logfile != NULL && log.loglevel != 0 )
	 {	if ( strcmp(logfile,"-")==0 || strcmp(logfile,"-2")==0 )
			log.fwlog=stderr;
		else if ( strcmp(logfile,"-1")==0 )
			log.fwlog=stdout;
		else
			log.fwlog=fopen(logfile,"wb");
		if ( log.fwlog==NULL )
		 {	fprintf(stderr,_("%s: error: unable to create log file '%s'.\n"),progbasename,logfile);
			return(1);
		 }
		if ( log.loglevel<0 )
			log.loglevel=+1;
	 }
	else if ( log.loglevel>0 )
		log.fwlog=stderr;
	else
	 {	log.fwlog=NULL;
		log.loglevel=0;
	 }

	if ( prio > 0 )
	 {	i=setpriority(PRIO_PROCESS,0,prio-128);
		if ( i<0 )
		 {	fprintf(stderr,_("%s: warning: unable to set scheduling priority to %d.\n"),progbasename,prio-128);	}
	 }

	if ( ctrlport != NULL ) 
	 {	sock=remote_control_port_bind(ctrlport,&pctrlport,1,1);
		if ( sock<0 )
		 {	fprintf(stderr,_("%s: error: unable to create or bind control socket to %s.\n"),progbasename,ctrlport);
			return(1);
		 }
 	 }
	else
	 {	sock=-1;
		pctrlport=NULL;
	 }

	p.fallback_to_die=1;

	nrshell=nrhost;
	rshells=(remoteshell *)xmalloc(sizeof(remoteshell)*nrshell);
	hsck=-1;
	for ( i=0 ; i<nrhost ; i++ )
	 {	
		if ( rhosts[i].hostspec==NULL && (rhosts[i].num_processes==PEXEC_MNP_AUTO || rhosts[i].num_processes==PEXEC_MNP_MANAGED) && hsck<0 )
		 {	if ( hypcport==NULL )
				hypcport=PEXEC_DEFAULT_HYPERVISOR_SOCKET;
			hsck=remote_control_port_connect(hypcport);
			if ( hsck<0 && rhosts[i].num_processes==PEXEC_MNP_MANAGED )
			 {	fprintf(stderr,_("%s: error: unable to connect hypervisor socket '%s'.\n"),progbasename,hypcport);
				return(1);
			 }
			else if ( hsck<0 )
			 {	rhosts[i].num_processes=PEXEC_MNP_NCPU;
				hsck=-1;
			 }
			else
				rhosts[i].num_processes=0;

			/* num_processes can only be PEXEC_MNP_NCPU if negative: */
			if ( rhosts[i].num_processes < 0 )
				rhosts[i].num_processes=get_number_of_cpus();
		 }


		r=remote_shell_init(p.rsh,p.rshargs,pexec_self,ctrlport,timeout,&rhosts[i],&rshells[i],prio);
		if ( r ) 
		 {	fprintf(stderr,_("%s: error: unable to connect or initialize remote shell '%s' and/or pexec daemon '%s' to the host '%s' (reason code=%d).\n"),progbasename,p.rsh,pexec_self,rhosts[i].hostspec,r);
			return(1);
		 }

		/* fprintf(stderr,"[%s:%d]\n",rhosts[i].hostspec,rhosts[i].num_processes); */
	 }

	if ( pctrlport != NULL )
		bind_variable_export(ctrlenv,pctrlport);

	if ( sock>=0 )	listen(sock,256);

	p.log=&log;

	/* well, it's turned out to be necessary: */
	signal(SIGPIPE,SIG_IGN);

	ret=pexec_do_parallelized_execution(&p,params,nparam,rshells,nrshell,sock,hsck);

	for ( i=0,rs=rshells ; i<nrshell ; i++,rs++ )
	 {	/* ignore local "remote shells" or what */
		if ( rs->pid < 0 )
			continue;

		/* send the exit command via the tunnel to stop child pexec: */
		hprintf(rs->fhsend,"exit\n");	
		/* wait for the remote shell process to die: */
		waitpid(rs->pid,&status,0);	
	 }

	if ( rshells != NULL )	free(rshells);
 
	if ( log.fwlog != NULL && fileno(log.fwlog) != 1 && fileno(log.fwlog) != 2 )	
		fclose(log.fwlog);
	if ( p.fwout != NULL && fileno(p.fwout) != 1 && fileno(p.fwout) != 2 )
		fclose(p.fwout);
	if ( p.fwerr != NULL && fileno(p.fwerr) != 1 && fileno(p.fwerr) != 2 )
		fclose(p.fwerr);

	if ( sock>=0 )
	 {	close(sock);
		if ( pctrlport != NULL && is_unix_socket_name(pctrlport) )
			unlink(pctrlport);
	 }
	if ( hsck>=0 )
	 {	hprintf(hsck,"close\n");
		fdwait(hsck,0);
		i=close(hsck);
		/* fprintf(stderr,"pexec: default: hypervisor socket: close()=%d\n",i); */
	 }

	if ( pnames != NULL )		free(pnames);
	if ( list   != NULL )		free(list);

	if ( params != NULL )		free(params);

	if ( pctrlport != NULL )	free(pctrlport);

	free(p.rshcmd);
	free(p.rshargs);

	for ( i=0 ; i<nrhost ; i++ )
	 {	if ( rhosts[i].hostspec != NULL )
			free(rhosts[i].hostspec);
	 }
	free(rhosts);

	return(ret);
  }

 else if ( mode == PEXEC_MODE_DAEMON )
  {	int	fhrecv,fhsend;
	int	num_processes;
	int	sock,hsck;
	char	*pctrlport;

	num_processes=-1;
	ctrlport=NULL;
	prio=0;

	for ( i=zeroarg ; i<argc ; i++ )
	 {	if ( strcmp(argv[i],"-n")==0 || strcmp(argv[i],"--number")==0 )
		 {	if ( (++i)>=argc )
			 {	fprint_err_invarg1(argv[i-1]);return(1);	}
			else if ( strcmp(argv[i],PEXEC_ISTR_AUTO)==0 )
				num_processes=PEXEC_MNP_AUTO;
			else if ( strcmp(argv[i],PEXEC_ISTR_MANAGED)==0 )
				num_processes=PEXEC_MNP_MANAGED;
			else if ( strcmp(argv[i],PEXEC_ISTR_NCPU)==0 )
				num_processes=PEXEC_MNP_NCPU;
			else if ( sscanf(argv[i],"%d",&num_processes)<1 || num_processes <= 0 )
			 {	fprint_err_invarg1(argv[i-1]);
				return(1);
			 }
		 }
		/* -y, --bind: control socket binding port or path */
		else if ( strcmp(argv[i],"-y")==0 || strcmp(argv[i],"--bind")==0 )
		 {	if ( (++i)>=argc )
			 {	fprint_err_invarg1(argv[i-1]);return(1);	}
			ctrlport=argv[i];
		 }
		/* -C, --control: hypervisor socket binding port or path */
		else if ( strcmp(argv[i],"-C")==0 || strcmp(argv[i],"--control")==0 )
		 {	if ( (++i)>=argc )
			 {	fprint_err_invarg1(argv[i-1]);return(1);	}
			hypcport=argv[i];
		 }
		/* -z, --nice: priority of pexec and the (initial) priority of the processes */
		else if ( strcmp(argv[i],"-z")==0 || strcmp(argv[i],"--nice")==0 )
		 {	if ( (++i)>=argc )
			 {	fprint_err_invarg1(argv[i-1]);return(1);	}
			else if ( sscanf(argv[i],"%d",&prio)<1 )
			 {	fprint_err_invarg1(argv[i-1]);
				return(1);
			 }
			prio+=128;
		 }
		else
		 {	fprint_err_invarg0(argv[i]);	
			return(1);
		 }
	 }

	if ( num_processes==PEXEC_MNP_AUTO || num_processes==PEXEC_MNP_MANAGED )
	 {	if ( hypcport==NULL )
			hypcport=PEXEC_DEFAULT_HYPERVISOR_SOCKET;
		hsck=remote_control_port_connect(hypcport);
		if ( hsck<0 && num_processes==PEXEC_MNP_MANAGED )
		 {	fprintf(stderr,_("%s: error: unable to connect hypervisor socket '%s'.\n"),progbasename,hypcport);
			return(1);
		 }
		else if ( hsck<0 )
		 {	num_processes=PEXEC_MNP_NCPU;
			hsck=-1;
		 }
		else
			num_processes=0;
	 }
	else
		hsck=-1;

	/* num_processes can only be PEXEC_MNP_NCPU if negative: */
	if ( num_processes < 0 )
		num_processes=get_number_of_cpus();

	if ( ctrlport != NULL ) 
	 {	sock=remote_control_port_bind(ctrlport,&pctrlport,1,1);
		if ( sock<0 )
		 {	fprintf(stderr,_("%s: error: unable to create or bind control socket to %s.\n"),progbasename,ctrlport);
			return(1);
		 }
 	 }
	else
	 {	sock=-1;
		pctrlport=NULL;
	 }

	fhrecv=fileno(stdin);	/* always 0 */
	fhsend=fileno(stdout);	/* always 1 */

	if ( pctrlport != NULL )
		bind_variable_export(ctrlenv,pctrlport);

	if ( sock>=0 )	listen(sock,256);

	if ( prio > 0 )
	 {	i=setpriority(PRIO_PROCESS,0,prio-128);
		/*
		if ( i<0 )
		 {	fprintf(stderr,_("%s: warning: unable to set scheduling priority to %d.\n"),progbasename,prio-128);	}
		*/
	 }

	/* well, it's turned out to be necessary: */
	signal(SIGPIPE,SIG_IGN);

	pexec_daemon_main_loop(fhrecv,fhsend,num_processes,sock,hsck);

	if ( sock>=0 )
	 {	close(sock);
		if ( pctrlport != NULL && is_unix_socket_name(pctrlport) )
			unlink(pctrlport);
	 }
	if ( hsck>=0 )
	 {	hprintf(hsck,"close\n");
		fdwait(hsck,0);
		i=close(hsck);
	 }

	if ( pctrlport != NULL )	free(pctrlport);

	return(0);
  }

 else if ( mode == PEXEC_MODE_REMOTECONTROL )
  {	char	*ctrlport;
	int	sock;
	int	task;
	char	*file,*name;
	FILE	*fr,*fw;
	int	azeroarg;

	ctrlport=NULL;

	task=PEXEC_REMOTE_STATUS;	/* default */
	name=NULL;
	file=NULL;
	is_shell_commands=0;
	azeroarg=-1;

	for ( i=zeroarg ; i<argc ; i++ )
	 {	if ( strcmp(argv[i],"-p")==0 || strcmp(argv[i],"--connect")==0 )
		 {	if ( (++i)>=argc )
			 {	fprint_err_invarg1(argv[i-1]);return(1);	}
			ctrlport=argv[i];
		 }
		else if ( strcmp(argv[i],"-t")==0 || strcmp(argv[i],"--status")==0 )
			task=PEXEC_REMOTE_STATUS;
		else if ( strcmp(argv[i],"-l")==0 || strcmp(argv[i],"--lock")==0 || strcmp(argv[i],"--mutex-lock")==0 )
		 {	if ( (++i)>=argc )
			 {	fprint_err_invarg1(argv[i-1]);return(1);	}
			task=PEXEC_REMOTE_LOCK;
			name=argv[i];
		 }
		else if ( strcmp(argv[i],"-u")==0 || strcmp(argv[i],"--unlock")==0 || strcmp(argv[i],"--mutex-unlock")==0 )
		 {	if ( (++i)>=argc )
			 {	fprint_err_invarg1(argv[i-1]);return(1);	}
			task=PEXEC_REMOTE_UNLOCK;
			name=argv[i];
		 }
		else if ( strcmp(argv[i],"-m")==0 || strcmp(argv[i],"--mutex")==0 )
		 {	if ( (++i)>=argc )
			 {	fprint_err_invarg1(argv[i-1]);return(1);	}
			name=argv[i];
		 }
		else if ( strcmp(argv[i],"-d")==0 || strcmp(argv[i],"--dump")==0 )
		 {	if ( (++i)>=argc )
			 {	fprint_err_invarg1(argv[i-1]);return(1);	}
			task=PEXEC_REMOTE_DUMP;
			file=argv[i];
		 }
		else if ( strcmp(argv[i],"-s")==0 || strcmp(argv[i],"--save")==0 )
		 {	if ( (++i)>=argc )
			 {	fprint_err_invarg1(argv[i-1]);return(1);	}
			task=PEXEC_REMOTE_SAVE;
			file=argv[i];
		 }
		else if ( strcmp(argv[i],"-c")==0 || strcmp(argv[i],"--shell-command")==0 )
			is_shell_commands=1;
		else if ( strcmp(argv[i],"-a")==0 || strcmp(argv[i],"--atomic")==0 )
			task=PEXEC_REMOTE_ATOMIC;
		else if ( strcmp(argv[i],"--")==0 && task==PEXEC_REMOTE_ATOMIC )
		 {	azeroarg=i+1;
			break;
		 }
		else if ( argv[i][0]=='-' )
		 {	fprint_err_invarg0(argv[i]);	
			return(1);
		 }
		else if ( task==PEXEC_REMOTE_ATOMIC )
		 {	azeroarg=i;
			break;
		 }
		else
		 {	fprint_err_invarg0(argv[i]);	
			return(1);
		 }
	 }

	if ( ctrlport==NULL || strcmp(ctrlport,PEXEC_ISTR_AUTO)==0 )
		ctrlport=getenv(ctrlenv);
	if ( ctrlport==NULL )
	 {	fprintf(stderr,_("%s: error: connection port has not been defined or set in the environment.\n"),progbasename);
		return(1);
	 }

	sock=remote_control_port_connect(ctrlport);
	if ( sock<0 )
	 {	fprintf(stderr,_("%s: error: unable to connect to '%s'.\n"),progbasename,ctrlport);
		return(1);
	 }

	switch ( task )
	 {   case PEXEC_REMOTE_STATUS:
		ret=remote_status(sock,stdout);
		break;	
	     case PEXEC_REMOTE_LOCK:
		ret=remote_lock(sock,NULL,name);
		break;
	     case PEXEC_REMOTE_UNLOCK:
		ret=remote_unlock(sock,NULL,name);
		break;
	     case PEXEC_REMOTE_DUMP:
		if ( file != NULL && (fr=fopen(file,"rb")) != NULL && name != NULL )
		 {	ret=remote_copy(sock,name,fr,stdout);
			fclose(fr);
		 }
		else
		 {	fprintf(stderr,_("%s: error: invalid or unspecified input file (%s).\n"),progbasename,(file==NULL?"-":file)); 
			ret=1;
		 }
		break;
	     case PEXEC_REMOTE_SAVE:
		if ( file != NULL && (fw=fopen(file,"wb")) != NULL && name != NULL )
		 {	ret=remote_copy(sock,name,stdin,fw);
			fclose(fw);
		 }
		else
		 {	fprintf(stderr,_("%s: error: invalid or unspecified output file (%s).\n"),progbasename,(file==NULL?"-":file)); 
			ret=1;
		 }
		break;
	     case PEXEC_REMOTE_ATOMIC:
		if ( name != NULL )	remote_lock(sock,NULL,name);
		if ( azeroarg<argc )
			ret=remote_atomic_execute(is_shell_commands,p.shell,argc-azeroarg,argv+azeroarg,sock);
		else
			ret=0;
		if ( name != NULL )	remote_unlock(sock,NULL,name);
		break;
	     default:
		ret=1;
		break;
	 }

	remote_disconnect(sock);

	close(sock);

	return(ret);
  } 
 else if ( mode == PEXEC_MODE_HYPERVISOR )
  {	int	startstop;
	int	num_processes,use_load,use_fifo;
	int	sock;

	hypcport=NULL;
	startstop=-1;
	num_processes=-1;
	use_load=-1;
	use_fifo=0;
	
	for ( i=zeroarg ; i<argc ; i++ )
	 {	if ( strcmp(argv[i],"-n")==0 || strcmp(argv[i],"--number")==0 )
		 {	if ( (++i)>=argc )
			 {	fprint_err_invarg1(argv[i-1]);return(1);	}
			else if ( strcmp(argv[i],PEXEC_ISTR_AUTO)==0 || strcmp(argv[i],PEXEC_ISTR_NCPU)==0 )
				num_processes=PEXEC_MNP_NCPU;
			else if ( sscanf(argv[i],"%d",&num_processes)<1 || num_processes <= 0 )
			 {	fprint_err_invarg1(argv[i-1]);
				return(1);
			 }
		 }
		/* -l, --load: use load */
		else if ( strcmp(argv[i],"-l")==0 || strcmp(argv[i],"--load")==0 || strcmp(argv[i],"--use-load")==0 )
		 {	if ( (++i)>=argc )
			 {	fprint_err_invarg1(argv[i-1]);return(1);	}
			else if ( strcmp(argv[i],"0")==0 || strcmp(argv[i], "1m")==0 || strcmp(argv[i], "1min")==0 )
				use_load=0;
			else if ( strcmp(argv[i],"1")==0 || strcmp(argv[i], "5m")==0 || strcmp(argv[i], "5min")==0 )
				use_load=1;
			else if ( strcmp(argv[i],"2")==0 || strcmp(argv[i],"15m")==0 || strcmp(argv[i],"15min")==0 )
				use_load=2;
			else
			 {	fprint_err_invarg1(argv[i-1]);
				return(1);
			 }
		 }
		/* -C, --control: hypervisor socket binding port or path */
		else if ( strcmp(argv[i],"-C")==0 || strcmp(argv[i],"--control")==0 )
		 {	if ( (++i)>=argc )
			 {	fprint_err_invarg1(argv[i-1]);return(1);	}
			hypcport=argv[i];
		 }
		/* -L, --log: name of output log file or '-' or '-[12]' */
		else if ( strcmp(argv[i],"-L")==0 || strcmp(argv[i],"--log")==0 )
		 {	if ( (++i)>=argc )
			 {	fprint_err_invarg1(argv[i-1]);return(1);	}
			logfile=argv[i];
		 }
		/* -V, --verbose: increase log level by one: */
		else if ( strcmp(argv[i],"-V")==0 || strcmp(argv[i],"--verbose")==0 )
			log.loglevel++;
		/* -W, --log-level: log level, some non-negative integer */
		else if ( strcmp(argv[i],"-W")==0 || strcmp(argv[i],"--log-level")==0 )
		 {	if ( (++i)>=argc )
			 {	fprint_err_invarg1(argv[i-1]);return(1);	}
			else if ( sscanf(argv[i],"%d",&log.loglevel)<1 || log.loglevel<0 )
			 {	fprint_err_invarg1(argv[i-1]);
				return(1);
			 }
		 } 
		/* -f, --fifo: first in first out queue processing: */
		else if ( strcmp(argv[i],"-f")==0 || strcmp(argv[i],"--fifo")==0 )
			use_fifo=1;
		else if ( strcmp(argv[i],"-s")==0 || strcmp(argv[i],"--lifo")==0 || strcmp(argv[i],"--stack")==0 )
			use_fifo=0;
		/* start: start in daemon mode (i.e. detach from terminal): */ 
		else if ( strcmp(argv[i],"start")==0 )
			startstop=0;
		/* stop: stop running but detached pexec hypervisor: */
		else if ( strcmp(argv[i],"stop")==0 )
			startstop=1;
		else
		 {	fprint_err_invarg0(argv[i]);	
			return(1);
		 }
	 }

	if ( num_processes <= 0 )
		num_processes=get_number_of_cpus();

	if ( hypcport==NULL ) 
		hypcport=PEXEC_DEFAULT_HYPERVISOR_SOCKET;

	/* explicit stop: send stop request to the hypervisor daemon: */
	if ( startstop>0 )
	 {	sock=remote_control_port_connect(hypcport);
		if ( sock<0 )
		 {	fprintf(stderr,_("%s: error: unable to connect to hypervisor control socket '%s' in order to stop the service.\n"),progbasename,hypcport);
			return(1);
		 }
		pexec_hypervisor_stop(sock);
		close(sock);
		return(0);
	 }

	sock=remote_control_port_bind(hypcport,NULL,0,1);

	if ( sock<0 )
	 {	fprintf(stderr,_("%s: error: unable to create or bind control socket to %s.\n"),progbasename,hypcport);
		return(1);
	 }

	/* explicit start: put the process into background: */
	if ( ! (startstop<0) )	
		background(0,0);

	if ( sock>=0 )	listen(sock,256);

	/* well, it's turned out to be necessary: */
	signal(SIGPIPE,SIG_IGN);

	pexec_hypervisor_main_loop(sock,num_processes,use_load,use_fifo);

	if ( is_unix_socket_name(hypcport) )
		unlink(hypcport);

	return(0);
  }
 else
  {	fprintf(stderr,_("%s: error: internal: invalid mode code %d.\n"),progbasename,mode);
	return(2);
  }

}

/*****************************************************************************/
                                    
