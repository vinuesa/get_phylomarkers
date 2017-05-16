/*****************************************************************************/
/* pexec.h								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Simple program for shell command execution in parallel		     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Copyright (C) 2007, 2008-2009; Pal, A. (apal@szofi.elte.hu)		     */
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

#ifndef	__PEXEC_H_INCLUDED
#define	__PEXEC_H_INCLUDED	1

/*****************************************************************************/

#define			PEXEC_VERSION			"1.0rc8"
#define			PEXEC_LAST_MODIFICATION		"2009.07.02"

/*****************************************************************************/

#define			PEXEC_DEFAULT_SHELL		"/bin/sh"
#define			PEXEC_DEFAULT_RSH		"/usr/bin/ssh"
#define			PEXEC_DEFAULT_NULLFILE		"/dev/null"
#define			PEXEC_DEFAULT_ENVVARIABLE	"PEXEC_REMOTE_PORT"
#define			PEXEC_DEFAULT_HYPERVISOR_SOCKET	"/tmp/pexec.sock"
#define			PEXEC_LOAD_CHECK_PERIOD		2

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define			PEXEC_MODE_DEFAULT		0
#define			PEXEC_MODE_DAEMON		1
#define			PEXEC_MODE_REMOTECONTROL	2	
#define			PEXEC_MODE_HYPERVISOR		3

#define			PEXEC_REMOTE_STATUS		1
#define			PEXEC_REMOTE_LOCK		2
#define			PEXEC_REMOTE_UNLOCK		3
#define			PEXEC_REMOTE_DUMP		4
#define			PEXEC_REMOTE_SAVE		5
#define			PEXEC_REMOTE_ATOMIC		6

#define			PEXEC_MNP_AUTO			-1
#define			PEXEC_MNP_MANAGED		-2
#define			PEXEC_MNP_NCPU			-3

/*****************************************************************************/

typedef struct	child	child;

typedef struct
 {	int		pid;		/* pid (non-zero for real remote)    */
	int		fhsend,fhrecv;	/* stdin/stdout of the remote shell  */
	int		num_processes;	/* max num of processes or 0 (see V) */
	int		achild;		/* active associated children 	     */
	int		estatus;	/* execution status (if managed)     */
	linebuffer	lrsh;		/* line buffer for fhrecv	     */
 } remoteshell;

struct child
 {	child		*prev,*next;	/* list chain pointers		     */
	int		id;		/* internal child identifier	     */
	int		pid;		/* PID of the child		     */
	remoteshell	*rs;		/* remote shell pointer		     */
	int		fdstdout;	/* standard output pipe if collected */
	linebuffer 	lout;		/* stdout line buffer if reformatted */
	int		fdstderr;	/* standard error pipe, if collected */
	linebuffer 	lerr;		/* stderr line buffer if reformatted */
 };

typedef struct
 {	int	is_shell;
	char	**argv;		/* concat'ed with spaces if is_shell is true */
	int	argc;
 } command;

typedef struct
 {	char	*hostspec;
	int	num_processes;
 } remotehost;

typedef struct
 {	char	*name;
	int	no_touch_std;
	command	c;
	int	id;
	int	status;
 } parameter;

typedef struct	client	client;

struct client
 {	client		*prev,*next;	/* linked list pointers		     */
	int		peer;		/* socket of peer (server)	     */
	linebuffer	lcli;		/* main line buffer 		     */
	int		status;		/* some status, for arbitarary purp. */
 };

typedef struct
 {	int		qid;		/* queue id */
	void		*data;		/* qid >= 0: *rsh, otherwise *client */
 } pendingclient;

typedef struct	imutex	imutex;

struct imutex
 {	imutex		*prev,*next;
	char		*name;
	int		state;
	pendingclient	*pclients;
	int		npclient;
 };

typedef struct	dqueue	dqueue;

struct dqueue
 {	dqueue		*prev,*next;
	client		*qclient;
	int		id;
 };

typedef struct
 {	FILE	*fwlog;
	int	loglevel;
 } logdata;

typedef struct
 {	char	*in,
		*out,
		*err;
	FILE	*fwout;
	char	*formatout;
	FILE	*fwerr;
	char	*formaterr;
	int	omit_newlines;
	char	*envvarname;
	char	*shell;
	char	*rsh,*rshcmd,**rshargs;
	int	fallback_to_die;	
	logdata	*log;
 } paralleldata;

typedef struct
 {	child	*childlist;
	int	achild;
	int	sock;			/* remote control listen socket */
	int	hsck;			/* hypervisor connection socket */
	client	*clientlist;
	imutex	*imutexlist;
	dqueue	*dqueuelist;
	int	iqueue;
	int	nfinished,npending,nparam;
	time_t	t0;
 } parallelstatus;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

typedef struct	request	request;

struct request
 {	client	*cl;
 };

typedef struct
 {	client	*clientlist;		/* all clients connected	     */
	int	num_processes;		/* max number of processes	     */
	int	use_load;		/* 0, 1 or 2: 1, 5 or 15 min loadvg  */
	int	use_fifo;		/* use FIFO queuing insead of LIFO   */
	request	*requests;		/* list of requests		     */
	int	nrequest;		/* total requests		     */
	int	nrunning;		/* total running jobs		     */
 } hypervisorstatus;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

typedef struct
 {	int	signal;
	int	pid;
	int	exitstatus;
	int	exitsignal;
 } signalinfo;

/*****************************************************************************/

#endif

/*****************************************************************************/
                                                     
