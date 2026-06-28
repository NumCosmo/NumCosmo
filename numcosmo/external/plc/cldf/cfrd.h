#include <errno.h>
#include <stdio.h>
//#define _GNU_SOURCE
#include <string.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <poll.h>
#include <unistd.h>
#include <fcntl.h>
#include <limits.h>
#include "errorlist.h"
#include "io.h"

void cfrd_pipe_err(int twinpip[2],error **err);
int cfrd_pollout_err(int fd, error **err);
int cfrd_pollin_err(int fd, long delay, error **err);
void cfrd_send(int fd,char* mess,int nb,error **err);
int cfrd_read_err(int fd, char* buf,  error **err);
int cfrd_startup(int pipin[2],int pipout[2],char *arg0, char *arg1, error **err);
