# include "cfrd.h"

void cfrd_pipe_err(int twinpip[2],error **err) {
  testErrorRetVA(pipe(twinpip) == -1,-98689,"Cannot create unnamed pipe (%s)",*err,__LINE__,,strerror(errno));
}

int cfrd_pollout_err(int fd, error **err) {
  struct pollfd fds;  
  int pollerr;
  
  fds.fd     = fd;
  fds.events = POLLOUT;
  pollerr  = poll(&fds,1,0);
  //_DEBUGHERE_("","");
  //_DEBUGHERE_("->%d",pollerr);
  testErrorRetVA(pollerr==0,-241342,"TIMEOUT ! >%.3f s",*err,__LINE__,pollerr,0);
  testErrorRetVA(pollerr<0,-241342,"Error on poll (%s)",*err,__LINE__,pollerr,strerror(errno));
  testErrorRetVA(fds.revents & POLLHUP,-241342,"Hangup ! (%d)",*err,__LINE__,pollerr,fds.revents);
  testErrorRetVA(fds.revents & POLLERR,-241342,"Poll err ! (%d)",*err,__LINE__,pollerr,fds.revents);
  return pollerr;
  
}
int cfrd_pollin_err(int fd, long delay, error **err) {
  struct pollfd fds;  
  int pollerr;
  
  if (delay<0) {
    return -1;
  }
  //_DEBUGHERE_("%d",delay);
  fds.fd     = fd;
  fds.events = POLLIN;
  pollerr  = poll(&fds,1,delay/1000);
  //_DEBUGHERE_("","");
  //_DEBUGHERE_("->%d",pollerr);
  testErrorRetVA(pollerr==0,-241342,"TIMEOUT ! >%.3f s",*err,__LINE__,pollerr,delay/1e6);
  testErrorRetVA(pollerr<0,-241342,"Error on poll (%s)",*err,__LINE__,pollerr,strerror(errno));
  testErrorRetVA(fds.revents & POLLHUP,-241342,"Hangup ! (%d)",*err,__LINE__,pollerr,fds.revents);
  testErrorRetVA(fds.revents & POLLERR,-241342,"Poll err ! (%d)",*err,__LINE__,pollerr,fds.revents);
  return pollerr;
  
}
void cfrd_send(int fd,char* mess,int nb,error **err) {
  cfrd_pollout_err(fd, err);
  forwardError(*err,__LINE__,);
  if (nb<0) {
    nb = strlen(mess);
  }
  //_DEBUGHERE_("SEND %d '%s' (%d %d)",fd,mess,nb,strlen(mess));
  testErrorRetVA(write(fd,mess,nb)==-1,-241342,"Error writing (%s)",*err,__LINE__,,strerror(errno));
}

int cfrd_read_err(int fd, char* buf,  error **err) {
  int res;
  int i;

  for(i=0;i<500;i++) {
    //cfrd_pollin_err(fd, delay, err);
    //forwardError(*err,__LINE__,);
    res = read(fd,&(buf[i]),1);
    testErrorRetVA(res==-1,-241342,"Error on read(%s)",*err,_LINE__,-1,strerror(errno));
    if (buf[i]=='\n') {
      buf[i] = '\0';
      return i;
    }     
  }
  return 0;
}

int cfrd_startup(int pipin[2],int pipout[2],char *arg0, char *arg1, error **err) {
  int childpid;

  cfrd_pipe_err(pipin,err); //pfd[0] is for reading, //pfd[1] is for writing
  forwardError(*err,__LINE__,-1);
  cfrd_pipe_err(pipout,err); //pfd[0] is for reading, //pfd[1] is for writing
  forwardError(*err,__LINE__,-1);

  childpid=fork();
  testErrorRet(childpid == -1,-241342,"fork failed",*err,__LINE__,-1);

  /* this is the child */
  if (childpid == 0) {
    close(pipin[1]);
    close(pipout[0]);
    
    /* Duplicate stdin*/
    dup2(pipin[0],0);
    close(pipin[0]);   
    
    /* Duplicate stdout */
    dup2(pipout[1],1);
    close(pipout[1]);
    write(1,"test\n",5);
    execlp(arg0,arg0,arg1,(char*) NULL);
    // never returns from here, but we never know
    _DEBUGHERE_("bad oh bad%s","");
    exit(-1);
  } /* end of child process case */

  close(pipin[0]);
  close(pipout[1]);  
  return childpid;
}
