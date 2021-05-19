#include <glib.h>                                                                                                                        
#include <numcosmo/numcosmo.h>
#include <numcosmo/misc/CNearTree.h>
#include <numcosmo/misc/CVector.h>


gint
main(){
CNearTreeHandle treehandle;
int bReturn;
double v[3];
double * vBest;
void * vvBest;
double vSearch[3];
double dRad = 3.;
v[0] = 1.;
v[1]=2.;
v[2]=3.;
bReturn=!CNearTreeInsert(treehandle,&v[0],NULL);

if(!CNearTreeNearestNeighbor (treehandle,dRad,&vvBest,NULL,vSearch))
{vBest= (double *)vvBest;}
printf("%f", vSearch[1]);
return 0;
}