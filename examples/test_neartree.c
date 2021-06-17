/*
 *  main.c
 *  NearTree
 *
 *  Copyright 2001, 2008 Larry Andrews.  All rights reserved
 *  Revised 12 Dec 2008 for sourceforge release -- H. J. Bernstein
 */

/**********************************************************************
 *                                                                    *
 * YOU MAY REDISTRIBUTE NearTree UNDER THE TERMS OF THE LGPL          *
 *                                                                    *
 **********************************************************************/

/************************* LGPL NOTICES *******************************
 *                                                                    *
 * This library is free software; you can redistribute it and/or      *
 * modify it under the terms of the GNU Lesser General Public         *
 * License as published by the Free Software Foundation; either       *
 * version 2.1 of the License, or (at your option) any later version. *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  *
 * Lesser General Public License for more details.                    *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License along with this library; if not, write to the Free         *
 * Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,    *
 * MA  02110-1301  USA                                                *
 *                                                                    *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <numcosmo/misc/CVector.h>
#include <numcosmo/numcosmo.h>
#include <numcosmo/misc/rhrand.h>
CRHrand rhr;
#include <numcosmo/misc/CNearTree.h>

int
main (int argc, char **argv)
{
  CNearTreeHandle treehandle;
  double v[3], vSearch[3];
  double *vBest;
  void *vvBest;
  CVectorHandle vReturn;
  CVectorHandle dDs;
  CVectorHandle stIndices;
  CVectorHandle oReturn;
  double xdist;
  long i, j, k;
  const long lMaxRow = 10;
  double dRad        = 100000.;
  void *qualquer;
  double qualquer2[3];
  size_t indp;
  double w[3], n[3], l[3];
  NcmVector *ncmvec;
  double *pointervec;
  
/**  Aqui eu crio os vetores, que no caso seram nosso sample entao nao precisa criar **/
  w[0]       = 1.;
  w[1]       = 2;
  w[2]       = 3.;
  n[0]       = 23.;
  n[1]       = 1212;
  n[2]       = 2133.;
  l[0]       = 231.;
  l[1]       = 22;
  l[2]       = 43.;
  v[0]       = 11.;
  v[1]       = 12;
  v[2]       = 13.;
  vSearch[0] = 1., vSearch[1] = 2.;
  vSearch[2] = 3.;
  
/** Eu crio um objeto do neartree e crio alguns objetos CVector que eu vou precisar pras contas. Mudar tamanho deles pro tamanho da dimensao, da tree tambem ***/
  CNearTreeCreate (&treehandle, 3, CNEARTREE_TYPE_DOUBLE);
  CVectorCreate (&vReturn, sizeof (void *), 10);
  CVectorCreate (&oReturn, sizeof (void *), 10);
  CVectorCreate (&dDs, sizeof (double), 10);
  CVectorCreate (&stIndices, sizeof (size_t), 10);
  
/***Depois de criar a tree, a gente precisa inserir todos nossos vetores nela. Aqui vai entrar um FOR,
 *  ai a gente precisa copiar os vetores da numcosmo pra um vetor normal, e passar o pointer do primeiro termo de cada vetor
 *  pra inserir na arvore. No meu caso, os vetores v,n e l. Nao colocar na arvore o vetor que queremos procurar a distancia****/
  CNearTreeInsert (treehandle, &v[0], NULL);
  CNearTreeInsert (treehandle, &n[0], NULL);
  CNearTreeInsert (treehandle, &l[0], NULL);
  /*** Essa funcao acha todos os vetores a um distancia dRad do nosso. A gente pode por o tamanho maximo do sample e depois pegar so os primeiros 20% de vetores. ***/
  CNearTreeFindInSphere (treehandle, dRad, vReturn, oReturn, vSearch, 0);
  {
    size_t index, jndex, kndex, tndex;
    void *prevpoint;
    void *foundpoint;
    void *localmetrics;
    void *localindices;
    int redo;
    
    /**Um print que veio do exemplo, mas bom pra saber como usar as coisas.
     *  A gente n consegue pegar ou tirar elementos da arvore com a documentacao do neartree,
     *  precisa usar funcoes do CVector. Eu tentei fazer contas de outro jeito mas eu nao consegui
     *  remover vetores da arvore que eu ja tinha inserido, por meu plano era inserir todos, e ir procurando o vizinho mais proximo,
     *  e cada vez que achar um, remover da arvore. Usei CVectoreRemoveElement na estrutura da arvore, mas nao consegui**/
    fprintf (stdout, " Returned %lu items within %g of [%g,%g,%g]\n",
             (long unsigned int) vReturn->size, dRad, vSearch[0], vSearch[1], vSearch[2]);
    
    for (index = 0; index < CVectorSize (vReturn); index++) /*Nessa parte da pra por index < 0.2 * , ai cortaria o itnervalo, mas ele dai corta os menores 20%*/
    {
      CVectorGetElement (vReturn, &foundpoint, index);
      xdist = CNearTreeDist (treehandle, (void CNEARTREE_FAR *) foundpoint,
                             (void CNEARTREE_FAR *) vSearch);
      CNearTreeSortIn (dDs, stIndices, xdist, index, CVectorSize (vReturn));
    }
    
    /**Essa parte eu nao consegui entender muito bem. Quando ele calcula todos os vetores do raio, ele guarda em vReturn. E aqui pra baixo ele
     *  da um jeito de acessar todos eles e deixar em ordem do mais proximo pro mais longe**/
    CVectorGetElementptr (dDs, &localmetrics, 0);
    CVectorGetElementptr (stIndices, &localindices, 0);
    
    redo = 1;
    
    while (redo == 1)
    {
      redo = 0;
      
      for (index = 1; index < CVectorSize (stIndices); index++)
      {
        CVectorGetElement (stIndices, &jndex, index - 1);
        CVectorGetElement (stIndices, &kndex, index);
        CVectorGetElement (vReturn, &prevpoint, jndex);
        CVectorGetElement (vReturn, &foundpoint, kndex);
        
        if (fabs (((double *) localmetrics)[jndex] - ((double *) localmetrics)[kndex])
            <= DBL_EPSILON * fabs (((double *) localmetrics)[jndex] + ((double *) localmetrics)[kndex]))
        {
          if ((((double *) foundpoint)[0] < ((double *) prevpoint)[0])
              || ((((double *) foundpoint)[0] == ((double *) prevpoint)[0])
                  && (((double *) foundpoint)[1] < ((double *) prevpoint)[1]))
              || ((((double *) foundpoint)[0] == ((double *) prevpoint)[0])
                  && (((double *) foundpoint)[1] == ((double *) prevpoint)[1])
                  && (((double *) foundpoint)[1] < ((double *) prevpoint)[1]))
             )
          {
            tndex                                = ((size_t *) localindices)[index - 1];
            ((size_t *) localindices)[index - 1] = ((size_t *) localindices)[index];
            ((size_t *) localindices)[index]     = tndex;
            CVectorSetFlags (stIndices, 0);
            redo = 1;
          }
        }
      }
    }
    
    /**Quando acaba de organizar, ele imprime todos os vetores. Ai eh so pegar os 20% primeiros vetores e copiar pra um vetor da NumCosmo
     *  , e copiar pra um GPtrArray**/
    for (index = 0; index < CVectorSize (stIndices); index++)
    {
      CVectorGetElement (stIndices, &jndex, index);
      CVectorGetElement (vReturn, &foundpoint, jndex);
      xdist = CNearTreeDist (treehandle, (void CNEARTREE_FAR *) foundpoint, (void CNEARTREE_FAR *) vSearch);
      fprintf (stdout, "\t [%g,%g,%g] DISTANCIA %g\n ", ((double *) foundpoint)[0],
               ((double *) foundpoint)[1],
               ((double *) foundpoint)[2], xdist);
      pointervec = &((double *) foundpoint)[0];
      ncmvec     = ncm_vector_new_data_dup (pointervec, 3, 1);
    }
  }
  
  fprintf (stdout, " -------------------------------------------------------------\n");
  /* for i */
  
  CNearTreeFree (&treehandle);
}

