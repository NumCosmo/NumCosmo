/*
 *  CNearTree.c
 *  NearTree
 *
 *  Based on TNear.h C++ Template
 *  Copyright 2001 Larry Andrews.  All rights reserved
 *
 *  C Version created by Herbert J. Bernstein on 11/29/08
 *  with permission from Larry Andrews.
 *  Copyright 2008 Larry Andrews and Herbert J. Bernstein. 
 *  All rights reserved.
 *
 *  Revised 30 May 2009, release with full containerization of C++
 *                       version and KNear/Far in C++ and C, LCA + HJB
 */

/**********************************************************************
 *                                                                    *
 * YOU MAY REDISTRIBUTE THE CNearTree API UNDER THE TERMS OF THE LGPL *
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

/* Notices from original C++ template:
 
 Nearest Neighbor algorithm after Kalantari and McDonald,
 (IEEE Transactions on Software Engineering, v. SE-9, pp.
 631-634,1983)
 modified to use recursion instead of a double-linked tree
 and simplified so that it does a bit less checking for
 things like is the distance to the right less than the
 distance to the left; it was found that these checks little
 to no difference.
 
 copyright by Larry Andrews, 2001
 may be freely distributed or used as long as this copyright notice
 is included
 */


#ifdef __cplusplus

extern "C" {
    
#endif
    
#ifndef USE_LOCAL_HEADERS
#include <CNearTree.h>
#else
#include "CNearTree.h"
#endif
#include <math.h>
    
#ifdef CNEARTREE_SAFE_TRIANG
#define TRIANG(a,b,c) (  (((b)+(c))-(a) >= 0) \
                      || ((b)-((a)-(c)) >= 0) \
                      || ((c)-((a)-(b)) >= 0))    
#else
#define TRIANG(a,b,c) (  (((b)+(c))-(a) >= 0))
#endif
    
#ifndef CNEARTREE_DIMSAMPLES
#define CNEARTREE_DIMSAMPLES 4
#endif
    
#define max2(x,y) ((x)>=(y)?(x):(y)) 
    
    
    /*
     =======================================================================
     double CNearTreeDistsq(void CNEARTREE_FAR * coord1, 
     void CNEARTREE_FAR * coord2,
     size_t treedim, long treetype)
     
     function to return the square of the Euclidean distance between two 
     coordinate vectors.
     
     treedim -- the dimension of the vectors
     treetype -- a long integer flag for type of the vectors
     CNEARTREE_TYPE_DOUBLE for double
     CNEARTREE_TYPE_INTEGER for integer
     CNEARTREE_TYPE_STRING for strings
     
     =======================================================================
     */
    
    double CNearTreeDistsq(void CNEARTREE_FAR * coord1, 
                           void CNEARTREE_FAR * coord2,
                           size_t treedim, long treetype) {
        size_t index;
        double distsq;
        
        if (treetype == CNEARTREE_TYPE_DOUBLE) {
            
            double CNEARTREE_FAR * dcoord1;
            double CNEARTREE_FAR * dcoord2;
            
            dcoord1 = (double CNEARTREE_FAR *)coord1;
            dcoord2 = (double CNEARTREE_FAR *)coord2;
            
            distsq = (dcoord1[0]-dcoord2[0])*(dcoord1[0]-dcoord2[0]);
            
            for (index=1; index < treedim; index++) {
                distsq += (dcoord1[index]-dcoord2[index])
                *(dcoord1[index]-dcoord2[index]);
            }
            
            return distsq;
            
        } else if (treetype == CNEARTREE_TYPE_INTEGER) {
            
            int CNEARTREE_FAR * icoord1;
            int CNEARTREE_FAR * icoord2;
            
            icoord1 = (int CNEARTREE_FAR *)coord1;
            icoord2 = (int CNEARTREE_FAR *)coord2;
            
            distsq = ((double)(icoord1[0]-icoord2[0]))*((double)(icoord1[0]-icoord2[0]));
            
            for (index=1; index < treedim; index++) {
                distsq += ((double)(icoord1[index]-icoord2[index]))
                *((double)(icoord1[index]-icoord2[index]));
            }
            
            return distsq;
            
            
        } else if (treetype == CNEARTREE_TYPE_STRING) {
            
            char CNEARTREE_FAR * string1;
            char CNEARTREE_FAR * string2;
            char CNEARTREE_FAR * stringx;
            char c1, c2;
            
            string1 = (char CNEARTREE_FAR *)coord1;
            string2 = (char CNEARTREE_FAR *)coord2;

            distsq = 0;
            
            c1=c2=0;
            stringx = NULL;
            for (index=0; index < treedim; index++) {
                c1 = string1[index];
                c2 = string2[index];
                if (c1 == 0) {
                    if (c2 == 0) break;
                    stringx = (char CNEARTREE_FAR *)coord2;
                    break;
                } else if (c2 == 0) {
                    stringx = (char CNEARTREE_FAR *)coord1;
                    break;
                }
                if (c1 != c2) distsq++;
            }
            if (index < treedim && stringx) {
                for (; index < treedim; index++) {
                    if (stringx[index] == 0) break;
                    distsq++;
                }
            }
            return distsq*distsq;
        
        } else return -1.0;
        
    }
    
    /*
     =======================================================================
     int CNearTreeSetNorm(const CNearTreeHandle treehandle, int treenorm);
     
     function to set the norm to use for this tree
     
     treenorm should be one of CNEARTREE_NORM_L1 for an L-1 norm
     CNEARTREE_NORM_L2 for an L-2 norm
     CNEARTREE_NORM_LINF for an L-infinity norm
     CNEARTREE_NORM_L2LAZY for an L-2 norm close in,
                           but an L-1 norm for pruning.
     CNEARTREE_NORM_SPHERE for a sphere-based norm
     CNEARTREE_NORM_HSPHERE for a hemisphere-based norm
     CNEARTREE_NORM_HAMMING for a Hamming distance norm
     
     the function returns CNEARTREE_BAD_ARGUMENT for an invalid argument
     CNEARTREE_SUCCESS (0) otherwise
     =======================================================================
     */
    
    int CNearTreeSetNorm(const CNearTreeHandle treehandle, int treenorm) {
        
        if (!treehandle ||
            (treenorm != CNEARTREE_NORM_L1
             && treenorm != CNEARTREE_NORM_L2
             && treenorm != CNEARTREE_NORM_LINF
             && treenorm != CNEARTREE_NORM_L2LAZY
             && treenorm != CNEARTREE_NORM_SPHERE
             && treenorm != CNEARTREE_NORM_HSPHERE
             && treenorm != CNEARTREE_NORM_HAMMING)
            || (treehandle->m_iflags & CNEARTREE_NORM)!=CNEARTREE_NORM_UNKNOWN ) return CNEARTREE_BAD_ARGUMENT;
        treehandle->m_iflags &= ~CNEARTREE_NORM_UNKNOWN;
        treehandle->m_iflags |= treenorm;
        return CNEARTREE_SUCCESS;
    }

    /*
     =======================================================================
     double CNearTreeDist(const CNearTreeHandle treehandle, 
     void CNEARTREE_FAR * coord1, 
     void CNEARTREE_FAR * coord2))
     
     function to return the distance (L1, L2 or L-infinity) between two 
     coordinate vectors according to the parameters of the given tree
     For L2LAZY tree this returns the L2 norm distance.
     
     =======================================================================
     */
    
    double CNearTreeDist(const CNearTreeHandle treehandle, 
                         void CNEARTREE_FAR * coord1,
                         void CNEARTREE_FAR * coord2)  {
        size_t index;
        double dist, distsq;
        
        size_t treedim;
        long treetype;
        long treenorm;
        
        double CNEARTREE_FAR * dcoord1 = NULL;
        double CNEARTREE_FAR * dcoord2 = NULL;
        int CNEARTREE_FAR * icoord1 = NULL;
        int CNEARTREE_FAR * icoord2 = NULL;
        char CNEARTREE_FAR * svalue1 = NULL;
        char CNEARTREE_FAR * svalue2 = NULL;
        
        
        treedim = treehandle->m_szdimension;
        treetype = treehandle->m_iflags&CNEARTREE_TYPE;
        treenorm = treehandle->m_iflags&CNEARTREE_NORM;
        
        if (!treehandle || !coord1 || !coord2) return -1.;
        
        if (treenorm == CNEARTREE_NORM_UNKNOWN || treenorm == 0) {
            treehandle->m_iflags &= ~CNEARTREE_NORM;
            if (treetype == CNEARTREE_TYPE_STRING) {
                treehandle->m_iflags |= CNEARTREE_NORM_HAMMING;
                treenorm = CNEARTREE_NORM_HAMMING;
            } else {
                treehandle->m_iflags |= CNEARTREE_NORM_L2;
                treenorm = CNEARTREE_NORM_L2;
            }
            
        }
        
        
        if (treetype == CNEARTREE_TYPE_DOUBLE) {
            
            dcoord1 = (double CNEARTREE_FAR *)coord1;
            dcoord2 = (double CNEARTREE_FAR *)coord2;
        } else if (treetype == CNEARTREE_TYPE_INTEGER) {
            icoord1 = (int CNEARTREE_FAR *)coord1;
            icoord2 = (int CNEARTREE_FAR *)coord2;
        } else if (treetype == CNEARTREE_TYPE_STRING) {
            svalue1 = (char CNEARTREE_FAR *)coord1;
            svalue2 = (char CNEARTREE_FAR *)coord2;
        } else return -1.0;
        
        if (treedim == 1) {
            if (treetype == CNEARTREE_TYPE_DOUBLE) {
                return fabs(dcoord1[0]-dcoord2[0]);
            } else if (treetype == CNEARTREE_TYPE_INTEGER) {
                return fabs((double)(icoord1[0]-icoord2[0]));
            } else if (treetype == CNEARTREE_TYPE_STRING) {
                return (svalue1[0]==svalue2[0])?0.:1.;
            } else {
                return -1.0;
            }
        }
        
        switch (treenorm) {
            case CNEARTREE_NORM_HAMMING:
                if (treetype == CNEARTREE_TYPE_STRING) {
                    char c1, c2;
                    char CNEARTREE_FAR * stringx;
                    dist = 0.;
                    stringx = NULL;
                    for (index=0; index < treedim; index++){
                        c1 = svalue1[index];
                        c2 = svalue2[index];
                        if (c1 == 0) {
                            if (c2 == 0) break;
                            stringx = (char CNEARTREE_FAR *)coord2;
                            break;
                        } else if (c2 == 0) {
                            stringx = (char CNEARTREE_FAR *)coord1;
                            break;
                        }
                        if (c1 != c2) dist++;
                    }
                    if (index < treedim && stringx) {
                        for (; index < treedim; index++) {
                            if (stringx[index] == 0) break;
                            dist++;
                        }                        
                    }
                } else return -1.0;
                return dist;
            case CNEARTREE_NORM_L1:
                if (treetype == CNEARTREE_TYPE_DOUBLE) {
                    dist= fabs(dcoord1[0]-dcoord2[0]);
                    for (index=1; index < treedim; index++) {
                        dist += fabs(dcoord1[index]-dcoord2[index]);
                    }
                } else if (treetype == CNEARTREE_TYPE_INTEGER) {
                    dist = fabs((double)(icoord1[0]-icoord2[0]));
                    for (index=1; index < treedim; index++) {
                        dist += fabs((double)(dcoord1[index]-dcoord2[index]));
                    }
                } else return -1.0;
                return dist;
                                
            case CNEARTREE_NORM_LINF:
                if (treetype == CNEARTREE_TYPE_DOUBLE) {
                    dist= fabs(dcoord1[0]-dcoord2[0]);
                    for (index=1; index < treedim; index++) {
                        dist = max2(dist,fabs(dcoord1[index]-dcoord2[index]));
                    }
                } else if (treetype == CNEARTREE_TYPE_INTEGER) {
                    dist = fabs((double)(icoord1[0]-icoord2[0]));
                    for (index=1; index < treedim; index++) {
                        dist = max2(dist,fabs((double)(dcoord1[index]-dcoord2[index])));
                    }
                } else return -1.0;
                return dist;
                
            case CNEARTREE_NORM_SPHERE:
            case CNEARTREE_NORM_HSPHERE:
                if (treetype == CNEARTREE_TYPE_DOUBLE) {
                    double dot, cosangle, angle, norm1, norm2, norm;
                    dot = norm1 = norm2 = 0.;
                    for (index=0; index < treedim; index++) {
                        dot += dcoord1[index]*dcoord2[index];
                        norm1 += dcoord1[index]*dcoord1[index];
                        norm2 += dcoord2[index]*dcoord2[index];
                    }
                    norm1 = sqrt(norm1);
                    norm2 = sqrt(norm2);
                    if (norm1 <= DBL_MIN) {
                        dist = norm2;
                    } else if (norm2 <= DBL_MIN) {
                        dist = norm1;
                    } else {
                        cosangle = dot/(norm1*norm2);
                        if (cosangle > 1.) cosangle = 1.;
                        if (cosangle < -1.) cosangle = -1.;
                        if (treenorm == CNEARTREE_NORM_HSPHERE && cosangle < 0.) cosangle = -cosangle;
                        
                        angle = atan2(sqrt(1.-cosangle*cosangle),cosangle);
                        norm = (norm1<norm2?norm1:norm2);
                        dist = norm*fabs(angle)+fabs(norm1-norm2);
                    }
                    
                } else if (treetype == CNEARTREE_TYPE_INTEGER) {
                    double dot, cosangle, angle, norm1, norm2, norm;
                    dot = norm1 = norm2 = 0.;
                    for (index=0; index < treedim; index++) {
                        dot += ((double)icoord1[index])*((double)icoord2[index]);
                        norm1 += ((double)icoord1[index])*((double)icoord1[index]);
                        norm2 += ((double)icoord2[index])*((double)icoord2[index]);
                    }
                    norm1 = sqrt(norm1);
                    norm2 = sqrt(norm2);
                    if (norm1 <= DBL_MIN) {
                        dist = norm2;
                    } else if (norm2 <= DBL_MIN) {
                        dist = norm1;
                    } else {
                        cosangle = dot/(norm1*norm2);
                        if (cosangle > 1.) cosangle = 1.;
                        if (cosangle < -1.) cosangle = -1.;
                        
                        angle = atan2(sqrt(1.-cosangle*cosangle),cosangle);
                        norm = (norm1<norm2?norm1:norm2);
                        dist = norm*fabs(angle)+fabs(norm1-norm2);
                    }
                } else return -1.0;
                return dist;
                
            case CNEARTREE_NORM_L2:
            case CNEARTREE_NORM_L2LAZY:
            default:
                if (treetype == CNEARTREE_TYPE_DOUBLE) {
                    distsq = (dcoord1[0]-dcoord2[0])*(dcoord1[0]-dcoord2[0]);
                    
                    for (index=1; index < treedim; index++) {
                        distsq += (dcoord1[index]-dcoord2[index])
                        *(dcoord1[index]-dcoord2[index]);
                    }
                } else if (treetype == CNEARTREE_TYPE_INTEGER) {
                    distsq = ((double)(icoord1[0]-icoord2[0]))*((double)(icoord1[0]-icoord2[0]));
                    
                    for (index=1; index < treedim; index++) {
                        distsq += ((double)(icoord1[index]-icoord2[index]))
                        *((double)(icoord1[index]-icoord2[index]));
                    }
                } else return -1;
                return sqrt(distsq);
        }
        
    }
    
    
    /*
     =======================================================================
     int CNearTreeNodeCreate ( const CNearTreeHandle treehandle,  
     CNearTreeNodeHandle CNEARTREE_FAR * treenodehandle)     
     
     Create a CNearTreeNode
     
     returns a pointer to the newly allocated block of memory as a 
     CNearTreeNodeHandle in *treenodehandle
     
     flags are inherited from the treehandle  
     
     creates an empty tree with no right or left node and with the dMax-below
     set to negative values so that any match found will be stored since it will
     greater than the negative value
     
     =======================================================================
     */
    
    int CNearTreeNodeCreate (const CNearTreeHandle treehandle,
                             CNearTreeNodeHandle CNEARTREE_FAR * treenodehandle) {
        
        long treenorm, treetype, treexflags;
        
        if (!treehandle || !treenodehandle) return CNEARTREE_BAD_ARGUMENT;
        
        treenorm = (treehandle->m_iflags) & CNEARTREE_NORM;
        
        if (!treenorm) treenorm = CNEARTREE_NORM_UNKNOWN;
        
        treetype = (treehandle->m_iflags) & CNEARTREE_TYPE;
        
        if ( (treetype != CNEARTREE_TYPE_DOUBLE) 
              && (treetype != CNEARTREE_TYPE_INTEGER)
              && (treetype != CNEARTREE_TYPE_STRING) ) return CNEARTREE_BAD_ARGUMENT;
        
        treexflags =  (treehandle->m_iflags) & CNEARTREE_XFLAGS;
        
        *treenodehandle = (CNearTreeNodeHandle)CNEARTREE_MALLOC(sizeof(CNearTreeNode));
        if (!(*treenodehandle)) {
            return CNEARTREE_MALLOC_FAILED;
        }
        
        (*treenodehandle)->m_indexLeft  = 0;
        (*treenodehandle)->m_indexRight = 0;
        (*treenodehandle)->m_dMaxLeft     = -1.;        /* negative for an empty branch */
        (*treenodehandle)->m_dMaxRight    = -1.;        /* negative for an empty branch */
        (*treenodehandle)->m_iflags       = treetype;   /* no data, no children */
        (*treenodehandle)->m_iflags      |= treenorm;
        (*treenodehandle)->m_iflags      |= treexflags;
        (*treenodehandle)->m_pLeftBranch  = NULL;
        (*treenodehandle)->m_pRightBranch = NULL;
        (*treenodehandle)->m_iTreeSize    = 0;
#ifdef CNEARTREE_INSTRUMENTED
        (*treenodehandle)->m_Height       = 0;
#endif
        
        return CNEARTREE_SUCCESS;
    }
    
    
    /*
     =======================================================================
     int CNearTreeCreate ( CNearTreeHandle CNEARTREE_FAR * treehandle, 
     size_t treedim, long treetype)
     
     Create a CNearTree
     
     returns a pointer to the newly allocated block of memory as a 
     CNearTreeHandle in *treehandle
     
     
     treedim -- the dimension of the vectors
     treetype -- double or integer flag for type of the vectors ored with norm
     and ored with execution flags
     CNEARTREE_TYPE_DOUBLE for double
     CNEARTREE_TYPE_INTEGER for integer
     CNEARTREE_TYPE_STRING for strings
     ored with
     CNEARTREE_NORM_L1        for the sum of the absolute values
     CNEARTREE_NORM_L2        for the square root of the sum of the squares
     CNEARTREE_NORM_L2LAZY    for the square root of the sum of the squares
                              but with the tree as L1
     CNEARTREE_NORM_LINF      for the max
     CNEARTREE_NORM_SPHERE    for norm as spherical angular distance
     CNEARTREE_NORM_HSPHERE   for norm as hemispherical angular distance
     CNEARTREE_NORM_HAMMING   for norm as string hamming distance
     ored with
     CNTF_NOPREPRUNE    0x10000L    flag to suppress all search prepruning
     CNTF_FORCEPREPRUNE 0x20000L    flag to force search prepruning
     CNTF_NOFLIP        0x40000L    flag to suppress flips on insert
     CNTF_FORCEFLIP     0x80000L    flag to force flips on insert
     CNTF_NODEFER      0x100000L    flag to prevent deferred insert 
     
     
     creates an empty tree with no right or left node and with the dMax-below
     set to negative values so that any match found will be stored since it will
     greater than the negative value
     
     =======================================================================
     */
    
    int CNearTreeCreate ( CNearTreeHandle CNEARTREE_FAR * treehandle, 
                         size_t treedim, long treetype) {
        
        long treenorm;
        
        long treeflip;
        
        long treedefer;
        
        long treepreprune;
        
        if (!treehandle) return CNEARTREE_BAD_ARGUMENT;
        
        treenorm = treetype & CNEARTREE_NORM;
        
        treeflip = treetype & CNEARTREE_FLIP;
        
        treedefer = treetype & CNTF_NODEFER;
        
        treepreprune = treetype & (CNTF_NOPREPRUNE|CNTF_FORCEPREPRUNE);
        
        treeflip = treetype & (CNTF_NOFLIP|CNTF_FORCEFLIP);
                
        if (!treenorm) treenorm = CNEARTREE_NORM_UNKNOWN;
                
        treetype &= CNEARTREE_TYPE;
        
        if (treedim == 0 
            || (treetype != CNEARTREE_TYPE_DOUBLE 
                && treetype != CNEARTREE_TYPE_INTEGER
                && treetype != CNEARTREE_TYPE_STRING)) return CNEARTREE_BAD_ARGUMENT;
        
        *treehandle = (CNearTreeHandle)CNEARTREE_MALLOC(sizeof(CNearTree));
        if (!(*treehandle)) {
            return CNEARTREE_MALLOC_FAILED;
        }
        
        (*treehandle)->m_iflags        = treetype;  /* no data, no children         */
        (*treehandle)->m_iflags       |= treenorm;  /* record the chosen norm       */
        (*treehandle)->m_iflags       |= treeflip;  /* record whether to flip       */
        (*treehandle)->m_iflags       |= treedefer; /* record whether to defer      */
        (*treehandle)->m_iflags       |= treeflip;  /* record whether to flip       */
        (*treehandle)->m_iflags       |= treepreprune;
                                                    /* record whether to preprune   */
        (*treehandle)->m_szdimension   = treedim;   /* number of ints or doubles    */
        (*treehandle)->m_szsize        = 0;         /* number of nodes in the tree  */
        (*treehandle)->m_szdepth       = 0;         /* depth of the tree            */
        if( CNearTreeNodeCreate(*treehandle,&((*treehandle)->m_ptTree))) {
            CNEARTREE_FREE(*treehandle);
            return CNEARTREE_MALLOC_FAILED;
        }
        (*treehandle)->m_DelayedIndices = NULL;
        (*treehandle)->m_ObjectStore    = NULL;
        (*treehandle)->m_CoordStore     = NULL;
        (*treehandle)->m_DiamEstimate   = 0.;
        (*treehandle)->m_SumSpacings    = 0.;
        (*treehandle)->m_SumSpacingsSq  = 0.;
        (*treehandle)->m_DimEstimate    = 0;
        (*treehandle)->m_DimEstimateEsd = 0;
#ifdef CNEARTREE_INSTRUMENTED
        (*treehandle)->m_NodeVisits     = 0;
#endif
        CRHrandSrandom(&((*treehandle)->m_rhr),0);

         
        return CNEARTREE_SUCCESS;
        
    }
    
    /*
     =======================================================================
     int CNearTreeNodeFree ( CNearTreeNodeHandle CNEARTREE_FAR * treenodehandle )
     
     Free a CNearTreeNode
     
     recursively frees the NearTreeNode tree with the handle *trenodeehandle
     and nulls the treenodehandle.
     
     note that the objects referenced are not freed.
     =======================================================================
     */
    
    int CNearTreeNodeFree ( CNearTreeNodeHandle CNEARTREE_FAR * treenodehandle ) {
        
        int errorcode;
        
        if (!treenodehandle) return CNEARTREE_BAD_ARGUMENT;
        
        errorcode = CNEARTREE_SUCCESS;
        
        if ((*treenodehandle)->m_pLeftBranch) {
            errorcode |= CNearTreeNodeFree(&((*treenodehandle)->m_pLeftBranch));
        }
        
        if ((*treenodehandle)->m_pRightBranch) {
            errorcode |= CNearTreeNodeFree(&((*treenodehandle)->m_pRightBranch));
        }
        
        CNEARTREE_FREE(*treenodehandle);
        
        *treenodehandle = NULL;
        
        return errorcode;
    }
    /*
     =======================================================================
     int CNearTreeFree ( CNearTreeHandle CNEARTREE_FAR * treehandle )
     
     Free a CNearTree
     
     Frees the NearTree with the handle *treehandle
     and nulls the treehandle.
     
     note that the objects referenced are not freed.
     =======================================================================
     */
    
    int CNearTreeFree ( CNearTreeHandle CNEARTREE_FAR * treehandle ) {
        
        int errorcode, errorcodev;
        
        if (!treehandle) return CNEARTREE_BAD_ARGUMENT;
        
        errorcode = CNEARTREE_SUCCESS;
        
        if ((*treehandle)->m_ptTree) {
            errorcode |= CNearTreeNodeFree(&((*treehandle)->m_ptTree));
        }
        
        if ((*treehandle)->m_DelayedIndices) {
            errorcodev = CVectorFree(&((*treehandle)->m_DelayedIndices));
            if (errorcodev) errorcode |= CNEARTREE_FREE_FAILED ;
        }
        
        if ((*treehandle)->m_ObjectStore) {
            errorcodev = CVectorFree(&((*treehandle)->m_ObjectStore));
            if (errorcodev) errorcode |= CNEARTREE_FREE_FAILED ;
        }
        
        if ((*treehandle)->m_CoordStore) {
            errorcodev = CVectorFree(&((*treehandle)->m_CoordStore));
            if (errorcodev) errorcode |= CNEARTREE_FREE_FAILED ;
        }
        
        CNEARTREE_FREE(*treehandle);
        
        *treehandle = NULL;
        
        return errorcode;
    }
    
    /*
     =======================================================================
     int CNearTreeClear ( const CNearTreeHandle treehandle )
     
     Clear a CNearTree
     
     Clears the NearTree with the handle *treehandle
     
     note that the objects referenced are not freed.
     =======================================================================
     */
    
    int CNearTreeClear ( const CNearTreeHandle treehandle ) {
        
        int errorcode, errorcodev;
        
        if (!treehandle) return CNEARTREE_BAD_ARGUMENT;
        
        errorcode = CNEARTREE_SUCCESS;
        errorcodev = 0;
        
        if (treehandle->m_ptTree) {
            errorcode |= CNearTreeNodeFree(&(treehandle->m_ptTree));
            treehandle->m_ptTree = NULL;
        }
        
        errorcode |= CNearTreeNodeCreate(treehandle,&((treehandle)->m_ptTree));
        
        if (treehandle->m_DelayedIndices) {
            errorcodev = CVectorFree(&(treehandle->m_DelayedIndices));
            if (errorcodev) errorcode |= CNEARTREE_FREE_FAILED ;
            if (!errorcodev) treehandle->m_DelayedIndices = NULL;
        }
        
        if (treehandle->m_ObjectStore) {
            errorcodev = CVectorFree(&(treehandle->m_ObjectStore));
            if (errorcodev) errorcode |= CNEARTREE_FREE_FAILED ;
            if (!errorcodev) treehandle->m_ObjectStore = NULL;
        }
        
        if (treehandle->m_CoordStore) {
            errorcodev = CVectorFree(&(treehandle->m_CoordStore));
            if (errorcodev) errorcode |= CNEARTREE_FREE_FAILED ;
            if (!errorcodev) treehandle->m_CoordStore = NULL;
        }
        
        if (!errorcode) {
            (treehandle)->m_szsize        = 0;         /* number of nodes in the tree  */
            (treehandle)->m_szdepth       = 0;         /* depth of in the tree         */
            (treehandle)->m_DiamEstimate  = 0.;
            (treehandle)->m_SumSpacings   = 0.;
            (treehandle)->m_SumSpacingsSq = 0.;
            (treehandle)->m_DimEstimate   = 0;
            (treehandle)->m_DimEstimateEsd= 0;
#ifdef CNEARTREE_INSTRUMENTED
            (treehandle)->m_NodeVisits    = 0;
#endif            
        }
        return errorcode;
    }
    
    /*
     =======================================================================
     int CNearTreeZeroIfEmpty (const CNearTreeHandle treehandle)
     
     Test for an empty CNearTree, returning 0 in that case
     
     =======================================================================
     */
    
    int CNearTreeZeroIfEmpty ( const CNearTreeHandle treehandle )
    {
        return (treehandle == NULL
                || treehandle->m_CoordStore == NULL
                || CVectorSize(treehandle->m_CoordStore) == 0)?0:1;
    }
    
    /*
     =======================================================================
     int CNearTreeGetSize (const CNearTreeHandle treehandle, size_t CNEARTREE_FAR * size)
     
     Return the number of objects in the tree or delay queue in size
     
     =======================================================================
     */
    
    int CNearTreeGetSize ( const CNearTreeHandle treehandle, size_t CNEARTREE_FAR * size )
    {
        if (!treehandle || !size ) return CNEARTREE_BAD_ARGUMENT;
        
        *size = treehandle->m_szsize;
        
        if (treehandle->m_DelayedIndices) *size += CVectorSize(treehandle->m_DelayedIndices);
        
        return CNEARTREE_SUCCESS;
        
    }
    
    /*
     =======================================================================
     int CNearTreeGetDelayedSize (const CNearTreeHandle treehandle, size_t CNEARTREE_FAR * size)
     
     Return the number of objects in the delay queue tree in size
     This is a deprecated alternate name for CNearTreeGetDeferredSize
     
     =======================================================================
     */
    
    int CNearTreeGetDelayedSize ( const CNearTreeHandle treehandle, size_t CNEARTREE_FAR * size )
    {
        if (!treehandle || !size ) return CNEARTREE_BAD_ARGUMENT;
        
        *size = 0;
        
        if (treehandle->m_DelayedIndices) *size = CVectorSize(treehandle->m_DelayedIndices);
        
        return CNEARTREE_SUCCESS;
        
    }
    
    /*
     =======================================================================
     int CNearTreeGetDeferredSize (const CNearTreeHandle treehandle, size_t CNEARTREE_FAR * size)
     
     Return the number of objects in the delay queue tree in size
     
     =======================================================================
     */
    
    int CNearTreeGetDeferredSize ( const CNearTreeHandle treehandle, size_t CNEARTREE_FAR * size )
    {
        if (!treehandle || !size ) return CNEARTREE_BAD_ARGUMENT;
        
        *size = 0;
        
        if (treehandle->m_DelayedIndices) *size = CVectorSize(treehandle->m_DelayedIndices);
        
        return CNEARTREE_SUCCESS;
        
    }
    
    
    /*
     =======================================================================
     int CNearTreeGetDepth (const CNearTreeHandle treehandle, size_t CNEARTREE_FAR * depth)
     
     Return the depth of the tree in depth
     
     =======================================================================
     */
    
    int CNearTreeGetDepth ( const CNearTreeHandle treehandle, size_t CNEARTREE_FAR * depth )
    {
        if (!treehandle || !depth ) return CNEARTREE_BAD_ARGUMENT;
        
        *depth = treehandle->m_szdepth;
        
        return CNEARTREE_SUCCESS;
        
    }
    
    /*
     =======================================================================
     int CNearTreeGetHeight (const CNearTreeHandle treehandle, size_t CNEARTREE_FAR * height)
     
     Return the height of the tree in height, of instrumented, otherwise the depth
     
     =======================================================================
     */
    
    
    int CNearTreeGetHeight ( const CNearTreeHandle treehandle, size_t CNEARTREE_FAR * height )
    {
        if (!treehandle || !height ) return CNEARTREE_BAD_ARGUMENT;
        
#ifdef CNEARTREE_INSTRUMENTED

        *height = 0;
        
        if (treehandle->m_ptTree) *height = treehandle->m_ptTree->m_Height;
        
        *height = treehandle->m_szdepth;
        
        return CNEARTREE_SUCCESS;

#else
        return CNearTreeGetDepth (treehandle, height);
#endif
        
    }
    
    /*
     =======================================================================
     
     int CNearTreeGetFlags(const CNearTreeHandle treehandle, 
     long CNEARTREE_FAR *flags, 
     const long mask)
     
     int CNearTreeSetFlags(const CNearTreeHandle treehandle,
     const long flags,
     const long mask)
     
     functions to get or set the execution control flags:
     
     CNTF_NOPREPRUNE    to suppress all search prepruning
     CNTF_FORCEPREPRUNE to force search prepruning      
     CNTF_NOFLIP        flag to suppress flips on insert
     CNTF_FORCEFLIP     flag to force flips on insert
     CNTF_NODEFER       flag to prevent deferred insert
     
     The desired flags may be ored.  If mask is non-zero it is
     used as a bit mask, so a call with a flags of zero and a
     mask equal to a particular flag, clears that flag.
     =======================================================================
     */
    
    int CNearTreeGetFlags(const CNearTreeHandle treehandle, 
                          long CNEARTREE_FAR *flags, 
                          const long mask) {
        if (!flags) {
            if (mask) {
                *flags = (treehandle->m_iflags)&mask;
            } else {
                *flags = treehandle->m_iflags;
            }
            return CNEARTREE_SUCCESS;
        } else {
            return CNEARTREE_BAD_ARGUMENT;
        }
    }
    
    int CNearTreeSetFlags(const CNearTreeHandle treehandle,
                          const long flags,
                          const long mask){
        if (mask) {
            treehandle->m_iflags = (flags&mask) | ((treehandle->m_iflags)&~mask);
        } else {
            treehandle->m_iflags = flags;
        }
        return CNEARTREE_SUCCESS;
    }
    
    
    /*
     =======================================================================
      int CNearTreeGetMeanSpacing ( const CNearTreeHandle treehandle, 
                                    double CNEARTREE_FAR * spacing  )
    
      Get an estimate of the spacing of points
      
    
     =======================================================================
     */
     
    int CNearTreeGetMeanSpacing ( const CNearTreeHandle treehandle, 
                                 double CNEARTREE_FAR * spacing  )
    {
        if ( !treehandle || !spacing ) return CNEARTREE_BAD_ARGUMENT;
        
        *spacing = treehandle->m_SumSpacings/(double)((1+CVectorSize(treehandle->m_DelayedIndices)));
        
        return CNEARTREE_SUCCESS;
    }
    
     /*
     =======================================================================
      int CNearTreeGetVarSpacing ( const CNearTreeHandle treehandle, 
                                   double CNEARTREE_FAR * varspacing  )
    
      Get an estimate of variance of the spacing of points
      
    
     =======================================================================
      */
    int CNearTreeGetVarSpacing ( const CNearTreeHandle treehandle, 
                                 double CNEARTREE_FAR * varspacing  )

    {
        double meanSpacing;
        
        if ( !treehandle || !varspacing ) return CNEARTREE_BAD_ARGUMENT;

        meanSpacing = 
          treehandle->m_SumSpacings/(double)((1+CVectorSize(treehandle->m_DelayedIndices)));

        *varspacing = treehandle->m_SumSpacingsSq/
          (double)((1+CVectorSize(treehandle->m_DelayedIndices)))-meanSpacing*meanSpacing;

        return CNEARTREE_SUCCESS;

    }
    
    /*
     =======================================================================
     int CNearTreeCount(const CNearTreeHandle treehandle, 
                        size_t CNEARTREE_FAR * count)
     =======================================================================
     */
    
     int CNearTreeCount(const CNearTreeHandle treehandle, 
                        size_t CNEARTREE_FAR * count){
         
         *count = 0;
         
         return (CNearTreeNodeCount(treehandle->m_ptTree,count));
                
     }
    
    
#ifndef CNEARTREE_INSTRUMENTED
     /*
     =======================================================================
      int CNearTreeGetNodeVisits (  const CNearTreeHandle treehandle,
                                    size_t CNEARTREE_FAR * visits)
    
      Get the number of visits to nodes
      Dummy version to return 0
      
    
     =======================================================================
      */
    int CNearTreeGetNodeVisits ( const CNearTreeHandle treehandle,
                                 size_t CNEARTREE_FAR * visits)
    {
        if ( !treehandle || !visits ) return CNEARTREE_BAD_ARGUMENT;
        
        *visits = 0;
        
        return CNEARTREE_SUCCESS;
        
    }
    
#endif
    
#ifdef CNEARTREE_INSTRUMENTED
     /*
     =======================================================================
      int CNearTreeGetNodeVisits ( const CNearTreeHandle treehandle,
                                   size_t CNEARTREE_FAR * visits)
      
      Get the number of visits to nodes
      
    
     =======================================================================
      */
    int CNearTreeGetNodeVisits ( const CNearTreeHandle treehandle,
                                size_t CNEARTREE_FAR * visits)
    {
        if ( !treehandle || !visits ) return CNEARTREE_BAD_ARGUMENT;
        
        *visits = treehandle->m_NodeVisits;
        
        return CNEARTREE_SUCCESS;

    }
     /*
     =======================================================================
      int CNearTreeSetNodeVisits (  const CNearTreeHandle treehandle,
                                    const size_t visits )
    
      Set the number of visits to nodes
      
    
     =======================================================================
      */
     int CNearTreeSetNodeVisits (  const CNearTreeHandle treehandle,
                                  const size_t visits )
    {
        if ( !treehandle ) return CNEARTREE_BAD_ARGUMENT;
        
        treehandle->m_NodeVisits = visits;
        
        return CNEARTREE_SUCCESS;
        
    }
#endif    
    
     /*
     =======================================================================
      int CNearTreeGetDiamEstimate (  const CNearTreeHandle treehandle,
                                      double CNEARTREE_FAR * diamest )
    
      Get an estimate of the diameter
      
    
     =======================================================================
      */
    int CNearTreeGetDiamEstimate ( const CNearTreeHandle treehandle,
                                   double CNEARTREE_FAR * diamest )
    {
        if ( !treehandle || !diamest ) return CNEARTREE_BAD_ARGUMENT;

        if ( treehandle->m_DelayedIndices ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }        
        
        if (treehandle->m_DiamEstimate<=0.) {
            
            CNearTreeNodeHandle pt = treehandle->m_ptTree;
            
            treehandle->m_DiamEstimate = 0.;
            
            if (pt) {
                
                if ((pt->m_iflags&CNEARTREE_FLAG_LEFT_DATA)
                    &&(pt->m_iflags&CNEARTREE_FLAG_RIGHT_DATA)){
                    
                    CNTM_DistL2(treehandle->m_DiamEstimate,treehandle, 
                                CVectorElementAt(treehandle->m_CoordStore,pt->m_indexLeft),
                                CVectorElementAt(treehandle->m_CoordStore,pt->m_indexRight));
                }
                
                if ((pt->m_iflags&CNEARTREE_FLAG_LEFT_CHILD) &&
                    pt->m_dMaxLeft > treehandle->m_DiamEstimate)
                    treehandle->m_DiamEstimate = pt->m_dMaxLeft;
                
                if ((pt->m_iflags&CNEARTREE_FLAG_RIGHT_CHILD) &&
                    pt->m_dMaxRight > treehandle->m_DiamEstimate)
                    treehandle->m_DiamEstimate = pt->m_dMaxRight;                
            
            }
            
        }

        *diamest = treehandle->m_DiamEstimate;
        
        return CNEARTREE_SUCCESS;
        
    }

    
     /*
     =======================================================================
      int CNearTreeGetDimEstimateEsd ( const CNearTreeHandle treehandle,
                                 double CNEARTREE_FAR * dimestesd )
    
      Get the current best estimate of the dimension esd
    
     =======================================================================
      */
     int CNearTreeGetDimEstimateEsd ( const CNearTreeHandle treehandle,
                              double CNEARTREE_FAR * dimestesd )
    {
        double dimest;
        
        int errorcode;
        
        if ( !treehandle || !dimestesd ) return CNEARTREE_BAD_ARGUMENT;
        
        if ( treehandle->m_DelayedIndices ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }        

        if ( treehandle->m_DimEstimate == 0.) {
            
            if (!(errorcode=CNearTreeGetDimEstimate(treehandle,&dimest,0.1))) {
                
                *dimestesd = DBL_MAX;
                
                return errorcode;

                
            }

        }
        
        *dimestesd = treehandle->m_DimEstimateEsd;
        
        return CNEARTREE_SUCCESS;
    } 
    
    
     /*
     =======================================================================
      int CNearTreeGetDimEstimate ( const CNearTreeHandle treehandle,
                                       double CNEARTREE_FAR * diamest,
                                       const double DimEstimateEsd )
    
      Get an estimate of the dimension of the collection of points
      in the tree, to within the specified esd
      
    
     =======================================================================
      */
    int CNearTreeGetDimEstimate ( const CNearTreeHandle treehandle,
                                 double CNEARTREE_FAR * dimest,
                                 const double DimEstimateEsd )
    {
        size_t estsize;
        double estd;
        double estdim = 0.;
        double estdimsq = 0.;
        size_t trials;
        size_t n, ii;
        size_t probe_index;
        double testlim;
        double meanSpacing;
        double pointdensity;
        double targetradius;
        double rat;
        double shrinkfactor;
        double dummy;
        int bResult;
        size_t elsize;
        long goodtrials;
        
        size_t poplarge, popsmall, poptrial;
        CVectorHandle sampledisklarge;
        
        if ( !treehandle || !dimest ) return CNEARTREE_BAD_ARGUMENT;
        
        if ( treehandle->m_DelayedIndices ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }
        
        elsize = sizeof(double);
        if ((treehandle->m_iflags)&CNEARTREE_TYPE_DOUBLE) {
            elsize = sizeof(double);
        } else if ((treehandle->m_iflags)&CNEARTREE_TYPE_INTEGER) {
            elsize = sizeof(int);
        } else elsize = sizeof(char);

        if ( treehandle->m_DimEstimate == DBL_MAX ) {
            *dimest = 0.;
            return CNEARTREE_NOT_FOUND;
        }
        if ( treehandle->m_DimEstimate > 0. 
            && (treehandle->m_DimEstimateEsd <= DimEstimateEsd || 
                DimEstimateEsd <= 0.) ) {
            *dimest = treehandle->m_DimEstimate;
            return CNEARTREE_SUCCESS;
        }
        
        estsize = CVectorSize(treehandle->m_CoordStore);
                
        testlim = (DimEstimateEsd<=0.)?0.01:(DimEstimateEsd*DimEstimateEsd);
        meanSpacing = treehandle->m_SumSpacings/(double)((1+treehandle->m_szsize));
        probe_index = 0;
        
        
        /*  Do not try to get a dimension extimate with fewer than
         32 points or a diameter less than DBL_EPSILON*/
        
        if (estsize < 32 || (double)treehandle->m_DiamEstimate < DBL_EPSILON) {
            treehandle->m_DimEstimate = treehandle->m_DimEstimateEsd = DBL_MAX;
            *dimest = 0.;
           return CNEARTREE_NOT_FOUND;
        }
        
        /* Estimate the number of points per unit distance
         and a target radius that would produce 4096 points
         in dimension 1.  If this would bring us beyond
         the diameter/1.1, reduce to that size.  */
        pointdensity = ((double)estsize)/((double)treehandle->m_DiamEstimate);
        targetradius = 4096./pointdensity;
        if (targetradius <  meanSpacing*10.) targetradius = meanSpacing*10.;
        if (targetradius > (double)treehandle->m_DiamEstimate/1.1)
            targetradius = (double)treehandle->m_DiamEstimate/1.1;
        
        /*  Now try to find a smaller adjusted target radius that will
         contain a reasonable number of points*/
        
        shrinkfactor = 4.;
        n = (size_t)(((double)estsize-1u) * (CRHrandUrand(&(treehandle->m_rhr))));
        dummy = CRHrandUrand(&(treehandle->m_rhr)); dummy=CRHrandUrand(&(treehandle->m_rhr));
        
        if (CVectorCreate(&(sampledisklarge), (treehandle->m_szdimension)*elsize,10)) {
                      return CNEARTREE_MALLOC_FAILED;
        }
                      
        bResult=!CNearTreeFindInSphere( treehandle, (double)targetradius/shrinkfactor,
                                       sampledisklarge, NULL, 
                                       CVectorElementAt(treehandle->m_CoordStore,n), 1);
        poptrial = CVectorSize(sampledisklarge);
        do { 
            shrinkfactor = shrinkfactor/1.1;
            popsmall=poptrial;
            bResult=!CNearTreeFindInSphere( treehandle, (double)targetradius/shrinkfactor,
                                           sampledisklarge, NULL, 
                                           CVectorElementAt(treehandle->m_CoordStore,n), 1);
            poptrial = CVectorSize(sampledisklarge);
            n = (size_t)(((double)estsize-1u) * (CRHrandUrand(&(treehandle->m_rhr))));
            dummy = CRHrandUrand(&(treehandle->m_rhr)); dummy=CRHrandUrand(&(treehandle->m_rhr));
        } while (poptrial < 256 && shrinkfactor > 1. && poptrial <= popsmall+10);
        
        targetradius /= shrinkfactor;
        targetradius *= 1.1;
        
        goodtrials = 0;
        trials = (size_t)sqrt(0.5+(double)estsize);
        if (trials < 10) trials = 10;
                
        for(ii = 0; ii < trials; ii++) {
            n = (size_t)(((double)estsize-1u) * (CRHrandUrand(&(treehandle->m_rhr))));
            dummy=CRHrandUrand(&(treehandle->m_rhr)); dummy=CRHrandUrand(&(treehandle->m_rhr));

            bResult=!CNearTreeFindInSphere( treehandle, (double)targetradius,
                                           sampledisklarge, NULL, 
                                           CVectorElementAt(treehandle->m_CoordStore,n), 1);
            poplarge = CVectorSize(sampledisklarge);
            if (poplarge > 0) {
                bResult=!CNearTreeFindInSphere( treehandle, (double)targetradius/1.1,
                                           sampledisklarge, NULL, 
                                           CVectorElementAt(treehandle->m_CoordStore,n), 1);
                popsmall = CVectorSize(sampledisklarge);
                if (popsmall > 0 && popsmall< poplarge) {
                    rat = (double)poplarge/(double)popsmall;
                    estd = log(rat)/log(1.1);
                    estdim += estd;
                    estdimsq += estd*estd;
                    goodtrials++;
                    if (goodtrials > (1L+(long)(trials))/2 && 
                        fabs(estdimsq/((double)goodtrials) - estdim*estdim/((double)(goodtrials*goodtrials))) <= testlim) break;                    
                }
            }
        }
        if (goodtrials < 1) {
            treehandle->m_DimEstimate = treehandle->m_DimEstimateEsd = DBL_MAX;
            *dimest = 0.;
            CVectorFree(&sampledisklarge);
            return CNEARTREE_NOT_FOUND;
        }
        treehandle->m_DimEstimate = estdim/((double)goodtrials);
        treehandle->m_DimEstimateEsd = sqrt(fabs(estdimsq/((double)goodtrials) 
                -  treehandle->m_DimEstimate*treehandle->m_DimEstimate));
        if (treehandle->m_DimEstimate + 3.*treehandle->m_DimEstimateEsd< 0.) {
            treehandle->m_DimEstimate = treehandle->m_DimEstimateEsd = DBL_MAX;
            *dimest = 0.;
            CVectorFree(&sampledisklarge);
            return CNEARTREE_NOT_FOUND;
        }
        *dimest = treehandle->m_DimEstimate;
        return CNEARTREE_SUCCESS;
    }
    
    
    
    
    
    /*
     =======================================================================
     int CNearTreeImmediateInsert ( const CNearTreeHandle treehandle, 
     const void CNEARTREE_FAR * coord, 
     const void * obj )
     
     Function to insert some "point" as an object into a CNearTree for
     later searching
     
     coord is a coordinate vector for an object, obj, to be inserted into a
     Neartree.  A static copy of the coordinates and a pointer to the
     object are inserted
     
     Three possibilities exist: put the datum into the left
     position (first test),into the right position, or else
     into a node descending from the nearer of those positions
     when they are both already used.
     
     return 0 for success, nonzero for an error
     
     =======================================================================
     */
    
    int CNearTreeImmediateInsert( const CNearTreeHandle treehandle,
                        const void CNEARTREE_FAR * coord, 
                        const void CNEARTREE_FAR * obj ) {
        
        /* do a bit of precomputing if it is possible so that we can
         reduce the number of calls to sqrt */
        
        size_t depth;
        
        size_t index;
        
        size_t elsize;
        
        int errorcode;
        
        if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        
        if ( !(treehandle->m_ptTree) ) return CNEARTREE_BAD_ARGUMENT;
                
        if ( treehandle->m_ObjectStore == NULL) {
            if (CVectorCreate(&(treehandle->m_ObjectStore),sizeof(void *),10)) {
                return CNEARTREE_MALLOC_FAILED;
            }
        }
        
        
        if ( treehandle->m_CoordStore == NULL) {
            elsize = sizeof(double);
            if ((treehandle->m_iflags)&CNEARTREE_TYPE_DOUBLE) {
                elsize = sizeof(double);
            } else if ((treehandle->m_iflags)&CNEARTREE_TYPE_INTEGER) {
                elsize = sizeof(int);
            } else elsize = sizeof(char);
            if (CVectorCreate(&(treehandle->m_CoordStore),
                              (treehandle->m_szdimension)*elsize,10)) {
                return CNEARTREE_MALLOC_FAILED;
            }
        }
        
        index = CVectorSize(treehandle->m_ObjectStore);
        if (CVectorAddElement(treehandle->m_ObjectStore,(void CNEARTREE_FAR *)&obj)) return CNEARTREE_CVECTOR_FAILED;
        if (CVectorAddElement(treehandle->m_CoordStore,(void CNEARTREE_FAR *)coord)) return CNEARTREE_CVECTOR_FAILED;
        
        
        depth = 0;
        
        
        if ((treehandle->m_iflags& CNTF_FORCEFLIP)||!(treehandle->m_iflags& CNTF_NOFLIP)) {
            
            errorcode = 
            CNearTreeNodeInsert_Flip( treehandle, treehandle->m_ptTree, index, &depth);
            
        } else {
            
            errorcode = 
            CNearTreeNodeInsert( treehandle, treehandle->m_ptTree, index, &depth);

        }
        if (errorcode == 0) {
            
            if (depth > treehandle->m_szdepth) treehandle->m_szdepth = depth;
            
            (treehandle->m_szsize)++;
            
        }
        
        return errorcode;
        
    }
    
    /* int CNearTreeNodeCount(const CNearTreeNodeHandle treenodehandle, size_t CNEARTREE_FAR * count)
     */
    int CNearTreeNodeCount(const CNearTreeNodeHandle treenodehandle, size_t CNEARTREE_FAR * count){
        
        int errorcode=0;
        
        if ( (treenodehandle->m_iflags)&CNEARTREE_FLAG_RIGHT_DATA ) {
            (*count)++;
        }
        
        if ( (treenodehandle->m_iflags)&CNEARTREE_FLAG_LEFT_DATA ) {
            (*count)++;
        }
        
        if ( (treenodehandle->m_iflags)&CNEARTREE_FLAG_RIGHT_CHILD ) {
            errorcode |= CNearTreeNodeCount(treenodehandle->m_pRightBranch, count);
        }
        
        if ( (treenodehandle->m_iflags)&CNEARTREE_FLAG_LEFT_CHILD ) {
            errorcode |= CNearTreeNodeCount(treenodehandle->m_pLeftBranch,count);
        }
        
        return errorcode;
        
    }
  /*  
    =======================================================================
    int CNearTreeNodeReInsert_Flip ( const CNearTreeHandle treehandle,
                                     const CNearTreeNodeHandle treenodehandle,
                                     const CNearTreeNodeHandle pntn,
                                      size_t CNEARTREE_FAR * depth)
    
    Function to reinsert the elements from a detached a node into a CNearTree 
    for later searching
   
    treehandle is the handle of the overall neartree being used
   
    treenodehandle is the handle of the node in the existing tree at which to start
    the insertion
        
    pntn is the handle of the previously detached node from which to 
    extract the nodes for insertion
        
    depth is used to keep track of the depth at which the insertion is done
                
    return 0 for success, nonzero for an error
                    
    =======================================================================
  */
                    
    int CNearTreeNodeReInsert_Flip( const CNearTreeHandle treehandle,
                            const CNearTreeNodeHandle treenodehandle,
                            const CNearTreeNodeHandle pntn,
                            size_t CNEARTREE_FAR * depth) {
        
        size_t tempdepth1, tempdepth2, tempdepth3, tempdepth4;
        
        double SumSpacings, SumSpacingsSq;
        
        int errorcode;
        
        errorcode = 0;
        
        tempdepth1 = tempdepth2 = tempdepth3 = tempdepth4 = *depth;
        
        SumSpacings = treehandle->m_SumSpacings;
        
        SumSpacingsSq = treehandle->m_SumSpacingsSq;
        
        if ( (pntn->m_iflags)&CNEARTREE_FLAG_RIGHT_DATA ) {
            errorcode |= CNearTreeNodeInsert_Flip(treehandle,treenodehandle,
                                         pntn->m_indexRight,
                                         &tempdepth1);
        }

        if ( (pntn->m_iflags)&CNEARTREE_FLAG_LEFT_DATA ) {
             errorcode |= CNearTreeNodeInsert_Flip(treehandle,treenodehandle,
                                         pntn->m_indexLeft,
                                         &tempdepth2);
        }
        
        if ( (pntn->m_iflags)&CNEARTREE_FLAG_RIGHT_CHILD ) {
             errorcode |= CNearTreeNodeReInsert_Flip(treehandle, treenodehandle, pntn->m_pRightBranch,
                                       &tempdepth3);
        }

        if ( (pntn->m_iflags)&CNEARTREE_FLAG_LEFT_CHILD ) {
             errorcode |= CNearTreeNodeReInsert_Flip(treehandle, treenodehandle, pntn->m_pLeftBranch,
                                       &tempdepth4);
        }
        
        *depth = tempdepth1>*depth?tempdepth1:*depth;
        *depth = tempdepth2>*depth?tempdepth2:*depth;
        *depth = tempdepth3>*depth?tempdepth3:*depth;
        *depth = tempdepth4>*depth?tempdepth4:*depth;
        
        treehandle->m_SumSpacings = SumSpacings;
        treehandle->m_SumSpacingsSq = SumSpacingsSq;
        
        return errorcode;
    }
        
    
    /*
     =======================================================================
     int CNearTreeNodeInsert ( const CNearTreeHandle treehandle,
     const const CNearTreeNodeHandle treenodehandle, 
     size_t index;
     size_t CNEARTREE_FAR * depth)
     
     Function to insert some "point" as an object into a node in a CNearTree for
     later searching
     
     treenodehandle is the handle of the node at which to start the insertion
     
     index is the index of the object and coordinates to add from 
     treehandle->m_ObjectStore and treehandle->CoordStore into a NearTree.
     A static copy of the coordinates and a pointer to the
     object are inserted
     
     Three possibilities exist: put the datum into the left
     position (first test),into the right position, or else
     into a node descending from the nearer of those positions
     when they are both already used.
     
     depth is used to keep track of the depth at which the insertion is done
     
     return 0 for success, nonzero for an error
     
     =======================================================================
     =======================================================================
     int CNearTreeNodeInsert_Flip ( const CNearTreeHandle treehandle,
     const const CNearTreeNodeHandle treenodehandle, 
     size_t index;
     size_t CNEARTREE_FAR * depth)
     
     Function to insert some "point" as an object into a node in a CNearTree for
     later searching
     
     treenodehandle is the handle of the node at which to start the insertion
     
     index is the index of the object and coordinates to add from 
     treehandle->m_ObjectStore and treehandle->CoordStore into a NearTree.
     A static copy of the coordinates and a pointer to the
     object are inserted
     
     Three possibilities exist: put the datum into the left
     position (first test),into the right position, or else
     into a node descending from the nearer of those positions
     when they are both already used.
          
     depth is used to keep track of the depth at which the insertion is done
     
     return 0 for success, nonzero for an error
     
     =======================================================================
     */
    
    int CNearTreeNodeInsert( const CNearTreeHandle treehandle,
                                 const CNearTreeNodeHandle treenodehandle,
                                 const size_t index,
                                 size_t CNEARTREE_FAR * depth) {
        
        /* do a bit of precomputing if it is possible so that we can
         reduce the number of calls to sqrt */
        
        double dTempRight =  0.;
        double dTempLeft  =  0.;
        double dTempLeftRight = 0;
        int errorcode = 0;
        void CNEARTREE_FAR * coord;
        void CNEARTREE_FAR * coordLeft;
        void CNEARTREE_FAR * coordRight;
        const void CNEARTREE_FAR * obj;
        size_t n = index;
        (*depth)++;
        
        if ( !treehandle || !treenodehandle 
            || index+1 >  CVectorSize(treehandle->m_ObjectStore)) return CNEARTREE_BAD_ARGUMENT;
        
        obj = CVectorElementAt(treehandle->m_ObjectStore,index);
        coord = CVectorElementAt(treehandle->m_CoordStore,index);
        coordLeft = NULL;
        coordRight = NULL;
        (treenodehandle->m_iTreeSize)++;
        
        
        if ( !((treenodehandle->m_iflags)&CNEARTREE_FLAG_LEFT_DATA) ) {
            treenodehandle->m_indexLeft = n;
            treenodehandle->m_dMaxLeft = -1.;
            treenodehandle->m_iflags |= CNEARTREE_FLAG_LEFT_DATA;
#ifdef CNEARTREE_INSTRUMENTED
            treenodehandle->m_Height = 1;
#endif
            return CNEARTREE_SUCCESS;
        }
        
        coordLeft = CVectorElementAt(treehandle->m_CoordStore,treenodehandle->m_indexLeft);
        dTempLeft = CNearTreeDist(treehandle, coord, coordLeft);
         
        if (  !((treenodehandle->m_iflags)&CNEARTREE_FLAG_RIGHT_DATA)  ){ 
            treenodehandle->m_indexRight = n;
            treenodehandle->m_dMaxRight = -1.;
            treenodehandle->m_iflags |= CNEARTREE_FLAG_RIGHT_DATA;
            treehandle->m_SumSpacings += dTempLeft;
            treehandle->m_SumSpacingsSq += dTempLeft*dTempLeft;
            return CNEARTREE_SUCCESS;
        }
        
        coordRight = CVectorElementAt(treehandle->m_CoordStore,treenodehandle->m_indexRight);
        dTempRight = CNearTreeDist(treehandle, coord, coordRight);
        dTempLeftRight = CNearTreeDist(treehandle, coordLeft, coordRight);
        
        if ( dTempLeft > dTempRight ) {
            if (  !((treenodehandle->m_iflags)&CNEARTREE_FLAG_RIGHT_CHILD) ) {
                if ( (errorcode = CNearTreeNodeCreate(treehandle, &(treenodehandle->m_pRightBranch)))) return errorcode;
                treenodehandle->m_iflags |= CNEARTREE_FLAG_RIGHT_CHILD;
                treenodehandle->m_dMaxRight = dTempRight;
            } else if ( treenodehandle->m_dMaxRight < dTempRight ) {
                treenodehandle->m_dMaxRight = dTempRight;
            }
            
            /* If the left branch is empty, we are going to put something here */
            
            if (!((treenodehandle->m_pRightBranch->m_iflags) & CNEARTREE_FLAG_LEFT_DATA )) {
                treenodehandle->m_pRightBranch->m_indexLeft = n;
                treenodehandle->m_pRightBranch->m_dMaxLeft = -1.;
                treenodehandle->m_pRightBranch->m_iflags |= CNEARTREE_FLAG_LEFT_DATA;
                treenodehandle->m_pRightBranch->m_iTreeSize = 1;
                treehandle->m_SumSpacings += dTempRight;
                treehandle->m_SumSpacingsSq += dTempRight*dTempRight;
#ifdef CNEARTREE_INSTRUMENTED
                treenodehandle->m_pRightBranch->m_Height = 1;
                if (treenodehandle->m_Height < 2) treenodehandle->m_Height=2;
#endif
                (*depth)++;
                return CNEARTREE_SUCCESS;
            }
            if ( (errorcode=CNearTreeNodeInsert(treehandle, treenodehandle->m_pRightBranch, n, depth))) return errorcode;
#ifdef CNEARTREE_INSTRUMENTED
            treenodehandle->m_Height = 1+(treenodehandle->m_pRightBranch->m_Height);
            if ((treenodehandle->m_iflags & CNEARTREE_FLAG_LEFT_CHILD)
                && (treenodehandle->m_pLeftBranch->m_Height) >= treenodehandle->m_Height) treenodehandle->m_Height = 1+(treenodehandle->m_pLeftBranch->m_Height);
#endif
            return CNEARTREE_SUCCESS;
        } else { /* ((double)(t - *m_tLeft) <= (double)(t - *m_tRight) ) */
            
            if (  !((treenodehandle->m_iflags)&CNEARTREE_FLAG_LEFT_CHILD) ) {
                if ( (errorcode = CNearTreeNodeCreate(treehandle, &(treenodehandle->m_pLeftBranch)))) return errorcode;
                treenodehandle->m_iflags |= CNEARTREE_FLAG_LEFT_CHILD;
                treenodehandle->m_dMaxLeft = dTempLeft;
            } else if ( treenodehandle->m_dMaxLeft < dTempLeft ) {
                treenodehandle->m_dMaxLeft = dTempLeft;
            }
            
            /* If the left branch is empty, we are going to put something here */
            
            if (!((treenodehandle->m_pLeftBranch->m_iflags) & CNEARTREE_FLAG_LEFT_DATA )) {
                treenodehandle->m_pLeftBranch->m_indexLeft = n;
                treenodehandle->m_pLeftBranch->m_dMaxLeft = -1.;
                treenodehandle->m_pLeftBranch->m_iflags |= CNEARTREE_FLAG_LEFT_DATA;
                treenodehandle->m_pLeftBranch->m_iTreeSize = 1;
                treehandle->m_SumSpacings += dTempLeft;
                treehandle->m_SumSpacingsSq += dTempLeft*dTempLeft;
#ifdef CNEARTREE_INSTRUMENTED
                treenodehandle->m_pLeftBranch->m_Height = 1;
                if (treenodehandle->m_Height < 2) treenodehandle->m_Height=2;
#endif
                (*depth)++;
                return CNEARTREE_SUCCESS;
            }
            if ( (errorcode=CNearTreeNodeInsert(treehandle, treenodehandle->m_pLeftBranch, n, depth))) return errorcode;
#ifdef CNEARTREE_INSTRUMENTED
            treenodehandle->m_Height = 1+(treenodehandle->m_pLeftBranch->m_Height);
            if ((treenodehandle->m_iflags & CNEARTREE_FLAG_RIGHT_CHILD)
                && (treenodehandle->m_pRightBranch->m_Height) >= treenodehandle->m_Height) treenodehandle->m_Height = 1+(treenodehandle->m_pRightBranch->m_Height);
#endif
            return CNEARTREE_SUCCESS;
        }
        
    }
    
    
    int CNearTreeNodeInsert_Flip( const CNearTreeHandle treehandle,
                            const CNearTreeNodeHandle treenodehandle,
                            const size_t index,
                            size_t CNEARTREE_FAR * depth) {
        
        /* do a bit of precomputing if it is possible so that we can
         reduce the number of calls to sqrt */
        
        double dTempRight =  0.;
        double dTempLeft  =  0.;
        double dTempLeftRight = 0;
        int errorcode = 0;
        void CNEARTREE_FAR * coord;
        void CNEARTREE_FAR * coordLeft;
        void CNEARTREE_FAR * coordRight;
        const void CNEARTREE_FAR * obj;
        size_t n = index;
        (*depth)++;
        
        if ( !treehandle || !treenodehandle 
            || index+1 >  CVectorSize(treehandle->m_ObjectStore)) return CNEARTREE_BAD_ARGUMENT;
        
        obj = CVectorElementAt(treehandle->m_ObjectStore,index);
        coord = CVectorElementAt(treehandle->m_CoordStore,index);
        coordLeft = NULL;
        coordRight = NULL;
        (treenodehandle->m_iTreeSize)++;            
        
        if ( !((treenodehandle->m_iflags)&CNEARTREE_FLAG_LEFT_DATA) ) {
            treenodehandle->m_indexLeft = n;
            treenodehandle->m_dMaxLeft = -1.;
            treenodehandle->m_iflags |= CNEARTREE_FLAG_LEFT_DATA;
#ifdef CNEARTREE_INSTRUMENTED
            treenodehandle->m_Height = 1;
#endif
            return CNEARTREE_SUCCESS;
        }
        
        coordLeft = CVectorElementAt(treehandle->m_CoordStore,treenodehandle->m_indexLeft);
        dTempLeft = CNearTreeDist(treehandle, coord, coordLeft);
        
        
        if (  !((treenodehandle->m_iflags)&CNEARTREE_FLAG_RIGHT_DATA)  ){ 
            treenodehandle->m_indexRight = n;
            treenodehandle->m_dMaxRight = -1.;
            treenodehandle->m_iflags |= CNEARTREE_FLAG_RIGHT_DATA;
            treehandle->m_SumSpacings += dTempLeft;
            treehandle->m_SumSpacingsSq += dTempLeft*dTempLeft;
            return CNEARTREE_SUCCESS;
        }
        
        coordRight = CVectorElementAt(treehandle->m_CoordStore,treenodehandle->m_indexRight);
        dTempRight = CNearTreeDist(treehandle, coord, coordRight);
        dTempLeftRight = CNearTreeDist(treehandle, coordLeft, coordRight);
        
        if ( dTempLeft > dTempRight ) {
            if (  !((treenodehandle->m_iflags)&CNEARTREE_FLAG_RIGHT_CHILD) ) {
                if ( (errorcode = CNearTreeNodeCreate(treehandle, &(treenodehandle->m_pRightBranch)))) return errorcode;
                treenodehandle->m_iflags |= CNEARTREE_FLAG_RIGHT_CHILD;
                treenodehandle->m_dMaxRight = dTempRight;
            } else if ( treenodehandle->m_dMaxRight < dTempRight ) {
                treenodehandle->m_dMaxRight = dTempRight;
            }
            
            /* If the left branch is empty, we are going to put something here */

            if (!((treenodehandle->m_pRightBranch->m_iflags) & CNEARTREE_FLAG_LEFT_DATA )) {
                treenodehandle->m_pRightBranch->m_indexLeft = n;
                treenodehandle->m_pRightBranch->m_dMaxLeft = -1.;
                treenodehandle->m_pRightBranch->m_iflags |= CNEARTREE_FLAG_LEFT_DATA;
                treenodehandle->m_pRightBranch->m_iTreeSize = 1;
                treehandle->m_SumSpacings += dTempRight;
                treehandle->m_SumSpacingsSq += dTempRight*dTempRight;
#ifdef CNEARTREE_INSTRUMENTED
                treenodehandle->m_pRightBranch->m_Height = 1;
                if (treenodehandle->m_Height < 2) treenodehandle->m_Height=2;
#endif
                (*depth)++;
                /* See if it would be better to put the new node at this level and drop the current
                   Right node down one level */
                if (dTempRight > dTempLeftRight) {
                    treenodehandle->m_pRightBranch->m_indexLeft = treenodehandle->m_indexRight;
                    treenodehandle->m_indexRight = n;
                }
                return CNEARTREE_SUCCESS;
            }
            if ( (errorcode=CNearTreeNodeInsert(treehandle, treenodehandle->m_pRightBranch, n, depth))) return errorcode;
#ifdef CNEARTREE_INSTRUMENTED
            treenodehandle->m_Height = 1+(treenodehandle->m_pRightBranch->m_Height);
            if ((treenodehandle->m_iflags & CNEARTREE_FLAG_LEFT_CHILD)
                && (treenodehandle->m_pLeftBranch->m_Height) >= treenodehandle->m_Height) treenodehandle->m_Height = 1+(treenodehandle->m_pLeftBranch->m_Height);
#endif
            return CNEARTREE_SUCCESS;
        } else { /* ((double)(t - *m_tLeft) <= (double)(t - *m_tRight) ) */
            
            if (  !((treenodehandle->m_iflags)&CNEARTREE_FLAG_LEFT_CHILD) ) {
                if ( (errorcode = CNearTreeNodeCreate(treehandle, &(treenodehandle->m_pLeftBranch)))) return errorcode;
                treenodehandle->m_iflags |= CNEARTREE_FLAG_LEFT_CHILD;
                treenodehandle->m_dMaxLeft = dTempLeft;
            } else if ( treenodehandle->m_dMaxLeft < dTempLeft ) {
                treenodehandle->m_dMaxLeft = dTempLeft;
            }
            
            /* If the left branch is empty, we are going to put something here */
            
            if (!((treenodehandle->m_pLeftBranch->m_iflags) & CNEARTREE_FLAG_LEFT_DATA )) {
                treenodehandle->m_pLeftBranch->m_indexLeft = n;
                treenodehandle->m_pLeftBranch->m_dMaxLeft = -1.;
                treenodehandle->m_pLeftBranch->m_iflags |= CNEARTREE_FLAG_LEFT_DATA;
                treenodehandle->m_pLeftBranch->m_iTreeSize = 1;
                treehandle->m_SumSpacings += dTempLeft;
                treehandle->m_SumSpacingsSq += dTempLeft*dTempLeft;
#ifdef CNEARTREE_INSTRUMENTED
                treenodehandle->m_pLeftBranch->m_Height = 1;
                if (treenodehandle->m_Height < 2) treenodehandle->m_Height=2;
#endif
                (*depth)++;
                /* See if it would be better to put the new node at this level and drop the current
                 Left node down one level */
                if (dTempLeft > dTempLeftRight ) {
                    treenodehandle->m_pLeftBranch->m_indexLeft = treenodehandle->m_indexLeft;
                    treenodehandle->m_indexLeft = n;
                }
                return CNEARTREE_SUCCESS;
            }
            if ( (errorcode=CNearTreeNodeInsert(treehandle, treenodehandle->m_pLeftBranch, n, depth))) return errorcode;
#ifdef CNEARTREE_INSTRUMENTED
            treenodehandle->m_Height = 1+(treenodehandle->m_pLeftBranch->m_Height);
            if ((treenodehandle->m_iflags & CNEARTREE_FLAG_RIGHT_CHILD)
                && (treenodehandle->m_pRightBranch->m_Height) >= treenodehandle->m_Height) treenodehandle->m_Height = 1+(treenodehandle->m_pRightBranch->m_Height);
#endif
            return CNEARTREE_SUCCESS;
        }
        
    }
    
    
    /*
     =======================================================================
     int CNearTreeInsert ( const CNearTreeHandle treehandle, 
     const void CNEARTREE_FAR * coord, 
     const void * obj )
     
     Function to queue some "point" as an object for future insertion
     into a CNearTree for later searching
     
     coord is a coordinate vector for an object, obj, to be inserted into a
     Neartree.  A static copy of the coordinates and a pointer to the
     object are queued for insertion.  The exact order of insertion
     is not predetermined.  It will be partially randomized to help to
     balance the tree.
     
     The insertions will be completed by a call to 
     CNearTreeCompleteDelayedInsert(const CNearTreeHandle treehandle) 
     or by execution of any search.
     
     return 0 for success, nonzero for an error
     
     =======================================================================
     */
    
    int CNearTreeInsert( const CNearTreeHandle treehandle,
                               const void CNEARTREE_FAR * coord, 
                               const void CNEARTREE_FAR * obj ) {
        
        int treetype;
        
        size_t treedim;
        
        size_t index;
        
        int errorcode;
        
        errorcode =  CNEARTREE_SUCCESS;
        
        if (!treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
                
        treetype = (treehandle->m_iflags) & CNEARTREE_TYPE;
        
        treedim = treehandle->m_szdimension;
        
        if ( treehandle->m_ObjectStore == NULL) {
            if (CVectorCreate(&(treehandle->m_ObjectStore),sizeof(void *),10)) {
                return CNEARTREE_MALLOC_FAILED;
            }
        }
        
        if ( treehandle->m_CoordStore == NULL) {
            if (CVectorCreate(&(treehandle->m_CoordStore),
                              treedim*((treetype&CNEARTREE_TYPE_DOUBLE)?sizeof(double):sizeof(int)),10)) {
                return CNEARTREE_MALLOC_FAILED;
            }
        }
        
        index = CVectorSize(treehandle->m_ObjectStore);
        if (CVectorAddElement(treehandle->m_ObjectStore,(void CNEARTREE_FAR *)&obj)) return CNEARTREE_CVECTOR_FAILED;
        if (CVectorAddElement(treehandle->m_CoordStore,(void CNEARTREE_FAR *)coord)) return CNEARTREE_CVECTOR_FAILED;
        
        if (treehandle->m_DelayedIndices==NULL) {
            if (CVectorCreate(&(treehandle->m_DelayedIndices),sizeof(size_t),10)) {
                return CNEARTREE_MALLOC_FAILED;
            }
        }
        if (CVectorAddElement(treehandle->m_DelayedIndices,&index)) return CNEARTREE_CVECTOR_FAILED;
        
        if (treehandle->m_iflags & CNTF_NODEFER && treehandle->m_szdepth < 100) {
            errorcode = CNearTreeCompleteDelayedInsert(treehandle);
        }
        
        treehandle->m_DimEstimate = 0;
        
        treehandle->m_DimEstimateEsd= 0;
        
        return errorcode;
        
    }
    
    
    /*
     =======================================================================
     int CNearTreeDelayedInsert ( const CNearTreeHandle treehandle, 
     const void CNEARTREE_FAR * coord, 
     const void * obj )
     
     *** THIS IS A DEPRECATED ALTLERNATE CALL TO CNearTreeInsert ***
     
     Function to queue some "point" as an object for future insertion
     into a CNearTree for later searching
     
     coord is a coordinate vector for an object, obj, to be inserted into a
     Neartree.  A static copy of the coordinates and a pointer to the
     object are queued for insertion.  The exact order of insertion
     is not predetermined.  It will be partially randomized to help to
     balance the tree.
     
     The insertions will be completed by a call to 
     CNearTreeCompleteDelayedInsert(const CNearTreeHandle treehandle) 
     or by execution of any search.
     
     return 0 for success, nonzero for an error
     
     =======================================================================
     */
    
    int CNearTreeDelayedInsert( const CNearTreeHandle treehandle,
                        const void CNEARTREE_FAR * coord, 
                        const void CNEARTREE_FAR * obj ) {
        return CNearTreeInsert(treehandle,coord,obj);
    }
        
    /*
     =======================================================================
     int CNearTreeCompleteDelayedInsert ( const CNearTreeHandle treehandle )
     
     Function to dequeue the "points" queued as an objects for future insertion
     into a CNearTree for later searching
     
     return 0 for success, nonzero for an error
     
     =======================================================================
     */
    
    int CNearTreeCompleteDelayedInsert ( const CNearTreeHandle treehandle ) {
        
        size_t nqueued;
        
        size_t ntarget;
        
        size_t nrandom;
        
        int npass;
        
        size_t ielement, kelement=0, oelement;
        
        int errorcode;
        
        size_t depth;
        
        size_t errsize;
        
        size_t added;
        
        size_t dummyindex;
        
        size_t index;
        
        double dummyrand;
        
        if (!treehandle) return CNEARTREE_BAD_ARGUMENT;
        
        if (treehandle->m_DelayedIndices == NULL ) 
            return CNEARTREE_SUCCESS;
        
        if (treehandle->m_ObjectStore == NULL || treehandle->m_CoordStore == NULL ||
            CVectorSize(treehandle->m_ObjectStore) != CVectorSize(treehandle->m_CoordStore))
            return CNEARTREE_BAD_ARGUMENT;
        
        errorcode = 0;
        
        nqueued = CVectorSize(treehandle->m_DelayedIndices);
        
        ntarget = nqueued + treehandle->m_szsize;
        
        dummyindex = CVectorSize(treehandle->m_ObjectStore);
        
        nrandom = (size_t)sqrt((double)nqueued);
        
        if ((treehandle->m_iflags& CNTF_FORCEFLIP)||!(treehandle->m_iflags& CNTF_NOFLIP)){
            
            for (ielement = 0; ielement < nrandom; ielement++) {
                
                kelement = (int)(CRHrandUrand(&(treehandle->m_rhr))*((double)(nqueued)));
                
                dummyrand = CRHrandUrand(&(treehandle->m_rhr)) + CRHrandUrand(&(treehandle->m_rhr));
                
                oelement = kelement;
                do {
                    if (kelement >= nqueued) {
                        kelement = 0;
                        CVectorSetSize(treehandle->m_DelayedIndices,oelement);
                        nqueued = oelement;
                    }
                    if (CVectorGetElement(treehandle->m_DelayedIndices,&index,kelement)) return  CNEARTREE_BAD_ARGUMENT;
                    kelement ++;
                } while (index == dummyindex);
                if (kelement == nqueued) {
                    nqueued--;
                    CVectorSetSize(treehandle->m_DelayedIndices,nqueued);
                }
                kelement--;
                if (CVectorSetElement(treehandle->m_DelayedIndices,&dummyindex,kelement)) errorcode |= CNEARTREE_CVECTOR_FAILED;
                depth = 0;
                errorcode |= CNearTreeNodeInsert_Flip(treehandle,treehandle->m_ptTree,index,&depth);
                if (depth > treehandle->m_szdepth) treehandle->m_szdepth = depth;
                (treehandle->m_szsize)++;
            }
            
            npass = 0;
            
            while (treehandle->m_szsize < ntarget) {
                npass++;
                errsize = (treehandle->m_ptTree->m_iTreeSize)>>(-1+treehandle->m_szdepth);
                if ( errsize < 1 ) {

                    kelement = (int)(CRHrandUrand(&(treehandle->m_rhr))*((double)(nqueued)));
                    
                    dummyrand = CRHrandUrand(&(treehandle->m_rhr)) + CRHrandUrand(&(treehandle->m_rhr));
                    
                    oelement = kelement;
                    do {
                        if (kelement >= nqueued) {
                            kelement = 0;
                            CVectorSetSize(treehandle->m_DelayedIndices,oelement);
                            nqueued = oelement;
                        }
                        if (CVectorGetElement(treehandle->m_DelayedIndices,&index,kelement)) return  CNEARTREE_BAD_ARGUMENT;
                        kelement ++;
                    } while (index == dummyindex);
                    if (kelement == nqueued) {
                        nqueued--;
                        CVectorSetSize(treehandle->m_DelayedIndices,nqueued);
                    }
                    kelement--;
                    if (CVectorSetElement(treehandle->m_DelayedIndices,&dummyindex,kelement)) errorcode |= CNEARTREE_CVECTOR_FAILED;
                    depth = 0;
                    errorcode |= CNearTreeNodeInsert_Flip(treehandle,treehandle->m_ptTree,index,&depth);
                    if (depth > treehandle->m_szdepth) treehandle->m_szdepth = depth;
                    (treehandle->m_szsize)++;
                    
                } else {
                    
                    added=0;
                    
                    for (ielement = 0; ielement < nqueued; ielement++) {
                        if (CVectorGetElement(treehandle->m_DelayedIndices,&index,ielement)) return CNEARTREE_CVECTOR_FAILED;
                        if (index == dummyindex) continue;
                        depth = 0;
                        errorcode |= CNearTreeNodeInsert_Flip(treehandle,treehandle->m_ptTree,index,&depth);
                        if (CVectorSetElement(treehandle->m_DelayedIndices,&dummyindex,ielement)) errorcode |= CNEARTREE_CVECTOR_FAILED;
                        if (depth > treehandle->m_szdepth) treehandle->m_szdepth = depth;
                        (treehandle->m_szsize)++;
                        added++;
                        if (added > 100 || added > 8*npass) {
                          errsize = (treehandle->m_ptTree->m_iTreeSize)>>(-1+treehandle->m_szdepth);
                          if ( errsize < 1 ) break;  
                        }
                    }
                }
            }
            
        } else {
            
            for (ielement = 0; ielement < nrandom; ielement++) {
                
                kelement = (int)(CRHrandUrand(&(treehandle->m_rhr))*((double)(nqueued)));
                
                dummyrand = CRHrandUrand(&(treehandle->m_rhr)) + CRHrandUrand(&(treehandle->m_rhr));
                
                oelement = kelement;
                do {
                    if (kelement >= nqueued) {
                        kelement = 0;
                        CVectorSetSize(treehandle->m_DelayedIndices,oelement);
                        nqueued = oelement;
                    }
                    if (CVectorGetElement(treehandle->m_DelayedIndices,&index,kelement)) return  CNEARTREE_BAD_ARGUMENT;
                    kelement ++;
                } while (index == dummyindex);
                if (kelement == nqueued) {
                    nqueued--;
                    CVectorSetSize(treehandle->m_DelayedIndices,nqueued);
                }
                kelement--;
                if (CVectorSetElement(treehandle->m_DelayedIndices,&dummyindex,kelement)) errorcode |= CNEARTREE_CVECTOR_FAILED;
                depth = 0;
                errorcode |= CNearTreeNodeInsert(treehandle,treehandle->m_ptTree,index,&depth);
                if (depth > treehandle->m_szdepth) treehandle->m_szdepth = depth;
                (treehandle->m_szsize)++;
            }
            
            npass = 0;
            
            while (treehandle->m_szsize < ntarget) {
                npass++;
                errsize = (treehandle->m_ptTree->m_iTreeSize)>>(-1+treehandle->m_szdepth);
                if ( errsize < 1 ) {
                    kelement = (int)(CRHrandUrand(&(treehandle->m_rhr))*((double)(nqueued)));
                    
                    dummyrand = CRHrandUrand(&(treehandle->m_rhr)) + CRHrandUrand(&(treehandle->m_rhr));
                    
                    oelement = kelement;
                    do {
                        if (kelement >= nqueued) {
                            kelement = 0;
                            CVectorSetSize(treehandle->m_DelayedIndices,oelement);
                            nqueued = oelement;
                        }
                        if (CVectorGetElement(treehandle->m_DelayedIndices,&index,ielement)) return  CNEARTREE_BAD_ARGUMENT;
                        kelement ++;
                    } while (index == dummyindex);
                    if (kelement == nqueued) {
                        nqueued--;
                        CVectorSetSize(treehandle->m_DelayedIndices,nqueued);
                    }
                    kelement--;
                    if (CVectorSetElement(treehandle->m_DelayedIndices,&dummyindex,kelement)) errorcode |= CNEARTREE_CVECTOR_FAILED;
                    depth = 0;
                    errorcode |= CNearTreeNodeInsert(treehandle,treehandle->m_ptTree,index,&depth);
                    if (depth > treehandle->m_szdepth) treehandle->m_szdepth = depth;
                    (treehandle->m_szsize)++;
                    
                } else {
                    
                    added = 0;
                    
                    for (ielement = 0; ielement < nqueued; ielement++) {
                        if (CVectorGetElement(treehandle->m_DelayedIndices,&index,ielement)) return CNEARTREE_CVECTOR_FAILED;
                        if (index == dummyindex) continue;
                        depth = 0;
                        errorcode |= CNearTreeNodeInsert(treehandle,treehandle->m_ptTree,index,&depth);
                        if (CVectorSetElement(treehandle->m_DelayedIndices,&dummyindex,ielement)) errorcode |= CNEARTREE_CVECTOR_FAILED;
                        if (depth > treehandle->m_szdepth) treehandle->m_szdepth = depth;
                        (treehandle->m_szsize)++;
                        added++;
                        if (added > 100 || added > 8*npass) {
                            errsize = (treehandle->m_ptTree->m_iTreeSize)>>(-1+treehandle->m_szdepth);
                            if ( errsize < 1 ) break;
                        }
                    }
                }
            }
        }
        if (CVectorFree(&(treehandle->m_DelayedIndices))) errorcode |= CNEARTREE_CVECTOR_FAILED;
        
        return (errorcode != 0)?errorcode:CNEARTREE_SUCCESS;
    }
    
    /*
     =======================================================================
     int CNearTreeNearestNeighbor ( const CNearTreeHandle treehandle, 
     const double dRadius,  
     void CNEARTREE_FAR * CNEARTREE_FAR * coordClosest,
     void CNEARTREE_FAR * CNEARTREE_FAR * objClosest,   
     const void CNEARTREE_FAR * coord )
     
     Function to search a Neartree for the object closest to some probe point, coord. This function
     is only here so that the function CNearTreeNearest can be called without having dRadius const
     
     dRadius is the maximum search radius - any point farther than dRadius from the probe
     point will be ignored
     coordClosest is a pointer to the coordinate vector of the nearest point
     objClosest is the address into which a pointer to the object associated with coordClosest
     will be stored
     coord  is the probe point
     the return value is true only if a point was found
     
     This version does a balanced search

     =======================================================================
     =======================================================================

     int CNearTreeLeftNearestNeighbor ( const CNearTreeHandle treehandle, 
     const double dRadius,  
     void CNEARTREE_FAR * CNEARTREE_FAR * coordClosest,
     void CNEARTREE_FAR * CNEARTREE_FAR * objClosest,   
     const void CNEARTREE_FAR * coord )
     
     Function to search a Neartree for the object closest to some probe point, coord. This function
     is only here so that the function CNearTreeNearest can be called without having dRadius const
     
     dRadius is the maximum search radius - any point farther than dRadius from the probe
     point will be ignored
     coordClosest is a pointer to the coordinate vector of the nearest point
     objClosest is the address into which a pointer to the object associated with coordClosest
     will be stored
     coord  is the probe point
     the return value is true only if a point was found
     
     This version does a left search first
     
     =======================================================================
     */
    int CNearTreeNearestNeighbor (const CNearTreeHandle treehandle, 
                                  const double dRadius,  
                                  void CNEARTREE_FAR * CNEARTREE_FAR * coordClosest,
                                  void CNEARTREE_FAR * CNEARTREE_FAR * objClosest, 
                                  const void CNEARTREE_FAR * coord ) {
        
        double dSearchRadius = dRadius;
        if (!treehandle || ! coord ) return CNEARTREE_BAD_ARGUMENT;
        if ( treehandle->m_DelayedIndices ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }        
        if (!(treehandle->m_ptTree->m_iflags&CNEARTREE_FLAG_LEFT_DATA)) 
            return CNEARTREE_NOT_FOUND;
        return ( CNearTreeNearest ( treehandle, &dSearchRadius, coordClosest, objClosest, coord ) );
    }
    
    int CNearTreeLeftNearestNeighbor (const CNearTreeHandle treehandle, 
                                  const double dRadius,  
                                  void CNEARTREE_FAR * CNEARTREE_FAR * coordClosest,
                                  void CNEARTREE_FAR * CNEARTREE_FAR * objClosest, 
                                  const void CNEARTREE_FAR * coord ) {
        
        double dSearchRadius = dRadius;
        if (!treehandle || ! coord ) return CNEARTREE_BAD_ARGUMENT;
        if ( treehandle->m_DelayedIndices ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }        
        if (!(treehandle->m_ptTree->m_iflags&CNEARTREE_FLAG_LEFT_DATA)) 
            return CNEARTREE_NOT_FOUND;
        return ( CNearTreeLeftNearest ( treehandle, &dSearchRadius, coordClosest, objClosest, coord ) );
    }
    
    
    /*
     =======================================================================
     int CNearTreeFarthestNeighbor ( const CNearTreeHandle treehandle, 
     void CNEARTREE_FAR * CNEARTREE_FAR * coordFarthest,
     void CNEARTREE_FAR * CNEARTREE_FAR * objFarthest,   
     const void CNEARTREE_FAR * coord )
     
     Function to search a Neartree for the object farthest some probe point, coord.
     
     coordClosest is a pointer to the coordinate vector of the nearest point
     objClosest is the address into which a pointer to the object associated with coordClosest
     will be stored
     coord  is the probe point
     the return value is 0 only if a point was found
     
     =======================================================================
     */
    int CNearTreeFarthestNeighbor (const CNearTreeHandle treehandle, 
                                   void CNEARTREE_FAR *  CNEARTREE_FAR * coordFarthest,
                                   void CNEARTREE_FAR * CNEARTREE_FAR * objFarthest,   
                                   const void CNEARTREE_FAR * coord ) {
        double dSearchRadius = DBL_MIN;
        if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        if ( treehandle->m_DelayedIndices ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }        
        if (!(treehandle->m_ptTree->m_iflags&CNEARTREE_FLAG_LEFT_DATA)) 
            return CNEARTREE_NOT_FOUND;
        return ( CNearTreeFindFarthest ( treehandle, &dSearchRadius, coordFarthest, objFarthest, coord ) );
    }
    
    /*
     =======================================================================
     int CNearTreeFindInSphere ( const CNearTreeHandle treehandle, 
     const double dRadius,
     CVectorHandle coordclosest,
     CVectorHandle objClosest,
     const void CNEARTREE_FAR * coord,
     int resetcount)
     
     Function to search a Neartree for the set of objects closer to some probe point, coord,
     than dRadius. This is only here so that objClosest can be cleared before starting the work.
     
     dRadius is the maximum search radius - any point farther than dRadius from the probe
     point will be ignored
     coordClosest is a vector of pointers to coordinate tuples of nearest points
     objClosest is a vector of objects and is the returned set of nearest points
     to the probe point that can be found in the Neartree
     coord  is the probe point
     resetcount should be non-zero to clear coordclosest and objClosest on entry
     return value is 0 if points were found
     
     =======================================================================
     */
    
    int CNearTreeFindInSphereL2LAZY ( const CNearTreeHandle treehandle,
                                     const double dRadius,
                                     CVectorHandle coordClosest,
                                     CVectorHandle objClosest,
                                     const void CNEARTREE_FAR * coord,
                                     int resetcount) {
        double dDR, dDL;
        int nopoints;
        CVectorHandle sStack;
        CVectorHandle coords;
        CVectorHandle objs;
        void CNEARTREE_FAR * xcoord;
        void CNEARTREE_FAR * xobj;
        CNearTreeNodeHandle pt;
        double dist, drat;
        enum  { left, right, end } eDir;
        
        eDir = left; /* examine the left nodes first */
        
        nopoints = 1;
        dist = 0.;
        
        dDR = dDL = DBL_MAX;
        
        if (dRadius < 0.) return 1;
        
        if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        
        if ( treehandle->m_DelayedIndices ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }
        
        if (!(treehandle->m_iflags & CNEARTREE_NORM_L2LAZY)) return CNEARTREE_BAD_ARGUMENT;
        
        drat = sqrt((double)(treehandle->m_szdimension));
        
        pt = treehandle->m_ptTree;
        
        coords=treehandle->m_CoordStore;
        
        objs=treehandle->m_ObjectStore;
        
        if (!pt) return CNEARTREE_BAD_ARGUMENT;
        
        if (!(pt->m_iflags&CNEARTREE_FLAG_LEFT_DATA)) return CNEARTREE_NOT_FOUND;
        
#ifdef CNEARTREE_INSTRUMENTED
        (treehandle->m_NodeVisits)++;
#endif                                               
        
        CVectorCreate(&sStack,sizeof(CNearTreeNodeHandle),10);
        
        /* clear the contents of the return vector so that things don't accidentally accumulate */
        if (resetcount) {
            if (coordClosest) CVectorClear( coordClosest );
            if (objClosest) CVectorClear( objClosest );
        }
        
        while (!(eDir == end && CVectorSize(sStack) == 0)) {
            
            if ( eDir == right ) {
                dDR = DBL_MAX;
                if ((pt->m_iflags)&CNEARTREE_FLAG_RIGHT_DATA) {
                    CNTM_DistL1(dDR,treehandle, 
                                (void CNEARTREE_FAR *)coord, 
                                CVectorElementAt(coords,pt->m_indexRight));
                    if (dDR <= dRadius ){
                        nopoints = 0;
                        if (coordClosest) {
                            xcoord = CVectorElementAt(coords,pt->m_indexRight);
                            CVectorAddElement(coordClosest,&xcoord);
                        }
                        if (objClosest) {
                            xobj = CVectorElementAt(objs,pt->m_indexRight);
                            CVectorAddElement(objClosest,&xobj);
                        }
                    } else if (dDR <= dRadius*drat) {
                        CNTM_DistL2(dist,treehandle, 
                                    (void CNEARTREE_FAR *)coord,
                                    CVectorElementAt(coords,pt->m_indexRight));
                        if (dist <= dRadius) {
                            nopoints = 0;
                            if (coordClosest) {
                                xcoord = CVectorElementAt(coords,pt->m_indexRight);
                                CVectorAddElement(coordClosest,&xcoord);
                            }
                            if (objClosest) {
                                xobj = CVectorElementAt(objs,pt->m_indexRight);
                                CVectorAddElement(objClosest,&xobj);
                            }
                        }
                    }
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_RIGHT_CHILD)&& 
                    (TRIANG(dDR,pt->m_dMaxRight,dRadius))){
                    /* we did the left and now we finished the right, go down */
                    pt = pt->m_pRightBranch;
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif                                                           
                    eDir = left;
                } else {
                    eDir = end;
                }
            }
            if ( eDir == left ) {
                if ((pt->m_iflags)&CNEARTREE_FLAG_LEFT_DATA) {
                    CNTM_DistL1(dDL,treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexLeft));
                    if (dDL <= dRadius) {
                         nopoints = 0;
                        if (coordClosest) {
                            xcoord = CVectorElementAt(coords,pt->m_indexLeft);
                            CVectorAddElement(coordClosest,&xcoord);
                        }
                        if (objClosest) {
                            xobj = CVectorElementAt(objs,pt->m_indexLeft);
                            CVectorAddElement(objClosest,&xobj);
                        }
                    } else if ( dDL <= dRadius*drat ) {
                        CNTM_DistL2(dist,treehandle,
                                    (void CNEARTREE_FAR *)coord,
                                    CVectorElementAt(coords,pt->m_indexLeft))
                        if (dist <= dRadius){ 
                            nopoints = 0;
                            if (coordClosest) {
                                xcoord = CVectorElementAt(coords,pt->m_indexLeft);
                                CVectorAddElement(coordClosest,&xcoord);
                            }
                            if (objClosest) {
                                xobj = CVectorElementAt(objs,pt->m_indexLeft);
                                CVectorAddElement(objClosest,&xobj);
                            }
                        }
                    }
                }
                if (pt->m_iflags&CNEARTREE_FLAG_RIGHT_DATA) {
                    CVectorAddElement(sStack,&pt);
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_LEFT_CHILD)&&
                    (TRIANG(dDL,pt->m_dMaxLeft,dRadius))){
                    pt = pt->m_pLeftBranch;
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif                                                           
                } else {
                    eDir = end;
                }
            }
            if ( eDir == end && CVectorSize(sStack) != 0 ) {
                CVectorGetElement(sStack,&pt,CVectorSize(sStack)-1);
                CVectorRemoveElement(sStack,CVectorSize(sStack)-1);
                eDir = right;
            }
            
        }    
        
        CVectorFree(&sStack);
        return  nopoints?CNEARTREE_NOT_FOUND:CNEARTREE_SUCCESS;
    }     
    
    
    int CNearTreeFindInSphere ( const CNearTreeHandle treehandle,
                               const double dRadius,
                               CVectorHandle coordClosest,
                               CVectorHandle objClosest,
                               const void CNEARTREE_FAR * coord,
                               int resetcount) {
        double dDR, dDL;
        int nopoints;
        CVectorHandle sStack;
        CVectorHandle coords;
        CVectorHandle objs;
        void CNEARTREE_FAR * xcoord;
        void CNEARTREE_FAR * xobj;
        CNearTreeNodeHandle pt;
        enum  { left, right, end } eDir;
        
        if ((treehandle->m_iflags & CNEARTREE_NORM_L2LAZY)) 
            return CNearTreeFindInSphereL2LAZY(treehandle,dRadius,coordClosest, objClosest,coord,resetcount);

        eDir = left; /* examine the left nodes first */
        
        nopoints = 1;
        
        dDR = dDL = DBL_MAX;
        
        if (dRadius < 0.) return 1;
        
        if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        
        if ( treehandle->m_DelayedIndices ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }
        
        pt = treehandle->m_ptTree;
        
        coords=treehandle->m_CoordStore;
        
        objs=treehandle->m_ObjectStore;
        
        if (!pt) return CNEARTREE_BAD_ARGUMENT;
        
        if (!(pt->m_iflags&CNEARTREE_FLAG_LEFT_DATA)) return CNEARTREE_NOT_FOUND;
        
#ifdef CNEARTREE_INSTRUMENTED
        (treehandle->m_NodeVisits)++;
#endif                                               
        
        CVectorCreate(&sStack,sizeof(CNearTreeNodeHandle),10);
        
        /* clear the contents of the return vector so that things don't accidentally accumulate */
        if (resetcount) {
            if (coordClosest) CVectorClear( coordClosest );
            if (objClosest) CVectorClear( objClosest );
        }
        
        while (!(eDir == end && CVectorSize(sStack) == 0)) {
            
            if ( eDir == right ) {
                dDR = DBL_MAX;
                if ((pt->m_iflags)&CNEARTREE_FLAG_RIGHT_DATA) {
                    dDR = CNearTreeDist(treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexRight));
                    if (dDR <= dRadius ) {
                        nopoints = 0;
                        if (coordClosest) {
                            xcoord = CVectorElementAt(coords,pt->m_indexRight);
                            CVectorAddElement(coordClosest,&xcoord);
                        }
                        if (objClosest) {
                            xobj = CVectorElementAt(objs,pt->m_indexRight);
                            CVectorAddElement(objClosest,&xobj);
                        }
                    }
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_RIGHT_CHILD)&& 
                    (TRIANG(dDR,pt->m_dMaxRight,dRadius))){
                    /* we did the left and now we finished the right, go down */
                    pt = pt->m_pRightBranch;
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif                                                           
                    eDir = left;
                } else {
                    eDir = end;
                }
            }
            if ( eDir == left ) {
                if ((pt->m_iflags)&CNEARTREE_FLAG_LEFT_DATA) {
                    dDL = CNearTreeDist(treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexLeft));
                    if (dDL <= dRadius ) {
                        nopoints = 0;
                        if (coordClosest) {
                            xcoord = CVectorElementAt(coords,pt->m_indexLeft);
                            CVectorAddElement(coordClosest,&xcoord);
                        }
                        if (objClosest) {
                            xobj = CVectorElementAt(objs,pt->m_indexLeft);
                            CVectorAddElement(objClosest,&xobj);
                        }
                    }
                }
                if (pt->m_iflags&CNEARTREE_FLAG_RIGHT_DATA) {
                    CVectorAddElement(sStack,&pt);
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_LEFT_CHILD)&&
                    (TRIANG(dDL,pt->m_dMaxLeft,dRadius))){
                    pt = pt->m_pLeftBranch;
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif                                                           
                } else {
                    eDir = end;
                }
            }
            if ( eDir == end && CVectorSize(sStack) != 0 ) {
                CVectorGetElement(sStack,&pt,CVectorSize(sStack)-1);
                CVectorRemoveElement(sStack,CVectorSize(sStack)-1);
                eDir = right;
            }
            
        }    
        
        CVectorFree(&sStack);
        return  nopoints?CNEARTREE_NOT_FOUND:CNEARTREE_SUCCESS;
    }     
    
    /*
     =======================================================================
     int CNearTreeFindTreeInSphere ( const CNearTreeHandle treehandle, 
     const double dRadius,
     CNearTreeHandle foundClosest,
     const void CNEARTREE_FAR * coord,
     int resetcount)
     
     Function to search a Neartree for the set of objects closer to some probe point, coord,
     than dRadius. This is only here so that objClosest can be cleared before starting the work.
     
     dRadius is the maximum search radius - any point farther than dRadius from the probe
     point will be ignored
     foundClosest is an existing CNearTree to which the found points will be added
     coord  is the probe point
     resetcount should be non-zero to clear found closest on entry
     return value is 0 if points were found
     
     =======================================================================
     */
    int CNearTreeFindTreeInSphere ( const CNearTreeHandle treehandle,
                                   const double dRadius,
                                   CNearTreeHandle foundClosest,
                                   const void CNEARTREE_FAR * coord,
                                   int resetcount) {
        
        CVectorHandle objs;
        CVectorHandle coords;
        size_t ii;
        int errorcode;
        
        if ( !treehandle || !coord || !foundClosest ) return CNEARTREE_BAD_ARGUMENT;
        
        if ( treehandle->m_DelayedIndices ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }
        
        /* clear the contents of the return tree so that things don't accidentally accumulate */
        if (resetcount) {
            CNearTreeClear( foundClosest );
        }
        
        if (CVectorCreate(&objs,sizeof(void *),10)) {
            return CNEARTREE_MALLOC_FAILED;
        }
        if (CVectorCreate(&coords,sizeof(void *),10)) {
            return CNEARTREE_MALLOC_FAILED;
        }
        
        CNearTreeFindInSphere(treehandle, dRadius, coords, objs, coord, 0);
        
        for (ii=0; ii<CVectorSize(coords);ii++) {
            if ((errorcode=CNearTreeInsert(foundClosest,*((void CNEARTREE_FAR * CNEARTREE_FAR *)CVectorElementAt(coords,ii)),CVectorElementAt(objs,ii))),errorcode) {
                CVectorFree(&coords);
                CVectorFree(&objs);
                return errorcode;
            }
        }
        
        errorcode = 0;
        if (CVectorSize(coords) == 0) errorcode=CNEARTREE_NOT_FOUND;
        
        CVectorFree(&coords);
        CVectorFree(&objs);
        
        return errorcode;
    }
    
    
    /*
     =======================================================================
     int CNearTreeFindOutSphere ( const CNearTreeHandle treehandle, 
     const double dRadius,
     CVectorHandle coordOutside,
     CVectorHandle objOutside,
     const void CNEARTREE_FAR * coord,
     int resetcount)
     
     Function to search a Neartree for the set of objects further from some probe point, coord,
     than dRadius. 
     
     dRadius is the maximum search radius - any point closer than dRadius from the probe
     point will be ignored
     
     coordOutside is a vector of pointers to coordinate tuples and is the 
     returned set of distant points from the probe point that can be found in the Neartree
     
     objOutside is a vector of objects and is the returned set of distant points
     from the probe point that can be found in the Neartree
     
     coord  is the probe point
     
     resetcount should be non-zero to clear objClosest on entry
     
     return value is 0 if points were found
     
     =======================================================================
     */
    
    int CNearTreeFindOutSphereL2LAZY ( const CNearTreeHandle treehandle,
                                const double dRadius,
                                CVectorHandle coordOutside,
                                CVectorHandle objOutside,
                                const void * coord,
                                int resetcount){
        double dDR, dDL;
        int nopoints;
        CVectorHandle sStack;
        CVectorHandle coords;
        CVectorHandle objs;
        void CNEARTREE_FAR * xcoord;
        void CNEARTREE_FAR * xobj;
        CNearTreeNodeHandle pt;
        double dist, drat;
        enum  { left, right, end } eDir;
        
        eDir = left; /* examine the left nodes first */
        
        nopoints = 1;
        dist = 0.;
        
        dDR = dDL = DBL_MAX;
        
        if (dRadius < 0.) return 1;
        
        if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        
        if ( treehandle->m_DelayedIndices ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }
        
        drat = sqrt((double)(treehandle->m_szdimension));
        
        pt = treehandle->m_ptTree;
        
        coords=treehandle->m_CoordStore;
        
        objs=treehandle->m_ObjectStore;
        
        if (!pt) return CNEARTREE_BAD_ARGUMENT;
        
        if (!(pt->m_iflags&CNEARTREE_FLAG_LEFT_DATA)) return CNEARTREE_NOT_FOUND;
        
#ifdef CNEARTREE_INSTRUMENTED
        (treehandle->m_NodeVisits)++;
#endif                                       
        
        CVectorCreate(&sStack,sizeof(CNearTreeNodeHandle),10);
        
        /* clear the contents of the return vector so that things don't accidentally accumulate */
        if (resetcount) {
            if (coordOutside) CVectorClear( coordOutside );
            if (objOutside) CVectorClear( objOutside );
        }
        
        while (!(eDir == end && CVectorSize(sStack) == 0)) {
            
            if ( eDir == right ) {
                dDR = DBL_MAX;
                if ((pt->m_iflags)&CNEARTREE_FLAG_RIGHT_DATA) {
                    CNTM_DistL1(dDR,treehandle, 
                                (void CNEARTREE_FAR *)coord, 
                                CVectorElementAt(coords,pt->m_indexRight));
                    if (dDR >= dRadius*drat ) {
                        nopoints = 0;
                        if (coordOutside) {
                            xcoord = CVectorElementAt(coords,pt->m_indexRight);
                            CVectorAddElement(coordOutside,&xcoord);
                        }
                        if (objOutside) {
                            xobj = CVectorElementAt(objs,pt->m_indexRight);
                            CVectorAddElement(objOutside,&xobj);
                        }
                    } else if (dDR >= dRadius ) {
                        CNTM_DistL2(dist,treehandle, 
                                    (void CNEARTREE_FAR *)coord, 
                                    CVectorElementAt(coords,pt->m_indexRight));
                        if (dist >= dRadius) {
                            nopoints = 0;
                            if (coordOutside) {
                                xcoord = CVectorElementAt(coords,pt->m_indexRight);
                                CVectorAddElement(coordOutside,&xcoord);
                            }
                            if (objOutside) {
                                xobj = CVectorElementAt(objs,pt->m_indexRight);
                                CVectorAddElement(objOutside,&xobj);
                            }
                        }
                    }
                    
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_RIGHT_CHILD)&& 
                    (TRIANG(dRadius,dDR,pt->m_dMaxRight))){
                    /* we did the left and now we finished the right, go down */
                    pt = pt->m_pRightBranch;
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif                                                           
                    eDir = left;
                } else {
                    eDir = end;
                }
            }
            if ( eDir == left ) {
                if ((pt->m_iflags)&CNEARTREE_FLAG_LEFT_DATA) {
                    CNTM_DistL1(dDL,treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexLeft));
                    if (dDL >= dRadius*drat ) {
                        nopoints = 0;
                        if (coordOutside) {
                            xcoord = CVectorElementAt(coords,pt->m_indexLeft);
                            CVectorAddElement(coordOutside,&xcoord);
                        }
                        if (objOutside) {
                            xobj = CVectorElementAt(objs,pt->m_indexLeft);
                            CVectorAddElement(objOutside,&xobj);
                        }
                    } else if (dDL >= dRadius) {
                        CNTM_DistL2(dist,treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexLeft));
                        if (dist >= dRadius)  {
                            nopoints = 0;
                            if (coordOutside) {
                                xcoord = CVectorElementAt(coords,pt->m_indexLeft);
                                CVectorAddElement(coordOutside,&xcoord);
                            }
                            if (objOutside) {
                                xobj = CVectorElementAt(objs,pt->m_indexLeft);
                                CVectorAddElement(objOutside,&xobj);
                            }                            
                        }
                    }
                }
                if (pt->m_iflags&CNEARTREE_FLAG_RIGHT_DATA) {
                    CVectorAddElement(sStack,&pt);
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_LEFT_CHILD)&&
                    (TRIANG(dRadius,dDL,pt->m_dMaxLeft))){
                    pt = pt->m_pLeftBranch;
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif                                                           
                } else {
                    eDir = end;
                }
            }
            if ( eDir == end && CVectorSize(sStack) != 0 ) {
                CVectorGetElement(sStack,&pt,CVectorSize(sStack)-1);
                CVectorRemoveElement(sStack,CVectorSize(sStack)-1);
                eDir = right;
            }
            
        }    
        
        CVectorFree(&sStack);
        return  nopoints?CNEARTREE_NOT_FOUND:CNEARTREE_SUCCESS;
    }     
    
    
    int CNearTreeFindOutSphere ( const CNearTreeHandle treehandle,
                                const double dRadius,
                                CVectorHandle coordOutside,
                                CVectorHandle objOutside,
                                const void * coord,
                                int resetcount){
        double dDR, dDL;
        int nopoints;
        CVectorHandle sStack;
        CVectorHandle coords;
        CVectorHandle objs;
        void CNEARTREE_FAR * xcoord;
        void CNEARTREE_FAR * xobj;
        CNearTreeNodeHandle pt;
        enum  { left, right, end } eDir;
        
        eDir = left; /* examine the left nodes first */
        
        nopoints = 1;
        
        dDR = dDL = DBL_MAX;
        
        if (dRadius < 0.) return 1;
        
        if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        
        if ((treehandle->m_iflags & CNEARTREE_NORM_L2LAZY)) 
            return CNearTreeFindOutSphereL2LAZY(treehandle,dRadius,coordOutside,objOutside,coord,resetcount);
        
        if ( treehandle->m_DelayedIndices ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }
        
        pt = treehandle->m_ptTree;
        
        coords=treehandle->m_CoordStore;
        
        objs=treehandle->m_ObjectStore;
        
        if (!pt) return CNEARTREE_BAD_ARGUMENT;
        
        if (!(pt->m_iflags&CNEARTREE_FLAG_LEFT_DATA)) return CNEARTREE_NOT_FOUND;

#ifdef CNEARTREE_INSTRUMENTED
        (treehandle->m_NodeVisits)++;
#endif                                       
        
        CVectorCreate(&sStack,sizeof(CNearTreeNodeHandle),10);
        
        /* clear the contents of the return vector so that things don't accidentally accumulate */
        if (resetcount) {
            if (coordOutside) CVectorClear( coordOutside );
            if (objOutside) CVectorClear( objOutside );
        }
        
        while (!(eDir == end && CVectorSize(sStack) == 0)) {
            
            if ( eDir == right ) {
                dDR = DBL_MAX;
                if ((pt->m_iflags)&CNEARTREE_FLAG_RIGHT_DATA) {
                    dDR = CNearTreeDist(treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexRight));
                    if (dDR >= dRadius ) {
                        nopoints = 0;
                        if (coordOutside) {
                            xcoord = CVectorElementAt(coords,pt->m_indexRight);
                            CVectorAddElement(coordOutside,&xcoord);
                        }
                        if (objOutside) {
                            xobj = CVectorElementAt(objs,pt->m_indexRight);
                            CVectorAddElement(objOutside,&xobj);
                        }
                    }
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_RIGHT_CHILD)&& 
                    (TRIANG(dRadius,dDR,pt->m_dMaxRight))){
                    /* we did the left and now we finished the right, go down */
                    pt = pt->m_pRightBranch;
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif                                                           
                    eDir = left;
                } else {
                    eDir = end;
                }
            }
            if ( eDir == left ) {
                if ((pt->m_iflags)&CNEARTREE_FLAG_LEFT_DATA) {
                    dDL = CNearTreeDist(treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexLeft));
                    if (dDL >= dRadius ) {
                        nopoints = 0;
                        if (coordOutside) {
                            xcoord = CVectorElementAt(coords,pt->m_indexLeft);
                            CVectorAddElement(coordOutside,&xcoord);
                        }
                        if (objOutside) {
                            xobj = CVectorElementAt(objs,pt->m_indexLeft);
                            CVectorAddElement(objOutside,&xobj);
                        }
                    }
                }
                if (pt->m_iflags&CNEARTREE_FLAG_RIGHT_DATA) {
                    CVectorAddElement(sStack,&pt);
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_LEFT_CHILD)&&
                    (TRIANG(dRadius,dDL,pt->m_dMaxLeft))){
                    pt = pt->m_pLeftBranch;
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif                                                           
                } else {
                    eDir = end;
                }
            }
            if ( eDir == end && CVectorSize(sStack) != 0 ) {
                CVectorGetElement(sStack,&pt,CVectorSize(sStack)-1);
                CVectorRemoveElement(sStack,CVectorSize(sStack)-1);
                eDir = right;
            }
            
        }    
        
        CVectorFree(&sStack);
        return  nopoints?CNEARTREE_NOT_FOUND:CNEARTREE_SUCCESS;
    }     
    
    /*
     =======================================================================
     int CNearTreeFindTreeOutSphere ( const CNearTreeHandle treehandle, 
     const double dRadius,
     CNearTreeHandle foundOutside,
     const void CNEARTREE_FAR * coord,
     int resetcount)
     
     Function to search a Neartree for the set of objects further from some probe point, coord,
     than dRadius.
     
     dRadius is the maximum search radius - any point closer than dRadius from the probe
     point will be ignored
     foundOutside is an existing CNearTree to which the found points will be added
     coord  is the probe point
     resetcount should be non-zero to clear found closest on entry
     return value is 0 if points were found
     
     =======================================================================
     */
    int CNearTreeFindTreeOutSphere ( const CNearTreeHandle treehandle,
                                    const double dRadius,
                                    CNearTreeHandle foundOutside,
                                    const void CNEARTREE_FAR * coord,
                                    int resetcount){
        
        CVectorHandle objs;
        CVectorHandle coords;
        size_t ii;
        int errorcode;
        
        if ( !treehandle || !coord || !foundOutside ) return CNEARTREE_BAD_ARGUMENT;

        if ( treehandle->m_DelayedIndices ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }
        
        /* clear the contents of the return tree so that things don't accidentally accumulate */
        if (resetcount) {
            CNearTreeClear( foundOutside );
        }
        
        if (CVectorCreate(&objs,sizeof(void *),10)) {
            return CNEARTREE_MALLOC_FAILED;
        }
        if (CVectorCreate(&coords,sizeof(void *),10)) {
            return CNEARTREE_MALLOC_FAILED;
        }
        
        CNearTreeFindOutSphere(treehandle, dRadius, coords, objs, coord, 0);
        
        for (ii=0; ii<CVectorSize(coords);ii++) {
            if ((errorcode=CNearTreeInsert(foundOutside,*((void CNEARTREE_FAR * CNEARTREE_FAR *)CVectorElementAt(coords,ii)),CVectorElementAt(objs,ii))),errorcode) {
                CVectorFree(&coords);
                CVectorFree(&objs);
                return errorcode;
            }
        }
        errorcode = 0;
        if (CVectorSize(coords) == 0) errorcode=CNEARTREE_NOT_FOUND;
        
        CVectorFree(&coords);
        CVectorFree(&objs);
        
        return errorcode;
        
        
    }
    
    /*
     =======================================================================
     int CNearTreeFindInAnnulus ( const CNearTreeHandle treehandle, 
     const double dRadiusInner,
     const double dRadiusOuter,
     CVectorHandle coordInRing,
     CVectorHandle objInRing,
     const void CNEARTREE_FAR * coord,
     int resetcount)
     
     Function to search a Neartree for the set of objects closer to some probe point, coord,
     than dRadiusOuter and further than dRadiusInner. 
     
     dRadiusInner is the minimum search radius - any point closer than dRadius from the probe
     point will be ignored
     dRadiusOuter is the maximum search radius - any point further than dRadius from the probe
     point will be ignored
     coordInRing is a vector of pointers to coordinate tuples of nearest points
     objInRing is a vector of objects and is the returned set of nearest points
     to the probe point that can be found in the Neartree
     coord  is the probe point
     resetcount should be non-zero to clear coordInRing and objInRing on entry
     return value is 0 if points were found
     
     =======================================================================
     */
    
    int CNearTreeFindInAnnulusL2LAZY ( const CNearTreeHandle treehandle,
                                const double dRadiusInner,
                                const double dRadiusOuter,
                                CVectorHandle coordInRing,
                                CVectorHandle objInRing,
                                const void CNEARTREE_FAR * coord,
                                int resetcount) {
        double dDR, dDL;
        int nopoints;
        CVectorHandle sStack;
        CVectorHandle coords;
        CVectorHandle objs;
        void CNEARTREE_FAR * xcoord;
        void CNEARTREE_FAR * xobj;
        CNearTreeNodeHandle pt;
        double dist, drat;
        enum  { left, right, end } eDir;
        
        eDir = left; /* examine the left nodes first */
        
        nopoints = 1;
        dist = 0.;
        
        dDR = dDL = DBL_MAX;
        
        if (dRadiusInner < 0. || dRadiusOuter < 0. || dRadiusInner > dRadiusOuter) return CNEARTREE_BAD_ARGUMENT;
        
        if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        
        if ( treehandle->m_DelayedIndices ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }
        
        drat = sqrt((double)(treehandle->m_szdimension));
        
        pt = treehandle->m_ptTree;
        
        coords=treehandle->m_CoordStore;
        
        objs=treehandle->m_ObjectStore;
        
        if (!pt) return CNEARTREE_BAD_ARGUMENT;
        
        if (!(pt->m_iflags&CNEARTREE_FLAG_LEFT_DATA)) return CNEARTREE_NOT_FOUND;
        
#ifdef CNEARTREE_INSTRUMENTED
        (treehandle->m_NodeVisits)++;
#endif                                               
        
        CVectorCreate(&sStack,sizeof(CNearTreeNodeHandle),10);
        
        /* clear the contents of the return vector so that things don't accidentally accumulate */
        if (resetcount) {
            if (coordInRing) CVectorClear( coordInRing );
            if (objInRing) CVectorClear( objInRing );
        }
        
        while (!(eDir == end && CVectorSize(sStack) == 0)) {
            
            if ( eDir == right ) {
                dDR = DBL_MAX;
                if ((pt->m_iflags)&CNEARTREE_FLAG_RIGHT_DATA) {
                    CNTM_DistL1(dDR,treehandle, 
                                (void CNEARTREE_FAR *)coord, 
                                CVectorElementAt(coords,pt->m_indexRight));                    
                    if (dDR <= dRadiusOuter && dDR >= dRadiusInner*drat ) {
                        nopoints = 0;
                        if (coordInRing) {
                            xcoord = CVectorElementAt(coords,pt->m_indexRight);
                            CVectorAddElement(coordInRing,&xcoord);
                        }
                        if (objInRing) {
                            xobj = CVectorElementAt(objs,pt->m_indexRight);
                            CVectorAddElement(objInRing,&xobj);
                        }
                    } else if (dDR <= dRadiusOuter*drat && dDR >= dRadiusInner ) {
                        CNTM_DistL2(dist,treehandle, 
                                    (void CNEARTREE_FAR *)coord,
                                    CVectorElementAt(coords,pt->m_indexRight));
                        if (dist <= dRadiusOuter && dist >= dRadiusInner ) {
                            nopoints = 0;
                            if (coordInRing) {
                                xcoord = CVectorElementAt(coords,pt->m_indexRight);
                                CVectorAddElement(coordInRing,&xcoord);
                            }
                            if (objInRing) {
                                xobj = CVectorElementAt(objs,pt->m_indexRight);
                                CVectorAddElement(objInRing,&xobj);
                            }
                            
                        }
                    }
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_RIGHT_CHILD)&& 
                    (TRIANG(dDR,pt->m_dMaxRight,dRadiusOuter))&&
                    (TRIANG(dRadiusInner,dDR,pt->m_dMaxRight))){
                    /* we did the left and now we finished the right, go down */
                    pt = pt->m_pRightBranch;
                    
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif                                       
                    eDir = left;
                } else {
                    eDir = end;
                }
            }
            if ( eDir == left ) {
                if ((pt->m_iflags)&CNEARTREE_FLAG_LEFT_DATA) {
                    dDL = CNearTreeDist(treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexLeft));
                    if (dDL <= dRadiusOuter && dDL >= dRadiusInner*drat) {
                        nopoints = 0;
                        if (coordInRing) {
                            xcoord = CVectorElementAt(coords,pt->m_indexLeft);
                            CVectorAddElement(coordInRing,&xcoord);
                        }
                        if (objInRing) {
                            xobj = CVectorElementAt(objs,pt->m_indexLeft);
                            CVectorAddElement(objInRing,&xobj);
                        }
                    } else if (dDL <= dRadiusOuter*drat && dDL >= dRadiusInner ) {
                        CNTM_DistL2(dist,treehandle, 
                                    (void CNEARTREE_FAR *)coord,
                                    CVectorElementAt(coords,pt->m_indexLeft));
                        if (dist <= dRadiusOuter && dist >= dRadiusInner ) {
                            nopoints = 0;
                            if (coordInRing) {
                                xcoord = CVectorElementAt(coords,pt->m_indexLeft);
                                CVectorAddElement(coordInRing,&xcoord);
                            }
                            if (objInRing) {
                                xobj = CVectorElementAt(objs,pt->m_indexLeft);
                                CVectorAddElement(objInRing,&xobj);
                            }
                            
                        }
                    }
                    
                }
                if (pt->m_iflags&CNEARTREE_FLAG_RIGHT_DATA) {
                    CVectorAddElement(sStack,&pt);
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_LEFT_CHILD)&&
                    (TRIANG(dDL,pt->m_dMaxLeft,dRadiusOuter))&&
                    (TRIANG(dRadiusInner,dDL,pt->m_dMaxLeft))){
                    pt = pt->m_pLeftBranch;
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif                                                           
                } else {
                    eDir = end;
                }
            }
            if ( eDir == end && CVectorSize(sStack) != 0 ) {
                CVectorGetElement(sStack,&pt,CVectorSize(sStack)-1);
                CVectorRemoveElement(sStack,CVectorSize(sStack)-1);
                eDir = right;
            }
            
        }    
        
        CVectorFree(&sStack);
        return  nopoints?CNEARTREE_NOT_FOUND:CNEARTREE_SUCCESS;
    }     
    
    
    int CNearTreeFindInAnnulus ( const CNearTreeHandle treehandle,
                                const double dRadiusInner,
                                const double dRadiusOuter,
                                CVectorHandle coordInRing,
                                CVectorHandle objInRing,
                                const void CNEARTREE_FAR * coord,
                                int resetcount) {
        double dDR, dDL;
        int nopoints;
        CVectorHandle sStack;
        CVectorHandle coords;
        CVectorHandle objs;
        void CNEARTREE_FAR * xcoord;
        void CNEARTREE_FAR * xobj;
        CNearTreeNodeHandle pt;
        enum  { left, right, end } eDir;
        
        eDir = left; /* examine the left nodes first */
        
        nopoints = 1;
        
        dDR = dDL = DBL_MAX;
        
        if (dRadiusInner < 0. || dRadiusOuter < 0. || dRadiusInner > dRadiusOuter) return CNEARTREE_BAD_ARGUMENT;
        
        if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        
        if ((treehandle->m_iflags & CNEARTREE_NORM_L2LAZY)) 
            return CNearTreeFindInAnnulusL2LAZY(treehandle,dRadiusInner,dRadiusOuter,coordInRing,objInRing,coord,resetcount);
        
        if ( treehandle->m_DelayedIndices ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }
        
        pt = treehandle->m_ptTree;
        
        coords=treehandle->m_CoordStore;
        
        objs=treehandle->m_ObjectStore;
        
        if (!pt) return CNEARTREE_BAD_ARGUMENT;
        
        if (!(pt->m_iflags&CNEARTREE_FLAG_LEFT_DATA)) return CNEARTREE_NOT_FOUND;
        
#ifdef CNEARTREE_INSTRUMENTED
        (treehandle->m_NodeVisits)++;
#endif                                               
        
        CVectorCreate(&sStack,sizeof(CNearTreeNodeHandle),10);
        
        /* clear the contents of the return vector so that things don't accidentally accumulate */
        if (resetcount) {
            if (coordInRing) CVectorClear( coordInRing );
            if (objInRing) CVectorClear( objInRing );
        }
        
        while (!(eDir == end && CVectorSize(sStack) == 0)) {
            
            if ( eDir == right ) {
                dDR = DBL_MAX;
                if ((pt->m_iflags)&CNEARTREE_FLAG_RIGHT_DATA) {
                    dDR = CNearTreeDist(treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexRight));
                    if (dDR <= dRadiusOuter && dDR >= dRadiusInner ) {
                        nopoints = 0;
                        if (coordInRing) {
                            xcoord = CVectorElementAt(coords,pt->m_indexRight);
                            CVectorAddElement(coordInRing,&xcoord);
                        }
                        if (objInRing) {
                            xobj = CVectorElementAt(objs,pt->m_indexRight);
                            CVectorAddElement(objInRing,&xobj);
                        }
                    }
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_RIGHT_CHILD)&& 
                    (TRIANG(dDR,pt->m_dMaxRight,dRadiusOuter))&&
                    (TRIANG(dRadiusInner,dDR,pt->m_dMaxRight))){
                    /* we did the left and now we finished the right, go down */
                    pt = pt->m_pRightBranch;

#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif                                       
                    eDir = left;
                } else {
                    eDir = end;
                }
            }
            if ( eDir == left ) {
                if ((pt->m_iflags)&CNEARTREE_FLAG_LEFT_DATA) {
                    dDL = CNearTreeDist(treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexLeft));
                    if (dDL <= dRadiusOuter && dDL >= dRadiusInner) {
                        nopoints = 0;
                        if (coordInRing) {
                            xcoord = CVectorElementAt(coords,pt->m_indexLeft);
                            CVectorAddElement(coordInRing,&xcoord);
                        }
                        if (objInRing) {
                            xobj = CVectorElementAt(objs,pt->m_indexLeft);
                            CVectorAddElement(objInRing,&xobj);
                        }
                    }
                }
                if (pt->m_iflags&CNEARTREE_FLAG_RIGHT_DATA) {
                    CVectorAddElement(sStack,&pt);
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_LEFT_CHILD)&&
                    (TRIANG(dDL,pt->m_dMaxLeft,dRadiusOuter))&&
                    (TRIANG(dRadiusInner,dDL,pt->m_dMaxLeft))){
                    pt = pt->m_pLeftBranch;
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif                                                           
                } else {
                    eDir = end;
                }
            }
            if ( eDir == end && CVectorSize(sStack) != 0 ) {
                CVectorGetElement(sStack,&pt,CVectorSize(sStack)-1);
                CVectorRemoveElement(sStack,CVectorSize(sStack)-1);
                eDir = right;
            }
            
        }    
        
        CVectorFree(&sStack);
        return  nopoints?CNEARTREE_NOT_FOUND:CNEARTREE_SUCCESS;
    }     
    
    /*
     =======================================================================
     int CNearTreeFindTreeInAnnulus ( const CNearTreeHandle treehandle, 
     const double dRadiusInner,
     const double dRadiusOuter,
     CNearTreeHandle foundInRing,
     const void CNEARTREE_FAR * coord,
     int resetcount)
     
     Function to search a Neartree for the set of objects closer to some probe point, coord,
     than dRadiusOuter and further than dRadiusInner. 
     
     dRadiusInner is the minimum search radius - any point closer than dRadius from the probe
     point will be ignored
     dRadiusOuter is the maximum search radius - any point further than dRadius from the probe
     point will be ignored
     foundInRing is an existing CNearTree to which the found points will be added
     coord  is the probe point
     resetcount should be non-zero to clear found InRing on entry
     return value is 0 if points were found
     
     =======================================================================
     */
    int CNearTreeFindTreeInAnnulus ( const CNearTreeHandle treehandle,
                                    const double dRadiusInner,
                                    const double dRadiusOuter,
                                    CNearTreeHandle foundInRing,
                                    const void CNEARTREE_FAR * coord,
                                    int resetcount) {
        
        CVectorHandle objs;
        CVectorHandle coords;
        size_t ii;
        int errorcode;
        
        if ( !treehandle || !coord || !foundInRing ) return CNEARTREE_BAD_ARGUMENT;
        
        if ( treehandle->m_DelayedIndices ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }
        
        /* clear the contents of the return tree so that things don't accidentally accumulate */
        if (resetcount) {
            CNearTreeClear( foundInRing );
        }
        
        if (CVectorCreate(&objs,sizeof(void *),10)) {
            return CNEARTREE_MALLOC_FAILED;
        }
        if (CVectorCreate(&coords,sizeof(void *),10)) {
            return CNEARTREE_MALLOC_FAILED;
        }
        
        CNearTreeFindInAnnulus(treehandle, dRadiusInner, dRadiusOuter, coords, objs, coord, 0);
        
        for (ii=0; ii<CVectorSize(coords);ii++) {
            if ((errorcode=CNearTreeInsert(foundInRing,*((void CNEARTREE_FAR * CNEARTREE_FAR *)CVectorElementAt(coords,ii)),CVectorElementAt(objs,ii))),errorcode) {
                CVectorFree(&coords);
                CVectorFree(&objs);
                return errorcode;
            }
        }      
        errorcode = 0;
        if (CVectorSize(coords) == 0) errorcode=CNEARTREE_NOT_FOUND;
        
        CVectorFree(&coords);
        CVectorFree(&objs);
        
        return errorcode;
        
    }
    
    /*  int CNearTreeSortIn(CVectorHandle metrics,
     CVectorHandle indices,
     double metric,
     size_t index,
     size_t k);
     
     CNearTreeSortIn inserts a new metric and index into the vectors
     metrics and indices, sorted on non-decreasing metric,
     with the size of the vectors capped at k, or uncapped if k = 0;
     
     */
    
    int CNearTreeSortIn(CVectorHandle metrics,
                        CVectorHandle indices,
                        double metric,
                        size_t index,
                        size_t k) {
        
        return CNearTreeSortIn2(metrics, indices, metric, index, k, 0, 0, 0.);
        
    }

    
    /*  int CNearTreeSortIn2(CVectorHandle metrics,
     CVectorHandle indices, 
     double metric, 
     size_t index,
     size_t k
     int shell,
     size_t sSizeCur,
     double dRadiusCur);
     
     CNearTreeSortIn2 inserts a new metric and index into the vectors
     metrics and indices, sorted on non-decreasing metric,
     with the size of the vectors capped at k, or uncapped if k = 0;
     
     if shell is non-zero only the results at the minimal distance are
     retained
     
     */ 
    
    int CNearTreeSortIn2(CVectorHandle metrics,
                         CVectorHandle indices,
                         double metric,
                         size_t index,
                         size_t k,
                         int shell,
                         size_t sSizeCur,
                         double dRadiusCur) {
        
        double localmetric=metric;
        double tempmetric;
        void * localmetrics;
        size_t localindex=index;
        size_t tempindex;
        void * localindices;
        size_t low, mid, high;
        
        if (!metrics
            || !indices 
            || CVectorSize(metrics) != CVectorSize(indices)) return CNEARTREE_BAD_ARGUMENT;
        
        if (CVectorSize(metrics) <= sSizeCur) {
            CVectorAddElement(metrics,&localmetric);
            CVectorAddElement(indices,&localindex);
            return CNEARTREE_SUCCESS;
        }
        
        low = 0;
        if (shell && sSizeCur > 0) low = sSizeCur;
        high = CVectorSize(metrics)-1;
        CVectorGetElementptr(metrics, &localmetrics,0);
        CVectorGetElementptr(indices, &localindices,0);
        
        if (( shell && localmetric < ((double *)localmetrics)[low] )
            || (localmetric < ((double *)localmetrics)[low] && k == 1)){
            CVectorSetSize(metrics,sSizeCur);
            CVectorSetSize(indices,sSizeCur);
            high = 0;
            CVectorAddElement(metrics,&localmetric);
            CVectorAddElement(indices,&localindex);
        
        } else if (localmetric <= ((double *)localmetrics)[low]) {
            
            tempmetric =  ((double *)localmetrics)[high];
            tempindex =  ((size_t *)localindices)[high];
            for (mid = high; mid > low; mid--) {
                ((double *)localmetrics)[mid] = ((double *)localmetrics)[mid-1];
                ((size_t *)localindices)[mid] = ((size_t *)localindices)[mid-1];
            }
            ((double *)localmetrics)[low] = localmetric;
            ((size_t *)localindices)[low] = localindex;
            CVectorSetFlags(metrics,0);
            CVectorSetFlags(indices,0);
            if (k==0 || CVectorSize(metrics)-sSizeCur < k) {
                CVectorAddElement(metrics,&tempmetric);
                CVectorAddElement(indices,&tempindex); 
            }
            
        } else if (localmetric >= ((double *)localmetrics)[high]) {
            CVectorSetFlags(metrics,0);
            CVectorSetFlags(indices,0);
            if ((!shell && (k==0 || CVectorSize(metrics) < k))
                || (shell
                    && localmetric == ((double *)localmetrics)[high]
                    && (k==0 || CVectorSize(metrics)-sSizeCur < k))){
                CVectorAddElement(metrics,&localmetric);
                CVectorAddElement(indices,&localindex); 
            }
            
            
        } else {
            
            /* ((double *)localmetrics)[low] < localmetric < ((double *)localmetrics)[high]*/
            
            while (low < high-1) {
                mid = (low+high)/2;
                if (localmetric == ((double *)localmetrics)[mid]) {
                    low = mid;
                    break;
                }
                if (localmetric < ((double *)localmetrics)[mid]) {
                    high = mid;
                } else {
                    low = mid;
                }      
            }
            
            /* Insert the new item just above low */
            
            if (low < high) {
                tempmetric =  ((double *)localmetrics)[high];
                tempindex =  ((size_t *)localindices)[high];
                for (mid = high; mid > low+1; mid--) {
                    ((double *)localmetrics)[mid] = ((double *)localmetrics)[mid-1];
                    ((size_t *)localindices)[mid] = ((size_t *)localindices)[mid-1];
                }
                ((double *)localmetrics)[low+1] = localmetric;
                ((size_t *)localindices)[low+1] = localindex;
            } else {
                tempmetric = localmetric;
                tempindex = localindex;
            }
            CVectorSetFlags(metrics,0);
            CVectorSetFlags(indices,0);
            if (k==0 || CVectorSize(metrics) < k) { 
                CVectorAddElement(metrics,&tempmetric);
                CVectorAddElement(indices,&tempindex); 
            }
            
        }
        
        return CNEARTREE_SUCCESS;
        
    }
    
    /*
     =======================================================================
     int CNearTreeFindKNearest ( const CNearTreeHandle treehandle,
     const size_t k,
     const double dRadius,
     CVectorHandle coordClosest,
     CVectorHandle objClosest,
     const void * coord,
     int resetcount);
     
     Function to search a Neartree for the set of objects closer to some probe point, coord,
     than dRadius. 
     
     k is the maximum number of closest neighbors to return.  Finds this many if passible.
     
     dRadius is the maximum search radius - any point further than dRadius from the probe
     point will be ignored
     
     coordClosest is a vector of pointers to coordinate tuples and is the 
     returned set of farthest points from the probe point that can be found in the Neartree
     
     objClosest is a vector of objects and is the returned set of farthest points
     from the probe point that can be found in the Neartree
     
     coord  is the probe point
     
     resetcount should be non-zero to clear objClosest on entry
     
     return value is 0 if points were found
     
     =======================================================================
     */
    
    int CNearTreeFindKNearest (const CNearTreeHandle treehandle,
                               const size_t k,
                               const double dRadius,
                               CVectorHandle coordClosest,
                               CVectorHandle objClosest,
                               const void * coord,
                               int resetcount){
        
        CVectorHandle distances;
        CVectorHandle indices;
        CVectorHandle coords;
        CVectorHandle objs;
        void CNEARTREE_FAR * xcoord;
        void CNEARTREE_FAR * xobj;
        size_t ii, index;
        int lresult;
        
        coords=treehandle->m_CoordStore;
        objs=treehandle->m_ObjectStore;

        /* prepare coordClosest and objClosest */
        
        if (resetcount) {
            if (coordClosest) CVectorClear( coordClosest );
            if (objClosest) CVectorClear( objClosest );
        }
        
        if (!CNearTreeZeroIfEmpty(treehandle)) return CNEARTREE_NOT_FOUND;
        
        if (CVectorCreate(&distances,sizeof(double),k) ||
            CVectorCreate(&indices,sizeof(size_t),k))
            
            return CNEARTREE_MALLOC_FAILED;

        
        if ((treehandle->m_iflags&CNTF_SKNN)) {
            lresult = CNearTreeFindKNearest_Sphere ( treehandle,
                                                 k,
                                                 dRadius,
                                                 distances,
                                                 indices,
                                                 coord);
        } else {
            lresult = CNearTreeFindKNearest_Annular ( treehandle,
                                                  k,
                                                  dRadius,
                                                  distances,
                                                  indices,
                                                  coord);
        }
        
        /* Then transfer out the results to coordClosest and objClosest */
        
        for (ii = 0; ii < CVectorSize(indices) && ii < k; ii++) {
            CVectorGetElement(indices,(void *)&index,ii);
            if (coordClosest)  {
                xcoord = CVectorElementAt(coords,index);
                CVectorAddElement(coordClosest,&xcoord);
            }
            if (objClosest) {
                xobj = CVectorElementAt(objs,index);
                CVectorAddElement(objClosest,&xobj);
            }
            
        }
        
        return lresult;
        
    }

#include <stdio.h>
    
    int CNearTreeFindKNearest_Annular (const CNearTreeHandle treehandle,
                                       const size_t k,
                                       const double dRadius,
                                       CVectorHandle distances,
                                       CVectorHandle indices,
                                       const void * coord){
        
        double dRadiusInner = 0;
        double dRadiusOuterSave;
        double dRadiusOuter;
        double radlist[CNEARTREE_DIMSAMPLES];
        double dimlist[CNEARTREE_DIMSAMPLES-1];
        double dimest = 1.;
        double foundatrad[CNEARTREE_DIMSAMPLES];
        double drat;
        int numrad;
        int shell, closed;
        int error;
        long lFound;
        int l2lazy;
        size_t size;
        
        drat = 1.;        
        l2lazy = 0;
        if ((treehandle->m_iflags&CNEARTREE_NORM_L2LAZY)) {
            /* compute the worst case ratio of the L1 norm to the L2 norm */
            drat = sqrt((double)(treehandle->m_szdimension));
            l2lazy = 1;
        }
        
        
        dRadiusOuter = treehandle->m_SumSpacings/sqrt((double)(1+CVectorSize(treehandle->m_ObjectStore)));
        if (dRadiusOuter <= dRadiusInner) dRadiusOuter = dRadiusInner+1.;
        if (dRadiusOuter > dRadius) dRadiusOuter = dRadius;
        numrad = 0;
        shell = 1;
        closed = 1;

        /*fprintf (stderr,
                 "Entering CNearTreeFindKNearest_Annular, dRadius %g, dRadiusInner %g, dRadiusOuter %g\n",
                 dRadius, dRadiusInner, dRadiusOuter);*/
        
        /* First find the nearest k inner shell */
        do {
            dRadiusOuterSave = dRadiusOuter;
            size = CVectorSize(distances);
            /*fprintf(stderr,"k %d, dRadius %g, dRadiusInner, dRadiusOuter %g %g\n", (int)k,
                    (double)dRadius, (double)(dRadiusInner), (double)(dRadiusOuter));*/
            error = CNearTreeFindKNearInAnnulus(treehandle,k-size,shell,closed,
                                                dRadiusInner,
                                                dRadiusOuter,
                                                distances,
                                                indices,
                                                coord);
            lFound = CVectorSize(distances)-size;
            if (lFound > 0) {
                CVectorGetElement(distances,&dRadiusOuter,CVectorSize(distances)-1);
                /*fprintf (stderr,
                         "Ran CNearTreeFindKNearInAnnulus, lFound %ld , dRadiusInner %g, dRadiusOuter %g\n",
                         (long)lFound, dRadiusInner, dRadiusOuter);*/
                break;
            }
            dRadiusOuter = dRadiusOuterSave+(dRadiusOuterSave-dRadiusInner)*2.;
            dRadiusInner = dRadiusOuterSave;
            if (dRadiusOuter <= dRadiusInner) dRadiusOuter = dRadiusInner+1.;
            if (dRadiusOuter > dRadius) dRadiusOuter = dRadius;
        } while (lFound <= 0 && dRadiusOuterSave < dRadius);

        if (lFound < 1) return (0L);
        while ( CVectorSize(distances) < k && dRadiusOuter < dRadius) {
            if (numrad < CNEARTREE_DIMSAMPLES) {
                foundatrad[numrad] = (double)CVectorSize(distances);
                radlist[numrad++] = dRadiusOuter;
                if (numrad > 1) {
                    foundatrad[numrad-1]+= foundatrad[numrad-2];
                    dimlist[numrad-2]
                    = log(foundatrad[numrad-1]-foundatrad[numrad-2])
                    /(log(radlist[numrad-1]-radlist[numrad-2])+1.e-38);
                }
                shell = 1;
                closed= 0;
                dRadiusInner = dRadiusOuter;
                dRadiusOuter = dRadiusInner
                +treehandle->m_SumSpacings
                /sqrt((double)(1+CVectorSize(treehandle->m_ObjectStore)));
                if (dRadiusOuter <= dRadiusInner) dRadiusOuter = dRadiusInner+1.;
                if (dRadiusOuter > dRadius) dRadiusOuter = dRadius;
                if (numrad == CNEARTREE_DIMSAMPLES) {
                    int ii;
                    dimest = 0.;
                    for (ii=0; ii < numrad-1; ii++) {
                        dimest += dimlist[ii];
                    }
                    dimest = dimest/((double) (numrad-1));
                }
            } else {
                shell = 0;
                closed = 0;
                dRadiusInner = dRadiusOuter;
                dRadiusOuter = dRadiusInner*pow(((double)k)/((double)CVectorSize(distances)),1./(3.*dimest));
                if (dRadiusOuter <= dRadiusInner) dRadiusOuter = dRadiusInner+1.;
                if (dRadiusOuter > dRadius) dRadiusOuter = dRadius;
            }
            /*fprintf(stderr,"dRadiusInner, dRadiusOuter %g %g\n",
             (double)(dRadiusInner), (double)(dRadiusOuter)); */
            /* Add any point from the next annular shell */
            do {
                dRadiusOuterSave = dRadiusOuter;
                size = CVectorSize(distances);
                /*fprintf(stderr,"k %d, dRadius %g, dRadiusInner, dRadiusOuter %g %g\n", (int)k,
                        (double)dRadius, (double)(dRadiusInner), (double)(dRadiusOuter));*/
                error = CNearTreeFindKNearInAnnulus(treehandle,k-size,shell,closed,
                                                    dRadiusInner,
                                                    dRadiusOuter,
                                                    distances,
                                                    indices,
                                                    coord);
                lFound = CVectorSize(distances) - size;
                if (lFound > 0) {
                    CVectorGetElement(distances,&dRadiusOuter,CVectorSize(distances)-1);
                    /*fprintf (stderr,
                             "Ran CNearTreeFindKNearInAnnulus, lFound %ld , dRadiusInner %g, dRadiusOuter %g\n",
                             (long)lFound, dRadiusInner, dRadiusOuter);*/
                    break;
                }
                dRadiusOuter = dRadiusOuterSave+(dRadiusOuterSave-dRadiusInner)*1.1;
                dRadiusInner = dRadiusOuterSave;
                if (dRadiusOuter <= dRadiusInner) dRadiusOuter = dRadiusInner+1.;
                if (dRadiusOuter > dRadius) dRadiusOuter = dRadius;
                
            } while (lFound <= 0 && dRadiusOuterSave < dRadius);
            /* fprintf(stderr,"found %ld points\n",lFound); */
            if (dRadiusOuter >= dRadius-1.e-36) break;
        }
        
        return error;

    }

    int CNearTreeFindKNearest_Sphere (const CNearTreeHandle treehandle,
                               const size_t k,
                               const double dRadius,
                               CVectorHandle distances,
                               CVectorHandle indices,
                               const void * coord){
        double dDR, dDL, dTarget, dist, drat;
        int nopoints, l2lazy;
        CVectorHandle sStack;
        CVectorHandle coords;
        CVectorHandle objs;
        void CNEARTREE_FAR * xcoord;
        void CNEARTREE_FAR * xobj;
        CNearTreeNodeHandle pt;
        size_t ii, index;
        enum  { left, right, end } eDir;
        
        eDir = left; /* examine the left nodes first */
        
        nopoints = 1;
        
        dDR = dDL = DBL_MAX;
        
        dTarget = dRadius;
        
        drat = 1.;
        
        l2lazy = 0;
        
        if ((treehandle->m_iflags&CNEARTREE_NORM_L2LAZY)) {
            
            drat = sqrt((double)(treehandle->m_szdimension));
            l2lazy = 1;
        }
        
        if (dRadius < 0.) return CNEARTREE_BAD_ARGUMENT;
        
        if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        
        if ( treehandle->m_DelayedIndices ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }
        
        pt = treehandle->m_ptTree;
        
        coords=treehandle->m_CoordStore;
        
        objs=treehandle->m_ObjectStore;
        
        if (!pt) return CNEARTREE_BAD_ARGUMENT;
        
        if (!(pt->m_iflags&CNEARTREE_FLAG_LEFT_DATA)) return CNEARTREE_NOT_FOUND;
        
#ifdef CNEARTREE_INSTRUMENTED
        (treehandle->m_NodeVisits)++;
#endif                                               
        
        if (CVectorCreate(&sStack,sizeof(CNearTreeNodeHandle),10))
            return CNEARTREE_MALLOC_FAILED;
        
        
        while (!(eDir == end && CVectorSize(sStack) == 0)) {
            
            if ( eDir == right ) {
                dDR = DBL_MAX;
                if ((pt->m_iflags)&CNEARTREE_FLAG_RIGHT_DATA) {
                    dist = dDR = CNearTreeDist(treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexRight));
                    if (dDR <= dTarget*drat ) {
                        if (l2lazy) {
                            CNTM_DistL2(dist,treehandle, 
                                    (void CNEARTREE_FAR *)coord,
                                    CVectorElementAt(coords,pt->m_indexRight));
                        }
                        if (!l2lazy || dist <= dTarget) {
                            nopoints = 0;
                            CNearTreeSortIn2(distances,indices,dist,pt->m_indexRight,k,0,0,0.);
                            if (CVectorSize(distances)==k) {
                                CVectorGetElement(distances,&dTarget,k-1);
                            }
                        }
                    }
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_RIGHT_CHILD)&& 
                    (TRIANG(dDR,pt->m_dMaxRight,dTarget))){
                    /* we did the left and now we finished the right, go down */
                    pt = pt->m_pRightBranch;
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif                                                           
                    eDir = left;
                } else {
                    eDir = end;
                }
            }
            if ( eDir == left ) {
                if ((pt->m_iflags)&CNEARTREE_FLAG_LEFT_DATA) {
                    dist = dDL = CNearTreeDist(treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexLeft));
                    if (dDL <= dTarget*drat ) {
                        if (l2lazy) {
                            CNTM_DistL2(dist,treehandle, 
                                        (void CNEARTREE_FAR *)coord,
                                        CVectorElementAt(coords,pt->m_indexLeft));
                        } 
                        if (!l2lazy || dist <= dTarget) {
                            nopoints = 0;
                            CNearTreeSortIn2(distances,indices,dist,pt->m_indexLeft,k,0,0,0.);
                            if (CVectorSize(distances)==k) {
                                CVectorGetElement(distances,&dTarget,k-1);
                        }
                        }
                    }
                }
                if (pt->m_iflags&CNEARTREE_FLAG_RIGHT_DATA) {
                    CVectorAddElement(sStack,&pt);
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_LEFT_CHILD)&&
                    (TRIANG(dDL,pt->m_dMaxLeft,dTarget))){
                    pt = pt->m_pLeftBranch;
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif                                                           
                } else {
                    eDir = end;
                }
            }
            if ( eDir == end && CVectorSize(sStack) != 0 ) {
                CVectorGetElement(sStack,&pt,CVectorSize(sStack)-1);
                CVectorRemoveElement(sStack,CVectorSize(sStack)-1);
                eDir = right;
            }
            
        }    
        
        CVectorFree(&sStack);
        
        return  nopoints?CNEARTREE_NOT_FOUND:CNEARTREE_SUCCESS;
    }     

    /*
     =======================================================================
     int CNearTreeFindKNearInAnnulus ( const CNearTreeHandle treehandle,
     const size_t k,
     const int shell,
     const int closed,
     const double dRadiusInner,
     double dRadiusOuter,
     CVectorHandle distInRing,
     CVectorHandle indexInRing,
     const void * coord);
     
     Function to search a Neartree for the set of objects closer to some probe point, coord,
     than dRadius.
     
     k is the maximum number of closest neighbors to return.  Finds this many if passible.
     
     shell if non-zero, the search only returns hits in the nearest thin shell
     
     closed if non-zero, point at dRadiusInner will be included
     
     dRadiusInner is the minimum radius - any point closer than dRadiusInner from the probe point will be ignored, if closed is zero, any point at dRadiusInner will
          also be ignored
     
     dRadiusOuter is the maximum search radius - any point further than dRadius from the probe point will be ignored,
     
     distInRing is an existing vector of the distances to the k closest points
          to the probe point in the annulus that can be found in the Neartree.  The caller
          should have created this vector with CVectorCreate(&distanceInRing,sizeof(double),k)
     
     indexInRing is an existing vector of indices of the coordinates and objects of 
          k closest points from the probe point in the annulus
          that can be found in the Neartree
     
     coord  is the probe point
     
     return value is 0 if points were found
     
     =======================================================================
     */
    
    int CNearTreeFindKNearInAnnulus (const CNearTreeHandle treehandle,
                                     const size_t k,
                                     const int shell,
                                     const int closed,
                                     const double dRadiusInner,
                                     double dRadiusOuter,
                                     CVectorHandle distInRing,
                                     CVectorHandle indexInRing,
                                     const void * coord){
        double dDR, dDL, dTarget, dist, drat;
        double dRadiusCur;
        int nopoints, l2lazy;
        CVectorHandle sStack;
        CVectorHandle coords;
        CVectorHandle objs;
        void CNEARTREE_FAR * xdist;
        void CNEARTREE_FAR * xcoord;
        void CNEARTREE_FAR * xobj;
        CNearTreeNodeHandle pt;
        size_t ii, index, size;
        enum  { left, right, end } eDir;
        
        eDir = left; /* examine the left nodes first */
        
        nopoints = 1;
        
        dDR = dDL = DBL_MAX;
        
        dTarget = dRadiusOuter;
        
        drat = 1.;
        
        l2lazy = 0;
        
        if ((treehandle->m_iflags&CNEARTREE_NORM_L2LAZY)) {
            
            /* compute the worst case ratio of the L1 norm to the L2 norm */
            drat = sqrt((double)(treehandle->m_szdimension));
            l2lazy = 1;
            
        }
        
        if (dRadiusInner < 0.
            || dRadiusOuter < 0.
            || dRadiusInner > dRadiusOuter) return CNEARTREE_BAD_ARGUMENT;
        
        if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        
        if ( treehandle->m_DelayedIndices ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }
        
        dRadiusCur = dRadiusInner;
        size = CVectorSize(distInRing);
        if (size > 0) {
            CVectorGetElement(distInRing,&dRadiusCur,size-1);
        }
        
        pt = treehandle->m_ptTree;
        
        coords=treehandle->m_CoordStore;
        
        objs=treehandle->m_ObjectStore;
        
        if (!pt) return CNEARTREE_BAD_ARGUMENT;
        
        if (!(pt->m_iflags&CNEARTREE_FLAG_LEFT_DATA)) return CNEARTREE_NOT_FOUND;
        
#ifdef CNEARTREE_INSTRUMENTED
        (treehandle->m_NodeVisits)++;
#endif
        
        if (CVectorCreate(&sStack,sizeof(CNearTreeNodeHandle),10))
            return CNEARTREE_MALLOC_FAILED;
        
        while (!(eDir == end && CVectorSize(sStack) == 0)) {
            
            if ( eDir == right ) {
                dDR = DBL_MAX;
                if ((pt->m_iflags)&CNEARTREE_FLAG_RIGHT_DATA) {
                    dist = dDR = CNearTreeDist(treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexRight));
                    if (dDR <= dRadiusOuter*drat
                        && (dDR > dRadiusInner
                           || ((closed||l2lazy) && dDR == dRadiusInner )  )) {
                        if (l2lazy) {
                            CNTM_DistL2(dist,treehandle,
                                        (void CNEARTREE_FAR *)coord,
                                        CVectorElementAt(coords,pt->m_indexRight));
                        }
                        if (!l2lazy
                            || (dist <= dRadiusOuter
                                && (dist > dRadiusInner
                                    || (closed && dist == dRadiusInner ) ))) {
                                nopoints = 0;
                                CNearTreeSortIn2(distInRing,indexInRing,dist,pt->m_indexRight,k,shell, size, dRadiusCur);
                                if (CVectorSize(distInRing)==k) {
                                    CVectorGetElement(distInRing,&dRadiusOuter,k-1);
                                }
                            }
                    }
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_RIGHT_CHILD)&&
                    (TRIANG(dDR,pt->m_dMaxRight,dRadiusOuter))&&
                    (TRIANG(dRadiusInner,dDR,pt->m_dMaxRight))){
                    /* we did the left and now we finished the right, go down */
                    pt = pt->m_pRightBranch;
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif
                    eDir = left;
                } else {
                    eDir = end;
                }
            }
            if ( eDir == left ) {
                if ((pt->m_iflags)&CNEARTREE_FLAG_LEFT_DATA) {
                    dist = dDL = CNearTreeDist(treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexLeft));
                    if (dDL <= dRadiusOuter*drat
                        && (dDL > dRadiusInner
                         || ((closed||l2lazy) && dDL == dRadiusInner )  )) {
                        if (l2lazy) {
                            CNTM_DistL2(dist,treehandle,
                                        (void CNEARTREE_FAR *)coord,
                                        CVectorElementAt(coords,pt->m_indexLeft));
                        }
                        if (!l2lazy
                            || (dist <= dRadiusOuter
                            && (dist > dRadiusInner || (closed && dist == dRadiusInner ) ))) {
                            nopoints = 0;
                            CNearTreeSortIn2(distInRing,indexInRing,dist,pt->m_indexLeft,k,shell,size,dRadiusCur);
                            if (CVectorSize(distInRing)==k) {
                                CVectorGetElement(distInRing,&dRadiusOuter,k-1);
                            }
                        }
                    }
                }
                if (pt->m_iflags&CNEARTREE_FLAG_RIGHT_DATA) {
                    CVectorAddElement(sStack,&pt);
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_LEFT_CHILD)&&
                    (TRIANG(dDL,pt->m_dMaxLeft,dRadiusOuter))&&
                    (TRIANG(dRadiusInner,dDL,pt->m_dMaxRight))){
                    pt = pt->m_pLeftBranch;
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif
                } else {
                    eDir = end;
                }
            }
            if ( eDir == end && CVectorSize(sStack) != 0 ) {
                CVectorGetElement(sStack,&pt,CVectorSize(sStack)-1);
                CVectorRemoveElement(sStack,CVectorSize(sStack)-1);
                eDir = right;
            }
            
        }
        
        CVectorFree(&sStack);
        
        return  nopoints?CNEARTREE_NOT_FOUND:CNEARTREE_SUCCESS;
    }

    
    
    /*
     =======================================================================
     int CNearTreeFindKTreeNearest ( const CNearTreeHandle treehandle,
     const size_t k,
     const double dRadius,
     CNearTreeHandle foundNearest,
     const void CNEARTREE_FAR * coord,
     int resetcount);
     
     Function to search a Neartree for the set of k objects closer to some probe point, coord,
     than dRadius.
     
     k is the maximum number of closest neighbors to return.  Finds this many if passible.
     dRadius is the minimum search radius - any point closer than dRadius to the probe
     point will be ignored
     foundClosest is an existing CNearTree to which the found points will be added
     coord  is the probe point
     resetcount should be non-zero to clear foundClosest on entry
     return value is 0 if points were found
     
     =======================================================================
     */
    int CNearTreeFindKTreeNearest (const CNearTreeHandle treehandle,
                                   const size_t k,
                                   const double dRadius,
                                   CNearTreeHandle foundClosest,
                                   const void CNEARTREE_FAR * coord,
                                   int resetcount){
        
        CVectorHandle objs;
        CVectorHandle coords;
        size_t ii;
        int errorcode;
        
        if ( !treehandle || !coord || !foundClosest ) return CNEARTREE_BAD_ARGUMENT;
        
        if ( treehandle->m_DelayedIndices ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }
        
        /* clear the contents of the return tree so that things don't accidentally accumulate */
        if (resetcount) {
            CNearTreeClear( foundClosest );
        }
        
        if (CVectorCreate(&objs,sizeof(void *),10)) {
            return CNEARTREE_MALLOC_FAILED;
        }
        if (CVectorCreate(&coords,sizeof(void *),10)) {
            return CNEARTREE_MALLOC_FAILED;
        }
        
        CNearTreeFindKNearest(treehandle, k, dRadius, coords, objs, coord, 0);
        
        for (ii=0; ii<CVectorSize(coords);ii++) {
            if ((errorcode=CNearTreeInsert(foundClosest,*((void CNEARTREE_FAR * CNEARTREE_FAR *)CVectorElementAt(coords,ii)),CVectorElementAt(objs,ii))),errorcode) {
                CVectorFree(&coords);
                CVectorFree(&objs);
                return errorcode;
            }
        }      
        errorcode = 0;
        if (CVectorSize(coords) == 0) errorcode=CNEARTREE_NOT_FOUND;
        
        CVectorFree(&coords);
        CVectorFree(&objs);
        
        return errorcode;
        
    }
    
    
    /*
     =======================================================================
     int CNearTreeFindKFarthest ( const CNearTreeHandle treehandle,
     const size_t k,
     const double dRadius,
     CVectorHandle coordFarthest,
     CVectorHandle objFarthest,
     const void * coord,
     int resetcount);
     
     Function to search a Neartree for the set of objects farther from some probe point, coord,
     than dRadius. 
     
     k is the maximum number of farthest neighbors to return.  Finds this many if passible.
     
     dRadius is the maximum search radius - any point farther than dRadius from the probe
     point will be ignored
     
     coordFarthest is a vector of pointers to coordinate tuples and is the 
     returned set of farthest points from the probe point that can be found in the Neartree
     
     objFarthest is a vector of objects and is the returned set of farthest points
     from the probe point that can be found in the Neartree
     
     coord  is the probe point
     
     resetcount should be non-zero to clear objFarthest on entry
     
     return value is 0 if points were found
     
     =======================================================================
     */
    
    int CNearTreeFindKFarthest ( const CNearTreeHandle treehandle,
                                const size_t k,
                                const double dRadius,
                                CVectorHandle coordFarthest,
                                CVectorHandle objFarthest,
                                const void * coord,
                                int resetcount){
        double dDR, dDL, dTarget, dist, drat;
        int nopoints, l2lazy;
        CVectorHandle sStack;
        CVectorHandle coords;
        CVectorHandle objs;
        CVectorHandle dDs;
        CVectorHandle stIndices;
        void CNEARTREE_FAR * xcoord;
        void CNEARTREE_FAR * xobj;
        CNearTreeNodeHandle pt;
        size_t ii, index;
        enum  { left, right, end } eDir;
        
        eDir = left; /* examine the left nodes first */
        
        nopoints = 1;
        
        dDR = dDL = DBL_MAX;
        
        dTarget = dRadius;
        
        drat = 1.;
        
        l2lazy = 0;
        
        if ((treehandle->m_iflags&CNEARTREE_NORM_L2LAZY)) {
            
            drat = sqrt((double)(treehandle->m_szdimension));
            l2lazy = 1;
        }
        
        
        if (dRadius < 0.) return 1;
        
        if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        
        if ( treehandle->m_DelayedIndices ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }
        
        pt = treehandle->m_ptTree;
        
        coords=treehandle->m_CoordStore;
        
        objs=treehandle->m_ObjectStore;
        
        if (!pt) return CNEARTREE_BAD_ARGUMENT;
        
        if (!(pt->m_iflags&CNEARTREE_FLAG_LEFT_DATA)) return CNEARTREE_NOT_FOUND;
        
#ifdef CNEARTREE_INSTRUMENTED
        (treehandle->m_NodeVisits)++;
#endif                                               
        
        if (CVectorCreate(&sStack,sizeof(CNearTreeNodeHandle),10) ||
            
            CVectorCreate(&dDs,sizeof(double),k) ||
            
            CVectorCreate(&stIndices,sizeof(size_t),k)) return CNEARTREE_MALLOC_FAILED;
        
        /* clear the contents of the return vector so that things don't accidentally accumulate */
        if (resetcount) {
            if (coordFarthest) CVectorClear( coordFarthest );
            if (objFarthest) CVectorClear( objFarthest );
        }
        
        while (!(eDir == end && CVectorSize(sStack) == 0)) {
            
            if ( eDir == right ) {
                dDR = DBL_MAX;
                if ((pt->m_iflags)&CNEARTREE_FLAG_RIGHT_DATA) {
                    dist = dDR = CNearTreeDist(treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexRight));
                    if (dDR >= dTarget/drat ) {
                        if (l2lazy) {
                            CNTM_DistL2(dist,treehandle, 
                                        (void CNEARTREE_FAR *)coord,
                                        CVectorElementAt(coords,pt->m_indexRight));
                        }
                        if (!l2lazy || dist >= dTarget ) {
                            nopoints = 0;
                            CNearTreeSortIn2(dDs,stIndices,-dist,pt->m_indexRight,k,0,0,0.);
                            if (CVectorSize(dDs)==k) {
                                CVectorGetElement(dDs,&dTarget,k-1);
                                dTarget = -dTarget;
                            }
                        }
                    }
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_RIGHT_CHILD)&& 
                    (TRIANG(dTarget,dDR,pt->m_dMaxRight))){
                    /* we did the left and now we finished the right, go down */
                    pt = pt->m_pRightBranch;
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif                                                           
                    eDir = left;
                } else {
                    eDir = end;
                }
            }
            if ( eDir == left ) {
                if ((pt->m_iflags)&CNEARTREE_FLAG_LEFT_DATA) {
                    dist = dDL = CNearTreeDist(treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexLeft));
                    if (dDL >= dTarget/drat ) {
                        if (l2lazy) {
                            CNTM_DistL2(dist,treehandle, 
                                        (void CNEARTREE_FAR *)coord,
                                        CVectorElementAt(coords,pt->m_indexLeft));
                        }   
                        if (!l2lazy || dist >= dTarget ) {
                            nopoints = 0;
                            CNearTreeSortIn2(dDs,stIndices,-dist,pt->m_indexLeft,k,0,0,0.);
                            if (CVectorSize(dDs)==k) {
                                CVectorGetElement(dDs,&dTarget,k-1);
                                dTarget = -dTarget;
                            }
                        }
                    }
                }
                if (pt->m_iflags&CNEARTREE_FLAG_RIGHT_DATA) {
                    CVectorAddElement(sStack,&pt);
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_LEFT_CHILD)&&
                    (TRIANG(dTarget,dDL,pt->m_dMaxLeft))){
                    pt = pt->m_pLeftBranch;
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif                                                           
                } else {
                    eDir = end;
                }
            }
            if ( eDir == end && CVectorSize(sStack) != 0 ) {
                CVectorGetElement(sStack,&pt,CVectorSize(sStack)-1);
                CVectorRemoveElement(sStack,CVectorSize(sStack)-1);
                eDir = right;
            }
            
        }    
        
        CVectorFree(&sStack);
        
        if (!nopoints) {
            for (ii = 0; ii < CVectorSize(stIndices); ii++) {
                CVectorGetElement(stIndices,&index,ii);
                if (coordFarthest)  {
                    xcoord = CVectorElementAt(coords,index);
                    CVectorAddElement(coordFarthest,&xcoord);
                }
                if (objFarthest) {
                    xobj = CVectorElementAt(objs,index);
                    CVectorAddElement(objFarthest,&xobj);
                }
                
            }
        }
        
        CVectorFree(&dDs);
        CVectorFree(&stIndices);
        
        return  nopoints?CNEARTREE_NOT_FOUND:CNEARTREE_SUCCESS;
    }     
    
    /*
     =======================================================================
     int CNearTreeFindKTreeFarthest ( const CNearTreeHandle treehandle,
     const size_t k,
     const double dRadius,
     CNearTreeHandle foundFarthest,
     const void CNEARTREE_FAR * coord,
     int resetcount);
     
     Function to search a Neartree for the set of k objects farther from some probe point, coord,
     than dRadius.
     
     k is the maximum number of farthest neighbors to return.  Finds this many if passible.
     dRadius is the minimum search radius - any point farther than dRadius from the probe
     point will be ignored
     foundFarthest is an existing CNearTree to which the found points will be added
     coord  is the probe point
     resetcount should be non-zero to clear foundFarthest on entry
     return value is 0 if points were found
     
     =======================================================================
     */
    int CNearTreeFindKTreeFarthest ( const CNearTreeHandle treehandle,
                                    const size_t k,
                                    const double dRadius,
                                    CNearTreeHandle foundFarthest,
                                    const void CNEARTREE_FAR * coord,
                                    int resetcount) {
        
        CVectorHandle objs;
        CVectorHandle coords;
        size_t ii;
        int errorcode;
        
        if ( !treehandle || !coord || !foundFarthest ) return CNEARTREE_BAD_ARGUMENT;
        
        if ( treehandle->m_DelayedIndices ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }
        
        /* clear the contents of the return tree so that things don't accidentally accumulate */
        if (resetcount) {
            CNearTreeClear( foundFarthest );
        }
        
        if (CVectorCreate(&objs,sizeof(void *),10)) {
            return CNEARTREE_MALLOC_FAILED;
        }
        if (CVectorCreate(&coords,sizeof(void *),10)) {
            return CNEARTREE_MALLOC_FAILED;
        }
        
        CNearTreeFindKFarthest(treehandle, k, dRadius, coords, objs, coord, 0);
        
        for (ii=0; ii<CVectorSize(coords);ii++) {
            if ((errorcode=CNearTreeInsert(foundFarthest,*((void CNEARTREE_FAR * CNEARTREE_FAR *)CVectorElementAt(coords,ii)),CVectorElementAt(objs,ii))),errorcode) {
                CVectorFree(&coords);
                CVectorFree(&objs);
                return errorcode;
            }
        }      
        errorcode = 0;
        if (CVectorSize(coords) == 0) errorcode=CNEARTREE_NOT_FOUND;
        
        CVectorFree(&coords);
        CVectorFree(&objs);
        
        return errorcode;
        
    }
    
    
    
    
    /*
     =======================================================================
     int CNearTreeNearest ( const CNearTreeHandle treehandle, 
     double CNEARTREE_FAR *dRadius,  
     void CNEARTREE_FAR * CNEARTREE_FAR * coordClosest,
     void CNEARTREE_FAR * CNEARTREE_FAR * objClosest,
     const void CNEARTREE_FAR * coord )
     
     Function to search a Neartree for the object closest to some probe point, coord.
     This function is called by CNearTreeNearestNeighbor.
     
     *dRadius is the smallest currently known distance of an object from the probe point
     and is modified as the search progresses.
     coordClosest is a pointer to the coordinate vector of the nearest point
     objClosest is a pointer to a pointer to hold the corresponding object or is NULL
     coord  is the probe point
     the return value is 0 only if a point was found within dRadius
     
     This version search down whichever branch seems shortest first
     
     =======================================================================
     =======================================================================
     int CNearTreeLeftNearest ( const CNearTreeHandle treehandle, 
     double CNEARTREE_FAR *dRadius,  
     void CNEARTREE_FAR * CNEARTREE_FAR * coordClosest,
     void CNEARTREE_FAR * CNEARTREE_FAR * objClosest,
     const void CNEARTREE_FAR * coord )
     
     Function to search a Neartree for the object closest to some probe point, coord.
     This function is called by CNearTreeNearestNeighbor.
     
     *dRadius is the smallest currently known distance of an object from the probe point
     and is modified as the search progresses.
     coordClosest is a pointer to the coordinate vector of the nearest point
     objClosest is a pointer to a pointer to hold the corresponding object or is NULL
     coord  is the probe point
     the return value is 0 only if a point was found within dRadius
     
     This version searches down the left branch first
     
     =======================================================================
     */
    
    
    int CNearTreeNearestL2LAZY ( const CNearTreeHandle treehandle, 
                                double CNEARTREE_FAR * dRadius,  
                                void CNEARTREE_FAR * CNEARTREE_FAR * coordClosest,
                                void CNEARTREE_FAR * CNEARTREE_FAR * objClosest,
                                const void CNEARTREE_FAR * coord ) {
        double   dDR, dDL, dDRC, dDLC, dDRCsq, dDLCsq;
        CVectorHandle sStack;
        CVectorHandle coords;
        CVectorHandle objs;
        CNearTreeNodeHandle pt;
        void CNEARTREE_FAR * pobjClosest;
        void CNEARTREE_FAR * pcoordClosest;
        double drat;
        
        pobjClosest = NULL;
        pcoordClosest = NULL;
        
        dDR = dDL = dDRC = dDLC = -1.;
        
        drat = sqrt((double)(treehandle->m_szdimension));
        
        if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        
        if ( treehandle->m_DelayedIndices ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }
        
        pt = treehandle->m_ptTree;
        
        if (!pt) return CNEARTREE_BAD_ARGUMENT;
        
        if (!(pt->m_iflags&CNEARTREE_FLAG_LEFT_DATA)) return CNEARTREE_NOT_FOUND;
        
        coords=treehandle->m_CoordStore;
        
        objs=treehandle->m_ObjectStore;        
        
        CVectorCreate(&sStack,sizeof(CNearTreeNodeHandle),10);
        
#ifdef CNEARTREE_INSTRUMENTED
        (treehandle->m_NodeVisits)++;
#endif          
        if (!(pt->m_iflags&(CNEARTREE_FLAG_LEFT_DATA|CNEARTREE_FLAG_RIGHT_DATA))) {
            
            CVectorFree(&sStack);
            
            return CNEARTREE_NOT_FOUND;
            
        }
        
        while ((pt->m_iflags&(CNEARTREE_FLAG_LEFT_DATA|CNEARTREE_FLAG_RIGHT_DATA))
               || CVectorSize(sStack) != 0) {
            
            if (!(pt->m_iflags&(CNEARTREE_FLAG_LEFT_DATA|CNEARTREE_FLAG_RIGHT_DATA))) {                    
                if (CVectorSize(sStack) !=0 ) {
                    CVectorGetElement(sStack,&pt,CVectorSize(sStack)-1);
                    CVectorRemoveElement(sStack,CVectorSize(sStack)-1);
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif    
                    continue;
                }
                break;
            }
            
            if ((pt->m_iflags&(CNEARTREE_FLAG_LEFT_DATA))) {
                CNTM_DistL1(dDL,treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexLeft));
                dDLC = dDL;
                if (dDL <= *dRadius) {
                    *dRadius = dDL;
                    pobjClosest = CVectorElementAt(objs,pt->m_indexLeft);
                    pcoordClosest = CVectorElementAt(coords,pt->m_indexLeft);                    
                } else if (dDL <= (*dRadius)*drat) {
                    CNTM_DistL2sq(dDLCsq,treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexLeft));
                    if (dDLCsq <= (*dRadius)*(*dRadius)) {
                        *dRadius = dDLC = sqrt(dDLCsq);
                        pobjClosest = CVectorElementAt(objs,pt->m_indexLeft);
                        pcoordClosest = CVectorElementAt(coords,pt->m_indexLeft);                    
                    }
                }
            }
            if ((pt->m_iflags&(CNEARTREE_FLAG_RIGHT_DATA))) {
#ifdef CNEARTREE_INSTRUMENTED
                (treehandle->m_NodeVisits)++;
#endif    
                CNTM_DistL1(dDR,treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexLeft));
                dDRC = dDR;
                if (dDR <= *dRadius) {
                    *dRadius = dDR;
                    pobjClosest = CVectorElementAt(objs,pt->m_indexRight);
                    pcoordClosest = CVectorElementAt(coords,pt->m_indexRight);                    
                } else if (dDR <= (*dRadius)*drat) {
                    CNTM_DistL2sq(dDRCsq,treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexRight));
                    if (dDRCsq <= (*dRadius)*(*dRadius)) {
                        *dRadius = dDRC = sqrt(dDRCsq);
                        pobjClosest = CVectorElementAt(objs,pt->m_indexRight);
                        pcoordClosest = CVectorElementAt(coords,pt->m_indexRight);                    
                    }
                }
            }
            

            /*
             See if both branches are populated.  In that case, save one branch
             on the stack, and process the other one based on which one seems
             smaller, but useful first]
             */
            
            if ((pt->m_iflags&CNEARTREE_FLAG_RIGHT_CHILD) && (pt->m_iflags&CNEARTREE_FLAG_LEFT_CHILD)){
                if (dDLC+pt->m_dMaxLeft < dDRC+pt->m_dMaxRight ) {
                    if ( TRIANG(dDLC,pt->m_dMaxLeft,*dRadius)) {
                        if ( TRIANG(dDRC,pt->m_dMaxRight,*dRadius)) {
                            CVectorAddElement(sStack,&(pt->m_pRightBranch));
                        }
                        pt = pt->m_pLeftBranch;
#ifdef CNEARTREE_INSTRUMENTED
                        (treehandle->m_NodeVisits)++;
#endif                                
                        continue;
                    }
                    /* If we are here, the left branch was not useful
                     Fall through to use the right
                     */
                }
                
                /* We come here either because pursuing the left branch was not useful
                 or the right branch looks shorter
                 */
                if ( TRIANG(dDRC,pt->m_dMaxRight,*dRadius)) {
                    if ( TRIANG(dDLC,pt->m_dMaxLeft,*dRadius)) {
                        CVectorAddElement(sStack,&(pt->m_pLeftBranch));
                    }
                    pt = pt->m_pRightBranch;
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif                
                    continue;
                } 
                
            }
            
            /* Only one branch is viable, try them one at a time
             */
            if ( (pt->m_iflags&CNEARTREE_FLAG_LEFT_CHILD) && TRIANG(dDLC,pt->m_dMaxLeft,*dRadius)) {
                pt = pt->m_pLeftBranch;
#ifdef CNEARTREE_INSTRUMENTED
                (treehandle->m_NodeVisits)++;
#endif                            
                continue;
            }
            
            if (  (pt->m_iflags&CNEARTREE_FLAG_RIGHT_CHILD) && TRIANG(dDRC,pt->m_dMaxRight,*dRadius)) {
                pt = pt->m_pRightBranch;
#ifdef CNEARTREE_INSTRUMENTED
                (treehandle->m_NodeVisits)++;
#endif                            
                continue;
            }
            if ( CVectorSize(sStack) != 0 ) {
                CVectorGetElement(sStack,&pt,CVectorSize(sStack)-1);
                CVectorRemoveElement(sStack,CVectorSize(sStack)-1);
#ifdef CNEARTREE_INSTRUMENTED
                (treehandle->m_NodeVisits)++;
#endif                            
                continue;
            }
            break;
        }
        
        CVectorFree(&sStack);
        if (coordClosest) *coordClosest = pcoordClosest;
        if (objClosest) *objClosest = pobjClosest;
        return  pcoordClosest?CNEARTREE_SUCCESS:CNEARTREE_NOT_FOUND;
    }   /* NearestL2LAZY*/

    
    
    int CNearTreeNearest ( const CNearTreeHandle treehandle, 
                          double CNEARTREE_FAR * dRadius,  
                          void CNEARTREE_FAR * CNEARTREE_FAR * coordClosest,
                          void CNEARTREE_FAR * CNEARTREE_FAR * objClosest,
                          const void CNEARTREE_FAR * coord ) {
        double   dDR, dDL;
        CVectorHandle sStack;
        CVectorHandle coords;
        CVectorHandle objs;
        CNearTreeNodeHandle pt;
        void CNEARTREE_FAR * pobjClosest;
        void CNEARTREE_FAR * pcoordClosest;
        
        pobjClosest = NULL;
        pcoordClosest = NULL;
        
        if ((treehandle->m_iflags&CNEARTREE_NORM_L2LAZY)) 
            return CNearTreeNearestL2LAZY(treehandle,dRadius,coordClosest,objClosest,coord);
        
        dDR = dDL = -1.;
        
        if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        
        if ( treehandle->m_DelayedIndices ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }
        
        pt = treehandle->m_ptTree;
        
        if (!pt) return CNEARTREE_BAD_ARGUMENT;
        
        if (!(pt->m_iflags&CNEARTREE_FLAG_LEFT_DATA)) return CNEARTREE_NOT_FOUND;
        
        coords=treehandle->m_CoordStore;
        
        objs=treehandle->m_ObjectStore;        
        
        CVectorCreate(&sStack,sizeof(CNearTreeNodeHandle),10);
        
#ifdef CNEARTREE_INSTRUMENTED
        (treehandle->m_NodeVisits)++;
#endif          
        if (!(pt->m_iflags&(CNEARTREE_FLAG_LEFT_DATA|CNEARTREE_FLAG_RIGHT_DATA))) {
         
           CVectorFree(&sStack);
         
           return CNEARTREE_NOT_FOUND;

        }
        
        while ((pt->m_iflags&(CNEARTREE_FLAG_LEFT_DATA|CNEARTREE_FLAG_RIGHT_DATA))
            || CVectorSize(sStack) != 0) {
                
                if (!(pt->m_iflags&(CNEARTREE_FLAG_LEFT_DATA|CNEARTREE_FLAG_RIGHT_DATA))) {                    
                    if (CVectorSize(sStack) !=0 ) {
                        CVectorGetElement(sStack,&pt,CVectorSize(sStack)-1);
                        CVectorRemoveElement(sStack,CVectorSize(sStack)-1);
#ifdef CNEARTREE_INSTRUMENTED
                        (treehandle->m_NodeVisits)++;
#endif    
                        continue;
                    }
                    break;
                }
                
                if ((pt->m_iflags&(CNEARTREE_FLAG_LEFT_DATA))) {
                    dDL = CNearTreeDist(treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexLeft));
                    if (dDL <= *dRadius ) {
                        *dRadius = dDL;
                        pobjClosest = CVectorElementAt(objs,pt->m_indexLeft);
                        pcoordClosest = CVectorElementAt(coords,pt->m_indexLeft);
                    }                    
                }
                if ((pt->m_iflags&(CNEARTREE_FLAG_RIGHT_DATA))) {
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif                      
                    dDR = CNearTreeDist(treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexRight));
                    if (dDR <= *dRadius ) {
                        *dRadius = dDR;
                        pobjClosest = CVectorElementAt(objs,pt->m_indexRight);
                        pcoordClosest = CVectorElementAt(coords,pt->m_indexRight);
                    }                    
                }
                
                /*
                 See if both branches are populated.  In that case, save one branch
                 on the stack, and process the other one based on which one seems
                 smaller, but useful first]
                 */
                
                if ((pt->m_iflags&CNEARTREE_FLAG_RIGHT_CHILD) && (pt->m_iflags&CNEARTREE_FLAG_LEFT_CHILD)){
                    if (dDL+pt->m_dMaxLeft < dDR+pt->m_dMaxRight ) {
                        if ( TRIANG(dDL,pt->m_dMaxLeft,*dRadius)) {
                            if ( TRIANG(dDR,pt->m_dMaxRight,*dRadius)) {
                                CVectorAddElement(sStack,&(pt->m_pRightBranch));
                            }
                            pt = pt->m_pLeftBranch;
#ifdef CNEARTREE_INSTRUMENTED
                            (treehandle->m_NodeVisits)++;
#endif                                
                            continue;
                        }
                        /* If we are here, the left branch was not useful
                         Fall through to use the right
                         */
                    }
                    
                    /* We come here either because pursuing the left branch was not useful
                     of the right branch look shorter
                     */
                    if ( TRIANG(dDR,pt->m_dMaxRight,*dRadius)) {
                        if ( TRIANG(dDL,pt->m_dMaxLeft,*dRadius)) {
                            CVectorAddElement(sStack,&(pt->m_pLeftBranch));
                        }
                        pt = pt->m_pRightBranch;
#ifdef CNEARTREE_INSTRUMENTED
                        (treehandle->m_NodeVisits)++;
#endif                
                        continue;
                    } 
                    
                }
                
                /* Only one branch is viable, try them one at a time
                 */
                if ( (pt->m_iflags&CNEARTREE_FLAG_LEFT_CHILD) && TRIANG(dDL,pt->m_dMaxLeft,*dRadius)) {
                    pt = pt->m_pLeftBranch;
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif                            
                    continue;
                }
                
                if (  (pt->m_iflags&CNEARTREE_FLAG_RIGHT_CHILD) && TRIANG(dDR,pt->m_dMaxRight,*dRadius)) {
                    pt = pt->m_pRightBranch;
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif                            
                    continue;
                }
                if ( CVectorSize(sStack) != 0 ) {
                    CVectorGetElement(sStack,&pt,CVectorSize(sStack)-1);
                    CVectorRemoveElement(sStack,CVectorSize(sStack)-1);
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif                            
                    continue;
                }
                    break;
        }
        
        CVectorFree(&sStack);
        if (coordClosest) *coordClosest = pcoordClosest;
        if (objClosest) *objClosest = pobjClosest;
        return  pcoordClosest?CNEARTREE_SUCCESS:CNEARTREE_NOT_FOUND;
    }   /* Nearest */
    
    
         
    int CNearTreeLeftNearest ( const CNearTreeHandle treehandle, 
                               double CNEARTREE_FAR * dRadius,  
                               void CNEARTREE_FAR * CNEARTREE_FAR * coordClosest,
                               void CNEARTREE_FAR * CNEARTREE_FAR * objClosest,
                               const void CNEARTREE_FAR * coord ) {
         double   dDR, dDL;
         CVectorHandle sStack;
         CVectorHandle coords;
         CVectorHandle objs;
         CNearTreeNodeHandle pt;
         void CNEARTREE_FAR * pobjClosest;
         void CNEARTREE_FAR * pcoordClosest;
         enum  { left, right, end } eDir;
         
         eDir = left; /* examine the left nodes first */
         
         pobjClosest = NULL;
         pcoordClosest = NULL;
         
         if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
         
         if ( treehandle->m_DelayedIndices ) {
             if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
         }
         
         pt = treehandle->m_ptTree;
         
         if (!pt) return CNEARTREE_BAD_ARGUMENT;
         
         if (!(pt->m_iflags&CNEARTREE_FLAG_LEFT_DATA)) return CNEARTREE_NOT_FOUND;
         
         coords=treehandle->m_CoordStore;
         
         objs=treehandle->m_ObjectStore;        
         
         CVectorCreate(&sStack,sizeof(CNearTreeNodeHandle),10);
         
#ifdef CNEARTREE_INSTRUMENTED
         (treehandle->m_NodeVisits)++;
#endif                            
         
         while (!(eDir == end && CVectorSize(sStack) == 0)) {
             
             if ( eDir == right ) {
                 dDR = CNearTreeDist(treehandle, (void CNEARTREE_FAR *)coord,  CVectorElementAt(coords,pt->m_indexRight));
                 if (dDR <= *dRadius ) {
                     *dRadius = dDR;
                     pobjClosest = CVectorElementAt(objs,pt->m_indexRight);
                     pcoordClosest = CVectorElementAt(coords,pt->m_indexRight);
                 }
                 if ((pt->m_iflags&CNEARTREE_FLAG_RIGHT_CHILD)&& 
                     (TRIANG(dDR,pt->m_dMaxRight,*dRadius))) {
                     /* we did the left and now we finished the right, go down */
                     pt = pt->m_pRightBranch;
#ifdef CNEARTREE_INSTRUMENTED
                     (treehandle->m_NodeVisits)++;
#endif                                      
                     eDir = left;
                 } else {
                     eDir = end;
                 }
             }
             if ( eDir == left ) {
                 dDL = CNearTreeDist(treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexLeft));
                 if (dDL <= *dRadius ) {
                     *dRadius = dDL;
                     pobjClosest = CVectorElementAt(objs,pt->m_indexLeft);
                     pcoordClosest = CVectorElementAt(coords,pt->m_indexLeft);
                 }
                 if (pt->m_iflags&CNEARTREE_FLAG_RIGHT_DATA) {
                     CVectorAddElement(sStack,&pt);
                 }
                 if ((pt->m_iflags&CNEARTREE_FLAG_LEFT_CHILD)&&
                     (TRIANG(dDL,pt->m_dMaxLeft,*dRadius))) {
                     pt = pt->m_pLeftBranch;
#ifdef CNEARTREE_INSTRUMENTED
                     (treehandle->m_NodeVisits)++;
#endif                                       
                 } else {
                     eDir = end;
                 }
             }
             if ( eDir == end && CVectorSize(sStack) != 0 ) {
                 CVectorGetElement(sStack,&pt,CVectorSize(sStack)-1);
                 CVectorRemoveElement(sStack,CVectorSize(sStack)-1);
                 eDir = right;
#ifdef CNEARTREE_INSTRUMENTED
                 (treehandle->m_NodeVisits)++;
#endif    
             }
             
         }    
         
         CVectorFree(&sStack);
         if (coordClosest) *coordClosest = pcoordClosest;
         if (objClosest) *objClosest = pobjClosest;
         return  pcoordClosest?CNEARTREE_SUCCESS:CNEARTREE_NOT_FOUND;
     }   /* LeftNearest */
          
    
    /*
     =======================================================================
     int CNearTreeFindFarthest ( const CNearTreeHandle treehandle,
     double CNEARTREE_FAR * dRadius,  
     void CNEARTREE_FAR * CNEARTREE_FAR * coordFarthest,
     void CNEARTREE_FAR * CNEARTREE_FAR * objFarthest,
     const void CNEARTREE_FAR * coord )
     
     Function to search a Neartree for the object farthest from some probe point, coord.
     This function is called by CNearTreeFarthestNeighbor.
     
     *dRadius is the largest currently known distance of an object from the probe point
     and is modified as the search progresses.
     coordFarthest is a pointer to the returned farthest point
     from the probe point that can be found in the Neartree
     objFarthest is a pointer to a pointer to hold the corresponding object or is NULL
     coord  is the probe point
     the return value is 0 only if a point was found.
     
     =======================================================================
     */
    int CNearTreeFindFarthest ( const CNearTreeHandle treehandle,
                               double CNEARTREE_FAR * dRadius,  
                               void CNEARTREE_FAR * CNEARTREE_FAR * coordFarthest,
                               void CNEARTREE_FAR * CNEARTREE_FAR * objFarthest,
                               const void CNEARTREE_FAR * coord ) {
        double   dDR, dDL;
        CVectorHandle sStack;
        CVectorHandle coords;
        CVectorHandle objs;
        
        CNearTreeNodeHandle pt;
        void CNEARTREE_FAR * pobjFarthest;
        void CNEARTREE_FAR * pcoordFarthest;
        enum  { left, right, end } eDir;
        
        eDir = left; /* examine the left nodes first */
        
        pobjFarthest = NULL;
        pcoordFarthest = NULL;
        
        if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        
        if ( treehandle->m_DelayedIndices ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }
        
        pt = treehandle->m_ptTree;
        
        if (!pt) return CNEARTREE_BAD_ARGUMENT;
        
        if (!(pt->m_iflags&CNEARTREE_FLAG_LEFT_DATA)) return CNEARTREE_NOT_FOUND;
        
#ifdef CNEARTREE_INSTRUMENTED
         (treehandle->m_NodeVisits)++;
#endif                                                
        
        coords=treehandle->m_CoordStore;
        
        objs=treehandle->m_ObjectStore;        
        
        CVectorCreate(&sStack,sizeof(CNearTreeNodeHandle),10);
        
        while (!(eDir == end && CVectorSize(sStack) == 0)) {
            
            if ( eDir == right ) {
                dDR = CNearTreeDist(treehandle, (void CNEARTREE_FAR *)coord,  CVectorElementAt(coords,pt->m_indexRight));
                if (dDR >= *dRadius ) {
                    *dRadius = dDR;
                    pobjFarthest = CVectorElementAt(objs,pt->m_indexRight);
                    pcoordFarthest = CVectorElementAt(coords,pt->m_indexRight);
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_RIGHT_CHILD)&& 
                    (TRIANG(*dRadius,dDR,pt->m_dMaxRight))) {
                    /* we did the left and now we finished the right, go down */
                    pt = pt->m_pRightBranch;
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif                                       
                    eDir = left;
                } else {
                    eDir = end;
                }
            }
            if ( eDir == left ) {
                dDL = CNearTreeDist(treehandle, (void CNEARTREE_FAR *)coord, CVectorElementAt(coords,pt->m_indexLeft));
                if (dDL >= *dRadius ) {
                    *dRadius = dDL;
                    pobjFarthest = CVectorElementAt(objs,pt->m_indexLeft);
                    pcoordFarthest = CVectorElementAt(coords,pt->m_indexLeft);
                }
                if (pt->m_iflags&CNEARTREE_FLAG_RIGHT_DATA) {
                    CVectorAddElement(sStack,&pt);
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_LEFT_CHILD)&&
                    (TRIANG(*dRadius,dDL,pt->m_dMaxLeft))) {
                    pt = pt->m_pLeftBranch;
#ifdef CNEARTREE_INSTRUMENTED
                    (treehandle->m_NodeVisits)++;
#endif                                       
                } else {
                    eDir = end;
                }
            }
            if ( eDir == end && CVectorSize(sStack) != 0 ) {
                CVectorGetElement(sStack,&pt,CVectorSize(sStack)-1);
                CVectorRemoveElement(sStack,CVectorSize(sStack)-1);
                eDir = right;
            }
            
        }    
        
        CVectorFree(&sStack);
        if (coordFarthest) *coordFarthest = pcoordFarthest;
        if (objFarthest) *objFarthest = pobjFarthest;
        return  pcoordFarthest?CNEARTREE_SUCCESS:CNEARTREE_NOT_FOUND;
    }
    
    /*
     =======================================================================
     int CNearTreeObjects ( const CNearTreeHandle treehandle, CVectorHandle CNEARTREE_FAR * vectorhandle)
     
     Function to return the vector of objects in the tree.  This vector
     is not guaranteed to be in the same order as the order of insertion
     
     vectorhandle -- a pointer to a CVectorHandle
     
     =======================================================================
     */
    
    int CNearTreeObjects ( const CNearTreeHandle treehandle, CVectorHandle CNEARTREE_FAR * vectorhandle) {
        if ( !treehandle || !vectorhandle ) return CNEARTREE_BAD_ARGUMENT;
        
        *vectorhandle = treehandle->m_ObjectStore;
        
        return CNEARTREE_SUCCESS;
        
    }
    
    
    /*
     =======================================================================
     int CNearTreeCoords ( const CNearTreeHandle treehandle, CVectorHandle CNEARTREE_FAR * vectorhandle)
     
     Function to return the vector of coordinates in the tree.  This vector
     is not guaranteed to be in the same order as the order of insertion
     
     vectorhandle -- a pointer to a CVectorHandle
     
     =======================================================================
     */
    
    int CNearTreeCoords ( const CNearTreeHandle treehandle, CVectorHandle CNEARTREE_FAR * vectorhandle) {
        if ( !treehandle || !vectorhandle ) return CNEARTREE_BAD_ARGUMENT;
        
        *vectorhandle = treehandle->m_CoordStore;
        
        return CNEARTREE_SUCCESS;
        
    }
    
    
#ifdef __cplusplus
    
}

#endif

