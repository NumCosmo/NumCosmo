/*
 *  CNearTree.h
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

#ifndef CNEARTREE_H_INCLUDED
#define CNEARTREE_H_INCLUDED

#ifdef __cplusplus

extern "C" {
    
#endif
    
#ifdef CNEARTREE_USE_FAR
#include <malloc.h>
#define CNEARTREE_FAR __far
#define CNEARTREE_MALLOC _fmalloc
#define CNEARTREE_FREE _ffree
#define CNEARTREE_MEMSET _fmemset
#define CNEARTREE_MEMMOVE _fmemmove
#else
#include <stdlib.h>
#define CNEARTREE_FAR
#define CNEARTREE_MALLOC malloc
#define CNEARTREE_FREE free
#define CNEARTREE_MEMSET memset
#define CNEARTREE_MEMMOVE memmove
#endif
    
#include <limits.h>
#include <float.h>
#include <math.h>
#ifndef USE_LOCAL_HEADERS
#include <CVector.h>
#else
#include "CVector.h"
#endif
#ifndef CVECTOR_FAR
#define CVECTOR_FAR CNEARTREE_FAR
#endif
    
    
    
#ifndef USE_LOCAL_HEADERS
#include <rhrand.h>
#else
#include "rhrand.h"
#endif
    
    
    /* function returns */
#define CNEARTREE_SUCCESS       0
#define CNEARTREE_MALLOC_FAILED 1
#define CNEARTREE_BAD_ARGUMENT  2
#define CNEARTREE_NOT_FOUND     4
#define CNEARTREE_FREE_FAILED   8
#define CNEARTREE_CVECTOR_FAILED   16
    
    
    /* flags 
     0 for n data, n children
     */
    
#define CNEARTREE_FLAG_LEFT_DATA   1L          /* 0x0001 */
#define CNEARTREE_FLAG_RIGHT_DATA  2L          /* 0x0002 */
#define CNEARTREE_FLAG_LEFT_CHILD  4L          /* 0x0004 */
#define CNEARTREE_FLAG_RIGHT_CHILD 8L          /* 0x0008 */
    
#define CNEARTREE_DATA_OR_CHILDREN 15L         /* 0x000F */
    
#define CNEARTREE_TYPE_DOUBLE      16L         /* 0x0010 */
#define CNEARTREE_TYPE_INTEGER     32L         /* 0x0020 */
#define CNEARTREE_TYPE_STRING      64L         /* 0x0040 */
    
#define CNEARTREE_TYPE            112L         /* 0x0070 */
    
#define CNEARTREE_NORM_UNKNOWN    128L         /* 0x0080 */
#define CNEARTREE_NORM_L1         256L         /* 0x0100 */
#define CNEARTREE_NORM_L2         512L         /* 0x0200 */
#define CNEARTREE_NORM_LINF      1024L         /* 0x0400 */
#define CNEARTREE_NORM_SPHERE    2048L         /* 0x0800 */
#define CNEARTREE_NORM_HSPHERE   4096L         /* 0x1000 */
#define CNEARTREE_NORM_HAMMING   8192L         /* 0x2000 */
#define CNEARTREE_NORM_L2LAZY   16384L         /* 0x4000 */
    
#define CNEARTREE_NORM          32640L         /* 0x7F80 */
    
    
/*  Execution Control Flags */
    
#define CNTF_NOPREPRUNE    0x10000L     /*flag to supress all search prepruning */
#define CNTF_FORCEPREPRUNE 0x20000L     /*flag to force search prepruning       */
#define CNTF_NOFLIP        0x40000L     /*flag to suppress flips on insert      */
#define CNTF_FORCEFLIP     0x80000L     /*flag to force flips on insert         */
#define CNTF_NODEFER      0x100000L     /*flag to prevent deferred insert       */
#define CNTF_SKNN         0x200000L     /*flag to use spherical KNN             */

    
#define CNEARTREE_XFLAGS  0x3F0000L     /*mask for execution flags */
    
#ifdef CNEARTREE_FORCEPREPRUNE
  #define CNFT_FLAGDEFAULTPRUNE  CNTF_FORCEPREPRUNE
  #ifdef CNEARTREE_NOPREPRUNE
    #error "CNEARTREE_NOPREPRUNE conflicts with  CNEARTREE_FORCEPREPRUNE"
  #endif
#else
  #ifdef CNEARTREE_NOPREPRUNE
    #define  CNFT_FLAGDEFAULTPRUNE  CNTF_NOPREPRUNE
  #else
    #define  CNFT_FLAGDEFAULTPRUNE  0
  #endif
#endif
    
#ifdef CNEARTREE_FORCEFLIP
  #define CNFT_FLAGDEFAULTFLIP  CNTF_FORCEFLIP
  #ifdef CNEARTREE_NOFLIP
    #error "CNEARTREE_NOFLIP conflicts with  CNEARTREE_FORCEFLIP"
  #endif
#else
  #ifdef CNEARTREE_NOFLIP
    #define CNFT_FLAGDEFAULTFLIP CNTF_NOFLIP
  #else
    #define CNFT_FLAGDEFAULTFLIP 0
  #endif
#endif
    
#ifdef CNEARTREE_NODEFER
    #define CNFT_FLAGDEFAULTDEFER CNTF_NODEFER
#else
    #define CNFT_FLAGDEFAULTDEFER 0
#endif
    
#define CNTF_FLAGSDEFAULT  (CNFT_FLAGDEFAULTPRUNE|CNFT_FLAGDEFAULTFLIP|CNFT_FLAGDEFAULTDEFER)
    
    
#define CNEARTREE_FLIP          0     /* Not used, defined for compatibility */
#define CNEARTREE_DEFER_ALL     0     /* Not used, defined for compatibility */
    
    
    
    typedef struct _CNearTreeNode {
        size_t           m_indexLeft;     /* index of left coords in m_CoordStore  
                                             and of left object in m_ObjectStore   */
        size_t           m_indexRight;    /* index of right coords in m_CoordStore 
                                             and of right object in m_ObjectStore  */
        double           m_dMaxLeft;      /* longest distance from the left object
                                             to anything below it in the tree      */
        double           m_dMaxRight;     /* longest distance from the right object 
                                             to anything below it in the tree      */
        struct _CNearTreeNode CNEARTREE_FAR * m_pLeftBranch;  
                                          /* tree descending from the left object  */
        struct _CNearTreeNode CNEARTREE_FAR * m_pRightBranch; 
                                          /* tree descending from the right object */
        long              m_iflags;       /* flags                   */
        size_t            m_iTreeSize;    /* size of this node tree  */
#ifdef CNEARTREE_INSTRUMENTED
        size_t            m_Height;        /* height of this node     */
#endif

    } CNearTreeNode;
    
    
    typedef CNearTreeNode CNEARTREE_FAR * CNearTreeNodeHandle;   
    
    typedef struct {
        CNearTreeNodeHandle m_ptTree;     /* pointer to the actual tree                  */
        size_t           m_szdimension;   /* dimension of the coordinates                */
        size_t           m_szsize;        /* size of this tree                           */
        size_t           m_szdepth;       /* depth of this tree                          */
        long             m_iflags;        /* flags                                       */
        CVectorHandle    m_ObjectStore;   /* all inserted objects                        */
        CVectorHandle    m_CoordStore;    /* all inserted coordinates                    */
        CVectorHandle    m_DelayedIndices;/* objects queued for insertion                */

        CRHrand          m_rhr;           /* random number generator                     */
        double           m_DiamEstimate;  /* estimated diameter */
        double           m_SumSpacings;   /* sum of spacings at time of insertion */
        double           m_SumSpacingsSq; /* sum of spacings squared at time of insertion */
        double           m_DimEstimate;   /* estimated dimension */
        double           m_DimEstimateEsd;/* estimated dimension estimated standard deviation */
#ifdef CNEARTREE_INSTRUMENTED
        size_t           m_NodeVisits;    /* number of node visits */
#endif
    } CNearTree;
    
    typedef CNearTree     CNEARTREE_FAR * CNearTreeHandle;
    
/* Distance Macros for L2LAZY */
    
    
#define CNTM_Abs(x) (((x)<0.)?(-(x)):(x))
#define CNTM_DistL2sq(distsq,treehandle,coord1,coord2) {     \
        size_t index;                                     \
        size_t treedim; \
        long treetype;  \
        double CNEARTREE_FAR * dcoord1 = NULL; \
        double CNEARTREE_FAR * dcoord2 = NULL; \
        int CNEARTREE_FAR * icoord1 = NULL; \
        int CNEARTREE_FAR * icoord2 = NULL; \
        treedim = (treehandle)->m_szdimension; \
        treetype = (treehandle)->m_iflags&CNEARTREE_TYPE; \
        if (treetype == CNEARTREE_TYPE_DOUBLE) { \
            dcoord1 = (double CNEARTREE_FAR *)(coord1);\
            dcoord2 = (double CNEARTREE_FAR *)(coord2);\
        } else if (treetype == CNEARTREE_TYPE_INTEGER) { \
            icoord1 = (int CNEARTREE_FAR *)(coord1); \
            icoord2 = (int CNEARTREE_FAR *)(coord2); \
        } \
        distsq = 0.; \
        if (treetype == CNEARTREE_TYPE_DOUBLE) { \
            distsq = (dcoord1[0]-dcoord2[0])*(dcoord1[0]-dcoord2[0]); \
            \
            for (index=1; index < treedim; index++) { \
                distsq += (dcoord1[index]-dcoord2[index]) \
                *(dcoord1[index]-dcoord2[index]); \
            } \
        } else if (treetype == CNEARTREE_TYPE_INTEGER) { \
            distsq = ((double)(icoord1[0]-icoord2[0]))*((double)(icoord1[0]-icoord2[0])); \
            \
            for (index=1; index < treedim; index++) { \
                distsq += ((double)(icoord1[index]-icoord2[index])) \
                *((double)(icoord1[index]-icoord2[index])); \
            } \
        }  \
        }

#define CNTM_DistL2(dist,treehandle,coord1,coord2) {     \
        size_t index;                                     \
        double distsq=0.;                                    \
        size_t treedim; \
        long treetype;  \
        double CNEARTREE_FAR * dcoord1 = NULL; \
        double CNEARTREE_FAR * dcoord2 = NULL; \
        int CNEARTREE_FAR * icoord1 = NULL; \
        int CNEARTREE_FAR * icoord2 = NULL; \
        treedim = (treehandle)->m_szdimension; \
        treetype = (treehandle)->m_iflags&CNEARTREE_TYPE; \
        if (treetype == CNEARTREE_TYPE_DOUBLE) { \
            dcoord1 = (double CNEARTREE_FAR *)(coord1);\
            dcoord2 = (double CNEARTREE_FAR *)(coord2);\
        } else if (treetype == CNEARTREE_TYPE_INTEGER) { \
            icoord1 = (int CNEARTREE_FAR *)(coord1); \
            icoord2 = (int CNEARTREE_FAR *)(coord2); \
        } \
        if (treedim == 1) { \
            if (treetype == CNEARTREE_TYPE_DOUBLE) { \
                (dist) = CNTM_Abs(dcoord1[0]-dcoord2[0]); \
            } else if (treetype == CNEARTREE_TYPE_INTEGER) { \
                (dist) = CNTM_Abs((double)(icoord1[0]-icoord2[0])); \
            } \
        } else { \
            if (treetype == CNEARTREE_TYPE_DOUBLE) { \
                distsq = (dcoord1[0]-dcoord2[0])*(dcoord1[0]-dcoord2[0]); \
                \
                for (index=1; index < treedim; index++) { \
                    distsq += (dcoord1[index]-dcoord2[index]) \
                    *(dcoord1[index]-dcoord2[index]); \
                } \
            } else if (treetype == CNEARTREE_TYPE_INTEGER) { \
                distsq = ((double)(icoord1[0]-icoord2[0]))*((double)(icoord1[0]-icoord2[0])); \
                \
                for (index=1; index < treedim; index++) { \
                    distsq += ((double)(icoord1[index]-icoord2[index])) \
                    *((double)(icoord1[index]-icoord2[index])); \
                } \
            }  \
            (dist) = sqrt(distsq);\
        }\
        }
  
#define CNTM_DistL1(dist,treehandle,coord1,coord2) {\
        size_t index;   \
        size_t treedim; \
        long treetype;  \
        double CNEARTREE_FAR * dcoord1 = NULL; \
        double CNEARTREE_FAR * dcoord2 = NULL; \
        int CNEARTREE_FAR * icoord1 = NULL; \
        int CNEARTREE_FAR * icoord2 = NULL; \
        \
        treedim = treehandle->m_szdimension; \
        treetype = treehandle->m_iflags&CNEARTREE_TYPE; \
        (dist) = 0.; \
        \
        if (treetype == CNEARTREE_TYPE_DOUBLE) { \
            dcoord1 = (double CNEARTREE_FAR *)(coord1);\
            dcoord2 = (double CNEARTREE_FAR *)(coord2);\
        } else if (treetype == CNEARTREE_TYPE_INTEGER) { \
            icoord1 = (int CNEARTREE_FAR *)(coord1); \
            icoord2 = (int CNEARTREE_FAR *)(coord2); \
        } \
       if (treetype == CNEARTREE_TYPE_DOUBLE) { \
            (dist)= fabs(dcoord1[0]-dcoord2[0]); \
            for (index=1; index < treedim; index++) { \
                (dist) += CNTM_Abs(dcoord1[index]-dcoord2[index]); \
            } \
        } else if (treetype == CNEARTREE_TYPE_INTEGER) { \
            (dist) = fabs((double)(icoord1[0]-icoord2[0])); \
            for (index=1; index < treedim; index++) { \
                (dist) += CNTM_Abs((double)(dcoord1[index]-dcoord2[index])); \
            } \
        }   \
    }

    
    /*
     =======================================================================
     double CNearTreeDistsq(void CNEARTREE_FAR * coord1, 
     void CNEARTREE_FAR * coord2,  
     size_t treedim, 
     int treetype)
     
     function to return the square of the Euclidean distance between two 
     coordinate vectors.  
     
     THIS FUNCTION IS DEPRECATED
     
     treedim -- the dimension of the vectors
     
     treetype -- and integer flag for type of the vectors
     CNEARTREE_TYPE_DOUBLE for double
     CNEARTREE_TYPE_INTEGER for integer
     
     =======================================================================
     */
    
    double CNearTreeDistsq(void CNEARTREE_FAR * coord1,
                           void CNEARTREE_FAR * coord2, 
                           size_t treedim, 
                           long treetype);

     /*
     =======================================================================
     double CNearTreeDist(const CNearTreeHandle treehandle, 
     void CNEARTREE_FAR * coord1, 
     void CNEARTREE_FAR * coord2)
     
     function to return the distance (L1, L2 or L-infinity) between two 
     coordinate vectors according to the parameters of the given tree  
     For L2LAZY tree this returns the L2 norm distance.
     
     =======================================================================
     */
    
    double CNearTreeDist(const CNearTreeHandle treehandle, 
                              void CNEARTREE_FAR * coord1,
                              void CNEARTREE_FAR * coord2);
    /*
     =======================================================================
     int CNearTreeSetNorm(const CNearTreeHandle treehandle, int treenorm);
     
     function to set the norm to use used for this tree
     
     treenorm should be one of CNEARTREE_NORM_L1 for an L-1 norm
     CNEARTREE_NORM_L2 for an L-2 norm
     CNEARTREE_NORM_LINF for an L-infinity norm
     CNEARTREE_NORM_LAZY for an L-1 norm used for L-2 searching
     
     the function returns CNEARTREE_BAD_ARGUMENT for an invalid argument
     CNEARTREE_SUCCESS (0) otherwise
     =======================================================================
     */
    
    int CNearTreeSetNorm(const CNearTreeHandle treehandle, int treenorm);
    
    
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
    
    int CNearTreeNodeCreate ( const CNearTreeHandle treehandle,  
                             CNearTreeNodeHandle CNEARTREE_FAR * treenodehandle);
    
    
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
     CNEARTREE_NORM_LINF      for the max
     CNEARTREE_NORM_SPHERE    for norm as spherical angular distance
     CNEARTREE_NORM_HSPHERE   for norm as hemispherical angular distance
     CNEARTREE_NORM_HAMMING   for norm as string hamming distance
     ored with
     CNTF_NOPREPRUNE    0x10000L    flag to supress all search prepruning
     CNTF_FORCEPREPRUNE 0x20000L    flag to force search prepruning
     CNTF_NOFLIP        0x40000L    flag to suppress flips on insert
     CNTF_FORCEFLIP     0x80000L    flag to force flips on insert
     CNTF_NODEFER      0x100000L    flag to prevent deferred insert 
     
     creates an empty tree with no right or left node and with the dMax-below
     set to negative values so that any match found will be stored since it will
     greater than the negative value
     
     =======================================================================
     */
    
    int CNearTreeCreate(CNearTreeHandle CNEARTREE_FAR * treehandle, 
                        size_t treedim, 
                        long treetype);
    
    /*
     =======================================================================
     int CNearTreeFree ( CNearTreeHandle CNEARTREE_FAR * treehandle )
     
     Free a CNearTree
     
     recursively frees the NearTree with the handle *treehandle
     and nulls the treehandle.
     
     note that the objects referenced are not freed.
     =======================================================================
     */
    
    int CNearTreeFree(CNearTreeHandle CNEARTREE_FAR * treehandle);
    
    /*
     =======================================================================
     int CNearTreeClear ( const CNearTreeHandle treehandle )
     
     Clear a CNearTree
     
     Clears the NearTree with the handle *treehandle
     
     note that the objects referenced are not freed.
     =======================================================================
     */
    
    int CNearTreeClear ( const CNearTreeHandle treehandle );
    
    
    /*
     =======================================================================
     size_t CNearTreeSize (const CNearTreeHandle treehandle)
     
     Macro to get Tree Size with no error checking
     
     =======================================================================
     */
    
#define CNearTreeSize(treehandle) \
((treehandle)?(((treehandle)->m_CoordStore)?(CVectorSize((treehandle)->m_CoordStore)):0):0)
    
    /*
     =======================================================================
     int CNearTreeGetSize (const CNearTreeHandle treehandle, size_t CNEARTREE_FAR * size)
     
     Return the number of objects in the tree in size
     
     =======================================================================
     */
    
    int CNearTreeGetSize (const CNearTreeHandle treehandle, size_t CNEARTREE_FAR * size);
    
    /*
     =======================================================================
     int CNearTreeGetDelayedSize (const CNearTreeHandle treehandle, size_t CNEARTREE_FAR * size)
     
     Return the number of objects in the delay queue tree in size
     This is a deprecated alternate name for CNearTreeGetDeferredSize
     
     =======================================================================
     */
    
    int CNearTreeGetDelayedSize (const CNearTreeHandle treehandle, size_t CNEARTREE_FAR * size);

    /*
     =======================================================================
     int CNearTreeGetDeferredSize (const CNearTreeHandle treehandle, size_t CNEARTREE_FAR * size)
     
     Return the number of objects in the delay queue tree in size
     
     =======================================================================
     */
    
    int CNearTreeGetDeferredSize ( const CNearTreeHandle treehandle, size_t CNEARTREE_FAR * size );
    
    /*
     =======================================================================
     int CNearTreeGetTotalSize (const CNearTreeHandle treehandle, size_t CNEARTREE_FAR * size)
     
     Return the number of objects in both in the tree and in the delay queue tree in size
     Identical to CNearTreeGetSize, retained to support older code
     
     =======================================================================
     */
    
#define CNearTreeGetTotalSize(treehandle,size) CNearTreeGetSize(treehandle,size)
    /*
     =======================================================================
     int CNearTreeGetDepth (const CNearTreeHandle treehandle, size_t CNEARTREE_FAR * depth)
     
     Return the depth of the tree in depth
     
     =======================================================================
     */
    
    int CNearTreeGetDepth (const CNearTreeHandle treehandle, size_t CNEARTREE_FAR * depth);

    /*
     =======================================================================
     int CNearTreeGetHeight (const CNearTreeHandle treehandle, size_t CNEARTREE_FAR * height)
     
     Return the height of the tree in height, of instrumented, otherwise the depth
     
     =======================================================================
     */
    
    
    int CNearTreeGetHeight ( const CNearTreeHandle treehandle, size_t CNEARTREE_FAR * height );

    
    /*
     =======================================================================

     int CNearTreeGetFlags(const CNearTreeHandle treehandle, 
     long CNEARTREE_FAR *flags, 
     const long mask)
     
     int CNearTreeSetFlags(const CNearTreeHandle treehandle,
     const long flags,
     const long mask)
     
     functions to get or set the execution control flags:
     
        CNTF_NOPREPRUNE    to supress all search prepruning
        CNTF_FORCEPREPRUNE to force search prepruning      
        CNTF_NOFLIP        flag to suppress flips on insert
        CNTF_FORCEFLIP     flag to force flips on insert
        CNTF_NODEFER       flag to prevent deferred insert
    
     The desired flags may be ored.  If mask is no zero it is
     used as a bit mask, so a call with a flags of zero and a
     mask equal to a particular flag, clears that flag.
     =======================================================================
     */
    
    int CNearTreeGetFlags(const CNearTreeHandle treehandle, 
                          long CNEARTREE_FAR *flags, 
                          const long mask);
    
    int CNearTreeSetFlags(const CNearTreeHandle treehandle,
                          const long flags,
                          const long mask);
    
    /*
     =======================================================================
     int CNearTreeGetMeanSpacing ( const CNearTreeHandle treehandle, 
     double CNEARTREE_FAR * spacing  )
     
     Get an estimate of the spacing of points
     
     
     =======================================================================
     */
    
    int CNearTreeGetMeanSpacing ( const CNearTreeHandle treehandle, 
                                 double CNEARTREE_FAR * spacing  );
    
    /*
     =======================================================================
     int CNearTreeGetVarSpacing ( const CNearTreeHandle treehandle, 
     double CNEARTREE_FAR * varspacing  )
     
     Get an estimate of variance of the spacing of points
     
     
     =======================================================================
     */
    int CNearTreeGetVarSpacing ( const CNearTreeHandle treehandle, 
                                double CNEARTREE_FAR * varspacing  );
    
    /*
     =======================================================================
     int CNearTreeCount(const CNearTreeHandle treehandle, 
     size_t CNEARTREE_FAR * count)
     =======================================================================
     */
    
    int CNearTreeCount(const CNearTreeHandle treehandle, 
                       size_t CNEARTREE_FAR * count);
                        
    
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
                                size_t CNEARTREE_FAR * visits);
    
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
                                size_t CNEARTREE_FAR * visits);
    /*
     =======================================================================
     int CNearTreeSetNodeVisits (  const CNearTreeHandle treehandle,
     const size_t visits )
     
     Set the number of visits to nodes
     
     
     =======================================================================
     */
    int CNearTreeSetNodeVisits (  const CNearTreeHandle treehandle,
                                const size_t visits );
#endif    
    
    /*
     =======================================================================
     int CNearTreeGetDiamEstimate (  const CNearTreeHandle treehandle,
     double CNEARTREE_FAR * diamest )
     
     Get an estimate of the diameter
     
     
     =======================================================================
     */
    int CNearTreeGetDiamEstimate ( const CNearTreeHandle treehandle,
                                  double CNEARTREE_FAR * diamest );
    
    /*
     =======================================================================
     int CNearTreeGetDimEstimateEsd ( const CNearTreeHandle treehandle,
     double CNEARTREE_FAR * dimestesd )
     
     Get the current best estimate of the dimension esd
     
     =======================================================================
     */
    int CNearTreeGetDimEstimateEsd ( const CNearTreeHandle treehandle,
                                    double CNEARTREE_FAR * dimestesd );
    
    
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
                                 const double DimEstimateEsd );
    
    /*
     =======================================================================
    int CNearTreeNodeCount(const CNearTreeNodeHandle treenodehandle, 
                           size_t CNEARTREE_FAR * count)
     =======================================================================
     */
    
    int CNearTreeNodeCount(const CNearTreeNodeHandle treenodehandle, 
                           size_t CNEARTREE_FAR * count);
    /*
     =======================================================================
     int CNearTreeImmediateInsert ( const CNearTreeHandle treehandle, 
     const void CNEARTREE_FAR * coord, 
     const void CNEARTREE_FAR * obj )
     
     Function to insert some "point" as an object into a CNearTree for
     later searching, but immediately, not into the queue used for
     normal insertions
     
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
                        const void CNEARTREE_FAR * obj );
    
    /*  
     =======================================================================
     int CNearTreeNodeReInsert_Flip ( const CNearTreeHandle treehandle,
     const const CNearTreeNodeHandle treenodehandle, 
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
                                   size_t CNEARTREE_FAR * depth);
        
        
    
    
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
                            size_t CNEARTREE_FAR * depth );
    int CNearTreeNodeInsert_Flip( const CNearTreeHandle treehandle,
                            const CNearTreeNodeHandle treenodehandle,
                            const size_t index, 
                            size_t CNEARTREE_FAR * depth );    
    
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
                               const void * obj );
    
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
                               const void CNEARTREE_FAR * obj );
        
    /*
     =======================================================================
     int CNearTreeCompleteDelayedInsert ( const CNearTreeHandle treehandle )
     
     Function to dequeue the "points" queued as an objects for future insertion
     into a CNearTree for later searching
     
     return 0 for success, nonzero for an error
     
     =======================================================================
     */
    
    int CNearTreeCompleteDelayedInsert( const CNearTreeHandle treehandle );
    
    /*
     =======================================================================
     int CNearTreeZeroIfEmpty (const CNearTreeHandle treehandle)
     
     Test for an empty CNearTree, returning 0 in that case
     
     =======================================================================
     */
    
    int CNearTreeZeroIfEmpty (const CNearTreeHandle treehandle);
    
    /*
     =======================================================================
     int CNearTreeNearestNeighbor ( const CNearTreeHandle treehandle, 
     const double dRadius,  
     void CNEARTREE_FAR *  CNEARTREE_FAR * coordClosest,
     void CNEARTREE_FAR * CNEARTREE_FAR * objClosest,   
     const void CNEARTREE_FAR * coord )
     
     Function to search a Neartree for the object closest to some probe point, coord. This function
     is only here so that the function CNearTreeNearest can be called without having dRadius const
     
     dRadius is the maximum search radius - any point farther than dRadius from the probe
     point will be ignored
     
     coordClosest is the coordinate vector into which the coordinates of the nearest point
     will be stored
     
     objClosest is the address into which a pointer to the object associated with coordClosest
     will be stored
     
     coord  is the probe point
     
     the return value is true only if a point was found
     
     this version searches down the shortest branch first
     
     =======================================================================
 
     =======================================================================
     int CNearTreeLeftNearestNeighbor ( const CNearTreeHandle treehandle, 
     const double dRadius,  
     void CNEARTREE_FAR *  CNEARTREE_FAR * coordClosest,
     void CNEARTREE_FAR * CNEARTREE_FAR * objClosest,   
     const void CNEARTREE_FAR * coord )
     
     Function to search a Neartree for the object closest to some probe point, coord. This function
     is only here so that the function CNearTreeNearest can be called without having dRadius const
     
     dRadius is the maximum search radius - any point farther than dRadius from the probe
     point will be ignored
     
     coordClosest is the coordinate vector into which the coordinates of the nearest point
     will be stored
     
     objClosest is the address into which a pointer to the object associated with coordClosest
     will be stored
     
     coord  is the probe point
     
     the return value is true only if a point was found
     
     this version searches down the left branch first
     
     =======================================================================
     */
    
    int CNearTreeNearestNeighbor (const CNearTreeHandle treehandle, 
                                  const double dRadius,  
                                  void CNEARTREE_FAR * CNEARTREE_FAR * coordClosest,
                                  void CNEARTREE_FAR * CNEARTREE_FAR * objClosest, 
                                  const void CNEARTREE_FAR * coord );
    int CNearTreeLeftNearestNeighbor (const CNearTreeHandle treehandle, 
                                  const double dRadius,  
                                  void CNEARTREE_FAR * CNEARTREE_FAR * coordClosest,
                                  void CNEARTREE_FAR * CNEARTREE_FAR * objClosest, 
                                  const void CNEARTREE_FAR * coord );
    /*
     =======================================================================
     int CNearTreeFarthestNeighbor ( const CNearTreeHandle treehandle, 
     void CNEARTREE_FAR *  CNEARTREE_FAR * coordFarthest,
     void CNEARTREE_FAR * CNEARTREE_FAR * objFarthest,   
     const void CNEARTREE_FAR * coord )
     
     Function to search a Neartree for the object farthest some probe point, coord.
     
     coordClosest is the coordinate vector into which the coordinates of the nearest point
     will be stored
     
     objClosest is the address into which a pointer to the object associated with coordClosest
     will be stored
     
     coord  is the probe point
     
     the return value is true only if a point was found
     
     =======================================================================
     */
    
    int CNearTreeFarthestNeighbor (const CNearTreeHandle treehandle, 
                                   void CNEARTREE_FAR * CNEARTREE_FAR * coordFarthest,
                                   void CNEARTREE_FAR * CNEARTREE_FAR * objFarthest,   
                                   const void CNEARTREE_FAR * coord );
    /*
     =======================================================================
     int CNearTreeFindInSphere ( const CNearTreeHandle treehandle, 
     const double CNEARTREE_FAR * dRadius,
     CVectorHandle coordClosest,
     CVectorHandle objClosest,
     const void CNEARTREE_FAR * coord,
     int resetcount)
     
     Function to search a Neartree for the set of objects closer to some probe point, coord,
     than dRadius. This is only here so that objClosest can be cleared before starting the work.
     
     dRadius is the maximum search radius - any point farther than dRadius from the probe
     point will be ignored
     
     coordClosest is a vector of pointers to coordinate tuples and is the 
     returned set of nearest points to the probe point that can be found in the Neartree
     
     objClosest is a vector of objects and is the returned set of nearest points
     to the probe point that can be found in the Neartree
     
     coord  is the probe point
     
     resetcount should be non-zero to clear objClosest on entry
     
     return value is 0 if points were found
     
     =======================================================================
     */
    
    int CNearTreeFindInSphere ( const CNearTreeHandle treehandle,
                               const double dRadius,
                               CVectorHandle coordClosest,
                               CVectorHandle objClosest,
                               const void * coord,
                               int resetcount);
    
    /*
     =======================================================================
     int CNearTreeFindTreeInSphere ( const CNearTreeHandle treehandle, 
     const double CNEARTREE_FAR * dRadius,
     CNearTreeHandle foundClosest,
     const void CNEARTREE_FAR * coord,
     int resetcount)
     
     Function to search a Neartree for the set of objects closer to some probe point, coord,
     than dRadius.
     
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
                                   int resetcount);
    
    /*
     =======================================================================
     int CNearTreeFindOutSphere ( const CNearTreeHandle treehandle, 
     const double CNEARTREE_FAR * dRadius,
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
    
    int CNearTreeFindOutSphere ( const CNearTreeHandle treehandle,
                                const double dRadius,
                                CVectorHandle coordOutside,
                                CVectorHandle objOutside,
                                const void * coord,
                                int resetcount);
    
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
                                    int resetcount);
    
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
    int CNearTreeFindInAnnulus ( const CNearTreeHandle treehandle,
                                const double dRadiusInner,
                                const double dRadiusOuter,
                                CVectorHandle coordInRing,
                                CVectorHandle objInRing,
                                const void CNEARTREE_FAR * coord,
                                int resetcount);
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
                                    int resetcount);
    
    
    /*
     =======================================================================

     int CNearTreeSortIn(CVectorHandle metrics, 
     CVectorHandle indices, 
     double metric, 
     size_t index, 
     size_t k);
     
     CNearTreeSortIn inserts a new metric and index into the vectors
     metrics and indices, sorted on non-decreasing metric,
     with the size of the vectors capped at k, or uncapped if k = 0;

     =======================================================================
     */

    int CNearTreeSortIn(CVectorHandle metrics,
                        CVectorHandle indices,
                        double metric,
                        size_t index,
                        size_t k);
    /*
     =======================================================================
     
     int CNearTreeSortIn2(CVectorHandle metrics,
     CVectorHandle indices,
     double metric,
     size_t index,
     size_t k,
     int shell,
     size_t sSizeCur,
     double dRadiusCur);
     
     CNearTreeSortIn2 inserts a new metric and index into the vectors
     metrics and indices, sorted on non-decreasing metric,
     with the size of the vectors capped at k, or uncapped if k = 0;
     
     =======================================================================
     */
    
    int CNearTreeSortIn2(CVectorHandle metrics,
                         CVectorHandle indices,
                         double metric,
                         size_t index,
                         size_t k,
                         int shell,
                         size_t sSizeCur,
                         double dRadiusCur);
    
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
                               int resetcount);
   
    int CNearTreeFindKNearest_Annular (const CNearTreeHandle treehandle,
                                       const size_t k,
                                       const double dRadius,
                                       CVectorHandle distances,
                                       CVectorHandle indices,
                                       const void * coord);

    int CNearTreeFindKNearest_Sphere (const CNearTreeHandle treehandle,
                                      const size_t k,
                                      const double dRadius,
                                      CVectorHandle distances,
                                      CVectorHandle indices,
                                      const void * coord);

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
                                     const void * coord);
    
    
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
     
     coordClosest is the returned closest point
     to the probe point that can be found in the Neartree
     
     objClosest is a pointer to a pointer to hold the corresponding object or is NULL
     
     coord  is the probe point
     
     the return value is 0 only if a point was found within dRadius
     
     this version search down the shortest branch first
     
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
     
     coordClosest is the returned closest point
     to the probe point that can be found in the Neartree
     
     objClosest is a pointer to a pointer to hold the corresponding object or is NULL
     
     coord  is the probe point
     
     the return value is 0 only if a point was found within dRadius
     
     this version searches down the left branch first
     
     =======================================================================
     
     */
    int CNearTreeNearest ( const CNearTreeHandle treehandle, 
                          double CNEARTREE_FAR * dRadius,  
                          void CNEARTREE_FAR * CNEARTREE_FAR * coordClosest,
                          void CNEARTREE_FAR * CNEARTREE_FAR * objClosest,
                          const void CNEARTREE_FAR * coord );

    int CNearTreeLeftNearest ( const CNearTreeHandle treehandle, 
                          double CNEARTREE_FAR * dRadius,  
                          void CNEARTREE_FAR * CNEARTREE_FAR * coordClosest,
                          void CNEARTREE_FAR * CNEARTREE_FAR * objClosest,
                          const void CNEARTREE_FAR * coord );
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
     
     coordFarthest is the returned farthest point
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
                               const void CNEARTREE_FAR * coord );
    
    
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
     
     dRadius is the maximum search radius - any point closer than dRadius to the probe
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
    
    int CNearTreeFindKNearest ( const CNearTreeHandle treehandle,
                               const size_t k,
                               const double dRadius,
                               CVectorHandle coordClosest,
                               CVectorHandle objClosest,
                               const void * coord,
                               int resetcount);
    
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
    int CNearTreeFindKTreeNearest ( const CNearTreeHandle treehandle,
                                   const size_t k,
                                   const double dRadius,
                                   CNearTreeHandle foundClosest,
                                   const void CNEARTREE_FAR * coord,
                                   int resetcount);
    
    
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
                                int resetcount);
    
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
                                    int resetcount);
    
    
    
    
    /*
     =======================================================================
     The following macro is provided here to ensure operation with older
     versions of CVector
     =======================================================================
     */
    
#ifndef CVectorElementAt
    /* CVectorElementAt -- return the element at the given index as a void pointer without checking
     and without protection against relocation */
    
#define CVectorElementAt(vectorhandle,index) ((void CNEARTREE_FAR *)(((char *)((vectorhandle)->array))+(index)*(vectorhandle)->elementsize))
#endif
    
    
    /*
     =======================================================================
     int CNearTreeObjects ( const CNearTreeHandle treehandle, CVectorHandle CNEARTREE_FAR * vectorhandle)
     
     Function to return the vector of objects in the tree.  This vector
     is not guaranteed to be in the same order as the order of insertion
     
     vectorhandle -- a pointer to a CVectorHandle
     
     =======================================================================
     */
    
    int CNearTreeObjects ( const CNearTreeHandle treehandle, CVectorHandle CNEARTREE_FAR * vectorhandle);    
    
    /*
     =======================================================================
     int CNearTreeCoords ( const CNearTreeHandle treehandle, CVectorHandle CNEARTREE_FAR * vectorhandle)
     
     Function to return the vector of coordinates in the tree.  This vector
     is not guaranteed to be in the same order as the order of insertion
     
     vectorhandle -- a pointer to a CVectorHandle
     
     =======================================================================
     */
    
    int CNearTreeCoords ( const CNearTreeHandle treehandle, CVectorHandle CNEARTREE_FAR * vectorhandle);    
    
    /*
     =======================================================================
     void CNEARTREE_FAR * CNearTreeObjectAt ( const CNearTreeHandle treehandle, size_t index)
     
     Function to return the an object pointer at the given index
     
     implemented as a macro
     
     =======================================================================
     */
    
#define CNearTreeObjectAt(treehandle,index) CVectorElementAt(treehandle->m_ObjectStore,index)
    
    /*
     =======================================================================
     void CNEARTREE_FAR * CNearTreeCoordAt ( const CNearTreeHandle treehandle, size_t index)
     
     Function to return the a coordinate pointer at the given index
     
     implemented as a macro
     
     =======================================================================
     */
    
#define CNearTreeCoordAt(treehandle,index) CVectorElementAt(treehandle->m_CoordStore,index)
    
    
    
#ifdef __cplusplus
    
}

#endif


#endif
