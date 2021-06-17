/*
 *  CVector.h
 *  CVector
 *
 *  Created by Herbert J. Bernstein on 11/28/08.
 *  Copyright 2008 Herbert J. Bernstein. All rights reserved.
 *
 */

/**********************************************************************
 *                                                                    *
 * YOU MAY REDISTRIBUTE THE CVector API UNDER THE TERMS OF THE LGPL   *
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

#ifndef CVECTOR_H_INCLUDED
#define CVECTOR_H_INCLUDED

#ifdef __cplusplus

extern "C" {
    
#endif
    
    /* CVector flags */
    
#define CVECTOR_FLAGS_NO_RELOCATION  1
#define CVECTOR_FLAGS_NO_RESIZE      2
    
    
    /* CVector error return values */
    
#define CVECTOR_MALLOC_FAILED -1
#define CVECTOR_BAD_ARGUMENT  -2
#define CVECTOR_NOT_FOUND     -4
#define CVECTOR_NO_RELOCATION -8
#define CVECTOR_NO_RESIZE     -16
    
#ifdef CVECTOR_USE_FAR
#include <malloc.h>
#include <string.h>
#define CVECTOR_FAR __far
#define CVECTOR_MALLOC _fmalloc
#define CVECTOR_FREE _ffree
#define CVECTOR_MEMSET _fmemset
#define CVECTOR_MEMMOVE _fmemmove
#else
#include <stdlib.h>
#include <string.h>
#define CVECTOR_FAR
#define CVECTOR_MALLOC malloc
#define CVECTOR_FREE free
#define CVECTOR_MEMSET memset
#define CVECTOR_MEMMOVE memmove
#endif
    
    typedef struct {
        size_t size;             /* size of the vector      */
        size_t capacity;         /* capacity of the vector  */
        size_t elementsize;      /* size of an element      */
        void CVECTOR_FAR * array;        /* the array of elements   */
        unsigned int flags;      /* flags                   */
    } CVector;
    
    typedef CVector CVECTOR_FAR * CVectorHandle;
    
    
    /*  CVectorAddElement -- add an element to a CVector */
    
    int CVectorAddElement(const CVectorHandle vectorhandle, const void CVECTOR_FAR * element);

    /* CVectorCapacity -- macro to return the CVector capacity */
    
#define CVectorCapacity(vectorhandle)  (vectorhandle)->capacity
    
    /* CVectorClear -- clear a generic vector */
    
    int CVectorClear(const CVectorHandle vectorhandle);

    /* CVectorCreate -- create a CVector */
    
    int CVectorCreate(CVectorHandle CVECTOR_FAR * vectorhandle, const size_t elementsize, const size_t capacity);
    
    /* CVectorElementAt -- return the element at the given index as a void pointer without checking
       and without protection against relocation */
        
#define CVectorElementAt(vectorhandle,index) ((void CVECTOR_FAR *)(((char *)((vectorhandle)->array))+(index)*(vectorhandle)->elementsize))

    /* CVectorFree -- remove a CVector */
    
    int CVectorFree(CVectorHandle CVECTOR_FAR * vectorhandle);
    
    /* CVectorGetCapacity - function to return the CVector capacity */
    
    int CVectorGetCapacity(const CVectorHandle vectorhandle, size_t CVECTOR_FAR * capacity);
    
    /* CVectorGetElement -- get a copy of an element from a CVector */
    
    int CVectorGetElement(const CVectorHandle vectorhandle, void CVECTOR_FAR * element, const size_t index);
    
    /* CVectorGetElementptr -- get a pointer to an element from a CVector */
    
    int CVectorGetElementptr(const CVectorHandle vectorhandle, void CVECTOR_FAR ** elementptr, const size_t index);
        
    /* CVectorGetFlags - function to return the CVector flags */
    
    int CVectorGetFlags(const CVectorHandle vectorhandle, unsigned int CVECTOR_FAR * flags);
    
    /* CVectorGetSize - function to return the CVector size */
    
    int CVectorGetSize(const CVectorHandle vectorhandle, size_t CVECTOR_FAR * size);
        
    /* CVectorRemoveElement -- remove an element from a generic vector */
    
    int CVectorRemoveElement(const CVectorHandle vectorhandle, const size_t index);

    /* CVectorSetCapacity - function to set the CVector capacity */
    
    int CVectorSetCapacity(const CVectorHandle vectorhandle, const size_t capacity);       
        
    /* CVectorSetElement -- set a copy of an element into a CVector */
    
    int CVectorSetElement(const CVectorHandle vectorhandle, const void CVECTOR_FAR * element, const size_t index);

    /* CVectorSetFags - function to set the CVector flags */
    
    int CVectorSetFlags(const CVectorHandle vectorhandle, const unsigned int flags);
    
    /* CVectorSetSize - function to set the CVector size */
    
    int CVectorSetSize(const CVectorHandle vectorhandle, const size_t size);

    /* CVectorSize -- macro to return the CVector size */
    
#define CVectorSize(vectorhandle)  (vectorhandle)->size
    


#ifdef __cplusplus
    
}

#endif


#endif
