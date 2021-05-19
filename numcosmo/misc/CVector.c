/*
 *  CVector.c
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

#ifdef __cplusplus

extern "C" {
    
#endif
    
#include "CVector.h"



    
    /* CVectorCreate -- create a generic vector */
    
    int CVectorCreate(CVectorHandle CVECTOR_FAR * vectorhandle, const size_t elementsize, const size_t capacity) {
        
        size_t cap = capacity;
        
        if ((vectorhandle==NULL)) { return CVECTOR_BAD_ARGUMENT; }
        
        *vectorhandle = (CVectorHandle)CVECTOR_MALLOC(sizeof(CVector));
        
        if ((*vectorhandle==NULL)) {
            return CVECTOR_MALLOC_FAILED;
        }
        
        (*vectorhandle)->size = 0;
        (*vectorhandle)->flags = 0;
        (*vectorhandle)->capacity = 0;
        (*vectorhandle)->elementsize = elementsize;
        if (!cap) { cap = 10; }
        (*vectorhandle)->array = (void CVECTOR_FAR *)CVECTOR_MALLOC(cap*elementsize);
        if ((*vectorhandle)->array) {
            (*vectorhandle)->capacity = cap;
            return 0;
        }
        CVECTOR_FREE(*vectorhandle);
        *vectorhandle = NULL;
        return CVECTOR_MALLOC_FAILED;
    }
    
    
    /*  CVectorAddElement -- add an element to a generic vector 
     equivalent to vector::push_back  */
    
    
    int CVectorAddElement(const CVectorHandle vectorhandle, const void CVECTOR_FAR * element) {
        
        size_t newcap;
        
        int errorcode;
        
        if ((vectorhandle==NULL)) { return CVECTOR_BAD_ARGUMENT; }
        
        if ( (vectorhandle->flags&CVECTOR_FLAGS_NO_RESIZE) ) { return CVECTOR_NO_RESIZE; }
        
        if (vectorhandle->size >= vectorhandle->capacity) {
            
            newcap = vectorhandle->capacity*2;
            
            if (newcap < 1) { newcap = 1; }
            
            errorcode = CVectorSetCapacity (vectorhandle, newcap);
            
            if (errorcode != 0) {
                
                newcap = vectorhandle->capacity+(size_t)(vectorhandle->capacity>>2);
                
                if (newcap < 1) { newcap = 1; }
                
                errorcode = CVectorSetCapacity (vectorhandle, newcap);
                
                if (errorcode != 0) { return errorcode; }
            }
            
        }
        
        CVECTOR_MEMMOVE(((char CVECTOR_FAR *)(vectorhandle->array))+vectorhandle->size*vectorhandle->elementsize,
                (const char CVECTOR_FAR *)element, vectorhandle->elementsize);
        vectorhandle->size ++;
        return 0;
        
    }
    
    /* CVectorGetElement -- get a copy of an element from a generic vector */
    
    int CVectorGetElement(const CVectorHandle vectorhandle, void CVECTOR_FAR * element, const size_t index) {
        
        if ((vectorhandle==NULL)) { return CVECTOR_BAD_ARGUMENT; }
        
        if (index >= 0 && index < vectorhandle->size) {
            
            CVECTOR_MEMMOVE((char *)element,((char *)(vectorhandle->array))+index*vectorhandle->elementsize,
                    vectorhandle->elementsize);
            
            return 0;
            
        } else {
            
            return CVECTOR_NOT_FOUND;
        }
        
    }
    
    
    /* CVectorGetElementptr -- get a pointer to an element from a generic vector */
    
    int CVectorGetElementptr(const CVectorHandle vectorhandle, void CVECTOR_FAR ** elementptr, const size_t index) {
        
        if ((vectorhandle==NULL)) { return CVECTOR_BAD_ARGUMENT; }
        
        if (index >= 0 && index < vectorhandle->size) {
            
            *elementptr = (void CVECTOR_FAR*)(((char *)(vectorhandle->array))+index*vectorhandle->elementsize);
            
            vectorhandle->flags |= CVECTOR_FLAGS_NO_RELOCATION;
            
            return 0;
            
        } else {
            
            return CVECTOR_NOT_FOUND;
        }
        
    }
    
    
    /* CVectorSetElement -- set a copy of an element into a generic vector */
    
    int CVectorSetElement(const CVectorHandle vectorhandle, const void CVECTOR_FAR * element, const size_t index) {
        
        size_t newcap;
        
        int errorcode;
        
        if ((vectorhandle==NULL) ) { return CVECTOR_BAD_ARGUMENT; }
        
        if (index >= vectorhandle->capacity) {
            
            newcap = index+vectorhandle->capacity+1;
            
            errorcode = CVectorSetCapacity(vectorhandle, newcap);
            
            if (errorcode != 0 ) {
                
                newcap = index*1.2;
                if (newcap < index+128) { newcap = index+128; }
                
                errorcode = CVectorSetCapacity(vectorhandle, newcap);
                
                if (errorcode != 0) { return errorcode; }
                
            }
        }
        
        
        if (index >= 0 && index < vectorhandle->capacity) {
            
            CVECTOR_MEMMOVE(((char *)(vectorhandle->array))+index*vectorhandle->elementsize,(char *)element,
                    vectorhandle->elementsize);
            
            if (index >= vectorhandle->size) { vectorhandle->size = index+1; }
            return 0;
            
        } else {
            
            return CVECTOR_NOT_FOUND;
        }
        
    }
    
    /* CVectorRemoveElement -- remove an element from a generic vector */
    
    /* keeps elements 0 -- index-1, discards element index
     moves elements index+1 through vectorhandle->size-1
     into element index through vectorhandle->size-2
     
     i.e. moves characters (index+1)*(vectorhandle->elementsize)
     through (vectorhandle->size)*(vectorhandle->elementsize)-1
     to index*(vectorhandle->elementsize)
     */
    
    int CVectorRemoveElement(const CVectorHandle vectorhandle, const size_t index) {
        
        if ((vectorhandle==NULL)) { return CVECTOR_BAD_ARGUMENT; }
        
        if ((vectorhandle->flags&CVECTOR_FLAGS_NO_RELOCATION)) { return CVECTOR_NO_RELOCATION; }
        
        if (index >= vectorhandle->size || index < 0 ) { return CVECTOR_NOT_FOUND; }
        
        if (index == vectorhandle->size-1) {
            vectorhandle->size--;
            return 0;
        }
        
        CVECTOR_MEMMOVE((char *)vectorhandle->array+index*(vectorhandle->elementsize),
                (char *)vectorhandle->array+(index+1)*(vectorhandle->elementsize),
                (vectorhandle->size-1-index)*(vectorhandle->elementsize));
        
        vectorhandle->size--;
        return 0;
    }
    
    /* CVectorClear -- clear a generic vector */
    
    int CVectorClear(const CVectorHandle vectorhandle) {
        
        if ((vectorhandle==NULL)) { return CVECTOR_BAD_ARGUMENT; }
        
        if (vectorhandle->flags & CVECTOR_FLAGS_NO_RESIZE) { return CVECTOR_NO_RESIZE; }
        
        vectorhandle->size = 0;
        
        return 0;
        
    }
    
    /* CVectorFree -- remove a generic vector */
    
    int CVectorFree(CVectorHandle CVECTOR_FAR * vectorhandle) {
        
        if ((vectorhandle==NULL)) { return CVECTOR_BAD_ARGUMENT; }
        
        if (*vectorhandle) {
            
            if ((*vectorhandle)->flags & CVECTOR_FLAGS_NO_RESIZE) { return CVECTOR_NO_RESIZE; }
            
            if ((*vectorhandle)->flags & CVECTOR_FLAGS_NO_RELOCATION) { return CVECTOR_NO_RELOCATION; }
            
            if ((*vectorhandle)->array) {
                
                CVECTOR_FREE((*vectorhandle)->array);
                
            }
            
            CVECTOR_FREE(*vectorhandle);
            
        }
        
        *vectorhandle = 0;
        
        return 0;
        
    }
    
    /* CVectorGetCapacity - function to return the CVector capacity */
    
    int CVectorGetCapacity(const CVectorHandle vectorhandle, size_t CVECTOR_FAR * capacity) {
        
        if ((vectorhandle==NULL)||!(capacity)) { return CVECTOR_BAD_ARGUMENT; }
        
        *capacity = vectorhandle->capacity;
        
        return 0;
    }
    
    /* CVectorGetSize - function to return the CVector size */
    
    int CVectorGetSize(const CVectorHandle vectorhandle, size_t CVECTOR_FAR * size) {
        
        if ((vectorhandle==NULL)||!(size)) { return CVECTOR_BAD_ARGUMENT; }
        
        *size = vectorhandle->size;
        
        return 0;
    }
    
    /* CVectorGetFlags - function to return the CVector flags */
    
    int CVectorGetFlags(const CVectorHandle vectorhandle, unsigned int CVECTOR_FAR * flags) {
        
        if ((vectorhandle==NULL)||!(flags)) { return CVECTOR_BAD_ARGUMENT; }
        
        *flags = vectorhandle->flags;
        
        return 0;
    }
    
    /* CVectorSetCapacity - function to set the CVector capacity */
    
    int CVectorSetCapacity(const CVectorHandle vectorhandle, const size_t capacity) {
        
        void CVECTOR_FAR * temparray;
        
        if ((vectorhandle==NULL) || capacity < vectorhandle->size) { return CVECTOR_BAD_ARGUMENT; }
        
        if (capacity == vectorhandle->capacity) { return 0; }
        
        if (vectorhandle->flags&CVECTOR_FLAGS_NO_RELOCATION) { return CVECTOR_NO_RELOCATION; }
        
        if (capacity) {
            
            temparray = CVECTOR_MALLOC(capacity*vectorhandle->elementsize);
            if (!temparray)   { return CVECTOR_MALLOC_FAILED; }
            
            if (vectorhandle->size) {   
              CVECTOR_MEMMOVE((char *)temparray, (char *)vectorhandle->array, vectorhandle->size*vectorhandle->elementsize); 
            }
            CVECTOR_FREE(vectorhandle->array);
            
        } else {
            temparray = NULL;
        }
        
        vectorhandle->array = temparray;
        vectorhandle->capacity = capacity;
        
        return 0;
    }
    
    /* CVectorSetSize - function to set the CVector size */
    
    int CVectorSetSize(const CVectorHandle vectorhandle, const size_t size) {
        
        int errorcode;
        
        if ((vectorhandle==NULL) ) { return CVECTOR_BAD_ARGUMENT; }
        
        if (size == vectorhandle->size) { return 0; }
        
        if ((vectorhandle->flags & CVECTOR_FLAGS_NO_RESIZE)) { return CVECTOR_NO_RESIZE; }
        
        if ( size > vectorhandle->capacity ) {
            
            if ((vectorhandle->flags & CVECTOR_FLAGS_NO_RELOCATION) ) { return CVECTOR_NO_RELOCATION; }
            
            errorcode = CVectorSetCapacity(vectorhandle,size);
            
            if (errorcode != 0) { return errorcode; }
            
        }
        
        if ( size <= vectorhandle->capacity ) {
            
            if (size > vectorhandle->size) {
                
                CVECTOR_MEMSET(((char *)vectorhandle->array)+(vectorhandle->size)*(vectorhandle->elementsize),
                       0, (vectorhandle->size-size)*(vectorhandle->elementsize));
            }
            
            vectorhandle->size = size;
            
        }
        
        
        return 0;
    }
    
    
    /* CVectorSetFags - function to set the CVector flags */
    
    int CVectorSetFlags(const CVectorHandle vectorhandle, const unsigned int flags) {
        
        if ((vectorhandle==NULL) ) { return CVECTOR_BAD_ARGUMENT; }
        
        vectorhandle->flags = flags;
        
        return 0;
        
    }
    
#ifdef __cplusplus
    
}

#endif


