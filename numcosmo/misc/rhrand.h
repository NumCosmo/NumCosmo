/*
 *  rhrand.h
 *  RHrand
 *  
 *  Based on:
 *  Random Number generator by Rob Harrison, derived from 
 *  "one in J.M.Hammersley and D.C. Handscomb, "Monte Carlo 
 *  Methods," Methuen & Co., London and Wiley & Sons, New 
 *  York, 1964, p47".  See also, D. E. Knuth "The Art of 
 *  Computer Programming", Volume 2, "Seminumerical
 *  Alogorithms, Third Edition, Addison-Wesley, Reading MA,
 *  1997.
 *
 *  This adaptation by H. J. Bernstein
 *  Copyright 2009 Rob Harrison, Larry Andrews and Herbert J. Bernstein
 *
 */

/**********************************************************************
 *                                                                    *
 * YOU MAY REDISTRIBUTE THE rhrand API UNDER THE TERMS OF THE LGPL    *
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

#ifndef RHRAND_H_INCLUDED
#define RHRAND_H_INCLUDED

#ifdef __cplusplus

#include <cmath>

class RHrand {
    
private:
    double randomNumberBuffer[56];
    int indx;
    int jndx;
    int kndx;
    double dTemp;
    
public:
    static const int RHRAND_MAX = 32767;
    
    RHrand( void ) {                 /* default constructor */
        srandom(0);
    };
    
    
    RHrand( const int iseed ) {            /* constructor with seed */
        srandom(iseed);
    };
    
    ~RHrand(void ) {                /* nothing to do for destructor */ 
    };
    
    void srandom ( const int iseed ) {
        jndx = iseed;
        if (jndx < 0) jndx = -jndx;
        for( indx=0; indx<(int)(sizeof(randomNumberBuffer)/sizeof(randomNumberBuffer[0])); ++indx )
        {
            jndx = (jndx*2349 + 14867)%32767;
            randomNumberBuffer[indx] = fabs((double)(jndx/32767.0));
        }
        indx = 55;
        kndx = 54;
        jndx = 31;
        return;
    };
    
    double urand( void ) {
        indx = indx%55 + 1;
        jndx = jndx%55 + 1;
        kndx = kndx%55 + 1;
        randomNumberBuffer[indx-1] = modf( randomNumberBuffer[jndx-1]+randomNumberBuffer[kndx-1], &dTemp );
        return  randomNumberBuffer[indx-1];
    };
    
    int random( void ) {
        return (int)(urand()*((double)RHRAND_MAX));
    };
    
};

#ifndef RHRAND_NOCCODE
extern "C" {
#endif
    
#endif

#ifndef RHRAND_NOCCODE

#include <math.h>
    
    typedef struct CRHrand_ {
        double buffer[56];
        int indx;
        int jndx;
        int kndx;
        double dTemp;
    } CRHrand;
    
    typedef CRHrand * CRHrandHandle;
    
#define CRHRAND_MAX 32767

#define CRHrandSrandom(randhandle,iseed) { \
        (randhandle)->jndx = iseed; \
        if ((randhandle)->jndx < 0) (randhandle)->jndx = -(randhandle)->jndx; \
        for((randhandle)->indx=0; (randhandle)->indx<(sizeof((randhandle)->buffer)/sizeof((randhandle)->buffer[0])); ++(randhandle)->indx ) \
        { \
            (randhandle)->jndx = ((randhandle)->jndx*2349 + 14867)%32767; \
            ((randhandle)->buffer)[(randhandle)->indx] = fabs((double)(((randhandle)->jndx)/32767.0)); \
        } \
        (randhandle)->indx = 55; \
        (randhandle)->kndx = 54; \
        (randhandle)->jndx = 31; \
    }
#define CRHrandUrand(randhandle) ( \
        (randhandle)->indx = (randhandle)->indx%55 + 1, \
        (randhandle)->jndx = (randhandle)->jndx%55 + 1, \
        (randhandle)->kndx = (randhandle)->kndx%55 + 1, \
        ((randhandle)->buffer)[(randhandle)->indx-1] \
        = modf( ((randhandle)->buffer)[(randhandle)->jndx-1]+ \
               ((randhandle)->buffer)[(randhandle)->kndx-1], &(randhandle)->dTemp ) , \
        ((randhandle)->buffer)[(randhandle)->indx-1]  ) 

#define CRHrandRandom(randhandle) ((int)(CRHrandUrand(randhandle)*(double)CRHRAND_MAX))
    
    
    
#ifdef __cplusplus
    
}

#endif
#endif


#endif


