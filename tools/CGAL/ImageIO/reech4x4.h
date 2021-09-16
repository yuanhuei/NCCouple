// Copyright (c) 2005-2008 ASCLEPIOS Project, INRIA Sophia-Antipolis (France)
// All rights reserved.
//
// This file is part of the ImageIO Library, and as been adapted for
// CGAL (www.cgal.org).
// You can redistribute it and/or  modify it under the terms of the
// GNU Lesser General Public License as published by the Free Software Foundation;
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// These files are provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: https://github.com/CGAL/cgal/blob/releases/CGAL-4.14.3/CGAL_ImageIO/include/CGAL/ImageIO/reech4x4.h $
// $Id: reech4x4.h 0698f79 2017-10-20T23:34:14+02:00 Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     :  ASCLEPIOS Project (INRIA Sophia-Antipolis), Laurent Rineau

/*************************************************************************
 * reech4x4.h -
 *
 * $Id: reech4x4.h 0698f79 2017-10-20T23:34:14+02:00 Sébastien Loriot
 *
 * Copyright©INRIA 1999
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 *
 *
 * ADDITIONS, CHANGES
 *	
 *	
 *	
 *
 */


/* CAUTION
   DO NOT EDIT THIS FILE,
   UNLESS YOU HAVE A VERY GOOD REASON 
 */

#ifndef _reech4x4_h_
#define _reech4x4_h_









extern void Reech3DTriLin4x4_u8 ( void* theBuf, /* buffer to be resampled */
			     int *theDim,  /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat   /* transformation matrix */
			     );
extern void Reech3DTriLin4x4gb_u8 ( void* theBuf, /* buffer to be resampled */
			     int *theDim, /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat,   /* transformation matrix */
			     float gain,
			     float bias );
extern void Reech3DNearest4x4_u8 ( void* theBuf, /* buffer to be resampled */
			      int *theDim,  /* dimensions of this buffer */
			      void* resBuf, /* result buffer */
			      int *resDim,  /* dimensions of this buffer */
			      double* mat   /* transformation matrix */
			      );
extern void Reech2DTriLin4x4_u8 ( void* theBuf, /* buffer to be resampled */
			     int *theDim,  /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat   /* transformation matrix */
			     );
extern void Reech2DTriLin4x4gb_u8 ( void* theBuf, /* buffer to be resampled */
			     int *theDim, /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat,   /* transformation matrix */
			     float gain,
			     float bias );
extern void Reech2DNearest4x4_u8 ( void* theBuf, /* buffer to be resampled */
			      int *theDim,  /* dimensions of this buffer */
			      void* resBuf, /* result buffer */
			      int *resDim,  /* dimensions of this buffer */
			      double* mat   /* transformation matrix */
			      );






extern void Reech3DTriLin4x4_s8 ( void* theBuf, /* buffer to be resampled */
			     int *theDim,  /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat   /* transformation matrix */
			     );
extern void Reech3DTriLin4x4gb_s8 ( void* theBuf, /* buffer to be resampled */
			     int *theDim, /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat,   /* transformation matrix */
			     float gain,
			     float bias );
extern void Reech3DNearest4x4_s8 ( void* theBuf, /* buffer to be resampled */
			      int *theDim,  /* dimensions of this buffer */
			      void* resBuf, /* result buffer */
			      int *resDim,  /* dimensions of this buffer */
			      double* mat   /* transformation matrix */
			      );
extern void Reech2DTriLin4x4_s8 ( void* theBuf, /* buffer to be resampled */
			     int *theDim,  /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat   /* transformation matrix */
			     );
extern void Reech2DTriLin4x4gb_s8 ( void* theBuf, /* buffer to be resampled */
			     int *theDim, /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat,   /* transformation matrix */
			     float gain,
			     float bias );
extern void Reech2DNearest4x4_s8 ( void* theBuf, /* buffer to be resampled */
			      int *theDim,  /* dimensions of this buffer */
			      void* resBuf, /* result buffer */
			      int *resDim,  /* dimensions of this buffer */
			      double* mat   /* transformation matrix */
			      );






extern void Reech3DTriLin4x4_u16 ( void* theBuf, /* buffer to be resampled */
			     int *theDim,  /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat   /* transformation matrix */
			     );
extern void Reech3DTriLin4x4gb_u16 ( void* theBuf, /* buffer to be resampled */
			     int *theDim, /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat,   /* transformation matrix */
			     float gain,
			     float bias );
extern void Reech3DNearest4x4_u16 ( void* theBuf, /* buffer to be resampled */
			      int *theDim,  /* dimensions of this buffer */
			      void* resBuf, /* result buffer */
			      int *resDim,  /* dimensions of this buffer */
			      double* mat   /* transformation matrix */
			      );
extern void Reech2DTriLin4x4_u16 ( void* theBuf, /* buffer to be resampled */
			     int *theDim,  /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat   /* transformation matrix */
			     );
extern void Reech2DTriLin4x4gb_u16 ( void* theBuf, /* buffer to be resampled */
			     int *theDim, /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat,   /* transformation matrix */
			     float gain,
			     float bias );
extern void Reech2DNearest4x4_u16 ( void* theBuf, /* buffer to be resampled */
			      int *theDim,  /* dimensions of this buffer */
			      void* resBuf, /* result buffer */
			      int *resDim,  /* dimensions of this buffer */
			      double* mat   /* transformation matrix */
			      );






extern void Reech3DTriLin4x4_s16 ( void* theBuf, /* buffer to be resampled */
			     int *theDim,  /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat   /* transformation matrix */
			     );
extern void Reech3DTriLin4x4gb_s16 ( void* theBuf, /* buffer to be resampled */
			     int *theDim, /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat,   /* transformation matrix */
			     float gain,
			     float bias );
extern void Reech3DNearest4x4_s16 ( void* theBuf, /* buffer to be resampled */
			      int *theDim,  /* dimensions of this buffer */
			      void* resBuf, /* result buffer */
			      int *resDim,  /* dimensions of this buffer */
			      double* mat   /* transformation matrix */
			      );
extern void Reech2DTriLin4x4_s16 ( void* theBuf, /* buffer to be resampled */
			     int *theDim,  /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat   /* transformation matrix */
			     );
extern void Reech2DTriLin4x4gb_s16 ( void* theBuf, /* buffer to be resampled */
			     int *theDim, /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat,   /* transformation matrix */
			     float gain,
			     float bias );
extern void Reech2DNearest4x4_s16 ( void* theBuf, /* buffer to be resampled */
			      int *theDim,  /* dimensions of this buffer */
			      void* resBuf, /* result buffer */
			      int *resDim,  /* dimensions of this buffer */
			      double* mat   /* transformation matrix */
			      );






extern void Reech3DTriLin4x4_r32 ( void* theBuf, /* buffer to be resampled */
			     int *theDim,  /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat   /* transformation matrix */
			     );
extern void Reech3DTriLin4x4gb_r32 ( void* theBuf, /* buffer to be resampled */
			     int *theDim, /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat,   /* transformation matrix */
			     float gain,
			     float bias );
extern void Reech3DNearest4x4_r32 ( void* theBuf, /* buffer to be resampled */
			      int *theDim,  /* dimensions of this buffer */
			      void* resBuf, /* result buffer */
			      int *resDim,  /* dimensions of this buffer */
			      double* mat   /* transformation matrix */
			      );
extern void Reech2DTriLin4x4_r32 ( void* theBuf, /* buffer to be resampled */
			     int *theDim,  /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat   /* transformation matrix */
			     );
extern void Reech2DTriLin4x4gb_r32 ( void* theBuf, /* buffer to be resampled */
			     int *theDim, /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat,   /* transformation matrix */
			     float gain,
			     float bias );
extern void Reech2DNearest4x4_r32 ( void* theBuf, /* buffer to be resampled */
			      int *theDim,  /* dimensions of this buffer */
			      void* resBuf, /* result buffer */
			      int *resDim,  /* dimensions of this buffer */
			      double* mat   /* transformation matrix */
			      );







extern void Reech4x4_verbose ( );
extern void Reech4x4_noverbose ( );

#ifdef CGAL_HEADER_ONLY
#include <CGAL/ImageIO/reech4x4_impl.h>
#endif // CGAL_HEADER_ONLY


#endif
