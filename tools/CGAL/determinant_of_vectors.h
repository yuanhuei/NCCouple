// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: https://github.com/CGAL/cgal/blob/releases/CGAL-4.14.3/NewKernel_d/include/CGAL/determinant_of_vectors.h $
// $Id: determinant_of_vectors.h 0698f79 2017-10-20T23:34:14+02:00 Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Marc Glisse

#ifndef CGAL_DETVEC_H
#define CGAL_DETVEC_H
#include <CGAL/determinant.h>
#include <CGAL/predicates/sign_of_determinant.h>

namespace CGAL {
  // TODO: determine whether it is better to pass them by lines or columns.

  template <class NT, class Vector> inline
  NT determinant_of_vectors(Vector const&a, Vector const&b){
    return determinant<NT>(a[0],a[1],b[0],b[1]);
  }
  template <class NT, class Vector> inline
  typename Sgn<NT>::result_type
  sign_of_determinant_of_vectors(Vector const&a, Vector const&b){
    return sign_of_determinant<NT>(a[0],a[1],b[0],b[1]);
  }

  template <class NT, class Vector>
  NT determinant_of_vectors(Vector const&a, Vector const&b,
      Vector const&c){
    return determinant<NT>(a[0],a[1],a[2],b[0],b[1],b[2],c[0],c[1],c[2]);
  }
  template <class NT, class Vector>
  typename Sgn<NT>::result_type
  sign_of_determinant_of_vectors(Vector const&a, Vector const&b,
      Vector const&c){
    return sign_of_determinant<NT>(a[0],a[1],a[2],b[0],b[1],b[2],c[0],c[1],c[2]);
  }

  template <class NT, class Vector>
  NT determinant_of_vectors(Vector const&a, Vector const&b,
      Vector const&c, Vector const&d){
    return determinant<NT>(
	a[0],a[1],a[2],a[3],
	b[0],b[1],b[2],b[3],
	c[0],c[1],c[2],c[3],
	d[0],d[1],d[2],d[3]);
  }
  template <class NT, class Vector>
  typename Sgn<NT>::result_type
  sign_of_determinant_of_vectors(Vector const&a, Vector const&b,
      Vector const&c, Vector const&d){
    return sign_of_determinant<NT>(
	a[0],a[1],a[2],a[3],
	b[0],b[1],b[2],b[3],
	c[0],c[1],c[2],c[3],
	d[0],d[1],d[2],d[3]);
  }

  template <class NT, class Vector>
  NT determinant_of_vectors(Vector const&a, Vector const&b,
      Vector const&c, Vector const&d, Vector const&e){
    return determinant<NT>(
	a[0],a[1],a[2],a[3],a[4],
	b[0],b[1],b[2],b[3],b[4],
	c[0],c[1],c[2],c[3],c[4],
	d[0],d[1],d[2],d[3],d[4],
	e[0],e[1],e[2],e[3],e[4]);
  }
  template <class NT, class Vector>
  typename Sgn<NT>::result_type
  sign_of_determinant_of_vectors(Vector const&a, Vector const&b,
      Vector const&c, Vector const&d, Vector const&e){
    return sign_of_determinant<NT>(
	a[0],a[1],a[2],a[3],a[4],
	b[0],b[1],b[2],b[3],b[4],
	c[0],c[1],c[2],c[3],c[4],
	d[0],d[1],d[2],d[3],d[4],
	e[0],e[1],e[2],e[3],e[4]);
  }

  template <class NT, class Vector>
  NT determinant_of_vectors(Vector const&a, Vector const&b,
      Vector const&c, Vector const&d, Vector const&e, Vector const&f){
    return determinant<NT>(
	a[0],a[1],a[2],a[3],a[4],a[5],
	b[0],b[1],b[2],b[3],b[4],b[5],
	c[0],c[1],c[2],c[3],c[4],c[5],
	d[0],d[1],d[2],d[3],d[4],d[5],
	e[0],e[1],e[2],e[3],e[4],e[5],
	f[0],f[1],f[2],f[3],f[4],f[5]);
  }
  template <class NT, class Vector>
  typename Sgn<NT>::result_type
  sign_of_determinant_of_vectors(Vector const&a, Vector const&b,
      Vector const&c, Vector const&d, Vector const&e, Vector const&f){
    return sign_of_determinant<NT>(
	a[0],a[1],a[2],a[3],a[4],a[5],
	b[0],b[1],b[2],b[3],b[4],b[5],
	c[0],c[1],c[2],c[3],c[4],c[5],
	d[0],d[1],d[2],d[3],d[4],d[5],
	e[0],e[1],e[2],e[3],e[4],e[5],
	f[0],f[1],f[2],f[3],f[4],f[5]);
  }

}
#endif
