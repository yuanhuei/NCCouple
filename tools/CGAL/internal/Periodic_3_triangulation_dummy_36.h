// Copyright (c) 2009   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: https://github.com/CGAL/cgal/blob/releases/CGAL-4.14.3/Periodic_3_triangulation_3/include/CGAL/internal/Periodic_3_triangulation_dummy_36.h $
// $Id: Periodic_3_triangulation_dummy_36.h 99617c3 2017-12-07T16:58:47+00:00 Andreas Fabri
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Manuel Caroli <Manuel.Caroli@sophia.inria.fr>

#ifdef CGAL_INCLUDE_FROM_PERIODIC_3_TRIANGULATION_3_H

template < class GT, class TDS >
inline std::vector<typename Periodic_3_triangulation_3<GT,TDS>::Vertex_handle >
Periodic_3_triangulation_3<GT,TDS>::insert_dummy_points() {

    
static const unsigned V[216][4] = {
  { 3, 9, 4, 12 },
  { 28, 8, 3, 6 },
  { 0, 11, 9, 3 },
  { 11, 19, 13, 10 },
  { 7, 15, 12, 13 },
  { 6, 12, 7, 15 },
  { 2, 8, 6, 17 },
  { 3, 12, 4, 6 },
  { 25, 24, 16, 22 },
  { 3, 14, 12, 6 },
  { 2, 1, 10, 4 },
  { 8, 7, 5, 13 },
  { 27, 5, 28, 8 },
  { 29, 21, 27, 18 },
  { 26, 25, 17, 23 },
  { 21, 30, 27, 22 },
  { 29, 6, 3, 4 },
  { 17, 8, 6, 14 },
  { 22, 30, 31, 25 },
  { 6, 2, 17, 0 },
  { 11, 19, 14, 13 },
  { 2, 4, 10, 5 },
  { 12, 21, 23, 18 },
  { 16, 14, 17, 22 },
  { 31, 2, 34, 33 },
  { 25, 31, 23, 26 },
  { 25, 22, 17, 23 },
  { 14, 5, 11, 13 },
  { 2, 16, 17, 10 },
  { 9, 24, 16, 10 },
  { 19, 18, 13, 10 },
  { 17, 2, 11, 0 },
  { 21, 32, 27, 30 },
  { 4, 7, 13, 5 },
  { 11, 2, 10, 5 },
  { 0, 2, 11, 5 },
  { 6, 15, 7, 0 },
  { 8, 5, 14, 13 },
  { 30, 8, 31, 2 },
  { 31, 8, 6, 2 },
  { 19, 28, 20, 22 },
  { 13, 15, 21, 16 },
  { 10, 18, 25, 19 },
  { 9, 26, 15, 24 },
  { 12, 23, 20, 18 },
  { 30, 8, 28, 31 },
  { 16, 8, 14, 13 },
  { 15, 26, 17, 23 },
  { 4, 12, 9, 10 },
  { 6, 17, 15, 0 },
  { 10, 1, 9, 4 },
  { 14, 19, 20, 22 },
  { 23, 14, 20, 22 },
  { 16, 1, 9, 10 },
  { 22, 28, 23, 31 },
  { 4, 3, 1, 9 },
  { 0, 3, 9, 1 },
  { 17, 11, 15, 0 },
  { 19, 28, 34, 20 },
  { 9, 12, 18, 10 },
  { 7, 6, 4, 12 },
  { 33, 4, 2, 5 },
  { 33, 35, 1, 4 },
  { 0, 5, 11, 3 },
  { 4, 7, 12, 13 },
  { 7, 16, 13, 8 },
  { 32, 34, 26, 31 },
  { 26, 23, 15, 21 },
  { 7, 1, 16, 8 },
  { 3, 8, 14, 6 },
  { 4, 13, 10, 5 },
  { 16, 1, 10, 2 },
  { 13, 5, 11, 10 },
  { 13, 21, 18, 19 },
  { 9, 20, 18, 12 },
  { 32, 0, 6, 7 },
  { 3, 14, 9, 12 },
  { 7, 15, 1, 0 },
  { 27, 29, 4, 7 },
  { 28, 8, 6, 31 },
  { 15, 9, 1, 0 },
  { 7, 15, 13, 16 },
  { 4, 13, 12, 10 },
  { 15, 11, 9, 0 },
  { 2, 17, 11, 10 },
  { 7, 1, 15, 16 },
  { 8, 1, 16, 2 },
  { 15, 1, 9, 16 },
  { 3, 11, 9, 14 },
  { 14, 5, 3, 11 },
  { 8, 5, 3, 14 },
  { 19, 34, 28, 27 },
  { 23, 14, 22, 17 },
  { 11, 19, 26, 20 },
  { 18, 29, 35, 33 },
  { 12, 13, 18, 10 },
  { 23, 21, 29, 18 },
  { 17, 14, 6, 12 },
  { 15, 17, 6, 12 },
  { 16, 8, 17, 14 },
  { 2, 8, 17, 16 },
  { 24, 26, 15, 21 },
  { 35, 3, 1, 4 },
  { 12, 14, 23, 17 },
  { 24, 30, 32, 21 },
  { 9, 11, 26, 20 },
  { 29, 6, 4, 7 },
  { 33, 30, 24, 25 },
  { 21, 24, 30, 22 },
  { 18, 29, 20, 35 },
  { 35, 34, 3, 28 },
  { 24, 21, 15, 16 },
  { 20, 28, 23, 22 },
  { 25, 33, 30, 31 },
  { 29, 23, 28, 31 },
  { 12, 17, 23, 15 },
  { 9, 11, 15, 26 },
  { 12, 14, 20, 23 },
  { 29, 28, 3, 6 },
  { 24, 35, 30, 33 },
  { 9, 20, 26, 24 },
  { 15, 11, 17, 26 },
  { 9, 20, 24, 18 },
  { 27, 7, 5, 8 },
  { 19, 27, 28, 22 },
  { 9, 18, 24, 10 },
  { 35, 20, 34, 28 },
  { 12, 23, 21, 15 },
  { 22, 30, 28, 31 },
  { 29, 32, 6, 7 },
  { 25, 16, 17, 22 },
  { 24, 21, 16, 22 },
  { 20, 11, 19, 14 },
  { 9, 11, 20, 14 },
  { 9, 14, 20, 12 },
  { 10, 25, 11, 19 },
  { 35, 32, 26, 24 },
  { 12, 21, 18, 13 },
  { 12, 15, 21, 13 },
  { 22, 31, 23, 25 },
  { 9, 24, 15, 16 },
  { 35, 32, 30, 1 },
  { 11, 25, 17, 26 },
  { 32, 35, 0, 1 },
  { 13, 22, 21, 19 },
  { 22, 13, 14, 19 },
  { 16, 13, 14, 22 },
  { 13, 16, 21, 22 },
  { 34, 2, 0, 5 },
  { 28, 5, 3, 8 },
  { 24, 30, 22, 25 },
  { 21, 23, 29, 32 },
  { 33, 35, 30, 1 },
  { 19, 33, 34, 27 },
  { 27, 8, 28, 30 },
  { 27, 30, 7, 8 },
  { 19, 21, 27, 22 },
  { 18, 33, 25, 19 },
  { 10, 18, 24, 25 },
  { 10, 24, 16, 25 },
  { 10, 16, 17, 25 },
  { 10, 17, 11, 25 },
  { 27, 30, 28, 22 },
  { 18, 20, 29, 23 },
  { 32, 23, 31, 26 },
  { 11, 19, 25, 26 },
  { 33, 2, 30, 31 },
  { 30, 1, 7, 8 },
  { 27, 7, 4, 5 },
  { 18, 21, 27, 19 },
  { 24, 32, 30, 35 },
  { 30, 1, 2, 33 },
  { 33, 31, 25, 34 },
  { 34, 5, 28, 27 },
  { 32, 31, 6, 0 },
  { 29, 23, 31, 32 },
  { 18, 35, 24, 33 },
  { 35, 34, 26, 32 },
  { 29, 3, 35, 4 },
  { 29, 35, 3, 28 },
  { 33, 1, 2, 4 },
  { 34, 5, 0, 3 },
  { 35, 3, 0, 1 },
  { 19, 34, 25, 26 },
  { 30, 1, 8, 2 },
  { 33, 2, 34, 5 },
  { 28, 5, 34, 3 },
  { 19, 34, 26, 20 },
  { 31, 2, 0, 34 },
  { 18, 33, 24, 25 },
  { 30, 32, 27, 7 },
  { 30, 32, 7, 1 },
  { 21, 29, 27, 32 },
  { 32, 0, 7, 1 },
  { 32, 29, 27, 7 },
  { 29, 31, 28, 6 },
  { 31, 2, 6, 0 },
  { 29, 20, 28, 23 },
  { 20, 34, 26, 35 },
  { 32, 0, 34, 31 },
  { 18, 35, 20, 24 },
  { 20, 35, 26, 24 },
  { 18, 27, 33, 19 },
  { 33, 29, 35, 4 },
  { 33, 5, 34, 27 },
  { 19, 33, 25, 34 },
  { 34, 31, 25, 26 },
  { 24, 32, 26, 21 },
  { 21, 23, 32, 26 },
  { 29, 31, 6, 32 },
  { 33, 27, 29, 4 },
  { 27, 4, 33, 5 },
  { 18, 27, 29, 33 },
  { 29, 20, 35, 28 },
  { 34, 3, 0, 35 },
  { 35, 0, 34, 32 }
};

  
static const unsigned N[216][4] = {
   { 48, 7, 76, 55 },
   { 69, 118, 79, 149 },
   { 88, 56, 63, 83 },
   { 30, 72, 135, 20 },
   { 138, 64, 81, 5 },
   { 4, 36, 98, 60 },
   { 17, 19, 100, 39 },
   { 60, 16, 9, 0 },
   { 131, 130, 150, 159 },
   { 97, 7, 69, 76 },
   { 50, 21, 180, 71 },
   { 33, 37, 65, 123 },
   { 149, 154, 123, 173 },
   { 169, 212, 96, 192 },
   { 26, 47, 25, 142 },
   { 162, 156, 108, 32 },
   { 7, 178, 106, 118 },
   { 69, 97, 99, 6 },
   { 113, 139, 150, 128 },
   { 31, 49, 196, 6 },
   { 145, 27, 3, 132 },
   { 70, 34, 61, 10 },
   { 96, 44, 137, 127 },
   { 92, 130, 146, 99 },
   { 185, 172, 166, 188 },
   { 164, 14, 206, 139 },
   { 92, 14, 139, 130 },
   { 72, 20, 37, 89 },
   { 160, 84, 71, 100 },
   { 159, 53, 125, 140 },
   { 95, 3, 42, 73 },
   { 35, 57, 19, 84 },
   { 190, 15, 104, 192 },
   { 11, 70, 168, 64 },
   { 21, 72, 35, 84 },
   { 34, 63, 148, 31 },
   { 77, 75, 49, 5 },
   { 27, 46, 11, 90 },
   { 39, 166, 184, 45 },
   { 6, 196, 38, 79 },
   { 112, 51, 124, 58 },
   { 111, 147, 81, 138 },
   { 157, 135, 30, 158 },
   { 101, 140, 120, 116 },
   { 163, 74, 22, 117 },
   { 79, 128, 38, 154 },
   { 37, 146, 65, 99 },
   { 14, 115, 67, 121 },
   { 59, 50, 82, 0 },
   { 57, 36, 19, 98 },
   { 55, 48, 10, 53 },
   { 40, 52, 145, 132 },
   { 51, 112, 92, 117 },
   { 50, 29, 71, 87 },
   { 114, 139, 128, 112 },
   { 56, 50, 0, 102 },
   { 55, 80, 182, 2 },
   { 83, 49, 31, 121 },
   { 126, 187, 40, 91 },
   { 95, 125, 48, 74 },
   { 7, 64, 5, 106 },
   { 21, 185, 211, 180 },
   { 102, 180, 203, 152 },
   { 89, 2, 181, 35 },
   { 4, 82, 33, 60 },
   { 46, 11, 68, 81 },
   { 206, 164, 199, 177 },
   { 127, 101, 208, 47 },
   { 86, 65, 167, 85 },
   { 17, 9, 1, 90 },
   { 72, 21, 33, 82 },
   { 10, 28, 86, 53 },
   { 34, 3, 70, 27 },
   { 169, 30, 144, 137 },
   { 44, 59, 134, 122 },
   { 36, 129, 193, 174 },
   { 134, 0, 9, 88 },
   { 80, 193, 36, 85 },
   { 106, 168, 194, 210 },
   { 39, 195, 45, 1 },
   { 56, 77, 83, 87 },
   { 41, 65, 85, 4 },
   { 95, 48, 70, 64 },
   { 2, 80, 57, 116 },
   { 161, 34, 28, 31 },
   { 87, 81, 68, 77 },
   { 71, 100, 184, 68 },
   { 53, 140, 85, 80 },
   { 133, 76, 89, 2 },
   { 63, 88, 27, 90 },
   { 89, 69, 37, 149 },
   { 173, 124, 153, 58 },
   { 23, 26, 103, 52 },
   { 187, 105, 132, 165 },
   { 203, 176, 212, 109 },
   { 30, 59, 82, 137 },
   { 13, 163, 22, 151 },
   { 9, 98, 103, 17 },
   { 97, 5, 115, 49 },
   { 17, 23, 46, 100 },
   { 99, 28, 86, 6 },
   { 67, 111, 207, 43 },
   { 55, 62, 178, 182 },
   { 92, 115, 97, 117 },
   { 32, 207, 108, 170 },
   { 93, 120, 133, 116 },
   { 60, 78, 129, 16 },
   { 150, 189, 113, 119 },
   { 150, 15, 131, 104 },
   { 213, 200, 94, 163 },
   { 186, 179, 126, 214 },
   { 41, 140, 131, 101 },
   { 54, 52, 40, 197 },
   { 166, 18, 172, 107 },
   { 54, 195, 175, 197 },
   { 47, 127, 98, 103 },
   { 121, 43, 105, 83 },
   { 52, 44, 103, 134 },
   { 1, 16, 195, 179 },
   { 152, 107, 176, 170 },
   { 201, 43, 122, 105 },
   { 142, 47, 116, 57 },
   { 200, 125, 74, 120 },
   { 11, 12, 155, 168 },
   { 162, 40, 156, 91 },
   { 158, 29, 59, 122 },
   { 58, 110, 213, 198 },
   { 67, 138, 115, 22 },
   { 45, 54, 18, 162 },
   { 75, 106, 194, 209 },
   { 23, 26, 8, 160 },
   { 147, 8, 108, 111 },
   { 20, 51, 133, 93 },
   { 132, 134, 88, 105 },
   { 117, 74, 76, 133 },
   { 165, 3, 42, 161 },
   { 207, 201, 170, 177 },
   { 73, 95, 138, 22 },
   { 41, 137, 4, 127 },
   { 25, 26, 18, 54 },
   { 111, 87, 29, 43 },
   { 191, 152, 143, 170 },
   { 14, 121, 165, 161 },
   { 182, 193, 141, 215 },
   { 156, 73, 145, 147 },
   { 20, 51, 144, 146 },
   { 145, 23, 147, 46 },
   { 131, 144, 146, 41 },
   { 35, 181, 185, 188 },
   { 90, 1, 12, 186 },
   { 18, 8, 107, 108 },
   { 175, 192, 208, 96 },
   { 141, 171, 62, 119 },
   { 204, 91, 202, 205 },
   { 45, 162, 155, 12 },
   { 167, 123, 154, 190 },
   { 15, 124, 144, 169 },
   { 205, 42, 202, 189 },
   { 189, 159, 42, 125 },
   { 8, 160, 158, 29 },
   { 130, 161, 159, 28 },
   { 142, 135, 160, 84 },
   { 128, 124, 15, 154 },
   { 197, 96, 44, 109 },
   { 25, 66, 208, 175 },
   { 183, 142, 93, 135 },
   { 38, 113, 24, 171 },
   { 68, 155, 184, 191 },
   { 33, 211, 123, 78 },
   { 156, 202, 73, 13 },
   { 141, 119, 136, 104 },
   { 180, 166, 152, 184 },
   { 206, 205, 24, 113 },
   { 12, 91, 204, 186 },
   { 196, 75, 199, 209 },
   { 164, 209, 151, 114 },
   { 119, 189, 94, 200 },
   { 66, 136, 215, 198 },
   { 102, 203, 16, 179 },
   { 110, 118, 213, 178 },
   { 10, 61, 62, 171 },
   { 63, 214, 186, 148 },
   { 56, 143, 102, 214 },
   { 206, 165, 187, 205 },
   { 86, 38, 171, 167 },
   { 148, 204, 61, 24 },
   { 181, 110, 149, 173 },
   { 198, 93, 58, 183 },
   { 148, 199, 24, 196 },
   { 107, 158, 157, 176 },
   { 194, 155, 191, 32 },
   { 193, 167, 141, 190 },
   { 194, 32, 151, 13 },
   { 77, 191, 143, 75 },
   { 78, 190, 129, 192 },
   { 79, 118, 209, 114 },
   { 19, 174, 188, 39 },
   { 112, 114, 163, 213 },
   { 177, 201, 126, 187 },
   { 188, 66, 174, 215 },
   { 201, 122, 176, 109 },
   { 136, 120, 200, 198 },
   { 153, 157, 169, 212 },
   { 178, 62, 210, 94 },
   { 173, 153, 211, 185 },
   { 172, 183, 153, 157 },
   { 25, 183, 66, 172 },
   { 208, 101, 104, 136 },
   { 164, 207, 67, 151 },
   { 174, 129, 175, 195 },
   { 78, 203, 211, 212 },
   { 61, 204, 168, 210 },
   { 210, 94, 202, 13 },
   { 126, 179, 197, 109 },
   { 182, 215, 110, 181 },
   { 199, 177, 143, 214 }
};

static const int O[216][4] = {
  { 0,0,0,0 },
  { 0,1,5,5 },
  { 4,0,4,4 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 2,0,4,0 },
  { 0,0,0,0 },
  { 2,2,0,0 },
  { 4,0,4,4 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 0,1,0,1 },
  { 0,0,0,0 },
  { 6,2,0,4 },
  { 0,0,0,0 },
  { 0,1,1,1 },
  { 0,0,4,0 },
  { 0,0,0,2 },
  { 4,2,0,6 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 0,3,2,2 },
  { 2,0,4,6 },
  { 2,0,0,4 },
  { 0,0,0,0 },
  { 2,0,0,2 },
  { 2,2,0,2 },
  { 0,0,0,0 },
  { 0,2,2,6 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 4,0,0,0 },
  { 0,0,0,2 },
  { 0,0,0,0 },
  { 0,1,0,3 },
  { 0,1,5,3 },
  { 0,0,4,0 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 2,2,0,2 },
  { 0,0,0,0 },
  { 0,1,0,0 },
  { 0,0,0,0 },
  { 4,6,0,4 },
  { 0,0,0,0 },
  { 4,0,4,6 },
  { 0,0,0,0 },
  { 0,0,4,0 },
  { 4,0,4,0 },
  { 0,2,2,2 },
  { 0,0,4,0 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 0,2,4,6 },
  { 0,0,0,4 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 0,1,1,1 },
  { 0,0,1,1 },
  { 4,0,0,4 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 4,2,6,0 },
  { 2,0,0,0 },
  { 0,2,0,0 },
  { 4,0,0,4 },
  { 0,0,0,0 },
  { 0,2,2,2 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 0,3,1,1 },
  { 4,0,4,4 },
  { 0,0,2,2 },
  { 0,0,1,1 },
  { 0,1,5,0 },
  { 0,2,2,2 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 4,2,6,6 },
  { 2,0,2,2 },
  { 0,2,0,0 },
  { 0,2,0,2 },
  { 0,2,2,0 },
  { 4,0,4,0 },
  { 0,0,4,0 },
  { 0,0,4,0 },
  { 0,0,0,0 },
  { 4,0,0,0 },
  { 0,0,4,4 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 0,0,4,4 },
  { 4,0,4,4 },
  { 0,0,0,0 },
  { 2,0,0,0 },
  { 2,2,0,0 },
  { 0,1,1,1 },
  { 4,0,4,0 },
  { 2,0,0,0 },
  { 4,0,4,4 },
  { 0,1,1,1 },
  { 2,0,2,2 },
  { 0,2,0,0 },
  { 0,0,0,0 },
  { 4,0,5,0 },
  { 2,0,0,0 },
  { 4,0,4,0 },
  { 2,2,0,0 },
  { 4,4,0,0 },
  { 4,0,4,4 },
  { 6,2,4,6 },
  { 4,0,4,4 },
  { 4,0,5,5 },
  { 2,2,0,2 },
  { 0,0,0,0 },
  { 4,2,0,6 },
  { 0,0,0,0 },
  { 0,1,1,1 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 4,4,0,0 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 0,0,1,1 },
  { 2,0,0,0 },
  { 2,0,0,0 },
  { 4,0,0,0 },
  { 4,0,4,0 },
  { 4,0,4,4 },
  { 0,0,0,0 },
  { 2,0,2,2 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 0,0,4,2 },
  { 2,2,0,0 },
  { 2,0,0,3 },
  { 2,2,0,6 },
  { 0,2,3,3 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 0,1,5,1 },
  { 0,1,5,1 },
  { 2,0,0,2 },
  { 0,0,0,0 },
  { 2,2,0,3 },
  { 0,0,0,0 },
  { 0,1,0,0 },
  { 0,0,1,1 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 2,2,0,2 },
  { 2,0,0,2 },
  { 2,0,2,2 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 4,4,0,6 },
  { 0,0,0,4 },
  { 2,3,0,0 },
  { 0,3,1,1 },
  { 0,1,1,1 },
  { 0,0,0,0 },
  { 2,0,0,2 },
  { 0,3,3,2 },
  { 2,0,2,2 },
  { 0,1,0,0 },
  { 4,0,5,7 },
  { 4,4,0,4 },
  { 0,0,0,0 },
  { 6,2,6,4 },
  { 0,1,0,1 },
  { 4,4,5,0 },
  { 0,1,1,1 },
  { 0,1,5,5 },
  { 0,1,1,1 },
  { 0,0,0,4 },
  { 0,3,1,3 },
  { 0,1,0,1 },
  { 0,1,0,5 },
  { 0,0,4,4 },
  { 0,3,7,2 },
  { 0,0,0,0 },
  { 0,0,0,1 },
  { 0,0,1,3 },
  { 0,0,0,0 },
  { 0,3,1,3 },
  { 0,0,0,1 },
  { 4,0,0,5 },
  { 0,3,5,7 },
  { 4,4,0,4 },
  { 4,0,4,4 },
  { 4,7,2,0 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 0,0,0,0 },
  { 0,0,0,1 },
  { 0,1,0,0 },
  { 0,0,0,0 },
  { 2,0,2,6 },
  { 2,0,2,0 },
  { 0,0,0,2 },
  { 4,0,5,4 },
  { 0,0,0,1 },
  { 0,1,0,1 },
  { 0,0,0,0 },
  { 4,4,4,0 },
  { 0,5,5,4 },
  { 6,7,2,4 }
};
 
  clear();

  Vertex_handle vertices[36];
  Cell_handle cells[216];

  // Initialise vertices:
  for (int i=0; i<4; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
        // Initialise virtual vertices out of the domain for debugging
        vertices[9*i+3*j+k] = _tds.create_vertex();
        Point p(k*(1.0/3.0) + i*(1.0/6.0),
                j*(1.0/3.0) + i*(1.0/6.0), i*(1.0/4.0) );
        p = Point((p.x() > FT(0.9375) ? (std::max)( p.x()-1, FT(0) ) : p.x()),
                  (p.y() > FT(0.9375) ? (std::max)( p.y()-1, FT(0) ) : p.y()), p.z());
        p = Point((domain().xmax()-domain().xmin())*p.x(),
                  (domain().xmax()-domain().xmin())*p.y(),
                  (domain().xmax()-domain().xmin())*p.z());
        p = Point(p.x() + domain().xmin(),
                  p.y() + domain().ymin(),
                  p.z() + domain().zmin());
        vertices[9*i+3*j+k]->set_point(p);
      }
    }
  }

  // Create cells:
  for (int i=0; i<216; i++) {
    cells[i] = _tds.create_cell();
    cells[i]->set_vertices(vertices[V[i][0]], vertices[V[i][1]], vertices[V[i][2]], vertices[V[i][3]]);
  }
  for (int i=0; i<216; i++) {
    cells[i]->set_neighbors(cells[N[i][0]], cells[N[i][1]], cells[N[i][2]], cells[N[i][3]]);
    set_offsets(cells[i], O[i][0], O[i][1], O[i][2], O[i][3]);
  }



  vertices[0]->set_cell(cells[2]);
  vertices[1]->set_cell(cells[10]);
  vertices[2]->set_cell(cells[6]);
  vertices[3]->set_cell(cells[0]);
  vertices[4]->set_cell(cells[0]);
  vertices[5]->set_cell(cells[11]);
  vertices[6]->set_cell(cells[1]);
  vertices[7]->set_cell(cells[4]);
  vertices[8]->set_cell(cells[1]);
  vertices[9]->set_cell(cells[0]);
  vertices[10]->set_cell(cells[3]);
  vertices[11]->set_cell(cells[2]);
  vertices[12]->set_cell(cells[0]);
  vertices[13]->set_cell(cells[3]);
  vertices[14]->set_cell(cells[9]);
  vertices[15]->set_cell(cells[4]);
  vertices[16]->set_cell(cells[8]);
  vertices[17]->set_cell(cells[6]);
  vertices[18]->set_cell(cells[13]);
  vertices[19]->set_cell(cells[3]);
  vertices[20]->set_cell(cells[40]);
  vertices[21]->set_cell(cells[13]);
  vertices[22]->set_cell(cells[8]);
  vertices[23]->set_cell(cells[14]);
  vertices[24]->set_cell(cells[8]);
  vertices[25]->set_cell(cells[8]);
  vertices[26]->set_cell(cells[14]);
  vertices[27]->set_cell(cells[12]);
  vertices[28]->set_cell(cells[1]);
  vertices[29]->set_cell(cells[13]);
  vertices[30]->set_cell(cells[15]);
  vertices[31]->set_cell(cells[18]);
  vertices[32]->set_cell(cells[32]);
  vertices[33]->set_cell(cells[24]);
  vertices[34]->set_cell(cells[24]);
  vertices[35]->set_cell(cells[62]);

  _tds.set_dimension(3);
  _cover = make_array(1,1,1);

  std::vector<Vertex_handle> ret_vector(36);
  for (int i=0; i<36; i++) {
        ret_vector[i] = vertices[i];
  }

  return ret_vector;
}

#endif // CGAL_INCLUDE_FROM_PERIODIC_3_TRIANGULATION_3_H
