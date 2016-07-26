/******************************************************************************
 *  File: CubicMeshGen.cpp
 *
 *  This file is part of isostuffer
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
#include "CubicMeshGen.h"

REAL CubicMeshGen::dx = 0.;
REAL CubicMeshGen::dy = 0.;
REAL CubicMeshGen::dz = 0.;

//int CubicMeshGen::NX = 85;
//int CubicMeshGen::NY = 2;
//int CubicMeshGen::NZ = 142;
//
//REAL CubicMeshGen::SIZE = 0.0035;

//int CubicMeshGen::NX = 66;
//int CubicMeshGen::NY = 2;
//int CubicMeshGen::NZ = 66;
//
//REAL CubicMeshGen::SIZE = 0.011434/2.;

//int CubicMeshGen::NX = 140;
//int CubicMeshGen::NY = 2;
//int CubicMeshGen::NZ = 188;
//
//REAL CubicMeshGen::SIZE = 1.213/188.;

int CubicMeshGen::NX = 5;
int CubicMeshGen::NY = 240;
int CubicMeshGen::NZ = 1;

REAL CubicMeshGen::SIZE = 0.003;
