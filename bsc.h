/*************************************************************************\

Copyright 2014 Zhejiang University.
All Rights Reserved.

Permission to use, copy, modify and distribute this software and its
documentation for educational, research and non-profit purposes, without
fee, and without a written agreement is hereby granted, provided that the
above copyright notice appear in all
copies.

The authors may be contacted via:

EMail:   tang_m@zju.edu.cn


\**************************************************************************/
/*************************************************************************\

Copyright 2010 The University of North Carolina at Chapel Hill.
All Rights Reserved.

Permission to use, copy, modify and distribute this software and its
documentation for educational, research and non-profit purposes, without
fee, and without a written agreement is hereby granted, provided that the
above copyright notice and the following three paragraphs appear in all
copies.

IN NO EVENT SHALL THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL BE
LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE
USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY
OF NORTH CAROLINA HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH
DAMAGES.

THE UNIVERSITY OF NORTH CAROLINA SPECIFICALLY DISCLAIM ANY
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
NORTH CAROLINA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

The authors may be contacted via:

US Mail:       GAMMA Research Group at UNC
Department of Computer Science
Sitterson Hall, CB #3175
University of N. Carolina
Chapel Hill, NC 27599-3175

Phone:            (919)962-1749

EMail:              geom@cs.unc.edu; tang_m@zju.edu.cn


\**************************************************************************/

namespace rootparity
{
	namespace   // unnamed namespace for local functions
	{
		/**
		 * @brief Detect collisions between a vertex and a triangular face.
		 *
		 * Looks for collisions between the vertex a0 and the face (b0, c0, d0)
		 * as they move towards a1 and (b1, c1, d1).
		 *
		 * @param  a0  Start position of the vertex.
		 * @param  b0  Start position of the first vertex of the face.
		 * @param  c0  Start position of the second vertex of the face.
		 * @param  d0  Start position of the third vertex of the face.
		 * @param  a1  End position of the vertex.
		 * @param  b1  End position of the first vertex of the face.
		 * @param  c1  End position of the second vertex of the face.
		 * @param  d1  End position of the third vertex of the face.
		 *
		 * @returns  True if the vertex and face collide.
		 */
		bool
		Intersect_VF_robust(
				const Vec3d &a0, const Vec3d &b0, const Vec3d &c0, const Vec3d &d0,
				const Vec3d &a1, const Vec3d &b1, const Vec3d &c1, const Vec3d &d1);

		/**
		 * @brief Detect collisions between two edges as they move.
		 *
		 * Looks for collisions between edges (a0, b0) and (c0, d0) as they
		 * move towards (a1, b1) and (c1, d1).
		 *
		 * @param  a0  Start position of the first edge's first vertex.
		 * @param  b0  Start position of the first edge's second vertex.
		 * @param  c0  Start position of the second edge's first vertex.
		 * @param  d0  Start position of the second edge's second vertex.
		 * @param  a1  End position of the first edge's first vertex.
		 * @param  b1  End position of the first edge's second vertex.
		 * @param  c1  End position of the second edge's first vertex.
		 * @param  d1  End position of the second edge's second vertex.
		 *
		 * @returns True if the edges collide.
		 */
		bool
		Intersect_EE_robust(
				const Vec3d &a0, const Vec3d &b0, const Vec3d &c0, const Vec3d &d0,
				const Vec3d &a1, const Vec3d &b1, const Vec3d &c1, const Vec3d &d1);
	}
}  // namespace rootparity
