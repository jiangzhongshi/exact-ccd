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

#pragma warning(disable: 4996)
#include "bsc.h"
#include "rootparitycollisiontest.h"
#include <cstdlib>
#include <cfloat>

namespace bsc
{
	using namespace rootparity;

	template<unsigned int N, class T>
	inline void make_vector( const Vec<N,double>& in, Vec<N,T>& out )
	{
		for(int i = 0; i < N; i++){
			out[i] = T(in[i]);
		}
	}

	//					| a  b |
	//  return		| c  d | = a*d - b*c
	template <class T>
	inline T det2x2(const T &a, const T &b, const T &c, const T &d)
	{
		return a*d-b*c;
	}

	// from polynomial decomposition theorem
	template <class T>
	inline bool bezierDecomposition(const T &k0, const T &k1, const T &k2, const T &k3,
						const T &j0, const T &j1, const T &j2,
						T &m0, T &m1,
						T &n0, T &n1)
	{
		T A = (j1-j2)*T(2.0);
		T B = j0-j2;
		T C = k2*T(3.0) - k3*T(2.0) - k0;
		T D = k1*T(3.0) - k0*T(2.0) -k3;
		T E = j2-j0;
		T F = (j1-j0)*T(2.0);

		T tt = det2x2(A, B, E, F);
		 if ( tt.is_certainly_zero() ) {
			printf("@@@det = 0.\n");
			return false;
		}

	/*
		m0 = det2x2(A, B, C, D)/tt;
		m1 = det2x2(F, E, D, C)/tt;
		n0 = k0-m0*j0;
		n1 = k3-m1*j2;
	*/
		m0 = det2x2(A, B, C, D);
		m1 = det2x2(F, E, D, C);
		n0 = k0*tt-m0*j0;
		n1 = k3*tt-m1*j2;

		return true;
	}

	template <class T>
	inline T _evaluateBezier2(const T &p0, const T &p1, const T &p2, const T &t, const T &s)
	{
		T s2 = s*s;
		T t2 = t*t;

		return p0*s2+p1*T(2.0)*s*t+p2*t2;
	}

	template <class T>
	inline T evaluateBezier1(const T &p0, const T &p1, const T &t)
	{
		T s = T(1.0)-t;
		return p0*s + p1*t;
	}

	template <class T>
	inline T evaluateBezier2(const T &p0, const T &p1, const T &p2, const T &t)
	{
		T s = T(1.0)-t;
		return _evaluateBezier2(p0, p1, p2, t, s);
	}

	template <class T>
	inline T _evaluateBezier(const T &p0, const T &p1, const T &p2, const T &p3, const T &t, const T &s)
	{
		T s2 = s*s;
		T s3 = s2*s;
		T t2 = t*t;
		T t3 = t2*t;

		return p0*s3+p1*T(3.0)*s2*t+p2*T(3.0)*s*t2+p3*t3;
	}

	template <class T>
	inline T evaluateBezier(const T &p0, const T &p1, const T &p2, const T &p3, const T &t)
	{
		T s = T(1.0)-t;
		return _evaluateBezier(p0, p1, p2, p3, t, s);
	}

	template <class T>
	inline T _evaluateBezier4(const T &p0, const T &p1, const T &p2, const T &p3, const T &p4, const T &t, const T &s)
	{
		T s2 = s*s;
		T s3 = s2*s;
		T s4 = s2*s2;
		T t2 = t*t;
		T t3 = t2*t;
		T t4 = t2*t2;

		return p0*s4+p1*4*s3*t+p2*6*s2*t2+p3*4*s*t3+p4*t4;
	}

	template <class T>
	inline T evaluateBezier4(const T &p0, const T &p1, const T &p2, const T &p3, const T &p4, const T &t)
	{
		T s = T(1.0)-t;
		return _evaluateBezier4(p0, p1, p2, p3, p4, t, s);
	}

	template <class T>
	inline Vec<3, T> norm(const Vec<3, T> &p1, const Vec<3, T> &p2, const Vec<3, T> &p3)
	{
		return cross(p2-p1, p3-p1);
	}

	template <class T>
	inline void getBezier4(
			const Vec<3, T> &a0, const Vec<3, T> &b0, const Vec<3, T> &c0, const Vec<3, T> &d0,
			const Vec<3, T> &a1, const Vec<3, T> &b1, const Vec<3, T> &c1, const Vec<3, T> &d1,
			T &l0, T &l1, T &l2, T &l3, T &l4, int which, bool ee_test)
	{
		Vec<3, T> n0 = norm(a0, b0, c0);
		Vec<3, T> n1 = norm(a1, b1, c1);
		Vec<3, T> deltaN = norm(a1-a0, b1-b0, c1-c0);
		Vec<3, T> nX = (n0+n1-deltaN)*T(0.5);

		Vec<3, T> m0, m1, deltaM, mX;

		if (which == 0) { // (bt-pt) x (ct-pt) . nt
			m0 = norm(d0, b0, c0);
			m1 = norm(d1, b1, c1);
			deltaM = norm(d1-d0, b1-b0, c1-c0);
		} else
		if (which == 1) { // ct-pt x at-pt . nt
			m0 = norm(d0, c0, a0);
			m1 = norm(d1, c1, a1);
			deltaM = norm(d1-d0, c1-c0, a1-a0);
		} else
		if (which == 2) {// at-pt x bt-pt .nt
			m0 = norm(d0, a0, b0);
			m1 = norm(d1, a1,b1);
			deltaM = norm(d1-d0, a1-a0, b1-b0);
		} else
			printf("@@@Imposible be here!");

		mX = (m0+m1-deltaM)*T(0.5);

		l0 = dot(m0, n0)*T(6.0);
		l1 = (dot(m0, nX) + dot(mX, n0))*T(3.0);
		l2 = dot(m0, n1) + dot(mX, nX)*T(4.0) + dot(m1, n0);
		l3 = (dot(mX, n1) + dot(m1, nX))*T(3.0);
		l4 = dot(m1, n1)*T(6.0);

		if (!ee_test && which != 0) {
			l0 = -l0, l1 = -l1, l2 = -l2, l3 = -l3, l4 = -l4;
		}

		if (which ==2 && ee_test) {
			l0 = -l0, l1 = -l1, l2 = -l2, l3 = -l3, l4 = -l4;
		}
	}

	inline bool
	DNF_Culling(
			const Vec3d &a0, const Vec3d &b0, const Vec3d &c0, const Vec3d &d0,
			const Vec3d &a1, const Vec3d &b1, const Vec3d &c1, const Vec3d &d1)
	{
		Vec3d n0 = norm(a0, b0, c0);
		Vec3d n1 = norm(a1, b1, c1);
		Vec3d delta = norm(a1-a0, b1-b0, c1-c0);
		Vec3d nX = (n0+n1-delta)*0.5;

		Vec3d pa0 = d0-a0;
		Vec3d pa1 = d1-a1;

		double A = dot(n0, pa0);
		double B = dot(n1, pa1);
		double C = dot(nX, pa0);
		double D = dot(nX, pa1);
		double E = dot(n1, pa0);
		double F = dot(n0, pa1);

		double p0 = A;
		double p1 = C*2.0+F;
		double p2 = D*2.0+E;
		double p3 = B;

		double e1, e2; // conservative bounds, can be calculated on the fly

		e1 = DBL_EPSILON*100;
		e2 = DBL_EPSILON*100;

		if (p0 > e1 && p1 > e2 && p2 > e2 && p3 > e1)
			return false;

		if (p0 < -e1 && p1 < -e2 && p2 < -e2 && p3 < -e1)
			return false;

		return true;
	}

	template <class T>
	inline bool
	getBezier(
			const Vec<3, T> &a0, const Vec<3, T> &b0, const Vec<3, T> &c0, const Vec<3, T> &d0,
			const Vec<3, T> &a1, const Vec<3, T> &b1, const Vec<3, T> &c1, const Vec<3, T> &d1,
			T &p0, T &p1, T &p2, T &p3)
	{
		Vec<3, T> n0 = norm(a0, b0, c0);
		Vec<3, T> n1 = norm(a1, b1, c1);
		Vec<3, T> delta = norm(a1-a0, b1-b0, c1-c0);
		Vec<3, T> nX = (n0+n1-delta)*T(0.5);

		Vec<3, T> pa0 = d0-a0;
		Vec<3, T> pa1 = d1-a1;

		T A = dot(n0, pa0);
		T B = dot(n1, pa1);
		T C = dot(nX, pa0);
		T D = dot(nX, pa1);
		T E = dot(n1, pa0);
		T F = dot(n0, pa1);

		p0 = A*T(3.0);
		p1 = C*T(2.0)+F;
		p2 = D*T(2.0)+E;
		p3 = B*T(3.0);

		//if (p0 > 0 && p1 > 0 && p2 > 0 && p3 > 0)
		//if (sign(p0)>0 && sign(p1)>0 && sign(p2)>0 && sign(p3)>0)
		if (p0.is_certainly_positive() && p1.is_certainly_positive() &&
			p2.is_certainly_positive() && p3.is_certainly_positive())
			return false;

		//if (p0 < 0 && p1 < 0 && p2 < 0 && p3 < 0)
		//if (sign(p0)<0 && sign(p1)<0 && sign(p2)<0 && sign(p3)<0)
		if (p0.is_certainly_negative() && p1.is_certainly_negative() &&
			p2.is_certainly_negative() && p3.is_certainly_negative())
			return false;

		return true;
	}


	template <class T>
	class bcrv {
	public:
		T k0, k1, k2, k3;
		T kk0, kk1, kk2;
		int ct;

	public:
		bcrv(const T& k0, const T& k1, const T& k2, const T& k3,
			const T& kk0, const T& kk1, const T& kk2, int ct)
		{
			this->k0 = k0;
			this->k1 = k1;
			this->k2 = k2;
			this->k3 = k3;
			this->kk0 = kk0;
			this->kk1 = kk1;
			this->kk2 = kk2;
			this->ct = ct;
		}
	};

	template <class T>
	inline bool diffSign(const T& a, const T& b)
	{
		//return (sign(a) < 0 && sign(b) > 0) || (sign(a) > 0 && sign(b) < 0);
		return !a.same_sign(b);
	}

	template <class T>
	inline bool sameSign(const T& a, const T& b)
	{
		//return !diffSign(a, b);
		return a.same_sign(b);
	}

	template <class T>
	inline T lineRoot(T a, T b)
	{
		return a/(a-b);
	}


	template <class T>
	inline int
	bezierClassification(const T& k0, const T& k1, const T& k2, const T& k3,
						 T &inflexion, T &kk0, T &kk1, T &kk2)
	{
		//if (sign(k0) > 0 && sign(k1) < 0 && sign(k2) < 0 && sign(k3) < 0)
		if (k0.is_certainly_positive() &&
			k1.is_certainly_negative() &&
			k2.is_certainly_negative() &&
			k3.is_certainly_negative())
			return 0;

		//if (sign(k0) < 0 && sign(k1) > 0 && sign(k2) > 0 && sign(k3) > 0)
		if (k0.is_certainly_negative() &&
			k1.is_certainly_positive() &&
			k2.is_certainly_positive() &&
			k3.is_certainly_positive())
			return 0;

		//if (sign(k3) > 0 && sign(k1) < 0 && sign(k2) < 0 && sign(k0) < 0)
		if (k3.is_certainly_positive() &&
			k1.is_certainly_negative() &&
			k2.is_certainly_negative() &&
			k0.is_certainly_negative())
			return 0;

		//if (sign(k3) < 0 && sign(k1) > 0 && sign(k2) > 0 && sign(k0) > 0)
		if (k3.is_certainly_negative() &&
			k1.is_certainly_positive() &&
			k2.is_certainly_positive() &&
			k0.is_certainly_positive())
			return 0;

		// f'' = 6*(k2-2*k1+k0)*B^0_1 + 6*(k3-2*k2+k1)*B^1_1
		T a = k2-k1*T(2.0)+k0;
		T b = k3-k2*T(2.0)+k1;

		if (diffSign(a, b)) {
			//inflexion = lineRoot(a, b);
			return 2; // 1 inflexion
		}

		// f = 3*(k1-k0) B^2_0 + 3*(k2-k1)*B^2_1 + 3*(k3-k2)*B^2_2
		kk0 = k1-k0;
		kk1 = k2-k1;
		kk2 = k3-k2;

		if (diffSign(kk0, kk2))
			return 1; // no inflexion, 1 extreme
		else
			return 0; // no inflexion, no extreme
	}

	template <class T>
	inline int
	coplanarTest(bcrv<T> &c)
	{
		if (c.ct == 0) {// we only need to make sure sign(k0) != sign(k3)
			if (diffSign(c.k0, c.k3))
				return 1;
			else
				return 0;
		} else {
			if (diffSign(c.k0, c.k3))
				return 1;

			if ((c.kk0+c.kk2-c.kk1*T(2.0)).is_certainly_zero()) {
				printf("@@@degenerated ...\n");

				//T t = lineRoot(c.kk0, c.kk2);
				//T fk = evaluateBezier(c.k0, c.k1, c.k2, c.k3, t);
				T fk = _evaluateBezier(c.k0, c.k1, c.k2, c.k3, c.kk0, -c.kk2);
				if (c.kk0.is_certainly_negative())
					fk = -fk;

				if (sameSign(c.k0, fk))
					return 0;
				else
					return 2;
			}

			T s0, s1;
			T t0, t1;

			if (false == bezierDecomposition(c.k0, c.k1, c.k2, c.k3, c.kk0, c.kk1, c.kk2, s0, s1, t0, t1))
			{
				printf("@@@degenerated ...\n");
				//T t = lineRoot(c.kk0, c.kk2);
				//T fk = evaluateBezier(c.k0, c.k1, c.k2, c.k3, t);
				T fk = _evaluateBezier(c.k0, c.k1, c.k2, c.k3, c.kk0, -c.kk2);
				if (c.kk0.is_certainly_negative())
					fk = -fk;

				if (sameSign(c.k0, fk))
					return 0;
				else
					return 2;
			}

			if (sameSign(t0, t1)) {
				if (sameSign(c.k0, t0))
					return 0;
				else
					return 2;
			}

			//T t = lineRoot(t0, t1);
			//T fk = evaluateBezier2(c.kk0, c.kk1, c.kk2, t);
			T fk = _evaluateBezier2(c.k0, c.k1, c.k2, t0, -t1);

			if (sameSign(fk, c.kk0)) {
				if (sameSign(c.k0, t1))
					return 0;
				else
					return 2;
			} else {
				if (sameSign(c.k0, t0))
					return 0;
				else
					return 2;
			}
		}
	}



	template <class T>
	inline bool getSimplifyed(
			const T& k0, const T& k1, const T& k2, const T& k3,
			const T& l0, const T& l1, const T& l2, const T& l3, const T& l4,
			T& j0, T& j1, T& j2)
	{
		T kk0 = k0*T(4.0);
		T kk1 = k0+k1*T(3.0);
		T kk2 = (k1+k2)*T(2.0);
		T kk3 = k2*T(3.0)+k3;
		T kk4 = k3*T(4.0);

		T s0 = (l1*kk0 - l0*kk1)*T(12.0);
		T s1 = (l2*kk0 - l0*kk2)*T(6.0);
		T s2 = (l3*kk0 - l0*kk3)*T(4.0);
		T s3 = (l4*kk0 - l0*kk4)*T(3.0);

		j0 = (s1*k0-s0*k1)*T(6.0);
		j1 = (s2*k0-s0*k2)*T(3.0);
		j2 = (s3*k0-s0*k3)*T(2.0);

		return true;
	}


	template <class T>
	inline bool getSigns(const T& t0, const T& t1, bcrv<T> &c, T &lt0, T &lt1)
	{
		if (sameSign(t0, t1)) {
			lt0 = t0;
			lt1 = t0;
			return true;
		}

		if ((c.ct == 0) ||
			(c.ct == 1 && diffSign(c.k0, c.k3))) {
			//T t = lineRoot(t0, t1);
			//T ft = evaluateBezier(c.k0, c.k1, c.k2, c.k3, t);
			T ft = _evaluateBezier(c.k0, c.k1, c.k2, c.k3, t0, -t1);
			if (t0.is_certainly_negative())
				ft = -ft;

			if (sameSign(ft, c.k0)) {
				lt0 = t1;
				lt1 = t1;
			}
			else {
				lt0 = t0;
				lt1 = t0;
			}
			return true;
		}

		if (c.ct == 1) {
//				T t = lineRoot(t0, t1);
//				T ft = evaluateBezier(c.k0, c.k1, c.k2, c.k3, t);
			T ft = _evaluateBezier(c.k0, c.k1, c.k2, c.k3, t0, -t1);
			if (t0.is_certainly_negative())
				ft = -ft;

			if (diffSign(ft, c.k0)) {
				lt0 = t0;
				lt1 = t1;
				return true;
			}

			//T fk = evaluateBezier2(c.kk0, c.kk1, c.kk2, t);
			T fk = _evaluateBezier2(c.kk0, c.kk1, c.kk2, t0, -t1);

			if (sameSign(fk, c.kk0))
				lt0 = lt1 = t1;
			else
				lt0 = lt1 = t0;

			return true;
		}

		printf("Impossible to be here!\n");
		return false;
	}




	template <class T>
	inline bool insideTest(
			const Vec<3, T> &a0, const Vec<3, T> &b0, const Vec<3, T> &c0, const Vec<3, T> &d0,
			const Vec<3, T> &a1, const Vec<3, T> &b1, const Vec<3, T> &c1, const Vec<3, T> &d1,
			bcrv<T> &c, bool ee_test)
	{
		T l0, l1, l2, l3, l4;
		T j0, j1, j2;
		T s0, s1; // for L(t)
		T t0, t1; // for K(t)

		T lt0, lt1, kt0, kt1; // for signs of lt and kt
		bool bt0 = true, bt1 = true;

		for (int i=0; i<3; i++) {
			getBezier4(a0, b0, c0, d0, a1, b1, c1, d1, l0, l1, l2, l3, l4, i, ee_test);
			getSimplifyed(c.k0, c.k1, c.k2, c.k3, l0, l1, l2, l3, l4, j0, j1, j2);

			if ((j0+j2-j1*T(2.0)).is_certainly_zero()) {// degenerate j0, j1, j2
				//printf("@@@Degenerated...\n");
				getSigns(j0, j2, c, lt0, lt1);
				//printf("lt0=%lf, lt1=%lf\n", lt0, lt1);

				if (c.ct == 0) {
					//if (lt0 < 0)
					if (lt0.is_certainly_negative())
						return false;
				} else {
					//if (lt0 < 0)
					if (lt0.is_certainly_negative())
						bt0 = false;

					//if (lt1 < 0)
					if (lt1.is_certainly_negative())
						bt1 = false;

					if (!bt0 && !bt1)
						return false;
				}

				continue;
			}

			if (false == bezierDecomposition(c.k0, c.k1, c.k2, c.k3, j0, j1, j2, s0, s1, t0, t1)) {
				getSigns(j0, j2, c, lt0, lt1);

				if (c.ct == 0) {
					//if (lt0 < 0)
					if (lt0.is_certainly_negative())
						return false;
				} else {
					//if (lt0 < 0)
					if (lt0.is_certainly_negative())
						bt0 = false;

					//if (lt1 < 0)
					if (lt1.is_certainly_negative())
						bt1 = false;

					if (!bt0 && !bt1)
						return false;
				}

				continue;
			}

			getSigns(t0, t1, c, lt0, lt1);
			getSigns(s0, s1, c, kt0, kt1);

			if (c.ct == 0) {
				if (sameSign(lt0, kt0))
					return false;

				continue;
			}

			// kill an possiblity
			if (sameSign(lt0, kt0))
				bt0 = false;

			// kill an possiblity
			if (sameSign(lt1, kt1))
				bt1 = false;

			//if no possiblity left, return false ...
			if (!bt0 && !bt1)
				return false;
		}

		return true;
	}

	template <class T>
	inline bool
	Intersect_robust(
			const Vec<3, T> &a0, const Vec<3, T> &b0, const Vec<3, T> &c0, const Vec<3, T> &d0,
			const Vec<3, T> &a1, const Vec<3, T> &b1, const Vec<3, T> &c1, const Vec<3, T> &d1,
			bool ee_test)
	{
		T::begin_special_arithmetic();

		T k0, k1, k2, k3;

		if (!getBezier(a0, b0, c0, d0, a1, b1, c1, d1, k0, k1, k2, k3)) {
			T::end_special_arithmetic();
			return false;
		}

		T inflexion;
		T kk0, kk1, kk2;
		int ct = bezierClassification(k0, k1, k2, k3, inflexion, kk0, kk1, kk2);

		if (ct == 2) {
			T t = inflexion;
			Vec<3, T> at = a0*(T(1.0)-t)+a1*t;
			Vec<3, T> bt = b0*(T(1.0)-t)+b1*t;
			Vec<3, T> ct = c0*(T(1.0)-t)+c1*t;
			Vec<3, T> dt = d0*(T(1.0)-t)+d1*t;

			bool ret1 = Intersect_robust(a0, b0, c0, d0, at, bt, ct, dt, ee_test);
			bool ret2 = Intersect_robust(at, bt, ct, dt, a1, b1, c1, d1, ee_test);

			T::end_special_arithmetic();
			return ret1 || ret2;
		}

		bcrv<T> c(k0, k1, k2, k3, kk0, kk1, kk2, ct);

		if (!coplanarTest(c)) {
			T::end_special_arithmetic();
			return false;
		}

		if (!insideTest(a0, b0, c0, d0, a1, b1, c1, d1, c, ee_test)) {
			T::end_special_arithmetic();
			return false;
		}

		T::end_special_arithmetic();
		return true;
	}

	bool
	Intersect_VF_robust(
			const Vec3d &a0, const Vec3d &b0, const Vec3d &c0, const Vec3d &d0,
			const Vec3d &a1, const Vec3d &b1, const Vec3d &c1, const Vec3d &d1)
	{
		//DNF culling with conservative bound
		if (!DNF_Culling(a0, b0, c0, d0, a1, b1, c1, d1))
			return false;

		Vec3Interval ia0, ia1, ib0, ib1, ic0, ic1, id0, id1;
		make_vector(a0, ia0);
		make_vector(b0, ib0);
		make_vector(c0, ic0);
		make_vector(d0, id0);
		make_vector(a1, ia1);
		make_vector(b1, ib1);
		make_vector(c1, ic1);
		make_vector(d1, id1);
		bool ret = Intersect_robust(ia0, ib0, ic0, id0, ia1, ib1, ic1, id1, false);

		if (!ret) {
			Vec3e ea0, ea1, eb0, eb1, ec0, ec1, ed0, ed1;
			make_vector(a0, ea0);
			make_vector(b0, eb0);
			make_vector(c0, ec0);
			make_vector(d0, ed0);
			make_vector(a1, ea1);
			make_vector(b1, eb1);
			make_vector(c1, ec1);
			make_vector(d1, ed1);
			ret = Intersect_robust(ea0, eb0, ec0, ed0, ea1, eb1, ec1, ed1, false);
		}

		return ret;
	}

	bool
	Intersect_EE_robust(
			const Vec3d &a0, const Vec3d &b0, const Vec3d &c0, const Vec3d &d0,
			const Vec3d &a1, const Vec3d &b1, const Vec3d &c1, const Vec3d &d1)
	{
		//DNF culling with conservative bound
		if (!DNF_Culling(a0, b0, c0, d0, a1, b1, c1, d1))
			return false;

		Vec3Interval ia0, ia1, ib0, ib1, ic0, ic1, id0, id1;
		make_vector(a0, ia0);
		make_vector(b0, ib0);
		make_vector(c0, ic0);
		make_vector(d0, id0);
		make_vector(a1, ia1);
		make_vector(b1, ib1);
		make_vector(c1, ic1);
		make_vector(d1, id1);
		bool ret = Intersect_robust(ia0, ib0, ic0, id0, ia1, ib1, ic1, id1, true);

		if (!ret) {
			Vec3e ea0, ea1, eb0, eb1, ec0, ec1, ed0, ed1;
			make_vector(a0, ea0);
			make_vector(b0, eb0);
			make_vector(c0, ec0);
			make_vector(d0, ed0);
			make_vector(a1, ea1);
			make_vector(b1, eb1);
			make_vector(c1, ec1);
			make_vector(d1, ed1);
			ret = Intersect_robust(ea0, eb0, ec0, ed0, ea1, eb1, ec1, ed1, true);
		}

		return ret;
	}
}  // namespace rootparity
