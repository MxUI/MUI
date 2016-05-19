/* $Id$ */

#ifndef MUI_POINT_H_
#define MUI_POINT_H_

#include<array>
#include<cmath>
#include<cassert>

namespace mui {

using uint = unsigned int;

// error-checking heplers
template<typename T1, typename T2> struct same_type { static bool const yes = false; };
template<typename T> struct same_type<T,T> { static bool const yes = true; };

// Reference:
// Expression template: http://en.wikipedia.org/wiki/Expression_templates
// CRTP: http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern

/*---------------------------------------------------------------------------
                              Interface
---------------------------------------------------------------------------*/

template<class VEX, typename SCALAR, uint D>
struct vexpr {
	using TYPE_ = SCALAR;
	static const uint D_ = D;

	// the only work in constructor is type checking
	inline vexpr() {
		static_assert( same_type<TYPE_, typename VEX::TYPE_>::yes, "point element type mismatch" );
		static_assert( D_ == VEX::D_, "point dimensionality mismatch" );
	}

	// dereferencing using static polymorphism
	inline SCALAR operator[] (uint i) const {
		assert( i < D );
		return static_cast<VEX const &>(*this)[i];
	}

	inline operator VEX      & ()       { return static_cast<VEX      &>(*this); }
	inline operator VEX const& () const { return static_cast<VEX const&>(*this); }

	inline uint d() const { return D_; }
};

/*---------------------------------------------------------------------------
                                 Container
---------------------------------------------------------------------------*/

template<typename SCALAR, uint D>
struct point : public vexpr<point<SCALAR,D>, SCALAR, D> {
protected:
	SCALAR x[D];
public:
	using TYPE_ = SCALAR;
	static const uint D_ = D;

	// default constructor
	inline point() {}
	// construct from scalar constant
	explicit inline point(float  const s) { for(uint i = 0 ; i < D ; i++) x[i] = s; }
	explicit inline point(double const s) { for(uint i = 0 ; i < D ; i++) x[i] = s; }
	explicit inline point(int    const s) { for(uint i = 0 ; i < D ; i++) x[i] = s; }
	explicit inline point(uint   const s) { for(uint i = 0 ; i < D ; i++) x[i] = s; }
	explicit inline point(long   const s) { for(uint i = 0 ; i < D ; i++) x[i] = s; }
	explicit inline point(ulong  const s) { for(uint i = 0 ; i < D ; i++) x[i] = s; }
	// construct from C-array
	explicit inline point(float  const *ps) { for(uint i = 0 ; i < D ; i++) x[i] = ps[i]; }
	explicit inline point(double const *ps) { for(uint i = 0 ; i < D ; i++) x[i] = ps[i]; }
	explicit inline point(int    const *ps) { for(uint i = 0 ; i < D ; i++) x[i] = ps[i]; }
	explicit inline point(uint   const *ps) { for(uint i = 0 ; i < D ; i++) x[i] = ps[i]; }
	explicit inline point(long   const *ps) { for(uint i = 0 ; i < D ; i++) x[i] = ps[i]; }
	explicit inline point(ulong  const *ps) { for(uint i = 0 ; i < D ; i++) x[i] = ps[i]; }
	// construct from parameter pack
	// 'head' differentiate it from constructing from vector expression
	template<typename ...T> inline point(SCALAR const head, T const ... tail ) {
		std::array<TYPE_,D_> s( { head, static_cast<SCALAR>(tail)... } );
		for(uint i = 0 ; i < D ; i++) x[i] = s[i];
	}
	// construct from any vector expression
	template<class E> inline point( const vexpr<E,SCALAR,D> &u ) {
		for(uint i = 0 ; i < D ; i++) x[i] = u[i];
	}

	// point must be assignable, while other expressions may not
	inline SCALAR      & operator [] (uint i)       { assert( i < D ); return x[i]; }
	inline SCALAR const& operator [] (uint i) const { assert( i < D ); return x[i]; }

	// STL-style direct data accessor
	inline SCALAR      * data()       { return x; }
	inline SCALAR const* data() const { return x; }

	// assign from any vector expression
	template<class E> inline point & operator += ( const vexpr<E,SCALAR,D> &u ) {
		for(uint i = 0 ; i < D ; i++) x[i] += u[i];
		return *this;
	}
	template<class E> inline point & operator -= ( const vexpr<E,SCALAR,D> &u ) {
		for(uint i = 0 ; i < D ; i++) x[i] -= u[i];
		return *this;
	}
	template<class E> inline point & operator *= ( const vexpr<E,SCALAR,D> &u ) {
		for(uint i = 0 ; i < D ; i++) x[i] *= u[i];
		return *this;
	}
	template<class E> inline point & operator /= ( const vexpr<E,SCALAR,D> &u ) {
		for(uint i = 0 ; i < D ; i++) x[i] /= u[i];
		return *this;
	}
	// conventional vector-scalar operators
	inline point & operator += ( SCALAR const u ) {
		for(uint i = 0 ; i < D ; i++) x[i] += u;
		return *this;
	}
	inline point & operator -= ( SCALAR const u ) {
		for(uint i = 0 ; i < D ; i++) x[i] -= u;
		return *this;
	}
	inline point & operator *= ( SCALAR const u ) {
		for(uint i = 0 ; i < D ; i++) x[i] *= u;
		return *this;
	}
	inline point & operator /= ( SCALAR const u ) {
		for(uint i = 0 ; i < D ; i++) x[i] /= u;
		return *this;
	}
};

/*---------------------------------------------------------------------------
                         Arithmetic Functors
---------------------------------------------------------------------------*/

template<class E1, class E2, typename SCALAR, uint D>
struct vexpr_add: public vexpr<vexpr_add<E1,E2,SCALAR,D>, SCALAR, D> {
	inline vexpr_add( vexpr<E1,SCALAR,D> const& u, vexpr<E2,SCALAR,D> const& v ) : u_(u), v_(v) {}
	inline SCALAR operator [] (uint i) const { return u_[i] + v_[i]; }
protected:
	E1 const& u_;
	E2 const& v_;
};

template<class E1, class E2, typename SCALAR, uint D>
struct vexpr_sub: public vexpr<vexpr_sub<E1,E2,SCALAR,D>, SCALAR, D> {
	inline vexpr_sub( vexpr<E1,SCALAR,D> const& u, vexpr<E2,SCALAR,D> const& v ) : u_(u), v_(v) {}
	inline SCALAR operator [] (uint i) const { return u_[i] - v_[i]; }
protected:
	E1 const& u_;
	E2 const& v_;
};

template<class E, typename SCALAR, uint D>
struct vexpr_neg: public vexpr<vexpr_neg<E,SCALAR,D>, SCALAR, D> {
	inline vexpr_neg( vexpr<E,SCALAR,D> const& u ) : u_(u) {}
	inline SCALAR operator [] (uint i) const { return -u_[i]; }
protected:
	E const& u_;
};

template<class E1, class E2, typename SCALAR, uint D>
struct vexpr_mul: public vexpr<vexpr_mul<E1,E2,SCALAR,D>, SCALAR, D> {
	inline vexpr_mul( vexpr<E1,SCALAR,D> const& u, vexpr<E2,SCALAR,D> const& v ) : u_(u), v_(v) {}
	inline SCALAR operator [] (uint i) const { return u_[i] * v_[i]; }
protected:
	E1 const& u_;
	E2 const& v_;
};

template<class E, typename SCALAR, uint D>
struct vexpr_scale: public vexpr<vexpr_scale<E,SCALAR,D>, SCALAR, D> {
	inline vexpr_scale( vexpr<E,SCALAR,D> const& u, SCALAR const a ) : u_(u), a_(a) {}
	inline SCALAR operator [] (uint i) const { return u_[i] * a_; }
protected:
	E      const& u_;
	SCALAR const  a_;
};

template<class E1, class E2, typename SCALAR, uint D>
struct vexpr_div: public vexpr<vexpr_div<E1,E2,SCALAR,D>, SCALAR, D> {
	inline vexpr_div( vexpr<E1,SCALAR,D> const& u, vexpr<E2,SCALAR,D> const& v ) : u_(u), v_(v) {}
	inline SCALAR operator [] (uint i) const { return u_[i] / v_[i]; }
protected:
	E1 const& u_;
	E2 const& v_;
};

template<class E, typename SCALAR, uint D>
struct vexpr_rscale: public vexpr<vexpr_rscale<E,SCALAR,D>, SCALAR, D> {
	inline vexpr_rscale( vexpr<E,SCALAR,D> const& u, SCALAR const a ) : u_(u), a_(a) {}
	inline SCALAR operator [] (uint i) const { return u_[i] / a_; }
protected:
	E      const& u_;
	SCALAR const  a_;
};

template<class E, typename SCALAR, uint D>
struct vexpr_rcp: public vexpr<vexpr_rcp<E,SCALAR,D>, SCALAR, D> {
	inline vexpr_rcp( vexpr<E,SCALAR,D> const& u ) : u_(u) {}
	inline SCALAR operator [] (uint i) const { return SCALAR(1)/u_[i]; }
protected:
	E const& u_;
};

template<class E, typename SCALAR, uint D>
struct vexpr_scale_rcp: public vexpr<vexpr_scale_rcp<E,SCALAR,D>, SCALAR, D> {
	inline vexpr_scale_rcp( SCALAR const a, vexpr<E,SCALAR,D> const& u ) : a_(a), u_(u) {}
	inline SCALAR operator [] (uint i) const { return a_ / u_[i]; }
protected:
	SCALAR const  a_;
	E      const& u_;
};

template<class E, class OP, typename SCALAR, uint D>
struct vexpr_apply1: public vexpr<vexpr_apply1<E,OP,SCALAR,D>, SCALAR, D> {
	inline vexpr_apply1( vexpr<E,SCALAR,D> const& u, OP const& op ) : u_(u), o_(op) {}
	inline SCALAR operator [] (uint i) const { return o_( u_[i] ); }
protected:
	E  const& u_;
	OP const& o_;
};

template<class E1, class E2, class OP, typename SCALAR, uint D>
struct vexpr_apply2: public vexpr<vexpr_apply2<E1,E2,OP,SCALAR,D>, SCALAR, D> {
	inline vexpr_apply2( vexpr<E1,SCALAR,D> const& u, vexpr<E2,SCALAR,D> const& v, OP const& op ) : u_(u), v_(v), o_(op) {}
	inline SCALAR operator [] (uint i) const { return o_( u_[i], v_[i] ); }
protected:
	E1 const& u_;
	E2 const& v_;
	OP const& o_;
};

/*---------------------------------------------------------------------------
                         Operator Overloads
---------------------------------------------------------------------------*/

template<class E1, class E2, typename SCALAR, uint D> inline
vexpr_add<E1, E2, SCALAR, D> operator + ( vexpr<E1,SCALAR,D> const &u, vexpr<E2,SCALAR,D> const &v ) {
	return vexpr_add<E1, E2, SCALAR, D>( u, v );
}

template<class E1, class E2, typename SCALAR, uint D> inline
vexpr_sub<E1, E2, SCALAR, D> operator - ( vexpr<E1,SCALAR,D> const &u, vexpr<E2,SCALAR,D> const &v ) {
	return vexpr_sub<E1, E2, SCALAR, D>( u, v );
}

template<class E, typename SCALAR, uint D> inline
vexpr_neg<E, SCALAR, D> operator - ( vexpr<E,SCALAR,D> const &u ) {
	return vexpr_neg<E, SCALAR, D>( u );
}

template<class E1, class E2, typename SCALAR, uint D> inline
vexpr_mul<E1, E2, SCALAR, D> operator * ( vexpr<E1,SCALAR,D> const &u, vexpr<E2,SCALAR,D> const &v ) {
	return vexpr_mul<E1, E2, SCALAR, D>( u, v );
}

template<class E, typename SCALAR, uint D> inline
vexpr_scale<E, SCALAR, D> operator * ( vexpr<E,SCALAR,D> const &u, SCALAR const a ) {
	return vexpr_scale<E, SCALAR, D>( u, a );
}

template<class E, typename SCALAR, uint D> inline
vexpr_scale<E, SCALAR, D> operator * ( SCALAR const a, vexpr<E,SCALAR,D> const &u ) {
	return vexpr_scale<E, SCALAR, D>( u, a );
}

template<class E, typename SCALAR, uint D> inline
vexpr_scale<E, SCALAR, D> operator / ( vexpr<E,SCALAR,D> const &u, SCALAR const a ) {
	return vexpr_rscale<E, SCALAR, D>( u, a );
}

template<class E1, class E2, typename SCALAR, uint D> inline
vexpr_div<E1, E2, SCALAR, D> operator / ( vexpr<E1,SCALAR,D> const &u, vexpr<E2,SCALAR,D> const &v ) {
	return vexpr_div<E1, E2, SCALAR, D>( u, v );
}

template<class E, typename SCALAR, uint D> inline
vexpr_scale_rcp<E, SCALAR, D> operator / ( SCALAR const a, vexpr<E,SCALAR,D> const &u ) {
	return vexpr_scale_rcp<E, SCALAR, D>( a, u );
}

/*---------------------------------------------------------------------------
                         Math functions
---------------------------------------------------------------------------*/

template<class E1, class E2, typename SCALAR> inline
point<SCALAR,3U> cross( vexpr<E1,SCALAR,3U> const &u, vexpr<E2,SCALAR,3U> const &v ) {
        point<SCALAR,3U> x;
        for(uint i=0;i<3;i++)  x[i] = u[(i+1U)%3U] * v[(i+2U)%3U] - u[(i+2U)%3U] * v[(i+1U)%3U];
        return x;
}

// generic reduction template
template<class E, class OP, typename SCALAR, uint D> inline
SCALAR reduce( vexpr<E,SCALAR,D> const &u, OP const & op ) {
	SCALAR core( u[0] );
	for(uint i = 1 ; i < D ; i++) core = op( core, u[i] );
	return core;
}

// biggest element within a vector
template<class E, typename SCALAR, uint D> inline
SCALAR max( vexpr<E,SCALAR,D> const &u ) {
	return reduce( u, [](SCALAR a, SCALAR b){return a>b?a:b;} );
}

// smallest element within a vector
template<class E, typename SCALAR, uint D> inline
SCALAR min( vexpr<E,SCALAR,D> const &u ) {
	return reduce( u, [](SCALAR a, SCALAR b){return a<b?a:b;} );
}

// smallest element within a vector
template<class E, typename SCALAR, uint D> inline
SCALAR sum( vexpr<E,SCALAR,D> const &u ) {
	return reduce( u, [](SCALAR a, SCALAR b){return a+b;} );
}

// smallest element within a vector
template<class E, typename SCALAR, uint D> inline
SCALAR mean( vexpr<E,SCALAR,D> const &u ) {
	return sum(u) / double(D);
}

// inner product
template<class E1, class E2, typename SCALAR, uint D> inline
SCALAR dot( vexpr<E1,SCALAR,D> const &u, vexpr<E2,SCALAR,D> const &v ) {
	return sum( u * v );
}

// square of L2 norm
template<class E, typename SCALAR, uint D> inline
SCALAR normsq( vexpr<E,SCALAR,D> const &u ) {
	return sum( u * u );
}

// L2 norm
template<class E, typename SCALAR, uint D> inline
SCALAR norm( vexpr<E,SCALAR,D> const &u ) {
	return std::sqrt( normsq(u) );
}

// element-wise arbitrary function applied for each element
template<class E, class OP, typename SCALAR, uint D> inline
vexpr_apply1<E, OP, SCALAR, D> apply( vexpr<E,SCALAR,D> const &u, OP const& op ) {
	return vexpr_apply1<E, OP, SCALAR, D>( u, op );
}

// element-wise arbitrary function applied element-wisely between 2 vectors
template<class E1, class E2, class OP, typename SCALAR, uint D> inline
vexpr_apply2<E1, E2, OP, SCALAR, D> apply( vexpr<E1,SCALAR,D> const &u, vexpr<E2,SCALAR,D> const &v, OP const& op ) {
	return vexpr_apply2<E1, E2, OP, SCALAR, D>( u, v, op );
}

}

#endif /* POINT_H_ */
