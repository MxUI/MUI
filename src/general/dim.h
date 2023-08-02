/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library                       *
*                                                                            *
* Copyright (C) 2019 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis  *
*                                                                            *
* This software is jointly licensed under the Apache License, Version 2.0    *
* and the GNU General Public License version 3, you may use it according     *
* to either.                                                                 *
*                                                                            *
* ** Apache License, version 2.0 **                                          *
*                                                                            *
* Licensed under the Apache License, Version 2.0 (the "License");            *
* you may not use this file except in compliance with the License.           *
* You may obtain a copy of the License at                                    *
*                                                                            *
* http://www.apache.org/licenses/LICENSE-2.0                                 *
*                                                                            *
* Unless required by applicable law or agreed to in writing, software        *
* distributed under the License is distributed on an "AS IS" BASIS,          *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
* See the License for the specific language governing permissions and        *
* limitations under the License.                                             *
*                                                                            *
* ** GNU General Public License, version 3 **                                *
*                                                                            *
* This program is free software: you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by       *
* the Free Software Foundation, either version 3 of the License, or          *
* (at your option) any later version.                                        *
*                                                                            *
* This program is distributed in the hope that it will be useful,            *
* but WITHOUT ANY WARRANTY; without even the implied warranty of             *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
* GNU General Public License for more details.                               *
*                                                                            *
* You should have received a copy of the GNU General Public License          *
* along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
*****************************************************************************/

/**
 * @file dim.h
 * @author Y. H. Tang
 * @date 20 March 2014
 * @brief File providing data specialisms at different dimensionalities.
 */

#ifndef DIM_H_
#define DIM_H_

#include "util.h"

namespace mui {

namespace dim {

/******************************************************************************
                                Documentation
******************************************************************************/

// using namespace dim;
// using namespace dim::mechanical;
// using namespace dim::electrical;
// using namespace dim::magnetic;
// using namespace dim::optical;
// using namespace dim::radioactive;
// using namespace dim::chemical;

/******************************************************************************
                                Base class
******************************************************************************/

// MUI convention:
// [M]ass
// [L]ength
// [T]ime
// [T]emperature
// [A]mount of substance
// [E]lectrical current
// [L]uminos intensity
// amount of [I]nformation

template<int... MLTTAELI>
struct dim
{
	double a; // prefactor
	inline dim() : a(0) {}
	inline dim( double _a_ ) : a(_a_) {}
	inline dim( const dim &other ) : a(other.a) {}
	template<int... ANOTHER_MLTTAELI> dim( const dim<ANOTHER_MLTTAELI...> &error ) {
		error.Cannot_Assign_Quantities_With_Incompatible_Dimensionalities;
	}

	// copy-assign
	inline dim& operator = ( const dim &other ) { a = other.a; return *this; }
	template<int... ANOTHER_MLTTAELI>
		class Cannot_Assign_Quantities_With_Incompatible_Dimensionalities
		operator = ( const dim<ANOTHER_MLTTAELI...> &other );

	// get numeric value
	inline operator double () { return a; }

	// conversion
	inline dim convert_to( const dim &other ) {
		return dim( a / other.a );
	}
	template<int... ANOTHER_MLTTAELI>
		class Cannot_Assign_Quantities_With_Incompatible_Dimensionalities
		convert_to( const dim<ANOTHER_MLTTAELI...> &other );

	// serialization
	friend inline ostream& operator << ( ostream &out, dim q ) {
		out << q.a;
		return out;
	}
	friend inline istream& operator >> ( istream &in, dim q ) {
		in >> q.a;
		return in;
	}
};

/******************************************************************************
                                Arithmetic
******************************************************************************/

// addition
template<int... MLTTAELI1, int... MLTTAELI2> inline
class Cannot_Add_Values_With_Mismatching_Dimensionality operator + ( const dim<MLTTAELI1...> &d1, const dim<MLTTAELI2...> &d2 );

template<int... MLTTAELI> inline
dim<MLTTAELI...> operator + ( const dim<MLTTAELI...> &d1, const dim<MLTTAELI...> &d2 ) {
	return dim<MLTTAELI...>( d1.a + d2.a );
}

// negate
template<int... MLTTAELI> inline
dim<MLTTAELI...> operator - ( const dim<MLTTAELI...> &d ) {
	return dim<MLTTAELI...>( -d.a );
}

// substraction
template<int... MLTTAELI1, int... MLTTAELI2> inline
class Cannot_Subtract_Values_With_Mismatching_Dimensionality operator - ( const dim<MLTTAELI1...> &d1, const dim<MLTTAELI2...> &d2 );

template<int... MLTTAELI> inline
dim<MLTTAELI...> operator - ( const dim<MLTTAELI...> &d1, const dim<MLTTAELI...> &d2 ) {
	return dim<MLTTAELI...>( d1.a - d2.a );
}

// multiplication
// case 0: all dimensions cancel, result is scalar
template<int... MLTTAELI> inline
double operator * ( const dim<MLTTAELI...> &u1, const dim<-MLTTAELI...> &u2 ) {
	return u1.a * u2.a;
}
// case 1: result in something with new dimensionality
template<int... MLTTAELI1, int... MLTTAELI2> inline
dim<(MLTTAELI1+MLTTAELI2)...> operator * ( const dim<MLTTAELI1...> &u1, const dim<MLTTAELI2...> &u2 ) {
	dim<(MLTTAELI1+MLTTAELI2)...> u;
	u.a = u1.a * u2.a;
	return u;
}

// multiply by constant
template<int... MLTTAELI> inline
dim<MLTTAELI...> operator * ( const dim<MLTTAELI...> &u, const double f ) {
	return dim<MLTTAELI...>( u.a * f );
}
template<int... MLTTAELI> inline
dim<MLTTAELI...> operator * ( const double f, const dim<MLTTAELI...> &u ) {
	return dim<MLTTAELI...>( u.a * f );
}

// division
// case 0: all dimensions cancel, result is scalar
template<int... MLTTAELI> inline
double operator / ( const dim<MLTTAELI...> &u1, const dim<MLTTAELI...> &u2 ) {
	return u1.a / u2.a;
}
// case 1: result in something with new dimensionality
template<int... MLTTAELI1, int... MLTTAELI2> inline
dim<(MLTTAELI1-MLTTAELI2)...> operator / ( const dim<MLTTAELI1...> &u1, const dim<MLTTAELI2...> &u2 ) {
	dim<(MLTTAELI1-MLTTAELI2)...> u;
	u.a = u1.a / u2.a;
	return u;
}

// divide by constant
template<int... MLTTAELI> inline
dim<MLTTAELI...> operator / ( const dim<MLTTAELI...> &u, const double f ) {
	return u * (1.0/f);
}

// inversion
template<int... MLTTAELI> inline
dim<-MLTTAELI...> operator / ( const double f, const dim<MLTTAELI...> &u ) {
	return dim<-MLTTAELI...>( f / u.a );
}

/******************************************************************************
                              Dimensions
******************************************************************************/

// Fundamental dimensions
//                      M L T T A E L I
using mass        = dim<1,0,0,0,0,0,0,0>;
using length      = dim<0,1,0,0,0,0,0,0>;
using time        = dim<0,0,1,0,0,0,0,0>;
using temperature = dim<0,0,0,1,0,0,0,0>;
using amount      = dim<0,0,0,0,1,0,0,0>;
using current     = dim<0,0,0,0,0,1,0,0>;
using luminos     = dim<0,0,0,0,0,0,1,0>;
using information = dim<0,0,0,0,0,0,0,1>;

// Derived dimensions - General
using angle       = decltype( length() / length() );
using area        = decltype( length() * length() );
using frequency   = decltype( 1.0 / time() );
using velocity    = decltype( length() / time() );
using acceleration= decltype( velocity() / time() );
using force       = decltype( mass() * acceleration() );
using pressure    = decltype( force() / area() );
using energy      = decltype( force() * length() );
using solid_angle = decltype( area() / area() );
using power       = decltype( energy() / time() );

// Derived dimensions - By discipline
namespace mechanical {
}
namespace electrical {
	using charge      = decltype( current() * time() );
	using voltage     = decltype( power() / current() );
	using capacitance = decltype( charge() / voltage() );
	using resistance  = decltype( voltage() / current() );
	using conductance = decltype( current() / voltage() );
}
namespace magnetic {
	using electrical::voltage;
	using flux        = decltype( voltage() * time() );
	using strength    = decltype( flux() / area() );
	using inductance  = decltype( flux() / current() );
}
namespace optical {
	using flux        = decltype( luminos() * solid_angle() );
	using illuminance = decltype( flux() / area() );
}
namespace radioactive {
	using activity    = decltype( 1.0 / time() );
	using dose        = decltype( energy() / mass() );
}
namespace chemical {
	using catativity  = decltype( amount() / time() );
}

/******************************************************************************
                                 Units
******************************************************************************/

// Example
inline angle operator "" _deg ( long double u ) {
	return angle( u * ( PI / 360.0 ) );
}
inline angle operator "" _deg ( unsigned long long u ) {
	return angle( u * ( PI / 360.0 ) );
}

// Short-hand macro
#define make_unit(dimension,suffix,conversion) \
		inline dimension operator "" _##suffix ( long double u ) { \
			return dimension( conversion ); \
		} \
		inline dimension operator "" _##suffix ( unsigned long long u ) { \
			return dimension( conversion ); \
		}

// Fundamental units
make_unit(mass,ton,u*1e3)
make_unit(mass,kg,u)
make_unit(mass,g,u*1e-3)
make_unit(mass,mg,u*1e-6)
make_unit(length,km,u*1e3)
make_unit(length,m,u)
make_unit(length,cm,u*1e-2)
make_unit(length,mm,u*1e-3)
make_unit(length,um,u*1e-6)
make_unit(length,nm,u*1e-9)
make_unit(length,Angstrom,u*1e-10)
make_unit(time,day,u*86400.0)
make_unit(time,hr,u*3600.0)
make_unit(time,min,u*60.0)
make_unit(time,s,u)
make_unit(time,ms,u*1e-3)
make_unit(time,us,u*1e-6)
make_unit(time,ns,u*1e-9)
make_unit(time,ps,u*1e-12)
make_unit(time,fs,u*1e-15)
make_unit(temperature,K,u)
make_unit(temperature,C,u+273.15)
make_unit(amount,mol,u*6.02214178999999989284864e23)
make_unit(current,Amp,u)
make_unit(current,mAmp,u*1e-3)
make_unit(luminos,cd,u)
make_unit(information,bit,u*0.125)
make_unit(information,nibble,u*0.5)
make_unit(information,byte,u)
make_unit(information,KB,u*1e3)
make_unit(information,MB,u*1e6)
make_unit(information,GB,u*1e9)
make_unit(information,TB,u*1e12)
make_unit(information,PB,u*1e15)

// Derived dimensions - General
make_unit(angle,rad,u)
make_unit(frequency,Hz,u)
make_unit(velocity,kmph,u*(1000.0/3600.0))
make_unit(acceleration,G,u*9.80665)
make_unit(force,kN,u*1e3)
make_unit(force,N,u)
make_unit(force,mN,u*1e-3)
make_unit(force,uN,u*1e-6)
make_unit(force,nN,u*1e-9)
make_unit(force,pN,u*1e-12)
make_unit(pressure,GPa,u*1e9)
make_unit(pressure,MPa,u*1e6)
make_unit(pressure,kPa,u*1e3)
make_unit(pressure,Pa,u)
make_unit(energy,GJ,u*1e9)
make_unit(energy,MJ,u*1e6)
make_unit(energy,kJ,u*1e3)
make_unit(energy,J,u)
make_unit(energy,mJ,u*1e-3)
make_unit(energy,uJ,u*1e-6)
make_unit(energy,nJ,u*1e-9)
make_unit(energy,pJ,u*1e-12)
make_unit(energy,eV,u*1.602176565e-19)
make_unit(solid_angle,sr,u)
make_unit(power,MW,u*1e6)
make_unit(power,kW,u*1e3)
make_unit(power,W,u)
make_unit(power,mW,u*1e-3)
make_unit(power,uW,u*1e-6)

// Derived dimensions - By discipline
namespace mechanical {
}
namespace electrical {
	make_unit(charge,C,u)
	make_unit(charge,e,u*1.60217648700000002946104e-19)
	make_unit(voltage,kV,u*1e3)
	make_unit(voltage,V,u)
	make_unit(voltage,mV,u*1e-3)
	make_unit(capacitance,F,u)
	make_unit(capacitance,mF,u*1e-3)
	make_unit(capacitance,uF,u*1e-6)
	make_unit(capacitance,nF,u*1e-9)
	make_unit(capacitance,pF,u*1e-12)
	make_unit(resistance,Mohm,u*1e6)
	make_unit(resistance,kohm,u*1e3)
	make_unit(resistance,ohm,u)
	make_unit(conductance,S,u)
}
namespace magnetic {
	make_unit(flux,Wb,u)
	make_unit(strength,T,u)
	make_unit(inductance,H,u)
}
namespace optical {
	make_unit(flux,lm,u)
	make_unit(illuminance,lux,u)
}
namespace radioactive {
	make_unit(activity,Bq,u)
	make_unit(dose,Gy,u)
	make_unit(dose,Sv,u)
}
namespace chemical {
	make_unit(catativity,kat,u)
}

#undef make_unit

}

}

#endif /* DIM_H_ */
