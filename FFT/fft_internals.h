#ifndef FFT_INTERNALS_H
#define FFT_INTERNALS_H

namespace FFT_INTERNALS
{
	/**
	* here follow the definitions for all supported datatypes, It is quite important that the * operator returns an accurate result (easy for float and double).
	* The fixed point version(s) return a value with twice the size of the input operant to avoid loss of precision. The / operator is overloaded but not used (untested).
	*/

	struct f32_complex {
		using	ACCTYPE = f32_complex;
		using	DATATYPE = float;

		DATATYPE re;
		DATATYPE im;

		f32_complex operator +(const f32_complex& a) const { f32_complex ret; ret.re = re + a.re; ret.im = im + a.im; return ret; }
		f32_complex operator -(const f32_complex& a) const { f32_complex ret; ret.re = re - a.re; ret.im = im - a.im; return ret; }
		f32_complex operator *(const f32_complex& a) const { f32_complex ret; ret.re = re * a.re - im * a.im; ret.im = im * a.re + re * a.im; return ret; }
		f32_complex operator /(const f32_complex& a) const { f32_complex ret; float absval = a.re * a.re + a.im * a.im; ret.re = (re * a.re + im * a.im) / absval; ret.im = (im * a.re - re * a.im) / absval; return ret; }
	};

	struct f64_complex {

		using	ACCTYPE = f64_complex;
		using	DATATYPE = double;

		DATATYPE re;
		DATATYPE im;

		f64_complex operator +(const f64_complex& a) const { f64_complex ret; ret.re = re + a.re; ret.im = im + a.im; return ret; }
		f64_complex operator -(const f64_complex& a) const { f64_complex ret; ret.re = re - a.re; ret.im = im - a.im; return ret; }
		f64_complex operator *(const f64_complex& a) const { f64_complex ret; ret.re = re * a.re - im * a.im; ret.im = im * a.re + re * a.im; return ret; }
		f64_complex operator /(const f64_complex& a) const { f64_complex ret; double absval = a.re * a.re + a.im * a.im; ret.re = (re * a.re + im * a.im) / absval; ret.im = (im * a.re - re * a.im) / absval; return ret; }
	};

	struct i64_complex;
	struct i32_complex;
	struct i16_complex;

	struct i64_complex {
		using ACCTYPE = i64_complex;
		using OFWTYPE = long long;
		using DATATYPE = long long;

		DATATYPE re;
		DATATYPE im;

		i64_complex operator +(const i64_complex& a) const { i64_complex ret; ret.re = re + a.re; ret.im = im + a.im; return ret; }
		i64_complex operator -(const i64_complex& a) const { i64_complex ret; ret.re = re - a.re; ret.im = im - a.im; return ret; }
		i64_complex operator *(const i64_complex& a) const { i64_complex ret; ret.re = re * a.re - im * a.im; ret.im = im * a.re + re * a.im; return ret; }
		i64_complex operator /(const i64_complex& a) const { i64_complex ret; long long absval = a.re * a.re + a.im * a.im; ret.re = (re * a.re + im * a.im) / absval; ret.im = (im * a.re - re * a.im) / absval; return ret; }//that is bullshit
		i64_complex operator << (const size_t shift) const { i64_complex retVal; retVal.re = static_cast<long long>(re) << shift; retVal.im = static_cast<long long>(im) << shift; return retVal; }
		i64_complex operator >> (const size_t shift) const { i64_complex retVal; retVal.re = static_cast<long long>(re) >> shift; retVal.im = static_cast<long long>(im) >> shift; return retVal; }
		operator i16_complex();
		operator i32_complex();
	};

	struct i32_complex {
		using ACCTYPE = i64_complex;
		using OFWTYPE = long long;
		using DATATYPE = int;

		DATATYPE re;
		DATATYPE im;

		i32_complex operator +(const i32_complex& a) const { i32_complex ret; ret.re = static_cast<int>(re) + a.re; ret.im = static_cast<int>(im) + a.im; return ret; }
		i32_complex operator -(const i32_complex& a) const { i32_complex ret; ret.re = static_cast<int>(re) - a.re; ret.im = static_cast<int>(im) - a.im; return ret; }
		i64_complex operator *(const i32_complex& a) const { i64_complex ret; ret.re = static_cast<long long>(re) * a.re - static_cast<long long>(im) * a.im; ret.im = static_cast<long long>(im) * a.re + static_cast<long long>(re) * a.im; return ret; }
		i64_complex operator /(const i32_complex& a) const { i64_complex ret; int absval = a.re * a.re + a.im * a.im; ret.re = (re * a.re + im * a.im) / absval; ret.im = (im * a.re - re * a.im) / absval; return ret; }//that is bullshit
		i64_complex operator << (const size_t shift) const { i64_complex retVal; retVal.re = static_cast<long long>(re) << shift; retVal.im = static_cast<long long>(im) << shift; return retVal; }
		i64_complex operator >> (const size_t shift) const { i64_complex retVal; retVal.re = static_cast<long long>(re) >> shift; retVal.im = static_cast<long long>(im) >> shift; return retVal; }
		operator i16_complex();
	};


	struct i16_complex {
		using ACCTYPE = i32_complex;
		using OFWTYPE = int;
		using DATATYPE = short;
		DATATYPE re;
		DATATYPE im;

		i32_complex operator +(const i16_complex& a) const { i32_complex ret; ret.re = static_cast<int>(re) + a.re; ret.im = static_cast<int>(im) + a.im; return ret; }
		i32_complex operator -(const i16_complex& a) const { i32_complex ret; ret.re = static_cast<int>(re) - a.re; ret.im = static_cast<int>(im) - a.im; return ret; }
		i32_complex operator *(const i16_complex& a) const { i32_complex ret; ret.re = static_cast<int>(re) * a.re - static_cast<int>(im) * a.im; ret.im = static_cast<int>(im) * a.re + static_cast<int>(re) * a.im; return ret; }
		i32_complex operator /(const i16_complex& a) const { i32_complex ret; short absval = a.re * a.re + a.im * a.im; ret.re = (re * a.re + im * a.im) / absval; ret.im = (im * a.re - re * a.im) / absval; return ret; }//that is bullshit}
		i32_complex operator << (const size_t shift) const { i32_complex retVal; retVal.re = static_cast<int>(re) << shift; retVal.im = static_cast<int>(im) << shift; return retVal; }
		i32_complex operator >> (const size_t shift) const { i32_complex retVal; retVal.re = static_cast<int>(re) >> shift; retVal.im = static_cast<int>(im) >> shift; return retVal; }
	};

	i64_complex::operator i16_complex() { i16_complex ret; ret.re = static_cast<i16_complex::DATATYPE>(re); ret.im = static_cast<i16_complex::DATATYPE>(im); return ret; }
	i64_complex::operator i32_complex() { i32_complex ret; ret.re = static_cast<i32_complex::DATATYPE>(re); ret.im = static_cast<i32_complex::DATATYPE>(im); return ret; }
	i32_complex::operator i16_complex() { i16_complex ret; ret.re = static_cast<i16_complex::DATATYPE>(re); ret.im = static_cast<i32_complex::DATATYPE>(im); return ret; }
}

#endif //FFT_INTERNALS
