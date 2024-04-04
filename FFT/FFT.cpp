// FFT.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <math.h>
#include <assert.h>
#include <iostream>
#include "plot.h"

typedef enum EN_ScalingMethod{SCALE_INPUT, SCALE_OUTPUT, SCALEINPOUT}EN_ScalingMethod;

//#define _USETWIDDLE_
#define _RECURRENCE_
//#define _C1x_

/**
* here follow the definitions for all supported datatypes, It is quite important that the * operator returns an accurate result (easy for float and double).
* The fixed point version(s) return a value with twice the size of the input operant to avoid loss of precision. The / operator is overloaded but not used (untested).
*/

struct FComplex {
	typedef FComplex	ACCTYPE;
	typedef float		DATATYPE;

	float re;
	float im;

	FComplex operator +(const FComplex &a){FComplex ret;ret.re=re+a.re;ret.im=im+a.im;return ret;}
	FComplex operator -(const FComplex &a){FComplex ret;ret.re=re-a.re;ret.im=im-a.im;return ret;}
	FComplex operator *(const FComplex &a){FComplex ret;ret.re=re*a.re-im*a.im;ret.im=im*a.re+re*a.im;return ret;}
	FComplex operator /(const FComplex &a){FComplex ret;float absval = a.re*a.re+a.im*a.im;ret.re=(re*a.re+im*a.im)/absval;ret.im=(im*a.re-re*a.im)/absval;return ret;}
	FComplex operator << (const size_t shift){(void)shift;return *this;}
};

struct DComplex {
	
	typedef DComplex	ACCTYPE;
	typedef double		DATATYPE;

	double re;
	double im;

	DComplex operator +(const DComplex &a){DComplex ret;ret.re=re+a.re;ret.im=im+a.im;return ret;}
	DComplex operator -(const DComplex &a){DComplex ret;ret.re=re-a.re;ret.im=im-a.im;return ret;}
	DComplex operator *(const DComplex &a){DComplex ret;ret.re=re*a.re-im*a.im;ret.im=im*a.re+re*a.im;return ret;}
	DComplex operator /(const DComplex &a){DComplex ret;double absval = a.re*a.re+a.im*a.im;ret.re=(re*a.re+im*a.im)/absval;ret.im=(im*a.re-re*a.im)/absval;return ret;}
	DComplex operator << (const size_t shift){(void)shift;return *this;}
};

struct i64Complex;
struct i32Complex;
struct i16Complex;

struct i64Complex {
	typedef i64Complex	ACCTYPE;
	typedef long long	OFWTYPE;
	typedef long long	DATATYPE;

	long long re;
	long long im;

	i64Complex operator +(const i64Complex &a){i64Complex ret;ret.re=re+a.re;ret.im=im+a.im;return ret;}
	i64Complex operator -(const i64Complex &a){i64Complex ret;ret.re=re-a.re;ret.im=im-a.im;return ret;}
	i64Complex operator *(const i64Complex &a){i64Complex ret;ret.re=re*a.re-im*a.im;ret.im=im*a.re+re*a.im;return ret;}
	i64Complex operator /(const i64Complex &a){i64Complex ret;long long absval = a.re*a.re+a.im*a.im;ret.re=(re*a.re+im*a.im)/absval;ret.im=(im*a.re-re*a.im)/absval;return ret;}//that is bullshit
	i64Complex operator << (const size_t shift){i64Complex retVal;retVal.re=static_cast<long long>(re)<<shift;retVal.im=static_cast<long long>(im)<<shift;return retVal;}
	operator i16Complex();
	operator i32Complex();
};

struct i32Complex {
	typedef i64Complex	ACCTYPE;
	typedef long long	OFWTYPE;
	typedef int			DATATYPE;

	int re;
	int im;

	i32Complex operator +(const i32Complex &a){i32Complex ret;ret.re=static_cast<int>(re)+a.re;ret.im=static_cast<int>(im)+a.im;return ret;}
	i32Complex operator -(const i32Complex &a){i32Complex ret;ret.re=static_cast<int>(re)-a.re;ret.im=static_cast<int>(im)-a.im;return ret;}
	i64Complex operator *(const i32Complex &a){i64Complex ret;ret.re=static_cast<long long>(re)*a.re-static_cast<long long>(im)*a.im;ret.im=static_cast<long long>(im)*a.re+static_cast<long long>(re)*a.im;return ret;}
	i64Complex operator /(const i32Complex &a){i64Complex ret;int absval = a.re*a.re+a.im*a.im;ret.re=(re*a.re+im*a.im)/absval;ret.im=(im*a.re-re*a.im)/absval;return ret;}//that is bullshit
	i64Complex operator << (const size_t shift){i64Complex retVal;retVal.re=static_cast<long long>(re)<<shift;retVal.im=static_cast<long long>(im)<<shift;return retVal;}
	operator i16Complex();
};


struct i16Complex {
	typedef i32Complex	ACCTYPE;
	typedef int			OFWTYPE;
	typedef short		DATATYPE;
	short re;
	short im;

	i32Complex operator +(const i16Complex &a){i32Complex ret;ret.re=static_cast<int>(re)+a.re;ret.im=static_cast<int>(im)+a.im;return ret;}
	i32Complex operator -(const i16Complex &a){i32Complex ret;ret.re=static_cast<int>(re)-a.re;ret.im=static_cast<int>(im)-a.im;return ret;}
	i32Complex operator *(const i16Complex &a){i32Complex ret;ret.re=static_cast<int>(re)*a.re-static_cast<int>(im)*a.im;ret.im=static_cast<int>(im)*a.re+static_cast<int>(re)*a.im;return ret;}
	i32Complex operator /(const i16Complex &a){i32Complex ret;short absval = a.re*a.re+a.im*a.im;ret.re=(re*a.re+im*a.im)/absval;ret.im=(im*a.re-re*a.im)/absval;return ret;}//that is bullshit}
	i32Complex operator << (const size_t shift){i32Complex retVal;retVal.re=static_cast<int>(re)<<shift;retVal.im=static_cast<int>(im)<<shift;return retVal;}
};

 i64Complex::operator i16Complex() {i16Complex ret;ret.re=re>>48;ret.im=im>>48;return ret;}
 i64Complex::operator i32Complex() {i32Complex ret;ret.re=re>>32;ret.im=im>>32;return ret;}
 i32Complex::operator i16Complex() {i16Complex ret;ret.re=re>>16;ret.im=im>>16;return ret;}




/**
* fixed point implementation of the FFT algorithm. 
* The fixed point version(s) return a value with twice the size of the input operant to avoid loss of precision
*/

template<size_t loglen, typename TW, typename T, typename TFFT> class CFFT {
	private:
		CFFT();
		TW      ma_TwiddleFactors[1<<(loglen-1)];

		void    CreateTwiddleFactors(size_t length, TW twiddlefactors[], double scaling);

		TW      GetCosFactor(size_t stage, size_t stageshift);
		TW      GetSinFactor(size_t stage, size_t stageshift);


		template <size_t bitsize> size_t swap(size_t pointer)
		{
			const size_t size			= bitsize>>1;//for an example for 11 size is 5
			const size_t lowermask		= ~(static_cast<size_t>(static_cast<size_t>(-1)>>size)<<size);
			const size_t uppermask		= lowermask<<(bitsize-size);
			const size_t mask			= ~(static_cast<size_t>(static_cast<size_t>(-1)>>bitsize)<<bitsize);
			const size_t remaining		= (~(lowermask|uppermask))&mask;
			pointer						= (swap<size>(lowermask&pointer)<<(bitsize-size))|swap<size>((uppermask&pointer)>>(bitsize-size))|(remaining&pointer);
			return pointer;
		}

		template <> size_t swap<2>(size_t pointer)
		{
			return ((pointer<<1)|(pointer>>1))&3;
		}

		template <> size_t swap<1>(size_t pointer)
		{
			return pointer&1;
		}

		template <size_t bitsize> size_t bitreversal(size_t pointer)
		{
			const size_t size		= (bitsize)>>1;//for an example for 11 size is 5
			const size_t lowermask	= ~(static_cast<size_t>(static_cast<size_t>(-1)>>size)<<size);
			const size_t uppermask	= lowermask<<(bitsize-size);
			const size_t mask		= ~(static_cast<size_t>(static_cast<size_t>(-1)>>(bitsize))<<(bitsize));
			const size_t remaining	= (~(lowermask|uppermask))&mask;
			pointer					= (swap<size>(lowermask&pointer)<<(bitsize-size))|swap<size>((uppermask&pointer)>>(bitsize-size))|(remaining&pointer);
			return pointer;
		}

		void ReorderSamples(T in[], TFFT data[]);
		void ReorderSamples(TFFT in[], TFFT data[]);
		void ReorderRealSamples(T in[], TFFT data[]);
		void ReorderRealSamplesInPlace(TFFT in[]);
		void Convert2HalfDFT(TFFT transform[]);
		void Convert2HalfDFT(TFFT transform[], size_t exponent);
		void Convert2ComplexDFT(TFFT transform[]);
		void Convert2ComplexDFT(TFFT transform[], size_t exponent);
	public:
		CFFT(double scaling){CreateTwiddleFactors(sizeof(ma_TwiddleFactors)/sizeof(ma_TwiddleFactors[0]), ma_TwiddleFactors, scaling);}
		void CalculateFFT(T buffer[], TFFT ffttransform[]);
		void CalculateFFT(T buffer[], TFFT ffttransform[], EN_ScalingMethod rescale, size_t *pBlockExponent);
		void CalculateIFFT(TFFT buffer[], TFFT ffttransform[]);
		void CalculateIFFT(TFFT buffer[], TFFT ffttransform[], EN_ScalingMethod rescale, size_t *pBlockExponent);
		void CalculateRealFFT(T buffer[], TFFT ffttransform[]);
		void CalculateRealFFT(T buffer[], TFFT ffttransform[], EN_ScalingMethod rescale, size_t *pBlockExponent);
		void CalculateRealIFFT(TFFT buffer[], TFFT ffttransform[]);
		void CalculateRealIFFT(TFFT buffer[], TFFT ffttransform[], EN_ScalingMethod rescale, size_t *pBlockExponent);
};

template<size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::CreateTwiddleFactors(size_t length, TW twiddlefactors[], double scaling)
{
	static const double PI2 = 6.283185307179586476925286f;

	size_t len        = 1<<loglen;
	if (2*length>=len)
	{
		const size_t len2 = len>>1;
		for (size_t i=0;i<len2;i++)
		{
			twiddlefactors[i] = static_cast<TW>(floor(scaling*cos(PI2*static_cast<double>(i)/len)+.5));
		}
	}
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::ReorderSamples(T in[], TFFT data[])
{
	size_t length = 1<<loglen;
	for (size_t i=0;i<length;i++)
	{
		size_t reordered	= bitreversal<loglen>(i);
		data[reordered].re	= in[i];
		data[reordered].im	= 0;
	}
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::ReorderSamples(TFFT in[], TFFT data[])
{
	size_t length = 1<<loglen;
	for (size_t i=0;i<length;i++)
	{
		size_t reordered	= bitreversal<loglen>(i);
		data[reordered]		= in[i];
	}
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::ReorderRealSamples(T in[], TFFT data[])
{
	size_t length = 1<<(loglen-1);
	for (size_t i=0;i<length;i++)
	{
		size_t reordered	= bitreversal<loglen-1>(i);
		data[reordered].re	= in[(i<<1)];
		data[reordered].im	= in[(i<<1)+1];
	}
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::ReorderRealSamplesInPlace(TFFT in[])
{
	size_t length = 1<<(loglen-1);
	for (size_t i=0;i<length;i++)
	{
		size_t reordered	= bitreversal<loglen-1>(i);
		if (i>=reordered)
		{
			TFFT dummy;
			dummy = in[reordered];
			in[reordered] = in[i];
			in[i] = dummy;
		}
	}
}

template <size_t loglen, typename TW, typename T, typename TFFT>
TW CFFT<loglen, TW, T, TFFT>::GetCosFactor(size_t stage, size_t stageshift)
{
	const size_t len = 1<<loglen;
	//more than .5 is not possible
	//size_t k = len*stage/stagesize;
	size_t k = (len*stage)>>stageshift;
	return ma_TwiddleFactors[k];
}

template <size_t loglen, typename TW, typename T, typename TFFT>
TW CFFT<loglen, TW, T, TFFT>::GetSinFactor(size_t stage, size_t stageshift)
{
	const size_t len        = 1<<loglen;
	//const size_t len2        = 1<<(loglen-1);
	const size_t len4        = 1<<(loglen-2);
	const size_t mask        = ~(((-1)>>(loglen))<<(loglen-1));
	//more than .5 is not possible
	//size_t k = (len*stage/stagesize+(len>>2))&mask;
	size_t k = (((len*stage)>>stageshift)+(len>>2))&mask;
	return k<len4?ma_TwiddleFactors[k]:-ma_TwiddleFactors[k];
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::CalculateFFT(T buffer[], TFFT ffttransform[])
{
	static const size_t bitcorrection = sizeof(reinterpret_cast<TFFT*>(0)->re) * 8 - 1;
	static const size_t wscaling_real = (sizeof(reinterpret_cast<TFFT*>(0)->re) - sizeof(TW)) * 8;
	static const size_t wscaling_imag = (sizeof(reinterpret_cast<TFFT*>(0)->im) - sizeof(TW)) * 8;

	size_t ldlen	= loglen;
	size_t length	= 1<<loglen;
	ReorderSamples(buffer, ffttransform);

	for(size_t i=0; i<ldlen;i++)
	{
		size_t stagesize	= static_cast<size_t>(1)<<(i+1);
		size_t stagesize2	= static_cast<size_t>(1)<<i;
		size_t looplength	= length>>(i+1);
		for (size_t R = 0; R < stagesize2; R++)
		{
			TFFT weight;
			weight.re = GetCosFactor(R, (i + 1));
			weight.im = -GetSinFactor(R, (i + 1));//-sin(PI2*RT_ratio);

			//optimized by the compiler if sizeof(TW)==sizeof(TFFT::ELEMENT)
			weight.re <<= wscaling_real;
			weight.im <<= wscaling_imag;
			

			for (size_t L = 0; L < looplength; L++)
			{
				size_t butterflyix = L << (i + 1);
				TFFT temp_a	=	ffttransform[butterflyix+R];
				TFFT temp_b =	ffttransform[butterflyix+R+stagesize2];//load individual stages
#ifdef _C1x_
				auto dummy = weight*temp_b;
				ffttransform[butterflyix+R]				= static_cast<TFFT>((temp_a<<bitcorrection)+dummy);
				ffttransform[butterflyix+R+stagesize2]	= static_cast<TFFT>((temp_a<<bitcorrection)-dummy);
#else
				typename TFFT::ACCTYPE dummy			= weight*temp_b;
				ffttransform[butterflyix+R]				= static_cast<TFFT>((temp_a<<bitcorrection)+dummy);
				ffttransform[butterflyix+R+stagesize2]	= static_cast<TFFT>((temp_a<<bitcorrection)-dummy);
#endif
			}
		}
	}
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::CalculateIFFT(TFFT buffer[], TFFT ffttransform[])
{
	static const size_t bitcorrection = sizeof(reinterpret_cast<TFFT*>(0)->re) * 8 - 1;
	static const size_t wscaling_real = (sizeof(reinterpret_cast<TFFT*>(0)->re) - sizeof(TW)) * 8;
	static const size_t wscaling_imag = (sizeof(reinterpret_cast<TFFT*>(0)->im) - sizeof(TW)) * 8;
	size_t ldlen	= loglen;
	size_t length	= 1<<loglen;
	ReorderSamples(buffer,ffttransform);

	for(size_t i=0; i<ldlen;i++)
	{
		size_t stagesize	= static_cast<size_t>(1)<<(i+1);
		size_t stagesize2	= static_cast<size_t>(1)<<i;
		size_t looplength	= length>>(i+1);
		for (size_t R = 0; R < stagesize2; R++)
		{
			TFFT weight;
			weight.re = GetCosFactor(R, (i + 1));
			weight.im = GetSinFactor(R, (i + 1));//-sin(PI2*RT_ratio);

			//optimized by the compiler if sizeof(TW)==sizeof(TFFT::ELEMENT)
			weight.re <<= wscaling_real;
			weight.im <<= wscaling_imag;

			for (size_t L = 0; L < looplength; L++)
			{
				size_t butterflyix = L << (i + 1);
				TFFT temp_a	=	ffttransform[butterflyix+R];
				TFFT temp_b =	ffttransform[butterflyix+R+stagesize2];//load individual stages	
#ifdef _C1x_
				auto dummy = weight*temp_b;
				ffttransform[butterflyix+R]				= static_cast<TFFT>((temp_a<<bitcorrection)+dummy);
				ffttransform[butterflyix+R+stagesize2]	= static_cast<TFFT>((temp_a<<bitcorrection)-dummy);
#else
				typename TFFT::ACCTYPE dummy			= weight*temp_b;
				ffttransform[butterflyix+R]				= static_cast<TFFT>((temp_a<<bitcorrection)+dummy);
				ffttransform[butterflyix+R+stagesize2]	= static_cast<TFFT>((temp_a<<bitcorrection)-dummy);
#endif
			}
		}
	}
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::CalculateFFT(T buffer[], TFFT ffttransform[], EN_ScalingMethod rescale, size_t *pBlockExponent)
{
	static const size_t bitcorrection = sizeof(reinterpret_cast<TFFT*>(0)->re) * 8 - 1;
	static const size_t wscaling_real = (sizeof(reinterpret_cast<TFFT*>(0)->re) - sizeof(TW)) * 8;
	static const size_t wscaling_imag = (sizeof(reinterpret_cast<TFFT*>(0)->im) - sizeof(TW)) * 8;

	size_t ldlen	= loglen;
	size_t length	= 1<<loglen;
	size_t Exponent	= 0;
	*pBlockExponent = loglen;
	ReorderSamples(buffer, ffttransform);

	if (rescale==SCALE_INPUT||rescale==SCALEINPOUT)
	{
		T MaxVal		= abs(buffer[0]);
		for (size_t i=1;i<length;i++)
		{
			T dummy = abs(buffer[i]);
			if (dummy>MaxVal)
			{
				MaxVal=dummy;
			}
		}

		if (MaxVal!=0)
		{
			T mask			=	0;
			mask			= (~mask)<<(sizeof(T)*8-2);
			for (Exponent=0;!(MaxVal&mask);MaxVal<<=1){Exponent++;}
			Exponent		+= (sizeof(reinterpret_cast<TFFT*>(0)->re)-sizeof(T))*8;

			//scale to the max is done in the butterfly
			/*for (size_t i=0;i<length;i++)
			{
				ffttransform[i].re <<= Exponent;
			}
			*pBlockExponent = Exponent;*/
		}
	}

	
	for(size_t i=0; i<ldlen;i++)
	{
		size_t stagesize	= 1<<(i+1);
		size_t stagesize2	= 1<<i;
		size_t looplength	= length>>(i+1);
		typename TFFT::DATATYPE MaxOccVal = 0;
		for (size_t R = 0; R < stagesize2; R++)
		{
			TFFT weight;
			weight.re = GetCosFactor(R, (i + 1));
			weight.im = -GetSinFactor(R, (i + 1));//-sin(PI2*RT_ratio);

			//optimized by the compiler if sizeof(TW)==sizeof(TFFT::ELEMENT)
			weight.re <<= wscaling_real;
			weight.im <<= wscaling_imag;

			for (size_t L = 0; L < looplength; L++)
			{
				size_t butterflyix = L << (i + 1);
				TFFT temp_a	=	ffttransform[butterflyix+R];
				TFFT temp_b =	ffttransform[butterflyix+R+stagesize2];//load individual stages
				temp_a.re	<<= Exponent;
				temp_a.im	<<= Exponent;
				temp_b.re	<<= Exponent;
				temp_b.im	<<= Exponent;
#ifdef _C1x_
				auto dummy = weight*temp_b;
				ffttransform[butterflyix+R]				= static_cast<TFFT>((temp_a<<bitcorrection)+dummy);
				ffttransform[butterflyix+R+stagesize2]	= static_cast<TFFT>((temp_a<<bitcorrection)-dummy);
#else
				typename TFFT::ACCTYPE dummy			= weight*temp_b;
				ffttransform[butterflyix+R]				= static_cast<TFFT>((temp_a<<bitcorrection)+dummy);
				ffttransform[butterflyix+R+stagesize2]	= static_cast<TFFT>((temp_a<<bitcorrection)-dummy);
#endif
				MaxOccVal = MaxOccVal>=abs(ffttransform[butterflyix+R].re)?MaxOccVal:abs(ffttransform[butterflyix+R].re);
				MaxOccVal = MaxOccVal>=abs(ffttransform[butterflyix+R+stagesize2].re)?MaxOccVal:abs(ffttransform[butterflyix+R+stagesize2].re);
				MaxOccVal = MaxOccVal>=abs(ffttransform[butterflyix+R].im)?MaxOccVal:abs(ffttransform[butterflyix+R].im);
				MaxOccVal = MaxOccVal>=abs(ffttransform[butterflyix+R+stagesize2].im)?MaxOccVal:abs(ffttransform[butterflyix+R+stagesize2].im);

			}
		}
		*pBlockExponent -= Exponent;


		if (MaxOccVal != 0)
		{
			typename TFFT::DATATYPE mask = 0;
			mask = (~mask) << (sizeof(TFFT::DATATYPE) * 8 - 2);
			for (Exponent = 0; !(MaxOccVal & mask); MaxOccVal <<= 1) { Exponent++; }
		}
		MaxOccVal	= 0;
	}
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::CalculateIFFT(TFFT buffer[], TFFT ffttransform[], EN_ScalingMethod rescale, size_t *pBlockExponent)
{
	static const size_t bitcorrection = sizeof(reinterpret_cast<TFFT*>(0)->re) * 8 - 1;
	static const size_t wscaling_real = (sizeof(reinterpret_cast<TFFT*>(0)->re) - sizeof(TW)) * 8;
	static const size_t wscaling_imag = (sizeof(reinterpret_cast<TFFT*>(0)->im) - sizeof(TW)) * 8;

	size_t ldlen	= loglen;
	size_t length	= 1<<loglen;
	size_t Exponent	= 0;
	*pBlockExponent = loglen;
	ReorderSamples(buffer, ffttransform);

	if (rescale==SCALE_INPUT||rescale==SCALEINPOUT)
	{
		typename TFFT::DATATYPE MaxVal1		= abs(buffer[0].re);
		typename TFFT::DATATYPE MaxVal2		= abs(buffer[0].im);
		typename TFFT::DATATYPE MaxVal		= MaxVal1>MaxVal2?MaxVal1:MaxVal2;
		for (size_t i=1;i<length;i++)
		{
			typename TFFT::DATATYPE dummy1 = abs(buffer[i].re);
			typename TFFT::DATATYPE dummy2 = abs(buffer[i].im);
			typename TFFT::DATATYPE dummy  = dummy1>dummy2?dummy1:dummy2;
			if (dummy>MaxVal)
			{
				MaxVal=dummy;
			}
		}

		if (MaxVal!=0)
		{
			T mask			=	0;
			mask			= (~mask)<<(sizeof(T)*8-2);
			for (Exponent=0;!(MaxVal&mask);MaxVal<<=1){Exponent++;}
			Exponent		+= (sizeof(reinterpret_cast<TFFT*>(0)->re)-sizeof(T))*8;

			//scale to the max is done in the butterfly
			/*for (size_t i=0;i<length;i++)
			{
				ffttransform[i].re <<= Exponent;
			}
			*pBlockExponent = Exponent;*/
		}
	}

	
	for(size_t i=0; i<ldlen;i++)
	{
		size_t stagesize	= 1<<(i+1);
		size_t stagesize2	= 1<<i;
		size_t looplength	= length>>(i+1);
		typename TFFT::DATATYPE MaxOccVal = 0;
		for (size_t R = 0; R < stagesize2; R++)
		{
			TFFT weight;
			weight.re = GetCosFactor(R, (i + 1));
			weight.im = -GetSinFactor(R, (i + 1));//-sin(PI2*RT_ratio);

			//optimized by the compiler if sizeof(TW)==sizeof(TFFT::ELEMENT)
			weight.re <<= wscaling_real;
			weight.im <<= wscaling_imag;

			for (size_t L = 0; L < looplength; L++)
			{
				size_t butterflyix = L << (i + 1);
				TFFT temp_a	=	ffttransform[butterflyix+R];
				TFFT temp_b =	ffttransform[butterflyix+R+stagesize2];//load individual stages
				temp_a.re	<<= Exponent;
				temp_a.im	<<= Exponent;
				temp_b.re	<<= Exponent;
				temp_b.im	<<= Exponent;
				
#ifdef _C1x_
				auto dummy = weight*temp_b;
				ffttransform[butterflyix+R]				= static_cast<TFFT>((temp_a<<bitcorrection)+dummy);
				ffttransform[butterflyix+R+stagesize2]	= static_cast<TFFT>((temp_a<<bitcorrection)-dummy);
#else
				typename TFFT::ACCTYPE dummy			= weight*temp_b;
				ffttransform[butterflyix+R]				= static_cast<TFFT>((temp_a<<bitcorrection)+dummy);
				ffttransform[butterflyix+R+stagesize2]	= static_cast<TFFT>((temp_a<<bitcorrection)-dummy);
#endif
				MaxOccVal = MaxOccVal>=abs(ffttransform[butterflyix+R].re)?MaxOccVal:abs(ffttransform[butterflyix+R].re);
				MaxOccVal = MaxOccVal>=abs(ffttransform[butterflyix+R+stagesize2].re)?MaxOccVal:abs(ffttransform[butterflyix+R+stagesize2].re);
				MaxOccVal = MaxOccVal>=abs(ffttransform[butterflyix+R].im)?MaxOccVal:abs(ffttransform[butterflyix+R].im);
				MaxOccVal = MaxOccVal>=abs(ffttransform[butterflyix+R+stagesize2].im)?MaxOccVal:abs(ffttransform[butterflyix+R+stagesize2].im);

			}
		}
		*pBlockExponent -= Exponent;
		if (MaxOccVal != 0)
		{
			typename TFFT::DATATYPE mask = 0;
			mask = (~mask) << (sizeof(TFFT::DATATYPE) * 8 - 2);
			for (Exponent = 0; !(MaxOccVal & mask); MaxOccVal <<= 1) { Exponent++; }
		}
	}
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::CalculateRealFFT(T buffer[], TFFT ffttransform[])
{
	static const size_t bitcorrection = sizeof(reinterpret_cast<TFFT*>(0)->re) * 8 - 1;
	static const size_t wscaling_real = (sizeof(reinterpret_cast<TFFT*>(0)->re) - sizeof(TW)) * 8;
	static const size_t wscaling_imag = (sizeof(reinterpret_cast<TFFT*>(0)->im) - sizeof(TW)) * 8;

	size_t ldlen	= loglen-1;
	size_t length	= 1<<(loglen-1);
	ReorderRealSamples(buffer, ffttransform);

	for(size_t i=0; i<ldlen;i++)
	{
		size_t stagesize	= 1<<(i+1);
		size_t stagesize2	= 1<<i;
		size_t looplength	= length>>(i+1);
		for (size_t R = 0; R < stagesize2; R++)
		{
			TFFT weight;
			weight.re = GetCosFactor(R, (i + 1));
			weight.im = -GetSinFactor(R, (i + 1));//-sin(PI2*RT_ratio);

			//optimized by the compiler if sizeof(TW)==sizeof(TFFT::ELEMENT)
			weight.re <<= wscaling_real;
			weight.im <<= wscaling_imag;

			for (size_t L = 0; L < looplength; L++)
			{
				size_t butterflyix = L << (i + 1);
				TFFT temp_a	=	ffttransform[butterflyix+R];
				TFFT temp_b =	ffttransform[butterflyix+R+stagesize2];//load individual stages
#ifdef _C1x_
				auto dummy = weight*temp_b;
				ffttransform[butterflyix+R]				= static_cast<TFFT>((temp_a<<bitcorrection)+dummy);
				ffttransform[butterflyix+R+stagesize2]	= static_cast<TFFT>((temp_a<<bitcorrection)-dummy);
#else
				typename TFFT::ACCTYPE dummy			= weight*temp_b;
				ffttransform[butterflyix+R]				= static_cast<TFFT>((temp_a<<bitcorrection)+dummy);
				ffttransform[butterflyix+R+stagesize2]	= static_cast<TFFT>((temp_a<<bitcorrection)-dummy);
#endif
			}
		}
	}
	Convert2HalfDFT(ffttransform);
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::CalculateRealIFFT(TFFT buffer[], TFFT ffttransform[])
{
	static const size_t bitcorrection = sizeof(reinterpret_cast<TFFT*>(0)->re) * 8 - 1;
	static const size_t wscaling_real = (sizeof(reinterpret_cast<TFFT*>(0)->re) - sizeof(TW)) * 8;
	static const size_t wscaling_imag = (sizeof(reinterpret_cast<TFFT*>(0)->im) - sizeof(TW)) * 8;

	size_t ldlen	= loglen-1;
	size_t length	= 1<<(loglen-1);
	
	memcpy(ffttransform,buffer,length*sizeof(TFFT));
	Convert2ComplexDFT(ffttransform);
	ReorderRealSamplesInPlace(ffttransform);
	//ReorderRealSamples(buffer, ffttransform);


	for(size_t i=0; i<ldlen;i++)
	{
		size_t stagesize	= 1<<(i+1);
		size_t stagesize2	= 1<<i;
		size_t looplength	= length>>(i+1);
		for (size_t R = 0; R < stagesize2; R++)
		{
			TFFT weight;
			weight.re = GetCosFactor(R, (i + 1));
			weight.im = GetSinFactor(R, (i + 1));//-sin(PI2*RT_ratio);

			//optimized by the compiler if sizeof(TW)==sizeof(TFFT::ELEMENT)
			weight.re <<= wscaling_real;
			weight.im <<= wscaling_imag;

			for (size_t L = 0; L < looplength; L++)
			{
				size_t butterflyix = L << (i + 1);
				TFFT temp_a	=	ffttransform[butterflyix+R];
				TFFT temp_b =	ffttransform[butterflyix+R+stagesize2];//load individual stages
#ifdef _C1x_
				auto dummy = weight*temp_b;
				ffttransform[butterflyix+R]				= static_cast<TFFT>((temp_a<<bitcorrection)+dummy);
				ffttransform[butterflyix+R+stagesize2]	= static_cast<TFFT>((temp_a<<bitcorrection)-dummy);
#else
				typename TFFT::ACCTYPE dummy			= weight*temp_b;
				ffttransform[butterflyix+R]				= static_cast<TFFT>((temp_a<<bitcorrection)+dummy);
				ffttransform[butterflyix+R+stagesize2]	= static_cast<TFFT>((temp_a<<bitcorrection)-dummy);
#endif
			}
		}
	}
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::CalculateRealFFT(T buffer[], TFFT ffttransform[], EN_ScalingMethod rescale, size_t *pBlockExponent)
{
	static const size_t bitcorrection = sizeof(reinterpret_cast<TFFT*>(0)->re) * 8 - 1;
	static const size_t wscaling_real = (sizeof(reinterpret_cast<TFFT*>(0)->re) - sizeof(TW)) * 8;
	static const size_t wscaling_imag = (sizeof(reinterpret_cast<TFFT*>(0)->im) - sizeof(TW)) * 8;

	size_t ldlen	= loglen-1;
	size_t length	= 1<<(loglen-1);
	size_t Exponent	= 0;
	*pBlockExponent = loglen;
	ReorderRealSamples(buffer, ffttransform);

	if (rescale==SCALE_INPUT||rescale==SCALEINPOUT)
	{
		T MaxVal		= abs(buffer[0]);
		for (size_t i=1;i<length;i++)
		{
			T dummy = abs(buffer[i]);
			if (dummy>MaxVal)
			{
				MaxVal=dummy;
			}
		}

		if (MaxVal!=0)
		{
			T mask			=	0;
			mask			= (~mask)<<(sizeof(T)*8-2);
			for (Exponent=0;!(MaxVal&mask);MaxVal<<=1){Exponent++;}
			Exponent		+= (sizeof(reinterpret_cast<TFFT*>(0)->re)-sizeof(T))*8;

			//scale to the max is done in the butterfly
			/*for (size_t i=0;i<length;i++)
			{
				ffttransform[i].re <<= Exponent;
			}
			*pBlockExponent = Exponent;*/
		}
	}

	
	for(size_t i=0; i<ldlen;i++)
	{
		//size_t stagesize	= static_cast<size_t>(1)<<(i+1);
		size_t stagesize2	= static_cast<size_t>(1)<<i;
		size_t looplength	= length>>(i+1);
		typename TFFT::DATATYPE MaxOccVal = 0;
		for (size_t R = 0; R < stagesize2; R++)
		{
			TFFT weight;
			weight.re = GetCosFactor(R, (i + 1));
			weight.im = -GetSinFactor(R, (i + 1));//-sin(PI2*RT_ratio);

			//optimized by the compiler if sizeof(TW)==sizeof(TFFT::ELEMENT)
			weight.re <<= wscaling_real;
			weight.im <<= wscaling_imag;

			for (size_t L = 0; L < looplength; L++)
			{
				size_t butterflyix = L << (i + 1);
				TFFT temp_a	=	ffttransform[butterflyix+R];
				TFFT temp_b =	ffttransform[butterflyix+R+stagesize2];//load individual stages
				temp_a.re	<<= Exponent;
				temp_a.im	<<= Exponent;
				temp_b.re	<<= Exponent;
				temp_b.im	<<= Exponent;
				
#ifdef _C1x_
				auto dummy = weight*temp_b;
				ffttransform[butterflyix+R]				= static_cast<TFFT>((temp_a<<bitcorrection)+dummy);
				ffttransform[butterflyix+R+stagesize2]	= static_cast<TFFT>((temp_a<<bitcorrection)-dummy);
#else
				typename TFFT::ACCTYPE dummy			= weight*temp_b;
				ffttransform[butterflyix+R]				= static_cast<TFFT>((temp_a<<bitcorrection)+dummy);
				ffttransform[butterflyix+R+stagesize2]	= static_cast<TFFT>((temp_a<<bitcorrection)-dummy);
#endif
				MaxOccVal = MaxOccVal>=abs(ffttransform[butterflyix+R].re)?MaxOccVal:abs(ffttransform[butterflyix+R].re);
				MaxOccVal = MaxOccVal>=abs(ffttransform[butterflyix+R+stagesize2].re)?MaxOccVal:abs(ffttransform[butterflyix+R+stagesize2].re);
				MaxOccVal = MaxOccVal>=abs(ffttransform[butterflyix+R].im)?MaxOccVal:abs(ffttransform[butterflyix+R].im);
				MaxOccVal = MaxOccVal>=abs(ffttransform[butterflyix+R+stagesize2].im)?MaxOccVal:abs(ffttransform[butterflyix+R+stagesize2].im);

			}
		}
		*pBlockExponent -= Exponent;
		if (MaxOccVal != 0)
		{
			typename TFFT::DATATYPE mask = 0;
			mask = (~mask) << (sizeof(TFFT::DATATYPE) * 8 - 2);
			for (Exponent = 0; !(MaxOccVal & mask); MaxOccVal <<= 1) { Exponent++; }
		}
	}
	Convert2HalfDFT(ffttransform, Exponent);
	*pBlockExponent -= Exponent;
}


template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::CalculateRealIFFT(TFFT buffer[], TFFT ffttransform[], EN_ScalingMethod rescale, size_t *pBlockExponent)
{
	static const size_t bitcorrection = sizeof(reinterpret_cast<TFFT*>(0)->re) * 8 - 1;
	static const size_t wscaling_real = (sizeof(reinterpret_cast<TFFT*>(0)->re) - sizeof(TW)) * 8;
	static const size_t wscaling_imag = (sizeof(reinterpret_cast<TFFT*>(0)->im) - sizeof(TW)) * 8;

	size_t ldlen	= loglen-1;
	size_t length	= 1<<(loglen-1);
	size_t Exponent	= 0;
	*pBlockExponent = loglen;
	
	memcpy(ffttransform,buffer,length*sizeof(TFFT));


	if (rescale==SCALE_INPUT||rescale==SCALEINPOUT)
	{
		typename TFFT::DATATYPE MaxVal1		= abs(ffttransform[0].re);
		typename TFFT::DATATYPE MaxVal2		= abs(ffttransform[0].im);
		typename TFFT::DATATYPE MaxVal		= MaxVal1<MaxVal2?MaxVal2:MaxVal1;
		for (size_t i=1;i<length;i++)
		{
			typename TFFT::DATATYPE MaxVal1	= abs(ffttransform[i].re);
			typename TFFT::DATATYPE MaxVal2	= abs(ffttransform[i].im);
			typename TFFT::DATATYPE dummy 	= MaxVal1<MaxVal2?MaxVal2:MaxVal1;
			if (dummy>MaxVal)
			{
				MaxVal=dummy;
			}
		}

		if (MaxVal!=0)
		{
			typename TFFT::DATATYPE mask  =	0;
			mask					= (~mask)<<(sizeof(T)*8-2);
			for (Exponent=0;!(MaxVal&mask);MaxVal<<=1){Exponent++;}
			//Exponent				+= (sizeof(reinterpret_cast<TFFT*>(0)->re)-sizeof(T))*8;

			//scale to the max is done in the butterfly
			/*for (size_t i=0;i<length;i++)
			{
				ffttransform[i].re <<= Exponent;
			}
			*pBlockExponent = Exponent;*/
		}
	}
	Convert2ComplexDFT(ffttransform,Exponent);
	ReorderRealSamplesInPlace(ffttransform);
	*pBlockExponent -= Exponent;

	if (rescale==SCALEINPOUT)
	{
		typename TFFT::DATATYPE MaxVal1		= abs(ffttransform[0].re);
		typename TFFT::DATATYPE MaxVal2		= abs(ffttransform[0].im);
		typename TFFT::DATATYPE MaxVal		= MaxVal1<MaxVal2?MaxVal2:MaxVal1;
		for (size_t i=1;i<length;i++)
		{
			typename TFFT::DATATYPE MaxVal1	= abs(ffttransform[i].re);
			typename TFFT::DATATYPE MaxVal2	= abs(ffttransform[i].im);
			typename TFFT::DATATYPE dummy	= MaxVal1<MaxVal2?MaxVal2:MaxVal1;
			if (dummy>MaxVal)
			{
				MaxVal=dummy;
			}
		}

		if (MaxVal!=0)
		{
			typename TFFT::DATATYPE mask		=	0;
			mask					= (~mask)<<(sizeof(T)*8-2);
			for (Exponent=0;!(MaxVal&mask);MaxVal<<=1){Exponent++;}
			//Exponent				+= (sizeof(reinterpret_cast<TFFT*>(0)->re)-sizeof(T))*8;

			//scale to the max is done in the butterfly
			/*for (size_t i=0;i<length;i++)
			{
			ffttransform[i].re <<= Exponent;
			}
			*pBlockExponent = Exponent;*/
		}
	}
	else
	{
		Exponent=0;
	}

	
	for (size_t i = 0; i < ldlen; i++)
	{
		//size_t stagesize	= static_cast<size_t>(1)<<(i+1);
		size_t stagesize2 = static_cast<size_t>(1) << i;
		size_t looplength = length >> (i + 1);
		typename TFFT::DATATYPE MaxOccVal = 0;
		for (size_t R = 0; R < stagesize2; R++)
		{
			TFFT weight;
			weight.re = GetCosFactor(R, (i + 1));
			weight.im = GetSinFactor(R, (i + 1));//-sin(PI2*RT_ratio);

			//optimized by the compiler if sizeof(TW)==sizeof(TFFT::ELEMENT)
			weight.re <<= wscaling_real;
			weight.im <<= wscaling_imag;

			for (size_t L = 0; L < looplength; L++)
			{
				size_t butterflyix = L << (i + 1);
				TFFT temp_a = ffttransform[butterflyix + R];
				TFFT temp_b = ffttransform[butterflyix + R + stagesize2];//load individual stages
				temp_a.re <<= Exponent;
				temp_a.im <<= Exponent;
				temp_b.re <<= Exponent;
				temp_b.im <<= Exponent;

#ifdef _C1x_
				auto dummy = weight * temp_b;
				ffttransform[butterflyix + R] = static_cast<TFFT>((temp_a << bitcorrection) + dummy);
				ffttransform[butterflyix + R + stagesize2] = static_cast<TFFT>((temp_a << bitcorrection) - dummy);
#else
				typename TFFT::ACCTYPE dummy = weight * temp_b;
				ffttransform[butterflyix + R] = static_cast<TFFT>((temp_a << bitcorrection) + dummy);
				ffttransform[butterflyix + R + stagesize2] = static_cast<TFFT>((temp_a << bitcorrection) - dummy);
#endif
				MaxOccVal = MaxOccVal >= abs(ffttransform[butterflyix + R].re) ? MaxOccVal : abs(ffttransform[butterflyix + R].re);
				MaxOccVal = MaxOccVal >= abs(ffttransform[butterflyix + R + stagesize2].re) ? MaxOccVal : abs(ffttransform[butterflyix + R + stagesize2].re);
				MaxOccVal = MaxOccVal >= abs(ffttransform[butterflyix + R].im) ? MaxOccVal : abs(ffttransform[butterflyix + R].im);
				MaxOccVal = MaxOccVal >= abs(ffttransform[butterflyix + R + stagesize2].im) ? MaxOccVal : abs(ffttransform[butterflyix + R + stagesize2].im);

			}
		}
		*pBlockExponent -= Exponent;

		if (MaxOccVal != 0)
		{
			typename TFFT::DATATYPE mask = 0;
			mask = (~mask) << (sizeof(TFFT::DATATYPE) * 8 - 2);
			for (Exponent = 0; !(MaxOccVal & mask); MaxOccVal <<= 1) { Exponent++; }
		}
	}
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::Convert2HalfDFT(TFFT transform[])
{
	static const float PI  = 3.1415926535897932384626434f;

	size_t length		= 1<<(loglen-1);
	size_t redlength	= length;
	size_t halflength	= length>>1;
	for (size_t i=1;i<halflength;i++)
	{
		//here the FB case
		typename TFFT::OFWTYPE H1R		= (static_cast<typename TFFT::OFWTYPE>(transform[i].re)+transform[redlength-i].re)>>1;
		typename TFFT::OFWTYPE H1I		= (static_cast<typename TFFT::OFWTYPE>(transform[i].im)-transform[redlength-i].im)>>1;//because conj complex
		typename TFFT::OFWTYPE H2R		= (static_cast<typename TFFT::OFWTYPE>(transform[i].re)-transform[redlength-i].re)>>1;//that is the second term
		typename TFFT::OFWTYPE H2I		= (static_cast<typename TFFT::OFWTYPE>(transform[i].im)+transform[redlength-i].im)>>1;//adding, because conj complex

		//to be checked

		TW cosine						= GetCosFactor(i,loglen);//cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		TW sine							= GetSinFactor(i,loglen);//sin(PI*i/length);
		
		typename TFFT::OFWTYPE res1i	= (H2R*cosine+H2I*sine)>>(sizeof(TW)*8-1);
		typename TFFT::OFWTYPE res1r	= -((H2I*cosine-H2R*sine)>>(sizeof(TW)*8-1));//complex i is integrated

		typename TFFT::OFWTYPE resreal	= (H1R-res1r);//0.5f*(H1R-res1r);
		typename TFFT::OFWTYPE resimag	= (H1I-res1i);//0.5f*(H1I-res1i);->already scaled see shift operation above

		//here follows the FN/2-b case (copy from above)
		//->HR1 is still valid
		//double H1R		= transform[i].re+transform[redlength-i].re;
		//H1I is just inverted
		//double H1I		= -transform[i].im-transform[redlength-i].im;//because conj complex

		//use HR1 as above, invert H1I
		//use HR2 as above, invert H2I (formula on paper contain an error)
		//double H2R		= transform[i].re-transform[redlength-i].re;//that is the second term
		//double H2I		= transform[i].im+transform[redlength-i].im;//adding, because conj complex

		//double cosine	= cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		//double sine		= sin(PI*i/length);
		
		res1i							= (H2R*cosine+H2I*sine)>>(sizeof(TW)*8-1);
		res1r							= (H2I*cosine-H2R*sine)>>(sizeof(TW)*8-1);//complex i is integrated

		typename TFFT::OFWTYPE resreal2	= (H1R-res1r);//0.5f*(H1R-res1r);
		typename TFFT::OFWTYPE resimag2	= (-H1I-res1i);//0.5f*(-H1I-res1i);->already scaled see shift operation above

		transform[i].re					= static_cast<typename TFFT::DATATYPE>(resreal>>1);
		transform[i].im					= static_cast<typename TFFT::DATATYPE>(resimag>>1);
		transform[length-i].re			= static_cast<typename TFFT::DATATYPE>(resreal2>>1);
		transform[length-i].im			= static_cast<typename TFFT::DATATYPE>(resimag2>>1);
	}
	typename TFFT::OFWTYPE lval			= static_cast<typename TFFT::OFWTYPE>(transform[0].re)+transform[0].im;
	typename TFFT::OFWTYPE hval			= static_cast<typename TFFT::OFWTYPE>(transform[0].re)-transform[0].im;
	transform[0].re						= static_cast<typename TFFT::DATATYPE>(lval>>1);
	transform[0].im						= static_cast<typename TFFT::DATATYPE>(hval>>1);
	transform[halflength].im			= -transform[halflength].im;//treatment of the half index
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::Convert2ComplexDFT(TFFT transform[])
{
	static const float PI2 = 6.2831853071795864769252868f;
	static const float PI  = 3.1415926535897932384626434f;

	size_t length		= 1<<(loglen-1);
	size_t redlength	= length;
	size_t halflength	= length>>1;
	for (size_t i=1;i<halflength;i++)
	{
		//here the FB case
		typename TFFT::OFWTYPE H1R		= (static_cast<typename TFFT::OFWTYPE>(transform[i].re)+transform[redlength-i].re)>>1;
		typename TFFT::OFWTYPE H1I		= (static_cast<typename TFFT::OFWTYPE>(transform[i].im)-transform[redlength-i].im)>>1;//because conj complex
		typename TFFT::OFWTYPE H2R		= (static_cast<typename TFFT::OFWTYPE>(transform[i].re)-transform[redlength-i].re)>>1;//that is the second term
		typename TFFT::OFWTYPE H2I		= (static_cast<typename TFFT::OFWTYPE>(transform[i].im)+transform[redlength-i].im)>>1;//adding, because conj complex

		TW cosine						= GetCosFactor(i,loglen);//cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		TW sine							= GetSinFactor(i,loglen);//sin(PI*i/length);
		
		//float res1i	= H2R*cosine-H2I*sine;
		//float res1r	= -(H2I*cosine+H2R*sine);//complex i is integrated
		typename TFFT::OFWTYPE res1i	= (H2R*cosine-H2I*sine)>>(sizeof(TW)*8-1);
		typename TFFT::OFWTYPE res1r	= -((H2I*cosine+H2R*sine)>>(sizeof(TW)*8-1));//complex i is integrated


		//float resreal	= 0.5f*(H1R+res1r);
		//float resimag	= 0.5f*(H1I+res1i);
		typename TFFT::OFWTYPE resreal	= (H1R+res1r);//0.5f*(H1R+res1r);
		typename TFFT::OFWTYPE resimag	= (H1I+res1i);//0.5f*(H1I+res1i);->already scaled see shift operation above

		//here follows the FN/2-b case (copy from above)
		//->HR1 is still valid
		//double H1R		= transform[i].re+transform[redlength-i].re;
		//H1I is just inverted
		//double H1I		= -transform[i].im-transform[redlength-i].im;//because conj complex

		//use HR1 as above, invert H1I
		//use HR2 as above, invert H2I (formula on paper contain an error)
		//double H2R		= transform[i].re-transform[redlength-i].re;//that is the second term
		//double H2I		= transform[i].im+transform[redlength-i].im;//adding, because conj complex

		//double cosine	= cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		//double sine		= sin(PI*i/length);
		
		//res1i					= -H2R*cosine+H2I*sine;
		//res1r					= (H2I*cosine+H2R*sine);//complex i is integrated
		res1i					= (-H2R*cosine+H2I*sine)>>(sizeof(TW)*8-1);
		res1r					= (H2I*cosine+H2R*sine)>>(sizeof(TW)*8-1);//complex i is integrated

		//float resreal2			= 0.5f*(H1R+res1r);
		//float resimag2			= -0.5f*(H1I+res1i);
		typename TFFT::OFWTYPE resreal2	= (H1R+res1r);//0.5f*(H1R+res1r);
		typename TFFT::OFWTYPE resimag2	= -(H1I+res1i);//-0.5f*(H1I+res1i);->already scaled see shift operation above

		//transform[i].re			= resreal;
		//transform[i].im			= resimag;
		//transform[length-i].re	= resreal2;
		//transform[length-i].im	= resimag2;
		transform[i].re				= static_cast<typename TFFT::DATATYPE>(resreal);//>>1
		transform[i].im				= static_cast<typename TFFT::DATATYPE>(resimag);//>>1
		transform[length-i].re		= static_cast<typename TFFT::DATATYPE>(resreal2);//>>1
		transform[length-i].im		= static_cast<typename TFFT::DATATYPE>(resimag2);//>>1
	}
	//float lval					= transform[0].re+transform[0].im;
	//float hval					= transform[0].re-transform[0].im;
	typename TFFT::OFWTYPE lval		= static_cast<typename TFFT::OFWTYPE>(transform[0].re)+transform[0].im;
	typename TFFT::OFWTYPE hval		= static_cast<typename TFFT::OFWTYPE>(transform[0].re)-transform[0].im;
	//transform[0].re				= 0.5f*lval;
	//transform[0].im				= 0.5f*hval;
	transform[0].re					= static_cast<typename TFFT::DATATYPE>(lval>>1);
	transform[0].im					= static_cast<typename TFFT::DATATYPE>(hval>>1);
	transform[halflength].im		= -transform[halflength].im;//treatment of the half index
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::Convert2ComplexDFT(TFFT transform[], size_t exponent)
{
	static const float PI2 = 6.2831853071795864769252868f;
	static const float PI  = 3.1415926535897932384626434f;

	size_t length		= 1<<(loglen-1);
	size_t redlength	= length;
	size_t halflength	= length>>1;
	for (size_t i=1;i<halflength;i++)
	{
		//here the FB case
		typename TFFT::OFWTYPE H1R		= ((static_cast<typename TFFT::OFWTYPE>(transform[i].re)+transform[redlength-i].re)<<exponent)>>1;
		typename TFFT::OFWTYPE H1I		= ((static_cast<typename TFFT::OFWTYPE>(transform[i].im)-transform[redlength-i].im)<<exponent)>>1;//because conj complex
		typename TFFT::OFWTYPE H2R		= ((static_cast<typename TFFT::OFWTYPE>(transform[i].re)-transform[redlength-i].re)<<exponent)>>1;//that is the second term
		typename TFFT::OFWTYPE H2I		= ((static_cast<typename TFFT::OFWTYPE>(transform[i].im)+transform[redlength-i].im)<<exponent)>>1;//adding, because conj complex

		TW cosine						= GetCosFactor(i,loglen);//cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		TW sine							= GetSinFactor(i,loglen);//sin(PI*i/length);
		
		//float res1i	= H2R*cosine-H2I*sine;
		//float res1r	= -(H2I*cosine+H2R*sine);//complex i is integrated
		typename TFFT::OFWTYPE res1i	= (H2R*cosine-H2I*sine)>>(sizeof(TW)*8-1);
		typename TFFT::OFWTYPE res1r	= -((H2I*cosine+H2R*sine)>>(sizeof(TW)*8-1));//complex i is integrated


		//float resreal	= 0.5f*(H1R+res1r);
		//float resimag	= 0.5f*(H1I+res1i);
		typename TFFT::OFWTYPE resreal	= (H1R+res1r);//0.5f*(H1R+res1r);
		typename TFFT::OFWTYPE resimag	= (H1I+res1i);//0.5f*(H1I+res1i);->already scaled see shift operation above

		//here follows the FN/2-b case (copy from above)
		//->HR1 is still valid
		//double H1R		= transform[i].re+transform[redlength-i].re;
		//H1I is just inverted
		//double H1I		= -transform[i].im-transform[redlength-i].im;//because conj complex

		//use HR1 as above, invert H1I
		//use HR2 as above, invert H2I (formula on paper contain an error)
		//double H2R		= transform[i].re-transform[redlength-i].re;//that is the second term
		//double H2I		= transform[i].im+transform[redlength-i].im;//adding, because conj complex

		//double cosine	= cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		//double sine		= sin(PI*i/length);
		
		//res1i					= -H2R*cosine+H2I*sine;
		//res1r					= (H2I*cosine+H2R*sine);//complex i is integrated
		res1i					= (-H2R*cosine+H2I*sine)>>(sizeof(TW)*8-1);
		res1r					= (H2I*cosine+H2R*sine)>>(sizeof(TW)*8-1);//complex i is integrated

		//float resreal2			= 0.5f*(H1R+res1r);
		//float resimag2			= -0.5f*(H1I+res1i);
		typename TFFT::OFWTYPE resreal2	= (H1R+res1r);//0.5f*(H1R+res1r);
		typename TFFT::OFWTYPE resimag2	= -(H1I+res1i);//-0.5f*(H1I+res1i);->already scaled see shift operation above

		//transform[i].re			= resreal;
		//transform[i].im			= resimag;
		//transform[length-i].re	= resreal2;
		//transform[length-i].im	= resimag2;
		transform[i].re				= static_cast<typename TFFT::DATATYPE>(resreal);//>>1
		transform[i].im				= static_cast<typename TFFT::DATATYPE>(resimag);//>>1
		transform[length-i].re		= static_cast<typename TFFT::DATATYPE>(resreal2);//>>1
		transform[length-i].im		= static_cast<typename TFFT::DATATYPE>(resimag2);//>>1
	}
	//float lval					= transform[0].re+transform[0].im;
	//float hval					= transform[0].re-transform[0].im;
	typename TFFT::OFWTYPE lval;//		= ((static_cast<typename TFFT::OFWTYPE>(transform[0].re)+transform[0].im)>>1)<<(exponent);
	typename TFFT::OFWTYPE hval;//		= ((static_cast<typename TFFT::OFWTYPE>(transform[0].re)-transform[0].im)>>1)<<(exponent);
	if (exponent)
	{
		lval = ((static_cast<typename TFFT::OFWTYPE>(transform[0].re)+transform[0].im))<<(exponent-1);
		hval = ((static_cast<typename TFFT::OFWTYPE>(transform[0].re)-transform[0].im))<<(exponent-1);
	}
	else
	{
		lval = ((static_cast<typename TFFT::OFWTYPE>(transform[0].re)+transform[0].im))>>1;
		hval = ((static_cast<typename TFFT::OFWTYPE>(transform[0].re)-transform[0].im))>>1;
	}
	//transform[0].re				= 0.5f*lval;
	//transform[0].im				= 0.5f*hval;
	transform[0].re					= static_cast<typename TFFT::DATATYPE>(lval);
	transform[0].im					= static_cast<typename TFFT::DATATYPE>(hval);//>>1 moved to exponent stuff ^
	transform[halflength].im		= -(transform[halflength].im<<exponent);//treatment of the half index
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::Convert2HalfDFT(TFFT transform[], size_t exponent)
{
	static const float PI  = 3.1415926535897932384626434f;

	size_t length		= 1<<(loglen-1);
	size_t redlength	= length;
	size_t halflength	= length>>1;
	for (size_t i=1;i<halflength;i++)
	{
		//here the FB case
		typename TFFT::OFWTYPE H1R		= ((static_cast<typename TFFT::OFWTYPE>(transform[i].re)+transform[redlength-i].re)<<exponent)>>1;
		typename TFFT::OFWTYPE H1I		= ((static_cast<typename TFFT::OFWTYPE>(transform[i].im)-transform[redlength-i].im)<<exponent)>>1;//because conj complex
		typename TFFT::OFWTYPE H2R		= ((static_cast<typename TFFT::OFWTYPE>(transform[i].re)-transform[redlength-i].re)<<exponent)>>1;//that is the second term
		typename TFFT::OFWTYPE H2I		= ((static_cast<typename TFFT::OFWTYPE>(transform[i].im)+transform[redlength-i].im)<<exponent)>>1;//adding, because conj complex

		//to be checked

		TW cosine						= GetCosFactor(i,loglen);//cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		TW sine							= GetSinFactor(i,loglen);//sin(PI*i/length);
		
		typename TFFT::OFWTYPE res1i	= (H2R*cosine+H2I*sine)>>(sizeof(TW)*8-1);
		typename TFFT::OFWTYPE res1r	= -((H2I*cosine-H2R*sine)>>(sizeof(TW)*8-1));//complex i is integrated

		typename TFFT::OFWTYPE resreal	= (H1R-res1r);//0.5f*(H1R-res1r);
		typename TFFT::OFWTYPE resimag	= (H1I-res1i);//0.5f*(H1I-res1i);->already scaled see shift operation above

		//here follows the FN/2-b case (copy from above)
		//->HR1 is still valid
		//double H1R		= transform[i].re+transform[redlength-i].re;
		//H1I is just inverted
		//double H1I		= -transform[i].im-transform[redlength-i].im;//because conj complex

		//use HR1 as above, invert H1I
		//use HR2 as above, invert H2I (formula on paper contain an error)
		//double H2R		= transform[i].re-transform[redlength-i].re;//that is the second term
		//double H2I		= transform[i].im+transform[redlength-i].im;//adding, because conj complex

		//double cosine	= cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		//double sine		= sin(PI*i/length);
		
		res1i							= (H2R*cosine+H2I*sine)>>(sizeof(TW)*8-1);
		res1r							= (H2I*cosine-H2R*sine)>>(sizeof(TW)*8-1);//complex i is integrated

		typename TFFT::OFWTYPE resreal2	= (H1R-res1r);//0.5f*(H1R-res1r);
		typename TFFT::OFWTYPE resimag2 = (-H1I-res1i);//0.5f*(-H1I-res1i);->already scaled see shift operation above

		transform[i].re					= static_cast<typename TFFT::DATATYPE>(resreal>>1);
		transform[i].im					= static_cast<typename TFFT::DATATYPE>(resimag>>1);
		transform[length-i].re			= static_cast<typename TFFT::DATATYPE>(resreal2>>1);
		transform[length-i].im			= static_cast<typename TFFT::DATATYPE>(resimag2>>1);
	}
	//typename TFFT::OFWTYPE lval			= ((static_cast<typename TFFT::OFWTYPE>(transform[0].re)+transform[0].im)>>1)<<exponent;
	//typename TFFT::OFWTYPE hval			= ((static_cast<typename TFFT::OFWTYPE>(transform[0].re)-transform[0].im)>>1)<<exponent;

	typename TFFT::OFWTYPE lval;
	typename TFFT::OFWTYPE hval;

	if (exponent)
	{
		lval							= ((static_cast<typename TFFT::OFWTYPE>(transform[0].re)+transform[0].im))<<(exponent-1);
		hval							= ((static_cast<typename TFFT::OFWTYPE>(transform[0].re)-transform[0].im))<<(exponent-1);
	}
	else
	{
		lval							= ((static_cast<typename TFFT::OFWTYPE>(transform[0].re)+transform[0].im))>>1;
		hval							= ((static_cast<typename TFFT::OFWTYPE>(transform[0].re)-transform[0].im))>>1;
	}

	transform[0].re						= static_cast<typename TFFT::DATATYPE>(lval);//>>1 moved to exponent above
	transform[0].im						= static_cast<typename TFFT::DATATYPE>(hval);//>>1 moved to exponent above
	transform[halflength].re			= transform[halflength].re<<exponent;//treatment of the half index
	transform[halflength].im			= -(transform[halflength].im<<exponent);//treatment of the half index
}

/**
* floating point implementation of the FFT algorithm (single precision). 
* The floating point version(s) return a value with twice the size of the input operant to avoid loss of precision
*/

template<size_t loglen> class CFFT<loglen,float,float,FComplex> {
	private:
#ifdef _USETWIDDLE_
		float	ma_TwiddleFactors[1<<(loglen-1)];

		void    CreateTwiddleFactors(size_t length, float twiddlefactors[]);
		float   GetCosFactor(size_t stage, size_t stageshift);
		float   GetSinFactor(size_t stage, size_t stageshift);
#endif



		template <size_t bitsize> size_t swap(size_t pointer)
		{
			const size_t size			= bitsize>>1;//for an example for 11 size is 5
			const size_t lowermask		= ~(static_cast<size_t>(static_cast<size_t>(-1)>>size)<<size);
			const size_t uppermask		= lowermask<<(bitsize-size);
			const size_t mask			= ~(static_cast<size_t>(static_cast<size_t>(-1)>>bitsize)<<bitsize);
			const size_t remaining		= (~(lowermask|uppermask))&mask;
			pointer						= (swap<size>(lowermask&pointer)<<(bitsize-size))|swap<size>((uppermask&pointer)>>(bitsize-size))|(remaining&pointer);
			return pointer;
		}

		template <> size_t swap<2>(size_t pointer)
		{
			return ((pointer<<1)|(pointer>>1))&3;
		}

		template <> size_t swap<1>(size_t pointer)
		{
			return pointer&1;
		}

		template <size_t bitsize> size_t bitreversal(size_t pointer)
		{
			const size_t size		= (bitsize)>>1;//for an example for 11 size is 5
			const size_t lowermask	= ~(static_cast<size_t>(static_cast<size_t>(-1)>>size)<<size);
			const size_t uppermask	= lowermask<<(bitsize-size);
			const size_t mask		= ~(static_cast<size_t>(static_cast<size_t>(-1)>>(bitsize))<<(bitsize));
			const size_t remaining	= (~(lowermask|uppermask))&mask;
			pointer					= (swap<size>(lowermask&pointer)<<(bitsize-size))|swap<size>((uppermask&pointer)>>(bitsize-size))|(remaining&pointer);
			return pointer;
		}

		void ReorderSamples(float in[], FComplex data[]);
		void ReorderSamples(FComplex in[], FComplex data[]);
		void ReorderRealSamples(float in[], FComplex data[]);
		void ReorderRealSamplesInPlace(FComplex in[]);
		void Convert2HalfDFT(FComplex transform[]);
		void Convert2ComplexDFT(FComplex transform[]);
	public:
#ifdef _USETWIDDLE_
		CFFT(){CreateTwiddleFactors(sizeof(ma_TwiddleFactors)/sizeof(ma_TwiddleFactors[0]), ma_TwiddleFactors);}
#else
		CFFT(){};
#endif
		void CalculateFFT(float buffer[], FComplex ffttransform[]);
		void CalculateIFFT(FComplex buffer[], FComplex ffttransform[]);
		void CalculateRealFFT(float buffer[], FComplex ffttransform[]);
		void CalculateRealIFFT(FComplex buffer[], FComplex ffttransform[]);
};

#ifdef _USETWIDDLE_
template<size_t loglen>
void CFFT<loglen, float, float, FComplex>::CreateTwiddleFactors(size_t length, float twiddlefactors[])
{
	static const float PI2 = 6.283185307179586476925286f;

	size_t len        = 1<<loglen;
	if (2*length>=len)
	{
		const size_t len2 = len>>1;
		for (size_t i=0;i<len2;i++)
		{
			twiddlefactors[i] = static_cast<float>(cos(PI2*static_cast<float>(i)/len));
		}
	}
}
#endif

template <size_t loglen>
void CFFT<loglen, float, float, FComplex>::ReorderSamples(float in[], FComplex data[])
{
	size_t length = 1<<loglen;
	for (size_t i=0;i<length;i++)
	{
		size_t reordered	= bitreversal<loglen>(i);
		data[reordered].re	= in[i];
		data[reordered].im	= 0;
	}
}

template <size_t loglen>
void CFFT<loglen, float, float, FComplex>::ReorderSamples(FComplex in[], FComplex data[])
{
	size_t length = 1<<loglen;
	for (size_t i=0;i<length;i++)
	{
		size_t reordered	= bitreversal<loglen>(i);
		data[reordered].re	= in[i].re;
		data[reordered].im	= in[i].im;
	}
}

template <size_t loglen>
void CFFT<loglen, float, float, FComplex>::ReorderRealSamples(float in[], FComplex data[])
{
	size_t length = 1<<(loglen-1);
	for (size_t i=0;i<length;i++)
	{
		size_t reordered	= bitreversal<loglen-1>(i);
		data[reordered].re	= in[(i<<1)];
		data[reordered].im	= in[(i<<1)+1];
	}
}

template <size_t loglen>
void CFFT<loglen, float, float, FComplex>::ReorderRealSamplesInPlace(FComplex in[])
{
	size_t length = 1<<(loglen-1);
	for (size_t i=0;i<length;i++)
	{
		size_t reordered	= bitreversal<loglen-1>(i);
		if (i>=reordered)
		{
			FComplex dummy;
			dummy = in[reordered];
			in[reordered] = in[i];
			in[i] = dummy;
		}
	}
}

#ifdef _USETWIDDLE_
template <size_t loglen>
float CFFT<loglen, float, float, FComplex>::GetCosFactor(size_t stage, size_t stageshift)
{
	const size_t len = 1<<loglen;
	//more than .5 is not possible
	//size_t k = len*stage/stagesize;
	size_t k = (len*stage)>>stageshift;
	return ma_TwiddleFactors[k];
}

template <size_t loglen>
float CFFT<loglen, float, float, FComplex>::GetSinFactor(size_t stage, size_t stageshift)
{
	const size_t len        = 1<<loglen;
	//const size_t len2        = 1<<(loglen-1);
	const size_t len4        = 1<<(loglen-2);
	const size_t mask        = ~(((-1)>>(loglen))<<(loglen-1));
	//more than .5 is not possible
	//size_t k = (len*stage/stagesize+(len>>2))&mask;
	size_t k = (((len*stage)>>stageshift)+(len>>2))&mask;
	return k<len4?ma_TwiddleFactors[k]:-ma_TwiddleFactors[k];
}
#endif

template <size_t loglen>
void CFFT<loglen, float, float, FComplex>::CalculateFFT(float buffer[], FComplex ffttransform[])
{
	static const float PI2 = 6.2831853071795864769252868f;
	size_t ldlen	= loglen;
	size_t length	= 1<<loglen;
	ReorderSamples(buffer, ffttransform);

	for(size_t i=0; i<ldlen;i++)
	{
		size_t stagesize	= static_cast<size_t>(1)<<(i+1);
		size_t stagesize2	= static_cast<size_t>(1)<<i;
		size_t looplength	= length>>(i+1);

#ifndef _USETWIDDLE_
#ifdef _RECURRENCE_
			const float delta		= PI2/stagesize;
			const float delta_2		= 0.5f*delta;
			const float sin_delta	= sin(delta_2);
			const float alpha		= 2.0f*sin_delta*sin_delta;
			const float beta		= sin(delta);
			float cos_rec			= 1.0f;
			float sin_rec			= 0.0f;
#endif
#endif

		for (size_t R = 0; R < stagesize2; R++)
		{
			FComplex weight;
#ifdef _USETWIDDLE_
			weight.re = GetCosFactor(R, (i + 1));
			weight.im = -GetSinFactor(R, (i + 1));
#else
#ifndef _RECURRENCE_
			const  float RT_Ratio = static_cast<float>(R) / stagesize;
			weight.re = cos(PI2 * RT_Ratio);
			weight.im = -sin(PI2 * RT_Ratio);
#else
			float sin_old = sin_rec;
			float cos_old = cos_rec;
			cos_rec = cos_old - (alpha * cos_old + beta * sin_old);
			sin_rec = sin_old - (alpha * sin_old - beta * cos_old);
			weight.re = cos_old;
			weight.im = -sin_old;
#endif
#endif

			for (size_t L = 0; L < looplength; L++)
			{
				size_t butterflyix = L << (i + 1);
				FComplex temp_a	=	ffttransform[butterflyix+R];
				FComplex temp_b =	ffttransform[butterflyix+R+stagesize2];//load individual stages
#ifdef _C1x_
				auto dummy								= weight*temp_b;
				ffttransform[butterflyix+R]				= temp_a+dummy;
				ffttransform[butterflyix+R+stagesize2]	= temp_a-dummy;
#else
				FComplex::ACCTYPE dummy					= weight*temp_b;
				ffttransform[butterflyix+R]				= temp_a+dummy;
				ffttransform[butterflyix+R+stagesize2]	= temp_a-dummy;
#endif
			}
		}
	}
}

template <size_t loglen>
void CFFT<loglen, float, float, FComplex>::CalculateIFFT(FComplex buffer[], FComplex ffttransform[])
{
	static const float PI2 = 6.2831853071795864769252868f;
	size_t ldlen	= loglen;
	size_t length	= 1<<loglen;
	ReorderSamples(buffer, ffttransform);

	for(size_t i=0; i<ldlen;i++)
	{
		size_t stagesize	= static_cast<size_t>(1)<<(i+1);
		size_t stagesize2	= static_cast<size_t>(1)<<i;
		size_t looplength	= length>>(i+1);

#ifndef _USETWIDDLE_
#ifdef _RECURRENCE_
			const float delta		= PI2/stagesize;
			const float delta_2		= 0.5f*delta;
			const float sin_delta	= sin(delta_2);
			const float alpha		= 2.0f*sin_delta*sin_delta;
			const float beta		= sin(delta);
			float cos_rec = 1.0f;
			float sin_rec = 0.0f;
#endif
#endif

		for (size_t R = 0; R < stagesize2; R++)
		{
			FComplex weight;
#ifdef _USETWIDDLE_
			weight.re = GetCosFactor(R, (i + 1));
			weight.im = GetSinFactor(R, (i + 1));
#else
#ifndef _RECURRENCE_
			const  float RT_Ratio = static_cast<float>(R) / stagesize;
			weight.re = cos(PI2 * RT_Ratio);
			weight.im = sin(PI2 * RT_Ratio);
#else
			float sin_old = sin_rec;
			float cos_old = cos_rec;
			cos_rec = cos_old - (alpha * cos_old + beta * sin_old);
			sin_rec = sin_old - (alpha * sin_old - beta * cos_old);
			weight.re = cos_old;
			weight.im = sin_old;
#endif
#endif
			for (size_t L = 0; L < looplength; L++)
			{
				size_t butterflyix = L << (i + 1);
				FComplex temp_a	=	ffttransform[butterflyix+R];
				FComplex temp_b =	ffttransform[butterflyix+R+stagesize2];//load individual stages
				
#ifdef _C1x_
				auto dummy								= weight*temp_b;
				ffttransform[butterflyix+R]				= temp_a+dummy;
				ffttransform[butterflyix+R+stagesize2]	= temp_a-dummy;
#else
				FComplex::ACCTYPE dummy					= weight*temp_b;
				ffttransform[butterflyix+R]				= temp_a+dummy;
				ffttransform[butterflyix+R+stagesize2]	= temp_a-dummy;
#endif
			}
		}
	}

	const float scaling = 1.0f/length;
	for (size_t i=0;i<length;i++)
	{
		ffttransform[i].re *= scaling;
		ffttransform[i].im *= scaling;
	}
}

template <size_t loglen>
void CFFT<loglen, float, float, FComplex>::CalculateRealFFT(float buffer[], FComplex ffttransform[])
{
	size_t ldlen	= loglen-1;
	size_t length	= 1<<(loglen-1);
	ReorderRealSamples(buffer, ffttransform);

	for(size_t i=0; i<ldlen;i++)
	{
#ifndef _USETWIDDLE_
		size_t stagesize	= static_cast<size_t>(1)<<(i+1);
#endif
		size_t stagesize2	= static_cast<size_t>(1)<<i;
		size_t looplength	= length>>(i+1);
		for (size_t R = 0; R < stagesize2; R++)
		{
			FComplex weight;
#ifdef _USETWIDDLE_
			weight.re = GetCosFactor(R, (i + 1));
			weight.im = -GetSinFactor(R, (i + 1));
#else
			static const float PI2 = 6.2831853071795864769252868f;
			const  float RT_Ratio = static_cast<float>(R) / stagesize;
			weight.re = cos(PI2 * RT_Ratio);
			weight.im = -sin(PI2 * RT_Ratio);
#endif
			for (size_t L = 0; L < looplength; L++)
			{
				size_t butterflyix = L << (i + 1);
				FComplex temp_a	=	ffttransform[butterflyix+R];
				FComplex temp_b =	ffttransform[butterflyix+R+stagesize2];//load individual stages
				
#ifdef _C1x_
				auto dummy								= weight*temp_b;
				ffttransform[butterflyix+R]				= temp_a+dummy;
				ffttransform[butterflyix+R+stagesize2]	= temp_a-dummy;
#else
				FComplex::ACCTYPE dummy					= weight*temp_b;
				ffttransform[butterflyix+R]				= temp_a+dummy;
				ffttransform[butterflyix+R+stagesize2]	= temp_a-dummy;
#endif
			}
		}
	}
	Convert2HalfDFT(ffttransform);
}

template <size_t loglen>
void CFFT<loglen, float, float, FComplex>::CalculateRealIFFT(FComplex buffer[], FComplex ffttransform[])
{
	static const float PI2 = 6.2831853071795864769252868f;
	size_t ldlen	= loglen-1;
	size_t length	= 1<<(loglen-1);
	memcpy(ffttransform, buffer,length*sizeof(FComplex));
	Convert2ComplexDFT(ffttransform);
	ReorderRealSamplesInPlace(ffttransform);

	for(size_t i=0; i<ldlen;i++)
	{
		size_t stagesize	= static_cast<size_t>(1)<<(i+1);
		size_t stagesize2	= static_cast<size_t>(1)<<i;
		size_t looplength	= length>>(i+1);

#ifndef _USETWIDDLE_
#ifdef _RECURRENCE_
			const float delta		= PI2/stagesize;
			const float delta_2		= 0.5f*delta;
			const float sin_delta	= sin(delta_2);
			const float alpha		= 2.0f*sin_delta*sin_delta;
			const float beta		= sin(delta);
			float cos_rec = 1.0f;
			float sin_rec = 0.0f;
#endif
#endif

		for (size_t R = 0; R < stagesize2; R++)
		{
			FComplex weight;
#ifdef _USETWIDDLE_
			weight.re = GetCosFactor(R, (i + 1));
			weight.im = GetSinFactor(R, (i + 1));
#else
#ifndef _RECURRENCE_
			const  float RT_Ratio = static_cast<float>(R) / stagesize;
			weight.re = cos(PI2 * RT_Ratio);
			weight.im = sin(PI2 * RT_Ratio);
#else
			float sin_old = sin_rec;
			float cos_old = cos_rec;
			cos_rec = cos_old - (alpha * cos_old + beta * sin_old);
			sin_rec = sin_old - (alpha * sin_old - beta * cos_old);
			weight.re = cos_old;
			weight.im = sin_old;
#endif
#endif
			for (size_t L = 0; L < looplength; L++)
			{
				size_t butterflyix = L << (i + 1);
				FComplex temp_a	=	ffttransform[butterflyix+R];
				FComplex temp_b =	ffttransform[butterflyix+R+stagesize2];//load individual stages

#ifdef _C1x_
				auto dummy								= weight*temp_b;
				ffttransform[butterflyix+R]				= temp_a+dummy;
				ffttransform[butterflyix+R+stagesize2]	= temp_a-dummy;
#else
				FComplex::ACCTYPE dummy					= weight*temp_b;
				ffttransform[butterflyix+R]				= temp_a+dummy;
				ffttransform[butterflyix+R+stagesize2]	= temp_a-dummy;
#endif
			}
		}
	}

	const float scaling = 1.0f/length;
	for (size_t i=0;i<length;i++)
	{
		ffttransform[i].re *= scaling;
		ffttransform[i].im *= scaling;
	}
}

template <size_t loglen>
void CFFT<loglen, float, float, FComplex>::Convert2HalfDFT(FComplex transform[])
{
	static const float PI2 = 6.2831853071795864769252868f;
	static const float PI  = 3.1415926535897932384626434f;

	size_t length		= 1<<(loglen-1);
	size_t redlength	= length;
	size_t halflength	= length>>1;

#ifndef _USETWIDDLE_
#ifdef _RECURRENCE_
			const float	delta		= PI/length;
			const float	delta_2		= 0.5f*delta;
			const float	sin_delta	= sin(delta_2);
			const float	alpha		= 2.0f*sin_delta*sin_delta;
			const float	beta		= sin(delta);
			float		 sineiter	= 0;
			float		 cosineiter	= 1.0f;
#endif
#endif

	for (size_t i=1;i<halflength;i++)
	{
		//here the FB case
		float H1R		= transform[i].re+transform[redlength-i].re;
		float H1I		= transform[i].im-transform[redlength-i].im;//because conj complex
		float H2R		= transform[i].re-transform[redlength-i].re;//that is the second term
		float H2I		= transform[i].im+transform[redlength-i].im;//adding, because conj complex

#ifdef _USETWIDDLE_
		float cosine	= GetCosFactor(i,loglen);//cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		float sine		= GetSinFactor(i,loglen);//sin(PI*i/length);
#else
#ifdef _RECURRENCE_

		float cosine	= cosineiter;
		float sine		= sineiter;
		cosineiter		= cosine-(alpha*cosine+beta*sine);
		sineiter		= sineiter-(alpha*sine-beta*cosine);
		cosine			= cosineiter;
		sine			= sineiter;
#else
		float cosine	= cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		float sine		= sin(PI*i/length);
#endif
#endif
		
		float res1i	= H2R*cosine+H2I*sine;
		float res1r	= -(H2I*cosine-H2R*sine);//complex i is integrated

		float resreal	= 0.5f*(H1R-res1r);
		float resimag	= 0.5f*(H1I-res1i);

		//here follows the FN/2-b case (copy from above)
		//->HR1 is still valid
		//double H1R		= transform[i].re+transform[redlength-i].re;
		//H1I is just inverted
		//double H1I		= -transform[i].im-transform[redlength-i].im;//because conj complex

		//use HR1 as above, invert H1I
		//use HR2 as above, invert H2I (formula on paper contain an error)
		//double H2R		= transform[i].re-transform[redlength-i].re;//that is the second term
		//double H2I		= transform[i].im+transform[redlength-i].im;//adding, because conj complex

		//double cosine	= cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		//double sine		= sin(PI*i/length);
		
		res1i					= H2R*cosine+H2I*sine;
		res1r					= (H2I*cosine-H2R*sine);//complex i is integrated

		float resreal2			= 0.5f*(H1R-res1r);
		float resimag2			= 0.5f*(-H1I-res1i);

		transform[i].re			= resreal;
		transform[i].im			= resimag;
		transform[length-i].re	= resreal2;
		transform[length-i].im	= resimag2;
	}
	float lval					= transform[0].re+transform[0].im;
	float hval					= transform[0].re-transform[0].im;
	transform[0].re				= lval;
	transform[0].im				= hval;
	transform[halflength].im	= -transform[halflength].im;//treatment of the half index
}

template <size_t loglen>
void CFFT<loglen, float, float, FComplex>::Convert2ComplexDFT(FComplex transform[])
{
	static const float PI2 = 6.2831853071795864769252868f;
	static const float PI  = 3.1415926535897932384626434f;

	size_t length		= 1<<(loglen-1);
	size_t redlength	= length;
	size_t halflength	= length>>1;
#ifndef _USETWIDDLE_
#ifdef _RECURRENCE_
			const float	delta		= PI/length;
			const float	delta_2	= 0.5f*delta;
			const float	sin_delta	= sin(delta_2);
			const float	alpha		= 2.0f*sin_delta*sin_delta;
			const float	beta		= sin(delta);
			float		 sineiter	= 0.0f;
			float		 cosineiter	= 1.0f;
#endif
#endif
	for (size_t i=1;i<halflength;i++)
	{
		//here the FB case
		float H1R		= transform[i].re+transform[redlength-i].re;
		float H1I		= transform[i].im-transform[redlength-i].im;//because conj complex
		float H2R		= transform[i].re-transform[redlength-i].re;//that is the second term
		float H2I		= transform[i].im+transform[redlength-i].im;//adding, because conj complex

#ifdef _USETWIDDLE_
		float cosine	= GetCosFactor(i,loglen);//cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		float sine		= GetSinFactor(i,loglen);//sin(PI*i/length);
#else
#ifdef _RECURRENCE_

		float cosine	= cosineiter;
		float sine		= sineiter;
		cosineiter		= cosine-(alpha*cosine+beta*sine);
		sineiter		= sineiter-(alpha*sine-beta*cosine);
		cosine			= cosineiter;
		sine			= sineiter;
#else
		float cosine	= cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		float sine		= sin(PI*i/length);
#endif
#endif
		
		float res1i	= H2R*cosine-H2I*sine;
		float res1r	= -(H2I*cosine+H2R*sine);//complex i is integrated

		float resreal	= 0.5f*(H1R+res1r);
		float resimag	= 0.5f*(H1I+res1i);

		//here follows the FN/2-b case (copy from above)
		//->HR1 is still valid
		//double H1R		= transform[i].re+transform[redlength-i].re;
		//H1I is just inverted
		//double H1I		= -transform[i].im-transform[redlength-i].im;//because conj complex

		//use HR1 as above, invert H1I
		//use HR2 as above, invert H2I (formula on paper contain an error)
		//double H2R		= transform[i].re-transform[redlength-i].re;//that is the second term
		//double H2I		= transform[i].im+transform[redlength-i].im;//adding, because conj complex

		//double cosine	= cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		//double sine		= sin(PI*i/length);
		
		res1i					= -H2R*cosine+H2I*sine;
		res1r					= (H2I*cosine+H2R*sine);//complex i is integrated

		float resreal2			= 0.5f*(H1R+res1r);
		float resimag2			= -0.5f*(H1I+res1i);

		transform[i].re			= resreal;
		transform[i].im			= resimag;
		transform[length-i].re	= resreal2;
		transform[length-i].im	= resimag2;
	}
	float lval					= transform[0].re+transform[0].im;
	float hval					= transform[0].re-transform[0].im;
	transform[0].re				= 0.5f*lval;
	transform[0].im				= 0.5f*hval;
	transform[halflength].im	= -transform[halflength].im;//treatment of the half index
}


/**
* floating point implementation of the FFT algorithm (single precision). 
* The floating point version(s) return a value with twice the size of the input operant to avoid loss of precision
*/

template<size_t loglen> class CFFT<loglen,double,double,DComplex> {
	private:
#ifdef _USETWIDDLE_
		double	ma_TwiddleFactors[1<<(loglen-1)];

		void    CreateTwiddleFactors(size_t length, double twiddlefactors[]);

		double	GetCosFactor(size_t stage, size_t stageshift);
		double	GetSinFactor(size_t stage, size_t stageshift);
#endif

		template <size_t bitsize> size_t swap(size_t pointer)
		{
			const size_t size			= bitsize>>1;//for an example for 11 size is 5
			const size_t lowermask		= ~(static_cast<size_t>(static_cast<size_t>(-1)>>size)<<size);
			const size_t uppermask		= lowermask<<(bitsize-size);
			const size_t mask			= ~(static_cast<size_t>(static_cast<size_t>(-1)>>bitsize)<<bitsize);
			const size_t remaining		= (~(lowermask|uppermask))&mask;
			pointer						= (swap<size>(lowermask&pointer)<<(bitsize-size))|swap<size>((uppermask&pointer)>>(bitsize-size))|(remaining&pointer);
			return pointer;
		}

		template <> size_t swap<2>(size_t pointer)
		{
			return ((pointer<<1)|(pointer>>1))&3;
		}

		template <> size_t swap<1>(size_t pointer)
		{
			return pointer&1;
		}

		template <size_t bitsize> size_t bitreversal(size_t pointer)
		{
			const size_t size		= (bitsize)>>1;//for an example for 11 size is 5
			const size_t lowermask	= ~(static_cast<size_t>(static_cast<size_t>(-1)>>size)<<size);
			const size_t uppermask	= lowermask<<(bitsize-size);
			const size_t mask		= ~(static_cast<size_t>(static_cast<size_t>(-1)>>(bitsize))<<(bitsize));
			const size_t remaining	= (~(lowermask|uppermask))&mask;
			pointer					= (swap<size>(lowermask&pointer)<<(bitsize-size))|swap<size>((uppermask&pointer)>>(bitsize-size))|(remaining&pointer);
			return pointer;
		}

		void ReorderSamples(double in[], DComplex data[]);
		void ReorderSamples(DComplex in[], DComplex data[]);
		void ReorderRealSamples(double in[], DComplex data[]);
		void ReorderRealSamplesInPlace(DComplex in[]);
		void Convert2HalfDFT(DComplex ffttransform[]);
		void Convert2ComplexDFT(DComplex transform[]);
	public:
#ifdef _DEBUG_RFFT_
		void Convert2HalfDFT(DComplex ffttransform[], DComplex result[]);
#endif
#ifdef _USETWIDDLE_
		CFFT(){CreateTwiddleFactors(sizeof(ma_TwiddleFactors)/sizeof(ma_TwiddleFactors[0]), ma_TwiddleFactors);}
#else
		CFFT(){};
#endif
		void CalculateFFT(double buffer[], DComplex ffttransform[]);
		void CalculateIFFT(DComplex buffer[], DComplex ffttransform[]);
		void CalculateRealFFT(double buffer[], DComplex ffttransform[]);
		void CalculateRealIFFT(DComplex buffer[], DComplex ffttransform[]);
};

#ifdef _USETWIDDLE_
template<size_t loglen>
void CFFT<loglen, double, double, DComplex>::CreateTwiddleFactors(size_t length, double twiddlefactors[])
{
	static const double PI2 = 6.283185307179586476925286f;

	size_t len        = 1<<loglen;
	if (2*length>=len)
	{
		const size_t len2 = len>>1;
		for (size_t i=0;i<len2;i++)
		{
			twiddlefactors[i] = static_cast<double>(cos(PI2*static_cast<double>(i)/len));
		}
	}
}
#endif

template <size_t loglen>
void CFFT<loglen, double, double, DComplex>::ReorderSamples(double in[], DComplex data[])
{
	size_t length = 1<<loglen;
	for (size_t i=0;i<length;i++)
	{
		size_t reordered	= bitreversal<loglen>(i);
		data[reordered].re	= in[i];
		data[reordered].im	= 0;
	}
}

template <size_t loglen>
void CFFT<loglen, double, double, DComplex>::ReorderSamples(DComplex in[], DComplex data[])
{
	size_t length = 1<<loglen;
	for (size_t i=0;i<length;i++)
	{
		size_t reordered	= bitreversal<loglen>(i);
		data[reordered].re	= in[i].re;
		data[reordered].im	= in[i].im;
	}
}

template <size_t loglen>
void CFFT<loglen, double, double, DComplex>::ReorderRealSamples(double in[], DComplex data[])
{
	size_t length = 1<<(loglen-1);
	for (size_t i=0;i<length;i++)
	{
		size_t reordered	= bitreversal<loglen-1>(i);
		data[reordered].re	= in[(i<<1)];
		data[reordered].im	= in[(i<<1)+1];
	}
}

template <size_t loglen>
void CFFT<loglen, double, double, DComplex>::ReorderRealSamplesInPlace(DComplex in[])
{
	size_t length = 1<<(loglen-1);
	for (size_t i=0;i<length;i++)
	{
		size_t reordered	= bitreversal<loglen-1>(i);
		if (i>=reordered)
		{
			DComplex dummy;
			dummy = in[reordered];
			in[reordered] = in[i];
			in[i] = dummy;
		}
	}
}

#ifdef _USETWIDDLE_
template <size_t loglen>
double CFFT<loglen, double, double, DComplex>::GetCosFactor(size_t stage, size_t stageshift)
{
	const size_t len = 1<<loglen;
	//more than .5 is not possible
	//size_t k = len*stage/stagesize;
	size_t k = (len*stage)>>stageshift;
	return ma_TwiddleFactors[k];
}

template <size_t loglen>
double CFFT<loglen, double, double, DComplex>::GetSinFactor(size_t stage, size_t stageshift)
{
	const size_t len        = 1<<loglen;
	//const size_t len2        = 1<<(loglen-1);
	const size_t len4        = 1<<(loglen-2);
	const size_t mask        = ~(((-1)>>(loglen))<<(loglen-1));
	//more than .5 is not possible
	//size_t k = (len*stage/stagesize+(len>>2))&mask;
	size_t k = (((len*stage)>>stageshift)+(len>>2))&mask;
	return k<len4?ma_TwiddleFactors[k]:-ma_TwiddleFactors[k];
}
#endif

template <size_t loglen>
void CFFT<loglen, double, double, DComplex>::CalculateFFT(double buffer[], DComplex ffttransform[])
{
	static const double PI2 = 6.2831853071795864769252868;
	size_t ldlen	= loglen;
	size_t length	= 1<<loglen;
	ReorderSamples(buffer, ffttransform);

	for(size_t i=0; i<ldlen;i++)
	{
		size_t stagesize	= 1<<(i+1);
		size_t stagesize2	= 1<<i;
		size_t looplength	= length>>(i+1);

#ifndef _USETWIDDLE_
#ifdef _RECURRENCE_
			const double delta		= PI2/stagesize;
			const double delta_2	= 0.5*delta;
			const double sin_delta	= sin(delta_2);
			const double alpha		= 2.0*sin_delta*sin_delta;
			const double beta		= sin(delta);
			double cos_rec			= 1.0f;
			double sin_rec			= 0.0f;
#endif
#endif
		for (size_t R = 0; R < stagesize2; R++)
		{
			DComplex weight;
#ifdef _USETWIDDLE_
			weight.re = GetCosFactor(R, (i + 1));
			weight.im = -GetSinFactor(R, (i + 1));
#else
#ifndef _RECURRENCE_
			const  double RT_Ratio = static_cast<float>(R) / stagesize;
			weight.re = cos(PI2 * RT_Ratio);
			weight.im = -sin(PI2 * RT_Ratio);
#else
			double sin_old = sin_rec;
			double cos_old = cos_rec;
			cos_rec = cos_old - (alpha * cos_old + beta * sin_old);
			sin_rec = sin_old - (alpha * sin_old - beta * cos_old);
			weight.re = cos_old;
			weight.im = -sin_old;
#endif
#endif
			for (size_t L = 0; L < looplength; L++)
			{
				size_t butterflyix = L << (i + 1);
				DComplex temp_a	=	ffttransform[butterflyix+R];
				DComplex temp_b =	ffttransform[butterflyix+R+stagesize2];//load individual stages

#ifdef _C1x_
				auto dummy								= weight*temp_b;
				ffttransform[butterflyix+R]				= temp_a+dummy;
				ffttransform[butterflyix+R+stagesize2]	= temp_a-dummy;
#else
				DComplex::ACCTYPE dummy					= weight*temp_b;
				ffttransform[butterflyix+R]				= temp_a+dummy;
				ffttransform[butterflyix+R+stagesize2]	= temp_a-dummy;
#endif
			}
		}
	}
}

template <size_t loglen>
void CFFT<loglen, double, double, DComplex>::CalculateIFFT(DComplex buffer[], DComplex ffttransform[])
{
	static const double PI2 = 6.2831853071795864769252868;
	size_t ldlen	= loglen;
	size_t length	= 1<<loglen;
	ReorderSamples(buffer, ffttransform);

	for(size_t i=0; i<ldlen;i++)
	{
		size_t stagesize	= 1<<(i+1);
		size_t stagesize2	= 1<<i;
		size_t looplength	= length>>(i+1);

#ifndef _USETWIDDLE_
#ifdef _RECURRENCE_
			const double delta		= PI2/stagesize;
			const double delta_2	= 0.5*delta;
			const double sin_delta	= sin(delta_2);
			const double alpha		= 2.0*sin_delta*sin_delta;
			const double beta		= sin(delta);
			double cos_rec			= 1.0f;
			double sin_rec			= 0.0f;
#endif
#endif
		for (size_t R = 0; R < stagesize2; R++)
		{
			DComplex weight;
#ifdef _USETWIDDLE_
			weight.re = GetCosFactor(R, (i + 1));
			weight.im = GetSinFactor(R, (i + 1));
#else
#ifndef _RECURRENCE_
			const  double RT_Ratio = static_cast<float>(R) / stagesize;
			weight.re = cos(PI2 * RT_Ratio);
			weight.im = sin(PI2 * RT_Ratio);
#else
			double sin_old = sin_rec;
			double cos_old = cos_rec;
			cos_rec = cos_old - (alpha * cos_old + beta * sin_old);
			sin_rec = sin_old - (alpha * sin_old - beta * cos_old);
			weight.re = cos_old;
			weight.im = sin_old;
#endif
#endif
			for (size_t L = 0; L < looplength; L++)
			{
				size_t butterflyix = L << (i + 1);
				DComplex temp_a	=	ffttransform[butterflyix+R];
				DComplex temp_b =	ffttransform[butterflyix+R+stagesize2];//load individual stages
#ifdef _C1x_
				auto dummy								= weight*temp_b;
				ffttransform[butterflyix+R]				= temp_a+dummy;
				ffttransform[butterflyix+R+stagesize2]	= temp_a-dummy;
#else
				DComplex::ACCTYPE dummy					= weight*temp_b;
				ffttransform[butterflyix+R]				= temp_a+dummy;
				ffttransform[butterflyix+R+stagesize2]	= temp_a-dummy;
#endif
			}
		}
	}
	
	const double scaling = 1.0/length;
	for (size_t i=0;i<length;i++)
	{
		ffttransform[i].re *= scaling;
		ffttransform[i].im *= scaling;
	}
}

template <size_t loglen>
void CFFT<loglen, double, double, DComplex>::CalculateRealIFFT(DComplex buffer[], DComplex ffttransform[])
{
	static const double PI2	= 6.2831853071795864769252868;
	size_t ldlen			= loglen-1;
	size_t length			= 1<<(loglen-1);
	memcpy(ffttransform, buffer,length*sizeof(DComplex));
	Convert2ComplexDFT(ffttransform);
	ReorderRealSamplesInPlace(ffttransform);

	for(size_t i=0; i<ldlen;i++)
	{
		size_t stagesize	= static_cast<size_t>(1)<<static_cast<size_t>(i+1);
		size_t stagesize2	= static_cast<size_t>(1)<<i;
		size_t looplength	= length>>(i+1);

#ifndef _USETWIDDLE_
#ifdef _RECURRENCE_
			const double delta		= PI2/stagesize;
			const double delta_2	= 0.5*delta;
			const double sin_delta	= sin(delta_2);
			const double alpha		= 2.0*sin_delta*sin_delta;
			const double beta		= sin(delta);
			double cos_rec			= 1.0;
			double sin_rec			= 0.0;
#endif
#endif
		for (size_t R = 0; R < stagesize2; R++)
		{
				DComplex weight;
#ifdef _USETWIDDLE_
				weight.re = GetCosFactor(R, (i + 1));
				weight.im = GetSinFactor(R, (i + 1));
#else
#ifndef _RECURRENCE_
				const  double RT_Ratio = static_cast<float>(R) / stagesize;
				weight.re = cos(PI2 * RT_Ratio);
				weight.im = sin(PI2 * RT_Ratio);
#else
				double sin_old = sin_rec;
				double cos_old = cos_rec;
				cos_rec = cos_old - (alpha * cos_old + beta * sin_old);
				sin_rec = sin_old - (alpha * sin_old - beta * cos_old);
				weight.re = cos_old;
				weight.im = sin_old;
#endif
#endif
			for (size_t L = 0; L < looplength; L++)
			{
				size_t butterflyix = L << (i + 1);
				DComplex temp_a	=	ffttransform[butterflyix+R];
				DComplex temp_b =	ffttransform[butterflyix+R+stagesize2];//load individual stages
#ifdef _C1x_
				auto dummy								= weight*temp_b;
				ffttransform[butterflyix+R]				= temp_a+dummy;
				ffttransform[butterflyix+R+stagesize2]	= temp_a-dummy;
#else
				DComplex::ACCTYPE dummy					= weight*temp_b;
				ffttransform[butterflyix+R]				= temp_a+dummy;
				ffttransform[butterflyix+R+stagesize2]	= temp_a-dummy;
#endif
			}
		}
	}
	
	const double scaling = 1.0/length;
	for (size_t i=0;i<length;i++)
	{
		ffttransform[i].re *= scaling;
		ffttransform[i].im *= scaling;
	}
}

template <size_t loglen>
void CFFT<loglen, double, double, DComplex>::CalculateRealFFT(double buffer[], DComplex ffttransform[])
{
	size_t ldlen	= loglen-1;
	size_t length	= 1<<(loglen-1);
	ReorderRealSamples(buffer, ffttransform);

	for(size_t i=0; i<ldlen;i++)
	{
#ifndef _USETWIDDLE_
		size_t stagesize	= static_cast<size_t>(1)<<(i+1);
#endif
		size_t stagesize2	= static_cast<size_t>(1)<<i;
		size_t looplength	= length>>(i+1);
		for (size_t R = 0; R < stagesize2; R++)
		{
			DComplex weight;
#ifdef _USETWIDDLE_
			weight.re = GetCosFactor(R, (i + 1));
			weight.im = -GetSinFactor(R, (i + 1));
#else
			static const double PI2 = 6.2831853071795864769252868;
			const  double RT_Ratio = static_cast<float>(R) / stagesize;
			weight.re = cos(PI2 * RT_Ratio);
			weight.im = -sin(PI2 * RT_Ratio);
#endif
			for (size_t L = 0; L < looplength; L++)
			{
				size_t butterflyix = L << (i + 1);
				DComplex temp_a	=	ffttransform[butterflyix+R];
				DComplex temp_b =	ffttransform[butterflyix+R+stagesize2];//load individual stages		
#ifdef _C1x_
				auto dummy								= weight*temp_b;
				ffttransform[butterflyix+R]				= temp_a+dummy;
				ffttransform[butterflyix+R+stagesize2]	= temp_a-dummy;
#else
				DComplex::ACCTYPE dummy					= weight*temp_b;
				ffttransform[butterflyix+R]				= temp_a+dummy;
				ffttransform[butterflyix+R+stagesize2]	= temp_a-dummy;
#endif
			}
		}
	}
	Convert2HalfDFT(ffttransform);
}

#ifdef _DEBUG_RFFT_
template <size_t loglen>
void CFFT<loglen, double, double, DComplex>:: Convert2HalfDFT(DComplex transform[], DComplex result[])
{
	static const double PI2 = 6.2831853071795864769252868;
	static const double PI  = 3.1415926535897932384626434;

	size_t length	= 1<<(loglen-1);
	size_t redlength= length;
	for (size_t i=1;i<length;i++)
	{
		double H1R		= transform[i].re+transform[redlength-i].re;
		double H1I		= transform[i].im-transform[redlength-i].im;//because conj complex

		double H2R		= transform[i].re-transform[redlength-i].re;//that is the second term
		double H2I		= transform[i].im+transform[redlength-i].im;//adding, because conj complex

		double cosine	= cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		double sine		= sin(PI*i/length);
		
		double res1i	= H2R*cosine+H2I*sine;
		double res1r	= -(H2I*cosine-H2R*sine);//complex i is integrated

		double resreal	= 0.5*(H1R-res1r);
		double resimag	= 0.5*(H1I-res1i);

		result[i].re = resreal;
		result[i].im = resimag;
	}
	double lval = transform[0].re+transform[0].im;
	double hval = transform[0].re-transform[0].im;
	result[0].re = lval;
	result[0].im = hval;
}
#endif

template <size_t loglen>
void CFFT<loglen, double, double, DComplex>::Convert2HalfDFT(DComplex transform[])
{
	static const double PI2 = 6.2831853071795864769252868;
	static const double PI  = 3.1415926535897932384626434;

	size_t length		= 1<<(loglen-1);
	size_t redlength	= length;
	size_t halflength	= length>>1;

#ifndef _USETWIDDLE_
#ifdef _RECURRENCE_
			const double delta		= PI/length;
			const double delta_2	= 0.5*delta;
			const double sin_delta	= sin(delta_2);
			const double alpha		= 2.0*sin_delta*sin_delta;
			const double beta		= sin(delta);
			double		 sineiter	= 0;
			double		 cosineiter	= 1;
#endif
#endif

	for (size_t i=1;i<halflength;i++)
	{
		//here the FB case
		double H1R		= transform[i].re+transform[redlength-i].re;
		double H1I		= transform[i].im-transform[redlength-i].im;//because conj complex
		double H2R		= transform[i].re-transform[redlength-i].re;//that is the second term
		double H2I		= transform[i].im+transform[redlength-i].im;//adding, because conj complex

#ifdef _USETWIDDLE_
		double cosine	= GetCosFactor(i,loglen);//cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		double sine		= GetSinFactor(i,loglen);//sin(PI*i/length);
#else
#ifdef _RECURRENCE_

		double cosine	= cosineiter;
		double sine		= sineiter;
		cosineiter		= cosine-(alpha*cosine+beta*sine);
		sineiter		= sineiter-(alpha*sine-beta*cosine);
		cosine			= cosineiter;
		sine			= sineiter;
#else
		double cosine	= cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		double sine		= sin(PI*i/length);
#endif
#endif
		
		double res1i	= H2R*cosine+H2I*sine;
		double res1r	= -(H2I*cosine-H2R*sine);//complex i is integrated

		double resreal	= 0.5*(H1R-res1r);
		double resimag	= 0.5*(H1I-res1i);

		//here follows the FN/2-b case (copy from above)
		//->HR1 is still valid
		//double H1R		= transform[i].re+transform[redlength-i].re;
		//H1I is just inverted
		//double H1I		= -transform[i].im-transform[redlength-i].im;//because conj complex

		//use HR1 as above, invert H1I
		//use HR2 as above, invert H2I (formula on paper contain an error)
		//double H2R		= transform[i].re-transform[redlength-i].re;//that is the second term
		//double H2I		= transform[i].im+transform[redlength-i].im;//adding, because conj complex

		//double cosine	= cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		//double sine		= sin(PI*i/length);
		
		res1i					= H2R*cosine+H2I*sine;
		res1r					= (H2I*cosine-H2R*sine);//complex i is integrated

		double resreal2			= 0.5*(H1R-res1r);
		double resimag2			= 0.5*(-H1I-res1i);

		transform[i].re			= resreal;
		transform[i].im			= resimag;
		transform[length-i].re	= resreal2;
		transform[length-i].im	= resimag2;
	}
	double lval					= transform[0].re+transform[0].im;
	double hval					= transform[0].re-transform[0].im;
	transform[0].re				= lval;
	transform[0].im				= hval;
	transform[halflength].im	= -transform[halflength].im;//treatment of the half index
}

template <size_t loglen>
void CFFT<loglen, double, double, DComplex>::Convert2ComplexDFT(DComplex transform[])
{
	static const double PI2 = 6.2831853071795864769252868;
	static const double PI  = 3.1415926535897932384626434;

	size_t length		= 1<<(loglen-1);
	size_t redlength	= length;
	size_t halflength	= length>>1;

#ifndef _USETWIDDLE_
#ifdef _RECURRENCE_
	const double delta		= PI/length;
	const double delta_2	= 0.5*delta;
	const double sin_delta	= sin(delta_2);
	const double alpha		= 2.0*sin_delta*sin_delta;
	const double beta		= sin(delta);
	double		 sineiter	= 0.;
	double		 cosineiter	= 1.;
#endif
#endif

	for (size_t i=1;i<halflength;i++)
	{
		//here the FB case
		double H1R		= transform[i].re+transform[redlength-i].re;
		double H1I		= transform[i].im-transform[redlength-i].im;//because conj complex
		double H2R		= transform[i].re-transform[redlength-i].re;//that is the second term
		double H2I		= transform[i].im+transform[redlength-i].im;//adding, because conj complex

#ifdef _USETWIDDLE_
		double cosine	= GetCosFactor(i,loglen);//cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		double sine		= GetSinFactor(i,loglen);//sin(PI*i/length);
#else
#ifdef _RECURRENCE_
		double cosine	= cosineiter;
		double sine		= sineiter;
		cosineiter		= cosine-(alpha*cosine+beta*sine);
		sineiter		= sineiter-(alpha*sine-beta*cosine);
		cosine			= cosineiter;
		sine			= sineiter;
#else
		double cosine	= cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		double sine		= sin(PI*i/length);
#endif
#endif
		
		double res1i	= H2R*cosine-H2I*sine;
		double res1r	= -(H2I*cosine+H2R*sine);//complex i is integrated

		double resreal	= 0.5f*(H1R+res1r);
		double resimag	= 0.5f*(H1I+res1i);

		//here follows the FN/2-b case (copy from above)
		//->HR1 is still valid
		//double H1R		= transform[i].re+transform[redlength-i].re;
		//H1I is just inverted
		//double H1I		= -transform[i].im-transform[redlength-i].im;//because conj complex

		//use HR1 as above, invert H1I
		//use HR2 as above, invert H2I (formula on paper contain an error)
		//double H2R		= transform[i].re-transform[redlength-i].re;//that is the second term
		//double H2I		= transform[i].im+transform[redlength-i].im;//adding, because conj complex

		//double cosine	= cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		//double sine		= sin(PI*i/length);
		
		res1i					= -H2R*cosine+H2I*sine;
		res1r					= (H2I*cosine+H2R*sine);//complex i is integrated

		double resreal2			= 0.5f*(H1R+res1r);
		double resimag2			= -0.5f*(H1I+res1i);

		transform[i].re			= resreal;
		transform[i].im			= resimag;
		transform[length-i].re	= resreal2;
		transform[length-i].im	= resimag2;
	}

		double H1R		= transform[0].re+transform[0].re;
		double H1I		= transform[0].im-transform[0].im;//because conj complex
		double H2R		= transform[0].re-transform[0].re;//that is the second term
		double H2I		= transform[0].im+transform[0].im;//adding, because conj complex
		double cosine	= cos(PI*0/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		double sine		= sin(PI*0/length);

		double res1i	= H2R*cosine-H2I*sine;
		double res1r	= -(H2I*cosine+H2R*sine);//complex i is integrated

		double resreal	= 0.5f*(H1R+res1r);
		double resimag	= 0.5f*(H1I+res1i);


	double lval					= transform[0].re+transform[0].im;
	double hval					= transform[0].re-transform[0].im;
	transform[0].re				= 0.5*lval;
	transform[0].im				= 0.5*hval;
	transform[halflength].im	= -transform[halflength].im;//treatment of the half index
}



size_t Exponent;
static i32Complex TransformScaling[1024];
static i32Complex Transform[1024];
static i32Complex RealTransform[1024];

static int TwiddleFactors[512];
static int Signal[1024];

static FComplex fTransform[1024];
static FComplex fBackTransform[1024];
static float fTwiddle[512];
static float fSignal[1024];

static DComplex dfTransform[1024];
static DComplex dfBackTransform[1024];
static DComplex dfrTransform[512];
static double dfTwiddle[512];
static double dfSignal[1024];
static DComplex Reference[1024];
static DComplex RefSplitSignal[1024];
static DComplex result[1024];
static float resultf[1024];
static double resultd[1024];

void CalculateDFT(size_t length, DComplex *buffer, DComplex *result)
{
	static const double PI2=6.2831853071795864769252868;
	for (size_t i=0;i<length;i++)
	{
		DComplex accu;
		accu.re = 0;
		accu.im = 0;
		for (size_t j=0;j<length;j++)
		{
			double cosine = cos((PI2*j*i)/length);
			double sine = sin((PI2*j*i)/length);
			accu.re += buffer[j].re*cosine+buffer[j].im*sine;
			accu.im += buffer[j].im*cosine-buffer[j].re*sine;
		}
		result[i].re=accu.re;
		result[i].im=accu.im;
	}
}

int _tmain(int argc, _TCHAR* argv[])
{
	CFFT<10,int, int, i32Complex> m_FFT(2147483647);
	CFFT<10,int, int, i32Complex> m_FFT2(2147483647);
	CFFT<10,float, float, FComplex> mf_FFT;
	CFFT<10,double, double, DComplex> mdf_FFT;
	CFFT<10,double, double, DComplex> mrdf_FFT;

	for (unsigned int i=0;i<1024;i++)
	{
		fSignal[i]=1.f*cos(512.f*(6.283185307179586476925286f*static_cast<float>(i))/1024.f)+1.f*sin((2*6.283185307179586476925286f*static_cast<float>(i))/1024.f)-10e3f;
		dfSignal[i]=1e8*cos(512*(6.283185307179586476925286*i)/1024)+1e7*sin((6.283185307179586476925286*i)/1024)-50e3;
		Signal[i]=static_cast<int>(1e8*cos(512*(6.283185307179586476925286*i)/1024)+1e7*sin((6.283185307179586476925286*i)/1024)-0e3);
	}
	//Signal[0]=1E9;
	//dfSignal[0]=1e9;
	for (size_t i=0;i<1024;i++)
	{
		//RefSplitSignal[i].re = cos(512*(6.283185307179586476925286*(2*i))/1024)+sin(3*(6.283185307179586476925286*(2*i))/1024)+2;
		//RefSplitSignal[i].im = cos(512*(6.283185307179586476925286*(2*i+1))/1024)+sin(3*(6.283185307179586476925286*(2*i+1))/1024)+2;
		//RealTransform[i].re = 2e9*cos(512*(6.283185307179586476925286*(2*i))/1024)+4e8*sin(3*(6.283185307179586476925286*(2*i))/1024);
		//RealTransform[i].re = 2e9*cos(512*(6.283185307179586476925286*(2*i))/1024)+4e8*sin(3*(6.283185307179586476925286*(2*i))/1024);
		RefSplitSignal[i].re = (10*cos(6.283185307179586476925286*i/1024*512));
		RefSplitSignal[i].im = 0;
	}

	/*m_FFT.CalculateFFT(Signal,Transform);
	m_FFT.CalculateFFT(Signal,TransformScaling,SCALEINPOUT,&Exponent);
	mf_FFT.CalculateFFT(fSignal,fTransform);
	mdf_FFT.CalculateFFT(dfSignal,dfTransform);*/
	//mf_FFT.CalculateFFT(fSignal,fTransform);
	//mf_FFT.CalculateIFFT(fTransform,fBackTransform);
	mrdf_FFT.CalculateRealFFT(dfSignal,dfrTransform);
	mrdf_FFT.CalculateRealIFFT(dfrTransform,dfrTransform);
	//CalculateDFT(1024, RefSplitSignal, Reference);

	mdf_FFT.CalculateRealFFT(dfSignal,dfTransform);
	mdf_FFT.CalculateRealIFFT(dfTransform,reinterpret_cast<DComplex*>(resultd));

	CFFT<10, double, double, DComplex> dummy;
	mf_FFT.CalculateRealFFT(fSignal,fTransform);
	mf_FFT.CalculateRealIFFT(fTransform, reinterpret_cast<FComplex*>(resultf));
	std::vector<double> real;
	std::vector<double> imag;
	for (size_t i = 0; i < 1024; ++i)
	{
		real.push_back(fTransform[i].re);
		imag.push_back(fTransform[i].im);
	}
	//real[5] = 1;

	plot<double>(real, imag);
	//plot<double>(real);
	size_t Exponent1, Exponent2;
	m_FFT.CalculateRealFFT(Signal,TransformScaling,SCALEINPOUT,&Exponent1);
	m_FFT.CalculateRealIFFT(TransformScaling,Transform,SCALEINPOUT,&Exponent2);

	size_t exponent;
	m_FFT2.CalculateRealFFT(Signal,RealTransform,SCALEINPOUT,&exponent);

	DComplex test_signal[1024];
	DComplex test_result[1024];
	size_t index = 0;
	for (DComplex& a : test_signal)
	{
		a.re = fSignal[index++];
		a.im = 0;
	}
	CalculateDFT(1024, test_signal, test_result);

	//dummy.Convert2HalfDFT(Reference);
	return 0;
}
