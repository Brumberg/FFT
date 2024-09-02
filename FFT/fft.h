#ifndef FFT_H
#define FFT_H

#include "fft_internals.h"
#include <limits>

#define _USETWIDDLE_
#define _RECURRENCE_
#define _C1x_

enum class EN_ScalingMethod { SCALE_INPUT, SCALE_OUTPUT, SCALEINPOUT };

/**
* fixed point implementation of the FFT algorithm.
* The fixed point version(s) return a value with twice the size of the input operant to avoid loss of precision
*/

template<size_t loglen, typename TW, typename T, typename TFFT> class CFFT {
private:
	TW      ma_TwiddleFactors[1 << (loglen - 1)];

	void    CreateTwiddleFactors(size_t length, TW twiddlefactors[], double scaling);

	TW      GetCosFactor(size_t stage, size_t stageshift);
	TW      GetSinFactor(size_t stage, size_t stageshift);


	template <size_t bitsize> size_t swap(size_t pointer)
	{
		static const size_t size = bitsize >> 1;//for an example for 11 size is 5
		static const size_t lowermask = ~(static_cast<size_t>(static_cast<size_t>(-1) >> size) << size);
		static const size_t uppermask = lowermask << (bitsize - size);
		static const size_t mask = ~(static_cast<size_t>(static_cast<size_t>(-1) >> bitsize) << bitsize);
		static const size_t remaining = (~(lowermask | uppermask)) & mask;
		pointer = (swap<size>(lowermask & pointer) << (bitsize - size)) | swap<size>((uppermask & pointer) >> (bitsize - size)) | (remaining & pointer);
		return pointer;
	}

	template <> size_t swap<2>(size_t pointer)
	{
		return ((pointer << 1) | (pointer >> 1)) & 3;
	}

	template <> size_t swap<1>(size_t pointer)
	{
		return pointer & 1;
	}

	template <size_t bitsize> size_t bitreversal(size_t pointer)
	{
		static const size_t size = (bitsize) >> 1;//for an example for 11 size is 5
		static const size_t lowermask = ~(static_cast<size_t>(static_cast<size_t>(-1) >> size) << size);
		static const size_t uppermask = lowermask << (bitsize - size);
		static const size_t mask = ~(static_cast<size_t>(static_cast<size_t>(-1) >> (bitsize)) << (bitsize));
		static const size_t remaining = (~(lowermask | uppermask)) & mask;
		pointer = (swap<size>(lowermask & pointer) << (bitsize - size)) | swap<size>((uppermask & pointer) >> (bitsize - size)) | (remaining & pointer);
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
	CFFT() { CreateTwiddleFactors(sizeof(ma_TwiddleFactors) / sizeof(ma_TwiddleFactors[0]), ma_TwiddleFactors, std::numeric_limits<TW>::max()); }//other scaling factors require explitcit treatment of halflength, index 0 and length in the  Convert2HalfDFT and Convert2ComplexDFT functions
	void CalculateFFT(T buffer[], TFFT ffttransform[]);
	void CalculateFFT(T buffer[], TFFT ffttransform[], EN_ScalingMethod rescale, size_t* pBlockExponent);
	void CalculateIFFT(TFFT buffer[], TFFT ffttransform[]);
	void CalculateIFFT(TFFT buffer[], TFFT ffttransform[], EN_ScalingMethod rescale, size_t* pBlockExponent);
	void CalculateRealFFT(T buffer[], TFFT ffttransform[]);
	void CalculateRealFFT(T buffer[], TFFT ffttransform[], EN_ScalingMethod rescale, size_t* pBlockExponent);
	void CalculateRealIFFT(TFFT buffer[], TFFT ffttransform[]);
	void CalculateRealIFFT(TFFT buffer[], TFFT ffttransform[], EN_ScalingMethod rescale, size_t* pBlockExponent);
};

template<size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::CreateTwiddleFactors(size_t length, TW twiddlefactors[], double scaling)
{
	static const double PI2 = 6.283185307179586476925286;

	static const size_t len = 1 << loglen;
	if (2 * length >= len)
	{
		const size_t len2 = len >> 1;
		for (size_t i = 0; i < len2; i++)
		{
			twiddlefactors[i] = static_cast<TW>(floor(scaling * cos(PI2 * static_cast<double>(i) / len) + .5));
		}
	}
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::ReorderSamples(T in[], TFFT data[])
{
	static const size_t length = 1 << loglen;
	for (size_t i = 0; i < length; i++)
	{
		size_t reordered = bitreversal<loglen>(i);
		data[reordered].re = in[i];
		data[reordered].im = 0;
	}
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::ReorderSamples(TFFT in[], TFFT data[])
{
	static const size_t length = 1 << loglen;
	for (size_t i = 0; i < length; i++)
	{
		size_t reordered = bitreversal<loglen>(i);
		data[reordered] = in[i];
	}
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::ReorderRealSamples(T in[], TFFT data[])
{
	static const size_t length = 1 << (loglen - 1);
	for (size_t i = 0; i < length; i++)
	{
		size_t reordered = bitreversal<loglen - 1>(i);
		data[reordered].re = in[(i << 1)];
		data[reordered].im = in[(i << 1) + 1];
	}
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::ReorderRealSamplesInPlace(TFFT in[])
{
	static const size_t length = 1 << (loglen - 1);
	for (size_t i = 0; i < length; i++)
	{
		size_t reordered = bitreversal<loglen - 1>(i);
		if (i >= reordered)
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
	static const size_t len = 1 << loglen;
	//more than .5 is not possible
	//size_t k = len*stage/stagesize;
	size_t k = (len * stage) >> stageshift;
	return ma_TwiddleFactors[k];
}

template <size_t loglen, typename TW, typename T, typename TFFT>
TW CFFT<loglen, TW, T, TFFT>::GetSinFactor(size_t stage, size_t stageshift)
{
	static const size_t len = 1 << loglen;
	//const size_t len2        = 1<<(loglen-1);
	static const size_t len4 = 1 << (loglen - 2);
	static const size_t mask = ~(((-1) >> (loglen)) << (loglen - 1));
	//more than .5 is not possible
	//size_t k = (len*stage/stagesize+(len>>2))&mask;
	size_t k = (((len * stage) >> stageshift) + (len >> 2)) & mask;
	return k < len4 ? ma_TwiddleFactors[k] : -ma_TwiddleFactors[k];
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::CalculateFFT(T buffer[], TFFT ffttransform[])
{
	static const size_t bitcorrection = sizeof(reinterpret_cast<TFFT*>(0)->re) * 8 - 1;
	static const size_t wscaling_real = (sizeof(reinterpret_cast<TFFT*>(0)->re) - sizeof(TW)) * 8;
	static const size_t wscaling_imag = (sizeof(reinterpret_cast<TFFT*>(0)->im) - sizeof(TW)) * 8;

	static const size_t ldlen = loglen;
	static const size_t length = 1 << loglen;
	ReorderSamples(buffer, ffttransform);

	for (size_t i = 0; i < ldlen; i++)
	{
		size_t stagesize = static_cast<size_t>(1) << (i + 1);
		size_t stagesize2 = static_cast<size_t>(1) << i;
		size_t looplength = length >> (i + 1);
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
				TFFT temp_a = ffttransform[butterflyix + R];
				TFFT temp_b = ffttransform[butterflyix + R + stagesize2];//load individual stages

				static const size_t req_downshifts = (sizeof(TFFT::OFWTYPE) - sizeof(T)) * 8u;
#ifdef _C1x_
				auto dummy = weight * temp_b;
				auto temp_a_scaled = temp_a << bitcorrection;
				ffttransform[butterflyix + R] = static_cast<TFFT>((temp_a_scaled + dummy) >> req_downshifts);
				ffttransform[butterflyix + R + stagesize2] = static_cast<TFFT>((temp_a_scaled - dummy) >> req_downshifts);
#else
				typename TFFT::ACCTYPE dummy = weight * temp_b;
				typename TFFT::ACCTYPE temp_a_scaled = temp_a << bitcorrection;
				ffttransform[butterflyix + R] = static_cast<TFFT>((temp_a_scaled + dummy) >> req_downshifts);
				ffttransform[butterflyix + R + stagesize2] = static_cast<TFFT>((temp_a_scaled - dummy) >> req_downshifts);
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
	static const size_t ldlen = loglen;
	static const size_t length = 1 << loglen;
	ReorderSamples(buffer, ffttransform);

	for (size_t i = 0; i < ldlen; i++)
	{
		size_t stagesize = static_cast<size_t>(1) << (i + 1);
		size_t stagesize2 = static_cast<size_t>(1) << i;
		size_t looplength = length >> (i + 1);
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
				
				static const size_t req_downshifts = (sizeof(TFFT::OFWTYPE) - sizeof(T)) * 8u - 1u;//scaling does not apply
#ifdef _C1x_
				auto dummy = weight * temp_b;
				auto temp_a_scaled = temp_a << bitcorrection;
				ffttransform[butterflyix + R] = static_cast<TFFT>((temp_a_scaled + dummy) >> req_downshifts);
				ffttransform[butterflyix + R + stagesize2] = static_cast<TFFT>((temp_a_scaled - dummy) >> req_downshifts);
#else
				typename TFFT::ACCTYPE dummy = weight * temp_b;
				typename TFFT::ACCTYPE temp_a_scaled = temp_a << bitcorrection;
				ffttransform[butterflyix + R] = static_cast<TFFT>((temp_a_scaled + dummy) >> req_downshifts);
				ffttransform[butterflyix + R + stagesize2] = static_cast<TFFT>((temp_a_scaled - dummy) >> req_downshifts);
#endif
			}
		}
	}
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::CalculateFFT(T buffer[], TFFT ffttransform[], EN_ScalingMethod rescale, size_t* pBlockExponent)
{
	static const size_t bitcorrection = sizeof(reinterpret_cast<TFFT*>(0)->re) * 8 - 1;
	static const size_t wscaling_real = (sizeof(reinterpret_cast<TFFT*>(0)->re) - sizeof(TW)) * 8;
	static const size_t wscaling_imag = (sizeof(reinterpret_cast<TFFT*>(0)->im) - sizeof(TW)) * 8;

	static const size_t ldlen = loglen;
	static const size_t length = 1 << loglen;
	size_t Exponent = 0;
	*pBlockExponent = loglen;
	ReorderSamples(buffer, ffttransform);

	if ((rescale == EN_ScalingMethod::SCALE_INPUT) || (rescale == EN_ScalingMethod::SCALEINPOUT))
	{
		T MaxVal = abs(buffer[0]);
		for (size_t i = 1; i < length; i++)
		{
			T dummy = abs(buffer[i]);
			if (dummy > MaxVal)
			{
				MaxVal = dummy;
			}
		}

		if (MaxVal != 0)
		{
			T mask = 0;
			mask = (~mask) << (sizeof(T) * 8 - 2);
			for (Exponent = 0; !(MaxVal & mask); MaxVal <<= 1) { Exponent++; }
			Exponent += (sizeof(reinterpret_cast<TFFT*>(0)->re) - sizeof(T)) * 8;

			//scale to the max is done in the butterfly
			/*for (size_t i=0;i<length;i++)
			{
				ffttransform[i].re <<= Exponent;
			}
			*pBlockExponent = Exponent;*/
		}
	}


	for (size_t i = 0; i < ldlen; i++)
	{
		size_t stagesize = 1 << (i + 1);
		size_t stagesize2 = 1 << i;
		size_t looplength = length >> (i + 1);
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
				TFFT temp_a = ffttransform[butterflyix + R];
				TFFT temp_b = ffttransform[butterflyix + R + stagesize2];//load individual stages
				temp_a.re <<= Exponent;
				temp_a.im <<= Exponent;
				temp_b.re <<= Exponent;
				temp_b.im <<= Exponent;

				static const size_t req_downshifts = (sizeof(TFFT::OFWTYPE) - sizeof(T)) * 8u;

#ifdef _C1x_
				auto dummy = weight * temp_b;
				auto temp_a_scaled = temp_a << bitcorrection;
				ffttransform[butterflyix + R] = static_cast<TFFT>((temp_a_scaled + dummy) >> req_downshifts);
				ffttransform[butterflyix + R + stagesize2] = static_cast<TFFT>((temp_a_scaled - dummy) >> req_downshifts);
#else
				typename TFFT::ACCTYPE dummy = weight * temp_b;
				typename TFFT::ACCTYPE temp_a_scaled = temp_a << bitcorrection;
				ffttransform[butterflyix + R] = static_cast<TFFT>((temp_a_scaled + dummy) >> req_downshifts);
				ffttransform[butterflyix + R + stagesize2] = static_cast<TFFT>((temp_a_scaled - dummy) >> req_downshifts);
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
		MaxOccVal = 0;
	}
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::CalculateIFFT(TFFT buffer[], TFFT ffttransform[], EN_ScalingMethod rescale, size_t* pBlockExponent)
{
	static const size_t bitcorrection = sizeof(reinterpret_cast<TFFT*>(0)->re) * 8 - 1;
	static const size_t wscaling_real = (sizeof(reinterpret_cast<TFFT*>(0)->re) - sizeof(TW)) * 8;
	static const size_t wscaling_imag = (sizeof(reinterpret_cast<TFFT*>(0)->im) - sizeof(TW)) * 8;

	static const size_t ldlen = loglen;
	static const size_t length = 1 << loglen;
	size_t Exponent = 0;
	*pBlockExponent = loglen;
	ReorderSamples(buffer, ffttransform);

	if ((rescale == EN_ScalingMethod::SCALE_INPUT) || (rescale == EN_ScalingMethod::SCALEINPOUT))
	{
		typename TFFT::DATATYPE MaxVal1 = abs(buffer[0].re);
		typename TFFT::DATATYPE MaxVal2 = abs(buffer[0].im);
		typename TFFT::DATATYPE MaxVal = MaxVal1 > MaxVal2 ? MaxVal1 : MaxVal2;
		for (size_t i = 1; i < length; i++)
		{
			typename TFFT::DATATYPE dummy1 = abs(buffer[i].re);
			typename TFFT::DATATYPE dummy2 = abs(buffer[i].im);
			typename TFFT::DATATYPE dummy = dummy1 > dummy2 ? dummy1 : dummy2;
			if (dummy > MaxVal)
			{
				MaxVal = dummy;
			}
		}

		if (MaxVal != 0)
		{
			T mask = 0;
			mask = (~mask) << (sizeof(T) * 8 - 2);
			for (Exponent = 0; !(MaxVal & mask); MaxVal <<= 1) { Exponent++; }
			Exponent += (sizeof(reinterpret_cast<TFFT*>(0)->re) - sizeof(T)) * 8;

			//scale to the max is done in the butterfly
			/*for (size_t i=0;i<length;i++)
			{
				ffttransform[i].re <<= Exponent;
			}
			*pBlockExponent = Exponent;*/
		}
	}


	for (size_t i = 0; i < ldlen; i++)
	{
		size_t stagesize = 1 << (i + 1);
		size_t stagesize2 = 1 << i;
		size_t looplength = length >> (i + 1);
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
				TFFT temp_a = ffttransform[butterflyix + R];
				TFFT temp_b = ffttransform[butterflyix + R + stagesize2];//load individual stages
				temp_a.re <<= Exponent;
				temp_a.im <<= Exponent;
				temp_b.re <<= Exponent;
				temp_b.im <<= Exponent;

				static const size_t req_downshifts = (sizeof(TFFT::OFWTYPE) - sizeof(T)) * 8u;
#ifdef _C1x_
				auto dummy = weight * temp_b;
				auto temp_a_scaled = temp_a << bitcorrection;
				ffttransform[butterflyix + R] = static_cast<TFFT>((temp_a_scaled + dummy) >> req_downshifts);
				ffttransform[butterflyix + R + stagesize2] = static_cast<TFFT>((temp_a_scaled - dummy) >> req_downshifts);
#else
				typename TFFT::ACCTYPE dummy = weight * temp_b;
				typename TFFT::ACCTYPE temp_a_scaled = temp_a << bitcorrection;
				ffttransform[butterflyix + R] = static_cast<TFFT>((temp_a_scaled + dummy) >> req_downshifts);
				ffttransform[butterflyix + R + stagesize2] = static_cast<TFFT>((temp_a_scaled - dummy) >> req_downshifts);
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
void CFFT<loglen, TW, T, TFFT>::CalculateRealFFT(T buffer[], TFFT ffttransform[])
{
	static const size_t bitcorrection = sizeof(reinterpret_cast<TFFT*>(0)->re) * 8 - 1;
	static const size_t wscaling_real = (sizeof(reinterpret_cast<TFFT*>(0)->re) - sizeof(TW)) * 8;
	static const size_t wscaling_imag = (sizeof(reinterpret_cast<TFFT*>(0)->im) - sizeof(TW)) * 8;

	static const size_t ldlen = loglen - 1;
	static const size_t length = 1 << (loglen - 1);
	ReorderRealSamples(buffer, ffttransform);

	for (size_t i = 0; i < ldlen; i++)
	{
		size_t stagesize = static_cast<size_t>(1) << (i + 1);
		size_t stagesize2 = static_cast<size_t>(1) << i;
		size_t looplength = length >> (i + 1);
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
				TFFT temp_a = ffttransform[butterflyix + R];
				TFFT temp_b = ffttransform[butterflyix + R + stagesize2];//load individual stages

				static const size_t req_downshifts = (sizeof(TFFT::OFWTYPE) - sizeof(T)) * 8u;

#ifdef _C1x_
				auto dummy = weight * temp_b;
				auto temp_a_scaled = temp_a << bitcorrection;
				ffttransform[butterflyix + R] = static_cast<TFFT>((temp_a_scaled + dummy) >> req_downshifts);
				ffttransform[butterflyix + R + stagesize2] = static_cast<TFFT>((temp_a_scaled - dummy) >> req_downshifts);
#else
				typename TFFT::ACCTYPE dummy = weight * temp_b;
				typename TFFT::ACCTYPE temp_a_scaled = temp_a << bitcorrection;
				ffttransform[butterflyix + R] = static_cast<TFFT>((temp_a_scaled + dummy) >> req_downshifts);
				ffttransform[butterflyix + R + stagesize2] = static_cast<TFFT>((temp_a_scaled - dummy) >> req_downshifts);
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

	static const size_t ldlen = loglen - 1;
	static const size_t length = 1 << (loglen - 1);

	memcpy(ffttransform, buffer, length * sizeof(TFFT));
	Convert2ComplexDFT(ffttransform);
	ReorderRealSamplesInPlace(ffttransform);
	//ReorderRealSamples(buffer, ffttransform);


	for (size_t i = 0; i < ldlen; i++)
	{
		size_t stagesize = static_cast<size_t>(1) << (i + 1);
		size_t stagesize2 = static_cast<size_t>(1) << i;
		size_t looplength = length >> (i + 1);
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
				
				static const size_t req_downshifts = (sizeof(TFFT::OFWTYPE) - sizeof(T)) * 8u - 1u;

#ifdef _C1x_
				auto dummy = weight * temp_b;
				auto temp_a_scaled = temp_a << bitcorrection;
				ffttransform[butterflyix + R] = static_cast<TFFT>((temp_a_scaled + dummy) >> req_downshifts);
				ffttransform[butterflyix + R + stagesize2] = static_cast<TFFT>((temp_a_scaled - dummy) >> req_downshifts);
#else
				typename TFFT::ACCTYPE dummy = weight * temp_b;
				typename TFFT::ACCTYPE temp_a_scaled = temp_a << bitcorrection;
				ffttransform[butterflyix + R] = static_cast<TFFT>((temp_a_scaled + dummy) >> req_downshifts);
				ffttransform[butterflyix + R + stagesize2] = static_cast<TFFT>((temp_a_scaled - dummy) >> req_downshifts);
#endif
			}
		}
	}
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::CalculateRealFFT(T buffer[], TFFT ffttransform[], EN_ScalingMethod rescale, size_t* pBlockExponent)
{
	static const size_t bitcorrection = sizeof(reinterpret_cast<TFFT*>(0)->re) * 8 - 1;
	static const size_t wscaling_real = (sizeof(reinterpret_cast<TFFT*>(0)->re) - sizeof(TW)) * 8;
	static const size_t wscaling_imag = (sizeof(reinterpret_cast<TFFT*>(0)->im) - sizeof(TW)) * 8;

	static const size_t ldlen = loglen - 1;
	static const size_t length = 1 << (loglen - 1);
	size_t Exponent = 0;
	*pBlockExponent = loglen;
	ReorderRealSamples(buffer, ffttransform);

	if ((rescale == EN_ScalingMethod::SCALE_INPUT) || (rescale == EN_ScalingMethod::SCALEINPOUT))
	{
		T MaxVal = abs(buffer[0]);
		for (size_t i = 1; i < length; i++)
		{
			T dummy = abs(buffer[i]);
			if (dummy > MaxVal)
			{
				MaxVal = dummy;
			}
		}

		if (MaxVal != 0)
		{
			T mask = 0;
			mask = (~mask) << (sizeof(T) * 8 - 2);
			for (Exponent = 0; !(MaxVal & mask); MaxVal <<= 1) { Exponent++; }
			Exponent += (sizeof(reinterpret_cast<TFFT*>(0)->re) - sizeof(T)) * 8;
		}
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
			weight.im = -GetSinFactor(R, (i + 1));//-sin(PI2*RT_ratio);

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

				static const size_t req_downshifts = (sizeof(TFFT::OFWTYPE) - sizeof(T)) * 8u;

#ifdef _C1x_
				auto dummy = weight * temp_b;
				auto temp_a_scaled = temp_a << bitcorrection;
				ffttransform[butterflyix + R] = static_cast<TFFT>((temp_a_scaled + dummy) >> req_downshifts);
				ffttransform[butterflyix + R + stagesize2] = static_cast<TFFT>((temp_a_scaled - dummy) >> req_downshifts);
#else
				typename TFFT::ACCTYPE dummy = weight * temp_b;
				typename TFFT::ACCTYPE temp_a_scaled = temp_a << bitcorrection;
				ffttransform[butterflyix + R] = static_cast<TFFT>((temp_a_scaled + dummy) >> req_downshifts);
				ffttransform[butterflyix + R + stagesize2] = static_cast<TFFT>((temp_a_scaled - dummy) >> req_downshifts);
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
	Convert2HalfDFT(ffttransform, Exponent);
	*pBlockExponent -= Exponent;
}


template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::CalculateRealIFFT(TFFT buffer[], TFFT ffttransform[], EN_ScalingMethod rescale, size_t* pBlockExponent)
{
	static const size_t bitcorrection = sizeof(reinterpret_cast<TFFT*>(0)->re) * 8 - 1;
	static const size_t wscaling_real = (sizeof(reinterpret_cast<TFFT*>(0)->re) - sizeof(TW)) * 8;
	static const size_t wscaling_imag = (sizeof(reinterpret_cast<TFFT*>(0)->im) - sizeof(TW)) * 8;

	static const size_t ldlen = loglen - 1;
	static const size_t length = 1 << (loglen - 1);
	size_t Exponent = 0;
	*pBlockExponent = loglen;

	memcpy(ffttransform, buffer, length * sizeof(TFFT));

	if ((rescale == EN_ScalingMethod::SCALE_INPUT) || (rescale == EN_ScalingMethod::SCALEINPOUT))
	{
		typename TFFT::DATATYPE MaxVal1 = abs(ffttransform[0].re);
		typename TFFT::DATATYPE MaxVal2 = abs(ffttransform[0].im);
		typename TFFT::DATATYPE MaxVal = MaxVal1 < MaxVal2 ? MaxVal2 : MaxVal1;
		for (size_t i = 1; i < length; i++)
		{
			typename TFFT::DATATYPE MaxVal1 = abs(ffttransform[i].re);
			typename TFFT::DATATYPE MaxVal2 = abs(ffttransform[i].im);
			typename TFFT::DATATYPE dummy = MaxVal1 < MaxVal2 ? MaxVal2 : MaxVal1;
			if (dummy > MaxVal)
			{
				MaxVal = dummy;
			}
		}

		if (MaxVal != 0)
		{
			typename TFFT::DATATYPE mask = 0;
			mask = (~mask) << (sizeof(T) * 8 - 2);
			for (Exponent = 0; !(MaxVal & mask); MaxVal <<= 1) { Exponent++; }
		}
	}
	Convert2ComplexDFT(ffttransform, Exponent);
	ReorderRealSamplesInPlace(ffttransform);
	*pBlockExponent -= Exponent;

	if (rescale == EN_ScalingMethod::SCALEINPOUT)
	{
		typename TFFT::DATATYPE MaxVal1 = abs(ffttransform[0].re);
		typename TFFT::DATATYPE MaxVal2 = abs(ffttransform[0].im);
		typename TFFT::DATATYPE MaxVal = MaxVal1 < MaxVal2 ? MaxVal2 : MaxVal1;
		for (size_t i = 1; i < length; i++)
		{
			typename TFFT::DATATYPE MaxVal1 = abs(ffttransform[i].re);
			typename TFFT::DATATYPE MaxVal2 = abs(ffttransform[i].im);
			typename TFFT::DATATYPE dummy = MaxVal1 < MaxVal2 ? MaxVal2 : MaxVal1;
			if (dummy > MaxVal)
			{
				MaxVal = dummy;
			}
		}

		if (MaxVal != 0)
		{
			typename TFFT::DATATYPE mask = 0;
			mask = (~mask) << (sizeof(T) * 8 - 2);
			for (Exponent = 0; !(MaxVal & mask); MaxVal <<= 1) { Exponent++; }
		}
	}
	else
	{
		Exponent = 0;
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

				static const size_t req_downshifts = (sizeof(TFFT::OFWTYPE) - sizeof(T)) * 8u;

#ifdef _C1x_
				auto dummy = weight * temp_b;
				auto temp_a_scaled = temp_a << bitcorrection;
				ffttransform[butterflyix + R] = static_cast<TFFT>((temp_a_scaled + dummy) >> req_downshifts);
				ffttransform[butterflyix + R + stagesize2] = static_cast<TFFT>((temp_a_scaled - dummy) >> req_downshifts);
#else
				typename TFFT::ACCTYPE dummy = weight * temp_b;
				typename TFFT::ACCTYPE temp_a_scaled = temp_a << bitcorrection;
				ffttransform[butterflyix + R] = static_cast<TFFT>((temp_a_scaled + dummy) >> req_downshifts);
				ffttransform[butterflyix + R + stagesize2] = static_cast<TFFT>((temp_a_scaled - dummy) >> req_downshifts);
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
	static const float PI = 3.1415926535897932384626434f;

	static const size_t length = 1 << (loglen - 1);
	static const size_t redlength = length;
	static const size_t halflength = length >> 1;
	for (size_t i = 1; i < halflength; i++)
	{
		//here the FB case
		typename TFFT::OFWTYPE H1R = (static_cast<typename TFFT::OFWTYPE>(transform[i].re) + transform[redlength - i].re) >> 1;
		typename TFFT::OFWTYPE H1I = (static_cast<typename TFFT::OFWTYPE>(transform[i].im) - transform[redlength - i].im) >> 1;//because conj complex
		typename TFFT::OFWTYPE H2R = (static_cast<typename TFFT::OFWTYPE>(transform[i].re) - transform[redlength - i].re) >> 1;//that is the second term
		typename TFFT::OFWTYPE H2I = (static_cast<typename TFFT::OFWTYPE>(transform[i].im) + transform[redlength - i].im) >> 1;//adding, because conj complex

		//to be checked

		TW cosine = GetCosFactor(i, loglen);//cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		TW sine = GetSinFactor(i, loglen);//sin(PI*i/length);

		typename TFFT::OFWTYPE res1i = (H2R * cosine + H2I * sine) >> (sizeof(TW) * 8 - 1);
		typename TFFT::OFWTYPE res1r = -((H2I * cosine - H2R * sine) >> (sizeof(TW) * 8 - 1));//complex i is integrated

		typename TFFT::OFWTYPE resreal = (H1R - res1r);//0.5f*(H1R-res1r);
		typename TFFT::OFWTYPE resimag = (H1I - res1i);//0.5f*(H1I-res1i);->already scaled see shift operation above

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

		res1i = (H2R * cosine + H2I * sine) >> (sizeof(TW) * 8 - 1);
		res1r = (H2I * cosine - H2R * sine) >> (sizeof(TW) * 8 - 1);//complex i is integrated

		typename TFFT::OFWTYPE resreal2 = (H1R - res1r);//0.5f*(H1R-res1r);
		typename TFFT::OFWTYPE resimag2 = (-H1I - res1i);//0.5f*(-H1I-res1i);->already scaled see shift operation above

		transform[i].re = static_cast<typename TFFT::DATATYPE>(resreal >> 1);
		transform[i].im = static_cast<typename TFFT::DATATYPE>(resimag >> 1);
		transform[length - i].re = static_cast<typename TFFT::DATATYPE>(resreal2 >> 1);
		transform[length - i].im = static_cast<typename TFFT::DATATYPE>(resimag2 >> 1);
	}
	typename TFFT::OFWTYPE lval = static_cast<typename TFFT::OFWTYPE>(transform[0].re) + transform[0].im;
	typename TFFT::OFWTYPE hval = static_cast<typename TFFT::OFWTYPE>(transform[0].re) - transform[0].im;
	transform[0].re = static_cast<typename TFFT::DATATYPE>(lval >> 1);
	transform[0].im = static_cast<typename TFFT::DATATYPE>(hval >> 1);
	transform[halflength].re = transform[halflength].re/2;
	transform[halflength].im = -transform[halflength].im/2;//treatment of the half index
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::Convert2ComplexDFT(TFFT transform[])
{
	static const float PI2 = 6.2831853071795864769252868f;
	static const float PI = 3.1415926535897932384626434f;

	static const size_t length = 1 << (loglen - 1);
	static const size_t redlength = length;
	static const size_t halflength = length >> 1;
	for (size_t i = 1; i < halflength; i++)
	{
		//here the FB case
		typename TFFT::OFWTYPE H1R = (static_cast<typename TFFT::OFWTYPE>(transform[i].re) + transform[redlength - i].re) >> 1;
		typename TFFT::OFWTYPE H1I = (static_cast<typename TFFT::OFWTYPE>(transform[i].im) - transform[redlength - i].im) >> 1;//because conj complex
		typename TFFT::OFWTYPE H2R = (static_cast<typename TFFT::OFWTYPE>(transform[i].re) - transform[redlength - i].re) >> 1;//that is the second term
		typename TFFT::OFWTYPE H2I = (static_cast<typename TFFT::OFWTYPE>(transform[i].im) + transform[redlength - i].im) >> 1;//adding, because conj complex

		TW cosine = GetCosFactor(i, loglen);//cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		TW sine = GetSinFactor(i, loglen);//sin(PI*i/length);

		//float res1i	= H2R*cosine-H2I*sine;
		//float res1r	= -(H2I*cosine+H2R*sine);//complex i is integrated
		typename TFFT::OFWTYPE res1i = (H2R * cosine - H2I * sine) >> (sizeof(TW) * 8 - 1);
		typename TFFT::OFWTYPE res1r = -((H2I * cosine + H2R * sine) >> (sizeof(TW) * 8 - 1));//complex i is integrated


		//float resreal	= 0.5f*(H1R+res1r);
		//float resimag	= 0.5f*(H1I+res1i);
		typename TFFT::OFWTYPE resreal = (H1R + res1r);//0.5f*(H1R+res1r);
		typename TFFT::OFWTYPE resimag = (H1I + res1i);//0.5f*(H1I+res1i);->already scaled see shift operation above

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
		res1i = (-H2R * cosine + H2I * sine) >> (sizeof(TW) * 8 - 1);
		res1r = (H2I * cosine + H2R * sine) >> (sizeof(TW) * 8 - 1);//complex i is integrated

		//float resreal2			= 0.5f*(H1R+res1r);
		//float resimag2			= -0.5f*(H1I+res1i);
		typename TFFT::OFWTYPE resreal2 = (H1R + res1r);//0.5f*(H1R+res1r);
		typename TFFT::OFWTYPE resimag2 = -(H1I + res1i);//-0.5f*(H1I+res1i);->already scaled see shift operation above

		//transform[i].re			= resreal;
		//transform[i].im			= resimag;
		//transform[length-i].re	= resreal2;
		//transform[length-i].im	= resimag2;
		transform[i].re = static_cast<typename TFFT::DATATYPE>(resreal);//>>1
		transform[i].im = static_cast<typename TFFT::DATATYPE>(resimag);//>>1
		transform[length - i].re = static_cast<typename TFFT::DATATYPE>(resreal2);//>>1
		transform[length - i].im = static_cast<typename TFFT::DATATYPE>(resimag2);//>>1
	}
	//float lval					= transform[0].re+transform[0].im;
	//float hval					= transform[0].re-transform[0].im;
	typename TFFT::OFWTYPE lval = static_cast<typename TFFT::OFWTYPE>(transform[0].re) + transform[0].im;
	typename TFFT::OFWTYPE hval = static_cast<typename TFFT::OFWTYPE>(transform[0].re) - transform[0].im;
	//transform[0].re				= 0.5f*lval;
	//transform[0].im				= 0.5f*hval;
	transform[0].re = static_cast<typename TFFT::DATATYPE>(lval >> 1);
	transform[0].im = static_cast<typename TFFT::DATATYPE>(hval >> 1);
	transform[halflength].im = -transform[halflength].im;//treatment of the half index
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::Convert2ComplexDFT(TFFT transform[], size_t exponent)
{
	static const float PI2 = 6.2831853071795864769252868f;
	static const float PI = 3.1415926535897932384626434f;

	static const size_t length = 1 << (loglen - 1);
	static const size_t redlength = length;
	static const size_t halflength = length >> 1;
	for (size_t i = 1; i < halflength; i++)
	{
		//here the FB case
		typename TFFT::OFWTYPE H1R = ((static_cast<typename TFFT::OFWTYPE>(transform[i].re) + transform[redlength - i].re) << exponent) >> 1;
		typename TFFT::OFWTYPE H1I = ((static_cast<typename TFFT::OFWTYPE>(transform[i].im) - transform[redlength - i].im) << exponent) >> 1;//because conj complex
		typename TFFT::OFWTYPE H2R = ((static_cast<typename TFFT::OFWTYPE>(transform[i].re) - transform[redlength - i].re) << exponent) >> 1;//that is the second term
		typename TFFT::OFWTYPE H2I = ((static_cast<typename TFFT::OFWTYPE>(transform[i].im) + transform[redlength - i].im) << exponent) >> 1;//adding, because conj complex

		TW cosine = GetCosFactor(i, loglen);//cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		TW sine = GetSinFactor(i, loglen);//sin(PI*i/length);

		//float res1i	= H2R*cosine-H2I*sine;
		//float res1r	= -(H2I*cosine+H2R*sine);//complex i is integrated
		typename TFFT::OFWTYPE res1i = (H2R * cosine - H2I * sine) >> (sizeof(TW) * 8 - 1);
		typename TFFT::OFWTYPE res1r = -((H2I * cosine + H2R * sine) >> (sizeof(TW) * 8 - 1));//complex i is integrated


		//float resreal	= 0.5f*(H1R+res1r);
		//float resimag	= 0.5f*(H1I+res1i);
		typename TFFT::OFWTYPE resreal = (H1R + res1r);//0.5f*(H1R+res1r);
		typename TFFT::OFWTYPE resimag = (H1I + res1i);//0.5f*(H1I+res1i);->already scaled see shift operation above

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
		res1i = (-H2R * cosine + H2I * sine) >> (sizeof(TW) * 8 - 1);
		res1r = (H2I * cosine + H2R * sine) >> (sizeof(TW) * 8 - 1);//complex i is integrated

		//float resreal2			= 0.5f*(H1R+res1r);
		//float resimag2			= -0.5f*(H1I+res1i);
		typename TFFT::OFWTYPE resreal2 = (H1R + res1r);//0.5f*(H1R+res1r);
		typename TFFT::OFWTYPE resimag2 = -(H1I + res1i);//-0.5f*(H1I+res1i);->already scaled see shift operation above

		//transform[i].re			= resreal;
		//transform[i].im			= resimag;
		//transform[length-i].re	= resreal2;
		//transform[length-i].im	= resimag2;
		transform[i].re = static_cast<typename TFFT::DATATYPE>(resreal);//>>1
		transform[i].im = static_cast<typename TFFT::DATATYPE>(resimag);//>>1
		transform[length - i].re = static_cast<typename TFFT::DATATYPE>(resreal2);//>>1
		transform[length - i].im = static_cast<typename TFFT::DATATYPE>(resimag2);//>>1
	}
	//float lval					= transform[0].re+transform[0].im;
	//float hval					= transform[0].re-transform[0].im;
	typename TFFT::OFWTYPE lval;//		= ((static_cast<typename TFFT::OFWTYPE>(transform[0].re)+transform[0].im)>>1)<<(exponent);
	typename TFFT::OFWTYPE hval;//		= ((static_cast<typename TFFT::OFWTYPE>(transform[0].re)-transform[0].im)>>1)<<(exponent);
	if (exponent)
	{
		lval = ((static_cast<typename TFFT::OFWTYPE>(transform[0].re) + transform[0].im)) << (exponent - 1);
		hval = ((static_cast<typename TFFT::OFWTYPE>(transform[0].re) - transform[0].im)) << (exponent - 1);
	}
	else
	{
		lval = ((static_cast<typename TFFT::OFWTYPE>(transform[0].re) + transform[0].im)) >> 1;
		hval = ((static_cast<typename TFFT::OFWTYPE>(transform[0].re) - transform[0].im)) >> 1;
	}
	//transform[0].re				= 0.5f*lval;
	//transform[0].im				= 0.5f*hval;
	transform[0].re = static_cast<typename TFFT::DATATYPE>(lval);
	transform[0].im = static_cast<typename TFFT::DATATYPE>(hval);//>>1 moved to exponent stuff ^
	transform[halflength].re = (transform[halflength].re << exponent);//treatment of the half index
	transform[halflength].im = -(transform[halflength].im << exponent);//treatment of the half index
}

template <size_t loglen, typename TW, typename T, typename TFFT>
void CFFT<loglen, TW, T, TFFT>::Convert2HalfDFT(TFFT transform[], size_t exponent)
{
	static const float PI = 3.1415926535897932384626434f;

	static const size_t length = 1 << (loglen - 1);
	static const size_t redlength = length;
	static const size_t halflength = length >> 1;
	for (size_t i = 1; i < halflength; i++)
	{
		//here the FB case
		typename TFFT::OFWTYPE H1R = ((static_cast<typename TFFT::OFWTYPE>(transform[i].re) + transform[redlength - i].re) << exponent) >> 1;
		typename TFFT::OFWTYPE H1I = ((static_cast<typename TFFT::OFWTYPE>(transform[i].im) - transform[redlength - i].im) << exponent) >> 1;//because conj complex
		typename TFFT::OFWTYPE H2R = ((static_cast<typename TFFT::OFWTYPE>(transform[i].re) - transform[redlength - i].re) << exponent) >> 1;//that is the second term
		typename TFFT::OFWTYPE H2I = ((static_cast<typename TFFT::OFWTYPE>(transform[i].im) + transform[redlength - i].im) << exponent) >> 1;//adding, because conj complex

		//to be checked

		TW cosine = GetCosFactor(i, loglen);//cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		TW sine = GetSinFactor(i, loglen);//sin(PI*i/length);

		typename TFFT::OFWTYPE res1i = (H2R * cosine + H2I * sine) >> (sizeof(TW) * 8 - 1);
		typename TFFT::OFWTYPE res1r = -((H2I * cosine - H2R * sine) >> (sizeof(TW) * 8 - 1));//complex i is integrated

		typename TFFT::OFWTYPE resreal = (H1R - res1r);//0.5f*(H1R-res1r);
		typename TFFT::OFWTYPE resimag = (H1I - res1i);//0.5f*(H1I-res1i);->already scaled see shift operation above

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

		res1i = (H2R * cosine + H2I * sine) >> (sizeof(TW) * 8 - 1);
		res1r = (H2I * cosine - H2R * sine) >> (sizeof(TW) * 8 - 1);//complex i is integrated

		typename TFFT::OFWTYPE resreal2 = (H1R - res1r);//0.5f*(H1R-res1r);
		typename TFFT::OFWTYPE resimag2 = (-H1I - res1i);//0.5f*(-H1I-res1i);->already scaled see shift operation above

		transform[i].re = static_cast<typename TFFT::DATATYPE>(resreal >> 1);
		transform[i].im = static_cast<typename TFFT::DATATYPE>(resimag >> 1);
		transform[length - i].re = static_cast<typename TFFT::DATATYPE>(resreal2 >> 1);
		transform[length - i].im = static_cast<typename TFFT::DATATYPE>(resimag2 >> 1);
	}
	//typename TFFT::OFWTYPE lval			= ((static_cast<typename TFFT::OFWTYPE>(transform[0].re)+transform[0].im)>>1)<<exponent;
	//typename TFFT::OFWTYPE hval			= ((static_cast<typename TFFT::OFWTYPE>(transform[0].re)-transform[0].im)>>1)<<exponent;

	typename TFFT::OFWTYPE lval;
	typename TFFT::OFWTYPE hval;
	typename TFFT::OFWTYPE midr;
	typename TFFT::OFWTYPE midi;
	if (exponent)
	{
		lval = ((static_cast<typename TFFT::OFWTYPE>(transform[0].re) + transform[0].im)) << (exponent - 1);
		hval = ((static_cast<typename TFFT::OFWTYPE>(transform[0].re) - transform[0].im)) << (exponent - 1);
		midr = static_cast<typename TFFT::OFWTYPE>(transform[halflength].re) << (exponent - 1);
		midi = static_cast<typename TFFT::OFWTYPE>(transform[halflength].im) << (exponent - 1);
	}
	else
	{
		lval = ((static_cast<typename TFFT::OFWTYPE>(transform[0].re) + transform[0].im)) >> 1;
		hval = ((static_cast<typename TFFT::OFWTYPE>(transform[0].re) - transform[0].im)) >> 1;
		midr = static_cast<typename TFFT::OFWTYPE>(transform[halflength].re) >> 1;
		midi = static_cast<typename TFFT::OFWTYPE>(transform[halflength].im) >> 1;
	}

	transform[0].re = static_cast<typename TFFT::DATATYPE>(lval);//>>1 moved to exponent above
	transform[0].im = static_cast<typename TFFT::DATATYPE>(hval);//>>1 moved to exponent above
	transform[halflength].re = static_cast<typename TFFT::DATATYPE>(midr);//treatment of the half index
	transform[halflength].im = -static_cast<typename TFFT::DATATYPE>(midi);//treatment of the half index
}

/**
* floating point implementation of the FFT algorithm (single precision).
* The floating point version(s) return a value with twice the size of the input operant to avoid loss of precision
*/

template<size_t loglen> class CFFT<loglen, float, float, FFT_INTERNALS::f32_complex> {
private:
#ifdef _USETWIDDLE_
	float	ma_TwiddleFactors[1 << (loglen - 1)];

	void    CreateTwiddleFactors(size_t length, float twiddlefactors[]);
	float   GetCosFactor(size_t stage, size_t stageshift);
	float   GetSinFactor(size_t stage, size_t stageshift);
#endif



	template <size_t bitsize> size_t swap(size_t pointer)
	{
		static const size_t size = bitsize >> 1;//for an example for 11 size is 5
		static const size_t lowermask = ~(static_cast<size_t>(static_cast<size_t>(-1) >> size) << size);
		static const size_t uppermask = lowermask << (bitsize - size);
		static const size_t mask = ~(static_cast<size_t>(static_cast<size_t>(-1) >> bitsize) << bitsize);
		static const size_t remaining = (~(lowermask | uppermask)) & mask;
		pointer = (swap<size>(lowermask & pointer) << (bitsize - size)) | swap<size>((uppermask & pointer) >> (bitsize - size)) | (remaining & pointer);
		return pointer;
	}

	template <> size_t swap<2>(size_t pointer)
	{
		return ((pointer << 1) | (pointer >> 1)) & 3;
	}

	template <> size_t swap<1>(size_t pointer)
	{
		return pointer & 1;
	}

	template <size_t bitsize> size_t bitreversal(size_t pointer)
	{
		const size_t size = (bitsize) >> 1;//for an example for 11 size is 5
		const size_t lowermask = ~(static_cast<size_t>(static_cast<size_t>(-1) >> size) << size);
		const size_t uppermask = lowermask << (bitsize - size);
		const size_t mask = ~(static_cast<size_t>(static_cast<size_t>(-1) >> (bitsize)) << (bitsize));
		const size_t remaining = (~(lowermask | uppermask)) & mask;
		pointer = (swap<size>(lowermask & pointer) << (bitsize - size)) | swap<size>((uppermask & pointer) >> (bitsize - size)) | (remaining & pointer);
		return pointer;
	}

	void ReorderSamples(float in[], FFT_INTERNALS::f32_complex data[]);
	void ReorderSamples(FFT_INTERNALS::f32_complex in[], FFT_INTERNALS::f32_complex data[]);
	void ReorderRealSamples(float in[], FFT_INTERNALS::f32_complex data[]);
	void ReorderRealSamplesInPlace(FFT_INTERNALS::f32_complex in[]);
	void Convert2HalfDFT(FFT_INTERNALS::f32_complex transform[]);
	void Convert2ComplexDFT(FFT_INTERNALS::f32_complex transform[]);
public:
#ifdef _USETWIDDLE_
	CFFT() { CreateTwiddleFactors(sizeof(ma_TwiddleFactors) / sizeof(ma_TwiddleFactors[0]), ma_TwiddleFactors); }
#else
	CFFT() {};
#endif
	void CalculateFFT(float buffer[], FFT_INTERNALS::f32_complex ffttransform[]);
	void CalculateIFFT(FFT_INTERNALS::f32_complex buffer[], FFT_INTERNALS::f32_complex ffttransform[]);
	void CalculateRealFFT(float buffer[], FFT_INTERNALS::f32_complex ffttransform[]);
	void CalculateRealIFFT(FFT_INTERNALS::f32_complex buffer[], FFT_INTERNALS::f32_complex ffttransform[]);
};

#ifdef _USETWIDDLE_
template<size_t loglen>
void CFFT<loglen, float, float, FFT_INTERNALS::f32_complex>::CreateTwiddleFactors(size_t length, float twiddlefactors[])
{
	static const float PI2 = 6.283185307179586476925286f;

	static const size_t len = 1 << loglen;
	if (2 * length >= len)
	{
		const size_t len2 = len >> 1;
		for (size_t i = 0; i < len2; i++)
		{
			twiddlefactors[i] = static_cast<float>(cos(PI2 * static_cast<float>(i) / len));
		}
	}
}
#endif

template <size_t loglen>
void CFFT<loglen, float, float, FFT_INTERNALS::f32_complex>::ReorderSamples(float in[], FFT_INTERNALS::f32_complex data[])
{
	static const size_t length = 1 << loglen;
	for (size_t i = 0; i < length; i++)
	{
		size_t reordered = bitreversal<loglen>(i);
		data[reordered].re = in[i];
		data[reordered].im = 0;
	}
}

template <size_t loglen>
void CFFT<loglen, float, float, FFT_INTERNALS::f32_complex>::ReorderSamples(FFT_INTERNALS::f32_complex in[], FFT_INTERNALS::f32_complex data[])
{
	static const size_t length = 1 << loglen;
	for (size_t i = 0; i < length; i++)
	{
		size_t reordered = bitreversal<loglen>(i);
		data[reordered].re = in[i].re;
		data[reordered].im = in[i].im;
	}
}

template <size_t loglen>
void CFFT<loglen, float, float, FFT_INTERNALS::f32_complex>::ReorderRealSamples(float in[], FFT_INTERNALS::f32_complex data[])
{
	static const size_t length = 1 << (loglen - 1);
	for (size_t i = 0; i < length; i++)
	{
		size_t reordered = bitreversal<loglen - 1>(i);
		data[reordered].re = in[(i << 1)];
		data[reordered].im = in[(i << 1) + 1];
	}
}

template <size_t loglen>
void CFFT<loglen, float, float, FFT_INTERNALS::f32_complex>::ReorderRealSamplesInPlace(FFT_INTERNALS::f32_complex in[])
{
	static const size_t length = 1 << (loglen - 1);
	for (size_t i = 0; i < length; i++)
	{
		size_t reordered = bitreversal<loglen - 1>(i);
		if (i >= reordered)
		{
			FFT_INTERNALS::f32_complex dummy;
			dummy = in[reordered];
			in[reordered] = in[i];
			in[i] = dummy;
		}
	}
}

#ifdef _USETWIDDLE_
template <size_t loglen>
float CFFT<loglen, float, float, FFT_INTERNALS::f32_complex>::GetCosFactor(size_t stage, size_t stageshift)
{
	static const size_t len = 1 << loglen;
	//more than .5 is not possible
	//size_t k = len*stage/stagesize;
	size_t k = (len * stage) >> stageshift;
	return ma_TwiddleFactors[k];
}

template <size_t loglen>
float CFFT<loglen, float, float, FFT_INTERNALS::f32_complex>::GetSinFactor(size_t stage, size_t stageshift)
{
	static const size_t len = 1 << loglen;
	//const size_t len2        = 1<<(loglen-1);
	static const size_t len4 = 1 << (loglen - 2);
	static const size_t mask = ~(((-1) >> (loglen)) << (loglen - 1));
	//more than .5 is not possible
	//size_t k = (len*stage/stagesize+(len>>2))&mask;
	size_t k = (((len * stage) >> stageshift) + (len >> 2)) & mask;
	return k < len4 ? ma_TwiddleFactors[k] : -ma_TwiddleFactors[k];
}
#endif

template <size_t loglen>
void CFFT<loglen, float, float, FFT_INTERNALS::f32_complex>::CalculateFFT(float buffer[], FFT_INTERNALS::f32_complex ffttransform[])
{
	static const float PI2 = 6.2831853071795864769252868f;
	static const size_t ldlen = loglen;
	static const size_t length = 1 << loglen;
	ReorderSamples(buffer, ffttransform);

	for (size_t i = 0; i < ldlen; i++)
	{
		size_t stagesize = static_cast<size_t>(1) << (i + 1);
		size_t stagesize2 = static_cast<size_t>(1) << i;
		size_t looplength = length >> (i + 1);

#ifndef _USETWIDDLE_
#ifdef _RECURRENCE_
		const float delta = PI2 / stagesize;
		const float delta_2 = 0.5f * delta;
		const float sin_delta = sin(delta_2);
		const float alpha = 2.0f * sin_delta * sin_delta;
		const float beta = sin(delta);
		float cos_rec = 1.0f;
		float sin_rec = 0.0f;
#endif
#endif

		for (size_t R = 0; R < stagesize2; R++)
		{
			FFT_INTERNALS::f32_complex weight;
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
				FFT_INTERNALS::f32_complex temp_a = ffttransform[butterflyix + R];
				FFT_INTERNALS::f32_complex temp_b = ffttransform[butterflyix + R + stagesize2];//load individual stages
#ifdef _C1x_
				auto dummy = weight * temp_b;
				ffttransform[butterflyix + R] = temp_a + dummy;
				ffttransform[butterflyix + R + stagesize2] = temp_a - dummy;
#else
				FFT_INTERNALS::f32_complex::ACCTYPE dummy = weight * temp_b;
				ffttransform[butterflyix + R] = temp_a + dummy;
				ffttransform[butterflyix + R + stagesize2] = temp_a - dummy;
#endif
			}
		}
	}
}

template <size_t loglen>
void CFFT<loglen, float, float, FFT_INTERNALS::f32_complex>::CalculateIFFT(FFT_INTERNALS::f32_complex buffer[], FFT_INTERNALS::f32_complex ffttransform[])
{
	static const float PI2 = 6.2831853071795864769252868f;
	static const size_t ldlen = loglen;
	static const size_t length = 1 << loglen;
	ReorderSamples(buffer, ffttransform);

	for (size_t i = 0; i < ldlen; i++)
	{
		size_t stagesize = static_cast<size_t>(1) << (i + 1);
		size_t stagesize2 = static_cast<size_t>(1) << i;
		size_t looplength = length >> (i + 1);

#ifndef _USETWIDDLE_
#ifdef _RECURRENCE_
		const float delta = PI2 / stagesize;
		const float delta_2 = 0.5f * delta;
		const float sin_delta = sin(delta_2);
		const float alpha = 2.0f * sin_delta * sin_delta;
		const float beta = sin(delta);
		float cos_rec = 1.0f;
		float sin_rec = 0.0f;
#endif
#endif

		for (size_t R = 0; R < stagesize2; R++)
		{
			FFT_INTERNALS::f32_complex weight;
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
				FFT_INTERNALS::f32_complex temp_a = ffttransform[butterflyix + R];
				FFT_INTERNALS::f32_complex temp_b = ffttransform[butterflyix + R + stagesize2];//load individual stages

#ifdef _C1x_
				auto dummy = weight * temp_b;
				ffttransform[butterflyix + R] = temp_a + dummy;
				ffttransform[butterflyix + R + stagesize2] = temp_a - dummy;
#else
				FFT_INTERNALS::f32_complex::ACCTYPE dummy = weight * temp_b;
				ffttransform[butterflyix + R] = temp_a + dummy;
				ffttransform[butterflyix + R + stagesize2] = temp_a - dummy;
#endif
			}
		}
	}

	const float scaling = 1.0f / length;
	for (size_t i = 0; i < length; i++)
	{
		ffttransform[i].re *= scaling;
		ffttransform[i].im *= scaling;
	}
}

template <size_t loglen>
void CFFT<loglen, float, float, FFT_INTERNALS::f32_complex>::CalculateRealFFT(float buffer[], FFT_INTERNALS::f32_complex ffttransform[])
{
	static const size_t ldlen = loglen - 1;
	static const size_t length = 1 << (loglen - 1);
	ReorderRealSamples(buffer, ffttransform);

	for (size_t i = 0; i < ldlen; i++)
	{
#ifndef _USETWIDDLE_
		size_t stagesize = static_cast<size_t>(1) << (i + 1);
#endif
		size_t stagesize2 = static_cast<size_t>(1) << i;
		size_t looplength = length >> (i + 1);
		for (size_t R = 0; R < stagesize2; R++)
		{
			FFT_INTERNALS::f32_complex weight;
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
				FFT_INTERNALS::f32_complex temp_a = ffttransform[butterflyix + R];
				FFT_INTERNALS::f32_complex temp_b = ffttransform[butterflyix + R + stagesize2];//load individual stages

#ifdef _C1x_
				auto dummy = weight * temp_b;
				ffttransform[butterflyix + R] = temp_a + dummy;
				ffttransform[butterflyix + R + stagesize2] = temp_a - dummy;
#else
				FFT_INTERNALS::f32_complex::ACCTYPE dummy = weight * temp_b;
				ffttransform[butterflyix + R] = temp_a + dummy;
				ffttransform[butterflyix + R + stagesize2] = temp_a - dummy;
#endif
			}
		}
	}
	Convert2HalfDFT(ffttransform);
}

template <size_t loglen>
void CFFT<loglen, float, float, FFT_INTERNALS::f32_complex>::CalculateRealIFFT(FFT_INTERNALS::f32_complex buffer[], FFT_INTERNALS::f32_complex ffttransform[])
{
	static const float PI2 = 6.2831853071795864769252868f;
	static const size_t ldlen = loglen - 1;
	static const size_t length = 1 << (loglen - 1);
	memcpy(ffttransform, buffer, length * sizeof(FFT_INTERNALS::f32_complex));
	Convert2ComplexDFT(ffttransform);
	ReorderRealSamplesInPlace(ffttransform);

	for (size_t i = 0; i < ldlen; i++)
	{
		size_t stagesize = static_cast<size_t>(1) << (i + 1);
		size_t stagesize2 = static_cast<size_t>(1) << i;
		size_t looplength = length >> (i + 1);

#ifndef _USETWIDDLE_
#ifdef _RECURRENCE_
		const float delta = PI2 / stagesize;
		const float delta_2 = 0.5f * delta;
		const float sin_delta = sin(delta_2);
		const float alpha = 2.0f * sin_delta * sin_delta;
		const float beta = sin(delta);
		float cos_rec = 1.0f;
		float sin_rec = 0.0f;
#endif
#endif

		for (size_t R = 0; R < stagesize2; R++)
		{
			FFT_INTERNALS::f32_complex weight;
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
				FFT_INTERNALS::f32_complex temp_a = ffttransform[butterflyix + R];
				FFT_INTERNALS::f32_complex temp_b = ffttransform[butterflyix + R + stagesize2];//load individual stages

#ifdef _C1x_
				auto dummy = weight * temp_b;
				ffttransform[butterflyix + R] = temp_a + dummy;
				ffttransform[butterflyix + R + stagesize2] = temp_a - dummy;
#else
				FFT_INTERNALS::f32_complex::ACCTYPE dummy = weight * temp_b;
				ffttransform[butterflyix + R] = temp_a + dummy;
				ffttransform[butterflyix + R + stagesize2] = temp_a - dummy;
#endif
			}
		}
	}

	const float scaling = 1.0f / length;
	for (size_t i = 0; i < length; i++)
	{
		ffttransform[i].re *= scaling;
		ffttransform[i].im *= scaling;
	}
}

template <size_t loglen>
void CFFT<loglen, float, float, FFT_INTERNALS::f32_complex>::Convert2HalfDFT(FFT_INTERNALS::f32_complex transform[])
{
	static const float PI2 = 6.2831853071795864769252868f;
	static const float PI = 3.1415926535897932384626434f;

	size_t length = 1 << (loglen - 1);
	size_t redlength = length;
	size_t halflength = length >> 1;

#ifndef _USETWIDDLE_
#ifdef _RECURRENCE_
	const float	delta = PI / length;
	const float	delta_2 = 0.5f * delta;
	const float	sin_delta = sin(delta_2);
	const float	alpha = 2.0f * sin_delta * sin_delta;
	const float	beta = sin(delta);
	float		 sineiter = 0;
	float		 cosineiter = 1.0f;
#endif
#endif

	for (size_t i = 1; i < halflength; i++)
	{
		//here the FB case
		float H1R = transform[i].re + transform[redlength - i].re;
		float H1I = transform[i].im - transform[redlength - i].im;//because conj complex
		float H2R = transform[i].re - transform[redlength - i].re;//that is the second term
		float H2I = transform[i].im + transform[redlength - i].im;//adding, because conj complex

#ifdef _USETWIDDLE_
		float cosine = GetCosFactor(i, loglen);//cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		float sine = GetSinFactor(i, loglen);//sin(PI*i/length);
#else
#ifdef _RECURRENCE_

		float cosine = cosineiter;
		float sine = sineiter;
		cosineiter = cosine - (alpha * cosine + beta * sine);
		sineiter = sineiter - (alpha * sine - beta * cosine);
		cosine = cosineiter;
		sine = sineiter;
#else
		float cosine = cos(PI * i / length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		float sine = sin(PI * i / length);
#endif
#endif

		float res1i = H2R * cosine + H2I * sine;
		float res1r = -(H2I * cosine - H2R * sine);//complex i is integrated

		float resreal = 0.5f * (H1R - res1r);
		float resimag = 0.5f * (H1I - res1i);

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

		res1i = H2R * cosine + H2I * sine;
		res1r = (H2I * cosine - H2R * sine);//complex i is integrated

		float resreal2 = 0.5f * (H1R - res1r);
		float resimag2 = 0.5f * (-H1I - res1i);

		transform[i].re = resreal;
		transform[i].im = resimag;
		transform[length - i].re = resreal2;
		transform[length - i].im = resimag2;
	}
	float lval = transform[0].re + transform[0].im;
	float hval = transform[0].re - transform[0].im;
	transform[0].re = lval;
	transform[0].im = hval;
	transform[halflength].im = -transform[halflength].im;//treatment of the half index
}

template <size_t loglen>
void CFFT<loglen, float, float, FFT_INTERNALS::f32_complex>::Convert2ComplexDFT(FFT_INTERNALS::f32_complex transform[])
{
	static const float PI2 = 6.2831853071795864769252868f;
	static const float PI = 3.1415926535897932384626434f;

	size_t length = 1 << (loglen - 1);
	size_t redlength = length;
	size_t halflength = length >> 1;
#ifndef _USETWIDDLE_
#ifdef _RECURRENCE_
	const float	delta = PI / length;
	const float	delta_2 = 0.5f * delta;
	const float	sin_delta = sin(delta_2);
	const float	alpha = 2.0f * sin_delta * sin_delta;
	const float	beta = sin(delta);
	float		 sineiter = 0.0f;
	float		 cosineiter = 1.0f;
#endif
#endif
	for (size_t i = 1; i < halflength; i++)
	{
		//here the FB case
		float H1R = transform[i].re + transform[redlength - i].re;
		float H1I = transform[i].im - transform[redlength - i].im;//because conj complex
		float H2R = transform[i].re - transform[redlength - i].re;//that is the second term
		float H2I = transform[i].im + transform[redlength - i].im;//adding, because conj complex

#ifdef _USETWIDDLE_
		float cosine = GetCosFactor(i, loglen);//cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		float sine = GetSinFactor(i, loglen);//sin(PI*i/length);
#else
#ifdef _RECURRENCE_

		float cosine = cosineiter;
		float sine = sineiter;
		cosineiter = cosine - (alpha * cosine + beta * sine);
		sineiter = sineiter - (alpha * sine - beta * cosine);
		cosine = cosineiter;
		sine = sineiter;
#else
		float cosine = cos(PI * i / length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		float sine = sin(PI * i / length);
#endif
#endif

		float res1i = H2R * cosine - H2I * sine;
		float res1r = -(H2I * cosine + H2R * sine);//complex i is integrated

		float resreal = 0.5f * (H1R + res1r);
		float resimag = 0.5f * (H1I + res1i);

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

		res1i = -H2R * cosine + H2I * sine;
		res1r = (H2I * cosine + H2R * sine);//complex i is integrated

		float resreal2 = 0.5f * (H1R + res1r);
		float resimag2 = -0.5f * (H1I + res1i);

		transform[i].re = resreal;
		transform[i].im = resimag;
		transform[length - i].re = resreal2;
		transform[length - i].im = resimag2;
	}
	float lval = transform[0].re + transform[0].im;
	float hval = transform[0].re - transform[0].im;
	transform[0].re = 0.5f * lval;
	transform[0].im = 0.5f * hval;
	transform[halflength].im = -transform[halflength].im;//treatment of the half index
}


/**
* floating point implementation of the FFT algorithm (single precision).
* The floating point version(s) return a value with twice the size of the input operant to avoid loss of precision
*/

template<size_t loglen> class CFFT<loglen, double, double, FFT_INTERNALS::f64_complex> {
private:
#ifdef _USETWIDDLE_
	double	ma_TwiddleFactors[1 << (loglen - 1)];

	void    CreateTwiddleFactors(size_t length, double twiddlefactors[]);

	double	GetCosFactor(size_t stage, size_t stageshift);
	double	GetSinFactor(size_t stage, size_t stageshift);
#endif

	template <size_t bitsize> size_t swap(size_t pointer)
	{
		static const size_t size = bitsize >> 1;//for an example for 11 size is 5
		static const size_t lowermask = ~(static_cast<size_t>(static_cast<size_t>(-1) >> size) << size);
		static const size_t uppermask = lowermask << (bitsize - size);
		static const size_t mask = ~(static_cast<size_t>(static_cast<size_t>(-1) >> bitsize) << bitsize);
		static const size_t remaining = (~(lowermask | uppermask)) & mask;
		pointer = (swap<size>(lowermask & pointer) << (bitsize - size)) | swap<size>((uppermask & pointer) >> (bitsize - size)) | (remaining & pointer);
		return pointer;
	}

	template <> size_t swap<2>(size_t pointer)
	{
		return ((pointer << 1) | (pointer >> 1)) & 3;
	}

	template <> size_t swap<1>(size_t pointer)
	{
		return pointer & 1;
	}

	template <size_t bitsize> size_t bitreversal(size_t pointer)
	{
		static const size_t size = (bitsize) >> 1;//for an example for 11 size is 5
		static const size_t lowermask = ~(static_cast<size_t>(static_cast<size_t>(-1) >> size) << size);
		static const size_t uppermask = lowermask << (bitsize - size);
		static const size_t mask = ~(static_cast<size_t>(static_cast<size_t>(-1) >> (bitsize)) << (bitsize));
		static const size_t remaining = (~(lowermask | uppermask)) & mask;
		pointer = (swap<size>(lowermask & pointer) << (bitsize - size)) | swap<size>((uppermask & pointer) >> (bitsize - size)) | (remaining & pointer);
		return pointer;
	}

	void ReorderSamples(double in[], FFT_INTERNALS::f64_complex data[]);
	void ReorderSamples(FFT_INTERNALS::f64_complex in[], FFT_INTERNALS::f64_complex data[]);
	void ReorderRealSamples(double in[], FFT_INTERNALS::f64_complex data[]);
	void ReorderRealSamplesInPlace(FFT_INTERNALS::f64_complex in[]);
	void Convert2HalfDFT(FFT_INTERNALS::f64_complex ffttransform[]);
	void Convert2ComplexDFT(FFT_INTERNALS::f64_complex transform[]);
public:
#ifdef _DEBUG_RFFT_
	void Convert2HalfDFT(FFT_INTERNALS::f64_complex ffttransform[], FFT_INTERNALS::f64_complex result[]);
#endif
#ifdef _USETWIDDLE_
	CFFT() { CreateTwiddleFactors(sizeof(ma_TwiddleFactors) / sizeof(ma_TwiddleFactors[0]), ma_TwiddleFactors); }
#else
	CFFT() {};
#endif
	void CalculateFFT(double buffer[], FFT_INTERNALS::f64_complex ffttransform[]);
	void CalculateIFFT(FFT_INTERNALS::f64_complex buffer[], FFT_INTERNALS::f64_complex ffttransform[]);
	void CalculateRealFFT(double buffer[], FFT_INTERNALS::f64_complex ffttransform[]);
	void CalculateRealIFFT(FFT_INTERNALS::f64_complex buffer[], FFT_INTERNALS::f64_complex ffttransform[]);
};

#ifdef _USETWIDDLE_
template<size_t loglen>
void CFFT<loglen, double, double, FFT_INTERNALS::f64_complex>::CreateTwiddleFactors(size_t length, double twiddlefactors[])
{
	static const double PI2 = 6.283185307179586476925286;

	static const size_t len = 1 << loglen;
	if (2 * length >= len)
	{
		const size_t len2 = len >> 1;
		for (size_t i = 0; i < len2; i++)
		{
			twiddlefactors[i] = static_cast<double>(cos(PI2 * static_cast<double>(i) / len));
		}
	}
}
#endif

template <size_t loglen>
void CFFT<loglen, double, double, FFT_INTERNALS::f64_complex>::ReorderSamples(double in[], FFT_INTERNALS::f64_complex data[])
{
	static const size_t length = 1 << loglen;
	for (size_t i = 0; i < length; i++)
	{
		size_t reordered = bitreversal<loglen>(i);
		data[reordered].re = in[i];
		data[reordered].im = 0;
	}
}

template <size_t loglen>
void CFFT<loglen, double, double, FFT_INTERNALS::f64_complex>::ReorderSamples(FFT_INTERNALS::f64_complex in[], FFT_INTERNALS::f64_complex data[])
{
	static const size_t length = 1 << loglen;
	for (size_t i = 0; i < length; i++)
	{
		size_t reordered = bitreversal<loglen>(i);
		data[reordered].re = in[i].re;
		data[reordered].im = in[i].im;
	}
}

template <size_t loglen>
void CFFT<loglen, double, double, FFT_INTERNALS::f64_complex>::ReorderRealSamples(double in[], FFT_INTERNALS::f64_complex data[])
{
	static const size_t length = 1 << (loglen - 1);
	for (size_t i = 0; i < length; i++)
	{
		size_t reordered = bitreversal<loglen - 1>(i);
		data[reordered].re = in[(i << 1)];
		data[reordered].im = in[(i << 1) + 1];
	}
}

template <size_t loglen>
void CFFT<loglen, double, double, FFT_INTERNALS::f64_complex>::ReorderRealSamplesInPlace(FFT_INTERNALS::f64_complex in[])
{
	static const size_t length = 1 << (loglen - 1);
	for (size_t i = 0; i < length; i++)
	{
		size_t reordered = bitreversal<loglen - 1>(i);
		if (i >= reordered)
		{
			FFT_INTERNALS::f64_complex dummy;
			dummy = in[reordered];
			in[reordered] = in[i];
			in[i] = dummy;
		}
	}
}

#ifdef _USETWIDDLE_
template <size_t loglen>
double CFFT<loglen, double, double, FFT_INTERNALS::f64_complex>::GetCosFactor(size_t stage, size_t stageshift)
{
	static const size_t len = 1 << loglen;
	//more than .5 is not possible
	//size_t k = len*stage/stagesize;
	size_t k = (len * stage) >> stageshift;
	return ma_TwiddleFactors[k];
}

template <size_t loglen>
double CFFT<loglen, double, double, FFT_INTERNALS::f64_complex>::GetSinFactor(size_t stage, size_t stageshift)
{
	static const size_t len = 1 << loglen;
	//static const size_t len2        = 1<<(loglen-1);
	static const size_t len4 = 1 << (loglen - 2);
	static const size_t mask = ~(((-1) >> (loglen)) << (loglen - 1));
	//more than .5 is not possible
	//size_t k = (len*stage/stagesize+(len>>2))&mask;
	size_t k = (((len * stage) >> stageshift) + (len >> 2)) & mask;
	return k < len4 ? ma_TwiddleFactors[k] : -ma_TwiddleFactors[k];
}
#endif

template <size_t loglen>
void CFFT<loglen, double, double, FFT_INTERNALS::f64_complex>::CalculateFFT(double buffer[], FFT_INTERNALS::f64_complex ffttransform[])
{
	static const double PI2 = 6.2831853071795864769252868;
	static const size_t ldlen = loglen;
	static const size_t length = 1 << loglen;
	ReorderSamples(buffer, ffttransform);

	for (size_t i = 0; i < ldlen; i++)
	{
		size_t stagesize = static_cast<size_t>(1) << (i + 1);
		size_t stagesize2 = static_cast<size_t>(1) << i;
		size_t looplength = length >> (i + 1);

#ifndef _USETWIDDLE_
#ifdef _RECURRENCE_
		const double delta = PI2 / stagesize;
		const double delta_2 = 0.5 * delta;
		const double sin_delta = sin(delta_2);
		const double alpha = 2.0 * sin_delta * sin_delta;
		const double beta = sin(delta);
		double cos_rec = 1.0f;
		double sin_rec = 0.0f;
#endif
#endif
		for (size_t R = 0; R < stagesize2; R++)
		{
			FFT_INTERNALS::f64_complex weight;
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
				FFT_INTERNALS::f64_complex temp_a = ffttransform[butterflyix + R];
				FFT_INTERNALS::f64_complex temp_b = ffttransform[butterflyix + R + stagesize2];//load individual stages

#ifdef _C1x_
				auto dummy = weight * temp_b;
				ffttransform[butterflyix + R] = temp_a + dummy;
				ffttransform[butterflyix + R + stagesize2] = temp_a - dummy;
#else
				FFT_INTERNALS::f64_complex::ACCTYPE dummy = weight * temp_b;
				ffttransform[butterflyix + R] = temp_a + dummy;
				ffttransform[butterflyix + R + stagesize2] = temp_a - dummy;
#endif
			}
		}
	}
}

template <size_t loglen>
void CFFT<loglen, double, double, FFT_INTERNALS::f64_complex>::CalculateIFFT(FFT_INTERNALS::f64_complex buffer[], FFT_INTERNALS::f64_complex ffttransform[])
{
	static const double PI2 = 6.2831853071795864769252868;
	static const size_t ldlen = loglen;
	static const size_t length = 1 << loglen;
	ReorderSamples(buffer, ffttransform);

	for (size_t i = 0; i < ldlen; i++)
	{
		size_t stagesize = static_cast<size_t>(1) << (i + 1);
		size_t stagesize2 = static_cast<size_t>(1) << i;
		size_t looplength = length >> (i + 1);

#ifndef _USETWIDDLE_
#ifdef _RECURRENCE_
		const double delta = PI2 / stagesize;
		const double delta_2 = 0.5 * delta;
		const double sin_delta = sin(delta_2);
		const double alpha = 2.0 * sin_delta * sin_delta;
		const double beta = sin(delta);
		double cos_rec = 1.0f;
		double sin_rec = 0.0f;
#endif
#endif
		for (size_t R = 0; R < stagesize2; R++)
		{
			FFT_INTERNALS::f64_complex weight;
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
				FFT_INTERNALS::f64_complex temp_a = ffttransform[butterflyix + R];
				FFT_INTERNALS::f64_complex temp_b = ffttransform[butterflyix + R + stagesize2];//load individual stages
#ifdef _C1x_
				auto dummy = weight * temp_b;
				ffttransform[butterflyix + R] = temp_a + dummy;
				ffttransform[butterflyix + R + stagesize2] = temp_a - dummy;
#else
				FFT_INTERNALS::f64_complex::ACCTYPE dummy = weight * temp_b;
				ffttransform[butterflyix + R] = temp_a + dummy;
				ffttransform[butterflyix + R + stagesize2] = temp_a - dummy;
#endif
			}
		}
	}

	const double scaling = 1.0 / length;
	for (size_t i = 0; i < length; i++)
	{
		ffttransform[i].re *= scaling;
		ffttransform[i].im *= scaling;
	}
}

template <size_t loglen>
void CFFT<loglen, double, double, FFT_INTERNALS::f64_complex>::CalculateRealIFFT(FFT_INTERNALS::f64_complex buffer[], FFT_INTERNALS::f64_complex ffttransform[])
{
	static const double PI2 = 6.2831853071795864769252868;
	static const size_t ldlen = loglen - 1;
	static const size_t length = 1 << (loglen - 1);
	memcpy(ffttransform, buffer, length * sizeof(FFT_INTERNALS::f64_complex));
	Convert2ComplexDFT(ffttransform);
	ReorderRealSamplesInPlace(ffttransform);

	for (size_t i = 0; i < ldlen; i++)
	{
		size_t stagesize = static_cast<size_t>(1) << static_cast<size_t>(i + 1);
		size_t stagesize2 = static_cast<size_t>(1) << i;
		size_t looplength = length >> (i + 1);

#ifndef _USETWIDDLE_
#ifdef _RECURRENCE_
		const double delta = PI2 / stagesize;
		const double delta_2 = 0.5 * delta;
		const double sin_delta = sin(delta_2);
		const double alpha = 2.0 * sin_delta * sin_delta;
		const double beta = sin(delta);
		double cos_rec = 1.0;
		double sin_rec = 0.0;
#endif
#endif
		for (size_t R = 0; R < stagesize2; R++)
		{
			FFT_INTERNALS::f64_complex weight;
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
				FFT_INTERNALS::f64_complex temp_a = ffttransform[butterflyix + R];
				FFT_INTERNALS::f64_complex temp_b = ffttransform[butterflyix + R + stagesize2];//load individual stages
#ifdef _C1x_
				auto dummy = weight * temp_b;
				ffttransform[butterflyix + R] = temp_a + dummy;
				ffttransform[butterflyix + R + stagesize2] = temp_a - dummy;
#else
				FFT_INTERNALS::f64_complex::ACCTYPE dummy = weight * temp_b;
				ffttransform[butterflyix + R] = temp_a + dummy;
				ffttransform[butterflyix + R + stagesize2] = temp_a - dummy;
#endif
			}
		}
	}

	const double scaling = 1.0 / length;
	for (size_t i = 0; i < length; i++)
	{
		ffttransform[i].re *= scaling;
		ffttransform[i].im *= scaling;
	}
}

template <size_t loglen>
void CFFT<loglen, double, double, FFT_INTERNALS::f64_complex>::CalculateRealFFT(double buffer[], FFT_INTERNALS::f64_complex ffttransform[])
{
	size_t ldlen = loglen - 1;
	size_t length = 1 << (loglen - 1);
	ReorderRealSamples(buffer, ffttransform);

	for (size_t i = 0; i < ldlen; i++)
	{
#ifndef _USETWIDDLE_
		size_t stagesize = static_cast<size_t>(1) << (i + 1);
#endif
		size_t stagesize2 = static_cast<size_t>(1) << i;
		size_t looplength = length >> (i + 1);
		for (size_t R = 0; R < stagesize2; R++)
		{
			FFT_INTERNALS::f64_complex weight;
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
				FFT_INTERNALS::f64_complex temp_a = ffttransform[butterflyix + R];
				FFT_INTERNALS::f64_complex temp_b = ffttransform[butterflyix + R + stagesize2];//load individual stages		
#ifdef _C1x_
				auto dummy = weight * temp_b;
				ffttransform[butterflyix + R] = temp_a + dummy;
				ffttransform[butterflyix + R + stagesize2] = temp_a - dummy;
#else
				FFT_INTERNALS::f64_complex::ACCTYPE dummy = weight * temp_b;
				ffttransform[butterflyix + R] = temp_a + dummy;
				ffttransform[butterflyix + R + stagesize2] = temp_a - dummy;
#endif
			}
		}
	}
	Convert2HalfDFT(ffttransform);
}

#ifdef _DEBUG_RFFT_
template <size_t loglen>
void CFFT<loglen, double, double, FFT_INTERNALS::f64_complex>::Convert2HalfDFT(FFT_INTERNALS::f64_complex transform[], FFT_INTERNALS::f64_complex result[])
{
	static const double PI2 = 6.2831853071795864769252868;
	static const double PI = 3.1415926535897932384626434;

	size_t length = 1 << (loglen - 1);
	size_t redlength = length;
	for (size_t i = 1; i < length; i++)
	{
		double H1R = transform[i].re + transform[redlength - i].re;
		double H1I = transform[i].im - transform[redlength - i].im;//because conj complex

		double H2R = transform[i].re - transform[redlength - i].re;//that is the second term
		double H2I = transform[i].im + transform[redlength - i].im;//adding, because conj complex

		double cosine = cos(PI * i / length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		double sine = sin(PI * i / length);

		double res1i = H2R * cosine + H2I * sine;
		double res1r = -(H2I * cosine - H2R * sine);//complex i is integrated

		double resreal = 0.5 * (H1R - res1r);
		double resimag = 0.5 * (H1I - res1i);

		result[i].re = resreal;
		result[i].im = resimag;
	}
	double lval = transform[0].re + transform[0].im;
	double hval = transform[0].re - transform[0].im;
	result[0].re = lval;
	result[0].im = hval;
}
#endif

template <size_t loglen>
void CFFT<loglen, double, double, FFT_INTERNALS::f64_complex>::Convert2HalfDFT(FFT_INTERNALS::f64_complex transform[])
{
	static const double PI2 = 6.2831853071795864769252868;
	static const double PI = 3.1415926535897932384626434;

	size_t length = 1 << (loglen - 1);
	size_t redlength = length;
	size_t halflength = length >> 1;

#ifndef _USETWIDDLE_
#ifdef _RECURRENCE_
	const double delta = PI / length;
	const double delta_2 = 0.5 * delta;
	const double sin_delta = sin(delta_2);
	const double alpha = 2.0 * sin_delta * sin_delta;
	const double beta = sin(delta);
	double		 sineiter = 0;
	double		 cosineiter = 1;
#endif
#endif

	for (size_t i = 1; i < halflength; i++)
	{
		//here the FB case
		double H1R = transform[i].re + transform[redlength - i].re;
		double H1I = transform[i].im - transform[redlength - i].im;//because conj complex
		double H2R = transform[i].re - transform[redlength - i].re;//that is the second term
		double H2I = transform[i].im + transform[redlength - i].im;//adding, because conj complex

#ifdef _USETWIDDLE_
		double cosine = GetCosFactor(i, loglen);//cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		double sine = GetSinFactor(i, loglen);//sin(PI*i/length);
#else
#ifdef _RECURRENCE_

		double cosine = cosineiter;
		double sine = sineiter;
		cosineiter = cosine - (alpha * cosine + beta * sine);
		sineiter = sineiter - (alpha * sine - beta * cosine);
		cosine = cosineiter;
		sine = sineiter;
#else
		double cosine = cos(PI * i / length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		double sine = sin(PI * i / length);
#endif
#endif

		double res1i = H2R * cosine + H2I * sine;
		double res1r = -(H2I * cosine - H2R * sine);//complex i is integrated

		double resreal = 0.5 * (H1R - res1r);
		double resimag = 0.5 * (H1I - res1i);

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

		res1i = H2R * cosine + H2I * sine;
		res1r = (H2I * cosine - H2R * sine);//complex i is integrated

		double resreal2 = 0.5 * (H1R - res1r);
		double resimag2 = 0.5 * (-H1I - res1i);

		transform[i].re = resreal;
		transform[i].im = resimag;
		transform[length - i].re = resreal2;
		transform[length - i].im = resimag2;
	}
	double lval = transform[0].re + transform[0].im;
	double hval = transform[0].re - transform[0].im;
	transform[0].re = lval;
	transform[0].im = hval;
	transform[halflength].im = -transform[halflength].im;//treatment of the half index
}

template <size_t loglen>
void CFFT<loglen, double, double, FFT_INTERNALS::f64_complex>::Convert2ComplexDFT(FFT_INTERNALS::f64_complex transform[])
{
	static const double PI2 = 6.2831853071795864769252868;
	static const double PI = 3.1415926535897932384626434;

	size_t length = 1 << (loglen - 1);
	size_t redlength = length;
	size_t halflength = length >> 1;

#ifndef _USETWIDDLE_
#ifdef _RECURRENCE_
	const double delta = PI / length;
	const double delta_2 = 0.5 * delta;
	const double sin_delta = sin(delta_2);
	const double alpha = 2.0 * sin_delta * sin_delta;
	const double beta = sin(delta);
	double		 sineiter = 0.;
	double		 cosineiter = 1.;
#endif
#endif

	for (size_t i = 1; i < halflength; i++)
	{
		//here the FB case
		double H1R = transform[i].re + transform[redlength - i].re;
		double H1I = transform[i].im - transform[redlength - i].im;//because conj complex
		double H2R = transform[i].re - transform[redlength - i].re;//that is the second term
		double H2I = transform[i].im + transform[redlength - i].im;//adding, because conj complex

#ifdef _USETWIDDLE_
		double cosine = GetCosFactor(i, loglen);//cos(PI*i/length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		double sine = GetSinFactor(i, loglen);//sin(PI*i/length);
#else
#ifdef _RECURRENCE_
		double cosine = cosineiter;
		double sine = sineiter;
		cosineiter = cosine - (alpha * cosine + beta * sine);
		sineiter = sineiter - (alpha * sine - beta * cosine);
		cosine = cosineiter;
		sine = sineiter;
#else
		double cosine = cos(PI * i / length);//length is reduced length, therefore we compensate with pi instead of 2*pi
		double sine = sin(PI * i / length);
#endif
#endif

		double res1i = H2R * cosine - H2I * sine;
		double res1r = -(H2I * cosine + H2R * sine);//complex i is integrated

		double resreal = 0.5f * (H1R + res1r);
		double resimag = 0.5f * (H1I + res1i);

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

		res1i = -H2R * cosine + H2I * sine;
		res1r = (H2I * cosine + H2R * sine);//complex i is integrated

		double resreal2 = 0.5f * (H1R + res1r);
		double resimag2 = -0.5f * (H1I + res1i);

		transform[i].re = resreal;
		transform[i].im = resimag;
		transform[length - i].re = resreal2;
		transform[length - i].im = resimag2;
	}

	double H1R = transform[0].re + transform[0].re;
	double H1I = transform[0].im - transform[0].im;//because conj complex
	double H2R = transform[0].re - transform[0].re;//that is the second term
	double H2I = transform[0].im + transform[0].im;//adding, because conj complex
	double cosine = cos(PI * 0 / length);//length is reduced length, therefore we compensate with pi instead of 2*pi
	double sine = sin(PI * 0 / length);

	double res1i = H2R * cosine - H2I * sine;
	double res1r = -(H2I * cosine + H2R * sine);//complex i is integrated

	double resreal = 0.5f * (H1R + res1r);
	double resimag = 0.5f * (H1I + res1i);


	double lval = transform[0].re + transform[0].im;
	double hval = transform[0].re - transform[0].im;
	transform[0].re = 0.5 * lval;
	transform[0].im = 0.5 * hval;
	transform[halflength].im = -transform[halflength].im;//treatment of the half index
}

#endif
