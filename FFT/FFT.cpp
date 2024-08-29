// FFT.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <math.h>
#include <assert.h>
#include <iostream>
#include "plot.h"
#include "fft_internals.h"
#include "fft.h"

#ifdef _GTEST_
#include "gtest/gtest.h"
#include "gmock/gmock.h"

using sint32 = int;

#include "DSPArray.h"
#endif

size_t Exponent;
static FFT_INTERNALS::i32_complex TransformScaling[1024];
static FFT_INTERNALS::i32_complex Transform[1024];
static FFT_INTERNALS::i32_complex RealTransform[1024];

static int TwiddleFactors[512];
static int Signal[1024];

static FFT_INTERNALS::f32_complex fTransform[1024];
static FFT_INTERNALS::f32_complex fBackTransform[1024];
static float fTwiddle[512];
static float fSignal[1024];

static FFT_INTERNALS::f64_complex dfTransform[1024];
static FFT_INTERNALS::f64_complex dfBackTransform[1024];
static FFT_INTERNALS::f64_complex dfrTransform[512];
static double dfTwiddle[512];
static double dfSignal[1024];
static FFT_INTERNALS::f64_complex Reference[1024];
static FFT_INTERNALS::f64_complex RefSplitSignal[1024];
static FFT_INTERNALS::f64_complex result[1024];
static float resultf[1024];
static double resultd[1024];


#ifdef _GTEST_
template <typename CType>
void CalculateDFT(size_t length, CType* buffer, CType* result)
{
	static const double PI2 = 6.2831853071795864769252868;
	for (size_t i = 0; i < length; i++)
	{
		FFT_INTERNALS::f64_complex accu;
		accu.re = 0;
		accu.im = 0;
		for (size_t j = 0; j < length; j++)
		{
			double cosine = cos((PI2 * j * i) / length);
			double sine = sin((PI2 * j * i) / length);
			accu.re += buffer[j].re * cosine + buffer[j].im * sine;
			accu.im += buffer[j].im * cosine - buffer[j].re * sine;
		}
		result[i].re = static_cast<decltype(CType::re)>(accu.re);
		result[i].im = static_cast<decltype(CType::im)>(accu.im);
	}
}

template <typename T, size_t N>
void GenerateSignalPairHannSym(CDSPArray<T, N>& up, CDSPArray<T, N>& dn,
	double SamplingFrequency, double BurstFrequency,
	double ToF, double DeltaT,
	double BufferDelay, double SignalDuration,
	double Amplitude, double Offset)
{
	static const double pi = 3.14159265358979324f;
	up.resize(N);
	dn.resize(N);
	up.fill(0);
	dn.fill(0);
	ToF -= BufferDelay;
	const double SamplesDelay = std::trunc(((ToF - abs(DeltaT / 2.f)) * SamplingFrequency));
	double intshift = std::trunc(DeltaT * SamplingFrequency);
	double fractshift = (DeltaT * SamplingFrequency) - intshift;
	const double SamplesWave = std::trunc(SignalDuration * SamplingFrequency);
	const double arg = 2.f * pi * BurstFrequency / SamplingFrequency;
	const double arg2 = 2.f * pi / (SignalDuration * SamplingFrequency);
	const size_t carrierlength_lo = std::min(static_cast<size_t>(SamplesWave), (N < (static_cast<size_t>(SamplesDelay) + static_cast<size_t>(abs(intshift)))) ? 0u : (N - static_cast<size_t>(SamplesDelay) - static_cast<size_t>(abs(intshift))));
	const size_t carrierlength_hi = std::min(static_cast<size_t>(SamplesWave), (N < (static_cast<size_t>(SamplesDelay))) ? 0u : (N - static_cast<size_t>(SamplesDelay)));
	const size_t uiinitshift = static_cast<size_t>(abs(intshift));

	if (DeltaT >= 0)
	{
		size_t NSamplesDelay = static_cast<size_t>(SamplesDelay);
		typename CDSPArray<T, N>::iterator itemup = up.begin() + NSamplesDelay + 1u;
		typename CDSPArray<T, N>::iterator itemdn = dn.begin() + NSamplesDelay + uiinitshift + 1u;

		for (size_t i = 1; (i < carrierlength_hi); ++i)
		{
			const double fi = static_cast<double>(i);
			const double formula = (0.5f * (1.f - cos(arg2 * fi)) * sin(arg * fi)) * Amplitude + Offset;
			const T Tformula = static_cast<T>(formula);
			up.insert(itemup, Tformula);
			++itemup;
		}

		for (size_t i = 1; (i < carrierlength_lo); ++i)
		{
			const double fi = static_cast<double>(i);
			const double delta_approx = fi - fractshift;
			const double formula = ((0.5f * (1.f - cos(arg2 * (delta_approx))) * sin(arg * (delta_approx)))) * Amplitude + Offset;
			const T Tformula = static_cast<T>(formula);
			dn.insert(itemdn, Tformula);
			++itemdn;
		}
	}
	else
	{
		intshift = -intshift;
		fractshift = -fractshift;
		size_t NSamplesDelay = static_cast<size_t>(SamplesDelay);
		typename CDSPArray<T, N>::iterator itemup = up.begin() + NSamplesDelay + uiinitshift + 1u;
		typename CDSPArray<T, N>::iterator itemdn = dn.begin() + NSamplesDelay + 1u;

		for (size_t i = 1; (i < carrierlength_hi); ++i)
		{
			const double fi = static_cast<double>(i);
			const double formula = (0.5f * (1.f - cos(arg2 * fi)) * sin(arg * fi)) * Amplitude + Offset;
			const T Tformula = static_cast<T>(formula);
			dn.insert(itemdn, Tformula);
			++itemdn;
		}

		for (size_t i = 1; (i < carrierlength_lo); ++i)
		{
			const double fi = static_cast<double>(i);
			const double delta_approx = fi - fractshift;
			const double formula = ((0.5f * (1.f - cos(arg2 * (delta_approx))) * sin(arg * (delta_approx)))) * Amplitude + Offset;
			const T Tformula = static_cast<T>(formula);
			up.insert(itemup, Tformula);
			++itemup;
		}
	}
}

TEST(sample_test_case, float32_cos_sin)
{
	CFFT<10, float, float, FFT_INTERNALS::f32_complex> mf_FFT;

	static CDSPArray<FFT_INTERNALS::f32_complex, 512> spectrum;
	static CDSPArray<float, 1024> ref_signal;
	static CDSPArray<float, 1024> transformed_signal;
	for (size_t i = 0; i < 1024; i++)
	{
		ref_signal.push_back(1.f * cos(512.f * (6.283185307179586476925286f * static_cast<float>(i)) / 1024.f) + 1.f * sin((2 * 6.283185307179586476925286f * static_cast<float>(i)) / 1024.f) - 10e3f);
	}
	spectrum.resize(512);
	transformed_signal.resize(1024);
	mf_FFT.CalculateRealFFT(ref_signal.data(), spectrum.data());
	mf_FFT.CalculateRealIFFT(spectrum.data(), reinterpret_cast<FFT_INTERNALS::f32_complex*>(transformed_signal.data()));
	for (size_t i = 0; i < 1024; ++i)
	{
		ASSERT_FLOAT_EQ(ref_signal[i], transformed_signal[i]);
	}
	
	for (size_t i = 0; i < 512; ++i)
	{
		if (i == 0)
		{
			ASSERT_NEAR(spectrum[i].re, -1024 * 10e3, 1e-1);
			ASSERT_NEAR(spectrum[i].im, 1024, 1e-1);
		}
		else if (i == 2)
		{
			ASSERT_NEAR(spectrum[i].re, 0, 1e-4);
			ASSERT_NEAR(spectrum[i].im, -512, 1e-1);
		}
		else
		{
			ASSERT_NEAR(spectrum[i].re, 0, 5e-2);
			ASSERT_NEAR(spectrum[i].im, 0, 5e-2);
		}
	}
}

TEST(sample_test_case, float64_cos_sin)
{
	CFFT<10, double, double, FFT_INTERNALS::f64_complex> mf_FFT;

	static CDSPArray<FFT_INTERNALS::f64_complex, 512> spectrum;
	static CDSPArray<double, 1024> ref_signal;
	static CDSPArray<double, 1024> transformed_signal;
	for (size_t i = 0; i < 1024; i++)
	{
		ref_signal.push_back(1. * cos(512. * (6.283185307179586476925286 * static_cast<double>(i)) / 1024.) + 1. * sin((2. * 6.283185307179586476925286 * static_cast<double>(i)) / 1024.) - 10.e3);
	}
	spectrum.resize(512);
	transformed_signal.resize(1024);
	mf_FFT.CalculateRealFFT(ref_signal.data(), spectrum.data());
	mf_FFT.CalculateRealIFFT(spectrum.data(), reinterpret_cast<FFT_INTERNALS::f64_complex*>(transformed_signal.data()));
	for (size_t i = 0; i < 1024; ++i)
	{
		ASSERT_DOUBLE_EQ(ref_signal[i], transformed_signal[i]);
	}

	for (size_t i = 0; i < 512; ++i)
	{
		if (i == 0)
		{
			ASSERT_NEAR(spectrum[i].re, -1024 * 10e3, 1e-1);
			ASSERT_NEAR(spectrum[i].im, 1024, 1e-1);
		}
		else if (i == 2)
		{
			ASSERT_NEAR(spectrum[i].re, 0, 1e-4);
			ASSERT_NEAR(spectrum[i].im, -512, 1e-1);
		}
		else
		{
			ASSERT_NEAR(spectrum[i].re, 0, 5e-2);
			ASSERT_NEAR(spectrum[i].im, 0, 5e-2);
		}
	}
}

TEST(sample_test_case, float64_cos_sin_complex)
{
	CFFT<10, double, double, FFT_INTERNALS::f64_complex> mf_FFT;

	static CDSPArray<FFT_INTERNALS::f64_complex, 512> spectrum;
	static CDSPArray<double, 1024> ref_signal;
	static CDSPArray<double, 1024> transformed_signal;
	for (size_t i = 0; i < 1024; i++)
	{
		ref_signal.push_back(65536 * 2048 * cos(512. * (6.283185307179586476925286 * static_cast<double>(i)) / 1024.) + 65536 * 4096. * sin((2 * 6.283185307179586476925286 * static_cast<double>(i)) / 1024.) - 65536 * 10);
	}
	spectrum.resize(512);
	transformed_signal.resize(1024);
	mf_FFT.CalculateRealFFT(ref_signal.data(), spectrum.data());
	static const FFT_INTERNALS::f64_complex normalizer = { 1. / 1024, 0 };
	for (auto& i : spectrum)
	{
		i = i * normalizer;
	}
	mf_FFT.CalculateRealIFFT(spectrum.data(), reinterpret_cast<FFT_INTERNALS::f64_complex*>(transformed_signal.data()));
	for (size_t i = 0; i < 1024; ++i)
	{
		ASSERT_NEAR(ref_signal[i], transformed_signal[i]*1024, 1e-3);
	}

	for (size_t i = 0; i < 512; ++i)
	{
		if (i == 0)
		{
			ASSERT_NEAR(spectrum[i].re, -10*65536, 1e-1);
			ASSERT_NEAR(spectrum[i].im, 65536 * 2048, 1e-1);
		}
		else if (i == 2)
		{
			ASSERT_NEAR(spectrum[i].re, 0, 1e-4);
			ASSERT_NEAR(spectrum[i].im, -65536 * 4096./2, 1e-1);
		}
		else
		{
			ASSERT_NEAR(spectrum[i].re, 0, 5e-2);
			ASSERT_NEAR(spectrum[i].im, 0, 5e-2);
		}
	}
}

TEST(sample_test_case, int16_cos_sin)
{
	CFFT<10, short, short, FFT_INTERNALS::i16_complex> m_FFT(32767);

	static CDSPArray<FFT_INTERNALS::i16_complex, 512> spectrum;
	static CDSPArray<short, 1024> ref_signal;
	static CDSPArray<short, 1024> transformed_signal;
	for (size_t i = 0; i < 1024; i++)
	{
		const double v = 8192. * cos(512. * (6.283185307179586476925286 * static_cast<double>(i)) / 1024.) + 4096. * sin((2 * 6.283185307179586476925286 * static_cast<double>(i)) / 1024.) - 10;
		ref_signal.push_back(static_cast<short>(v));
	}
	spectrum.resize(512);
	transformed_signal.resize(1024);
	m_FFT.CalculateRealFFT(ref_signal.data(), spectrum.data());
	m_FFT.CalculateRealIFFT(spectrum.data(), reinterpret_cast<FFT_INTERNALS::i16_complex*>(transformed_signal.data()));
	for (size_t i = 0; i < 1024; ++i)
	{
		//ASSERT_NEAR(i16Signal[i], i16ftSignal[i], 30);
	}

	for (size_t i = 0; i < 512; ++i)
	{
		if (i == 0)
		{
			ASSERT_NEAR(spectrum[i].re, -10, 10);
			ASSERT_NEAR(spectrum[i].im, 8192, 10);
		}
		else if (i == 2)
		{
			ASSERT_NEAR(spectrum[i].re, 0, 10);
			ASSERT_NEAR(spectrum[i].im, -2048, 10);
		}
		else
		{
			ASSERT_NEAR(spectrum[i].re, 0, 10);
			ASSERT_NEAR(spectrum[i].im, 0, 10);
		}
	}
#if 0	
	std::vector<double> signal_ref;
	std::vector<double> signal_transformed;
	for (size_t i = 0; i < 1024; ++i)
	{
		signal_ref.push_back(i16Signal[i]);
	}
	for (size_t i = 0; i < 512; ++i)
	{
		signal_transformed.push_back(i16Transform[i].re);
		signal_transformed.push_back(i16Transform[i].im);
	}
	plot<double>(signal_ref, signal_transformed);
#endif
}

TEST(sample_test_case, int16_cos_sin_scaling)
{
	CFFT<10, short, short, FFT_INTERNALS::i16_complex> m_FFT(32767);

	static FFT_INTERNALS::i16_complex i16Transform[1024];
	static short i16Signal[1024];
	static short i16ftSignal[1024];
	for (size_t i = 0; i < 1024; i++)
	{
		double v = 8192. * cos(512. * (6.283185307179586476925286 * static_cast<double>(i)) / 1024.) + 4096. * sin((2 * 6.283185307179586476925286 * static_cast<double>(i)) / 1024.) - 10;
		i16Signal[i] = static_cast<short>(v);
	}
	size_t exponent;
	m_FFT.CalculateRealFFT(i16Signal, i16Transform, EN_ScalingMethod::SCALEINPOUT, &exponent);
	size_t exponent_back;
	m_FFT.CalculateRealIFFT(i16Transform, reinterpret_cast<FFT_INTERNALS::i16_complex*>(i16ftSignal), EN_ScalingMethod::SCALEINPOUT, &exponent_back);
#if 0
	for (size_t i = 0; i < 1024; ++i)
	{
		//ASSERT_NEAR(i16Signal[i], i16ftSignal[i], 30);
	}

	for (size_t i = 0; i < 1024; ++i)
	{
		if (i == 0)
		{
			ASSERT_NEAR(i16Transform[i].re, -10, 10);
			ASSERT_NEAR(i16Transform[i].im, 8192, 10);
		}
		else if (i == 2)
		{
			ASSERT_NEAR(i16Transform[i].re, 0, 10);
			ASSERT_NEAR(i16Transform[i].im, -2048, 10);
		}
		else
		{
			ASSERT_NEAR(i16Transform[i].re, 0, 10);
			ASSERT_NEAR(i16Transform[i].im, 0, 10);
		}
	}
#endif
#if 0	
	std::vector<double> signal_ref;
	std::vector<double> signal_transformed;
	for (size_t i = 0; i < 1024; ++i)
	{
		signal_ref.push_back(i16Signal[i]);
	}
	for (size_t i = 0; i < 1024; ++i)
	{
		signal_transformed.push_back(i16ftSignal[i]);
	}
	plot<double>(signal_ref, signal_transformed);
#endif
}

TEST(sample_test_case, int32_cos_sin_real)
{
	CFFT<10, int, int, FFT_INTERNALS::i32_complex> m_FFT(2147483647);

	static CDSPArray<FFT_INTERNALS::i32_complex, 512> spectrum;
	static CDSPArray<int, 1024> ref_signal;
	static CDSPArray<int, 1024> transformed_signal;
	for (size_t i = 0; i < 1024; i++)
	{
		const double v = 65536 * 2048 * cos(512. * (6.283185307179586476925286 * static_cast<double>(i)) / 1024.) + 65536 * 4096. * sin((2 * 6.283185307179586476925286 * static_cast<double>(i)) / 1024.) - 65536 * 10;
		ref_signal.push_back(static_cast<int>(v));
	}
	spectrum.resize(512);
	transformed_signal.resize(1024);
	m_FFT.CalculateRealFFT(ref_signal.data(), spectrum.data());
	m_FFT.CalculateRealIFFT(spectrum.data(), reinterpret_cast<FFT_INTERNALS::i32_complex*>(transformed_signal.data()));
	for (size_t i = 0; i < 1024; ++i)
	{
		ASSERT_NEAR(ref_signal[i], 2* transformed_signal[i], 2048);
	}

	for (size_t i = 0; i < 512; ++i)
	{
		if (i == 0)
		{
			ASSERT_NEAR(spectrum[i].re, -10*65536, 256);
			ASSERT_NEAR(spectrum[i].im, 65536 * 2048, 256);
		}
		else if (i == 2)
		{
			ASSERT_NEAR(spectrum[i].re, 0, 256);
			ASSERT_NEAR(spectrum[i].im, -65536 * 4096./2, 256);
		}
		else
		{
			ASSERT_NEAR(spectrum[i].re, 0, 256);
			ASSERT_NEAR(spectrum[i].im, 0, 256);
		}
	}
#if 1	
	std::vector<double> signal_ref;
	std::vector<double> signal_transformed;
	for (size_t i = 0; i < 30; ++i)
	{
		signal_ref.push_back(ref_signal[i]);
	}
	for (size_t i = 0; i < 30; ++i)
	{
		signal_transformed.push_back(2*transformed_signal[i]);
	}
	plot<double>(signal_ref, signal_transformed);
#endif
}

TEST(sample_test_case, int32_cos_sin_complex_high_freq)
{
	CFFT<10, int, int, FFT_INTERNALS::i32_complex> m_FFT(2147483647);

	static CDSPArray<FFT_INTERNALS::i32_complex, 1024> spectrum;
	static CDSPArray<int, 1024> ref_signal;
	static CDSPArray<FFT_INTERNALS::i32_complex, 1024> transformed_signal;
	for (size_t i = 0; i < 1024; i++)
	{
		const double v = 65536 * 32767 * cos(512. * (6.283185307179586476925286 * static_cast<double>(i)) / 1024.);// +65536 * 4096. * sin((2 * 6.283185307179586476925286 * static_cast<double>(i)) / 1024.) - 65536 * 10;
		int vv = static_cast<int>(v);
		ref_signal.push_back(vv);
	}
	spectrum.resize(1024);
	transformed_signal.resize(1024);
	m_FFT.CalculateFFT(ref_signal.data(), spectrum.data());
	m_FFT.CalculateIFFT(spectrum.data(), reinterpret_cast<FFT_INTERNALS::i32_complex*>(transformed_signal.data()));
	for (size_t i = 0; i < 1024; ++i)
	{
		ASSERT_NEAR(ref_signal[i], transformed_signal[i].re, 2048);
	}

	for (size_t i = 0; i < 1024; ++i)
	{
		if (i == 512)
		{
			ASSERT_NEAR(spectrum[i].re, 65536 * 32767, 512);
			ASSERT_NEAR(spectrum[i].im, 0, 512);
		}
		else
		{
			ASSERT_NEAR(spectrum[i].re, 0, 256);
			ASSERT_NEAR(spectrum[i].im, 0, 256);
		}
	}
#if 0	
	std::vector<double> signal_ref;
	std::vector<double> signal_transformed;
	for (size_t i = 0; i < 100; ++i)
	{
		signal_ref.push_back(ref_signal[i]);
	}
	for (size_t i = 0; i < 100; ++i)
	{
		signal_transformed.push_back(transformed_signal[i].re);
	}
	plot<double>(signal_ref, signal_transformed);
#endif
}

TEST(sample_test_case, int32_cos_sin_complex_mixed_freq)
{
	CFFT<10, int, int, FFT_INTERNALS::i32_complex> m_FFT(2147483647);

	static CDSPArray<FFT_INTERNALS::i32_complex, 1024> spectrum;
	static CDSPArray<int, 1024> ref_signal;
	static CDSPArray<FFT_INTERNALS::i32_complex, 1024> transformed_signal;
	for (size_t i = 0; i < 1024; i++)
	{
		const double v = 65536 * 2048 * cos(512. * (6.283185307179586476925286 * static_cast<double>(i)) / 1024.) + 65536 * 4096. * sin((2 * 6.283185307179586476925286 * static_cast<double>(i)) / 1024.) - 65536 * 10;
		int vv = static_cast<int>(v);
		ref_signal.push_back(vv);
	}
	spectrum.resize(1024);
	transformed_signal.resize(1024);
	m_FFT.CalculateFFT(ref_signal.data(), spectrum.data());
	m_FFT.CalculateIFFT(spectrum.data(), reinterpret_cast<FFT_INTERNALS::i32_complex*>(transformed_signal.data()));
	for (size_t i = 0; i < 1024; ++i)
	{
		ASSERT_NEAR(ref_signal[i], transformed_signal[i].re, 2048);
	}

	for (size_t i = 0; i < 1024; ++i)
	{
		if (i == 512)
		{
			ASSERT_NEAR(spectrum[i].re, 65536 * 2048, 512);
			ASSERT_NEAR(spectrum[i].im, 0, 512);
		}
		else if (i == 2)
		{
			ASSERT_NEAR(spectrum[i].re, 0, 512);
			ASSERT_NEAR(spectrum[i].im, -65536 * 4096/2, 512);
		}
		else if (i == 1024-2)
		{
			ASSERT_NEAR(spectrum[i].re, 0, 512);
			ASSERT_NEAR(spectrum[i].im, 65536 * 4096 / 2, 512);
		}
		else if (i == 0)
		{
			ASSERT_NEAR(spectrum[i].re, -65536*10, 512);
			ASSERT_NEAR(spectrum[i].im, 0, 512);
		}
		else
		{
			ASSERT_NEAR(spectrum[i].re, 0, 256);
			ASSERT_NEAR(spectrum[i].im, 0, 256);
		}
	}
#if 0	
	std::vector<double> signal_ref;
	std::vector<double> signal_transformed;
	for (size_t i = 0; i < 100; ++i)
	{
		signal_ref.push_back(ref_signal[i]);
	}
	for (size_t i = 0; i < 100; ++i)
	{
		signal_transformed.push_back(transformed_signal[i].re);
	}
	plot<double>(signal_ref, signal_transformed);
#endif
}

TEST(sample_test_case, float32_Hanning)
{
	CFFT<10, float, float, FFT_INTERNALS::f32_complex> mf_FFT;
	CDSPArray<float, 1024> up;
	CDSPArray<float, 1024> dn;
	static const double Fs = 16e6;
	static const double Fx = 1e6;
	static const double ToF = 10e-6;
	static const double DeltaT = 10e-9;
	static const double SignalDuration = 15e-6;
	static const double SignalAmplitude = 1.;
	static const double SignalOffset = 0.;

	GenerateSignalPairHannSym<float, 1024>(up, dn,
		Fs, Fx,
		ToF, DeltaT,
		0., SignalDuration,
		SignalAmplitude, SignalOffset);
	
	CDSPArray<FFT_INTERNALS::f32_complex, 1024> reference_complex_signal;
	CDSPArray<FFT_INTERNALS::f32_complex, 1024> reference_complex_spectrum;
	
	for (const auto& i : up)
	{
		FFT_INTERNALS::f32_complex t;
		t.re = i;
		t.im = 0;
		reference_complex_signal.push_back(t);
	}
	CDSPArray<FFT_INTERNALS::f32_complex, 1024> calculated_complex_spectrum;
	CDSPArray<FFT_INTERNALS::f32_complex, 1024> calculated_complex_signal;
	calculated_complex_spectrum.resize(1024);
	CalculateDFT(1024, reference_complex_signal.data(), reference_complex_spectrum.data());
	mf_FFT.CalculateFFT(up.data(), calculated_complex_spectrum.data());
	mf_FFT.CalculateIFFT(calculated_complex_spectrum.data(), calculated_complex_signal.data());
	for (size_t i = 0; i < 1024; ++i)
	{
		ASSERT_NEAR(reference_complex_spectrum[i].re, calculated_complex_spectrum[i].re, 1e-4);
		ASSERT_NEAR(reference_complex_spectrum[i].im, calculated_complex_spectrum[i].im, 1e-4);
	}
	for (size_t i = 0; i < 1024; ++i)
	{
		ASSERT_NEAR(reference_complex_signal[i].re, calculated_complex_signal[i].re, 1e-4);
		ASSERT_NEAR(reference_complex_signal[i].im, calculated_complex_signal[i].im, 1e-4);
	}

	mf_FFT.CalculateRealFFT(up.data(), calculated_complex_spectrum.data());
	mf_FFT.CalculateRealIFFT(calculated_complex_spectrum.data(), calculated_complex_signal.data());
	for (size_t i = 1; i < 512; ++i)
	{
		ASSERT_NEAR(reference_complex_spectrum[i].re, calculated_complex_spectrum[i].re, 1e-4);
		ASSERT_NEAR(reference_complex_spectrum[i].im, calculated_complex_spectrum[i].im, 1e-4);
	}
	for (size_t i = 0; i < 512; ++i)
	{
		ASSERT_NEAR(reference_complex_signal[2*i].re, calculated_complex_signal[i].re, 1e-4);
		ASSERT_NEAR(reference_complex_signal[2*i+1].re, calculated_complex_signal[i].im, 1e-4);
	}
}

TEST(sample_test_case, float64_Hanning)
{
	CFFT<10, double, double, FFT_INTERNALS::f64_complex> mf_FFT;
	CDSPArray<double, 1024> up;
	CDSPArray<double, 1024> dn;
	static const double Fs = 16e6;
	static const double Fx = 1e6;
	static const double ToF = 10e-6;
	static const double DeltaT = 10e-9;
	static const double SignalDuration = 15e-6;
	static const double SignalAmplitude = 1.;
	static const double SignalOffset = 0.;

	GenerateSignalPairHannSym<double, 1024>(up, dn,
		Fs, Fx,
		ToF, DeltaT,
		0., SignalDuration,
		SignalAmplitude, SignalOffset);

	CDSPArray<FFT_INTERNALS::f64_complex, 1024> reference_complex_signal;
	CDSPArray<FFT_INTERNALS::f64_complex, 1024> reference_complex_spectrum;

	for (const auto& i : up)
	{
		FFT_INTERNALS::f64_complex t;
		t.re = i;
		t.im = 0;
		reference_complex_signal.push_back(t);
	}
	CDSPArray<FFT_INTERNALS::f64_complex, 1024> calculated_complex_spectrum;
	CDSPArray<FFT_INTERNALS::f64_complex, 1024> calculated_complex_signal;
	calculated_complex_spectrum.resize(1024);
	CalculateDFT(1024, reference_complex_signal.data(), reference_complex_spectrum.data());
	mf_FFT.CalculateFFT(up.data(), calculated_complex_spectrum.data());
	mf_FFT.CalculateIFFT(calculated_complex_spectrum.data(), calculated_complex_signal.data());
	for (size_t i = 0; i < 1024; ++i)
	{
		ASSERT_NEAR(reference_complex_spectrum[i].re, calculated_complex_spectrum[i].re, 1e-4);
		ASSERT_NEAR(reference_complex_spectrum[i].im, calculated_complex_spectrum[i].im, 1e-4);
	}
	for (size_t i = 0; i < 1024; ++i)
	{
		ASSERT_NEAR(reference_complex_signal[i].re, calculated_complex_signal[i].re, 1e-4);
		ASSERT_NEAR(reference_complex_signal[i].im, calculated_complex_signal[i].im, 1e-4);
	}

	mf_FFT.CalculateRealFFT(up.data(), calculated_complex_spectrum.data());
	mf_FFT.CalculateRealIFFT(calculated_complex_spectrum.data(), calculated_complex_signal.data());
	for (size_t i = 1; i < 512; ++i)
	{
		ASSERT_NEAR(reference_complex_spectrum[i].re, calculated_complex_spectrum[i].re, 1e-4);
		ASSERT_NEAR(reference_complex_spectrum[i].im, calculated_complex_spectrum[i].im, 1e-4);
	}
	for (size_t i = 0; i < 512; ++i)
	{
		ASSERT_NEAR(reference_complex_signal[2 * i].re, calculated_complex_signal[i].re, 1e-4);
		ASSERT_NEAR(reference_complex_signal[2 * i + 1].re, calculated_complex_signal[i].im, 1e-4);
	}
}

int main(int argc, char** argv)
{
	//testing::InitGoogleTest(&argc, argv);
	testing::InitGoogleMock(&argc, argv);
	RUN_ALL_TESTS();
	//std::getchar();
}
#else
int _tmain(int argc, _TCHAR* argv[])
{
	CFFT<10,int, int, FFT_INTERNALS::i32_complex> m_FFT(2147483647);
	CFFT<10,int, int, FFT_INTERNALS::i32_complex> m_FFT2(2147483647);
	CFFT<10,float, float, FFT_INTERNALS::f32_complex> mf_FFT;
	CFFT<10,double, double, FFT_INTERNALS::f64_complex> mdf_FFT;
	CFFT<10,double, double, FFT_INTERNALS::f64_complex> mrdf_FFT;

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
	mdf_FFT.CalculateRealIFFT(dfTransform,reinterpret_cast<FFT_INTERNALS::f64_complex*>(resultd));

	CFFT<10, double, double, FFT_INTERNALS::f64_complex> dummy;
	mf_FFT.CalculateRealFFT(fSignal,fTransform);
	mf_FFT.CalculateRealIFFT(fTransform, reinterpret_cast<FFT_INTERNALS::f32_complex*>(resultf));
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
	m_FFT.CalculateRealFFT(Signal,TransformScaling, EN_ScalingMethod::SCALEINPOUT, &Exponent1);
	m_FFT.CalculateRealIFFT(TransformScaling,Transform, EN_ScalingMethod::SCALEINPOUT, &Exponent2);

	size_t exponent;
	m_FFT2.CalculateRealFFT(Signal,RealTransform, EN_ScalingMethod::SCALEINPOUT, &exponent);

	FFT_INTERNALS::f64_complex test_signal[1024];
	FFT_INTERNALS::f64_complex test_result[1024];
	size_t index = 0;
	for (FFT_INTERNALS::f64_complex& a : test_signal)
	{
		a.re = fSignal[index++];
		a.im = 0;
	}
	CalculateDFT(1024, test_signal, test_result);

	//dummy.Convert2HalfDFT(Reference);
	return 0;
}
#endif
