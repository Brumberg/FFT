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

template <size_t loglen, typename T, size_t len>
void DoScaling(size_t fft_scale, size_t ifft_scale, CDSPArray<T, len>& src)
{
	static const size_t base_scaling = loglen;
	const size_t req_scaling = fft_scale + ifft_scale;
	if (req_scaling > base_scaling)
	{
		const size_t scaling = req_scaling - base_scaling;
		for (auto& i : src)
		{
			i <<= scaling;
		}
	}
	else
	{
		const size_t scaling = base_scaling - req_scaling;
		for (auto& i : src)
		{
			i >>= scaling;
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
		ASSERT_NEAR(ref_signal[i], transformed_signal[i], 5e-5);
	}

	for (size_t i = 0; i < 512; ++i)
	{
		if (i == 0)
		{
			ASSERT_NEAR(spectrum[i].re, -1024 * 10e3, 1e-10);
			ASSERT_NEAR(spectrum[i].im, 1024, 1e-10);
		}
		else if (i == 2)
		{
			ASSERT_NEAR(spectrum[i].re, 0, 1e-10);
			ASSERT_NEAR(spectrum[i].im, -512, 1e-10);
		}
		else
		{
			ASSERT_NEAR(spectrum[i].re, 0, 1e-10);
			ASSERT_NEAR(spectrum[i].im, 0, 1e-10);
		}
	}

#if 0	
	std::vector<double> signal_ref;
	std::vector<double> signal_transformed;
	std::vector<double> signal_diff;
	for (size_t i = 0; i < 1024; ++i)
	{
		signal_ref.push_back(ref_signal[i]);
	}
	for (size_t i = 0; i < 1024; ++i)
	{
		signal_transformed.push_back(transformed_signal[i]);
	}
	plot<double>(signal_ref, signal_transformed);
	for (size_t i = 0; i < 1024; ++i)
	{
		signal_diff.push_back(signal_ref[i] - transformed_signal[i]);
	}
	plot<double>(signal_diff);
#endif
}

TEST(sample_test_case, float64_cos_sin_complex)
{
	CFFT<10, double, double, FFT_INTERNALS::f64_complex> mf_FFT;

	static CDSPArray<FFT_INTERNALS::f64_complex, 512> spectrum;
	static CDSPArray<double, 1024> ref_signal;
	static CDSPArray<double, 1024> transformed_signal;
	for (size_t i = 0; i < 1024; i++)
	{
		ref_signal.push_back(65536. * 2048. * cos(512. * (6.283185307179586476925286 * static_cast<double>(i)) / 1024.) + 65536. * 4096. * sin((2 * 6.283185307179586476925286 * static_cast<double>(i)) / 1024.) - 65536. * 10.);
	}
	spectrum.resize(512);
	transformed_signal.resize(1024);
	mf_FFT.CalculateRealFFT(ref_signal.data(), spectrum.data());
	static const FFT_INTERNALS::f64_complex normalizer = { 1. / 1024., 0. };
	for (auto& i : spectrum)
	{
		i = i * normalizer;
	}
	mf_FFT.CalculateRealIFFT(spectrum.data(), reinterpret_cast<FFT_INTERNALS::f64_complex*>(transformed_signal.data()));
	for (size_t i = 0; i < 1024; ++i)
	{
		ASSERT_NEAR(ref_signal[i], transformed_signal[i]*1024, 64);
	}

	for (size_t i = 0; i < 512; ++i)
	{
		if (i == 0)
		{
			ASSERT_NEAR(spectrum[i].re, -10.*65536., 50.);
			ASSERT_NEAR(spectrum[i].im, 65536. * 2048., 50.);
		}
		else if (i == 2)
		{
			ASSERT_NEAR(spectrum[i].re, 0., 5.);
			ASSERT_NEAR(spectrum[i].im, -65536. * 4096./2, 30.);
		}
		else
		{
			ASSERT_NEAR(spectrum[i].re, 0., 50.);
			ASSERT_NEAR(spectrum[i].im, 0., 50.);
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
	for (size_t i = 0; i < 150; ++i)
	{
		ASSERT_NEAR(ref_signal[i], 2*transformed_signal[i], 1300);
	}
	for (size_t i = 150; i < 1024-150; ++i)
	{
		ASSERT_NEAR(ref_signal[i], 2 * transformed_signal[i], 100);
	}
	for (size_t i = 1024 - 150; i < 1024; ++i)
	{
		ASSERT_NEAR(ref_signal[i], 2 * transformed_signal[i], 1300);
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
	std::vector<double> signal_diff;
	for (size_t i = 0; i < 1024; ++i)
	{
		signal_ref.push_back(ref_signal[i]);
	}
	for (size_t i = 0; i < 1024; ++i)
	{
		signal_transformed.push_back(2*transformed_signal[i]);
	}
	plot<double>(signal_ref, signal_transformed);
	for (size_t i = 0; i < 1024; ++i)
	{
		signal_diff.push_back(signal_ref[i] - 2 * transformed_signal[i]);
	}
	plot<double>(signal_diff);
#endif
}

TEST(sample_test_case, int16_cos_sin_scaling)
{
	CFFT<10, short, short, FFT_INTERNALS::i16_complex> m_FFT(32767);

	static CDSPArray <FFT_INTERNALS::i16_complex, 1024> i16Transform;
	static CDSPArray<short, 1024> ref_signal;
	static CDSPArray<short, 1024> transformed_signal;
	for (size_t i = 0; i < 1024; i++)
	{
		double v = 8192. * cos(512. * (6.283185307179586476925286 * static_cast<double>(i)) / 1024.) + 4096. * sin((2 * 6.283185307179586476925286 * static_cast<double>(i)) / 1024.) - 10;
		ref_signal.push_back(static_cast<short>(v));
	}
	size_t fft_scaling;
	i16Transform.resize(512);
	m_FFT.CalculateRealFFT(ref_signal.data(), i16Transform.data(), EN_ScalingMethod::SCALEINPOUT, &fft_scaling);
	size_t ifft_scaling;
	m_FFT.CalculateRealIFFT(reinterpret_cast<FFT_INTERNALS::i16_complex*>(i16Transform.data()), reinterpret_cast<FFT_INTERNALS::i16_complex*>(transformed_signal.data()), EN_ScalingMethod::SCALEINPOUT, &ifft_scaling);
	DoScaling<10>(fft_scaling, ifft_scaling, transformed_signal);
	const size_t exponent = fft_scaling;
	for (size_t i = 0; i < 1024; ++i)
	{
		ASSERT_NEAR(ref_signal[i], transformed_signal[i], 750);
	}

	for (size_t i = 0; i < 1024; ++i)
	{
		if (i == 0)
		{
			ASSERT_NEAR(i16Transform[i].re, -10<< (10 - exponent), 5);
			ASSERT_NEAR(i16Transform[i].im, 8192<< (10 - exponent), 5);
		}
		else if (i == 2)
		{
			ASSERT_NEAR(i16Transform[i].re, 0, 5);
			ASSERT_NEAR(i16Transform[i].im, -4096 << (10 - exponent - 1), 5);
		}
		else
		{
			ASSERT_NEAR(i16Transform[i].re, 0, 5);
			ASSERT_NEAR(i16Transform[i].im, 0, 5);
		}
	}

#if 0	
	std::vector<double> signal_ref;
	std::vector<double> signal_transformed;
	std::vector<double> signal_diff;

	for (size_t i = 0; i < 1024; ++i)
	{
		signal_ref.push_back(ref_signal[i]);
	}
	for (size_t i = 0; i < 1024; ++i)
	{
		signal_transformed.push_back(transformed_signal[i]);
	}
	plot<double>(signal_ref, signal_transformed);

	for (size_t i = 0; i < 1024; ++i)
	{
		signal_diff.push_back(ref_signal[i] - transformed_signal[i]);
	}
	plot<double>(signal_diff);
#endif
}

TEST(sample_test_case, int16_cos_sin_low_scaling)
{
	CFFT<10, short, short, FFT_INTERNALS::i16_complex> m_FFT(32767);

	static CDSPArray<FFT_INTERNALS::i16_complex, 512> spectrum;
	static CDSPArray<short, 1024> ref_signal;
	static CDSPArray<short, 1024> transformed_signal;

	for (size_t i = 0; i < 1024; i++)
	{
		double v = 120 * cos(512. * (6.283185307179586476925286 * static_cast<double>(i)) / 1024.) + 196. * sin((2 * 6.283185307179586476925286 * static_cast<double>(i)) / 1024.);
		ref_signal.push_back(static_cast<short>(v));
	}
	transformed_signal.resize(1024);
	size_t fft_scaling;
	m_FFT.CalculateRealFFT(ref_signal.data(), spectrum.data(), EN_ScalingMethod::SCALEINPOUT, &fft_scaling);
	size_t ifft_scaling;
	m_FFT.CalculateRealIFFT(spectrum.data(), reinterpret_cast<FFT_INTERNALS::i16_complex*>(transformed_signal.data()), EN_ScalingMethod::SCALEINPOUT, &ifft_scaling);
	DoScaling<10>(fft_scaling, ifft_scaling, transformed_signal);

	for (size_t i = 0; i < 1024; ++i)
	{
		ASSERT_NEAR(ref_signal[i], transformed_signal[i], 512);
	}

	for (size_t i = 0; i < 512; ++i)
	{
		if (i == 2)
		{
			ASSERT_NEAR(spectrum[i].re, 0 << (10 - fft_scaling - 1), 15);
			ASSERT_NEAR(spectrum[i].im, -196 << (10 - fft_scaling - 1), 80);
		}
		else if (i == 0)
		{
			ASSERT_NEAR(spectrum[i].re, 0 << (10 - fft_scaling), 15);
			ASSERT_NEAR(spectrum[i].im, 120 << (10 - fft_scaling), 58);
		}
		else
		{
			ASSERT_NEAR(spectrum[i].re, 0, 40);
			ASSERT_NEAR(spectrum[i].im, 0, 40);
		}
	}

#if 0	
	std::vector<double> signal_ref;
	std::vector<double> signal_transformed;
	std::vector<double> signal_diff;

	for (size_t i = 0; i < 1024; ++i)
	{
		signal_ref.push_back(ref_signal[i]);
	}
	for (size_t i = 0; i < 1024; ++i)
	{
		signal_transformed.push_back(transformed_signal[i]);
		}
	plot<double>(signal_ref, signal_transformed);

	for (size_t i = 0; i < 1024; ++i)
	{
		signal_diff.push_back(signal_ref[i] - signal_transformed[i]);
	}
	plot<double>(signal_diff);
#endif
}

TEST(sample_test_case, int16_cos_sin_mixed_freq)
{
	CFFT<10, short, short, FFT_INTERNALS::i16_complex> m_FFT(32767);

	static CDSPArray<FFT_INTERNALS::i16_complex, 512> spectrum;
	static CDSPArray<short, 1024> ref_signal;
	static CDSPArray<short, 1024> transformed_signal;

	for (size_t i = 0; i < 1024; i++)
	{
		double v =
			+196. * sin((70 * 6.283185307179586476925286 * static_cast<double>(i)) / 1024.)
			+ 517. * cos((257 * 6.283185307179586476925286 * static_cast<double>(i)) / 1024.)
			- 517. * sin((257 * 6.283185307179586476925286 * static_cast<double>(i)) / 1024.)
			- 2000. * cos((384 * 6.283185307179586476925286 * static_cast<double>(i)) / 1024.)
			+ 1000. * sin((384 * 6.283185307179586476925286 * static_cast<double>(i)) / 1024.);

		ref_signal.push_back(static_cast<short>(v));
	}
	transformed_signal.resize(1024);
	size_t fft_scaling;
	m_FFT.CalculateRealFFT(ref_signal.data(), spectrum.data(), EN_ScalingMethod::SCALEINPOUT, &fft_scaling);
	size_t ifft_scaling;
	m_FFT.CalculateRealIFFT(spectrum.data(), reinterpret_cast<FFT_INTERNALS::i16_complex*>(transformed_signal.data()), EN_ScalingMethod::SCALEINPOUT, &ifft_scaling);
	DoScaling<10>(fft_scaling, ifft_scaling, transformed_signal);

	for (size_t i = 0; i < 1024; ++i)
	{
		ASSERT_NEAR(ref_signal[i], transformed_signal[i], 512);
	}

	for (size_t i = 0; i < 512; ++i)
	{
		if (i == 70)
		{
			ASSERT_NEAR(spectrum[i].re, 0 << (10 - fft_scaling - 1), 15);
			ASSERT_NEAR(spectrum[i].im, -196 << (10 - fft_scaling - 1), 15);
		}
		else if (i == 257)
		{
			ASSERT_NEAR(spectrum[i].re, 517 << (10 - fft_scaling - 1), 15);
			ASSERT_NEAR(spectrum[i].im, 517 << (10 - fft_scaling - 1), 15);
		}
		else if (i == 384)
		{
			ASSERT_NEAR(spectrum[i].re, -2000 << (10 - fft_scaling - 1), 15);
			ASSERT_NEAR(spectrum[i].im, -1000 << (10 - fft_scaling - 1), 15);
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
	std::vector<double> signal_diff;

	for (size_t i = 0; i < 1024; ++i)
	{
		signal_ref.push_back(ref_signal[i]);
	}
	for (size_t i = 0; i < 1024; ++i)
	{
		signal_transformed.push_back(transformed_signal[i]);
	}
	plot<double>(signal_ref, signal_transformed);

	for (size_t i = 0; i < 1024; ++i)
	{
		signal_diff.push_back(signal_ref[i] - signal_transformed[i]);
	}
	plot<double>(signal_diff);
#endif
}

TEST(sample_test_case, int16_cos_sin_mixed_freq_midrange_check_single_bin)
{
	CFFT<10, short, short, FFT_INTERNALS::i16_complex> m_FFT(32767);

	static CDSPArray<FFT_INTERNALS::i16_complex, 512> spectrum;
	static CDSPArray<short, 1024> ref_signal;
	static CDSPArray<short, 1024> transformed_signal;

	for (size_t i = 0; i < 1024; i++)
	{
		double v =
			+517. * cos((256 * 6.283185307179586476925286 * static_cast<double>(i)) / 1024.)
			-517. * sin((256 * 6.283185307179586476925286 * static_cast<double>(i)) / 1024.);

		ref_signal.push_back(static_cast<short>(v));
	}
	transformed_signal.resize(1024);
	size_t fft_scaling;
	m_FFT.CalculateRealFFT(ref_signal.data(), spectrum.data(), EN_ScalingMethod::SCALEINPOUT, &fft_scaling);
	size_t ifft_scaling;
	m_FFT.CalculateRealIFFT(spectrum.data(), reinterpret_cast<FFT_INTERNALS::i16_complex*>(transformed_signal.data()), EN_ScalingMethod::SCALEINPOUT, &ifft_scaling);
	DoScaling<10>(fft_scaling, ifft_scaling, transformed_signal);

	for (size_t i = 0; i < 1024; ++i)
	{
		ASSERT_NEAR(ref_signal[i], transformed_signal[i], 512);
	}

	for (size_t i = 0; i < 512; ++i)
	{
		if (i == 256)
		{
			ASSERT_NEAR(spectrum[i].re, 517 << (10 - fft_scaling - 1), 15);
			ASSERT_NEAR(spectrum[i].im, 517 << (10 - fft_scaling - 1), 15);
		}
		else
		{
			ASSERT_NEAR(spectrum[i].re, 0, 5);
			ASSERT_NEAR(spectrum[i].im, 0, 5);
		}
	}

#if 0	
	std::vector<double> signal_ref;
	std::vector<double> signal_transformed;
	std::vector<double> signal_diff;

	for (size_t i = 0; i < 1024; ++i)
	{
		signal_ref.push_back(ref_signal[i]);
	}
	for (size_t i = 0; i < 1024; ++i)
	{
		signal_transformed.push_back(transformed_signal[i]);
	}
	plot<double>(signal_ref, signal_transformed);

	for (size_t i = 0; i < 1024; ++i)
	{
		signal_diff.push_back(signal_ref[i] - signal_transformed[i]);
	}
	plot<double>(signal_diff);
#endif
}

TEST(sample_test_case, int16_cos_sin_mixed_freq_midrange_check)
{
	CFFT<10, short, short, FFT_INTERNALS::i16_complex> m_FFT(32767);

	static CDSPArray<FFT_INTERNALS::i16_complex, 512> spectrum;
	static CDSPArray<short, 1024> ref_signal;
	static CDSPArray<short, 1024> transformed_signal;

	for (size_t i = 0; i < 1024; i++)
	{
		double v =
			+196. * sin((70 * 6.283185307179586476925286 * static_cast<double>(i)) / 1024.)
			+ 517. * cos((256 * 6.283185307179586476925286 * static_cast<double>(i)) / 1024.)
			- 517. * sin((256 * 6.283185307179586476925286 * static_cast<double>(i)) / 1024.)
			- 2000. * cos((384 * 6.283185307179586476925286 * static_cast<double>(i)) / 1024.)
			+ 1000. * sin((384 * 6.283185307179586476925286 * static_cast<double>(i)) / 1024.);
		ref_signal.push_back(static_cast<short>(v));
	}
	transformed_signal.resize(1024);
	size_t fft_scaling;
	m_FFT.CalculateRealFFT(ref_signal.data(), spectrum.data(), EN_ScalingMethod::SCALEINPOUT, &fft_scaling);
	size_t ifft_scaling;
	m_FFT.CalculateRealIFFT(spectrum.data(), reinterpret_cast<FFT_INTERNALS::i16_complex*>(transformed_signal.data()), EN_ScalingMethod::SCALEINPOUT, &ifft_scaling);
	DoScaling<10>(fft_scaling, ifft_scaling, transformed_signal);

	for (size_t i = 0; i < 1024; ++i)
	{
		ASSERT_NEAR(ref_signal[i], transformed_signal[i], 512);
	}

	for (size_t i = 0; i < 512; ++i)
	{
		if (i == 70)
		{
			ASSERT_NEAR(spectrum[i].re, 0 << (10 - fft_scaling - 1), 10);
			ASSERT_NEAR(spectrum[i].im, -196 << (10 - fft_scaling - 1), 64);
		}
		else if (i == 256)
		{
			ASSERT_NEAR(spectrum[i].re, 517 << (10 - fft_scaling - 1), 10);
			ASSERT_NEAR(spectrum[i].im, 517 << (10 - fft_scaling - 1), 10);
		}
		else if (i == 384)
		{
			ASSERT_NEAR(spectrum[i].re, -2000 << (10 - fft_scaling - 1), 10);
			ASSERT_NEAR(spectrum[i].im, -1000 << (10 - fft_scaling - 1), 10);
		}
		else
		{
			ASSERT_NEAR(spectrum[i].re, 0, 40);
			ASSERT_NEAR(spectrum[i].im, 0, 40);
		}
	}

#if 0	
	std::vector<double> signal_ref;
	std::vector<double> signal_transformed;
	std::vector<double> signal_diff;

	for (size_t i = 0; i < 1024; ++i)
	{
		signal_ref.push_back(ref_signal[i]);
	}
	for (size_t i = 0; i < 1024; ++i)
	{
		signal_transformed.push_back(transformed_signal[i]);
	}
	plot<double>(signal_ref, signal_transformed);

	for (size_t i = 0; i < 1024; ++i)
	{
		signal_diff.push_back(signal_ref[i] - signal_transformed[i]);
	}
	plot<double>(signal_diff);
#endif
}

TEST(sample_test_case, int16_cos_sin_impulse_response)
{
	CFFT<10, short, short, FFT_INTERNALS::i16_complex> m_FFT(32767);

	static CDSPArray<FFT_INTERNALS::i16_complex, 512> spectrum;
	static CDSPArray<short, 1024> ref_signal;
	static CDSPArray<short, 1024> transformed_signal;

	for (size_t i = 0; i < 1024; i++)
	{
		double v = +0;
		ref_signal.push_back(static_cast<short>(v));
	}
	ref_signal[0] = -32768;
	transformed_signal.resize(1024);
	size_t fft_scaling;
	m_FFT.CalculateRealFFT(ref_signal.data(), spectrum.data(), EN_ScalingMethod::SCALEINPOUT, &fft_scaling);
	size_t ifft_scaling;
	m_FFT.CalculateRealIFFT(spectrum.data(), reinterpret_cast<FFT_INTERNALS::i16_complex*>(transformed_signal.data()), EN_ScalingMethod::SCALEINPOUT, &ifft_scaling);
	DoScaling<10>(fft_scaling, ifft_scaling, transformed_signal);

	for (size_t i = 0; i < 1024; ++i)
	{
		ASSERT_NEAR(ref_signal[i], transformed_signal[i], 34);
	}

	for (size_t i = 0; i < 512; ++i)
	{
			ASSERT_NEAR(spectrum[i].re, -32768 >> (fft_scaling), 40);//because fft scaling is zero and the shift is in lower direction
			if (i == 0)
			{
				ASSERT_NEAR(spectrum[i].re, -32768 >> (fft_scaling), 40);//because fft scaling is zero and the shift is in lower direction
			}
			else
			{
				ASSERT_NEAR(spectrum[i].im, 0, 40);
			}
	}

#if 0	
	std::vector<double> signal_ref;
	std::vector<double> signal_transformed;
	std::vector<double> signal_diff;

	for (size_t i = 0; i < 1024; ++i)
	{
		signal_ref.push_back(ref_signal[i]);
	}
	for (size_t i = 0; i < 1024; ++i)
	{
		signal_transformed.push_back(transformed_signal[i]);
	}
	plot<double>(signal_ref, signal_transformed);

	for (size_t i = 0; i < 1024; ++i)
	{
		signal_diff.push_back(signal_ref[i] - signal_transformed[i]);
	}
	plot<double>(signal_diff);
#endif
}

TEST(sample_test_case, int32_cos_sin_impulse_response)
{
	CFFT<10, int, int, FFT_INTERNALS::i32_complex> m_FFT(32767*65536);

	static CDSPArray<FFT_INTERNALS::i32_complex, 512> spectrum;
	static CDSPArray<int, 1024> ref_signal;
	static CDSPArray<int, 1024> transformed_signal;

	for (size_t i = 0; i < 1024; i++)
	{
		double v = +0;
		ref_signal.push_back(static_cast<int>(v));
	}
	ref_signal[0] = -32767*32767;
	transformed_signal.resize(1024);
	size_t fft_scaling;
	m_FFT.CalculateRealFFT(ref_signal.data(), spectrum.data(), EN_ScalingMethod::SCALEINPOUT, &fft_scaling);
	size_t ifft_scaling;
	//spectrum[0].im *= -1;
	m_FFT.CalculateRealIFFT(spectrum.data(), reinterpret_cast<FFT_INTERNALS::i32_complex*>(transformed_signal.data()), EN_ScalingMethod::SCALEINPOUT, &ifft_scaling);
	DoScaling<10>(fft_scaling, ifft_scaling, transformed_signal);
#if 0
	for (size_t i = 0; i < 1024; ++i)
	{
		ASSERT_NEAR(ref_signal[i], transformed_signal[i], 34);
	}

	for (size_t i = 0; i < 512; ++i)
	{
		ASSERT_NEAR(spectrum[i].re, -32768*65536 >> (fft_scaling), 40);//because fft scaling is zero and the shift is in lower direction
		if (i == 0)
		{
			ASSERT_NEAR(spectrum[i].re, -32768*65536 >> (fft_scaling), 40);//because fft scaling is zero and the shift is in lower direction
		}
		else
		{
			ASSERT_NEAR(spectrum[i].im, 0, 40);
		}
	}
#endif
#if 0	
	std::vector<double> signal_ref;
	std::vector<double> signal_transformed;
	std::vector<double> signal_diff;

	for (size_t i = 0; i < 1024; ++i)
	{
		signal_ref.push_back(ref_signal[i]);
	}
	for (size_t i = 0; i < 1024; ++i)
	{
		signal_transformed.push_back(transformed_signal[i]);
	}
	plot<double>(signal_ref, signal_transformed);

	for (size_t i = 0; i < 1024; ++i)
	{
		signal_diff.push_back(signal_ref[i] - signal_transformed[i]);
	}
	plot<double>(signal_diff);
#endif
}

TEST(sample_test_case, float32_cos_sin_impulse_response)
{
	CFFT<10, float, float, FFT_INTERNALS::f32_complex> m_FFT;

	static CDSPArray<FFT_INTERNALS::f32_complex, 512> spectrum;
	static CDSPArray<float, 1024> ref_signal;
	static CDSPArray<float, 1024> transformed_signal;

	for (size_t i = 0; i < 1024; i++)
	{
		double v = +0;
		ref_signal.push_back(static_cast<float>(v));
	}
	ref_signal[0] = -32768 * 65536;
	transformed_signal.resize(1024);
	m_FFT.CalculateRealFFT(ref_signal.data(), spectrum.data());
	m_FFT.CalculateRealIFFT(spectrum.data(), reinterpret_cast<FFT_INTERNALS::f32_complex*>(transformed_signal.data()));
	//DoScaling<10>(fft_scaling, ifft_scaling, transformed_signal);

	for (size_t i = 0; i < 1024; ++i)
	{
		ASSERT_NEAR(ref_signal[i], transformed_signal[i], 500);
	}

	for (size_t i = 0; i < 512; ++i)
	{
		ASSERT_NEAR(spectrum[i].re, -32768 * 65536, 1000);//because fft scaling is zero and the shift is in lower direction
		if (i == 0)
		{
			ASSERT_NEAR(spectrum[i].re, -32768 * 65536, 1000);//because fft scaling is zero and the shift is in lower direction
		}
		else
		{
			ASSERT_NEAR(spectrum[i].im, 0, 1000);
		}
	}

#if 0	
	std::vector<double> signal_ref;
	std::vector<double> signal_transformed;
	std::vector<double> signal_diff;

	for (size_t i = 0; i < 1024; ++i)
	{
		signal_ref.push_back(ref_signal[i]);
	}
	for (size_t i = 0; i < 1024; ++i)
	{
		signal_transformed.push_back(transformed_signal[i]);
	}
	plot<double>(signal_ref, signal_transformed);

	for (size_t i = 0; i < 1024; ++i)
	{
		signal_diff.push_back(signal_ref[i] - signal_transformed[i]);
	}
	plot<double>(signal_diff);
#endif
}

TEST(sample_test_case, int32_DC_test)
{
	CFFT<10, int, int, FFT_INTERNALS::i32_complex> m_FFT(2147483647);

	static CDSPArray<FFT_INTERNALS::i32_complex, 512> spectrum;
	static CDSPArray<int, 1024> ref_signal;
	static CDSPArray<int, 1024> transformed_signal;
	for (size_t i = 0; i < 1024; i++)
	{
		const double v = 65536 * 16384;
		ref_signal.push_back(static_cast<int>(v));
	}
	spectrum.resize(512);
	transformed_signal.resize(1024);
	m_FFT.CalculateRealFFT(ref_signal.data(), spectrum.data());
	m_FFT.CalculateRealIFFT(spectrum.data(), reinterpret_cast<FFT_INTERNALS::i32_complex*>(transformed_signal.data()));
	for (size_t i = 0; i < 1024; ++i)
	{
		ASSERT_NEAR(ref_signal[i], 2 * transformed_signal[i], 2048);
	}

	for (size_t i = 0; i < 512; ++i)
	{
		if (i == 0)
		{
			ASSERT_NEAR(spectrum[i].re, 65536 * 16384, 256);
			ASSERT_NEAR(spectrum[i].im, 0, 256);
		}
		else
		{
			ASSERT_NEAR(spectrum[i].re, 0, 32);
			ASSERT_NEAR(spectrum[i].im, 0, 32);
		}
	}
#if 0	
	std::vector<double> signal_ref;
	std::vector<double> signal_transformed;
	for (size_t i = 0; i < 30; ++i)
	{
		signal_ref.push_back(ref_signal[i]);
	}
	for (size_t i = 0; i < 30; ++i)
	{
		signal_transformed.push_back(2 * transformed_signal[i]);
	}
	plot<double>(signal_ref, signal_transformed);
#endif
}

TEST(sample_test_case, int32_DC_test_with_scaling)
{
	CFFT<10, int, int, FFT_INTERNALS::i32_complex> m_FFT(2147483647);

	static CDSPArray<FFT_INTERNALS::i32_complex, 512> spectrum;
	static CDSPArray<int, 1024> ref_signal;
	static CDSPArray<int, 1024> transformed_signal;
	for (size_t i = 0; i < 1024; i++)
	{
		const double v = 65536 * 16384;
		ref_signal.push_back(static_cast<int>(v));
	}
	spectrum.resize(512);
	transformed_signal.resize(1024);
	m_FFT.CalculateRealFFT(ref_signal.data(), spectrum.data());
	m_FFT.CalculateRealIFFT(spectrum.data(), reinterpret_cast<FFT_INTERNALS::i32_complex*>(transformed_signal.data()));
	size_t fft_scaling;
	m_FFT.CalculateRealFFT(ref_signal.data(), spectrum.data(), EN_ScalingMethod::SCALEINPOUT, &fft_scaling);
	size_t ifft_scaling;
	m_FFT.CalculateRealIFFT(spectrum.data(), reinterpret_cast<FFT_INTERNALS::i32_complex*>(transformed_signal.data()), EN_ScalingMethod::SCALEINPOUT, &ifft_scaling);
	DoScaling<10>(fft_scaling, ifft_scaling, transformed_signal);
	for (size_t i = 0; i < 1024; ++i)
	{
		ASSERT_NEAR(ref_signal[i], transformed_signal[i], 2048);
	}
	
	for (size_t i = 0; i < 512; ++i)
	{
		if (i == 0)
		{
			ASSERT_NEAR(spectrum[i].re, ((65536 * 16384 - 1)) << (10 - fft_scaling), 256);
			ASSERT_NEAR(spectrum[i].im, 0, 256);
		}
		else
		{
			ASSERT_NEAR(spectrum[i].re, 0, 32);
			ASSERT_NEAR(spectrum[i].im, 0, 32);
		}
	}
#if 0	
	std::vector<double> signal_ref;
	std::vector<double> signal_transformed;
	for (size_t i = 0; i < 30; ++i)
	{
		signal_ref.push_back(ref_signal[i]);
	}
	for (size_t i = 0; i < 30; ++i)
	{
		signal_transformed.push_back(2 * transformed_signal[i]);
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
#if 0	
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
	return 0;
}
#endif
