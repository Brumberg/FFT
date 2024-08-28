#pragma once

#ifndef DSPARRAY_
#define DSPARRAY_
#include "carray.h"
#include <limits>
#ifdef _ADI_COMPILER
#include "mathematics.h"
#include <vector.h>
#endif
//#include <utility>

template <typename T, size_t N>
class CDSPArray : public CArray<T, N>
{
	public:
		enum class EN_Scaling{EN_Factor, EN_Shift};
		CDSPArray():CArray<T,N>() {};
		T max() const;
		T max(typename CArray<T, N>::const_iterator const begin, typename CArray<T, N>::const_iterator const end) const;
		T min() const;
		T min(typename CArray<T, N>::const_iterator const begin, typename CArray<T, N>::const_iterator const end) const;
		T absmax() const;
		T absmax(typename CArray<T, N>::const_iterator const begin, typename CArray<T, N>::const_iterator const end) const;
		
		typename CDSPArray<T, N>::iterator findmax();
		typename CDSPArray<T, N>::iterator findmax(typename CArray<T, N>::iterator const begin, typename CArray<T, N>::iterator const end);
		typename CDSPArray<T, N>::iterator findmin();
		typename CDSPArray<T, N>::iterator findmin(typename CArray<T, N>::iterator const begin, typename CArray<T, N>::iterator const end);
		typename CDSPArray<T, N>::iterator findabsmax();
		typename CDSPArray<T, N>::iterator findabsmax(typename CArray<T, N>::iterator const begin, typename CArray<T, N>::iterator const end);

		typename CDSPArray<T, N>::const_iterator findmax(typename CArray<T, N>::const_iterator const begin, typename CArray<T, N>::const_iterator const end) const;
		typename CDSPArray<T, N>::const_iterator findmin(typename CArray<T, N>::const_iterator const begin, typename CArray<T, N>::const_iterator const end) const;
		typename CDSPArray<T, N>::const_iterator findabsmax(typename CArray<T, N>::const_iterator const begin, typename CArray<T, N>::const_iterator const end) const;

		sint32	NormShifts() const;
		sint32	NormFactor() const;
		void	Normalize(sint32 factor, EN_Scaling scaling);
		void	Downscale(sint32 factor, EN_Scaling scaling);
};

template <typename T, size_t N>
T CDSPArray<T, N>::max() const
{
	auto const iter = std::max_element(this->cbegin(), this->cend());
	return *iter;
}

template <typename T, size_t N>
T CDSPArray<T, N>::max(typename CArray<T, N>::const_iterator const begin, typename CArray<T, N>::const_iterator const end) const
{
	auto const iter = std::max_element(begin, end);
	return *iter;
}

template <typename T, size_t N>
T CDSPArray<T, N>::min() const
{
	auto iter = std::min_element(this->begin(), this->end());
	return *iter;
}

template <typename T, size_t N>
T CDSPArray<T, N>::min(typename CArray<T, N>::const_iterator const begin, typename CArray<T, N>::const_iterator const end) const
{
	auto iter = std::min_element(begin, end);
	return *iter;
}

template <typename T, size_t N>
T CDSPArray<T, N>::absmax() const
{
	T absmax_ = 0;
	if (this->cbegin() != this->cend())
	{
		absmax_ = static_cast<T>(abs(*this->cbegin()));
		for (auto iter = this->cbegin(); iter != this->cend(); ++iter)
		{
			const T item = static_cast<T>(abs(*iter));
			absmax_ = (absmax_ < item) ? item : absmax_;
		}
	}
	else
	{
		//this is an exception
	}
	return absmax_;
}

template <typename T, size_t N>
T CDSPArray<T, N>::absmax(typename CArray<T, N>::const_iterator const begin, typename CArray<T, N>::const_iterator const end) const
{
	T absmax_ = 0;
	for (auto iter = begin; iter != end; ++iter)
	{
		const T item = static_cast<T>(abs(*iter));
		absmax_ = (absmax_ < item) ? item : absmax_;
	}
	return absmax_;
}

template <typename T, size_t N>
typename CDSPArray<T, N>::iterator CDSPArray<T, N>::findmax()
{
	return std::max_element(this->begin(), this->end());
}

template <typename T, size_t N>
typename CDSPArray<T, N>::iterator CDSPArray<T, N>::findmax(typename CArray<T, N>::iterator const begin, typename CArray<T, N>::iterator const end)
{
	return std::max_element(begin, end);
}

template <typename T, size_t N>
typename CDSPArray<T, N>::iterator CDSPArray<T, N>::findmin()
{
	return std::min_element(this->begin(), this->end());
}

template <typename T, size_t N>
typename CDSPArray<T, N>::iterator CDSPArray<T, N>::findmin(typename CArray<T, N>::iterator const begin, typename CArray<T, N>::iterator const end)
{
	return std::min_element(begin, end);
}

template <typename T, size_t N>
typename CDSPArray<T, N>::iterator CDSPArray<T, N>::findabsmax()
{
	T val = abs(*(this->begin()));
	typename CDSPArray<T, N>::iterator Ix = this->begin();
	for (auto& i : *this)
	{
		T absval = abs(i);
		if (val < absval)
		{
			val = absval;
			Ix = &i;
		}
	}
	return Ix;
}

template <typename T, size_t N>
typename CDSPArray<T, N>::iterator CDSPArray<T, N>::findabsmax(typename CArray<T, N>::iterator const begin, typename CArray<T, N>::iterator const end)
{
	T val = abs(*begin);
	typename CDSPArray<T, N>::iterator Ix = begin;
	for (; begin != end; ++begin)
	{
		T absval = abs(*begin);
		if (val < absval)
		{
			val = absval;
			Ix = begin;
		}
	}
	return Ix;
}

template <typename T, size_t N>
typename CDSPArray<T, N>::const_iterator CDSPArray<T, N>::findmax(typename CArray<T, N>::const_iterator const begin, typename CArray<T, N>::const_iterator const end) const
{
	return std::max_element(begin, end);
}

template <typename T, size_t N>
typename CDSPArray<T, N>::const_iterator CDSPArray<T, N>::findmin(typename CArray<T, N>::const_iterator const begin, typename CArray<T, N>::const_iterator const end) const
{
	return std::min_element(begin, end);
}

template <typename T, size_t N>
typename CDSPArray<T, N>::const_iterator CDSPArray<T, N>::findabsmax(typename CArray<T, N>::const_iterator const begin, typename CArray<T, N>::const_iterator const end) const
{
	T val = abs(*begin);
	typename CDSPArray<T, N>::const_iterator Ix = begin;
	for (; begin != end; ++begin)
	{
		T absval = abs(*begin);
		if (val < absval)
		{
			val = absval;
			Ix = begin;
		}
	}
	return Ix;
}

template <typename T, size_t N>
sint32 CDSPArray<T, N>::NormShifts() const
{
	T val = absmax();
	sint32 reqshifts = 0;
	if (val != 0)
	{
		//val is a signed type
		static const T maxrangedb = static_cast<T>((8u * sizeof(T)) - 2u);
		static const T maxrange = static_cast<T>(1)<<maxrangedb;
		while (val < maxrange)
		{
			val *= 2;
			reqshifts++;
		}
	}
	return reqshifts;
}

template <typename T, size_t N>
sint32 CDSPArray<T, N>::NormFactor() const
{
	const T val = absmax();
	sint32 normfactor;

	if (val>0)
	{
		normfactor = INT_MAX / val;
	}
	else
	{
		const size_t shifts = ((sizeof(sint32) * 8u)) - ((sizeof(T) * 8u));
		normfactor = (static_cast<sint32>(1) << shifts);
	}
	return normfactor;
}

template <typename T, size_t N>
void CDSPArray<T, N>::Normalize(sint32 factor, EN_Scaling scaling)
{
	static const sint32 downscaling = static_cast<sint32>(((sizeof(sint32) * 8u)) - ((sizeof(T) * 8u)));//lint !e778 BOL: evaluation to zero is OK
	if (scaling == EN_Scaling::EN_Factor)
	{
		const sint32 downscalingfactor = static_cast<sint32>(1) << downscaling;
		for (auto& i : *this)
		{
			const sint32 imr = (static_cast<sint32>(i) * factor);
			//i = static_cast<T>(imr >> downscaling);
			i = static_cast<T>(imr/downscalingfactor);
		}
	}
	else
	{
		for (auto& i : *this)
		{
			i = i<<factor;
		}
	}
}

template <typename T, size_t N>
void CDSPArray<T, N>::Downscale(sint32 factor, EN_Scaling scaling)
{
	if (scaling == EN_Scaling::EN_Factor)
	{
			for (auto& i : *this)
			{
				const sint32 imr = (static_cast<sint32>(i) / factor);
				i = imr;
			}
	}
	else
	{
		for (auto& i : *this)
		{
			i = i>>factor;
		}
	}
}


#ifdef _ADI_COMPILER
//specialized class
template <size_t N>
class CDSPArray<sint16, N> : public CArray<sint16, N>
{
	public:
		enum class EN_Scaling{EN_Factor, EN_Shift};
		CDSPArray():CArray<sint16,N>() {};
		sint16 max() const;
		sint16 max(typename CArray<sint16, N>::const_iterator const begin, typename CArray<sint16, N>::const_iterator const end) const;//lint !e1746
		sint16 min() const;
		sint16 min(typename CArray<sint16, N>::const_iterator const begin, typename CArray<sint16, N>::const_iterator const end) const;//lint !e1746
		sint16 absmax() const;
		sint16 absmax(typename CArray<sint16, N>::const_iterator const begin, typename CArray<sint16, N>::const_iterator const end) const;//lint !e1746
		typename CDSPArray<sint16, N>::iterator findmax();
		typename CDSPArray<sint16, N>::iterator findmax(typename CArray<sint16, N>::iterator const begin, typename CArray<sint16, N>::iterator const end);
		typename CDSPArray<sint16, N>::iterator findmin();
		typename CDSPArray<sint16, N>::iterator findmin(typename CArray<sint16, N>::iterator const begin, typename CArray<sint16, N>::iterator const end);
		typename CDSPArray<sint16, N>::iterator findabsmax();
		typename CDSPArray<sint16, N>::iterator findabsmax(typename CArray<sint16, N>::iterator const begin, typename CArray<sint16, N>::iterator const end);

		typename CDSPArray<sint16, N>::const_iterator findmax(typename CArray<sint16, N>::const_iterator const begin, typename CArray<sint16, N>::const_iterator const end) const;
		typename CDSPArray<sint16, N>::const_iterator findmin(typename CArray<sint16, N>::const_iterator const begin, typename CArray<sint16, N>::const_iterator const end) const;
		typename CDSPArray<sint16, N>::const_iterator findabsmax(typename CArray<sint16, N>::const_iterator const begin, typename CArray<sint16, N>::const_iterator const end) const;

		sint32	NormShifts() const;
		sint32	NormFactor() const;
		void	Normalize(sint32 factor, EN_Scaling scaling);
		void	Downscale(sint32 factor, EN_Scaling scaling);
};

template <size_t N>
sint16 CDSPArray<sint16, N>::max() const
{
	return vecmax_fr16(this->data(), static_cast<sint32>(this->size()));
}

template <size_t N>
sint16 CDSPArray<sint16, N>::max(typename CArray<sint16, N>::const_iterator const begin, typename CArray<sint16, N>::const_iterator const end) const
{
	return vecmax_fr16(const_cast<sint16*>(begin.operator->()) , static_cast<sint32>(begin.diff(end)));//lint !e9005
}

template <size_t N>
sint16 CDSPArray<sint16, N>::min() const
{
	return vecmin_fr16(this->data(), static_cast<sint32>(this->size()));
}

template <size_t N>
sint16 CDSPArray<sint16, N>::min(typename CArray<sint16, N>::const_iterator const begin, typename CArray<sint16, N>::const_iterator const end) const
{
	return vecmin_fr16(static_cast<sint16*>(begin.operator->()), static_cast<sint32>(begin.diff(end)));
}

template <size_t N>
sint16 CDSPArray<sint16, N>::absmax() const
{
	sint16 absmax_ = 0;
	for (auto iter = this->cbegin(); iter != this->cend(); ++iter)
	{
		const sint16 item = static_cast<sint16>(abs(*iter));
		absmax_ = (absmax_ < item) ? item : absmax_;
	}
	return absmax_;
}

template <size_t N>
sint16 CDSPArray<sint16, N>::absmax(typename CArray<sint16, N>::const_iterator const begin, typename CArray<sint16, N>::const_iterator const end) const
{
	sint16 absmax_ = 0;
	for (auto iter = begin; iter != end; ++iter)
	{
		const sint16 item = static_cast<sint16>(abs(*iter));
		absmax_ = (absmax_ < item) ? item : absmax_;
	}
	return absmax_;
}

template <size_t N>
typename CDSPArray<sint16, N>::iterator CDSPArray<sint16, N>::findmax()
{
	const typename CDSPArray<sint16, N>::iterator begin = this->begin();
	const typename CDSPArray<sint16, N>::iterator end = this->end();
	const size_t ix = static_cast<size_t>(vecmaxloc_fr16(this->data(), end.diff(begin)));
	return begin+ix;
}

template <size_t N>
typename CDSPArray<sint16, N>::iterator CDSPArray<sint16, N>::findmax(typename CArray<sint16, N>::iterator const begin, typename CArray<sint16, N>::iterator const end)
{
	const size_t ix = static_cast<size_t>(vecmaxloc_fr16(static_cast<sint16*>(begin.operator->()), end.diff(begin)));
	return begin+ix;
}

template <size_t N>
typename CDSPArray<sint16, N>::iterator CDSPArray<sint16, N>::findmin()
{
	const size_t ix = static_cast<size_t>(vecminloc_fr16(this->data(), this->end().diff(this->begin())));
	return this->begin()+ix;
}

template <size_t N>
typename CDSPArray<sint16, N>::iterator findmin(typename CArray<sint16, N>::iterator const begin, typename CArray<sint16, N>::iterator const end)
{
	const size_t ix = static_cast<size_t>(vecminloc_fr16(static_cast<sint16*>(begin.operator->()), end.diff(begin)));
	return begin + ix;
}

template <size_t N>
typename CDSPArray<sint16, N>::iterator CDSPArray<sint16, N>::findabsmax()
{
	typename CDSPArray<sint16, N>::iterator begin = this->begin();
	typename CDSPArray<sint16, N>::iterator const end = this->end();

	sint16 val = static_cast<sint16>(abs(*begin));
	typename CDSPArray<sint16, N>::iterator Ix = begin;
	for (; begin != end; ++begin)
	{
		const sint16 absval = static_cast<sint16>(abs(*begin));
		if (val < absval)
		{
			val = absval;
			Ix = begin;
		}
	}
	return Ix;
}

template <size_t N>
typename CDSPArray<sint16, N>::iterator CDSPArray<sint16, N>::findabsmax(typename CArray<sint16, N>::iterator const begin, typename CArray<sint16, N>::iterator const end)
{
	sint16 val = abs(*begin);
	typename CDSPArray<sint16, N>::iterator Ix = begin;
	for (; begin != end; ++begin)
	{
		sint16 absval = abs(*begin);
		if (val < absval)
		{
			val = absval;
			Ix = begin;
		}
	}
	return Ix;
}

template <size_t N>
typename CDSPArray<sint16, N>::const_iterator CDSPArray<sint16, N>::findmax(typename CArray<sint16, N>::const_iterator const begin, typename CArray<sint16, N>::const_iterator const end) const
{
	const size_t ix = static_cast<size_t>(vecmaxloc_fr16(const_cast<sint16*>(begin.operator->()), end.diff(begin)));//lint !e9005
	return begin + ix;
}

template <size_t N>
typename CDSPArray<sint16, N>::const_iterator CDSPArray<sint16, N>::findmin(typename CArray<sint16, N>::const_iterator const begin, typename CArray<sint16, N>::const_iterator const end) const
{
	const size_t ix = static_cast<size_t>(vecminloc_fr16(static_cast<sint16*>(begin.operator->()), end.diff(begin)));
	return begin + ix;
}

template <size_t N>
typename CDSPArray<sint16, N>::const_iterator CDSPArray<sint16, N>::findabsmax(typename CArray<sint16, N>::const_iterator const begin, typename CArray<sint16, N>::const_iterator const end) const
{
	sint16 val = abs(*begin);
	typename CDSPArray<sint16, N>::const_iterator Ix = begin;
	for (; begin != end; ++begin)
	{
		sint16 absval = abs(*begin);
		if (val < absval)
		{
			val = absval;
			Ix = begin;
		}
	}
	return Ix;
}

template <size_t N>
sint32 CDSPArray<sint16, N>::NormShifts() const
{
	const sint16 val = absmax();
	return norm_fr1x16(val);
}

template <size_t N>
sint32 CDSPArray<sint16, N>::NormFactor() const
{
	const sint16 val = absmax();
	sint32 normfactor;

	if (val>0)
	{
		normfactor = INT_MAX / val;
	}
	else
	{
		const uint32 shifts = ((sizeof(sint32) * 8u)) - ((sizeof(sint16) * 8u));
		normfactor = (static_cast<sint32>(1) << shifts);
	}
	return normfactor;
}

template <size_t N>
void CDSPArray<sint16, N>::Normalize(sint32 factor, EN_Scaling scaling)
{
	static const sint32 downscaling = static_cast<sint32>(((sizeof(sint32) * 8u)) - ((sizeof(sint16) * 8u)));
	if (scaling == EN_Scaling::EN_Factor)
	{
		const sint32 downscalingfactor = static_cast<sint32>(1) << downscaling;
		for (auto& i : *this)
		{
			const sint32 imr = (static_cast<sint32>(i) * factor);
			//i = static_cast<T>(imr >> downscaling);
			i = static_cast<sint16>(imr/downscalingfactor);
		}
	}
	else
	{
		vecsnorm(static_cast<size_t>(this->end().diff(this->begin())),factor,this->data());
	}
}

template <size_t N>
void CDSPArray<sint16, N>::Downscale(sint32 factor, EN_Scaling scaling)
{
	if (scaling == EN_Scaling::EN_Factor)
	{
		for (auto& i : *this)
		{
			const sint32 imr = (static_cast<sint32>(i) / factor);
			i = static_cast<sint16>(imr);
		}
	}
	else
	{
		for (auto& i : *this)
		{
			i = i>>factor;
		}
	}
}
#endif
#endif //DSPARRAY_
