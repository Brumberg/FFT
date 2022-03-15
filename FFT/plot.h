#pragma once
#include <vector>
#include <utility>
#include <boost/tuple/tuple.hpp>
#include "gnuplot-iostream/gnuplot-iostream.h"
#include <algorithm>
#undef max
#undef min

template <typename T>
void plot(std::vector<T> signal)
{
    std::vector<std::pair<double, T>> signal2plot;
    size_t k = 0;
    T absmaxval = 0;
    for (const auto& i : signal)
    {
        signal2plot.push_back(std::pair<double, T>(static_cast<double>(k), i));
        absmaxval = std::max(absmaxval, abs(i));
        ++k;
    }

    double absmaximum = static_cast<double>(absmaxval);
    double range = pow(10.f, floor(log10(absmaximum)));
    double scale = /*std::ceil*/(absmaximum / range);
    if (scale > 5.)
    {
        range = 10.f * range;
    }
    else if (scale > 2)
    {
        range = 5.f * range;
    }
    else
    {
        range = 2.f * range;
    }


    Gnuplot gp;
    gp << "clear;\n";
    gp << "set xrange [0:" << signal2plot.size() << "]\nset yrange [" << -range << ":" << range << "]\n";
    gp << "plot '-' with lines title 'real'\n";

    gp.send1d(signal2plot);
    gp << "set xlabel 'Time [Samples]' font 'Helvetica, 12'\n" <<
        "set ylabel 'Amplitude' font 'Helvetica, 12'" << std::endl;
    gp << "set title 'Visualization of signal pair' font 'Helvetica, 13'" << std::endl;
    gp << "set tics font 'Helvetica, 10'" << std::endl;

    gp << "set grid xtics lt 0 lw 1 lc 'gray'\n" <<
        "set grid ytics lt 0 lw 2 lc 'gray'" << std::endl;
    gp << "replot;" << std::endl;
}

template <typename T>
void plot(std::vector<T> signal1, std::vector<T> signal2)
{
    std::vector<std::pair<double, T>> signal2plot1;
    std::vector<std::pair<double, T>> signal2plot2;
    size_t k = 0;
    T absmaxval = 0;
    for (const auto& i : signal1)
    {
        signal2plot1.push_back(std::pair<double, T>(static_cast<double>(k), i));
        absmaxval = std::max(absmaxval, abs(i));
        ++k;
    }
    k = 0;
    for (const auto& i : signal2)
    {
        signal2plot2.push_back(std::pair<double, T>(static_cast<double>(k), i));
        absmaxval = std::max(absmaxval, abs(i));
        ++k;
    }

    double absmaximum = static_cast<double>(absmaxval);
    double range = pow(10.f, floor(log10(absmaximum)));
    double scale = /*std::ceil*/(absmaximum / range);
    if (scale > 5.)
    {
        range = 10.f * range;
    }
    else if (scale > 2)
    {
        range = 5.f * range;
    }
    else
    {
        range = 2.f * range;
    }


    Gnuplot gp;
    gp << "clear;\n";
    gp << "set xrange [0:" << signal2plot1.size() << "]\nset yrange [" << -range << ":" << range << "]\n";
    gp << "plot '-' with lines title 'up', '-' with lines title 'dn'\n";

    gp.send1d(signal2plot1);
    gp.send1d(signal2plot2);
    gp << "set xlabel 'Time [Samples]' font 'Helvetica, 12'\n" <<
        "set ylabel 'Amplitude' font 'Helvetica, 12'" << std::endl;
    gp << "set title 'Visualization of signal pair' font 'Helvetica, 13'" << std::endl;
    gp << "set tics font 'Helvetica, 10'" << std::endl;

    gp << "set grid xtics lt 0 lw 1 lc 'gray'\n" <<
        "set grid ytics lt 0 lw 2 lc 'gray'" << std::endl;
    gp << "replot;" << std::endl;
}