#ifndef B_SPLINE_CXX_H
#define B_SPLINE_CXX_H


namespace BSpline{
    // Template voodoo
    // Direct implemetation of expressions given on
    // https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node16.html
    // Computes the kth order B-spline via de Boor recursion.
    template <unsigned k>
    double BSpline(double x, const double *t)
    {
        if (*t <= x && x < *(t+k))
        {
            double a=0;
            double h = *(t+k-1)-*t;
            a += (fabs(h) > 1e-16) ? BSpline<k-1>(x, t) * (x - *t) / h : 0;
            h = (*(t+k) - *(t+1));
            a += (fabs(h) > 1e-16) ? BSpline<k-1>(x, (t+1)) * (*(t+k) - x) / h : 0;
            return a;
        }
        else
            return 0;
    }

    template <>
    double BSpline<1>(double x, const double *t)
    {
        if (*t <= x && x < *(t+1))
            return 1.;
        else
            return 0.;
    }

    template <unsigned k>
    double DBSpline(double x, const double* t)
    {
        if (*t <= x && x < *(t+k))
        {
            double a=0;
            double h = *(t+k-1)-*t;
            a += (fabs(h) > 1e-16) ? BSpline<k-1>(x, t) * (k - 1) / h : 0;
            h = (*(t+k) - *(t+1));
            a -= (fabs(h) > 1e-16) ? BSpline<k-1>(x, (t+1)) * (k - 1) / h : 0;
            return a;
        }
        else
            return 0;
    }

    template <>
    double DBSpline<1>(double x, const double *t)
    {
        return 0.;
    }
};

#endif /* end of include guard: B_SPLINE_CXX_H */
