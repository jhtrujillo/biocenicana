package org.cenicana.bio.utils;

import java.util.function.DoubleFunction;

/**
 * Mathematical utilities for GWAS, including Brent's optimizer and T-distribution.
 */
public class GwasMathUtils {

    /**
     * Brent's method for finding the minimum of a function in a given interval [a, b].
     */
    public static double brentOptimize(double a, double b, double tol, DoubleFunction<Double> func) {
        double x = a + 0.381966 * (b - a);
        double w = x;
        double v = x;
        double u, delta = 0, d = 0;
        double fx = func.apply(x);
        double fw = fx;
        double fv = fx;

        for (int iter = 0; iter < 100; iter++) {
            double middle = 0.5 * (a + b);
            if (Math.abs(x - middle) <= (2.0 * tol - 0.5 * (b - a))) break;

            if (Math.abs(delta) > tol) {
                double r = (x - w) * (fx - fv);
                double q = (x - v) * (fx - fw);
                double p = (x - v) * q - (x - w) * r;
                q = 2.0 * (q - r);
                if (q > 0) p = -p;
                q = Math.abs(q);
                double temp = delta;
                delta = d;
                if (Math.abs(p) < Math.abs(0.5 * q * temp) && p > q * (a - x) && p < q * (b - x)) {
                    d = p / q;
                    u = x + d;
                    if (u - a < 2.0 * tol || b - u < 2.0 * tol) d = Math.copySign(tol, middle - x);
                } else {
                    delta = (x >= middle) ? a - x : b - x;
                    d = 0.381966 * delta;
                }
            } else {
                delta = (x >= middle) ? a - x : b - x;
                d = 0.381966 * delta;
            }

            u = x + (Math.abs(d) >= tol ? d : Math.copySign(tol, d));
            double fu = func.apply(u);

            if (fu <= fx) {
                if (u >= x) a = x; else b = x;
                v = w; fv = fw; w = x; fw = fx; x = u; fx = fu;
            } else {
                if (u < x) a = u; else b = u;
                if (fu <= fw || w == x) {
                    v = w; fv = fw; w = u; fw = fu;
                } else if (fu <= fv || v == x || v == w) {
                    v = u; fv = fu;
                }
            }
        }
        return x;
    }

    /**
     * Cumulative distribution function (CDF) for the T-distribution via regularized incomplete beta.
     * Falls back to normal approximation for df > 1000.
     */
    public static double tCDF(double t, double df) {
        if (df <= 0) return 0.5;
        // Use normal approximation for large df
        if (df > 1000) return normalCDF(t);
        // I_x(a,b) where x = df/(df+t²), a = df/2, b = 1/2
        return 0.5 * (1.0 + Math.signum(t) * (1.0 - betaIncomplete(df / (df + t * t), df / 2.0, 0.5)));
    }

    /**
     * Regularized Incomplete Beta Function I_x(a, b) approximation.
     */
    public static double betaIncomplete(double x, double a, double b) {
        if (x < 0 || x > 1) return 0;
        if (x == 0) return 0;
        if (x == 1) return 1;

        // Use symmetry identity for numerical stability
        if (x > (a + 1.0) / (a + b + 2.0)) {
            return 1.0 - betaIncomplete(1.0 - x, b, a);
        }

        double eps = 1e-10;
        double fpmin = 1e-30;

        double qab = a + b;
        double qap = a + 1.0;
        double qam = a - 1.0;
        double c = 1.0;
        double d = 1.0 - qab * x / qap;
        if (Math.abs(d) < fpmin) d = fpmin;
        d = 1.0 / d;
        double h = d;

        for (int m = 1; m <= 100; m++) {
            int m2 = 2 * m;
            double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
            d = 1.0 + aa * d;
            if (Math.abs(d) < fpmin) d = fpmin;
            c = 1.0 + aa / c;
            if (Math.abs(c) < fpmin) c = fpmin;
            d = 1.0 / d;
            h *= d * c;
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
            d = 1.0 + aa * d;
            if (Math.abs(d) < fpmin) d = fpmin;
            c = 1.0 + aa / c;
            if (Math.abs(c) < fpmin) c = fpmin;
            d = 1.0 / d;
            double del = d * c;
            h *= del;
            if (Math.abs(del - 1.0) < eps) break;
        }

        double logBeta = logGamma(a) + logGamma(b) - logGamma(a + b);
        double factor = Math.exp(a * Math.log(x) + b * Math.log(1.0 - x) - logBeta) / a;
        return factor * h;
    }

    public static double logGamma(double x) {
        double[] coeff = {76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};
        double y = x;
        double tmp = x + 5.5;
        tmp -= (x + 0.5) * Math.log(tmp);
        double ser = 1.000000000190015;
        for (double c : coeff) ser += c / ++y;
        return -tmp + Math.log(2.5066282746310005 * ser / x);
    }

    public static double normalCDF(double x) {
        return 0.5 * (1.0 + erf(x / Math.sqrt(2.0)));
    }

    public static double erf(double z) {
        double t = 1.0 / (1.0 + 0.5 * Math.abs(z));
        double ans = 1.0 - t * Math.exp(-z * z - 1.26551223 +
                t * (1.00002368 +
                t * (0.37409196 +
                t * (0.09678418 +
                t * (-0.18628806 +
                t * (0.27886807 +
                t * (-1.13520398 +
                t * (1.48851587 +
                t * (-0.82215223 +
                t * 0.17087277)))))))));
        if (z >= 0) return ans;
        else return -ans;
    }
    public static double fCDF(double f, double df1, double df2) {
        if (f <= 0) return 0.0;
        double x = (df1 * f) / (df1 * f + df2);
        return betaIncomplete(x, df1 / 2.0, df2 / 2.0);
    }

    public static double calcMissing(double[] dosages) {
        int missing = 0;
        for (double d : dosages) if (Double.isNaN(d)) missing++;
        return (double) missing / dosages.length;
    }

    public static double calcMaf(double[] dosages, int ploidy) {
        double sum = 0;
        int count = 0;
        for (double d : dosages) {
            if (!Double.isNaN(d)) {
                sum += d;
                count++;
            }
        }
        if (count == 0) return 0;
        double freq = sum / (count * ploidy);
        return Math.min(freq, 1.0 - freq);
    }

    public static double calcGenoFreq(double[] dosages, int ploidy) {
        int n = dosages.length;
        int[] counts = new int[ploidy + 1];
        int total = 0;
        for (double d : dosages) {
            if (!Double.isNaN(d)) {
                int geno = (int) Math.round(d);
                if (geno >= 0 && geno <= ploidy) {
                    counts[geno]++;
                    total++;
                }
            }
        }
        if (total == 0) return 1.0;
        int max = 0;
        for (int c : counts) if (c > max) max = c;
        return (double) max / total;
    }

    public static double calculateR2(double[] x, double[] y) {
        int n = Math.min(x.length, y.length);
        double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, sumY2 = 0;
        int count = 0;
        for (int i = 0; i < n; i++) {
            if (Double.isNaN(x[i]) || Double.isNaN(y[i])) continue;
            sumX += x[i];
            sumY += y[i];
            sumXY += x[i] * y[i];
            sumX2 += x[i] * x[i];
            sumY2 += y[i] * y[i];
            count++;
        }
        if (count < 2) return 0;
        double num = count * sumXY - sumX * sumY;
        double den = Math.sqrt((count * sumX2 - sumX * sumX) * (count * sumY2 - sumY * sumY));
        if (den == 0) return 0;
        double r = num / den;
        return r * r;
    }
}

