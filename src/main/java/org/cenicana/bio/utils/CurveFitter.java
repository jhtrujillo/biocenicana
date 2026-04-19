package org.cenicana.bio.utils;

/**
 * Fits 5 decay models to binned LD data and selects the best by coefficient of determination (R²).
 * Models: Exponential, Power Law, Hyperbolic, Polynomial-2, Polynomial-3.
 */
public class CurveFitter {

    public String name;
    public String equation;
    public double r2;
    public double[] fitX;
    public double[] fitY;
    /** Half-decay distance computed analytically from the fitted model equation. -1 if not computable. */
    public long analyticalHalfDecay = -1;

    /** Entry point: given raw binned data, try all models and return the best. */
    public static CurveFitter selectBest(double[] xd, double[] yd, int numPoints) {
        int n = xd.length;
        if (n < 3) return fallback(xd, numPoints);

        double maxX = xd[n - 1];
        double meanY = 0;
        for (double y : yd) meanY += y;
        meanY /= n;
        double ssTot = 0;
        for (double y : yd) ssTot += (y - meanY) * (y - meanY);

        double[] evalX = new double[numPoints + 1];
        for (int i = 0; i <= numPoints; i++) evalX[i] = maxX * i / numPoints;

        CurveFitter best = null;

        // ── 1. Exponential: r² = a * exp(-b * d) ────────────────────────────
        {
            double sx = 0, sy = 0, sxy = 0, sxx = 0;
            int m = 0;
            for (int i = 0; i < n; i++) {
                if (yd[i] > 0) {
                    double x = xd[i], y = Math.log(yd[i]);
                    sx += x; sy += y; sxy += x * y; sxx += x * x; m++;
                }
            }
            if (m >= 2) {
                double den = m * sxx - sx * sx;
                double b = den != 0 ? Math.max(0, -(m * sxy - sx * sy) / den) : 0;
                double a = Math.exp((sy + b * sx) / m);
                double[] fy = new double[numPoints + 1];
                for (int i = 0; i <= numPoints; i++) fy[i] = clamp(a * Math.exp(-b * evalX[i]));
                double r2 = computeR2(xd, yd, ssTot, d -> a * Math.exp(-b * d));
                String eq = String.format(java.util.Locale.US, "r&sup2; = %.4f &times; e<sup>&minus;%.3e &times; d</sup>", a, b);
                // Analytical half-decay: a*exp(-b*d) = a/2 → d = ln(2)/b
                long hd = b > 0 ? Math.round(Math.log(2) / b) : -1L;
                best = keep(best, "Exponential", eq, r2, evalX, fy, hd);
            }
        }

        // ── 2. Power: r² = a * d^(-b) ────────────────────────────────────────
        {
            double sx = 0, sy = 0, sxy = 0, sxx = 0;
            int m = 0;
            for (int i = 0; i < n; i++) {
                if (xd[i] > 0 && yd[i] > 0) {
                    double x = Math.log(xd[i]), y = Math.log(yd[i]);
                    sx += x; sy += y; sxy += x * y; sxx += x * x; m++;
                }
            }
            if (m >= 2) {
                double den = m * sxx - sx * sx;
                double b = den != 0 ? Math.max(0, -(m * sxy - sx * sy) / den) : 0;
                double a = Math.exp((sy + b * sx) / m);
                double[] fy = new double[numPoints + 1];
                for (int i = 0; i <= numPoints; i++)
                    fy[i] = evalX[i] > 0 ? clamp(a * Math.pow(evalX[i], -b)) : clamp(a);
                double r2 = computeR2(xd, yd, ssTot, d -> d > 0 ? a * Math.pow(d, -b) : a);
                String eq = String.format(java.util.Locale.US, "r&sup2; = %.4f &times; d<sup>&minus;%.4f</sup>", a, b);
                // Analytical half-decay: a*d^(-b) = a/2 → d = 2^(1/b)
                long hd = b > 0 ? Math.round(Math.pow(2, 1.0 / b)) : -1L;
                best = keep(best, "Power Law", eq, r2, evalX, fy, hd);
            }
        }

        // ── 3. Hyperbolic: r² = c / (1 + b*d) ──────────────────────────────
        {
            double sx = 0, sy = 0, sxy = 0, sxx = 0;
            int m = 0;
            for (int i = 0; i < n; i++) {
                if (yd[i] > 0) {
                    double x = xd[i], y = 1.0 / yd[i];
                    sx += x; sy += y; sxy += x * y; sxx += x * x; m++;
                }
            }
            if (m >= 2) {
                double den = m * sxx - sx * sx;
                double slope = den != 0 ? (m * sxy - sx * sy) / den : 0;
                double intercept = (sy - slope * sx) / m;
                double c = intercept > 0 ? 1.0 / intercept : 1.0;
                double b = Math.max(0, slope * c);
                double[] fy = new double[numPoints + 1];
                for (int i = 0; i <= numPoints; i++) fy[i] = clamp(c / (1 + b * evalX[i]));
                double r2 = computeR2(xd, yd, ssTot, d -> c / (1 + b * d));
                String eq = String.format(java.util.Locale.US, "r&sup2; = %.4f / (1 + %.3e &times; d)", c, b);
                // Analytical half-decay: c/(1+b*d) = c/2 → d = 1/b
                long hd = b > 0 ? Math.round(1.0 / b) : -1L;
                best = keep(best, "Hyperbolic", eq, r2, evalX, fy, hd);
            }
        }

        // ── 4. Polynomial degree 2 (normalized) ─────────────────────────────
        {
            double[] t = normalize(xd, maxX);
            double[] cf = polyFit(t, yd, 2);
            if (cf != null) {
                double[] fy = new double[numPoints + 1];
                for (int i = 0; i <= numPoints; i++) {
                    double ti = evalX[i] / maxX;
                    fy[i] = clamp(cf[0] + cf[1] * ti + cf[2] * ti * ti);
                }
                double r2 = computeR2(xd, yd, ssTot, d -> {
                    double ti = d / maxX;
                    return cf[0] + cf[1] * ti + cf[2] * ti * ti;
                });
                double b1 = cf[1] / maxX, b2 = cf[2] / (maxX * maxX);
                String eq = String.format(java.util.Locale.US, "r&sup2; = %.4f + %.3e&thinsp;d + %.3e&thinsp;d&sup2;", cf[0], b1, b2);
                // Half-decay for poly-2: binary search on the fitted curve
                double r2max_p2 = cf[0]; // approx at d=0
                long hd = binarySearchHalfDecay(d -> { double ti=d/maxX; return cf[0]+cf[1]*ti+cf[2]*ti*ti; }, r2max_p2, maxX);
                best = keep(best, "Polynomial (degree 2)", eq, r2, evalX, fy, hd);
            }
        }

        // ── 5. Polynomial degree 3 (normalized) ─────────────────────────────
        {
            double[] t = normalize(xd, maxX);
            double[] cf = polyFit(t, yd, 3);
            if (cf != null) {
                double[] fy = new double[numPoints + 1];
                for (int i = 0; i <= numPoints; i++) {
                    double ti = evalX[i] / maxX;
                    fy[i] = clamp(cf[0] + cf[1]*ti + cf[2]*ti*ti + cf[3]*ti*ti*ti);
                }
                double r2 = computeR2(xd, yd, ssTot, d -> {
                    double ti = d / maxX;
                    return cf[0] + cf[1]*ti + cf[2]*ti*ti + cf[3]*ti*ti*ti;
                });
                double b1 = cf[1]/maxX, b2 = cf[2]/(maxX*maxX), b3 = cf[3]/(maxX*maxX*maxX);
                String eq = String.format(java.util.Locale.US, "r&sup2; = %.4f + %.3e&thinsp;d + %.3e&thinsp;d&sup2; + %.3e&thinsp;d&sup3;", cf[0], b1, b2, b3);
                double r2max_p3 = cf[0];
                long hd = binarySearchHalfDecay(d -> { double ti=d/maxX; return cf[0]+cf[1]*ti+cf[2]*ti*ti+cf[3]*ti*ti*ti; }, r2max_p3, maxX);
                best = keep(best, "Polynomial (degree 3)", eq, r2, evalX, fy, hd);
            }
        }

        return best != null ? best : fallback(xd, numPoints);
    }

    // ── Helpers ──────────────────────────────────────────────────────────────

    private static CurveFitter keep(CurveFitter best, String name, String eq, double r2, double[] fx, double[] fy, long analyticalHd) {
        if (best == null || r2 > best.r2) {
            CurveFitter c = new CurveFitter();
            c.name = name; c.equation = eq; c.r2 = r2; c.fitX = fx; c.fitY = fy;
            c.analyticalHalfDecay = analyticalHd;
            return c;
        }
        return best;
    }

    /** Binary search for the distance where predict(d) drops to target/2. */
    private static long binarySearchHalfDecay(Predict pred, double r2AtZero, double maxX) {
        double target = r2AtZero / 2.0;
        double lo = 0, hi = maxX;
        for (int i = 0; i < 60; i++) {
            double mid = (lo + hi) / 2.0;
            if (pred.predict(mid) > target) lo = mid; else hi = mid;
        }
        double result = (lo + hi) / 2.0;
        return result > 0 && result < maxX ? Math.round(result) : -1L;
    }

    @FunctionalInterface
    private interface Predict { double predict(double x); }

    private static double computeR2(double[] xd, double[] yd, double ssTot, Predict pred) {
        if (ssTot <= 0) return 0;
        double ssRes = 0;
        for (int i = 0; i < xd.length; i++) { double e = yd[i] - pred.predict(xd[i]); ssRes += e * e; }
        return 1.0 - ssRes / ssTot;
    }

    private static double[] normalize(double[] xd, double maxX) {
        double[] t = new double[xd.length];
        for (int i = 0; i < xd.length; i++) t[i] = xd[i] / maxX;
        return t;
    }

    private static double clamp(double v) { return Math.max(0, Math.min(1.0, v)); }

    /** Fit polynomial of given degree using normal equations with Gaussian elimination. */
    private static double[] polyFit(double[] x, double[] y, int deg) {
        int n = x.length, d = deg + 1;
        double[][] A = new double[d][d];
        double[] b = new double[d];
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                double s = 0; for (double xi : x) s += Math.pow(xi, i + j); A[i][j] = s;
            }
            double s = 0; for (int k = 0; k < n; k++) s += Math.pow(x[k], i) * y[k]; b[i] = s;
        }
        return gaussianElimination(A, b);
    }

    private static double[] gaussianElimination(double[][] A, double[] b) {
        int n = A.length;
        double[][] aug = new double[n][n + 1];
        for (int i = 0; i < n; i++) { for (int j = 0; j < n; j++) aug[i][j] = A[i][j]; aug[i][n] = b[i]; }
        for (int k = 0; k < n; k++) {
            int max = k;
            for (int i = k + 1; i < n; i++) if (Math.abs(aug[i][k]) > Math.abs(aug[max][k])) max = i;
            double[] tmp = aug[k]; aug[k] = aug[max]; aug[max] = tmp;
            if (Math.abs(aug[k][k]) < 1e-14) return null;
            for (int i = k + 1; i < n; i++) {
                double f = aug[i][k] / aug[k][k];
                for (int j = k; j <= n; j++) aug[i][j] -= f * aug[k][j];
            }
        }
        double[] sol = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            sol[i] = aug[i][n];
            for (int j = i + 1; j < n; j++) sol[i] -= aug[i][j] * sol[j];
            sol[i] /= aug[i][i];
        }
        return sol;
    }

    private static CurveFitter fallback(double[] xd, int numPoints) {
        CurveFitter c = new CurveFitter();
        c.name = "No fit"; c.equation = "N/A"; c.r2 = 0;
        c.fitX = new double[numPoints + 1]; c.fitY = new double[numPoints + 1];
        return c;
    }
}
