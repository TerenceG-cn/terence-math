package com.tce.math.utils;

/**
 *org.apache.commons.math3.util.FastMath
 */
public class FastMath {
    /** 0x7fffffffffffffffl */
    private static final long MASK_NON_SIGN_LONG = Long.MAX_VALUE;
    /**
     * 0x40000000 - 用来把double值分为两个已清除低位的部分
     * = 2^30.
     */
    private static final long HEX_40000000 = 0x40000000L; // 1073741824L


    /**
     * 求绝对值
     * @param x 求绝对值的double值
     * @return abs(x)
     */
    public static double abs(double x) {
        return Double.longBitsToDouble(MASK_NON_SIGN_LONG & Double.doubleToRawLongBits(x));
    }


    /**
     *arc cos
     * @param x
     * @return arc cosine of x
     */
    public static double acos(double x) {
        //x is NaN
        if (x != x) {
            return Double.NaN;
        }
        //no defined
        if (x > 1.0 || x < -1.0) {
            return Double.NaN;
        }

        if (x == -1.0) {
            return Math.PI;
        }

        if (x == 1.0) {
            return 0.0;
        }

        if (x == 0) {
            return Math.PI/2.0;
        }

        /* Compute acos(x) = atan(sqrt(1-x*x)/x) */

        /* Split x */
        double temp = x * HEX_40000000;
        final double xa = x + temp - temp;//high 32
        final double xb = x - xa;//low 32

        /* Square it */
        double ya = xa*xa;
        double yb = xa*xb*2.0 + xb*xb;

        /* Subtract from 1 */
        ya = -ya;
        yb = -yb;

        double za = 1.0 + ya;
        double zb = -(za - 1.0 - ya);

        temp = za + yb;
        zb += -(temp - za - yb);
        za = temp;

        /* Square root */
        double y = Math.sqrt(za);
        temp = y * HEX_40000000;
        ya = y + temp - temp;
        yb = y - ya;

        /* Extend precision of sqrt */
        yb += (za - ya*ya - 2*ya*yb - yb*yb) / (2.0*y);

        /* Contribution of zb to sqrt */
        yb += zb / (2.0*y);
        y = ya+yb;
        yb = -(y - ya - yb);

        // Compute ratio r = y/x
        double r = y/x;

        // Did r overflow?
        if (Double.isInfinite(r)) { // x is effectively zero
            return Math.PI/2; // so return the appropriate value
        }

        double ra = doubleHighPart(r);
        double rb = r - ra;

        rb += (y - ra*xa - ra*xb - rb*xa - rb*xb) / x;  // Correct for rounding in division
        rb += yb / x;  // Add in effect additional bits of sqrt.

        temp = ra + rb;
        rb = -(temp - ra - rb);
        ra = temp;

        return atan(ra, rb, x<0);
    }

    /**
     * TODO
     * @param x
     * @return
     */
    public static double atan(double x) {
        return atan(x, 0.0, false);
    }

    /** Internal helper function to compute arctangent.
     * @param xa number from which arctangent is requested
     * @param xb extra bits for x (may be 0.0)
     * @param leftPlane if true, result angle must be put in the left half plane
     * @return atan(xa + xb) (or angle shifted by {@code PI} if leftPlane is true)
     */
    private static double atan(double xa, double xb, boolean leftPlane) {
        if (xa == 0.0) { // Matches +/- 0.0; return correct sign
            return leftPlane ? copySign(Math.PI, xa) : xa;
        }

        final boolean negate;
        if (xa < 0) {
            // negative
            xa = -xa;
            xb = -xb;
            negate = true;
        } else {
            negate = false;
        }

        if (xa > 1.633123935319537E16) { // Very large input
            return (negate ^ leftPlane) ? (-Math.PI * F_1_2) : (Math.PI * F_1_2);
        }

        /* Estimate the closest tabulated arctan value, compute eps = xa-tangentTable */
        final int idx;
        if (xa < 1) {
            idx = (int) (((-1.7168146928204136 * xa * xa + 8.0) * xa) + 0.5);
        } else {
            final double oneOverXa = 1 / xa;
            idx = (int) (-((-1.7168146928204136 * oneOverXa * oneOverXa + 8.0) * oneOverXa) + 13.07);
        }

        final double ttA = TANGENT_TABLE_A[idx];
        final double ttB = TANGENT_TABLE_B[idx];

        double epsA = xa - ttA;
        double epsB = -(epsA - xa + ttA);
        epsB += xb - ttB;

        double temp = epsA + epsB;
        epsB = -(temp - epsA - epsB);
        epsA = temp;

        /* Compute eps = eps / (1.0 + xa*tangent) */
        temp = xa * HEX_40000000;
        double ya = xa + temp - temp;
        double yb = xb + xa - ya;
        xa = ya;
        xb += yb;

        //if (idx > 8 || idx == 0)
        if (idx == 0) {
            /* If the slope of the arctan is gentle enough (< 0.45), this approximation will suffice */
            //double denom = 1.0 / (1.0 + xa*tangentTableA[idx] + xb*tangentTableA[idx] + xa*tangentTableB[idx] + xb*tangentTableB[idx]);
            final double denom = 1d / (1d + (xa + xb) * (ttA + ttB));
            //double denom = 1.0 / (1.0 + xa*tangentTableA[idx]);
            ya = epsA * denom;
            yb = epsB * denom;
        } else {
            double temp2 = xa * ttA;
            double za = 1d + temp2;
            double zb = -(za - 1d - temp2);
            temp2 = xb * ttA + xa * ttB;
            temp = za + temp2;
            zb += -(temp - za - temp2);
            za = temp;

            zb += xb * ttB;
            ya = epsA / za;

            temp = ya * HEX_40000000;
            final double yaa = (ya + temp) - temp;
            final double yab = ya - yaa;

            temp = za * HEX_40000000;
            final double zaa = (za + temp) - temp;
            final double zab = za - zaa;

            /* Correct for rounding in division */
            yb = (epsA - yaa * zaa - yaa * zab - yab * zaa - yab * zab) / za;

            yb += -epsA * zb / za / za;
            yb += epsB / za;
        }


        epsA = ya;
        epsB = yb;

        /* Evaluate polynomial */
        final double epsA2 = epsA * epsA;

        /*
    yb = -0.09001346640161823;
    yb = yb * epsA2 + 0.11110718400605211;
    yb = yb * epsA2 + -0.1428571349122913;
    yb = yb * epsA2 + 0.19999999999273194;
    yb = yb * epsA2 + -0.33333333333333093;
    yb = yb * epsA2 * epsA;
         */

        yb = 0.07490822288864472;
        yb = yb * epsA2 - 0.09088450866185192;
        yb = yb * epsA2 + 0.11111095942313305;
        yb = yb * epsA2 - 0.1428571423679182;
        yb = yb * epsA2 + 0.19999999999923582;
        yb = yb * epsA2 - 0.33333333333333287;
        yb = yb * epsA2 * epsA;


        ya = epsA;

        temp = ya + yb;
        yb = -(temp - ya - yb);
        ya = temp;

        /* Add in effect of epsB.   atan'(x) = 1/(1+x^2) */
        yb += epsB / (1d + epsA * epsA);

        final double eighths = EIGHTHS[idx];

        //result = yb + eighths[idx] + ya;
        double za = eighths + ya;
        double zb = -(za - eighths - ya);
        temp = za + yb;
        zb += -(temp - za - yb);
        za = temp;

        double result = za + zb;

        if (leftPlane) {
            // Result is in the left plane
            final double resultb = -(result - za - zb);
            final double pia = 1.5707963267948966 * 2;
            final double pib = 6.123233995736766E-17 * 2;

            za = pia - result;
            zb = -(za - pia + result);
            zb += pib - resultb;

            result = za + zb;
        }


        if (negate ^ leftPlane) {
            result = -result;
        }

        return result;
    }
}
