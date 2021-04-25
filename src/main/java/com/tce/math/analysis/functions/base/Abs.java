package com.tce.math.analysis.functions.base;


import com.tce.math.analysis.UnivariateFunction;
import com.tce.math.utils.FastMath;

/**
 *绝对值
 */
public class Abs implements UnivariateFunction {

    @Override
    public double value(double x) {
        return FastMath.abs(x);
    }
}
