package com.tce.math.analysis.functions.arithmetic;

import com.tce.math.analysis.BivariateRealFunction;

/**
 * 实数加法
 */
public class Add implements BivariateRealFunction {

    @Override
    public double value(double x, double y) {
        return x+y;
    }
}
