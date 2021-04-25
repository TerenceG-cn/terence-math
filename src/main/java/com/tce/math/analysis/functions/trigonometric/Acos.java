package com.tce.math.analysis.functions.trigonometric;


import com.tce.math.analysis.DerivativeStructure;
import com.tce.math.analysis.DimensionMismatchException;
import com.tce.math.analysis.UnivariateDifferentiableFunction;
import com.tce.math.analysis.UnivariateFunction;
import com.tce.math.utils.FastMath;

public class Acos implements UnivariateDifferentiableFunction {
    @Override
    public DerivativeStructure value(DerivativeStructure t) throws DimensionMismatchException {
        return null;
    }

    @Override
    public double value(double x) {
        return FastMath.acos(x);
    }
}
