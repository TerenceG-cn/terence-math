package com.tce.math.analysis;

/**
 * 单变量可导函数接口
 */
public interface UnivariateDifferentiableFunction extends UnivariateFunction {
    DerivativeStructure value(DerivativeStructure t) throws DimensionMismatchException;
}
