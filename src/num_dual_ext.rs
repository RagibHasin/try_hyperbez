use std::f64;

use nalgebra::{allocator::Allocator, DefaultAllocator, Dim};
use num_dual::{Dual2Vec, DualNum, DualVec};

pub trait DualNumExt<T: Copy>: DualNum<T> + Copy {
    fn re_mut(&mut self) -> &mut T;
}

impl<T, F, D> DualNumExt<T> for DualVec<T, F, D>
where
    T: DualNum<F> + Copy,
    D: Dim,
    DualVec<T, F, D>: DualNum<T> + Copy,
    DefaultAllocator: Allocator<D>,
{
    fn re_mut(&mut self) -> &mut T {
        &mut self.re
    }
}

impl<T, F, D> DualNumExt<T> for Dual2Vec<T, F, D>
where
    T: DualNum<F> + Copy,
    D: Dim,
    Dual2Vec<T, F, D>: DualNum<T> + Copy,
    DefaultAllocator: Allocator<D, D> + Allocator<nalgebra::Const<1>, D>,
    num_dual::Derivative<T, F, nalgebra::Const<1>, D>: Copy,
{
    fn re_mut(&mut self) -> &mut T {
        &mut self.re
    }
}

pub fn norm_radians<D: DualNumExt<f64>>(mut theta: D) -> D {
    let re = theta.re_mut();
    *re = re.rem_euclid(f64::consts::TAU);
    if *re > f64::consts::PI {
        *re -= f64::consts::TAU;
    } else if *re < -f64::consts::PI {
        *re += f64::consts::TAU;
    }
    theta
}
