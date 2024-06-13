use std::iter;
use std::ops::Div;
use num_traits::NumAssign;


pub fn mean<T, I>(xs: I) -> T
where T: NumAssign + Copy + Div<f64, Output=T> + iter::Sum<T>,
      I: AsRef<[T]>,
{
    let xs = xs.as_ref();
    let len = xs.len();
    assert!(len > 0);
    xs.iter().cloned().sum::<T>() / len as f64
}


/// Perform cumulative sum.
pub fn cumsum<'a, T, I, R>(xs: I) -> R
where T: NumAssign + Copy + 'a,
      I: Iterator<Item=&'a T>,
      R: FromIterator<T>
{
    xs.scan(T::zero(), |acc, &x| {
            *acc += x;
            Some(*acc)
        })
        .collect()
}


/// Perform cumsum and then integrate it using trapezoidal rule. The result is stored
/// in the given contiguous array.
///
/// Similar to `scipy.integrate.cumtrapz`
///
/// Requirements: `ys.len() >= 2` and `ret.len() = ys.len() - 1`.
pub fn cumtrapz<T, In, Out>(ys: In, dx: T, mut ret: Out)
where T: NumAssign + Copy + Div<f64, Output=T>,
      In: AsRef<[T]>,
      Out: AsMut<[T]>,
{
    let ys  = ys.as_ref();
    let ret = ret.as_mut();
    assert!(ys.len() >= 2);
    assert!(ret.len() == ys.len() - 1);

    let len = ys.len();
    let mut sum = T::zero();

    for i in 1 .. len {
        sum += ys[i-1] + ys[i];
        ret[i - 1] = sum * dx / 2.0f64;
    }
}


/// Perform cumsum and then integrate it using trapezoidal rule then return the result
/// in a `Vec`.
///
/// Similar to `scipy.integrate.cumtrapz`
///
/// Requirements: `ys.len() >= 2`.
pub fn cumtrapz2<T, In>(ys: In, dx: T, init: bool) -> Vec<T>
where T: NumAssign + Copy + Div<f64, Output=T>,
      In: AsRef<[T]>,
{
    let len = ys.as_ref().len();
    assert!(len >= 2);

    if init {
        let mut ret = vec![T::zero(); len];
        cumtrapz(ys, dx, &mut ret[1..]);
        ret
    } else {
        let mut ret = vec![T::zero(); len - 1];
        cumtrapz(ys, dx, &mut ret);
        ret
    }
}


/// Calculate the self-correlate function. The result is stored in the given contiguous array.
///
/// `y[i] = SUM(x[0 .. n-i] * x[i .. n])`
///
/// where `xs.len() == n` and `ys.len() == n`
pub fn auto_correlation<T, In, Out>(xs: In, mut ys: Out)
where T: NumAssign + iter::Sum + Div<f64, Output=T> + Copy,
      In: AsRef<[T]>,
      Out: AsMut<[T]>,
{
    let xs = xs.as_ref();
    let ys = ys.as_mut();

    let len = xs.len();
    assert!(len >= 2);

    let xs_mean = mean(xs);
    let xs = xs.iter()      // shift to 0
        .cloned()
        .map(|x| x - xs_mean)
        .collect::<Vec<_>>();

    assert!(ys.len() == xs.len());

    for i in 0 .. len {
        // `y[i] = SUM(x[0 .. n-i] * x[i .. n])`
        ys[i] = iter::zip(xs[.. len-i].iter(), xs[i ..].iter())
            .map(|(&x1, &x2)| x1 * x2)
            .sum::<T>();
    }
}


/// Calculate the self-correlate function and then return the result in a new `Vec`
///
/// `y[i] = SUM(x[0 .. n-i] * x[i .. n])`
///
/// where `xs.len() == n` and `ys.len() == n`
pub fn auto_correlation2<T, In>(xs: In) -> Vec<T>
where T: NumAssign + iter::Sum + Div<f64, Output=T> + Copy,
      In: AsRef<[T]>
{
    let len = xs.as_ref().len();
    assert!(len >= 2);

    let mut ret = vec![T::zero(); len];
    auto_correlation(xs, &mut ret);

    ret
}


/// Calculate the frequencies for generic FFT.
pub fn fft_freq_1d<R>(n: usize, dx: f64) -> R
where R: FromIterator<f64>
{
    let n = n as i64;

    let positive_end = n/2;
    let negative_beg = 1 + n/2 - n;

    (0 ..= positive_end).chain(negative_beg ..= -1)
        .map(|x| x as f64 / (dx * n as f64))
        .collect()
}


/// Calculate the frequencies for real-to-complex FFT.
pub fn rfft_freq_1d<R>(n: usize, dx: f64) -> R
where R: FromIterator<f64>
{
    let end = n / 2;
    (0 ..= end)
        .map(|x| x as f64 / (dx * n as f64))
        .collect()
}


#[cfg(test)]
mod tests{
    use super::*;

    #[test]
    fn test_cumsum() {
        let ys = &[1i32, 2, 3, 4, 5];
        let expect = &[1i32, 3, 6, 10, 15];
        let ret: Vec<i32> = cumsum(ys.iter());
        assert_eq!(ret, expect);
    }


    #[test]
    fn test_cumtrapz() {
        let ys = &[1.0f64, 2.0, 3.0, 4.0, 5.0];
        let dx = 0.3f64;
        let ret = cumtrapz2(ys, dx, false);
        assert_eq!(ret[2], 2.25);
        assert_eq!(ret[3], 3.5999999999999996);
    }


    #[test]
    fn test_auto_correlation() {
        let xs = &[1.0, 2.0, 3.0, 4.0, 5.0];
        let ys = auto_correlation2(xs);
        assert_eq!(ys[0], 10.0);
        assert_eq!(ys[1], 4.0);
        assert_eq!(ys[4], -4.0);
    }


    #[test]
    fn test_fft_freq() {
        let freq: Vec<f64> = fft_freq_1d(4, 0.5);
        assert_eq!(freq[1],  0.5);
        assert_eq!(freq[3], -0.5);

        let freq: Vec<f64> = fft_freq_1d(5, 0.5);
        assert_eq!(freq[1],  0.4);
        assert_eq!(freq[2],  0.8);
        assert_eq!(freq[3], -0.8);

        let freq: Vec<f64> = rfft_freq_1d(4, 0.5);
        assert_eq!(freq.len(), 3);
        assert_eq!(freq[2], 1.0);

        let freq: Vec<f64> = rfft_freq_1d(5, 0.5);
        assert_eq!(freq.len(), 3);
        assert_eq!(freq[2], 0.8);
    }
}
