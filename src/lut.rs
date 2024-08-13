use core::f64;
use std::{
    fs,
    io::{self, Read, Seek},
    mem::{self, MaybeUninit},
    slice,
};

use crate::*;

// pub const THRESHOLD: f64 = 1e-8;
// pub const VERIFY_THRESHOLD: f64 = 1e-4;
// pub const N_ITER: usize = 10;
// pub const THETA_STEP: usize = 3;
// pub const KAPPAS: &[f64] = &[
//     -1008.000, -875.444, -760.047, -659.588, -572.134, -496.000, -429.722, -372.023, -321.794,
//     -278.067, -240.000, -206.861, -178.012, -152.897, -131.033, -112.000, -95.430, -81.006,
//     -68.449, -57.517, -48.000, -39.715, -32.503, -26.224, -20.758, -16.000, -11.858, -8.251,
//     -6.627, -5.112, -3.698, -2.379, -1.500, -1.000, -0.500, 0.000, 0.500, 1.000, 1.500, 2.379,
//     3.698, 5.112, 6.627, 8.251, 11.858, 16.000, 20.758, 26.224, 32.503, 39.715, 48.000, 57.517,
//     68.449, 81.006, 95.430, 112.000, 131.033, 152.897, 178.012, 206.861, 240.000, 278.067, 321.794,
//     372.023, 429.722, 496.000, 572.134, 659.588, 760.047, 875.444, 1008.000,
// ];
// pub const ABS_MAX_KAPPA: f64 = 1160.;

pub const THRESHOLD: f64 = 1e-8;
pub const VERIFY_THRESHOLD: f64 = 1e-3;
pub const N_ITER: usize = 7;
pub const THETA_STEP: usize = 10;
pub const KAPPAS: &[f64] = &[
    -48.000, -39.715, -32.503, -26.224, -20.758, -16.000, -11.858, -8.251, -6.627, -5.112, -3.698,
    -2.379, -1.500, -1.000, -0.500, 0.000, 0.500, 1.000, 1.500, 2.379, 3.698, 5.112, 6.627, 8.251,
    11.858, 16.000, 20.758, 26.224, 32.503, 39.715, 48.000,
];
pub const ABS_MAX_KAPPA: f64 = 57.517;

pub const N_THETA: usize = 360 / THETA_STEP;
pub const N_KAPPA: usize = KAPPAS.len();
pub const N_ELEM: usize = (N_THETA * N_KAPPA).pow(2);

pub fn ij_from_theta(theta: f64) -> Option<usize> {
    let theta = theta.to_degrees();
    (theta.abs() <= 180.).then(|| (((theta + 180.) / THETA_STEP as f64).round() as usize) % N_THETA)
}

pub fn kl_from_kappa(kappa: f64) -> Option<usize> {
    if kappa.abs() > ABS_MAX_KAPPA {
        return None;
    }
    for (i, r) in KAPPAS.windows(2).map(|r| (r[0] + r[1]) / 2.).enumerate() {
        if kappa < r {
            return Some(i);
        }
    }
    Some(N_KAPPA - 1)
}

pub fn theta_from_ij(i: usize) -> f64 {
    ((i * THETA_STEP) as f64 - 180.).to_radians()
}

pub fn kappa_from_kl(k: usize) -> f64 {
    KAPPAS[k]
}

pub fn index_for_query(theta0: f64, theta1: f64, kappa0: f64, kappa1: f64) -> Option<[usize; 4]> {
    Some([
        ij_from_theta(theta0)?,
        ij_from_theta(theta1)?,
        kl_from_kappa(kappa0)?,
        kl_from_kappa(kappa1)?,
    ])
}

pub fn index_from_result(theta0: f64, theta1: f64, kappa0: f64, kappa1: f64) -> Option<[usize; 4]> {
    Some([
        ij_from_theta(theta0)?,
        ij_from_theta(theta1)?,
        kl_from_kappa(kappa0)?,
        kl_from_kappa(kappa1)?,
    ])
}

pub fn result_from_index([i, j, k, l]: [usize; 4]) -> [f64; 4] {
    [
        theta_from_ij(i),
        theta_from_ij(j),
        kappa_from_kl(k),
        kappa_from_kl(l),
    ]
}

pub fn index_of_hb(params: &hb::HyperbezParams) -> Option<[usize; 4]> {
    let arg_uv = params.integrate(1.).angle();
    let theta1 = params.theta(1.);
    let kappa0 = params.kappa(0.);
    let kappa1 = params.kappa(1.);
    (arg_uv.is_finite() && theta1.is_finite() && kappa0.is_finite() && kappa1.is_finite())
        .then(|| index_from_result(-arg_uv, -arg_uv + theta1, kappa0, kappa1))
        .flatten()
}

pub fn clamp_kl(k: usize) -> Option<usize> {
    (k < N_KAPPA).then_some(k)
}

pub fn next_indices_of([i, j, k, l]: [usize; 4]) -> impl Iterator<Item = [usize; 4]> {
    (0b0001..=0b1111).filter_map(move |f| {
        Some([
            (i + ((f & 0b1000) >> 3)) % N_THETA,
            (j + ((f & 0b0100) >> 2)) % N_THETA,
            clamp_kl(k + ((f & 0b0010) >> 1))?,
            clamp_kl(l + (f & 0b0001))?,
        ])
    })
}

pub fn prev_indices_of([i, j, k, l]: [usize; 4]) -> impl Iterator<Item = [usize; 4]> {
    (0b0001..=0b1111).rev().filter_map(move |f| {
        Some([
            i.checked_sub((f & 0b1000) >> 3).unwrap_or(N_THETA - 1),
            j.checked_sub((f & 0b0100) >> 2).unwrap_or(N_THETA - 1),
            k.checked_sub((f & 0b0010) >> 1)?,
            l.checked_sub(f & 0b0001)?,
        ])
    })
}

pub fn roll_from_index([i, j, k, l]: [usize; 4]) -> usize {
    i * N_THETA * N_KAPPA * N_KAPPA + j * N_KAPPA * N_KAPPA + k * N_KAPPA + l
}

pub fn index_from_roll(roll: usize) -> [usize; 4] {
    [
        roll / N_THETA / N_KAPPA / N_KAPPA,
        roll / N_KAPPA / N_KAPPA % N_THETA,
        roll / N_KAPPA % N_KAPPA,
        roll % N_KAPPA,
    ]
}

pub fn solve_from_index_init(
    index: [usize; 4],
    init: [f64; 3],
) -> Result<([f64; 4], f64, usize), hb::SolveError> {
    let [theta0, theta1, kappa0, kappa1] = result_from_index(index);
    hb::solve_for_params_exact(theta0, theta1, kappa0, kappa1, init, THRESHOLD, N_ITER)
}

pub fn verify_for_index(index: [usize; 4], estimate: [f64; 4]) -> bool {
    let [theta0, theta1, kappa0, kappa1] = result_from_index(index);
    let [a, b, c, d] = estimate;
    let hb_params = hb::HyperbezParams::new(a, b, c, d);
    let arg_uv = -hb_params.integrate(1.).angle();

    [
        theta0 - arg_uv,
        theta1 - arg_uv - hb_params.theta(1.),
        kappa0 - hb_params.kappa(0.),
        kappa1 - hb_params.kappa(1.),
    ]
    .into_iter()
    .all(|delta| delta.abs() < VERIFY_THRESHOLD)
}

#[allow(unused)]
pub fn range_step_inclusive(start: f64, stop: f64, step: f64) -> impl Iterator<Item = f64> {
    (0..)
        .map(move |i| start + i as f64 * step)
        .take_while(move |i| *i <= stop)
}

#[allow(unused)]
pub fn range_step(start: f64, stop: f64, step: f64) -> impl Iterator<Item = f64> {
    (0..)
        .map(move |i| start + i as f64 * step)
        .take_while(move |i| *i < stop)
}

pub struct Lut {
    estimate: Box<[[f64; 4]]>,
    variance: Box<[f64]>,
}

impl Lut {
    pub fn new() -> Self {
        Lut {
            estimate: vec![[f64::NAN; 4]; N_ELEM].into_boxed_slice(),
            variance: vec![f64::INFINITY; N_ELEM].into_boxed_slice(),
        }
    }

    pub fn read_from_files(iter: usize) -> io::Result<Self> {
        let mut estimate = Box::new_uninit_slice(N_ELEM);
        let mut variance = Box::new_uninit_slice(N_ELEM);

        let mut f = fs::File::open(format!("lut_estimate_{iter}.bin"))?;
        f.read_exact(unsafe {
            MaybeUninit::slice_assume_init_mut(MaybeUninit::slice_as_bytes_mut(&mut estimate))
        })?;

        f = fs::File::open(format!("lut_variance_{iter}.bin"))?;
        f.read_exact(unsafe {
            MaybeUninit::slice_assume_init_mut(MaybeUninit::slice_as_bytes_mut(&mut variance))
        })?;

        Ok(Lut {
            estimate: unsafe { estimate.assume_init() },
            variance: unsafe { variance.assume_init() },
        })
    }

    pub fn write_to_files(&self, iter: usize) -> io::Result<()> {
        fs::write(format!("lut_estimate_{iter}.bin"), unsafe {
            slice::from_raw_parts(
                self.estimate.as_ptr() as *const u8,
                mem::size_of_val(self.estimate.as_ref()),
            )
        })?;
        fs::write(format!("lut_variance_{iter}.bin"), unsafe {
            slice::from_raw_parts(
                self.variance.as_ptr() as *const u8,
                mem::size_of_val(self.variance.as_ref()),
            )
        })
    }

    pub fn read_from_f32_file() -> io::Result<Self> {
        let buf = fs::read("lut_estimate_f32.bin")?;
        Ok(Lut::read_from_f32(&buf))
    }

    pub fn read_from_f32(buf: &[u8]) -> Self {
        let estimate = buf
            .chunks_exact(16)
            .map(|buf| {
                std::array::from_fn(|i| {
                    f32::from_be_bytes(buf[i * 2..i * 2 + 4].try_into().unwrap()) as f64
                })
            })
            .collect::<Box<_>>();

        Lut {
            estimate,
            variance: Box::new([]),
        }
    }

    pub fn write_f32_estimate_to_file(&self) -> io::Result<()> {
        let bytes = self
            .estimate
            .iter()
            .flat_map(|e| e.map(|x| x as f32).map(f32::to_be_bytes))
            .flatten()
            .collect::<Vec<_>>();
        fs::write("lut_estimate_f32.bin", bytes)
    }

    pub fn read_from_f16_file() -> io::Result<Self> {
        let buf = fs::read("lut_estimate_f16.bin")?;
        Ok(Lut::read_from_f16(&buf))
    }

    pub fn read_from_f16(buf: &[u8]) -> Self {
        let estimate = buf
            .chunks_exact(8)
            .map(|buf| {
                std::array::from_fn(|i| {
                    half::f16::from_be_bytes([buf[i * 2], buf[i * 2 + 1]]).to_f64()
                })
            })
            .collect::<Box<_>>();

        Lut {
            estimate,
            variance: Box::new([]),
        }
    }

    pub fn write_f16_estimate_to_file(&self) -> io::Result<()> {
        let bytes = self
            .estimate
            .iter()
            .flat_map(|e| e.map(half::f16::from_f64).map(half::f16::to_be_bytes))
            .flatten()
            .collect::<Vec<_>>();
        fs::write("lut_estimate_f16.bin", bytes)
    }

    fn populate(&mut self, index: [usize; 4], init: [f64; 3]) {
        let roll = roll_from_index(index);
        let existing_var = &mut self.variance[roll];
        if *existing_var <= THRESHOLD {
            return;
        }
        let Ok((solved, e, ..)) = solve_from_index_init(index, init) else {
            return;
        };
        if !existing_var.is_finite() && verify_for_index(index, solved) {
            self.estimate[roll] = solved;
            *existing_var = e;
        }
    }

    fn populate_adj(&mut self, index: [usize; 4], init: [f64; 3]) {
        for adj_index in prev_indices_of(index).chain(next_indices_of(index)) {
            self.populate(adj_index, init)
            // let opp_index = array::from_fn(|i| 2 * index[i] - adj_index[i]);
            // let opp_roll = roll_from_index(opp_index);
            // if self.variance[opp_roll].is_finite() {
            //     let opp_init = self.estimate[opp_roll];
            //     let res = result_from_index(index);
            //     let opp_res = result_from_index(opp_index);
            //     let adj_res = result_from_index(adj_index);

            //     let adj_init = todo!();
            //     self.populate(adj_index, adj_init)
            // }
        }
    }

    pub fn init(&mut self) {
        // let mut kappa0_min = f64::INFINITY;
        // let mut kappa0_max = f64::NEG_INFINITY;
        // let mut kappa1_min = f64::INFINITY;
        // let mut kappa1_max = f64::NEG_INFINITY;
        for a in range_step_inclusive(-20., 20., 1.) {
            for b in range_step_inclusive(-20., 20., 1.) {
                for c in range_step_inclusive(1., 20., 1.) {
                    let d_limit = hb::d_limit_rounded(a, b, c);
                    for d in range_step_inclusive(d_limit.start, d_limit.end, 0.1) {
                        let params = hb::HyperbezParams::new(a, b, c, d);
                        // let kappa0 = params.kappa(0.);
                        // let kappa1 = params.kappa(1.);
                        // kappa0_min = kappa0_min.min(kappa0);
                        // kappa0_max = kappa0_max.max(kappa0);
                        // kappa1_min = kappa1_min.min(kappa1);
                        // kappa1_max = kappa1_max.max(kappa1);

                        let init = [a, c, d];

                        let Some(index) = index_of_hb(&params) else {
                            continue;
                        };
                        // let abs_max_theta = hb::abs_max_theta(&params);
                        // if abs_max_theta.abs() >= f64::consts::TAU
                        // // || (abs_max_theta - params.theta(1.)).abs() >= f64::consts::TAU
                        // {
                        //     continue;
                        // }

                        self.populate(index, init);
                        self.populate_adj(index, init);
                        // print!("L");
                    }
                    // print!("K");
                }
                // print!("J");
            }
            // print!("I");
        }
        // dbg!(kappa0_min, kappa0_max, kappa1_min, kappa1_max);
    }

    pub fn iterate(&mut self) {
        for roll in 0..N_ELEM {
            // skip if roll is not filled
            if self.variance[roll].is_infinite() {
                continue;
            }

            let [a, _, c, d] = self.estimate[roll];
            self.populate_adj(index_from_roll(roll), [a, c, d]);
        }
    }

    pub fn verify(&self) -> [usize; 3] {
        let mut pass = 0;
        let mut fail = 0;
        let mut skip = 0;
        for roll in 0..N_ELEM {
            // skip if roll is not filled
            if self.variance[roll].is_infinite() {
                skip += 1;
                continue;
            }

            if verify_for_index(index_from_roll(roll), self.estimate[roll]) {
                pass += 1;
            } else {
                fail += 1;
            }
        }
        [pass, fail, skip]
    }

    pub fn print_status(&self) {
        let total = N_ELEM;
        let filled = self.variance.iter().filter(|v| v.is_finite()).count();
        // let unacceptable = self
        //     .variance
        //     .iter()
        //     .filter(|v| v.is_finite() && **v > THRESHOLD)
        //     .count();
        let max = self
            .variance
            .iter()
            .copied()
            .max_by(f64::total_cmp)
            .unwrap();
        let min = self
            .variance
            .iter()
            .copied()
            .min_by(f64::total_cmp)
            .unwrap();
        let empty = total - filled;
        dbg!(total, filled, /*unacceptable,*/ empty, max, min);
    }

    pub fn solve_for_params(
        &self,
        theta0: f64,
        theta1: f64,
        kappa0: f64,
        kappa1: f64,
        threshold: f64,
        n_iter: usize,
    ) -> Option<[Option<hb::HyperbezParams>; 4]> {
        let hb_from_params =
            |([a, b, c, d], ..): ([f64; 4], f64, usize)| hb::HyperbezParams::new(a, b, c, d);

        let mut params = [
            self.solve_for_params_exact(theta0, theta1, -kappa0, -kappa1, threshold, n_iter)?,
            self.solve_for_params_exact(theta0, theta1, -kappa0, kappa1, threshold, n_iter)?,
            self.solve_for_params_exact(theta0, theta1, kappa0, -kappa1, threshold, n_iter)?,
            self.solve_for_params_exact(theta0, theta1, kappa0, kappa1, threshold, n_iter)?,
        ]
        .map(|p| p.map(hb_from_params).ok());

        // TODO: figure out if these should be left as results

        let f_base = |p: &hb::HyperbezParams| p.integrate(1.).length_squared();
        params.sort_by(|a, b| {
            a.as_ref()
                .map(f_base)
                .zip(b.as_ref().map(f_base))
                .map_or(std::cmp::Ordering::Equal, |(a, b)| a.total_cmp(&b))
        });

        Some(params)
    }

    pub fn solve_for_params_exact(
        &self,
        theta0: f64,
        theta1: f64,
        kappa0: f64,
        kappa1: f64,
        threshold: f64,
        n_iter: usize,
    ) -> Option<Result<([f64; 4], f64, usize), hb::SolveError>> {
        let index = index_for_query(theta0, theta1, kappa0, kappa1)?;
        let [a, _, c, d] = self.estimate[roll_from_index(index)];
        // dbg!(index, init);
        Some(hb::solve_for_params_exact(
            theta0,
            theta1,
            kappa0,
            kappa1,
            [a, c, d],
            threshold,
            n_iter,
        ))
    }
}

impl Default for Lut {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(not(target_arch = "wasm32"))]
fn params_from_buffer(index: [usize; 4]) -> [f64; 4] {
    let roll = roll_from_index(index);
    let mut buf = [0u8; 8];
    let mut estimate_file = std::fs::File::open("lut_estimate_f16.bin").unwrap();
    estimate_file.seek(io::SeekFrom::Start(roll as _)).unwrap();
    estimate_file.read_exact(&mut buf).unwrap();
    std::array::from_fn(|i| half::f16::from_le_bytes([buf[i * 2], buf[i * 2 + 1]]).to_f64())
}

#[cfg(not(target_arch = "wasm32"))]
pub fn solve_for_params(
    theta0: f64,
    theta1: f64,
    kappa0: f64,
    kappa1: f64,
    threshold: f64,
    n_iter: usize,
) -> Option<[Option<hb::HyperbezParams>; 2]> {
    let index0 = index_for_query(theta0, theta1, -kappa0, -kappa1)?;
    let init0 = params_from_buffer(index0);
    dbg!(index0, init0);
    let params0 = hb::solve_for_params_exact(
        theta0,
        theta1,
        -kappa0,
        -kappa1,
        [init0[0], init0[2], init0[3]],
        threshold,
        n_iter,
    );

    let index1 = index_for_query(theta0, theta1, kappa0, -kappa1)?;
    let init1 = params_from_buffer(index1);
    dbg!(index1, init1);
    let params1 = hb::solve_for_params_exact(
        theta0,
        theta1,
        kappa0,
        -kappa1,
        [init1[0], init1[2], init1[3]],
        threshold,
        n_iter,
    );

    let hb_from_params =
        |([a, b, c, d], ..): ([f64; 4], f64, usize)| hb::HyperbezParams::new(a, b, c, d);

    // TODO: fihue out if these should be left as results
    let params0 = params0.map(hb_from_params).ok();
    let params1 = params1.map(hb_from_params).ok();

    let f_base = |p: &hb::HyperbezParams| p.integrate(1.).length_squared();

    Some(
        if params0.as_ref().map(f_base) >= params1.as_ref().map(f_base) {
            [params0, params1]
        } else {
            [params1, params0]
        },
    )
}

#[test]
fn test_theta_mappings() {
    for t in -400..400 {
        eprintln!(
            "{t} => {:?} => {:?}",
            ij_from_theta((t as f64).to_radians()),
            ij_from_theta((t as f64).to_radians())
        );
    }
}

#[test]
fn test_kappa_mappings() {
    for t in range_step(-1200., -500., 10.)
        .chain(range_step(-500., -100., 5.))
        .chain(range_step(-100., -20., 2.))
        .chain(range_step(-20., -4., 0.5))
        .chain(range_step(-4., 4., 0.1))
        .chain(range_step(4., 20., 0.5))
        .chain(range_step(20., 100., 2.))
        .chain(range_step(100., 500., 5.))
        .chain(range_step_inclusive(500., 1200., 10.))
    {
        eprintln!("{t:>7.1} => {:>2?}", kl_from_kappa(t));
    }
}

#[test]
fn test_half() {
    let i = 6.514942734853722e-7f64;
    let j = half::f16::from_f64(i);
    let k = j.to_be_bytes();
    let l = half::f16::from_be_bytes([k[0], k[1]]);
    dbg!(i, j, k, l);
}

#[test]
fn test_half2() {
    let i = [[2f64, 3., 4., 7.], [3., 5., 2., 7.]];
    let j = i
        .iter()
        .flat_map(|e| e.map(half::f16::from_f64).map(half::f16::to_be_bytes))
        .flatten()
        .collect::<Box<_>>();
    let k = j
        .chunks_exact(8)
        .map(|buf| std::array::from_fn(|i| half::f16::from_be_bytes([buf[i * 2], buf[i * 2 + 1]])))
        .collect::<Box<[[_; 4]]>>();
    let l = k
        .iter()
        .map(|buf| buf.map(half::f16::to_f64))
        .collect::<Box<[[_; 4]]>>();
    dbg!(i, j, k, l);
}

#[test]
fn test_solve_0_0_0_0_lutey() {
    let theta0 = (30.566393078928684f64).to_radians();
    let theta1 = (-26.72938643415364f64).to_radians();
    let kappa0 = 1.;
    let kappa1 = 0.7071067811865475;
    let threshold = 1e-4;
    let n_iter = 5;

    let lut = Lut::read_from_files(19).unwrap();
    let params = lut.solve_for_params(theta0, theta1, kappa0, kappa1, threshold, n_iter);
    dbg!(&params);

    let lut = Lut::read_from_f16_file().unwrap();
    let params = lut.solve_for_params(theta0, theta1, kappa0, kappa1, threshold, n_iter);
    dbg!(&params);
}

#[test]
fn test_make_0_0_0_0() {
    let params = hb::HyperbezParams::new(
        0.8574765235303419,
        -1.1866876695084798,
        -0.315043950143241,
        -0.15689405666050837,
    );
    dbg!(hb::DebugHbParams(&params));
}

#[test]
fn test_make_0_0_0_1() {
    let params = hb::HyperbezParams::new(
        0.9083586923675648,
        -1.2221115709196142,
        -0.2956153701871316,
        -0.16665328299490925,
    );
    dbg!(hb::DebugHbParams(&params));
}

// #[test]
// fn test_solve_0_0_0_0() {
//     let params = solve_for_params(
//         30.6f64.to_radians(),
//         (-26.7f64).to_radians(),
//         1.,
//         0.71,
//         1e-4,
//         5,
//     );
//     dbg!(&params);
// }

#[test]
fn test_solve_35_25_7_7() {
    let [a, b, c, d] = [-20., -20., 1., -1.5];
    let param = hb::HyperbezParams::new(a, b, c, d);
    let index = index_of_hb(&param).unwrap();
    let [theta0, theta1, kappa0, kappa1] = result_from_index(index);
    let theta0 = theta0.to_degrees();
    let theta1 = theta1.to_degrees();
    let (solved, err, ..) = solve_from_index_init(index, [a, c, d]).unwrap();
    let solved_param = hb::HyperbezParams::new(solved[0], solved[1], solved[2], solved[3]);
    let (solved_exact, err_exact, ..) = hb::solve_for_params_exact(
        (30.566393078928684f64).to_radians(),
        (-26.72938643415364f64).to_radians(),
        -1.0,
        -0.7071067811865475,
        [solved[0], solved[2], solved[3]],
        1e-8,
        10,
    )
    .unwrap();
    let solved_exact_param =
        hb::HyperbezParams::new(solved_exact[0], b, solved_exact[1], solved_exact[2]);
    dbg!(
        hb::DebugHbParams(&param),
        index,
        theta0,
        theta1,
        kappa0,
        kappa1,
        solved,
        err,
        hb::DebugHbParams(&solved_param),
        solved_exact,
        err_exact,
        hb::DebugHbParams(&solved_exact_param),
    );
}

#[test]
fn test_iter_0_0_0_0() {
    let [a, b, c, d] = [-20., -20., 1., -1.5];
    let params = hb::HyperbezParams::new(a, b, c, d);
    dbg!(hb::DebugHbParams(&params));
    for index in next_indices_of(index_of_hb(&params).unwrap()) {
        let [theta0, theta1, kappa0, kappa1] = result_from_index(index);
        dbg!((index, theta0, theta1, kappa0, kappa1));

        let Ok(solved) = solve_from_index_init(index, [a, c, d]) else {
            continue;
        };
        let roll = roll_from_index(index);
        dbg!(roll, solved);
    }
}
