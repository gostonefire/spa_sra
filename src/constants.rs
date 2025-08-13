pub const PI: f64 = std::f64::consts::PI;
pub const SUN_RADIUS: f64 = 0.26667;

pub const L_COUNT: usize = 6;
pub const B_COUNT: usize = 2;
pub const R_COUNT: usize = 5;
pub const Y_COUNT: usize = 63;

pub const L_MAX_SUBCOUNT: usize = 64;
pub const B_MAX_SUBCOUNT: usize = 5;
pub const R_MAX_SUBCOUNT: usize = 40;

pub const TERM_A: usize = 0;
pub const TERM_B: usize = 1;
pub const TERM_C: usize = 2;
pub const TERM_COUNT: usize = 3;

pub const TERM_X0: usize = 0;
pub const TERM_X1: usize = 1;
pub const TERM_X2: usize = 2;
pub const TERM_X3: usize = 3;
pub const TERM_X4: usize = 4;
pub const TERM_X_COUNT: usize = 5;

pub const TERM_PSI_A: usize = 0;
pub const TERM_PSI_B: usize = 1;
pub const TERM_EPS_C: usize = 2;
pub const TERM_EPS_D: usize = 3;
pub const TERM_PE_COUNT: usize = 4;

pub const JD_MINUS: usize = 0;
pub const JD_ZERO: usize = 1;
pub const JD_PLUS: usize = 2;
pub const JD_COUNT: usize = 3;

pub const SUN_TRANSIT: usize = 0;
pub const SUN_RISE: usize = 1;
pub const SUN_SET: usize = 2;
pub const SUN_COUNT: usize = 3;

pub const TERM_Y_COUNT: usize = TERM_X_COUNT;

pub const L_SUBCOUNT: [i64;L_COUNT] = [64,34,20,7,3,1];
pub const B_SUBCOUNT: [i64;B_COUNT] = [5,2];
pub const R_SUBCOUNT: [i64;R_COUNT] = [40,10,6,2,1];
