//! Utility functions for other applications (such as NREL's SAMPA).
//! They are all used internally by the main SPA algorithm, but was also made public in
//! the original C-code header file, hence they are also here denoted as utility functions.

use crate::constants::{PI, SUN_RADIUS};

/// Utility function from original C-code
pub fn rad2deg(radians: f64) -> f64 {
    180.0/PI * radians
}

/// Utility function from original C-code
pub fn deg2rad(degrees: f64) -> f64 {
    PI/180.0 * degrees
}

/// Utility function from original C-code
pub fn limit_degrees(mut degrees: f64) -> f64 {
    degrees /= 360.0;
    let limited = 360.0 * (degrees - degrees.floor());
    
    if limited < 0.0 {
        limited + 360.0
    } else {
        limited
    }
}

/// Utility function from original C-code
pub fn third_order_polynomial(a: f64, b: f64, c: f64, d: f64, x: f64) -> f64 {
    ((a * x + b) * x + c) * x + d
}

/// Utility function from original C-code
pub fn geocentric_right_ascension(lamda: f64, epsilon: f64, beta: f64) -> f64 {
    let lamda_rad: f64   = deg2rad(lamda);
    let epsilon_rad: f64 = deg2rad(epsilon);

    limit_degrees(rad2deg((lamda_rad.sin() * epsilon_rad.cos() -
        deg2rad(beta).tan() * epsilon_rad.sin()).atan2(lamda_rad.cos())))
}

/// Utility function from original C-code
pub fn geocentric_declination(beta: f64, epsilon: f64, lamda: f64) -> f64 {
    let beta_rad: f64    = deg2rad(beta);
    let epsilon_rad: f64 = deg2rad(epsilon);

    rad2deg((beta_rad.sin() * epsilon_rad.cos() +
        beta_rad.cos() * epsilon_rad.sin() * deg2rad(lamda).sin()).asin())
}

/// Utility function from original C-code
pub fn observer_hour_angle(nu: f64, longitude: f64, alpha_deg: f64) -> f64 {
    limit_degrees(nu + longitude - alpha_deg)
}

/// Utility function from original C-code
pub fn right_ascension_parallax_and_topocentric_dec(latitude: f64, elevation: f64, xi: f64, h: f64, delta: f64, delta_alpha: &mut f64, delta_prime: &mut f64) {
    let lat_rad: f64   = deg2rad(latitude);
    let xi_rad: f64    = deg2rad(xi);
    let h_rad: f64     = deg2rad(h);
    let delta_rad: f64 = deg2rad(delta);
    let u: f64 = (0.99664719 * lat_rad.tan()).atan();
    let y: f64 = 0.99664719 * u.sin() + elevation * lat_rad.sin() / 6378140.0;
    let x: f64 =              u.cos() + elevation * lat_rad.cos() / 6378140.0;

    let delta_alpha_rad: f64 = (             - x * xi_rad.sin()  * h_rad.sin())
        .atan2(        delta_rad.cos() - x * xi_rad.sin()  * h_rad.cos());

    *delta_prime = rad2deg(((delta_rad.sin() - y * xi_rad.sin()) * delta_alpha_rad.cos())
        .atan2(        delta_rad.cos() - x * xi_rad.sin()  * h_rad.cos()));

    *delta_alpha = rad2deg(delta_alpha_rad);
}

/// Utility function from original C-code
pub fn topocentric_right_ascension(alpha_deg: f64, delta_alpha: f64) -> f64 {
    alpha_deg + delta_alpha
}

/// Utility function from original C-code
pub fn topocentric_local_hour_angle(h: f64, delta_alpha: f64) -> f64 {
    h - delta_alpha
}

/// Utility function from original C-code
pub fn topocentric_elevation_angle(latitude: f64, delta_prime: f64, h_prime: f64) -> f64 {
    let lat_rad: f64         = deg2rad(latitude);
    let delta_prime_rad: f64 = deg2rad(delta_prime);

    rad2deg((lat_rad.sin() * delta_prime_rad.sin() + 
             lat_rad.cos() * delta_prime_rad.cos() * deg2rad(h_prime).cos()).asin())
}

/// Utility function from original C-code
pub fn atmospheric_refraction_correction(pressure: f64, temperature: f64, atmos_refract: f64, e0: f64) -> f64 {
    let mut del_e: f64 = 0.0;

    if e0 >= -1.0 * (SUN_RADIUS + atmos_refract) {
        del_e = (pressure / 1010.0) * (283.0 / (273.0 + temperature)) * 
            1.02 / (60.0 * deg2rad(e0 + 10.3/(e0 + 5.11)).tan());
    }

    del_e
}

/// Utility function from original C-code
pub fn topocentric_elevation_angle_corrected(e0: f64, delta_e: f64) -> f64 {
    e0 + delta_e
}

/// Utility function from original C-code
pub fn topocentric_zenith_angle(e: f64) -> f64 {
    90.0 - e
}

/// Utility function from original C-code
pub fn topocentric_azimuth_angle_astro(h_prime: f64, latitude: f64, delta_prime: f64) -> f64 {
    let h_prime_rad: f64 = deg2rad(h_prime);
    let lat_rad: f64     = deg2rad(latitude);

    limit_degrees(rad2deg(h_prime_rad.sin().atan2(h_prime_rad.cos() * lat_rad.sin() - deg2rad(delta_prime).tan() * lat_rad.cos())))
}

/// Utility function from original C-code
pub fn topocentric_azimuth_angle(azimuth_astro: f64) -> f64 {
    limit_degrees(azimuth_astro + 180.0)
}
