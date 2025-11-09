//! Defined error and its codes and messages

use std::fmt::{Display, Formatter};

pub(crate) const MESSAGES: [&str;18] = [
    "",
    "4-digit year outside valid range: -2000 to 6000",
    "2-digit month outside valid range: 1 to  12",
    "2-digit day outside valid range: 1 to  31",
    "Observer local hour outside valid range: 0 to  24",
    "Observer local minute outside valid range: 0 to  59",
    "Observer local second outside valid range: 0 to <60",
    "Delta t outside valid range: -8000 to 8000 seconds",
    "Observer time zone (negative west of Greenwich) outside valid range: -18 to 18 hours",
    "Observer longitude (negative west of Greenwich) outside valid range: -180 to 180 degrees",
    "Observer latitude (negative south of equator) outside valid range: -90 to 90 degrees",
    "Observer elevation [meters] outside valid range: -6500000 or higher meters",
    "Annual average local pressure [millibars] outside valid range: 0 to 5000 millibars",
    "Annual average local temperature [degrees Celsius] outside valid range: -273 to 6000 degrees Celsius",
    "Surface slope (measured from the horizontal plane) outside valid range: -360 to 360 degrees",
    "Surface azimuth rotation (measured from south to projection of surface normal on horizontal plane, negative east) outside valid range: -360 to 360 degrees",
    "Atmospheric refraction at sunrise and sunset (0.5667 deg is typical) outside valid range: -5 to 5 degrees",
    "Delta UT1 outside valid range: -1 to 1 second (exclusive)",
];

/// Spa related errors
/// * 1 -> 4-digit year outside valid range: -2000 to 6000
/// * 2 -> 2-digit month outside valid range: 1 to  12
/// * 3 -> 2-digit day outside valid range: 1 to  31
/// * 4 -> Observer local hour outside valid range: 0 to  24
/// * 5 -> Observer local minute outside valid range: 0 to  59
/// * 6 -> Observer local second outside valid range: 0 to <60
/// * 17 -> Delta UT1 outside valid range: -1 to 1 second (exclusive)
/// * 7 -> Delta t outside valid range: -8000 to 8000 seconds
/// * 8 -> Observer time zone (negative west of Greenwich) outside valid range: -18 to 18 hours
/// * 9 -> Observer longitude (negative west of Greenwich) outside valid range: -180 to 180 degrees
/// * 10 -> Observer latitude (negative south of equator) outside valid range: -90 to 90 degrees
/// * 11 -> Observer elevation \[meters\] outside valid range: -6500000 or higher meters
/// * 12 -> Annual average local pressure \[millibars\] outside valid range: 0 to 5000 millibars
/// * 13 -> Annual average local temperature \[degrees Celsius\] outside valid range: -273 to 6000 degrees Celsius
/// * 14 -> Surface slope (measured from the horizontal plane) outside valid range: -360 to 360 degrees
/// * 15 -> Surface azimuth rotation (measured from south to projection of surface normal on horizontal plane, negative east) outside valid range: -360 to 360 degrees
/// * 16 -> Atmospheric refraction at sunrise and sunset (0.5667 deg is typical) outside valid range: -5 to 5 degrees
#[derive(Debug)]
pub struct SpaError {
    /// Error code
    pub code: i64,
    /// Error message
    pub message: &'static str,
}

impl Display for SpaError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "SpaError: {} -> {}", self.code, self.message)
    }
}

impl std::error::Error for SpaError {}