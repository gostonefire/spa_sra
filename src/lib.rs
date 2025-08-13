#![warn(missing_docs)]
//! # SPA_SRA
//!
//! The SPA_SRA (Solar Position Algorithm for Solar Radiation Applications) calculates the solar zenith and
//! azimuth angles in the period from the year -2000 to 6000, with uncertainties of +/- 0.0003 degrees
//! based on the date, time, and location on Earth. (Reference: Reda, I.; Andreas, A., Solar Position Algorithm
//! for Solar Radiation Applications, Solar Energy. Vol. 76(5), 2004; pp. 577-589).
//!
//! It can also calculate the surface incidence angle for e.g. a solar panel.
//! The surface incidence angle is the angle between an incoming ray (like light or radar) and
//! a line perpendicular to the surface at the point where the ray hits.
//!
//! Further information on this algorithm is available in the following NREL technical report (pdf):
//! [Reda, I.; Andreas, A. (2003). Solar Position Algorithm for Solar Radiation Applications. 55 pp.; NREL Report No. TP-560-34302, Revised January 2008.]
//!
//! ## License and Acknowledgement
//!
//! SPA_SRA is a port with rust specific adjustments from the original C-code prepared by employees
//! of the Alliance for Sustainable Energy, LLC. Hence, the [MIT license] comes with an acknowledgement
//! and disclaimer related to that original code. For personal use it should be fine though.
//!
//! ## Usage overview
//!
//! The library has two main entry points:
//! * Using the [SpaBuilder] (which has an optional flavor if using the library with the chrono_0_4 feature)
//! * Using the [SpaData] struct and its sub-structs directly, which may be favourable for a more fine-grained integration.
//!
//! If choosing to use the chrono_0_4 feature then date, time and timezone is managed through the [chrono] crate.
//! If however the library is used without the chrono_0_4 feature, there are no external crates being
//! imported, and date, time and timezone are entered as simple numbers
//!
//! It is possible to use a hybrid between the [SpaBuilder] and the [SpaData] struct.
//!
//! However, if doing a hybrid when also using the chrono_0_4 feature it gets slightly more complicated.
//! Changing date and time directly in the [Input] struct doesn't automatically update the timezone
//! offset numeric field (which is used by the SPA-SRA algorithm) and can change between two hours
//! due to sunlight saving times.
//!
//! ## Example using the SpaData struct
//!
//! The code below uses the SpaData struct directly, instead of using the builder, to sweep
//! over one day and report a vector of incidence together with sunrise and sunset.
//! ```rust
//! use std::ops::Add;
//! use chrono::{DateTime, Local, TimeDelta, TimeZone};
//! use spa_sra::errors::SpaError;
//! use spa_sra::spa::{Function, Input, SpaData};
//!
//! fn main() {
//!
//!     match get_day_incidence() {
//!         Err(e) => { println!("{}", e) },
//!         Ok(incidence) => {
//!             println!("Sunrise: {}", incidence.0);
//!             println!("Sunset:  {}", incidence.1);
//!             println!("Incidence: {:?}", incidence.2);
//!         }
//!     }
//! }
//!
//! fn get_day_incidence() -> Result<(DateTime<Local>, DateTime<Local>, Vec<f64>), SpaError> {
//!     let mut result: Vec<f64> = Vec::new();
//!
//!     // First decide on a date, time is not relevant at this point
//!     let date_time = Local::now()
//!         .timezone()
//!         .with_ymd_and_hms(2025, 8, 12, 0, 0, 0)
//!         .unwrap();
//!
//!     // Create an Input instance with relevant parameters, and we are happy with the
//!     // defaults for atmospheric_refraction, delta_ut1 and delta_t.
//!     // We only need sunrise and sunset first, so we chose a specific function for that
//!     let mut input = Input::from_date_time(date_time);
//!     input.latitude = 56.198569;
//!     input.longitude = 15.637826;
//!     input.pressure = 1013.0;
//!     input.temperature = 10.0;
//!     input.elevation = 10.0;
//!     input.slope = 25.0;
//!     input.azm_rotation = 0.0;
//!     input.function = Function::SpaZaRts;
//!
//!     // Create and calculate a SpaData struct
//!     let mut spa = SpaData::new(input);
//!     spa.spa_calculate()?;
//!
//!     // Now we got sunrise and sunset for the chosen day
//!     let sunrise = spa.get_sunrise();
//!     let sunset = spa.get_sunset();
//!     let mut time_of_interest = sunrise;
//!
//!     // We are only interested in the incidence from now on so we set that function
//!     spa.input.function = Function::SpaZaInc;
//!
//!     // Loop through the day with a one minute incrementation and save the incidence to result
//!     while time_of_interest < sunset {
//!         spa.input.date_time(time_of_interest);
//!         spa.spa_calculate()?;
//!         result.push(spa.spa_za_inc.incidence);
//!
//!         time_of_interest = time_of_interest.add(TimeDelta::minutes(1));
//!     }
//!
//!     Ok((sunrise, sunset, result))
//! }
//! ```
//!
//! [MIT license]: https://github.com/gostonefire/spa_sra/blob/main/LICENSE
//! [Reda, I.; Andreas, A. (2003). Solar Position Algorithm for Solar Radiation Applications. 55 pp.; NREL Report No. TP-560-34302, Revised January 2008.]: http://www.nrel.gov/docs/fy08osti/34302.pdf
//! [chrono]: https://docs.rs/chrono/latest/chrono/

#[cfg(feature = "chrono_0_4")]
use chrono::{DateTime, Datelike, Offset, TimeZone, Timelike};
use crate::errors::{SpaError, MESSAGES};
use crate::spa::{Function, Input, SpaData, Std};

mod earth_periodic_terms;
mod constants;
mod nutation_obliquity_periodic_terms;
pub mod utils;
pub mod spa;
pub mod errors;

/// Builder for creating and calculating an operational SpaData struct
///
/// It comes with a couple of possible entrypoint to choose between:
/// * `SpaBuilder::new` to create one from scratch and not using the chrono_0_4 feature
/// * `SpaBuilder::from_input` to create one from an already created [SpaData] struct
/// * `SpaBuilder::from_date_time` to create a new timezone aware builder
///
pub struct SpaBuilder<T: Clone> {
    input: Input<T>,
}

/// Implementation of Std specific methods
///
impl SpaBuilder<Std> {
    /// Creates a new SpaBuilder.
    ///
    /// These fields will have quite valid defaults, at least within the neighborhood of 2025-07-10, those are:
    /// * [`SpaBuilder::atmospheric_refraction`] / [`Input::atmos_refract`]: 0.5667 (which is a typical value)
    /// * [`SpaBuilder::delta_ut1`] / [`Input::delta_ut1`]: 0.1 (doesn't change that often)
    /// * [`SpaBuilder::delta_t`] / [`Input::delta_t`]: 69.084 (doesn't change that often)
    ///
    /// # Example
    ///
    /// Usage which gives the same results as in the original C program test code:
    /// ```rust
    /// use spa_sra::{SpaBuilder};
    /// use spa_sra::errors::SpaError;
    /// use spa_sra::spa::{Function, SpaData, Std};
    ///
    /// fn main() {
    ///
    ///     match get_spa() {
    ///         Err(e) => { println!("{}", e) },
    ///         Ok(spa) => {
    ///             println!("Julian Day:    {:.6}", spa.spa_za.jd);
    ///             println!("L:             {:.6e} degrees", spa.spa_za.l);
    ///             println!("B:             {:.6e} degrees", spa.spa_za.b);
    ///             println!("R:             {:.6} AU", spa.spa_za.r);
    ///             println!("H:             {:.6} degrees", spa.spa_za.h);
    ///             println!("Delta Psi:     {:.6e} degrees", spa.spa_za.del_psi);
    ///             println!("Delta Epsilon: {:.6e} degrees", spa.spa_za.del_epsilon);
    ///             println!("Epsilon:       {:.6} degrees", spa.spa_za.epsilon);
    ///
    ///             println!("Zenith:        {:.6} degrees", spa.spa_za.zenith);
    ///             println!("Azimuth:       {:.6} degrees", spa.spa_za.azimuth);
    ///             println!("Incidence:     {:.6} degrees", spa.spa_za_inc.incidence);
    ///
    ///             let mut min: f64 = 60.0 * (spa.spa_za_rts.sunrise - (spa.spa_za_rts.sunrise as i64) as f64);
    ///             let mut sec: f64 = 60.0 * (min - (min as i64) as f64);
    ///             println!("Sunrise:       {:0>2}:{:0>2}:{:0>2} Local Time", spa.spa_za_rts.sunrise as i64, min as i64, sec as i64);
    ///
    ///             min = 60.0 * (spa.spa_za_rts.sunset - (spa.spa_za_rts.sunset as i64) as f64);
    ///             sec = 60.0 * (min - (min as i64) as f64);
    ///             println!("Sunset:        {:0>2}:{:0>2}:{:0>2} Local Time", spa.spa_za_rts.sunset as i64, min as i64, sec as i64);
    ///         }
    ///     }
    /// }
    ///
    /// fn get_spa() -> Result<SpaData<Std>, SpaError> {
    ///
    ///     let spa = SpaBuilder::new()
    ///         .date(2003, 10, 17)?
    ///         .time(12, 30, 30.0)?
    ///         .timezone(-7.0)?
    ///         .lat_long(39.742476, -105.1786)?
    ///         .pressure(820.0)?
    ///         .temperature(11.0)?
    ///         .atmospheric_refraction(0.5667)?
    ///         .elevation(1830.14)?
    ///         .slope(30.0)?
    ///         .azimuth_rotation(-10.0)?
    ///         .delta_ut1(0.0)?
    ///         .delta_t(67.0)?
    ///         .execute(Function::SpaAll)?;
    ///
    ///     Ok(spa)
    /// }
    /// ```
    ///
    pub fn new() -> Self {
        Self {
            input: Input::new_std(),
        }
    }

    /// Sets date
    /// 
    /// # Arguments
    /// 
    /// * 'year' - 4-digit year, valid range: -2000 to 6000
    /// * 'month' - 2-digit month, valid range: 1 to  12
    /// * 'day' - 2-digit day, valid range: 1 to  31
    pub fn date(mut self, year: i64, month: i64, day: i64) -> Result<Self, SpaError> {
        if year        < -2000   || year        > 6000   { return Err(SpaError{ code: 1, message: MESSAGES[1] }) };
        if month       < 1       || month       > 12     { return Err(SpaError{ code: 2, message: MESSAGES[2] }) };
        if day         < 1       || day         > 31     { return Err(SpaError{ code: 3, message: MESSAGES[3] }) };

        self.input.year = year;
        self.input.month = month;
        self.input.day = day;

        Ok(self)
    }

    /// Sets time
    ///
    /// # Arguments
    ///
    /// * 'hour' - Observer local hour, valid range: 0 to  24
    /// * 'minute' - Observer local minute, valid range: 0 to  59
    /// * 'second' - Observer local second, valid range: 0 to <60 (accepts fraction)
    pub fn time(mut self, hour: i64, minute: i64, second: f64) -> Result<Self, SpaError> {
        if hour        < 0       || hour        > 24     { return Err(SpaError{ code: 4, message: MESSAGES[4] }) };
        if minute      < 0       || minute      > 59     { return Err(SpaError{ code: 5, message: MESSAGES[5] }) };
        if second      < 0.0     || second      >=60.0   { return Err(SpaError{ code: 6, message: MESSAGES[6] }) };
        if hour        == 24     && minute      > 0      { return Err(SpaError{ code: 5, message: MESSAGES[5] }) };
        if hour        == 24     && second      > 0.0    { return Err(SpaError{ code: 6, message: MESSAGES[6] }) };

        self.input.hour = hour;
        self.input.minute = minute;
        self.input.second = second;

        Ok(self)
    }

    /// Sets observer time zone (negative west of Greenwich)
    ///
    /// # Arguments
    /// 
    /// * 'timezone' - valid range: -18 to 18 hours
    pub fn timezone(mut self, timezone: f64) -> Result<Self, SpaError> {
        if timezone.abs()      > 18.0       { return Err(SpaError{ code: 8, message: MESSAGES[8] }) };

        self.input.timezone = timezone;
        
        Ok(self)
    }
}

/// Implementation of chrono_0_4 specific methods
///
#[cfg(feature = "chrono_0_4")]
impl <T: TimeZone>SpaBuilder<T> {
    /// Creates a new timezone aware SpaBuilder from the given [DateTime]
    ///
    /// These fields will have quite valid defaults, at least within the neighborhood of 2025-07-10, those are:
    /// * [`SpaBuilder::atmospheric_refraction`] / [`Input::atmos_refract`]: 0.5667 (which is a typical value)
    /// * [`SpaBuilder::delta_ut1`] / [`Input::delta_ut1`]: 0.1 (doesn't change that often)
    /// * [`SpaBuilder::delta_t`] / [`Input::delta_t`]: 69.084 (doesn't change that often)
    ///
    /// This method is dependent on the feature "chrono_0_4" which will include the [chrono] crate.
    ///
    /// # Arguments
    ///
    /// * 'date_time' - a [DateTime] object including the time zone
    ///
    /// # Example
    ///
    /// Accepts defaults for atmospheric_refraction, delta_ut1 and delta_t.
    /// Parameters corresponds roughly for an existing solar park in Karlskrona, Sweden.
    /// Datetime is whatever the local time is at time of execution.
    /// ```rust
    /// use chrono::Local;
    /// use spa_sra::{SpaBuilder};
    /// use spa_sra::errors::SpaError;
    /// use spa_sra::spa::{Function, SpaData};
    ///
    ///
    /// fn main() {
    ///
    ///     match get_spa() {
    ///         Err(e) => { println!("{}", e) },
    ///         Ok(spa) => {
    ///             println!("Altitude:                 {:.6} degrees", spa.spa_za.e);
    ///             println!("Zenith:                   {:.6} degrees", spa.spa_za.zenith);
    ///             println!("Declination (geocentric): {:.6} degrees", spa.spa_za.delta);
    ///
    ///             println!("Azimuth:                  {:.6} degrees", spa.spa_za.azimuth);
    ///             println!("Incidence:                {:.6} degrees", spa.spa_za_inc.incidence);
    ///
    ///             println!("Sunrise:                  {} (Local Time)", spa.get_sunrise());
    ///             println!("Sunset:                   {} (Local Time)", spa.get_sunset());
    ///             println!("Sun transit:              {} (Local time)", spa.get_suntransit());
    ///         }
    ///     }
    /// }
    ///
    /// fn get_spa() -> Result<SpaData<Local>, SpaError> {
    ///
    ///     let spa = SpaBuilder::from_date_time(Local::now())
    ///         .lat_long(56.198569, 15.637826)?
    ///         .pressure(1013.0)?
    ///         .temperature(10.0)?
    ///         .elevation(10.0)?
    ///         .slope(25.0)?
    ///         .azimuth_rotation(0.0)?
    ///         .execute(Function::SpaAll)?;
    ///
    ///     Ok(spa)
    /// }
    /// ```
    ///
    pub fn from_date_time(date_time: DateTime<T>) -> Self {
        let mut input = Input::new_tz(date_time.timezone());
        input.year= date_time.year() as i64;
        input.month = date_time.month() as i64;
        input.day = date_time.day() as i64;
        input.hour = date_time.hour() as i64;
        input.minute = date_time.minute() as i64;
        input.second = date_time.second() as f64 + date_time.nanosecond() as f64 / 1_000_000_000f64;
        input.timezone = date_time.offset().fix().local_minus_utc() as f64 / 3600.0;

        Self {
            input,
        }
    }

    /// Sets all date, time and timezone fields from the given date_time parameter.
    ///
    /// This is useful if the builder was created using [SpaBuilder::from_input]
    ///
    /// This method is dependent on the feature "chrono_0_4" which will include the [chrono] crate.
    ///
    /// # Arguments
    ///
    /// * 'date_time' - a [DateTime] object including the time zone
    pub fn date_time(mut self, date_time: DateTime<T>) -> Self {
        self.input.year= date_time.year() as i64;
        self.input.month = date_time.month() as i64;
        self.input.day = date_time.day() as i64;
        self.input.hour = date_time.hour() as i64;
        self.input.minute = date_time.minute() as i64;
        self.input.second = date_time.second() as f64 + date_time.nanosecond() as f64 / 1_000_000_000f64;
        self.input.timezone = date_time.offset().fix().local_minus_utc() as f64 / 3600.0;
        self.input.tz = date_time.timezone();

        self
    }
}

/// Implementation of common methods for all T
///
impl<T: Clone> SpaBuilder<T> {
    /// Creates a new SpaBuilder from an existing [Input] struct.
    ///
    /// Useful to for instance re-run calculations after just changing one or a few settings
    ///
    /// # Arguments
    ///
    /// * 'input' - existing [Input] struct
    ///
    /// # Example
    ///
    /// Accepts defaults for atmospheric_refraction, delta_ut1 and delta_t.
    /// Parameters corresponds roughly for an existing solar park in Karlskrona, Sweden.
    /// Date, time and timezone first corresponds to sun transit for the day and then a new
    /// run is executed for later in the afternoon the same day.
    /// ```rust
    /// use spa_sra::SpaBuilder;
    /// use spa_sra::errors::SpaError;
    /// use spa_sra::spa::{Function, SpaData, Std};
    ///
    /// fn main() {
    ///
    ///     match get_spa() {
    ///         Err(e) => { println!("{}", e) },
    ///         Ok(spa) => {
    ///             println!("Altitude:                 {:.6} degrees", spa.spa_za.e);
    ///             println!("Azimuth:                  {:.6} degrees", spa.spa_za.azimuth);
    ///             println!("Incidence:                {:.6} degrees", spa.spa_za_inc.incidence);
    ///         }
    ///     }
    /// }
    ///
    /// fn get_spa() -> Result<SpaData<Std>, SpaError> {
    ///     let spa = SpaBuilder::new()
    ///         .date(2025, 8, 12)?
    ///         .time(13, 2, 28.2817)?
    ///         .timezone(2.0)?
    ///         .lat_long(56.198569, 15.637826)?
    ///         .pressure(1013.0)?
    ///         .temperature(10.0)?
    ///         .elevation(10.0)?
    ///         .slope(25.0)?
    ///         .azimuth_rotation(0.0)?
    ///         .execute(Function::SpaAll)?;
    ///
    ///     // Do something with the result here
    ///     // .
    ///     // .
    ///
    ///     // Now use the same input parameters again, but change the time
    ///     let spa = SpaBuilder::from_input(spa.input)
    ///         .time(17, 0, 0.0)?
    ///         .execute(Function::SpaAll)?;
    ///
    ///     Ok(spa)
    /// }
    /// ```
    ///
    pub fn from_input(input: Input<T>) -> Self {
        Self {
            input,
        }
    }

    /// Sets observer latitude (negative south of the equator) and longitude (negative west of Greenwich)
    ///
    /// # Arguments
    ///
    /// * 'latitude' - valid range: -90 to 90 degrees
    /// * 'longitude' - valid range: -180 to 180 degrees
    pub fn lat_long(mut self, latitude: f64, longitude: f64) -> Result<Self, SpaError> {
        if longitude.abs()     > 180.0      { return Err(SpaError{ code: 9, message: MESSAGES[9] }) };
        if latitude.abs()      > 90.0       { return Err(SpaError{ code: 10, message: MESSAGES[10] }) };

        self.input.latitude = latitude;
        self.input.longitude = longitude;

        Ok(self)
    }

    /// Sets annual average local pressure \[millibars\]
    ///
    /// # Arguments
    ///
    /// * 'pressure' - valid range: 0 to 5000 millibars
    pub fn pressure(mut self, pressure: f64) -> Result<Self, SpaError> {
        if pressure    < 0.0     || pressure    > 5000.0 { return Err(SpaError{ code: 12, message: MESSAGES[12] }) };

        self.input.pressure = pressure;

        Ok(self)
    }

    /// Sets annual average local temperature \[degrees Celsius\]
    ///
    /// # Arguments
    ///
    /// * 'temperature' - valid range: -273 to 6000 degrees Celsius
    pub fn temperature(mut self, temperature: f64) -> Result<Self, SpaError> {
        if temperature <= -273.0 || temperature > 6000.0 { return Err(SpaError{ code: 13, message: MESSAGES[13] }) };

        self.input.temperature = temperature;

        Ok(self)
    }

    /// Sets atmospheric refraction at sunrise and sunset (0.5667 deg is typical)
    ///
    /// # Arguments
    ///
    /// * 'atmos_refract' - valid range: -5 to 5 degrees
    pub fn atmospheric_refraction(mut self, atmos_refract: f64) -> Result<Self, SpaError> {
        if atmos_refract.abs() > 5.0        { return Err(SpaError{ code: 16, message: MESSAGES[16] }) };

        self.input.atmos_refract = atmos_refract;

        Ok(self)
    }

    /// Sets observer elevation \[meters\]
    ///
    /// # Arguments
    ///
    /// * 'elevation' - valid range: -6500000 or higher meters
    pub fn elevation(mut self, elevation: f64) -> Result<Self, SpaError> {
        if elevation           < -6500000.0 { return Err(SpaError{ code: 11, message: MESSAGES[11] }) };

        self.input.elevation = elevation;

        Ok(self)
    }

    /// Sets surface slope (measured from the horizontal plane).
    ///
    /// This value is used to calculate the surface incidence angle for e.g. a solar panel.
    /// The surface incidence angle is the angle between an incoming ray (like light or radar) and
    /// a line perpendicular to the surface at the point where the ray hits.
    ///
    /// No need to set this unless [Function::SpaZaInc] or [Function::SpaAll] will be used
    /// when executing.
    ///
    /// # Arguments
    ///
    /// * 'slope' - valid range: -360 to 360 degrees
    pub fn slope(mut self, slope: f64) -> Result<Self, SpaError> {
        if slope.abs()  > 360.0 { return Err(SpaError{ code: 14, message: MESSAGES[14] }) };

        self.input.slope = slope;

        Ok(self)
    }

    /// Sets surface azimuth rotation (measured from south to projection of
    /// surface normal on horizontal plane, negative east)
    ///
    /// This value is used to calculate the surface incidence angle for e.g. a solar panel.
    /// The surface incidence angle is the angle between an incoming ray (like light or radar) and
    /// a line perpendicular to the surface at the point where the ray hits.
    ///
    /// No need to set this unless [Function::SpaZaInc] or [Function::SpaAll] will be used
    /// when executing.
    ///
    /// # Arguments
    ///
    /// * 'azm_rotation' - -360 to 360 degrees
    pub fn azimuth_rotation(mut self, azm_rotation: f64) -> Result<Self, SpaError> {
        if azm_rotation.abs() > 360.0 { return Err(SpaError{ code: 15, message: MESSAGES[15] }) };

        self.input.azm_rotation = azm_rotation;

        Ok(self)
    }

    /// Fractional second difference between UTC and UT which is used
    /// to adjust UTC for earth's irregular rotation rate and is derived
    /// from observation only and is reported in this bulletin:
    /// <http://maia.usno.navy.mil/ser7/ser7.dat>,
    /// where delta_ut1 = DUT1
    ///
    /// # Arguments
    ///
    /// * 'delta_ut1' - valid range: -1 to 1 second (exclusive)
    pub fn delta_ut1(mut self, delta_ut1: f64) -> Result<Self, SpaError> {
        if delta_ut1   <= -1.0   || delta_ut1   >= 1.0   { return Err(SpaError{ code: 17, message: MESSAGES[17] }) };

        self.input.delta_ut1 = delta_ut1;

        Ok(self)
    }

    /// Difference between earth rotation time and terrestrial time
    /// It is derived from observation only and is reported in this
    /// bulletin: <http://maia.usno.navy.mil/ser7/ser7.dat>,
    /// where delta_t = 32.184 + (TAI-UTC) - DUT1
    ///
    /// # Arguments
    ///
    /// * 'delta_t' - valid range: -8000 to 8000 seconds
    pub fn delta_t(mut self, delta_t: f64) -> Result<Self, SpaError> {
        if delta_t.abs()       > 8000.0     { return Err(SpaError{ code: 7, message: MESSAGES[7] }) };

        self.input.delta_t = delta_t;

        Ok(self)
    }

    /// Builds and executes the spa calculations.
    ///
    /// It is important to call all builder functions before this since the default SpaData struct
    /// only has defaults for some parameters that doesn't change that very often (delta_ut1, delta_t and
    /// atmospheric_refraction). See [`SpaBuilder<Std>::new`] for more information.
    ///
    /// Also, if you are not interested in the incidence angle of solar radiation you can omit
    /// giving values for slope and azimuth_rotation.
    ///
    /// The argument to this function sets what output is expected:
    /// * Function::SpaZa    - calculate zenith and azimuth
    /// * Function::SpaZaInc - calculate zenith, azimuth, and incidence
    /// * Function::SpaZaRts - calculate zenith, azimuth, and sun rise/transit/set values
    /// * Function::SpaAll   - calculate all SPA output values
    ///
    /// # Arguments
    ///
    /// * 'function' - switch to choose functions for desired output (from enumeration [Function])
    pub fn execute(mut self, function: Function) -> Result<SpaData<T>, SpaError> {
        self.input.function = function;
        let mut spa_data = SpaData::new(self.input);
        spa_data.spa_calculate()?;

        Ok(spa_data)
    }
}
