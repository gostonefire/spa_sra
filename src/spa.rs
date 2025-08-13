//! This module hosts the main SPA algorithm and the [SpaData] struct together with its sub-structs. 

#[cfg(feature = "chrono_0_4")]
use chrono::{DateTime, TimeZone, Timelike, Datelike, Offset};
use std::marker::PhantomData;
use crate::constants::{B_COUNT, B_SUBCOUNT, JD_COUNT, JD_MINUS, JD_PLUS, JD_ZERO, L_COUNT, L_SUBCOUNT, R_COUNT, R_SUBCOUNT, SUN_COUNT, SUN_RADIUS, SUN_RISE, SUN_SET, SUN_TRANSIT, TERM_A, TERM_B, TERM_C, TERM_COUNT, TERM_EPS_C, TERM_EPS_D, TERM_PSI_A, TERM_PSI_B, TERM_X0, TERM_X1, TERM_X2, TERM_X3, TERM_X4, TERM_X_COUNT, TERM_Y_COUNT, Y_COUNT};
use crate::earth_periodic_terms::{B_TERMS, L_TERMS, R_TERMS};
use crate::errors::{SpaError, MESSAGES};
use crate::utils::{atmospheric_refraction_correction, deg2rad, geocentric_declination, geocentric_right_ascension, limit_degrees, observer_hour_angle, rad2deg, right_ascension_parallax_and_topocentric_dec, third_order_polynomial, topocentric_azimuth_angle, topocentric_azimuth_angle_astro, topocentric_elevation_angle, topocentric_elevation_angle_corrected, topocentric_local_hour_angle, topocentric_right_ascension, topocentric_zenith_angle};
use crate::nutation_obliquity_periodic_terms::{PE_TERMS, Y_TERMS};

/// Enumeration for function codes to select desired final outputs from SPA
///
#[derive(PartialEq, Clone, Default)]
pub enum Function {
    /// calculate zenith and azimuth
    SpaZa,
    /// calculate zenith, azimuth, and incidence
    SpaZaInc,
    /// calculate zenith, azimuth, and sun rise/transit/set values
    SpaZaRts,
    /// calculate all SPA output values
    #[default]
    SpaAll,
}

/// Marker struct when not using Chrono
#[derive(Clone)]
#[derive(Default)]
pub struct Std;

/// Input values for all [Function] functions
///
#[derive(Clone, Default)]
pub struct Input<T: Clone> {
    /// 4-digit year,      valid range: -2000 to 6000, error code: 1
    pub year: i64,
    /// 2-digit month,         valid range: 1 to  12,  error code: 2
    pub month: i64,
    /// 2-digit day,           valid range: 1 to  31,  error code: 3
    pub day: i64,
    /// Observer local hour,   valid range: 0 to  24,  error code: 4
    pub hour: i64,
    /// Observer local minute, valid range: 0 to  59,  error code: 5
    pub minute: i64,
    /// Observer local second, valid range: 0 to <60,  error code: 6
    pub second: f64,

    /// Fractional second difference between UTC and UT which is used
    /// to adjust UTC for earth's irregular rotation rate and is derived
    /// from observation only and is reported in this bulletin:
    /// <http://maia.usno.navy.mil/ser7/ser7.dat>,
    /// where delta_ut1 = DUT1
    ///
    /// valid range: -1 to 1 second (exclusive), error code 17
    pub delta_ut1: f64,

    /// Difference between earth rotation time and terrestrial time
    /// It is derived from observation only and is reported in this
    /// bulletin: <http://maia.usno.navy.mil/ser7/ser7.dat>,
    /// where delta_t = 32.184 + (TAI-UTC) - DUT1
    ///
    /// valid range: -8000 to 8000 seconds, error code: 7
    pub delta_t: f64,

    /// Observer time zone (negative west of Greenwich)
    ///
    /// valid range: -18   to   18 hours,   error code: 8
    pub timezone: f64,

    /// Observer longitude (negative west of Greenwich)
    ///
    /// valid range: -180  to  180 degrees, error code: 9
    pub longitude: f64,

    /// Observer latitude (negative south of equator)
    ///
    /// valid range: -90   to   90 degrees, error code: 10
    pub latitude: f64,

    /// Observer elevation \[meters\]
    ///
    /// valid range: -6500000 or higher meters,    error code: 11
    pub elevation: f64,

    /// Annual average local pressure \[millibars\]
    ///
    /// valid range:    0 to 5000 millibars,       error code: 12
    pub pressure: f64,

    /// Annual average local temperature \[degrees Celsius\]
    ///
    /// valid range: -273 to 6000 degrees Celsius, error code; 13
    pub temperature: f64,

    /// Surface slope (measured from the horizontal plane)
    ///
    /// valid range: -360 to 360 degrees, error code: 14
    pub slope: f64,

    /// Surface azimuth rotation (measured from south to projection of
    /// surface normal on horizontal plane, negative east)
    ///
    /// valid range: -360 to 360 degrees, error code: 15
    pub azm_rotation: f64,

    /// Atmospheric refraction at sunrise and sunset (0.5667 deg is typical)
    ///
    /// valid range: -5   to   5 degrees, error code: 16
    pub atmos_refract: f64,

    /// Switch to choose functions for desired output (from enumeration)
    pub function: Function,

    #[cfg(feature = "chrono_0_4")]
    pub(crate) tz: T,
    _phantom: PhantomData<T>,
}

impl Input<Std> {
    /// Creates a new standard [Input] struct where some values that doesn't change often are defaulted
    ///
    pub fn new_std() -> Self {
        Input {
            delta_ut1: 0.1,
            delta_t: 69.084,
            atmos_refract: 0.5667,
            #[cfg(feature = "chrono_0_4")]
            tz: Std,
            ..Default::default()
        }
    }
}

#[cfg(feature = "chrono_0_4")]
impl<T: TimeZone> Input<T> {
    /// Creates a new timezone aware [Input] struct where some values that doesn't change often are defaulted
    ///
    /// This method is dependent on the feature "chrono_0_4" which will include the [chrono] crate.
    ///
    /// # Arguments
    ///
    /// * 'tz' - a Chrono [TimeZone]
    pub fn new_tz(tz: T) -> Self {
        Input {
            year: 0,
            month: 0,
            day: 0,
            hour: 0,
            minute: 0,
            second: 0.0,
            delta_ut1: 0.1,
            delta_t: 69.084,
            timezone: 0.0,
            longitude: 0.0,
            latitude: 0.0,
            elevation: 0.0,
            pressure: 0.0,
            temperature: 0.0,
            slope: 0.0,
            azm_rotation: 0.0,
            atmos_refract: 0.5667,
            function: Default::default(),
            tz,
            _phantom: PhantomData,
        }
    }

    /// Creates a new timezone aware [Input] struct where some values that doesn't change often are defaulted.
    /// Also sets all date, time and timezone fields from the given date_time parameter.
    ///
    /// This method is dependent on the feature "chrono_0_4" which will include the [chrono] crate.
    ///
    /// # Arguments
    ///
    /// * 'date_time' - a [DateTime] object including the time zone
    pub fn from_date_time(date_time: DateTime<T>) -> Self {
        Input {
            year: date_time.year() as i64,
            month: date_time.month() as i64,
            day: date_time.day() as i64,
            hour: date_time.hour() as i64,
            minute: date_time.minute() as i64,
            second: date_time.second() as f64 + date_time.nanosecond() as f64 / 1_000_000_000f64,
            delta_ut1: 0.1,
            delta_t: 69.084,
            timezone: date_time.offset().fix().local_minus_utc() as f64 / 3600.0,
            longitude: 0.0,
            latitude: 0.0,
            elevation: 0.0,
            pressure: 0.0,
            temperature: 0.0,
            slope: 0.0,
            azm_rotation: 0.0,
            atmos_refract: 0.5667,
            function: Default::default(),
            tz: date_time.timezone(),
            _phantom: PhantomData,
        }
    }

    /// Sets all date, time and timezone fields from the given date_time parameter.
    ///
    /// This method is dependent on the feature "chrono_0_4" which will include the [chrono] crate.
    ///
    /// # Arguments
    ///
    /// * 'date_time' - a [DateTime] object including the time zone
    pub fn date_time(&mut self, date_time: DateTime<T>) {
        self.tz = date_time.timezone();
        self.year= date_time.year() as i64;
        self.month = date_time.month() as i64;
        self.day = date_time.day() as i64;
        self.hour = date_time.hour() as i64;
        self.minute = date_time.minute() as i64;
        self.second = date_time.second() as f64 + date_time.nanosecond() as f64 / 1_000_000_000f64;
        self.timezone = date_time.offset().fix().local_minus_utc() as f64 / 3600.0;
    }
}

impl<T: Clone> Input<T> {
    /// Validates inputs
    ///
    fn validate_inputs(&self) -> i64 {
        if self.year        < -2000   || self.year        > 6000   { return 1 };
        if self.month       < 1       || self.month       > 12     { return 2 };
        if self.day         < 1       || self.day         > 31     { return 3 };
        if self.hour        < 0       || self.hour        > 24     { return 4 };
        if self.minute      < 0       || self.minute      > 59     { return 5 };
        if self.second      < 0.0     || self.second      >=60.0   { return 6 };
        if self.pressure    < 0.0     || self.pressure    > 5000.0 { return 12 };
        if self.temperature <= -273.0 || self.temperature > 6000.0 { return 13 };
        if self.delta_ut1   <= -1.0   || self.delta_ut1   >= 1.0   { return 17 };
        if self.hour        == 24     && self.minute      > 0      { return 5 };
        if self.hour        == 24     && self.second      > 0.0    { return 6 };

        if self.delta_t.abs()       > 8000.0     { return 7 };
        if self.timezone.abs()      > 18.0       { return 8 };
        if self.longitude.abs()     > 180.0      { return 9 };
        if self.latitude.abs()      > 90.0       { return 10 };
        if self.atmos_refract.abs() > 5.0        { return 16 };
        if self.elevation           < -6500000.0 { return 11 };

        if self.function == Function::SpaZaInc || self.function == Function::SpaAll {
            if self.slope.abs()  > 360.0 { return 14 };
            if self.azm_rotation.abs() > 360.0 { return 15 };
        }

        0
    }
}

/// Intermediate and final output values for all available [Function] functions
///
#[derive(Clone, Default)]
pub struct SpaZa {
    //-----------------Intermediate OUTPUT VALUES--------------------

    /// Julian day
    pub jd: f64,
    /// Julian century
    pub jc: f64,

    /// Julian ephemeris day
    pub jde: f64,
    /// Julian ephemeris century
    pub jce: f64,
    /// Julian ephemeris millennium
    pub jme: f64,

    /// earth heliocentric longitude \[degrees\]
    pub l: f64,
    /// earth heliocentric latitude \[degrees\]
    pub b: f64,
    /// earth radius vector \[Astronomical Units, AU\]
    pub r: f64,

    /// geocentric longitude \[degrees\]
    pub theta: f64,
    /// geocentric latitude \[degrees\]
    pub beta: f64,

    /// mean elongation (moon-sun) \[degrees\]
    pub x0: f64,
    /// mean anomaly (sun) \[degrees\]
    pub x1: f64,
    /// mean anomaly (moon) \[degrees\]
    pub x2: f64,
    /// argument latitude (moon) \[degrees\]
    pub x3: f64,
    /// ascending longitude (moon) \[degrees\]
    pub x4: f64,

    /// nutation longitude \[degrees\]
    pub del_psi: f64,
    /// nutation obliquity \[degrees\]
    pub del_epsilon: f64,
    /// ecliptic mean obliquity \[arc seconds\]
    pub epsilon0: f64,
    /// ecliptic true obliquity  \[degrees\]
    pub epsilon: f64,

    /// aberration correction \[degrees\]
    pub del_tau: f64,
    /// apparent sun longitude \[degrees\]
    pub lamda: f64,
    /// Greenwich mean sidereal time \[degrees\]
    pub nu0: f64,
    /// Greenwich sidereal time \[degrees\]
    pub nu: f64,

    /// geocentric sun right ascension \[degrees\]
    pub alpha: f64,
    /// geocentric sun declination \[degrees\]
    pub delta: f64,

    /// observer hour angle \[degrees\]
    pub h: f64,
    /// sun equatorial horizontal parallax \[degrees\]
    pub xi: f64,
    /// sun right ascension parallax \[degrees\]
    pub del_alpha: f64,
    /// topocentric sun declination \[degrees\]
    pub delta_prime: f64,
    /// topocentric sun right ascension \[degrees\]
    pub alpha_prime: f64,
    /// topocentric local hour angle \[degrees\]
    pub h_prime: f64,

    /// topocentric elevation angle (uncorrected) \[degrees\]
    pub e0: f64,
    /// atmospheric refraction correction \[degrees\]
    pub del_e: f64,
    /// topocentric elevation angle (corrected) \[degrees\]
    pub e: f64,

    /// topocentric zenith angle \[degrees\]
    pub zenith: f64,
    /// topocentric azimuth angle (westward from south) \[for astronomers\]
    pub azimuth_astro: f64,
    /// topocentric azimuth angle (eastward from north) \[for navigators and solar radiation\]
    pub azimuth: f64,

}

impl SpaZa {
    /// Creates a new `SpaZa` struct with zeroed values
    ///
    fn new() -> Self {
        SpaZa {
            ..Default::default()
        }
    }
}

/// Final output values for [Function::SpaZaInc] and [Function::SpaAll]
///
#[derive(Clone)]
pub struct SpaZaInc {
    /// surface incidence angle \[degrees\]
    pub incidence: f64,
}

/// Intermediate and final output values for [Function::SpaZaRts] and [Function::SpaAll]
///
#[derive(Clone, Default)]
pub struct SpaZaRts {
    //-----------------Intermediate OUTPUT VALUES--------------------

    /// equation of time \[minutes\]
    pub eot: f64,
    /// sunrise hour angle \[degrees\]
    pub srha: f64,
    /// sunset hour angle \[degrees\]
    pub ssha: f64,
    /// sun transit altitude \[degrees\]
    pub sta: f64,

    //---------------------Final OUTPUT VALUES------------------------

    /// local sun transit time (or solar noon) \[fractional hour\]
    pub suntransit: f64,
    /// local sunrise time (+/- 30 seconds) \[fractional hour\]
    pub sunrise: f64,
    /// local sunset time (+/- 30 seconds) \[fractional hour\]
    pub sunset: f64,
}

impl SpaZaRts {
    /// Creates a new `SpaZaRts` struct with zeroed values
    ///
    fn new() -> Self {
        SpaZaRts {
            ..Default::default()
        }
    }
}

/// The main SpaData struct that holds both input data and structs for result data
/// after a Spa calculation
///
#[derive(Clone)]
pub struct SpaData<T: Clone> {
    /// Input values for all [Function] functions
    pub input: Input<T>,

    /// Intermediate and final output values for all available [Function] functions
    pub spa_za: SpaZa,

    /// Final output values for [Function::SpaZaInc] and [Function::SpaAll]
    pub spa_za_inc: SpaZaInc,

    /// Intermediate and final output values for [Function::SpaZaRts] and [Function::SpaAll]
    pub spa_za_rts: SpaZaRts,
}

/// Implementation of Chrono_0_4 specific methods
///
#[cfg(feature = "chrono_0_4")]
impl<T: TimeZone> SpaData<T> {
    /// Returns the sunrise as a [DateTime] object including nanoseconds
    ///
    /// This method is dependent on the feature "chrono_0_4" which will include the [chrono] crate.
    ///
    pub fn get_sunrise(&self) -> DateTime<T> {
        let time_comp = get_time_components(self.spa_za_rts.sunrise);

        self.input.tz.with_ymd_and_hms(self.input.year as i32, self.input.month as u32, self.input.day as u32,
                            time_comp.0, time_comp.1, time_comp.2).unwrap()
            .with_nanosecond(time_comp.3).unwrap()
    }

    /// Returns the sunset as a [DateTime] object including nanoseconds
    ///
    /// This method is dependent on the feature "chrono_0_4" which will include the [chrono] crate.
    ///
    pub fn get_sunset(&self) -> DateTime<T> {
        let time_comp = get_time_components(self.spa_za_rts.sunset);

        self.input.tz.with_ymd_and_hms(self.input.year as i32, self.input.month as u32, self.input.day as u32,
                            time_comp.0, time_comp.1, time_comp.2).unwrap()
            .with_nanosecond(time_comp.3).unwrap()
    }

    /// Returns the suntransit as a [DateTime] object including nanoseconds
    ///
    /// This method is dependent on the feature "chrono_0_4" which will include the [chrono] crate.
    ///
    pub fn get_suntransit(&self) -> DateTime<T> {
        let time_comp = get_time_components(self.spa_za_rts.suntransit);

        self.input.tz.with_ymd_and_hms(self.input.year as i32, self.input.month as u32, self.input.day as u32,
                            time_comp.0, time_comp.1, time_comp.2).unwrap()
            .with_nanosecond(time_comp.3).unwrap()
    }
}

/// Implementation of common methods for all T
///
impl<T: Clone> SpaData<T> {
    /// Returns a new [SpaData] struct
    ///
    /// # Arguments
    ///
    /// * 'input' - an instance of the [Input] struct
    pub fn new(input: Input<T>) -> Self {
        SpaData {
            input,
            spa_za: SpaZa::new(),
            spa_za_inc: SpaZaInc { incidence: 0.0 },
            spa_za_rts: SpaZaRts::new(),
        }
    }

    /// Calculate all SPA parameters and put into structure
    ///
    pub fn spa_calculate(&mut self) -> Result<(), SpaError> {
        let result: i64 = self.input.validate_inputs();
        if result == 0 {
            self.spa_za.jd = julian_day(self.input.year, self.input.month, self.input.day, self.input.hour,
                                        self.input.minute, self.input.second, self.input.delta_ut1, self.input.timezone);

            self.calculate_geocentric_sun_right_ascension_and_declination();

            self.spa_za.h  = observer_hour_angle(self.spa_za.nu, self.input.longitude, self.spa_za.alpha);
            self.spa_za.xi = sun_equatorial_horizontal_parallax(self.spa_za.r);

            right_ascension_parallax_and_topocentric_dec(self.input.latitude, self.input.elevation, self.spa_za.xi,
                                                         self.spa_za.h, self.spa_za.delta, &mut self.spa_za.del_alpha, &mut self.spa_za.delta_prime);

            self.spa_za.alpha_prime = topocentric_right_ascension(self.spa_za.alpha, self.spa_za.del_alpha);
            self.spa_za.h_prime     = topocentric_local_hour_angle(self.spa_za.h, self.spa_za.del_alpha);

            self.spa_za.e0      = topocentric_elevation_angle(self.input.latitude, self.spa_za.delta_prime, self.spa_za.h_prime);
            self.spa_za.del_e   = atmospheric_refraction_correction(self.input.pressure, self.input.temperature,
                                                                    self.input.atmos_refract, self.spa_za.e0);
            self.spa_za.e       = topocentric_elevation_angle_corrected(self.spa_za.e0, self.spa_za.del_e);

            self.spa_za.zenith        = topocentric_zenith_angle(self.spa_za.e);
            self.spa_za.azimuth_astro = topocentric_azimuth_angle_astro(self.spa_za.h_prime, self.input.latitude,
                                                                        self.spa_za.delta_prime);
            self.spa_za.azimuth       = topocentric_azimuth_angle(self.spa_za.azimuth_astro);

            if self.input.function == Function::SpaZaInc || self.input.function == Function::SpaAll {
                self.spa_za_inc.incidence  = surface_incidence_angle(self.spa_za.zenith, self.spa_za.azimuth_astro,
                                                                     self.input.azm_rotation, self.input.slope);
            }

            if self.input.function == Function::SpaZaRts || self.input.function == Function::SpaAll {
                self.calculate_eot_and_sun_rise_transit_set();
            }

            Ok(())

        } else {
            Err(SpaError{ code: result, message: MESSAGES[result as usize] })
        }
    }

    /// Calculate Equation of Time (EOT) and Sun Rise, Transit, & Set (RTS)
    ///
    fn calculate_eot_and_sun_rise_transit_set(&mut self) {
        let mut sun_rts: SpaData<T> = self.clone();

        let mut alpha: [f64;JD_COUNT] = [0.0; JD_COUNT];
        let mut delta: [f64;JD_COUNT] = [0.0; JD_COUNT];

        let mut m_rts: [f64;SUN_COUNT] = [0.0; SUN_COUNT];
        let mut nu_rts: [f64;SUN_COUNT] = [0.0; SUN_COUNT];
        let mut h_rts: [f64;SUN_COUNT] = [0.0; SUN_COUNT];

        let mut alpha_prime: [f64;SUN_COUNT] = [0.0; SUN_COUNT];
        let mut delta_prime: [f64;SUN_COUNT] = [0.0; SUN_COUNT];
        let mut h_prime: [f64;SUN_COUNT] = [0.0; SUN_COUNT];

        let h0_prime: f64 = -1.0 * (SUN_RADIUS + self.input.atmos_refract);

        let m: f64 = sun_mean_longitude(self.spa_za.jme);
        self.spa_za_rts.eot = eot(m, self.spa_za.alpha, self.spa_za.del_psi, self.spa_za.epsilon);

        sun_rts.input.hour = 0; sun_rts.input.minute = 0; sun_rts.input.second = 0.0;
        sun_rts.input.delta_ut1 = 0.0; sun_rts.input.timezone = 0.0;

        sun_rts.spa_za.jd = julian_day (sun_rts.input.year,   sun_rts.input.month,  sun_rts.input.day,       sun_rts.input.hour,
                                        sun_rts.input.minute, sun_rts.input.second, sun_rts.input.delta_ut1, sun_rts.input.timezone);

        sun_rts.calculate_geocentric_sun_right_ascension_and_declination();
        let nu: f64 = sun_rts.spa_za.nu;

        sun_rts.input.delta_t = 0.0;
        sun_rts.spa_za.jd -= 1.0;

        for i in 0..JD_COUNT {
            sun_rts.calculate_geocentric_sun_right_ascension_and_declination();
            alpha[i] = sun_rts.spa_za.alpha;
            delta[i] = sun_rts.spa_za.delta;
            sun_rts.spa_za.jd += 1.0;
        }

        m_rts[SUN_TRANSIT] = approx_sun_transit_time(alpha[JD_ZERO], self.input.longitude, nu);
        let h0: f64 = sun_hour_angle_at_rise_set(self.input.latitude, delta[JD_ZERO], h0_prime);

        if h0 >= 0.0 {
            approx_sun_rise_and_set(&mut m_rts, h0);

            for i in 0..SUN_COUNT {
                nu_rts[i]      = nu + 360.985647*m_rts[i];

                let n: f64     = m_rts[i] + self.input.delta_t / 86400.0;
                alpha_prime[i] = rts_alpha_delta_prime(&alpha, n);
                delta_prime[i] = rts_alpha_delta_prime(&delta, n);

                h_prime[i]     = limit_degrees180pm(nu_rts[i] + self.input.longitude - alpha_prime[i]);

                h_rts[i]       = rts_sun_altitude(self.input.latitude, delta_prime[i], h_prime[i]);
            }

            self.spa_za_rts.srha = h_prime[SUN_RISE];
            self.spa_za_rts.ssha = h_prime[SUN_SET];
            self.spa_za_rts.sta  = h_rts[SUN_TRANSIT];

            self.spa_za_rts.suntransit = dayfrac_to_local_hr(m_rts[SUN_TRANSIT] - h_prime[SUN_TRANSIT] / 360.0, self.input.timezone);

            self.spa_za_rts.sunrise = dayfrac_to_local_hr(sun_rise_and_set(&m_rts, &h_rts, &delta_prime,
                                                                           self.input.latitude, &h_prime, h0_prime, SUN_RISE), self.input.timezone);

            self.spa_za_rts.sunset  = dayfrac_to_local_hr(sun_rise_and_set(&m_rts, &h_rts, &delta_prime,
                                                                           self.input.latitude, &h_prime, h0_prime, SUN_SET), self.input.timezone);

        } else {
            self.spa_za_rts.srha       = -99999.0;
            self.spa_za_rts.ssha       = -99999.0;
            self.spa_za_rts.sta        = -99999.0;
            self.spa_za_rts.suntransit = -99999.0;
            self.spa_za_rts.sunrise    = -99999.0;
            self.spa_za_rts.sunset     = -99999.0;
        }
    }

    /// Calculate required SPA parameters to get the right ascension (alpha) and declination (delta)
    /// Note: JD must be already calculated and in structure
    ///
    /// # Arguments
    ///
    /// * 'spa' - the [SpaData] struct
    fn calculate_geocentric_sun_right_ascension_and_declination(&mut self) {
        let mut x: [f64;TERM_X_COUNT] = [0.0; TERM_X_COUNT];

        self.spa_za.jc = julian_century(self.spa_za.jd);

        self.spa_za.jde = julian_ephemeris_day(self.spa_za.jd, self.input.delta_t);
        self.spa_za.jce = julian_ephemeris_century(self.spa_za.jde);
        self.spa_za.jme = julian_ephemeris_millennium(self.spa_za.jce);

        self.spa_za.l = earth_heliocentric_longitude(self.spa_za.jme);
        self.spa_za.b = earth_heliocentric_latitude(self.spa_za.jme);
        self.spa_za.r = earth_radius_vector(self.spa_za.jme);

        self.spa_za.theta = geocentric_longitude(self.spa_za.l);
        self.spa_za.beta  = geocentric_latitude(self.spa_za.b);

        x[TERM_X0] = { self.spa_za.x0 = mean_elongation_moon_sun(self.spa_za.jce); self.spa_za.x0 };
        x[TERM_X1] = { self.spa_za.x1 = mean_anomaly_sun(self.spa_za.jce);         self.spa_za.x1 };
        x[TERM_X2] = { self.spa_za.x2 = mean_anomaly_moon(self.spa_za.jce);        self.spa_za.x2 };
        x[TERM_X3] = { self.spa_za.x3 = argument_latitude_moon(self.spa_za.jce);   self.spa_za.x3 };
        x[TERM_X4] = { self.spa_za.x4 = ascending_longitude_moon(self.spa_za.jce); self.spa_za.x4 };

        self.nutation_longitude_and_obliquity(&x);

        self.spa_za.epsilon0 = ecliptic_mean_obliquity(self.spa_za.jme);
        self.spa_za.epsilon  = ecliptic_true_obliquity(self.spa_za.del_epsilon, self.spa_za.epsilon0);

        self.spa_za.del_tau   = aberration_correction(self.spa_za.r);
        self.spa_za.lamda     = apparent_sun_longitude(self.spa_za.theta, self.spa_za.del_psi, self.spa_za.del_tau);
        self.spa_za.nu0       = greenwich_mean_sidereal_time(self.spa_za.jd, self.spa_za.jc);
        self.spa_za.nu        = greenwich_sidereal_time(self.spa_za.nu0, self.spa_za.del_psi, self.spa_za.epsilon);

        self.spa_za.alpha = geocentric_right_ascension(self.spa_za.lamda, self.spa_za.epsilon, self.spa_za.beta);
        self.spa_za.delta = geocentric_declination(self.spa_za.beta, self.spa_za.epsilon, self.spa_za.lamda);
    }

    fn nutation_longitude_and_obliquity(&mut self, x: &[f64]) {
        let mut sum_psi: f64 = 0.0;
        let mut sum_epsilon: f64 = 0.0;

        for i in 0..Y_COUNT {
            let xy_term_sum  = deg2rad(xy_term_summation(i, x));
            sum_psi     += (PE_TERMS[i][TERM_PSI_A] + self.spa_za.jce * PE_TERMS[i][TERM_PSI_B]) * xy_term_sum.sin();
            sum_epsilon += (PE_TERMS[i][TERM_EPS_C] + self.spa_za.jce * PE_TERMS[i][TERM_EPS_D]) * xy_term_sum.cos();
        }

        self.spa_za.del_psi = sum_psi / 36000000.0;
        self.spa_za.del_epsilon = sum_epsilon / 36000000.0;
    }
}

fn integer(value: f64) -> f64 {
    (value as i64) as f64
}

fn limit_degrees180pm(mut degrees: f64) -> f64 {
    degrees /= 360.0;
    let limited = 360.0 * (degrees - degrees.floor());

    if limited < -180.0 {
        limited + 360.0
    } else if limited > 180.0 {
        limited - 360.0
    } else {
        limited
    }
}

fn limit_degrees180(mut degrees: f64) -> f64 {
    degrees /= 180.0;
    let limited = 180.0 * (degrees - degrees.floor());

    if limited < 0.0 {
        limited + 180.0
    } else {
        limited
    }
}

fn limit_zero2one(value: f64) -> f64 {
    let limited = value - value.floor();

    if limited < 0.0 {
        limited + 1.0
    } else {
        limited
    }
}


fn limit_minutes(minutes: f64) -> f64 {
    let limited = minutes;

    if limited < -20.0 {
        limited + 1440.0
    } else if limited > 20.0 {
        limited - 1440.0
    } else {
        limited
    }
}

fn dayfrac_to_local_hr(dayfrac: f64, timezone: f64) -> f64 {
    24.0 * limit_zero2one(dayfrac + timezone / 24.0)
}

fn julian_day(mut year: i64, mut month: i64, day: i64, hour: i64, minute: i64, second: f64, dut1: f64, tz: f64) -> f64 {
    let day_decimal: f64 = day as f64 + (hour as f64 - tz + (minute as f64 + (second + dut1) / 60.0) / 60.0) / 24.0;

    if month < 3 {
        month += 12;
        year -= 1;
    }

    let mut julian_day: f64 = integer(365.25 * (year as f64 + 4716.0)) + integer(30.6001 * (month as f64 + 1.0)) + day_decimal - 1524.5;

    if julian_day > 2299160.0 {
        let a = integer(year as f64 / 100.0);
        julian_day += 2.0 - a + integer(a / 4.0);
    }

    julian_day
}

fn julian_century(jd: f64) -> f64 {
    (jd - 2451545.0) / 36525.0
}

fn julian_ephemeris_day(jd: f64, delta_t: f64) -> f64 {
    jd + delta_t / 86400.0
}

fn julian_ephemeris_century(jde: f64) -> f64 {
    (jde - 2451545.0) / 36525.0
}

fn julian_ephemeris_millennium(jce: f64) -> f64 {
    jce / 10.0
}

fn earth_periodic_term_summation(terms: &[[f64;TERM_COUNT]], count: i64, jme: f64) -> f64 {
    let mut sum:f64 = 0.0;
    for i in 0..count as usize {
        sum += terms[i][TERM_A] * (terms[i][TERM_B] + terms[i][TERM_C] * jme).cos();
    }

    sum
}

fn earth_values(term_sum: &[f64], count: usize, jme: f64) -> f64 {
    let mut sum:f64 = 0.0;
    for i in 0..count {
        sum += term_sum[i] * jme.powi(i as i32);
    }

    sum /= 1.0e8;

    sum
}

fn earth_heliocentric_longitude(jme: f64) -> f64 {
    let mut sum: [f64;L_COUNT] = [0.0; L_COUNT];
    for i in 0..L_COUNT {
        sum[i] = earth_periodic_term_summation(&L_TERMS[i], L_SUBCOUNT[i], jme);
    }

    limit_degrees(rad2deg(earth_values(&sum, L_COUNT, jme)))
}

fn earth_heliocentric_latitude(jme: f64) -> f64 {
    let mut sum: [f64;B_COUNT] = [0.0; B_COUNT];
    for i in 0..B_COUNT {
        sum[i] = earth_periodic_term_summation(&B_TERMS[i], B_SUBCOUNT[i], jme);
    }

    rad2deg(earth_values(&sum, B_COUNT, jme))
}

fn earth_radius_vector(jme: f64) -> f64 {
    let mut sum: [f64;R_COUNT] = [0.0; R_COUNT];
    for i in 0..R_COUNT {
        sum[i] = earth_periodic_term_summation(&R_TERMS[i], R_SUBCOUNT[i], jme);
    }

    earth_values(&sum, R_COUNT, jme)
}

fn geocentric_longitude(l: f64) -> f64 {
    let theta: f64 = l + 180.0;

    if theta >= 360.0 {
        theta - 360.0
    } else {
        theta
    }
}

fn geocentric_latitude(b: f64) -> f64 {
    -b
}

fn mean_elongation_moon_sun(jce: f64) -> f64 {
    third_order_polynomial(1.0 / 189474.0, -0.0019142, 445267.11148, 297.85036, jce)
}

fn mean_anomaly_sun(jce: f64) -> f64 {
    third_order_polynomial(-1.0 / 300000.0, -0.0001603, 35999.05034, 357.52772, jce)
}

fn mean_anomaly_moon(jce: f64) -> f64 {
    third_order_polynomial(1.0 / 56250.0, 0.0086972, 477198.867398, 134.96298, jce)
}

fn argument_latitude_moon(jce: f64) -> f64 {
    third_order_polynomial(1.0 / 327270.0, -0.0036825, 483202.017538, 93.27191, jce)
}

fn ascending_longitude_moon(jce: f64) -> f64 {
    third_order_polynomial(1.0 / 450000.0, 0.0020708, -1934.136261, 125.04452, jce)
}

fn xy_term_summation(i: usize, x: &[f64]) -> f64 {
    let mut sum:f64 = 0.0;
    for j in 0..TERM_Y_COUNT {
        sum += x[j] * Y_TERMS[i][j] as f64;
    }

    sum
}

fn ecliptic_mean_obliquity(jme: f64) -> f64 {
    let u: f64 = jme / 10.0;

    84381.448 + u * (-4680.93 + u * (-1.55 + u * (1999.25 + u * (-51.38 + u * (-249.67 +
        u * (  -39.05 + u * ( 7.12 + u * (  27.87 + u * (  5.79 + u * 2.45)))))))))
}

fn ecliptic_true_obliquity(delta_epsilon: f64, epsilon0: f64) -> f64 {
    delta_epsilon + epsilon0 / 3600.0
}

fn aberration_correction(r: f64) -> f64 {
    -20.4898 / (3600.0 * r)
}

fn apparent_sun_longitude(theta: f64, delta_psi: f64, delta_tau: f64) -> f64 {
    theta + delta_psi + delta_tau
}

fn greenwich_mean_sidereal_time(jd: f64, jc: f64) -> f64 {
    limit_degrees(280.46061837 + 360.98564736629 * (jd - 2451545.0) +
        jc * jc * (0.000387933 - jc/38710000.0))
}

fn greenwich_sidereal_time(nu0: f64, delta_psi: f64, epsilon: f64) -> f64 {
    nu0 + delta_psi * deg2rad(epsilon).cos()
}

fn sun_equatorial_horizontal_parallax(r: f64) -> f64 {
    8.794 / (3600.0 * r)
}

fn surface_incidence_angle(zenith: f64, azimuth_astro: f64, azm_rotation: f64, slope: f64) -> f64 {
    let zenith_rad: f64 = deg2rad(zenith);
    let slope_rad: f64  = deg2rad(slope);

    rad2deg((zenith_rad.cos() * slope_rad.cos()  +
        slope_rad.sin() * zenith_rad.sin() * deg2rad(azimuth_astro - azm_rotation).cos()).acos())
}

fn sun_mean_longitude(jme: f64) -> f64 {
    limit_degrees(280.4664567 + jme * (360007.6982779 + jme * (0.03032028 +
        jme * (1.0 / 49931.0 + jme * (-1.0 / 15300.0 + jme * (-1.0 / 2000000.0))))))
}

fn eot(m: f64, alpha: f64, del_psi: f64, epsilon: f64) -> f64 {
    limit_minutes(4.0 * (m - 0.0057183 - alpha + del_psi * deg2rad(epsilon).cos()))
}

fn approx_sun_transit_time(alpha_zero: f64, longitude: f64, nu: f64) -> f64 {
    (alpha_zero - longitude - nu) / 360.0
}

fn sun_hour_angle_at_rise_set(latitude: f64, delta_zero: f64, h0_prime: f64) -> f64 {
    let mut h0: f64             = -99999.0;
    let latitude_rad: f64   = deg2rad(latitude);
    let delta_zero_rad: f64 = deg2rad(delta_zero);
    let argument: f64       = (deg2rad(h0_prime).sin() - latitude_rad.sin() * delta_zero_rad.sin()) /
        (latitude_rad.cos() * delta_zero_rad.cos());

    if argument.abs() <= 1.0 {
        h0 = limit_degrees180(rad2deg(argument.acos()));
    }

    h0
}

fn approx_sun_rise_and_set(m_rts: &mut [f64], h0: f64) {
    let h0_dfrac: f64 = h0 / 360.0;

    m_rts[SUN_RISE]    = limit_zero2one(m_rts[SUN_TRANSIT] - h0_dfrac);
    m_rts[SUN_SET]     = limit_zero2one(m_rts[SUN_TRANSIT] + h0_dfrac);
    m_rts[SUN_TRANSIT] = limit_zero2one(m_rts[SUN_TRANSIT]);
}

fn rts_alpha_delta_prime(ad: &[f64], n: f64) -> f64 {
    let mut a: f64 = ad[JD_ZERO] - ad[JD_MINUS];
    let mut b: f64 = ad[JD_PLUS] - ad[JD_ZERO];

    if a.abs() >= 2.0 {
        a = limit_zero2one(a);
    }
    if b.abs() >= 2.0 {
        b = limit_zero2one(b);
    }

    ad[JD_ZERO] + n * (a + b + (b - a) * n) / 2.0
}

fn rts_sun_altitude(latitude: f64, delta_prime: f64, h_prime: f64) -> f64 {
    let latitude_rad: f64    = deg2rad(latitude);
    let delta_prime_rad: f64 = deg2rad(delta_prime);

    rad2deg((latitude_rad.sin() * delta_prime_rad.sin() +
        latitude_rad.cos() * delta_prime_rad.cos() * deg2rad(h_prime).cos()).asin())
}

fn sun_rise_and_set(m_rts: &[f64], h_rts: &[f64], delta_prime: &[f64], latitude: f64, h_prime: &[f64], h0_prime: f64, sun: usize) -> f64 {
    m_rts[sun] + (h_rts[sun] - h0_prime) /
        (360.0 * deg2rad(delta_prime[sun]).cos() * deg2rad(latitude).cos() * deg2rad(h_prime[sun]).sin())
}

/// Returns time components given an hour with fractions
///
/// # Arguments
///
/// * 'frac_time' - hour with fraction
#[cfg(feature = "chrono_0_4")]
fn get_time_components(frac_time: f64) -> (u32, u32, u32, u32) {
    let min: f64 = 60.0 * (frac_time - (frac_time as i64) as f64);
    let sec: f64 = 60.0 * (min - (min as i64) as f64);
    let nano: u32 = ((sec - (sec as i64) as f64) * 1_000_000_000f64) as u32;

        (frac_time as u32, min as u32, sec as u32, nano)
}


#[cfg(test)]
mod tests {
    use crate::spa::{Function, Input, SpaData};

    /// Checks spa calculation output against original test results
    ///
    fn check_test_data<T: Clone>(spa: &SpaData<T>) {
        assert_eq!(format!("{:.6}", spa.spa_za.jd), "2452930.312847",            "Julian Day:    2452930.312847");
        assert_eq!(format!("{:.6e}", spa.spa_za.l), "2.401826e1",                "L:             2.401826e+01 degrees");
        assert_eq!(format!("{:.6e}", spa.spa_za.b), "-1.011219e-4",              "B:             -1.011219e-04 degrees");
        assert_eq!(format!("{:.6}", spa.spa_za.r), "0.996542",                   "R:             0.996542 AU");
        assert_eq!(format!("{:.6}", spa.spa_za.h), "11.105902",                  "H:             11.105902 degrees");
        assert_eq!(format!("{:.6e}", spa.spa_za.del_psi), "-3.998404e-3",        "Delta Psi:     -3.998404e-03 degrees");
        assert_eq!(format!("{:.6e}", spa.spa_za.del_epsilon), "1.666568e-3",     "Delta Epsilon: 1.666568e-03 degrees");
        assert_eq!(format!("{:.6}", spa.spa_za.epsilon), "23.440465",            "Epsilon:       23.440465 degrees");
        assert_eq!(format!("{:.6}", spa.spa_za.zenith), "50.111622",             "Zenith:        50.111622 degrees");
        assert_eq!(format!("{:.6}", spa.spa_za.azimuth), "194.340241",           "Azimuth:       194.340241 degrees");
        assert_eq!(format!("{:.6}", spa.spa_za_inc.incidence), "25.187000",      "Incidence:     25.187000 degrees");

        let mut min: f64 = 60.0 * (spa.spa_za_rts.sunrise - (spa.spa_za_rts.sunrise as i64) as f64);
        let mut sec: f64 = 60.0 * (min - (min as i64) as f64);
        assert_eq!(format!("{:0>2}:{:0>2}:{:0>2}", spa.spa_za_rts.sunrise as i64, min as i64, sec as i64), "06:12:43",
                   "Sunrise: 06:12:43 Local Time");

        min = 60.0 * (spa.spa_za_rts.sunset - (spa.spa_za_rts.sunset as i64) as f64);
        sec = 60.0 * (min - (min as i64) as f64);
        assert_eq!(format!("{:0>2}:{:0>2}:{:0>2}", spa.spa_za_rts.sunset as i64, min as i64, sec as i64), "17:20:19",
                   "Sunset: 17:20:19 Local Time");
    }

    #[test]
    fn spa_data_new_std() {
        // enter required input values into SPA structure
        let mut input = Input::new_std();
        input.year          =  2003;
        input.month         =  10;
        input.day           =  17;
        input.hour          =  12;
        input.minute        =  30;
        input.second        =  30.0;
        input.timezone      = -7.0;
        input.delta_ut1     =  0.0;
        input.delta_t       =  67.0;
        input.longitude     = -105.1786;
        input.latitude      =  39.742476;
        input.elevation     =  1830.14;
        input.pressure      =  820.0;
        input.temperature   =  11.0;
        input.slope         =  30.0;
        input.azm_rotation  = -10.0;
        input.atmos_refract =  0.5667;
        input.function      = Function::SpaAll;

        let mut spa = SpaData::new(input);

        // call the SPA calculate function
        if let Err(e) = spa.spa_calculate() {
            panic!("{}", e);
        }

        // test and test results according original code
        check_test_data(&spa);
    }

    #[cfg(feature = "chrono_0_4")]
    #[test]
    fn spa_data_new_tz() {
        use chrono::{FixedOffset, TimeZone};

        // Create timezone of interest
        let tz = FixedOffset::east_opt(-7 * 3600).unwrap();

        // Get a new Input instance with the concrete type of FixedOffset and set parameters
        let mut input = Input::new_tz(tz);
        input.delta_ut1 = 0.0;
        input.delta_t = 67.0;
        input.longitude = -105.1786;
        input.latitude = 39.742476;
        input.elevation = 1830.14;
        input.pressure = 820.0;
        input.temperature = 11.0;
        input.slope = 30.0;
        input.azm_rotation = -10.0;
        input.atmos_refract = 0.5667;
        input.function = Function::SpaAll;

        // Create a DateTime instance from timezone
        let date_time = tz
            .with_ymd_and_hms(2003, 10, 17, 12, 30, 30)
            .unwrap();

        // Set date, time and timezone (numerical for the algorithm)
        input.date_time(date_time);

        // Get a new SpaData instance
        let mut spa = SpaData::new(input);

        // call the SPA calculate function
        if let Err(e) = spa.spa_calculate() {
            panic!("{}", e);
        }

        // test and test results according original code
        check_test_data(&spa);
    }

    #[cfg(feature = "chrono_0_4")]
    #[test]
    fn spa_data_from_date_time() {
        use chrono::{FixedOffset, TimeZone};

        // Create timezone of interest
        let date_time = FixedOffset::east_opt(-7 * 3600)
            .unwrap()
            .with_ymd_and_hms(2003, 10, 17, 12, 30, 30)
            .unwrap();

        // Get a new Input instance with the concrete type of FixedOffset and set parameters
        let mut input = Input::from_date_time(date_time);
        input.delta_ut1 = 0.0;
        input.delta_t = 67.0;
        input.longitude = -105.1786;
        input.latitude = 39.742476;
        input.elevation = 1830.14;
        input.pressure = 820.0;
        input.temperature = 11.0;
        input.slope = 30.0;
        input.azm_rotation = -10.0;
        input.atmos_refract = 0.5667;
        input.function = Function::SpaAll;

        // Get a new SpaData instance
        let mut spa = SpaData::new(input);

        // call the SPA calculate function
        if let Err(e) = spa.spa_calculate() {
            panic!("{}", e);
        }

        // test and test results according original code
        check_test_data(&spa);
    }
}
