use spa_sra::spa::{Function, SpaData};
use spa_sra::SpaBuilder;

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

/// Checks (not equal) spa calculation output against original test results
///
fn check_test_data_ne<T: Clone>(spa: &SpaData<T>) {
    assert_ne!(format!("{:.6}", spa.spa_za.jd), "2452930.312847",            "Julian Day:    2452930.312847");
    assert_ne!(format!("{:.6e}", spa.spa_za.l), "2.401826e1",                "L:             2.401826e+01 degrees");
    assert_ne!(format!("{:.6e}", spa.spa_za.b), "-1.011219e-4",              "B:             -1.011219e-04 degrees");
    assert_ne!(format!("{:.6}", spa.spa_za.r), "0.996542",                   "R:             0.996542 AU");
    assert_ne!(format!("{:.6}", spa.spa_za.h), "11.105902",                  "H:             11.105902 degrees");
    assert_ne!(format!("{:.6e}", spa.spa_za.del_psi), "-3.998404e-3",        "Delta Psi:     -3.998404e-03 degrees");
    assert_ne!(format!("{:.6e}", spa.spa_za.del_epsilon), "1.666568e-3",     "Delta Epsilon: 1.666568e-03 degrees");
    assert_ne!(format!("{:.6}", spa.spa_za.epsilon), "23.440465",            "Epsilon:       23.440465 degrees");
    assert_ne!(format!("{:.6}", spa.spa_za.zenith), "50.111622",             "Zenith:        50.111622 degrees");
    assert_ne!(format!("{:.6}", spa.spa_za.azimuth), "194.340241",           "Azimuth:       194.340241 degrees");
    assert_ne!(format!("{:.6}", spa.spa_za_inc.incidence), "25.187000",      "Incidence:     25.187000 degrees");

    let mut min: f64 = 60.0 * (spa.spa_za_rts.sunrise - (spa.spa_za_rts.sunrise as i64) as f64);
    let mut sec: f64 = 60.0 * (min - (min as i64) as f64);
    assert_ne!(format!("{:0>2}:{:0>2}:{:0>2}", spa.spa_za_rts.sunrise as i64, min as i64, sec as i64), "06:12:43",
               "Sunrise: 06:12:43 Local Time");

    min = 60.0 * (spa.spa_za_rts.sunset - (spa.spa_za_rts.sunset as i64) as f64);
    sec = 60.0 * (min - (min as i64) as f64);
    assert_ne!(format!("{:0>2}:{:0>2}:{:0>2}", spa.spa_za_rts.sunset as i64, min as i64, sec as i64), "17:20:19",
               "Sunset: 17:20:19 Local Time");
}
#[test]
fn build_and_calculate() {
    let spa = SpaBuilder::new()
        .date(2003, 10, 17).unwrap_or_else(|e| panic!("{}", e))
        .time(12, 30, 30.0).unwrap_or_else(|e| panic!("{}", e))
        .timezone(-7.0).unwrap_or_else(|e| panic!("{}", e))
        .lat_long(39.742476, -105.1786).unwrap_or_else(|e| panic!("{}", e))
        .pressure(820.0).unwrap_or_else(|e| panic!("{}", e))
        .temperature(11.0).unwrap_or_else(|e| panic!("{}", e))
        .atmospheric_refraction(0.5667).unwrap_or_else(|e| panic!("{}", e))
        .elevation(1830.14).unwrap_or_else(|e| panic!("{}", e))
        .slope(30.0).unwrap_or_else(|e| panic!("{}", e))
        .azimuth_rotation(-10.0).unwrap_or_else(|e| panic!("{}", e))
        .delta_ut1(0.0).unwrap_or_else(|e| panic!("{}", e))
        .delta_t(67.0).unwrap_or_else(|e| panic!("{}", e))
        .execute(Function::SpaAll).unwrap_or_else(|e| panic!("{}", e));

    // test and test results according original code
    check_test_data(&spa);
}

#[test]
fn from_input_std() {
    let spa = SpaBuilder::new()
        .date(2025, 8, 12).unwrap_or_else(|e| panic!("{}", e))
        .time(12, 30, 30.0).unwrap_or_else(|e| panic!("{}", e))
        .timezone(-7.0).unwrap_or_else(|e| panic!("{}", e))
        .lat_long(39.742476, -105.1786).unwrap_or_else(|e| panic!("{}", e))
        .pressure(820.0).unwrap_or_else(|e| panic!("{}", e))
        .temperature(11.0).unwrap_or_else(|e| panic!("{}", e))
        .atmospheric_refraction(0.5667).unwrap_or_else(|e| panic!("{}", e))
        .elevation(1830.14).unwrap_or_else(|e| panic!("{}", e))
        .slope(30.0).unwrap_or_else(|e| panic!("{}", e))
        .azimuth_rotation(-10.0).unwrap_or_else(|e| panic!("{}", e))
        .delta_ut1(0.0).unwrap_or_else(|e| panic!("{}", e))
        .delta_t(67.0).unwrap_or_else(|e| panic!("{}", e))
        .execute(Function::SpaAll).unwrap_or_else(|e| panic!("{}", e));

    // Should not result in original test results
    check_test_data_ne(&spa);

    let spa = SpaBuilder::from_input(spa.input)
        .date(2003, 10, 17).unwrap_or_else(|e| panic!("{}", e))
        .execute(Function::SpaAll).unwrap_or_else(|e| panic!("{}", e));

    // test and test results according original code
    check_test_data(&spa);
}

#[cfg(feature = "chrono_0_4")]
#[test]
fn feature_chrono_0_4() {
    use chrono::{FixedOffset, TimeZone, Timelike};

    let date_time = FixedOffset::east_opt(-7 * 3600)
        .unwrap()
        .with_ymd_and_hms(2003, 10, 17, 12, 30, 30)
        .unwrap();

    let spa = SpaBuilder::from_date_time(date_time)
        .lat_long(39.742476, -105.1786).unwrap_or_else(|e| panic!("{}", e))
        .pressure(820.0).unwrap_or_else(|e| panic!("{}", e))
        .temperature(11.0).unwrap_or_else(|e| panic!("{}", e))
        .atmospheric_refraction(0.5667).unwrap_or_else(|e| panic!("{}", e))
        .elevation(1830.14).unwrap_or_else(|e| panic!("{}", e))
        .slope(30.0).unwrap_or_else(|e| panic!("{}", e))
        .azimuth_rotation(-10.0).unwrap_or_else(|e| panic!("{}", e))
        .delta_ut1(0.0).unwrap_or_else(|e| panic!("{}", e))
        .delta_t(67.0).unwrap_or_else(|e| panic!("{}", e))
        .execute(Function::SpaAll).unwrap_or_else(|e| panic!("{}", e));

    // test and test results according original code
    check_test_data(&spa);

    let was = spa.get_sunrise();
    let should = FixedOffset::east_opt(-7 * 3600)
        .unwrap()
        .with_ymd_and_hms(2003, 10, 17, 06, 12, 43)
        .unwrap()
        .with_nanosecond(439793424)
        .unwrap();

    assert_eq!(was, should);
}

#[cfg(feature = "chrono_0_4")]
#[test]
fn from_input_tz() {
    use chrono::{FixedOffset, TimeZone};

    let date_time = FixedOffset::east_opt(-7 * 3600)
        .unwrap()
        .with_ymd_and_hms(2025, 8, 12, 12, 30, 30)
        .unwrap();

    let spa = SpaBuilder::from_date_time(date_time)
        .lat_long(39.742476, -105.1786).unwrap_or_else(|e| panic!("{}", e))
        .pressure(820.0).unwrap_or_else(|e| panic!("{}", e))
        .temperature(11.0).unwrap_or_else(|e| panic!("{}", e))
        .atmospheric_refraction(0.5667).unwrap_or_else(|e| panic!("{}", e))
        .elevation(1830.14).unwrap_or_else(|e| panic!("{}", e))
        .slope(30.0).unwrap_or_else(|e| panic!("{}", e))
        .azimuth_rotation(-10.0).unwrap_or_else(|e| panic!("{}", e))
        .delta_ut1(0.0).unwrap_or_else(|e| panic!("{}", e))
        .delta_t(67.0).unwrap_or_else(|e| panic!("{}", e))
        .execute(Function::SpaAll).unwrap_or_else(|e| panic!("{}", e));

    // Should not result in original test results
    check_test_data_ne(&spa);

    let date_time = FixedOffset::east_opt(-7 * 3600)
        .unwrap()
        .with_ymd_and_hms(2003, 10, 17, 12, 30, 30)
        .unwrap();

    let spa = SpaBuilder::from_input(spa.input)
        .date_time(date_time)
        .execute(Function::SpaAll).unwrap_or_else(|e| panic!("{}", e));

    // test and test results according original code
    check_test_data(&spa);
}

#[cfg(feature = "chrono_0_4")]
#[test]
fn hybrid_chrono_0_4() {
    use chrono::{FixedOffset, TimeZone};

    let date_time = FixedOffset::east_opt(2 * 3600)
        .unwrap()
        .with_ymd_and_hms(2025, 8, 12, 12, 30, 30)
        .unwrap();

    let mut spa = SpaBuilder::from_date_time(date_time)
        .lat_long(39.742476, -105.1786).unwrap_or_else(|e| panic!("{}", e))
        .pressure(820.0).unwrap_or_else(|e| panic!("{}", e))
        .temperature(11.0).unwrap_or_else(|e| panic!("{}", e))
        .atmospheric_refraction(0.5667).unwrap_or_else(|e| panic!("{}", e))
        .elevation(1830.14).unwrap_or_else(|e| panic!("{}", e))
        .slope(30.0).unwrap_or_else(|e| panic!("{}", e))
        .azimuth_rotation(-10.0).unwrap_or_else(|e| panic!("{}", e))
        .delta_ut1(0.0).unwrap_or_else(|e| panic!("{}", e))
        .delta_t(67.0).unwrap_or_else(|e| panic!("{}", e))
        .execute(Function::SpaAll).unwrap_or_else(|e| panic!("{}", e));

    // Should not result in original test results
    check_test_data_ne(&spa);

    // Create date time according original test
    let date_time = FixedOffset::east_opt(-7 * 3600)
        .unwrap()
        .with_ymd_and_hms(2003, 10, 17, 12, 30, 30)
        .unwrap();

    // Set date, time and timezone (numerical for the algorithm)
    spa.input.date_time(date_time);

    // Execute algorithm directly from the SpaData instance
    if let Err(e) = spa.spa_calculate() {
        panic!("{}", e);
    }

    // test and test results according original code
    check_test_data(&spa);
}
