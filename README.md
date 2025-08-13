# SPA_SRA

The SPA_SRA (Solar Position Algorithm for Solar Radiation Applications) calculates the solar zenith and
azimuth angles in the period from the year -2000 to 6000, with uncertainties of +/- 0.0003 degrees
based on the date, time, and location on Earth. (Reference: Reda, I.; Andreas, A., Solar Position Algorithm
for Solar Radiation Applications, Solar Energy. Vol. 76(5), 2004; pp. 577-589).

It can also calculate the surface incidence angle for e.g. a solar panel.
The surface incidence angle is the angle between an incoming ray (like light or radar) and
a line perpendicular to the surface at the point where the ray hits.

Further information on this algorithm is available in the following NREL technical report (pdf):
[Reda, I.; Andreas, A. (2003). Solar Position Algorithm for Solar Radiation Applications. 55 pp.; NREL Report No. TP-560-34302, Revised January 2008.](http://www.nrel.gov/docs/fy08osti/34302.pdf)

## License and Acknowledgement

SPA_SRA is a port with rust specific adjustments from the original C-code prepared by employees
of the Alliance for Sustainable Energy, LLC. Hence, the [MIT license](https://github.com/gostonefire/spa_sra/blob/main/LICENSE) comes with an acknowledgement
and disclaimer related to that original code. For personal use it should be fine though.
