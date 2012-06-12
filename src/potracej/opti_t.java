package potracej;

public class opti_t {
    double pen;	   /* penalty */
    dpoint_t c[];   /* curve parameters */
    double t, s;	   /* curve parameters */
    double alpha;	   /* curve parameter */

    public opti_t() {
        c = new dpoint_t[2];
        c[0] = new dpoint_t();
        c[1] = new dpoint_t();
    }
}
