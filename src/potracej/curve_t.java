package potracej;

public class curve_t {
    public int n;                    /* number of segments */
    public CurveTag tag[];                 /* tag[n]: POTRACE_CURVETO or POTRACE_CORNER */
    public dpoint_t c[][]; /* c[n][3]: control points. */

    public enum CurveTag {
        POTRACE_CURVETO, POTRACE_CORNER
    }
}
