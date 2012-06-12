package potracej;

import java.util.Arrays;

public class privcurve_t {

    public int n;            /* number of segments */
    public curve_t.CurveTag tag[];         /* tag[n]: POTRACE_CORNER or POTRACE_CURVETO */
    public dpoint_t c[][]; /* c[n][i]: control points.
  		       c[n][0] is unused for tag[n]=POTRACE_CORNER */
    /* the remainder of this structure is special to privcurve, and is
       used in EPS debug output and special EPS "short coding". These
       fields are valid only if "alphacurve" is set. */
    public int alphacurve;   /* have the following fields been initialized? */
    public dpoint_t vertex[]; /* for POTRACE_CORNER, this equals c[1] */
    public double[] alpha;    /* only for POTRACE_CURVETO */
    public double[] alpha0;   /* "uncropped" alpha parameter - for debug output only */
    public double[] beta;

    public privcurve_t(int n) {
        this.n = n;
        this.tag = new curve_t.CurveTag[n];
        this.c = new dpoint_t[n][];
        this.vertex = new dpoint_t[n];
        for(int i=0; i<n; i++) {
            c[i] = new dpoint_t[3];
            for(int k=0; k<3; k++) {
                c[i][k] = new dpoint_t();
            }
            vertex[i] = new dpoint_t();
        }
        alpha = new double[n];
        alpha0 = new double[n];
        beta = new double[n];
    }

    @Override
    public String toString() {
        return "privcurve_t{" +
                "alpha0=" + alpha0 +
                ", n=" + n +
                ", tag=" + (tag == null ? null : Arrays.asList(tag)) +
                ", c=" + (c == null ? null : Arrays.asList(c)) +
                ", alphacurve=" + alphacurve +
                ", vertex=" + (vertex == null ? null : Arrays.asList(vertex)) +
                ", alpha=" + alpha +
                ", beta=" + beta +
                '}';
    }
}
