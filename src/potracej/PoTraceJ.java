package potracej;

import java.text.DecimalFormat;
import java.util.ArrayList;

import static potracej.Bitmap.BM_GET;

/**
 * User: Sanny Sanoff (san@sysdate.com)
 * Date: 6/12/12
 *
 *
 */
public class PoTraceJ {
    param_t param;

    public PoTraceJ(param_t param) {
        this.param = param;
    }

    public path_t trace(Bitmap bitmap) {
        long l = System.currentTimeMillis();
        path_t paths = bm_to_pathlist(bitmap);
        l = System.currentTimeMillis() - l;
        //System.out.println("Decompose: "+l+" msec");
        l = System.currentTimeMillis();
        process_path(paths);
        l = System.currentTimeMillis() - l;
        //System.out.println("Curves: "+l+" msec");
        return paths;
    }

    public static String f(double d) {
        DecimalFormat decimalFormat = new DecimalFormat();
        decimalFormat.setMaximumFractionDigits(6);
        decimalFormat.setMinimumFractionDigits(6);
        return decimalFormat.format(d);
    }

    long times[] = new long[10];

    public void resetTimers() {
        for (int i = 0; i < times.length; i++) {
            times[i] = 0;
        }
    }

    private void process_path(path_t paths) {
        path_t p = paths;
        long ol;
        while (p != null) {
            long l = System.currentTimeMillis();
            calc_sums(p.priv);
            ol = l; l = System.currentTimeMillis();
            times[0] += ol - l;

            calc_lon(p.priv);
            ol = l; l = System.currentTimeMillis();
            times[1] += ol - l;

            bestpolygon(p.priv);
            ol = l; l = System.currentTimeMillis();
            times[2] += ol - l;

            adjust_vertices(p.priv);
            ol = l; l = System.currentTimeMillis();
            times[3] += ol - l;

            if (p.sign == '-') {
                reverse(p.priv.curve);
            }
            ol = l; l = System.currentTimeMillis();
            times[4] += ol - l;
//            privcurve_t curve = p.priv.curve;
//            for(int i=0; i< curve.n; i++) {
//                System.out.println(i+": "+curve.vertex[i].x+" "+curve.vertex[i].y);
//            }

            smooth(p.priv.curve, param.alphamax);
            ol = l; l = System.currentTimeMillis();
            times[5] += ol - l;

            if (param.opticurve != 0) {
                opticurve(p.priv, param.opttolerance);
                p.priv.fcurve = p.priv.ocurve;
            } else {
                p.priv.fcurve = p.priv.curve;
            }
            smooth(p.priv.curve, param.alphamax);
            ol = l; l = System.currentTimeMillis();
            times[6] += ol - l;

//            for(int i=0; i< curve.n; i++) {
//                System.out.println(i+": "+f(curve.c[i][0].x)+","+f(curve.c[i][0].y)+" "+f(curve.c[i][1].x)+","+f(curve.c[i][1].y)+" " + f(curve.c[i][2].x)+","+f(curve.c[i][2].y));
//                System.out.println(i+":                 v="+f(curve.vertex[i].x)+","+f(curve.vertex[i].y)+" "+f(curve.alpha[i])+" "+f(curve.beta[i])+" "+f(curve.alpha0[i])+" "+curve.tag[i]+" ");
//            }
            p.curve = new curve_t();
            privcurve_to_curve(p.priv.fcurve, p.curve);
            smooth(p.priv.curve, param.alphamax);
            ol = l; l = System.currentTimeMillis();
            times[7] += ol - l;

            p = p.next;
        }
    }

    public void printTimers() {
        for (int i = 0; i < times.length; i++) {
            long time = times[i];
            System.out.println("Timer "+i+" = "+(-time));
        }
    }

    private void privcurve_to_curve(privcurve_t pc, curve_t c) {
        c.n = pc.n;
        c.tag = pc.tag;
        c.c = pc.c;
    }

    private void opticurve(privpath_t pp, double opttolerance) {
        int m = pp.curve.n;
        int[] pt = new int[m + 1];     /* pt[m+1] */
        double[] pen = new double[m + 1]; /* pen[m+1] */
        int[] len = new int[m + 1];    /* len[m+1] */
        opti_t[] opt = new opti_t[m + 1]; /* opt[m+1] */
        for (int i = 0; i < opt.length; i++) {
            opt[i] = new opti_t();
        }
        int om;
        int i, j, r;
        dpoint_t p0;
        int i1;
        double area;
        double alpha;
        double[] s;
        double[] t;

        int[] convc = new int[m]; /* conv[m]: pre-computed convexities */
        double[] areac = new double[m + 1]; /* cumarea[m+1]: cache for fast area computation */

        /* pre-calculate convexity: +1 = right turn, -1 = left turn, 0 = corner */
        for (i = 0; i < m; i++) {
            if (pp.curve.tag[i] == curve_t.CurveTag.POTRACE_CURVETO) {
                convc[i] = sign(dpara(pp.curve.vertex[mod(i - 1, m)], pp.curve.vertex[i], pp.curve.vertex[mod(i + 1, m)]));
            } else {
                convc[i] = 0;
            }
        }

        /* pre-calculate areas */
        area = 0.0;
        areac[0] = 0.0;
        p0 = pp.curve.vertex[0];
        for (i = 0; i < m; i++) {
            i1 = mod(i + 1, m);
            if (pp.curve.tag[i1] == curve_t.CurveTag.POTRACE_CURVETO) {
                alpha = pp.curve.alpha[i1];
                area += 0.3 * alpha * (4 - alpha) * dpara(pp.curve.c[i][2], pp.curve.vertex[i1], pp.curve.c[i1][2]) / 2;
                area += dpara(p0, pp.curve.c[i][2], pp.curve.c[i1][2]) / 2;
            }
            areac[i + 1] = area;
        }

        pt[0] = -1;
        pen[0] = 0;
        len[0] = 0;

        /* Fixme: we always start from a fixed point -- should find the best
           curve cyclically */

        for (j = 1; j <= m; j++) {
            /* calculate best path from 0 to j */
            pt[j] = j - 1;
            pen[j] = pen[j - 1];
            len[j] = len[j - 1] + 1;

            for (i = j - 2; i >= 0; i--) {
                opti_t o = new opti_t();
                r = opti_penalty(pp, i, mod(j, m), o, opttolerance, convc, areac);
                if (r > 0) break;
                if (len[j] > len[i] + 1 || (len[j] == len[i] + 1 && pen[j] > pen[i] + o.pen)) {
                    pt[j] = i;
                    pen[j] = pen[i] + o.pen;
                    len[j] = len[i] + 1;
                    opt[j] = o;
                }
            }
        }
        om = len[m];
        pp.ocurve = new privcurve_t(om);
        s = new double[om];
        t = new double[om];

        j = m;
        for (i = om - 1; i >= 0; i--) {
            if (pt[j] == j - 1) {
                pp.ocurve.tag[i] = pp.curve.tag[mod(j, m)];
                pp.ocurve.c[i][0] = pp.curve.c[mod(j, m)][0];
                pp.ocurve.c[i][1] = pp.curve.c[mod(j, m)][1];
                pp.ocurve.c[i][2] = pp.curve.c[mod(j, m)][2];
                pp.ocurve.vertex[i] = pp.curve.vertex[mod(j, m)];
                pp.ocurve.alpha[i] = pp.curve.alpha[mod(j, m)];
                pp.ocurve.alpha0[i] = pp.curve.alpha0[mod(j, m)];
                pp.ocurve.beta[i] = pp.curve.beta[mod(j, m)];
                s[i] = t[i] = 1.0;
            } else {
                pp.ocurve.tag[i] = curve_t.CurveTag.POTRACE_CURVETO;
                pp.ocurve.c[i][0] = opt[j].c[0];
                pp.ocurve.c[i][1] = opt[j].c[1];
                pp.ocurve.c[i][2] = pp.curve.c[mod(j, m)][2];
                pp.ocurve.vertex[i] = interval(opt[j].s, pp.curve.c[mod(j, m)][2], pp.curve.vertex[mod(j, m)]);
                pp.ocurve.alpha[i] = opt[j].alpha;
                pp.ocurve.alpha0[i] = opt[j].alpha;
                s[i] = opt[j].s;
                t[i] = opt[j].t;
            }
            j = pt[j];
        }

        /* calculate beta parameters */
        for (i = 0; i < om; i++) {
            i1 = mod(i + 1, om);
            pp.ocurve.beta[i] = s[i] / (s[i] + t[i1]);
        }
        pp.ocurve.alphacurve = 1;

    }

    private void smooth(privcurve_t curve, double alphamax) {
        int m = curve.n;

        int i, j, k;
        double dd, denom, alpha;
        dpoint_t p2, p3, p4;

        /* examine each vertex and find its best fit */
        for (i = 0; i < m; i++) {
            j = mod(i + 1, m);
            k = mod(i + 2, m);
            p4 = interval(1 / 2.0, curve.vertex[k], curve.vertex[j]);

            denom = ddenom(curve.vertex[i], curve.vertex[k]);
            if (denom != 0.0) {
                dd = dpara(curve.vertex[i], curve.vertex[j], curve.vertex[k]) / denom;
                dd = Math.abs(dd);
                alpha = dd > 1 ? (1 - 1.0 / dd) : 0;
                alpha = alpha / 0.75;
            } else {
                alpha = 4 / 3.0;
            }
            curve.alpha0[j] = alpha;     /* remember "original" value of alpha */

            if (alpha > alphamax) {  /* pointed corner */
                curve.tag[j] = curve_t.CurveTag.POTRACE_CORNER;
                curve.c[j][1] = curve.vertex[j];
                curve.c[j][2] = p4;
            } else {
                if (alpha < 0.55) {
                    alpha = 0.55;
                } else if (alpha > 1) {
                    alpha = 1;
                }
                p2 = interval(.5 + .5 * alpha, curve.vertex[i], curve.vertex[j]);
                p3 = interval(.5 + .5 * alpha, curve.vertex[k], curve.vertex[j]);
                curve.tag[j] = curve_t.CurveTag.POTRACE_CURVETO;
                curve.c[j][0] = p2;
                curve.c[j][1] = p3;
                curve.c[j][2] = p4;
            }
            curve.alpha[j] = alpha;    /* store the "cropped" value of alpha */
            curve.beta[j] = 0.5;
        }
        curve.alphacurve = 1;
        return;
    }

    private void reverse(privcurve_t curve) {
        int m = curve.n;
        int i, j;
        dpoint_t tmp;

        for (i = 0, j = m - 1; i < j; i++, j--) {
            tmp = curve.vertex[i];
            curve.vertex[i] = curve.vertex[j];
            curve.vertex[j] = tmp;
        }
    }

    private void adjust_vertices(privpath_t pp) {
        int m = pp.m;
        int[] po = pp.po;
        int n = pp.pt.size();
        point_t pt[] = pp.pt.toArray(new point_t[n]);
        int x0 = pp.x0;
        int y0 = pp.y0;

        dpoint_t[] ctr = new dpoint_t[m];      /* ctr[m] */
        dpoint_t[] dir = new dpoint_t[m];      /* dir[m] */
        double q[][][] = new double[m][][];      /* q[m] */
        for (int i = 0; i < m; i++) {
            q[i] = new_quadform();
            ctr[i] = new dpoint_t();
            dir[i] = new dpoint_t();
        }
        double v[] = new double[3];
        double d;
        int i, j, k, l;
        dpoint_t s = new dpoint_t();
        int r;

        pp.curve = new privcurve_t(m);

        /* calculate "optimal" point-slope representation for each line
           segment */
        for (i = 0; i < m; i++) {
            j = po[mod(i + 1, m)];
            j = mod(j - po[i], n) + po[i];
            pointslope(pp, po[i], j, ctr[i], dir[i]);
        }

        /* represent each line segment as a singular quadratic form; the
           distance of a point (x,y) from the line segment will be
           (x,y,1)Q(x,y,1)^t, where Q=q[i]. */
        for (i = 0; i < m; i++) {
            d = sq(dir[i].x) + sq(dir[i].y);
            if (d == 0.0) {
                for (j = 0; j < 3; j++) {
                    for (k = 0; k < 3; k++) {
                        q[i][j][k] = 0;
                    }
                }
            } else {
                v[0] = dir[i].y;
                v[1] = -dir[i].x;
                v[2] = -v[1] * ctr[i].y - v[0] * ctr[i].x;
                for (l = 0; l < 3; l++) {
                    for (k = 0; k < 3; k++) {
                        q[i][l][k] = v[l] * v[k] / d;
                    }
                }
            }
        }

        /* now calculate the "intersections" of consecutive segments.
           Instead of using the actual intersection, we find the point
           within a given unit square which minimizes the square distance to
           the two lines. */
        for (i = 0; i < m; i++) {
            double[][] Q = new_quadform();
            dpoint_t w = new dpoint_t();
            double dx, dy;
            double det;
            double min, cand; /* minimum and candidate for minimum of quad. form */
            double xmin, ymin;    /* coordinates of minimum */
            int z;

            /* let s be the vertex, in coordinates relative to x0/y0 */
            s.x = pt[po[i]].x - x0;
            s.y = pt[po[i]].y - y0;

            /* intersect segments i-1 and i */

            j = mod(i - 1, m);

            /* add quadratic forms */
            for (l = 0; l < 3; l++) {
                for (k = 0; k < 3; k++) {
                    Q[l][k] = q[j][l][k] + q[i][l][k];
                }
            }

            while (true) {
                /* minimize the quadratic form Q on the unit square */
                /* find intersection */

                det = Q[0][0] * Q[1][1] - Q[0][1] * Q[1][0];
                if (det != 0.0) {
                    w.x = (-Q[0][2] * Q[1][1] + Q[1][2] * Q[0][1]) / det;
                    w.y = (Q[0][2] * Q[1][0] - Q[1][2] * Q[0][0]) / det;
                    break;
                }

                /* matrix is singular - lines are parallel. Add another,
               orthogonal axis, through the center of the unit square */
                if (Q[0][0] > Q[1][1]) {
                    v[0] = -Q[0][1];
                    v[1] = Q[0][0];
                } else if (Q[1][1] != 0) {
                    v[0] = -Q[1][1];
                    v[1] = Q[1][0];
                } else {
                    v[0] = 1;
                    v[1] = 0;
                }
                d = sq(v[0]) + sq(v[1]);
                v[2] = -v[1] * s.y - v[0] * s.x;
                for (l = 0; l < 3; l++) {
                    for (k = 0; k < 3; k++) {
                        Q[l][k] += v[l] * v[k] / d;
                    }
                }
            }
            dx = Math.abs(w.x - s.x);
            dy = Math.abs(w.y - s.y);
            if (dx <= .5 && dy <= .5) {
                pp.curve.vertex[i].x = w.x + x0;
                pp.curve.vertex[i].y = w.y + y0;
                continue;
            }

            /* the minimum was not in the unit square; now minimize quadratic
        on boundary of square */
            min = quadform(Q, s);
            xmin = s.x;
            ymin = s.y;

            if (Q[0][0] == 0.0) {
                //
            } else {
                for (z = 0; z < 2; z++) {   /* value of the y-coordinate */
                    w.y = s.y - 0.5 + z;
                    w.x = -(Q[0][1] * w.y + Q[0][2]) / Q[0][0];
                    dx = Math.abs(w.x - s.x);
                    cand = quadform(Q, w);
                    if (dx <= .5 && cand < min) {
                        min = cand;
                        xmin = w.x;
                        ymin = w.y;
                    }
                }
            }
            if (Q[1][1] == 0.0) {
                //
            } else {
                for (z = 0; z < 2; z++) {   /* value of the x-coordinate */
                    w.x = s.x - 0.5 + z;
                    w.y = -(Q[1][0] * w.x + Q[1][2]) / Q[1][1];
                    dy = Math.abs(w.y - s.y);
                    cand = quadform(Q, w);
                    if (dy <= .5 && cand < min) {
                        min = cand;
                        xmin = w.x;
                        ymin = w.y;
                    }
                }
            }
            /* check four corners */
            for (l = 0; l < 2; l++) {
                for (k = 0; k < 2; k++) {
                    w.x = s.x - 0.5 + l;
                    w.y = s.y - 0.5 + k;
                    cand = quadform(Q, w);
                    if (cand < min) {
                        min = cand;
                        xmin = w.x;
                        ymin = w.y;
                    }
                }
            }

            pp.curve.vertex[i].x = xmin + x0;
            pp.curve.vertex[i].y = ymin + y0;
            continue;
        }

//        for(i=0; i<m; i++) {
//            System.out.println("vertex["+i+"]="+pp.curve.vertex[i]);
//        }

    }

    private double[][] new_quadform() {
        double[][] retval = new double[3][];
        for (int i = 0; i < retval.length; i++) {
            retval[i] = new double[3];
        }
        return retval;
    }

    private void bestpolygon(privpath_t pp) {
        int i, j, m, k;
        int n = pp.pt.size();
        double[] pen = new double[n + 1]; /* pen[n+1]: penalty vector */
        int[] prev = new int[n + 1];   /* prev[n+1]: best path pointer vector */
        int[] clip0 = new int[n];  /* clip0[n]: longest segment pointer, non-cyclic */
        int[] clip1 = new int[n + 1];  /* clip1[n+1]: backwards segment pointer, non-cyclic */
        int[] seg0 = new int[n + 1];    /* seg0[m+1]: forward segment bounds, m<=n */
        int[] seg1 = new int[n + 1];   /* seg1[m+1]: backward segment bounds, m<=n */
        double thispen;
        double best;
        int c;

        /* calculate clipped paths */
        for (i = 0; i < n; i++) {
            c = mod(pp.lon[mod(i - 1, n)] - 1, n);
            if (c == i) {
                c = mod(i + 1, n);
            }
            if (c < i) {
                clip0[i] = n;
            } else {
                clip0[i] = c;
            }
        }

        /* calculate backwards path clipping, non-cyclic. j <= clip0[i] iff
           clip1[j] <= i, for i,j=0..n. */
        j = 1;
        for (i = 0; i < n; i++) {
            while (j <= clip0[i]) {
                clip1[j] = i;
                j++;
            }
        }

        /* calculate seg0[j] = longest path from 0 with j segments */
        i = 0;
        for (j = 0; i < n; j++) {
            seg0[j] = i;
            i = clip0[i];
        }
        seg0[j] = n;
        m = j;

        /* calculate seg1[j] = longest path to n with m-j segments */
        i = n;
        for (j = m; j > 0; j--) {
            seg1[j] = i;
            i = clip1[i];
        }
        seg1[0] = 0;

        /* now find the shortest path with m segments, based on penalty3 */
        /* note: the outer 2 loops jointly have at most n iterations, thus
           the worst-case behavior here is quadratic. In practice, it is
           close to linear since the inner loop tends to be short. */
        pen[0] = 0;
        for (j = 1; j <= m; j++) {
            for (i = seg1[j]; i <= seg0[j]; i++) {
                best = -1;
                int limit = clip1[i];
                for (k = seg0[j - 1]; k >= limit; k--) {
                    thispen = penalty3(pp, k, i) + pen[k];
                    if (best < 0 || thispen < best) {
                        prev[i] = k;
                        best = thispen;
                    }
                }
                pen[i] = best;
            }
        }

        pp.m = m;
        pp.po = new int[m];

        /* read off shortest path */
        for (i = n, j = m - 1; i > 0; j--) {
            i = prev[i];
            pp.po[j] = i;
        }

//        for(i=0; i<m; i++) {
//            System.out.println("po["+i+"]="+pp.po[i]);
//        }
    }

    public static final int INFTY = 10000000;    /* it suffices that this is longer than any
			                                    path; it need not be really infinite */

    public static double COS179 = -0.999847695156;	 /* the cosine of 179 degrees */

    private void calc_lon(privpath_t pp) {
        int n = pp.pt.size();
        point_t pt[] = pp.getPTArray();
        int i, j, k, k1;
        //int[] ct = new int[4];
        int dir;
        int iter = 0;
        point_t constraint0 = new point_t();
        point_t constraint1 = new point_t();
        point_t cur = new point_t();
        point_t off = new point_t();
        int[] pivk = new int[n];  /* pivk[n] */
        int[] nc = new int[n];    /* nc[n]: next corner */
        point_t dk = new point_t();  /* direction of k-k1 */
        int a, b, c, d;

        /* initialize the nc data structure. Point from each point to the
           furthest future point to which it is connected by a vertical or
           horizontal segment. We take advantage of the fact that there is
           always a direction change at 0 (due to the path decomposition
           algorithm). But even if this were not so, there is no harm, as
           in practice, correctness does not depend on the word "furthest"
           above.  */
        k = 0;
        for (i = n - 1; i >= 0; i--) {
            if (pt[i].x != pt[k].x && pt[i].y != pt[k].y) {
                k = i + 1;  /* necessarily i<n-1 in this case */
            }
            nc[i] = k;
            //System.out.println("NC["+i+"]="+k);
        }

        pp.lon = new int[n];

        /* determine pivot points: for each i, let pivk[i] be the furthest k
           such that all j with i<j<k lie on a line connecting i,k. */
        int ct;
        foundkloop:
        for (i = n - 1; i >= 0; i--) {
            ct = 0; // ct[0] = ct[1] = ct[2] = ct[3] = 0;

            /* keep track of "directions" that have occurred */
            dir = (3 + 3 * (pt[mod(i + 1, n)].x - pt[i].x) + (pt[mod(i + 1, n)].y - pt[i].y)) / 2;
            ct |= (1 << dir);
            // ct[dir]++;

            constraint0.x = 0;
            constraint0.y = 0;
            constraint1.x = 0;
            constraint1.y = 0;

            /* find the next k such that no straight line from i to k */
            k = nc[i];
            k1 = i;
constvioloop:
            while (true) {
                //iter++;
                int ptkx = pt[k].x;
                int ptky = pt[k].y;
                dir = (3 + 3 * sign(ptkx - pt[k1].x) + sign(ptky - pt[k1].y)) / 2;
                ct |= (1 << dir);
                //ct[dir]++;

                /* if all four "directions" have occurred, cut this path */
                //if (ct[0] != 0 && ct[1] != 0 && ct[2] != 0 && ct[3] != 0) {
                if (ct == 0xF) {
                    pivk[i] = k1;
                    continue foundkloop;
                }

                int curx = cur.x = ptkx - pt[i].x;
                int cury = cur.y = ptky - pt[i].y;

                /* see if current constraint is violated */
                if (xprod(constraint0, cur) < 0 || xprod(constraint1, cur) > 0) {
                    //System.out.println("o!");
                    break constvioloop;
                }

                /* else, update constraint */
                if (Math.abs(curx) <= 1 && Math.abs(cury) <= 1) {
                    /* no constraint */
                } else {
                    off.x = curx + ((cury >= 0 && (cury > 0 || curx < 0)) ? 1 : -1);
                    off.y = cury + ((curx <= 0 && (curx < 0 || cury < 0)) ? 1 : -1);
                    if (xprod(constraint0, off) >= 0) {
                        constraint0.x = off.x;
                        constraint0.y = off.y;
                    }
                    off.x = curx + ((cury <= 0 && (cury < 0 || curx < 0)) ? 1 : -1);
                    off.y = cury + ((curx >= 0 && (curx > 0 || cury < 0)) ? 1 : -1);
                    if (xprod(constraint1, off) <= 0) {
                        constraint1.x = off.x;
                        constraint1.y = off.y;
                    }
                }
                k1 = k;
                k = nc[k1];
                if (!cyclic(k, i, k1)) {
                    break;
                }
            }
            /* k1 was the last "corner" satisfying the current constraint, and
     k is the first one violating it. We now need to find the last
     point along k1..k which satisfied the constraint. */
            dk.x = sign(pt[k].x - pt[k1].x);
            dk.y = sign(pt[k].y - pt[k1].y);
            cur.x = pt[k1].x - pt[i].x;
            cur.y = pt[k1].y - pt[i].y;
            /* find largest integer j such that xprod(constraint[0], cur+j*dk)
          >= 0 and xprod(constraint[1], cur+j*dk) <= 0. Use bilinearity
          of xprod. */
            a = xprod(constraint0, cur);
            b = xprod(constraint0, dk);
            c = xprod(constraint1, cur);
            d = xprod(constraint1, dk);
            /* find largest integer j such that a+j*b>=0 and c+j*d<=0. This
      can be solved with integer arithmetic. */
            j = INFTY;
            if (b < 0) {
                j = floordiv(a, -b);
            }
            if (d > 0) {
                j = Math.min(j, floordiv(-c, d));
            }
            pivk[i] = mod(k1 + j, n);
            foundk:
            ;
        } /* for i */

        /* clean up: for each i, let lon[i] be the largest k such that for
           all i' with i<=i'<k, i'<k<=pivk[i']. */

        j = pivk[n - 1];
        pp.lon[n - 1] = j;
        for (i = n - 2; i >= 0; i--) {
            if (cyclic(i + 1, pivk[i], j)) {
                j = pivk[i];
            }
            pp.lon[i] = j;
        }
        for (i = n - 1; cyclic(mod(i + 1, n), j, pp.lon[i]); i--) {
            pp.lon[i] = j;
        }
    }

    private void calc_sums(privpath_t pp) {
        int i, x, y;
        int n = pp.pt.size();
        pp.sums = new sums_t[n + 1];
        for(i=0; i<n+1; i++) {
            pp.sums[i] = new sums_t();
        }

        /* origin */
        pp.x0 = pp.pt.get(0).x;
        pp.y0 = pp.pt.get(0).y;

        /* preparatory computation for later fast summing */
        pp.sums[0].x2 = pp.sums[0].xy = pp.sums[0].y2 = pp.sums[0].x = pp.sums[0].y = 0;
        for (i = 0; i < n; i++) {
            x = pp.pt.get(i).x - pp.x0;
            y = pp.pt.get(i).y - pp.y0;
            pp.sums[i + 1].x = pp.sums[i].x + x;
            pp.sums[i + 1].y = pp.sums[i].y + y;
            pp.sums[i + 1].x2 = pp.sums[i].x2 + x * x;
            pp.sums[i + 1].xy = pp.sums[i].xy + x * y;
            pp.sums[i + 1].y2 = pp.sums[i].y2 + y * y;
        }

//        for (i = 0; i < n; i++) {
//            System.out.println("Sums["+i+"]: "+pp.sums[i].x+" "+pp.sums[i].y+" "+pp.sums[i].x2+" "+pp.sums[i].xy+" "+pp.sums[i].y2);
//        }


    }

    private path_t bm_to_pathlist(Bitmap bm) {
        ArrayList<path_t> paths = new ArrayList<path_t>();
        //path_t plist = null;
        Bitmap bm1 = bm.dup();
        bm1.clearexcess();
        point_t pt = new point_t(0, bm1.h - 1);
        int count = 0;
        //System.out.println(bm1.toDebugString());
        while (findnext(bm1, pt) == 0) {
            char sign = bm.get(pt.x, pt.y) ? '+' : '-';
            path_t p = findpath(bm1, pt.x, pt.y + 1, sign, param.turnPolicy);
            xor_path(bm1, p);
            //System.out.println("PATH: area="+p.area);
            if (p.area <= param.turdsize) {
                // nothing
            } else {
                count++;
                for(int i=0; i<p.priv.pt.size(); i++) {
                    //System.out.println("point "+i+"="+p.priv.pt.get(i));
                }
                //System.out.println("===============================================================================");
                paths.add(p);
            }
        }
//        System.out.println("count="+count);
        return pathlist_to_tree(paths, bm1);
    }

    void print_path(ArrayList<path_t> paths) {
        int sum = 0;
        for (int i = 0; i < paths.size(); i++) {
            path_t p = paths.get(i);
            int n1 = p.area;
            int n2 = p.next != null ? p.next.area : 0;
            int n3 = (p.sibling != null ? p.sibling.area : 0);
            int n4 = p.childlist != null ? p.childlist.area : 0;
            System.out.print("P"+i+"="+n1+" n="+n2
                +" s="+n3
                +" cl="+n4
            );
            sum = ((sum + n1 + n2 + n3 + n4) << 1)| ((sum >> 31) & 1);
        }
        //System.out.println("\njsum="+sum);
    }

    private path_t pathlist_to_tree(ArrayList<path_t> paths, Bitmap bm) {
        int iter = 0;
        bm.clear();
        path_t p;
        if (paths.size() == 0) return null;
        for (int i = 0; i < paths.size(); i++) {
            p = paths.get(i);
            p.sibling = i < paths.size() - 1 ? paths.get(i+1) : null;
            p.next = p.sibling;
            p.childlist = null;
        }
        path_t heap = paths.get(0);
        path_t cur = null;
        path_t head = null;
        while (heap != null) {
            //System.out.println("step0: ............. cur="+(cur != null ? cur.area : 0)+" heap="+(heap != null ? heap.area : 0) + " head="+(head != null ? head.area : 0));
            /* unlink first sublist */
            cur = heap;
            heap = heap.childlist;
            cur.childlist = null;

            /* unlink first path */
            head = cur;
            cur = cur.next;
            head.next = null;

            //System.out.println("step0a: ............. cur="+(cur != null ? cur.area : 0)+" heap="+(heap != null ? heap.area : 0) + " head="+(head != null ? head.area : 0));
            xor_path(bm, head);
            bbox_t bbox = setbbox_path(head);


            path_t_Holder nexts = new path_t_Holder(null);
            path_t_Holder children = new path_t_Holder(null);

            /* now do insideness test for each element of cur; append it to
           head->childlist if it's inside head, else append it to
           head->next. */
            // hook_in=&head->childlist;
            //System.out.println("step1, begin loop....");
            p = cur;
            while (true) {
                if (p != null) {
                    cur = p.next;
                    p.next = null;
                } else break;
                //
                //System.out.println("step1, iter="+(iter++));
                //print_path(paths);
                if (p.priv.pt.get(0).y <= bbox.y0) {
                    nexts.addLast(p);
                    path_t headNextLast = find_last(head.next);
                    p.next = null;
                    if (headNextLast != null) {
                        headNextLast.next = p;
                    } else {
                        head.next = p;
                    }
                    //System.out.println("step1, E0iter="+(iter++));
                    //print_path(paths);
                    find_last(head).next = cur;
                    // ext = cur;
                    //p.next = cur;
//                    nexts.addLast(cur);
                    //System.out.println("step1, E0iter="+(iter++));
                    //print_path(paths);
                    break;
                }
                if (BM_GET(bm, p.priv.pt.get(0).x, p.priv.pt.get(0).y - 1)) {
//                    children.addLast(p);
                    path_t headChildListLast = find_last(head.childlist);
                    p.next = null;
                    if (headChildListLast != null) {
                        headChildListLast.next = p;
                    } else {
                        head.childlist = p;
                    }
                } else {
//                    nexts.addLast(p);
                    path_t headNextLast = find_last(head.next);
                    p.next = null;
                    if (headNextLast != null) {
                        headNextLast.next = p;
                    } else {
                        head.next = p;
                    }
                }
                //
                p = cur;
            }
            /* clear bm */
            bm.clear_with_bbox(bbox);
            /* now schedule head->childlist and head->next for further
           processing */
            //System.out.println("step1, E1iter="+(iter++));
            //print_path(paths);
//            head.childlist = children.getValue();
//            head.next = nexts.getValue();
            if (head.next != null) {
                head.next.childlist = heap;
                heap = head.next;
            }
            if (head.childlist != null) {
                head.childlist.childlist = heap;
                heap = head.childlist;
            }
            //System.out.println("step1, E2iter="+(iter++));
            //print_path(paths);
        }

        /* copy sibling structure from "next" to "sibling" component */
        p = paths.get(0);
        while (p != null) {
            path_t p1 = p.sibling;
            p.sibling = p.next;
            p = p1;
        }
        /* reconstruct a new linked list ("next") structure from tree
       ("childlist", "sibling") structure. This code is slightly messy,
       because we use a heap to make it tail recursive: the heap
       contains a list of childlists which still need to be
       processed. */
        heap = paths.get(0);
        if (heap != null) {
            heap.next = null;  /* heap is a linked list of childlists */
        }

        path_t plist = null;
        while (heap != null) {
            //System.out.println("step2, iter="+(iter++));
            //print_path(paths);
            path_t heap1 = heap.next;
            for (p = heap; p != null; p = p.sibling) {
                /* p is a positive path */
                /* append to linked list */
                path_t plistLast = find_last(plist);
                p.next = null;
                if (plistLast != null) {
                    plistLast.next = p;
                } else {
                    plist = p;
                }
//                System.out.println("step2-0, iter="+(iter++));
//                print_path(paths);
                /* go through its children */
                for (path_t p1 = p.childlist; p1 != null; p1 = p1.sibling) {
                    /* append to linked list */
                    plistLast = find_last(plist);
                    p1.next = null;
                    if (plistLast != null) {
                        plistLast.next = p1;
                    } else {
                        plist = p1;
                    }
                    /* append its childlist to heap, if non-empty */
                    if (p1.childlist != null) {
                        if (heap1 != null) {
                            find_last(heap1).next = p1.childlist;
                        } else {
                            heap1 = p1.childlist;
                        }
                    }
//                    System.out.println("step2a, iter="+(iter++));
//                    print_path(paths);
                }
            }
            heap = heap1;
        }

//        System.out.println("Final plist");
//        print_path(paths);
        return plist;
    }

    private path_t find_last(path_t lst) {
        path_t prev = null;
        while (lst != null) {
            prev = lst;
            lst = lst.next;
        }
        return prev;
    }

    private bbox_t setbbox_path(path_t p) {
        bbox_t bbox = new bbox_t();

        bbox.y0 = Integer.MAX_VALUE;
        bbox.y1 = 0;
        bbox.x0 = Integer.MAX_VALUE;
        bbox.x1 = 0;

        for (point_t pt : p.priv.pt) {
            int x = pt.x;
            int y = pt.y;

            if (x < bbox.x0) {
                bbox.x0 = x;
            }
            if (x > bbox.x1) {
                bbox.x1 = x;
            }
            if (y < bbox.y0) {
                bbox.y0 = y;
            }
            if (y > bbox.y1) {
                bbox.y1 = y;
            }
        }
        return bbox;
    }

    private void xor_path(Bitmap bm, path_t p) {
        if (p.priv.pt.size() == 0) return; /* a path of length 0 is silly, but legal */
        int y1 = p.priv.pt.get(p.priv.pt.size() - 1).y;
        int xa = p.priv.pt.get(0).x;
        for (point_t pt : p.priv.pt) {
            if (pt.y != y1) {
                /* efficiently invert the rectangle [x,xa] x [y,y1] */
                xor_to_ref(bm, pt.x, Math.min(pt.y, y1), xa);
                y1 = pt.y;
            }
        }

    }

    private void xor_to_ref(Bitmap bm, int x, int y, int xa) {
        if (x > xa) {
            for (int i = xa; i < x; i++) bm.toggle(i, y);
        } else {
            for (int i = x; i < xa; i++) bm.toggle(i, y);
        }
    }

    private path_t findpath(Bitmap bm, int x0, int y0, char sign, param_t.TurnPolicy turnpolicy) {
        int x, y, dirx, diry, size, area, tmp;
        boolean c, d;
        ArrayList<point_t> pt = new ArrayList<point_t>();

        x = x0;
        y = y0;
        dirx = 0;
        diry = -1;
        size = 0;
        area = 0;

        while (true) {
            pt.add(new point_t(x, y));

            /* move to next point */
            x += dirx;
            y += diry;
            area += x * diry;

            /* path complete? */
            if (x == x0 && y == y0) {
                break;
            }

            /* determine next direction */
            c = BM_GET(bm, x + (dirx + diry - 1) / 2, y + (diry - dirx - 1) / 2);
            d = BM_GET(bm, x + (dirx - diry - 1) / 2, y + (diry + dirx - 1) / 2);

            if (c && !d) {               /* ambiguous turn */
                if (turnpolicy == param_t.TurnPolicy.POTRACE_TURNPOLICY_RIGHT
                        || (turnpolicy == param_t.TurnPolicy.POTRACE_TURNPOLICY_BLACK && sign == '+')
                        || (turnpolicy == param_t.TurnPolicy.POTRACE_TURNPOLICY_WHITE && sign == '-')
                        || (turnpolicy == param_t.TurnPolicy.POTRACE_TURNPOLICY_RANDOM && detrand(x, y) != 0)
                        || (turnpolicy == param_t.TurnPolicy.POTRACE_TURNPOLICY_MAJORITY && majority(bm, x, y))
                        || (turnpolicy == param_t.TurnPolicy.POTRACE_TURNPOLICY_MINORITY && !majority(bm, x, y))) {
                    tmp = dirx;              /* right turn */
                    dirx = diry;
                    diry = -tmp;
                } else {
                    tmp = dirx;              /* left turn */
                    dirx = -diry;
                    diry = tmp;
                }
            } else if (c) {              /* right turn */
                tmp = dirx;
                dirx = diry;
                diry = -tmp;
            } else if (!d) {             /* left turn */
                tmp = dirx;
                dirx = -diry;
                diry = tmp;
            }
        } /* while this path */

        /* allocate new path object */
        return new path_t(pt, area, sign);

    }

    int z;
    static byte detrand_t[] = {
            /* non-linear sequence: constant term of inverse in GF(8),
          mod x^8+x^4+x^3+x+1 */
            0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1,
            0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
            1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1,
            0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0,
            0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0,
            0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0,
            0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1,
            1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0,
            0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1,
            1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
    };

    static int detrand(int x, int y) {

        if (true) throw new RuntimeException("Untested: detrand");
        /* 0x04b3e375 and 0x05a8ef93 are chosen to contain every possible
      5-bit sequence */
        int z = ((0x04b3e375 * x) ^ y) * 0x05a8ef93;
        z = detrand_t[z & 0xff] ^ detrand_t[(z >> 8) & 0xff] ^ detrand_t[(z >> 16) & 0xff] ^ detrand_t[(z >> 24) & 0xff];
        return z;
    }

    static boolean majority(Bitmap bm, int x, int y) {
        int i, a, ct;

        for (i = 2; i < 5; i++) { /* check at "radius" i */
            ct = 0;
            for (a = -i + 1; a <= i - 1; a++) {
                ct += BM_GET(bm, x + a, y + i - 1) ? 1 : -1;
                ct += BM_GET(bm, x + i - 1, y + a - 1) ? 1 : -1;
                ct += BM_GET(bm, x + a - 1, y - i) ? 1 : -1;
                ct += BM_GET(bm, x - i, y + a) ? 1 : -1;
            }
            if (ct > 0) {
                return true;
            } else if (ct < 0) {
                return false;
            }
        }
        return false;
    }


    private int findnext(Bitmap bm, point_t p) {
        int x0 = p.x;
        for (int y = p.y; y >= 0; y--) {
            for (int x = x0; x < bm.w; x++) {
                if (bm.get(x, y)) {
                    p.x = x;
                    p.y = y;
                    return 0;
                }
            }
            x0 = 0;
        }
        return 1;
    }

    static final int mod(int a, int n) {
        return a >= n ? a % n : a >= 0 ? a : n - 1 - (-1 - a) % n;
    }

    static int sign(int x) {
        return ((x) > 0 ? 1 : (x) < 0 ? -1 : 0);
    }

    static final int sign(double x) {
        return ((x) > 0 ? 1 : (x) < 0 ? -1 : 0);
    }

    static final int abs(int x) {
        return x > 0 ? x : -x;
    }

    /* calculate p1 x p2 */
    static final int xprod(point_t p1, point_t p2) {
        return p1.x * p2.y - p1.y * p2.x;
    }

    /* calculate (p1-p0)x(p3-p2) */
    static final double cprod(dpoint_t p0, dpoint_t p1, dpoint_t p2, dpoint_t p3) {
        double x1, y1, x2, y2;

        x1 = p1.x - p0.x;
        y1 = p1.y - p0.y;
        x2 = p3.x - p2.x;
        y2 = p3.y - p2.y;

        return x1 * y2 - x2 * y1;
    }

    /* calculate (p1-p0)*(p2-p0) */
    static double iprod(dpoint_t p0, dpoint_t p1, dpoint_t p2) {
        double x1, y1, x2, y2;

        x1 = p1.x - p0.x;
        y1 = p1.y - p0.y;
        x2 = p2.x - p0.x;
        y2 = p2.y - p0.y;

        return x1 * x2 + y1 * y2;
    }

    /* calculate (p1-p0)*(p3-p2) */
    static double iprod1(dpoint_t p0, dpoint_t p1, dpoint_t p2, dpoint_t p3) {
        double x1, y1, x2, y2;

        x1 = p1.x - p0.x;
        y1 = p1.y - p0.y;
        x2 = p3.x - p2.x;
        y2 = p3.y - p2.y;

        return x1 * x2 + y1 * y2;
    }

    /* calculate distance between two points */
    static double ddist(dpoint_t p, dpoint_t q) {
        return Math.sqrt(sq(p.x - q.x) + sq(p.y - q.y));
    }

    static double sq(double d) {
        return d * d;
    }


    /* calculate point of a bezier curve */
    static dpoint_t bezier(double t, dpoint_t p0, dpoint_t p1, dpoint_t p2, dpoint_t p3) {
        double s = 1 - t;
        dpoint_t res = new dpoint_t();

        /* Note: a good optimizing compiler (such as gcc-3) reduces the
       following to 16 multiplications, using common subexpression
       elimination. */

        res.x = s * s * s * p0.x + 3 * (s * s * t) * p1.x + 3 * (t * t * s) * p2.x + t * t * t * p3.x;
        res.y = s * s * s * p0.y + 3 * (s * s * t) * p1.y + 3 * (t * t * s) * p2.y + t * t * t * p3.y;

        return res;
    }

    /* calculate the point t in [0..1] on the (convex) bezier curve
       (p0,p1,p2,p3) which is tangent to q1-q0. Return -1.0 if there is no
       solution in [0..1]. */
    static double tangent(dpoint_t p0, dpoint_t p1, dpoint_t p2, dpoint_t p3, dpoint_t q0, dpoint_t q1) {
        double A, B, C;   /* (1-t)^2 A + 2(1-t)t B + t^2 C = 0 */
        double a, b, c;   /* a t^2 + b t + c = 0 */
        double d, s, r1, r2;

        A = cprod(p0, p1, q0, q1);
        B = cprod(p1, p2, q0, q1);
        C = cprod(p2, p3, q0, q1);

        a = A - 2 * B + C;
        b = -2 * A + 2 * B;
        c = A;

        d = b * b - 4 * a * c;

        if (a == 0 || d < 0) {
            return -1.0;
        }

        s = Math.sqrt(d);

        r1 = (-b + s) / (2 * a);
        r2 = (-b - s) / (2 * a);

        if (r1 >= 0 && r1 <= 1) {
            return r1;
        } else if (r2 >= 0 && r2 <= 1) {
            return r2;
        } else {
            return -1.0;
        }
    }

    /* return 1 if a <= b < c < a, in a cyclic sense (mod n) */
    static final boolean cyclic(int a, int b, int c) {
        if (a <= c) {
            return (a <= b && b < c);
        } else {
            return (a <= b || b < c);
        }
    }

    static final int floordiv(int a, int n) {
        return a >= 0 ? a / n : -1 - (-1 - a) / n;
    }

    static double penalty3(privpath_t pp, int i, int j) {
        int n = pp.pt.size();
        point_t pt[] = pp.getPTArray();
        sums_t[] sums = pp.sums;

        /* assume 0<=i<j<=n  */
        double x, y, x2, xy, y2;
        double k;
        double a, b, c, s;
        double px, py, ex, ey;

        int r = 0; /* rotations from i to j */

        if (j >= n) {
            j -= n;
            r = 1;
        }

        /* critical inner loop: the "if" gives a 4.6 percent speedup */
        if (r == 0) {
            x = sums[j + 1].x - sums[i].x;
            y = sums[j + 1].y - sums[i].y;
            x2 = sums[j + 1].x2 - sums[i].x2;
            xy = sums[j + 1].xy - sums[i].xy;
            y2 = sums[j + 1].y2 - sums[i].y2;
            k = j + 1 - i;
        } else {
            x = sums[j + 1].x - sums[i].x + sums[n].x;
            y = sums[j + 1].y - sums[i].y + sums[n].y;
            x2 = sums[j + 1].x2 - sums[i].x2 + sums[n].x2;
            xy = sums[j + 1].xy - sums[i].xy + sums[n].xy;
            y2 = sums[j + 1].y2 - sums[i].y2 + sums[n].y2;
            k = j + 1 - i + n;
        }

        px = (pt[i].x + pt[j].x) / 2.0 - pt[0].x;
        py = (pt[i].y + pt[j].y) / 2.0 - pt[0].y;
        ey = (pt[j].x - pt[i].x);
        ex = -(pt[j].y - pt[i].y);

        a = ((x2 - 2 * x * px) / k + px * px);
        b = ((xy - x * py - y * px) / k + px * py);
        c = ((y2 - 2 * y * py) / k + py * py);

        s = ex * ex * a + 2 * ex * ey * b + ey * ey * c;

        return Math.sqrt(s);
    }

    static double quadform(double[][] Q, dpoint_t w) {
        double v[] = new double[3];
        int i, j;
        double sum;

        v[0] = w.x;
        v[1] = w.y;
        v[2] = 1;
        sum = 0.0;

        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                sum += v[i] * Q[i][j] * v[j];
            }
        }
        return sum;
    }


    static void pointslope(privpath_t pp, int i, int j, dpoint_t ctr, dpoint_t dir) {
        /* assume i<j */

        int n = pp.pt.size();
        sums_t[] sums = pp.sums;

        double x, y, x2, xy, y2;
        double k;
        double a, b, c, lambda2, l;
        int r = 0; /* rotations from i to j */

        while (j >= n) {
            j -= n;
            r += 1;
        }
        while (i >= n) {
            i -= n;
            r -= 1;
        }
        while (j < 0) {
            j += n;
            r -= 1;
        }
        while (i < 0) {
            i += n;
            r += 1;
        }

        x = sums[j + 1].x - sums[i].x + r * sums[n].x;
        y = sums[j + 1].y - sums[i].y + r * sums[n].y;
        x2 = sums[j + 1].x2 - sums[i].x2 + r * sums[n].x2;
        xy = sums[j + 1].xy - sums[i].xy + r * sums[n].xy;
        y2 = sums[j + 1].y2 - sums[i].y2 + r * sums[n].y2;
        k = j + 1 - i + r * n;

        ctr.x = x / k;
        ctr.y = y / k;

        a = (x2 - (double) x * x / k) / k;
        b = (xy - (double) x * y / k) / k;
        c = (y2 - (double) y * y / k) / k;

        lambda2 = (a + c + Math.sqrt((a - c) * (a - c) + 4 * b * b)) / 2; /* larger e.value */

        /* now find e.vector for lambda2 */
        a -= lambda2;
        c -= lambda2;

        if (Math.abs(a) >= Math.abs(c)) {
            l = Math.sqrt(a * a + b * b);
            if (l != 0) {
                dir.x = -b / l;
                dir.y = a / l;
            }
        } else {
            l = Math.sqrt(c * c + b * b);
            if (l != 0) {
                dir.x = -c / l;
                dir.y = b / l;
            }
        }
        if (l == 0) {
            dir.x = dir.y = 0;   /* sometimes this can happen when k=4:
			      the two eigenvalues coincide */
        }
    }

    static dpoint_t interval(double lambda, dpoint_t a, dpoint_t b) {
        dpoint_t res = new dpoint_t();

        res.x = a.x + lambda * (b.x - a.x);
        res.y = a.y + lambda * (b.y - a.y);
        return res;
    }

    static double ddenom(dpoint_t p0, dpoint_t p2) {
        point_t r = dorth_infty(p0, p2);

        return r.y * (p2.x - p0.x) - r.x * (p2.y - p0.y);
    }

    static point_t dorth_infty(dpoint_t p0, dpoint_t p2) {
        point_t r = new point_t();

        r.y = (int) Math.signum(p2.x - p0.x);
        r.x = (int) -Math.signum(p2.y - p0.y);

        return r;
    }

    static double dpara(dpoint_t p0, dpoint_t p1, dpoint_t p2) {
        double x1, y1, x2, y2;

        x1 = p1.x - p0.x;
        y1 = p1.y - p0.y;
        x2 = p2.x - p0.x;
        y2 = p2.y - p0.y;

        return x1 * y2 - x2 * y1;
    }

    int opti_penalty(privpath_t pp, int i, int j, opti_t res, double opttolerance, int[] convc, double[] areac) {
        int m = pp.curve.n;
        int k, k1, k2, conv, i1;
        double area, alpha, d, d1, d2;
        dpoint_t p0, p1, p2, p3, pt;
        double A, R, A1, A2, A3, A4;
        double s, t;

        /* check convexity, corner-freeness, and maximum bend < 179 degrees */

        if (i == j) {  /* sanity - a full loop can never be an opticurve */
            return 1;
        }

        k = i;
        i1 = mod(i + 1, m);
        k1 = mod(k + 1, m);
        conv = convc[k1];
        if (conv == 0) {
            return 1;
        }
        d = ddist(pp.curve.vertex[i], pp.curve.vertex[i1]);
        for (k = k1; k != j; k = k1) {
            k1 = mod(k + 1, m);
            k2 = mod(k + 2, m);
            if (convc[k1] != conv) {
                return 1;
            }
            if (sign(cprod(pp.curve.vertex[i], pp.curve.vertex[i1], pp.curve.vertex[k1], pp.curve.vertex[k2])) != conv) {
                return 1;
            }
            if (iprod1(pp.curve.vertex[i], pp.curve.vertex[i1], pp.curve.vertex[k1], pp.curve.vertex[k2]) < d * ddist(pp.curve.vertex[k1], pp.curve.vertex[k2]) * COS179) {
                return 1;
            }
        }

        /* the curve we're working in: */
        p0 = pp.curve.c[mod(i, m)][2];
        p1 = pp.curve.vertex[mod(i + 1, m)];
        p2 = pp.curve.vertex[mod(j, m)];
        p3 = pp.curve.c[mod(j, m)][2];

        /* determine its area */
        area = areac[j] - areac[i];
        area -= dpara(pp.curve.vertex[0], pp.curve.c[i][2], pp.curve.c[j][2]) / 2;
        if (i >= j) {
            area += areac[m];
        }

        /* find intersection o of p0p1 and p2p3. Let t,s such that o =
      interval(t,p0,p1) = interval(s,p3,p2). Let A be the area of the
      triangle (p0,o,p3). */

        A1 = dpara(p0, p1, p2);
        A2 = dpara(p0, p1, p3);
        A3 = dpara(p0, p2, p3);
        /* A4 = dpara(p1, p2, p3); */
        A4 = A1 + A3 - A2;

        if (A2 == A1) {  /* this should never happen */
            return 1;
        }

        t = A3 / (A3 - A4);
        s = A2 / (A2 - A1);
        A = A2 * t / 2.0;

        if (A == 0.0) {  /* this should never happen */
            return 1;
        }

        R = area / A;     /* relative area */
        alpha = 2 - Math.sqrt(4 - R / 0.3);  /* overall alpha for p0-o-p3 curve */

        res.c[0] = interval(t * alpha, p0, p1);
        res.c[1] = interval(s * alpha, p3, p2);
        res.alpha = alpha;
        res.t = t;
        res.s = s;

        p1 = res.c[0];
        p2 = res.c[1];  /* the proposed curve is now (p0,p1,p2,p3) */

        res.pen = 0;

        /* calculate penalty */
        /* check tangency with edges */
        for (k = mod(i + 1, m); k != j; k = k1) {
            k1 = mod(k + 1, m);
            t = tangent(p0, p1, p2, p3, pp.curve.vertex[k], pp.curve.vertex[k1]);
            if (t < -.5) {
                return 1;
            }
            pt = bezier(t, p0, p1, p2, p3);
            d = ddist(pp.curve.vertex[k], pp.curve.vertex[k1]);
            if (d == 0.0) {  /* this should never happen */
                return 1;
            }
            d1 = dpara(pp.curve.vertex[k], pp.curve.vertex[k1], pt) / d;
            if (Math.abs(d1) > opttolerance) {
                return 1;
            }
            if (iprod(pp.curve.vertex[k], pp.curve.vertex[k1], pt) < 0 || iprod(pp.curve.vertex[k1], pp.curve.vertex[k], pt) < 0) {
                return 1;
            }
            res.pen += sq(d1);
        }

        /* check corners */
        for (k = i; k != j; k = k1) {
            k1 = mod(k + 1, m);
            t = tangent(p0, p1, p2, p3, pp.curve.c[k][2], pp.curve.c[k1][2]);
            if (t < -.5) {
                return 1;
            }
            pt = bezier(t, p0, p1, p2, p3);
            d = ddist(pp.curve.c[k][2], pp.curve.c[k1][2]);
            if (d == 0.0) {  /* this should never happen */
                return 1;
            }
            d1 = dpara(pp.curve.c[k][2], pp.curve.c[k1][2], pt) / d;
            d2 = dpara(pp.curve.c[k][2], pp.curve.c[k1][2], pp.curve.vertex[k1]) / d;
            d2 *= 0.75 * pp.curve.alpha[k1];
            if (d2 < 0) {
                d1 = -d1;
                d2 = -d2;
            }
            if (d1 < d2 - opttolerance) {
                return 1;
            }
            if (d1 < d2) {
                res.pen += sq(d1 - d2);
            }
        }

        return 0;
    }


}
