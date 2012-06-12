package potracej;

import java.util.ArrayList;
import java.util.Arrays;

public class privpath_t {
    public ArrayList<point_t> pt;     /* pt[len]: path as extracted from bitmap */
    public int lon[];        /* lon[len]: (i,lon[i]) = longest straight line from i */

    public int x0, y0;      /* origin for sums */
    public sums_t sums[];    /* sums[len+1]: cache for fast summing */

    public int m;           /* length of optimal polygon */
    public int po[];         /* po[m]: optimal polygon */

    public privcurve_t curve;   /* curve[m]: array of curve elements */
    public privcurve_t ocurve;  /* ocurve[om]: array of curve elements */
    public privcurve_t fcurve;  /* final curve: this points to either curve or
    		                    ocurve. Do not free this separately. */
    private point_t[] PTArray;

    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append("privpath_t{");
        sb.append("lon=" + intArrayToString(lon));
        sb.append(" y0=" + x0);
        sb.append(" y0=" + x0);
        sb.append(" m=" + m);
        sb.append(" po=" + intArrayToString(po));
        sb.append(" sums=" + Arrays.asList(sums));
        sb.append(" fcurve=" + fcurve);
        sb.append(" ocurve=" + ocurve);
        sb.append(" fcurve=" + fcurve);
        sb.append("}/*privpath_t*/");
        return sb.toString();
    }

    public static String intArrayToString(int[] arr) {
        StringBuffer sb = new StringBuffer();
        sb.append("{");
        for (int i = 0; i < arr.length; i++) {
            sb.append(arr[i]);
            sb.append(",");
        }
        sb.append("}");
        return sb.toString();
    }

    public point_t[] getPTArray() {
        if (PTArray == null) {
            PTArray = pt.toArray(new point_t[pt.size()]);
        }
        return PTArray;
    }
}
