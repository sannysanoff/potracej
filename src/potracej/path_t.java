package potracej;

import java.util.ArrayList;

public class path_t {
    public int area;                         /* area of the bitmap path */
    public char sign;                         /* '+' or '-', depending on orientation */
    public curve_t curve;            /* this path's vector data */

    public path_t next;      /* linked list structure */

    public path_t childlist; /* tree structure */
    public path_t sibling;   /* tree structure */

    privpath_t priv = new privpath_t();  /* private state */

    public path_t(ArrayList<point_t> pt, int area, char sign) {
        this.priv.pt = pt;
        this.area = area;
        this.sign = sign;
    }

    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append("path{("+priv.pt.size()+" points):");
        for (point_t point_t : priv.pt) {
            sb.append(point_t.toString()+" ");
        }
        sb.append(", priv="+priv+"}");
        //sb.append(", \nnext="+next+"}");
        return sb.toString();
    }

}
