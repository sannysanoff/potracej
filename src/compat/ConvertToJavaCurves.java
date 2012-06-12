package compat;

import potracej.curve_t;
import potracej.dpoint_t;
import potracej.path_t;

import java.awt.geom.Point2D;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Created with IntelliJ IDEA.
 * User: san
 * Date: 6/12/12
 * Time: 12:28 PM
 * To change this template use File | Settings | File Templates.
 */
public class ConvertToJavaCurves {

    public static String f(double d) {
        DecimalFormat decimalFormat = new DecimalFormat();
        decimalFormat.setMaximumFractionDigits(1);
        return decimalFormat.format(d);
    }

    public static int convert(path_t plist, HashSet<Point2D.Double> points, ArrayList<PathElement> result) {
        int nodeCount = 0;
        int i;
        path_t child;

        path_t node;
        for (node = plist; node != null; node = node.next) {
            curve_t curve = node.curve;
            //g_message("node->fm:%d\n", node->fm);
            if (curve.n == 0)
                continue;
            dpoint_t[] pt = curve.c[curve.n - 1];
            double x0 = 0.0;
            double y0 = 0.0;
            double x1 = 0.0;
            double y1 = 0.0;
            double x2 = pt[2].x;
            double y2 = pt[2].y;
            //Have we been here already?
            Point2D.Double point = new Point2D.Double(x2, y2);
            if (points.contains(point)) {
                continue;
            } else {
                points.add(point);
            }

            result.add(new PathElement(PathElement.Type.MOVE_TO, x2, y2, 0, 0, 0, 0));
            nodeCount++;

            for (i = 0; i < curve.n; i++) {
                pt = curve.c[i];
                x0 = pt[0].x;
                y0 = pt[0].y;
                x1 = pt[1].x;
                y1 = pt[1].y;
                x2 = pt[2].x;
                y2 = pt[2].y;
                switch (curve.tag[i]) {
                    case POTRACE_CORNER:
                        result.add(new PathElement(PathElement.Type.LINE_TO, x1, y1, 0, 0, 0, 0));
                        result.add(new PathElement(PathElement.Type.LINE_TO, x2, y2, 0, 0, 0, 0));
                        break;
                    case POTRACE_CURVETO:
                        result.add(new PathElement(PathElement.Type.CURVE_TO, x0, y0, x1, y1, x2, y2));
                        break;
                    default:
                        break;
                }
                nodeCount++;
            }
            result.add(new PathElement(PathElement.Type.CLOSE_PATH, 0, 0, 0, 0, 0, 0));

            for (child = node.childlist; child != null; child = child.sibling) {
                result.add(new PathElement(PathElement.Type.PUSH_PARENT, 0, 0, 0, 0, 0, 0));
                nodeCount += convert(child, points, result);
                result.add(new PathElement(PathElement.Type.PUSH_PARENT, 0, 0, 0, 0, 0, 0));
            }
        }

        return nodeCount;

    }

}
