package compat;

/**
* Created with IntelliJ IDEA.
* User: san
* Date: 6/12/12
* Time: 12:29 PM
* To change this template use File | Settings | File Templates.
*/
public class PathElement {

    public enum Type {
        MOVE_TO,
        LINE_TO,
        CURVE_TO,
        PUSH_PARENT,
        POP_PARENT,
        CLOSE_PATH
    }

    Type type;
    double p0x, p0y;
    double p1x, p1y;
    double p2x, p2y;

    public PathElement(Type type, double p0x, double p0y) {
        this(type, p0x, p0y, 0, 0, 0, 0);
    }

    public PathElement(Type type, double p0x, double p0y, double p1x, double p1y, double p2x, double p2y) {
        this.type = type;
        this.p0x = p0x;
        this.p0y = p0y;
        this.p1x = p1x;
        this.p1y = p1y;
        this.p2x = p2x;
        this.p2y = p2y;
    }

    public double getP0x() {
        return p0x;
    }

    public double getP0y() {
        return p0y;
    }

    public double getP1x() {
        return p1x;
    }

    public double getP1y() {
        return p1y;
    }

    public double getP2x() {
        return p2x;
    }

    public double getP2y() {
        return p2y;
    }

    public Type getType() {
        return type;
    }

    @Override
    public String toString() {
        switch (type) {
            case CLOSE_PATH: return "X";
            case PUSH_PARENT: return "P";
            case POP_PARENT: return "p";
            case MOVE_TO: return "M,"+p0x+","+p0y;
            case LINE_TO: return "L,"+p0x+","+p0y;
            case CURVE_TO: return "C,"+p0x+","+p0y+","+p1x+","+p1y+","+p2x+","+p2y;
            default:
                throw new RuntimeException("Sorry.");
        }
    }
}
