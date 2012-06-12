package potracej;

public class point_t {
    public int x;
    public int y;



    public point_t(point_t o) {
        this.x = o.x;
        this.y = o.y;
    }


    public point_t() {
    }
    public point_t(int x, int y) {
        this.x = x;
        this.y = y;
    }

    public String toString() {
        return "{"+x+","+y+"}";
    }
}
