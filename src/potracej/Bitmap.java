package potracej;

/**
 * Grayscale bitmap
 */
public class Bitmap {
    byte[] arr; // level of gray
    int w, h;

    public Bitmap(int w, int h) {
        arr = new byte[w * h];
        this.w = w;
        this.h = h;
    }

    public Bitmap(byte[] arr, int w, int h) {
        this.arr = arr;
        this.w = w;
        this.h = h;
    }

    public void put(int x, int y, byte b) {
        arr[x + w * y] = b;
    }

    public void put(int x, int y, int i) {
        arr[x + w * y] = (byte) (i & 0xFF);
    }

    public boolean get(int x, int y) {
        if (x < 0 || y < 0 || x >= w || y >= h) {
            return false;
        }
        return arr[x + w * y] != 0;
    }

    public Bitmap dup() {
        Bitmap bm = new Bitmap(w, h);
        System.arraycopy(arr, 0, bm.arr, 0, arr.length);
        return bm;
    }

    public void clearexcess() {
        // not needed
    }

    public static boolean BM_GET(Bitmap bm, int x, int y) {
        return bm.get(x, y);
    }

    public void toggle(int x, int y) {
        put(x, y, get(x, y) ? 0 : 1);
    }

    public void clear() {
        for (int i = 0; i < arr.length; i++) {
            arr[i] = 0;
        }
    }

    public void clear_with_bbox(bbox_t bbox) {
        for (int y = bbox.y0; y < bbox.y1; y++) {
            for (int i = bbox.x0; i < bbox.x1; i++) {
                put(i, y, 0);
            }
        }
    }

    public String toDebugString() {
        StringBuffer sb = new StringBuffer();
        sb.append("====================================================================================================================================================================================\n");
        for(int y=0; y<h; y++) {
            for(int x=0; x<w; x++) {
                if (get(x, y))
                    sb.append("#");
                else sb.append(".");
            }
            sb.append("\n");
        }
        return sb.toString();
    }

    public int getWidth() {
        return w;
    }

    public int getHeight() {
        return h;
    }
}
