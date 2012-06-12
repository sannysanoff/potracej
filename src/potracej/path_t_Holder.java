package potracej;

import java.util.ArrayList;

public class path_t_Holder {
    path_t value;
    path_t last;

    public path_t_Holder(path_t value) {
        this.value = value;
    }

    public void addLast(path_t path) {
        if (value == null) {
            value = last = path;
        } else {
            last.next = path;
            last = path;
        }
        path.next = null;
    }

    public path_t getValue() {
        return value;
    }

    public path_t getLast() {
        return last;
    }
}
