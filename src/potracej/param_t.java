package potracej;

public class param_t {

    enum TurnPolicy {
        POTRACE_TURNPOLICY_RIGHT,
        POTRACE_TURNPOLICY_BLACK,
        POTRACE_TURNPOLICY_WHITE,
        POTRACE_TURNPOLICY_RANDOM,
        POTRACE_TURNPOLICY_MAJORITY,
        POTRACE_TURNPOLICY_MINORITY

    }

    public TurnPolicy turnPolicy = TurnPolicy.POTRACE_TURNPOLICY_MINORITY;   /* resolves ambiguous turns in path decomposition */
    public int turdsize = 2;            /* area of largest path to be ignored */
    public double alphamax = 1;      /* corner threshold */
    public int opticurve = 1;       /* use curve optimization? */
    public double opttolerance = 0.2; /* curve optimization tolerance */
}
