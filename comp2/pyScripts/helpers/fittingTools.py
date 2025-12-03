
class GlobalChi2(object):
    def __init__(self, f1, f2, f3):
        self._f1 = f1
        self._f2 = f2
        self._f3 = f3

    def __call__(self, par):
        # parameter vector is first background (in common 1 and 2) and then is
        # signal (only in 2)

        # the zero-copy way to get a numpy array from a double *
        par_arr = np.frombuffer(par, dtype=np.float64, count=-1)

        p1 = par_arr[ipar1]
        p2 = par_arr[ipar2]
        p3 = par_arr[ipar3]

        return self._f1(p1) + self._f2(p2) + self._f3(p3)