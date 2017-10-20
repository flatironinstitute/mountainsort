#include "p_compute_amplitudes.h"
#include <QFile>
#include <diskreadmda32.h>
#include <mda.h>
#include <mda32.h>
#include "mlutil.h"

namespace P_compute_amplitudes {
}

bool p_compute_amplitudes(QString timeseries, QString event_times, QString amplitudes_out, P_compute_amplitudes_opts opts)
{
    DiskReadMda32 X(timeseries);
    Mda ET(event_times);
    bigint L = ET.totalSize();
    Mda32 A(1, L);
    for (bigint i = 0; i < L; i++) {
        double amp = 0;
        bigint time0 = ET.value(i);
        if (opts.central_channel) {
            amp = X.value(opts.central_channel - 1, time0);
        }
        else {
            double maxval = 0;
            for (bigint m = 0; m < X.N1(); m++) {
                double val = X.value(m, time0);
                if (qAbs(val) > qAbs(maxval))
                    maxval = val;
            }
            amp = maxval;
        }
        A.setValue(amp, i);
    }
    return A.write32(amplitudes_out);
}

namespace P_consolidate_clusters {
}
